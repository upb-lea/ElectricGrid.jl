using DrWatson
@quickactivate "MicroGridSimWithRL"

using PyCall
@pyinclude(srcdir("nodeconstructor.py"))

using BenchmarkTools
using Statistics
using DifferentialEquations
using LSODA
using ControlSystems
using JSON

include(srcdir("env.jl"))

global yout

function execute_env(env::SimEnv, actions::Matrix{Float64})
    RLBase.reset!(env)
    for i = 1:length(actions[1,:])
        env(actions[:,i])
    end
end

function timing_experiment_simulation(repeat::Int64=5, loops::Int64=10, num_nodes::Vector{Int64}=[5,10,25,50,100],
                                        t_end::Vector{Float64}=[0.001,0.01,0.1], ts::Float64=1e-4,
                                        methods::Vector{String}=["control", "env_without_agent", "lsoda"],
                                        parameter::Dict=Dict())
    
    if isempty(parameter)
        parameter["R_source"] = 0.4
        parameter["L_source"] = 2.3e-3
        parameter["C_source"] = 10e-6
        parameter["L_cabel"] = 2.3e-3
        parameter["R_cabel"] = 0.4
        parameter["R_load"] = 14
        parameter["V_dc"] = 300
    end

    limits = Dict("i_lim" => 20, "v_lim" => 600)
    ref = 200

    t_result_mean = zeros(length(methods), length(num_nodes), length(t_end))
    t_result_std = zeros(length(methods), length(num_nodes), length(t_end))
    t_num_samples = zeros(Int64, length(methods), length(num_nodes), length(t_end))

    timer = nothing

    for n = 1:length(methods)
        for k = 1:length(num_nodes)
            nc = py"NodeConstructor"(num_nodes[k], num_nodes[k], parameter)
            A, B, C, D = nc.get_sys()

            #println("Size of x: $(length(A[1,:]))")
            #println("Size of u: $(length(B[1,:]))")

            if methods[n] == "env_without_agent"
                norm_array = vcat( [limits[i] for j = 1:nc.num_source for i in ["i_lim", "v_lim"]], [limits["i_lim"] for i = 1:nc.num_connections] )
                global env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, v_dc=parameter["V_dc"], ts=rationalize(ts))
            elseif methods[n] == "lsoda"
                function f!(du, u, p, t)
                    du[:] = A * u + B * p
                end
            elseif methods[n] == "control"
                Ad = exp(A*ts)
                Bd = A \ (Ad - C) * B
                global sys_d = ss(Ad, Bd, C, D, ts)
            end

            for l = 1:length(t_end)
                timer = nothing
                t = collect(0:ts:t_end[l])
                
                if methods[n] == "env_without_agent"
                    RLBase.reset!(env)
                    #global actions = rand(Float64, ( length(B[1,:]), length(t)-1 )) .*2 .-1
                    global u = [0.7 for i = 1:length(t)]
                    global uu = [u for i = 1:length(B[1,:]) ]
                    global actions = mapreduce(permutedims, vcat, uu)
                    timer = @benchmark execute_env(env, actions) samples = repeat evals = loops seconds = 1000
                elseif methods[n] == "lsoda"
                    p = [230.0 for i = 1:length(B[1,:])]
                    tspan = (0.0,t_end[l])
                    u0 = [0.0 for i = 1:length(A[1,:])]
                    global problem = ODEProblem(f!,u0,tspan,p)
                    timer = @benchmark solve(problem,lsoda(), reltol=1e-6, abstol=1e-6) samples = repeat evals = loops seconds = 1000
                elseif methods[n] == "control"
                    global x0 = [0.0 for i = 1:length(A[1,:])]
                    #global u = rand(Float64, ( length(t) )) .*2 .-1
                    global u = [230.0 for i = 1:length(t)]
                    global uu = [u for i = 1:length(B[1,:]) ]
                    global uuu = mapreduce(permutedims, vcat, uu)
                    global ttt = t
                    #yout = lsim(sys,uuu,ttt,x0=x0)
                    timer = @benchmark lsim(sys_d,uuu,ttt,x0=x0) samples = repeat evals = loops seconds = 1000
                end

                if timer !== nothing
                    t_result_mean[n,k,l] = mean(timer).time / 1_000_000_000
                    t_result_std[n,k,l] = std(timer).time / 1_000_000_000
                    t_num_samples[n,k,l] = length(timer.times)
                end

                println("done for $(methods[n]) with $(num_nodes[k]) nodes and $(t_end[l]) steps")
            end
        end
    end
    result_dict = Dict()
    result_dict["methods"] = methods
    result_dict["times_mean"] = t_result_mean
    result_dict["times_std"] = t_result_std
    result_dict["num_samples"] = t_num_samples
    result_dict["t_end"] = t_end
    result_dict["num_grid_nodes"] = num_nodes
    result_dict["info"] = "Logs the mean and std of the execution time for all defined methods to simulate the power " *
                            "grid for different simulation times (t_end) and grid size (num_grid_nodes). num_grid_nodes " *
                            "thereby defines the number of sources and the number of loads " *
                            "(grid size = 2*num_grid_nodes). " *
                            "Each experiment is executed loops*repeats-times " *
                            "while the mean and std is calculated based on repeats"
    
    result_dict
end

r = timing_experiment_simulation(5, 5, [5,10,20], [0.001,0.01,0.03])

open(datadir("juliaresults.json"),"w") do f 
    write(f, JSON.json(r)) 
end