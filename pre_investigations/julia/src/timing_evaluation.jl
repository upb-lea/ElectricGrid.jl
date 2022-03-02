using PyCall
@pyinclude(srcdir("nodeconstructor.py"))

using BenchmarkTools
using Statistics
using DifferentialEquations
using LSODA
using ControlSystems
using JSON
using Plots
using Flux
using StableRNGs

include(srcdir("env.jl"))

global yout

global rngg = StableRNG(123)
global initt = glorot_uniform(rngg)

global create_actor(na, ns) = Chain(
    Dense(ns, 400, relu; init = initt),
    Dense(400, 300, relu; init = initt),
    Dense(300, na, tanh; init = initt),
) |> gpu

global create_critic(na, ns) = Chain(
    Dense(ns + na, 400, relu; init = initt),
    Dense(400, 300, relu; init = initt),
    Dense(300, 1; init = initt),
) |> gpu

function create_agent(na, ns)
    Agent(
        policy = DDPGPolicy(
            behavior_actor = NeuralNetworkApproximator(
                model = create_actor(na, ns),
                optimizer = ADAM(),
            ),
            behavior_critic = NeuralNetworkApproximator(
                model = create_critic(na, ns),
                optimizer = ADAM(),
            ),
            target_actor = NeuralNetworkApproximator(
                model = create_actor(na, ns),
                optimizer = ADAM(),
            ),
            target_critic = NeuralNetworkApproximator(
                model = create_critic(na, ns),
                optimizer = ADAM(),
            ),
            γ = 0.99f0,
            ρ = 0.995f0,
            na = na,
            batch_size = 64,
            start_steps = 0,
            start_policy = RandomPolicy(-1.0..1.0; rng = rngg),
            update_after = 1000,
            update_freq = 1,
            act_limit = 1.0,
            act_noise = 0.1,
            rng = rngg,
        ),
        trajectory = CircularArraySARTTrajectory(
            capacity = 100,
            state = Vector{Float32} => (ns,),
            action = Float32 => (na, ),
        ),
    )
end

function execute_env(env::SimEnv, actions::Matrix{Float64}, debug::Bool)
    output = Vector{Float64}()

    RLBase.reset!(env)
    for i = 1:length(actions[1,:])
        env(actions[:,i])
        if debug append!(output, env.state[2]*600) end
    end

    return output
end

function execute_env(env::SimEnv, agent::Agent, t_len::Int, debug::Bool)
    output = Vector{Float64}()

    RLBase.reset!(env)
    for i = 1:t_len
        action = agent(env)
        env(action)
        if debug append!(output, env.state[2]*600) end
    end

    return output
end

function timing_experiment_simulation(repeat::Int64=5, loops::Int64=10, num_nodes::Vector{Int64}=[5,10,25,50,100],
                                        t_end::Vector{Float64}=[0.001,0.01,0.1], num_cm::Int64=1, ts::Float64=1e-4,
                                        methods::Vector{String}=["control", "env_without_agent", "env_with_agent", "lsoda"],
                                        parameter::Dict=Dict(), debug::Bool=false)
    
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

    t_result_mean = zeros(length(t_end), length(num_nodes), length(methods))
    t_result_std = zeros(length(t_end), length(num_nodes), length(methods))
    t_num_samples = zeros(Int64, length(t_end), length(num_nodes), length(methods))


    for n = 1:length(methods)
        for k = 1:length(num_nodes)

            CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(num_nodes[k]) * ".json"))

            for l = 1:length(t_end)
                timer = nothing
                t = collect(0:ts:t_end[l])
                time_list = []

                for c = 1:num_cm
                    CM = reduce(hcat, CM_list[c])'
                    CM = convert(Matrix{Int64}, CM)

                    #nc = py"NodeConstructor"(num_nodes[k], num_nodes[k], parameter)
                    nc = py"NodeConstructor"(num_nodes[k], num_nodes[k], parameter, CM=CM)
                    A, B, C, D = nc.get_sys()
                    ns = length(A[1,:])
                    na = length(B[1,:])

                    #println("Size of x: $(ns)")
                    #println("Size of u: $(na)")

                    if debug
                        result = nothing
                    end

                    if methods[n] == "env_without_agent"
                        norm_array = vcat( [limits[i] for j = 1:nc.num_source for i in ["i_lim", "v_lim"]], [limits["i_lim"] for i = 1:nc.num_connections] )
                        global env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, v_dc=parameter["V_dc"], ts=rationalize(ts))
                    elseif methods[n] == "env_with_agent"
                        norm_array = vcat( [limits[i] for j = 1:nc.num_source for i in ["i_lim", "v_lim"]], [limits["i_lim"] for i = 1:nc.num_connections] )
                        global env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, v_dc=parameter["V_dc"], ts=rationalize(ts))
                        global agent = create_agent(na, ns)
                    elseif methods[n] == "lsoda"
                        function f!(du, u, p, t)
                            du[:] = A * u + B * p
                        end
                    elseif methods[n] == "control"
                        Ad = exp(A*ts)
                        Bd = A \ (Ad - C) * B
                        global sys_d = ss(Ad, Bd, C, D, ts)
                    end
                    
                    # --- #

                    if methods[n] == "env_without_agent"
                        RLBase.reset!(env)
                        #global actions = rand(Float64, ( na, length(t)-1 )) .*2 .-1
                        global u = [0.76666 for i = 1:length(t)]
                        global uu = [u for i = 1:na ]
                        global actions = mapreduce(permutedims, vcat, uu)
                        timer = @benchmark execute_env(env, actions, false) samples = repeat evals = loops seconds = 1000
                        if debug
                            result = execute_env(env, actions, true)
                        end
                    elseif methods[n] == "env_with_agent"
                        RLBase.reset!(env)
                        timer = @benchmark execute_env(env, agent, length(t), false) samples = repeat evals = loops seconds = 1000
                        if debug
                            result = execute_env(env, agent, length(t), true)
                        end
                    elseif methods[n] == "lsoda"
                        p = [230.0 for i = 1:na]
                        tspan = (0.0,t_end[l])
                        u0 = [0.0 for i = 1:ns]
                        global problem = ODEProblem(f!,u0,tspan,p)
                        timer = @benchmark solve(problem,lsoda(), reltol=1e-6, abstol=1e-6) samples = repeat evals = loops seconds = 1000
                        if debug
                            result = solve(problem,lsoda(), reltol=1e-6, abstol=1e-6)
                            result = [result.u[h][2] for h = 1:length(result)]
                        end
                    elseif methods[n] == "control"
                        global x0 = [0.0 for i = 1:ns]
                        #global u = rand(Float64, ( length(t) )) .*2 .-1
                        global u = [230.0 for i = 1:length(t)]
                        global uu = [u for i = 1:na ]
                        global uuu = mapreduce(permutedims, vcat, uu)
                        global ttt = t
                        #yout = lsim(sys,uuu,ttt,x0=x0)
                        timer = @benchmark lsim(sys_d,uuu,ttt,x0=x0) samples = repeat evals = loops seconds = 1000
                        if debug
                            result, _, _, _ = lsim(sys_d,uuu,ttt,x0=x0)
                            result = result'[:,2]
                        end
                    end

                    append!(time_list, timer.times)
                    #println("done for c=$(c) with time_list $(time_list)")

                    if debug
                        display(plot(result, title="$(methods[n]), $(num_nodes[k]), $(t_end[l]), $(c)"))
                    end
                end

                

                if timer !== nothing
                    t_result_mean[l,k,n] = mean(time_list) / 1_000_000_000
                    t_result_std[l,k,n] = std(time_list) / 1_000_000_000
                    t_num_samples[l,k,n] = length(time_list)
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