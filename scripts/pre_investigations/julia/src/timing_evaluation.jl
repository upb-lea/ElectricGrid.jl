using PyCall
@pyinclude(srcdir("nodeconstructor.py"))
@pyinclude(srcdir("nodeconstructorcable.py"))

using BenchmarkTools
using Statistics
using DifferentialEquations
using LSODA
using ControlSystems
using JSON
using Plots
using Flux
using StableRNGs
using CUDA

include(srcdir("env.jl"))
include(srcdir("custom_control.jl"))

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
    if debug
        output = zeros(length(env.A[1,:]),length(actions[1,:])+1)
    else
        output = 0.0
    end

    RLBase.reset!(env)
    for i = 1:length(actions[1,:])
        env(actions[:,i])
        if debug output[:,i+1] = env.state.*env.norm_array end
    end

    return output
end

function execute_env(env::SimEnv, agent::Agent, t_len::Int, debug::Bool)
    if debug
        output = zeros(length(env.Ad[1,:]),t_len+1)
    else
        output = 0.0
    end

    RLBase.reset!(env)
    for i = 1:t_len
        action = agent(env)
        env(action)
        if debug output[:,i+1] = env.state.*env.norm_array end
    end

    return output
end

function execute_env_cuda(env::SimEnv, agent::Agent, t_len::Int, debug::Bool)
    CUDA.@sync begin
        println("ho")
        if debug
            output = zeros(length(env.Ad[1,:]),t_len+1)
        else
            output = 0.0
        end

        RLBase.reset!(env)
        for i = 1:t_len
            action = CuArray(agent(env))
            env(action)
            if debug
                output[:,i+1] = Array(env.state).*env.norm_array
            end
        end
    end
    return output
end

function control_stepwise(sys_d,uuu,ttt,x0)
    tempx = x0
    for i = 1:length(uuu[1,:])-1
        xout = custom_lsim(sys_d,uuu[:,i:i+1],ttt[i:i+1],x0=tempx)
        tempx = xout[:,2]
    end
end

function cudalsim(sys_d,uuu,ttt,x0,stepwise)
    CUDA.@sync begin
        if stepwise
            tempx = x0
            for i = 1:length(uuu[1,:])-1
                xout = custom_lsim(sys_d,uuu[:,i:i+1],ttt[i:i+1],x0=tempx)
                tempx = xout[:,2]
            end
        else
            custom_lsim(sys_d,uuu,ttt,x0=x0)
        end
    end
end

function cudasolve(problem, tol, tss)
    CUDA.@sync begin
        solve(problem, BS3(), reltol=tol, abstol=tol, saveat=tss)
    end
end

function timing_experiment_simulation(repeat::Int64=5, loops::Int64=10, num_nodes::Vector{Int64}=[5,10,25,50,100],
                                        t_end::Vector{Float64}=[0.001,0.01,0.1], num_cm::Int64=1, ts::Float64=1e-4,
                                        methods::Vector{String}=["control", "env_without_agent", "env_with_agent", "lsoda"],
                                        parameter::Dict=Dict(), debug::Bool=false, stepwise::Bool=true, fixed_stepsize::Bool=false)
    
    if isempty(parameter)
        parameter["R_source"] = 0.4
        parameter["L_source"] = 2.3e-3
        parameter["C_source"] = 10e-6
        parameter["L_cable"] = 2.3e-3
        parameter["R_cable"] = 0.4
        parameter["R_load"] = 14
        parameter["V_dc"] = 300
    end

    limits = Dict("i_lim" => 20, "v_lim" => 600)
    ref = 200

    if debug
        t_result_mean = zeros(num_cm, length(t_end), length(num_nodes), length(methods))
        t_result_std = zeros(num_cm, length(t_end), length(num_nodes), length(methods))
        t_num_samples = zeros(Int64, num_cm, length(t_end), length(num_nodes), length(methods))
        #results = Array{Matrix}(undef, num_cm, length(t_end), length(num_nodes), length(methods))
        results = Dict()
    else
        t_result_mean = zeros(length(t_end), length(num_nodes), length(methods))
        t_result_std = zeros(length(t_end), length(num_nodes), length(methods))
        t_num_samples = zeros(Int64, length(t_end), length(num_nodes), length(methods))
    end


    for n = 1:length(methods)
        for k = 1:length(num_nodes)

            CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(num_nodes[k]) * ".json"))

            for l = 1:length(t_end)
                timer = nothing
                t = collect(0:ts:t_end[l])
                time_list = []

                for c = 1:num_cm
                    CM = reduce(hcat, CM_list[c])'
                    CM = convert(Matrix{Int}, CM)

                    #nc = py"NodeConstructor"(num_nodes[k], num_nodes[k], parameter, CM=CM)
                    nc = py"NodeConstructorCable"(num_nodes[k], num_nodes[k], CM=CM)
                    A, B, C, D = nc.get_sys()
                    ns = length(A[1,:])
                    na = length(B[1,:])

                    #println("Size of x: $(ns)")
                    #println("Size of u: $(na)")

                    if debug
                        result = nothing
                        resulttoplot = nothing
                    end

                    if methods[n] == "env_without_agent"
                        norm_array = vcat( [limits[i] for j = 1:nc.num_source for i in ["i_lim", "v_lim"]], [limits["i_lim"] for i = 1:nc.num_connections] )
                        global env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, v_dc=parameter["V_dc"], ts=rationalize(ts))
                    elseif methods[n] == "env_with_agent" || methods[n] == "env_with_agent_cuda"
                        norm_array = vcat( [limits[i] for j = 1:nc.num_source for i in ["i_lim", "v_lim"]], [limits["i_lim"] for i = 1:nc.num_connections] )
                        if methods[n] == "env_with_agent_cuda"
                            Ad = exp(A*ts)
                            Bd = A \ (Ad - C) * B
                            Ad = CuArray(Ad)
                            Bd = CuArray(Bd)
                            C = CuArray(C)
                            D = CUDA.zeros(Float64, size(Bd))
                            x0 = CuArray([ 0.0 for i = 1:length(A[1,:]) ])
                            global env = SimEnv(Ad=Ad, Bd=Bd, C=C, D=D, x0=x0, norm_array=norm_array, v_dc=parameter["V_dc"], ts=rationalize(ts))
                        else
                            global env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, v_dc=parameter["V_dc"], ts=rationalize(ts))
                        end
                        global agent = create_agent(na, ns)
                    elseif methods[n] == "expm1" || methods[n] == "expm2" || methods[n] == "expm4" || methods[n] == "expm8" || methods[n] == "expm16" || methods[n] == "expm32"
                        global tss = ts
                        global AA = A

                        if methods[n] == "expm1"
                            BLAS.set_num_threads(1)
                        elseif methods[n] == "expm2"
                            BLAS.set_num_threads(2)
                        elseif methods[n] == "expm4"
                            BLAS.set_num_threads(4)
                        elseif methods[n] == "expm8"
                            BLAS.set_num_threads(8)
                        elseif methods[n] == "expm16"
                            BLAS.set_num_threads(16)
                        else
                            BLAS.set_num_threads(32)
                        end
                    elseif methods[n] == "lsoda" || methods[n] == "DP5" || methods[n] == "BS3" || methods[n] == "BS3_32" || methods[n] == "BS3_CUDA" || methods[n] == "BS3_CUDA32"
                        global tss = ts

                        if methods[n] == "BS3_32" || methods[n] == "BS3_CUDA32"
                            global tss = Float32(tss)
                            A = Float32.(A)
                            B = Float32.(B)
                        end

                        if methods[n] == "BS3_CUDA" || methods[n] == "BS3_CUDA32"
                            A = CuArray(A)
                            B = CuArray(B)
                        end

                        function f!(du, u, p, t)
                            du[:] = A * u + B * p
                        end
                    elseif methods[n] == "control" || methods[n] == "control1" || methods[n] == "control16" || methods[n] == "control2" || methods[n] == "control4"
                        Ad = exp(A*ts)
                        Bd = A \ (Ad - C) * B
                        global sys_d = ss(Ad, Bd, C, D, ts)
                        if methods[n] == "control"
                            BLAS.set_num_threads(8)
                        elseif methods[n] == "control1"
                            BLAS.set_num_threads(1)
                        elseif methods[n] == "control2"
                            BLAS.set_num_threads(2)
                        elseif methods[n] == "control4"
                            BLAS.set_num_threads(4)
                        else
                            BLAS.set_num_threads(16)
                        end
                    elseif methods[n] == "control32bit"
                        BLAS.set_num_threads(8)
                        Ad = exp(A*ts)
                        Bd = A \ (Ad - C) * B
                        Ad = Float32.(Ad)
                        Bd = Float32.(Bd)
                        C = Float32.(C)
                        D = zeros(Float32, size(Bd))
                        global sys_d = HeteroStateSpace(Ad, Bd, C, D, ts)
                    elseif methods[n] == "controlCUDA"
                        print("a")
                        Ad = exp(A*ts)
                        Bd = A \ (Ad - C) * B
                        Ad = CuArray(Ad)
                        Bd = CuArray(Bd)
                        C = CuArray(C)
                        D = CUDA.zeros(Float64, size(Bd))
                        global sys_d = HeteroStateSpace(Ad, Bd, C, D, ts)
                    elseif methods[n] == "controlCUDA32"
                        Ad = exp(A*ts)
                        Bd = A \ (Ad - C) * B
                        Ad = CuArray(Float32.(Ad))
                        Bd = CuArray(Float32.(Bd))
                        C = CuArray(Float32.(C))
                        D = CUDA.zeros(Float32, size(Bd))
                        global sys_d = HeteroStateSpace(Ad, Bd, C, D, ts)
                    end
                    
                    # --- #

                    if methods[n] == "env_without_agent"
                        RLBase.reset!(env)
                        #global actions = rand(Float64, ( na, length(t)-1 )) .*2 .-1
                        global u = [0.76666 for i = 1:length(t)]
                        global uu = [u for i = 1:na ]
                        global actions = mapreduce(permutedims, vcat, uu)
                        execute_env(env, actions, false)
                        timer = @benchmark execute_env(env, actions, false) samples = repeat evals = loops seconds = 1000
                        if debug
                            result = execute_env(env, actions, true)
                            resulttoplot = result[2,:]
                        end
                    elseif methods[n] == "env_with_agent"
                        RLBase.reset!(env)
                        global tt = t
                        execute_env(env, agent, length(tt), false)
                        timer = @benchmark execute_env(env, agent, length(tt), false) samples = repeat evals = loops seconds = 1000
                        if debug
                            result = execute_env(env, agent, length(tt), true)
                            resulttoplot = result[2,:]
                        end
                    elseif methods[n] == "env_with_agent_cuda"
                        RLBase.reset!(env)
                        global tt = t
                        execute_env_cuda(env, agent, length(tt), false)
                        timer = @benchmark execute_env_cuda(env, agent, length(tt), false) samples = repeat evals = loops seconds = 1000
                        if debug
                            result = execute_env_cuda(env, agent, length(tt), true)
                            resulttoplot = result[2,:]
                        end
                    elseif methods[n] == "expm1" || methods[n] == "expm2" || methods[n] == "expm4" || methods[n] == "expm8" || methods[n] == "expm16" || methods[n] == "expm32"
                        exp(AA*tss)
                        timer = @benchmark exp(AA*tss) samples = repeat evals = loops seconds = 1000
                    elseif methods[n] == "lsoda" || methods[n] == "DP5" || methods[n] == "BS3" || methods[n] == "BS3_32" || methods[n] == "BS3_CUDA" || methods[n] == "BS3_CUDA32"
                        p = [230.0 for i = 1:na]
                        tspan = (0.0,t_end[l])
                        u0 = [0.0 for i = 1:ns]

                        if methods[n] == "BS3_32" || methods[n] == "BS3_CUDA32"
                            p = Float32.(p)
                            tspan = (0.0f0,Float32(t_end[l]))
                            u0 = Float32.(u0)
                            global tol = Float32(1e-6)
                        else
                            global tol = 1e-6
                        end

                        if methods[n] == "BS3_CUDA" || methods[n] == "BS3_CUDA32"
                            p = CuArray(p)
                            u0 = CuArray(u0)
                        end

                        global problem = ODEProblem(f!,u0,tspan,p)
                        if methods[n] == "lsoda"
                            if fixed_stepsize
                                solve(problem, lsoda(), adaptive=false, dt=tss, dtmin=tss, dtmax=tss)
                                timer = @benchmark solve(problem, lsoda(), adaptive=false, dt=tss, dtmin=tss, dtmax=tss) samples = repeat evals = loops seconds = 1000
                            else
                                solve(problem, lsoda(), reltol=tol, abstol=tol, saveat=tss)
                                timer = @benchmark solve(problem, lsoda(), reltol=tol, abstol=tol, saveat=tss) samples = repeat evals = loops seconds = 1000
                            end
                        elseif methods[n] == "DP5"
                            if fixed_stepsize
                                solve(problem, DP5(), adaptive=false, dt=tss, dtmin=tss, dtmax=tss)
                                timer = @benchmark solve(problem, DP5(), adaptive=false, dt=tss, dtmin=tss, dtmax=tss) samples = repeat evals = loops seconds = 1000
                            else
                                solve(problem, DP5(), reltol=tol, abstol=tol, saveat=tss)
                                timer = @benchmark solve(problem, DP5(), reltol=tol, abstol=tol, saveat=tss) samples = repeat evals = loops seconds = 1000
                            end
                        elseif methods[n] == "BS3" || methods[n] == "BS3_32"
                            if fixed_stepsize
                                solve(problem, BS3(), adaptive=false, dt=tss, dtmin=tss, dtmax=tss)
                                timer = @benchmark solve(problem, BS3(), adaptive=false, dt=tss, dtmin=tss, dtmax=tss) samples = repeat evals = loops seconds = 1000
                            else
                                solve(problem, BS3(), reltol=tol, abstol=tol, saveat=tss)
                                timer = @benchmark solve(problem, BS3(), reltol=tol, abstol=tol, saveat=tss) samples = repeat evals = loops seconds = 1000
                            end
                        else
                            cudasolve(problem, tol, tss)
                            timer = @benchmark cudasolve(problem, tol, tss) samples = repeat evals = loops seconds = 1000
                        end
                        if debug
                            if methods[n] == "lsoda"
                                if fixed_stepsize
                                    result = solve(problem, lsoda(), adaptive=false, dt=tss, dtmin=tss, dtmax=tss)
                                else
                                    result = solve(problem, lsoda(), reltol=tol, abstol=tol, saveat=tss)
                                end
                            elseif methods[n] == "DP5"
                                if fixed_stepsize
                                    result = solve(problem, DP5(), adaptive=false, dt=tss, dtmin=tss, dtmax=tss)
                                else
                                    result = solve(problem, DP5(), reltol=tol, abstol=tol, saveat=tss)
                                end
                            else
                                if fixed_stepsize
                                    result = solve(problem, BS3(), adaptive=false, dt=tss, dtmin=tss, dtmax=tss)
                                else
                                    result = solve(problem, BS3(), reltol=tol, abstol=tol, saveat=tss)
                                end
                            end
                            if methods[n] == "BS3_CUDA" || methods[n] == "BS3_CUDA32"
                                u_temp = []
                                for i = 1:length(result.u)
                                    append!(u_temp,[Array(result.u[i])])
                                end
                                resulttoplot = [u_temp[h][2] for h = 1:length(u_temp)]
                            else
                                resulttoplot = [result.u[h][2] for h = 1:length(result)]
                            end
                        end
                    elseif methods[n] == "control" || methods[n] == "control1" || methods[n] == "control16" || methods[n] == "control2" || methods[n] == "control4"
                        global x0 = [0.0 for i = 1:ns]
                        #global u = rand(Float64, ( length(t) )) .*2 .-1
                        global u = [230.0 for i = 1:length(t)]
                        global uu = [u for i = 1:na ]
                        global uuu = mapreduce(permutedims, vcat, uu)
                        global ttt = t
                        #yout = lsim(sys,uuu,ttt,x0=x0)
                        if stepwise
                            control_stepwise(sys_d,uuu,ttt,x0)
                            timer = @benchmark control_stepwise(sys_d,uuu,ttt,x0) samples = repeat evals = loops seconds = 1000
                        else
                            custom_lsim(sys_d,uuu,ttt,x0=x0)
                            timer = @benchmark custom_lsim(sys_d,uuu,ttt,x0=x0) samples = repeat evals = loops seconds = 1000
                        end
                        if debug
                            result, _, _, _ = lsim(sys_d,uuu,ttt,x0=x0)
                            resulttoplot = result'[:,2]
                        end
                    elseif methods[n] == "control32bit"
                        global x0 = [0.0 for i = 1:ns]
                        global u = [230.0 for i = 1:length(t)]
                        global uu = [u for i = 1:na ]
                        global uuu = mapreduce(permutedims, vcat, uu)
                        global ttt = t
                        global x0 = Float32.(x0)
                        global uuu = Float32.(uuu)

                        if stepwise
                            control_stepwise(sys_d,uuu,ttt,x0)
                            timer = @benchmark control_stepwise(sys_d,uuu,ttt,x0) samples = repeat evals = loops seconds = 1000
                        else
                            custom_lsim(sys_d,uuu,ttt,x0=x0)
                            timer = @benchmark custom_lsim(sys_d,uuu,ttt,x0=x0) samples = repeat evals = loops seconds = 1000
                        end
                        if debug
                            result, _, _, _ = lsim(sys_d,uuu,ttt,x0=x0)
                            resulttoplot = result'[:,2]
                        end
                    elseif methods[n] == "controlCUDA" || methods[n] == "controlCUDA32"
                        global x0 = [0.0 for i = 1:ns]
                        global u = [230.0 for i = 1:length(t)]
                        global uu = [u for i = 1:na ]
                        global uuu = mapreduce(permutedims, vcat, uu)
                        global ttt = t

                        if methods[n] == "controlCUDA"
                            global x0 = CuArray(x0)
                            global uuu = CuArray(uuu)
                            #global ttt = CuArray(ttt)
                        else
                            global x0 = CuArray(Float32.(x0))
                            global uuu = CuArray(Float32.(uuu))
                        end
                        
                        cudalsim(sys_d,uuu,ttt,x0,stepwise)
                        global stepwise_g = stepwise
                        timer = @benchmark cudalsim(sys_d,uuu,ttt,x0,stepwise_g) samples = repeat evals = loops seconds = 1000
                        if debug
                            result, _, _, _ = lsim(sys_d,uuu,ttt,x0=x0)
                            result = Matrix(result)
                            resulttoplot = result'[:,2]
                        end
                    end

                    if debug
                        t_result_mean[c,l,k,n] = mean(timer.times) / 1_000_000_000
                        t_result_std[c,l,k,n] = std(timer.times) / 1_000_000_000
                        t_num_samples[c,l,k,n] = length(timer.times)
                    else
                        append!(time_list, timer.times)
                        #println("done for c=$(c) with time_list $(time_list)")
                    end

                    if debug
                        results["$(methods[n])_c$(c)_nodes$(num_nodes[k])_tend$(t_end[l])"] = result'
                        display(plot(resulttoplot, title="$(methods[n]), $(num_nodes[k]), $(t_end[l]), $(c)"))
                    end
                end

                

                if timer !== nothing && !debug
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
    result_dict["name"] = "julia_"
    if debug result_dict["name"] *= "debug_" end
    result_dict["name"] *= "$(join(methods,","))_$(join(num_nodes,","))_$(join(t_end,","))_cm=$(num_cm)_loops=$(loops)_repeat=$(repeat)"
    if debug result_dict["y_debug"] = results end
    
    result_dict
end