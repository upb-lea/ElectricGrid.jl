using ReinforcementLearning
using IntervalSets
using LinearAlgebra
using ControlSystems
using CUDA
using DataStructures

include("./custom_control.jl")
include("./nodeconstructor.jl")


mutable struct SimEnv <: AbstractEnv
    nc
    sys_d
    action_space
    state_space
    done
    featurize
    prepare_action
    reward_function
    x0
    x
    t0
    t
    ts
    state
    maxsteps
    steps
    state_ids
    v_dc
    norm_array
    convert_state_to_cpu
    reward
    action
    action_ids
    action_delay_buffer
end

function SimEnv(; maxsteps = 500, ts = 1/10_000, action_space = nothing, state_space = nothing, prepare_action = nothing, featurize = nothing, reward_function = nothing, CM = nothing, num_sources = nothing, num_loads = nothing, parameters = nothing, x0 = nothing, t0 = 0.0, state_ids = nothing, v_dc = nothing, norm_array = nothing, convert_state_to_cpu = true, use_gpu = false, reward = nothing, action = nothing, action_ids = nothing, action_delay = 0)
    
    if !(isnothing(CM) || isnothing(num_sources) || isnothing(num_loads))

        if isnothing(parameters)
            nc = NodeConstructor(num_sources = num_sources, num_loads = num_loads, CM = CM)
        else
            nc = NodeConstructor(num_sources = num_sources, num_loads = num_loads, CM = CM, parameters = parameters)
        end

        A, B, C, D = get_sys(nc)
        Ad = exp(A*ts)
        Bd = A \ (Ad - C) * B

        if use_gpu
            Ad = CuArray(A)
            Bd = CuArray(B)
            C = CuArray(C)
            if isa(D, Array)
                D = CuArray(D)
            end
        end

        sys_d = HeteroStateSpace(Ad, Bd, C, D, Float64(ts))

    else
        # Construct standard env with 2 sources, 1 load
        println("INFO: Three phase electric power grid with 2 sources and 1 load is created! Parameters are drawn randomly! To change, please define parameters (see nodeconstructor)")
        CM = [ 0. 0. 1.
               0. 0. 2
              -1. -2. 0.]

        if isnothing(parameters)
            nc = NodeConstructor(num_sources = 2, num_loads = 1, CM = CM)
        else
            nc = NodeConstructor(num_sources = 2, num_loads = 1, CM = CM, parameters = parameters)
        end

        A, B, C, D = get_sys(nc)
        Ad = exp(A*ts)
        Bd = A \ (Ad - C) * B
        sys_d = HeteroStateSpace(Ad, Bd, C, D, Float64(ts))

    end

    if isnothing(action_space)
        action_space = Space([ -1.0..1.0 for i = 1:length(sys_d.B[1,:]) ], )
    end

    if isnothing(featurize)
        featurize = function(x0 = nothing, t0 = nothing; env = nothing, name = nothing)
            if !isnothing(name)
                return env.state
            elseif isnothing(env)
                return x0
            else
                return env.state
            end
        end
    end

    if isnothing(prepare_action)
        prepare_action = function(env) 
            env.action
        end
    end

    if isnothing(reward_function)
        reward_function = function(env, name = nothing) 
            return 0.0
        end
    end

    if isnothing(x0)
        x0 = [ 0.0 for i = 1:length(sys_d.A[1,:]) ]
    end

    x = x0
    t = t0

    state = featurize(x0,t0)

    if isnothing(state_space)
        state_space = Space([ -1.0..1.0 for i = 1:length(state) ], )
    end

    if use_gpu
        if isa(x0, Array)
            x0 = CuArray(x0)
        end
        if !(convert_state_to_cpu) && isa(state, Array)
            state = CuArray(state)
        end
    end

    if isnothing(state_ids)
        if isnothing(nc)
            state_ids = []
            println("WARNING: No state_ids array specified - observing states with DataHook not possible")
        else
            state_ids = get_state_ids(nc)
        end
    end

    if isnothing(action_ids)
        if isnothing(nc)
            action_ids = []
            println("WARNING: No state_ids array specified - observing states with DataHook not possible")
        else
            action_ids = get_action_ids(nc)
        end
    end
    # TODO: take vdc per source from the nc.parameters[] something like:
    # v_dc = [env.nc.parameters["sources"][n]["vdc"] for n = 1:env.nc.num_sources]
    if isnothing(v_dc)
        println("INFO: v_dc = 350V will get applied to all actions")
        v_dc = 350 * ones(length(action_space))
    elseif isa(v_dc, Number)
        println("INFO: v_dc = $(v_dc)V will get applied to all actions")
        v_dc = v_dc * ones(length(action_space))
    end
    
    #TODO: norm_array from parameters Dict
    if isnothing(norm_array)
        if isnothing(nc)
            println("INFO: norm_array set to ones - if neccessary please define norm_array in env initialization")
            norm_array = ones(length(sys_d.A[1,:]))
        else
            println("INFO: Generating standard norm_array from nodeconstructor")
            states = get_state_ids(nc)
            norm_array = []
            println("WARNING: limits set to fixed value - define in nc.parameters")
            for state_name in states
                if occursin("_i", state_name)
                    #push!(norm_array, limits["i_lim"])
                    push!(norm_array, 20.0)
                elseif occursin("_u", state_name)
                    #push!(norm_array, limits["v_lim"])
                    push!(norm_array, 600.0)
                end
            end
        end
    end

    if isnothing(reward)
        reward = 0.0
    end

    if isnothing(action)
        action = zeros(length(action_space))
    end

    if action_delay == 0
        action_delay_buffer = nothing
    else
        action_delay_buffer = CircularBuffer{Vector{Float64}}(action_delay)
        fill!(action_delay_buffer, zeros(length(action_space)))
    end

    SimEnv(nc, sys_d, action_space, state_space, false, featurize, prepare_action, reward_function, x0, x, t0, t, ts, state, maxsteps, 0, state_ids, v_dc, norm_array, convert_state_to_cpu, reward, action, action_ids, action_delay_buffer)
end

RLBase.action_space(env::SimEnv) = env.action_space
RLBase.state_space(env::SimEnv) = env.state_space
RLBase.reward(env::SimEnv) =  env.reward

function RLBase.reward(env::SimEnv, name::String)
    return env.reward_function(env, name)
end

RLBase.is_terminated(env::SimEnv) = env.done
RLBase.state(env::SimEnv) = env.state

function RLBase.state(env::SimEnv, name::String)
    return env.featurize(;env = env, name = name)
end

function RLBase.reset!(env::SimEnv)
    env.state = env.convert_state_to_cpu ? Array(env.featurize(env.x0, env.t0)) : env.featurize(env.x0, env.t0)
    env.x = env.x0
    if !isnothing(env.action_delay_buffer)
        empty!(env.action_delay_buffer)
        fill!(env.action_delay_buffer, zeros(length(env.action_space)))
    end
    env.t = env.t0
    env.steps = 0
    env.reward = 0.0
    env.done = false
    nothing
end

function (env::SimEnv)(action)
    env.steps += 1
    # i_1, v_1 , i_L = 
    # why tt??
    tt = [env.t, env.t + env.ts]

    env.t = tt[2]

    if !isnothing(env.action_delay_buffer)
        env.action = env.action_delay_buffer[1]
        push!(env.action_delay_buffer, action)
    else
        env.action = action
    end
    
    env.action = env.action .* env.v_dc/2

    env.action = env.prepare_action(env)
    
    if env.sys_d.A isa CuArray && !(env.action isa CuArray)
        if env.action isa Array
            env.action = CuArray(env.action)
        else
            env.action = CuArray([env.action])
        end
    end
    u = [env.action env.action]

    xout_d = custom_lsim(env.sys_d, u, tt, x0=env.x)
    #xout_d = [env.x env.x]

    env.x = xout_d[:,2]
    #env.x = xout_d'[2,:]
    if env.convert_state_to_cpu
        env.state = Array(xout_d)'[2,:] ./ env.norm_array
    else
        env.state = xout_d'[2,:] ./ env.norm_array
    end

    env.state = env.featurize(; env = env)

    # reward
    # loss_error = 1e-1
    # hardcoded values - change later
    # use functions outside this reward function - normalised
    # P_load = (env.norm_array[end] * env.state[end])^2 / 14
    # push!(PLoad, P_load)
    # P_diff = -abs(P_required - P_load)   

    # env.reward = -sqrt((P_source - (P_R + P_load + loss_error))^2)
    # Power constraint
    env.reward = env.reward_function(env)

    env.done = env.steps >= env.maxsteps
end