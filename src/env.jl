using ReinforcementLearning
using IntervalSets
using LinearAlgebra
using ControlSystems
using CUDA

include("./custom_control.jl")
include("./nodeconstructor.jl")


mutable struct SimEnv <: AbstractEnv
    nc
    action_space
    state_space
    done
    x0
    x
    maxsteps
    steps
    t
    ts
    sys_d
    state_ids
    prepare_data
    featurize
    rewardfunction
    norm_array
    v_dc
    reward
    convert_state_to_cpu
    state
    action
end

function SimEnv(; maxsteps = 500, ts = 1/10_000, action_space = nothing, state_space = nothing, prepare_data = nothing, featurize = nothing, rewardfunction = nothing, sys_d = nothing, Ad = nothing, Bd = nothing, A = nothing, B = nothing, C = nothing, D = nothing, CM = nothing, num_sources = nothing, num_loads = nothing, parameters = nothing, x0 = nothing, t = 0.0)
    
    # TODO: Shall we store nc in env to make it accessible later on to check parameters? (@janstenner)
    if !(isnothing(sys_d))
        # take sys_d... what to do here?
    elseif !(isnothing(Ad) || isnothing(Bd) || isnothing(C) || isnothing(D))

        sys_d = HeteroStateSpace(Ad, Bd, C, D, Float64(ts))

    elseif !(isnothing(A) || isnothing(B) || isnothing(C) || isnothing(D))

        Ad = exp(A*ts)
        Bd = A \ (Ad - C) * B
        sys_d = HeteroStateSpace(Ad, Bd, C, D, Float64(ts))

    elseif !(isnothing(CM) || isnothing(num_sources) || isnothing(num_loads))

        nc = NodeConstructor(num_sources = num_sources, num_loads = num_loads, CM = CM)

        A, B, C, D = get_sys(nc)
        Ad = exp(A*ts)
        Bd = A \ (Ad - C) * B
        sys_d = HeteroStateSpace(Ad, Bd, C, D, Float64(ts))

    else
        # Construct standard env with 2 sources, 1 load
        println("Three phase electric power grid with 2 sources and 1 load is created! Parameters are drwan randomly! To change, please define parameters (see nodeconstructor)")
        CM = [ 0. 0. 1.
               0. 0. 2
              -1. -2. 0.]

        nc = NodeConstructor(num_sources = 2, num_loads = 1, CM = CM)

        A, B, C, D = get_sys(nc)
        Ad = exp(A*ts)
        Bd = A \ (Ad - C) * B
        sys_d = HeteroStateSpace(Ad, Bd, C, D, Float64(ts))

    end

    A
    B
    C
    D = 0
    action_space::Space{Vector{ClosedInterval{Float64}}} = Space([ -1.0..1.0 for i = 1:length(B[1,:]) ], )
    state_space::Space{Vector{ClosedInterval{Float64}}} = Space([ -1.0..1.0 for i = 1:length(A[1,:]) ], )
    done::Bool = false
    x0 = [ 0.0 for i = 1:length(A[1,:]) ]
    x = x0
    maxsteps::Int = 300
    steps::Int = 0
    t = 0
    ts = 1/10_000
    Ad = exp(A*ts)
    Bd = A \ (Ad - C) * B
    sys_d = HeteroStateSpace(Ad, Bd, C, D, Float64(ts))
    state_ids::Vector{String}
    rewardfunction
    featurize = [ 0.0 for i = 1:length(state_space) ]
    norm_array::Vector{Float64} = [ 600.0 for i = 1:length(A[1,:]) ]
    v_dc::Float64 = 300
    reward::Float64 = 0
    convert_state_to_cpu::Bool = true
    state = convert_state_to_cpu ? Array(featurize(env)) : 
end

RLBase.action_space(env::SimEnv) = env.action_space
RLBase.state_space(env::SimEnv) = env.state_space
RLBase.reward(env::SimEnv) =  env.reward 


RLBase.is_terminated(env::SimEnv) = env.done
RLBase.state(env::SimEnv) = env.state

function RLBase.reset!(env::SimEnv)
    env.state = env.convert_state_to_cpu ? Array(env.x0) : env.x0
    env.x = env.x0
    if !(isnothing(featurize))
        featurize(env)
    end
    env.t = 0
    env.steps = 0
    env.reward = 0
    env.done = false
    nothing
end

function (env::SimEnv)(action)
    env.steps += 1
    # i_1, v_1 , i_L = 
    # why tt??
    tt = [env.t, env.t + env.ts]

    env.t = tt[2]

    env.action = action
    prepare_data(env)

    env.action *= env.v_dc
    if env.Ad isa CuArray && !(env.action isa CuArray)
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

    featurize(env)

    # reward
    # loss_error = 1e-1
    # hardcoded values - change later
    # use functions outside this reward function - normalised
    # P_load = (env.norm_array[end] * env.state[end])^2 / 14
    # push!(PLoad, P_load)
    # P_diff = -abs(P_required - P_load)   

    # env.reward = -sqrt((P_source - (P_R + P_load + loss_error))^2)
    # Power constraint
    env.reward = env.rewardfunction(env)

    env.done = env.steps >= env.maxsteps
end