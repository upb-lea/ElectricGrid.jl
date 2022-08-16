using ReinforcementLearning
using IntervalSets
using LinearAlgebra
using ControlSystems
using CUDA

include(srcdir("custom_control.jl"))

# required power at the load
P_required = 1000 # W

# required Voltage at transformer [v_1]
V_required = 230  # V

PLoad = []
# --- RL ENV ---

Base.@kwdef mutable struct SimEnv <: AbstractEnv
    A
    B
    C
    D = 0
    action_space::Space{Vector{ClosedInterval{Float64}}} = Space([ -1.0..1.0 for i = 1:length(B[1,:]) ], )
    observation_space::Space{Vector{ClosedInterval{Float64}}} = Space([ -1.0..1.0 for i = 1:length(A[1,:]) ], )
    done::Bool = false
    x0 = [ 0.0 for i = 1:length(A[1,:]) ]
    x = x0
    maxsteps::Int = 300
    steps::Int = 0
    t::Rational = 0
    ts::Rational = 1//10_000
    sys_d = HeteroStateSpace(exp(A*ts), A \ (Ad - C) * B, C, D, Float64(ts))
    state_ids::Vector{String}
    norm_array::Vector{Float64} = [ 600.0 for i = 1:length(A[1,:]) ]
    v_dc::Float64 = 300
    reward::Float64 = 0
    convert_state_to_cpu::Bool = true
    state = convert_state_to_cpu ? Array(x0) : x0
end

RLBase.action_space(env::SimEnv) = env.action_space
RLBase.state_space(env::SimEnv) = env.observation_space
RLBase.reward(env::SimEnv) =  env.reward 


RLBase.is_terminated(env::SimEnv) = env.done
RLBase.state(env::SimEnv) = env.state

function RLBase.reset!(env::SimEnv)
    env.state = env.convert_state_to_cpu ? Array(env.x0) : env.x0
    env.x = env.x0
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

    action *= env.v_dc
    if env.Ad isa CuArray && !(action isa CuArray)
        if action isa Array
            action = CuArray(action)
        else
            action = CuArray([action])
        end
    end
    u = [action action]

    xout_d = custom_lsim(env.sys_d, u, tt, x0=env.x)
    #xout_d = [env.x env.x]

    env.x = xout_d[:,2]
    #env.x = xout_d'[2,:]
    if env.convert_state_to_cpu
        env.state = Matrix(xout_d)'[2,:] ./ env.norm_array
    else
        env.state = xout_d'[2,:] ./ env.norm_array
    end

    # reward
    # loss_error = 1e-1
    # hardcoded values - change later
    # use functions outside this reward function - normalised
    # P_load = (env.norm_array[end] * env.state[end])^2 / 14
    # push!(PLoad, P_load)
    # P_diff = -abs(P_required - P_load)   

    # env.reward = -sqrt((P_source - (P_R + P_load + loss_error))^2)
    # Power constraint
    env.reward = 1#reward_func("Power_exp", env)

    env.done = env.steps >= env.maxsteps
end