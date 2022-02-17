using ReinforcementLearning
using IntervalSets
using LinearAlgebra
using ControlSystems

# --- RL ENV ---

Base.@kwdef mutable struct SimEnv <: AbstractEnv
    action_space::Space{Vector{ClosedInterval{Float64}}} = Space([0.0..240.0, 0.0..240.0], )
    observation_space::Space{Vector{ClosedInterval{Float64}}} = Space([ typemin(Float64)..typemax(Float64),
                                                                        typemin(Float64)..typemax(Float64),
                                                                        typemin(Float64)..typemax(Float64),
                                                                        typemin(Float64)..typemax(Float64),
                                                                        typemin(Float64)..typemax(Float64),
                                                                        typemin(Float64)..typemax(Float64) ], )
    done::Bool = false
    x0::Vector{Float64} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    x::Vector{Float64} = x0
    state::Vector{Float64} = x0
    maxsteps::Int = 1000
    steps::Int = 0
    t::Rational = 0
    ts::Rational = 1//10_000
    A = 0
    B = 0
    C = 0
    D = 0
    Ad::AbstractMatrix = exp(A*ts)
    Bd::AbstractMatrix = A \ (Ad - C) * B
    sys_d::StateSpace = ss(Ad, Bd, C, 0, Float64(ts))
end


RLBase.action_space(env::SimEnv) = env.action_space
RLBase.state_space(env::SimEnv) = env.observation_space
RLBase.reward(env::SimEnv) = abs(env.state[2] - 150.0)
RLBase.is_terminated(env::SimEnv) = env.done
RLBase.state(env::SimEnv) = env.state

function RLBase.reset!(env::SimEnv) where {A,T}
    env.state = x0
    env.t = 0
    env.steps = 0
    env.done = false
    nothing
end

function (env::SimEnv)(action)
    env.steps += 1

    tt = [env.t, env.t + env.ts]

    env.t = tt[2]

    u = [action action]

    yout_d, tout_d, xout_d, uout_d = lsim(env.sys_d, u, tt, x0=env.x)

    env.x = xout_d'[2,:]
    env.state = yout_d'[2,:]

    env.done = env.steps >= env.maxsteps
end