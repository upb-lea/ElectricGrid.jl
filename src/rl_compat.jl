#=
Minimal vendored subset of the ReinforcementLearning.jl v0.10-era API
(ReinforcementLearningBase v0.9 / ReinforcementLearningCore v0.8, MIT license,
https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl).

ElectricGrid was written against RL.jl v0.10, whose dependency stack (CUDA 3.x,
Flux 0.12) cannot be installed on Julia >= 1.10. ElectricGrid only uses a small,
stable slice of that API (the env interface, Agent/NamedPolicy wrappers, the
SART trajectory, stages/hooks and the run loop), so that slice is vendored here
with the same semantics instead of depending on RL.jl.

Inside ElectricGrid, `RLBase`, `RLCore` and `ReinforcementLearningCore` are all
aliases for this module (see ElectricGrid.jl).
=#
module RLCompat

using CircularArrayBuffers: CircularArrayBuffer
using IntervalSets: ClosedInterval, leftendpoint, rightendpoint
using Random
import StatsBase: sample
import Flux
import CUDA
import CUDA: device
using MacroTools: @forward

export AbstractEnv, AbstractPolicy, AbstractHook, AbstractStage,
    AbstractTrajectory, AbstractApproximator
export PreExperimentStage, PostExperimentStage, PreEpisodeStage,
    PostEpisodeStage, PreActStage, PostActStage
export PRE_EXPERIMENT_STAGE, POST_EXPERIMENT_STAGE, PRE_EPISODE_STAGE,
    POST_EPISODE_STAGE, PRE_ACT_STAGE, POST_ACT_STAGE
export Space
export Agent, NamedPolicy, RandomPolicy, NeuralNetworkApproximator
export Trajectory, CircularArraySARTTrajectory, BatchSampler, SART, SARTS
export StopAfterEpisode, StopAfterStep
export send_to_device, send_to_host
export legal_action_space_mask

#####
# Abstract types and stages
#####

abstract type AbstractEnv end
abstract type AbstractPolicy end
abstract type AbstractHook end
abstract type AbstractTrajectory end
abstract type AbstractApproximator end

abstract type AbstractStage end
struct PreExperimentStage <: AbstractStage end
struct PostExperimentStage <: AbstractStage end
struct PreEpisodeStage <: AbstractStage end
struct PostEpisodeStage <: AbstractStage end
struct PreActStage <: AbstractStage end
struct PostActStage <: AbstractStage end

const PRE_EXPERIMENT_STAGE = PreExperimentStage()
const POST_EXPERIMENT_STAGE = PostExperimentStage()
const PRE_EPISODE_STAGE = PreEpisodeStage()
const POST_EPISODE_STAGE = PostEpisodeStage()
const PRE_ACT_STAGE = PreActStage()
const POST_ACT_STAGE = PostActStage()

# Default no-op stage callbacks, as in old RLCore.
(p::AbstractPolicy)(::AbstractStage, ::AbstractEnv) = nothing
(p::AbstractPolicy)(::PreActStage, ::AbstractEnv, action) = nothing
(h::AbstractHook)(args...) = nothing

struct EmptyHook <: AbstractHook end

#####
# Env interface (extended by ElectricGrid via `RLBase.state(env::...) = ...` etc.)
#####

function state end
function state_space end
function action_space end
function reward end
function reset! end
function is_terminated end
function update! end
function DynamicStyle end
function legal_action_space_mask end

struct Simultaneous end
struct Sequential end

#####
# Space (old RLBase.Space)
#####

struct Space{T}
    s::T
end

Base.length(s::Space) = length(s.s)
Base.getindex(s::Space, i...) = getindex(s.s, i...)
Base.iterate(s::Space, args...) = iterate(s.s, args...)
Base.in(x, s::Space) = all(((xi, si),) -> xi in si, zip(x, s.s))
Random.rand(rng::AbstractRNG, s::Space) = map(x -> _rand_in(rng, x), s.s)

# uniform sampling from a ClosedInterval without pirating Random.rand
_rand_in(rng::AbstractRNG, i::ClosedInterval) =
    leftendpoint(i) + rand(rng) * (rightendpoint(i) - leftendpoint(i))
_rand_in(rng::AbstractRNG, x) = rand(rng, x)

#####
# Trajectories (old RLCore CircularArraySARTTrajectory)
#####

const SART = (:state, :action, :reward, :terminal)
const SARTS = (:state, :action, :reward, :terminal, :next_state)

struct Trajectory{T} <: AbstractTrajectory
    traces::T
end

Base.getindex(t::Trajectory, s::Symbol) = getindex(t.traces, s)
Base.haskey(t::Trajectory, s::Symbol) = haskey(t.traces, s)
Base.keys(t::Trajectory) = keys(t.traces)

const CircularArraySARTTrajectory =
    Trajectory{<:NamedTuple{SART,<:Tuple{Vararg{CircularArrayBuffer}}}}

_scalar_eltype(T::Type) = T <: AbstractArray ? eltype(T) : T
_buffer_from_spec((T, shape)::Pair, capacity::Int) =
    CircularArrayBuffer{_scalar_eltype(T)}(shape..., capacity)

function CircularArraySARTTrajectory(;
        capacity::Int,
        state = Int => (),
        action = Int => (),
        reward = Float32 => (),
        terminal = Bool => (),
    )
    # As in old RLCore: state/action traces hold one extra (dummy) frame, so
    # that `next_state` of the last transition is available.
    Trajectory((
        state = _buffer_from_spec(state, capacity + 1),
        action = _buffer_from_spec(action, capacity + 1),
        reward = _buffer_from_spec(reward, capacity),
        terminal = _buffer_from_spec(terminal, capacity),
    ))
end

# number of stored frames = size along the last (frame) dimension
_nframes(cb::CircularArrayBuffer) = size(cb, ndims(cb))

Base.length(t::CircularArraySARTTrajectory) = _nframes(t[:terminal])

#####
# Batch sampling (old RLCore BatchSampler)
#####

struct BatchSampler{traces}
    batch_size::Int
end

select_last_dim(x::AbstractArray{T,N}, inds) where {T,N} =
    @views x[ntuple(_ -> (:), N - 1)..., inds]

function sample(rng::AbstractRNG, t::CircularArraySARTTrajectory, s::BatchSampler{SARTS})
    inds = rand(rng, 1:length(t), s.batch_size)
    batch = NamedTuple{SARTS}((
        select_last_dim(t[:state], inds),
        select_last_dim(t[:action], inds),
        select_last_dim(t[:reward], inds),
        select_last_dim(t[:terminal], inds),
        select_last_dim(t[:state], inds .+ 1),
    ))
    inds, batch
end

#####
# Device helpers (old RLCore extended CUDA.device the same way)
#####

device(x::AbstractArray) = Val(:cpu)
function device(x)
    ps = Flux.params(x)
    isempty(ps) ? Val(:cpu) : device(first(ps))
end

send_to_device(d) = x -> send_to_device(d, x)
send_to_device(::Val{:cpu}, x) = Flux.cpu(x)
send_to_device(::CUDA.CuDevice, x) = Flux.gpu(x)
send_to_host(x) = Flux.cpu(x)

#####
# NeuralNetworkApproximator (old RLCore)
#####

Base.@kwdef struct NeuralNetworkApproximator{M,O} <: AbstractApproximator
    model::M
    optimizer::O = nothing
end

(app::NeuralNetworkApproximator)(args...; kwargs...) = app.model(args...; kwargs...)

@forward NeuralNetworkApproximator.model Flux.testmode!,
Flux.trainmode!,
Flux.params,
device

Flux.functor(x::NeuralNetworkApproximator) =
    (model = x.model,), y -> NeuralNetworkApproximator(y.model, x.optimizer)

update!(app::NeuralNetworkApproximator, gs) =
    Flux.Optimise.update!(app.optimizer, Flux.params(app), gs)

Base.copyto!(dest::NeuralNetworkApproximator, src::NeuralNetworkApproximator) =
    Flux.loadparams!(dest.model, Flux.params(src))

#####
# NamedPolicy (old RLCore)
#####

# no P<:AbstractPolicy constraint: ElectricGrid wraps `nothing` when a setup
# has no classically controlled sources (multi_controller.jl SetupAgents)
struct NamedPolicy{P,N} <: AbstractPolicy
    name::N
    policy::P
end

Base.nameof(p::NamedPolicy) = p.name

# explicit signature (not Vararg) to avoid ambiguity with the generic no-ops
update!(p::NamedPolicy, t::AbstractTrajectory, e::AbstractEnv, s::AbstractStage) =
    update!(p.policy, t, e, s)
update!(p::NamedPolicy, x) = update!(p.policy, x)

(p::NamedPolicy)(env::AbstractEnv) = p.policy(env)
(p::NamedPolicy)(stage::AbstractStage, env::AbstractEnv) = p.policy(stage, env)
(p::NamedPolicy)(stage::PreActStage, env::AbstractEnv, action) =
    p.policy(stage, env, action)

#####
# RandomPolicy (old RLCore)
#####

struct RandomPolicy{S,R<:AbstractRNG} <: AbstractPolicy
    action_space::S
    rng::R
end

RandomPolicy(s; rng = Random.GLOBAL_RNG) = RandomPolicy(s, rng)

(p::RandomPolicy)(env::AbstractEnv) = _rand_in(p.rng, p.action_space)

#####
# Agent (old RLCore)
#####

Base.@kwdef struct Agent{P<:AbstractPolicy,T<:AbstractTrajectory} <: AbstractPolicy
    policy::P
    trajectory::T
end

Flux.functor(x::Agent) = (policy = x.policy,), y -> Agent(y.policy, x.trajectory)

(agent::Agent)(env::AbstractEnv) = agent.policy(env)

function (agent::Agent)(stage::AbstractStage, env::AbstractEnv)
    update!(agent.trajectory, agent.policy, env, stage)
    update!(agent.policy, agent.trajectory, env, stage)
end

function (agent::Agent)(stage::PreActStage, env::AbstractEnv, action)
    update!(agent.trajectory, agent.policy, env, stage, action)
    update!(agent.policy, agent.trajectory, env, stage)
end

# Default updates: no-ops, as in old RLCore (policies like CustomDDPGPolicy and
# the PreEpisode/PreAct trajectory methods below override the relevant cases;
# ElectricGrid overrides PostAct/PostEpisode for ElectricGridEnv).
update!(::AbstractPolicy, ::AbstractTrajectory, ::AbstractEnv, ::AbstractStage) = nothing
update!(::AbstractTrajectory, ::AbstractPolicy, ::AbstractEnv, ::AbstractStage) = nothing

function update!(
        trajectory::AbstractTrajectory,
        ::AbstractPolicy,
        ::AbstractEnv,
        ::PreEpisodeStage,
    )
    if length(trajectory) > 0
        pop!(trajectory[:state])
        pop!(trajectory[:action])
        if haskey(trajectory, :legal_actions_mask)
            pop!(trajectory[:legal_actions_mask])
        end
    end
end

function update!(
        trajectory::AbstractTrajectory,
        policy::AbstractPolicy,
        env::AbstractEnv,
        ::PreActStage,
        action,
    )
    s = policy isa NamedPolicy ? state(env, nameof(policy)) : state(env)
    push!(trajectory[:state], s)
    push!(trajectory[:action], action)
    if haskey(trajectory, :legal_actions_mask)
        lasm =
            policy isa NamedPolicy ? legal_action_space_mask(env, nameof(policy)) :
            legal_action_space_mask(env)
        push!(trajectory[:legal_actions_mask], lasm)
    end
end

# NOTE: the generic PostActStage/PostEpisodeStage trajectory updates are not
# vendored; ElectricGrid overrides both for `env::ElectricGridEnv`
# (see multi_controller.jl and agent_ddpg.jl).

#####
# Stop conditions (old RLCore, without progress meters)
#####

mutable struct StopAfterEpisode
    episode::Int
    cur::Int
end

StopAfterEpisode(episode; cur = 0, is_show_progress = false) =
    StopAfterEpisode(episode, cur)

function (s::StopAfterEpisode)(agent, env)
    is_terminated(env) && (s.cur += 1)
    s.cur >= s.episode
end

mutable struct StopAfterStep
    step::Int
    cur::Int
end

StopAfterStep(step; cur = 1, is_show_progress = false) = StopAfterStep(step, cur)

function (s::StopAfterStep)(args...)
    res = s.cur >= s.step
    s.cur += 1
    res
end

#####
# run loop (old RLCore _run; extends Base.run so the unqualified `run(...)`
# calls in render.jl and `RLBase.run(...)` in the tests keep working)
#####

function Base.run(
        policy::AbstractPolicy,
        env::AbstractEnv,
        stop_condition = StopAfterEpisode(1),
        hook = EmptyHook(),
    )
    hook(PRE_EXPERIMENT_STAGE, policy, env)
    policy(PRE_EXPERIMENT_STAGE, env)
    is_stop = false
    while !is_stop
        reset!(env)
        policy(PRE_EPISODE_STAGE, env)
        hook(PRE_EPISODE_STAGE, policy, env)

        while !is_terminated(env) # one episode
            action = policy(env)

            policy(PRE_ACT_STAGE, env, action)
            hook(PRE_ACT_STAGE, policy, env, action)

            env(action)

            policy(POST_ACT_STAGE, env)
            hook(POST_ACT_STAGE, policy, env)

            if stop_condition(policy, env)
                is_stop = true
                break
            end
        end # end of an episode

        if is_terminated(env)
            policy(POST_EPISODE_STAGE, env)  # let the policy see the last observation
            hook(POST_EPISODE_STAGE, policy, env)
        end
    end
    hook(POST_EXPERIMENT_STAGE, policy, env)
    hook
end

end # module RLCompat
