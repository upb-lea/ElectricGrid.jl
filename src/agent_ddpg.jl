using Random
using Setfield: @set
using Zygote: ignore
import CUDA: device
using MacroTools: @forward

#functions to make the training flag work for the whole RL framework
(agent::Agent)(env::ElectricGridEnv, training::Bool) = agent.policy(env, training)
(p::NamedPolicy)(env::ElectricGridEnv, training::Bool) = p.policy(env, p.name, training)
(policy::AbstractPolicy)(env, training::Bool) = policy(env)
(policy::AbstractPolicy)(env, name::Union{Nothing, String}, training::Bool) = policy(env, name)
(policy::AbstractPolicy)(stage::AbstractStage, env::AbstractEnv, training::Bool) = policy(stage, env)
(policy::AbstractPolicy)(stage::AbstractStage, env::AbstractEnv, action, training::Bool) = policy(stage, env, action)


#re-definition needed to use ElectricGrid-internal action-space functions
function RLBase.update!(
        trajectory::AbstractTrajectory,
        policy::AbstractPolicy,
        env::ElectricGridEnv,
        ::PostEpisodeStage,
    )
    # Note that for trajectories like `CircularArraySARTTrajectory`, data are
    # stored in a SARSA format, which means we still need to generate a dummy
    # action at the end of an episode.

    #s = policy isa NamedPolicy ? state(env, nameof(policy)) : state(env)
    #A = policy isa NamedPolicy ? action_space(env, nameof(policy)) : action_space(env)
    #a = rand(A)

    #just copy the last values to achieve that - otherwise the osc module broke the commented out code from above

    push!(trajectory[:state], trajectory[:state][:,end])
    push!(trajectory[:action], trajectory[:action][:,end])
    if haskey(trajectory, :legal_actions_mask)
        lasm =
            policy isa NamedPolicy ? legal_action_space_mask(env, nameof(policy)) :
            legal_action_space_mask(env)
        push!(trajectory[:legal_actions_mask], lasm)
    end
end


Base.@kwdef struct CustomNeuralNetworkApproximator{M,O} <: AbstractApproximator
    model::M
    optimizer::O = nothing
end

# some model may accept multiple inputs
(app::CustomNeuralNetworkApproximator)(args...; kwargs...) = app.model(args...; kwargs...)

@forward CustomNeuralNetworkApproximator.model Flux.testmode!,
Flux.trainmode!,
Flux.params,
device

functor(x::CustomNeuralNetworkApproximator) =
    (model=x.model,), y -> CustomNeuralNetworkApproximator(y.model, x.optimizer)

RLBase.update!(app::CustomNeuralNetworkApproximator, gs) =
    Flux.Optimise.update!(app.optimizer, Flux.params(app), gs)

Base.copyto!(dest::CustomNeuralNetworkApproximator, src::CustomNeuralNetworkApproximator) =
    Flux.loadparams!(dest.model, Flux.params(src))

mutable struct CustomDDPGPolicy{
        BA<:CustomNeuralNetworkApproximator,
        BC<:CustomNeuralNetworkApproximator,
        TA<:CustomNeuralNetworkApproximator,
        TC<:CustomNeuralNetworkApproximator,
        P,
        R<:AbstractRNG,
    } <: AbstractPolicy

    behavior_actor::BA
    behavior_critic::BC
    target_actor::TA
    target_critic::TC
    γ::Float32
    ρ::Float32
    na::Int
    batch_size::Int
    start_steps::Int
    start_policy::P
    update_after::Int
    update_freq::Int
    act_limit::Float64
    act_noise::Float64
    update_step::Int
    rng::R
    # for logging
    actor_loss::Float32
    critic_loss::Float32
end

Flux.functor(x::CustomDDPGPolicy) = (
    ba = x.behavior_actor,
    bc = x.behavior_critic,
    ta = x.target_actor,
    tc = x.target_critic,
),
y -> begin
    x = @set x.behavior_actor = y.ba
    x = @set x.behavior_critic = y.bc
    x = @set x.target_actor = y.ta
    x = @set x.target_critic = y.tc
    x
end

function CustomDDPGPolicy(;
        behavior_actor,
        behavior_critic,
        target_actor,
        target_critic,
        start_policy,
        γ = 0.99f0,
        ρ = 0.995f0,
        na = 1,
        batch_size = 32,
        start_steps = 10000,
        update_after = 1000,
        update_freq = 50,
        act_limit = 1.0,
        act_noise = 0.1,
        update_step = 0,
        rng = Random.GLOBAL_RNG,
    )
    copyto!(behavior_actor, target_actor)  # force sync
    copyto!(behavior_critic, target_critic)  # force sync
    CustomDDPGPolicy(
        behavior_actor,
        behavior_critic,
        target_actor,
        target_critic,
        γ,
        ρ,
        na,
        batch_size,
        start_steps,
        start_policy,
        update_after,
        update_freq,
        act_limit,
        act_noise,
        update_step,
        rng,
        0.0f0,
        0.0f0,
    )
end

function (::CustomDDPGPolicy)(::AbstractStage, ::AbstractEnv)
    nothing
end

function (p::CustomDDPGPolicy)(env::ElectricGridEnv, name::Union{Nothing, String} = nothing, training::Bool = false)
    p.update_step += 1

    if p.update_step <= p.start_steps
        p.start_policy(env)
    else
        D = device(p.behavior_actor)
        s = isnothing(name) ? state(env) : state(env, name)
        s = Flux.unsqueeze(s, ndims(s) + 1)
        actions = p.behavior_actor(send_to_device(D, s)) |> vec |> send_to_host
        if training
            c = clamp.(actions .+ randn(p.rng, p.na) .* repeat([p.act_noise], p.na), -p.act_limit, p.act_limit)
        else
            c = clamp.(actions, -p.act_limit, p.act_limit)
        end
        #p.na == 1 && return c[1]
        c
    end
end

function RLBase.update!(
        p::CustomDDPGPolicy,
        traj::CircularArraySARTTrajectory,
        ::AbstractEnv,
        ::PreActStage,
    )
    length(traj) > p.update_after || return
    p.update_step % p.update_freq == 0 || return
    inds, batch = sample(p.rng, traj, BatchSampler{SARTS}(p.batch_size))
    RLBase.update!(p, batch)
end

function RLBase.update!(p::CustomDDPGPolicy, batch::NamedTuple{SARTS})
    s, a, r, t, s′ = send_to_device(device(p), batch)

    A = p.behavior_actor
    C = p.behavior_critic
    Aₜ = p.target_actor
    Cₜ = p.target_critic

    γ = p.γ
    ρ = p.ρ


    # !!! we have several assumptions here, need revisit when we have more complex environments
    # state is vector
    # action is scalar
    a′ = Aₜ(s′)
    qₜ = Cₜ(vcat(s′, a′)) |> vec
    y = r .+ γ .* (1 .- t) .* qₜ
    a = Flux.unsqueeze(a, ndims(a)+1)

    gs1 = gradient(Flux.params(C)) do
        q = C(vcat(s, a)) |> vec
        loss = mean((y .- q) .^ 2)
        ignore() do
            p.critic_loss = loss
        end
        loss
    end

    RLBase.update!(C, gs1)

    gs2 = gradient(Flux.params(A)) do
        loss = -mean(C(vcat(s, A(s))))
        ignore() do
            p.actor_loss = loss
        end
        loss
    end

    RLBase.update!(A, gs2)

    # polyak averaging
    for (dest, src) in zip(Flux.params([Aₜ, Cₜ]), Flux.params([A, C]))
        dest .= ρ .* dest .+ (1 - ρ) .* src
    end
end


# also in a sep src

global rngg = StableRNG(123)
global initt = Flux.glorot_uniform(rngg)

global CreateActor(na, ns) = Chain(
    Dense(ns, 20, relu; init = initt),
    Dense(20, 10, relu; init = initt),
    Dense(10, na, tanh; init = initt),
)

global CreateCritic(na, ns) = Chain(
    Dense(ns + na, 20, relu; init = initt),
    Dense(20, 10, relu; init = initt),
    Dense(10, 1; init = initt),
)

function CreateAgentDdpg(;na, ns, batch_size = 32, use_gpu = true)
    Agent(
        policy = CustomDDPGPolicy(
            behavior_actor = CustomNeuralNetworkApproximator(
                model = use_gpu ? CreateActor(na, ns) |> gpu : CreateActor(na, ns),
                optimizer = Flux.ADAM(),
            ),
            behavior_critic = CustomNeuralNetworkApproximator(
                model = use_gpu ? CreateCritic(na, ns) |> gpu : CreateCritic(na, ns),
                optimizer = Flux.ADAM(),
            ),
            target_actor = CustomNeuralNetworkApproximator(
                model = use_gpu ? CreateActor(na, ns) |> gpu : CreateActor(na, ns),
                optimizer = Flux.ADAM(),
            ),
            target_critic = CustomNeuralNetworkApproximator(
                model = use_gpu ? CreateCritic(na, ns) |> gpu : CreateCritic(na, ns),
                optimizer = Flux.ADAM(),
            ),
            γ = 0.999f0,
            ρ = 0.895f0,
            na = na,
            batch_size = batch_size,
            start_steps = -1,
            start_policy = RandomPolicy(-1.0..1.0; rng = rngg),
            update_after = 50, #1000
            update_freq = 10,
            act_limit = 1.0,
            act_noise = 0.032,
            rng = rngg,
        ),
        trajectory = CircularArraySARTTrajectory(
            capacity = 20_000,
            state = Vector{Float32} => (ns,),
            action = Float32 => (na, ),
        ),
    )
end
