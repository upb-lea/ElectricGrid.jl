Base.@kwdef struct DareNeuralNetworkApproximator{M,O} <: AbstractApproximator
    model::M
    optimizer::O = nothing
end

# some model may accept multiple inputs
(app::DareNeuralNetworkApproximator)(args...; kwargs...) = app.model(args...; kwargs...)

@forward DareNeuralNetworkApproximator.model Flux.testmode!,
Flux.trainmode!,
Flux.params,
device

functor(x::DareNeuralNetworkApproximator) =
    (model=x.model,), y -> DareNeuralNetworkApproximator(y.model, x.optimizer)

RLBase.update!(app::DareNeuralNetworkApproximator, gs) =
    Flux.Optimise.update!(app.optimizer, Flux.params(app), gs)

Base.copyto!(dest::DareNeuralNetworkApproximator, src::DareNeuralNetworkApproximator) =
    Flux.loadparams!(dest.model, Flux.params(src))



# also in a sep src

global rngg = StableRNG(123)
global initt = Flux.glorot_uniform(rngg)

global create_actor(na, ns) = Chain(
    Dense(ns, 40, relu; init = initt),
    Dense(40, 30, relu; init = initt),
    Dense(30, na, tanh; init = initt),
)

global create_critic(na, ns) = Chain(
    Dense(ns + na, 100, relu; init = initt),
    Dense(100, 100, relu; init = initt),
    Dense(100, 100, relu; init = initt),
    Dense(100, 1; init = initt),
)

function create_agent_ddpg(;na, ns, batch_size = 32, use_gpu = true)
    Agent(
        policy = DDPGPolicy(
            behavior_actor = DareNeuralNetworkApproximator(
                model = use_gpu ? create_actor(na, ns) |> gpu : create_actor(na, ns),
                optimizer = Flux.ADAM(),
            ),
            behavior_critic = DareNeuralNetworkApproximator(
                model = use_gpu ? create_critic(na, ns) |> gpu : create_critic(na, ns),
                optimizer = Flux.ADAM(),
            ),
            target_actor = DareNeuralNetworkApproximator(
                model = use_gpu ? create_actor(na, ns) |> gpu : create_actor(na, ns),
                optimizer = Flux.ADAM(),
            ),
            target_critic = DareNeuralNetworkApproximator(
                model = use_gpu ? create_critic(na, ns) |> gpu : create_critic(na, ns),
                optimizer = Flux.ADAM(),
            ),
            γ = 0.99f0,
            ρ = 0.995f0,
            na = na,
            batch_size = batch_size,
            start_steps = 0,
            start_policy = RandomPolicy(-1.0..1.0; rng = rngg),
            update_after = 50, #1000 
            update_freq = 10,
            act_limit = 1.0,
            act_noise = 0.1,
            rng = rngg,
        ),
        trajectory = CircularArraySARTTrajectory(
            capacity = 800,
            state = Vector{Float32} => (ns,),
            action = Float32 => (na, ),
        ),
    )
end

function (::DDPGPolicy)(::AbstractStage, ::AbstractEnv)
    nothing
end

function (p::DDPGPolicy)(env::SimEnv, name::Any = nothing)
    p.update_step += 1

    if p.update_step <= p.start_steps
        p.start_policy(env)
    else
        D = device(p.behavior_actor)
        s = isnothing(name) ? state(env) : state(env, name)
        s = Flux.unsqueeze(s, ndims(s) + 1)
        actions = p.behavior_actor(send_to_device(D, s)) |> vec |> send_to_host
        c = clamp.(actions .+ randn(p.rng, p.na) .* repeat([p.act_noise], p.na), -p.act_limit, p.act_limit)
        p.na == 1 && return c[1]
        c
    end
end