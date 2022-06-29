using ReinforcementLearning
using Flux
using StableRNGs
using IntervalSets

# also in a sep src
global rngg = StableRNG(123)
global initt = Flux.glorot_uniform(rngg)

global create_actor(na, ns) = Chain(
    Dense(ns, 40, relu; init = initt),
    Dense(40, 30, relu; init = initt),
    Dense(30, na, tanh; init = initt),
)

global create_critic(na, ns) = Chain(
    Dense(ns + na, 40, relu; init = initt),
    Dense(40, 30, relu; init = initt),
    Dense(30, 1; init = initt),
)

function create_agent(na, ns, use_gpu = true)
    Agent(
        policy = DDPGPolicy(
            behavior_actor = NeuralNetworkApproximator(
                model = use_gpu ? create_actor(na, ns) |> gpu : create_actor(na, ns),
                optimizer = Flux.ADAM(),
            ),
            behavior_critic = NeuralNetworkApproximator(
                model = use_gpu ? create_critic(na, ns) |> gpu : create_critic(na, ns),
                optimizer = Flux.ADAM(),
            ),
            target_actor = NeuralNetworkApproximator(
                model = use_gpu ? create_actor(na, ns) |> gpu : create_actor(na, ns),
                optimizer = Flux.ADAM(),
            ),
            target_critic = NeuralNetworkApproximator(
                model = use_gpu ? create_critic(na, ns) |> gpu : create_critic(na, ns),
                optimizer = Flux.ADAM(),
            ),
            γ = 0.99f0,
            ρ = 0.995f0,
            na = na,
            batch_size = 64,
            start_steps = 0,
            start_policy = RandomPolicy(-1.0..1.0; rng = rngg),
            update_after = 150, #1000 
            update_freq = 1,
            act_limit = 1.0,
            act_noise = 0.1,
            rng = rngg,
        ),
        trajectory = CircularArraySARTTrajectory(
            capacity = 10000,
            state = Vector{Float32} => (ns,),
            action = Float32 => (na, ),
        ),
    )
end