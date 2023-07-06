using ElectricGrid
using ReinforcementLearning
using Flux
using Flux.Losses
using StableRNGs
using IntervalSets

# include agetn_td3 here
src_dir = joinpath(dirname(pathof(ElectricGrid)))
include(src_dir * "/agent_td3.jl")

CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]

R_load, L_load, _, _ = ParallelLoadImpedance(50e3, 0.95, 230)

parameters = Dict{Any, Any}(
                    "source" => Any[
                                    Dict{Any, Any}(
                                        "pwr" => 200e3,
                                        "control_type" => "RL",
                                        "mode" => "my_agent",
                                        "fltr" => "L",
                                        #"L1" => 0.0008,
                                        ),
                                    Dict{Any, Any}(
                                        "pwr" => 200e3,
                                        "fltr" => "LC",
                                        "control_type" => "classic",
                                        "mode" => "Droop",),
                                    ],
                    "grid" => Dict{Any, Any}(
                        "phase" => 3,
                        "ramp_end" => 0.04,)
    )


function reference(t)
    if t < 0.04
        return [0.0, 0.0, 0.0]
    end

    θ = 2*pi*50*t
    θph = [θ; θ - 120π/180; θ + 120π/180]
    return +10 * cos.(θph) 
end

featurize_ddpg = function(state, env, name)
    if name == "my_agent"
        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        state = vcat(state, reference(env.t)/norm_ref)
    end
end

function reward_function(env, name = nothing)
    if name == "classic"
        return 0        
    else
        state_to_control_1 = env.state[findfirst(x -> x == "source1_i_L1_a", env.state_ids)]
        state_to_control_2 = env.state[findfirst(x -> x == "source1_i_L1_b", env.state_ids)]
        state_to_control_3 = env.state[findfirst(x -> x == "source1_i_L1_c", env.state_ids)]

        state_to_control = [state_to_control_1, state_to_control_2, state_to_control_3]

        if any(abs.(state_to_control).>1)
            return -1
        else

            refs = reference(env.t)
            norm_ref = env.nc.parameters["source"][1]["i_limit"]          
            r = 1-1/3*(sum((abs.(refs/norm_ref - state_to_control)/2).^0.5))
            return r 
        end
    end

end

env = ElectricGridEnv(
    #CM =  CM,
    parameters = parameters,
    t_end = 1,
    reward_function = reward_function,
    featurize = featurize_ddpg,
    action_delay = 0,
    verbosity = 0)

seed = 123

rng = StableRNG(seed)
init = glorot_uniform(rng)
create_actor() = Chain(
    Dense(ns, 32, relu; init = init),
    Dense(32, 64, relu; init = init),
    Dense(64, na, tanh; init = init)
)

create_critic_model() = Chain(
    Dense(ns + na, 64, relu; init = init),
    Dense(64, 64, relu; init = init),
    Dense(64, 1; init = init)
)

struct TD3Critic
    critic_1::Flux.Chain
    critic_2::Flux.Chain
end


create_critic() = TD3Critic(
    create_critic_model(), 
    create_critic_model()
    )

ns = length(env.state)
na = length(env.agent_dict["my_agent"]["action_ids"])

agent = Agent(
    policy = TD3Policy(
        behavior_actor = NeuralNetworkApproximator(
            model = create_actor(),
            optimizer = ADAM(0.001),
        ),
        behavior_critic = NeuralNetworkApproximator(
            model = create_critic(),
            optimizer = ADAM(0.001),
        ),
        target_actor = NeuralNetworkApproximator(
            model = create_actor(),
            optimizer = ADAM(0.001),
        ),
        target_critic = NeuralNetworkApproximator(
            model = create_critic(),
            optimizer = ADAM(0.001),
        ),
        γ = 0.99f0,
        ρ = 0.99f0,
        batch_size = 64,
        start_steps = 1000,
        start_policy = RandomPolicy(-1.0..1.0; rng = rng),
        update_after = 1000,
        update_freq = 1,
        policy_freq = 2,
        target_act_limit = 1.0,
        target_act_noise = 0.1,
        act_limit = 1.0,
        act_noise = 0.1,
        rng = rng,
    ),

    trajectory = CircularArraySARTTrajectory(
            capacity = 10_000_000,
            state = Vector{Float32} => (ns,),
            action = Float32 => (),
    ),
)



run(agent, env, 1_000, verbosity = 1)



agent = CreateAgentDdpg(na = length(env.agent_dict["my_agent"]["action_ids"]),
    ns = length(state(env, "my_agent")),
    use_gpu = false)

my_custom_agents = Dict("my_agent" => agent)

controllers = SetupAgents(env, my_custom_agents)


function learn()
    steps_total = 1_500_000
    steps_loop = 50_000

    Learn(controllers, env, steps = steps_loop)
    while length(controllers.hook.df[!,"reward"]) <= steps_total

        println("Steps so far: $(length(controllers.hook.df[!,"reward"]))")
        Learn(controllers, env, steps = steps_loop, hook = learnhook)

    end
end


learn()

plot_rewardresults(controllers = controllers)

hook = DataHook(collect_state_ids = env.state_ids,
                collect_action_ids = env.action_ids)

Simulate(controllers, env, hook=hook)


RenderHookResults(hook = hook,
                    states_to_plot  = env.state_ids,
                    actions_to_plot = env.action_ids,
                    plot_reward=true)
