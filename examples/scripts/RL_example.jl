using Dare
using ReinforcementLearning
using StableRNGs
using Flux
using Flux.Losses
using IntervalSets

R_load, L_load, X, Z = Parallel_Load_Impedance(100e3, 1, 230)

# define grid using CM
CM = [0. 1.
    -1. 0.]

# Set parameters accoring graphic above
parameters = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("pwr" => 200e3, "control_type" => "RL", "mode" => "user_def", "fltr" => "LC"),
                    ],
    "load"   => Any[
                    Dict{Any, Any}("impedance" => "R", "R" => R_load,"v_limit"=>1e4, "i_limit"=>1e4),
                    ],
    "grid" => Dict{Any, Any}("phase" => 3)
)


function reference(t)
    return [1]
end


featurize_ddpg = function(state, env, name)

    norm_ref = env.nc.parameters["source"][1]["i_limit"]
    state = vcat(state, reference(env.t)/norm_ref)

end

env = SimEnv(
    CM = CM,
    parameters = parameters,
    t_end = 0.1,
    featurize = featurize_ddpg,
    #reward_function = reward_function,
    action_delay = 0)

#agent = create_agent_ddpg(na = length(env.agent_dict["user_def"]["action_ids"]), ns = length(state(env, "user_def")), use_gpu = false)


rng = StableRNG(1)
init = glorot_uniform(rng)

ns = length(env.agent_dict["user_def"]["state_ids"])
na = length(env.agent_dict["user_def"]["action_ids"])

create_actor() = Chain(
    Dense(ns, 30, relu; init = init),
    Dense(30, 30, relu; init = init),
    Dense(30, 1, tanh; init = init),
) |> gpu

create_critic() = Chain(
    Dense(ns + na, 30, relu; init = init),
    Dense(30, 30, relu; init = init),
    Dense(30, 1; init = init),
) |> gpu

agent = Agent(
    policy = DDPGPolicy(
        behavior_actor = NeuralNetworkApproximator(
            model = create_actor(),
            optimizer = ADAM(),
        ),
        behavior_critic = NeuralNetworkApproximator(
            model = create_critic(),
            optimizer = ADAM(),
        ),
        target_actor = NeuralNetworkApproximator(
            model = create_actor(),
            optimizer = ADAM(),
        ),
        target_critic = NeuralNetworkApproximator(
            model = create_critic(),
            optimizer = ADAM(),
        ),
        γ = 0.99f0,
        ρ = 0.995f0,
        na = 1,
        batch_size = 64,
        start_steps = 1000,
        start_policy = RandomPolicy(-1.0..1.0; rng = rng),
        update_after = 1000,
        update_freq = 1,
        act_limit = 1.0,
        act_noise = 0.1,
        rng = rng,
    ),
    trajectory = CircularArraySARTTrajectory(
        capacity = 10000,
        state = Vector{Float32} => (ns,),
        action = Float32 => (na, ),
    ),
)

controllers = setup_agents(env, Dict("user_def" => agent))


#run(ma["dare_ddpg_1"]["policy"], env)

learn(controllers, env, num_episodes = 1)

#learn(agent, env)
