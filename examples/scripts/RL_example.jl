# In this script we set up an experiment with one source that is controlled by RL that is connected to one load. Note that this is a single phase setup. The goal of the RL agent will be to control the input voltage of the source such that the corresponding current will match a reference signal which is a constant of 10 A.

using ElectricGrid
using PlotlyJS

R_load, L_load, X, Z = ParallelLoadImpedance(100e3, 1, 230)

# define grid using CM
CM = [0. 1.
    -1. 0.]

# Set parameters accoring graphic above
parameters = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("pwr" => 200e3, "control_type" => "RL", "mode" => "my_ddpg", "fltr" => "L"),
                    ],
    "load"   => Any[
                    Dict{Any, Any}("impedance" => "R", "R" => R_load,"v_limit"=>1e4, "i_limit"=>1e4),
                    ],
    "grid" => Dict{Any, Any}("phase" => 1)
)


# The reference function just returns the 10.0 A and the featurize_ddpg adds it to the controller's state (normalized).

function reference(t)
    return 10.0
end


featurize_ddpg = function(state, env, name)
    if name == "my_ddpg"
        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        state = vcat(state, reference(env.t)/norm_ref)
    end
end


# The reward function computes the error and takes the square root.

function reward_function(env, name = nothing)
    index_1 = findfirst(x -> x == "source1_i_L1", env.state_ids)
    state_to_control = env.state[index_1]

    if any(abs.(state_to_control).>1)
        return -1
    else

            refs = reference(env.t)
            norm_ref = env.nc.parameters["source"][1]["i_limit"]
            r = 1-((abs(refs/norm_ref - state_to_control)/2)^0.5)
            return r
        end
    end
end

env = ElectricGridEnv(
    CM = CM,
    parameters = parameters,
    t_end = 0.1,
    featurize = featurize_ddpg,
    reward_function = reward_function,
    action_delay = 0)

#agent = CreateAgentDdpg(na = length(env.agent_dict["my_ddpg"]["action_ids"]), ns = length(state(env, "my_ddpg")), use_gpu = false)


rng = StableRNG(1)
init = glorot_uniform(rng)

ns = length(state(env, "my_ddpg"))#length(env.agent_dict["my_ddpg"]["state_ids"])
na = length(env.agent_dict["my_ddpg"]["action_ids"])

CreateActor() = Chain(
    Dense(ns, 30, relu; init = init),
    Dense(30, 30, relu; init = init),
    Dense(30, 1, tanh; init = init),
) |> gpu

CreateCritic() = Chain(
    Dense(ns + na, 30, relu; init = init),
    Dense(30, 30, relu; init = init),
    Dense(30, 1; init = init),
) |> gpu

agent = Agent(
    policy = DDPGPolicy(
        behavior_actor = NeuralNetworkApproximator(
            model = CreateActor(),
            optimizer = ADAM(),
        ),
        behavior_critic = NeuralNetworkApproximator(
            model = CreateCritic(),
            optimizer = ADAM(),
        ),
        target_actor = NeuralNetworkApproximator(
            model = CreateActor(),
            optimizer = ADAM(),
        ),
        target_critic = NeuralNetworkApproximator(
            model = CreateCritic(),
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

controllers = SetupAgents(env, Dict("my_ddpg" => agent))

#run(ma["ElectricGrid_ddpg_1"]["policy"], env)

Learn(controllers, env, num_episodes = 80)

#Learn(agent, env)


states_to_plot = ["source1_i_L1"]
actions_to_plot = ["source1_u"]

hook = DataHook(collect_state_ids = states_to_plot, collect_action_ids = actions_to_plot)

Simulate(controllers, env, hook=hook)

RenderHookResults(hook = hook,
                  states_to_plot  = states_to_plot,
                  actions_to_plot = actions_to_plot)
