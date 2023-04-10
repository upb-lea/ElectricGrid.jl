using JEG
using ReinforcementLearning

"""
This scipt contains the content of the RL_Complex_DEMO.ipynb notebook.
For comments and more documentation see refering notebookin examples/notebooks
"""
CM = [0.0   1.0  0
     -1.0   0.0  2.0
     0  -2.0  0.0]

parameters =
Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "control_type" => "RL",
                        "mode" => "my_ddpg",
                        "fltr" => "L"),
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "fltr" => "LC",
                        "control_type" =>
                        "RL", "mode" => "my_ddpg"),
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "fltr" => "L",
                        "control_type" =>
                        "RL", "mode" => "my_ddpg"),
                    ],
    "grid" => Dict{Any, Any}(
        "phase" => 1,
        "ramp_end" => 0.04,)
)

function reference(t)

    return [-10, 230, -15]
end

featurize_ddpg = function(state, env, name)
    if name == "my_ddpg"
        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        state = vcat(state, reference(env.t)/norm_ref)
    end
end

function reward_function(env, name = nothing)
    if name == "my_ddpg"
        state_to_control_1 = env.state[findfirst(x -> x == "source1_i_L1", env.state_ids)]
        state_to_control_2 = env.state[findfirst(x -> x == "source2_v_C_filt", env.state_ids)]
        state_to_control_3 = env.state[findfirst(x -> x == "source3_i_L1", env.state_ids)]

        state_to_control = [state_to_control_1, state_to_control_2, state_to_control_3]

        if any(abs.(state_to_control).>1)
            return -1
        else

            refs = reference(env.t)
            norm_ref = env.nc.parameters["source"][1]["i_limit"]
            r = 1-1/3*(sum((abs.(refs/norm_ref - state_to_control)/2).^0.5))
            return r
        end
    else
        return 1
    end

end

env = ElectricGridEnv(
    CM = CM,
    parameters = parameters,
    t_end = 1,
    featurize = featurize_ddpg,
    reward_function = reward_function,
    action_delay = 0)

    agent = CreateAgentDdpg(na = length(env.agent_dict["my_ddpg"]["action_ids"]),
    ns = length(state(env, "my_ddpg")),
    use_gpu = false)

my_custom_agents = Dict("my_ddpg" => agent)

controllers = SetupAgents(env, my_custom_agents)


agent = CreateAgentDdpg(na = length(env.agent_dict["my_ddpg"]["action_ids"]),
                          ns = length(state(env, "my_ddpg")),
                          use_gpu = false)

my_custom_agents = Dict("my_ddpg" => agent)

controllers = SetupAgents(env, my_custom_agents)

Learn(controllers, env, num_episodes = 20)

hook = DataHook(collect_state_ids = env.state_ids,
                collect_action_ids = env.action_ids)

hook = Simulate(controllers, env, hook=hook)


RenderHookResults(hook = hook,
                    states_to_plot  = env.state_ids,
                    actions_to_plot = env.action_ids,
                    plot_reward=true)
