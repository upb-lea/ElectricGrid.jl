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
        r = 1-((abs.(refs/norm_ref - state_to_control)/2).^0.5)
        return r 
    end
end


# Now we can set up the env and the controllers and define a training functions.

env = ElectricGridEnv(
    CM = CM, 
    parameters = parameters, 
    t_end = 0.025, 
    featurize = featurize_ddpg, 
    reward_function = reward_function, 
    action_delay = 0);


agent = CreateAgentDdpg(na = length(env.agent_dict["my_ddpg"]["action_ids"]),
                          ns = length(state(env, "my_ddpg")),
                          use_gpu = false);

controllers = SetupAgents(env, Dict("my_ddpg" => agent))


# We start training with 30_000 steps with the default action noise (0.032).

learnhook = DataHook()
Learn(controllers, env, steps = 30_000, hook = learnhook);


# Here we define a second training function with an action noise scheduler enabled.

function learn2()
    num_steps = 10_000

    an_scheduler_loops = 7


    for j in 1:1
        an = 0.001 * exp10.(collect(LinRange(0.0, -13, an_scheduler_loops)))
        for i in 1:an_scheduler_loops
            controllers.agents["my_ddpg"]["policy"].policy.policy.act_noise = an[i]
            println("Steps so far: $(length(controllers.hook.df[!,"reward"]))")
            println("next action noise level: $(an[i])")
            Learn(controllers, env, steps = num_steps, hook = learnhook)
        end
    end
end


# Run it!

learn2()


# Run an additional 20_000 steps of training with action noise set to 0.0.

controllers.agents["my_ddpg"]["policy"].policy.policy.act_noise = 0.0
Learn(controllers, env, steps = 20_000, hook = learnhook);


# We now plot the learning rate per time step.

p = plot(learnhook.df[!, :reward])
display(p)


# Now we also want to run a simulation with the fully trained agent and plot the results.

states_to_plot = ["source1_i_L1"]
actions_to_plot = ["source1_u"]

hook = DataHook(collect_state_ids = states_to_plot, collect_action_ids = actions_to_plot)

Simulate(controllers, env, hook=hook)

RenderHookResults(hook = hook,
                  states_to_plot  = states_to_plot,
                  actions_to_plot = actions_to_plot)
