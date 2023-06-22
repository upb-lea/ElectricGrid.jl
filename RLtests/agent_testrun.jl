using ElectricGrid

CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]

R_load, L_load, _, _ = ParallelLoadImpedance(50e3, 0.95, 230)

parameters =
Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "control_type" => "RL",
                        "mode" => "my_ddpg",
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
    if name == "ElectricGrid_ddpg_1"

        #state = state[findall(x -> split(x, "_")[2] == "i" , env.agent_dict["ElectricGrid_ddpg_1"]["state_ids"])]

        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        state = vcat(state, reference(env.t)/norm_ref)

        # θ = 2*pi*50*env.t

        # state_to_control_1 = env.state[findfirst(x -> x == "source1_i_L1_a", env.state_ids)]
        # state_to_control_2 = env.state[findfirst(x -> x == "source1_i_L1_b", env.state_ids)]
        # state_to_control_3 = env.state[findfirst(x -> x == "source1_i_L1_c", env.state_ids)]

        # state_to_control = [state_to_control_1, state_to_control_2, state_to_control_3]

        # Il_dq0 = DQ0Transform(state_to_control, θ)
        # state = vcat(state, Il_dq0)
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


agent = CreateAgentDdpg(na = length(env.agent_dict["my_ddpg"]["action_ids"]),
    ns = length(state(env, "my_ddpg")),
    use_gpu = false)

my_custom_agents = Dict("my_ddpg" => agent)

controllers = SetupAgents(env, my_custom_agents)

learnhook = DataHook()

function learn()
    steps_total = 1_500_000
    steps_loop = 50_000

    while length(controllers.hook.df[!,"reward"]) <= steps_total

        println("Steps so far: $(length(controllers.hook.df[!,"reward"]))")
        Learn(controllers, env, steps = steps_loop, hook = learnhook)

    end
end


learn()

plot_rewardresults()

hook = DataHook(collect_state_ids = env.state_ids,
                collect_action_ids = env.action_ids)

Simulate(controllers, env, hook=hook)


RenderHookResults(hook = hook,
                    states_to_plot  = env.state_ids,
                    actions_to_plot = env.action_ids,
                    plot_reward=true)
