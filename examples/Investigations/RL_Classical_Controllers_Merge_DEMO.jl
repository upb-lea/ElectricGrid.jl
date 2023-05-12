using ElectricGrid
#using ReinforcementLearning

"""
This scipt contains the content of the RL_Classical_Controllers_Merge_DEMO.ipynb notebook.
For comments and more documentation see refering notebookin examples/notebooks
"""

println("...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]

R_load, L_load, _, _ = ParallelLoadImpedance(50e3, 0.99, 230)

parameters =
Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "control_type" => "RL",
                        "fltr" => "LCL",),
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "fltr" => "LCL",
                        "control_type" => "classic",
                        "mode" => "Swing",),
                    ],
    "load"   => Any[
                    Dict{Any, Any}("impedance" => "RL",
                                    "R" => R_load,
                                    "L" => L_load,
                                    "v_limit" => 10e3),
                                    ],
    "cable"   => Any[
                    Dict{Any, Any}("R" => 0.1,
                                    "L" => 0.0025e-3,
                                    "C" => 0.001e-4),
                    Dict{Any, Any}("R" => 0.1,
                                    "L" => 0.0025e-3,
                                    "C" => 0.001e-4),
                    ],
    "grid" => Dict{Any, Any}(
        "phase" => 3,
        "ramp_end" => 0.04,)
)



function reward_function(env, name = nothing)

    if name == "my_ddpg"

        p_ref = 0e3
        q_ref = 0e3

        pq0_ref = [p_ref; q_ref; 0]

        v_a_meas = env.x[findfirst(x -> x == "source1_v_C_cables_a", env.state_ids)]
        v_b_meas = env.x[findfirst(x -> x == "source1_v_C_cables_b", env.state_ids)]
        v_c_meas  = env.x[findfirst(x -> x == "source1_v_C_cables_c", env.state_ids)]
        i_a_meas = env.state[findfirst(x -> x == "source1_i_L2_a", env.state_ids)]
        i_b_meas = env.state[findfirst(x -> x == "source1_i_L2_b", env.state_ids)]
        i_c_meas  = env.state[findfirst(x -> x == "source1_i_L2_c", env.state_ids)]

        v_abc_meas = [v_a_meas; v_b_meas; v_c_meas]
        i_abc_meas = [i_a_meas; i_b_meas; i_c_meas]

        i_abc_ref = InvClarkeTransform(Inv_p_q_v(ClarkeTransform(v_abc_meas), pq0_ref))

        i_abc_ref = clamp.(i_abc_ref/env.nc.parameters["source"][1]["i_limit"], -0.95, 0.95)

        i_abc_ref= [0; 0; 0]

        #println(i_abc_ref)

        if any((abs.(env.state)).>=1)
            return -1
        end

        r = 1 - 1/3*(sum((abs.(i_abc_ref - i_abc_meas)/2).^0.5))

        return r * (1 - agent.policy.γ)/2

    else
        return 0
    end

end

env = ElectricGridEnv(
    CM =  CM,
    parameters = parameters,
    t_end = 1,
    reward_function = reward_function,
    action_delay = 0,
    verbosity = 0)


controllers = SetupAgents(env)

function learn()
    while true
        Learn(controllers, env, num_episodes = 20)
    end
end

hook = DataHook(collect_state_ids = env.state_ids,
                collect_action_ids = env.action_ids)

Simulate(controllers, env, hook=hook)


RenderHookResults(hook = hook,
                    states_to_plot  = env.state_ids,
                    actions_to_plot = env.action_ids,
                    plot_reward=true)

println("...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")
