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
                        "fltr" => "L",),
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



function reward_function(env, name = nothing)

    if name == "classic"
        return 0
    else
        p_ref = 1e3
        q_ref = 0e3

        pq0_ref = [p_ref; q_ref; 0]

        v_a_meas = env.x[findfirst(x -> x == "source1_v_C_cables_a", env.state_ids)]
        v_b_meas = env.x[findfirst(x -> x == "source1_v_C_cables_b", env.state_ids)]
        v_c_meas  = env.x[findfirst(x -> x == "source1_v_C_cables_c", env.state_ids)]
        i_a_meas = env.state[findfirst(x -> x == "source1_i_L1_a", env.state_ids)]
        i_b_meas = env.state[findfirst(x -> x == "source1_i_L1_b", env.state_ids)]
        i_c_meas  = env.state[findfirst(x -> x == "source1_i_L1_c", env.state_ids)]

        v_abc_meas = [v_a_meas; v_b_meas; v_c_meas]
        i_abc_meas = [i_a_meas; i_b_meas; i_c_meas]

        i_abc_ref = InvClarkeTransform(Inv_p_q_v(ClarkeTransform(v_abc_meas), pq0_ref))

        i_abc_ref = clamp.(i_abc_ref/env.nc.parameters["source"][1]["i_limit"], -0.95, 0.95)

        #i_abc_ref= [0; 0; 0]

        #println(i_abc_ref)

        # if any((abs.(env.state)).>=1)
        #     return -1
        # end

        r = 1 - 1/3*(sum((abs.(i_abc_ref - i_abc_meas)/2).^0.5))

        r += (0.2 - clamp(maximum(abs.(env.state)), 0.2, 1.0))

        return r * (1 - 0.99)/2
    end

end

function reference(t)
    θ = 2*pi*50*t
    θph = [θ; θ - 120π/180; θ + 120π/180]
    return +10 * cos.(θph) 
end

featurize_ddpg = function(state, env, name)
    if name != "classic"
        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        state = vcat(state, reference(env.t)/norm_ref)
    end
end

function reward_function2(env, name = nothing)
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
    reward_function = reward_function2,
    featurize = featurize_ddpg,
    action_delay = 0,
    verbosity = 0)


controllers = SetupAgents(env)

# Let it learn till 19_000 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function learn()
    while true
        Learn(controllers, env, num_episodes = 50)
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
