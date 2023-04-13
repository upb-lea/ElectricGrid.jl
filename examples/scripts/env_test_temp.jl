using JEG

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
                        "fltr" => "L",
                        "i_limit"      => 2000.),
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "fltr" => "LC",
                        "control_type" =>
                        "RL", "mode" => "my_ddpg"),
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "fltr" => "L",
                        "control_type" =>
                        "RL", "mode" => "my_ddpg",
                        "i_limit"      => 2000.),
                    ],
        #"load"   => Any[
        #    Dict{Any, Any}("impedance" => "RLC", "R" => R_load, "v_limit" => 1e4, "i_limit" => 1e4)
        #    ],
    "grid" => Dict{Any, Any}(
        "phase" => 1,
        "ramp_end" => 0.04,)
)

function reference(t)

    #return [-1., -2.]
    return [-10., 230., -15.]
end

featurize_ddpg = function(state, env, name)
    if name == "my_ddpg"
        #norm_ref = env.nc.parameters["source"][1]["i_limit"]
        #state = vcat(state, reference(env.t)/norm_ref)

        refs = reference(env.t)
        refs[1] = refs[1] / env.nc.parameters["source"][1]["i_limit"]
        refs[2] = refs[2] / env.nc.parameters["source"][2]["v_limit"]
        refs[3] = refs[3] / env.nc.parameters["source"][3]["i_limit"]

        state = vcat(state, refs)
    end
end

function reward_function(env, name=nothing)
    if name == "my_ddpg"
        #println("Inside reward")
        state_to_control_1 = env.state[findfirst(x -> x == "source1_i_L1", env.state_ids)]
        state_to_control_2 = env.state[findfirst(x -> x == "source2_v_C_filt", env.state_ids)]
        state_to_control_3 = env.state[findfirst(x -> x == "source3_i_L1", env.state_ids)]

        state_to_control = [state_to_control_1, state_to_control_2, state_to_control_3]
        #state_to_control = [state_to_control_1, state_to_control_3]

        if any(abs.(state_to_control).>1)
            return -1
        else
            refs = reference(env.t)
            refs[1] = refs[1] / env.nc.parameters["source"][1]["i_limit"]
            refs[2] = refs[2] / env.nc.parameters["source"][2]["v_limit"]
            refs[3] = refs[3] / env.nc.parameters["source"][3]["i_limit"]

            r = 1-1/3*(sum((abs.(refs - state_to_control)/2).^0.5))

            #refs = reference(env.t)
            #norm_ref = env.nc.parameters["source"][1]["i_limit"]
            #r = 1-1/3*(sum((abs.(refs/norm_ref - state_to_control)/2).^0.5))

            #println(r)
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
    action_delay = 0,
    verbosity = 2);