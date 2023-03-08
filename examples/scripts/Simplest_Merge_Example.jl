using JEG
using ReinforcementLearning
using IntervalSets
using Statistics

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Configuration

#-------------------------------------------------------------------------------
# Time simulation

Timestep = 1e-4  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 0.1   # total run time, seconds
num_eps  = 100    # number of episodes to run

#-------------------------------------------------------------------------------
# Connectivity Matrix

#-------------------------------------------------------------------------------
# Parameters

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
=#
# Filter
R_1 = 400e-3;
L_1 = 2.3e-3;
R_c = 7e-3;
C_1 = 1e-5;

S_load = 100e3
pf_load = 0.98
v_rms = 230
R_load, L_load, X, Z = ParallelLoadImpedance(S_load, pf_load, v_rms)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3, "control_type" => "RL", "fltr" => "L"),
                        Dict{Any, Any}("pwr" => 200e3, "fltr" => "LC", "control_type" => "RL"),
                        Dict{Any, Any}("pwr" => 200e3, "fltr" => "LC", "control_type" => "classic", "mode" => 1),
                        #Dict{Any, Any}("pwr" => 200e3, "control_type" => "RL", "mode" => "sac", "fltr" => "L"),
                        #Dict{Any, Any}("pwr" => 200e3, "control_type" => "RL", "mode" => "ddpg2", "fltr" => "L"),
                        ],
         "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => R_load, "L" => L_load,"v_limit"=>1e4, "i_limit"=>10e4),
                        ],
        #"cable"   => Any[
        #                Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-8, "i_limit" => 10e8,"v_limit"=>1000, ),
         #               ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.04, "phase" => 3 )
    )


#agents = Dict("ddpg" => RL.DDPG(params1), "sac" => RL.DDPG(params1), )

#_______________________________________________________________________________
# Defining the environment

function reference(t)
    θ = 2*pi*50*t
    θph = [θ; θ - 120π/180; θ + 120π/180]
    #i = [10 * cos.(2*pi*50*t .- 2/3*pi*(i-1)) for i = 1:3]
    return [-10, -10, -10]#.* cos.(θph)#, -1, -1]#[2, 2, 2]# * cos.(θph)
end

function reward(env, name = nothing)

    index_1 = findfirst(x -> x == "source1_i_L1_a", env.state_ids)
    index_2 = findfirst(x -> x == "source1_i_L1_b", env.state_ids)
    index_3 = findfirst(x -> x == "source1_i_L1_c", env.state_ids)

    u_l1 = env.state[index_1]
    u_l2 = env.state[index_2]
    u_l3 = env.state[index_3]

    # better:
    #u = env.state[findall(x -> x in env.state_ids, env.state_ids_RL)]
    # problem, u_cable is state as well

    u = [u_l1, u_l2, u_l3]

    if any(abs.(u).>1)
        return -1
    else

        refs = reference(env.t)
        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        #r = 1-(mean(abs.(refs/norm_ref - u)/2))                 # TODO: replace by correct entry of norm array
        r = 1-(mean((abs.(refs/norm_ref - u)/2).^0.5))
        return r #* (1-0.999f0)
    end
end




env = SimEnv(ts = Timestep, num_sources = 3, num_loads = 1, parameters = parameters, t_end = t_end, verbosity = 2, reward_function = reward, action_delay = 0)

#_______________________________________________________________________________
# Setting up data hooks

hook = data_hook(collect_sources  = [1],
                plot_rewards = true)

#_______________________________________________________________________________
# Running the Time Simulation


ma = setup_agents(env)

learn(ma, env, num_episodes = num_eps)

hook = data_hook(collect_state_ids = env.state_ids,
                collect_action_ids = env.action_ids)


hook = simulate(ma, env, hook=hook)

#RLBase.run(ma, env, StopAfterEpisode(num_eps), hook);
#_______________________________________________________________________________
# Plotting


plot_hook_results(hook = hook,
                    #episode = hook.bestepisode,
                    #episode = num_eps,
                    #states_to_plot  = ["source1_i_L1_a", "source2_i_L1_a", "source2_v_C_filt_a"],
                    states_to_plot  = env.state_ids,
                    actions_to_plot = env.action_ids,
                    plot_reward=true)

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")
