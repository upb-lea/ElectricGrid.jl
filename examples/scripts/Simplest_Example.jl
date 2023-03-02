using Dare

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Configuration 

#-------------------------------------------------------------------------------
# Time simulation

Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 0.8     # total run time, seconds
num_eps  = 1       # number of episodes to run

#-------------------------------------------------------------------------------
# Connectivity Matrix

CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]

#= CM = [ 0. 1.
        -1. 0.] =#

#-------------------------------------------------------------------------------
# Parameters

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
=#

R_load, L_load, _, _ = Parallel_Load_Impedance(100e3, 0.99, 230)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3),
                        Dict{Any, Any}("pwr" => 100e3),
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => R_load, "L" => L_load),
                        ],
        "cable"   => Any[
                        Dict{Any, Any}("R" => 0.208, "L" => 0.00025, "C" => 0.4e-3, "i_limit" => 10e4,),
                        Dict{Any, Any}("R" => 0.208, "L" => 0.00025, "C" => 0.4e-3, "i_limit" => 10e4,),
                        ],
        #= "cable"   => Any[
                        Dict{Any, Any}("R" => 1e-5, "L" => 1e-6, "C" => 1e-7, "i_limit" => 10e4,),
                        ], =#
        "grid" => Dict{Any, Any}("ramp_end" => 0.04)
    )
#_______________________________________________________________________________
# Defining the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2, action_delay = 1)

#num_sources = 1
#num_loads = 0
#p_load_total, q_load_total, s_load_total, s_source_total = CheckPowerBalance(parameters,  num_sources, num_loads, CM)
#_______________________________________________________________________________
# initialising the agents 

Multi_Agent = setup_agents(env)
Source = Multi_Agent.agents["classic"]["policy"].policy.Source

#_______________________________________________________________________________
# running the time simulation 

hook = simulate(Multi_Agent, env, num_episodes = num_eps)

#_______________________________________________________________________________
# Plotting

plot_hook_results(hook = hook, 
                    states_to_plot  = ["cable1_i_L_a"], 
                    actions_to_plot = [],  
                    power_p         = [1 2], 
                    power_q         = [], 
                    vrms            = [1 2], 
                    irms            = [],
                    freq            = [1 2],
                    angles          = [1 2],
                    i_sat           = [],
                    v_sat           = [],
                    i_err_t         = [],
                    v_err_t         = [])

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")