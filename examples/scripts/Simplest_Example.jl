using Dare

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Configuration 

#-------------------------------------------------------------------------------
# Time simulation

Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 1.0     # total run time, seconds

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
                        Dict{Any, Any}("pwr" => 200e3, "mode" => 4),
                        Dict{Any, Any}("pwr" => 100e3, "mode" => 4),
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => 2.64, "L" => 0.006),
                        ],
        "cable"   => Any[
                        Dict{Any, Any}("R" => 0.208, "L" => 0.00025, "C" => 0.4e-3, "i_limit" => 10e4,),
                        Dict{Any, Any}("R" => 0.208, "L" => 0.00025, "C" => 0.4e-3, "i_limit" => 10e4,),
                        ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.04)
    )
#_______________________________________________________________________________
# Defining the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2)

#_______________________________________________________________________________
# Setting up data hooks

hook = DataHook(collect_sources  = [1 2],
                vrms             = [1 2], 
                irms             = [1 2], 
                power_pq         = [1 2],
                freq             = [1 2],
                angles           = [1 2])

#_______________________________________________________________________________
# initialising the agents 

Multi_Agent = setup_agents(env)
Source = Multi_Agent.agents["classic"]["policy"].policy.Source

#_______________________________________________________________________________
# running the time simulation 

RLBase.run(Multi_Agent, env, StopAfterEpisode(num_eps), hook);  
#_______________________________________________________________________________
# Plotting

plot_hook_results(hook = hook, 
                    states_to_plot  = [], 
                    actions_to_plot = [],  
                    power_p         = [1 2], 
                    power_q         = [], 
                    vrms            = [1 2], 
                    irms            = [],
                    freq            = [1 2],
                    angles          = [1 2])

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")