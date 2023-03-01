using Dare

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Configuration 

#-------------------------------------------------------------------------------
# Time simulation

t_end    = 0.8     # total run time, seconds

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

R_load, L_load, _, _ = Parallel_Load_Impedance(200e3, 0.99, 230)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3),
                        Dict{Any, Any}("pwr" => 100e3),
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", 
                                        "R" => R_load, 
                                        "L" => L_load),
                        ],
        "cable"   => Any[
                        Dict{Any, Any}("R" => 0.1, 
                                        "L" => 0.25e-3, 
                                        "C" => 0.1e-4),
                        Dict{Any, Any}("R" => 0.1, 
                                        "L" => 0.25e-3, 
                                        "C" => 0.1e-4),
                        ],
    )
#_______________________________________________________________________________
# Defining the environment

env = SimEnv(CM = CM, parameters = parameters, t_end = t_end, verbosity = 2)

#_______________________________________________________________________________
# initialising the agents 

Multi_Agent = setup_agents(env)
Source = Multi_Agent.agents["classic"]["policy"].policy.Source

#_______________________________________________________________________________
# running the time simulation 

hook = simulate(Multi_Agent, env)

#_______________________________________________________________________________
# Plotting

plot_hook_results(hook = hook, 
                    states_to_plot  = [], 
                    actions_to_plot = [],  
                    power_p         = [1 2], 
                    power_q         = [1 2], 
                    v_mag           = [1 2], 
                    i_mag           = [],
                    freq            = [1 2],
                    angles          = [1 2])

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")