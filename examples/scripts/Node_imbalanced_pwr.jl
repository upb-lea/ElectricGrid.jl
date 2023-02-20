# Copy of Simplest_Example.jl
# invoke layout_cables() by specifying the cables in parameters #L832 from nodeconstructor
# Verify if the optimizer fails due to the load_pwr being too large for the source

using Dare

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Configuration 

#-------------------------------------------------------------------------------
# Time simulation

Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 0.2     # total run time, seconds

#---------------------------------------------------------------------------
# Connectivity Matrix

#  CM = [ 0. 0. 1.
#         0. 0. 2.
#         -1. -2. 0.] 

CM = [ 0. 1.
        -1. 0.]

#-------------------------------------------------------------------------------
# Parameters

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
=#

source_pwr  =   200e3 

# Toggle theses lines to see the optimizer behaviour
load_pwr    =   170e3 
# load_pwr    =   190e3
# load_pwr    =   210e3       # Somehow the optimizer doesn't fail for this value : [warning: The apparent power drawn from the loads exceeds the apparent power ...]
# load_pwr    =   220e3       # The optimizer fails for this value


R, L_C, _, _ = Parallel_Load_Impedance(load_pwr, 0.95, 230; fsys = 50)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => source_pwr, "fltr" => "L", "R1" => 1e-3,"L1" => 1e-6, "mode" => 1, "v_δ_set" => 0.0, "pu" => 1.0),
                        # Dict{Any, Any}("pwr" => 200e3, "fltr" => "L", "R1" => 1e-3,"L1" => 1e-6, "mode" => 1, "v_δ_set" => 0.0, "pu" => 1.0),
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => R, "L" => L_C),
                        ],
        # invoke layout cables by not specifying the cables in parameters 
        # "cable"   => Any[
        #                 Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4, "i_limit" => 10e4,),
        #                 ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.02) #, "vrms" => 230.0) #, "f" => 50.0, "L" => 1e-6, "R" => 1e-3, "C" => 1e-4, "i_limit" => 10e4,)
    )
#_______________________________________________________________________________
# Defining the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2)

#_______________________________________________________________________________
# Setting up data hooks

hook = DataHook(collect_sources  = [1 ],
                vrms             = [1 ], 
                irms             = [1 ], 
                power_pq         = [1 ],
                freq             = [1 ],
                angles           = [1 ])
#_______________________________________________________________________________
# Running the Time Simulation

Power_System_Dynamics(env, hook)

#_______________________________________________________________________________
# Plotting

plot_hook_results(hook = hook, 
                    states_to_plot  = ["cable1_i_L_a"], 
                    actions_to_plot = [],  
                    power_p         = [1 ], 
                    power_q         = [], 
                    vrms            = [1 ], 
                    irms            = [],
                    freq            = [1 ],
                    angles          = [1 ])


print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")


# Time domain : p, q
# Freq domain : P, Q