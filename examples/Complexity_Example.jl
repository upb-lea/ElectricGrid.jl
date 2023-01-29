using Dare
using Graphs

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Configuration 

#-------------------------------------------------------------------------------
# Time simulation

Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 0.2     # total run time, seconds

#-------------------------------------------------------------------------------
# Connectivity Matrix

num_nodes = 20 # The number of nodes

CM, num_cables  = SmallWorld(num_nodes)
#CM, num_cables = Barabasi_Albert(num_nodes)

#-------------------------------------------------------------------------------
# Parameters

parameters = Dict{Any, Any}()

parameters["source"] = Source_Setup(num_nodes, random = 1)
parameters["cable"] = Cable_Length_Setup(num_cables, random = 1)

#_______________________________________________________________________________
# Defining the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2)

#_______________________________________________________________________________
# Setting up data hooks

hook = DataHook(collect_vrms_ids = [1 2], 
                collect_irms_ids = [1 2], 
                collect_pq_ids   = [1 2], #collecting p and q for sources 1, 2
                collect_freq     = [1 2],
                collect_sources  = [1 2])

#_______________________________________________________________________________
# Running the Time Simulation

Power_System_Dynamics(env, hook)

#_______________________________________________________________________________
# Plotting

plot_hook_results(hook = hook, 
                    states_to_plot  = [], 
                    actions_to_plot = [],  
                    p_to_plot       = [1 2], 
                    q_to_plot       = [1 2], 
                    vrms_to_plot    = [1 2], 
                    irms_to_plot    = [],
                    freq_to_plot    = [])

drawGraph(CM, parameters, Layout = 1)

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")
