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

num_nodes = 4
num_sources = 2

CM, num_cables  = SmallWorld(num_nodes, p = 0.2, Z = 2, num_sources = num_sources)
#CM, num_cables = Barabasi_Albert(num_nodes)

#-------------------------------------------------------------------------------
# Parameters

parameters = Dict{Any, Any}()

parameters["source"], total_gen = Source_Setup(num_sources, random = 0, mode = 4)
parameters["load"] = Load_Setup(num_nodes - num_sources, total_gen, random = 0)
parameters["cable"] = Cable_Length_Setup(num_cables, random = 0)
parameters["grid"] = Dict("v_rms" => 230, "ramp_end" => 0.04, "process_start" => 0.05)

#_______________________________________________________________________________
# Defining the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2)

#_______________________________________________________________________________
# Setting up data hooks

hook = DataHook(collect_vrms_ids = 1:num_sources, 
                collect_irms_ids = 1:num_sources, 
                collect_pq_ids   = 1:num_sources,
                collect_freq     = 1:num_sources,
                collect_sources  = 1:num_sources)

#_______________________________________________________________________________
# Running the Time Simulation

Power_System_Dynamics(env, hook)

#_______________________________________________________________________________
# Plotting

# Spring Layout (Layout = 3) is better for Barabasi-Albert
# Circular Layout (Layout = 1) is better for SmallWolrd
drawGraph(CM, parameters, Layout = 3)

plot_hook_results(hook = hook, 
                    states_to_plot  = [], 
                    actions_to_plot = [],  
                    p_to_plot       = 1:num_sources, 
                    q_to_plot       = [], 
                    vrms_to_plot    = 1:num_sources, 
                    irms_to_plot    = [],
                    freq_to_plot    = [])

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")
