# Copy of Simplest_Example.jl
# invoke layout_cables() by specifying the cables in parameters #L832 from nodeconstructor
#  - 2 sources - 1 load  (mode 1)
# compare P, Q, p, q 
# Start with putting every source in "mode" => 1 or "Swing". This is open loop - no control
# All sources given filters that are "L" with small "L1" and "R1".
# Start with only sources no loads.
# All sources set to 1 p.u. peak
# All sources set to delta = 0 degrees at first. Then vary angles.
# Test for Vrms = 230 V, 100 V to make sure it works for different voltages
# Add in loads and try different topologies.

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


R, L_C, _, _ = Parallel_Load_Impedance(190e3, 0.95, 230; fsys = 50)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3, "fltr" => "L", "R1" => 1e-3,"L1" => 1e-6, "mode" => 1, "v_δ_set" => 0.0, "pu" => 1.0),
                        # Dict{Any, Any}("pwr" => 200e3, "fltr" => "L", "R1" => 1e-3,"L1" => 1e-6, "mode" => 1, "v_δ_set" => 0.0, "pu" => 1.0),
                        ],
        # No load
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => 2.64, "L" => 1e-3),
                        ], 
        # "cable"   => Any[
        #                 Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4, "i_limit" => 10e4,),
        #                 ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.0) #, "vrms" => 230.0) #, "f" => 50.0, "L" => 1e-6, "R" => 1e-3, "C" => 1e-4, "i_limit" => 10e4,)
    )
#_______________________________________________________________________________
# Defining the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2)

#_______________________________________________________________________________
# Setting up data hooks

hook = DataHook(collect_vrms_ids = [1 ], 
                collect_irms_ids = [1 ], 
                collect_pq_ids   = [1 ], #collecting p and q for sources 1, 2
                collect_freq     = [1 ],
                collect_sources  = [1 ])
#_______________________________________________________________________________
# Running the Time Simulation

Power_System_Dynamics(env, hook)

#_______________________________________________________________________________
# Plotting

plot_hook_results(hook = hook, 
                    states_to_plot  = [], 
                    actions_to_plot = [],  
                    p_to_plot       = [1], 
                    q_to_plot       = [1], 
                    vrms_to_plot    = [1], 
                    irms_to_plot    = [],
                    freq_to_plot    = [])

# steady state values from Source p , q from data_hook
for i in env.nc.num_sources
    println("Source $i : p = ", hook.df[!,"source$(i)_p"][end], " q = ", hook.df[!,"source$(i)_q"][end])
end

p_test = 0.0;
q_test = 0.0;

for i in env.nc.num_sources
    p_test += hook.df[!,"source$(i)_p"][end]
    q_test += hook.df[!,"source$(i)_q"][end]
    @show p_test, q_test
end



num_source = length(parameters["source"])
num_load   = length(parameters["load"])

total_P_load, total_Q_load, s_load_total, total_S_source = CheckPowerBalance(parameters, num_source, num_load, CM)

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")


# Time domain : p, q
# Freq domain : P, Q