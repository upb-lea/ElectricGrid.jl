using JEG

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Parameters 

#-------------------------------------------------------------------------------
# Time simulation

Timestep = 1e-4  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 5e-4     # total run time, seconds

#-------------------------------------------------------------------------------
# Connectivity Matrix

CM = [ 0. 1.
        -1. 0.]

#-------------------------------------------------------------------------------
# Sources

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3),
                        #Dict{Any, Any}("pwr" => 200e3)
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => 2.64, "L" => 0.006) 
                        ]
    )
#_______________________________________________________________________________
# Defining the environment

env = ElectricGridEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, action_delay = 3)

#_______________________________________________________________________________
# Setting up data hooks

hook = data_hook(collect_sources  = [1])

#_______________________________________________________________________________
# Running the Time Simulation

Power_System_Dynamics(env, hook)

#_______________________________________________________________________________
# Plotting
#state_data = hook.df[!,"source1_i_L1_a"][1:5,:]'
#println(state_data)

#action_data = hook.df[!,"source1_u_a"][1:5,:]'
#println(action_data)


#plot_hook_results(hook = hook, 
#                    states_to_plot  = [], 
#                    actions_to_plot = [])

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")
