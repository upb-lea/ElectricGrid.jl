
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

# CM = [ 0. 0. 1.
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

R_load, L_load, _, _ = Parallel_Load_Impedance(220e3, 0.95, 230)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3, "mode" => 1),
                        # Dict{Any, Any}("pwr" => 100e3, "mode" => 4),
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => 1.5, "L" => 1e-3),
                        ],
        # "cable"   => Any[
        #                 Dict{Any, Any}("R" => 0.208, "L" => 0.00025, "C" => 0.4e-3, "i_limit" => 10e4,),
                        # ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.04)
    )
#_______________________________________________________________________________
# Defining the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2, action_delay = 1)

Dare.optimizer_status
print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")

num_source = 1
num_load = 1


params = copy(env.nc.parameters)
pop!(params, "cable") 

params

status = Dict{String, Any}()

# function get_optim_status()
layout_cabels(CM, num_source, num_load, parameters, 2)
Dare.optimizer_status

