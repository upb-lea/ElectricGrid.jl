using ElectricGrid
using PlotlyJS

#= Description:

    The simplest example to get a classically controlled network running.

=#

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Configuration

#-------------------------------------------------------------------------------
# Time simulation

t_end    = 0.08     # total run time, seconds

#-------------------------------------------------------------------------------
# Connectivity Matrix

CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]

#-------------------------------------------------------------------------------
# Parameters

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
=#

R_load, L_load, _, _ = ParallelLoadImpedance(50e3, 0.99, 230)
R_load₂, C_load, _, _ = ParallelLoadImpedance(50e3, -0.99, 230)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 400e3,
                                        "mode" => "VSG",
                                        "fltr" => "LC",
                                        "v_δ_set" => 4),
                        Dict{Any, Any}("pwr" => 100e3,
                                        "fltr" => "LCL",
                                        "mode" => "PQ",
                                        "p_set" => 50e3,
                                        "q_set" => 10e3,
                                        "v_pu_set" => 1.0,
                                        "v_δ_set" => 1),
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "RLC",
                                        "R" => R_load + R_load₂,
                                        "L" => L_load,
                                        "C" => C_load,
                                        "v_limit" => 10e3),
                        ],
        "cable"   => Any[
                        Dict{Any, Any}("R" => 0.1,
                                        "L" => 0.0025e-3,
                                        "C" => 0.001e-4),
                        Dict{Any, Any}("R" => 0.1,
                                        "L" => 0.0025e-3,
                                        "C" => 0.001e-4),
                        ],
    )
#_______________________________________________________________________________
# Defining the environment

env = ElectricGridEnv(CM = CM, parameters = parameters, t_end = t_end, verbosity = 2)

#_______________________________________________________________________________
# initialising the agents

Multi_Agent = SetupAgents(env)
Source = Multi_Agent.agents["classic"]["policy"].policy.Source

#_______________________________________________________________________________
# running the time simulation

hook = Simulate(Multi_Agent, env)

#_______________________________________________________________________________
# Plotting

p = RenderHookResults(hook = hook,
                    states_to_plot  = ["source1_v_C_filt_a", "source1_v_C_filt_b", "source1_v_C_filt_c"],
                    actions_to_plot = [],
                    power_p_inv     = [1 2],
                    power_q_inv     = [1 2],
                    power_p_poc     = [1 2],
                    power_q_poc     = [1 2],
                    v_mag_inv       = [1 2],
                    v_mag_poc       = [1 2],
                    i_mag_inv       = [1 2],
                    i_mag_poc       = [1 2],
                    i_sat           = [],
                    i_dq            = [],
                    v_dq            = [],
                    freq            = [1 2],
                    angles          = [1 2],
                    return_plot     = true)


print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")
