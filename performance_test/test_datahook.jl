using ElectricGrid
using PlotlyJS
using TimerOutputs
using ReinforcementLearning

reset_timer!(to)

t_end    = 1.0     # total run time, seconds
CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]

R_load, L_load, _, _ = ParallelLoadImpedance(50e3, 0.99, 230)
R_load₂, C_load, _, _ = ParallelLoadImpedance(50e3, -0.99, 230)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 400e3,
                                        "control_type" => "RL",
                                        "fltr" => "LC",
                                        "v_δ_set" => 4),
                        Dict{Any, Any}("pwr" => 100e3,
                                        "fltr" => "LCL",
                                        "control_type" => "RL",
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

env = ElectricGridEnv(CM = CM, parameters = parameters, t_end = t_end, verbosity = 0)
Multi_Agent = SetupAgents(env)

all_sources = collect(1:env.nc.num_sources)
all_loads = collect(1:env.nc.num_loads)
all_cables = collect(1:env.nc.num_connections)

hook = DataHook3(collect_sources  = all_sources,
                collect_cables   = all_cables,
                #collect_loads    = all_loads
                )

disable_timer!(to)
    Simulate(Multi_Agent, env, hook = hook)
enable_timer!(to) 

Simulate(Multi_Agent, env, hook = hook)


show(to)