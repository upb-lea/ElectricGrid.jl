using ElectricGrid
using LinearAlgebra
using DataFrames
using PlotlyJS
"""
This scipt contains parts of the DC_linkModels.ipynb notebook.
For comments and more documentation see refering notebookin examples/notebooks.
"""

# constructing a simple grid
CM = [ 0. 1.
      -1. 0.]

# create suitable load
R_load, L_load, X, Z = ParallelLoadImpedance(5e3, .95, 230)

# configer the other parameters
parameters = Dict{Any,Any}(
    "source" => Any[
        Dict{Any,Any}("source_type" => "battery", "pwr" => 15000.0, "fltr" => "L",
            "L1" => 0.001, "C" => 2e-8, "L2" => 0.001, "R_C" => 0.01, "R2" => 0.05,
            "Q" => 100 * 26 * 3600, "mode" => 1, "parallel" => 48, "serial" => 200,
            "i_bat_limit" => 500, "i_limit" => 1e3, "v_limit" => 1e3), # p_set maybe 15k
            ],
    "load" => Any[
        Dict{Any,Any}("impedance" => "RL", "R" => R_load, "L" => L_load, "i_limit" => 500, "v_limit" => 1000)
    ],
    "grid" => Dict{Any,Any}("fs" => 1e4, "phase" => 3, "v_rms" => 230, "f_grid" => 50, "ramp_end" => 0.9)
)

# create env
env = ElectricGridEnv(num_sources=1, num_loads=1, CM=CM, parameters=parameters, t_end=2, verbosity=2)

#create agent
Multi_Agent = SetupAgents(env)

# define information we like to collect
hook = DataHook(collect_state_ids = [],
                collect_action_ids = [],
                collect_vdc_ids = [1],
                collect_soc_ids = [1],
                collect_idc_ids = [1],)

# run experiment
hook = Simulate(Multi_Agent, env, hook=hook)

# render results
p = RenderHookResults(hook = hook,
                    states_to_plot  = [],
                    vdc_to_plot     = [1],
                    soc_to_plot     = [1],
                    idc_to_plot     = [1],
                    return_plot     = false)
