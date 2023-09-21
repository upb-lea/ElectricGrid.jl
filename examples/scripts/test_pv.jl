using ElectricGrid
using LinearAlgebra
using DataFrames
using PlotlyJS

CM = [ 0.  1.
      -1.  0.]

R_load, L_load, X, Z = ParallelLoadImpedance(1e3, 1., 230)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("source_type" => "pv",
                        "fltr"=>"L", "L1"=>0.001, "C"=>2e-8, "L2"=>0.001, "R_C"=>0.01,
                        "R2"=> 0.05, "mode" => 1, "parallel" => 20, "serial" => 20,
                        "i_limit" => 1e4, "v_limit" => 1e5, "module_N_cell" => 36),
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "R", "R" => 1e100, "i_limit" => 10e4, "v_limit" => 10e4),
                        ],
        "cable"   => Any[
                        Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4, "i_limit" => 1e4, "v_limit" => 1e6,),
                        ],
        "grid" => Dict{Any, Any}("fs"=>1e4, "phase"=>3, "v_rms"=>230, "f_grid" => 50, "ramp_end" => 0.9)
    )


env = ElectricGridEnv(num_sources=1, num_loads=1, CM = CM, parameters = parameters, t_end=2, verbosity=2)

Multi_Agent = SetupAgents(env)

hook = DataHook(collect_state_ids = ["source1_i_L1_a", "source1_v_C_cables_a"],
                collect_action_ids = [],
                collect_vdc_ids = [1],
                collect_soc_ids = [1],
                collect_idc_ids = [1],)

hook = Simulate(Multi_Agent, env, hook=hook)

p = RenderHookResults(hook = hook,
                    states_to_plot  = ["source1_i_L1_a", "source1_v_C_cables_a"],
                    vdc_to_plot     = [1],
                    idc_to_plot     = [1],
                    return_plot     = false);
