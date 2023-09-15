using ElectricGrid
using LinearAlgebra
using DataFrames
using PlotlyJS

# CM = [ 0.  1.  2.
#       -1.  0.  0.
#       -2.  0.  0.]

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

hook = DataHook(collect_state_ids = ["source1_i_L1_a", "source1_v_C_cables_a"],#env.state_ids,
                collect_action_ids = [],#env.action_ids,
                collect_vdc_ids = [1],
                collect_soc_ids = [1],
                collect_idc_ids = [1],)

hook = Simulate(Multi_Agent, env, hook=hook)

p = RenderHookResults(hook = hook,
                    states_to_plot  = ["source1_i_L1_a", "source1_v_C_cables_a"],#GetStateIds(env.nc),
                    actions_to_plot = [],
                    power_p_inv     = [],
                    power_q_inv     = [],
                    power_p_poc     = [],
                    power_q_poc     = [],
                    v_mag_inv       = [],
                    v_mag_poc       = [],
                    i_mag_inv       = [],
                    i_mag_poc       = [],
                    freq            = [],
                    angles          = [],
                    vdc_to_plot     = [1],
                    soc_to_plot     = [],
                    idc_to_plot     = [1],
                    return_plot     = false);


# for val in range(-1,1,step=0.001)
#     # println(val)
#     x = env([val, -val/2, -val/2])
#     # println(env.v_dc)

#     append!(states, x)
#     append!(A, env.v_dc[1])
# end

# t = collect(1:size(A)[1])

# Plot([scatter(x=collect(1:size(A)[1]), y=A),
#     #   scatter(x=collect(1:size(B)[1]), y=states[1:6:end]),
#     #   scatter(x=collect(1:size(A)[1]), y=collect(range(-1,1,step=0.0001)))
#       ])
