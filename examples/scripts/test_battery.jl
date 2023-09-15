using ElectricGrid
using LinearAlgebra
using DataFrames
using PlotlyJS

# CM = [0.0 1.0 2.0 0.0
#     -1.0 0.0 0.0 3.0
#     -2.0 0.0 0.0 4.0
#     0.0 -3.0 -4.0 0]

CM = [ 0. 1.
      -1. 0.]

R_load, L_load, X, Z = ParallelLoadImpedance(5e3, .95, 230)


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
    # "cable" => Any[
    #     Dict{Any,Any}("Cb" => 4.0e-9, "Lb" => 0.000264, "Rb" => 0.002, "C" => 0.4e-2, "i_limit" => 1.0e13, "v_limit" => 374.05948724768365, "len" => 1.0, "L" => 0.25e-5, "R" => 0.1e-1),
    #     # Dict{Any,Any}("Cb" => 4.0e-9, "Lb" => 0.000264, "Rb" => 0.002, "C" => 0.4e-4, "i_limit" => 1.0e13, "v_limit" => 374.05948724768365, "len" => 1.0, "L" => 0.25e-5, "R" => 0.1e-1),
    #     # Dict{Any,Any}("Cb" => 4.0e-9, "Lb" => 0.000264, "Rb" => 0.002, "C" => 0.4e-4, "i_limit" => 1.0e13, "v_limit" => 374.05948724768365, "len" => 1.0, "L" => 0.25e-5, "R" => 0.1e-1),
    #     # Dict{Any,Any}("Cb" => 4.0e-9, "Lb" => 0.000264, "Rb" => 0.002, "C" => 0.4e-4, "i_limit" => 1.0e13, "v_limit" => 374.05948724768365, "len" => 1.0, "L" => 0.25e-5, "R" => 0.1e-1),
    # ],
    "grid" => Dict{Any,Any}("fs" => 1e4, "phase" => 3, "v_rms" => 230, "f_grid" => 50, "ramp_end" => 0.9)
)


states = []
states_all = []
Charging = []
Discharging = []
V_dc_char = []
V_dc_dis = []
env = ElectricGridEnv(num_sources=1, num_loads=1, CM=CM, parameters=parameters, t_end=2, verbosity=2)

Multi_Agent = SetupAgents(env)

hook = DataHook(collect_state_ids = [],#env.state_ids,
                collect_action_ids = [],#env.action_ids,
                collect_vdc_ids = [1],
                collect_soc_ids = [1],
                collect_idc_ids = [1],)

hook = Simulate(Multi_Agent, env, hook=hook)

p = RenderHookResults(hook = hook,
                    states_to_plot  = [],#GetStateIds(env.nc),
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
                    soc_to_plot     = [1],
                    idc_to_plot     = [1],
                    return_plot     = false)

# println(env.action_space)

# for val in range(0, 10, step=0.001)
#     x = env([10])
#     # println(env.v_dc)
#     append!(states, x[env.vdc_link_voltages.sources[1].state_idx])
#     append!(states_all, x)
#     append!(Charging, env.vdc_link_voltages.sources[1].BatteryBlock.SOC)
#     append!(V_dc_char, env.v_dc[1])
# end

# for val in range(0, 10, step=0.001)
#     x = env([-1])
#     # println(env.v_dc)
#     append!(states, x[env.vdc_link_voltages.sources[1].state_idx])
#     append!(Discharging, env.vdc_link_voltages.sources[1].BatteryBlock.SOC)
#     append!(V_dc_dis, env.v_dc[1])
# end

# for val in range(0, 450, step=0.001)
#     x = env([50])
#     # println(env.v_dc)
#     append!(states, x[env.vdc_link_voltages.sources[1].state_idx])
#     append!(Charging, env.vdc_link_voltages.sources[1].BatteryBlock.SOC)
#     append!(V_dc_char, env.v_dc[1])
# end

# for val in range(0, 450, step=0.001)
#     x = env([-50])
#     # println(env.v_dc)
#     append!(states, x[env.vdc_link_voltages.sources[1].state_idx])
#     append!(Discharging, env.vdc_link_voltages.sources[1].BatteryBlock.SOC)
#     append!(V_dc_dis, env.v_dc[1])
# end

# Plot([scatter(x=collect(1:size(Charging)[1]), y=Charging),
#     scatter(x=collect(1:size(V_dc_char)[1]), y=V_dc_char),
#     # scatter(x=collect(1:size(Discharging)[1]), y=Discharging),
#     # scatter(x=collect(1:size(V_dc_dis)[1]), y=V_dc_dis)
# ])

# p1 = Plot([scatter(x=collect(1:size(states_all)[1]), y=states_all)])

# plot_data = zeros((floor(Int, size(states_all)[1]/env.nc.num_spp), env.nc.num_spp))
# for i in 1:env.nc.num_spp
#     plot_data[:,i] = states_all[i:env.nc.num_spp:end]
# end
# state_names = GetStateIds(env.nc)


# plot_data[]
# p2 = Plot([scatter(x=collect(1:size(plot_data)[1]), y=plot_data[:,i], name=state_names[i]) for i in 1:env.nc.num_spp])
