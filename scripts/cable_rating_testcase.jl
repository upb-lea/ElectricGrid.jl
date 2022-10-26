
using Dare

#=
 File to prepare all functions needed in #24
 Most of the developed function will have to be shifted to the NodeConstructor or Power_System_Theory.jl

Based on documentation can be found in https://github.com/upb-lea/dare/issues/24 
=#


#################################
# 1 Generate random power system
#################################
CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

parameters = Dict{Any, Any}(
    # generate 2 sources, each 10 kVA with autodesigned filter
    "source" => Any[
                    Dict{Any, Any}("pwr"=>10000.0, "fltr"=>"LC", "source_type" => "ideal"),
                    Dict{Any, Any}("pwr"=>10000.0, "fltr"=>"L", "source_type" => "ideal")
                    ],
    # TODO: what if active load? 
    "load"   => Any[
                    Dict{Any, Any}("R"=>14, "impedance"=>"R")
                    #Dict{Any, Any}("C"=>0.381, "L"=>0.2, "R"=>14, "impedance"=>"RLC")
                    ],
    # desing the cables underrated
    "cable"  => Any[
                    Dict{Any, Any}("C"=>4.0e-7, "L"=>0.000264, "R"=>0.722),
                    Dict{Any, Any}("C"=>4.0e-7, "L"=>0.000264, "R"=>0.722)
                    ],
    "grid"   => Dict{Any, Any}("fs"=>10000.0, "phase"=>1, "v_rms"=>230)
)


env = SimEnv(num_sources = 2, num_loads = 1, CM = CM, parameters = parameters, maxsteps=600)

#################################
# 2 Calulate characteristic impedance Z_0 of cables
# Location: Should be done in NC (?) on the fly while layouting the grid
# TODO: Add Z_0 to parameter_dict of the cables?
#################################

Z_0 = Dict{Any, Any}(
        "Z_01" => sqrt(parameters["cable"][1]["L"]/parameters["cable"][1]["C"]),
        "Z_02" => sqrt(parameters["cable"][2]["L"]/parameters["cable"][2]["C"])
        )

#################################
# 3 Calulate SIL as P_max of cables from first design
# Location: Should be used in NC (?); implement in Power_System_Theory?
# TODO: Add P_max (in the end!!) to parameter_dict of the cables?
#################################
# Take for V = parameters["grid"] * sqrt(3) (-> phase neutral to phase-phase) ~ 400 V
# Add 5 % safety margin
P_max = Dict{Any, Any}(
    "P_max1" => 1.05 * parameters["grid"]["v_rms"] * sqrt(3) / Z_0["Z_01"],
    "P_max2" => 1.05 * parameters["grid"]["v_rms"] * sqrt(3) / Z_0["Z_02"]
    )

#################################
# 4 Solve power flow equation
# First: Calulate bus admittance matrix Y
#################################

# per hand - to be automated based on CM
# which cable number belongs to which source?
omega = 2*π*50  # TODO: env.nc.parameters["fg"]??
G_RL = env.nc.parameters["cable"][1]["R"] / (env.nc.parameters["cable"][1]["R"]^2 - omega^2 *env.nc.parameters["cable"][1]["L"]^2)
B_RL = (omega *env.nc.parameters["cable"][1]["L"])/ (env.nc.parameters["cable"][1]["R"]^2 - omega^2 *env.nc.parameters["cable"][1]["L"]^2)

Y_11 = env.nc.parameters["cable"][1]["C"] / 2 # C is splitted to 2 capacitors
Y_12 = 0
Y_13 = G_RL + B_RL*im
Y_22 = env.nc.parameters["cable"][1]["C"] / 2 # C is splitted to 2 capacitors
Y_21 = 0
Y_23 = G_RL + B_RL*im
Y_32 = G_RL + B_RL*im
Y_31 = G_RL + B_RL*im
Y_33 = env.nc.parameters["cable"][1]["C"] # 2 times!

Y_bus = [[Y_11, Y_12, Y_13],
         [Y_21, Y_22, Y_23],
         [Y_31, Y_32, Y_33]]