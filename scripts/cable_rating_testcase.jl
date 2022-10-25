
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


#################################
# 2 Calulate characteristic impedance Z_0 of cables
# Should be done in NC (?) on the fly while layouting the grid
# TODO: Add Z_0 to parameter_dict of the cables?
#################################

Z_0 = Dict{Any, Any}(
        "Z_01" => sqrt(parameters["cable"][1]["Lb"]/parameters["cable"][1]["Cb"]),
        "Z_02" => sqrt(parameters["cable"][2]["Lb"]/parameters["cable"][2]["Cb"])
        )

# Take for V = parameters["grid"] * sqrt(3) (-> phase neutral to phase-phase) ~ 400 V
# Add 5 % safety margin
P_max = Dict{Any, Any}(
    "P_max1" => 1.05 * parameters["grid"]["v_rms"] * sqrt(3) / Z_0["Z_01"],
    "P_max2" => 1.05 * parameters["grid"]["v_rms"] * sqrt(3) / Z_0["Z_02"]
    )

env = SimEnv(num_sources = 2, num_loads = 1, CM = CM, parameters = parameters, maxsteps=600)

