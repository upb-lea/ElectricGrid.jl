
using JEG

#=
 File to prepare all functions needed in #24
 Most of the developed function will have to be shifted to the NodeConstructor or power_system_theory.jl

Based on documentation can be found in https://github.com/upb-lea/JEG/issues/24 
=#


#################################
# 1 Generate random power system
#################################
CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

parameters = Dict{Any, Any}(
    # generate 2 sources, each 10 kVA with autodesigned filter
    "source" => Any[
                    Dict{Any, Any}("pwr"=>10000.0, "fltr"=>"LC", "source_type" => "ideal"),
                    Dict{Any, Any}("pwr"=>10000.0, "fltr"=>"L", "source_type" => "ideal"),
                    Dict{Any, Any}("pwr"=>10000.0, "fltr"=>"L", "source_type" => "ideal")
                    ],
    # TODO: what if active load? 
    "load"   => Any[
                    Dict{Any, Any}("R"=>14, "impedance"=>"R")
                    Dict{Any, Any}("C"=>0.381, "L"=>0.2, "R"=>14, "impedance"=>"RLC")
                    ],
    # desing the cables underrated
    #"cable"  => Any[
    #                Dict{Any, Any}("C"=>4.0e-7, "L"=>0.000264, "R"=>0.722),
    #                Dict{Any, Any}("C"=>4.0e-7, "L"=>0.000264, "R"=>0.722)
    #                ],
    "grid"   => Dict{Any, Any}("fs"=>10000.0, "phase"=>1, "v_rms"=>230)
)


env = SimEnv(num_sources = 5, num_loads = 4, CM = nothing, parameters = parameters, maxsteps=600)

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
# Location: Should be used in NC (?); implement in power_system_theory?
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
#=
cable_idx = 1
G_RL_1 = env.nc.parameters["cable"][cable_idx]["R"] / (env.nc.parameters["cable"][cable_idx]["R"]^2 - omega^2 *env.nc.parameters["cable"][1]["L"]^2)
B_RL_1 = (omega * env.nc.parameters["cable"][cable_idx]["L"])/ (env.nc.parameters["cable"][cable_idx]["R"]^2 - omega^2 *env.nc.parameters["cable"][1]["L"]^2)

B_shunt_1 = omega * env.nc.parameters["cable"][cable_idx]["C"] / 2 

cable_idx = 2
G_RL_2 = env.nc.parameters["cable"][cable_idx]["R"] / (env.nc.parameters["cable"][cable_idx]["R"]^2 - omega^2 *env.nc.parameters["cable"][1]["L"]^2)
B_RL_2 = (omega *env.nc.parameters["cable"][cable_idx]["L"])/ (env.nc.parameters["cable"][cable_idx]["R"]^2 - omega^2 *env.nc.parameters["cable"][1]["L"]^2)

B_shunt_2 = omega * env.nc.parameters["cable"][cable_idx]["C"] / 2 


# Load is always parrallel connection of devices
env.nc.parameters["load"]

Y_11 = G_RL_1 + B_RL_1*im  + B_shunt_1*im # C is splitted to 2 capacitors
Y_12 = 0
Y_13 = G_RL_1 + B_RL_1*im 
Y_22 = G_RL_2 + B_RL_2*im + B_shunt_2*im # C is splitted to 2 capacitors
Y_21 = 0
Y_23 = G_RL_2 + B_RL_2*im
Y_32 = G_RL_2 + B_RL_2*im
Y_31 = G_RL_1 + B_RL_1*im
Y_33 = G_RL_1 + B_RL_1*im + G_RL_2 + B_RL_2*im  + B_shunt_1*im + + B_shunt_2*im 

Y_bus_fixed = [[Y_11, Y_12, Y_13],
               [Y_21, Y_22, Y_23],
               [Y_31, Y_32, Y_33]]

G_bus = real(Y_bus)
B_bus = imag(Y_bus)
=#

################################# 
# Automated Y_bus calculation, to be shifted to the NC
# 2 sources, 1 load - in CM, first all sources, then all loads are listed
# 1 Bus per component -> 3 here. 
#       S1  S2 L1
#  S1  [ 0. 0. 1.
#  S2    0. 0. 2
#  L1   -1. -2. 0.]



num_buses = env.nc.tot_ele  # Y -> num_buses x num_buses
Y_bus = zeros(Complex{Float64}, num_buses, num_buses)

for col in 1:num_buses
    for row in 1:num_buses

        if env.nc.CM[row, col] != 0  # we have a cable connected
            # CM index defines the number of the cable
            cable_idx = abs(Int(env.nc.CM[row, col]))
            #println(cable_idx)
            #println(row)
            #println(col)
            G_RL = env.nc.parameters["cable"][cable_idx]["R"] / (env.nc.parameters["cable"][cable_idx]["R"]^2 + omega^2 *env.nc.parameters["cable"][1]["L"]^2)
            B_RL = (omega *env.nc.parameters["cable"][cable_idx]["L"])/ (env.nc.parameters["cable"][cable_idx]["R"]^2 + omega^2 *env.nc.parameters["cable"][1]["L"]^2)
            Y_bus[row, col] = -G_RL - im*B_RL
        
        elseif row == col  # diagonal elements
            # Go through all col elements of that row to find the connected cable to that bus (non zero elements in CM[:, row])
            println(env.nc.CM[row, :])
            cable_idxs = filter(n -> n !=0, env.nc.CM[row, :])
            G = 0
            B = 0
            for idx in cable_idxs
                # add all RL 
                idx = abs(Int(idx))
                G += env.nc.parameters["cable"][idx]["R"] / (env.nc.parameters["cable"][idx]["R"]^2 + omega^2 *env.nc.parameters["cable"][1]["L"]^2)
                B += (omega *env.nc.parameters["cable"][idx]["L"])/ (env.nc.parameters["cable"][idx]["R"]^2 + omega^2 *env.nc.parameters["cable"][1]["L"]^2)
                # and add all shunt C connected to that bus since diagonal element
                B += omega * env.nc.parameters["cable"][idx]["C"] / 2 
            end
            println(G)
            println(B)
            println(G + im*B)
            Y_bus[row, col] = G + im*B
        end
    end
end
