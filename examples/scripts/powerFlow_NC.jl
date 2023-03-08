#using JEG   warum geht das nicht?
using DrWatson
@quickactivate "JEG"

using ReinforcementLearning
using PlotlyJS

include(srcdir("node_constructor.jl"))
include(srcdir("electric_grid_env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("ClassicalControl.jl"))
include(srcdir("power_system_theory.jl"))
include(srcdir("MultiAgentGridController.jl"))


CM = [  0   0   0   1
        0   0   0   2
        0   0   0   3
        -1  -2  -3  0  ] # user specified



num_source = 3 # user
num_load = 1 # user

# example to calculate RC load for 1000 VA per phase
#R_load, C_load, _, Z = ParallelLoadImpedance(3000, -0.9, 230)
R_load, L_load, _, Z = ParallelLoadImpedance(3000, 0.9, 230)

#===
#TODO: shift this for every load to the node_constructor
for (index, load) in enumerate(parameters["load"])
    # example for RL load
    print(load)
    if load["impedance"] == "RL"
        if !haskey(load, "pwr")
            # parallel R||L
            load["Z"] = 1im*parameters["grid"]["fg"]*2*pi*load["R"]*load["L"]/(load["R"]+1im*parameters["grid"]["fg"]*2*pi*load["L"])
        end
    elseif load["impedance"] == "RC"
        load["Z"] = load["R"]/(1+1im*parameters["grid"]["fg"]*2*pi*load["C"]*load["R"])
    end
    load["pwr"] = parameters["grid"]["v_rms"]^2 / abs(load["Z"]) * parameters["grid"]["phase"]
end
print(parameters["load"])
===#





# define same setting for env and check if 
parameters_nc = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("fltr"=>"LC", "pwr"=>1000000, "v_pu_set" => 1.0, "mode" => 1, "control_type" => "classic"),
                    Dict{Any, Any}("fltr"=>"LC", "pwr"=>1000000, "v_pu_set" => 1.0, "mode" => 1, "control_type" => "classic",
                                    "p_set"=>500, "q_set"=>500),
                    Dict{Any, Any}("fltr"=>"LC", "pwr"=>1000000, "v_pu_set" => 1.0, "mode" => 1, "control_type" => "classic"),
                    #Dict{Any, Any}("fltr"=>"LC", "pwr"=>1000, "v_pu_set" => 1.0, "mode" => "PQ Control”"),
                    #Dict{Any, Any}("fltr"=>"LC", "pwr"=>1000, "v_pu_set" => 1.0, "mode" => "Voltage Control”"),
                    ],
    "load"   => Any[
                    #Dict{Any, Any}("R"=>10, "L" => 0.16, "impedance"=>"RL")#, "pwr"=>1000)
                    Dict{Any, Any}("R"=>R_load, "L"=>L_load, "impedance"=>"RL") # ->1kVA
                    #Dict{Any, Any}("R"=>3526.058777777777777776, "L"=>23.0003863038165653084, "impedance"=>"RL") # ->1000kVA
                    #Dict{Any, Any}("R"=>R_load, "C"=>C_load, "impedance"=>"RC")
                    ],
    "grid"   => Dict{Any, Any}("fs"=>50.0, "phase"=>1, "v_rms"=>230, "f_grid" => 50),
    "cable" => Any[
                    Dict{Any, Any}("len"=>1),
                    Dict{Any, Any}("len"=>2),
                    Dict{Any, Any}("len"=>3)
                    ]
)
env = ElectricGridEnv(ts = 1e-4, use_gpu = false, CM = CM, num_sources = num_source, num_loads = num_load, parameters = parameters_nc, maxsteps = 100, action_delay = 1)

println("########################################################################")
println()
println("Check if the limits are fed correctly to the env")
println()
println("Current limit in env results in parameters:")
println(env.norm_array[10])
println("Current limit in NC results in parameters of cable 1:")
println(env.nc.parameters["cable"][1]["i_limit"])
println()
println()

