using ElectricGrid
using DifferentialEquations
using PlotlyJS
#using ReinforcementLearning
CM = [0. 1.
     -1. 0.]

S_source = 1e4

S_load = 1e2
pf_load = 1
v_rms = 230
R_load, L_load, X, Z = ParallelLoadImpedance(S_load, pf_load, v_rms)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => S_source, "control_type" => "classic", "mode" => "Swing", "fltr" => "LC", "i_limit" => 1e4, "v_limit" => 1e4,),
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "RLC", "R" => R_load, "v_limit" => 1e4, "i_limit" => 1e4)
                        ],
        "cable"   => Any[
                       Dict{Any, Any}("R" => 1e-3, "L" => x->1e-4, "C" => 1e-4),
                       ],
        "grid" => Dict{Any, Any}("fs"=>1e4, "phase"=>3, "v_rms"=>230, "f_grid" => 50, "ramp_end"=>0.00)
    )


env = ElectricGridEnv(CM = CM, parameters = parameters, verbosity = 2)
#env = ElectricGridEnv(num_sources = 2, num_loads = 1)
sol = []
for i = 1:1000
    env([0.2, 0.2, 0.2])
    append!(sol,env.x[2])
end
t_t = collect(env.t0:env.ts:env.t)

plot(t_t,sol)