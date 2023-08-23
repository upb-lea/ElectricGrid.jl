using ElectricGrid
using BenchmarkTools
using PlotlyJS
using DifferentialEquations

CM = [0. 1.
     -1. 0.]

S_source = 1e4

S_load = 1e2
pf_load = 1
v_rms = 230
R_load, L_load, X, Z = ParallelLoadImpedance(S_load, pf_load, v_rms)

value = 5
function NonlinearInductance(drop, length)
    # drop, is to how much it decreses with i to infinity
    # length, length is the point, where its half decresed
    return x->(1-(drop/100)*1/(1+exp(-(abs(x)-length)/2)))
end

l = NonlinearInductance(95,1)

parameters = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("pwr" => S_source, "control_type" => "classic", "mode" => "Swing", "fltr" => "LC", "i_limit" => 1e4, "v_limit" => 1e4,),
                    ],
    "load"   => Any[
                    Dict{Any, Any}("impedance" => "R", "R" => R_load, "v_limit" => 1e4, "i_limit" => 1e4)
                    ],
    "cable"   => Any[
                    Dict{Any, Any}("R" => 1e-3, "L" => x->1e-4, "C" => 1e-4),
                    ],
    "grid" => Dict{Any, Any}("fs"=>1e4, "phase"=>3, "v_rms"=>230, "f_grid" => 50, "ramp_end"=>0.00)
)

env1 = ElectricGridEnv(CM = CM, parameters = parameters, verbosity = 2);

sol_1 = []
function test1()
    for i = 1:1000
    env1([0.5, 0.5, 0.5])
    append!(sol_1,env1.x[value])
    end
end
test1()
t_t = collect(env1.t0:env1.ts:env1.t)
p1 = scatter(x=t_t,y=sol_1,mode="lines",name="nonlinear")

parameters2 = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("pwr" => S_source, "control_type" => "classic", "mode" => "Swing", "fltr" => "LC", "i_limit" => 1e4, "v_limit" => 1e4,),
                    ],
    "load"   => Any[
                    Dict{Any, Any}("impedance" => "R", "R" => R_load, "v_limit" => 1e4, "i_limit" => 1e4)
                    ],
    "cable"   => Any[
                    Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4),
                    ],
    "grid" => Dict{Any, Any}("fs"=>1e4, "phase"=>3, "v_rms"=>230, "f_grid" => 50, "ramp_end"=>0.00)
)

env2 = ElectricGridEnv(CM = CM, parameters = parameters2, verbosity = 2)

sol_2 = []
function test2()
    for i = 1:1000
        env2([0.5, 0.5, 0.5])
        append!(sol_2,env2.x[value])
    end
end
test2()
t_t = collect(env2.t0:env2.ts:env2.t)
p2 = scatter(x=t_t,y=sol_2,mode="lines",name="linear")

plot([p1,p2])