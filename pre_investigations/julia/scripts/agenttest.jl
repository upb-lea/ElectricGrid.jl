using DrWatson
@quickactivate "MicroGridSimWithRL"

# using DifferentialEquations
# using Sundials
# using Plots
# using LinearAlgebra
# using ControlSystems
# using BenchmarkTools
using ReinforcementLearning
# using IntervalSets


include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))

CM = [ 0.  1.
        -1.  0.]

parameters = Dict()
parameters["source"] = [Dict("fltr" => "LCL", "R" => 0.4, "L1" => 2.3e-3, "L2" => 2.3e-3,"C" => 10e-6)]
parameters["cable"] = [Dict("R" => 0.722, "L" => 0.955e-3, "C" => 8e-09)]
parameters["load"] = [Dict("impedance" => "R", "R" => 14)]

Grid_FC = NodeConstructor(num_source=1, num_loads=1, CM=CM, parameters=parameters)

#draw_graph(Grid_FC)   ---   not yet implemented

A, B, C, D = get_sys(Grid_FC)

env = SimEnv(A=A, B=B, C=C)

