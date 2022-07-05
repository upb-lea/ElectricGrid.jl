using DrWatson
@quickactivate "MicroGridSimWithRL"

using DifferentialEquations
using Sundials
using Plots
using LinearAlgebra
using ControlSystems
using BenchmarkTools
using ReinforcementLearning
using IntervalSets


include(srcdir("variables.jl"))
include(srcdir("env.jl"))

i_lim = 60.0
v_lim = 700.0
v_dc=230.0
env = SimEnv(A=A, B=B, C=C, norm_array=[i_lim, v_lim, i_lim, i_lim, v_lim, i_lim], v_dc=v_dc)

# --- Constant Input ---

RLBase.reset!(env)

output = Vector{Float64}()
for i = 1:50
    env([230.0 / 300.0, 230.0 / 300.0])
    append!(output, env.state[2])
end

plot(output)


# --- Sine Input ---

RLBase.reset!(env)

output = Vector{Float64}()
f0 = 50
V_eff = 230 * sqrt(2)
for i = 1:50
    local v_sin1 = V_eff * sin.(2*pi * f0 * i / 10_000)
    local v_sin2 = V_eff * sin.(2*pi * f0 * i / 10_000 + 1)

    uuu = [v_sin1, v_sin2]
    uuu /= 300.0
    env(u[:,i])
    append!(output, env.state[2])
end

plot(output)