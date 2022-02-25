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


env = SimEnv(A=A, B=B, C=C)

