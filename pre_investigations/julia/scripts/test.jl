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


# --- SIM ---

#sys_d = ss(Ad,Bd,C,0, ts)
#yout_d, tout_d, xout_d, uout_d = lsim(sys_d,u,t,x0=x0)

env = SimEnv(A=A, B=B, C=C)

RLBase.test_runnable!(env)

# --- PLOT ---

#plot(yout_d'[:,2])