using Plots
using DifferentialEquations

include("Complexity.jl")

#_______________________________________________________________________________
# Parameters - Time simulation

t_final = 220000.0
Timestep = 10 # time step in μs

#_______________________________________________________________________________
# Environment Calcs

N_s = (1/(Timestep*1e-6)) # time intervals
μ_s = 1/N_s # time step
N = convert(Int64, floor(t_final/μ_s)) + 1

t = 0:μ_s:t_final # time

#_______________________________________________________________________________
# Machine Parameters

μ_m = μ_s #sampling timestep

#_______________________________________________________________________________
# Dynamical System

dim = 2 # dimensions of Poincare section

N = convert(Int64, floor(t_final/μ_s)) + 1

#Driven_Duffing!
u0 = [1.0, 0.0, 0.0] # initial conditions
tspan = (0.0, t_final)

#Parameters
#p = [m, δ, α, β, γ, ω]
p = [1, 0.25, -1, 1, 0.4, 1]
#p = [1, 0.02, 1, 5, 8, 0.5]

prob = ODEProblem(Driven_Duffing!, u0, tspan, p);

#_______________________________________________________________________________
#%% Starting time simulation
sol = solve(prob)

periods = floor(tspan[2]*p[6]/2π)
section = (2π/(p[6]))*(0:1:periods)
Poincare = sol(section)

#_______________________________________________________________________________
#%% Plots

traj = plot(sol, idxs = (1, 2))
display(traj)

#traj = plot(sol, idxs = (1, 2))
traj = plot(Poincare[1, :], Poincare[2, :], 
            seriestype = :scatter,
            markercolor = :black,
            markersize = 1,
            markerstrokewidth = 0,
            legend = false)
display(traj)
