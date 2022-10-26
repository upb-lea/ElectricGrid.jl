using Plots

include("Complexity.jl")

#_______________________________________________________________________________
# Parameters - Time simulation

ω = 1
t_final = 15e3*2π/ω
μ_s = 2π/ω # time step in s

#_______________________________________________________________________________
# Environment Calcs

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
p = [1, 0.25, -1, 1, 0.4, ω] #ω = 1
#p = [1, 0.02, 1, 5, 8, ω] #ω = 0.5

prob = ODEProblem(Driven_Duffing!, u0, tspan, p);

#_______________________________________________________________________________
#%% Starting time simulation
sol = solve(prob)

periods = floor(tspan[2]*p[6]/2π)
section = (2π/(p[6]))*(0:1:periods)
Poincare = sol(section)

x = Poincare[1:2, :]

#_______________________________________________________________________________
#%% Emergence

D = 14

ϵ = Array{Float64, 1}(undef, dim) # Instrument resolution
ep = sum(x[1:1, :])/N
ϵ[:] = [0.5 for i in 1:dim]

x_range = Array{Float64, 2}(undef, dim, 2) # State space limits
x_range[1, 1] = maximum(Poincare[1, :])
x_range[1, 2] = minimum(Poincare[1, :])
x_range[2, 1] = maximum(Poincare[2, :])
x_range[2, 2] = minimum(Poincare[2, :])

Deus = ϵ_Machine(N, D, ϵ[1:2], x_range[1:2, :], μ_m, μ_s, δ = 0.05)

Cranking(Deus, Poincare[1:2, :], μ_s)

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
