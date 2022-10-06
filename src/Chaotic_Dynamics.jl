using Plots
using DifferentialEquations

include("Complexity.jl")

# Lorenz!
u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 100.0)
p = [10.0, 28.0, 8/3] 

#Driven_Duffing!
u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 220000.0)
#p = [m, δ, α, β, γ, ω]
p = [1, 0.25, -1, 1, 0.4, 1]
p = [1, 0.02, 1, 5, 8, 0.5]
#p = [1, 0.3, -1, 1, 0.5, 1.2]

prob = ODEProblem(Driven_Duffing!, u0, tspan, p);
sol = solve(prob);

periods = floor(tspan[2]*p[6]/2π)
section = (2π/(p[6]))*(0:1:periods)
Poincare = sol(section)

#traj = plot(sol, idxs = (1, 2))
traj = plot(Poincare[1, :], Poincare[2, :], 
            seriestype = :scatter,
            markercolor = :black,
            markersize = 1,
            markerstrokewidth = 0,
            legend = false)
display(traj)
