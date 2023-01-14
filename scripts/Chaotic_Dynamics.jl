using DrWatson
@quickactivate "dare"

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

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
u0 = [-1.0, -0.3, 0.0] # initial conditions
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

D = 10

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

u = reduce(hcat, sol.u)
df = DataFrame(u1 = u[1,:], u2 = u[2,:])

plot(df, x = :u1, y = :u2)
#display(traj)

#traj = plot(sol, idxs = (1, 2))

colrs = Array{Symbol, 1}(undef, Deus.k)
for i in 1:Deus.k
    
    if Deus.s[i] == 0
        colrs[i] = :red
    elseif Deus.s[i] == 1
        colrs[i] = :blue
    elseif Deus.s[i] == 2
        colrs[i] = :green
    end
end


#= traj = plot(Deus.x[1, 500:end], Deus.x[2, 500:end], 
            seriestype = :scatter,
            markercolor = colrs[500:end],
            markersize = 1,
            markerstrokewidth = 0,
            legend = false)
display(traj) =#

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
