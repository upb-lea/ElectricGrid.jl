#using DrWatson
#@quickactivate "dare"

using Plots
using DifferentialEquations
using VoronoiCells
using GeometryBasics

include("Complexity.jl")

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#_______________________________________________________________________________
# Parameters - Time simulation

Timestep = 10 # time step in μs
t_final = 1.3 # 0.75 # time in seconds, total simulation run time
fsys = 2000 # Hz, fundamental frequency of system
fsys = 1/(10e-6) 

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
dim = 2 # dimensions of state space
x = Array{Float64, 2}(undef, dim, N) # State space
λ = Array{Float64, 1}(undef, dim) # Lyapunov exponents

# Misiurewicz point:
r = 3.9277370017867516
# Accumulation board - onset of chaos
#r = 3.5699456718695445
r = 3.56

x[1, 1] = 0.4
x[2, 1] = 0.41

#_______________________________________________________________________________
#%% Starting time simulation
println("\nHere we go.\n")

@time begin

    println("Progress : 0.0 %")

    for i in 1:N-1

        global λ

        # Progress Bar
        if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
            flush(stdout)
            println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
        end

        x[1, i + 1], λ[1] = Logistic_Map(x[1, i], λ[1], r = r)
        x[2, i + 1], λ[2] = Logistic_Map(x[2, i], λ[2], r = r)

    end

    println("Progress : 100.0 %\n")

    λ = λ./N

    D = 14

    ϵ = Array{Float64, 1}(undef, dim) # Instrument resolution
    ep = sum(x[1:1, :])/N
    ϵ[:] = [0.5 for i in 1:dim]

    x_range = Array{Float64, 2}(undef, dim, 2) # State space limits
    x_range[:, 1] = [1.0 for i in 1:dim] # maximum
    x_range[:, 2] = [0.0 for i in 1:dim] # minimum

    Deus = ϵ_Machine(N, D, ϵ[1:1], x_range[1:1, :], μ_m, μ_s, δ = 0.05)

    Cranking(Deus, x[1:1, :], μ_s)
end

T_plot_end = 30

N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/N_s)*N_s))
N_range = 1:N_plot_end + 2
Nm_plot_end = convert(Int64, round((T_plot_end/fsys  - μ_m)*1/μ_m))
Nm_range = 1:Nm_plot_end + 2

p1 = plot(t[N_range], x[1, N_range], label = "x")
#p1 = plot!(Deus.t_m[Nm_range], Deus.s[Nm_range])
p1 = plot!(Deus.t_m[Nm_range], Deus.x_m[1, Nm_range], label = "xm")
#plot!(N_range, x[1, N_range])
#plot!(Nm_range, Deus.s[Nm_range])
#plot!(Nm_range, Deus.x_m[1, Nm_range])

#display(p1)

println("Cμ_t[end] = ", Deus.Cμ_t[end])
println("Hα[end] = ", Deus.Hα[end])
println("Period = ", 2^(Deus.Hα[end]))

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")

