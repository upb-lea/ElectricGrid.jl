using DrWatson
@quickactivate "dare"

using PlotlyJS
using DataFrames
using Distributions
using Base: llvmcall

include(srcdir("Kernel_Machine.jl"))
include(srcdir("Dif_Map.jl"))
include(srcdir("Dynamical_Systems.jl"))

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

#-------------------------------------------------------------------------------
# Damped Pendulum Dynamics

μ_s = 0.005 # time step size (seconds)
t_final = 30      # max simulation time (seconds)

# Other input constants 
m = 1 # mass (kg)
g = 9.81 # gravity (m/s^2)
L = 1 # length (m)
b = 0.1 # damping value (kg/m^2-s)

# initial angle (rad) is drawn uniformly in that range
θ₀_range = (π/180)*[30, 150]

# same for initial angular velocity (rad/s)
ω₀_range = [1, 2]

#-------------------------------------------------------------------------------
# Machine parameters

sampling = 5 # machine retain 1 in every so many samples

history = 1 # backward trajectory [seconds]
future = 1 # forward trajectory [seconds]

scale = 1 #bandwidth, i.e. Kernel Scale

#-------------------------------------------------------------------------------
# Damped Pendulum Dynamics - time simulation

tₛ = 0:μ_s:t_final # system time

u0 = [rand(Uniform(θ₀_range[1], θ₀_range[2])), rand(Uniform(ω₀_range[1], ω₀_range[2]))] # initial conditions
tspan = (0.0, t_final)

#Parameters
p = [m, g, L, b]

prob = ODEProblem(Pendulum!, u0, tspan, p)

sol = solve(prob, saveat = μ_s, wrap = Val(true))
u = reduce(hcat, sol.u)
#-------------------------------------------------------------------------------
# Machine Perspective

μ_m = sampling*μ_s # machine time step
npast = convert(Int, round(history/μ_m)) #Past Series sample size
nfuture = convert(Int, round(future/μ_m)) #Future series sample size

tₘ = 0:μ_m:t_final # machine time
N = length(tₘ)# number of samples
xₘ = sol(tₘ) # machine samples of the continuous system
uₘ = reduce(hcat, xₘ.u)

# positions - no velocity, not a 1 to 1 mapping of state space
# => require history for reconstruction
x = L * sin.(uₘ[1,:])
y = L .- L.*cos.(uₘ[1,:])
tₘ = xₘ.t

#-------------------------------------------------------------------------------
# Emergence

println("\n1. Generating gram matrices")
#Gx, Gy, index_map = series_Gxy([x y], scale, npast, nfuture)

#-------------------------------------------------------------------------------
# Plots

# State space (θ, ω)
df = DataFrame(x = u[1,:], y = u[2,:])
plot_θ_ω = plot(df, x = :θ, y = :ω)

# (x, y) positions from machine perspective

df = DataFrame(x = x, y = y, t = tₘ./30)
plot_x_y = plot(df, 
                x = :x, y = :y,
                marker = attr(
                    color = :t,
                    colorscale = "Cividis",
                    size = 7),
                mode = "markers")