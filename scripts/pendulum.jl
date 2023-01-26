using Dare
using Distributions
using DifferentialEquations
using CSV
using DataFrames
using PlotlyJS

#= using DrWatson
@quickactivate "dare"

using CSV
using DataFrames
using DifferentialEquations
using Distributions
using LinearAlgebra
using JuMP
import Ipopt
using KrylovKit #GenericSchur is another option
using NearestNeighbors
using NonNegLeastSquares
using PlotlyJS =#

#= include(srcdir("Kernel_Machine.jl"))
include(srcdir("Machine_Dynamics.jl"))
include(srcdir("Dif_Map.jl"))
include(srcdir("Dynamical_Systems.jl")) =#

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#-------------------------------------------------------------------------------
# Damped Pendulum Dynamics

μ_s = 0.005 # time step size (seconds)
t_final = 30 # max simulation time (seconds)

# Other input constants 
m = 1 # mass (kg)
g = 9.81 # gravity (m/s^2)
L = 1 # length (m)
b = 0.1 # damping value (kg/m^2-s)

# initial angle (rad) is drawn uniformly in this range
θ₀_range = (π/180)*[30, 150]

# same for initial angular velocity (rad/s)
ω₀_range = [1, 2]

#-------------------------------------------------------------------------------
# Machine parameters

sampling = 5 # machine retain 1 in every so many samples

history = 5 # backward trajectory [seconds]
future = 5 # forward trajectory [seconds]

scale = 1 #bandwidth, i.e. Kernel Scale

#-------------------------------------------------------------------------------
# Damped Pendulum Dynamics - time simulation

tₛ = 0:μ_s:t_final # system time

u0 = [rand(Uniform(θ₀_range[1], θ₀_range[2])), rand(Uniform(ω₀_range[1], ω₀_range[2]))] # initial conditions
u0 = [π/4, 0.0]
tspan = (0.0, t_final)

#Parameters
p = [m, g, L, b]
prob = ODEProblem(Pendulum!, u0, tspan, p)

sol = solve(prob, saveat = μ_s, wrap = Val(false))
u = reduce(hcat, sol.u)

#-------------------------------------------------------------------------------
# Machine Perspective

μ_m = sampling*μ_s # machine time step
npast = convert(Int, round(history/μ_m)) #Past Series sample size
nfuture = convert(Int, round(future/μ_m)) #Future series sample size

window_size = npast + nfuture

tₘ = 0:μ_m:t_final # machine time
xₘ = sol(tₘ) # machine samples of the continuous system
uₘ = reduce(hcat, xₘ.u)

# positions - no velocity, not a 1 to 1 mapping of state space
# => require history for reconstruction
x = L * sin.(uₘ[1,:])
y = L .- L.*cos.(uₘ[1,:])
tₘ = xₘ.t

data = [vec(x), vec(y)]

#******* Writing to CSV
#= dateframe_x_y_t = DataFrame(x = x, y = y, t = tₘ)
CSV.write("pendulum.csv", dateframe_x_y_t) =#
#*******

#******* or read from CSV
#= df = CSV.read(joinpath(pwd(), "pendulum.csv"), DataFrame)
data = [df.x[1:end], df.y[1:end]] =#
#*******

N = length(data[1])# number of samples
#-------------------------------------------------------------------------------
# Emergence - Pattern Discovery

@time begin

    println("\n1. Generating gram matrices")

    Gx, Gy, index_map = series_Gxy(data, scale, npast, nfuture)

    #= println("Gx = ")
    display(Gx)
    println("Gy = ")
    display(Gy) =#

    # Compute the state similarity matrix. See the paper
    # Embedding to get the similarity matrix between conditional distributions
    println("\n2. Computing Gs")
    Gs = Embed_States(Gx, Gy)

    #= println("Gs = ")
    display(Gs) =#

    # Compute a spectral basis for representing the causal states.
    # Find a reduced dimension embedding and extract the significant coordinates"
    println("\n3. Projection")
    eigenvalues, basis, coords, info = Spectral_Basis(Gs, num_basis = 15, scaled = true)

    #= println("eigenvalues = ")
    display(eigenvalues)
    println("basis = ")
    display(basis)
    println("coords = ")
    display(coords) =#

    # This is the forward operator in state space. It is built from consecutive
    # indices in the index map. Data series formed by multiple contiguous time
    # blocks are supported, as well as the handling of NaN values
    println("\n4. Forward Shift Operator")
    #= eigenvalues[2] = 5
    eigenvalues[3] = 6
    coords[:,2] = collect(1:6)
    coords[:,3] = collect(100*(3:8)) =#
    shift_op = Shift_Operator(coords, eigenvalues, index_map = index_map)

    #= println("shift_op = ")
    display(shift_op) =#

    # This is the expectation operator, using its default function that predicts
    # the first entry in the future sequence from the current state distribution. 
    # You can specify other functions, see the documentation
    println("\n5. Expectation Operator")
    expect_op = Expectation_Operator(coords, index_map, data)

    #= println("expect_op = ")
    display(expect_op) =#

    # Start from the last known point (represented by its coordinates) and
    # evolve the state for nfuture+1 points.
    println("\n6. Prediction")
    #= pred, dist = Predict(2*nfuture, coords[end - nfuture, :], shift_op, expect_op, return_dist = 2)
    final_dist = dist[:, end] =#

    pred, dist = Predict(2*nfuture, coords[end - nfuture, :], shift_op, expect_op, return_dist = 2, knn_convexity = 1, coords = coords)
    final_dist = dist[:, end]
end
#-------------------------------------------------------------------------------
# Plots

# Exact State space (θ, ω)
df_θ_ω = DataFrame(θ = u[1,:], ω = u[2,:])
plot_θ_ω = plot(df_θ_ω, x = :θ, y = :ω,
                Layout(
                    title = attr(
                        text = "Damped Pendulum: Exact State Space",  
                    ),
                    title_x = 0.5,
                    xaxis_title = "θ [rad]",
                    yaxis_title = "ω [rad/s]",)
                )

#display(plot_θ_ω)    

# (x, y) positions from machine perspective
nans = Array{Float64, 1}(undef, N - 2*nfuture)
nans= fill!(nans, NaN)

x̂ = vec([nans; pred[1]])
ŷ = vec([nans; pred[2]])

df_x_y_t = DataFrame(x = x, y = y, x̂ = x̂, ŷ = ŷ, t = tₘ)
plot_x_y = plot(df_x_y_t, 
                x = :x, y = :y,
                Layout(
                    title = attr(
                        text = "Damped Pendulum: Cartesian Coordinates",  
                    ),
                    title_x = 0.5,
                    xaxis_title = "x [m]",
                    yaxis_title = "y [m]",),
                marker = attr(
                    color = :t,
                    colorscale = "Cividis",
                    size = 7),
                mode = "markers")

#display(plot_x_y)

trace_x = scatter(df_x_y_t, x = :t, y = :y, name = "y")
trace_y = scatter(df_x_y_t, x = :t, y = :x, name = "x")
trace_ŷ = scatter(df_x_y_t, x = :t, y = :ŷ, name = "y'")
trace_x̂ = scatter(df_x_y_t, x = :t, y = :x̂, name = "x'")

plot_x_t = plot([trace_x, trace_y, trace_x̂, trace_ŷ],
                Layout(
                    title = attr(
                        text = "Damped Pendulum: Evolution in Time",     
                        ),
                    title_x = 0.5,
                    xaxis_title = "t [s]",
                    yaxis_title = "y,x [m]",
                    ),  
                )
display(plot_x_t)

N₁ = length(coords[:,2])
N₂ = length(dist[2,:])
nans = Array{Float64, 1}(undef, N₁ - N₂)
nans = fill!(nans, NaN)
Φ₁ = vec([nans; dist[2, :]])
Φ₂ = vec([nans; dist[3, :]])
Φ₃ = vec([nans; dist[4, :]])

df_Ψ_Φ = DataFrame(Ψ₁ = coords[:,2], Ψ₂ = coords[:,3], Ψ₃ = coords[:,4], Φ₁ = Φ₁, Φ₂ = Φ₂, Φ₃ = Φ₃)

#= trace_Ψ = scatter(df_Ψ_Φ, x = :Ψ₁, y = :Ψ₂, name = "Ψ")
trace_Φ = scatter(df_Ψ_Φ, x = :Φ₁, y = :Φ₂, name = "Φ")

plot_Ψ_Φ = plot([trace_Ψ, trace_Φ],
                Layout(
                    title = attr(
                        text = "Damped Pendulum: Reconstructed State Space",  
                    ),
                    title_x = 0.5,
                    xaxis_title = "Ψ₁",
                    yaxis_title = "Ψ₂",),
                )

display(plot_Ψ_Φ) =#

df_Ψ_Φ = DataFrame(Ψ₁ = coords[:,2], Ψ₂ = coords[:,3], Ψ₃ = coords[:,4], Φ₁ = Φ₁, Φ₂ = Φ₂, Φ₃ = Φ₃)

trace_Ψ = scatter3d(df_Ψ_Φ, x = :Ψ₁, y = :Ψ₂, z = :Ψ₃, name = "Ψ", mode = "lines")
trace_Φ = scatter3d(df_Ψ_Φ, x = :Φ₁, y = :Φ₂, z = :Φ₃, name = "Φ", mode = "lines")

plot_Ψ_Φ_3d = plot([trace_Ψ, trace_Φ],
                Layout(
                    title = attr(
                        text = "Damped Pendulum: Reconstructed State Space",  
                    ),
                    title_x = 0.5,
                    scene = attr(
                        xaxis_title = "Ψ₁",
                        yaxis_title = "Ψ₂",
                        zaxis_title = "Ψ₃",  
                    ),
                    ),
                )

display(plot_Ψ_Φ_3d)

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")