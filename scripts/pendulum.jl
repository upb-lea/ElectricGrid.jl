using DrWatson
@quickactivate "dare"

using PlotlyJS
using DataFrames
using Distributions
using CSV

include(srcdir("Kernel_Machine.jl"))
include(srcdir("Machine_Dynamics.jl"))
include(srcdir("Dif_Map.jl"))
include(srcdir("Dynamical_Systems.jl"))

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
u0 = [(π/180)*50, 1.5]
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
N = length(tₘ)# number of samples
xₘ = sol(tₘ) # machine samples of the continuous system
uₘ = reduce(hcat, xₘ.u)

# positions - no velocity, not a 1 to 1 mapping of state space
# => require history for reconstruction
x = L * sin.(uₘ[1,:])
y = L .- L.*cos.(uₘ[1,:])
tₘ = xₘ.t

data = [vec(x), vec(y)]

#-------------------------------------------------------------------------------
# Emergence

#= #******* for debugging
data = [vec(x[1:7]), vec(y[1:7])]
npast = 1
nfuture = 1
#******* =#

#******* or read from CSV
df = CSV.read(joinpath(pwd(), "pendulum.csv"), DataFrame)
data = [df.x[1:end], df.y[1:end]]
N = length(data[1])
#*******

println("\n1. Generating gram matrices")

Gx, Gy, index_map = series_Gxy(data, scale, npast, nfuture)

#= println("Gx = ")
display(Gx)
println("Gy = ")
display(Gy) =#

# Compute the state similarity matrix. See the paper
# Embedding to get the similarity matrix between conditional distributions
println("\n2. Computing Gs")
Gs = embed_states(Gx, Gy)

#= println("Gs = ")
display(Gs) =#

# Compute a spectral basis for representing the causal states. See the paper
# Find a reduced dimension embedding and extract the first two coordinates"
println("\n3. Projection")
eigenvalues, basis, coords, info = spectral_basis(Gs, num_basis = 3, scaled = true)

#= println("eigenvalues = ")
display(eigenvalues)
println("basis = ")
display(basis)
println("coords = ")
display(coords) =#

# This is the forward operator in state space. It is built from consecutive
# indices in the index map. Data series formed by multiple contiguous time
# blocks are supported, as well as the handling of NaN values
println("\n4. Shift Operator")
eigenvalues, basis, coords, info = spectral_basis(Gs, num_basis = 3, scaled = true)
#= eigenvalues[2] = 5
eigenvalues[3] = 6
coords[:,2] = collect(1:6)
coords[:,3] = collect(100*(3:8)) =#
shift_op = shift_operator(coords, eigenvalues, index_map = index_map)

println("shift_op = ")
display(shift_op)

# This is the expectation operator, using its default function that predicts
# the first entry in the future sequence from the current state distribution. 
# You can specify other functions, see the documentation
println("\n5. Expectation Operator")
expect_op = expectation_operator(coords, index_map, data)

#= println("expect_op = ")
display(expect_op) =#

# Start from the last known point (represented by its coordinates) and
# evolve the state for nfuture+1 points.
println("\n6. Prediction")
pred, dist = predict(2*nfuture + 1, coords[end - nfuture, :], shift_op, expect_op, return_dist = 2)
final_dist = dist[:, end]

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

display(plot_θ_ω)    

# (x, y) positions from machine perspective
nans = Array{Float64, 1}(undef, N - 2*nfuture - 1)
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

display(plot_x_y)

trace_x = scatter(df_x_y_t, x = :t, y = :y, name = "y")
trace_y = scatter(df_x_y_t, x = :t, y = :x, name = "x")
trace_x̂ = scatter(df_x_y_t, x = :t, y = :x̂, name = "x'")
trace_ŷ = scatter(df_x_y_t, x = :t, y = :ŷ, name = "y'")

plot_x_t = plot([trace_x, trace_y, trace_x̂, trace_ŷ],
                Layout(
                    title = attr(
                        text = "Damped Pendulum: Cartesian Coordinates",  
                        
                        ),
                    title_x = 0.5,
                    xaxis_title = "t [s]",
                    yaxis_title = "y,x [m]",
                    ),  
                )
display(plot_x_t)

df_Ψ₁_Ψ₂ = DataFrame(Ψ₁ = coords[:,2], Ψ₂ = coords[:,3])
plot_Ψ₁_Ψ₂ = plot(df_Ψ₁_Ψ₂, x = :Ψ₁, y = :Ψ₂,
                Layout(
                    title = attr(
                        text = "Damped Pendulum: Reconstructed State Space",  
                    ),
                    title_x = 0.5,
                    xaxis_title = "Ψ₁ [k₁*rad]",
                    yaxis_title = "Ψ₂ [k₂*rad/s]",),
                    line = attr(
                        width = 2,
                        color = "firebrick"),
                    mode = "line",
                )

display(plot_Ψ₁_Ψ₂)

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")