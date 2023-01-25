#using Dare

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
# Machine parameters

sampling = 1 # machine retain 1 in every so many samples

history = 22*12 # backward trajectory [months]
future = 22*12 # forward trajectory [months]

scale = 1 #bandwidth, i.e. Kernel Scale

#-------------------------------------------------------------------------------
# Machine Perspective

npast = convert(Int, round(history/sampling)) #Past Series sample size
nfuture = convert(Int, round(future/sampling)) #Future series sample size

window_size = npast + nfuture

#******* or read from CSV
df = CSV.read(joinpath(pwd(), "SN_m_tot_V2.0.csv"), DataFrame, header = false)
year = df.Column1
month = df.Column2
data = df.Column3 # sunspots
#*******

N = length(data)# number of samples
#-------------------------------------------------------------------------------
# Emergence - Pattern Discovery

println("\n1. Generating gram matrices")

Gx, Gy, index_map = series_Gxy([data], scale, npast, nfuture)

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
eigenvalues, basis, coords, info = Spectral_Basis(Gs, num_basis = 10, scaled = true)

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

pred, dist = Predict(2*nfuture, coords[end - nfuture, :], shift_op, expect_op, return_dist = 2, knn_convexity = 4, coords = coords)
final_dist = dist[:, end]

#-------------------------------------------------------------------------------
# Plots

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
nans= fill!(nans, NaN)
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