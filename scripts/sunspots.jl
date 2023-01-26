using Dare

using DataFrames
using Distributions
using DifferentialEquations
using DSP
using CSV
using PlotlyJS

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

df = CSV.read(joinpath(pwd(), "SN_m_tot_V2.0.csv"), DataFrame, header = false)
year = df.Column1
month = df.Column2
data = df.Column4[1:sampling:end] # sunspots

N = length(data)# number of samples
t = 1:1:N

#-------------------------------------------------------------------------------
# Filtering




#-------------------------------------------------------------------------------
# Emergence - Pattern Discovery

#= #= 1. Generating gram matrices =#
Gx, Gy, index_map = series_Gxy([data], scale, npast, nfuture)

#= 2. Computing Gs
    Compute the state similarity matrix.
    Embedding to get the similarity matrix between conditional distributions 
=#
Gs = Embed_States(Gx, Gy)

#= 3. Projection
    Compute a spectral basis for representing the causal states.
    Find a reduced dimension embedding and extract the significant coordinates
=#
eigenvalues, basis, coords, info = Spectral_Basis(Gs, num_basis = 10, scaled = true)

#= 4. Forward Shift Operator
    This is the forward operator in state space. It is built from consecutive
    indices in the index map. Data series formed by multiple contiguous time
    blocks are supported, as well as the handling of NaN values 
=#
shift_op = Shift_Operator(coords, eigenvalues, index_map = index_map)

#= 5. Expectation Operator
    This is the expectation operator, using its default function that predicts
    the first entry in the future sequence from the current state distribution. 
    You can specify other functions, see the documentation 
=#
expect_op = Expectation_Operator(coords, index_map, data)

#= 6. Prediction
    Start from the last known point (represented by its coordinates) and
    evolve the state for nfuture+1 points. 
=#
pred, dist = Predict(2*nfuture, coords[end - nfuture, :], shift_op, expect_op, return_dist = 2, knn_convexity = 4, coords = coords)
final_dist = dist[:, end] =#

#-------------------------------------------------------------------------------
# Plots

# (x, y) positions from machine perspective
nans = Array{Float64, 1}(undef, N - 2*nfuture)
nans = fill!(nans, NaN)

x̂ = vec([nans; pred[1]])
ŷ = vec([nans; pred[2]])

df_sunspots = DataFrame(months = t, sunspots = data)
plot_sun_t = plot(df_sunspots, 
                x = :months, y = :sunspots,
                Layout(
                    title = attr(
                        text = "Sunspots monthly series from SILSO",  
                    ),
                    title_x = 0.5,
                    xaxis_title = "Months",
                    yaxis_title = "Number of Sunspots",),
                mode = "lines")

display(plot_sun_t)



print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")