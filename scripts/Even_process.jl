using DrWatson
@quickactivate "dare"

using PlotlyJS
using DataFrames

include(srcdir("Kernel_Machine.jl"))
include(srcdir("Dif_Map.jl"))

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

npast = 10 #Past Series size
nfuture = 5 #Future series size

N = 5000 #Number of training samples

scale = 1 #bandwidth

#-------------------------------------------------------------

window_size = npast + nfuture

series_length = N + window_size - 1

series = zeros(series_length)

state = rand(0:1)

for i in 1:series_length

    global state

    if state == 1

        series[i] = 1
        state = 0

    else

        series[i] = rand(0:1)
        state = series[i]
    end
end

println("\n1. Generating gram matrices")
Gx, Gy, index_map = series_Gxy([series], scale, npast, nfuture)

#= println("Gx = ")
display(Gx)
println("Gy = ")
display(Gy) =#

# Compute the state similarity matrix. See the paper
println("\n2. Computing Gs")
Gs = embed_states(Gx, Gy)

#= println("Gs = ")
display(Gs) =#

# Compute a spectral basis for representing the causal states. See the paper
println("\n3. Projection")
eigenvalues, basis, coords, info = spectral_basis(Gs, num_basis = 2, scaled = false)

#= println("eigenvalues = ")
display(eigenvalues)
println("basis = ")
display(basis)
println("coords = ")
display(coords) =#

df = DataFrame(Ψ₁ = coords[:,2])
p = plot(df, x = :Ψ₁, kind = "histogram", nbinsx = 100, histnorm = "probability density")

display(p)

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")