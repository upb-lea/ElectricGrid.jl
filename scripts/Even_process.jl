using Dare

using Distributions
using DifferentialEquations
using CSV
using DataFrames
using PlotlyJS

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........")

npast = 10 #Past Series size
nfuture = 10 #Future series size

N = 8000 #Number of training samples

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

@time begin

    println("\n1. Generating gram matrices")
    Gx, Gy, index_map = series_Gxy([series], scale, npast, nfuture)

    # Compute the state similarity matrix.
    println("\n2. Computing Gs")
    Gs = Embed_States(Gx, Gy)

    # Compute a spectral basis for representing the causal states.
    println("\n3. Projection")
    eigenvalues, basis, coords = Spectral_Basis(Gs, num_basis = 2, scaled = false)

end

df = DataFrame(Ψ₁ = coords[:,2])
p = plot(df, x = :Ψ₁, kind = "histogram", nbinsx = 100, histnorm = "probability density")

display(p)

print("\n...........o0o----ooo0o0ooo~~~   END   ~~~ooo0o0ooo----o0o...........\n")