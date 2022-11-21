using SpecialFunctions
using LinearAlgebra
using Plots
using StatsBase
using LinearAlgebra
using AbstractAlgebra

A = Array{Float64, 2}(undef, 3, 3)
C = Array{Float64, 2}(undef, 1, 3)

function charpoly_coef(λ)

    α = Array{Float64, 1}(undef, length(λ))

    α[1] = sum(λ)
    α[end] = prod(λ)

    return α
end

A = [0 1 0;
    0 0 1;
    0 2 -1]

C = [1 0 0]

λ1 = -10
λ2 = -10
λ3 = -10

B = matrix(ZZ, A)
Zx, x = ZZ["x"]
f = charpoly(Zx, B)

eigvals(A)

λ = [10 10 10]