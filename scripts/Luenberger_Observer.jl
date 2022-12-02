using SpecialFunctions
using LinearAlgebra
using Plots
using StatsBase
using LinearAlgebra
using Combinatorics

function charpoly_coef(λ)

    # given the roots, this function finds the coefficients

    n = length(λ)
    ds = 1:n

    λ = -1*λ

    α = Array{Float64, 1}(undef, length(λ))

    for i in 1:length(λ)

        x = combinations(ds, n-i+1)
        y = collect(x)

        α[i] = 0.0

        for j in eachindex(y)

            α[i] = α[i] + prod(λ[y[j]])
        end
    end

    return α
end

function Observability(C, A; n = size(A,1))

    O = Array{Float64, 2}(undef, n*size(C,1), size(A,2))

    iter = size(C, 1)
    dim = 1:iter

    O[dim, :] = C

    for i in 2:n

        O[dim .+ iter, :] = O[dim, :]*A
        dim = dim .+ iter
    end

    return O, rank(O)
end

function Ackermann_Gain_Matrix(A, C, λ; n = size(A, 1))

    #= Theory
        For a single-output, observable system (A, C) and the desired closed-loop 
        eigenvalues being the roots of the characteristic polynomial αd_A the 
        feedback of the constant gain matrix can be found.
    =#

    αd_A = Array{Float64, 2}(undef, size(A,1), size(A,2))

    α = charpoly_coef(λ)

    αd_A = α[1]*I

    for i in 2:n
        αd_A = αd_A + α[i]*(A^(i-1))
    end

    αd_A = αd_A .+ A^n

    O, _ = Observability(C, A)

    output = zeros(n)
    output[end] = 1.0

    K = αd_A*inv(O)*output

    return K
end

function Mulit_Gain_Matrix(A, C, λ, p; n = size(A, 1))

    #given parameter vectors 

    v = Array{Float64, 2}(undef, size(A,1), size(A,2))

    for i in 1:n

        v[:, i] = -transpose(p[:, i])*C*inv(λ[i]*I - A)
        #v[:, i] = v[:, i]/sum(v[:, i])

    end

    K = -inv(transpose(v))*transpose(p)

    return K, v
end

A = Array{Float64, 2}(undef, 3, 3)
C = Array{Float64, 2}(undef, 1, 3)

A = [0 1 0;
    0 0 1;
    0.0 2 -1]

C = [1.0 0 0]

O, r = Observability(C, A)

eigvals(A)

λ = [-10 -10 -10]
p = [2.0 0 1;
     0 2 0]

α = charpoly_coef(λ)

#K, v = Mulit_Gain_Matrix(A, C, λ, p; n = size(A, 1))

K = Ackermann_Gain_Matrix(A, C, λ)

#= 
i = 3

#= v = Array{Float64, 2}(undef, size(A,1), size(A,2))

v[:, i] = -transpose(p[:, i])*C*inv(λ[i]*I - A)
v[:, 1] = [0.2; 0.6; 0.2]
v[:, 2] = [0.07; 0.2; 0.73]
v[:, 3] = [0.03; 0.15; 0.03] =#

for i in 1:3
    p[:, i] = -inv(C*transpose(C))*C*(λ[i]*I - A)*v[:, i]
end
p =#

#K = -inv(transpose(v))*transpose(p)

#= A = [0 20.6;
    1 0]

C = [0.0 1]

λ = [-10 -10]

K = Ackermann_Gain_Matrix(A, C, λ) =#