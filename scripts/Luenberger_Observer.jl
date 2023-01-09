using LinearAlgebra
using StatsBase
using Combinatorics

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

function Multi_Gain_Matrix_par(A, C, λ, p)

    n = size(A, 1)
    v = Array{Float64, 2}(undef, size(A,1), size(A,2))

    for i in 1:n

        v[:, i] = -transpose(p[:, i])*C*inv(I*λ[i] - A)
    end

    K = inv(transpose(v))*transpose(p)

    return K, v
end

function Multi_Gain_Matrix_vec(A, C, λ; v = Matrix(I, size(A,1), size(A,2)))

    n = size(A, 1)
    p = Array{Float64, 2}(undef, size(C,1), size(C,2))

    for i in 1:n
        p[:, i] = -inv(C*transpose(C))*C*(λ[i]*I - transpose(A))*v[:, i]
    end  

    K = inv(transpose(v))*transpose(p)

    return K, p
end

function Feed_Gain_Matrix_par(A, B, λ, p)

    n = size(A, 1)
    v = Array{Float64, 2}(undef, size(A,1), size(A,2))

    for i in 1:n

        v[:, i] = -inv(I*λ[i] - A)*B*p[:, i]
    end

    F = p*inv(v)

    return F, v
end

function Feed_Gain_Matrix_vec(A, B, λ; v = Matrix(I, size(A,1), size(A,2)))

    n = size(A, 1)
    p = Array{Float64, 2}(undef, size(B,2), size(B,1))

    for i in 1:n

        p[:, i] = -inv(transpose(B)*B)*transpose(B)*(λ[i]*I - A)*v[:, i]
    end

    F = p*inv(v)

    return F, p
end

function Gram_Schmidt(V)

    # The Gram-Schmidt process return the vectors as an orthonormal set

    N = size(V,1) # number of vectors in V
    R = similar(V, Float64)
    R = fill!(R, 0)

    for i in 1:N

        for j in i:-1:1

            R[:,i] = R[:,i] + dot(V[:, i], R[:, j])*R[:, j]

        end

        R[:,i] = normalize(V[:,i] - R[:,i])

    end

    return R
end

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

A = [1. -1 0;
     -1 2 -1;
     0. -1 1;]

C = [0. 1. 0.;
     0. 0. 1.]

p = [2. 0. 1.;
     0. 2. 0.]

λ = [-2 -2 -10]

_, r = Observability(C, A)

K, _ = Multi_Gain_Matrix_par(A, C, λ, p)

println()
println("K = ", round.(K, digits = 3))
println("v = ", round.(v, digits = 3))
println("p = ", round.(p, digits = 3))
println("λ = ", round.(eigvals(A - K*C), digits = 3))

A = [0 1 0;
    0 0 1;
    0.0 2 -1]

C = [1.0 0 0]

λ = [-10. -10 -10]

K = Ackermann_Gain_Matrix(A, C, λ)

#= println()
println("K = ", round.(K, digits = 3))
println("λ = ", real(round.(eigvals(A - K*C), digits = 3))) =#

A = [0. 1 0 ;
     0 0 1 ;
     0. 2 -1];

B = [0. 1.;
     1. 1.;
     0. 0.]

p = [1. 0. 1.;
     0. 1. 1.]

λ = [-5. -1 -3]

F, v = Feed_Gain_Matrix_par(A, B, λ, p)

#= println()
println("F = ", round.(F, digits = 3))
println("v = ", round.(v, digits = 3))
println("p = ", round.(p, digits = 3))
println("λ = ", real(round.(eigvals(A - B*F), digits = 3))) =#

A = [23.0    1.0  262.0    8.0;
    -1.0   23.0   -8.0  262.0;
    -20.0   -1.0   40.0    1.0;
    1.0  -20.0   -1.0   40.0]

C = [10.0   0.0  -36.0   -1.0;
    -0.0  10.0    1.0  -36.0]

p = [2.5 8. 1. 1.5;
     0.5 2. 9. 1.2]

λ = [-0.0 -0.0 -0.01 -0.01]

K, _ = Multi_Gain_Matrix_par(A, C, λ, p)

#= println()
println("K = ", round.(K, digits = 3))
println("v = ", round.(v, digits = 3))
println("p = ", round.(p, digits = 3))
println("λ = ", real.(round.(eigvals(A - K*C), digits = 3))) =#

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")