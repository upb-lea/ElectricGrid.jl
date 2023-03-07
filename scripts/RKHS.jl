using Dare
using LinearAlgebra

# number of samples
m = 1000

# regularisation parameter
ε = 10e-4

# data
x = Array{Int64, 1}(undef, m) 
y = Array{Int64, 1}(undef, m)

# feature functions
ϕ = [1 0; 0 1]

# feature matrices
Υ = Array{Float64, 2}(undef, 2, m)
Φ = Array{Float64, 2}(undef, 2, m)

# Gram matrices
Gˣ = Array{Float64, 2}(undef, m, m) 
Gʸ = Array{Float64, 2}(undef, m, m)

# kernel mean
μx = Array{Float64, 1}(undef, 2) 
μy = Array{Float64, 1}(undef, 2)

# kernel (uncentered) cross-covariance
Cxx = Array{Float64, 2}(undef, 2, 2) 
Cyy = Array{Float64, 2}(undef, 2, 2)
Cyx = Array{Float64, 2}(undef, 2, 2) 
Cxy = Array{Float64, 2}(undef, 2, 2) 

# kernel (centered) cross-covariance
Cxxₕ = Array{Float64, 2}(undef, 2, 2) 
Cyyₕ = Array{Float64, 2}(undef, 2, 2)
Cyxₕ = Array{Float64, 2}(undef, 2, 2) 
Cxyₕ = Array{Float64, 2}(undef, 2, 2) 

# kernel conditional embedding
Cy_x = Array{Float64, 2}(undef, 2, 2) 
Cx_y = Array{Float64, 2}(undef, 2, 2) 

# (centered)kernel conditional embedding
Cy_xₕ = Array{Float64, 2}(undef, 2, 2) 
Cx_yₕ = Array{Float64, 2}(undef, 2, 2) 

# kernel conditional mean
μy_x = Array{Float64, 2}(undef, 2, 2) 
μx_y = Array{Float64, 2}(undef, 2, 2)

# (centered) kernel conditional mean
μy_xₕ = Array{Float64, 2}(undef, 2, 2) 
μx_yₕ = Array{Float64, 2}(undef, 2, 2)

# centring matrix
H = I - (1/m)*ones(m)*transpose(ones(m))

# weight vector
ωx = Array{Float64, 2}(undef, m, 2)
ωy = Array{Float64, 2}(undef, m, 2)

# centered weight vector
ωxₕ = Array{Float64, 2}(undef, m, 2)
ωyₕ = Array{Float64, 2}(undef, m, 2)

# kernel matrix
kx = Array{Float64, 2}(undef, m, 2)
ky = Array{Float64, 2}(undef, m, 2)

for i in 1:m

    state_x = rand()
    state_y = rand()
    state_y = state_x

    if state_x < 1/3
        x[i] = 1
    else
        x[i] = 2
    end

    if state_y < 1/2
        y[i] = 1
    else
        y[i] = 2
    end

    Υ[:,i] = ϕ[x[i],:]
    Φ[:,i] = ϕ[y[i],:]
end

μx = (1/m)*sum(Υ, dims = 2)
μy = (1/m)*sum(Φ, dims = 2)

Cxx = (1/m)*Υ*transpose(Υ)
Cyy = (1/m)*Φ*transpose(Φ)
Cyx = (1/m)*Φ*transpose(Υ)
Cxy = (1/m)*Υ*transpose(Φ)

Cxxₕ = (1/m)*Υ*H*transpose(Υ) # (1/m)*(Υ - μx*transpose(ones(m)))*transpose((Υ - μx*transpose(ones(m))))
Cyyₕ = (1/m)*Φ*H*transpose(Φ) # (1/m)*(Φ - μy*transpose(ones(m)))*transpose((Φ - μy*transpose(ones(m))))
Cyxₕ = (1/m)*Φ*H*transpose(Υ) # (1/m)*(Φ - μy*transpose(ones(m)))*transpose((Υ - μx*transpose(ones(m))))
Cxyₕ = (1/m)*Υ*H*transpose(Φ) # (1/m)*(Υ - μx*transpose(ones(m)))*transpose((Φ - μy*transpose(ones(m))))

Gˣ = transpose(Υ)*Υ
Gʸ = transpose(Φ)*Φ

Cy_x = Φ*inv(Gˣ + ε*I)*transpose(Υ) # Cyx*inv(Cxx)
Cx_y = Υ*inv(Gʸ + ε*I)*transpose(Φ) # Cxy*inv(Cyy)

Cy_xₕ = Cyxₕ*inv(Cxxₕ + ε*I) # Φ*inv(H*Gˣ + m*ε*I)*H*transpose(Υ)
Cx_yₕ = Cxyₕ*inv(Cyyₕ + ε*I) # Υ*inv(H*Gʸ + m*ε*I)*H*transpose(Φ)

kx = transpose(Υ)*ϕ
ky = transpose(Φ)*ϕ

ωx = inv(Gˣ + ε*I)*kx
ωy = inv(Gʸ + ε*I)*ky

ωxₕ = inv(H*Gˣ + m*ε*I)*H*kx
ωyₕ = inv(H*Gʸ + m*ε*I)*H*ky

μy_x = Φ*ωx
μx_y = Υ*ωy

μy_xₕ = Φ*ωxₕ
μx_yₕ = Υ*ωyₕ

@show round.(Cy_x - Cyx*inv(Cxx), digits = 5)
@show round.(μy_x[:,1] - Cy_x*ϕ[1,:], digits = 5)
@show round.(μy_x[:,2] - Cy_x*ϕ[2,:], digits = 5)

@show round.(Cx_y - Cxy*inv(Cyy))
@show round.(μx_y[:,1] - Cxy*inv(Cyy)*ϕ[1,:], digits = 5)
@show round.(μx_y[:,2] - Cxy*inv(Cyy)*ϕ[2,:], digits = 5)

println("\nCxx = ")
display(Cxx)
println("\nCxy = ")
display(Cxy)
println("\nCx_y = ")
display(Cx_y)
println("\nμx_y = ")
display(μx_y)
println("\nCyy = ")
display(Cyy)
println("\nCyx = ")
display(Cyx)
println("\nCy_x = ")
display(Cx_y)
println("\nμy_x = ")
display(μy_x)

hsic = norm(Cyx - μy*transpose(μx))
hsicₕ = norm(Cyxₕ)

println("\nhsic = ",hsic)
println("hsicₕ = ",hsicₕ)


