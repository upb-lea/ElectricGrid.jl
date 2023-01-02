using StatsBase
using Plots

#https://github.com/cantaro86/Financial-Models-Numerical-Methods/blob/master/6.1%20Ornstein-Uhlenbeck%20process%20and%20applications.ipynb

N = 100 # time steps
sources = 1 # number of paths/sources
T = 10 # final time

κ = 1.000 # mean reversion parameter
σ = 0.00 # Brownian motion scale (standard deviation) - sqrt(diffusion) 
γ = 10.5 #0.5 # asymptotoic mean

#Note: κ and σ can't be negative

X₀ = -γ

std_asy = sqrt(σ^2/(2*κ)) # asymptotic standard deviation

dt = T/N
t = 0:dt:T

X = Array{Float64, 2}(undef, N+1, sources)
Y = Array{Float64, 2}(undef, N+1, sources)
Z = Array{Float64, 2}(undef, N+1, sources)
V = Array{Float64, 2}(undef, N+1, sources)

X[1, :] = fill!(X[1, :], X₀)
Y[1, :] = fill!(Y[1, :], X₀)
Z[1, :] = fill!(Z[1, :], 0)
V[1, :] = fill!(V[1, :], X₀)

ingr = fill!(ingr, 0.)

dW = randn(N + 1, sources)

# Euler Maruyama
for i in 1:N

    global κ, γ, σ

    X[i + 1, :] = X[i, :] + κ*(γ .- X[i, :])*dt + σ*sqrt(dt)*dW[i, :]
end

# simplified

std_dt = sqrt(σ^2/(2*κ) * (1 - exp(-2*κ*dt)))

for i in 1:N

    global κ, γ, σ

    Y[i + 1, :] = γ .+ exp(-κ*dt)*(Y[i, :] .- γ) + std_dt * dW[i, :]
end

for i in 1:N

    Z[i + 1, :] = Z[i, :] .+ σ*sqrt(dt)*dW[i, :]
end

p = plot(t, X[:, 1],
title = "Ornstein-Uhlenbeck",
label = "Euler Maruyama",
xlabel = "Time [s]")

p = plot!(t, Y[:, 1],
label = "simple")

#= p = plot!(t, Z[:, 1],
label = "Wiener") =#
#= 
p = plot!(t, W[:, 1],
label = "wiki") =#

for i in 2:sources

    global p
    p = plot!(t, X[:, i],
    label = false)
    p = plot!(t, Y[:, i],
    label = false)
    #= p = plot!(t, Z[:, i],
    label = false) =#
end

display(p)

#= T = 100
α = 1 # mean reversion parameter
β = 1 # Brownian motion scale (standard deviation)
γ = 1 # asymptotoic mean

X₀ = γ

t = 1:T

dW = randn(T)
dW[1] = 0

W_cs = Array{Float64, 1}(undef, T)
W_cs[1] = 0

OU = Array{Float64, 1}(undef, T)
OU[1] = X₀

dW_cs = 0.
ingral = 0.

for t in 1:T

    global dW_cs, X₀, γ

    dW_cs += dW[t]
    W_cs[t] = dW_cs

    expₐ = exp(-α*t)

end

p5 = plot(t, W_cs,
title = "Wiener",
label = false,
xlabel = "Time [s]",
ylabel = "Current [A]")

display(p5) =#