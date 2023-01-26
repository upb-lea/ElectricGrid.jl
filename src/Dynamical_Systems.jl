#_______________________________________________________________________________
# Dynamic systems

function Logistic_Map(x, λ; r = 3.9277370017867516)

    # Misiurewicz point:
    # r = 3.9277370017867516
    # Accumulation board - onset of chaos
    #r = 3.5699456718695445

    x = r*x*(1 - x)

    if λ == 0 
        λ = log.(2, abs.(r .- 2*r*x[1]))
    else
        λ = λ .+ log.(2, abs.(r .- 2*r*x[1]))
    end

    return x, λ
end

function Lorenz!(du, u, p, t)

    #= parameters
        p[1] = σ - Prandtl
        p[2] = ρ - Raleigh
        p[3] = β - geometric aspect ratio

        u0 = [1.0, 0.0, 0.0]
        tspan = (0.0, 100.0)
        p = [10.0, 28.0, 8/3]
    =#

    x, y, z = u
    σ, ρ, β = p

    du[1] = σ*(y - x)
    du[2] = ρ*x- y - x*z 
    du[3] = x*y - β*z

    return du
end

function Driven_Duffing!(du, u, p, t)

    #= Theory
        The driven duffing system, a forced harmonic ocillator, is an example 
        of a nonautomonous system. That is, it has a time dependence. It is also 
        an example of a three-dimensional system. Similarly, an nth-order time-
        dependent equation is a special case of an (n+1)-dimensional system. 
        By this trick, we may remove any time dependence by adding an extra 
        dimension to the system. A more physical motivation is that we need 
        to know three numbers u[1], u[2], and u[3], to predict that future,
        given the present. It is a 2nd order differential equaion

        m = mass
        δ = controls the amount of damping
        α = controls the linear stiffness
        β = controls the amount of non-linearity in the restoring force; 
            if zero the Duffing equation is a damped and driven simple harmonic oscillator
        γ = amplitude of the periodic driving force
        ω = angular frequency

        m*ddx + δ*dx + α*x + β*x^3 = γ*cos(ω*t)

        p = [m, δ, α, β, γ, ω]

        p = [1, 0.25, -1, 1, 0.4, 1]
        p = [1, 0.3, -1, 1, 0.5, 1.2] # Chaotic
    =#

    x, y, z = u
    m, δ, α, β, γ, ω = p

    du[1] = y
    du[2] = (1/m)*(-δ*y - α*x - β*x^3 + γ*cos(z))
    du[3] = ω

    return du
end

function Pendulum!(du, u, p, t)

    θ, ω = u
    m, g, L, b = p

    du[1] = ω
    du[2] = -b/m*ω - (g/L)*sin(θ)

    return du
end