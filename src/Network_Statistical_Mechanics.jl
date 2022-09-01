#using DrWatson
#@quickactivate "dare"

using Plots

include("Complexity.jl")

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#_______________________________________________________________________________
# Parameters - Time simulation

Timestep = 10 # time step in μs
t_final = 0.1 #0.75 # time in seconds, total simulation run time
fsys = 2000 # Hz, fundamental frequency of system
fsys = 1/(10e-6) 

#_______________________________________________________________________________
# Environment Calcs
N_s = (1/(Timestep*1e-6)) # time intervals
μ_s = 1/N_s # time step
N = convert(Int64, floor(t_final/μ_s)) + 1

t = 0:μ_s:t_final # time

#_______________________________________________________________________________
# Machine Parameters

#= Every other one
    D = 5
    Timestep = 10
    t_final = 2
    δ = 0.05
=#
#= Logistic Map
    r = 3.9277370017867516
    D = 14
    Timestep = 10
    t_final = 1.3
    δ = 0.05
    ϵ[:] = 0.5
    λ = 0.7898

    r = 3.7
    D = 16
    Timestep = 10
    t_final = 2
    δ = 0.05
    ϵ[:] = 0.5
    λ = 0.5118
=#
#= CCA
    T = 1    T = 2   T = 3      T = 4
    Cμ = 2   Cμ = 3  Cμ = 8.5   Cμ = 0 

    dim = 1
    1:200:3801
    L = 5000
    D = 150
    Timestep = 10
    t_final = 0.005
    δ = 0.05
    ϵ[:] = 0.25
    0.0 <= x_range <= 4.0
=#
#= periodic
    fsys = 2000
    t_final = 0.05
    D = 4
    δ = 0.05
    dim = 1
    ϵ[:] = [0.5
    x_range[:,1] = [1.0
    x_range[:,2] = [0.0
    L = 1
=#
#= Weiss's Even Process
    D = 9
    Timestep = 10
    t_final = 2.0
    δ = 0.05
    ϵ[:] = 0.5
=#

dim = 1 # dimensions of state space
#hl = 0.125 # lower frequency - harmonic order
#hf = 1 # upper frequency - harmonic order
#μ_m = 0.5/(hf*fsys) #sampling timestep
μ_m = μ_s #sampling timestep

#= 
    ϵ = 0.5
    P = 2        P = 4        P = 8          P = 16
    hl = 1       hl = 0.5     hl = 0.25      hl = 0.125       
    D0           D0 + 3       D0 + 7         D0 + 15
=#

#D0 = 4*hf

#_______________________________________________________________________________
## Generating Connectivity Matrix

L = 1 # number of agents

Z = 2 # coordination number, if Z = L, then we have a fully connected graph. If Z = 2, p = 0 then circular
p = 0.4 # the fraction of random connections

#Net = SmallWorld(L, Z, p) # create small world graph
#Net = Barabasi_Albert(L, m = 1) # create scale-free network
Net, L = SquareLattice(L, torus = true) # create Neighbour connected graph

#_______________________________________________________________________________
# Dynamical System

x = Array{Float64, 2}(undef, L, N) # State space
x[:, 1] = initialise_CCA(x[:, 1])

#= The Logistic Map

    1. Discrete quadratic nonlineariy
    2. 0 <= r <= 4
    3. 0 <= s[1] <= 1

    x[i+1] = r*x[i]*(1 - x[i])

    The control parameter r governs the degree of nonlinearity. At a Misiurewicz
    parameter value the chaotic behaviour is governed by an absolutely continuous
    invariant measure. The consequence is that the statistical properties are
    particularly well-behaved. These parameter values are determined by the condition
    that the iterates f^(n)(xc) of the map's maximum xc = 1/2 are asymptotically
    stable.

    We can choose a partitioning {[0, xc],(xc,1)}

    #er = 0.01
    #r = er + 3.0 # period 2
    #r = er + 3.449490 # period 4
    #r = er + 3.544090 # period 8
    #r = er + 3.564407 # period 16?
    #r = er + 3.568750 # period 32?
    #r = er + 3.56969 # period 64?
    #r = er + 3.56989 # period 128?
=#

# Misiurewicz point:
r = 3.9277370017867516
# Accumulation board - onset of chaos
#r = 3.5699456718695445
#r = 3.567
#r = 4.0

#x[:,1] = 1/2 # initial conditions

x[:, 1] .= 0.4
logdFdx = log.(2, abs.(r .- 2*r*x[:,1]))

#_______________________________________________________________________________
#%% Starting time simulation
println("\nHere we go.\n")
even = 0
@time begin

    println("Progress : 0.0 %")

    for i in 1:N-1

        global even, logdFdx

        # Progress Bar
        if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
            flush(stdout)
            println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
        end

        #= CCA #
        x[:, i + 1] = CCA_dynamics(Net, x[:, i]; T = 3)
        =#

        # Logistic Map
        x[1, i + 1] = r*x[1, i]*(1 - x[1, i])
        logdFdx = logdFdx .+ log.(2, abs.(r .- 2*r*x[:,i]))
        #

        #= periodic
        x[1, i + 1] = 1*sin(2π*hf*fsys*t[i+1]) + 1*sin(2π*hl*fsys*t[i+1])
        =#

        #= Every other One
        if i == 1
            Deus.k = N
        end
        if i%2 == 0
            Deus.s[i] = 1
        else
            b = rand()
            if b > 0.5
                Deus.s[i] = 0
            else
                Deus.s[i] = 1
            end
        end
        =#

        #= Weiss's Even Process
        if i == 1

            Deus.s[i] = 0
            even = 0
            Deus.k = N
        elseif i == 2

            if Deus.s[i-1] == 0

                b = rand()
                if b > 0.5
                    Deus.s[i] = 0
                    even = 0
                else
                    Deus.s[i] = 1
                    even = even + 1
                end
            else

                Deus.s[i] = 1
                even = even + 1
            end

        elseif Deus.s[i-1] == 0

            b = rand()
            if b > 0.5
                Deus.s[i] = 0
                even = 0
            else
                Deus.s[i] = 1
                even = even + 1
            end
        elseif Deus.s[i-1] == 1 && even%2 == 1

            Deus.s[i] = 1
            even = even + 1
        elseif Deus.s[i-1] == 1 && even%2 == 0

            b = rand()
            if b > 0.5
                Deus.s[i] = 0
                even = 0
            else
                Deus.s[i] = 1
                even = even + 1
            end
        end
        =#
        
    end

    println("Progress : 100.0 %\n")

    D = 32 #convert(Int, ceil(D0 + 2/hl - 1)) + 120

    ϵ = Array{Float64, 1}(undef, dim) # Instrument resolution
    ep = sum(x[1:1, :])/N
    ϵ[:] = [0.5 for i in 1:dim]

    x_range = Array{Float64, 2}(undef, dim, 2) # State space limits
    x_range[:, 1] = [1.0 for i in 1:dim] # maximum
    x_range[:, 2] = [0.0 for i in 1:dim] # minimum

    Deus = ϵ_Machine(N, D, ϵ, x_range, μ_m, μ_s, δ = 0.05)

    Cranking(Deus, x[1:1, D:end], μ_s)
end

T_plot_end = 30

N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/N_s)*N_s))
N_range = 1:N_plot_end + 2
Nm_plot_end = convert(Int64, round((T_plot_end/fsys  - μ_m)*1/μ_m))
Nm_range = 1:Nm_plot_end + 2

p1 = plot(t[N_range], x[1, N_range], label = "x")
#p1 = plot!(Deus.t_m[Nm_range], Deus.s[Nm_range])
p1 = plot!(Deus.t_m[Nm_range], Deus.x_m[1, Nm_range], label = "xm")
#plot!(N_range, x[1, N_range])
#plot!(Nm_range, Deus.s[Nm_range])
#plot!(Nm_range, Deus.x_m[1, Nm_range])

display(p1)
#

#_______________________________________________________________________________
## Draw Graph
#Draw_Plots(Deus)
#Draw_Graph(Net, x, N, run = 0)
#All_Plots(Deus)

#=
    Rlr = eigvals(Deus.Tα[1])
    Rlr0 = maximum(real.(Rlr))
    RVl = eigvecs(transpose(Deus.Tα[1]))
    _, ilr0 = findmax(real.(Rlr))

    Stat_dist_R = abs.(RVl[1:end, ilr0])/abs(sum(RVl[1:end, ilr0]))
    a = 0
    h = (1/(1-a))*log(2, Rlr0)

    println("hl = ", 2/2^(Deus.Hα[end]))
    println("Zα[1]/L = ", log(2, Deus.Zα[1])/Deus.Tree.L)
    println("logdFdx/N = ", logdFdx/N)
    println("logdFdx/N = ", Deus.hα[end])
    rintln("Zα[1]/L = ", log(Deus.Zα[1] - Deus.Tree.Nodes[1].morph_branches)/Deus.Tree.D)
    g = 0
    for i in 1:length(Deus.Tree.Nodes)
        global g
        if Deus.Tree.Nodes[i].generation == Deus.Tree.L
            g = g + Deus.Tree.Nodes[i].morph_branches
        end

    end
    println("Zα[1]/L = ", log(2, g)/Deus.Tree.D)
=#

println("Cμ_t[end] = ", Deus.Cμ_t[end])
println("Hα[end] = ", Deus.Hα[end])
println("Period = ", 2^(Deus.Hα[end]))

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")

