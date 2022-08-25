#using DrWatson
#@quickactivate "dare"

include("Complexity.jl")

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#_______________________________________________________________________________
# Parameters - Time simulation

Timestep = 10 # time step in μs
t_final = 0.005 #0.75 # time in seconds, total simulation run time
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
    t_final = 10
    δ = 0.005
    ϵ[:] = 0.5
=#
#= Logistic Map
    D = 12 - or 16
    Timestep = 10
    t_final = 0.75 - or 3.75
    δ = 0.05
    ϵ[:] = 0.5
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
D = 150 #convert(Int, ceil(D0 + 2/hl - 1)) + 120

x_range = Array{Float64, 2}(undef, dim, 2) # State space limits
ϵ = Array{Float64, 1}(undef, dim) # Instrument resolution

ϵ[:] = [0.25 for i in 1:dim]

x_range[:,1] = [4.0 for i in 1:dim] # maximum
x_range[:,2] = [0.0 for i in 1:dim] # minimum

Deus = ϵ_Machine(N, D, ϵ, x_range, μ_m, μ_s)

#_______________________________________________________________________________
## Generating Connectivity Matrix

L = 5000 # number of agents

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

    This partition is generating and the resulting binary sequences completely
    capture the statistical properties of the Logistic map. In other words, there is
    a one-to-one mapping between infinite binary sequences and almost all points on
    the attractor.

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

#x[1] = 1/2 # initial conditions

#x[1] = 0

#_______________________________________________________________________________
#%% Starting time simulation
println("\nHere we go.\n")

@time begin

    println("Progress : 0.0 %")

    for i in 1:N-1

        # Progress Bar
        if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
            flush(stdout)
            println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
        end

        # CCA #
        #Sampling(Deus, x[1:200:3801, i], i, μ_s)   
        #Sampling(Deus, x[1:1, i], i, μ_s)  
        x[:, i + 1] = CCA_dynamics(Net, x[:, i]; T = 3)
        #

        #= Map
        x[1, i + 1] = r*x[1, i]*(1 - x[1, i])
        =#

        #= periodic
        x[1, i + 1] = 1*sin(2π*hf*fsys*t[i+1]) + 1*sin(2π*hl*fsys*t[i+1])
        =#

        #= Every other One
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
        
    end
    #Sampling(Deus, x[:, N], N, μ_s) 
    println("Progress : 100.0 %\n")
end

Cranking(Deus, x[1:1, :], μ_s)

T_plot_end = 30

N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/N_s)*N_s))
N_range = 1:N_plot_end + 2
Nm_plot_end = convert(Int64, round((T_plot_end/fsys  - μ_m)*1/μ_m))
Nm_range = 1:Nm_plot_end + 2

fig = Figure()
p = Axis(fig[1, 1])
lines!(p, t[N_range], x[1, N_range])
#lines!(p, Deus.t_m[Nm_range], Deus.s[Nm_range])
lines!(p, Deus.t_m[Nm_range], Deus.x_m[1, Nm_range])
#lines!(p, N_range, x[1, N_range])
#lines!(p, Nm_range, Deus.s[Nm_range])
#lines!(p, Nm_range, Deus.x_m[1, Nm_range])

#display(fig)
#

#_______________________________________________________________________________
## Draw Graph
Draw_Plots(Deus)
#Draw_Graph(Net, x, N, run = 0)

println("Cμ_t[end] = ", Deus.Cμ_t[end])
println("Hα[end] = ", Deus.Hα[end])
println("Period = ", 2^(Deus.Hα[end]))
println("hl = ", 2/2^(Deus.Hα[end]))

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")

