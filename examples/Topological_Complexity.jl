using DrWatson
@quickactivate "dare"

include(srcdir("Complexity.jl"))

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")
#=
function Draw_Graph(Deus::ϵ_Machine, Net, nodekey, N)

    cm = SimpleGraph(Net) # create graph object from Connecitivity Matrix

    x = LinRange(0, 10, 100)
    y = sin.(x)

    fig = Figure(resolution=(1600, 800))

    fig[1, 1:2] = title = Label(fig, "Topological Complexity", textsize = 30)
    title.tellwidth = false

    Cμ_t_ax = Axis(fig[2, 2], title = "Time Dependent Complexity - Knowledge Relaxation", xlabel = "Time [seconds]", ylabel = "Cμ [bits]")
    lines!(Cμ_t_ax, Deus.t_m, Deus.Cμ_t, color = :red)

    fig[2, 1] = graph_ax = Axis(fig)
    
    nodecolours = Array{Symbol, 1}(undef, L)
    nodecolours = fill!(nodecolours, :blue)
       
    #set_theme!(resolution = (800, 800))
    p = graphplot!(graph_ax, cm;
            edge_color = [RGBAf(0,0,0,0) for i in 1:ne(cm)], #RGBAf(0,0,0,0)
            node_size = 20,
            node_marker = :rect,
            fontsize = 1,
            node_color = nodecolours,
            #curve_distance =-.5,
            #curve_distance_usage = true
            )
    
    #Circular_Layout(p, graph_ax, L)
    Plane_Layout(p, graph_ax, L)
    #p.layout = Spring(Ptype = Float32)
    hidedecorations!(graph_ax)
    #hidespines!(ax)
    #ax.aspect = DataAspect()
    display(fig)

    for t in 1:N

        for n in 1:L

            if nodekey[n, t] == 0
                nodecolours[n] = :green4   
            elseif nodekey[n, t] == 1
                nodecolours[n] = :dodgerblue
            elseif nodekey[n, t] == 2
                nodecolours[n] = :royalblue4
            elseif nodekey[n, t] == 3
                nodecolours[n] = :navyblue
            else
                nodecolours[n, t] = :blue
            end
        end

        p.node_color[][:] = nodecolours[:]
        p.node_color[] = p.node_color[]

        sleep(0.05)

    end
    
    return nothing
end
=#
#_______________________________________________________________________________
# Parameters - Time simulation
fsys = 2000 # Hz
Timestep = 10 # time step in μs
t_final = 0.005 # time in seconds, total simulation run time

#_______________________________________________________________________________
# Environment Calcs
N_s = (1/(Timestep*1e-6)) # time intervals
μ_s = 1/N_s # time step
N = convert(Int64, floor(t_final/μ_s)) + 1

t = 0:μ_s:t_final # time

#_______________________________________________________________________________
# Machine Parameters
D = 160 # 20 for map, 10 or 20 for every other 1, ϵ = 0.5 # 15 for map gives trouble
δ = 0.05 # 0.025 for map, 0.002 for every other, ϵ = 0.5 # 0.025 for map gives trouble

#= Every other one
    D = 4
    Timestep = 10
    t_final = 10
    δ = 0.005
=#
#= Logistic Map
    D = 12
    Timestep = 10
    t_final = 10
    δ = 0.035
=#
#= CCA
    T = 1    T = 2   T = 3      T = 4
    Cμ = 2   Cμ = 3  Cμ = 8.5   Cμ = 0 

    dim = 20
    1:200:3801
    L = 5000
    D = 10
    Timestep = 10
    t_final = 0.005
    δ = 0.05
=#

dim = 20 # dimensions of state space
#μ_m = 0.25/150 #sampling timestep
μ_m = μ_s #sampling timestep
t_m = 0:μ_m:t_final

x_range = Array{Float64, 2}(undef, dim, 2) # State space limits
ϵ = Array{Float64, 1}(undef, dim) # Instrument resolution

ϵ[:] = [0.25 for i in 1:dim]

x_range[:,1] = [4.0 for i in 1:dim] # maximum

x_range[:,2] = [0.0 for i in 1:dim] # minimum
#_range[1,2] = -1.1 # minimum

Deus = ϵ_Machine(N, D, δ, ϵ, x_range, μ_m, μ_s)

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
    =#

# Misiurewicz point:
#r = 3.9277370017867516

#xc = 1/2

#x[1] = xc # initial conditions
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

        #= CCA =#
        Sampling(Deus, x[1:200:3801, i], i, μ_s)   
        x[:, i + 1] = CCA_dynamics(Net, x[:, i]; T = 3)
        #= =#
        #= Map =
            Sampling(Deus, x[:, i], i, μ_s) 
            x[1, i + 1] = r*x[1, i]*(1 - x[1, i])
        = =#

        #= periodic
            x[1, i + 1] = sin(2π*fsys*t[i+1]) + 0.5*sin(2π*3*fsys*t[i+1])
        =#

        #= Every other One
            if i%2 == 0
                Deus.s[i] = 1
            else
                b = randn()
                if b > 0
                    Deus.s[i] = 0
                else
                    Deus.s[i] = 1
                end
            end
        =#
        
    end

    println("Progress : 100.0 %\n")
end

Cranking(Deus)

#=

    T_plot_end = 5

    N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/N_s)*N_s))
    N_range = 1:N_plot_end+2
    Nm_plot_end = convert(Int64, round((T_plot_end/fsys  - μ_m)*1/μ_m))
    Nm_range = 1:Nm_plot_end+2

    #p = plot(t[N_range], x[1, N_range])
    #p = plot!(t_m[Nm_range], Deus.s[Nm_range])
    #p = plot!(Deus.t_m[Nm_range], Deus.x_m[1, Nm_range])

    display(p)
=#

#_______________________________________________________________________________
## Draw Graph

#Draw_Graph(Deus, Net, x, N)

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
println(Deus.Cμ_t[end])
