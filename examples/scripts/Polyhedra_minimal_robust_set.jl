# https://juliapolyhedra.github.io/Polyhedra.jl/stable/generated/Minimal%20Robust%20Positively%20Invariant%20Set/


A = [1 1; 0 1] - [1; 1] * [1.17 1.03]

using Polyhedra
Wv = vrep([[x, y] for x in [-1.0, 1.0] for y in [-1.0, 1.0]])

using GLPK
using JuMP
lib = DefaultLibrary{Float64}(GLPK.Optimizer)
W = polyhedron(Wv, lib)


function Fs(s::Integer, verbose=1)
    @assert s ≥ 1
    F = W
    A_W = W
    for i in 1:(s-1)
        A_W = A * A_W
        F += A_W
        if verbose ≥ 1
            println("Number of points after adding A^$i * W: ", npoints(F))
        end
        removevredundancy!(F)
        if verbose ≥ 1
            println("Number of points after removing redundant ones: ", npoints(F))
        end
    end
    return F
end

F_1 = Fs(1)
F_2 = Fs(2)
F_3 = Fs(3)

using Plots
p = plot()
for i in 10:-1:1
    plot!(Fs(i, 0))
end
plot!()
display(p)
#=
using PlotlyJS
layout = Layout(
            plot_bgcolor="#f1f3f7",
            title = "DoubleIntegrator",
            xaxis_title = "x1",
            yaxis_title = "x2",
            width = 1200,
            height = 850,
            margin=attr(l=100, r=200, b=80, t=100, pad=10)
        )

traces = []

for i in 10:-1:1
    push!(traces, Fs(i, 0))
end

p = PlotlyBase.Plot(traces, layout)#, config = PlotConfig(scrollZoom=true))
display(p)
=#
