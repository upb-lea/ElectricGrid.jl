
using Polyhedra
using Plots
import QHull
import GLPK

lib1 = QHull.Library()
lib2 = DefaultLibrary{Float64}(GLPK.Optimizer)

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

x = 10*rand(1, 100000, 1);

v = Vector{MixedMatVRep{Float64, Matrix{Float64}}}(undef, size(x,1))
p1 = Vector{QHull.Polyhedron}(undef, size(x,1))
p2 = Vector{DefaultPolyhedron{Float64, MixedMatHRep{Float64, Matrix{Float64}}, MixedMatVRep{Float64, Matrix{Float64}}}}(undef, size(x,1))

err = Vector{Float64}(undef, size(x,1))

for i in axes(x, 1)
    v[i] = vrep(x[i, :, :])    
end

t1 = @elapsed begin

    for i in axes(x, 1)
        p1[i] = polyhedron(v[i], lib1)
        removevredundancy!(p1[i])
    end

end

t2 = @elapsed begin

    for i in axes(x, 1)
        p2[i] = polyhedron(v[i], lib2)
        removevredundancy!(p2[i])
    end

end

c = 1
for i in axes(x, 1)
    if size(p1[i].vrep.V, 1) == size(p2[i].vrep.V, 1)
        err[i] = maximum(p1[i].vrep.V .- p2[i].vrep.V)
    end
end
println(maximum(err))

println("Qhull =", t1)
println("GLPK =", t2)

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")