using PolygonArea

using LazySets

using VoronoiCells
using GeometryBasics
using Polyhedra

using Plots

points = [Point(10*rand(), 10*rand()) for _ in 1:7]
boundaries = Rectangle(Point2(0, 0), Point2(10, 10))
tess = voronoicells(points, boundaries);

p = plot(tess, legend = :topleft)
p = scatter!(tess.Generators, label = "Generators")
display(p)

function Polygon_Union()

    return polygon
end

#= temp = Vector{Point2{Float64}}(undef, 1)
temp[1] = tess.Cells[1][1]
polygon1 = [tess.Cells[1]; temp]
temp[1] = tess.Cells[2][1]
polygon2 = [tess.Cells[2]; temp]
points = tess.Generators

p = plot(polygon1, label = "polygon1", colorbar=false)
p = plot!(polygon2, label = "polygon2")

t = reduce(vcat, transpose.(tess.Cells[1]))
temp = [t[i,:] for i in axes(t,1)]
Lz_poly1 = VPolygon(temp)
t = reduce(vcat, transpose.(tess.Cells[2]))
temp = [t[i,:] for i in axes(t,1)]
Lz_poly2 = VPolygon(temp)
t = reduce(vcat, transpose.(tess.Cells[3]))
temp = [t[i,:] for i in axes(t,1)]
Lz_poly3 = VPolygon(temp)

t = reduce(vcat, transpose.(points))
inside = Vector{Bool}(undef, size(t,1))

partition = UnionSet(Lz_poly1, Lz_poly2)

for i in axes(t,1)
    inside[i] = ∈(element(Singleton(t[i,:])), partition)
end

#p = plot(Lz_poly1, label = "polygon1", colorbar=false)
#p = plot!(Lz_poly2, label = "polygon2")
p = plot(partition, label = "partition", colorbar=false)
p = scatter!(points, marker_z = inside, label = "inside")
display(p)


#reduce(vcat, transpose.(points))

a = rectangle((0.25, 0.25), (2.0, 2.0))
b = rectangle((0.0, 0.0), (1.0, 1.0))
c = rectangle((0.5, 0.5), (1.5, 1.5))

x = b ∪ c
y = a ∩ x

plot(a, fill =:green)
plot!(b, fill =:blue)
plot!(c, fill =:red)

plot(y, fill =:green) =#