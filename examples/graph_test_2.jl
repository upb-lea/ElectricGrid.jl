using GraphRecipes
using Plots

n = 15
A = Float64[ rand() < 0.5 ? 0 : rand() for i=1:n, j=1:n]
for i=1:n
    A[i, 1:i-1] = A[1:i-1, i]
    A[i, i] = 0
end

CM_net = graphplot(A,
        markersize = 0.2,
        node_weights = 1:n,
        markercolor = range(colorant"blue", stop=colorant"blue", length=n),
        names = 1:n,
        fontsize = 10,
        linecolor = :red
        )

display(CM_net)