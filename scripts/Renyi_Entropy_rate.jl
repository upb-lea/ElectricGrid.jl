using LinearAlgebra
using SpecialFunctions

T = Deus.Tα[end]#[0.2 0.5 0.3; 0.1 0.7 0.2; 0.1 0.5 0.4] #

a = 0.5

if a == 0
    R = Float64.(T .!= 0) # the connectivity / adjacency matrix of T
else
    R = T.^a
end

Rlr = eigvals(R)
Rlr0 = maximum(real.(Rlr))
RVl = eigvecs(transpose(R))
_, ilr0 = findmax(real.(Rlr))

Stat_dist_R = abs.(RVl[1:end, ilr0])/abs(sum(RVl[1:end, ilr0]))

h = (1/(1-a))*log(2, Rlr0)

lr = eigvals(T)
Vl = eigvecs(transpose(T))

index = findfirst( ==(1.0), real(round.(lr, digits = 1)))

Stat_dist = abs.(Vl[1:end, index])/abs(sum(Vl[1:end, index]))

hg = 0
hμ = 0
for v in 1:size(T,2)
    global hg, hμ

    for vd in 1:size(T,2)
  
        if T[v,vd] != 0

            local p, pe
            pe = Stat_dist[v]*T[v,vd]
            p = pe*log(2, T[v,vd])
            hμ = hμ - p

            if v == 4 && vd == 2 && false
                hg = hg + Stat_dist_R[v]*Deus.CM[v][3, 1]^a
                hg = hg + Stat_dist_R[v]*Deus.CM[v][3, 2]^a
            else
                hg = hg + Stat_dist_R[v]*T[v, vd]^a
            end
        end
    end

end

hg = ((1 - a)^-1)*log(2, hg)

println("α = ", a, "\n")

println("hμ = ", hμ,"\n") # Shannon Entropy rate

println("hg = ", hg) # Closed form weighted expression of Renyi entropy rate
println("h = ", h) # Renyi entropy rate in terms of the dominant eigenvalue of the perturbed matrix - spectral radius
