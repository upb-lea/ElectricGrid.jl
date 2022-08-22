using LinearAlgebra
using SpecialFunctions

T = [0.2 0.5 0.3; 0.1 0.7 0.2; 0.1 0.5 0.4]

a = 0.0

Q = Float64.(T .!= 0) # the connectivity / adjacency matrix of T

if a == 0
    R = Q
else
    R = T.^a
end

Rlr = eigvals(R)
Rlr0 = maximum(real.(Rlr))

h = (1/(1-a))*log(2, Rlr0)

lr = eigvals(T)
Vl = eigvecs(transpose(T))

index = findfirst( ==(1.0), real(round.(lr, digits = 1)))

Stat_dist = abs.(Vl[1:end, index])/abs(sum(Vl[1:end, index]))

he = 0
ha = 0
hu = 0
hg = 0
for v in 1:size(T,2)
    global he, ha, hu, hg

    he = 0
    for vd in 1:size(T,2)
  
        if T[v,vd] != 0
            he = he + (T[v,vd]^a)
            hg = hg + ((Stat_dist[v])/(sum(Stat_dist.^a)))^a*T[v, vd]

            hu = hu - Stat_dist[v]*T[v,vd]*log(2, T[v,vd])
        end
    end

    ha = ha + ((1 - a)^-1)*Stat_dist[v]*log(2, he)
end

hg = ((1 - a)^-1)*log(2, hg)

println("hu = ", hu)
println("\nha = ", ha)
#println("hg = ", hg)
println("h = ", h)
