
"""
Simple basic example policy which has internal time which is increased every call by ts
and returns as action a three-phase sinewave (3 sinwaves shifted by 120°) with a currently fixed amplitide
"""
Base.@kwdef mutable struct Classical_Policy <: AbstractPolicy

    n_actions = 1
    action_space::Space{Vector{ClosedInterval{Float64}}} = Space([ -1.0..1.0 for i = 1:n_actions], )
    t = 0.0
    ts = 1e-4
    Source_Indices = [1]

end

function (p::Classical_Policy)(env::ElectricGridEnv, name::Union{String, Nothing} = nothing)
    p.t = p.t + p.ts
    #u = [230 * sin.(50*2*pi*p.t .+ 2/3*pi*(i-1)) for i = 1:length(p.action_space)]
    u = [.5 for i = 1:length(p.action_space)]
    #u = [0.5, 0, -0.5, 0.7, 0.2, -0.2]
    # u = [1 * sin.(50*2*pi*p.t .- 2/3*pi*(i-1)) for i = 1:3]
    return u
    #return [u[1], u[1], u[2], u[2], u[3], u[3]]   # order for 2 sources depening on GetActionIds(env.nc) 
end