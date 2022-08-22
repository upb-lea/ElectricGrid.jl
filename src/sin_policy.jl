
"""
Simple basic example policy which has internal time which is increased every call by ts
and returns as action a three-phase sinewave (3 sinwaves shifted by 120°) with a currently fixed amplitide
"""
Base.@kwdef mutable struct sin_policy <: AbstractPolicy

    n_actions = 1
    action_space::Space{Vector{ClosedInterval{Float64}}} = Space([ -1.0..1.0 for i = 1:n_actions], )
    t = 0.0
    ts = 1e-4


end

function (p::sin_policy)(env)
    p.t = p.t + p.ts
    u = [sqrt(2)*230 * sin.((50*2*π)*p.t .+ 2/3*π*(i-1)) for i = 1:length(p.action_space)]
    return u    
end