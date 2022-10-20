mutable struct dare_setup
    env_cuda
    agent_cuda
    num_sources
    num_loads
    CM
    parameters
    ts
    t_end
    env
    controller
    hook
end

function create_setup(;num_sources, num_loads, CM=nothing, parameters=nothing, ts=0.0001, t_end=0.05)
    
end