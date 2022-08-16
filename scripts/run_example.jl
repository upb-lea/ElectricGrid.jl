using ReinforcementLearning

include("../src/nodeconstructor.jl")
include("../src/env.jl")
include("../src/agent_ddpg.jl")
include("../src/data_hook.jl")


env_cuda = false
agent_cuda = false

CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

parameters = Dict()

nc = NodeConstructor(num_sources = 2, num_loads = 1, CM = CM)

A, B, C, D = get_sys(nc)

limits = Dict("i_lim" => 20, "v_lim" => 600)

states = get_state_ids(nc)
norm_array = []
for state_name in states
    if startswith(state_name, "i")
        push!(norm_array, limits["i_lim"])
    elseif startswith(state_name, "u")
        push!(norm_array, limits["v_lim"])
    end
end

ns = length(A[1,:])
na = length(B[1,:])

# time step
ts = 1e-5

V_source = 300

x0 = [ 0.0 for i = 1:length(A[1,:]) ]

if env_cuda
    A = CuArray(A)
    B = CuArray(B)
    C = CuArray(C)
    x0 = CuArray(x0)
end

env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, state_ids = states, x0=x0, v_dc=V_source, ts=rationalize(ts), convert_state_to_cpu=true)
agent = create_agent_ddpg(na = na, ns = ns, use_gpu = agent_cuda)

hook = DataHook(state_ids = ["u_f1", "u_f2", "u_l1"])

run(agent, env, StopAfterEpisode(10), hook)