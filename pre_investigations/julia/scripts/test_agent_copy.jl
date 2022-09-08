using DrWatson
@quickactivate "MicroGridSimWithRL"

using Flux
using CUDA
using ReinforcementLearning
using JSON
using BenchmarkTools

# include src
include(srcdir("agent_copy.jl"))
include(srcdir("nodeconstructor.jl"))  
include(srcdir("env.jl"))


CM = [ 0.  1.
            -1.  0.]

num_nodes = 50

CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(num_nodes) * ".json"))

CM = reduce(hcat, CM_list[1])'
CM = convert(Matrix{Int}, CM)

parameters = Dict()
nc = NodeConstructor(num_sources=num_nodes, num_loads=num_nodes, CM=CM)

#draw_graph(Grid_FC)   ---   not yet implemented

A, B, C, D = get_sys(nc)

limits = Dict("i_lim" => 20, "v_lim" => 600)

norm_array = vcat([limits[i] for j = 1:nc.num_sources for i in ["i_lim", "v_lim"]], [limits["i_lim"] for i = 1:nc.num_connections] )
norm_array = vcat( norm_array, [limits["v_lim"] for i = 1:nc.num_loads] )

states = get_states(nc)
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
Ad = exp(A*ts)
Bd = A \ (Ad - C) * B

env = SimEnv(A=A, B=B, C=C, Ad=Ad, Bd=Bd, norm_array=norm_array, x0=x0, v_dc=V_source, ts=rationalize(ts), convert_state_to_cpu=true)
gpu(env)
agent = create_agent(na, ns, 64)
gpu(agent.policy)

No_Episode = 100

run(agent,
            env,
            StopAfterEpisode(20)
            )