using DrWatson
@quickactivate "dare"

using ReinforcementLearning
using PlotlyJS

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))


function reward(env)
    #implement your reward function here

    P_required = 466 # W

    u_l1_index = findfirst(x -> x == "u_l1", env.state_ids)

    u_l1 = env.state[u_l1_index]

    P_load = (env.norm_array[u_l1_index] * u_l1)^2 / 14
    
    P_diff = -abs(P_required - P_load) 
    reward = exp(P_diff/130) - 1

    return reward
end

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

env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, state_ids = states, rewardfunction = reward, x0=x0, v_dc=V_source, ts=ts, convert_state_to_cpu=true)
agent = create_agent_ddpg(na = na, ns = ns, use_gpu = agent_cuda)

hook = DataHook(state_ids = ["u_f1", "u_f2", "u_l1"])

run(agent, env, StopAfterEpisode(10), hook)





# PLOT

for episode in 10:10

    df_to_plot = a = hook.df[(hook.df.episode .== episode), :]

    layout = Layout(
            plot_bgcolor="#f1f3f7",
            title = "Example Plot<br><sub>Episode " * string(episode) * "</sub>",
            xaxis_title = "Time in Seconds",
            yaxis_title = "State values",
            autosize = false,
            width = 1000,
            height = 650,
            margin=attr(l=100, r=10, b=80, t=100, pad=10)
        )

    trace1 = scatter(df_to_plot, x = :time, y = :u_f1, mode="lines", name = "U f1")
    trace2 = scatter(df_to_plot, x = :time, y = :u_f2, mode="lines", name = "U f2")
    trace3 = scatter(df_to_plot, x = :time, y = :u_l1, mode="lines", name = "U l1")

    p = plot([trace1, trace2, trace3], layout)
    display(p)

end