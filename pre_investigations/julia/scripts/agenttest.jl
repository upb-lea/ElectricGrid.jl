using DrWatson
@quickactivate "MicroGridSimWithRL"

using PProf, Profile
# using DifferentialEquations
# using Sundials
# using LinearAlgebra
# using ControlSystems
using BenchmarkTools
using ReinforcementLearning
using Flux
using StableRNGs
using IntervalSets
using TimerOutputs
using JSON
using PlotlyJS

import ReinforcementLearning.update!


# include(srcdir("collect_timing_results.jl"))
include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent.jl"))
include(srcdir("run_timed.jl"))

# collect_results = Dict()


# reset_timer!(timer::TimerOutput)

no_nodes = []
overall_run = []
inside_run = []
policy_update = []
grad_Actor = []
grad_Critic=[]
data_transfer=[]
update_Critic=[]
update_Actor=[]
env_calc = []
prepare_data = []

Timing_dict = Dict()

# global timer = TimerOutput()

function train_ddpg(timer::TimerOutput, 
                    env_cuda::Bool, 
                    agent_cuda::Bool, 
                    num_nodes::Int, 
                    No_Episodes::Int,
                    batch)
    
    env_cuda = env_cuda
    agent_cuda = agent_cuda

    num_nodes = num_nodes
    println(num_nodes)

    CM = [ 0.  1.
            -1.  0.]

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

    if env_cuda
        A = CuArray(A)
        B = CuArray(B)
        C = CuArray(C)
        Ad = CuArray(Ad)
        Bd = CuArray(Bd)
        x0 = CuArray(x0)
    end


    env = SimEnv(A=A, B=B, C=C, Ad=Ad, Bd=Bd, norm_array=norm_array, x0=x0, v_dc=V_source, ts=rationalize(ts), convert_state_to_cpu=true)
    agent = create_agent(na, ns, batch, agent_cuda)

    hook = TotalRewardPerEpisode()
    
    # time parts of the algorithm
    @timeit timer "Overall run" begin
        run(
            agent,
            env,
            timer,
            StopAfterEpisode(No_Episodes),
            hook
        )
        end
    
    nothing
end

P_required = 466 # W
V_required = 230 # V
PLoad = []
Pdiff = []

function reward_func(method::String, env::SimEnv)

    i_1, u_1, i_c1, u_l1 = Array(env.state)

    P_load = (env.norm_array[end] * u_l1)^2 / 14
    
    if method == "Power_exp"
        # push!(PLoad, P_load)
        P_diff = -abs(P_required - P_load) 
        # push!(Pdiff, P_diff)
        reward = exp(P_diff/130) - 1
    
    elseif method == "Power"
        reward = -abs(P_required - P_load) / (600 * 20)

    elseif method == "Voltage"
        reward = -((V_required - u_l1)/ 600) ^2
    end
    return reward
end

function collect_results!(timer::TimerOutput, node::Int)
    append!(no_nodes, node)
    append!(overall_run, TimerOutputs.time(timer["Overall run"]))
    append!(inside_run, TimerOutputs.time(timer["Overall run"]["inside run"]))
    append!(policy_update, TimerOutputs.time(timer["Overall run"]["inside run"]["Agent policy update"]))
    append!(env_calc, TimerOutputs.time(timer["Overall run"]["inside run"]["Env calculation"]))
    append!(grad_Actor, TimerOutputs.time(timer["Overall run"]["inside run"]["Agent policy update"]["gradients for Actor"]))
    append!(grad_Critic, TimerOutputs.time(timer["Overall run"]["inside run"]["Agent policy update"]["gradients for Critic"]))
    append!(data_transfer, TimerOutputs.time(timer["Overall run"]["inside run"]["Agent policy update"]["policy - transfer of data to device"]))
    append!(update_Critic, TimerOutputs.time(timer["Overall run"]["inside run"]["Agent policy update"]["update Critic network"]))
    append!(update_Actor, TimerOutputs.time(timer["Overall run"]["inside run"]["Agent policy update"]["update Actor network"]))
    append!(prepare_data, TimerOutputs.time(timer["Overall run"]["inside run"]["Agent prepare data"]))
end



# collect_results!(timer, 1)

env_cuda = true
agent_cuda = true

batch_sizes = [16, 32]
nodes = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
No_Episodes = 20


train_ddpg(reset_timer!(TimerOutput()), env_cuda, agent_cuda, 1, 10, 64)

for batch in batch_sizes
    global no_nodes = []
    global overall_run = []
    global inside_run = []
    global policy_update = []
    global grad_Actor = []
    global grad_Critic=[]
    global data_transfer=[]
    global update_Critic=[]
    global update_Actor=[]
    global env_calc = []
    global prepare_data = []

    for i = 1:length(nodes)
        to = TimerOutput()
        reset_timer!(to)
        train_ddpg(to, env_cuda, agent_cuda, nodes[i], No_Episodes, batch)
        collect_results!(to, nodes[i])
    end



    name = @ntuple env_cuda agent_cuda No_Episodes batch

    layout = Layout(
        plot_bgcolor="#f1f3f7",
        title = "Training time<br><sub>Env: " * (env_cuda ? "GPU" : "CPU") * ", Agent: " * (agent_cuda ? "GPU" : "CPU") * ", Batch Size: " * string(batch) * "</sub>",
        xaxis_title = "Network Size",
        yaxis_title = "Time in Seconds",
        autosize = false,
        width = 1000,
        height = 650,
        margin=attr(l=20, r=10, b=10, t=100, pad=10)
    )

    trace1 = scatter(x = nodes, y = overall_run ./1e9, mode="lines", name = "Overall Run")
    trace2 = scatter(x = nodes, y = env_calc ./1e9, mode="lines", name = "Env Calculation")
    trace3 = scatter(x = nodes, y = policy_update ./1e9, mode="lines", name = "Policy Update")

    p = plot([trace1, trace2, trace3], layout)
    display(p)
    savefig(p, plotsdir(savename("overall_run", name, "png")))



    layout["title"] = "Update Actor and Critic Networks<br><sub>Env: " * (env_cuda ? "GPU" : "CPU") * ", Agent: " * (agent_cuda ? "GPU" : "CPU") * ", Batch Size: " * string(batch) * "</sub>"

    trace1 = scatter(x = nodes, y = update_Actor ./1e9, mode="lines", name = "Update (Actor)")
    trace2 = scatter(x = nodes, y = update_Critic ./1e9, mode="lines", name = "Update (Critic)")

    p = plot([trace1, trace2], layout)
    display(p)
    savefig(p, plotsdir(savename("update_networks", name, "png")))
    


    layout["title"] = "Gradient Computations<br><sub>Env: " * (env_cuda ? "GPU" : "CPU") * ", Agent: " * (agent_cuda ? "GPU" : "CPU") * ", Batch Size: " * string(batch) * "</sub>"
    
    trace1 = scatter(x = nodes, y = grad_Actor ./1e9, mode="lines", name = "Gradient (Actor)")
    trace2 = scatter(x = nodes, y = grad_Critic ./1e9, mode="lines", name = "Gradient (Critic)")

    p = plot([trace1, trace2], layout)
    display(p)
    savefig(p, plotsdir(savename("gradients", name, "png")))



    layout["title"] = "Data transfer: CPU -> GPU<br><sub>Env: " * (env_cuda ? "GPU" : "CPU") * ", Agent: " * (agent_cuda ? "GPU" : "CPU") * ", Batch Size: " * string(batch) * "</sub>"

    trace1 = scatter(x = nodes, y = data_transfer ./1e9, mode="lines", name = "Gradient (Actor)")

    p = plot([trace1], layout)
    display(p)
    savefig(p, plotsdir(savename("data_transfer", name, "png")))



    Timing_dict[batch] = Dict()

    Timing_dict[batch]["no_nodes"] = nodes
    Timing_dict[batch]["overall_run"] = overall_run
    Timing_dict[batch]["inside_run"] = inside_run
    Timing_dict[batch]["policy_update"] = policy_update
    Timing_dict[batch]["env_calc"] = env_calc
    Timing_dict[batch]["prepare_data"] = prepare_data
end