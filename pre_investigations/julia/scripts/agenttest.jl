using DrWatson
@quickactivate "MicroGridSimWithRL"

using PProf, Profile
# using DifferentialEquations
# using Sundials
using Plots
# using LinearAlgebra
# using ControlSystems
using BenchmarkTools
using ReinforcementLearning
using Flux
using StableRNGs
using IntervalSets
using TimerOutputs
using JSON

# include(srcdir("collect_timing_results.jl"))
include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent.jl"))
include(srcdir("run_timed.jl"))

# collect_results = Dict()

global timer = TimerOutput()
# reset_timer!(timer::TimerOutput)

no_nodes = []
overall_run = []
inside_run = []
policy_update = []
env_calc = []
prepare_data = []


function train_ddpg(env_cuda::Bool, agent_cuda::Bool, num_nodes::Int, 
                    No_Episodes::Int)

    env_cuda = env_cuda
    agent_cuda = agent_cuda

    num_nodes = num_nodes

    CM = [ 0.  1.
            -1.  0.]

    CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(num_nodes) * ".json"))

    CM = reduce(hcat, CM_list[1])'
    CM = convert(Matrix{Int}, CM)

    parameters = Dict()
    # LC filter
    # parameters["source"] = [Dict("fltr" => "LC", "R" => 0.4, "L1" => 2.3e-3, "C" => 10e-6)]
    # parameters["cable"] = [Dict("R" => 0.722, "L" => 0.955e-3, "C" => 8e-09)]
    # parameters["load"] = [Dict("impedance" => "R", "R" => 14)]


    #nc = NodeConstructor(num_sources=1, num_loads=1, CM=CM, parameters=parameters)
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
    agent = create_agent(na, ns, agent_cuda)

    hook = TotalRewardPerEpisode()

    # No_Episodes = 50

    # run once to precompile
    run(agent, env, StopAfterEpisode(1), hook)
    
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
    timer
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
    append!(prepare_data, TimerOutputs.time(timer["Overall run"]["inside run"]["Agent prepare data"]))
end




# collect_results!(timer, 1)

nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30]

for i = 1:length(nodes)
    local_timer = train_ddpg(true, true, nodes[i], 5)
    collect_results!(local_timer, i)
end
# show(timer)



# @benchmark run(
#     agent,
#     env,
#     StopAfterEpisode(No_Episodes),
#     hook
# ) seconds = 120 evals = 2

# Plots.plot(nodes, overall_run)
Plots.plot(nodes, [overall_run, inside_run, policy_update],
    title = "Overall training time - on CPU",
    ylabel = "Time [ns]",
    xlabel = "No. of Nodes",
    label = ["overall run" "inside run" "policy update"],
    legend =:outertopright)

# plot(hook.rewards, 
#     title = "Total reward per episode",
#     xlabel = "Episodes",
#     ylabel = "Rewards",
#     legend = false)

# plot(PLoad,
#     title = "Power @ Load",
#     ylabel = "Power in Watts",
#     xlabel = "Time steps for 150 episodes [150 x 300]",
#     legend = false)

# plot(PLoad[end-300 : end],
#     title = "Power @ Load for the last episode",
#     ylabel = "Power in Watts",
#     xlabel = "Time steps",
#     legend = false
#     )

# Plots.savefig(p, plotsdir())

Timing_dict = Dict()
Timing_dict["no_nodes"] = nodes
Timing_dict["overall_run"] = overall_run
Timing_dict["inside_run"] = inside_run
Timing_dict["policy_update"] = policy_update
Timing_dict["env_calc"] = env_calc
Timing_dict["prepare_data"] = prepare_data