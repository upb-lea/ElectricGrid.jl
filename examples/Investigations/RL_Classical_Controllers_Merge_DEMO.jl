using ElectricGrid
#using ReinforcementLearning

"""
This scipt contains the content of the RL_Classical_Controllers_Merge_DEMO.ipynb notebook.
For comments and more documentation see refering notebookin examples/notebooks
"""

println("...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]

R_load, L_load, _, _ = ParallelLoadImpedance(50e3, 0.95, 230)

parameters =
Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "control_type" => "RL",
                        "fltr" => "L",
                        #"L1" => 0.0008,
                        ),
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "fltr" => "LC",
                        "control_type" => "classic",
                        "mode" => "Droop",),
                    ],
    "grid" => Dict{Any, Any}(
        "phase" => 3,
        "ramp_end" => 0.04,)
)



function reward_function(env, name = nothing)

    if name == "classic"
        return 0
    else
        p_ref = 1e3
        q_ref = 0e3

        pq0_ref = [p_ref; q_ref; 0]

        v_a_meas = env.x[findfirst(x -> x == "source1_v_C_cables_a", env.state_ids)]
        v_b_meas = env.x[findfirst(x -> x == "source1_v_C_cables_b", env.state_ids)]
        v_c_meas  = env.x[findfirst(x -> x == "source1_v_C_cables_c", env.state_ids)]
        i_a_meas = env.state[findfirst(x -> x == "source1_i_L1_a", env.state_ids)]
        i_b_meas = env.state[findfirst(x -> x == "source1_i_L1_b", env.state_ids)]
        i_c_meas  = env.state[findfirst(x -> x == "source1_i_L1_c", env.state_ids)]

        v_abc_meas = [v_a_meas; v_b_meas; v_c_meas]
        i_abc_meas = [i_a_meas; i_b_meas; i_c_meas]

        i_abc_ref = InvClarkeTransform(Inv_p_q_v(ClarkeTransform(v_abc_meas), pq0_ref))

        i_abc_ref = clamp.(i_abc_ref/env.nc.parameters["source"][1]["i_limit"], -0.95, 0.95)

        #i_abc_ref= [0; 0; 0]

        #println(i_abc_ref)

        # if any((abs.(env.state)).>=1)
        #     return -1
        # end

        r = 1 - 1/3*(sum((abs.(i_abc_ref - i_abc_meas)/2).^0.5))

        r += (0.2 - clamp(maximum(abs.(env.state)), 0.2, 1.0))

        return r * (1 - 0.99)/2
    end

end

function reference(t)
    if t < 0.04
        return [0.0, 0.0, 0.0]
    end

    θ = 2*pi*50*t
    θph = [θ; θ - 120π/180; θ + 120π/180]
    return +10 * cos.(θph) 
end

featurize_ddpg = function(state, env, name)
    if name == "ElectricGrid_ddpg_1"

        #state = state[findall(x -> split(x, "_")[2] == "i" , env.agent_dict["ElectricGrid_ddpg_1"]["state_ids"])]

        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        state = vcat(state, reference(env.t)/norm_ref)

        # θ = 2*pi*50*env.t

        # state_to_control_1 = env.state[findfirst(x -> x == "source1_i_L1_a", env.state_ids)]
        # state_to_control_2 = env.state[findfirst(x -> x == "source1_i_L1_b", env.state_ids)]
        # state_to_control_3 = env.state[findfirst(x -> x == "source1_i_L1_c", env.state_ids)]

        # state_to_control = [state_to_control_1, state_to_control_2, state_to_control_3]

        # Il_dq0 = DQ0Transform(state_to_control, θ)
        # state = vcat(state, Il_dq0)
    end
end

featurize_ddpg2 = function(state, env, name)
    if name == "ElectricGrid_ddpg_1"

        state = state[findall(x -> split(x, "_")[2] == "i" , env.agent_dict["ElectricGrid_ddpg_1"]["state_ids"])]

        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        state = vcat(state, reference(env.t)/norm_ref)
    end
end

function reward_function2(env, name = nothing)
    if name == "classic"
        return 0
        
    else
        state_to_control_1 = env.state[findfirst(x -> x == "source1_i_L1_a", env.state_ids)]
        state_to_control_2 = env.state[findfirst(x -> x == "source1_i_L1_b", env.state_ids)]
        state_to_control_3 = env.state[findfirst(x -> x == "source1_i_L1_c", env.state_ids)]

        state_to_control = [state_to_control_1, state_to_control_2, state_to_control_3]

        if any(abs.(state_to_control).>1)
            return -1
        else

            refs = reference(env.t)
            norm_ref = env.nc.parameters["source"][1]["i_limit"]          
            r = 1-1/3*(sum((abs.(refs/norm_ref - state_to_control)/2).^0.5))
            return r 
        end
    end

end

env = ElectricGridEnv(
    #CM =  CM,
    parameters = parameters,
    t_end = 1,
    reward_function = reward_function2,
    featurize = featurize_ddpg,
    action_delay = 0,
    verbosity = 0)


controllers = SetupAgents(env)

using FileIO, JLD2
#FileIO.save("controllers.jld2","controllers",controllers)
#controllers = FileIO.load("controllers.jld2","controllers");

learnhook = DataHook()

# Let it learn till 19_000 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function learn()
    num_steps = 50_000

    an_scheduler_loops = 20

    while true
        if length(controllers.hook.df[!,"reward"]) <= 1_500_000
            println("Steps so far: $(length(controllers.hook.df[!,"reward"]))")
            Learn(controllers, env, steps = num_steps, hook = learnhook)
        else
            an = 0.01 * exp10.(collect(LinRange(0.0, -10, an_scheduler_loops)))
            for i in 1:an_scheduler_loops
                controllers.agents["ElectricGrid_ddpg_1"]["policy"].policy.policy.act_noise = an[i]
                println("Steps so far: $(length(controllers.hook.df[!,"reward"]))")
                println("next action noise level: $(an[i])")
                Learn(controllers, env, steps = num_steps, hook = learnhook)
            end
        end
    end
end

function learn5()
    num_steps = 3_000_000

    global rewardresults1 = [[], [], [], [], []]
    rewardresults1 = convert(Vector{Vector{Float64}}, rewardresults1)
    for i in 1:5
        learnmore = true

        env = ElectricGridEnv(
            #CM =  CM,
            parameters = parameters,
            t_end = 1,
            reward_function = reward_function2,
            featurize = featurize_ddpg,
            action_delay = 0,
            verbosity = 0)

        controllers = SetupAgents(env)
        learnhook = DataHook()

        Learn(controllers, env, steps = num_steps, hook = learnhook)

        println(length(learnhook.df[!,"reward"]))

        rewardresults1[i] = convert(Vector{Float64}, learnhook.df[!,"reward"])
    end

    # global rewardresults2 = [[], [], [], [], []]
    # rewardresults2 = convert(Vector{Vector{Float64}}, rewardresults2)
    # for i in 1:5
    #     learnmore = true

    #     env = ElectricGridEnv(
    #         #CM =  CM,
    #         parameters = parameters,
    #         t_end = 1,
    #         reward_function = reward_function2,
    #         featurize = featurize_ddpg2,
    #         action_delay = 0,
    #         verbosity = 0)

    #     controllers = SetupAgents(env)
    #     learnhook = DataHook()

    #     Learn(controllers, env, steps = num_steps, hook = learnhook)

    #     println(length(learnhook.df[!,"reward"]))

    #     rewardresults2[i] = convert(Vector{Float64}, learnhook.df[!,"reward"])
    # end
end

#FileIO.save("rewardresults1.jld2","rewardresults1",rewardresults1)
#FileIO.save("rewardresults2.jld2","rewardresults2",rewardresults2)
#rewardresults1 = FileIO.load("rewardresults1.jld2","rewardresults1")
#rewardresults2 = FileIO.load("rewardresults2.jld2","rewardresults2")

#good in rewardresults1: 2, 3, 4
#good in rewardresults2: 2, 3

using Statistics
using PlotlyJS
function plot_rewardresults(n = 0)
    global plotresults = Float64[]
    global plotresults_std = Float64[]
    if n > 0
        xx = [i*1000+1 for i in 1:3000]
        for i in 1:3000
            temp = []
            if n == 1
                for j = 2:4
                    append!(temp,rewardresults1[j][(i-1)*1000+1:i*1000])
                end
            else
                for j = 2:3
                    append!(temp,rewardresults2[j][(i-1)*1000+1:i*1000])
                end
            end

            push!(plotresults, mean(temp))
            push!(plotresults_std, std(temp))
            #push!(plotresults_std, std([mean(temp[1:1000]),mean(temp[1001:2000]),mean(temp[2001:3000])]))
        end
    else
        rewardresults = convert(Vector{Float64}, controllers.hook.df[!,"reward"])
        global batch_size = Int(floor(length(rewardresults)/3000))
        xx = [i*batch_size+1 for i in 1:3000]
        for i in 1:3000
            temp = []
            append!(temp,rewardresults[(i-1)*batch_size+1:i*batch_size])
            push!(plotresults, mean(temp))
            push!(plotresults_std, std(temp))
        end
    end

    rl = Layout(
        plot_bgcolor = "white",
        font=attr(
            family="Arial",
            size=16,
            color="black"
        ),
        yaxis=attr(
            title="reward",
            showline=true,
            linewidth=2,
            linecolor="black",
            showgrid=true,
            gridwidth=1,
            gridcolor="LightGrey",
            ),
        xaxis=attr(
            title="time steps",
            showline=true,
            linewidth=2,
            linecolor="black",
            showgrid=true,
            gridwidth=1,
            gridcolor="LightGrey",
            ),
    )

    plot([
        scatter(x=xx, y=plotresults.+plotresults_std, mode="lines", line=attr(width=0.0)),
        scatter(x=xx, y=plotresults.-plotresults_std, mode="none", fillcolor="rgba(111, 120, 219, 0.3)", fill="tonexty", showlegend=false),
        scatter(x=xx, y=plotresults, mode="lines", line_color="indigo", showlegend=false),
    ], rl)
end

hook = DataHook(collect_state_ids = env.state_ids,
                collect_action_ids = env.action_ids)

Simulate(controllers, env, hook=hook)


RenderHookResults(hook = hook,
                    states_to_plot  = env.state_ids,
                    actions_to_plot = env.action_ids,
                    plot_reward=true)

#RenderHookResults(hook = hook, states_to_plot  = ["source1_i_L1_a", "source1_i_L1_b", "source1_i_L1_c"], actions_to_plot = [])
#RenderHookResults(hook = hook, states_to_plot  = [], actions_to_plot = ["source1_u_a", "source1_u_b", "source1_u_c"])

println("...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")



function plottest()
    a = [1.0, 2.0, 2.0, 3.0]
    b = [0.3, 0.1, 0.4, 0.1]

    plot([
        scatter(y=a.+b, mode="lines", line=attr(width=0.0), showlegend=false),
        scatter(y=a.-b, mode="none", fillcolor="rgba(111, 120, 219, 0.3)", fill="tonexty", showlegend=false),
        scatter(y=a, mode="lines", line_color="indigo"),
    ])
end