using ElectricGrid
using LinearAlgebra

mutable struct PLL

    total_phases::Int64 # always has to be 3 when using pll 3ph

    # initialise matrices
    θpll::Array{Float64} # two past terms
    fpll::Array{Float64} # four past terms
    pll_err_t::Array{Float64}
    pll_err::Array{Float64}

    Ki::Float64
    Kp::Float64

    ts::Float64
    fsys::Float64

    #=     function PLL(total_phases::Int64, θpll::Array{Float64},
                fpll::Array{Float64}, pll_err_t::Array{Float64}, pll_err::Array{Float64}, Ki::Float64, Kp::Float64, ts::Float64)

        new(total_phases, θpll, fpll, pll_err_t, pll_err, Ki, Kp, ts)
    end
    =#
    function PLL(fsys, ts; ωn = fsys + 20, ξ = 0.35)

        total_phases = 3 # always has to be 3 when using pll 3ph

        # initialise matrices
        θpll = zeros(total_phases, 2) # two past terms
        fpll = fsys*ones(total_phases, 4) # four past terms
        pll_err_t = zeros(total_phases)
        pll_err = zeros(total_phases, 3)

        # tuning
        Ki = (ωn)^2 # tuning -> to struct
        Kp = (ξ)*2*sqrt(Ki) # tuning

        #PLL(total_phases, θpll, fpll, pll_err_t, pll_err, Ki, Kp, ts)
        new(total_phases, θpll, fpll, pll_err_t, pll_err, Ki, Kp, ts, fsys)
    end

    function (pll::PLL)(V_filt_cap_new)

        θ = pll.θpll[1, end]

        v_αβγ = ClarkeTransform(V_filt_cap_new)

        if norm(v_αβγ) != 0
            v_αβγ = sqrt(3)*v_αβγ./norm(v_αβγ)
        end

        err_new = v_αβγ[2]*cos(θ) - v_αβγ[1]*sin(θ)

        f_new, err_t_new, err_int =
        PIController(err_new, pll.pll_err[1, :], pll.pll_err_t[1], pll.Kp, pll.Ki, pll.ts, bias = pll.fsys, max_t_err = 0.00015)

        θ = ThirdOrderIntegrator(θ, pll.ts, 2π*[pll.fpll[1, 2:end]; f_new])

        roll(pll.fpll)
        pll.fpll[1, end] = f_new[1]
        pll.fpll[2, end] = f_new[1]
        pll.fpll[3, end] = f_new[1]

        roll(pll.θpll)
        pll.θpll[1, end] = (θ.%(2π))
        pll.θpll[2, end] = ((θ - 120π/180).%(2π))
        pll.θpll[3, end] = ((θ + 120π/180).%(2π))

        pll.pll_err_t[1] = err_t_new[1]
        pll.pll_err[1, :] = err_int

        return pll.fpll, pll.θpll
    end
end

println("...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

CM = [ 0. 1.
        -1. 0.]

R_load, L_load, _, _ = ParallelLoadImpedance(50e3, 0.95, 230)

parameters =
Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "control_type" => "RL",
                        "fltr" => "LC",
                        "osc" => true,
                        #"L1" => 0.0008,
                        ),
                    ],
    "load"   => Any[
        Dict{Any, Any}(
            "impedance" => "RL",
            "R" => R_load,
            "L" => L_load,
            "v_limit" => 1e4,
            "i_limit" => 1e4)
        ],
    "grid" => Dict{Any, Any}(
        "phase" => 3,
        "ramp_end" => 0.04,)
)

#-------------------------------------------------------------------------------
total_num_RL = 1
total_phases = 3 # always has to be 3 when using pll 3ph

RL_pll = PLL(50, 0.0001)

function reference(env)

    if env.steps == 0
        RL_pll.θpll = zeros(total_phases, 2) # two past terms
        RL_pll.fpll = RL_pll.fsys*ones(total_phases, 4) # four past terms
        RL_pll.pll_err_t = zeros(total_phases)
        RL_pll.pll_err = zeros(total_phases, 3)
    end

    #V_filt_cap_new = env.state_abc_v
    vC_1 = env.state[findfirst(x -> x == "source1_v_C_cables_a", env.state_ids)]
    vC_2 = env.state[findfirst(x -> x == "source1_v_C_cables_b", env.state_ids)]
    vC_3 = env.state[findfirst(x -> x == "source1_v_C_cables_c", env.state_ids)]

    V_filt_cap_new = [vC_1, vC_2, vC_3]

    fpll, θpll = RL_pll(V_filt_cap_new)


    if env.t < 0.04
        return [0.0, 0.0, 0.0]
    end

    θph = θpll[ :, end]

    return +10 * cos.(θph .+ π/2)
end

global ref
featurize_ddpg = function(state, env, name)
    if name == "ElectricGrid_ddpg_1"

        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        global ref = reference(env)
        state = vcat(state, ref/norm_ref)

    end
end

function reward_function(env, name = nothing)
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
            refs = ref
            norm_ref = env.nc.parameters["source"][1]["i_limit"]
            r = 1-1/3*(sum((abs.(refs/norm_ref - state_to_control)/2).^0.5))
            #push!(ref_a, refs[1]/norm_ref)
            return r
        end
    end

end

env = ElectricGridEnv(
    CM =  CM,
    parameters = parameters,
    t_end = 1,
    reward_function = reward_function,
    featurize = featurize_ddpg,
    action_delay = 0,
    verbosity = 0)


controllers = SetupAgents(env)

learnhook = DataHook()

function learn1()
    steps_total = 1_500_000

    steps_loop = 50_000

    Learn(controllers, env, steps = steps_loop)

    while length(controllers.hook.df[!,"reward"]) <= steps_total

        println("Steps so far: $(length(controllers.hook.df[!,"reward"]))")
        Learn(controllers, env, steps = steps_loop, hook = learnhook)

    end

end

# second training phase with action noise scheduler
function learn2()
    num_steps = 50_000

    an_scheduler_loops = 20


    for j in 1:10
        an = 0.01 * exp10.(collect(LinRange(0.0, -10, an_scheduler_loops)))
        for i in 1:an_scheduler_loops
            controllers.agents["ElectricGrid_ddpg_1"]["policy"].policy.policy.act_noise = an[i]
            println("Steps so far: $(length(controllers.hook.df[!,"reward"]))")
            println("next action noise level: $(an[i])")
            Learn(controllers, env, steps = num_steps, hook = learnhook)
        end
    end
end


#learn1()

#plot_rewardresults(convert(Vector{Float64}, controllers.hook.df[!,"reward"]))


hook = DataHook(collect_state_ids = env.state_ids,
                collect_action_ids = env.action_ids)

#Simulate(controllers, env, hook=hook);

using PlotlyJS
#Plot([scatter(y=hook.df[!,"source1_v_C_cables_a"], name="source1_v_C_cables_a"), scatter(y=hook.df[!,"source1_i_L1_a"], name="source1_i_L1_a"), scatter(y=ref_a.*700, name="ref_a"), scatter(y=fpll_a[1:end], name="fpll_a"), scatter(y=θpll_a[1:end], name="θpll_a")])

#RenderHookResults(hook = hook, states_to_plot  = env.state_ids, actions_to_plot = env.action_ids, plot_reward=true)

println("...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")
