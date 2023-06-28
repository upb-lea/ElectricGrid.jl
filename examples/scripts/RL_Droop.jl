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

#-------------------------------------------------------------------------------
total_num_RL = 1
total_phases = 3 # always has to be 3 when using pll 3ph

# initialise matrices
θpll = zeros(total_num_RL, total_phases, 2) # two past terms
fpll = 50*ones(total_num_RL, total_phases, 4) # four past terms
pll_err_t = zeros(total_num_RL, total_phases)
pll_err = zeros(total_num_RL, total_phases, 3)

global θpll, fpll, pll_err_t, pll_err

# fsys = system frequency = 50 Hz
# ts = step interval
# V_filt_cap_new = instantaneous voltage over capacitor at poc
# phase = 1
# fpll[num_RL_source, phase, end = newest]
# θpll[num_RL_source, phase, end = newest]
# pll_err_t[num_source, phase]
# pll_err[num_source, phase, :]
#-------------------------------------------------------------------------------
RL_pll = PLL(50, 0.0001)

#= for i in 1:1000
    v1 = cos(50*2π*i)
    v2 = cos(50*2π*i + 120) # for checking pll,check 120 angles + or -
    v3 = cos(50*2π*i - 120)
    fpll, θpll = RL_pll([v1, v2, v3])
end =#


function reference(env)

    if env.t < 0.04
        return [0.0, 0.0, 0.0]
    end

    #V_filt_cap_new = env.state_abc_v
    vC_1 = env.state[findfirst(x -> x == "source1_v_C_cables_a", env.state_ids)]
    vC_2 = env.state[findfirst(x -> x == "source1_v_C_cables_b", env.state_ids)]
    vC_3 = env.state[findfirst(x -> x == "source1_v_C_cables_c", env.state_ids)]

    V_filt_cap_new = [vC_1, vC_2, vC_3]

    #= θ_fudge = 2*pi*50*env.t
    θ_f = [θ; θ - 120π/180; θ + 120π/180]
    V_filt_cap_new =cos.(θ_f)
    =#
    fpll, θpll = RL_pll(V_filt_cap_new) #TODO we do not know the env here

    θph = θpll[ :, end]

    return +10 * cos.(θph .+ π/2)
end

featurize_ddpg = function(state, env, name)
    if name == "ElectricGrid_ddpg_1"

        #state = state[findall(x -> split(x, "_")[2] == "i" , env.agent_dict["ElectricGrid_ddpg_1"]["state_ids"])]

        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        state = vcat(state, reference(env)/norm_ref)

        # θ = 2*pi*50*env.t

        # state_to_control_1 = env.state[findfirst(x -> x == "source1_i_L1_a", env.state_ids)]
        # state_to_control_2 = env.state[findfirst(x -> x == "source1_i_L1_b", env.state_ids)]
        # state_to_control_3 = env.state[findfirst(x -> x == "source1_i_L1_c", env.state_ids)]

        # state_to_control = [state_to_control_1, state_to_control_2, state_to_control_3]

        # Il_dq0 = DQ0Transform(state_to_control, θ)
        # state = vcat(state, Il_dq0)
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

            refs = reference(env)
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

hook = DataHook(collect_state_ids = env.state_ids,
                collect_action_ids = env.action_ids)

#Simulate(controllers, env, hook=hook);


#RenderHookResults(hook = hook, states_to_plot  = env.state_ids, actions_to_plot = env.action_ids, plot_reward=true)

println("...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")
