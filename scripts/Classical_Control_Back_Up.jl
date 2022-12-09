include("./env.jl")

mutable struct Classical_Controls

    #---------------------------------------------------------------------------
    # Physical Electrical Parameters

    Vdc::Vector{Float64} # DC link voltage
    Vrms::Vector{Float64} # nominal output voltage
    S::Vector{Float64} # rated nominal apparent power
    P::Vector{Float64} # rated nominal active power
    Q::Vector{Float64} # rated nominal active power

    # IGBT / Switch Ratings
    i_max::Vector{Float64}
    v_max::Vector{Float64}

    # Filter values
    Lf::Vector{Float64} # Filter values
    Cf::Vector{Float64}
    Rf::Vector{Float64}

    #---------------------------------------------------------------------------
    # Measurements

    T_eval::Int64 # number of periods to save to memory (for averaging, e.g. RMS)
    T_sp_rms::Float64
    V_ph::Array{Float64}
    I_ph::Array{Float64}

    p_q_inst::Array{Float64}
    p_inst::Array{Float64}
    Pm::Array{Float64}
    Qm::Array{Float64}

    #---------------------------------------------------------------------------
    # General System & Control

    Modes::Dict{String, Int64}
    Source_Modes::Vector{String}

    num_sources::Int64
    phases::Int64

    f_cntr::Float64 # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz
    fsys::Float64
    θsys::Float64
    ts::Float64 # s, sampling timestep
    N::Int64
    delay::Int64
    steps::Int64

    V_poc_loc::Matrix{Int64} # the position in the state vector where the POC Voltage is measured
    I_poc_loc::Matrix{Int64}
    I_inv_loc::Matrix{Int64}

    #---------------------------------------------------------------------------
    # Phase Locked Loops

    pll_err::Array{Float64}
    pll_err_t::Matrix{Float64}
    vd::Array{Float64}
    qvd::Array{Float64}
    fpll::Array{Float64}
    θpll::Array{Float64}

    #---------------------------------------------------------------------------
    # Interface (e.g. filters)

    V_filt_poc::Array{Float64}
    V_filt_inv::Array{Float64}
    I_filt_poc::Array{Float64}
    I_filt_inv::Array{Float64}
    p_q_filt::Array{Float64}

    #---------------------------------------------------------------------------
    # Current Controller

    Gi_cl::Array{TransferFunction} # Closed Loop transfer function

    I_dq0::Array{Float64} # DQ0 of I_filt_inv
    I_ref_dq0::Array{Float64}
    I_ref::Array{Float64}

    I_err::Array{Float64}
    I_err_t::Matrix{Float64}
    I_kp::Vector{Float64}
    I_ki::Vector{Float64}

    Vd_abc_new::Array{Float64} # output action

    #---------------------------------------------------------------------------
    # Voltage Controller

    Gv_cl::Array{TransferFunction} # Closed Loop transfer function

    V_δ_set::Matrix{Float64} # set points when also in swing mode
    V_pu_set::Matrix{Float64}

    V_dq0::Array{Float64} # DQ0 of V_filt_poc
    V_ref_dq0::Array{Float64}
    V_ref::Array{Float64}

    V_err::Array{Float64}
    V_err_t::Matrix{Float64}
    V_kp::Vector{Float64}
    V_ki::Vector{Float64}

    I_lim::Matrix{Float64} # interim state before passing through limiter

    #---------------------------------------------------------------------------
    # Droop (Classical, Synchronverter, and VSG) Mode

    Δfmax::Float64 # The drop (increase) in frequency that causes a 100% increase (decrease) in power
    ΔEmax::Float64 # The drop (increase) in rms voltage that causes a 100% increase (decrease) in reactive power (from nominal)
    τv::Float64  # time constant of the voltage loop
    τf::Float64  # time constant of the frequency-droop loop

    D::Matrix{Float64} # Droop coefficients
    ω_droop::Array{Float64}
    θ_droop::Array{Float64}

    #---------------------------------------------------------------------------
    # PQ Mode

    pq0_set::Matrix{Float64}

    #---------------------------------------------------------------------------
    # Synchronverter Mode

    J_sync::Vector{Float64} # Virtual Mass Moment of Inertia
    K_sync::Vector{Float64} # Reactive droop integrator gain
    ΔT_t::Vector{Float64}

    α_sync::Matrix{Float64}
    ω_sync::Matrix{Float64}
    θ_sync::Matrix{Float64}
    Δω_sync::Matrix{Float64}
    eq::Matrix{Float64}
    Mfif::Matrix{Float64}

    function Classical_Controls(Vdc::Vector{Float64}, Vrms::Vector{Float64},
        S::Vector{Float64}, P::Vector{Float64}, Q::Vector{Float64},
        i_max::Vector{Float64}, v_max::Vector{Float64},
        Lf::Vector{Float64}, Cf::Vector{Float64}, Rf::Vector{Float64},
        T_eval::Int64, T_sp_rms::Float64, V_ph::Array{Float64}, I_ph::Array{Float64}, 
        p_q_inst::Array{Float64}, p_inst::Array{Float64}, Pm::Array{Float64}, Qm::Array{Float64},
        Modes::Dict{String, Int64}, Source_Modes::Vector{String},
        num_sources::Int64, phases::Int64,
        f_cntr::Float64, fsys::Float64, θsys::Float64,
        ts::Float64, N::Int64, delay::Int64, steps::Int64,
        V_poc_loc::Matrix{Int64}, I_poc_loc::Matrix{Int64}, I_inv_loc::Matrix{Int64},
        pll_err::Array{Float64}, pll_err_t::Matrix{Float64},
        vd::Array{Float64}, qvd::Array{Float64},
        fpll::Array{Float64}, θpll::Array{Float64},
        V_filt_poc::Array{Float64}, V_filt_inv::Array{Float64},
        I_filt_poc::Array{Float64}, I_filt_inv::Array{Float64},
        p_q_filt::Array{Float64},
        Gi_cl::Array{TransferFunction},
        I_dq0::Array{Float64}, I_ref_dq0::Array{Float64}, I_ref::Array{Float64},
        I_err::Array{Float64}, I_err_t::Matrix{Float64},
        I_kp::Vector{Float64}, I_ki::Vector{Float64},
        Vd_abc_new::Array{Float64},
        Gv_cl::Array{TransferFunction}, V_δ_set::Matrix{Float64}, V_pu_set::Matrix{Float64},
        V_dq0::Array{Float64}, V_ref_dq0::Array{Float64}, V_ref::Array{Float64},
        V_err::Array{Float64}, V_err_t::Matrix{Float64},
        V_kp::Vector{Float64}, V_ki::Vector{Float64},
        I_lim::Matrix{Float64},
        Δfmax::Float64, ΔEmax::Float64, τv::Float64, τf::Float64,
        D::Matrix{Float64}, ω_droop::Array{Float64}, θ_droop::Array{Float64},
        pq0_set::Matrix{Float64},
        J_sync::Vector{Float64}, K_sync::Vector{Float64}, ΔT_t::Vector{Float64},
        α_sync::Matrix{Float64}, ω_sync::Matrix{Float64}, θ_sync::Matrix{Float64},
        Δω_sync::Matrix{Float64}, eq::Matrix{Float64}, Mfif::Matrix{Float64})

        new(Vdc, Vrms,
        S, P, Q,
        i_max, v_max,
        Lf, Cf, Rf,
        T_eval, T_sp_rms, V_ph, I_ph, 
        p_q_inst, p_inst, Pm, Qm,
        Modes, Source_Modes,
        num_sources, phases,
        f_cntr, fsys, θsys,
        ts, N, delay, steps,
        V_poc_loc, I_poc_loc, I_inv_loc,
        pll_err, pll_err_t,
        vd, qvd,
        fpll, θpll,
        V_filt_poc, V_filt_inv,
        I_filt_poc, I_filt_inv,
        p_q_filt,
        Gi_cl,
        I_dq0, I_ref_dq0, I_ref,
        I_err, I_err_t,
        I_kp, I_ki,
        Vd_abc_new,
        Gv_cl,V_δ_set, V_pu_set,
        V_dq0, V_ref_dq0, V_ref,
        V_err, V_err_t,
        V_kp, V_ki,
        I_lim,
        Δfmax, ΔEmax, τv, τf,
        D, ω_droop, θ_droop,
        pq0_set,
        J_sync, K_sync, ΔT_t,
        α_sync, ω_sync, θ_sync,
        Δω_sync, eq, Mfif)
    end

    function Classical_Controls(t_final, f_cntr, num_sources; delay = 1, phases = 3)

        #---------------------------------------------------------------------------
        # Physical Electrical Parameters

        Vdc = Array{Float64, 1}(undef, num_sources)
        Vdc = fill!(Vdc, 800.0)
        Vrms = Array{Float64, 1}(undef, num_sources)
        Vrms = fill!(Vrms, 230.0)

        S = Array{Float64, 1}(undef, num_sources)
        S = fill!(S, 50.0e3)
        P = Array{Float64, 1}(undef, num_sources)
        P = fill!(P, 40.0e3)
        Q = Array{Float64, 1}(undef, num_sources)
        Q = fill!(Q, 30.0e3)

        i_max = Array{Float64, 1}(undef, num_sources)
        i_max = fill!(i_max, sqrt(2)*S[1]/(Vrms[1]*3))
        v_max = Array{Float64, 1}(undef, num_sources)
        v_max = fill!(v_max, 1200/(2*sqrt(2)))

        Lf = Array{Float64, 1}(undef, num_sources)
        Lf = fill!(Lf, 0)
        Cf = Array{Float64, 1}(undef, num_sources)
        Cf = fill!(Cf, 0)
        Rf = Array{Float64, 1}(undef, num_sources)
        Rf = fill!(Rf, 0.4)

        #---------------------------------------------------------------------------
        # Measurements

        t_final = convert(Float64, t_final)
        ts = 1/f_cntr
        t = 0:ts:t_final

        steps = 0
        fsys = 50.0
        θsys = 0.0

        T_eval = 1 #number of periods to average over
        N = length(t)

        T_sp_rms = 10*fsys #samples in a second for rms calcs, x*fsys = x samples in a cycle

        # RMS Phasors
        V_ph = Array{Float64, 4}(undef, num_sources, phases, 3, N)
        V_ph = fill!(V_ph, 0)
        I_ph = Array{Float64, 4}(undef, num_sources, phases, 3, N)
        I_ph = fill!(I_ph, 0)

        # Instantaneous Real, Imaginary, and Zero powers
        p_q_inst = Array{Float64, 3}(undef, num_sources, phases, N)
        p_q_inst = fill!(p_q_inst, 0)

        p_inst = Array{Float64, 3}(undef, num_sources, phases, N) # instantaneous power at PCC
        p_inst = fill!(p_inst, 0)

        Pm = Array{Float64, 3}(undef, num_sources, phases+1, N) # 4th column is total
        Pm = fill!(Pm, 0)
        Qm = Array{Float64, 3}(undef, num_sources, phases+1, N) # 4th column is total
        Qm = fill!(Qm, 0)

        #---------------------------------------------------------------------------
        # General System & Control

        Modes = Dict("Swing Mode" => 1, "Voltage Control Mode" => 2, "PQ Control Mode" => 3,
        "Droop Control Mode" => 4, "Synchronverter Mode" => 5, "Self-Synchronverter Mode" => 6)

        Source_Modes = Array{String, 1}(undef, num_sources)
        Source_Modes = fill!(Source_Modes, "Voltage Control Mode")

        V_poc_loc = Array{Int64, 2}(undef, phases, num_sources)
        V_poc_loc = fill!(V_poc_loc, 1)
        I_poc_loc = Array{Int64, 2}(undef, phases, num_sources)
        I_poc_loc = fill!(I_poc_loc, 1)
        I_inv_loc = Array{Int64, 2}(undef, phases, num_sources)
        I_inv_loc = fill!(I_inv_loc, 1)

        #---------------------------------------------------------------------------
        # Phase Locked Loops

        vd = Array{Float64, 3}(undef, num_sources, phases, N) #3 for 3 phases
        vd = fill!(vd, 0.0)
        qvd = Array{Float64, 3}(undef, num_sources, phases, N) #3 for 3 phases
        qvd = fill!(qvd, 0.0)

        fpll = Array{Float64, 3}(undef, num_sources, phases, N) #3 for 3 phases
        fpll = fill!(fpll, fsys)

        θpll = Array{Float64, 3}(undef, num_sources, phases, N)
        θpll = fill!(θpll, 0)

        pll_err = Array{Float64, 3}(undef, num_sources, phases, 3) # 3 phases and 3rd order integration
        pll_err = fill!(pll_err, 0)
        pll_err_t = Array{Float64, 2}(undef, num_sources, phases) # PLL total integrated error
        pll_err_t = fill!(pll_err_t, 0)

        #---------------------------------------------------------------------------
        # Interface (e.g. Digital Signal Processing filters)

        V_filt_poc = Array{Float64, 3}(undef, num_sources, phases, N)
        V_filt_poc = fill!(V_filt_poc, 0)
        V_filt_inv = Array{Float64, 3}(undef, num_sources, phases, N)
        V_filt_inv = fill!(V_filt_inv, 0)

        I_filt_poc = Array{Float64, 3}(undef, num_sources, phases, N)
        I_filt_poc = fill!(I_filt_poc, 0)
        I_filt_inv = Array{Float64, 3}(undef, num_sources, phases, N)
        I_filt_inv = fill!(I_filt_inv, 0)

        p_q_filt = Array{Float64, 3}(undef, num_sources, phases, N)
        p_q_filt  = fill!(p_q_filt, 0)

        #---------------------------------------------------------------------------
        # Current Controller

        Gi_cl = Array{TransferFunction, 1}(undef, num_sources)

        I_dq0 = Array{Float64, 3}(undef, num_sources, phases, N)
        I_dq0 = fill!(I_dq0, 0)
        I_ref_dq0 = Array{Float64, 3}(undef, num_sources, phases, N)
        I_ref_dq0 = fill!(I_ref_dq0, 0)
        I_ref = Array{Float64, 3}(undef, num_sources, phases, N)
        I_ref = fill!(I_ref, 0)

        # Current Integrations
        I_err = Array{Float64, 3}(undef, num_sources, phases, 3) # 3 phases and 3rd order integration
        I_err = fill!(I_err, 0)
        I_err_t = Array{Float64, 2}(undef, num_sources, phases)
        I_err_t = fill!(I_err_t, 0)
        I_kp = Array{Float64, 1}(undef, num_sources)
        I_kp = fill!(I_kp, 0.5)
        I_ki = Array{Float64, 1}(undef, num_sources)
        I_ki = fill!(I_ki, 25)

        Vd_abc_new = Array{Float64, 3}(undef, num_sources, phases, N + delay)
        Vd_abc_new = fill!(Vd_abc_new, 0)

        #---------------------------------------------------------------------------
        # Voltage Controller

        V_δ_set = Array{Float64, 2}(undef, num_sources, phases)
        V_δ_set = fill!(V_δ_set, 0.0)
        V_pu_set = Array{Float64, 2}(undef, num_sources, phases)
        V_pu_set = fill!(V_pu_set, 1.0)

        Gv_cl = Array{TransferFunction, 1}(undef, num_sources)

        V_dq0 = Array{Float64, 3}(undef, num_sources, phases, N)
        V_dq0 = fill!(V_dq0, 0)
        V_ref_dq0 = Array{Float64, 3}(undef, num_sources, phases, N)
        V_ref_dq0 = fill!(V_ref_dq0, 0)

        V_ref = Array{Float64, 3}(undef, num_sources, phases, N)
        V_ref = fill!(V_ref, 0)

        # Voltage Integrations
        V_err = Array{Float64, 3}(undef, num_sources, phases, 3) # 3 phases and 3rd order integration
        V_err = fill!(V_err, 0)
        V_err_t = Array{Float64, 2}(undef, num_sources, phases)
        V_err_t = fill!(V_err_t, 0)
        V_kp = Array{Float64, 1}(undef, num_sources)
        V_kp = fill!(V_kp, 0.01)
        V_ki = Array{Float64, 1}(undef, num_sources)
        V_ki = fill!(V_ki, 20)

        I_lim = Array{Float64, 2}(undef, num_sources, 3)
        I_lim = fill!(I_lim, 0)

        #---------------------------------------------------------------------------
        # Droop (Classical, Synchronverter, and VSG) Mode

        #=
            Typical values for the frequency droop are a 100% increase in power for a
            frequency decrease between 3% and 5% (from nominal values)
        =#

        Δfmax = 0.03*fsys/100 # Hz # The drop in frequency, Hz, which will cause a 100% increase in active power
        ΔEmax = 0.05*Vrms[1]/100 # V # The drop in rms voltage, which will cause a 100% decrease in reactive power
        τv = 0.02 # time constant of the voltage loop
        τf = 0.02

        D = Array{Float64, 2}(undef, num_sources, 2)
        D[:,1] = fill!(D[:,1], 2π*Δfmax/P[1])
        D[:,2] = fill!(D[:,2], ΔEmax/Q[1])

        ω_droop = Array{Float64, 3}(undef, num_sources, phases, N) #3 for 3 phases
        ω_droop = fill!(ω_droop, fsys*2π)

        θ_droop = Array{Float64, 3}(undef, num_sources, phases, N)
        θ_droop = fill!(θ_droop, 0)

        #---------------------------------------------------------------------------
        # PQ Mode

        pq0_set = Array{Float64, 2}(undef, num_sources, 3) # real, imaginary, and zero power set points
        pq0_set = fill!(pq0_set, 0)

        #---------------------------------------------------------------------------
        # Synchronverter Mode

        J_sync = Array{Float64, 1}(undef, num_sources)
        J_sync = fill!(J_sync, 0)
        K_sync = Array{Float64, 1}(undef, num_sources)
        K_sync = fill!(K_sync, 0)
        ΔT_t = Array{Float64, 1}(undef, num_sources)
        ΔT_t = fill!(ΔT_t, 0)

        α_sync = Array{Float64, 2}(undef, num_sources, N)
        α_sync = fill!(α_sync, 0)
        ω_sync = Array{Float64, 2}(undef, num_sources, N)
        ω_sync = fill!(ω_sync, fsys*2π)
        θ_sync = Array{Float64, 2}(undef, num_sources, N)
        θ_sync = fill!(θ_sync, 0)
        Δω_sync = Array{Float64, 2}(undef, num_sources, N)
        Δω_sync = fill!(Δω_sync, 0)
        eq = Array{Float64, 2}(undef, num_sources, N)
        eq = fill!(eq, 0)
        Mfif = Array{Float64, 2}(undef, num_sources, N)
        Mfif = fill!(Mfif, 0)

        Classical_Controls(Vdc, Vrms,
        S, P, Q,
        i_max, v_max,
        Lf, Cf, Rf,
        T_eval, T_sp_rms, V_ph, I_ph, 
        p_q_inst, p_inst, Pm, Qm,
        Modes, Source_Modes,
        num_sources, phases,
        f_cntr, fsys, θsys,
        ts, N, delay, steps,
        V_poc_loc, I_poc_loc, I_inv_loc,
        pll_err, pll_err_t,
        vd, qvd,
        fpll, θpll,
        V_filt_poc, V_filt_inv,
        I_filt_poc, I_filt_inv,
        p_q_filt,
        Gi_cl,
        I_dq0, I_ref_dq0, I_ref,
        I_err, I_err_t,
        I_kp, I_ki,
        Vd_abc_new,
        Gv_cl,V_δ_set, V_pu_set,
        V_dq0, V_ref_dq0, V_ref,
        V_err, V_err_t,
        V_kp, V_ki,
        I_lim,
        Δfmax, ΔEmax, τv, τf,
        D, ω_droop, θ_droop,
        pq0_set,
        J_sync, K_sync, ΔT_t,
        α_sync, ω_sync, θ_sync,
        Δω_sync, eq, Mfif)
    end
end

Base.@kwdef mutable struct Classical_Policy <: AbstractPolicy

    n_actions = 1
    t_final = 0.04
    fs = 10e3
    num_sources = 1
    delay = 0
    phases = 3

    action_space::Space{Vector{ClosedInterval{Float64}}} = Space([ -1.0..1.0 for i = 1:n_actions], )
    Source::Classical_Controls = Classical_Controls(t_final, fs, num_sources, delay = 0, phases = phases)

    state_ids = nothing
    action_ids = nothing
end

function (Animo::Classical_Policy)(env::SimEnv, name::Union{String, Nothing} = nothing)

    if isnothing(name)
        Action = Classical_Control(Animo, env)
    else
        Action = Classical_Control(Animo, env, name)
    end

    return Action    
end

function (Animo::Classical_Policy)(::PostEpisodeStage, ::AbstractEnv)

    Source = Animo.Source
    t = 0:Source.ts:Animo.t_final

    θt = (2*π*Source.fsys*t).%(2π)

    Source.steps = 0

    Source.V_ph = fill!(Source.V_ph, 0)
    Source.I_ph = fill!(Source.I_ph, 0)

    Source.p_q_inst = fill!(Source.p_q_inst, 0)

    Source.p_inst = fill!(Source.p_inst, 0)

    Source.Pm = fill!(Source.Pm, 0)
    Source.Qm = fill!(Source.Qm, 0)

    Source.vd = fill!(Source.vd, 0.0)
    Source.qvd = fill!(Source.qvd, 0.0)

    Source.fpll = fill!(Source.fpll, Source.fsys)

    Source.θpll = fill!(Source.θpll, 0)

    Source.pll_err = fill!(Source.pll_err, 0)
    Source.pll_err_t = fill!(Source.pll_err_t, 0)

    Source.V_filt_poc = fill!(Source.V_filt_poc, 0)
    Source.V_filt_inv = fill!(Source.V_filt_inv, 0)

    Source.I_filt_poc = fill!(Source.I_filt_poc, 0)
    Source.I_filt_inv = fill!(Source.I_filt_inv, 0)

    Source.p_q_filt  = fill!(Source.p_q_filt, 0)

    Source.I_dq0 = fill!(Source.I_dq0, 0)
    Source.I_ref_dq0 = fill!(Source.I_ref_dq0, 0)
    Source.I_ref = fill!(Source.I_ref, 0)

    Source.I_err = fill!(Source.I_err, 0) #
    Source.I_err_t = fill!(Source.I_err_t, 0) #

    Source.Vd_abc_new = fill!(Source.Vd_abc_new, 0)

    Source.V_dq0 = fill!(Source.V_dq0, 0)
    Source.V_ref_dq0 = fill!(Source.V_ref_dq0, 0)

    Source.V_ref = fill!(Source.V_ref, 0)

    Source.V_err = fill!(Source.V_err, 0) # 
    Source.V_err_t = fill!(Source.V_err_t, 0) #

    Source.I_lim = fill!(Source.I_lim, 0)

    Source.ω_droop = fill!(Source.ω_droop, Source.fsys*2π)

    Source.θ_droop = fill!(Source.θ_droop, 0)

    Source.α_sync = fill!(Source.α_sync, 0)
    Source.ω_sync = fill!(Source.ω_sync, Source.fsys*2π)
    Source.θ_sync = fill!(Source.θ_sync, 0)
    Source.Δω_sync = fill!(Source.Δω_sync, 0)
    Source.eq = fill!(Source.eq, 0)
    Source.Mfif = fill!(Source.Mfif, 0)

end

function Classical_Control(Animo, env, name = nothing)
    
    Source_Interface(env, Animo, name)
    Source = Animo.Source

    for s in 1:Source.num_sources

        if Source.Source_Modes[s] == "Swing Mode"

            Swing_Mode(Source, s)
        elseif Source.Source_Modes[s] == "Voltage Control Mode"

            Voltage_Control_Mode(Source, s)
        elseif Source.Source_Modes[s] == "PQ Control Mode"

            PQ_Control_Mode(Source, s, Source.pq0_set[s, :])
        elseif Source.Source_Modes[s] == "Droop Control Mode"

            Droop_Control_Mode(Source, s)
        elseif Source.Source_Modes[s] == "Synchronverter Mode"

            Synchronverter_Mode(Source, s, pq0_ref = Source.pq0_set[s, :])
        elseif Source.Source_Modes[s] == "Self-Synchronverter Mode"

            Synchronverter_Mode(Source, s, pq0_ref = Source.pq0_set[s, :])
        end

    end

    Action = Env_Interface(env, Animo, name)
    Measurements(Source)

    return Action
end

function Source_Interface(env, Animo, name = nothing)

    Source = Animo.Source

    i = env.steps + 1
    Source.steps = i
    ω = 2*π*Source.fsys
    Source.θsys = (Source.θsys + Source.ts*ω)%(2*π)

    if !isnothing(name)
        state = RLBase.state(env, name)
    end

    for num_source in 1:Source.num_sources

        if isnothing(name)

            #-- Legacy Code

            Source.V_filt_poc[num_source, :, i] = env.x[Source.V_poc_loc[: , num_source]]
            Source.I_filt_poc[num_source, :, i] = env.x[Source.I_poc_loc[: , num_source]]
            Source.I_filt_inv[num_source, :, i] = env.x[Source.I_inv_loc[: , num_source]]
            
            Source.p_q_filt[num_source, :, i] =  p_q_theory(Source.V_filt_poc[num_source, :, i], Source.I_filt_poc[num_source, :, i]) 
        else
            env_ns = 2
            env_ns = num_source
            Source.V_filt_poc[num_source, :, i] = state[findall(contains(string(env_ns)*"_u_C_cables"), Animo.state_ids)]
            if isnothing(findfirst(contains("L2"), Animo.state_ids))
                Source.I_filt_poc[num_source, :, i] = state[findall(contains(string(env_ns)*"_i_L1"), Animo.state_ids)]
            else
                Source.I_filt_poc[num_source, :, i] = state[findall(contains(string(env_ns)*"_i_L2"), Animo.state_ids)]
            end
            Source.I_filt_inv[num_source, :, i] = state[findall(contains(string(env_ns)*"_i_L1"), Animo.state_ids)]
            
            Source.p_q_filt[num_source, :, i] =  p_q_theory(Source.V_filt_poc[num_source, :, i], Source.I_filt_poc[num_source, :, i])
        end

    end

    Animo.Source = Source

    return nothing
end

function Env_Interface(env, Animo, name = nothing)

    Source = Animo.Source

    i = Source.steps

    if isnothing(name)

        # -- Legacy Code
        
        _, B, _, _ = get_sys(env.nc)
        Action = zeros(length(B[1,:]))

        for ph in 1:3
            
            for s in 1:Source.num_sources

                Action[s + Source.num_sources*(ph - 1)] = Source.Vd_abc_new[s, ph, i]
            end
        end
    else
        #Action = Source.Vd_abc_new[:, :, i][:]
        #TODO: when the new action_ids format is released: switch
        letterdict = Dict("a" => 1, "b" => 2, "c" => 3)
        source_indices = unique([parse(Int64, SubString(split(x, "_")[1], 7)) for x in Animo.state_ids], dims=1)
        Action = [Source.Vd_abc_new[findfirst(y -> y == parse(Int64, SubString(split(x, "_")[1], 7)), source_indices),
                                    letterdict[split(x, "_")[3]],
                                    i] for x in Animo.action_ids]
    end

    return Action
end

function Collect_IDs(env, Source::Classical_Controls)

    A, _, _, _ = get_sys(env.nc)
    num_fltr_LCL = env.nc.num_fltr_LCL
    num_fltr_LC = env.nc.num_fltr_LC
    num_fltr_L = env.nc.num_fltr_L

    ns = length(A[1,:])/3 # get num of inputs
    np = 3 # number of phases

    x = 0
    for s in 1:Source.num_sources
        if s <= num_fltr_LCL

            x += 1
            Source.I_inv_loc[1:np, s] = Int.([(x + i*ns) for i in 0:(np-1)])
            x += 2
            Source.I_poc_loc[1:np, s] = Int.([(x + i*ns) for i in 0:(np-1)])
            x += 1
            Source.V_poc_loc[1:np, s] = Int.([(x + i*ns) for i in 0:(np-1)])
        
        elseif s <= num_fltr_LCL + num_fltr_LC

            x += 1
            Source.I_inv_loc[1:np, s] = Int.([(x + i*ns) for i in 0:(np-1)])
            Source.I_poc_loc[1:np, s] = Int.([(x + i*ns) for i in 0:(np-1)])
            x += 2
            Source.V_poc_loc[1:np, s] = Int.([(x + i*ns) for i in 0:(np-1)])
        
        elseif s <= num_fltr_LCL + num_fltr_LC + num_fltr_L

            x += 1
            Source.I_inv_loc[1:np, s] = Int.([(x + i*ns) for i in 0:(np-1)])
            Source.I_poc_loc[1:np, s] = Int.([(x + i*ns) for i in 0:(np-1)])
            x += 1
            Source.V_poc_loc[1:np, s] = Int.([(x + i*ns) for i in 0:(np-1)])

        end
    end

    return nothing
end

function Ramp(final, μ, i; t_end = 0.02, ramp = 0)

    if μ*i < t_end && ramp == 1

        x_out = final.*(μ*i/t_end)
    else
        x_out = final
    end

    return x_out
end

function D_Ramp(D, μ, i; t_end = 0.02, ramp = 0)

    if μ*i < t_end && ramp == 1

        Dout = D.*(μ*i/t_end)

    else
        Dout = D
    end

    return Dout
end

function Swing_Mode(Source::Classical_Controls, num_source; ramp = 0, t_end = 0.08)

    i = Source.steps
    
    δ = Source.V_δ_set[num_source, 1]
    pu = Source.V_pu_set[num_source, 1]

    θt = Source.θsys
    θph = [θt + δ; θt + δ - 120π/180; θt + δ + 120π/180]
    Vrms = Ramp(pu*Source.Vrms[num_source], Source.ts, i; t_end = t_end, ramp = ramp)
    Source.V_ref[num_source, :, i] = sqrt(2)*(Vrms)*cos.(θph)

    Source.Vd_abc_new[num_source, :, i] = Source.V_ref[num_source, :, i]/Source.Vdc[num_source]
    Phase_Locked_Loop_3ph(Source, num_source)

    return nothing
end

function Voltage_Control_Mode(Source::Classical_Controls, num_source; ramp = 0, t_end = 0.04)

    i = Source.steps

    δ = Source.V_δ_set[num_source, 1]
    pu = Source.V_pu_set[num_source, 1]

    θt = Source.θsys# + π/2
    θph = [θt + δ; θt + δ - 120π/180; θt + δ + 120π/180]
    Vrms = Ramp(pu*Source.Vrms[num_source], Source.ts, i; t_end = t_end, ramp = ramp)
    Source.V_ref[num_source, :, i] = sqrt(2)*Vrms*cos.(θph)
    
    Voltage_Controller(Source, num_source, θt)
    Current_Controller(Source, num_source, θt)

    Phase_Locked_Loop_3ph(Source, num_source)

    return nothing
end

function Droop_Control_Mode(Source::Classical_Controls, num_source; ramp = 0, t_end = 0.04)

    i = Source.steps

    Vrms = Ramp(Source.Vrms[num_source], Source.ts, i; t_end = t_end, ramp = ramp)
    Dout = D_Ramp([2π*Source.Δfmax; Source.ΔEmax], Source.ts, i, t_end = t_end, ramp = ramp)

    Source.D[num_source, 1] = Dout[1]/Source.P[num_source]
    Source.D[num_source, 2] = Dout[2]/Source.Q[num_source]

    Droop_Control(Source, num_source, Vrms = Vrms)
    θt = Source.θ_droop[num_source, 1, i]

    Phase_Locked_Loop_3ph(Source, num_source)
    Voltage_Controller(Source, num_source, θt)
    Current_Controller(Source, num_source, θt)

    return nothing
end

function PQ_Control_Mode(Source::Classical_Controls, num_source, pq0)

    i = Source.steps

    Phase_Locked_Loop_3ph(Source, num_source)
    #Phase_Locked_Loop_1ph(Source, num_source, ph = 1)
    #Phase_Locked_Loop_1ph(Source, num_source, ph = 2)
    #Phase_Locked_Loop_1ph(Source, num_source, ph = 3)
    θt = Source.θpll[num_source, 1 , i]

    if i*Source.ts > 2/Source.fsys
        PQ_Control(pq0_ref = pq0, Source, num_source, θt)
    else

        PQ_Control(pq0_ref = [0; 0; 0], Source, num_source, θt)
    end

    Current_Controller(Source, num_source, θt)

    return nothing
end

function Synchronverter_Mode(Source::Classical_Controls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]], ramp = 0, t_end = 0.04)

    i = Source.steps

    Vrms = Ramp(Source.Vrms[num_source], Source.ts, i; t_end = t_end, ramp = ramp)
    Dout = D_Ramp([2π*Source.Δfmax; Source.ΔEmax], Source.ts, i, t_end = t_end, ramp = ramp)

    Source.D[num_source, 1] = Dout[1]/Source.P[num_source]
    Source.D[num_source, 2] = Dout[2]/Source.Q[num_source]

    # Synchronverter parameters
    Dp_sync = 1/(Source.D[num_source, 1]*(2*π)*Source.fsys)
    Source.J_sync[num_source] = Source.τf*Dp_sync

    Dq_sync = sqrt(2)/(Source.D[num_source, 2])
    Source.K_sync[num_source] = Source.τv*Source.fsys*2π*Dq_sync

    if i*Source.ts > 0/Source.fsys
        Synchronverter_Control(Source, num_source; pq0_ref = pq0_ref, Vrms = Vrms)
        θ_S = Source.θ_sync[num_source, i]
        Voltage_Controller(Source, num_source, θ_S)
        Current_Controller(Source, num_source, θ_S)
    end

    Phase_Locked_Loop_3ph(Source, num_source)

    return nothing
end

function Self_Synchronverter_Mode(Source::Classical_Controls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]], SQ = 1, SP = 1)

    i = Source.steps

    if i*Source.ts > 3/Source.fsys
        Self_Synchronverter_Control(Source, num_source, pq0_ref = pq0_ref, SQ = SQ, SP = SP)
        θ_S = Source.θ_sync[num_source, i]
        Voltage_Controller(Kp = 0.01, Ki = 20, Source, num_source, θ_S)
        Current_Controller(Kp = 0.5, Ki = 25, Source, num_source, θ_S)
    end

    return nothing
end

function Phase_Locked_Loop_3ph(Source::Classical_Controls, num_source; ωn = 70, ξ = 0.7)

    #= A robost 3 phase phase locked loop

        The synchronizing unit provides the frequency and the phase of the fundamental
        component of the grid voltage as the references. The negative impact of a
        synchronizing unit are well known. Moreover, because PLLs are inherently
        nonlinear and so are the inverter controller and the power system, it is
        difficult and time-consuming to tune the PLL parameters to achieve
        satisfactory performace. A slow synchronizing unit could directly affect control
        performance and degrade system stability but a complex PLL is computationally
        expensive, which adds significant burden to the controller. Hence, the PLL
        needs to be done quickly and accurately in order to maintain synchronism,
        which makes the design of the controller and the PLL very challenging because
        the PLL is often not fast enough with acceptable accuracy and it also takes time
        for the power and voltage controllers to track the references provided by the
        PLL unit as well.

        The PLL circuit tracks continuously the fundamental frequency of the measured
        system voltages. The appropriate desing of the PLL should allow proper operation
        under distorted and unbalanced voltage waveforms. An interesting desing of
        a PLL circuit that is almost insensitive to imbalances and distortions in the
        voltage waveforms is implemented below. Theis synchronizing circuit determines
        automatically the system frequency and the phase angle of the fundamental
        positive sequence component of a three-phase generic input signal. This circuit
        has proven to be very effecive even under highly distorted system voltages.
        It has been used successfully in several power electronic devices, like
        active filters and FACTS devices. The algorithm is based on a fictitious
        instantaneous active power expression. Note that the fictitious currents in
        a three-phase system sum to zero. In fact the p_3Φ is not related to any
        instantaneous active power of the power system, although it could be
        considered as a variable in the PLL circuit with a dimension of power. This
        is why is ti called fictitious power. As no current is measured fro the power
        circuit, one may find it difficult to understand how the algorithm works. The
        fictitious current feedback signals ia(wt) = sin(wt) and ic(wt) = sin(wt + 2pi/3)
        are built up by the PLL circuit just calculating the time integral of the output
        ω of the PI controller. The PLL can reach a stable point of operation only if
        the input p_3Φ of the PI controller has, in steady state, a zero average
        value. Moreover, the control circuit should minimize oscillation in p_3Φ at
        low frequencies. The oscillating portion of the fictitious power is not well
        attenuated by the PI controller and may bring instability to the PLL control
        circuit. These constraints are only found if ω equals the system frequency,
        and the current ia(wt) becomes orthogonal to the fundamental positive-sequence
        component of the measured three-phase voltages va, vb, and vc.

        However if the point where ia(wt) lags va by 90 degrees is reached, this is still
        an unstable point of operation. At this point, an eventual disturbance that
        slightly increases the system frequency will increase the frequency of vab and vcb
        and will make the positive sequence voltage phasor rotate faster that the current
        phasor corresponding to the generated feedback signal ia(wt). Hence, the displacement
        angle between va and ia, becomes greater than 90 degrees. This leads to a negative
        average input fictitious power and, consequently, decreases the signal ω,
        making the phase angle between va and ia even greater. This positive feedback
        characterizes an unstable point of operation. Thus the PLL has only one stable
        point of operation, that is ia leading by 90 degrees the phase voltage va.
        Now if the same disturbance is verified, the displacement angle between the
        voltage and current phasors will be reduced so that the average power will
        be positive. This will increase the output signal ω and force the current
        phasor to rotate faster, maintaining the orthogonality (leading currents)
        between the positive sequence voltage and fictitious current signals. To be
        clear, ia(wt) = sin(wt) leads by 90 degrees the fundamental positve
        sequence component of the measured system voltages. Thus, a fictitious
        auxiliary current ia_2 = sin(wt - 90deg) will be in phase with V+1.

        err_new = 1*((v_abc[1] - v_abc[2])*cos(-θ)
        + (v_abc[3] - v_abc[2])*cos(-θ - 2π/3)) # this is magic

        v_αβγ = Clarke_Transform(v_abc)
        if norm(v_αβγ) != 0
            v_αβγ = v_αβγ./norm(v_αβγ)
        end
        err_new = v_αβγ[2]*cos(θ) - v_αβγ[1]*sin(θ) # this also works

        Both of these use the same ideas, but the 2nd was implemented below with it's 
        corresponding PI tuning.

        Well-designed PLL systems should meet the following criteria:
        1. ξ ~ 0.7 for optimum transient response (ITAE sense)
        2. narrow bandwidth (low ωn) for improved noise rejection,
            in order to produce a pureley sinusoidal output signal even in 
            the presence of input harmonics

        The PLL "lock range" is defined as the maximum initial frequency deviation
        between the reference input and Voltage Controlled Oscillator (VCO) output,
        which will still cause the PLL to get locked in a single beat. The "lock 
        range" can be shown to be approximately equal to the natural frequency ωn.

        Thus, a narrow-bandwidth PLL may fail in getting locked to the fundamental frequency
        component of a given input signal during the start-up transient if following 
        conditions are met simultaneously.
        1. The input signal contains higher-order harmonics or sub-harmonic components.
        2. One of these harmonics components has a frequency that is close to the initial 
        PI output (or equivalently, to the center frequency)
        3. The difference between the initial PI output and the target fundamental frequency
        is larger than the lock range.
        It is, however, very difficult to predict the actual behaviour of the PLL under 
        the above conditions because it depends on the relative amplitude of the harmonic 
        components. For instance, subharmonic oscillations at the reference input can 
        cause the PLL to lock at the lower subharmonic frequency even if the relative 
        magnitude is very low.

        The settling time for the controller is approximately ξ*ωn
        
    =#

    i = Source.steps

    Ki = ωn^2 # tuning
    Kp = ξ*2*sqrt(Ki) # tuning

    v_abc = Source.V_filt_poc[num_source, :, i]

    range, cnt_end = Integrator_Prep(i)
    f = Source.fpll[num_source, 1, range]
    θ = Source.θpll[num_source, 1, cnt_end]

    err_t = Source.pll_err_t[num_source, 1]
    err = Source.pll_err[num_source, 1, :]

    v_αβγ = Clarke_Transform(v_abc)
    if norm(v_αβγ) != 0
        v_αβγ = sqrt(3)*v_αβγ./norm(v_αβγ)
    end
    err_new = v_αβγ[2]*cos(θ) - v_αβγ[1]*sin(θ)

    f_new, err_t_new, err_int =
    PI_Controller(err_new, err, err_t, Kp, Ki, Source.ts, bias = Source.fsys)

    θ = Third_Order_Integrator(θ, Source.ts, 2π*[f; f_new])

    Source.fpll[num_source, :, i] = [f_new; f_new; f_new]
    Source.θpll[num_source, :, i] = [θ; θ - 120π/180; θ + 120π/180].%(2π)
    Source.pll_err_t[num_source, :] = [err_t_new; 0; 0]
    Source.pll_err[num_source, :, :] = [transpose(err_int); transpose(err_int); transpose(err_int)]

    return nothing
end

function Phase_Locked_Loop_1ph(Source::Classical_Controls, num_source; Kp = 0.001, Ki = 1, ph = 1, k_sogi = 0.8)

    i = Source.steps

    range_1, cnt_end_1 = Integrator_Prep(i)
    range_2, _ = Integrator_Prep(i-1)

    v_ph = Source.V_filt_poc[num_source, ph, i]
    v_ph_r = Source.V_filt_poc[num_source, ph, range_1]
    f_1 = Source.fpll[num_source, ph, range_1]
    f_2 = Source.fpll[num_source, ph, range_2]
    θ = Source.θpll[num_source, ph, cnt_end_1]
    ω_1 = 2π*f_1
    ω_2 = 2π*f_2

    vd_1 = Source.vd[num_source, ph, range_1]
    qvd_1 = Source.qvd[num_source, ph, range_1]
    vd_2 = Source.vd[num_source, ph, range_2]
    qvd_2 = Source.qvd[num_source, ph, range_2]

    err_t = Source.pll_err_t[num_source, ph]
    err = Source.pll_err[num_source, ph, :]

    #= k_sogi
    The level of filtering can be set from gain k_sogi as follows: - if k decreases
    the bandpass of the filter becomes narrower resulting in heavy filtering, but
    in the same time the dynamic response of the system will become slower
    =#

    #--- Second Order Generalised Integrator
    l = length(v_ph_r) - length(vd_2) + 1

    vd_sogi_old = (k_sogi*(v_ph_r[l:end] - vd_2) - qvd_2).*ω_2
    vd_sogi_new = (k_sogi*(v_ph - vd_1[end]) - qvd_1[end])*ω_1[end]
    vd_sogi = [vd_sogi_old; vd_sogi_new]
    vd_new = Third_Order_Integrator(vd_1[end], Source.ts, vd_sogi)

    qvd_int_old = vd_2.*ω_2
    qvd_int_new = vd_1[end]*ω_1[end]
    qvd_int = [qvd_int_old; qvd_int_new]
    qvd_new = Third_Order_Integrator(qvd_1[end], Source.ts, qvd_int)
    #----

    α_β_0 = [vd_new; qvd_new; 0]
    d_q_0 = Park_Transform(α_β_0, θ - π)

    vq = d_q_0[1]

    err_new = 0 - vq

    ω_new, err_t_new, err_int =
    PI_Controller(err_new, err, err_t, Kp, Ki, Source.ts, bias = 2*π*Source.fsys)

    Source.θpll[num_source, ph, i] =
    (Third_Order_Integrator(θ, Source.ts, [ω_1; ω_new]))%(2*π)

    Source.pll_err[num_source, ph, :] = err_int
    Source.vd[num_source, ph, i] = vd_new
    Source.qvd[num_source, ph, i] = qvd_new

    Source.pll_err_t[num_source, ph] = err_t_new[end]
    Source.fpll[num_source, ph, i] = ω_new[end]/(2*π)

    return nothing
end

function PQ_Control(Source::Classical_Controls, num_source, θ; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0], Kb = 1)

    i = Source.steps

    Kp = Source.V_kp[num_source]
    Ki = Source.V_ki[num_source]

    V_err = Source.V_err[num_source, :, :]
    V_err_t = Source.V_err_t[num_source, :]

    if norm(pq0_ref) > Source.S[num_source]
        pq0_ref = pq0_ref.*(Source.S[num_source]/norm(pq0_ref))
    end

    V_αβγ = Clarke_Transform(Source.V_filt_poc[num_source, :, i])
    I_αβγ = Clarke_Transform(Source.I_filt_poc[num_source, :, i])

    #-------------------------------------------------------------

    I_αβγ_ref = Inv_p_q_v(V_αβγ, pq0_ref)

    I_dq0_ref = Park_Transform(I_αβγ_ref, θ)
    I_dq0 = Park_Transform(I_αβγ, θ)

    #-------------------------------------------------------------

    V_αβγ_ref = Inv_p_q_i(I_αβγ, pq0_ref)

    V_dq0_ref = Park_Transform(V_αβγ_ref, θ)
    V_dq0 = Park_Transform(V_αβγ, θ)

    if sqrt(2/3)*norm(V_dq0_ref) > Source.Vdc[num_source]/2
        V_dq0_ref = V_dq0_ref.*((Source.Vdc[num_source]/2)/(sqrt(2/3)*norm(V_dq0_ref) ))
    end

    #-------------------------------------------------------------

    Source.V_ref[num_source, :, i] = Inv_DQ0_transform(V_dq0_ref, θ)

    Source.V_ref_dq0[num_source, :, i] = V_dq0_ref
    Source.V_dq0[num_source, :, i] = V_dq0

    if i > 1
        # Including Anti-windup - Back-calculation
        I_err_new = I_dq0_ref .- I_dq0 .+ Kb*(Source.I_ref_dq0[num_source, :, i - 1] .- Source.I_lim[num_source,:])
    else
        I_err_new = I_dq0_ref .- I_dq0
    end

    Source.I_lim[num_source,:], Source.V_err_t[num_source, :], 
    Source.V_err[num_source, :, :] =
    PI_Controller(I_err_new, V_err, V_err_t, Kp, Ki, Source.ts)

    Ip_ref = sqrt(2/3)*norm(Source.I_lim[num_source, :]) # peak set point

    if Ip_ref > Source.i_max[num_source]
        Source.I_ref_dq0[num_source, :, i] = Source.I_lim[num_source,:]*Source.i_max[num_source]/Ip_ref
    else
        Source.I_ref_dq0[num_source, :, i] = Source.I_lim[num_source,:]
    end

    return nothing
end

#= function PQ_Control(Source::Classical_Controls, num_source, θ; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0], Kb = 1)

    i = Source.steps

    Kp = Source.V_kp[num_source]
    Ki = Source.V_ki[num_source]

    V_err = Source.V_err[num_source, :, :]
    V_err_t = Source.V_err_t[num_source, :]

    if norm(pq0_ref) > Source.S[num_source]
        pq0_ref = pq0_ref.*(Source.S[num_source]/norm(pq0_ref))
    end

    V_αβγ = Clarke_Transform(Source.V_filt_poc[num_source, :, end])
    I_αβγ = Clarke_Transform(Source.I_filt_poc[num_source, :, end])

    #-------------------------------------------------------------

    I_αβγ_ref = Inv_p_q_v(V_αβγ, pq0_ref)

    I_dq0_ref = Park_Transform(I_αβγ_ref, θ)
    I_dq0 = Park_Transform(I_αβγ, θ)

    #-------------------------------------------------------------

    V_αβγ_ref = Inv_p_q_i(I_αβγ, pq0_ref)

    V_dq0_ref = Park_Transform(V_αβγ_ref, θ)
    V_dq0 = Park_Transform(V_αβγ, θ)

    if sqrt(2/3)*norm(V_dq0_ref) > Source.Vdc[num_source]/2
        V_dq0_ref = V_dq0_ref.*((Source.Vdc[num_source]/2)/(sqrt(2/3)*norm(V_dq0_ref) ))
    end

    #-------------------------------------------------------------

    Source.V_ref[num_source, :] = Inv_DQ0_transform(V_dq0_ref, θ)

    Source.V_ref_dq0[num_source, :] = V_dq0_ref
    Source.V_dq0[num_source, :] = V_dq0

    if i > 1
        # Including Anti-windup - Back-calculation
        I_err_new = I_dq0_ref .- I_dq0 .+ Kb*(Source.I_ref_dq0[num_source, :] .- Source.I_lim[num_source,:])
    else
        I_err_new = I_dq0_ref .- I_dq0
    end

    Source.I_lim[num_source,:], Source.V_err_t[num_source, :], 
    Source.V_err[num_source, :, :] =
    PI_Controller(I_err_new, V_err, V_err_t, Kp, Ki, Source.ts)

    Ip_ref = sqrt(2/3)*norm(Source.I_lim[num_source, :]) # peak set point

    if Ip_ref > Source.i_max[num_source]
        Source.I_ref_dq0[num_source, :] = Source.I_lim[num_source,:]*Source.i_max[num_source]/Ip_ref
    else
        Source.I_ref_dq0[num_source, :] = Source.I_lim[num_source, :]
    end

    return nothing
end =#

#= function PQ_Control(Source::Classical_Controls, num_source, θ; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0], Kb = 1)

    i = Source.steps

    Kp = Source.I_kp[num_source]
    Ki = Source.I_ki[num_source]

    I_err = Source.I_err[num_source, :, :]
    I_err_t = Source.I_err_t[num_source, :]

    if norm(pq0_ref) > Source.S[num_source]
        pq0_ref = pq0_ref.*(Source.S[num_source]/norm(pq0_ref))
    end

    V_αβγ = Clarke_Transform(Source.V_filt_poc[num_source, :, i])
    I_αβγ = Clarke_Transform(Source.I_filt_poc[num_source, :, i])

    #-------------------------------------------------------------

    I_αβγ_ref = Inv_p_q_v(V_αβγ, pq0_ref)

    I_dq0_ref = Park_Transform(I_αβγ_ref, θ)
    I_dq0 = Park_Transform(I_αβγ, θ)

    #I_dq0_ref = [50; -100; 0]

    I_err_new = I_dq0_ref .- I_dq0

    s_dq0_avg, Source.I_err_t[num_source, :], Source.I_err[num_source, :, :] =
    PI_Controller(I_err_new, I_err, I_err_t, Kp, Ki, Source.ts)

    Source.Vd_abc_new[num_source, :, i + Source.delay] = 0.5*Inv_DQ0_transform(s_dq0_avg, θ)
    
    Source.I_ref_dq0[num_source, :,i] = I_dq0_ref
    Source.I_dq0[num_source, :, i] = I_dq0
    Source.I_ref[num_source, :, i] = Inv_DQ0_transform(Source.I_ref_dq0[num_source, :,i], θ)

    return nothing
end =#

function Droop_Control(Source::Classical_Controls, num_source; Vrms = Source.Vrms[num_source])

    #=
    The droop control method has been referred to as the independent, autonomous,
    and wireless control due to elimination of intercommunication links between
    the converters.

    Droop control works best when switching ripples and high frequency harmonics
    are neglegtc, in which case the VSC can be modeled as an AC source.

    The droop coefficients, Dp and Dq, can be adjusted either heuristically or by
    tuning algorithms (e.g., particle swarm optimisation). In the former approach,
    Dp and Dq are determined based on the converter power rating and the maximum
    allowable voltage and frequency deviations.
    =#

    i = Source.steps

    range, range_end = Integrator_Prep(i)

    μ = Source.ts
    ω = Source.ω_droop[num_source, 1, range]
    θ = Source.θ_droop[num_source, 1, range_end]
    p_q = Source.p_q_filt[num_source, :, i]

    ω_new = Source.fsys*2*π - p_q[1]*Source.D[num_source, 1]
    Source.ω_droop[num_source, :, i] = [ω_new; ω_new; ω_new]

    Source.θ_droop[num_source, 1, i] = Third_Order_Integrator(θ, μ, [ω; ω_new])%(2*π)
    Source.θ_droop[num_source, 2, i] = (Source.θ_droop[num_source, 1, i] - 120*π/180)%(2*π)
    Source.θ_droop[num_source, 3, i] = (Source.θ_droop[num_source, 1, i] + 120*π/180)%(2*π)

    E_new = Vrms - p_q[2]*Source.D[num_source, 2]

    Source.V_ref[num_source, 1, i] = sqrt(2)*E_new*cos(Source.θ_droop[num_source, 1, i])
    Source.V_ref[num_source, 2, i] = sqrt(2)*E_new*cos(Source.θ_droop[num_source, 2, i])
    Source.V_ref[num_source, 3, i] = sqrt(2)*E_new*cos(Source.θ_droop[num_source, 3, i])

    return nothing
end

function PI_Controller(Error_new, Error_Hist, Error_t, Kp, Ki, μ; bias = zeros(length(Error_new)))

    d = length(Error_new)

    if d == 1
        Error_Hist = transpose(Error_Hist)
    end

    Err_t_new = Array{Float64, 1}(undef, d)
    Err_t_new = fill!(Err_t_new, 0)

    Err_d = [Error_Hist Error_new]

    Err_int = Err_d[:, 2:end]

    for j in 1:d
        Err_t_new[j] = Third_Order_Integrator(Error_t[j], μ, Err_int[j,:]) # integration
        #Err_t_new[j] = Error_t[j] + μ*Err_int[j, end] # integration
    end

    Action = bias + Kp.*Error_new .+ Ki.*Err_t_new

    return Action, Err_t_new, Err_int
end

function Current_Controller(Source::Classical_Controls, num_source, θ)

    #=
        When a grid-connected inverter is controlled as a current supply, the output
        voltage is mainted by the grid and the inverter only regulates current exchanged
        with the grid. Some simplified synchronizing methods can be adopted and no
        extra effort is needed to design the synchronizing unit. Because of the
        simplified control structure and the reduced demand on the synchronization
        unit, it is well known that it is much easier to control a grid-connected inverter
        as a current supply than to control it as a voltage supply. However, when an
        inverter is controlled as a current supply, it causes undesirable problems.
        For example, the controller needs to be changed when the inverter is disconnected
        from the grid to operate in the standalone mode or when the grid is weak because,
        it does not have the capability of regulating the voltage. A current-controlled
        inverter may also continue injecting currents into the grid when there is
        a fault on the grid, which might cause excessively high voltage. Moreover,
        a current-controlled inverter is difficult to take part in the regulation of
        the grid frequency and voltage, which is a must when the penetration of
        renewable energy exceeds a certain level.
    =#

    # --- For Debugging
    #Source.I_ref_dq0[num_source, :,i] = DQ0_transform(Source.I_ref[num_source, :, i], θ)
    # --- =#
    i = Source.steps

    Kp = Source.I_kp[num_source]
    Ki = Source.I_ki[num_source]

    Source.I_ref[num_source, :, i] = Inv_DQ0_transform(Source.I_ref_dq0[num_source, :,i], θ)
    Source.I_dq0[num_source, :, i] = DQ0_transform(Source.I_filt_inv[num_source, :, i], θ)

    I_dq0 = Source.I_dq0[num_source, :, i]
    I_ref_dq0 = Source.I_ref_dq0[num_source, :,i]
    I_err = Source.I_err[num_source, :, :]
    I_err_t = Source.I_err_t[num_source, :]

    I_err_new = I_ref_dq0 .- I_dq0

    s_dq0_avg, Source.I_err_t[num_source, :], Source.I_err[num_source, :, :] =
    PI_Controller(I_err_new, I_err, I_err_t, Kp, Ki, Source.ts)

    Source.Vd_abc_new[num_source, :, i + Source.delay] = 0.5*Inv_DQ0_transform(s_dq0_avg, θ)
    #Source.Vd_abc_new[num_source, :, i + Source.delay] = 0.5*Source.Vdc[num_source].*Inv_DQ0_transform(s_dq0_avg, θ)
    #= Theory:
        The switching functions s_abc(t) is generated by comparing the normalized
        voltage value with a triangular modulation carrier. The output of the
        comparator can be directly referred to as the switching function. Through
        geometric interpretation of this procedure it becomes clear that the time
        average of the switching function corresponds to the reference signal, as
        long as the phasor components of the reference signal can be assumed to be
        slowly varying. That is, the above equation holds when averaging over one
        pulse period.
    =#

    return nothing
end

function Voltage_Controller(Source::Classical_Controls, num_source, θ; Kb = 1)

    i = Source.steps

    Kp = Source.V_kp[num_source]
    Ki = Source.V_ki[num_source]

    Source.V_ref_dq0[num_source, :,i] = DQ0_transform(Source.V_ref[num_source, :, i], θ)
    Source.V_dq0[num_source, :, i] = DQ0_transform(Source.V_filt_poc[num_source, :, i], θ)
    V_dq0 = Source.V_dq0[num_source, :, i]
    V_ref_dq0 = Source.V_ref_dq0[num_source, :, i]

    if sqrt(2/3)*norm(Source.V_ref_dq0[num_source, :,i]) > Source.v_max[num_source]/2 #Source.Vdc[num_source]
        Source.V_ref_dq0[num_source, :,i] = Source.V_ref_dq0[num_source, :,i].*((Source.Vdc[num_source]/2)/(sqrt(2/3)*norm(V_ref_dq0)))
        V_ref_dq0 = Source.V_ref_dq0[num_source, :, i] # v_max Vdc
    end

    V_err = Source.V_err[num_source, :, :]
    V_err_t = Source.V_err_t[num_source, :]

    if i > 1
        # Including Anti-windup - Back-calculation
        V_err_new = V_ref_dq0 .- V_dq0 .+ Kb*(Source.I_ref_dq0[num_source, :, i - 1] .- Source.I_lim[num_source,:])
    else
        V_err_new = V_ref_dq0 .- V_dq0
    end

    Source.I_lim[num_source,:], Source.V_err_t[num_source, :], Source.V_err[num_source, :, :] =
    PI_Controller(V_err_new, V_err, V_err_t, Kp, Ki, Source.ts)

    # ---- Limiting Output (Saturation)
    Ip_ref = sqrt(2/3)*norm(Source.I_lim[num_source,:]) # peak set point

    if Ip_ref > Source.i_max[num_source]
        Source.I_ref_dq0[num_source, :, i] = Source.I_lim[num_source,:]*Source.i_max[num_source]/Ip_ref
    else
        Source.I_ref_dq0[num_source, :, i] = Source.I_lim[num_source,:]
    end

    return nothing
end

function Synchronverter_Control(Source::Classical_Controls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0], Vrms = Source.Vrms[num_source])

    #=
        A synchronverter is an inverter that mimics synchronous generators, which
        offers a mechanism for power systems to control grid-connected renewable energy
        and facilitates smart grid integration. Similar to other grid-connected inverters,
        it typically needs a dedicated synchronizing unit, e.g. a phase-locked loop,
        to provide a the phase, frequency, and amplitude of the grid voltages as
        references. However, in this implementation the synchronizing unit has been
        removed. It can automatically synchronize itself with the grid before connection
        and tracks the grid frequency after connection. All the function of the original
        synchronverter, such as frequency and voltage regulation, real power, and
        reactive power control, are maintained.

        Synchronverters are grid-friendly inverters that mimic synchronous generatons.
        A synchronverter includes the mathematical model of a synchronous machines and
        behaves in the same way, mathematically, as a synchronous generator to provide
        a voltage supply. Its controller is in principle a power controller with integrated
        capability of voltage and frequency regulation so it is able to achieve real
        power control, reactive power control, frequency regulation,a nd voltage regulation.
        Because of the embedded mathematical model, a utility is able to control a
        synchronverter in the same way as controlling synchronous generators, which
        consideratbly facilitates the grid connection of renewable energy and smart
        grid integration. Since a synchronous machine is inherently able to synchronize
        with the grid, it should be possible to integrate synchronization function
        into the power controller and make a synchronverter to synchronize with the
        grid without a dedicated PLL. This would remove the slow element in the closed
        loop system consisteing of the PLL, the inverter controller and the power
        system, and removes a major nonlinear element that affects the speed and
        accuracy of synchronization.

        A synchronverter will have all the good and bad properties of an SG, which is
        nonlinear system. For example, the undesirable phenomena, such as loss of
        stability due to underexcitation as well as hunting (oscillations around the
        synchronous frequency), could occur in a synchronverter. An advantage is that
        we can choose the parameters, such as intertia, friction coefficient, field
        inductance, and mutual inductances. (The energy that would be lost in the
        virtual mechanical friction is not lost in reality; it is directed back to
        the DC bus). Moreover, we can (and do) choose to have no magnetic saturation
        and no eddy currents. If we want, we can choose parameter values that are
        impossible in a real SG, and we can also vary the parameters while the system
        is operating. Synchronverters can alos be operated as synchronous motors based
        on the same mathematical derivation.

        In this model, the number of pairs of poles for each phase is 1 and hence the
        mechanical speed of the machine is the same as the electrical speed of the
        electromagnetic field.

        Similarly to the control of a synchronous generator, the controller of a
        synchronverter has two channels: one for the real power and the other for
        reactive power. The real power is controlled by a frequency droop loop,
        using the virtual mechanical friction coefficient Dp as the feedback gain.
        The loop regulated the virtual speed, ω, of the synchronous machine and creates
        the phase angle, θ, for the control signal, e.

        The reactive power is controlled by a voltage droop loop, using the voltage
        droop coefficient Dq. This loop regulates the field excitation Mf*if, which
        is proportional to the amplitude of the voltage generated. Hence, the frequency
        control, voltage control, real power control, and reactive power control are
        all integrated in one compact controller with only four parameters (Dp, Dq,
        K, and J)

        For grid-connected applications, a synchronization unit is needed to provide the
        grid information for the synchronverter to synchronize with the grid before
        connection and for the synchronverter to deliver the desired real and reactive
        powers after connection. If P and Q are controlled to be zero, then the
        generated voltage e is the same as the grid voltage vg. This condition is not
        common in the normal operation of an SG, but when it is satisfied, the SG
        can be connected to or disconnected from the grid without causing large
        transient dynamics. The same applies to a synchronverter to synchronize with
        the grid before connection. However, by making some changes to the core of
        the model the synchronverter can be connected to the grid safely and operated
        without the need of a dedicated PLL.

        When the VSG is operated in P-droop mode the real power is not the same as Pset,
        but deviated from the value. This is because preference is given to maintaining
        a minimum frequency deviation of the SVG from the nominal frequency. The same
        holds for Q-droop mode, where the reactive power is not the same as Qset, but
        deviated from the Qset. In this case preference is given to voltage deviation.

        Typical values for the frequency droop are a 100% increase in power for a
        frequency decrease between 3% and 5%.

        τv = 0.002 # s
        τf = 0.002 # s
        Dp = 0.2026
        Dq = 117.88
    =#
    i = Source.steps

    J = Source.J_sync[num_source]
    K = Source.K_sync[num_source]
    
    Dp = 1/(Source.D[num_source, 1]*(2*π)*Source.fsys)
    Dq = sqrt(2)/(Source.D[num_source, 2])

    range, range_end = Integrator_Prep(i)
    eq = Source.eq[num_source, range]
    Mfif = Source.Mfif[num_source, range]

    α = Source.α_sync[num_source, range]
    ω = Source.ω_sync[num_source, range]
    θ = Source.θ_sync[num_source, range_end] # only phase a

    ωsys = Source.fsys*2π # nominal grid frequency
    Vn = sqrt(2)*Vrms # nominal peak POC voltage
    Vg = sqrt(2/3)*norm(DQ0_transform(Source.V_filt_poc[num_source, :, i], 0)) # peak measured voltage

    μ = Source.ts

    i_abc = Source.I_filt_poc[num_source, :, i]

    sin_θ = sin.([θ; θ - 120*π/180; θ + 120*π/180])
    cos_θ = cos.([θ; θ - 120*π/180; θ + 120*π/180])

    #---- Integrate eq_new to find Mfif_new

    Q = ω[end]*Mfif[end]*dot(i_abc, sin_θ) # Reactive Power
    eq_new = (1/K)*(pq0_ref[2] + Dq*(Vn - Vg) - Q)

    Mfif_new = Third_Order_Integrator(Mfif[end], μ, [eq; eq_new])
    Mfif_new = Mfif[end] + μ*eq_new

    Source.Mfif[num_source, i] = Mfif_new
    Source.eq[num_source, i] = eq_new
    #----

    #---- Integrate α_new to find ω_new
    Tm = pq0_ref[1]/ωsys # Virtual machine Torque

    ω_err_new = ω[end] - ωsys
    ΔT_new = Dp*ω_err_new

    Te_new = Mfif_new*dot(i_abc, cos_θ) # New Electrical Torque

    α_new = (1/J)*(Tm - Te_new - ΔT_new) # New Angular Acceleration

    ω_new = Third_Order_Integrator(ω[end], μ, [α; α_new])
    ω_new = ω[end] + μ*α_new

    Source.ω_sync[num_source, i] = ω_new
    Source.α_sync[num_source, i] = α_new
    #----

    #---- Integrate ω_new to find θ_new
    θ_new = Third_Order_Integrator(θ, μ, [ω; ω_new])%(2*π)

    Source.θ_sync[num_source, i] = θ_new
    #----

    #----
    cos_θ_new = cos.([θ_new; θ_new - 120*π/180; θ_new + 120*π/180])
    e = ω_new*Mfif_new*cos_θ_new # three phase generated voltage
    #----

    Source.V_ref[num_source, 1, i] = e[1]
    Source.V_ref[num_source, 2, i] = e[2]
    Source.V_ref[num_source, 3, i] = e[3]

    return nothing
end

#= function Synchronverter_Control(Source::Classical_Controls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0], Vrms = Source.Vrms[num_source])

    #=
        A synchronverter is an inverter that mimics synchronous generators, which
        offers a mechanism for power systems to control grid-connected renewable energy
        and facilitates smart grid integration. Similar to other grid-connected inverters,
        it typically needs a dedicated synchronizing unit, e.g. a phase-locked loop,
        to provide a the phase, frequency, and amplitude of the grid voltages as
        references. However, in this implementation the synchronizing unit has been
        removed. It can automatically synchronize itself with the grid before connection
        and tracks the grid frequency after connection. All the function of the original
        synchronverter, such as frequency and voltage regulation, real power, and
        reactive power control, are maintained.

        Synchronverters are grid-friendly inverters that mimic synchronous generatons.
        A synchronverter includes the mathematical model of a synchronous machines and
        behaves in the same way, mathematically, as a synchronous generator to provide
        a voltage supply. Its controller is in principle a power controller with integrated
        capability of voltage and frequency regulation so it is able to achieve real
        power control, reactive power control, frequency regulation,a nd voltage regulation.
        Because of the embedded mathematical model, a utility is able to control a
        synchronverter in the same way as controlling synchronous generators, which
        consideratbly facilitates the grid connection of renewable energy and smart
        grid integration. Since a synchronous machine is inherently able to synchronize
        with the grid, it should be possible to integrate synchronization function
        into the power controller and make a synchronverter to synchronize with the
        grid without a dedicated PLL. This would remove the slow element in the closed
        loop system consisteing of the PLL, the inverter controller and the power
        system, and removes a major nonlinear element that affects the speed and
        accuracy of synchronization.

        A synchronverter will have all the good and bad properties of an SG, which is
        nonlinear system. For example, the undesirable phenomena, such as loss of
        stability due to underexcitation as well as hunting (oscillations around the
        synchronous frequency), could occur in a synchronverter. An advantage is that
        we can choose the parameters, such as intertia, friction coefficient, field
        inductance, and mutual inductances. (The energy that would be lost in the
        virtual mechanical friction is not lost in reality; it is directed back to
        the DC bus). Moreover, we can (and do) choose to have no magnetic saturation
        and no eddy currents. If we want, we can choose parameter values that are
        impossible in a real SG, and we can also vary the parameters while the system
        is operating. Synchronverters can alos be operated as synchronous motors based
        on the same mathematical derivation.

        In this model, the number of pairs of poles for each phase is 1 and hence the
        mechanical speed of the machine is the same as the electrical speed of the
        electromagnetic field.

        Similarly to the control of a synchronous generator, the controller of a
        synchronverter has two channels: one for the real power and the other for
        reactive power. The real power is controlled by a frequency droop loop,
        using the virtual mechanical friction coefficient Dp as the feedback gain.
        The loop regulated the virtual speed, ω, of the synchronous machine and creates
        the phase angle, θ, for the control signal, e.

        The reactive power is controlled by a voltage droop loop, using the voltage
        droop coefficient Dq. This loop regulates the field excitation Mf*if, which
        is proportional to the amplitude of the voltage generated. Hence, the frequency
        control, voltage control, real power control, and reactive power control are
        all integrated in one compact controller with only four parameters (Dp, Dq,
        K, and J)

        For grid-connected applications, a synchronization unit is needed to provide the
        grid information for the synchronverter to synchronize with the grid before
        connection and for the synchronverter to deliver the desired real and reactive
        powers after connection. If P and Q are controlled to be zero, then the
        generated voltage e is the same as the grid voltage vg. This condition is not
        common in the normal operation of an SG, but when it is satisfied, the SG
        can be connected to or disconnected from the grid without causing large
        transient dynamics. The same applies to a synchronverter to synchronize with
        the grid before connection. However, by making some changes to the core of
        the model the synchronverter can be connected to the grid safely and operated
        without the need of a dedicated PLL.

        When the VSG is operated in P-droop mode the real power is not the same as Pset,
        but deviated from the value. This is because preference is given to maintaining
        a minimum frequency deviation of the SVG from the nominal frequency. The same
        holds for Q-droop mode, where the reactive power is not the same as Qset, but
        deviated from the Qset. In this case preference is given to voltage deviation.

        Typical values for the frequency droop are a 100% increase in power for a
        frequency decrease between 3% and 5%.

        τv = 0.002 # s
        τf = 0.002 # s
        Dp = 0.2026
        Dq = 117.88
    =#

    J = Source.J_sync[num_source]
    K = Source.K_sync[num_source]
    
    Dp = 1/((2π)*(Source.fsys)*Source.D[num_source, 1])
    Dq = 1/(Source.D[num_source, 2])

    eq = Source.eq[num_source, :]
    Mfif = Source.Mfif[num_source]

    α = Source.α_sync[num_source, :]
    ω = Source.ω_sync[num_source, :]
    θ = Source.θ_sync[num_source] # only phase a

    ωsys = Source.fsys*2π # nominal grid frequency
    Vn = sqrt(2)*Vrms # nominal peak POC voltage
    Vg = sqrt(2/3)*norm(DQ0_transform(Source.V_filt_poc[num_source, :, end], 0)) # peak measured voltage

    μ = Source.ts

    i_abc = Source.I_filt_poc[num_source, :, end]

    sin_θ = sin.([θ; θ - 120*π/180; θ + 120*π/180])
    cos_θ = cos.([θ; θ - 120*π/180; θ + 120*π/180])

    #---- Integrate eq_new to find Mfif_new

    #Q = ω[end]*Mfif[end]*dot(i_abc, sin_θ) # Reactive Power
    Q = Source.p_q_filt[num_source, 2]
    eq_new = (1/K)*(pq0_ref[2] + Dq*(Vn - Vg) - Q)

    Mfif_new = Third_Order_Integrator(Mfif, μ, [eq[2:end]; eq_new])

    Source.Mfif[num_source] = Mfif_new
    Source.eq[num_source, :] = [eq[2:end]; eq_new]
    #----

    #---- Integrate α_new to find ω_new
    Tm = pq0_ref[1]/ωsys # Virtual machine Torque

    ω_err_new = ω[end] - ωsys
    ΔT_new = Dp*ω_err_new

    #Te_new = Mfif_new*dot(i_abc, cos_θ) # New Electrical Torque
    Te_new = Source.p_q_filt[num_source, 1]/ω[end] # New Electrical Torque

    α_new = (1/J)*(Tm - Te_new - ΔT_new) # New Angular Acceleration

    ω_new = Third_Order_Integrator(ω[end], μ, [α[2:end]; α_new])

    Source.ω_sync[num_source, :] = [ω[2:end]; ω_new]
    Source.α_sync[num_source, :] = [α[2:end]; α_new]
    #----

    #---- Integrate ω_new to find θ_new
    θ_new = Third_Order_Integrator(θ, μ, [ω[2:end]; ω_new])%(2*π)

    Source.θ_sync[num_source] = θ_new
    #----

    #----
    cos_θ_new = cos.([θ_new; θ_new - 120*π/180; θ_new + 120*π/180])
    e = ω_new*Mfif_new*cos_θ_new # three phase generated voltage
    #----

    Source.V_ref[num_source, 1] = e[1]
    Source.V_ref[num_source, 2] = e[2]
    Source.V_ref[num_source, 3] = e[3]

    if Source.steps*Source.ts >= 3.0 - Source.ts/2 && Source.steps*Source.ts <= 3.0 + Source.ts/2
        println("")
        println("--------------------------------------------------------")

        v_inv = 0.5*Source.Vdc[num_source].*Source.Vd_abc_new[num_source, :]
        p_q_inv = p_q_theory(v_inv, Source.I_filt_inv[num_source, :, end])

        Ip_ref = sqrt(2/3)*norm(DQ0_transform(Source.I_filt_inv[num_source, :, end], 0)) 

        println("")
        println("mode = ", mode)
        println("num_source = ", num_source)

        println("")
        println("Ip_ref = ", Ip_ref)
        println("i_max = ", Source.i_max[num_source])

        println("")
        println("ω_new = ", ω_new)
        println("ω[end] = ", ω[end])
        println("ω_err_new = ", ω_err_new)
        println("eq_new = ", eq_new)
        println("(Vn - Vg) = ", (Vn - Vg))
        println("(Vn) = ", (Vn))
        println("(Vg) = ", (Vg))
        println("(1/K)*(Dq) = ", (1/K)*(Dq*(Vn - Vg)))

        println("")
        println("p_filt = ", Source.p_q_filt[num_source, 1])
        println("Pset = ", pq0_ref[1])

        println("")
        println("q_filt = ", Source.p_q_filt[num_source, 2])
        println("Qset = ", pq0_ref[2])

        println("V_ph[ns,  1, 2] = ", Source.V_ph[num_source, 1, 2]/230)
        println("I_ph[ns,  1, 2] = ", sqrt(2)*Source.I_ph[num_source, 1, 2])

        println("")
        println("--------------------------------------------------------")
    end

    return nothing
end =#

function Self_Synchronverter_Control(Source::Classical_Controls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0], SQ = 1, SP = 1, SC = 0, Kp = 0.0005, Ki = 0.01)


    ##### REMEMBER TO CHANGE COS TO SINE AND VICE VERSA - DQ0 CONSISTENT

    i = Source.steps
    # SQ = 1 for Q-droop mode, 0 for Q-set mode
    # SP = 1 for P-droop mode, 0 for P-set mode
    # SC = 1 for Self-synchronisation mode

    J = Source.J_sync[num_source]
    K = Source.K_sync[num_source]
    Dq = sqrt(2)/(Source.D[num_source, 2])
    Dp = 1/(Source.D[num_source, 1]*(2*π)*Source.fsys)

    range, range_end = Integrator_Prep(i)
    Δω = Source.Δω_sync[num_source, range]
    ΔT_t = Source.ΔT_t[num_source]
    eq = Source.eq[num_source, range]
    Mfif = Source.Mfif[num_source, range]

    α = Source.α_sync[num_source, range]
    ω = Source.ω_sync[num_source, range]
    θ = Source.θ_sync[num_source, range_end] # only phase a

    ωsys = Source.fsys*2π # nominal grid frequency
    Vn = sqrt(2)*Source.Vrms[num_source] # nominal peak POC voltage
    Vg_abc = Source.V_filt_poc[num_source, :, i]
    Vg = sqrt(2/3)*norm(DQ0_transform(Vg_abc, 0)) # peak measured voltage

    μ = Source.ts
    pq0 = Source.p_q_filt[num_source, :, i]
    i_abc = Source.I_filt_poc[num_source, :, i]

    sin_θ = sin.([θ; θ - 120*π/180; θ + 120*π/180])
    cos_θ = cos.([θ; θ - 120*π/180; θ + 120*π/180])

    #---- Self-synchronisation

    if SC == 1
        SP = 0
        SQ = 0
    end

    i_c = i_abc
    #----

    #---- Integrate eq_new to find Mfif_new

    if SC == 1
        Q = -ω[end]*Mfif[end]*dot(i_c, cos_θ) # Reactive Power
    else
        Q = pq0[2]
    end

    eq_new = (1/K)*(pq0_ref[2] + SQ*Dq*(Vn - Vg) - Q)

    Mfif_new = Third_Order_Integrator(Mfif[end], μ, [eq; eq_new])
    Mfif_new = Mfif[end] + μ*eq_new

    Source.Mfif[num_source, i] = Mfif_new
    Source.eq[num_source, i] = eq_new
    #----

    #---- PI control to find Δω
    #=
    When Sp = 0, is turned ON, ΔT is controlled to be 0 in the steady state via
    the PI controller.
    =#
    if SP == 0
        ω_ref = Δω .+ ωsys
        ω_err = ω .- ω_ref
        ΔT = Dp*ω_err

        ΔT_t = Third_Order_Integrator(ΔT_t, μ, ΔT) # integration
        Δω_new = Kp*ΔT[end] + Ki*ΔT_t

        Source.Δω_sync[num_source, i] = Δω_new
        Source.ΔT_t[num_source] = ΔT_t
    else
        Δω_new = 0
        Source.Δω_sync[num_source, i] = 0
        Source.ΔT_t[num_source] = 0
    end
    #----

    #---- Integrate α_new to find ω_new
    Tm = pq0_ref[1]/ωsys # Virtual machine Torque

    ω_ref_new = Δω_new + ωsys
    ω_err_new = ω[end] - ω_ref_new
    ΔT_new = Dp*ω_err_new
    Te_new = Mfif_new*dot(i_c, sin_θ) # New Electrical Torque
    α_new = (1/J)*(Tm - Te_new - ΔT_new) # New Angular Acceleration

    ω_new = Third_Order_Integrator(ω[end], μ, [α; α_new])
    ω_new = ω[end] + μ*α_new

    Source.ω_sync[num_source, i] = ω_new
    Source.α_sync[num_source, i] = α_new
    #----

    #---- Integrate ω_new to find θ_new
    θ_new = Third_Order_Integrator(θ, μ, [ω; ω_new])%(2*π)

    Source.θ_sync[num_source, i] = θ_new
    #----

    #----
    sin_θ_new = sin.([θ_new; θ_new - 120*π/180; θ_new + 120*π/180])
    e = ω_new*Mfif_new*sin_θ_new # three phase generated voltage
    #----

    Source.V_ref[num_source, 1, i] = e[1]
    Source.V_ref[num_source, 2, i] = e[2]
    Source.V_ref[num_source, 3, i] = e[3]

    return nothing
end

function Butterworth_LPF(fc, x, y, μ)

    # 2nd Order Low Pass Butterworth Filter

    n = size(y,2)
    k = size(x,2)

    dy = size(y,1)
    dx = size(x,1)

    # Zero Padding
    if n < 2
        yd = Array{Float64, 2}(undef, dy, 2-n)
        yd = fill!(yd,0)
        y = [y yd]
        n = 2
    end
    if k < 3
        xd = Array{Float64, 2}(undef, dx, 3-k)
        xd = fill!(xd,0)
        x = [x xd]
        k = 3
    end

    #apply pre-warping transformation
    ωa = (2/μ)*tan(fc*2π*μ/2)

    #=
        Finding Coefficients of Transfer/Pulse function G(z)
        G(z) = (top_z_2*z^-2 + top_z_0*z^-1 + top_z_0)/(bot_z_2*z^-2 + bot_z_0*z^-1 + 1)

        G(Z) = Y(z)/X(z) -> output/input
        Y(z)*(bot_z_2*z^-2 + bot_z_1*z^-1 + 1) = X(z)*(top_z_2*z^-2 + top_z_1*z^-1 + top_z_0)
    =#

    ωc = ωa*μ/2
    a = ωc^2
    b = sqrt(2)*ωc
    top_z_0 = a/(a + b + 1)
    top_z_1 = 2*a/(a + b + 1)
    top_z_2 = top_z_0
    bot_z_2 = (a - b + 1)/(a + b + 1)
    bot_z_1 = (2*a - 2)/(a + b + 1)

    y_new = -1*bot_z_1.*y[:,n] .- bot_z_2*y[:,n-1] .+ top_z_0*x[:,k] .+ top_z_1*x[:,k-1] .+ top_z_2*x[:,k-2]

    return y_new
end

function First_Order_LPF(fc, x, y, μ)

    # 2nd Order Low Pass Butterworth Filter

    k = length(x)

    # Zero Padding
    if k < 2
        x = [0 x]
        k = 2
    end

    #apply pre-warping transformation
    ωa = (2/μ)*tan(fc*2π*μ/2)

    #=
        Finding Coefficients of Transfer/Pulse function G(z)
        G(z) = (ωd*z^-1 + ωd)/((ωd - 1)*z^-1 + ωd + 1)
    =#
    ωd = ωa*μ/2
    α = ωd/(ωd + 1)
    β = (ωd - 1)/(ωd + 1)

    y_new = -β*y + α*x[k] + α*x[k-1]

    return y_new
end

function Third_Order_Integrator(y_i, μ, u)

    if length(u) > 2
        y_next = y_i + (μ/12)*(23*u[3] - 16*u[2] + 5*u[1]) #3rd Order
    elseif length(u) == 2
        y_next = y_i + (μ/2)*(3*u[2] - u[1]) #2nd Order Integration
    elseif length(u) == 1
        y_next = y_i + (μ)*(u[1]) #Forward Euler
    else
        y_next = 0
    end

    return y_next
end

function Integrator_Prep(i; Order = 3)

    cnt_start = i - Order - 1
    cnt_end = i - 1
    if cnt_end < 1
        cnt_start = 1
        cnt_end = 1
    end
    if cnt_start < 1
        cnt_start = 1
    end
    range = cnt_start:cnt_end

    return range, cnt_end
end

#-------------------------------------------------------------------------------

function Measurements(Source::Classical_Controls)

    i = Source.steps

    Source.num_sources

    i_sp_rms = convert(Int64, round((1/(Source.ts*Source.T_sp_rms))))

    i_start = i - convert(Int64, round(Source.T_eval/(Source.fsys*Source.ts)))

    i_range = i_start:i_sp_rms:i

    for ns in 1:Source.num_sources

        V_poc = Source.V_filt_poc[ns, :, i]
        I_poc = Source.I_filt_poc[ns, :, i]

        Source.p_q_inst[ns, :, i] = p_q_theory(V_poc, I_poc) # real and imaginary powers

        Source.p_inst[ns, :, i] = V_poc.*I_poc

        # Phasors
        if i_range[1] >= 1 #waiting for one evaluation cycle to pass

                θ = Source.θpll[ns, 1, i_range]

                # Voltages
                v_signals = transpose(Source.V_filt_poc[ns, :, i_range])
                Source.V_ph[ns,  :, :, i] = RMS(θ, v_signals)
                Source.V_ph[ns,  :, 3, i] = Source.V_ph[ns,  :, 3, i]# .+ π/2

                # Currents
                i_signals = transpose(Source.I_filt_poc[ns, :, i_range])
                Source.I_ph[ns, :, :, i] = RMS(θ, i_signals)
                Source.I_ph[ns, :, 3, i] = Source.I_ph[ns, :, 3, i]# .+ π/2

                # Per phase (and total) Active and Reactiv Powers
                Source.Pm[ns, 1:3, i] = (Source.V_ph[ns, :, 2, i].*Source.I_ph[ns, :, 2, i]).*cos.(Source.V_ph[ns, :, 3, i] .- Source.I_ph[ns, :, 3, i])
                Source.Pm[ns, 4, i] = sum(Source.Pm[ns, 1:3, i])
                Source.Qm[ns, 1:3, i] = (Source.V_ph[ns, :, 2, i].*Source.I_ph[ns, :, 2, i]).*sin.(Source.V_ph[ns, :, 3, i] .- Source.I_ph[ns, :, 3, i])
                Source.Qm[ns, 4, i] = sum(Source.Qm[ns, 1:3, i])
        end
    end

    return nothing
end

function Current_PI_LoopShaping(Source::Classical_Controls, num_source)

    #=
        The current controller is designed for a short circuit
        The voltage controller is designed for an open circuit
        - assuming that cable capacitances can be neglected

        Controller Definitions
        ωcp - gain crossover frequency, is the frequency where the phase margin is measured,
        which is a 0-dB Gain crossing frequency
        ωcg - phase crossover frequency, is the frequency where the gain margin is measure,
        which is a -π phase crossing frequency
        pm - phase margin, the phase margin is the difference between the phase of the response
        and -π when the loop gain is 1
        gm - gain margin, the amount of gain variance required to make the loop gain unity
        at the frequency ωcg where the phase angle is -π. In other words, the gain margin
        is 1/g if g is the gain at the -π phase frequency. Negative gain margins indicate
        that stability is lost by decreasing the gain, while positive gain margins indicate
        that stability is lost by decreasing the gain, while positive gain margins indicate
        that stability is lost by increasing the gain.

        Dead Time in Digital Control Loops
        If the control scheme is implemented on a microcontroller or microprocessor,
        then a certain time is required to process the control algorithm. Therefore,
        a measured value can affect the voltage reference only after this time period
        has passed. In an appropriate manner, all these processes are synchronised with
        the clock cycle given by the pulse width modulation or vector modulation. This
        way, the digital control loop introduces a dead time of one sampling step.
        Together with the ZoH which samples the signals at the PoC, a total dead time of
        1.5 sampling steps of the current control loop results.
    =#

    #--------------------------------------
    # Current Controller

    # Short Circuit State Space
    A = [-Source.Rf[num_source]/Source.Lf[num_source]]
    B = [1/(Source.Lf[num_source])]
    C = [1]
    D = [0]
    sys = ss(A, B, C, D) # continuous
    tf_sys = tf(sys)

    Ts = 1/Source.f_cntr
    dly = Source.delay*Ts
    ZoH = tf([1], [Ts/2, 1]) # Transfer function of the sample and hold process
    Pade = tf([-dly/2, 1], [dly/2, 1]) # Pure first-order delay approximation
    PWM_gain = Source.Vdc[num_source]/2

    #SC = tf([1/(1*Lf)], [1, Rf/Lf]) # = tf(sys_sc)
    Gsc_ol = minreal(tf_sys*Pade*PWM_gain*ZoH) # Full transfer function of plant

    min_fp = 300 # Hz, minimum allowable gain cross-over frequency
    max_i = convert(Int64, floor(1/(min_fp*Ts)))

    #=
        Decreasing the cross-over frequency, making the controller slower until it is stable.
        The result is the fastest stable controller.
    =#
    for i in 6:max_i

        ωp = 2π/(i*Ts) # gain cross-over frequency
        pm = 60 # degrees, phase margin
        Gpi_i, _, _ = loopshapingPI(Gsc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)

        Gi_cl = minreal(Gsc_ol*Gpi_i/(1 + Gsc_ol*Gpi_i)) # closed loop transfer function

        if any(real(ControlSystems.poles(Gi_cl)) .> 0) == false

            # all the poles are on the left side
            ωp = 2π/((i)*Ts) # gain cross-over frequency
            pm = 60 # degrees, phase margin
            Gpi_i, kp_i, ki_i = loopshapingPI(Gsc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
            Gi_cl = minreal(Gsc_ol*Gpi_i/(1 + Gsc_ol*Gpi_i))

            Source.I_kp[num_source] = kp_i
            Source.I_ki[num_source] = ki_i
            Source.Gi_cl[num_source] = Gi_cl

            break

        end

        if i == max_i

            ωp = 2π/((i)*Ts) # gain cross-over frequency
            pm = 60 # degrees, phase margin
            Gpi_i, kp_i, ki_i = loopshapingPI(Gsc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
            Gi_cl = minreal(Gsc_ol*Gpi_i/(1 + Gsc_ol*Gpi_i))

            Source.I_kp[num_source] = kp_i
            Source.I_ki[num_source] = ki_i
            Source.Gi_cl[num_source] = Gi_cl
            
            println("\nError. PI Current Controller with Positive Poles.")
            println("Suggestion: Decrease Simulation Time Step")
            println("Source = ", num_source,"\n")
        end
    end

    return nothing
end

function Voltage_PI_LoopShaping(Source::Classical_Controls, num_source)

    #=
        The current controller is designed for a short circuit
        The voltage controller is designed for an open circuit
        - assuming that cable capacitances can be neglected

        Controller Definitions
        ωcp - gain crossover frequency, is the frequency where the phase margin is measured,
        which is a 0-dB Gain crossing frequency
        ωcg - phase crossover frequency, is the frequency where the gain margin is measure,
        which is a -π phase crossing frequency
        pm - phase margin, the phase margin is the difference between the phase of the response
        and -π when the loop gain is 1
        gm - gain margin, the amount of gain variance required to make the loop gain unity
        at the frequency ωcg where the phase angle is -π. In other words, the gain margin
        is 1/g if g is the gain at the -π phase frequency. Negative gain margins indicate
        that stability is lost by decreasing the gain, while positive gain margins indicate
        that stability is lost by decreasing the gain, while positive gain margins indicate
        that stability is lost by increasing the gain.

        Dead Time in Digital Control Loops
        If the control scheme is implemented on a microcontroller or microprocessor,
        then a certain time is required to process the control algorithm. Therefore,
        a measured value can affect the voltage reference only after this time period
        has passed. In an appropriate manner, all these processes are synchronised with
        the clock cycle given by the pulse width modulation or vector modulation. This
        way, the digital control loop introduces a dead time of one sampling step.
        Together with the ZoH which samples the signals at the PoC, a total dead time of
        1.5 sampling steps of the current control loop results.
    =#

    Goc_ol = minreal(Source.Gi_cl[num_source]*tf([1], [Source.Cf[num_source], 0]))

    Ts = 1/Source.f_cntr
    min_fp = 300 # Hz, minimum allowable gain cross-over frequency
    max_i = convert(Int64, floor(1/(min_fp*Ts)))

    #=
        Increasing the cross-over frequency, making the controller faster, until it is stable.
        The result is the slowest stable controller.
    =#
    for i in max_i:-1:1

        ωp = 2π/(i*Ts) # gain cross-over frequency
        pm = 60 # degrees, phase margin
        Gpi_v, _, _ = loopshapingPI(Goc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)

        Gv_cl = minreal(Goc_ol*Gpi_v/(1 + Goc_ol*Gpi_v)) # closed loop transfer function

        if any(real(ControlSystems.poles(Gv_cl)) .> 0) == false

            # all the poles are on the left side
            ωp = 2π/((i)*Ts) # gain cross-over frequency
            pm = 60 # degrees, phase margin
            Gpi_v, kp_v, ki_v = loopshapingPI(Goc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
            Gv_cl = minreal(Goc_ol*Gpi_v/(1 + Goc_ol*Gpi_v))

            Source.V_kp[num_source] = kp_v
            Source.V_ki[num_source] = ki_v
            Source.Gv_cl[num_source] = Gv_cl

            break
        end

        if i == 1

            ωp = 2π/((i)*Ts) # gain cross-over frequency
            pm = 60 # degrees, phase margin
            Gpi_v, kp_v, ki_v = loopshapingPI(Goc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
            Gv_cl = minreal(Goc_ol*Gpi_v/(1 + Goc_ol*Gpi_v))

            Source.V_kp[num_source] = kp_v
            Source.V_ki[num_source] = ki_v
            Source.Gv_cl[num_source] = Gv_cl

            println("\nError. PI Voltage Controller with Positive Poles.")
            println("Suggestion: Decrease Simulation Time Step")
            println("Source = ", num_source,"\n")
        end
    end

    return nothing
end

function Filter_Design(Sr, fs; Vrms = 230, Vdc = 800, ΔILf_ILf = 0.15, ΔVCf_VCf = 0.01537)

    #=
    The filtering capacitors C should be chosen such that the resonant frequency
    1/sqrt(Ls*C) is approximately sqrt(ωn * ωs), where the ωn is the angular
    frequency of the gird voltage, and ωs is the angular switching frequency
    used to turn on/off the switches)
    =#
    Ir = ΔILf_ILf
    Vr = ΔVCf_VCf

    #____________________________________________________________
    # Inductor Design
    Vorms = Vrms*1.05
    Vop = Vorms*sqrt(2)

    Zl = 3*Vorms*Vorms/Sr

    Iorms = Vorms/Zl
    Iop = Iorms*sqrt(2)

    ΔIlfmax = Ir*Iop

    Lf = Vdc/(4*fs*ΔIlfmax)

    #____________________________________________________________
    # Capacitor Design
    Vorms = Vrms*0.95
    Vop = Vorms*sqrt(2)

    Zl = 3*Vorms*Vorms/Sr

    Iorms = Vorms/Zl
    Iop = Iorms*sqrt(2)
    Ir = Vdc/(4*fs*Lf*Iop)
    ΔIlfmax = Ir*Iop
    ΔVcfmax = Vr*Vop

    Cf = ΔIlfmax/(8*fs*ΔVcfmax)

    fc = 1/(2*π*sqrt(Lf*Cf))

    return Lf, Cf, fc
end

function Source_Initialiser(env, Animo, modes; pf = 0.8)

    if isa(Animo, NamedPolicy)

        Animo.policy.state_ids
        Animo.policy.action_ids

        Source = Animo.policy.Source
        
        Mode_Keys = [k[1] for k in sort(collect(Source.Modes), by = x -> x[2])]

        source_indices = unique([parse(Int64, SubString(split(x, "_")[1], 7)) for x in Animo.policy.state_ids], dims=1)

        e = 1
        for ns in source_indices

            Srated = env.nc.parameters["source"][ns]["pwr"]
            Source.Rf[e] = env.nc.parameters["source"][ns]["R1"]
            Source.Lf[e] = env.nc.parameters["source"][ns]["L1"]
            Source.Cf[e] = env.nc.parameters["source"][ns]["C"]
            Source.Vdc[e] = env.nc.parameters["source"][ns]["vdc"]
            Source.Vrms[e] = env.nc.parameters["grid"]["v_rms"]

            Source.S[e] = Srated
            Source.P[e] = pf*Srated
            Source.Q[e] = sqrt(Srated^2 - Source.P[e]^2)

            if typeof(modes[e]) == Int64
                Source.Source_Modes[e] = Mode_Keys[modes[e]]
            else
                Source.Source_Modes[e] = modes[e]
            end

            Source.pq0_set[e, :] = [Source.P[e]; Source.Q[e]; 0]

            Source.i_max[e] = 1.15*sqrt(2)*Source.S[e]/(3*Source.Vrms[e])
            Source.v_max[e] = 1.5*Source.Vdc[e]/(sqrt(2))

            Current_PI_LoopShaping(Source, e)
            Voltage_PI_LoopShaping(Source, e)

            e += 1
        end
    else
        #-- Legacy Code
        Source = Animo.Source

        Mode_Keys = [k[1] for k in sort(collect(Source.Modes), by = x->x[2])]

        for ns in 1:Animo.num_sources

            Srated = env.nc.parameters["source"][ns]["pwr"]
            Source.Rf[ns] = env.nc.parameters["source"][ns]["R1"]
            Source.Lf[ns] = env.nc.parameters["source"][ns]["L1"]
            Source.Cf[ns] = env.nc.parameters["source"][ns]["C"]
            Source.Vdc[ns] = env.nc.parameters["source"][ns]["vdc"]
            Source.Vrms[ns] = env.nc.parameters["grid"]["v_rms"]

            Source.S[ns] = Srated
            Source.P[ns] = pf*Srated
            Source.Q[ns] = sqrt(Srated^2 - Source.P[ns]^2)

            if typeof(modes[ns]) == Int64
                Source.Source_Modes[ns] = Mode_Keys[modes[ns]]
            else
                Source.Source_Modes[ns] = modes[ns]
            end

            Source.pq0_set[ns, :] = [Source.P[ns]; Source.Q[ns]; 0]

            Source.i_max[ns] = 1.15*sqrt(2)*Source.S[ns]/(3*Source.Vrms[ns])
            Source.v_max[ns] = 1.5*Source.Vdc[ns]/(sqrt(2))

            Current_PI_LoopShaping(Source, ns)
            Voltage_PI_LoopShaping(Source, ns)
        end

        # find the indices in the state vector that correspond to the inverters
        Collect_IDs(env, Source)

    end

    return nothing
end

#= function Source_Initialiser(env, Animo, modes; pf = 0.8)

    if isa(Animo, NamedPolicy)

        Source = Animo.policy.Source
        
        Mode_Keys = [k[1] for k in sort(collect(Source.Modes), by = x -> x[2])]

        source_indices = unique([parse(Int64, SubString(split(x, "_")[1], 7)) for x in Animo.policy.state_ids], dims=1)

        e = 1
        for ns in source_indices

            Srated = env.nc.parameters["source"][ns]["pwr"]
            Source.Rf[e] = env.nc.parameters["source"][ns]["R1"]
            Source.Lf[e] = env.nc.parameters["source"][ns]["L1"]
            Source.Cf[e] = env.nc.parameters["source"][ns]["C"]
            Source.Vdc[e] = env.nc.parameters["source"][ns]["vdc"]
            Source.Vrms[e] = env.nc.parameters["grid"]["v_rms"]

            Source.S[e] = Srated
            Source.P[e] = pf*Srated
            Source.Q[e] = sqrt(Srated^2 - Source.P[e]^2)

            if typeof(modes[e]) == Int64
                Source.Source_Modes[e] = Mode_Keys[modes[e]]
            else
                Source.Source_Modes[e] = modes[e]
            end

            Source.pq0_set[e, :] = [Source.P[e]; Source.Q[e]; 0]

            Source.i_max[e] = 1.15*sqrt(2)*Source.S[e]/(3*Source.Vrms[e])
            Source.v_max[e] = 1.5*Source.Vdc[e]/2

            Current_PI_LoopShaping(Source, e)
            Voltage_PI_LoopShaping(Source, e)

            e += 1
        end
    else
    
        Println("\nError.\nClassical policy is not a named policy.\n")

    end

    return nothing
end =#

function Observer_Initialiser(Source::Classical_Controls, num_source)

    # Predictive Deadbeat Reduced-Order Observer initialisation

    #= Theory:
        Regarding the selection of observer type, there are three alternatives which are 
        the prediction-type observer, the current observer, and the reduced-order observer. 
        In the prediction-type observer, the estimated states are determined based on the
        past measurement or outputs. This means that the control signal does not utilize the 
        most current information on outputs, which leads to the inaccuracy of estimated values 
        and might cause control performance degradation. On the other hand, the reduced-order 
        observer can solve the drawback of the prediction-type observer by using the measured 
        states to estimate remaining unmeasurable states at the current time. However, if the 
        measurement variables are noisy, the imprecise measured states may influence directly 
        the feedback control inputs. For these reasons, the current full-state observer is 
        employed for the three phase LCL-filtered inverter, which yields a precise estimation 
        capability even under harmonically distorted grid voltage condition. However, in general 
        if we do have full knowledge about some states and creating a full-state observer is not 
        required. The reduced-order observer only estimates the unmeasured states. The amplification 
        of (measurement) noise in the system can be mitigated by introducing an auxiliary estimate.

        In continuous time we have seen, that the estimation error decreases asymptotically, i.e.,
        it only vanishes completely for t->∞. In contrast, we can choose K in the discrete-time 
        case such that the eigenvalues of (A - KCA) are at the origin. Finding a K with this 
        property delivers a deadbeat observer, i.e., an observer which requires the minimum 
        amount of time steps to reduce the estimation error to zero. For the n-dimensional system 
        the deadbeat observer requires n steps to estimation error to zero. 
    =#

    ns = num_source

    if Source.filter_type[ns] == "LCL"

        A = Array{Float64, 2}(undef, 3, 3)

        A = [-(Source.Rf_L1[ns] + Source.Rf_C[ns])/Source.Lf_1[ns] Source.Rf_C[ns]/Source.Lf_1[ns] -1/Source.Lf_1[ns];
             Source.Rf_C[ns]/Source.Lf_2[ns] -(Source.Rf_L2[ns] + Source.Rf_C[ns])/Source.Lf_2[ns] +1/Source.Lf_2[ns];
             1/Source.Cf[ns] -1/Source.Cf[ns] 0.0]

        B = [1/Source.Lf_1[ns]; 0; 0]
        D = [0; -1/Source.Lf_1[ns]; 0]

        #------------------------------------------------------------------------------------------------

        Ad = exp(A*Source.ts)
        Bd = inv(A)*(Ad - I)*B
        Dd = inv(A)*(Ad - I)*D

        Source.Ad[ns, :, :] = Ad[2:3, 2:3]
        Source.Bd[ns, :, :] = [Ad[2:3, 1] Bd[2:3] Dd[2:3]]
        Source.Dd[ns, :] = [Ad[1, 1] Bd[1] Dd[1]]

        Source.Cd[ns, :] = Ad[1, 2:3]

        C = transpose(Ad[1, 2:3])
        _, r = Observability(C, Source.Ad[ns, :, :])

        if r != size(Source.Ad[ns, :, :], 1)
            println("\nERROR: The system is not observable. The rank of 'O' is not equal to $(size(Source.Ad[ns, :, :], 1)).")
        end

        A1_1 = Source.Ad[ns, 1, 1]
        A1_2 = Source.Ad[ns, 1, 2]
        A2_1 = Source.Ad[ns, 2, 1]
        A2_2 = Source.Ad[ns, 2, 2]

        C1 = C[1]
        C2 = C[2]

        Source.Ko[ns, 1] = -(C2*A1_1^2 - A1_2*C1*A1_1 + A1_2*A2_1*C2 - A1_2*A2_2*C1)/(A1_2*C1^2 - A2_1*C2^2 - A1_1*C1*C2 + A2_2*C1*C2)
        Source.Ko[ns, 2] = (C1*A2_2^2 - A2_1*C2*A2_2 - A1_1*A2_1*C2 + A1_2*A2_1*C1)/(A1_2*C1^2 - A2_1*C2^2 - A1_1*C1*C2 + A2_2*C1*C2)
        
        #------------------------------------------------------------------------------------------------

        #= 
        Source.Ad[ns, :, :] = exp(A*Source.ts)
        
        Source.Bd[ns, :] = inv(A)*(Source.Ad[ns, :, :] - I)*B
        Source.Dd[ns, :] = inv(A)*(Source.Ad[ns, :, :] - I)*D

        Source.Cd[ns, :] = [1. 0. 0.] # measuring the inverter current

        C = transpose(Source.Cd[ns, :])*Source.Ad[ns, :, :]
        O, r = Observability(C, Source.Ad[ns, :, :])

        A1_1 = Source.Ad[ns, 1, 1]
        A1_2 = Source.Ad[ns, 1, 2]
        A1_3 = Source.Ad[ns, 1, 3]
        A2_1 = Source.Ad[ns, 2, 1]
        A2_2 = Source.Ad[ns, 2, 2]
        A2_3 = Source.Ad[ns, 2, 3]
        A3_1 = Source.Ad[ns, 3, 1]
        A3_2 = Source.Ad[ns, 3, 2]
        A3_3 = Source.Ad[ns, 3, 3]

        C1 = C[1]
        C2 = C[2]
        C3 = C[3]

        Source.Ko[ns, 1] = -(- A1_1^3*A2_2*C2*C3 + A1_1^3*A2_3*C2^2 - A1_1^3*A3_2*C3^2 + A1_1^3*A3_3*C2*C3 + A1_1^2*A1_2*A2_1*C2*C3 + A1_1^2*A1_2*A2_2*C1*C3 - 2*A1_1^2*A1_2*A2_3*C1*C2 + A1_1^2*A1_2*A3_1*C3^2 - A1_1^2*A1_2*A3_3*C1*C3 - A1_1^2*A1_3*A2_1*C2^2 + A1_1^2*A1_3*A2_2*C1*C2 - A1_1^2*A1_3*A3_1*C2*C3 + 2*A1_1^2*A1_3*A3_2*C1*C3 - A1_1^2*A1_3*A3_3*C1*C2 - A1_1*A1_2^2*A2_1*C1*C3 + A1_1*A1_2^2*A2_3*C1^2 + A1_1*A1_2*A1_3*A2_1*C1*C2 - A1_1*A1_2*A1_3*A2_2*C1^2 - A1_1*A1_2*A1_3*A3_1*C1*C3 + A1_1*A1_2*A1_3*A3_3*C1^2 - A1_1*A1_2*A2_1*A2_2*C2*C3 + A1_1*A1_2*A2_1*A2_3*C2^2 - 2*A1_1*A1_2*A2_1*A3_2*C3^2 + 2*A1_1*A1_2*A2_1*A3_3*C2*C3 + A1_1*A1_2*A2_2^2*C1*C3 - A1_1*A1_2*A2_2*A2_3*C1*C2 + A1_1*A1_2*A2_2*A3_1*C3^2 - A1_1*A1_2*A2_2*A3_3*C1*C3 - A1_1*A1_2*A2_3*A3_1*C2*C3 + 2*A1_1*A1_2*A2_3*A3_2*C1*C3 - A1_1*A1_2*A2_3*A3_3*C1*C2 + A1_1*A1_3^2*A3_1*C1*C2 - A1_1*A1_3^2*A3_2*C1^2 + A1_1*A1_3*A2_1*A3_2*C2*C3 - A1_1*A1_3*A2_1*A3_3*C2^2 - 2*A1_1*A1_3*A2_2*A3_1*C2*C3 + A1_1*A1_3*A2_2*A3_2*C1*C3 + A1_1*A1_3*A2_2*A3_3*C1*C2 + 2*A1_1*A1_3*A2_3*A3_1*C2^2 - 2*A1_1*A1_3*A2_3*A3_2*C1*C2 - A1_1*A1_3*A3_1*A3_2*C3^2 + A1_1*A1_3*A3_1*A3_3*C2*C3 + A1_1*A1_3*A3_2*A3_3*C1*C3 - A1_1*A1_3*A3_3^2*C1*C2 + A1_2^2*A2_1^2*C2*C3 - A1_2^2*A2_1*A2_2*C1*C3 - A1_2^2*A2_1*A2_3*C1*C2 + A1_2^2*A2_1*A3_1*C3^2 - A1_2^2*A2_1*A3_3*C1*C3 + A1_2^2*A2_2*A2_3*C1^2 - A1_2^2*A2_3*A3_1*C1*C3 + A1_2^2*A2_3*A3_3*C1^2 - A1_2*A1_3*A2_1^2*C2^2 + 2*A1_2*A1_3*A2_1*A2_2*C1*C2 - A1_2*A1_3*A2_2^2*C1^2 + A1_2*A1_3*A3_1^2*C3^2 - 2*A1_2*A1_3*A3_1*A3_3*C1*C3 + A1_2*A1_3*A3_3^2*C1^2 - A1_2*A2_1*A2_2*A3_2*C3^2 + A1_2*A2_1*A2_2*A3_3*C2*C3 + A1_2*A2_1*A2_3*A3_2*C2*C3 - A1_2*A2_1*A2_3*A3_3*C2^2 + A1_2*A2_2^2*A3_1*C3^2 - A1_2*A2_2^2*A3_3*C1*C3 - 2*A1_2*A2_2*A2_3*A3_1*C2*C3 + A1_2*A2_2*A2_3*A3_2*C1*C3 + A1_2*A2_2*A2_3*A3_3*C1*C2 + A1_2*A2_3^2*A3_1*C2^2 - A1_2*A2_3^2*A3_2*C1*C2 - A1_3^2*A2_1*A3_1*C2^2 + A1_3^2*A2_1*A3_2*C1*C2 + A1_3^2*A2_2*A3_1*C1*C2 - A1_3^2*A2_2*A3_2*C1^2 - A1_3^2*A3_1^2*C2*C3 + A1_3^2*A3_1*A3_2*C1*C3 + A1_3^2*A3_1*A3_3*C1*C2 - A1_3^2*A3_2*A3_3*C1^2 - A1_3*A2_1*A3_2^2*C3^2 + 2*A1_3*A2_1*A3_2*A3_3*C2*C3 - A1_3*A2_1*A3_3^2*C2^2 + A1_3*A2_2*A3_1*A3_2*C3^2 - A1_3*A2_2*A3_1*A3_3*C2*C3 - A1_3*A2_2*A3_2*A3_3*C1*C3 + A1_3*A2_2*A3_3^2*C1*C2 - A1_3*A2_3*A3_1*A3_2*C2*C3 + A1_3*A2_3*A3_1*A3_3*C2^2 + A1_3*A2_3*A3_2^2*C1*C3 - A1_3*A2_3*A3_2*A3_3*C1*C2)/(A1_1^2*A2_2*C1*C2*C3 - A1_1^2*A2_3*C1*C2^2 + A1_1^2*A3_2*C1*C3^2 - A1_1^2*A3_3*C1*C2*C3 - A1_1*A1_2*A2_1*C1*C2*C3 - A1_1*A1_2*A2_2*C1^2*C3 + 2*A1_1*A1_2*A2_3*C1^2*C2 - A1_1*A1_2*A3_1*C1*C3^2 + A1_1*A1_2*A3_3*C1^2*C3 + A1_1*A1_3*A2_1*C1*C2^2 - A1_1*A1_3*A2_2*C1^2*C2 + A1_1*A1_3*A3_1*C1*C2*C3 - 2*A1_1*A1_3*A3_2*C1^2*C3 + A1_1*A1_3*A3_3*C1^2*C2 + A1_1*A2_1*A2_2*C2^2*C3 - A1_1*A2_1*A2_3*C2^3 + A1_1*A2_1*A3_2*C2*C3^2 - A1_1*A2_1*A3_3*C2^2*C3 - A1_1*A2_2^2*C1*C2*C3 + A1_1*A2_2*A2_3*C1*C2^2 + A1_1*A2_2*A3_1*C2*C3^2 - A1_1*A2_2*A3_2*C1*C3^2 - A1_1*A2_3*A3_1*C2^2*C3 + A1_1*A2_3*A3_3*C1*C2^2 + A1_1*A3_1*A3_2*C3^3 - A1_1*A3_1*A3_3*C2*C3^2 - A1_1*A3_2*A3_3*C1*C3^2 + A1_1*A3_3^2*C1*C2*C3 + A1_2^2*A2_1*C1^2*C3 - A1_2^2*A2_3*C1^3 - A1_2*A1_3*A2_1*C1^2*C2 + A1_2*A1_3*A2_2*C1^3 + A1_2*A1_3*A3_1*C1^2*C3 - A1_2*A1_3*A3_3*C1^3 - A1_2*A2_1^2*C2^2*C3 + A1_2*A2_1*A2_2*C1*C2*C3 + A1_2*A2_1*A2_3*C1*C2^2 - 2*A1_2*A2_1*A3_1*C2*C3^2 + 2*A1_2*A2_1*A3_2*C1*C3^2 - A1_2*A2_2*A2_3*C1^2*C2 - A1_2*A2_2*A3_1*C1*C3^2 + A1_2*A2_2*A3_3*C1^2*C3 + 3*A1_2*A2_3*A3_1*C1*C2*C3 - 2*A1_2*A2_3*A3_2*C1^2*C3 - A1_2*A2_3*A3_3*C1^2*C2 - A1_2*A3_1^2*C3^3 + 2*A1_2*A3_1*A3_3*C1*C3^2 - A1_2*A3_3^2*C1^2*C3 - A1_3^2*A3_1*C1^2*C2 + A1_3^2*A3_2*C1^3 + A1_3*A2_1^2*C2^3 - 2*A1_3*A2_1*A2_2*C1*C2^2 + 2*A1_3*A2_1*A3_1*C2^2*C3 - 3*A1_3*A2_1*A3_2*C1*C2*C3 + A1_3*A2_1*A3_3*C1*C2^2 + A1_3*A2_2^2*C1^2*C2 + A1_3*A2_2*A3_2*C1^2*C3 - A1_3*A2_2*A3_3*C1^2*C2 - 2*A1_3*A2_3*A3_1*C1*C2^2 + 2*A1_3*A2_3*A3_2*C1^2*C2 + A1_3*A3_1^2*C2*C3^2 - A1_3*A3_1*A3_2*C1*C3^2 - A1_3*A3_1*A3_3*C1*C2*C3 + A1_3*A3_2*A3_3*C1^2*C3 + A2_1*A2_2*A3_2*C2*C3^2 - A2_1*A2_2*A3_3*C2^2*C3 - A2_1*A2_3*A3_2*C2^2*C3 + A2_1*A2_3*A3_3*C2^3 + A2_1*A3_2^2*C3^3 - 2*A2_1*A3_2*A3_3*C2*C3^2 + A2_1*A3_3^2*C2^2*C3 - A2_2^2*A3_1*C2*C3^2 + A2_2^2*A3_3*C1*C2*C3 + 2*A2_2*A2_3*A3_1*C2^2*C3 - A2_2*A2_3*A3_2*C1*C2*C3 - A2_2*A2_3*A3_3*C1*C2^2 - A2_2*A3_1*A3_2*C3^3 + A2_2*A3_1*A3_3*C2*C3^2 + A2_2*A3_2*A3_3*C1*C3^2 - A2_2*A3_3^2*C1*C2*C3 - A2_3^2*A3_1*C2^3 + A2_3^2*A3_2*C1*C2^2 + A2_3*A3_1*A3_2*C2*C3^2 - A2_3*A3_1*A3_3*C2^2*C3 - A2_3*A3_2^2*C1*C3^2 + A2_3*A3_2*A3_3*C1*C2*C3)
        Source.Ko[ns, 2] = (A1_1^2*A2_1*A2_2*C2*C3 - A1_1^2*A2_1*A2_3*C2^2 + A1_1^2*A2_1*A3_2*C3^2 - A1_1^2*A2_1*A3_3*C2*C3 - A1_1*A1_2*A2_1^2*C2*C3 - A1_1*A1_2*A2_1*A2_2*C1*C3 + 2*A1_1*A1_2*A2_1*A2_3*C1*C2 - A1_1*A1_2*A2_1*A3_1*C3^2 + A1_1*A1_2*A2_1*A3_3*C1*C3 + A1_1*A1_3*A2_1^2*C2^2 - A1_1*A1_3*A2_1*A2_2*C1*C2 + A1_1*A1_3*A2_1*A3_1*C2*C3 - 2*A1_1*A1_3*A2_1*A3_2*C1*C3 + A1_1*A1_3*A2_1*A3_3*C1*C2 + A1_1*A2_1*A2_2^2*C2*C3 - A1_1*A2_1*A2_2*A2_3*C2^2 + A1_1*A2_1*A2_2*A3_2*C3^2 - A1_1*A2_1*A2_2*A3_3*C2*C3 - A1_1*A2_2^3*C1*C3 + A1_1*A2_2^2*A2_3*C1*C2 + A1_1*A2_2*A2_3*A3_1*C2*C3 - 2*A1_1*A2_2*A2_3*A3_2*C1*C3 + A1_1*A2_2*A2_3*A3_3*C1*C2 - A1_1*A2_3^2*A3_1*C2^2 + A1_1*A2_3^2*A3_2*C1*C2 + A1_1*A2_3*A3_1*A3_2*C3^2 - A1_1*A2_3*A3_1*A3_3*C2*C3 - A1_1*A2_3*A3_2*A3_3*C1*C3 + A1_1*A2_3*A3_3^2*C1*C2 + A1_2^2*A2_1^2*C1*C3 - A1_2^2*A2_1*A2_3*C1^2 - A1_2*A1_3*A2_1^2*C1*C2 + A1_2*A1_3*A2_1*A2_2*C1^2 + A1_2*A1_3*A2_1*A3_1*C1*C3 - A1_2*A1_3*A2_1*A3_3*C1^2 - A1_2*A2_1^2*A2_2*C2*C3 + A1_2*A2_1^2*A3_2*C3^2 - A1_2*A2_1^2*A3_3*C2*C3 + A1_2*A2_1*A2_2^2*C1*C3 + A1_2*A2_1*A2_2*A2_3*C1*C2 - 2*A1_2*A2_1*A2_2*A3_1*C3^2 + 2*A1_2*A2_1*A2_2*A3_3*C1*C3 - A1_2*A2_2^2*A2_3*C1^2 + A1_2*A2_2*A2_3*A3_1*C1*C3 - A1_2*A2_2*A2_3*A3_3*C1^2 + A1_2*A2_3^2*A3_1*C1*C2 - A1_2*A2_3^2*A3_2*C1^2 - A1_2*A2_3*A3_1^2*C3^2 + 2*A1_2*A2_3*A3_1*A3_3*C1*C3 - A1_2*A2_3*A3_3^2*C1^2 - A1_3^2*A2_1*A3_1*C1*C2 + A1_3^2*A2_1*A3_2*C1^2 + A1_3*A2_1^2*A2_2*C2^2 - A1_3*A2_1^2*A3_2*C2*C3 + A1_3*A2_1^2*A3_3*C2^2 - 2*A1_3*A2_1*A2_2^2*C1*C2 + 2*A1_3*A2_1*A2_2*A3_1*C2*C3 - A1_3*A2_1*A2_2*A3_2*C1*C3 - A1_3*A2_1*A2_2*A3_3*C1*C2 + A1_3*A2_2^3*C1^2 - 2*A1_3*A2_2*A2_3*A3_1*C1*C2 + 2*A1_3*A2_2*A2_3*A3_2*C1^2 + A1_3*A2_3*A3_1^2*C2*C3 - A1_3*A2_3*A3_1*A3_2*C1*C3 - A1_3*A2_3*A3_1*A3_3*C1*C2 + A1_3*A2_3*A3_2*A3_3*C1^2 + A2_1*A2_2^2*A3_2*C3^2 - A2_1*A2_2^2*A3_3*C2*C3 - A2_1*A2_2*A2_3*A3_2*C2*C3 + A2_1*A2_2*A2_3*A3_3*C2^2 + A2_1*A2_3*A3_2^2*C3^2 - 2*A2_1*A2_3*A3_2*A3_3*C2*C3 + A2_1*A2_3*A3_3^2*C2^2 - A2_2^3*A3_1*C3^2 + A2_2^3*A3_3*C1*C3 + 2*A2_2^2*A2_3*A3_1*C2*C3 - A2_2^2*A2_3*A3_2*C1*C3 - A2_2^2*A2_3*A3_3*C1*C2 - A2_2*A2_3^2*A3_1*C2^2 + A2_2*A2_3^2*A3_2*C1*C2 - A2_2*A2_3*A3_1*A3_2*C3^2 + A2_2*A2_3*A3_1*A3_3*C2*C3 + A2_2*A2_3*A3_2*A3_3*C1*C3 - A2_2*A2_3*A3_3^2*C1*C2 + A2_3^2*A3_1*A3_2*C2*C3 - A2_3^2*A3_1*A3_3*C2^2 - A2_3^2*A3_2^2*C1*C3 + A2_3^2*A3_2*A3_3*C1*C2)/(A1_1^2*A2_2*C1*C2*C3 - A1_1^2*A2_3*C1*C2^2 + A1_1^2*A3_2*C1*C3^2 - A1_1^2*A3_3*C1*C2*C3 - A1_1*A1_2*A2_1*C1*C2*C3 - A1_1*A1_2*A2_2*C1^2*C3 + 2*A1_1*A1_2*A2_3*C1^2*C2 - A1_1*A1_2*A3_1*C1*C3^2 + A1_1*A1_2*A3_3*C1^2*C3 + A1_1*A1_3*A2_1*C1*C2^2 - A1_1*A1_3*A2_2*C1^2*C2 + A1_1*A1_3*A3_1*C1*C2*C3 - 2*A1_1*A1_3*A3_2*C1^2*C3 + A1_1*A1_3*A3_3*C1^2*C2 + A1_1*A2_1*A2_2*C2^2*C3 - A1_1*A2_1*A2_3*C2^3 + A1_1*A2_1*A3_2*C2*C3^2 - A1_1*A2_1*A3_3*C2^2*C3 - A1_1*A2_2^2*C1*C2*C3 + A1_1*A2_2*A2_3*C1*C2^2 + A1_1*A2_2*A3_1*C2*C3^2 - A1_1*A2_2*A3_2*C1*C3^2 - A1_1*A2_3*A3_1*C2^2*C3 + A1_1*A2_3*A3_3*C1*C2^2 + A1_1*A3_1*A3_2*C3^3 - A1_1*A3_1*A3_3*C2*C3^2 - A1_1*A3_2*A3_3*C1*C3^2 + A1_1*A3_3^2*C1*C2*C3 + A1_2^2*A2_1*C1^2*C3 - A1_2^2*A2_3*C1^3 - A1_2*A1_3*A2_1*C1^2*C2 + A1_2*A1_3*A2_2*C1^3 + A1_2*A1_3*A3_1*C1^2*C3 - A1_2*A1_3*A3_3*C1^3 - A1_2*A2_1^2*C2^2*C3 + A1_2*A2_1*A2_2*C1*C2*C3 + A1_2*A2_1*A2_3*C1*C2^2 - 2*A1_2*A2_1*A3_1*C2*C3^2 + 2*A1_2*A2_1*A3_2*C1*C3^2 - A1_2*A2_2*A2_3*C1^2*C2 - A1_2*A2_2*A3_1*C1*C3^2 + A1_2*A2_2*A3_3*C1^2*C3 + 3*A1_2*A2_3*A3_1*C1*C2*C3 - 2*A1_2*A2_3*A3_2*C1^2*C3 - A1_2*A2_3*A3_3*C1^2*C2 - A1_2*A3_1^2*C3^3 + 2*A1_2*A3_1*A3_3*C1*C3^2 - A1_2*A3_3^2*C1^2*C3 - A1_3^2*A3_1*C1^2*C2 + A1_3^2*A3_2*C1^3 + A1_3*A2_1^2*C2^3 - 2*A1_3*A2_1*A2_2*C1*C2^2 + 2*A1_3*A2_1*A3_1*C2^2*C3 - 3*A1_3*A2_1*A3_2*C1*C2*C3 + A1_3*A2_1*A3_3*C1*C2^2 + A1_3*A2_2^2*C1^2*C2 + A1_3*A2_2*A3_2*C1^2*C3 - A1_3*A2_2*A3_3*C1^2*C2 - 2*A1_3*A2_3*A3_1*C1*C2^2 + 2*A1_3*A2_3*A3_2*C1^2*C2 + A1_3*A3_1^2*C2*C3^2 - A1_3*A3_1*A3_2*C1*C3^2 - A1_3*A3_1*A3_3*C1*C2*C3 + A1_3*A3_2*A3_3*C1^2*C3 + A2_1*A2_2*A3_2*C2*C3^2 - A2_1*A2_2*A3_3*C2^2*C3 - A2_1*A2_3*A3_2*C2^2*C3 + A2_1*A2_3*A3_3*C2^3 + A2_1*A3_2^2*C3^3 - 2*A2_1*A3_2*A3_3*C2*C3^2 + A2_1*A3_3^2*C2^2*C3 - A2_2^2*A3_1*C2*C3^2 + A2_2^2*A3_3*C1*C2*C3 + 2*A2_2*A2_3*A3_1*C2^2*C3 - A2_2*A2_3*A3_2*C1*C2*C3 - A2_2*A2_3*A3_3*C1*C2^2 - A2_2*A3_1*A3_2*C3^3 + A2_2*A3_1*A3_3*C2*C3^2 + A2_2*A3_2*A3_3*C1*C3^2 - A2_2*A3_3^2*C1*C2*C3 - A2_3^2*A3_1*C2^3 + A2_3^2*A3_2*C1*C2^2 + A2_3*A3_1*A3_2*C2*C3^2 - A2_3*A3_1*A3_3*C2^2*C3 - A2_3*A3_2^2*C1*C3^2 + A2_3*A3_2*A3_3*C1*C2*C3)
        Source.Ko[ns, 3] = -(- A1_1^2*A2_2*A3_1*C2*C3 + A1_1^2*A2_3*A3_1*C2^2 - A1_1^2*A3_1*A3_2*C3^2 + A1_1^2*A3_1*A3_3*C2*C3 + A1_1*A1_2*A2_1*A3_1*C2*C3 + A1_1*A1_2*A2_2*A3_1*C1*C3 - 2*A1_1*A1_2*A2_3*A3_1*C1*C2 + A1_1*A1_2*A3_1^2*C3^2 - A1_1*A1_2*A3_1*A3_3*C1*C3 - A1_1*A1_3*A2_1*A3_1*C2^2 + A1_1*A1_3*A2_2*A3_1*C1*C2 - A1_1*A1_3*A3_1^2*C2*C3 + 2*A1_1*A1_3*A3_1*A3_2*C1*C3 - A1_1*A1_3*A3_1*A3_3*C1*C2 - A1_1*A2_1*A2_2*A3_2*C2*C3 + A1_1*A2_1*A2_3*A3_2*C2^2 - A1_1*A2_1*A3_2^2*C3^2 + A1_1*A2_1*A3_2*A3_3*C2*C3 + A1_1*A2_2^2*A3_2*C1*C3 - A1_1*A2_2*A2_3*A3_2*C1*C2 - A1_1*A2_2*A3_1*A3_3*C2*C3 + A1_1*A2_2*A3_2*A3_3*C1*C3 + A1_1*A2_3*A3_1*A3_3*C2^2 + A1_1*A2_3*A3_2^2*C1*C3 - 2*A1_1*A2_3*A3_2*A3_3*C1*C2 - A1_1*A3_1*A3_2*A3_3*C3^2 + A1_1*A3_1*A3_3^2*C2*C3 + A1_1*A3_2*A3_3^2*C1*C3 - A1_1*A3_3^3*C1*C2 - A1_2^2*A2_1*A3_1*C1*C3 + A1_2^2*A2_3*A3_1*C1^2 + A1_2*A1_3*A2_1*A3_1*C1*C2 - A1_2*A1_3*A2_2*A3_1*C1^2 - A1_2*A1_3*A3_1^2*C1*C3 + A1_2*A1_3*A3_1*A3_3*C1^2 + A1_2*A2_1^2*A3_2*C2*C3 - A1_2*A2_1*A2_2*A3_2*C1*C3 - A1_2*A2_1*A2_3*A3_2*C1*C2 + 2*A1_2*A2_1*A3_1*A3_3*C2*C3 - 2*A1_2*A2_1*A3_2*A3_3*C1*C3 + A1_2*A2_2*A2_3*A3_2*C1^2 + A1_2*A2_2*A3_1^2*C3^2 - A1_2*A2_2*A3_1*A3_3*C1*C3 - A1_2*A2_3*A3_1^2*C2*C3 - A1_2*A2_3*A3_1*A3_3*C1*C2 + 2*A1_2*A2_3*A3_2*A3_3*C1^2 + A1_2*A3_1^2*A3_3*C3^2 - 2*A1_2*A3_1*A3_3^2*C1*C3 + A1_2*A3_3^3*C1^2 + A1_3^2*A3_1^2*C1*C2 - A1_3^2*A3_1*A3_2*C1^2 - A1_3*A2_1^2*A3_2*C2^2 + 2*A1_3*A2_1*A2_2*A3_2*C1*C2 - 2*A1_3*A2_1*A3_1*A3_3*C2^2 + A1_3*A2_1*A3_2^2*C1*C3 + A1_3*A2_1*A3_2*A3_3*C1*C2 - A1_3*A2_2^2*A3_2*C1^2 - A1_3*A2_2*A3_1^2*C2*C3 + 2*A1_3*A2_2*A3_1*A3_3*C1*C2 - A1_3*A2_2*A3_2*A3_3*C1^2 + A1_3*A2_3*A3_1^2*C2^2 - A1_3*A2_3*A3_2^2*C1^2 - A1_3*A3_1^2*A3_3*C2*C3 + A1_3*A3_1*A3_2*A3_3*C1*C3 + A1_3*A3_1*A3_3^2*C1*C2 - A1_3*A3_2*A3_3^2*C1^2 - A2_1*A2_2*A3_2^2*C3^2 + A2_1*A2_2*A3_2*A3_3*C2*C3 + A2_1*A2_3*A3_2^2*C2*C3 - A2_1*A2_3*A3_2*A3_3*C2^2 - A2_1*A3_2^2*A3_3*C3^2 + 2*A2_1*A3_2*A3_3^2*C2*C3 - A2_1*A3_3^3*C2^2 + A2_2^2*A3_1*A3_2*C3^2 - A2_2^2*A3_2*A3_3*C1*C3 - 2*A2_2*A2_3*A3_1*A3_2*C2*C3 + A2_2*A2_3*A3_2^2*C1*C3 + A2_2*A2_3*A3_2*A3_3*C1*C2 + A2_2*A3_1*A3_2*A3_3*C3^2 - A2_2*A3_1*A3_3^2*C2*C3 - A2_2*A3_2*A3_3^2*C1*C3 + A2_2*A3_3^3*C1*C2 + A2_3^2*A3_1*A3_2*C2^2 - A2_3^2*A3_2^2*C1*C2 - A2_3*A3_1*A3_2*A3_3*C2*C3 + A2_3*A3_1*A3_3^2*C2^2 + A2_3*A3_2^2*A3_3*C1*C3 - A2_3*A3_2*A3_3^2*C1*C2)/(A1_1^2*A2_2*C1*C2*C3 - A1_1^2*A2_3*C1*C2^2 + A1_1^2*A3_2*C1*C3^2 - A1_1^2*A3_3*C1*C2*C3 - A1_1*A1_2*A2_1*C1*C2*C3 - A1_1*A1_2*A2_2*C1^2*C3 + 2*A1_1*A1_2*A2_3*C1^2*C2 - A1_1*A1_2*A3_1*C1*C3^2 + A1_1*A1_2*A3_3*C1^2*C3 + A1_1*A1_3*A2_1*C1*C2^2 - A1_1*A1_3*A2_2*C1^2*C2 + A1_1*A1_3*A3_1*C1*C2*C3 - 2*A1_1*A1_3*A3_2*C1^2*C3 + A1_1*A1_3*A3_3*C1^2*C2 + A1_1*A2_1*A2_2*C2^2*C3 - A1_1*A2_1*A2_3*C2^3 + A1_1*A2_1*A3_2*C2*C3^2 - A1_1*A2_1*A3_3*C2^2*C3 - A1_1*A2_2^2*C1*C2*C3 + A1_1*A2_2*A2_3*C1*C2^2 + A1_1*A2_2*A3_1*C2*C3^2 - A1_1*A2_2*A3_2*C1*C3^2 - A1_1*A2_3*A3_1*C2^2*C3 + A1_1*A2_3*A3_3*C1*C2^2 + A1_1*A3_1*A3_2*C3^3 - A1_1*A3_1*A3_3*C2*C3^2 - A1_1*A3_2*A3_3*C1*C3^2 + A1_1*A3_3^2*C1*C2*C3 + A1_2^2*A2_1*C1^2*C3 - A1_2^2*A2_3*C1^3 - A1_2*A1_3*A2_1*C1^2*C2 + A1_2*A1_3*A2_2*C1^3 + A1_2*A1_3*A3_1*C1^2*C3 - A1_2*A1_3*A3_3*C1^3 - A1_2*A2_1^2*C2^2*C3 + A1_2*A2_1*A2_2*C1*C2*C3 + A1_2*A2_1*A2_3*C1*C2^2 - 2*A1_2*A2_1*A3_1*C2*C3^2 + 2*A1_2*A2_1*A3_2*C1*C3^2 - A1_2*A2_2*A2_3*C1^2*C2 - A1_2*A2_2*A3_1*C1*C3^2 + A1_2*A2_2*A3_3*C1^2*C3 + 3*A1_2*A2_3*A3_1*C1*C2*C3 - 2*A1_2*A2_3*A3_2*C1^2*C3 - A1_2*A2_3*A3_3*C1^2*C2 - A1_2*A3_1^2*C3^3 + 2*A1_2*A3_1*A3_3*C1*C3^2 - A1_2*A3_3^2*C1^2*C3 - A1_3^2*A3_1*C1^2*C2 + A1_3^2*A3_2*C1^3 + A1_3*A2_1^2*C2^3 - 2*A1_3*A2_1*A2_2*C1*C2^2 + 2*A1_3*A2_1*A3_1*C2^2*C3 - 3*A1_3*A2_1*A3_2*C1*C2*C3 + A1_3*A2_1*A3_3*C1*C2^2 + A1_3*A2_2^2*C1^2*C2 + A1_3*A2_2*A3_2*C1^2*C3 - A1_3*A2_2*A3_3*C1^2*C2 - 2*A1_3*A2_3*A3_1*C1*C2^2 + 2*A1_3*A2_3*A3_2*C1^2*C2 + A1_3*A3_1^2*C2*C3^2 - A1_3*A3_1*A3_2*C1*C3^2 - A1_3*A3_1*A3_3*C1*C2*C3 + A1_3*A3_2*A3_3*C1^2*C3 + A2_1*A2_2*A3_2*C2*C3^2 - A2_1*A2_2*A3_3*C2^2*C3 - A2_1*A2_3*A3_2*C2^2*C3 + A2_1*A2_3*A3_3*C2^3 + A2_1*A3_2^2*C3^3 - 2*A2_1*A3_2*A3_3*C2*C3^2 + A2_1*A3_3^2*C2^2*C3 - A2_2^2*A3_1*C2*C3^2 + A2_2^2*A3_3*C1*C2*C3 + 2*A2_2*A2_3*A3_1*C2^2*C3 - A2_2*A2_3*A3_2*C1*C2*C3 - A2_2*A2_3*A3_3*C1*C2^2 - A2_2*A3_1*A3_2*C3^3 + A2_2*A3_1*A3_3*C2*C3^2 + A2_2*A3_2*A3_3*C1*C3^2 - A2_2*A3_3^2*C1*C2*C3 - A2_3^2*A3_1*C2^3 + A2_3^2*A3_2*C1*C2^2 + A2_3*A3_1*A3_2*C2*C3^2 - A2_3*A3_1*A3_3*C2^2*C3 - A2_3*A3_2^2*C1*C3^2 + A2_3*A3_2*A3_3*C1*C2*C3)

        if r != size(Source.Ad[ns, :, :], 1)
            println("\nERROR: The system is not observable. The rank of 'O' is not equal to $(size(Source.Ad[ns, :, :], 1)).")
        end =#

        #= println("ns = ", ns)
        println("Ad = ", Source.Ad[ns, :, :])
        println("Bd = ", Source.Bd[ns, :, :])
        println("Dd = ", Source.Dd[ns, :])
        println("Ko = ", Source.Ko[ns, :])
        println("r = ", r) =#
    end

    return nothing
end

function Observability(C, A; n = size(A,1))

    O = Array{Float64, 2}(undef, n*size(C,1), size(A,2))

    iter = size(C, 1)
    dim = 1:iter

    O[dim, :] = C

    for i in 2:n

        O[dim .+ iter, :] = O[dim, :]*A
        dim = dim .+ iter
    end

    return O, rank(O)
end

function Luenberger_Observer(Source::Classical_Controls, num_source)

    ns = num_source

    if Source.filter_type[ns] == "LCL" #&& 1 == 2

        A = Source.Ad[ns, :, :]
        B = Source.Bd[ns, :, :]
        C = transpose(Source.Cd[ns, :])
        D = transpose(Source.Dd[ns, :])

        K = Source.Ko[ns, 1:2]

        #----------------------------------------------------------------------

        #= A = Source.Ad[ns, :, :]
        B = Source.Bd[ns, :]
        D = Source.Dd[ns, :]
        C = transpose(Source.Cd[ns, :])

        K = Source.Ko[ns, 1:3] =#

        for ph in 1:1#Source.phases

            y = Source.I_filt_inv[ns, ph, end]

            yₚ = Source.I_filt_inv[ns, ph, end - 1]
            vₚ = (Source.Vdc[ns]/2)*Source.Vd_abc_new[ns, ph, end - Source.action_delay - 1]
            eₚ = Source.V_filt_poc[ns, ph, end - 1]

            uₚ = [yₚ; vₚ; eₚ]

            wp = Source.wp[ns, ph, 1:2]

            wp = (A - K*C)*wp + (B - K*D)*uₚ + K*y

            xp = wp

            Source.wp[ns, ph, 1:2] = wp

            #= Source.I_filt_poc[ns, ph, end] = xp[1]
            Source.V_filt_cap[ns, ph, end] = xp[2] =#

            if ph == 1 && ns == 1
                Source.debug[1] = Source.I_filt_inv[ns, ph, end] 
                Source.debug[2] = xp[1]
                Source.debug[3] = xp[2]
            end
           
            #----------------------------------------------------------------------

            #= y = Source.I_filt_inv[ns, ph, end] 
            vₚ = (Source.Vdc[ns]/2)*Source.Vd_abc_new[ns, ph, end - Source.action_delay - 1]
            eₚ = Source.V_filt_poc[ns, ph, end-1]

            wp = Source.wp[ns, ph, :]

            wp = (A - K*C*A)*wp + (B - K*C*B)*vₚ + (D - K*C*D)*eₚ + K*y

            Source.wp[ns, ph, :] = wp

            #= Source.I_filt_inv[ns, ph, end] = wp[1]
            Source.I_filt_poc[ns, ph, end] = wp[2]
            Source.V_filt_cap[ns, ph, end] = wp[3] =#

            if ph == 1 && ns == 1
                Source.debug[1] = wp[1]
                Source.debug[2] = wp[2]
                Source.debug[3] = wp[3]
            end =#

            #----------------------------------------------------------------------

            Vbase_LN_rms = Source.Vrms[ns]
            Ibase_rms = Source.S[ns]/(3*Source.Vrms[ns])
            I_err = abs(Source.I_filt_poc[ns, ph, end] - Source.debug[2])#/(Ibase_rms*sqrt(2))
            V_err = abs(Source.V_filt_cap[ns, ph, end] - Source.debug[3])#/(Vbase_LN_rms*sqrt(2))

            Source.debug[7] = Source.I_filt_poc[ns, ph, end] - Source.debug[1]
            Source.debug[8] = Source.I_filt_poc[ns, ph, end] - Source.debug[2]
            Source.debug[9] = Source.V_filt_cap[ns, ph, end] - Source.debug[3]
        
            Source.debug[10] = sqrt(I_err^2 + V_err^2)
            
        end

    end

    return nothing
end