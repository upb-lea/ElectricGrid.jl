
"""
    ClassicalControls

# Description
Mutable struct containing all of the variables and properties necessary to define any
classical controller
"""
mutable struct ClassicalControls

    #---------------------------------------------------------------------------
    # Physical Electrical Parameters

    Vdc::Vector{Float64} # DC link voltage
    Vrms::Vector{Float64} # nominal output voltage
    S::Vector{Float64} # rated nominal apparent power
    P::Vector{Float64} # rated nominal active power
    Q::Vector{Float64} # rated nominal active power
    pf::Vector{Float64} # power factor

    # IGBT / Switch Ratings
    i_max::Vector{Float64}
    v_max::Vector{Float64}

    # Filter values
    filter_type::Vector{String}
    Lf_1::Vector{Float64} # Filter values
    Lf_2::Vector{Float64} # Filter values
    Cf::Vector{Float64}
    Rf_L1::Vector{Float64} #inductor parasitic resistance
    Rf_L2::Vector{Float64} #inductor parasitic resistance
    Rf_C::Vector{Float64} #capacitor parasitic resistance

    #---------------------------------------------------------------------------
    # Measurements

    T_eval::Int64 # number of periods to save to memory (for averaging, e.g. RMS)
    T_sp_rms::Float64
    V_ph::Array{Float64}
    I_ph::Array{Float64}

    p_q_inst::Matrix{Float64}
    p_inst::Matrix{Float64}
    Pm::Matrix{Float64}
    Qm::Matrix{Float64}

    debug::Vector{Float64}

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

    f_avg::Matrix{Float64}
    θ_avg::Vector{Float64}

    N::Int64
    steps::Int64
    action_delay::Int64
    ramp_end::Float64
    process_start::Float64

    V_cable_loc::Matrix{Int64} # the position in the state vector where the POC Voltage is measured
    V_cap_loc::Matrix{Int64}
    I_poc_loc::Matrix{Int64}
    I_inv_loc::Matrix{Int64}

    Action_loc::Vector{Vector{Int64}}

    grid_forming::Vector{Int64}
    grid_following::Vector{Int64}

    f_source::Array{Float64}
    θ_source::Array{Float64}

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
    V_filt_cap::Array{Float64}
    I_filt_poc::Array{Float64}
    I_filt_inv::Array{Float64}
    p_q_inv::Matrix{Float64}
    p_q_poc::Matrix{Float64}

    #---------------------------------------------------------------------------
    # Current Controller

    Gi_cl::Array{TransferFunction} # Closed Loop transfer function

    I_dq0::Matrix{Float64} # DQ0 of I_filt_inv
    I_ref_dq0::Matrix{Float64}
    I_ref::Matrix{Float64}

    I_err::Array{Float64}
    I_err_t::Matrix{Float64}
    I_kp::Vector{Float64}
    I_ki::Vector{Float64}

    s_dq0_avg::Matrix{Float64}
    Vd_abc_new::Array{Float64} # output action

    s_lim::Matrix{Float64} # interim state before passing through limiter

    #---------------------------------------------------------------------------
    # Voltage Controller

    Gv_cl::Array{TransferFunction} # Closed Loop transfer function

    V_δ_set::Matrix{Float64} # set points when also in swing mode
    V_pu_set::Matrix{Float64}

    V_dq0::Matrix{Float64} # DQ0 of V_filt_cap
    V_ref_dq0::Matrix{Float64}
    V_ref::Matrix{Float64}

    V_err::Array{Float64}
    V_err_t::Matrix{Float64}
    V_kp::Vector{Float64}
    V_ki::Vector{Float64}

    I_lim::Matrix{Float64} # interim state before passing through limiter

    #---------------------------------------------------------------------------
    # Droop (Classical, Synchronverter, and VSG) Mode

    Δfmax::Float64 # The drop (increase) in frequency that causes a 100% increase (decrease) in power
    ΔEmax::Float64 # The drop (increase) in rms voltage that causes a 100% increase (decrease) in reactive power (from nominal)
    τv::Vector{Float64}  # time constant of the voltage loop
    τf::Vector{Float64}  # time constant of the frequency-droop loop

    D::Matrix{Float64} # Droop coefficients
    ω_droop::Array{Float64}
    θ_droop::Matrix{Float64}

    #---------------------------------------------------------------------------
    # PQ Mode (and PV Mode)

    pq0_set::Matrix{Float64}
    V_pre_dq0::Array{Float64} # V_dq0 is the filtered value of this
    V_dq0_inv::Matrix{Float64}  # V_dq0 filtered

    #---------------------------------------------------------------------------
    # (Self-)Synchronverter Mode

    J_sync::Vector{Float64} # Virtual Mass Moment of Inertia
    K_sync::Vector{Float64} # Reactive droop integrator gain
    ΔT_t::Vector{Float64}

    α_sync::Matrix{Float64}
    ω_sync::Matrix{Float64}
    θ_sync::Vector{Float64}
    Δω_sync::Matrix{Float64}
    eq::Matrix{Float64}
    Mfif::Vector{Float64}

    ΔT_err::Matrix{Float64}
    ΔT_err_t::Vector{Float64}
    ω_set::Vector{Float64}

    #---------------------------------------------------------------------------
    # Observer

    Observer::Vector{Bool}

    Ad_DQ::Array{Float64}
    Bd_DQ::Array{Float64}
    Cd_DQ::Array{Float64}
    Dd_DQ::Array{Float64}
    Ko_DQ::Array{Float64}
    xp_DQ::Matrix{Float64}

    Ad_0::Array{Float64}
    Bd_0::Array{Float64}
    Cd_0::Matrix{Float64}
    Dd_0::Matrix{Float64}
    Ko_0::Matrix{Float64}
    xp_0::Matrix{Float64}

    #---------------------------------------------------------------------------
    # Stochastic processes - Ornstein-Uhlenbeck

    κ::Vector{Float64} # mean reversion parameter
    σ::Vector{Float64} # Brownian motion scale (standard deviation) - sqrt(diffusion)
    γ::Vector{Float64} # asymptotoic mean
    k::Vector{Float64} # interpolation degree
    c_diff::Vector{Vector{Float64}} # interpolation degree

    X::Vector{Vector{Float64}}
    X₀::Vector{Float64}
    rol::Vector{Int64}
    cnt::Vector{Int64}

    function ClassicalControls(Vdc::Vector{Float64}, Vrms::Vector{Float64},
        S::Vector{Float64}, P::Vector{Float64}, Q::Vector{Float64}, pf::Vector{Float64},
        i_max::Vector{Float64}, v_max::Vector{Float64}, filter_type::Vector{String},
        Lf_1::Vector{Float64}, Lf_2::Vector{Float64}, Cf::Vector{Float64},
        Rf_L1::Vector{Float64}, Rf_L2::Vector{Float64}, Rf_C::Vector{Float64},
        T_eval::Int64, T_sp_rms::Float64, V_ph::Array{Float64}, I_ph::Array{Float64},
        p_q_inst::Matrix{Float64}, p_inst::Matrix{Float64}, Pm::Matrix{Float64}, Qm::Matrix{Float64},
        debug::Vector{Float64}, Modes::Dict{String, Int64}, Source_Modes::Vector{String},
        num_sources::Int64, phases::Int64,
        f_cntr::Float64, fsys::Float64, θsys::Float64,
        ts::Float64, f_avg::Matrix{Float64}, θ_avg::Vector{Float64},
        N::Int64, steps::Int64, action_delay::Int64, ramp_end::Float64, process_start::Float64,
        V_cable_loc::Matrix{Int64}, V_cap_loc::Matrix{Int64},
        I_poc_loc::Matrix{Int64}, I_inv_loc::Matrix{Int64},
        Action_loc::Vector{Vector{Int64}}, grid_forming::Vector{Int64}, grid_following::Vector{Int64},
        f_source::Array{Float64}, θ_source::Array{Float64},
        pll_err::Array{Float64}, pll_err_t::Matrix{Float64},
        vd::Array{Float64}, qvd::Array{Float64},
        fpll::Array{Float64}, θpll::Array{Float64},
        V_filt_poc::Array{Float64}, V_filt_cap::Array{Float64},
        I_filt_poc::Array{Float64}, I_filt_inv::Array{Float64},
        p_q_inv::Matrix{Float64}, p_q_poc::Matrix{Float64},
        Gi_cl::Array{TransferFunction},
        I_dq0::Matrix{Float64}, I_ref_dq0::Matrix{Float64}, I_ref::Matrix{Float64},
        I_err::Array{Float64}, I_err_t::Matrix{Float64},
        I_kp::Vector{Float64}, I_ki::Vector{Float64},
        s_dq0_avg::Matrix{Float64},Vd_abc_new::Array{Float64}, s_lim::Matrix{Float64},
        Gv_cl::Array{TransferFunction}, V_δ_set::Matrix{Float64}, V_pu_set::Matrix{Float64},
        V_dq0::Matrix{Float64}, V_ref_dq0::Matrix{Float64}, V_ref::Matrix{Float64},
        V_err::Array{Float64}, V_err_t::Matrix{Float64},
        V_kp::Vector{Float64}, V_ki::Vector{Float64},
        I_lim::Matrix{Float64},
        Δfmax::Float64, ΔEmax::Float64, τv::Vector{Float64}, τf::Vector{Float64},
        D::Matrix{Float64}, ω_droop::Array{Float64}, θ_droop::Matrix{Float64},
        pq0_set::Matrix{Float64}, V_pre_dq0::Array{Float64}, V_dq0_inv::Matrix{Float64},
        J_sync::Vector{Float64}, K_sync::Vector{Float64}, ΔT_t::Vector{Float64},
        α_sync::Matrix{Float64}, ω_sync::Matrix{Float64}, θ_sync::Vector{Float64},
        Δω_sync::Matrix{Float64}, eq::Matrix{Float64}, Mfif::Vector{Float64},
        ΔT_err::Matrix{Float64}, ΔT_err_t::Vector{Float64}, ω_set::Vector{Float64},
        Observer::Vector{Bool},
        Ad_DQ::Array{Float64}, Bd_DQ::Array{Float64}, Cd_DQ::Array{Float64}, Dd_DQ::Array{Float64},
        Ko_DQ::Array{Float64}, xp_DQ::Matrix{Float64},
        Ad_0::Array{Float64}, Bd_0::Array{Float64}, Cd_0::Matrix{Float64}, Dd_0::Matrix{Float64},
        Ko_0::Matrix{Float64}, xp_0::Matrix{Float64},
        κ::Vector{Float64}, σ::Vector{Float64}, γ::Vector{Float64}, k::Vector{Float64}, c_diff::Vector{Vector{Float64}},
        X::Vector{Vector{Float64}}, X₀::Vector{Float64}, rol::Vector{Int64}, cnt::Vector{Int64})

        new(Vdc, Vrms,
        S, P, Q, pf,
        i_max, v_max, filter_type,
        Lf_1, Lf_2, Cf,
        Rf_L1, Rf_L2, Rf_C,
        T_eval, T_sp_rms, V_ph, I_ph,
        p_q_inst, p_inst, Pm, Qm,
        debug, Modes, Source_Modes,
        num_sources, phases,
        f_cntr, fsys, θsys,
        ts, f_avg, θ_avg,
        N, steps, action_delay, ramp_end, process_start,
        V_cable_loc, V_cap_loc,
        I_poc_loc, I_inv_loc,
        Action_loc, grid_forming, grid_following,
        f_source, θ_source,
        pll_err, pll_err_t,
        vd, qvd,
        fpll, θpll,
        V_filt_poc, V_filt_cap,
        I_filt_poc, I_filt_inv,
        p_q_inv, p_q_poc,
        Gi_cl,
        I_dq0, I_ref_dq0, I_ref,
        I_err, I_err_t,
        I_kp, I_ki,
        s_dq0_avg, Vd_abc_new, s_lim,
        Gv_cl,V_δ_set, V_pu_set,
        V_dq0, V_ref_dq0, V_ref,
        V_err, V_err_t,
        V_kp, V_ki,
        I_lim,
        Δfmax, ΔEmax, τv, τf,
        D, ω_droop, θ_droop,
        pq0_set, V_pre_dq0, V_dq0_inv,
        J_sync, K_sync, ΔT_t,
        α_sync, ω_sync, θ_sync,
        Δω_sync, eq, Mfif,
        ΔT_err, ΔT_err_t, ω_set,
        Observer,
        Ad_DQ, Bd_DQ, Cd_DQ, Dd_DQ,
        Ko_DQ, xp_DQ,
        Ad_0, Bd_0, Cd_0, Dd_0,
        Ko_0, xp_0,
        κ, σ, γ, k, c_diff,
        X, X₀, rol, cnt)
    end

    function ClassicalControls(f_cntr, num_sources; phases = 3, action_delay = 1, fsys = 50.0)

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
        pf = Array{Float64, 1}(undef, num_sources)
        pf = fill!(pf, 0.8)

        i_max = Array{Float64, 1}(undef, num_sources)
        i_max = fill!(i_max, sqrt(2)*S[1]/(Vrms[1]*3))
        v_max = Array{Float64, 1}(undef, num_sources)
        v_max = fill!(v_max, 1200/(2*sqrt(2)))

        filter_type = Array{String, 1}(undef, num_sources)
        filter_type = fill!(filter_type, "LC")
        Lf_1 = Array{Float64, 1}(undef, num_sources)
        Lf_1 = fill!(Lf_1, 0)
        Lf_2 = Array{Float64, 1}(undef, num_sources)
        Lf_2 = fill!(Lf_2, 0)
        Cf = Array{Float64, 1}(undef, num_sources)
        Cf = fill!(Cf, 0)
        Rf_L1 = Array{Float64, 1}(undef, num_sources)
        Rf_L1 = fill!(Rf_L1, 0.4)
        Rf_L2 = Array{Float64, 1}(undef, num_sources)
        Rf_L2 = fill!(Rf_L2, 0.4)
        Rf_C = Array{Float64, 1}(undef, num_sources)
        Rf_C = fill!(Rf_C, 0.04)

        #---------------------------------------------------------------------------
        # Measurements

        ts = 1/f_cntr

        steps = 0
        θsys = 0.0
        ramp_end = 2/fsys
        process_start = 4/fsys

        T_eval = 1 #number of periods to average over (for rms calcs)
        N = 2#convert(Int64, round(T_eval/(fsys*ts))) + 1

        T_sp_rms = 5*fsys #samples in a second for rms calcs, x*fsys = x samples in a cycle

        # RMS Phasors
        V_ph = Array{Float64, 3}(undef, num_sources, phases, 3)
        V_ph = fill!(V_ph, 0)
        I_ph = Array{Float64, 3}(undef, num_sources, phases, 3)
        I_ph = fill!(I_ph, 0)

        # Instantaneous Real, Imaginary, and Zero powers
        p_q_inst = Array{Float64, 2}(undef, num_sources, phases)
        p_q_inst = fill!(p_q_inst, 0)

        p_inst = Array{Float64, 2}(undef, num_sources, phases) # instantaneous power at PCC
        p_inst = fill!(p_inst, 0)

        Pm = Array{Float64, 2}(undef, num_sources, phases+1) # 4th column is total
        Pm = fill!(Pm, 0)
        Qm = Array{Float64, 2}(undef, num_sources, phases+1) # 4th column is total
        Qm = fill!(Qm, 0)

        debug = Array{Float64, 1}(undef, 20) # put anything in here that needs to be debugged - for plotting
        debug = fill!(debug, 0)

        #---------------------------------------------------------------------------
        # General System & Control

        Modes = Dict(
                    "Swing" => 1,
                    "PQ" => 2,
                    "Droop" => 3,
                    "Synchronverter" => 4,
                    "Semi-Droop" => 5,
                    "Semi-Synchronverter" => 6,
                    "Voltage" => 7,
                    "Step" => 8,
                    "PV" => 9,
                    "Not Used 3" => 10)

        Source_Modes = Array{String, 1}(undef, num_sources)
        Source_Modes = fill!(Source_Modes, "Synchronverter")

        V_cable_loc = Array{Int64, 2}(undef, phases, num_sources)
        V_cable_loc = fill!(V_cable_loc, 1)
        V_cap_loc = Array{Int64, 2}(undef, phases, num_sources)
        V_cap_loc = fill!(V_cap_loc, 1)
        I_poc_loc = Array{Int64, 2}(undef, phases, num_sources)
        I_poc_loc = fill!(I_poc_loc, 1)
        I_inv_loc = Array{Int64, 2}(undef, phases, num_sources)
        I_inv_loc = fill!(I_inv_loc, 1)

        Action_loc = Vector{Vector{Int64}}(undef, 0)

        Order = 3 #integration Order

        f_avg = Array{Float64, 2}(undef, phases, Order)
        f_avg = fill!(f_avg, fsys)
        θ_avg = Array{Float64, 1}(undef, phases)
        θ_avg = fill!(θ_avg, 0.0)

        grid_forming = Array{Int64, 1}()
        grid_following = Array{Int64, 1}()

        f_source = Array{Float64, 2}(undef, num_sources, phases)
        f_source = fill!(f_source, 0.0)

        θ_source = Array{Float64, 3}(undef, num_sources, phases, N)
        θ_source = fill!(θ_source, 0.0)

        #---------------------------------------------------------------------------
        # Phase Locked Loops

        vd = Array{Float64, 3}(undef, num_sources, phases, Order) #3 for 3 phases *
        vd = fill!(vd, 0.0)
        qvd = Array{Float64, 3}(undef, num_sources, phases, Order) #3 for 3 phases *
        qvd = fill!(qvd, 0.0)

        fpll = Array{Float64, 3}(undef, num_sources, phases, Order + 1) #3 for 3 phases
        fpll = fill!(fpll, 0.0)

        θpll = Array{Float64, 3}(undef, num_sources, phases, N)
        θpll = fill!(θpll, 0.0)

        pll_err = Array{Float64, 3}(undef, num_sources, phases, Order) # 3 phases and 3rd order integration
        pll_err = fill!(pll_err, 0)
        pll_err_t = Array{Float64, 2}(undef, num_sources, phases) # PLL total integrated error
        pll_err_t = fill!(pll_err_t, 0)

        #---------------------------------------------------------------------------
        # Interface (e.g. Digital Signal Processing filters)

        V_filt_poc = Array{Float64, 3}(undef, num_sources, phases, N)
        V_filt_poc = fill!(V_filt_poc, 0)
        V_filt_cap = Array{Float64, 3}(undef, num_sources, phases, N)
        V_filt_cap = fill!(V_filt_cap, 0)

        I_filt_poc = Array{Float64, 3}(undef, num_sources, phases, N)
        I_filt_poc = fill!(I_filt_poc, 0.0)
        I_filt_inv = Array{Float64, 3}(undef, num_sources, phases, N)
        I_filt_inv = fill!(I_filt_inv, 0.0)

        p_q_inv = Array{Float64, 2}(undef, num_sources, phases)
        p_q_inv  = fill!(p_q_inv, 0)

        p_q_poc = Array{Float64, 2}(undef, num_sources, phases)
        p_q_poc  = fill!(p_q_poc, 0)

        #---------------------------------------------------------------------------
        # Current Controller

        Gi_cl = Array{TransferFunction, 1}(undef, num_sources)

        I_dq0 = Array{Float64, 2}(undef, num_sources, 3)
        I_dq0 = fill!(I_dq0, 0)
        I_ref_dq0 = Array{Float64, 2}(undef, num_sources, 3)
        I_ref_dq0 = fill!(I_ref_dq0, 0)
        I_ref = Array{Float64, 2}(undef, num_sources, phases)
        I_ref = fill!(I_ref, 0)

        # Current Integrations
        I_err = Array{Float64, 3}(undef, num_sources, phases, Order) # 3 phases and 3rd order integration
        I_err = fill!(I_err, 0)
        I_err_t = Array{Float64, 2}(undef, num_sources, phases)
        I_err_t = fill!(I_err_t, 0)
        I_kp = Array{Float64, 1}(undef, num_sources)
        I_kp = fill!(I_kp, 0.5)
        I_ki = Array{Float64, 1}(undef, num_sources)
        I_ki = fill!(I_ki, 25)

        s_dq0_avg = Array{Float64, 2}(undef, num_sources, phases)
        s_dq0_avg = fill!(s_dq0_avg, 0)

        Vd_abc_new = Array{Float64, 3}(undef, num_sources, phases, action_delay + 3) # dim 3 history terms
        Vd_abc_new = fill!(Vd_abc_new, 0)

        s_lim = Array{Float64, 2}(undef, num_sources, phases)
        s_lim = fill!(s_lim, 0.0)

        #---------------------------------------------------------------------------
        # Voltage Controller

        V_δ_set = Array{Float64, 2}(undef, num_sources, phases)
        V_δ_set = fill!(V_δ_set, 0.0)
        V_pu_set = Array{Float64, 2}(undef, num_sources, phases)
        V_pu_set = fill!(V_pu_set, 1.0)

        Gv_cl = Array{TransferFunction, 1}(undef, num_sources)

        V_dq0 = Array{Float64, 2}(undef, num_sources, 3)
        V_dq0 = fill!(V_dq0, 0)
        V_ref_dq0 = Array{Float64, 2}(undef, num_sources, 3)
        V_ref_dq0 = fill!(V_ref_dq0, 0)

        V_ref = Array{Float64, 2}(undef, num_sources, phases)
        V_ref = fill!(V_ref, 0)

        # Voltage Integrations
        V_err = Array{Float64, 3}(undef, num_sources, phases, Order) # 3 phases and 3rd order integration
        V_err = fill!(V_err, 0)
        V_err_t = Array{Float64, 2}(undef, num_sources, phases)
        V_err_t = fill!(V_err_t, 0)
        V_kp = Array{Float64, 1}(undef, num_sources)
        V_kp = fill!(V_kp, 0.01)
        V_ki = Array{Float64, 1}(undef, num_sources)
        V_ki = fill!(V_ki, 20)

        I_lim = Array{Float64, 2}(undef, num_sources, phases)
        I_lim = fill!(I_lim, 0)

        #---------------------------------------------------------------------------
        # Droop (Classical, Synchronverter, and VSG) Mode

        #=
            Typical values for the frequency droop are a 100% increase in power for a
            frequency decrease between 3% and 5% (from nominal values)
        =#

        Δfmax = 0.5/100 # Hz # The drop in frequency, Hz, which will cause a 100% increase in active power
        ΔEmax = 5/100 # V # The drop in rms voltage, which will cause a 100% decrease in reactive power

        τv = Array{Float64, 1}(undef, num_sources)
        τv = fill!(τv, 0.002 ) # time constant of the voltage loop # 0.02
        τf = Array{Float64, 1}(undef, num_sources)
        τf = fill!(τf, 0.002 ) # time constant of the frequency loop # 0.002

        D = Array{Float64, 2}(undef, num_sources, 2)
        D[:,1] = fill!(D[:,1], 2π*Δfmax/P[1])
        D[:,2] = fill!(D[:,2], ΔEmax/Q[1])

        ω_droop = Array{Float64, 3}(undef, num_sources, phases, Order) #3 for 3 phases
        ω_droop = fill!(ω_droop, fsys*2π)

        θ_droop = Array{Float64, 2}(undef, num_sources, phases)
        θ_droop = fill!(θ_droop, -ts*π*fsys)

        #---------------------------------------------------------------------------
        # PQ Mode

        pq0_set = Array{Float64, 2}(undef, num_sources, 3) # real, imaginary, and zero power set points
        pq0_set = fill!(pq0_set, 0)

        V_pre_dq0 = Array{Float64, 3}(undef, num_sources, 3, 3) # 3 for dq0, and 3 for history
        V_pre_dq0 = fill!(V_pre_dq0, 0)

        V_dq0_inv = Array{Float64, 2}(undef, num_sources, 3)
        V_dq0_inv = fill!(V_dq0_inv, 0)

        #---------------------------------------------------------------------------
        # Synchronverter Mode

        J_sync = Array{Float64, 1}(undef, num_sources)
        J_sync = fill!(J_sync, 0)
        K_sync = Array{Float64, 1}(undef, num_sources)
        K_sync = fill!(K_sync, 0)
        ΔT_t = Array{Float64, 1}(undef, num_sources)
        ΔT_t = fill!(ΔT_t, 0)

        α_sync = Array{Float64, 2}(undef, num_sources, Order)
        α_sync = fill!(α_sync, 0)
        ω_sync = Array{Float64, 2}(undef, num_sources, Order)
        ω_sync = fill!(ω_sync, fsys*2π)
        θ_sync = Array{Float64, 1}(undef, num_sources)
        θ_sync = fill!(θ_sync, -ts*π*fsys)
        Δω_sync = Array{Float64, 2}(undef, num_sources, Order)
        Δω_sync = fill!(Δω_sync, 0)
        eq = Array{Float64, 2}(undef, num_sources, Order)
        eq = fill!(eq, 0)
        Mfif = Array{Float64, 1}(undef, num_sources)
        Mfif = fill!(Mfif, 0)

        ΔT_err = Array{Float64, 2}(undef, num_sources, Order)
        ΔT_err = fill!(ΔT_err, 0)
        ΔT_err_t = Array{Float64, 1}(undef, num_sources)
        ΔT_err_t = fill!(ΔT_err_t, 0)
        ω_set = Array{Float64, 1}(undef, num_sources)
        ω_set = fill!(ω_set, 0)

        #---------------------------------------------------------------------------
        # Observers

        Observer = Array{Bool, 1}(undef, num_sources)
        Observer = fill!(Observer, false)

        Ad_0 = Array{Float64, 3}(undef, num_sources, 2, 2)
        Bd_0 = Array{Float64, 3}(undef, num_sources, 2, 3)
        Cd_0 = Array{Float64, 2}(undef, num_sources, 2)
        Dd_0 = Array{Float64, 2}(undef, num_sources, 3)
        Ko_0 = Array{Float64, 2}(undef, num_sources, 2)
        xp_0 = Array{Float64, 2}(undef, num_sources, 2)
        xp_0 = fill!(xp_0, 0)

        Ad_DQ = Array{Float64, 3}(undef, num_sources, 4, 4) #4, 4 #6, 6
        Bd_DQ = Array{Float64, 3}(undef, num_sources, 4, 6) #4, 6 #6, 2
        Cd_DQ = Array{Float64, 3}(undef, num_sources, 2, 4) #2, 4 #2, 6
        Dd_DQ = Array{Float64, 3}(undef, num_sources, 2, 6) #2, 6 #6, 2
        Ko_DQ = Array{Float64, 3}(undef, num_sources, 4, 2) #4, 2 #6, 2
        xp_DQ = Array{Float64, 2}(undef, num_sources, 4) #4 #6
        xp_DQ = fill!(xp_DQ, 0)

        #---------------------------------------------------------------------------
        # Stochastic processes - Ornstein-Uhlenbeck

        κ = Array{Float64, 1}(undef, num_sources) # mean reversion parameter
        σ = Array{Float64, 1}(undef, num_sources) # Brownian motion scale (standard deviation) - sqrt(diffusion)
        γ = Array{Float64, 1}(undef, num_sources) # asymptotoic mean
        k = Array{Float64, 1}(undef, num_sources) # interpolation degree
        c_diff = Vector{Vector{Float64}}(undef, num_sources) # polynomial coefficients

        X = Vector{Vector{Float64}}(undef, num_sources)
        X₀ = Array{Float64, 1}(undef, num_sources)
        rol = Array{Int64, 1}(undef, num_sources)
        cnt = Array{Int64, 1}(undef, num_sources)
        cnt = fill!(cnt, 0)

        ClassicalControls(Vdc, Vrms,
        S, P, Q, pf,
        i_max, v_max, filter_type,
        Lf_1, Lf_2, Cf,
        Rf_L1, Rf_L2, Rf_C,
        T_eval, T_sp_rms, V_ph, I_ph,
        p_q_inst, p_inst, Pm, Qm,
        debug, Modes, Source_Modes,
        num_sources, phases,
        f_cntr, fsys, θsys,
        ts, f_avg, θ_avg,
        N, steps, action_delay, ramp_end, process_start,
        V_cable_loc, V_cap_loc,
        I_poc_loc, I_inv_loc,
        Action_loc, grid_forming, grid_following,
        f_source, θ_source,
        pll_err, pll_err_t,
        vd, qvd,
        fpll, θpll,
        V_filt_poc, V_filt_cap,
        I_filt_poc, I_filt_inv,
        p_q_inv, p_q_poc,
        Gi_cl,
        I_dq0, I_ref_dq0, I_ref,
        I_err, I_err_t,
        I_kp, I_ki,
        s_dq0_avg, Vd_abc_new, s_lim,
        Gv_cl,V_δ_set, V_pu_set,
        V_dq0, V_ref_dq0, V_ref,
        V_err, V_err_t,
        V_kp, V_ki,
        I_lim,
        Δfmax, ΔEmax, τv, τf,
        D, ω_droop, θ_droop,
        pq0_set, V_pre_dq0, V_dq0_inv,
        J_sync, K_sync, ΔT_t,
        α_sync, ω_sync, θ_sync,
        Δω_sync, eq, Mfif,
        ΔT_err, ΔT_err_t, ω_set,
        Observer,
        Ad_DQ, Bd_DQ, Cd_DQ, Dd_DQ,
        Ko_DQ, xp_DQ,
        Ad_0, Bd_0, Cd_0, Dd_0,
        Ko_0, xp_0,
        κ, σ, γ, k, c_diff,
        X, X₀, rol, cnt)
    end
end

"""
    ClassicalPolicy()

# Description
The policy which is called when a classical controller is required.
"""
Base.@kwdef mutable struct ClassicalPolicy <: AbstractPolicy

    action_space::Space{Vector{ClosedInterval{Float64}}}
    Source::ClassicalControls

    state_ids::Vector{String}
    action_ids::Vector{String}
    Source_Indices::Vector{Int64}

    function ClassicalPolicy(action_space, Source, state_ids, action_ids, Source_Indices)
        new(action_space, Source, state_ids, action_ids, Source_Indices)
    end

    function ClassicalPolicy(env::ElectricGridEnv)

        Source_Indices = Array{Int64, 1}(undef, 0)
        Modes = Array{Any, 1}(undef, 0)

        for ns in 1:env.nc.num_sources

            if env.nc.parameters["source"][ns]["control_type"] == "classic"
                Source_Indices = [Source_Indices; Int(ns)]
                Modes = [Modes; env.nc.parameters["source"][ns]["mode"]]
            end
        end

        state_ids = GetStateIds(env.nc)
        action_ids = GetActionIds(env.nc)

        ssa = "source".*string.(Source_Indices)
        state_ids_classic = filter(x -> !isempty(findall(y -> y == split(x, "_")[1], ssa)), state_ids)
        action_ids_classic = filter(x -> !isempty(findall(y -> y == split(x, "_")[1], ssa)), action_ids)

        if env.action_delay_buffer !== nothing
            action_delay = length(env.action_delay_buffer)
        else
            action_delay = 0
        end

        if length(Source_Indices) > 0

            Source = ClassicalControls(1/env.ts,
                                        length(Source_Indices),
                                        phases = env.nc.parameters["grid"]["phase"],
                                        action_delay = action_delay,
                                        fsys = convert(Float64, env.nc.parameters["grid"]["f_grid"]))

            SourceInitialiser(env, Source, Modes, Source_Indices)

            #------------------------------------

            for s in axes(Source_Indices, 1)

                s_idx = string(Source_Indices[s])

                Source.V_cable_loc[:, s]  = findall(contains("source"*s_idx*"_v_C_cable"), state_ids_classic)
                Source.I_inv_loc[:, s] = findall(contains("source"*s_idx*"_i_L1"), state_ids_classic)

                if Source.filter_type[s] == "LC"

                    Source.V_cap_loc[:, s]  = findall(contains("source"*s_idx*"_v_C_filt"), state_ids_classic)
                    Source.I_poc_loc[:, s] = findall(contains("source"*s_idx*"_v_C_filt"), state_ids)

                elseif Source.filter_type[s] == "LCL"

                    Source.V_cap_loc[:, s]  = findall(contains("source"*s_idx*"_v_C_filt"), state_ids_classic)
                    Source.I_poc_loc[:, s] = findall(contains("source"*s_idx*"_i_L2"), state_ids_classic)

                elseif Source.filter_type[s] == "L"

                    Source.I_poc_loc[:, s] = findall(contains("source"*s_idx*"_v_C_cable"), state_ids)
                end
            end

            letterdict = Dict("a" => 1, "b" => 2, "c" => 3)

            Source.Action_loc = [[findfirst(y -> y == parse(Int64, SubString(split(x, "_")[1], 7)),
            Source_Indices), letterdict[split(x, "_")[3]]] for x in action_ids_classic]

            #------------------------------------

            animo = ClassicalPolicy(Space([-1.0..1.0 for i in 1:length(action_ids_classic)]), Source,
            state_ids_classic, action_ids_classic, Source_Indices)

            return animo

        else

            return nothing
        end
    end
end

function (Animo::ClassicalPolicy)(env::ElectricGridEnv, name::Union{String, Nothing})

    Action = ClassicalControl(Animo, env)

    return Action
end

function ResetPolicy(Animo::ClassicalPolicy)
    Source = Animo.Source

    Source.steps = 0
    Source.θsys = 0.0

    Source.f_avg = fill!(Source.f_avg, Source.fsys)
    Source.θ_avg = fill!(Source.θ_avg, 0.0)

    Source.vd = fill!(Source.vd, 0.0)
    Source.qvd = fill!(Source.qvd, 0.0)

    Source.f_source = fill!(Source.f_source, 0.0)
    Source.θ_source = fill!(Source.θ_source, 0.0)

    Source.fpll = fill!(Source.fpll, 0.0)
    Source.θpll = fill!(Source.θpll, 0.0)

    Source.pll_err = fill!(Source.pll_err, 0)
    Source.pll_err_t = fill!(Source.pll_err_t, 0)

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

    Source.V_pre_dq0 = fill!(Source.V_pre_dq0, 0)
    Source.V_dq0_inv = fill!(Source.V_dq0_inv, 0)

    Source.I_lim = fill!(Source.I_lim, 0)
    Source.s_lim = fill!(Source.s_lim, 0)

    Source.ω_droop = fill!(Source.ω_droop, Source.fsys*2π)

    Source.θ_droop = fill!(Source.θ_droop, -Source.ts*π*Source.fsys)

    Source.α_sync = fill!(Source.α_sync, 0)
    Source.ω_sync = fill!(Source.ω_sync, Source.fsys*2π)
    Source.θ_sync = fill!(Source.θ_sync, -Source.ts*π*Source.fsys)
    Source.Δω_sync = fill!(Source.Δω_sync, 0)
    Source.eq = fill!(Source.eq, 0)
    Source.Mfif = fill!(Source.Mfif, 0)

    Source.ΔT_err = fill!(Source.ΔT_err, 0)
    Source.ΔT_err_t = fill!(Source.ΔT_err_t, 0)
    Source.ω_set = fill!(Source.ω_set, 0)

    Source.xp_DQ = fill!(Source.xp_DQ, 0)
    Source.xp_0= fill!(Source.xp_0, 0)

    Source.V_filt_poc = fill!(Source.V_filt_poc, 0.0)
    Source.V_filt_cap = fill!(Source.V_filt_cap, 0.0)

    Source.I_filt_poc = fill!(Source.I_filt_poc, 0.0)
    Source.I_filt_inv = fill!(Source.I_filt_inv, 0.0)

    for ns in 1:Source.num_sources # initial conditions for stochastic process

        Source.X[ns] = fill!(Source.X[ns], Source.X₀[ns])

        n = convert(Int, round(Source.k[ns] + 1))
        if n > 0

            Δt = Source.rol[ns]*Source.ts
            t_data = Δt*collect(0:1:n)
            coef = DividedDiff(t_data, Source.X[ns])
            Source.c_diff[ns] = cat(coef, dims = 1)
        end
    end

    Source.debug = fill!(Source.debug, 0)

    return nothing
end

"""
    ClassicalControl(ClassicalPolicy, ElectricGridEnv)

# Description
Loops through all of the sources, obtaining measurments from the environment and
calculating the actions.
"""
function ClassicalControl(Animo::ClassicalPolicy, env::ElectricGridEnv)

    Source = Animo.Source
    SourceInterface(Animo, env)

    ramp_end = Source.ramp_end

    OrnsteinUhlenbeck(Source)

    Threads.@threads for ns in 1:Source.num_sources

        if Source.Source_Modes[ns] == "Swing"

            SwingMode(Source, ns, t_end = ramp_end)
        elseif Source.Source_Modes[ns] == "Voltage"

            VoltageControlMode(Source, ns, t_end = ramp_end)
        elseif Source.Source_Modes[ns] == "PQ"

            PQControlMode(Source, ns, Source.pq0_set[ns, :])
        elseif Source.Source_Modes[ns] == "PV"

            PVControlMode(Source, ns, Source.pq0_set[ns, :])
        elseif Source.Source_Modes[ns] == "Droop" || Source.Source_Modes[ns] == "Semi-Droop"

            DroopControlMode(Source, ns, t_end = ramp_end)
        elseif Source.Source_Modes[ns] == "Synchronverter"

            SynchronverterMode(Source, ns, pq0_ref = Source.pq0_set[ns, :], mode = 1, t_end = ramp_end)
        elseif Source.Source_Modes[ns] == "Semi-Synchronverter"

            SynchronverterMode(Source, ns, pq0_ref = Source.pq0_set[ns, :], mode = 2, t_end = ramp_end)
        elseif Source.Source_Modes[ns] == "Step"

            StepMode(Source, ns, t_end = ramp_end)
        elseif Source.Source_Modes[ns] == "Not Used 2"

            SynchronverterMode(Source, ns, pq0_ref = Source.pq0_set[ns, :], mode = 2, t_end = ramp_end)
        elseif Source.Source_Modes[ns] == "Not Used 3"

            SynchronverterMode(Source, ns, pq0_ref = Source.pq0_set[ns, :], mode = 2, t_end = ramp_end)
        end
    end

    Action = EnvInterface(Source)
    Measurements(Source)

    return Action
end

function roll(matrix, ns = nothing)
    if isnothing(ns)
        @views matrix[:, 1:end-1] = matrix[:, 2:end]
    else
        @views matrix[ns, :, 1:end-1] = matrix[ns, :, 2:end]
    end
end

function update_value(matrix, ns, vector)
    @views matrix[ns, :, end] = vector
end

"""
    SourceInterface(ClassicalPolicy, ElectricGridEnv)

# Description
"Measures" or "Observes" the relevant quantities necessary for control.
"""
function SourceInterface(Animo::ClassicalPolicy, env::ElectricGridEnv)

    Source = Animo.Source
    Source.steps = env.steps + 1
    ω = 2π*Source.fsys
    Source.θsys = (Source.θsys + Source.ts*ω)%(2π)

    state = env.x[findall(x -> x in Animo.state_ids, env.state_ids)]

    Threads.@threads for ns in 1:Source.num_sources

        #= # for single phase measurements
        Source.θ_source[ns, :, 1:end-1] = Source.θ_source[ns, :, 2:end]
        Source.V_filt_cap[ns, :, 1:end-1] = Source.V_filt_cap[ns, :, 2:end]
        Source.I_filt_poc[ns, :, 1:end-1] = Source.I_filt_poc[ns, :, 2:end]
        =#

        #----------------------------------------------------------------------------------
        # Rolling Windows
        roll(Source.Vd_abc_new, ns)
        roll(Source.V_filt_poc, ns)
        roll(Source.I_filt_inv, ns)

        #----------------------------------------------------------------------------------
        # New values from environment
        update_value(Source.I_filt_inv, ns, state[Source.I_inv_loc[:, ns]])

        if Source.filter_type[ns] == "LCL"

            if Source.Observer[ns]
                update_value(Source.V_filt_poc, ns, state[Source.V_cable_loc[:, ns]])
                LuenbergerObserver(Source, ns)
            else

                update_value(Source.V_filt_poc, ns, state[Source.V_cable_loc[:, ns]])

                update_value(Source.I_filt_poc, ns, state[Source.I_poc_loc[:, ns]])
                update_value(Source.V_filt_cap, ns, state[Source.V_cap_loc[:, ns]])
            end

            @views icap = Source.I_filt_inv[ns, :, end] .- state[Source.I_poc_loc[:, ns]]
            update_value(Source.V_filt_cap, ns, Source.V_filt_cap[ns, :, end] .+ Source.Rf_C[ns]*icap)

        elseif Source.filter_type[ns] == "LC"

            update_value(Source.V_filt_poc, ns, Source.V_filt_cap[ns, :, end])
            update_value(Source.I_filt_poc, ns, Source.I_filt_inv[ns, :, end] .- env.y[Source.I_poc_loc[:, ns]])

            @views icap = env.y[Source.I_poc_loc[:, ns]]
            update_value(Source.V_filt_cap, ns, state[Source.V_cable_loc[:, ns]] .+ Source.Rf_C[ns]*icap)

        elseif Source.filter_type[ns] == "L"

            update_value(Source.V_filt_cap, ns, state[Source.V_cable_loc[:, ns]])

            update_value(Source.V_filt_poc, ns, Source.V_filt_cap[ns, :, end])
            update_value(Source.I_filt_poc, ns, Source.I_filt_inv[ns, :, end] .- env.y[Source.I_poc_loc[:, ns]])
        end

        Source.p_q_inv[ns, :] =  pqTheory((Source.Vdc[ns]/2)*Source.Vd_abc_new[ns, :, end], Source.I_filt_inv[ns, :, end])

        if Source.filter_type[ns] != "L"

            Source.p_q_poc[ns, :] =  pqTheory(Source.V_filt_cap[ns, :, end], Source.I_filt_poc[ns, :, end])
        else

            Source.p_q_poc[ns, :] =  Source.p_q_inv[ns, :]
        end
    end

    return nothing
end

"""
    EnvInterface(ClassicalControls)

# Description
Passes the actions back to the environment.
"""
function EnvInterface(Source::ClassicalControls)

    Action = [Source.Vd_abc_new[Source.Action_loc[x][1], Source.Action_loc[x][2], end] for x in 1:length(Source.Action_loc)]

    return Action
end

"""
    Ramp(final, μ, i; t_end = 0.02)

# Description
Ramps up a signal.

# Arguments
- `final`: the final value.
- `μ`: time step size [s].
- `i`: time step.

# Keyword Arguments
- `t_end`: when the final value should be reached [s].
"""
function Ramp(final, μ, i; t_end = 0.02)

    if μ*i < t_end

        x_out = final.*(μ*i/t_end)
    else
        x_out = final
    end

    return x_out
end

"""
    SwingMode(Source::ClassicalControls, num_source; t_end = 0.04)

# Description
Open loop control. Produces 3 phase sinusoidal signals.
"""
function SwingMode(Source::ClassicalControls, num_source; t_end = 0.04)

    θ = Source.θsys + Source.V_δ_set[num_source, 1] - 0.5*Source.ts*2π*Source.fsys
    θph = [θ; θ - 120π/180; θ + 120π/180]

    Vrms = Ramp(Source.V_pu_set[num_source, 1]*Source.Vrms[num_source], Source.ts, Source.steps; t_end = t_end)
    Source.V_ref[num_source, :] = sqrt(2)*(Vrms)*cos.(θph)

    # debug for frequency step
    if Source.steps*Source.ts >= 2.0 && 1 == 5

        f = 60 + 0.0005*60
        ω = 2π*f
        Source.θ_source[num_source, 1, end] = (Source.θ_source[num_source, 1, end] + Source.ts*ω)%(2π)

        θ = Source.θ_source[num_source, 1, end]
        θph = [θ; θ - 120π/180; θ + 120π/180]

        Vrms = Ramp(Source.V_pu_set[num_source, 1]*Source.Vrms[num_source], Source.ts, Source.steps; t_end = t_end)
        Source.V_ref[num_source, :] = sqrt(2)*(Vrms)*cos.(θph)
    end

    # debug for voltage step
    if Source.steps*Source.ts >= 0.1 && 1 == 5

        Vdc = 150

        Source.V_ref[num_source, :] = 0.85*sqrt(2)*(Vrms)*cos.([θ; θ - 110π/180; θ + 110π/180]) .+ Vdc*[1; 1; 1]
    end

    Source.Vd_abc_new[num_source, :, end] = 2*Source.V_ref[num_source, :]/Source.Vdc[num_source]

    PhaseLockedLoop3ph(Source, num_source)

    Source.f_source[num_source, :, end] = Source.fsys*[1 1 1]
    Source.θ_source[num_source, :, end] = θph

    return nothing
end

"""
    StepMode(Source::ClassicalControls, num_source; t_end = 0.04)

# Description
Open loop control. Produces 3 stepped signals.
"""
function StepMode(Source::ClassicalControls, num_source; t_end = 0.04)

    Vrms = Ramp(Source.V_pu_set[num_source, 1]*Source.Vrms[num_source], Source.ts, Source.steps; t_end = t_end)
    Source.V_ref[num_source, :] = (Vrms)*[1; 1; 1]

    Source.Vd_abc_new[num_source, :, end] = 2*Source.V_ref[num_source, :]/Source.Vdc[num_source]

    Source.f_source[num_source, :, end] = [0 0 0]
    Source.θ_source[num_source, :, end] = [0; 0; 0]

    return nothing
end

"""
    VoltageControlMode(Source::ClassicalControls, num_source; t_end = 0.04)

# Description
Closed loop voltage control with an inner current loop. Produces 3 phase sinusoidal
voltage signals over the filter capacitor.
"""
function VoltageControlMode(Source::ClassicalControls, num_source; t_end = 0.04)

    ω = 2π*Source.fsys
    θ = Source.θsys + Source.V_δ_set[num_source, 1] - 0.5*Source.ts*ω
    θph = [θ; θ - 120π/180; θ + 120π/180]

    Vrms = Ramp(Source.V_pu_set[num_source, 1]*Source.Vrms[num_source], Source.ts, Source.steps; t_end = t_end)
    Source.V_ref[num_source, :] = sqrt(2)*Vrms*cos.([θ; θ - 120π/180; θ + 120π/180])

    VoltageController(Source, num_source, θ, ω)
    CurrentController(Source, num_source, θ, ω)

    PhaseLockedLoop3ph(Source, num_source)

    Source.f_source[num_source, :, end] = Source.fsys*[1 1 1]
    Source.θ_source[num_source, :, end] = θph

    return nothing
end

"""
    DroopControlMode(Source::ClassicalControls, num_source; t_end = 0.04)

# Description
Wrapper for simple grid forming control.
"""
function DroopControlMode(Source::ClassicalControls, num_source; t_end = 0.04)

    pu = Source.V_pu_set[num_source, 1]

    Vrms = Ramp(pu*Source.Vrms[num_source], Source.ts, Source.steps; t_end = t_end)

    DroopControl(Source, num_source, Vrms = Vrms)
    θ = Source.θ_droop[num_source, 1]
    ω = Source.ω_droop[num_source, 1, end]

    VoltageController(Source, num_source, θ, ω)
    CurrentController(Source, num_source, θ, ω)

    PhaseLockedLoop3ph(Source, num_source)

    Source.f_source[num_source, :, end] = (ω/2π)*[1 1 1]
    Source.θ_source[num_source, :, end] = [θ; θ - 120π/180; θ + 120π/180]

    return nothing
end

"""
    PQControlMode(Source::ClassicalControls, num_source, pq0)

# Description
Wrapper for simple grid following control. A controllable load on the real and imaginary power.
"""
function PQControlMode(Source::ClassicalControls, num_source, pq0)

    if norm(pq0) > Source.S[num_source]
        pq0 = pq0.*(Source.S[num_source]/norm(pq0))
    end

    PhaseLockedLoop3ph(Source, num_source)
    #PhaseLockedLoop1ph(Source, num_source, ph = 1)
    #PhaseLockedLoop1ph(Source, num_source, ph = 2)
    #PhaseLockedLoop1ph(Source, num_source, ph = 3)
    θ = Source.θpll[num_source, 1, end]
    ω = 2π*Source.fpll[num_source, 1, end]

    Filtering(Source, num_source, θ)
    #Source.V_dq0_inv[num_source, :] = DQ0Transform(Source.V_filt_poc[num_source, :, end], θ) # for when filter becomes unstable

    if Source.steps*Source.ts > Source.process_start

        PQControl(pq0_ref = pq0, Source, num_source, θ)
        CurrentController(Source, num_source, θ, ω)
    else

        PQControl(pq0_ref = [0.0; 0.0; 0.0], Source, num_source, θ)
        Source.I_ref_dq0[num_source, :] = [0.0; 0.0; 0.0]
        CurrentController(Source, num_source, θ, ω)
    end

    Source.f_source[num_source, :, end] = Source.fpll[num_source, :, end]
    Source.θ_source[num_source, :, end] = Source.θpll[num_source, :, end]

    return nothing
end

"""
    PVControlMode(Source::ClassicalControls, num_source, pq0)

# Description
Wrapper for more elaborate grid following control. A controllable load on the real power and voltage magnitude.
"""
function PVControlMode(Source::ClassicalControls, num_source, pq0)

    PhaseLockedLoop3ph(Source, num_source)
    #PhaseLockedLoop1ph(Source, num_source, ph = 1)
    #PhaseLockedLoop1ph(Source, num_source, ph = 2)
    #PhaseLockedLoop1ph(Source, num_source, ph = 3)
    θ = Source.θpll[num_source, 1, end] # positive phase sequence angle
    ω = 2π*Source.fpll[num_source, 1, end]

    Filtering(Source, num_source, θ)

    if Source.steps*Source.ts > 4/Source.fsys

        pq0_ref = PVControl(pq0_ref = pq0, Source, num_source)
        PQControl(pq0_ref = pq0_ref, Source, num_source, θ)
        CurrentController(Source, num_source, θ, ω)

    else

        PQControl(pq0_ref = [0.0; 0.0; 0.0], Source, num_source, θ)
        Source.I_ref_dq0[num_source, :] = [0.0; 0.0; 0.0]
        CurrentController(Source, num_source, θ, ω)
    end

    Source.f_source[num_source, :, end] = Source.fpll[num_source, :, end]
    Source.θ_source[num_source, :, end] = Source.θpll[num_source, :, end]

    return nothing
end

"""
    SynchronverterMode(Source::ClassicalControls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]], t_end = 0.04, mode = 2)

# Description
Wrapper for enhanced grid forming control.
"""
function SynchronverterMode(Source::ClassicalControls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]], t_end = 0.04, mode = 2)

    pu = Source.V_pu_set[num_source, 1]

    if norm(pq0_ref) > Source.S[num_source]
        pq0_ref = pq0_ref.*(Source.S[num_source]/norm(pq0_ref))
    end

    Vrms = Ramp(pu*Source.Vrms[num_source], Source.ts, Source.steps; t_end = t_end)

    SynchronverterControl(Source, num_source, pq0_ref = pq0_ref, Vrms = Vrms, mode = mode)

    VoltageController(Source, num_source, Source.θ_sync[num_source], Source.ω_sync[num_source, end])
    CurrentController(Source, num_source, Source.θ_sync[num_source], Source.ω_sync[num_source, end])

    PhaseLockedLoop3ph(Source, num_source)

    ω = Source.ω_sync[num_source, end]
    θ = Source.θ_sync[num_source]

    Source.f_source[num_source, :, end] = (ω/2π)*[1 1 1]
    Source.θ_source[num_source, :, end] = [θ; θ - 120π/180; θ + 120π/180]

    return nothing
end

"""
    PhaseLockedLoop3ph(Source::ClassicalControls, num_source; ωn = 70, ξ = 0.35)

# Description
Tuned 3 phase Phase Locked loop.
"""
function PhaseLockedLoop3ph(Source::ClassicalControls, num_source; ωn = Source.fsys + 20, ξ = 0.35)

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

        v_αβγ = ClarkeTransform(v_abc)
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

    Ki = ωn^2 # tuning
    Kp = ξ*2*sqrt(Ki) # tuning

    v_abc = Source.V_filt_cap[num_source, :, end]

    f = Source.fpll[num_source, 1, :]
    θ = Source.θpll[num_source, 1, end]

    err_t = Source.pll_err_t[num_source, 1]
    err = Source.pll_err[num_source, 1, :]

    v_αβγ = ClarkeTransform(v_abc)

    if norm(v_αβγ) != 0
        v_αβγ = sqrt(3)*v_αβγ./norm(v_αβγ)
    end

    err_new = v_αβγ[2]*cos(θ) - v_αβγ[1]*sin(θ)

    f_new, err_t_new, err_int =
    PIController(err_new, err, err_t, Kp, Ki, Source.ts, bias = Source.fsys, max_t_err = 0.00015)

    θ = ThirdOrderIntegrator(θ, Source.ts, 2π*[f[2:end]; f_new])

    Source.fpll[num_source, 1, :] = [f[2:end]; f_new[1]]
    Source.fpll[num_source, 2, :] = [f[2:end]; f_new[1]]
    Source.fpll[num_source, 3, :] = [f[2:end]; f_new[1]]

    Source.θpll[num_source, 1, :] = [Source.θpll[num_source, 1, 2:end]; (θ.%(2π))]
    Source.θpll[num_source, 2, :] = [Source.θpll[num_source, 2, 2:end]; ((θ - 120π/180).%(2π))]
    Source.θpll[num_source, 3, :] = [Source.θpll[num_source, 3, 2:end]; ((θ + 120π/180).%(2π))]

    Source.pll_err_t[num_source, :] = [err_t_new; 0; 0]
    Source.pll_err[num_source, :, :] = [transpose(err_int); transpose(err_int); transpose(err_int)]

    return nothing
end

function PhaseLockedLoop1ph(Source::ClassicalControls, num_source; Kp = 0.001, Ki = 1, ph = 1, k_sogi = 0.8)

    i = Source.steps

    range_1, cnt_end_1 = Integrator_Prep(i)
    range_2, _ = Integrator_Prep(i-1)

    v_ph = Source.V_filt_cap[num_source, ph, end]
    v_ph_r = Source.V_filt_cap[num_source, ph, range_1]
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
    vd_new = ThirdOrderIntegrator(vd_1[end], Source.ts, vd_sogi)

    qvd_int_old = vd_2.*ω_2
    qvd_int_new = vd_1[end]*ω_1[end]
    qvd_int = [qvd_int_old; qvd_int_new]
    qvd_new = ThirdOrderIntegrator(qvd_1[end], Source.ts, qvd_int)
    #----

    α_β_0 = [vd_new; qvd_new; 0]
    d_q_0 = ParkTransform(α_β_0, θ - π)

    vq = d_q_0[1]

    err_new = 0 - vq

    ω_new, err_t_new, err_int =
    PIController(err_new, err, err_t, Kp, Ki, Source.ts, bias = 2*π*Source.fsys)

    Source.θpll[num_source, ph, i] =
    (ThirdOrderIntegrator(θ, Source.ts, [ω_1; ω_new]))%(2*π)

    Source.pll_err[num_source, ph, :] = err_int
    Source.vd[num_source, ph, i] = vd_new
    Source.qvd[num_source, ph, i] = qvd_new

    Source.pll_err_t[num_source, ph] = err_t_new[end]
    Source.fpll[num_source, ph, i] = ω_new[end]/(2*π)

    return nothing
end

function SynchronverterControl(Source::ClassicalControls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0], mode = 2, Vrms = Source.Vrms[num_source])

    #= Modes:
            "Synchronverter Modes" - grid forming with power balancing via virtual motor (advanced controllable source/load)
        1 -> "Full-Synchronverter Mode" - droop control on real and imaginary powers
        2 -> "Infinite-Synchronverter Mode" - droop characteristic on real power, and active control on voltage
        *3 -> "Self-Synchronverter Mode" - active control on real and imaginary powers
        *4 -> "Semi-Synchronverter Mode" - droop characteristic on imaginary power, and active control on real power
        *5 -> "Null-Synchronverter Mode" - active control on real power and voltage

        * problems
    =#

    if mode == 1 # no changes
        S1 = 0 # PI = 0
        S2 = 1 # Dq = 1
    elseif mode == 2
        S1 = 0 # PI = 1
        S2 = 1 # Dq = 1 && 88
    elseif mode == 3
        S1 = 1 # PI = 0
        S2 = 0 # Dq = 1 && 88
    elseif mode == 4
        S1 = 1 # PI = 1
        S2 = 1 # Dq = 1
    elseif mode == 5
        S1 = 1 # PI = 1
        S2 = 1 # Dq = 1 && 88
    end

    Dq = S2*Source.D[num_source, 2]

    eq = Source.eq[num_source, :]
    α = Source.α_sync[num_source, :]
    ω = Source.ω_sync[num_source, :]
    θ = Source.θ_sync[num_source]

    ωsys = Source.fsys*2π # nominal grid frequency
    Vn = sqrt(2)*Vrms # nominal peak POC voltage
    Vg = sqrt(2)*ClarkeMag(Source.V_filt_cap[num_source, :, end]) # peak measured voltage

    #---- Integrate eq_new to find Mfif_new

    if mode == 2 || mode == 5
        eq_new = (1/Source.K_sync[num_source])*(Dq*(Vn - Vg))
    else
        eq_new = (1/Source.K_sync[num_source])*(pq0_ref[2] + Dq*(Vn - Vg) - Source.p_q_poc[num_source, 2])
    end

    Mfif_new = ThirdOrderIntegrator(Source.Mfif[num_source], Source.ts, [eq[2:end]; eq_new])

    Source.Mfif[num_source] = Mfif_new
    Source.eq[num_source, :] = [eq[2:end]; eq_new]
    #----

    #---- Integrate α_new to find ω_new
    Tm = pq0_ref[1]/ω[end] # Virtual machine Torque

    ω_err_new = ω[end] - ωsys - Source.ω_set[num_source]
    ΔT = Source.D[num_source, 1]*ω_err_new

    #~~~ PI Controller

    if S1 == 1 && 1 == 2

        Kp = 0.0001
        Ki = 0.001
        ω_set, ΔT_err_t, Source.ΔT_err[num_source, :] =
        PIController([-ΔT], Source.ΔT_err[num_source, :], Source.ΔT_err_t[num_source],
        Kp, Ki, Source.ts)

        Source.ω_set[num_source] = ω_set[1]
        Source.ΔT_err_t[num_source] = ΔT_err_t[1]
    else

        Source.ω_set[num_source] = 0
        Source.ΔT_err_t[num_source] = 0
        Source.ΔT_err[num_source, :] = fill!(Source.ΔT_err[num_source, :], 0.0)
    end
    #~~~

    Te_new = Source.p_q_poc[num_source, 1]/ω[end] # New Electrical Torque

    α_new = (1/Source.J_sync[num_source])*(Tm - Te_new - ΔT) # New Angular Acceleration
    Source.α_sync[num_source, :] = [α[2:end]; α_new]

    ω_new = ThirdOrderIntegrator(ω[end], Source.ts, Source.α_sync[num_source, :])

    Source.ω_sync[num_source, :] = [ω[2:end]; ω_new]

    #----

    #---- Integrate ω_new to find θ_new
    θ_new = ThirdOrderIntegrator(θ, Source.ts, Source.ω_sync[num_source, :])%(2π)

    Source.θ_sync[num_source] = θ_new
    #----

    #----
    cos_θ_new = cos.([θ_new; θ_new - 120*π/180; θ_new + 120*π/180])
    Source.V_ref[num_source, :] = ω_new*Mfif_new*cos_θ_new # three phase generated voltage
    #----

    return nothing
end

function PQControl(Source::ClassicalControls, num_source, θ; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0])

    #-------------------------------------------------------------

    V_αβγ = InvParkTransform(Source.V_dq0_inv[num_source, :], θ)
    I_αβγ_ref = Inv_p_q_v(V_αβγ, pq0_ref)
    Source.I_ref_dq0[num_source, :] = ParkTransform(I_αβγ_ref, θ)

    Ip_ref = sqrt(2/3)*norm(Source.I_ref_dq0[num_source, :]) # peak set point

    if Ip_ref > 0.98*Source.i_max[num_source]

        Source.I_ref_dq0[num_source, :] = Source.I_ref_dq0[num_source, :]*(0.98*Source.i_max[num_source]/(sqrt(2/3)*norm(Source.I_ref_dq0[num_source, :])))
    end

    #-------------------------------------------------------------
    I_αβγ = ClarkeTransform(Source.I_filt_inv[num_source, :, end])

    V_αβγ_ref = Inv_p_q_i(I_αβγ, pq0_ref)

    V_dq0_ref = ParkTransform(V_αβγ_ref, θ)

    if sqrt(2/3)*norm(V_dq0_ref) > Source.v_max[num_source]
        V_dq0_ref = V_dq0_ref.*((Source.v_max[num_source])/(sqrt(2/3)*norm(V_dq0_ref) ))
    end

    Source.V_ref[num_source, :] = InvDQ0Transform(V_dq0_ref, θ)

    Source.V_ref_dq0[num_source, :] = V_dq0_ref

    #-------------------------------------------------------------

    return nothing
end

function PVControl(Source::ClassicalControls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0])

    Vn = sqrt(2)*Source.V_pu_set[num_source, 1]*Source.Vrms[num_source] #peak
    #Vg = sqrt(2/3)*norm(Source.V_dq0_inv[num_source, :]) #peak
    Vg = sqrt(2/3)*norm(DQ0Transform(Source.V_filt_cap[num_source, :, end], 0))

    Kp = 200000*Source.V_kp[num_source]
    Ki = 5000*Source.V_ki[num_source]

    V_err = Source.V_err[num_source, :, 1]
    V_err_t = Source.V_err_t[num_source, 1]

    V_err_new = Vn - Vg

    q_ref, V_err_t, Source.V_err[num_source, :, 1] =
    PIController(V_err_new, V_err, V_err_t, Kp, Ki, Source.ts, bias = 0)

    pq0_ref[2] = q_ref[1]
    Source.V_err_t[num_source, 1] = V_err_t[1]

    return pq0_ref
end

function DroopControl(Source::ClassicalControls, num_source; Vrms = Source.Vrms[num_source])

    ω = Source.ω_droop[num_source, 1, :]
    θ = Source.θ_droop[num_source, 1]
    p_q = Source.p_q_poc[num_source, :]

    ωsys = Source.fsys*2π

    ω_new = ωsys - p_q[1]/(ωsys*Source.D[num_source, 1])
    Source.ω_droop[num_source, 1, 1:2] = ω[2:end]
    Source.ω_droop[num_source, :, end] = [ω_new; ω_new; ω_new]

    θ_new = ThirdOrderIntegrator(θ, Source.ts, [ω[2:end]; ω_new])%(2π)

    Source.θ_droop[num_source, :] = [θ_new; θ_new - 120*π/180; θ_new + 120*π/180].%(2π)

    e = sqrt(2)*Vrms - p_q[2]/Source.D[num_source, 2]

    Source.V_ref[num_source, :] = e*cos.(Source.θ_droop[num_source, :])

    return nothing
end

"""
    CurrentController(Source::ClassicalControls, num_source, θ, ω; Kb = 1)

# Description
Inner current control with anti-windup.
"""
function CurrentController(Source::ClassicalControls, num_source, θ, ω; Kb = 0.5)

    Kp = Source.I_kp[num_source]
    Ki = Source.I_ki[num_source]

    Source.I_ref[num_source, :] = InvDQ0Transform(Source.I_ref_dq0[num_source, :], θ)
    Source.I_dq0[num_source, :] = DQ0Transform(Source.I_filt_inv[num_source, :, end], θ)
    Source.V_dq0[num_source, :] = DQ0Transform(Source.V_filt_cap[num_source, :, end], θ)
    V_dq0 = Source.V_dq0[num_source, :]

    I_dq0 = Source.I_dq0[num_source, :]
    I_ref_dq0 = Source.I_ref_dq0[num_source, :]
    I_err = Source.I_err[num_source, :, :]
    I_err_t = Source.I_err_t[num_source, :]

    if Source.steps > 1
        # Including Anti-windup - Back-calculation
        I_err_new = I_ref_dq0 .- I_dq0 .+ Kb*(Source.Vdc[num_source]/2)*(Source.s_dq0_avg[num_source, :] .- Source.s_lim[num_source, :])
    else
        I_err_new = I_ref_dq0 .- I_dq0
    end

    max_I_err = Source.i_max[num_source]/sqrt(2) # maximum rms error

    if norm(I_err_new) > max_I_err/sqrt(1/3)

        I_err_new = I_err_new.*(max_I_err/(sqrt(1/3)*norm(I_err_new)))
    end

    Source.s_lim[num_source, :], Source.I_err_t[num_source, :], Source.I_err[num_source, :, :] =
    PIController(I_err_new, I_err, I_err_t, Kp, Ki, Source.ts, max_t_err = 0.3*sqrt(3))

    # cross-coupling / feedforward
    Source.s_lim[num_source, 1] = Source.s_lim[num_source, 1] -
    (Source.Lf_1[num_source]*ω*I_dq0[2] - V_dq0[1])*2/Source.Vdc[num_source]
    Source.s_lim[num_source, 2] = Source.s_lim[num_source, 2] +
    (Source.Lf_1[num_source]*ω*I_dq0[1] + V_dq0[2])*2/Source.Vdc[num_source]

    # ---- Limiting Output (Saturation)
    Vp_ref = (Source.Vdc[num_source]/2)*sqrt(2)*ClarkeMag(Source.s_lim[num_source,:]) # peak set point

    if Vp_ref > Source.v_max[num_source]
        Source.s_dq0_avg[num_source, :]  = Source.s_lim[num_source, :]*Source.v_max[num_source]/Vp_ref
    else
        Source.s_dq0_avg[num_source, :]  = Source.s_lim[num_source, :]
    end

    Source.Vd_abc_new[num_source, :, end] = InvDQ0Transform(Source.s_dq0_avg[num_source, :] , θ)

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

"""
    VoltageController(Source::ClassicalControls, num_source, θ, ω; Kb = 1)

# Description
Outer voltage control with anti-windup.
"""
function VoltageController(Source::ClassicalControls, num_source, θ, ω; Kb = 1)

    Kp = Source.V_kp[num_source]
    Ki = Source.V_ki[num_source]

    Source.V_ref_dq0[num_source, :] = DQ0Transform(Source.V_ref[num_source, :], θ)
    Source.V_dq0[num_source, :] = DQ0Transform(Source.V_filt_cap[num_source, :, end], θ)

    I_dq0_poc = DQ0Transform(Source.I_filt_poc[num_source, :, end], θ)
    V_dq0 = Source.V_dq0[num_source, :]
    V_ref_dq0 = Source.V_ref_dq0[num_source, :]

    V_err = Source.V_err[num_source, :, :]
    V_err_t = Source.V_err_t[num_source, :]

    if Source.steps > 1
        # Including Anti-windup - Back-calculation
        V_err_new = V_ref_dq0 .- V_dq0 .+ Kb*(Source.I_ref_dq0[num_source, :] .- Source.I_lim[num_source, :])
    else
        V_err_new = V_ref_dq0 .- V_dq0
    end

    max_V_err = Source.v_max[num_source]/sqrt(2) # maximum rms error

    if norm(V_err_new) > max_V_err/sqrt(1/3)

        V_err_new = V_err_new.*(max_V_err/(sqrt(1/3)*norm(V_err_new)))
    end

    Source.I_lim[num_source, :], Source.V_err_t[num_source, :], Source.V_err[num_source, :, :] =
    PIController(V_err_new, V_err, V_err_t, Kp, Ki, Source.ts, max_t_err = 3*sqrt(3))

    # cross-coupling / feedforward
    Source.I_lim[num_source, 1] = Source.I_lim[num_source, 1] + I_dq0_poc[1]
    - Source.Cf[num_source]*ω*V_dq0[2]
    Source.I_lim[num_source, 2] = Source.I_lim[num_source, 2] + I_dq0_poc[2]
    + Source.Cf[num_source]*ω*V_dq0[1]

    # ---- Limiting Output (Saturation)
    Ip_ref = sqrt(2)*ClarkeMag(Source.I_lim[num_source,:]) # peak set point

    if Ip_ref > 0.98*Source.i_max[num_source]
        Source.I_ref_dq0[num_source, :] = Source.I_lim[num_source, :]*0.98*Source.i_max[num_source]/Ip_ref
    else
        Source.I_ref_dq0[num_source, :] = Source.I_lim[num_source, :]
    end

    return nothing
end

"""
    PIController(Error_new, Error_Hist, Error_t, Kp, Ki, μ; bias = zeros(length(Error_new)))

# Description
Generic PI controller.
"""
function PIController(Error_new, Error_Hist, Error_t, Kp, Ki, μ; bias = zeros(length(Error_new)), max_t_err = nothing)

    d = length(Error_new)

    if d == 1
        Error_Hist = transpose(Error_Hist)
    end

    Err_t_new = Array{Float64, 1}(undef, d)
    Err_t_new = fill!(Err_t_new, 0)

    Err_d = [Error_Hist Error_new]

    Err_int = Err_d[:, 2:end]

    for j in 1:d
        Err_t_new[j] = ThirdOrderIntegrator(Error_t[j], μ, Err_int[j,:]) # integration
        #Err_t_new[j] = Error_t[j] + μ*Err_int[j, end] # integration
    end

    if !isnothing(max_t_err)

        if norm(Err_t_new) > max_t_err
            Err_t_new = Err_t_new.*(max_t_err/norm(Err_t_new))
        end
    end

    Action = bias + Kp.*Error_new .+ Ki.*Err_t_new

    return Action, Err_t_new, Err_int
end

"""
    ButterworthLPF(fc, x, y, μ) - Low Pass Filter
"""
function ButterworthLPF(fc, x, y, μ)

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

"""
    FirstOrderLPF(fc, x, y, μ) - Low Pass Filter
"""
function FirstOrderLPF(fc, x, y, μ)

    k = size(x,2)

    #= # Zero Padding
    if k < 2
        x = [0 x]
        k = 2
    end =#

    #apply pre-warping transformation
    ωa = (2/μ)*tan(fc*2π*μ/2)

    #=
        Finding Coefficients of Transfer/Pulse function G(z)
        G(z) = (ωd*z^-1 + ωd)/((ωd - 1)*z^-1 + ωd + 1)
    =#
    ωd = ωa*μ/2
    α = ωd/(ωd + 1)
    β = (ωd - 1)/(ωd + 1)

    if size(x, 1) == 1
        y_new = -β*y + α*x[k] + α*x[k-1]
    else
        y_new = -β*y .+ α*x[:, k] .+ α*x[:, k-1]
    end

    return y_new
end

"""
    ThirdOrderIntegrator(y_i, μ, u)
"""
function ThirdOrderIntegrator(y_i, μ, u)

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

"""
    Filtering(Source::ClassicalControls, num_source, θ)

# Description
First order Low pass filter on voltage DQ0 components.
"""
function Filtering(Source::ClassicalControls, num_source, θ)

    V_inv = (Source.Vdc[num_source]/2)*Source.Vd_abc_new[num_source, :, end]

    roll(Source.V_pre_dq0, num_source)
    Source.V_pre_dq0[num_source, :, end] = DQ0Transform(V_inv, θ)

    Source.V_dq0_inv[num_source, :] = FirstOrderLPF(100, Source.V_pre_dq0[num_source, :, :],
    Source.V_dq0_inv[num_source, :], Source.ts)

    return nothing
end

"""
    LuenbergerObserver(Source::ClassicalControls, num_source)

# Description
Discrete time approximate deadbeat Luenberger observer operating in DQ0 frame
"""
function LuenbergerObserver(Source::ClassicalControls, num_source)

    ns = num_source

    f_sys = Source.f_source[ns, 1, end]
    ω = 2π*f_sys
    θ = Source.θ_source[ns, 1, end]

    if Source.filter_type[ns] == "LCL"

        I_poc_DQ0 = [0.; 0.; 0.]
        V_cap_DQ0 = [0.; 0.; 0.]

        y_DQ0 = DQ0Transform(Source.I_filt_inv[ns, :, end], θ + 0.5*Source.ts*ω)

        vₚ_DQ0 = DQ0Transform((Source.Vdc[ns]/2)*Source.Vd_abc_new[ns, :, end - Source.action_delay - 1], θ)
        yₚ_DQ0 = DQ0Transform(Source.I_filt_inv[ns, :, end - 1], θ - 0.5*Source.ts*ω)
        eₚ_DQ0 = DQ0Transform(Source.V_filt_poc[ns, :, end - 1], θ - 0.5*Source.ts*ω)

        #----------------------------------------------------------------------
        # Zero component

        A = Source.Ad_0[ns, :, :]
        B = Source.Bd_0[ns, :, :]
        C = transpose(Source.Cd_0[ns, :])
        D = transpose(Source.Dd_0[ns, :])

        K = Source.Ko_0[ns, :]

        uₚ = [yₚ_DQ0[3]; vₚ_DQ0[3]; eₚ_DQ0[3]]

        Source.xp_0[ns, :] = (A - K*C)*Source.xp_0[ns, :] + (B - K*D)*uₚ + K*y_DQ0[3]

        I_poc_DQ0[3] = Source.xp_0[ns, 1]
        V_cap_DQ0[3] = Source.xp_0[ns, 2]

        #----------------------------------------------------------------------
        # DQ components

        A = Source.Ad_DQ[ns, :, :]
        B = Source.Bd_DQ[ns, :, :]
        C = Source.Cd_DQ[ns, :, :]
        D = Source.Dd_DQ[ns, :, :]

        K = Source.Ko_DQ[ns, :, :]

        uₚ = [yₚ_DQ0[1:2]; vₚ_DQ0[1:2]; eₚ_DQ0[1:2]]

        Source.xp_DQ[ns, :] = (A - K*C)*Source.xp_DQ[ns, :] + (B - K*D)*uₚ + K*y_DQ0[1:2]

        I_poc_DQ0[1:2] = Source.xp_DQ[ns, 1:2]
        V_cap_DQ0[1:2] = Source.xp_DQ[ns, 3:4]

        #----------------------------------------------------------------------

        Source.I_filt_poc[ns, :, end] = InvDQ0Transform(I_poc_DQ0, θ + 0.5*Source.ts*ω)
        Source.V_filt_cap[ns, :, end] = InvDQ0Transform(V_cap_DQ0, θ + 0.5*Source.ts*ω)

        if any(isnan.(Source.I_filt_poc[ns, :, end]))

            Source.I_filt_poc[ns, :, end] = fill!(Source.I_filt_poc[ns, :, end], 0.0)
            Source.xp_0[ns, 1] = 0.0
            Source.xp_DQ[ns, 1:2] = fill!(Source.xp_DQ[ns, 1:2], 0.0)
        end

        if any(isnan.(Source.V_filt_cap[ns, :, end]))

            Source.V_filt_cap[ns, :, end] = fill!(Source.V_filt_cap[ns, :, end], 0.0)
            Source.xp_0[ns, 2] = 0.0
            Source.xp_DQ[ns, 3:4] = fill!(Source.xp_DQ[ns, 3:4], 0.0)
        end

    end

    return nothing
end

"""
    NewtonInterpolation(coef, x_data, x)

# Description
Performs a Newton interpolation. Think of x as the point in time, and x_data as the points in time
where we know what values the function takes.
"""
function NewtonInterpolation(coef, x_data, x)

    #= Theory:
        When constructing interpolating polynomials, there is a tradeoff between
        having a better fit and having a smooth well-behaved fitting function. The
        more data points that are used in the interpolation, the higher the degree
        of the resulting polynomial, and therefore the greater oscillation it will
        exhibit between the data points. Therefore, a high-degree interpolation may be
        a poor predictor of the function between points, although the accuracy at the
        data points will be "perfect".
    =#

    n = length(x_data)
    p = coef[n]

    for k in 1:n
        p = coef[n - k + 1] + (x - x_data[n - k + 1])*p
    end

    return p
end

"""
    DividedDiff(x, y)

# Description
Calculates the coefficients required for Newton interpolation using a divided differences algorithm.
"""
function DividedDiff(x, y)

    n = length(y)

    coef = zeros(n, n)

    coef[:, 1] .= y

    for j in 2:n
        for i in 1:(n - j + 1)
            coef[i, j] = (coef[i + 1, j - 1] - coef[i, j - 1]) / (x[i + j - 1] - x[i])
        end
    end

    #= ind = findall(x->x==0, vec(coef))
    coef = deleteat!(vec(coef), ind) =#

    return vec(coef[1,:])
end

"""
    OrnsteinUhlenbeck(Source::ClassicalControls)

# Description
Produces an Ornstein Uhlenbeck process.
"""
function OrnsteinUhlenbeck(Source::ClassicalControls)

    if Source.steps*Source.ts >= Source.process_start

        for ns in 1:Source.num_sources

            κ = Source.κ[ns] # mean reversion parameter
            γ = Source.γ[ns] # asymptotoic mean
            σ = Source.σ[ns] # Brownian motion scale i.e. ∝ diffusion parameter

            if σ != 0

                Source.cnt[ns] += 1
                Δt = Source.rol[ns]*Source.ts

                if Source.cnt[ns] == Source.rol[ns]

                    Source.cnt[ns] = 0

                    if Source.k[ns] > 0
                        Source.X[ns][1:end-1] = Source.X[ns][2:end]
                    end

                    # Euler Maruyama
                    #Source.X[ns][end] = Source.X[ns][end] + κ*(γ .- Source.X[ns][end])*Δt + σ*sqrt(Δt)*randn()
                    #std_asy = sqrt(σ^2/(2*κ)) # asymptotic standard deviation

                    std_dt = sqrt(σ^2/(2*κ) * (1 - exp(-2*κ*Δt)))
                    Source.X[ns][end] = γ .+ exp(-κ*Δt)*(Source.X[ns][end] .- γ) + std_dt*randn()

                    if Source.k[ns] > 0
                        t_data = Δt*collect(0:1:Source.k[ns])
                        Source.c_diff[ns] = DividedDiff(t_data, Source.X[ns])
                    end
                end

                if Source.k[ns] > 0
                    t_data = Δt*collect(0:1:Source.k[ns])
                    Pset = NewtonInterpolation(Source.c_diff[ns], t_data, (Source.k[ns] - 1)*Δt + Source.ts*Source.cnt[ns])
                else
                    Pset = Source.X[ns][end]
                end

                Sset = Pset/Source.pf[ns]
                Source.pq0_set[ns, 1] = Pset
                Source.pq0_set[ns, 2] = sqrt(Sset^2 - Pset^2)*sign(Pset*Source.pf[ns])

                if Sset > Source.S[ns]

                    pq0_ref = Source.pq0_set[ns, :]
                    pq0_ref = pq0_ref.*(Source.S[ns]/norm(pq0_ref))

                    if Source.X[ns][end] > abs(Source.pq0_set[ns, 1])
                        Source.X[ns][end] = Source.pq0_set[ns, 1]
                    end
                end

            end

        end
    end

    return nothing
end

"""
    Measurements(Source::ClassicalControls)

# Description
Calculates RMS quantities and Active and Reactive power for every classical source.
"""
function Measurements(Source::ClassicalControls)

    #= i = Source.steps

    i_sp_rms = convert(Int64, round((1/(Source.ts*Source.T_sp_rms))))

    i_start = i - convert(Int64, round(Source.T_eval/(Source.fsys*Source.ts)))

    i_range = 1:i_sp_rms:Source.N

    for ns in 1:Source.num_sources

        V_poc = Source.V_filt_cap[ns, :, end]
        I_poc = Source.I_filt_poc[ns, :, end]

        Source.p_q_inst[ns, :] = pqTheory(V_poc, I_poc) # real and imaginary powers

        Source.p_inst[ns, :] = V_poc.*I_poc

        # Phasors
        if i_start >= 1 #waiting for one evaluation cycle to pass

            θ = Source.θ_source[ns, 1, i_range]

            # Voltages
            v_signals = transpose(Source.V_filt_cap[ns, :, i_range])
            Source.V_ph[ns,  :, :] = RMS(θ, v_signals)
            #Source.V_ph[ns,  :, 3] = Source.V_ph[ns,  :, 3] .+ π/2

            # Currents
            i_signals = transpose(Source.I_filt_poc[ns, :, i_range])
            Source.I_ph[ns, :, :] = RMS(θ, i_signals)
            #Source.I_ph[ns, :, 3] = Source.I_ph[ns, :, 3] .+ π/2

            # Per phase (and total) Active and Reactiv Powers
            Source.Pm[ns, 1:3] = (Source.V_ph[ns, :, 2].*Source.I_ph[ns, :, 2]).*cos.(Source.V_ph[ns, :, 3] .- Source.I_ph[ns, :, 3])
            Source.Pm[ns, 4] = sum(Source.Pm[ns, 1:3])
            Source.Qm[ns, 1:3] = (Source.V_ph[ns, :, 2].*Source.I_ph[ns, :, 2]).*sin.(Source.V_ph[ns, :, 3] .- Source.I_ph[ns, :, 3])
            Source.Qm[ns, 4] = sum(Source.Qm[ns, 1:3])
        end
    end =#

    # global metrics

    roll(Source.f_avg)
    roll(Source.θ_avg)

    grid_swing = findfirst(x->x == "Swing", Source.Source_Modes)
    grid_voltage = findfirst(x->x == "Voltage", Source.Source_Modes)

    if isnothing(grid_swing) &&  isnothing(grid_voltage)
        Source.f_avg[:, end] = sum(Source.f_source[Source.grid_forming , : , end], dims = 1)./length(Source.grid_forming)
    else
        Source.f_avg[:, end] = Source.fsys*[1 1 1]
    end

    for i in 1:Source.phases
        Source.θ_avg[i, end] = ThirdOrderIntegrator(Source.θ_avg[i, end], Source.ts, 2π*Source.f_avg[i, :]) # integration
        Source.θ_avg[i, end] = (Source.θ_avg[i, end] - (i - 1)*120π/180)%(2π)
    end

    return nothing
end

"""
    CurrentPILoopShaping(Source::ClassicalControls, num_source)

# Description
Tuning of proportional and integral gain for inner current controller.
"""
function CurrentPILoopShaping(Source::ClassicalControls, num_source)

    #= Theory:
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

    pole_flag = 0

    #--------------------------------------
    # Current Controller

    # Short Circuit State Space
    A = [-Source.Rf_L1[num_source]/Source.Lf_1[num_source]]
    B = [1/(Source.Lf_1[num_source])]
    C = [1]
    D = [0]
    sys = ss(A, B, C, D) # continuous
    tf_sys = tf(sys) |> ss

    Ts = Source.ts
    dly = Source.action_delay*Ts
    ZoH = tf([1], [Ts/2, 1]) |> ss # Transfer function of the sample and hold process
    Pade = tf([-dly/2, 1], [dly/2, 1]) |> ss # Pure first-order delay approximation
    PWM_gain = Source.Vdc[num_source]/2

    #SC = tf([1/(1*Lf_1)], [1, Rf_L1/Lf_1]) # = tf(sys_sc)
    Gsc_ol = minreal(tf_sys*Pade*PWM_gain*ZoH) |> ss # Full transfer function of plant

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

        Gi_cl = minreal(Gsc_ol*Gpi_i/(1 + Gsc_ol*Gpi_i)) |> ss # closed loop transfer function

        if any(real(ControlSystemsBase.poles(Gi_cl)) .> 0) == false

            # all the poles are on the left side
            ωp = 2π/((i)*Ts) # gain cross-over frequency
            pm = 60 # degrees, phase margin
            Gpi_i, kp_i, ki_i = loopshapingPI(Gsc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
            Gi_cl = minreal(Gsc_ol*Gpi_i/(1 + Gsc_ol*Gpi_i)) |> ss

            Source.I_kp[num_source] = kp_i
            Source.I_ki[num_source] = ki_i
            Source.Gi_cl[num_source] = Gi_cl

            break

        end

        if i == max_i

            pole_flag = 1

            ωp = 2π/((i)*Ts) # gain cross-over frequency
            pm = 60 # degrees, phase margin
            Gpi_i, kp_i, ki_i = loopshapingPI(Gsc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
            Gi_cl = minreal(Gsc_ol*Gpi_i/(1 + Gsc_ol*Gpi_i))

            Source.I_kp[num_source] = kp_i
            Source.I_ki[num_source] = ki_i
            Source.Gi_cl[num_source] = Gi_cl
        end
    end

    return pole_flag
end

"""
    VoltagePILoopShaping(Source::ClassicalControls, num_source)

# Description
Tuning of proportional and integral gain for outer voltage controller.
"""
function VoltagePILoopShaping(Source::ClassicalControls, num_source)

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

    pole_flag = 0

    Goc_ol = minreal(Source.Gi_cl[num_source]*tf([1], [Source.Cf[num_source], 0])) |> ss

    Ts = Source.ts
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

        Gv_cl = minreal(Goc_ol*Gpi_v/(1 + Goc_ol*Gpi_v)) |> ss # closed loop transfer function

        if any(real(ControlSystemsBase.poles(Gv_cl)) .> 0) == false

            # all the poles are on the left side
            ωp = 2π/((i)*Ts) # gain cross-over frequency
            pm = 60 # degrees, phase margin
            Gpi_v, kp_v, ki_v = loopshapingPI(Goc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
            Gv_cl = minreal(Goc_ol*Gpi_v/(1 + Goc_ol*Gpi_v)) |> ss

            Source.V_kp[num_source] = kp_v
            Source.V_ki[num_source] = ki_v
            Source.Gv_cl[num_source] = Gv_cl

            break
        end

        if i == 1

            ωp = 2π/((i)*Ts) # gain cross-over frequency
            pm = 60 # degrees, phase margin
            Gpi_v, kp_v, ki_v = loopshapingPI(Goc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
            Gv_cl = minreal(Goc_ol*Gpi_v/(1 + Goc_ol*Gpi_v)) |> ss

            Source.V_kp[num_source] = kp_v
            Source.V_ki[num_source] = ki_v
            Source.Gv_cl[num_source] = Gv_cl

            pole_flag = 1
        end
    end

    return pole_flag
end

"""
    SourceInitialiser(env, Source, modes, source_indices)

# Description
Initialises all of the sources with droop coefficients and proportional and integral gains.
"""
function SourceInitialiser(env, Source, modes, source_indices; seed = false)

    if seed
        Random.seed!(1)
    end

    # logging
    count_V_K = 0
    count_I_K = 0
    count_Dp = 0
    count_Dq = 0
    count_L_fltr = 0
    count_I_poles = 0
    count_V_poles = 0

    Mode_Keys = [k[1] for k in sort(collect(Source.Modes), by = x -> x[2])]

    Source.ramp_end = env.nc.parameters["grid"]["ramp_end"]
    Source.process_start = env.nc.parameters["grid"]["process_start"]
    Source.Δfmax = env.nc.parameters["grid"]["Δfmax"]/100
    Source.ΔEmax = env.nc.parameters["grid"]["ΔEmax"]/100

    ωsys = Source.fsys*2π

    Source.ω_sync = fill!(Source.ω_sync, ωsys)
    Source.ω_droop = fill!(Source.ω_droop, ωsys)

    e = 1

    for ns in source_indices

        Source.S[e] = env.nc.parameters["source"][ns]["pwr"]
        Source.pf[e] = env.nc.parameters["source"][ns]["pf"]

        Source.filter_type[e] = env.nc.parameters["source"][ns]["fltr"]
        Source.Rf_L1[e] = env.nc.parameters["source"][ns]["R1"]
        Source.Lf_1[e] = env.nc.parameters["source"][ns]["L1"]

        Source.Vdc[e] = env.nc.parameters["source"][ns]["vdc"]
        Source.i_max[e] = env.nc.parameters["source"][ns]["i_limit"]
        Source.v_max[e] = env.nc.parameters["source"][ns]["v_limit"]/2

        Source.Vrms[e] = env.nc.parameters["grid"]["v_rms"]
        Source.V_pu_set[e, 1] = env.nc.parameters["source"][ns]["v_pu_set"]
        Source.V_δ_set[e, 1] = env.nc.parameters["source"][ns]["v_δ_set"]*π/180

        Source.τv[e] = env.nc.parameters["source"][ns]["τv"] # time constant of the voltage loop
        Source.τf[e] = env.nc.parameters["source"][ns]["τf"] # time constant of the frequency loop

        Source.pq0_set[e, 1] = env.nc.parameters["source"][ns]["p_set"] # W, Real Power
        Source.pq0_set[e, 2] = env.nc.parameters["source"][ns]["q_set"] # VAi, Imaginary Power

        if Source.filter_type[e] == "LCL"
            Source.Rf_L2[e] = env.nc.parameters["source"][ns]["R2"]
            Source.Lf_2[e] = env.nc.parameters["source"][ns]["L2"]
        end

        if haskey(env.nc.parameters["source"][ns], "Observer") && Source.filter_type[e] == "LCL"
            Source.Observer[e] = env.nc.parameters["source"][ns]["Observer"]
        else
            Source.Observer[e] = false
        end

        Source.P[e] = Source.pf[e]*Source.S[e]
        Source.Q[e] = sqrt(Source.S[e]^2 - Source.P[e]^2)

        if typeof(modes[e]) == Int64
            Source.Source_Modes[e] = Mode_Keys[modes[e]]
        else
            if modes[e] == "VSG"
                Source.Source_Modes[e] = "Synchronverter"
            else
                Source.Source_Modes[e] = modes[e]
            end
        end

        if haskey(env.nc.parameters["source"][ns], "R_C")

            Source.Rf_C[e] = env.nc.parameters["source"][ns]["R_C"]
            Source.Cf[e] = env.nc.parameters["source"][ns]["C"]
        else

            Source.Cf[e] = Source.S[e]/(3*Source.fsys*2π*Source.Vrms[e]*Source.Vrms[e])
            Source.Rf_C[e] = Source.Cf[e]

            if (Source.Source_Modes[e] != "Swing"
                && Source.Source_Modes[e] != "PQ"
                && Source.Source_Modes[e] != "PV"
                && Source.Source_Modes[e] != "Step")
                count_L_fltr += 1
            end
        end

        if (Source.Source_Modes[e] == "Semi-Synchronverter" ||
            Source.Source_Modes[e] == "Synchronverter" ||
            Source.Source_Modes[e] == "Semi-Droop" ||
            Source.Source_Modes[e] == "Droop")

            if !haskey(env.nc.parameters["source"][ns], "Dp")
                Source.D[e, 1] = Source.S[e]/(ωsys*ωsys*Source.Δfmax)
                count_Dp += 1
            else
                Source.D[e, 1] = env.nc.parameters["source"][ns]["Dp"]
            end

            if !haskey(env.nc.parameters["source"][ns], "Dq")
                if Source.Source_Modes[e] == "Semi-Droop"
                    Source.D[e, 2] = 0.0
                else
                    Source.D[e, 2] = Source.S[e]/(Source.V_pu_set[e, 1]*Source.Vrms[e]*sqrt(2)*Source.ΔEmax)
                end
                count_Dq += 1
            else
                Source.D[e, 2] = env.nc.parameters["source"][ns]["Dq"]
            end
        end

        Source.J_sync[e] = Source.τf[e]*Source.D[e, 1] #total mass moment of inertia of the rotating masses, kg*m^2
        Source.K_sync[e] = Source.τv[e]*ωsys*Source.D[e, 2]

        if Source.Source_Modes[e] == "Synchronverter" || Source.Source_Modes[e] == "PQ"
            Source.σ[e] = abs(env.nc.parameters["source"][ns]["σ"]) # Brownian motion scale (standard deviation) - sqrt(diffusion)
        else
            Source.σ[e] = 0.0
        end

        Source.κ[e] = abs(env.nc.parameters["source"][ns]["κ"]) # mean reversion parameter
        Source.γ[e] = env.nc.parameters["source"][ns]["γ"] # asymptotic mean
        Source.k[e] = env.nc.parameters["source"][ns]["k"] # interpolation degree

        if Source.k[e] - 1 > env.nc.parameters["source"][ns]["Δt"]/Source.ts
            Source.k[e] = floor(env.nc.parameters["source"][ns]["Δt"]/ts) - 1
        end

        n = convert(Int, round(Source.k[e] + 1))
        X₀ = env.nc.parameters["source"][ns]["X₀"]*ones(n)

        Source.X[e] = cat(X₀, dims = 1)
        Source.X₀[e] = Source.X[e][end]

        Source.rol[e] = convert(Int64, round(env.nc.parameters["source"][ns]["Δt"]*env.nc.parameters["grid"]["fs"]))

        if n > 0

            Δt = Source.rol[e]*Source.ts
            t_data = Δt*collect(0:1:n)
            coef = DividedDiff(t_data, Source.X[e])
            Source.c_diff[e] = cat(coef, dims = 1)
        end

        if Source.Source_Modes[e] != "Swing" && Source.Source_Modes[e] != "Step"
            if !haskey(env.nc.parameters["source"][ns], "I_kp") && !haskey(env.nc.parameters["source"][ns], "I_ki")

                count_I_poles += CurrentPILoopShaping(Source, e)
                count_I_K += 1
            else

                Source.I_kp[e] = env.nc.parameters["source"][ns]["I_kp"]
                Source.I_ki[e] = env.nc.parameters["source"][ns]["I_ki"]
            end
        end

        if (Source.Source_Modes[e] != "Swing"
            && Source.Source_Modes[e] != "Step"
            && Source.Source_Modes[e] != "PQ"
            && Source.Source_Modes[e] != "PV")

            if !haskey(env.nc.parameters["source"][ns], "V_kp") && !haskey(env.nc.parameters["source"][ns], "V_ki")

                count_V_poles += VoltagePILoopShaping(Source, e)
                count_V_K += 1
            else

                Source.V_kp[e] = env.nc.parameters["source"][ns]["V_kp"]
                Source.V_ki[e] = env.nc.parameters["source"][ns]["V_ki"]
            end
        end

        if Source.Observer[e]

            ObserverInitialiser(Source, e)

            if isnan(Source.Ko_DQ[e, 1, 1])
                Source.Observer[e] = false
            end
        end

        e += 1
    end

    Source.grid_forming = union(findall(x->x == "Synchronverter", Source.Source_Modes),
                        findall(x->x == "Droop", Source.Source_Modes),
                        findall(x->x == "Swing", Source.Source_Modes),
                        findall(x->x == "Voltage", Source.Source_Modes),
                        findall(x->x == "Semi-Synchronverter", Source.Source_Modes),
                        findall(x->x == "Semi-Droop", Source.Source_Modes),
                        )

    Source.grid_following = union(findall(x->x == "PQ", Source.Source_Modes),
                        findall(x->x == "PV", Source.Source_Modes),
                        )

    #-------------------------------------------------------------------------------
    # Logging

    if env.verbosity > 0

        if count_L_fltr == 1
            @warn "$(count_L_fltr) source with an 'L' filter is being controlled. 'LCL' or 'LC' filters are preferred for grid forming sources."
        elseif count_L_fltr > 1
            @warn "$(count_L_fltr) source 'L' filters are being controlled. 'LCL' or 'LC' filters are preferred for grid forming sources."
        end

        if count_V_poles == 1
            @warn "$(count_V_poles) Voltage Controller with Positive Poles.
            Suggestion: Decrease simulation time step or choose different filter values."
        elseif count_L_fltr > 1
            @warn "$(count_V_poles) Voltage Controllers with Positive Poles.
            Suggestion: Decrease simulation time step or choose different filter values."
        end

        if count_I_poles == 1
            @warn "$(count_I_poles) Current Controller with Positive Poles.
            Suggestion: Decrease simulation time step or choose different filter values."
        elseif count_I_poles > 1
            @warn "$(count_I_poles) Current Controllers with Positive Poles.
            Suggestion: Decrease simulation time step or choose different filter values."
        end

    end

    mode_count = Array{Int64, 1}(undef, length(Mode_Keys))
    mode_count = fill!(mode_count, 0)

    if env.verbosity > 1

        if Source.num_sources == 1
            @info "$(Source.num_sources) 'classically' controlled source has been initialised."
        else
            @info "$(Source.num_sources) 'classically' controlled sources have been initialised."
        end

        for modes in eachindex(Mode_Keys)

            mode_count[modes] = length(findall(x->x == Mode_Keys[modes], Source.Source_Modes))

            if mode_count[modes] != 0
                if mode_count[modes] == 1
                    @info "$(mode_count[modes]) source has been set up in $(Mode_Keys[modes]) mode."
                else
                    @info "$(mode_count[modes]) sources have been set up in $(Mode_Keys[modes]) mode."
                end
            end

        end

        if (count_V_K == sum(mode_count[3:6]) && count_I_K == sum(mode_count[2:6])
            && count_Dp == mode_count[3] + mode_count[4] + mode_count[5] + mode_count[6]
            && count_Dq == mode_count[3] + mode_count[4] + mode_count[5] + mode_count[6])

            @info "All 'classically' controlled sources have been automatically set up with droop coeficients, and proportional and integral gains."

        else

            if count_V_K == 1
                @info "$(count_V_K) source has automatically calculated proportional and integral gains for its voltage control loop."
            elseif count_V_K > 1
                @info "$(count_V_K) sources have automatically calculated proportional and integral gains for their voltage control loops."
            end

            if count_I_K == 1
                @info "$(count_I_K) source has automatically calculated proportional and integral gains for its current control loop."
            elseif count_I_K > 1
                @info "$(count_I_K) sources have automatically calculated proportional and integral gains for their current control loops."
            end

            if count_Dp == 1
                @info "$(count_Dp) source has an automatically calculated frequency droop coefficient."
            elseif count_Dp > 1
                @info "$(count_Dp) sources have automatically calculated frequency droop coefficients."
            end

            if count_Dq == 1
                @info "$(count_Dq) source has an automatically calculated voltage droop coefficient."
            elseif count_Dq > 1
                @info "$(count_Dq) sources have automatically calculated voltage droop coefficients."
            end
        end

        if length(findall(Source.Observer)) == 1
            @info "$(length(findall(Source.Observer))) source has been set up with a Luenberger discrete Observer."
        elseif length(findall(Source.Observer)) > 1
            @info "$(length(findall(Source.Observer))) sources have been set up with Luenberger discrete Observers."
        end

        processes = length(findall(x->x > 0.0, convert.(Float64, Source.σ)))
        if processes == 1
            @info "$(processes) stochastic process will start at $(Source.process_start) [s]."
        elseif processes > 1
            @info "$(processes) stochastic processes will start at $(Source.process_start) [s]."
        end
    end

    return nothing
end

"""
    ObserverInitialiser(Source::ClassicalControls, ns)

# Description
Initialises the observers
"""
function ObserverInitialiser(Source::ClassicalControls, ns)

    # Predictive Approximate Deadbeat Reduced-Order Observer

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

    if Source.filter_type[ns] == "LCL"

        ω = 2π*Source.fsys

        # Finding the phase/zero-sequence system dynamics

        A = [-(Source.Rf_L1[ns] + Source.Rf_C[ns])/Source.Lf_1[ns] Source.Rf_C[ns]/Source.Lf_1[ns] -1/Source.Lf_1[ns];
        Source.Rf_C[ns]/Source.Lf_2[ns] -(Source.Rf_L2[ns] + Source.Rf_C[ns])/Source.Lf_2[ns] 1/Source.Lf_2[ns];
        1/Source.Cf[ns] -1/Source.Cf[ns] 0.0]

        B = [1/Source.Lf_1[ns]; 0; 0]
        D = [0; -1/Source.Lf_2[ns]; 0]

        #------------------------------------------------------------------------------------------------
        # DQ sequence observer
        # Converting to discrete form

        # Converting to DQ frame - under balanced symmetrical conditions the zero sequence variables can be removed
        # if zero sequence is non-trivial then because the components are independent of DQ they can be observed independently

        A_DQ = [A ω*I;
                -ω*I A]

        B_DQ = [B zeros(3,1);
               zeros(3,1) B]

        D_DQ = [D zeros(3,1);
               zeros(3,1) D]

        # reorganising matrices
        A_DQ = SwitchRows!(A_DQ, 2, 4)
        A_DQ = SwitchRows!(A_DQ, 3, 4)
        A_DQ = SwitchRows!(A_DQ, 4, 5)

        B_DQ = SwitchRows!(B_DQ, 2, 4)

        D_DQ = SwitchRows!(D_DQ, 2, 3)
        D_DQ = SwitchRows!(D_DQ, 4, 5)

        Ad = exp(A_DQ*Source.ts)
        Bd = A_DQ \ (Ad - I)*B_DQ
        Dd = A_DQ \ (Ad - I)*D_DQ

        Source.Ad_DQ[ns, :, :] = Ad[3:end, 3:end]
        Source.Bd_DQ[ns, :, :] = [Ad[3:end, 1:2] Bd[3:end, :] Dd[3:end, :]]
        Source.Dd_DQ[ns, :, :] = [Ad[1:2, 1:2] Bd[1:2, 1:2] Dd[1:2, 1:2]]

        Source.Cd_DQ[ns, :, :] = Ad[1:2, 3:end] #Source.Cd_DQ[ns, :, :] = [1 0 0 0 0 0; 0 0 0 1 0 0]

        _, r = Observability(Source.Cd_DQ[ns, :, :], Source.Ad_DQ[ns, :, :])

        if r != size(Source.Ad_DQ[ns, :, :], 1)
            @warn ("The DQ sequence states of the inverter are not observable.
            The rank of 'O' is not equal to $(size(Source.Ad_DQ[ns, :, :], 1)).")
        end

        #------------------------------------------------------------------------------------------------
        # Solving (A - K*C)^4 = [0]

        λ = [0.003; 0.003; 0.0; 0.0]

        p = [2.0 1.0 1.0 1.5;
             0.0 2.0 0.5 1.0]

        Source.Ko_DQ[ns, :, :],   = MultiGainMatrixPar(Source.Ad_DQ[ns, :, :], Source.Cd_DQ[ns, :, :], λ, p)

        if !isnan(Source.Ko_DQ[ns, 1, 1])
            λₒ = round.(eigvals(Source.Ad_DQ[ns, :, :] - Source.Ko_DQ[ns, :, :]*Source.Cd_DQ[ns, :, :]), digits = 3)

            if maximum(abs.(sort(λ) .- λₒ)) != 0
                Source.Ko_DQ[ns, :, :] = fill!(Source.Ko_DQ[ns, :, :], NaN)
            end
        end

        #(Source.Ad_DQ[ns, :, :] - Source.Ko_DQ[ns, :, :]*Source.Cd_DQ[ns, :, :])^4

        #------------------------------------------------------------------------------------------------
        # zero sequence observer

        # Converting to discrete form

        Ad = exp(A*Source.ts)
        Bd = inv(A)*(Ad - I)*B
        Dd = inv(A)*(Ad - I)*D

        Source.Ad_0[ns, :, :] = Ad[2:3, 2:3]
        Source.Bd_0[ns, :, :] = [Ad[2:3, 1] Bd[2:3] Dd[2:3]]
        Source.Dd_0[ns, :] = [Ad[1, 1]; Bd[1]; Dd[1]]

        Source.Cd_0[ns, :] = Ad[1, 2:3]

        C = transpose(Source.Cd_0[ns, :])
        _, r = Observability(C, Source.Ad_0[ns, :, :])

        if r != size(Source.Ad_0[ns, :, :], 1)
            @warn ("The 0 sequence states of the inverter are not observable.
            The rank of 'O' is not equal to $(size(Source.Ad_0[ns, :, :], 1)).")
        end

        Source.Ko_0[ns, :] = AckermannGainMatrix(Source.Ad_0[ns, :, :], C, [0; 0])

        #------------------------------------------------------------------------------------------------

    end

    return nothing
end

"""
    O, rank(O) = Observability(C, A)

# Description
Finds the observability matrix.

# Return Values
- `O`: Observability matrix
- `rank(O)`: rank of O
"""
function Observability(C, A)

    n = size(A,1)

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

"""
    SwitchRows!(A, row_1, row_2)

# Description
Switches rows and columns of a matrix
"""
function SwitchRows!(A, row_1, row_2)

    num_rows = size(A, 1)
    num_cols = size(A, 2)

    for i in 1:num_rows

        if i <= num_cols && row_1 <= num_rows && row_2 <= num_rows

            temp = A[row_1, i]
            A[row_1, i] = A[row_2, i]
            A[row_2, i] = temp
        end
    end

    for i in 1:num_cols

        if i <= num_rows && row_1 <= num_cols && row_2 <= num_cols

            temp = A[i, row_1]
            A[i, row_1] = A[i, row_2]
            A[i, row_2] = temp
        end
    end

    return A
end

"""
    α = CharpolyCoef(λ)

# Description
given the roots, this function finds the coefficients

# Return Values
- `α`: vector of length(λ)
"""
function CharpolyCoef(λ)

    # given the roots, this function finds the coefficients

    n = length(λ)

    λ = -1*λ

    α = Array{Float64, 1}(undef, length(λ))

    for i in 1:length(λ)

        x = combinations(1:n, n-i+1)
        y = collect(x)

        α[i] = 0.0

        for j in eachindex(y)

            α[i] = α[i] + prod(λ[y[j]])
        end
    end

    return α
end

"""
    K = AckermannGainMatrix(λ)

# Description
Finds the Ackermann Gain Matrix given the chosen eigenvalues.
"""
function AckermannGainMatrix(A, C, λ)

    #= Theory
        For a single-output, observable system (A, C) and the desired closed-loop
        eigenvalues being the roots of the characteristic polynomial αd_A the
        feedback of the constant gain matrix can be found.
    =#

    n = size(A, 1)

    αd_A = Array{Float64, 2}(undef, size(A,1), size(A,2))

    α = CharpolyCoef(λ)

    αd_A = α[1]*I

    for i in 2:n
        αd_A = αd_A + α[i]*(A^(i-1))
    end

    αd_A = αd_A .+ A^n

    O, _ = Observability(C, A)

    output = zeros(n)
    output[end] = 1.0

    K = αd_A*inv(O)*output

    return K
end

"""
    K, v = MultiGainMatrixPar(A, C, λ, p)

# Description
Finds the Gain Matrix given the chosen eigenvalues for a multi input system.
"""
function MultiGainMatrixPar(A, C, λ, p)

    n = size(A, 1)
    v = Array{Float64, 2}(undef, size(A,1), size(A,2))
    K = Array{Float64, 2}(undef, size(A,1), size(C,1))

    for i in 1:n

        v[:, i] = -transpose(p[:, i])*C*inv(I*λ[i] - A)
    end

    if round(det(v), digits = 9) != 0
        K = inv(transpose(v))*transpose(p)
    else
        K = fill!(K, NaN)
    end

    return K, v
end

function FeedGainMatrixPar(A, B, λ, p)

    n = size(A, 1)
    v = Array{Float64, 2}(undef, size(A,1), size(A,2))

    for i in 1:n

        v[:, i] = -inv(I*λ[i] - A)*B*p[:, i]
    end

    F = p*inv(v)

    return F, v
end
