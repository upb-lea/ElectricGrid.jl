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
    N::Int64
    steps::Int64

    V_cable_loc::Matrix{Int64} # the position in the state vector where the POC Voltage is measured
    V_cap_loc::Matrix{Int64}
    I_poc_loc::Matrix{Int64}
    I_inv_loc::Matrix{Int64}

    Action_loc::Vector{Vector{Int64}}

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
    p_q_filt::Matrix{Float64}

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
    Vd_abc_new::Matrix{Float64} # output action

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

    Ad::Array{Float64}
    Bd::Matrix{Float64}
    Cd::Matrix{Float64}
    Dd::Matrix{Float64}
    Ko::Matrix{Float64}

    function Classical_Controls(Vdc::Vector{Float64}, Vrms::Vector{Float64},
        S::Vector{Float64}, P::Vector{Float64}, Q::Vector{Float64},
        i_max::Vector{Float64}, v_max::Vector{Float64}, filter_type::Vector{String},
        Lf_1::Vector{Float64}, Lf_2::Vector{Float64}, Cf::Vector{Float64}, 
        Rf_L1::Vector{Float64}, Rf_L2::Vector{Float64}, Rf_C::Vector{Float64},
        T_eval::Int64, T_sp_rms::Float64, V_ph::Array{Float64}, I_ph::Array{Float64}, 
        p_q_inst::Matrix{Float64}, p_inst::Matrix{Float64}, Pm::Matrix{Float64}, Qm::Matrix{Float64},
        debug::Vector{Float64}, Modes::Dict{String, Int64}, Source_Modes::Vector{String},
        num_sources::Int64, phases::Int64,
        f_cntr::Float64, fsys::Float64, θsys::Float64,
        ts::Float64, N::Int64, steps::Int64,
        V_cable_loc::Matrix{Int64}, V_cap_loc::Matrix{Int64}, 
        I_poc_loc::Matrix{Int64}, I_inv_loc::Matrix{Int64},
        Action_loc::Vector{Vector{Int64}},
        pll_err::Array{Float64}, pll_err_t::Matrix{Float64},
        vd::Array{Float64}, qvd::Array{Float64},
        fpll::Array{Float64}, θpll::Array{Float64},
        V_filt_poc::Array{Float64}, V_filt_cap::Array{Float64}, 
        I_filt_poc::Array{Float64}, I_filt_inv::Array{Float64},
        p_q_filt::Matrix{Float64},
        Gi_cl::Array{TransferFunction},
        I_dq0::Matrix{Float64}, I_ref_dq0::Matrix{Float64}, I_ref::Matrix{Float64},
        I_err::Array{Float64}, I_err_t::Matrix{Float64},
        I_kp::Vector{Float64}, I_ki::Vector{Float64},
        s_dq0_avg::Matrix{Float64},Vd_abc_new::Matrix{Float64}, s_lim::Matrix{Float64},
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
        Ad::Array{Float64}, Bd::Matrix{Float64}, Cd::Matrix{Float64}, Dd::Matrix{Float64}, 
        Ko::Matrix{Float64})

        new(Vdc, Vrms,
        S, P, Q,
        i_max, v_max, filter_type,
        Lf_1, Lf_2, Cf, 
        Rf_L1, Rf_L2, Rf_C,
        T_eval, T_sp_rms, V_ph, I_ph, 
        p_q_inst, p_inst, Pm, Qm,
        debug, Modes, Source_Modes,
        num_sources, phases,
        f_cntr, fsys, θsys,
        ts, N, steps,
        V_cable_loc, V_cap_loc,
        I_poc_loc, I_inv_loc,
        Action_loc,
        pll_err, pll_err_t,
        vd, qvd,
        fpll, θpll,
        V_filt_poc, V_filt_cap, 
        I_filt_poc, I_filt_inv,
        p_q_filt,
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
        Ad, Bd, Cd, Dd,
        Ko)
    end

    function Classical_Controls(f_cntr, num_sources; phases = 3)

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
        fsys = 50.0
        θsys = 0.0

        T_eval = 1 #number of periods to average over (for rms calcs)
        N = convert(Int64, round(T_eval/(fsys*ts))) + 1

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

        debug = Array{Float64, 1}(undef, 6) # put anything in here that needs to be debugged - for plotting
        debug = fill!(debug, 0)

        #---------------------------------------------------------------------------
        # General System & Control

        Modes = Dict("Swing" => 1, 
        "Voltage Control" => 2, 
        "PQ Control" => 3,
        "PV Control" => 4,
        "Droop Control" => 5, 
        "Full-Synchronverter" => 6, 
        "Semi-Synchronverter" => 7,
        "Not Used 1" => 8,
        "Not Used 2" => 9,
        "Not Used 3" => 10)

        Source_Modes = Array{String, 1}(undef, num_sources)
        Source_Modes = fill!(Source_Modes, "Voltage Control")

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
        I_filt_poc = fill!(I_filt_poc, 0)
        I_filt_inv = Array{Float64, 3}(undef, num_sources, phases, N)
        I_filt_inv = fill!(I_filt_inv, 0)

        p_q_filt = Array{Float64, 2}(undef, num_sources, phases)
        p_q_filt  = fill!(p_q_filt, 0)

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

        Vd_abc_new = Array{Float64, 2}(undef, num_sources, phases)
        Vd_abc_new = fill!(Vd_abc_new, 0)

        s_lim = Array{Float64, 2}(undef, num_sources, phases)
        s_lim = fill!(s_lim, 0)

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
        θ_droop = fill!(θ_droop, 0)

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
        θ_sync = fill!(θ_sync, 0)
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

        Ad = Array{Float64, 3}(undef, num_sources, 3, 3)
        Bd = Array{Float64, 2}(undef, num_sources, 3)
        Cd = Array{Float64, 2}(undef, num_sources, 3)
        Dd = Array{Float64, 2}(undef, num_sources, 3)
        Ko = Array{Float64, 2}(undef, num_sources, 2)

        Classical_Controls(Vdc, Vrms,
        S, P, Q,
        i_max, v_max, filter_type,
        Lf_1, Lf_2, Cf, 
        Rf_L1, Rf_L2, Rf_C,
        T_eval, T_sp_rms, V_ph, I_ph, 
        p_q_inst, p_inst, Pm, Qm,
        debug, Modes, Source_Modes,
        num_sources, phases,
        f_cntr, fsys, θsys,
        ts, N, steps,
        V_cable_loc, V_cap_loc,
        I_poc_loc, I_inv_loc,
        Action_loc,
        pll_err, pll_err_t,
        vd, qvd,
        fpll, θpll,
        V_filt_poc, V_filt_cap,
        I_filt_poc, I_filt_inv,
        p_q_filt,
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
        Ad, Bd, Cd, Dd, 
        Ko)
    end
end

Base.@kwdef mutable struct Classical_Policy <: AbstractPolicy

    action_space::Space{Vector{ClosedInterval{Float64}}}
    Source::Classical_Controls

    state_ids::Vector{String}
    action_ids::Vector{String}
    Source_Indices::Vector{Int64}

    function Classical_Policy(action_space, Source, state_ids, action_ids, Source_Indices)
        new(action_space, Source, state_ids, action_ids, Source_Indices)
    end

    function Classical_Policy(env)

        Source_Indices = Array{Int64, 1}(undef, 0)
        Modes = Array{Any, 1}(undef, 0)

        for ns in 1:env.nc.num_sources

            if env.nc.parameters["source"][ns]["control_type"] == "classic"
                Source_Indices = [Source_Indices; Int(ns)]
                Modes = [Modes; env.nc.parameters["source"][ns]["mode"]]
            end
        end

        state_ids = get_state_ids(env.nc)
        action_ids = get_action_ids(env.nc)

        ssa = "source".*string.(Source_Indices)
        state_ids_classic = filter(x -> !isempty(findall(y -> y == split(x, "_")[1], ssa)), state_ids)     
        action_ids_classic = filter(x -> !isempty(findall(y -> y == split(x, "_")[1], ssa)), action_ids)

        Source = Classical_Controls(1/env.ts, length(Source_Indices), phases = env.nc.parameters["grid"]["phase"])

        Source_Initialiser(env, Source, Modes, Source_Indices)

        #------------------------------------

        for s in axes(Source_Indices, 1)

            s_idx = string(Source_Indices[s])
    
            Source.V_cable_loc[:, s]  = findall(contains(s_idx*"_v_C_cable"), state_ids_classic)
            Source.I_inv_loc[:, s] = findall(contains(s_idx*"_i_L1"), state_ids_classic)

            if Source.filter_type[s] == "LC"

                Source.V_cap_loc[:, s]  = findall(contains(s_idx*"_v_C_filt"), state_ids_classic)
                Source.I_poc_loc[:, s] = findall(contains(s_idx*"_v_C_filt"), state_ids)

            elseif Source.filter_type[s] == "LCL"

                Source.V_cap_loc[:, s]  = findall(contains(s_idx*"_v_C_filt"), state_ids_classic)
                Source.I_poc_loc[:, s] = findall(contains(s_idx*"_i_L2"), state_ids_classic)
            
            elseif Source.filter_type[s] == "L"

                Source.I_poc_loc[:, s] = findall(contains(s_idx*"_v_C_cable"), state_ids)
            end
        end

        letterdict = Dict("a" => 1, "b" => 2, "c" => 3)

        Source.Action_loc = [[findfirst(y -> y == parse(Int64, SubString(split(x, "_")[1], 7)), 
        Source_Indices), letterdict[split(x, "_")[3]]] for x in action_ids_classic]

        #------------------------------------

        animo = Classical_Policy(Space([-1.0..1.0 for i in 1:length(action_ids_classic)]), Source,
        state_ids_classic, action_ids_classic, Source_Indices)

        return animo
    end
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

    Source.steps = 0

    Source.vd = fill!(Source.vd, 0.0)
    Source.qvd = fill!(Source.qvd, 0.0)

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

    Source.ω_droop = fill!(Source.ω_droop, Source.fsys*2π)

    Source.θ_droop = fill!(Source.θ_droop, 0)

    Source.α_sync = fill!(Source.α_sync, 0)
    Source.ω_sync = fill!(Source.ω_sync, Source.fsys*2π)
    Source.θ_sync = fill!(Source.θ_sync, 0)
    Source.Δω_sync = fill!(Source.Δω_sync, 0)
    Source.eq = fill!(Source.eq, 0)
    Source.Mfif = fill!(Source.Mfif, 0)

    Source.ΔT_err = fill!(Source.ΔT_err, 0)
    Source.ΔT_err_t = fill!(Source.ΔT_err_t, 0)
    Source.ω_set = fill!(Source.ω_set, 0)

    # fill these with initial conditions
    Source.V_filt_poc = fill!(Source.V_filt_poc, 0)
    Source.V_filt_cap = fill!(Source.V_filt_cap, 0)

    Source.I_filt_poc = fill!(Source.I_filt_poc, 0)
    Source.I_filt_inv = fill!(Source.I_filt_inv, 0)

    return nothing
end

function Classical_Control(Animo, env, name = nothing)
    
    Source = Animo.Source
    Source_Interface(env, Source, name)

    ramp_end = 2/Source.fsys

    for s in 1:Source.num_sources

        if Source.Source_Modes[s] == "Swing"

            Swing_Mode(Source, s, t_end = ramp_end)
        elseif Source.Source_Modes[s] == "Voltage Control"

            Voltage_Control_Mode(Source, s, t_end = ramp_end)
        elseif Source.Source_Modes[s] == "PQ Control"

            PQ_Control_Mode(Source, s, Source.pq0_set[s, :])
        elseif Source.Source_Modes[s] == "PV Control"

            PV_Control_Mode(Source, s, Source.pq0_set[s, :])
        elseif Source.Source_Modes[s] == "Droop Control"

            Droop_Control_Mode(Source, s, t_end = ramp_end)
        elseif Source.Source_Modes[s] == "Full-Synchronverter"

            Synchronverter_Mode(Source, s, pq0_ref = Source.pq0_set[s, :], mode = 1, t_end = ramp_end)
        elseif Source.Source_Modes[s] == "Semi-Synchronverter"

            Synchronverter_Mode(Source, s, pq0_ref = Source.pq0_set[s, :], mode = 2, t_end = ramp_end)
        elseif Source.Source_Modes[s] == "Not Used 1"

            Synchronverter_Mode(Source, s, pq0_ref = Source.pq0_set[s, :], mode = 2, t_end = ramp_end)
        elseif Source.Source_Modes[s] == "Not Used 2"

            Synchronverter_Mode(Source, s, pq0_ref = Source.pq0_set[s, :], mode = 2, t_end = ramp_end)
        elseif Source.Source_Modes[s] == "Not Used 3"

            Synchronverter_Mode(Source, s, pq0_ref = Source.pq0_set[s, :], mode = 2, t_end = ramp_end)
        end
    end

    Action = Env_Interface(Source)
    Measurements(Source)

    return Action
end

function Source_Interface(env, Source::Classical_Controls, name = nothing)

    Source.steps = env.steps + 1
    ω = 2π*Source.fsys
    Source.θsys = (Source.θsys + Source.ts*ω)%(2π)

    observer = true

    if !isnothing(name)
        state = RLBase.state(env, name)
    end

    for ns in 1:Source.num_sources

        Source.V_filt_cap[ns, :, 1:end-1] = Source.V_filt_cap[ns, :, 2:end]
        Source.V_filt_poc[ns, :, 1:end-1] = Source.V_filt_poc[ns, :, 2:end]

        Source.I_filt_inv[ns, :, 1:end-1] = Source.I_filt_inv[ns, :, 2:end]  
        Source.I_filt_poc[ns, :, 1:end-1] = Source.I_filt_poc[ns, :, 2:end]

        Source.I_filt_inv[ns, :, end] = state[Source.I_inv_loc[:, ns]]

        if Source.filter_type[ns] == "LCL"

            if observer == true || 1 == 1

                Source.V_filt_poc[ns, :, end] = state[Source.V_cable_loc[:, ns]]
                Luenberger_Observer(Source, ns)
                
            else

                Source.V_filt_cap[ns, :, end] = state[Source.V_cap_loc[:, ns]]

                Source.V_filt_poc[ns, :, end] = state[Source.V_cable_loc[:, ns]]
                Source.I_filt_poc[ns, :, end] = state[Source.I_poc_loc[:, ns]]
            end

        elseif Source.filter_type[ns] == "LC"

            Source.V_filt_cap[ns, :, end] = state[Source.V_cable_loc[:, ns]]

            Source.V_filt_poc[ns, :, end] = Source.V_filt_cap[ns, :, end]
            Source.I_filt_poc[ns, :, end] = Source.I_filt_inv[ns, :, end] .- env.y[Source.I_poc_loc[:, ns]]

        elseif Source.filter_type[ns] == "L"

            Source.V_filt_cap[ns, :, end] = state[Source.V_cable_loc[:, ns]]

            Source.V_filt_poc[ns, :, end] = Source.V_filt_cap[ns, :, end] 
            Source.I_filt_poc[ns, :, end] = Source.I_filt_inv[ns, :, end] .- env.y[Source.I_poc_loc[:, ns]]
        end

        Source.p_q_filt[ns, :] =  p_q_theory((Source.Vdc[ns]/2)*Source.Vd_abc_new[ns, :], Source.I_filt_inv[ns, :, end])
    end

    return nothing
end

function Env_Interface(Source::Classical_Controls)

    Action = [Source.Vd_abc_new[Source.Action_loc[x][1], Source.Action_loc[x][2]] for x in 1:length(Source.Action_loc)]

    return Action
end

function Ramp(final, μ, i; t_end = 0.02)

    if μ*i < t_end

        x_out = final.*(μ*i/t_end)
    else
        x_out = final
    end

    return x_out
end

function Swing_Mode(Source::Classical_Controls, num_source; t_end = 0.04)
    
    δ = Source.V_δ_set[num_source, 1]
    pu = Source.V_pu_set[num_source, 1]

    θt = Source.θsys
    θph = [θt + δ; θt + δ - 120π/180; θt + δ + 120π/180]
    Vrms = Ramp(pu*Source.Vrms[num_source], Source.ts, Source.steps; t_end = t_end)
    Source.V_ref[num_source, :] = sqrt(2)*(Vrms)*cos.(θph)

    Source.Vd_abc_new[num_source, :] = 2*Source.V_ref[num_source, :]/Source.Vdc[num_source]

    Phase_Locked_Loop_3ph(Source, num_source)

    return nothing
end

function Voltage_Control_Mode(Source::Classical_Controls, num_source; t_end = 0.04)

    δ = Source.V_δ_set[num_source, 1]
    pu = Source.V_pu_set[num_source, 1]

    θt = Source.θsys
    θph = [θt + δ; θt + δ - 120π/180; θt + δ + 120π/180]
    Vrms = Ramp(pu*Source.Vrms[num_source], Source.ts, Source.steps; t_end = t_end)
    Source.V_ref[num_source, :] = sqrt(2)*Vrms*cos.(θph)
    ω = 2π*Source.fsys
    
    Voltage_Controller(Source, num_source, θt, ω)
    Current_Controller(Source, num_source, θt, ω)

    Phase_Locked_Loop_3ph(Source, num_source)

    return nothing
end

function Droop_Control_Mode(Source::Classical_Controls, num_source; t_end = 0.04)

    pu = Source.V_pu_set[num_source, 1]

    Vrms = Ramp(pu*Source.Vrms[num_source], Source.ts, Source.steps; t_end = t_end)

    Droop_Control(Source, num_source, Vrms = Vrms)
    θt = Source.θ_droop[num_source, 1]
    ω = Source.ω_droop[num_source, 1, end]

    Voltage_Controller(Source, num_source, θt, ω)
    Current_Controller(Source, num_source, θt, ω)

    Phase_Locked_Loop_3ph(Source, num_source)

    return nothing
end

function PQ_Control_Mode(Source::Classical_Controls, num_source, pq0)

    if norm(pq0) > Source.S[num_source]
        pq0 = pq0.*(Source.S[num_source]/norm(pq0))
    end

    Phase_Locked_Loop_3ph(Source, num_source)
    #Phase_Locked_Loop_1ph(Source, num_source, ph = 1)
    #Phase_Locked_Loop_1ph(Source, num_source, ph = 2)
    #Phase_Locked_Loop_1ph(Source, num_source, ph = 3)
    θ = Source.θpll[num_source, 1, end]
    ω = 2π*Source.fpll[num_source, 1, end]

    Filtering(Source, num_source, θ)

    if Source.steps*Source.ts > 4/Source.fsys

        PQ_Control(pq0_ref = pq0, Source, num_source, θ)
        Current_Controller(Source, num_source, θ, ω)
    else

        PQ_Control(pq0_ref = [0.0; 0.0; 0.0], Source, num_source, θ)
        Source.I_ref_dq0[num_source, :] = [0.0; 0.0; 0.0]
        Current_Controller(Source, num_source, θ, ω)
    end

    return nothing
end

function PV_Control_Mode(Source::Classical_Controls, num_source, pq0)

    Phase_Locked_Loop_3ph(Source, num_source)
    #Phase_Locked_Loop_1ph(Source, num_source, ph = 1)
    #Phase_Locked_Loop_1ph(Source, num_source, ph = 2)
    #Phase_Locked_Loop_1ph(Source, num_source, ph = 3)
    θ = Source.θpll[num_source, 1, end] # phase a
    ω = 2π*Source.fpll[num_source, 1, end]

    Filtering(Source, num_source, θ)

    if Source.steps*Source.ts > 4/Source.fsys

        pq0_ref = PV_Control(pq0_ref = pq0, Source, num_source)
        PQ_Control(pq0_ref = pq0_ref, Source, num_source, θ)
        Current_Controller(Source, num_source, θ, ω)
    else

        PQ_Control(pq0_ref = [0.0; 0.0; 0.0], Source, num_source, θ)
        Source.I_ref_dq0[num_source, :] = [0.0; 0.0; 0.0]
        Current_Controller(Source, num_source, θ, ω)
    end

    return nothing
end

function Synchronverter_Mode(Source::Classical_Controls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]], t_end = 0.04, mode = 2)

    pu = Source.V_pu_set[num_source, 1]

    Vrms = Ramp(pu*Source.Vrms[num_source], Source.ts, Source.steps; t_end = t_end)

    Synchronverter_Control(Source, num_source, pq0_ref = pq0_ref, Vrms = Vrms, mode = mode)

    Voltage_Controller(Source, num_source, Source.θ_sync[num_source], Source.ω_sync[num_source, end])
    Current_Controller(Source, num_source, Source.θ_sync[num_source], Source.ω_sync[num_source, end])

    Phase_Locked_Loop_3ph(Source, num_source)

    return nothing
end

function Phase_Locked_Loop_3ph(Source::Classical_Controls, num_source; ωn = 70, ξ = 0.35)

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

    Ki = ωn^2 # tuning
    Kp = ξ*2*sqrt(Ki) # tuning

    v_abc = Source.V_filt_cap[num_source, :, end]

    f = Source.fpll[num_source, 1, :]
    θ = Source.θpll[num_source, 1, end]

    err_t = Source.pll_err_t[num_source, 1]
    err = Source.pll_err[num_source, 1, :]

    v_αβγ = Clarke_Transform(v_abc)
    if norm(v_αβγ) != 0
        v_αβγ = sqrt(3)*v_αβγ./norm(v_αβγ)
    end
    err_new = v_αβγ[2]*cos(θ) - v_αβγ[1]*sin(θ)

    f_new, err_t_new, err_int =
    PI_Controller(err_new, err, err_t, Kp, Ki, Source.ts, bias = Source.fsys)

    θ = Third_Order_Integrator(θ, Source.ts, 2π*[f[2:end]; f_new])

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

function Phase_Locked_Loop_1ph(Source::Classical_Controls, num_source; Kp = 0.001, Ki = 1, ph = 1, k_sogi = 0.8)

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

function PQ_Control(Source::Classical_Controls, num_source, θ; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0])

    #-------------------------------------------------------------

    V_αβγ = Inv_Park_Transform(Source.V_dq0_inv[num_source, :], θ)
    I_αβγ_ref = Inv_p_q_v(V_αβγ, pq0_ref)
    Source.I_ref_dq0[num_source, :] = Park_Transform(I_αβγ_ref, θ)

    #-------------------------------------------------------------
    I_αβγ = Clarke_Transform(Source.I_filt_inv[num_source, :, end])

    V_αβγ_ref = Inv_p_q_i(I_αβγ, pq0_ref)

    V_dq0_ref = Park_Transform(V_αβγ_ref, θ)

    if sqrt(2/3)*norm(V_dq0_ref) > Source.v_max[num_source]
        V_dq0_ref = V_dq0_ref.*((Source.v_max[num_source])/(sqrt(2/3)*norm(V_dq0_ref) ))
    end

    Source.V_ref[num_source, :] = Inv_DQ0_transform(V_dq0_ref, θ)

    Source.V_ref_dq0[num_source, :] = V_dq0_ref

    #-------------------------------------------------------------

    return nothing
end

function PV_Control(Source::Classical_Controls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0])
    
    Vn = sqrt(2)*Source.V_pu_set[num_source, 1]*Source.Vrms[num_source] #peak
    Vg = sqrt(2/3)*norm(Source.V_dq0_inv[num_source, :]) #peak

    Kp = Source.V_kp[num_source]
    Ki = 250*Source.V_ki[num_source]

    V_err = Source.V_err[num_source, :, 1]
    V_err_t = Source.V_err_t[num_source, 1]

    V_err_new = Vn - Vg

    q_ref, V_err_t, Source.V_err[num_source, :, 1] =
    PI_Controller(V_err_new, V_err, V_err_t, Kp, Ki, Source.ts, bias = 0)

    pq0_ref[2] = q_ref[1]
    Source.V_err_t[num_source, 1] = V_err_t[1]

    return pq0_ref
end

function Droop_Control(Source::Classical_Controls, num_source; Vrms = Source.Vrms[num_source])

    #= Theory
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

    ω = Source.ω_droop[num_source, 1, :]
    θ = Source.θ_droop[num_source, 1]
    p_q = Source.p_q_filt[num_source, :]

    ω_new = Source.fsys*2*π - p_q[1]*Source.D[num_source, 1]
    Source.ω_droop[num_source, :, end] = [ω_new; ω_new; ω_new]

    Source.θ_droop[num_source, 1] = Third_Order_Integrator(θ, Source.ts, [ω[2:end]; ω_new])%(2*π)
    Source.θ_droop[num_source, 2] = (Source.θ_droop[num_source, 1] - 120*π/180)%(2*π)
    Source.θ_droop[num_source, 3] = (Source.θ_droop[num_source, 1] + 120*π/180)%(2*π)

    E_new = Vrms - p_q[2]*Source.D[num_source, 2]

    Source.V_ref[num_source, 1] = sqrt(2)*E_new*cos(Source.θ_droop[num_source, 1])
    Source.V_ref[num_source, 2] = sqrt(2)*E_new*cos(Source.θ_droop[num_source, 2])
    Source.V_ref[num_source, 3] = sqrt(2)*E_new*cos(Source.θ_droop[num_source, 3])

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

function Current_Controller(Source::Classical_Controls, num_source, θ, ω; Kb = 1)

    #= Theory:
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

    Kp = Source.I_kp[num_source]
    Ki = Source.I_ki[num_source]

    Source.I_ref[num_source, :] = Inv_DQ0_transform(Source.I_ref_dq0[num_source, :], θ)
    Source.I_dq0[num_source, :] = DQ0_transform(Source.I_filt_inv[num_source, :, end], θ)
    Source.V_dq0[num_source, :] = DQ0_transform(Source.V_filt_cap[num_source, :, end], θ)
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

    Source.s_lim[num_source, :], Source.I_err_t[num_source, :], Source.I_err[num_source, :, :] =
    PI_Controller(I_err_new, I_err, I_err_t, Kp, Ki, Source.ts)

    # cross-coupling / feedforward
    Source.s_lim[num_source, 1] = Source.s_lim[num_source, 1] - 
    (Source.Lf_1[num_source]*ω*I_dq0[2] - V_dq0[1])*2/Source.Vdc[num_source]
    Source.s_lim[num_source, 2] = Source.s_lim[num_source, 2] + 
    (Source.Lf_1[num_source]*ω*I_dq0[1] + V_dq0[2])*2/Source.Vdc[num_source]

    # ---- Limiting Output (Saturation)
    Vp_ref = (Source.Vdc[num_source]/2)*sqrt(2/3)*norm(Source.s_lim[num_source,:]) # peak set point

    if Vp_ref > Source.v_max[num_source]
        Source.s_dq0_avg[num_source, :]  = Source.s_lim[num_source, :]*Source.v_max[num_source]/Vp_ref
    else
        Source.s_dq0_avg[num_source, :]  = Source.s_lim[num_source, :]
    end

    Source.Vd_abc_new[num_source, :] = Inv_DQ0_transform(Source.s_dq0_avg[num_source, :] , θ)

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

function Voltage_Controller(Source::Classical_Controls, num_source, θ, ω; Kb = 1)

    Kp = Source.V_kp[num_source]
    Ki = Source.V_ki[num_source]

    Source.V_ref_dq0[num_source, :] = DQ0_transform(Source.V_ref[num_source, :], θ)
    Source.V_dq0[num_source, :] = DQ0_transform(Source.V_filt_cap[num_source, :, end], θ)

    I_dq0_poc = DQ0_transform(Source.I_filt_poc[num_source, :, end], θ)
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

    Source.I_lim[num_source, :], Source.V_err_t[num_source, :], Source.V_err[num_source, :, :] =
    PI_Controller(V_err_new, V_err, V_err_t, Kp, Ki, Source.ts)

    # cross-coupling / feedforward
    Source.I_lim[num_source, 1] = Source.I_lim[num_source, 1] + I_dq0_poc[1] 
    - Source.Cf[num_source]*ω*V_dq0[2] 
    Source.I_lim[num_source, 2] = Source.I_lim[num_source, 2] + I_dq0_poc[2] 
    + Source.Cf[num_source]*ω*V_dq0[1] 

    # ---- Limiting Output (Saturation)
    Ip_ref = sqrt(2/3)*norm(Source.I_lim[num_source,:]) # peak set point

    if Ip_ref > Source.i_max[num_source]
        Source.I_ref_dq0[num_source, :] = Source.I_lim[num_source, :]*Source.i_max[num_source]/Ip_ref
    else
        Source.I_ref_dq0[num_source, :] = Source.I_lim[num_source, :]
    end

    return nothing
end

function Synchronverter_Control(Source::Classical_Controls, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0], mode = 2, Vrms = Source.Vrms[num_source])

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
    Vg = sqrt(2/3)*norm(DQ0_transform(Source.V_filt_cap[num_source, :, end], 0)) # peak measured voltage

    #---- Integrate eq_new to find Mfif_new

    if mode == 2 || mode == 5
        eq_new = (1/Source.K_sync[num_source])*(Dq*(Vn - Vg))
    else
        eq_new = (1/Source.K_sync[num_source])*(pq0_ref[2] + Dq*(Vn - Vg) - Source.p_q_filt[num_source, 2])
    end

    Mfif_new = Third_Order_Integrator(Source.Mfif[num_source], Source.ts, [eq[2:end]; eq_new])

    Source.Mfif[num_source] = Mfif_new
    Source.eq[num_source, :] = [eq[2:end]; eq_new]
    #----

    #---- Integrate α_new to find ω_new
    Tm = pq0_ref[1]/ωsys # Virtual machine Torque

    ω_err_new = ω[end] - ωsys - Source.ω_set[num_source]
    ΔT = Source.D[num_source, 1]*ω_err_new

    #~~~ PI Controller

    if S1 == 1 && 1 == 2

        Kp = 0.0001
        Ki = 0.001
        ω_set, ΔT_err_t, Source.ΔT_err[num_source, :] =
        PI_Controller([-ΔT], Source.ΔT_err[num_source, :], Source.ΔT_err_t[num_source], 
        Kp, Ki, Source.ts) 

        Source.ω_set[num_source] = ω_set[1]
        Source.ΔT_err_t[num_source] = ΔT_err_t[1]
    else

        Source.ω_set[num_source] = 0
        Source.ΔT_err_t[num_source] = 0
        Source.ΔT_err[num_source, :] = fill!(Source.ΔT_err[num_source, :], 0.0)
    end
    #~~~

    Te_new = Source.p_q_filt[num_source, 1]/ω[end] # New Electrical Torque

    α_new = (1/Source.J_sync[num_source])*(Tm - Te_new - ΔT) # New Angular Acceleration

    ω_new = Third_Order_Integrator(ω[end], Source.ts, [α[2:end]; α_new])

    Source.ω_sync[num_source, :] = [ω[2:end]; ω_new]
    Source.α_sync[num_source, :] = [α[2:end]; α_new]

    #----

    #---- Integrate ω_new to find θ_new
    θ_new = Third_Order_Integrator(θ, Source.ts, [ω[2:end]; ω_new])%(2π)

    Source.θ_sync[num_source] = θ_new
    #----

    #----
    cos_θ_new = cos.([θ_new; θ_new - 120*π/180; θ_new + 120*π/180])
    Source.V_ref[num_source, :] = ω_new*Mfif_new*cos_θ_new # three phase generated voltage
    #----

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

function Filtering(Source::Classical_Controls, num_source, θ)

    V_inv = (Source.Vdc[num_source]/2)*Source.Vd_abc_new[num_source, :]

    Source.V_pre_dq0[num_source, :, 1:end-1] = Source.V_pre_dq0[num_source, :, 2:end]
    Source.V_pre_dq0[num_source, :, end] = DQ0_transform(V_inv, θ)

    Source.V_dq0_inv[num_source, :] = First_Order_LPF(500, Source.V_pre_dq0[num_source, :, :], 
    Source.V_dq0_inv[num_source, :], Source.ts)

    return nothing
end

function Luenberger_Observer(Source::Classical_Controls, num_source)

    ns = num_source

    if Source.filter_type[ns] == "LCL"

        A = Source.Ad[ns, 2:3, 2:3] #A22
        A11 = Source.Ad[ns, 1, 1]
        B = [Source.Ad[ns, 2:3, 1] Source.Bd[ns, 2:3] Source.Dd[ns, 2:3]]
        B1 = Source.Bd[ns, 1]
        C = transpose(Source.Ad[ns, 1, 2:3]) #A12
        D1 = Source.Dd[ns, 1]

        K = Source.Ko[ns, :]

        for ph in 1:Source.phases

            yₚ = Source.I_filt_inv[ns, ph, end - 1]
            vₚ = (Source.Vdc[ns]/2)*Source.Vd_abc_new[ns, ph]
            eₚ = Source.V_filt_poc[ns, ph, end - 1]

            u = [yₚ; vₚ; eₚ]

            xpₚ = [Source.V_filt_cap[ns, ph, end-1];
                  Source.I_filt_poc[ns, ph, end-1]]

            wpₚ = xpₚ - yₚ*K

            wp = (A - K*C*A)*wpₚ + (A*K - K*A11 - K*C*A*K)*yₚ - K*(B1*vₚ + D1*eₚ + C*B*u)

            xp = wp + K*Source.I_filt_inv[ns, ph, end]

            Source.V_filt_cap[ns, ph, end] = xp[1]
            Source.I_filt_poc[ns, ph, end] = xp[2]
        end

    end

    return nothing
end

#-------------------------------------------------------------------------------

function Measurements(Source::Classical_Controls)

    i = Source.steps

    Source.num_sources

    i_sp_rms = convert(Int64, round((1/(Source.ts*Source.T_sp_rms))))

    i_start = i - convert(Int64, round(Source.T_eval/(Source.fsys*Source.ts)))

    i_range = 1:i_sp_rms:Source.N

    for ns in 1:Source.num_sources

        V_poc = Source.V_filt_cap[ns, :, end]
        I_poc = Source.I_filt_poc[ns, :, end]

        Source.p_q_inst[ns, :] = p_q_theory(V_poc, I_poc) # real and imaginary powers

        Source.p_inst[ns, :] = V_poc.*I_poc

        # Phasors
        if i_start >= 1 #waiting for one evaluation cycle to pass

            θ = Source.θpll[ns, 1, i_range]

            # Voltages
            v_signals = transpose(Source.V_filt_cap[ns, :, i_range])
            Source.V_ph[ns,  :, :] = RMS(θ, v_signals)
            Source.V_ph[ns,  :, 3] = Source.V_ph[ns,  :, 3]# .+ π/2

            # Currents
            i_signals = transpose(Source.I_filt_poc[ns, :, i_range])
            Source.I_ph[ns, :, :] = RMS(θ, i_signals)
            Source.I_ph[ns, :, 3] = Source.I_ph[ns, :, 3]# .+ π/2

            # Per phase (and total) Active and Reactiv Powers
            Source.Pm[ns, 1:3] = (Source.V_ph[ns, :, 2].*Source.I_ph[ns, :, 2]).*cos.(Source.V_ph[ns, :, 3] .- Source.I_ph[ns, :, 3])
            Source.Pm[ns, 4] = sum(Source.Pm[ns, 1:3])
            Source.Qm[ns, 1:3] = (Source.V_ph[ns, :, 2].*Source.I_ph[ns, :, 2]).*sin.(Source.V_ph[ns, :, 3] .- Source.I_ph[ns, :, 3])
            Source.Qm[ns, 4] = sum(Source.Qm[ns, 1:3])
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
    A = [-Source.Rf_L1[num_source]/Source.Lf_1[num_source]]
    B = [1/(Source.Lf_1[num_source])]
    C = [1]
    D = [0]
    sys = ss(A, B, C, D) # continuous
    tf_sys = tf(sys)

    Ts = Source.ts
    dly = 1*Ts #To do: action_delay from env
    ZoH = tf([1], [Ts/2, 1]) # Transfer function of the sample and hold process
    Pade = tf([-dly/2, 1], [dly/2, 1]) # Pure first-order delay approximation
    PWM_gain = Source.Vdc[num_source]/2

    #SC = tf([1/(1*Lf_1)], [1, Rf_L1/Lf_1]) # = tf(sys_sc)
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
            
            println("\nWARNING: PI Current Controller with Positive Poles.")
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

            println("\nWARNING: PI Voltage Controller with Positive Poles.")
            println("Suggestion: Decrease Simulation Time Step")
            println("Source = ", num_source,"\n")
        end
    end

    return nothing
end

function Source_Initialiser(env, Source, modes, source_indices; pf = 0.8)

    Mode_Keys = [k[1] for k in sort(collect(Source.Modes), by = x -> x[2])]

    e = 1

    Source.fsys = env.nc.parameters["grid"]["f_grid"]
    Source.Δfmax = env.nc.parameters["grid"]["Δfmax"]
    Source.ΔEmax = env.nc.parameters["grid"]["ΔEmax"]

    for ns in source_indices

        Srated = env.nc.parameters["source"][ns]["pwr"]
        Source.filter_type[e] = env.nc.parameters["source"][ns]["fltr"]
        Source.Rf_L1[e] = env.nc.parameters["source"][ns]["R1"]
        Source.Lf_1[e] = env.nc.parameters["source"][ns]["L1"]
        Source.Vdc[e] = env.nc.parameters["source"][ns]["vdc"]
        Source.Vrms[e] = env.nc.parameters["grid"]["v_rms"]

        if haskey(env.nc.parameters["source"][1], "R_C")
            Source.Rf_C[e] = env.nc.parameters["source"][ns]["R_C"]
            Source.Cf[e] = env.nc.parameters["source"][ns]["C"]
        else
            Source.Cf[e] = 0.00001e-6
            Source.Rf_C[e] = 1*Source.Cf[e]
        end

        if Source.filter_type[e] == "LCL"
            Source.Rf_L2[e] = env.nc.parameters["source"][ns]["R2"]
            Source.Lf_2[e] = env.nc.parameters["source"][ns]["L2"]
        end

        Source.S[e] = Srated
        Source.P[e] = pf*Srated
        Source.Q[e] = sqrt(Srated^2 - Source.P[e]^2)

        if typeof(modes[e]) == Int64
            Source.Source_Modes[e] = Mode_Keys[modes[e]]
        else
            Source.Source_Modes[e] = modes[e]
        end

        #Droop (and synchronverter) parameters
        Source.τv[e] = env.nc.parameters["source"][ns]["τv"] # time constant of the voltage loop
        Source.τf[e] = env.nc.parameters["source"][ns]["τf"] # time constant of the frequency loop

        if Source.Source_Modes[e] == "Semi-Synchronverter" || Source.Source_Modes[e] == "Full-Synchronverter"

            Source.D[e, 1] = Source.S[e]/((2π)*(Source.fsys)*Source.fsys*2π*Source.Δfmax)
            Source.D[e, 2] = Source.S[e]/(Source.Vrms[e]*sqrt(2)*Source.ΔEmax)
        else

            Source.D[e, 1] = 2π*Source.fsys*2π*Source.Δfmax/Source.P[e]
            Source.D[e, 2] = Source.Vrms[e]*Source.ΔEmax/Source.Q[e]
        end

        Source.J_sync[e] = Source.τf[e]*Source.D[e, 1] #total mass moment of inertia of the rotating masses, kg*m^2
        Source.K_sync[e] = Source.τv[e]*Source.fsys*2π*Source.D[e, 2]

        Source.pq0_set[e, 1] = env.nc.parameters["source"][ns]["p_set"] # W, Real Power
        Source.pq0_set[e, 2] = env.nc.parameters["source"][ns]["q_set"] # VAi, Imaginary Power

        Source.V_pu_set[e, 1] = env.nc.parameters["source"][ns]["v_pu_set"]
        Source.V_δ_set[e, 1] = env.nc.parameters["source"][ns]["v_δ_set"]*π/180

        Source.i_max[e] = env.nc.parameters["source"][ns]["i_limit"] # should there be a 3? TODO
        Source.v_max[e] = env.nc.parameters["source"][ns]["v_limit"]/2

        Current_PI_LoopShaping(Source, e)
        Voltage_PI_LoopShaping(Source, e)

        Observer_Initialiser(Source, e)

        e += 1
    end

    return nothing
end

function Observer_Initialiser(Source::Classical_Controls, num_source)

    # Current Deadbeat Reduced-Order Observer initialisation

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

        Source.Bd[ns, :] = [1/Source.Lf_1[ns]; 0; 0]
        Source.Dd[ns, :] = [0; -1/Source.Lf_1[ns]; 0]

        Source.Ad[ns, :, :] = exp(A*Source.ts)
        
        Source.Bd[ns, :] = inv(A)*(Source.Ad[ns, :, :] - I)*Source.Bd[ns, :]
        Source.Dd[ns, :] = inv(A)*(Source.Ad[ns, :, :] - I)*Source.Dd[ns, :]

        Source.Cd[ns, :] = [1. 0. 0.] # measuring the inverter current

        #A11 = Source.Ad[ns, 1, 1]
        #A12 = Source.Ad[ns, 1, 2:3]
        #A21 = Source.Ad[ns, 2:3, 1]
        A22 = Source.Ad[ns, 2:3, 2:3]

        C = transpose(Source.Ad[ns, 1, 2:3])*A22
        _, r = Observability(C, A22)

        if r != size(A22, 1)
            println("\nERROR: The system is not observable. The rank of 'O' is not equal to $(size(A22, 1)).")
        end

        Source.Ko[ns, 1] = -(C[2]*A22[1, 1]^2 - A22[1, 2]*C[1]*A22[1, 1] + A22[1, 2]*A22[2, 1]*C[2] - A22[1, 2]*A22[2, 2]*C[1])/(A22[1, 2]*C[1]^2 - A22[2, 1]*C[2]^2 - A22[1, 1]*C[1]*C[2] + A22[2, 2]*C[1]*C[2])
        Source.Ko[ns, 2] = (C[1]*A22[2, 2]^2 - A22[2, 1]*C[2]*A22[2, 2] - A22[1, 1]*A22[2, 1]*C[2] + A22[1, 2]*A22[2, 1]*C[1])/(A22[1, 2]*C[1]^2 - A22[2, 1]*C[2]^2 - A22[1, 1]*C[1]*C[2] + A22[2, 2]*C[1]*C[2])
        
        #= println("Ad = ", Source.Ad[ns, :, :])
        println("Bd = ", Source.Bd[ns, :])
        println("Dd = ", Source.Dd[ns, :])
        println("Ko = ", Source.Ko[ns, :])
        println("O, r = ", O, " ", r) =#
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

