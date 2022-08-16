mutable struct Source_Controller <: AbstractPolicy

    #=
        This object contains the functions and parameters that are relevant for the
        control of a three-phase half bridge DC to AC Converter.
    =#

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
    # General System & Control

    Modes::Dict{String, Int64}
    Source_Modes::Vector{String}

    num_sources::Int64

    f_cntr::Float64 # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz
    fsys::Float64
    θsys::Float64
    μ_cntr::Float64 # s, sampling timestep
    N_cntr::Int64
    delay::Int64
    steps::Int64

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

    Vd_abc_new::Array{Float64}

    #---------------------------------------------------------------------------
    # Voltage Controller

    Gv_cl::Array{TransferFunction} # Closed Loop transfer function

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
    PQ_err::Array{Float64}
    PQ_err_t::Matrix{Float64}

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

    function Source_Controller(Vdc::Vector{Float64}, Vrms::Vector{Float64},
        S::Vector{Float64}, P::Vector{Float64}, Q::Vector{Float64},
        i_max::Vector{Float64}, v_max::Vector{Float64},
        Lf::Vector{Float64}, Cf::Vector{Float64}, Rf::Vector{Float64},
        Modes::Dict{String, Int64}, Source_Modes::Vector{String},
        num_sources::Int64,
        f_cntr::Float64, fsys::Float64, θsys::Float64,
        μ_cntr::Float64, N_cntr::Int64, delay::Int64, steps::Int64,
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
        Gv_cl::Array{TransferFunction},
        V_dq0::Array{Float64}, V_ref_dq0::Array{Float64}, V_ref::Array{Float64},
        V_err::Array{Float64}, V_err_t::Matrix{Float64},
        V_kp::Vector{Float64}, V_ki::Vector{Float64},
        I_lim::Matrix{Float64},
        Δfmax::Float64, ΔEmax::Float64, τv::Float64, τf::Float64,
        D::Matrix{Float64}, ω_droop::Array{Float64}, θ_droop::Array{Float64},
        pq0_set::Matrix{Float64}, PQ_err::Array{Float64}, PQ_err_t::Matrix{Float64},
        J_sync::Vector{Float64}, K_sync::Vector{Float64}, ΔT_t::Vector{Float64},
        α_sync::Matrix{Float64}, ω_sync::Matrix{Float64}, θ_sync::Matrix{Float64},
        Δω_sync::Matrix{Float64}, eq::Matrix{Float64}, Mfif::Matrix{Float64})

        new(Vdc, Vrms,
        S, P, Q,
        i_max, v_max,
        Lf, Cf, Rf,
        Modes, Source_Modes,
        num_sources,
        f_cntr, fsys, θsys,
        μ_cntr, N_cntr, delay, steps,
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
        Gv_cl,
        V_dq0, V_ref_dq0, V_ref,
        V_err, V_err_t,
        V_kp, V_ki,
        I_lim,
        Δfmax, ΔEmax, τv, τf,
        D, ω_droop, θ_droop,
        pq0_set, PQ_err, PQ_err_t,
        J_sync, K_sync, ΔT_t,
        α_sync, ω_sync, θ_sync,
        Δω_sync, eq, Mfif)
    end

    function Source_Controller(t_final, f_cntr, num_sources; delay = 1)

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
        Rf = fill!(Rf, 0.001)

        #---------------------------------------------------------------------------
        # General System & Control

        Modes = Dict("Swing Mode" => 1, "Voltage Control Mode" => 2, "PQ Control Mode" => 3,
        "Droop Control Mode" => 4, "Synchronverter Mode" => 5)

        Source_Modes = Array{String, 1}(undef, num_sources)
        Source_Modes = fill!(Source_Modes, "Voltage Control Mode")

        steps = 0
        fsys = 50.0
        θsys = 0.0

        μ_cntr = 1/f_cntr
        N_cntr = convert(Int64, floor(t_final/μ_cntr)) + 1

        t_final = convert(Float64, t_final)
        t = 0:μ_cntr:t_final
        θt = (2*π*fsys*t).%(2π)

        #---------------------------------------------------------------------------
        # Phase Locked Loops

        vd = Array{Float64, 3}(undef, num_sources, 3, N_cntr) #3 for 3 phases
        vd = fill!(vd, 00.0)
        qvd = Array{Float64, 3}(undef, num_sources, 3, N_cntr) #3 for 3 phases
        qvd = fill!(qvd, 00.0)

        fpll = Array{Float64, 3}(undef, num_sources, 3, N_cntr) #3 for 3 phases
        fpll = fill!(fpll, fsys)

        θpll = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        θpll = fill!(θpll, 0)

        pll_err = Array{Float64, 3}(undef, num_sources, 3, 3) # 3 phases and 3rd order integration
        pll_err = fill!(pll_err, 0)
        pll_err_t = Array{Float64, 2}(undef, num_sources, 3) # PLL total integrated error
        pll_err_t = fill!(pll_err_t, 0)

        #---------------------------------------------------------------------------
        # Interface (e.g. filters)

        V_filt_poc = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        V_filt_poc = fill!(V_filt_poc, 0)
        V_filt_inv = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        V_filt_inv = fill!(V_filt_inv, 0)

        I_filt_poc = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        I_filt_poc = fill!(I_filt_poc, 0)
        I_filt_inv = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        I_filt_inv = fill!(I_filt_inv, 0)

        p_q_filt = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        p_q_filt  = fill!(p_q_filt, 0)

        #---------------------------------------------------------------------------
        # Current Controller

        Gi_cl = Array{TransferFunction, 1}(undef, num_sources)

        I_dq0 = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        I_dq0 = fill!(I_dq0, 0)
        I_ref_dq0 = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        I_ref_dq0 = fill!(I_ref_dq0, 0)
        I_ref = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        I_ref = fill!(I_ref, 0)

        # Current Integrations
        I_err = Array{Float64, 3}(undef, num_sources, 3, 3) # 3 phases and 3rd order integration
        I_err = fill!(I_err, 0)
        I_err_t = Array{Float64, 2}(undef, num_sources, 3)
        I_err_t = fill!(I_err_t, 0)
        I_kp = Array{Float64, 1}(undef, num_sources)
        I_kp = fill!(I_kp, 0.5)
        I_ki = Array{Float64, 1}(undef, num_sources)
        I_ki = fill!(I_ki, 25)

        Vd_abc_new = Array{Float64, 3}(undef, num_sources, 3, N_cntr + delay)
        Vd_abc_new = fill!(Vd_abc_new, 0)

        #---------------------------------------------------------------------------
        # Voltage Controller

        Gv_cl = Array{TransferFunction, 1}(undef, num_sources)

        V_dq0 = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        V_dq0 = fill!(V_dq0, 0)
        V_ref_dq0 = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        V_ref_dq0 = fill!(V_ref_dq0, 0)

        V_ref = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        V_ref_p = sqrt(2)*Vrms[1]
        for j in 1:num_sources
            V_ref[j, 1, :] = V_ref_p*cos.(θt)
            V_ref[j, 2, :] = V_ref_p*cos.(θt .- 120*π/180)
            V_ref[j, 3, :] = V_ref_p*cos.(θt .+ 120*π/180)
        end

        # Voltage Integrations
        V_err = Array{Float64, 3}(undef, num_sources, 3, 3) # 3 phases and 3rd order integration
        V_err = fill!(V_err, 0)
        V_err_t = Array{Float64, 2}(undef, num_sources, 3)
        V_err_t = fill!(V_err_t, 0)
        V_kp = Array{Float64, 1}(undef, num_sources)
        V_kp = fill!(V_kp, 0.01)
        V_ki = Array{Float64, 1}(undef, num_sources)
        V_ki = fill!(V_ki, 20)

        I_lim = Array{Float64, 2}(undef, num_sources, 3)
        I_lim = fill!(I_lim, 0)

        #---------------------------------------------------------------------------
        # Droop (Classical, Synchronverter, and VSG) Mode

        Δfmax = 0.25 # Hz
        ΔEmax = 10.0 # V
        τv = 0.002 # time constant of the voltage loop
        τf = 0.002

        D = Array{Float64, 2}(undef, num_sources, 2)
        D[:,1] = fill!(D[:,1], 2π*Δfmax/P[1])
        D[:,2] = fill!(D[:,2], ΔEmax/Q[1])

        ω_droop = Array{Float64, 3}(undef, num_sources, 3, N_cntr) #3 for 3 phases
        ω_droop = fill!(ω_droop, fsys*2π)

        θ_droop = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        θ_droop = fill!(θ_droop, 0)

        #---------------------------------------------------------------------------
        # PQ Mode

        pq0_set = Array{Float64, 2}(undef, num_sources, 3) # real, imaginary, and zero power set points
        pq0_set = fill!(pq0_set, 0)
        PQ_err = Array{Float64, 3}(undef, num_sources, 3, 3) # 3 phases and 3rd order integration
        PQ_err = fill!(PQ_err, 0)
        PQ_err_t = Array{Float64, 2}(undef, num_sources, 3) # PLL total integrated error
        PQ_err_t = fill!(PQ_err_t, 0)

        #---------------------------------------------------------------------------
        # Synchronverter Mode

        J_sync = Array{Float64, 1}(undef, num_sources)
        J_sync = fill!(J_sync, 0)
        K_sync = Array{Float64, 1}(undef, num_sources)
        K_sync = fill!(K_sync, 0)
        ΔT_t = Array{Float64, 1}(undef, num_sources)
        ΔT_t = fill!(ΔT_t, 0)

        α_sync = Array{Float64, 2}(undef, num_sources, N_cntr)
        α_sync = fill!(α_sync, 0)
        ω_sync = Array{Float64, 2}(undef, num_sources, N_cntr)
        ω_sync = fill!(ω_sync, fsys*2π)
        θ_sync = Array{Float64, 2}(undef, num_sources, N_cntr)
        θ_sync = fill!(θ_sync, 0)
        Δω_sync = Array{Float64, 2}(undef, num_sources, N_cntr)
        Δω_sync = fill!(Δω_sync, 0)
        eq = Array{Float64, 2}(undef, num_sources, N_cntr)
        eq = fill!(eq, 0)
        Mfif = Array{Float64, 2}(undef, num_sources, N_cntr)
        Mfif = fill!(Mfif, 0)

        Source_Controller(Vdc, Vrms,
        S, P, Q,
        i_max, v_max,
        Lf, Cf, Rf,
        Modes, Source_Modes,
        num_sources,
        f_cntr, fsys, θsys,
        μ_cntr, N_cntr, delay, steps,
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
        Gv_cl,
        V_dq0, V_ref_dq0, V_ref,
        V_err, V_err_t,
        V_kp, V_ki,
        I_lim,
        Δfmax, ΔEmax, τv, τf,
        D, ω_droop, θ_droop,
        pq0_set, PQ_err, PQ_err_t,
        J_sync, K_sync, ΔT_t,
        α_sync, ω_sync, θ_sync,
        Δω_sync, eq, Mfif)
    end
end

mutable struct Environment

    steps::Int64

    μ::Float64 # s, sampling timestep
    N::Int64 # number of samples
    t_final::Float64 # total simulation time

    fsys::Float64 # Mostly for plotting
    fs::Vector{Float64}
    θs::Vector{Float64}
    T_sp_rms::Float64

    V_ph::Array{Float64}
    I_ph::Array{Float64}

    p_q_inst::Array{Float64}
    p_inst::Array{Float64}
    P::Array{Float64}
    Q::Array{Float64}

    A::Matrix{Float64}
    B::Vector{Float64}
    C::Matrix{Float64}
    D::Vector{Float64}

    Ad::Matrix{Float64}
    Bd::Vector{Float64}

    x::Matrix{Float64}
    y::Matrix{Float64}
    u::Matrix{Float64}

    V_poc_loc::Matrix{Int64}
    I_poc_loc::Matrix{Int64}
    I_inv_loc::Matrix{Int64}

    V_load_loc::Matrix{Int64}
    I_load_loc::Matrix{Int64}

    num_sources::Int64
    num_loads::Int64

    function Environment(steps::Int64, μ::Float64, N::Int64, t_final::Float64,
        fsys::Float64, fs::Vector{Float64}, θs::Vector{Float64},
        T_sp_rms::Float64,
        V_ph::Array{Float64}, I_ph::Array{Float64},
        p_q_inst::Array{Float64}, p_inst::Array{Float64}, P::Array{Float64}, Q::Array{Float64},
        A::Matrix{Float64}, B::Vector{Float64}, C::Matrix{Float64}, D::Vector{Float64},
        Ad::Matrix{Float64}, Bd::Vector{Float64},
        x::Matrix{Float64}, y::Matrix{Float64}, u::Matrix{Float64},
        V_poc_loc::Matrix{Int64}, I_poc_loc::Matrix{Int64}, I_inv_loc::Matrix{Int64},
        V_load_loc::Matrix{Int64}, I_load_loc::Matrix{Int64},
        num_sources::Int64, num_loads::Int64)

        new(steps, μ, N, t_final,
        fsys, fs, θs,
        T_sp_rms,
        V_ph, I_ph,
        p_q_inst, p_inst, P, Q,
        A, B, C, D,
        Ad, Bd,
        x, y, u,
        V_poc_loc, I_poc_loc, I_inv_loc,
        V_load_loc, I_load_loc,
        num_sources, num_loads)
    end

    function Environment(t_final, μ, A, B, C, D, num_sources, num_loads)

        steps = 1
        N = convert(Int64, floor(t_final/μ)) + 1
        num_nodes = num_sources + num_loads

        t_final = convert(Float64, t_final)
        fsys = 50.0
        T_sp_rms = 5*fsys #samples in a second for rms calcs, x*fsys = x samples in a cycle

        fs = Array{Float64, 1}(undef, N)
        fs = fill!(fs, fsys)

        θs = Array{Float64, 1}(undef, N)
        θs = fill!(θs, 0)

        # RMS Phasors
        V_ph = Array{Float64, 4}(undef, num_nodes, 3, 3, N)
        V_ph = fill!(V_ph, 0)
        I_ph = Array{Float64, 4}(undef, num_nodes, 3, 3, N)
        I_ph = fill!(I_ph, 0)

        # Instantaneous Real, Imaginary, and Zero powers
        p_q_inst = Array{Float64, 3}(undef, num_nodes, 3, N)
        p_q_inst = fill!(p_q_inst, 0)

        p_inst = Array{Float64, 3}(undef, num_nodes, 3, N) # instantaneous power at PCC
        p_inst = fill!(p_inst, 0)

        P = Array{Float64, 3}(undef, num_nodes, 4, N) # 4th column is total
        P = fill!(P, 0)
        Q = Array{Float64, 3}(undef, num_nodes, 4, N) # 4th column is total
        Q = fill!(Q, 0)

        Ad = exp(A*μ)
        Bd = inv(A)*(Ad - Matrix(I, size(A, 1), size(A, 1)))*B

        x = Array{Float64, 2}(undef, size(A, 1), N)
        x = fill!(x, 0)

        y = Array{Float64, 2}(undef, size(A, 1), N)
        y = fill!(y, 0)

        u = Array{Float64, 2}(undef, size(A, 1), N+10)
        u = fill!(u, 0)

        V_poc_loc = Array{Int64, 2}(undef, 3, num_sources)
        V_poc_loc = fill!(V_poc_loc, 1)
        I_poc_loc = Array{Int64, 2}(undef, 3, num_sources)
        I_poc_loc = fill!(I_poc_loc, 1)
        I_inv_loc = Array{Int64, 2}(undef, 3, num_sources)
        I_inv_loc = fill!(I_inv_loc, 1)

        V_load_loc = Array{Int64, 2}(undef, 3, num_loads)
        V_load_loc = fill!(V_load_loc, 1)
        I_load_loc = Array{Int64, 2}(undef, 3, num_loads)
        I_load_loc = fill!(I_load_loc, 1)

        Environment(steps, μ, N, t_final,
        fsys, fs, θs,
        T_sp_rms,
        V_ph, I_ph,
        p_q_inst, p_inst, P, Q,
        A, B, C, D,
        Ad, Bd,
        x, y, u,
        V_poc_loc, I_poc_loc, I_inv_loc,
        V_load_loc, I_load_loc,
        num_sources, num_loads)
    end
end

function Classical_Policy(Source::Source_Controller, Env::Environment)

    Source_Interface(Env, Source, active = 0, fc = 2000)

    for s in 1:Source.num_sources

        if Source.Modes[Source.Source_Modes[s]] == 1

            Swing_Mode(Source, s)
        elseif Source.Modes[Source.Source_Modes[s]] == 2

            Voltage_Control_Mode(Source, s)
        elseif Source.Modes[Source.Source_Modes[s]] == 3

            PQ_Control_Mode(Source, s, Source.pq0_set[s,:])
        elseif Source.Modes[Source.Source_Modes[s]] == 4

            Droop_Control_Mode(Source, s)
        elseif Source.Modes[Source.Source_Modes[s]] == 5

            Synchronverter_Mode(Source, s)
        end
    end

    Action = Env_Interface(Env, Source)

    return Action
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

function Swing_Mode(Source::Source_Controller, num_source; δ = 0, pu = 1, ramp = 0, t_end = 0.04)

    i = Source.steps

    θt = Source.θsys
    θph = [θt + δ; θt + δ - 120π/180; θt + δ + 120π/180]
    Vrms = V_Ramp(pu*Source.Vrms[num_source], Source.μ_cntr, i; t_end = t_end, ramp = ramp)
    Source.V_ref[num_source, :, i] = sqrt(2)*Vrms*sin.(θph)

    Source.Vd_abc_new[num_source, :, i] = Source.V_ref[num_source, :, i]
    Phase_Locked_Loop_3ph(Source, num_source)

    return nothing
end

function Voltage_Control_Mode(Source::Source_Controller, num_source; δ = 0, pu = 1, ramp = 0, t_end = 0.04)

    i = Source.steps

    θt = Source.θsys# + π/2
    θph = [θt + δ; θt + δ - 120π/180; θt + δ + 120π/180]
    Vrms = Ramp(pu*Source.Vrms[num_source], Source.μ_cntr, i; t_end = t_end, ramp = ramp)
    Source.V_ref[num_source, :, i] = sqrt(2)*Vrms*cos.(θph)
    
    Phase_Locked_Loop_3ph(Source, num_source)
    Voltage_Controller(Source, num_source, θt)
    Current_Controller(Source, num_source, θt)

    return nothing
end

function Droop_Control_Mode(Source::Source_Controller, num_source; ramp = 0, t_end = 0.04)

    i = Source.steps

    Vrms = Ramp(Source.Vrms[num_source], Source.μ_cntr, i; t_end = t_end, ramp = ramp)
    Dout = D_Ramp([2π*Source.Δfmax; Source.ΔEmax], Source.μ_cntr, i, t_end = t_end, ramp = ramp)

    Source.D[num_source, 1] = Dout[1]/Source.P[num_source]
    Source.D[num_source, 2] = Dout[2]/Source.Q[num_source]

    Droop_Control(Source, num_source, Vrms = Vrms)
    θt = Source.θ_droop[num_source, 1, i]

    Phase_Locked_Loop_3ph(Source, num_source)
    Voltage_Controller(Source, num_source, θt)
    Current_Controller(Source, num_source, θt)

    return nothing
end

function PQ_Control_Mode(Source::Source_Controller, num_source, pq0)

    i = Source.steps

    Phase_Locked_Loop_3ph(Source, num_source)
    θt = Source.θpll[num_source, 1 , i]

    if i*Source.μ_cntr > 2/Source.fsys
        PQ_Control(pq0_ref = pq0, Source, num_source, θt)
    else
        Voltage_Controller(Source, num_source, θt)
    end

    Phase_Locked_Loop_3ph(Source, num_source)
    Current_Controller(Source, num_source, θt)

    return nothing
end

function Synchronverter_Mode(Source::Source_Controller, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]], ramp = 0, t_end = 0.04)

    i = Source.steps

    Vrms = Ramp(Source.Vrms[num_source], Source.μ_cntr, i; t_end = t_end, ramp = ramp)
    Dout = D_Ramp([2π*Source.Δfmax; Source.ΔEmax], Source.μ_cntr, i, t_end = t_end, ramp = ramp)

    Source.D[num_source, 1] = Dout[1]/Source.P[num_source]
    Source.D[num_source, 2] = Dout[2]/Source.Q[num_source]

    # Synchronverter parameters
    Dp_sync = 1/(Source.D[num_source, 1]*(2*π)*Source.fsys)
    Source.J_sync[num_source] = Source.τf*Dp_sync

    Dq_sync = sqrt(2)/(Source.D[num_source, 2])
    Source.K_sync[num_source] = Source.τv*Source.fsys*2π*Dq_sync

    if i*Source.μ_cntr > 3/Source.fsys
        Synchronverter_Control(Source, num_source; pq0_ref = pq0_ref)
        θ_S = Source.θ_sync[num_source, i]
        Voltage_Controller(Source, num_source, θ_S)
        Current_Controller(Source, num_source, θ_S)
    end

    return nothing
end

function Self_Synchronverter_Mode(Source::Source_Controller, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]], SQ = 1, SP = 1)

    i = Source.steps

    if i*Source.μ_cntr > 3/Source.fsys
        Self_Synchronverter_Control(Source, num_source, pq0_ref = pq0_ref, SQ = SQ, SP = SP)
        θ_S = Source.θ_sync[num_source, i]
        Voltage_Controller(Kp = 0.01, Ki = 20, Source, num_source, θ_S)
        Current_Controller(Kp = 0.5, Ki = 25, Source, num_source, θ_S)
    end

    return nothing
end

function Phase_Locked_Loop_3ph(Source::Source_Controller, num_source; Kp = 0.02, Ki = 0.01)

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

        v_αβγ = Clarke_Transform(v_abc)
        err_new = v_αβγ[2]*cos(θ[n]) - v_αβγ[1]*sin(θ[n]) # this also works
    =#

    i = Source.steps

    v_abc = Source.V_filt_poc[num_source, :, i]

    range, cnt_end = Integrator_Prep(i)
    f = Source.fpll[num_source, 1, range]
    θ = Source.θpll[num_source, 1, cnt_end]

    err_t = Source.pll_err_t[num_source, 1]
    err = Source.pll_err[num_source, 1, :]

    err_new = 1*((v_abc[1] - v_abc[2])*cos(-θ)
    + (v_abc[3] - v_abc[2])*cos(-θ - 2π/3)) # this is magic

    f_new, err_t_new, err_int =
    PI_Controller(err_new, err, err_t, Kp, Ki, Source.μ_cntr, bias = Source.fsys)

    θ = Third_Order_Integrator(θ, Source.μ_cntr, 2π*[f; f_new])

    Source.fpll[num_source, :, i] = [f_new; f_new; f_new]
    Source.θpll[num_source, :, i] = [θ; θ - 120π/180; θ + 120π/180].%(2π)
    Source.pll_err_t[num_source, :] = [err_t_new; 0; 0]
    Source.pll_err[num_source, :, :] = [transpose(err_int); transpose(err_int); transpose(err_int)]

    return nothing
end

function Phase_Locked_Loop_1ph(Source::Source_Controller, num_source; Kp = 0.001, Ki = 1, ph = 1, k_sogi = 0.8)

    i = Source.steps

    range_1, cnt_end_1 = Integrator_Prep(i)
    range_2, cnt_end_2 = Integrator_Prep(i-1)

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
    vd_new = Third_Order_Integrator(vd_1[end], Source.μ_cntr, vd_sogi)

    qvd_int_old = vd_2.*ω_2
    qvd_int_new = vd_1[end]*ω_1[end]
    qvd_int = [qvd_int_old; qvd_int_new]
    qvd_new = Third_Order_Integrator(qvd_1[end], Source.μ_cntr, qvd_int)
    #----

    α_β_0 = [vd_new; qvd_new; 0]
    d_q_0 = Park_Transform(α_β_0, θ - π)

    vq = d_q_0[1]

    err_new = 0 - vq

    ω_new, err_t_new, err_int =
    PI_Controller(err_new, err, err_t, Kp, Ki, Source.μ_cntr, bias = 2*π*Source.fsys)

    Source.θpll[num_source, ph, i] =
    (Third_Order_Integrator(θ, Source.μ_cntr, [ω_1; ω_new]))%(2*π)

    Source.pll_err[num_source, ph, :] = err_int
    Source.vd[num_source, ph, i] = vd_new
    Source.qvd[num_source, ph, i] = qvd_new

    Source.pll_err_t[num_source, ph] = err_t_new[end]
    Source.fpll[num_source, ph, i] = ω_new[end]/(2*π)

    return nothing
end

function PQ_Control(Source::Source_Controller, num_source, θ; Kp = 0.1, Ki = 10, pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0])

    i = Source.steps

    i_abc = Source.I_filt_poc[num_source, :, i]
    v_abc = Source.V_filt_poc[num_source, :, i]
    PQ_err = Source.PQ_err[num_source, :, :]
    PQ_err_t = Source.PQ_err_t[num_source, :]

    V_αβγ = Clarke_Transform(v_abc)
    I_αβγ = Clarke_Transform(i_abc)

    I_αβγ_ref = Inv_p_q_v(V_αβγ, pq0_ref)

    I_ref_dq0 = Park_Transform(I_αβγ_ref, θ)

    I_dq0 = Park_Transform(I_αβγ, θ)

    PQ_err_new = (I_ref_dq0 - I_dq0)

    Source.I_ref_dq0[num_source, :, i], Source.PQ_err_t[num_source, :],
    Source.PQ_err[num_source, :, :] =
    PI_Controller(PQ_err_new, PQ_err, PQ_err_t, Kp, Ki, Source.μ_cntr)

    return nothing
end

function Droop_Control(Source::Source_Controller, num_source; Vrms = Source.Vrms[num_source])

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

    μ = Source.μ_cntr
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

function Current_Controller(Source::Source_Controller, num_source, θ)

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

    # --- to be removed
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
    PI_Controller(I_err_new, I_err, I_err_t, Kp, Ki, Source.μ_cntr)
    # Kp/1000, Ki/100

    Source.Vd_abc_new[num_source, :, i + Source.delay] = 0.5*Source.Vdc[num_source].*Inv_DQ0_transform(s_dq0_avg, θ)
    #=
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

function Voltage_Controller(Source::Source_Controller, num_source, θ; Kb = 1)

    i = Source.steps

    Kp = Source.V_kp[num_source]
    Ki = Source.V_ki[num_source]

    Source.V_ref_dq0[num_source, :,i] = DQ0_transform(Source.V_ref[num_source, :, i], θ)
    Source.V_dq0[num_source, :, i] = DQ0_transform(Source.V_filt_poc[num_source, :, i], θ)

    V_dq0 = Source.V_dq0[num_source, :, i]
    V_ref_dq0 = Source.V_ref_dq0[num_source, :, i]
    V_err = Source.V_err[num_source, :, :]
    V_err_t = Source.V_err_t[num_source, :]

    if i > 1
        # Including Anti-windup - Back-calculation
        V_err_new = V_ref_dq0 .- V_dq0 .+ Kb*(Source.I_ref_dq0[num_source, :, i - 1] .- Source.I_lim[num_source,:])
    else
        V_err_new = V_ref_dq0 .- V_dq0
    end

    Source.I_lim[num_source,:], Source.V_err_t[num_source, :], Source.V_err[num_source, :, :] =
    PI_Controller(V_err_new, V_err, V_err_t, Kp, Ki, Source.μ_cntr)

    # ---- Limiting Output (Saturation)
    Ip_ref = sqrt(2/3)*norm(Source.I_lim[num_source,:]) # peak set point

    if Ip_ref > Source.i_max[num_source]
        Source.I_ref_dq0[num_source, :, i] = Source.I_lim[num_source,:]*Source.i_max[num_source]/Ip_ref
    else
        Source.I_ref_dq0[num_source, :, i] = Source.I_lim[num_source,:]
    end

    return nothing
end

function Synchronverter_Control(Source::Source_Controller, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0], Vrms = Source.Vrms[num_source])

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
    Vn = sqrt(2)*Vrms # nominal peak POC voltage
    Vg = sqrt(2/3)*norm(DQ0_transform(Source.V_filt_poc[num_source, :, i], 0)) # peak measured voltage

    μ = Source.μ_cntr

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

function Self_Synchronverter_Control(Source::Source_Controller, num_source; pq0_ref = [Source.P[num_source]; Source.Q[num_source]; 0], SQ = 1, SP = 1, SC = 0, Kp = 0.0005, Ki = 0.01)


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

    μ = Source.μ_cntr
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

function Source_Interface(Env::Environment, Source::Source_Controller; active = 1, fc = 1000)

    i = Env.steps
    Source.steps = i
    ω = 2*π*Source.fsys
    Source.θsys = (Source.θsys + Source.μ_cntr*ω)%(2*π)

    for num_source in 1:Env.num_sources

        if active == 1

            x_start = i - 2
            if x_start < 1
                x_start = 1
            end
            x_range = x_start:i

            y_start = i - 2
            y_end = i - 1
            if y_start < 1
                y_start = 1
                y_end = 1
            end
            y_range = y_start:y_end

            V_poc = Env.x[Env.V_poc_loc[: , num_source], x_range]
            V_inv = Env.y[Env.I_inv_loc[: , num_source], x_range] .+ V_poc
            I_poc = Env.x[Env.I_poc_loc[: , num_source], x_range]
            I_inv = Env.x[Env.I_inv_loc[: , num_source], x_range]
            P_Q = Env.p_q_inst[num_source, :, x_range]

            Source.V_filt_poc[num_source, :, i] = Butterworth_LPF(fc, V_poc, Source.V_filt_poc[num_source, :, y_range], Source.μ_cntr)
            Source.V_filt_inv[num_source, :, i] = Butterworth_LPF(fc, V_inv, Source.V_filt_inv[num_source, :, y_range], Source.μ_cntr)
            Source.I_filt_poc[num_source, :, i] = Butterworth_LPF(fc, I_poc, Source.I_filt_poc[num_source, :, y_range], Source.μ_cntr)
            Source.I_filt_inv[num_source, :, i] = Butterworth_LPF(fc, I_inv, Source.I_filt_inv[num_source, :, y_range], Source.μ_cntr)
            Source.p_q_filt[num_source, :, i] = p_q_theory(Source.V_filt_poc[num_source, :, i], Source.I_filt_poc[num_source, :, i])
        else
            Source.V_filt_poc[num_source, :, i] = Env.x[Env.V_poc_loc[: , num_source], i]
            Source.V_filt_inv[num_source, :, i] = Env.y[Env.I_inv_loc[: , num_source], i] .+ Source.V_filt_poc[num_source, :, i]
            Source.I_filt_poc[num_source, :, i] = Env.x[Env.I_poc_loc[: , num_source], i]
            Source.I_filt_inv[num_source, :, i] = Env.x[Env.I_inv_loc[: , num_source], i]
            Source.p_q_filt[num_source, :, i] =  p_q_theory(Source.V_filt_poc[num_source, :, i], Source.I_filt_poc[num_source, :, i])
        end

    end

    return nothing
end

function Env_Interface(Env::Environment, Source::Source_Controller)

    i = Source.steps

    Action = zeros(length(Env.B), 1)

    for s in 1:Source.num_sources

        # Inverter Voltages - Control Actions
        #_______________________________________________________
        Action[Env.I_inv_loc[1, s]] = Source.Vd_abc_new[s, 1, i]
        Action[Env.I_inv_loc[2, s]] = Source.Vd_abc_new[s, 2, i]
        Action[Env.I_inv_loc[3, s]] = Source.Vd_abc_new[s, 3, i]
    end

    return Action
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

function Measurements(Env::Environment)

    i = Env.steps

    v_loc = [3 11 19; 4 12 20]
    i_loc = [5 13 21; 6 14 22]

    t_final = (Env.N - 1)*Env.μ
    t = 0:Env.μ:t_final # bad code

    num_nodes = size(v_loc, 1)

    T_eval = 1 #number of periods to average over
    i_sp_rms = convert(Int64, round((1/(Env.μ*Env.T_sp_rms))))

    i_start = i - convert(Int64, round(T_eval/(Env.fsys*Env.μ)))

    i_range = i_start:i_sp_rms:i

    for n in 1:num_nodes

        V_poc = [Env.x[v_loc[n,1], i]; Env.x[v_loc[n,2], i]; Env.x[v_loc[n,3], i]] # a, b, c components
        I_poc = [Env.x[i_loc[n,1], i]; Env.x[i_loc[n,2], i]; Env.x[i_loc[n,3], i]]

        Env.p_q_inst[n, :, i] = p_q_theory(V_poc, I_poc)

        Env.p_inst[n, :, i] = V_poc.*I_poc

        # Phasors
        if i_range[1] >= 1 #waiting for one evaluation cycle to pass

                # Voltages
                v_signals = [Env.x[v_loc[n,1], i_range] Env.x[v_loc[n,2], i_range] Env.x[v_loc[n,3], i_range]]

                # Currents
                i_signals = [Env.x[i_loc[n,1], i_range] Env.x[i_loc[n,2], i_range] Env.x[i_loc[n,3], i_range]]

                Env.V_ph[n,  :, :, i] = RMS(Env.θs[i_range], v_signals, Env.T_sp_rms)

                Env.I_ph[n, :, :, i] = RMS(Env.θs[i_range], i_signals, Env.T_sp_rms)

                Env.P[n, 1:3, i] = (Env.V_ph[n, :, 2, i].*Env.I_ph[n, :, 2, i]).*cos.(Env.V_ph[n, :, 3, i] .- Env.I_ph[n, :, 3, i])
                Env.P[n, 4, i] = sum(Env.P[n, 1:3, i])
                Env.Q[n, 1:3, i] = (Env.V_ph[n, :, 2, i].*Env.I_ph[n, :, 2, i]).*sin.(Env.V_ph[n, :, 3, i] .- Env.I_ph[n, :, 3, i])
                Env.Q[n, 4, i] = sum(Env.Q[n, 1:3, i])

        end
    end

    return nothing
end

function Evolution(Env::Environment, Action)

    i = Env.steps
    i_next = i + 1

    Env.u[:, i] = Action

    # Evolving System to next state
    #_______________________________________________________
    k = 1 # number of steps to evolve

    x0 = Env.x[:,i]

    xp = Array{Float64, 1}(undef, size(Env.Ad, 1))
    xp = fill!(xp, 0.0)

    if isempty(Env.Bd)
        xp = xp
    else
        for j in 0:(k-1)
            xp = xp + (Env.Ad^((k-1) - j))*(Env.u[:,i].*Env.Bd)
        end
    end

    Env.x[:, i_next] = (Env.Ad^k)*x0 + xp

    Env.y[:, i_next] = Env.C*Env.x[:, i_next] + Env.u[:,i].*Env.D # Currents and Voltages

    # Evolving System Frequency and Phase
    #_______________________________________________________
    i_start = i_next - 2
    if i_start < 1
        i_start = 1
    end
    i_range = i_start:i_next

    ω = 2*π*Env.fs[i_range]

    Env.θs[i_next] = Third_Order_Integrator(Env.θs[i], Env.μ, ω)%(2*π)

    Env.steps = i_next

    return nothing
end

function Simple_State_Space(u, Lf, Cf, LL, RL)

    a = [[0 -1/Lf 0];
        [1/Cf 0 -1/Cf];
        [0 1/LL -RL/LL]]
    a2 = [[0 -1/Lf 0];
        [1/Cf 0 -1/Cf];
        [0 1/LL -RL*2/LL]]
    A = [a zeros(3,3) zeros(3,3);
        zeros(3,3) a zeros(3,3);
        zeros(3,3) zeros(3,3) a]

    if u == 1
        A = [a2 zeros(3,3) zeros(3,3);
        zeros(3,3) a zeros(3,3);
        zeros(3,3) zeros(3,3) a]
    end

    b = [1/Lf; 0; -1/LL]
    B = [b; b; b]

    c = [[0 -1 0];
        [1 0 -1];
        [0 1 -RL]]
    C = [c zeros(3,3) zeros(3,3);
        zeros(3,3) c zeros(3,3);
        zeros(3,3) zeros(3,3) c]
    d = [1.0; 0.0; -1.0]
    D = [d; d; d]

    return A, B, C, D
end

function Two_Sources_One_Load(Source::Source_Controller, Vo_rms, SL1, pf1, SL2, pf2, Lt1, Lt2, Rt1, Rt2)

    Lf1 = Source.Lf[1]
    Cf1 = Source.Cf[1]

    Lf2 = Source.Lf[2]
    Cf2 = Source.Cf[2]

    RL1, LL, XL, ZL1 = Load_Impedance(SL1, pf1, Vo_rms)
    RL2, CL, XC, ZL2 = Load_Impedance(SL2, pf2, Vo_rms)

    a = [[0 0 -1/Lf1 0 0 0 0 0];
        [0 0 0 -1/Lf2 0 0 0 0];
        [1/Cf1 0 0 0 -1/Cf1 0 0 0];
        [0 1/Cf2 0 0 0 -1/Cf2 0 0];
        [0 0 1/Lt1 0 -Rt1/Lt1 0 -1/Lt1 0];
        [0 0 0 1/Lt2 0 -Rt2/Lt2 -1/Lt2 0];
        [0 0 0 0 1/CL 1/CL -1/(RL2*CL) -1/CL];
        [0 0 0 0 0 0 1/LL -RL1/LL]]

    c = [[0 0 -1 0 0 0 0 0];
        [0 0 0 -1 0 0 0 0];
        [1 0 0 0 -1 0 0 0];
        [0 1 0 0 0 -1 0 0];
        [0 0 1 0 -Rt1 0 -1 0];
        [0 0 0 1 0 -Rt2 -1 0];
        [0 0 0 0 1 1 -1/RL2 -1];
        [0 0 0 0 0 0 -1 RL1]]

    b = [1/Lf1; 1/Lf2; 0; 0; 0; 0; 0; 0]
    d = [1.0; 1.0; 0; 0; 0; 0; 0; 0]

    b2a = [1/Lf1; 0; 0; 0; 0; 0; 0; 0]
    b2b = [0; 1/Lf2; 0; 0; 0; 0; 0; 0]
    b2 = [b2a b2b]

    d2a = [1.0; 0; 0; 0; 0; 0; 0; 0]
    d2b = [0; 1.0; 0; 0; 0; 0; 0; 0]
    d2 = [d2a d2b]

    A = [a zeros(8,8) zeros(8,8);
        zeros(8,8) a zeros(8,8);
        zeros(8,8) zeros(8,8) a]

    B = [b; b; b]
    B2 = [b2; b2; b2]

    C = [c zeros(8,8) zeros(8,8);
        zeros(8,8) c zeros(8,8);
        zeros(8,8) zeros(8,8) c]

    D = [d; d; d]

    D2 = [d2; d2; d2]

    return A, B, C, D, B2, D2
end

function Source_Initialiser(Source::Source_Controller, mode; num_source = 1, Prated = 0, Qrated = 0, Srated = 0, pf = 0.8)

    if Prated == 0 && Qrated == 0 && Srated == 0
        Source.S[num_source] = 50e3
        Source.P[num_source] = 40e3
        Source.Q[num_source] = 30e3
    elseif Prated == 0 && Qrated == 0 && Srated != 0
        Source.S[num_source] = Srated
        Source.P[num_source] = pf*Srated
        Source.Q[num_source] = sqrt(Srated^2 - Source.P[num_source]^2)
    elseif Prated != 0 && Qrated != 0 && Srated == 0
        Source.P[num_source] = Prated
        Source.Q[num_source] = Qrated
        Source.S[num_source] = sqrt(Prated^2 + Qrated^2)
    elseif Prated != 0 && Qrated == 0 && Srated != 0
        Source.S[num_source] = Srated
        Source.P[num_source] = Prated
        Source.Q[num_source] = sqrt(Srated^2 - Prated^2)
    elseif Prated == 0 && Qrated != 0 && Srated != 0
        Source.S[num_source] = Srated
        Source.Q[num_source] = Qrated
        Source.P[num_source] = sqrt(Srated^2 - Qrated^2)
    elseif Prated != 0 && Qrated != 0 && Srated != 0
        Source.S[num_source] = Srated
        Source.P[num_source] = pf*Srated
        Source.Q[num_source] = sqrt(Srated^2 - Source.P[num_source]^2)
    elseif Prated != 0 && Qrated == 0 && Srated == 0
        Source.S[num_source] = Prated/pf
        Source.P[num_source] = Prated
        Source.Q[num_source] = sqrt(Source.S[num_source]^2 - Prated^2)
    end

    Source.Source_Modes[num_source] = mode

    Source.pq0_set[num_source, :] = [Source.P[num_source]; Source.Q[num_source]; 0]

    Source.i_max[num_source] = 1.15*sqrt(2)*Source.S[num_source]/(3*Source.Vrms[num_source])

    Filter_Design(Source, num_source, 0.15, 0.01537)

    Current_PI_LoopShaping(Source, num_source)
    Voltage_PI_LoopShaping(Source, num_source)

    return nothing
end

function Filter_Design(Source::Source_Controller, num_source, ΔILf_ILf, ΔVCf_VCf)

    #=
    The filtering capacitors C should be chosen such that the resonant frequency
    1/sqrt(Ls*C) is approximately sqrt(ωn * ωs), where the ωn is the angular
    frequency of the gird voltage, and ωs is the angular switching frequency
    used to turn on/off the switches)
    =#
    Ir = ΔILf_ILf
    Vr = ΔVCf_VCf

    Sr = Source.S[num_source]
    fs = Source.f_cntr
    Vdc = Source.Vdc[num_source]

    #____________________________________________________________
    # Inductor Design
    Vorms = Source.Vrms[num_source]*1.05
    Vop = Vorms*sqrt(2)

    Zl = 3*Vorms*Vorms/Sr

    Iorms = Vorms/Zl
    Iop = Iorms*sqrt(2)

    ΔIlfmax = Ir*Iop

    Lf = Vdc/(4*fs*ΔIlfmax)

    #____________________________________________________________
    # Capacitor Design
    Vorms = Source.Vrms[num_source]*0.95
    Vop = Vorms*sqrt(2)

    Zl = 3*Vorms*Vorms/Sr

    Iorms = Vorms/Zl
    Iop = Iorms*sqrt(2)
    Ir = Vdc/(4*fs*Lf*Iop)
    ΔIlfmax = Ir*Iop
    ΔVcfmax = Vr*Vop

    Cf = ΔIlfmax/(8*fs*ΔVcfmax)

    fc = 1/(2*π*sqrt(Lf*Cf))

    Source.Cf[num_source] = Cf
    Source.Lf[num_source] = Lf

    return Lf, Cf, fc
end

function Current_PI_LoopShaping(Source::Source_Controller, num_source)

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
    sys_sc = ss(A,B,C,D) # continuous

    Ts = 1/Source.f_cntr
    dly = Source.delay*Ts
    ZoH = tf([1], [Ts/2, 1]) # Transfer function of the sample and hold process
    Pade = tf([-dly/2, 1], [dly/2, 1]) # Pure first-order delay approximation
    PWM_gain = Source.Vdc[num_source]/2

    #SC = tf([1/(1*Lf)], [1, Rf/Lf]) # = tf(sys_sc)
    Gsc_ol = minreal(tf(sys_sc)*Pade*PWM_gain*ZoH) # Full transfer function of plant

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

        if any(real(poles(Gi_cl)) .> 0) == false

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
            println("\nError. PI Current Controller with Positive Poles.")
            println("Source = ", num_source,"\n")
        end
    end

    return nothing
end

function Voltage_PI_LoopShaping(Source::Source_Controller, num_source)

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

        if any(real(poles(Gv_cl)) .> 0) == false

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
            println("\nError. PI Voltage Controller with Positive Poles.")
            println("Source = ", num_source,"\n")
        end
    end

    return nothing
end
