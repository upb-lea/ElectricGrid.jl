mutable struct Source_Controller

    #=
        This object contains the functions and parameters that are relevant for the
        control of a three-phase half bridge DC to AC Converter.
    =#

    Vdc::Vector{Float64}
    Vrms::Vector{Float64}
    S::Vector{Float64} # rated nominal apparent power
    P::Vector{Float64} # rated nominal active power
    Q::Vector{Float64} # rated nominal active power

    f_cntr::Float64 # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz
    μ_cntr::Float64 # s, sampling timestep
    N_cntr::Int64
    delay::Int64

    fpll::Array{Float64}
    θpll::Array{Float64}

    V_filt_poc::Array{Float64}
    I_filt_poc::Array{Float64}
    I_filt_inv::Array{Float64}
    p_q_filt::Array{Float64}

    pll_err::Array{Float64}
    pll_err_t::Matrix{Float64}

    I_err::Array{Float64}
    I_err_t::Matrix{Float64}

    V_err::Array{Float64}
    V_err_t::Matrix{Float64}

    I_ref::Array{Float64}
    V_ref::Array{Float64}

    I_ref_dq0::Array{Float64}
    V_ref_dq0::Array{Float64}

    I_dq0::Array{Float64}
    V_dq0::Array{Float64}

    Vd_abc_new::Array{Float64}

    D::Matrix{Float64} # Droop coefficients
    ω_droop::Array{Float64}
    θ_droop::Array{Float64}

    function Source_Controller(Vdc::Vector{Float64}, Vrms::Vector{Float64},
        S::Vector{Float64}, P::Vector{Float64}, Q::Vector{Float64},
        f_cntr::Float64, μ_cntr::Float64, N_cntr::Int64, delay::Int64,
        fpll::Array{Float64}, θpll::Array{Float64},
        V_filt_poc::Array{Float64}, I_filt_poc::Array{Float64},
        I_filt_inv::Array{Float64}, p_q_filt::Array{Float64},
        pll_err::Array{Float64}, pll_err_t::Matrix{Float64},
        I_err::Array{Float64}, I_err_t::Matrix{Float64},
        V_err::Array{Float64}, V_err_t::Matrix{Float64},
        I_ref::Array{Float64}, V_ref::Array{Float64},
        I_ref_dq0::Array{Float64}, V_ref_dq0::Array{Float64},
        I_dq0::Array{Float64}, V_dq0::Array{Float64},
        Vd_abc_new::Array{Float64},
        D::Matrix{Float64}, ω_droop::Array{Float64}, θ_droop::Array{Float64})

        new(Vdc, Vrms,
        S, P, Q,
        f_cntr, μ_cntr, N_cntr, delay,
        fpll, θpll,
        V_filt_poc, I_filt_poc,
        I_filt_inv, p_q_filt,
        pll_err, pll_err_t,
        I_err, I_err_t,
        V_err, V_err_t,
        I_ref, V_ref,
        I_ref_dq0, V_ref_dq0,
        I_dq0, V_dq0,
        Vd_abc_new,
        D, ω_droop, θ_droop)
    end

    function Source_Controller(t_final, f_cntr, num_sources)

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

        μ_cntr = 1/f_cntr
        N_cntr = convert(Int64, floor(t_final/μ_cntr)) + 1
        delay = 1

        fpll = Array{Float64, 3}(undef, num_sources, 3, N_cntr) #3 for 3 phases
        fpll = fill!(fpll, 50.0)

        θpll = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        θpll = fill!(θpll, 0)

        V_filt_poc = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        V_filt_poc = fill!(V_filt_poc, 0)
        I_filt_poc = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        I_filt_poc = fill!(I_filt_poc, 0)
        I_filt_inv = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        I_filt_inv = fill!(I_filt_inv, 0)

        p_q_filt = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        p_q_filt  = fill!(p_q_filt, 0)

        pll_err = Array{Float64, 3}(undef, num_sources, 3, 3) # 3 phases and 3rd order integration
        pll_err = fill!(pll_err, 0)
        pll_err_t = Array{Float64, 2}(undef, num_sources, 3) # PLL total integrated error
        pll_err_t = fill!(pll_err_t, 0)

        # Current Integrations
        I_err = Array{Float64, 3}(undef, num_sources, 3, 3) # 3 phases and 3rd order integration
        I_err = fill!(I_err, 0)
        I_err_t = Array{Float64, 2}(undef, num_sources, 3)
        I_err_t = fill!(I_err_t, 0)

        # Voltage Integrations
        V_err = Array{Float64, 3}(undef, num_sources, 3, 3) # 3 phases and 3rd order integration
        V_err = fill!(V_err, 0)
        V_err_t = Array{Float64, 2}(undef, num_sources, 3)
        V_err_t = fill!(V_err_t, 0)

        I_ref = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        I_ref = fill!(I_ref, 0)
        V_ref = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        V_ref = fill!(V_ref, 0)

        I_ref_dq0 = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        I_ref_dq0 = fill!(I_ref_dq0, 0)
        V_ref_dq0 = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        V_ref_dq0 = fill!(V_ref_dq0, 0)

        I_dq0 = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        I_dq0 = fill!(I_dq0, 0)
        V_dq0 = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        V_dq0 = fill!(V_dq0, 0)

        Vd_abc_new = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        Vd_abc_new = fill!(Vd_abc_new, 0)

        D = Array{Float64, 2}(undef, num_sources, 2)
        D = fill!(D, 0)

        ω_droop = Array{Float64, 3}(undef, num_sources, 3, N_cntr) #3 for 3 phases
        ω_droop = fill!(ω_droop, 50.0*2*π)

        θ_droop = Array{Float64, 3}(undef, num_sources, 3, N_cntr)
        θ_droop = fill!(θ_droop, 5)

        Source_Controller(Vdc, Vrms,
        S, P, Q,
        f_cntr, μ_cntr, N_cntr, delay,
        fpll, θpll,
        V_filt_poc, I_filt_poc,
        I_filt_inv, p_q_filt,
        pll_err, pll_err_t,
        I_err, I_err_t,
        V_err, V_err_t,
        I_ref, V_ref,
        I_ref_dq0, V_ref_dq0,
        I_dq0, V_dq0,
        Vd_abc_new,
        D, ω_droop, θ_droop)
    end
end

mutable struct Environment

    μ::Float64 # s, sampling timestep
    N::Int64 # number of samples

    fsys::Float64 # Mostly for plotting
    fs::Vector{Float64}
    θs::Vector{Float64}
    Δfmax::Float64 # The maximum allowable angular frequency deviation
    ΔEmax::Float64 # The maximum allowable voltage deviation
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

    function Environment(μ::Float64, N::Int64,
        fsys::Float64, fs::Vector{Float64}, θs::Vector{Float64},
        Δfmax::Float64, ΔEmax::Float64, T_sp_rms::Float64,
        V_ph::Array{Float64}, I_ph::Array{Float64},
        p_q_inst::Array{Float64}, p_inst::Array{Float64}, P::Array{Float64}, Q::Array{Float64},
        A::Matrix{Float64}, B::Vector{Float64}, C::Matrix{Float64}, D::Vector{Float64},
        Ad::Matrix{Float64}, Bd::Vector{Float64},
        x::Matrix{Float64}, y::Matrix{Float64}, u::Matrix{Float64},
        V_poc_loc::Matrix{Int64}, I_poc_loc::Matrix{Int64}, I_inv_loc::Matrix{Int64},
        V_load_loc::Matrix{Int64}, I_load_loc::Matrix{Int64},
        num_sources::Int64, num_loads::Int64)

        new(μ, N,
        fsys, fs, θs,
        Δfmax, ΔEmax, T_sp_rms,
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

        N = convert(Int64, floor(t_final/μ)) + 1
        num_nodes = num_sources + num_loads

        fsys = 50.0
        T_sp_rms = 5*fsys #samples in a second for rms calcs, x*fsys = x samples in a cycle

        fs = Array{Float64, 1}(undef, N)
        fs = fill!(fs, fsys)

        θs = Array{Float64, 1}(undef, N)
        θs = fill!(θs, 0)

        Δfmax = 5.0 # Hz
        ΔEmax = 10.0 # V

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

        Environment(μ, N,
        fsys, fs, θs,
        Δfmax, ΔEmax, T_sp_rms,
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

function Phase_Locked_Loop_3ph(Source::Source_Controller, num_source, i, v_abc; Kp = 0.1, Ki = 10)

    #= A robost 3 phase phase locked loop
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
    =#

    cntr_cnt_start = i - 3
    cntr_cnt_end = i - 1
    if cntr_cnt_end < 1
        cntr_cnt_start = 1
        cntr_cnt_end = 1
    end
    if cntr_cnt_start < 1
        cntr_cnt_start = 1
    end
    cntr_range = cntr_cnt_start:cntr_cnt_end

    μ = Source.μ_cntr
    f = Source.fpll[num_source, 1, cntr_range]
    θ = Source.θpll[num_source, 1, cntr_range]
    err_t = Source.pll_err_t[num_source, 1]
    err = Source.pll_err[num_source, 1, :]

    n = length(f)

    θ_new = Array{Float64, 1}(undef, 3)
    θ_new = fill!(θ_new, 0)

    ω = 2*π*f

    err_new = 1*((v_abc[1] - v_abc[2])*cos(-θ[n])
    + (v_abc[3] - v_abc[2])*cos(-θ[n] - 2π/3)) # this is magic

    #v_αβγ = Clarke_Transform(v_abc)
    #err_new = v_αβγ[2]*cos(θ[n]) - v_αβγ[1]*sin(θ[n]) # this also works

    #Low pass filter
    cut_off = 10 #Hz
    rc = 1/(cut_off*2*π)
    α = μ/(rc + μ)
    α = 1

    err_new = err[n] + α*(err_new - err[n]) #err = p_3Φ [Watts]

    err_d = [err[:]; err_new]
    err_int = err_d[2:end]

    err_t_new = Third_Order_Integrator(err_t, μ, err_int) # integration

    f_new = 50 + (Kp*err_new + Ki*err_t_new)

    ω_new = 2*π*f_new

    ω_d = [ω; ω_new]
    if n > 2
        ω_int = ω_d[2:end]
    else
        ω_int = ω_d
    end

    θ_new[1] = Third_Order_Integrator(θ[n], μ, ω_int)%(2*π)
    θ_new[2] = (θ_new[1] - 120*π/180)%(2*π)
    θ_new[3] = (θ_new[1] + 120*π/180)%(2*π)

    Source.fpll[num_source, :, i] = [f_new; f_new; f_new]
    Source.θpll[num_source, :, i] = θ_new
    Source.pll_err_t[num_source, :] = [err_t_new; 0; 0]
    Source.pll_err[num_source, :, :] = [transpose(err_int); transpose(err_int); transpose(err_int)]

    return nothing
end

function PQ_Control(Source::Source_Controller, num_source, i, i_abc, v_abc, θ)

    V_αβγ = Clarke_Transform(v_abc)
    I_αβγ = Clarke_Transform(i_abc)

    pq0 = [V_αβγ[1] V_αβγ[2] 0; V_αβγ[2] -V_αβγ[1] 0; 0 0 V_αβγ[3]]*I_αβγ

    pq0_ref = [40e3; 10e3; 0]

    I_αβγ_ref = Inv_p_q_theory(V_αβγ, pq0_ref)
    I_abc_ref = Inv_Clarke_Transform(I_αβγ_ref)

    I_ref_dq0 = Park_Transform(I_αβγ_ref, θ)

    Source.I_ref_dq0[num_source, :, i] = Park_Transform(I_αβγ_ref, θ)
end

function Droop_Control(Source::Source_Controller, num_source, i, p_q, Env::Environment)

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

    cntr_cnt_start = i - 3
    cntr_cnt_end = i - 1
    if cntr_cnt_end < 1
        cntr_cnt_start = 1
        cntr_cnt_end = 1
    end
    if cntr_cnt_start < 1
        cntr_cnt_start = 1
    end
    cntr_range = cntr_cnt_start:cntr_cnt_end

    μ = Source.μ_cntr
    ω = Source.ω_droop[num_source, 1, cntr_range]
    θ = Source.θ_droop[num_source, 1, cntr_range]

    n = length(ω)

    θ_new = Array{Float64, 1}(undef, 3)
    θ_new = fill!(θ_new, 0)

    ω_new = 50*2*π - p_q[1]*Source.D[num_source, 1]/3
    Source.ω_droop[num_source, :, i] = [ω_new; ω_new; ω_new]

    ω_d = [ω; ω_new]
    if n > 2
        ω_int = ω_d[2:end]
    else
        ω_int = ω_d
    end

    Source.θ_droop[num_source, 1, i] = Third_Order_Integrator(θ[end], μ, ω_int)%(2*π)
    Source.θ_droop[num_source, 2, i] = (Source.θ_droop[num_source, 1, i] - 120*π/180)%(2*π)
    Source.θ_droop[num_source, 3, i] = (Source.θ_droop[num_source, 1, i] + 120*π/180)%(2*π)

    E_new = Source.Vrms[num_source] - p_q[2]*Source.D[num_source, 2]/3

    Source.V_ref[num_source, 1, i] = sqrt(2)*E_new*sin(Source.θ_droop[num_source, 1, i])
    Source.V_ref[num_source, 2, i] = sqrt(2)*E_new*sin(Source.θ_droop[num_source, 2, i])
    Source.V_ref[num_source, 3, i] = sqrt(2)*E_new*sin(Source.θ_droop[num_source, 3, i])
end

function Current_Controller(Source::Source_Controller, num_source, i, i_abc, θ; Ki = 10.0, Kp = 5)

    # --- to be removed
    #Source.I_ref_dq0[num_source, :,i] = DQ0_transform(Source.I_ref[num_source, :, i], θ)
    # --- =#
    Source.I_ref[num_source, :, i] = Inv_DQ0_transform(Source.I_ref_dq0[num_source, :,i], θ)

    Source.I_dq0[num_source, :, i] = DQ0_transform(i_abc, θ)

    I = Source.I_dq0[num_source, :, i]
    I_ref = Source.I_ref_dq0[num_source, :,i]
    I_err = Source.I_err[num_source, :, :]
    I_err_t = Source.I_err_t[num_source, :]
    μ = Source.μ_cntr

    n = size(I_err, 2)

    I_err_t_new = Array{Float64, 1}(undef, 3)
    I_err_t_new = fill!(I_err_t_new, 0)

    I_err_new = I_ref .- I

    I_err_d = [I_err I_err_new]

    I_err_int = I_err_d[:, 2:end]

    for i in 1:3
        I_err_t_new[i] = Third_Order_Integrator(I_err_t[i], μ, I_err_int[i,:]) # integration
        #I_err_t_new[i] = I_err_t[i] + μ*I_err_int[i, end] # integration
    end

    s_dq0_avg = Kp.*I_err_new/1000 .+ Ki.*I_err_t_new/100
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

    Source.I_err[num_source, :, :] = I_err_int
    Source.I_err_t[num_source, :] = I_err_t_new

    Source.Vd_abc_new[num_source, :, i] = 0.5*Source.Vdc[num_source].*Inv_DQ0_transform(s_dq0_avg, θ)

    return nothing
end

function Voltage_Controller(Source::Source_Controller, num_source, i, v_abc, θ; Ki = 100.0, Kp = 0.01)

    Rf_V = [Source.V_ref[num_source, 1, i];
            Source.V_ref[num_source, 2, i];
            Source.V_ref[num_source, 3, i]]

    Source.V_ref_dq0[num_source, :,i] = DQ0_transform(Rf_V, θ)
    Source.V_dq0[num_source, :, i] = DQ0_transform(v_abc, θ)

    V = Source.V_dq0[num_source, :, i]
    V_ref = Source.V_ref_dq0[num_source, :, i]
    V_err = Source.V_err[num_source, :, :]
    V_err_t = Source.V_err_t[num_source, :]
    μ = Source.μ_cntr

    n = size(V_err,2)

    V_err_t_new = Array{Float64, 1}(undef, 3)
    V_err_t_new = fill!(V_err_t_new, 0)

    V_err_new = V_ref .- V

    V_err_d = [V_err V_err_new]

    V_err_int = V_err_d[:, 2:end]

    for i in 1:3
        V_err_t_new[i] = Third_Order_Integrator(V_err_t[i], μ, V_err_int[i,:]) # integration
        #V_err_t_new[i] = V_err_t[i] + μ*V_err_int[i, end] # integration
    end

    I_new = [0; 0; 0] .+ Kp.*V_err_new .+ Ki.*V_err_t_new

    Source.I_ref_dq0[num_source, :, i] = I_new
    Source.V_err[num_source, :, :] = V_err_int
    Source.V_err_t[num_source, :] = V_err_t_new

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

    y_new = -1*bot_z_1.*y[:,n] .- bot_z_2*y[:,n-1] .+
            top_z_0*x[:,k] .+ top_z_1*x[:,k-1] .+ top_z_2*x[:,k-2]

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

#-------------------------------------------------------------------------------

function Measurements(Env::Environment, i)

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

function Evolution(Env::Environment, Source::Source_Controller, i)

    i_next = i + 1
    Vo_p = 230*sqrt(2)

    # Inverter Voltages - Control Actions
    #_______________________________________________________
    Env.u[1, i_next + Source.delay] = Source.Vd_abc_new[1, 1, i] # Vo_p*sin(Env.θs[i]) #
    Env.u[9, i_next + Source.delay] = Source.Vd_abc_new[1, 2, i] # Vo_p*sin(Env.θs[i] - 120*π/180) #
    Env.u[17, i_next + Source.delay] = Source.Vd_abc_new[1, 3, i] # Vo_p*sin(Env.θs[i] + 120*π/180) #

    # 2nd Source

    Env.u[2, i_next + Source.delay] = Source.Vd_abc_new[2, 1, i] # Vo_p*sin(Env.θs[i]) #
    Env.u[10, i_next + Source.delay] = Source.Vd_abc_new[2, 2, i] # Vo_p*sin(Env.θs[i] - 120*π/180) #
    Env.u[18, i_next + Source.delay] = Source.Vd_abc_new[2, 3, i] # Vo_p*sin(Env.θs[i] + 120*π/180) #

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

function Two_Sources_One_Load(Lf1, Lf2, Cf1, Cf2, LL, CL, RL1, RL2, Lt1, Lt2, Rt1, Rt2)

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

    A = [a zeros(8,8) zeros(8,8);
        zeros(8,8) a zeros(8,8);
        zeros(8,8) zeros(8,8) a]

    B = [b; b; b]

    C = [c zeros(8,8) zeros(8,8);
        zeros(8,8) c zeros(8,8);
        zeros(8,8) zeros(8,8) c]

    D = [d; d; d]

    return A, B, C, D
end

function Source_Initialiser(Source::Source_Controller, Env::Environment; num_source = 1, Prated = 0, Qrated = 0, Srated = 0)

    if Prated == 0 && Qrated == 0 && Srated == 0
        Source.S[num_source] = 50e3
        Source.P[num_source] = 40e3
        Source.Q[num_source] = 30e3
    elseif Prated == 0 && Qrated == 0 && Srated != 0
        Source.S[num_source] = Srated
        Source.P[num_source] = 0.8*Srated
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
        println("\nError. Too many arguments specified")
    end

    Source.D[num_source, 1] = 2*π*Env.Δfmax/Source.P[num_source]
    Source.D[num_source, 2] = Env.ΔEmax/Source.Q[num_source]
end
