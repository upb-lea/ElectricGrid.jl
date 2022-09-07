mutable struct Source_Controller

    #=
        This object contains the functions and parameters that are relevant for the
        control of a three-phase half bridge DC to AC Converter.
    =#

    Vdc::Vector{Float64}
    Vrms::Vector{Float64}
    S::Vector{Float64} # rated nominal apparent power
    
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

    function Source_Controller(Vdc::Vector{Float64}, Vrms::Vector{Float64},
        S::Vector{Float64},
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
        Vd_abc_new::Array{Float64})

        new(Vdc, Vrms,
        S,
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
        Vd_abc_new)
    end

    function Source_Controller(t_final, f_cntr, num_sources)

        Vdc = Array{Float64, 1}(undef, num_sources)
        Vdc = fill!(Vdc, 800.0)

        Vrms = Array{Float64, 1}(undef, num_sources)
        Vrms = fill!(Vrms, 230.0)

        S = Array{Float64, 1}(undef, num_sources)
        S = fill!(S, 50.0e3)

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

        Source_Controller(Vdc, Vrms,
        S,
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
        Vd_abc_new)
    end
end

mutable struct Environment

    μ::Float64 # s, sampling timestep
    N::Int64 # number of samples

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

    function Environment(μ::Float64, N::Int64,
        fsys::Float64, fs::Vector{Float64}, θs::Vector{Float64}, T_sp_rms::Float64,
        V_ph::Array{Float64}, I_ph::Array{Float64},
        p_q_inst::Array{Float64}, p_inst::Array{Float64}, P::Array{Float64}, Q::Array{Float64},
        A::Matrix{Float64}, B::Vector{Float64}, C::Matrix{Float64}, D::Vector{Float64},
        Ad::Matrix{Float64}, Bd::Vector{Float64},
        x::Matrix{Float64}, y::Matrix{Float64}, u::Matrix{Float64})

        new(μ, N,
        fsys, fs, θs, T_sp_rms, 
        V_ph, I_ph,
        p_q_inst, p_inst, P, Q,
        A, B, C, D,
        Ad, Bd,
        x, y, u)
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

        Environment(μ, N,
        fsys, fs, θs, T_sp_rms, 
        V_ph, I_ph,
        p_q_inst, p_inst, P, Q,
        A, B, C, D,
        Ad, Bd,
        x, y, u)
    end

end

function Phase_Locked_Loop_3ph(Source::Source_Controller, num_source, i, v_abc)

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
    f = Source.fpll[num_source, num_source, cntr_range]
    θ = Source.θpll[num_source, num_source, cntr_range]
    err_t = Source.pll_err_t[num_source, 1]
    err = Source.pll_err[num_source, 1, :]

    n = length(f)

    θ_new = Array{Float64, 1}(undef, 3)
    θ_new = fill!(θ_new, 0)

    Kp = 0.1
    Ki = 10
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

function Current_Controller(Source::Source_Controller, num_source, i, i_abc, θ; Ki = 100.0, Kp = 0.5)

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

    V_new = [0; 0; 0] + Kp*I_err_new .+ Ki*I_err_t_new

    Source.I_err[num_source, :, :] = I_err_int
    Source.I_err_t[num_source, :] = I_err_t_new

    Source.Vd_abc_new[num_source, :, i] = Inv_DQ0_transform(V_new, θ)

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

    I_new = [0; 0; 0] + Kp*V_err_new .+ Ki*V_err_t_new

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

function RMS(θ, t_signals, T_sp)

    # Calcutates the DC offset, RMS magnitude, and phase angle relative to the
    # frequency (θ = 2*π*t[:]) for a three phase system
    i = 1
    i_length = size(t_signals, 1)

    rms = Array{Float64, 2}(undef, 3, 3)
    rms = fill!(rms, 0)

    rm = sqrt(2)

    # Least squares method to find RMS of fundamental
    LS = [ones(i_length, 1) rm*cos.(θ) rm*sin.(θ)]
    LS_inv = pinv(LS)

    for ph in 1:3

        #println(t_signals[ph, i_range])
        AB = LS_inv*t_signals[:, ph]
        A0 = AB[1]
        A1 = AB[2]
        B1 = AB[3]

        a = sqrt(A1^2 + B1^2)

        if B1 < 0 && A1 > 0
            d = acos(B1/a)
        elseif B1 < 0 && A1 < 0
            d = -acos(B1/a)
        elseif B1 > 0 && A1 > 0
            d = acos(B1/a)
        else
            d = -acos(B1/a)
        end

        rms[ph,1] = A0
        rms[ph,2] = a
        rms[ph,3] = d
    end

    return rms
end

function Clarke_Transform(v_abc)
    #= Also known as the alpha-beta-(gamma), αβγ transformation. Implementing
    below is the the power invariant version, such that the matrix is unitary.
    This feature is very suitable when the focus is put on the analysis of
    instantaneous power in three-phase systems. The three-phase instantaneous
    active power has a clear and universally accepted physical meaning and is
    also valid during transient states.

    The transformation maps the three-phase instantaneous voltages/currents in the
    abc phases into the instantaneous voltages/currents on the αβγ axes.

    One advantage of applying the Clark transformation is to separate the zero-
    sequence components from the abc-phase components. The α and β axes make no
    contribution to the zero-sequence components. Typically no zero-sequence currents
    exist in a three-phase three-wire system. If the three-phase voltages are
    balanced in a four-wire system, no zero-sequence voltage is present.

    The transformation suggest an axis transformation. They are stationary axes
    and should not be confused with the concepts of voltage or current phasors.
    Here, instantaneous values of phase voltages and line currents referred to the
    stationary axes transformed into the αβ stationary axes, or vice versa. The
    a, b, and c axes aer spatially shifted by 120 degrees from each other, whereas
    the α and β axes are orthogonal, and the α axis is parallel to the a axis.
    The direction of the β axis was chosen in such a way that if voltage or current
    spatial vectors on the abc coordinates rotate in the abc sequence, they would
    rotate in the αβ sequenc on the αβ coordinates.

    =#
    return v_αβγ = sqrt(2/3)*[1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2; 1/sqrt(2) 1/sqrt(2) 1/sqrt(2)]*v_abc
end

function Inv_Clarke_Transform(v_αβγ)

    return v_abc = sqrt(2/3)*[1 0 1/sqrt(2);-1/2 sqrt(3)/2 1/sqrt(2);-1/2 -sqrt(3)/2 1/sqrt(2)]*v_αβγ
end

function Park_Transform(v_αβγ, θ)
    #=
    The αβγ to dq0 function performs a transformation of αβγ Clarke components
    in a fixed reference frame to dq0 Park components in a rotating reference
    frame. The rotating fame is aligned such that the α-axis and the d-axis are
    parrallel. This type of Park transform is also known as the cosine-based
    Park transformation.

    The position of the rotating frame is given by θ = ωt (where ω represents the
    frame rotation speed), the αβγ to dq0 transformation performs a -ωt rotation
    to the space vector Us = Uα + j*Uβ. The homopolar or zero-sequence component
    remains unchanged. That is, the angle, θ, in the Park transformation is the
    synchronous angular position. It can be made equal to the time integral of the
    fundamental angular frequency ω, determined by a PLL circuit.

    When using PLL circuits and the synchronous reference frames, special care
    should be taken to properly set up the angular position θ = ωt of the Park
    transformation in order to guarantee that the direct-axis current component
    id produces only real power (active power) with the system voltage. In this
    case, the quadrature axis iq produces only imaginary power (reactive power).
    This is the principal reason to set up θ = ωt in phase with the fundamental,
    positive-sequence voltage.
     =#
    return v_dq0 = [cos(θ) sin(θ) 0; -sin(θ) cos(θ) 0; 0 0 1]*v_αβγ
end

function Inv_Park_Transform(v_dq0, θ)
    return v_αβγ = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]*v_dq0
end

function DQ0_transform(v_abc,θ)
    return v_dq0 = Park_Transform(Clarke_Transform(v_abc), θ)
end

function Inv_DQ0_transform(v_dq0, θ)

    return v_abc = Inv_Clarke_Transform(Inv_Park_Transform(v_dq0, θ))
end

function p_q_theory(V_abc, I_abc)

    #=
        For a three-phase system with or without a neutral conductor in the steady-
        state or during transients, the three-phase instantaneous active power
        describes the total instantaneous energy flow per second between two
        subsystems.

        The p-q theory is defined in three-phase system with or without a
        neutral conductor. Three instantaneous powers - the instantaneous zero-sequence
        p0, the instantaneous real power p, and the instantaneous imaginary power q -
        are defined from the instantaneous phase voltages, and line current on the
        αβ0 axes as

        [p0; p; q] = [v0 0 0; 0 vα vβ; 0 vβ -vα]*[i0; iα; iβ]

        The imaginary power q does not contribute to the total energy flow between
        the source and the load and vice versa. The imaginary power q is a new quantity
        and needs a unit to distinguish this power from the traditional reactive power.
        The authors propose the units of "volt-ampere-imaginary" and the symbol [vai].
        The imaginary power q is proportional to the quantity of energy that is being
        exchanged between the phases of the system. It does not contribute to the
        energy transfer between the source and the load at any time. The term
        energy transfer is used here in its strict sense, meaning not only the energy
        delivered to the load, but also the energy oscillating between the source
        and the load.

        It is important to note that the conventional power theory defined reactive
        power as a component of the instantaneous (active) power, which has an average
        value equal to zero. Here, it is not so. The imaginary power means a sum
        of products of instantaneous three-phase voltage and current portions that
        do not contribute to the energy transfer between two subsystems at any time.

        For a linear circuit with sinusoidal voltages and currents, both the instantaneous
        powers aer constant. The real power is equal to the conventional definition of
        the three-phase active power, whereas the imaginary power q is equal to the
        conventional three-phase reactive power. If the load has inductive characteristics,
        the imaginary power has a positive value, and if the load is capacitive, its
        value is negative, in concordance with the most common difinition of reactive power.

        If however a three-phase balanced voltage with a three-phase balanced capacitive
        load is compared to an unbalanced load where the load is unbalanced with just
        one capacitor connected between two phases the difference to conventional theory is
        clearer. In the first case, there is no power flowing from the source to the
        load. Moreover, the imaginary power is constant and coincident with the
        conventional three-phase reactive power. In the second case, the real power
        has no constant part, but an oscillating part. The imaginary power has a constant
        value and an oscillating component. From the conventional power theory, it would
        be normal to expect only a reactive power (average imaginary power) and no real
        power at all. The reason for the difference is because the real power is not zero
        and the capacitor terminal voltage is verying as a sinusoidal wave and, therefore,
        it is being charged and discharged, justifying the energy flow given by p. In fact,
        if it is considered that a turbine is powering the generator, it will have to produce
        a positive torque when the capacitor is charging, or a negative torque when it
        is discharging. Of course, there will be torque ripples in its shaft. In the first
        example, with three balanced capacitors, one capacitor is discharging while
        the others are charging. In steady-stae conditions, there is not total
        three-phase energy flow from the source to the capacitors.

        The imaginary power corresponds to a power that is being exchanged among the
        three phases, without transferring any energy between the source and the load.
        The p-q theory has the prominent merit of allowing complete analysis and real
        time calculation of various powers and respective currents involved in three-
        phase circuits. However, this is not the main point. Knowing in real time the
        values of undesirable currents in a circuit allows us to eliminate them.
        For instance, if the oscillating powers are undesirable, by compensating the
        respective currents of the load and their correspondent currents then the
        compensated currents drawn from the network would be sinusoidal. This is one
        of the basic ideas of active filtering.

        The zero-sequence power has the same characteristics as the instantaneous power
        in a single-phase circuit. It has an average value and an oscillating component
        at twice the line frequency. The average value represents a unidirectional energy
        flow. It has the same characteristics as the conventional (average) active power.
        The oscillating component also transfers energy instantaneously. However, it is has
        an average value equal to zero, because it is oscillating. The analysis shows that,
        in principle, the average value of the zero-sequence power helps to increase the
        total energy transfer, and in this sense, it can be considered as a positive point.
        However, even for the simplest case of a zero-sequence component in the voltage
        and current, the zero-sequence power cannot produce constant power alone. In
        other words, p0 always consists of an oscillating component and a non-zero average
        component. The elimination of the oscillating component is accompanied by the
        elimination of the average component. The zero-sequence power exists only if there
        are zero-sequence voltage and current.

        In summary. The instantaneous powers that the p-q theory defines in the time
        domain are independent of the rms values of voltages and currents. This theory
        includes the conventional frequency-domain concepts of active and reactive
        power defined for three phase sinusoidal balanced systems as a particular case.
        Therefore, the p-q theory in the time domain is not contradictory but
        complementary to the conventional theories in the frequency domain.

        1. Zero-sequence components in the fundamental voltage and current and/or
        in the harmonics do not contribute to hte real power p or the imaginary power q
        2. The total instantaneous energy flow per unit time, that is, the three-phase
        instantaneous active power, even in a distorted and unbalanced system, is
        always equal to the sum of the real and the zero-sequence power (p3Φ = p + p9),
        and may contain several oscillating parts.
        3. The imaginary power q, independent of the presence of harmonic or unbalances,
        represents the energy quantity that is being exchanged between the phases of the
        system. This means that the imaginary power does not contribute to the energy transfer
        between the source and the load at any time. The term "energy transfer" is used
        here in a general manner, referring not only to the energy delivered to the load,
        but also to the energy oscillation between the source and load as well.
    =#

    V_αβγ = Clarke_Transform(V_abc)
    I_αβγ = Clarke_Transform(I_abc)

    return p = [V_αβγ[1] V_αβγ[2] 0; V_αβγ[2] -V_αβγ[1] 0; 0 0 V_αβγ[3]]*I_αβγ
end

function Filter_Design(Source::Source_Controller, num_source, ΔILf_ILf, ΔVCf_VCf)

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

    return Lf, Cf, fc
end

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
    Env.u[1, i_next + Source.delay] = Vo_p*sin(Env.θs[i])#Source.Vd_abc_new[1, 1, i]
    Env.u[9, i_next + Source.delay] = Vo_p*sin(Env.θs[i] - 120*π/180)#Source.Vd_abc_new[1, 2, i]
    Env.u[17, i_next + Source.delay] = Vo_p*sin(Env.θs[i] + 120*π/180)#Source.Vd_abc_new[1, 3, i]
    
    # 2nd Source
    on = 1
    
    Env.u[2, i_next + Source.delay] = on*Vo_p*sin(Env.θs[i])
    Env.u[10, i_next + Source.delay] = on*Vo_p*sin(Env.θs[i] - 120*π/180)
    Env.u[18, i_next + Source.delay] = on*Vo_p*sin(Env.θs[i] + 120*π/180)

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

function Load_Impedance(S, pf, vrms; fsys = 50)

    ω = 2*π*fsys

    if pf > 0

        P = pf*S #W, Active Power - average instance power
        Q = sqrt(S^2 - P^2) #VAr, Reactive Power

        Z = conj(3*vrms^2/(P + 1im*Q)) #Load Impedance

        R = real(Z) # Load Resistor
        X = imag(Z)/ω  # Load Reactor
        L_C = X/ω  # Impedance in series

    else
        P = -pf*S #W, Active Power - average instance power
        Q = -sqrt(S^2 - P^2) #VAr, Reactive Power

        Z = conj(3*vrms^2/(P + 1im*Q))

        Y = 1/Z # Impedance in parallel

        R = real(Y)^-1
        X = imag(Y)^-1

        L_C = 1/(X*ω)

    end

    return R, L_C, X, Z
end

function Two_Sources_One_Load(Lf1, Lf2, Cf1, Cf2, LL, CL, RL1, RL2, Lt1, Lt2, Rt1, Rt2)

    a = [[0 0 -1/Lf1 0 0 0 0 0];
        [0 0 0 -1/Lf2 0 0 0 0];
        [1/Cf1 0 0 0 -1/Cf1 0 0 0];
        [0 1/Cf2 0 0 0 -1/Cf2 0 0];
        [0 0 1/Lt1 0 -Rt1/Lt1 0 -1/Lt1 0];
        [0 0 0 1/Lt2 0 -Rt2/Lt2 -1/Lt2 0];
        [0 0 0 0 1/CL 1/CL -1/(RL2*CL) -1/CL];
        [0 0 0 0 0 0 -1/LL RL1/LL]]

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