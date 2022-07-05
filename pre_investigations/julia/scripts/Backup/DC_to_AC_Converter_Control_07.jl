using Plots
using LinearAlgebra
using FFTW

function Discrete_time(Ad, Bd, u, k, x0)

    xp = Array{Float64, 1}(undef, size(Ad,1))
    xp = fill!(xp, 0.0)

    if isempty(Bd) || u == 0
        xp = xp
    else
        for j in 0:(k-1)
            xp = xp + (Ad^((k-1) - j))*(u.*Bd)
        end
    end

    x = (Ad^k)*x0 + xp

    # Brute force - not discrete time matrices
    #dx_dt = A*x0+ u*B
    #x = x0 + μ.*dx_dt

    return x
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

function Phase_Locked_Loop(v_abc, t, μ, f, θ, err_t, err)

    #=
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

    n = size(f,2)

    θ_new = Array{Float64, 1}(undef, 3)
    θ_new = fill!(θ_new, 0)

    if n == 0
        f = 50
        θ = [0; 0; 0] #90*π/180
        n = 1
    end

    Kp = 0.01
    Ki = 1
    ω = 2*π*f[1]

    err_new = 1*((v_abc[1] - v_abc[2])*cos(-θ[1, n])
    + (v_abc[3] - v_abc[2])*cos(-θ[1, n] - 2π/3)) # this is magic

    #v_αβγ = Clarke_Transform(v_abc)
    #err_new = v_αβγ[2]*cos(θ[n]) - v_αβγ[1]*sin(θ[n]) # this also works

    #Low pass filter
    cut_off = 10 #Hz
    rc = 1/(cut_off*2*π)
    α = μ/(rc + μ)
    α = 1

    err_new = err[1, n] + α*(err_new - err[1, n]) #err = p_3Φ [Watts]

    err_d = [err[1,:]; err_new]
    err_int = err_d[2:end]

    err_t_new = Third_Order_Integrator(err_t[1], μ, err_int) # integration

    f_new = 50 + (Kp*err_new + Ki*err_t_new)

    ω_new = 2*π*f_new

    ω_d = [ω; ω_new]
    if n > 2
        ω_int = ω_d[2:end]
    else
        ω_int = ω_d
    end

    θ_new[1] = Third_Order_Integrator(θ[1, n], μ, ω_int)%(2*π)
    θ_new[2] = (θ_new[1] - 120*π/180)%(2*π)
    θ_new[3] = (θ_new[1] + 120*π/180)%(2*π)

    return [f_new; f_new; f_new], θ_new, [err_t_new; 0; 0], [transpose(err_int); transpose(err_int); transpose(err_int)]
end

function RMS(rms, θ, t_signals, T_sp)

    # Calcutates the DC offset, RMS magnitude, and phase angle relative to the
    # frequency (θ = 2*π*t[:]) for a three phase system
    i = 1
    i_length = size(t_signals,1)

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

        rms[ph,1,i] = A0
        rms[ph,2,i] = a
        rms[ph,3,i] = d
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

function Current_Controller(I, I_ref, I_err, I_err_t, μ)

    n = size(I_err,2)

    I_err_t_new = Array{Float64, 1}(undef, 3)
    I_err_t_new = fill!(I_err_t_new, 0)

    Ki = 50
    Kp = 1

    I_err_new = I_ref .- I

    #Low pass filter
    cut_off = 800 #Hz
    rc = 1/(cut_off*2*π)
    α = μ/(rc + μ)
    #α = 1
    I_err_new = I_err[:,n] .+ α*(I_err_new .- I_err[:,n])

    I_err_d = [I_err I_err_new]
    I_err_int = I_err_d[:, 2:end]

    for i in 1:3
        I_err_t_new[i] = Third_Order_Integrator(I_err_t[i], μ, I_err_int[i,:]) # integration
        #I_err_t_new[i] = I_err_t[i] + μ*I_err_int[i, end] # integration
    end

    V_new = [0; 0; 0] + Kp*I_err_new .+ Ki*I_err_t_new

    return V_new, I_err_int, I_err_t_new
end

function Voltage_Controller(V, V_ref, V_err, V_err_t, μ)

    n = size(V_err,2)

    V_err_t_new = Array{Float64, 1}(undef, 3)
    V_err_t_new = fill!(V_err_t_new, 0)

    Ki = 30
    Kp = 1.2

    V_err_new = V_ref .- V

    #Low pass filter
    fc = 1000 #Hz
    rc = 1/(fc*2*π)
    α = μ/(rc + μ)
    #α = 1
    V_err_new = V_err[:,n] .+ α*(V_err_new .- V_err[:,n])
    V_err_new[1] = First_Order_LPF(fc, V_err[1,(n-1):n], V_err_new[1], μ)
    V_err_new[2] = First_Order_LPF(fc, V_err[2,(n-1):n], V_err_new[2], μ)
    V_err_new[3] = First_Order_LPF(fc, V_err[3,(n-1):n], V_err_new[3], μ)

    V_err_d = [V_err V_err_new]
    V_err_int = V_err_d[:, 2:end]

    for i in 1:3
        V_err_t_new[i] = Third_Order_Integrator(V_err_t[i], μ, V_err_int[i,:]) # integration
        #V_err_t_new[i] = V_err_t[i] + μ*V_err_int[i, end] # integration
    end

    I_new = [0; 0; 0] + Kp*V_err_new .+ Ki*V_err_t_new

    return I_new, V_err_int, V_err_t_new
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
        x = [x 0]
        k = 2
    end

    #apply pre-warping transformation
    ωa = (2/μ)*tan(fc*2π*μ/2)

    #=
        Finding Coefficients of Transfer/Pulse function G(z)
        G(z) = (ωd*z^-1 + ωd)/((ωd - 1)*z^-1 + 1)
    =#

    ωd = ωa*μ/2

    y_new = (1 - ωd)*y + ωd*x[k] + ωd*x[k-1]

    return y_new
end

mutable struct Source_Controller

    #=
        This object contains the functions and parameters that are relevant for the
        control of a three-phase half bridge DC to AC Converter.
    =#

    Vdc::Float64
    f_cntr::Float64 # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz
    μ_cntr::Float64 # s, sampling timestep
    N_cntr::Int64

    fpll::Matrix{Float64}
    θpll::Matrix{Float64}

    V_filt_poc::Matrix{Float64}
    I_filt_poc::Matrix{Float64}
    I_filt_inv::Matrix{Float64}
    p_q_filt::Matrix{Float64}

    pll_err::Matrix{Float64}
    pll_err_t::Vector{Float64}

    I_err::Matrix{Float64}
    I_err_t::Vector{Float64}

    V_err::Matrix{Float64}
    V_err_t::Vector{Float64}

    function Source_Controller(Vdc::Float64,
        f_cntr::Float64, μ_cntr::Float64, N_cntr::Int64,
        fpll::Matrix{Float64}, θpll::Matrix{Float64},
        V_filt_poc::Matrix{Float64}, I_filt_poc::Matrix{Float64},
        I_filt_inv::Matrix{Float64}, p_q_filt::Matrix{Float64},
        pll_err::Matrix{Float64}, pll_err_t::Vector{Float64},
        I_err::Matrix{Float64}, I_err_t::Vector{Float64},
        V_err::Matrix{Float64}, V_err_t::Vector{Float64})

        new(Vdc,
        f_cntr, μ_cntr, N_cntr,
        fpll, θpll,
        V_filt_poc, I_filt_poc,
        I_filt_inv, p_q_filt,
        pll_err, pll_err_t,
        I_err, I_err_t,
        V_err, V_err_t)
    end

    function Source_Controller()

        Vdc = 800.0

        t_final = 1 # second

        f_cntr = 15e3 # Hz
        μ_cntr = 1/f_cntr
        N_cntr = convert(Int64, floor(t_final/μ_cntr)) + 1

        fpll = Array{Float64, 2}(undef, 3, N_cntr)
        fpll = fill!(fpll, 50.0)

        θpll = Array{Float64, 2}(undef, 3, N_cntr)
        θpll = fill!(θpll, 0)

        V_filt_poc = Array{Float64, 2}(undef, 3, N_cntr)
        V_filt_poc = fill!(V_filt_poc, 0)
        I_filt_poc = Array{Float64, 2}(undef, 3, N_cntr)
        I_filt_poc = fill!(I_filt_poc, 0)
        I_filt_inv = Array{Float64, 2}(undef, 3, N_cntr)
        I_filt_inv = fill!(I_filt_inv, 0)

        p_q_filt = Array{Float64, 2}(undef, 3, N_cntr)
        p_q_filt  = fill!(p_q_filt, 0)

        pll_err = Array{Float64, 2}(undef, 3, 3) # 3 phases and 3rd order integration
        pll_err = fill!(pll_err, 0)
        pll_err_t = Array{Float64, 1}(undef, 3) # PLL total integrated error
        pll_err_t = fill!(pll_err_t, 0)

        # Current Integrations
        I_err = Array{Float64, 2}(undef, 3, 3) # 3 phases and 3rd order integration
        I_err = fill!(I_err, 0)
        I_err_t = Array{Float64, 1}(undef, 3)
        I_err_t = fill!(I_err_t, 0)

        # Voltage Integrations
        V_err = Array{Float64, 2}(undef, 3, 3) # 3 phases and 3rd order integration
        V_err = fill!(V_err, 0)
        V_err_t = Array{Float64, 1}(undef, 3)
        V_err_t = fill!(V_err_t, 0)

        Source_Controller(Vdc,
        f_cntr, μ_cntr, N_cntr,
        fpll, θpll,
        V_filt_poc, I_filt_poc,
        I_filt_inv, p_q_filt,
        pll_err, pll_err_t,
        I_err, I_err_t,
        V_err, V_err_t)
    end
end

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

fsys = 50 # the frequency of the reference/control waveform
Vo_rms = 230 # rms output voltage
Vo_p = sqrt(2)*Vo_rms

Vdc = 800
#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 10 #time step in μs
t_final = 1 #time in seconds, total simulation run time

μ_cntr = 10e-6#6e-5 # s, 10e-4 #sampling timestep
f_cntr = 1/μ_cntr # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz
#f_cntr = 15e3 # Hz, Sampling/Switching frequency of controller ~ 15 kHz -> 50kHz

# Impedance Calculations
SL = 50e3 #VA, 3-ph Apparent Power
pf = 0.6 #power factor

ΔILf_ILf = 0.15

ΔVCf_VCf = 0.01537

#_______________________________________________________________________________
# Circuit Elements Calcs
PL = pf*SL #W, Active Power - average instance power
QL = sqrt(SL^2 - PL^2) #VAr, Reactive Power

RL = real(3*Vo_rms^2/(PL + 1im*QL)) # Load/Line Resistor
LL = -1*imag(3*Vo_rms^2/(PL + 1im*QL))/(2*π*fsys) # Load/Line Reactor
XL = (2*π*fsys)*LL
ZL = RL + 1im*XL #Load Impedance

# Reference currents
I_ref = Vo_rms/ZL
I_ref_p = abs(I_ref)*sqrt(2) # peak reference current
I_ref_rms = abs(I_ref) # peak reference current
I_ref_ang = angle(I_ref)

# Inv 1 Source Filter Calcs
io_p = Vo_p/abs(ZL) # the amplitude of the load current
ΔILf_max = ΔILf_ILf*I_ref_p
ΔVCf_max = ΔVCf_VCf*Vo_p

Lf = Vdc/(4*f_cntr*ΔILf_max) # Filter Inductor
Cf = ΔILf_max/(8*f_cntr*ΔVCf_max) # Filter Capacitor
fc = 1/(2*π*sqrt(Lf*Cf)) #Filter Cut-off frequency
#Lf = 0.0026 # Filter Inductor
#Cf = 7.71e-5 # Filter Capacitor
#_______________________________________________________________________________
# Simulation Calcs
Nps = (1/(Timestep*1e-6)) # time intervals, samples in a second + correction
μps = 1/Nps #corrected time step
t = 0:μps:t_final # time

N = length(t) # number of samples
#_______________________________________________________________________________
# State space representation
a = [[0 -1/Lf 0];
    [1/Cf 0 -1/Cf];
    [0 1/LL -RL/LL]]
a2 = [[0 -1/Lf 0];
    [1/Cf 0 -1/Cf];
    [0 1/LL -RL*2/LL]]
A = [a2 zeros(3,3) zeros(3,3);
    zeros(3,3) a zeros(3,3);
    zeros(3,3) zeros(3,3) a]

b = [1/Lf; 0; -1/LL]
B = [b; b; b]

c = [[0 -1 0];
    [1 0 -1];
    [0 1 -RL]]
C = [c zeros(3,3) zeros(3,3);
    zeros(3,3) c zeros(3,3);
    zeros(3,3) zeros(3,3) c]
d = [1; 0; -1]
D = [d; d; d]

# Discrete time matrices
Ad = exp(A*μps)
Bd = inv(A)*(Ad - Matrix(I, size(A,1), size(A,1)))*B

# state space variables
x = Array{Float64, 2}(undef, size(A,1), N)
x = fill!(x, 0)

# x = [iLfa; vCfa; iLLa; ...]

#_______________________________________________________________________________
# Control System Variables
# Control Calcs
μ_cntr = 1/f_cntr # s, 10e-4 #sampling timestep

#= Switching and Sampling Frequencies / Dead time in digital control loops

    Inverter controls are realized with the help of microcontrollers or DSP technologies.
    These controls rely on recurse control rules which are executed in a discrete-time
    manner (sampling), i.e. new reference values are provided only at a certain clock
    cycle. Since the controller can influence the inverter currents only via the switching
    elements, it is not reasonable from an economic point of view to execute the control
    cycle more often that the inverter's maximum switching frequency. Therefore, the
    sampling time μ_cntr is synchronised with the switching period μ_sw, and we can apply
    two strategies μ_cntr = μ_sw or μ_cntr = 0.5*μ_cntr. Either way, the transfer function
    of the sample and hold process is (despite the fact that no continuous-time reference
    value exists that could be samples, since the discrete-time reference falue is
    determined directly from a discrete-time controller):

    H(s) = (1/(s*μ_cntr))*(1 - exp(-s*μ_cntr))

    Hence the discrete-time reference value can be approximated either by a dead time
    delay element, whereas in bothe cases the significant time constant μ_cntr/2 occurs.

    When we assume a purely inductive load with a constant reverse voltage, a synchronisation
    of the current sampling with the PWM carrier causes the sampled current values to
    exactly represent the time-local avearge current values. The current ripple caused
    by the inverter pulsation is masked by this subtly sampling process making analog
    filtering reedundant. On the contrary, the applicatin of an analog pre-filter would
    cause the current sampling method to no longer work as originally intended. Although
    the harmonic pulsation components would be smoothed, the pre-filter would cause a
    phase shift in the fundamental component, leading to errors in the control loop.

    Dead time in digital control loops:
    If the control scheme is implemented on a microcontroller or microprocessor, then
    a certain time is required to process the control algorithm. Therefore, a measured value
    can affect the voltage reference only after this time period has passed. In an appropriate
    manner, all these processes are synchronised within the clock cycle given by the
    pulse width modulation or vector modulation. This way, the digital control introduces
    a dead time of one sampling step. Together with the discrete time application of the
    reference voltage for PWM a total dead time of 1.5 sampling steps of the control loop
    results.

    At varying DC input voltages the dead time related to the discrete-time processes
    causes a further problem. Both, for the vector- and PWM, the voltage reference
    v* initially needs to be referred to the input voltage vdc. There, the DC voltage
    value, which is known by the time when the reference is applied, is used. If the input
    voltage has been sampled synchronously together with the current measurements,
    then the voltage value at the previous sampling instant is given as:

    average(v*[k]) = v*[k]/(vdc[k - 1]/2)

    However, this normalized voltage reference is converted into a corresponding pulse
    sequence not before the next sampling interval. If the input voltage does not change
    or only slowly changes until that point, no problem occurs. When dealing with highly
    fluctuating input voltages, however, a voltage error is caused during the PWM
=#

cntr_cnt = 0
t_cntr = 0:μ_cntr:t_final # time
inter_length = 2 # number of samples to integrate over
i_sp = convert(Int64, round(μ_cntr/μps))
N_cntr = convert(Int64, floor(t_final/μ_cntr)) + 1

#_______________________________________________________________________________
# output variables

y = Array{Float64, 2}(undef, size(A, 1), N)
y = fill!(y, 0)
# y = [vLfa; iCfa; vLLa; ...]

p_inv = Array{Float64, 2}(undef, 3, N) # instantaneous power at PCC
p_net = Array{Float64, 2}(undef, 3, N) # instantaneous power at load/network

# initial conditions
x[1,1] = 0 #iLfa
x[2,1] = 0 #vCfa
x[3,1] = 0 #iLLa
x[4,1] = 0 #iLfb
x[5,1] = 0 #vCfb
x[6,1] = 0 #iLLb
x[7,1] = 0 #iLfc
x[8,1] = 0 #vCfc
x[9,1] = 0 #iLLc

# Inverter Initial Output Values and state placeholders
Vd_a = 0
Vd_b = 0
Vd_c = 0

#_______________________________________________________________________________
# Control Reference signals

Rf_Va = Array{Float64, 1}(undef, N)
Rf_Vb = Array{Float64, 1}(undef, N)
Rf_Vc = Array{Float64, 1}(undef, N)

Rf_Ia = Array{Float64, 1}(undef, N)
Rf_Ib = Array{Float64, 1}(undef, N)
Rf_Ic = Array{Float64, 1}(undef, N)

fs = Array{Float64, 1}(undef, N)
fs = fill!(fs, fsys)
Phase = Array{Float64, 1}(undef, N)
Phase = fill!(Phase, 0)

Vs_a = Array{Float64, 1}(undef, N)
Vs_b = Array{Float64, 1}(undef, N)
Vs_c = Array{Float64, 1}(undef, N)

Vd_a = Array{Float64, 1}(undef, N)
Vd_b = Array{Float64, 1}(undef, N)
Vd_c = Array{Float64, 1}(undef, N)

I_ref_dq0 = Array{Float64, 2}(undef, 3, N_cntr)

I_dq0 = Array{Float64, 2}(undef, 3, N_cntr)

V_ref_dq0 = Array{Float64, 2}(undef, 3, N_cntr)

V_dq0 = Array{Float64, 2}(undef, 3, N_cntr)
Vd_dq0_new = Array{Float64, 2}(undef, 3, N_cntr)

#_______________________________________________________________________________
# Control System Integrator errors

pll_err = Array{Float64, 2}(undef, 3, 3) #PLL errors - 3 for 3rd order integrations
pll_err = fill!(pll_err, 0)
pll_err_t = Array{Float64, 1}(undef, 3) #3 phases - PLL total integrated error
pll_err_t = fill!(pll_err_t, 0) #PLL total integrated error

# Current Integrations
I_err = Array{Float64, 2}(undef, 3, 3) # 3 phases and 3rd order integration
I_err = fill!(I_err, 0)
I_err_t = Array{Float64, 1}(undef, 3)
I_err_t = fill!(I_err_t, 0)

# Voltage Integrations
V_err = Array{Float64, 2}(undef, 3, 3) # 3 phases and 3rd order integration
V_err = fill!(V_err, 0)
V_err_t = Array{Float64, 1}(undef, 3)
V_err_t = fill!(V_err_t, 0)

#_______________________________________________________________________________
# Measurements

# RMS Phasors
Vinv_ph = Array{Float64, 3}(undef, 3, 3, N)
Vinv_ph = fill!(Vinv_ph, 0)
Iinv_ph = Array{Float64, 3}(undef, 3, 3, N)
Iinv_ph = fill!(Iinv_ph, 0)

# Instantaneous Real, Imaginary, and Zero powers
p_q_inst = Array{Float64, 2}(undef, 3, N)
p_q_inst = fill!(p_q_inst, 0)

# Phase Locked Loop Measurements
fpll = Array{Float64, 2}(undef, 3, N_cntr)
fpll = fill!(fpll, 50)
θpll = Array{Float64, 2}(undef, 3, N_cntr)
θpll = fill!(θpll, 0)

V_filt_poc = Array{Float64, 2}(undef, 3, N_cntr)
V_filt_poc = fill!(V_filt_poc, 0)
I_filt_poc = Array{Float64, 2}(undef, 3, N_cntr)
I_filt_poc = fill!(I_filt_poc, 0)
I_filt_inv = Array{Float64, 2}(undef, 3, N_cntr)
I_filt_inv = fill!(I_filt_inv, 0)

p_q_filt = Array{Float64, 2}(undef, 3, N_cntr)
p_q_filt  = fill!(p_q_filt, 0)

#%% Starting time simulation
println("\nHere we go.\n")

println("Progress : 0.0 %")

for i in 1:N-1

    #----
    if t[i] > 40
        fs[i] = 45 + 0.0*sin(2*π*10*t[i])
        Phase[i] = 10*π/180
    else
        fs[i] = 50
        Phase[i] = 0 #210*π/180#Phase[i] + 50*π/180
    end

    if t[i] < -0.1
        v_ref = Vo_p*t[i]/0.1
    else
        v_ref = Vo_p
    end
    #----

    if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
        flush(stdout)
        println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
    end

    i_next = i + 1

    global pll_err, pll_err_t, cntr_cnt
    global I_err, I_err_t
    global V_err, V_err_t
    global Vd_a, Vd_b, Vd_c

    global Vo_p, Phase

    # Control Reference Values
    #___________________________________________________________________________

    Rf_Va[i] = v_ref*sin(2*π*fs[i]*t[i] + Phase[i])
    Rf_Vb[i] = v_ref*sin(2*pi*fs[i]*t[i] - 120*π/180 + Phase[i])
    Rf_Vc[i] = v_ref*sin(2*pi*fs[i]*t[i] + 120*π/180 + Phase[i])

    # --- to be removed
    Rf_Ia[i] = I_ref_p*sin(2*π*fs[i]*t[i] + I_ref_ang)
    Rf_Ib[i] = I_ref_p*sin(2*pi*fs[i]*t[i] - 120*π/180 + I_ref_ang)
    Rf_Ic[i] = I_ref_p*sin(2*pi*fs[i]*t[i] + 120*π/180 + I_ref_ang)
    Rf_I = [Rf_Ia[i]; Rf_Ib[i]; Rf_Ic[i]]
    # ----

    # Environment Measurments
    #___________________________________________________________________________
    # System State and Powers
    V_poc = [x[2, i]; x[5, i]; x[8, i]] # a, b, c components
    I_poc = [x[3, i]; x[6, i]; x[9, i]]
    I_inv = [x[1, i]; x[4, i]; x[7, i]]

    # Instantaneous Power Calculations
    p_q_inst[:, i] = p_q_theory(V_poc, I_poc)
    #s_inst = p_q_inst[1] + 1im*p_q_inst[2]

    p_inv[1, i] = x[2, i]*x[3, i]
    p_inv[2, i] = x[5, i]*x[6, i]
    p_inv[3, i] = x[8, i]*x[9, i]

    p_net[1, i] = Vs_a[i]*x[3, i]
    p_net[2, i] = Vs_b[i]*x[6, i]
    p_net[3, i] = Vs_c[i]*x[9, i]

    # Phasor Outputs

    T_eval = 1 #number of periods to average over
    T_sp_rms = 10*fsys #samples in a second, x*fsys = x samples in a cycle
    i_start = i - convert(Int64, round(T_eval/(fsys*μps)))
    local i_sp_rms = convert(Int64, round((1/(μps*T_sp_rms))))

    i_range = i_start:i_sp_rms:i

    if i_range[1] >= 1 #waiting for one evaluation cycle to pass
        global T_sp_rms

        # Voltages
        v_signals = [x[2, i_range] x[5, i_range] x[8, i_range]]
        #v_signals = [Rf_a[i_range] Rf_b[i_range] Rf_c[i_range]]

        # Currents
        i_signals = [x[3, i_range] x[6, i_range] x[9, i_range]]

        Vinv_ph[ :, :, i_next] = RMS(Vinv_ph[ :, :, i], 2*pi*fs[i_range].*t[i_range],
        v_signals, T_sp_rms)

        Iinv_ph[ :, :, i_next] = RMS(Iinv_ph[ :, :, i], 2*pi*fs[i_range].*t[i_range],
        i_signals, T_sp_rms)
    end

    # System Dynamics
    #___________________________________________________________________________
    # Electrical
    if i > 1
        u = [Vd_a[i-1]; 0; Vs_a[i-1]; Vd_b[i-1]; 0; Vs_b[i-1]; Vd_c[i-1]; 0; Vs_c[i-1]]
    else
        u = [0; 0; 0; 0; 0; 0; 0; 0; 0]
    end

    x[:, i_next] = Discrete_time(Ad, Bd, u, 1, x[:,i])

    #  Outputs
    y[:, i_next] = C*x[:, i_next] + u.*D # Currents and Voltages

    #___________________________________________________________________________
    # Control
    if floor(i*μps/μ_cntr) - cntr_cnt > 0 && i >= i_sp

        cntr_cnt_start = cntr_cnt - inter_length
        if cntr_cnt_start < 1
            cntr_cnt_start = 1
        end
        cntr_range = cntr_cnt_start:cntr_cnt
        cntr_cnt_next = cntr_cnt + 1

        # Measurements and Filters
        #_______________________________________________________________________

        x_start = i - 2*i_sp
        if x_start < 1
            x_start = 1
        end
        x_range = x_start:i_sp:i

        y_start = cntr_cnt - 1
        y_end = cntr_cnt
        if y_start < 1
            y_start = 1
            y_end = 1
        end
        y_range = y_start:y_end

        V_poc = transpose([x[2, x_range] x[5, x_range] x[8, x_range]]) # a, b, c components
        I_poc = transpose([x[3, x_range] x[6, x_range] x[9, x_range]])
        I_inv = transpose([x[1, x_range] x[4, x_range] x[7, x_range]])

        V_filt_poc[:,cntr_cnt_next] = Butterworth_LPF(800, V_poc, V_filt_poc[:,y_range], μ_cntr)
        I_filt_poc[:,cntr_cnt_next] = Butterworth_LPF(800, I_poc, I_filt_poc[:,y_range], μ_cntr)
        I_filt_inv[:,cntr_cnt_next] = Butterworth_LPF(1000, I_inv, I_filt_inv[:,y_range], μ_cntr)

        p_q_filt[:, cntr_cnt_next] = p_q_theory(V_filt_poc[:,cntr_cnt_next], I_filt_poc[:,cntr_cnt_next])

        # Phase Locked Loop
        #_______________________________________________________________________
        fpll[:, cntr_cnt_next], θpll[:, cntr_cnt_next], pll_err_t, pll_err =
        Phase_Locked_Loop(V_filt_poc[:, cntr_cnt_next], t[i], μ_cntr,
        fpll[:, cntr_range], θpll[:, cntr_range], pll_err_t, pll_err)

        # DQ0 Transform
        #_______________________________________________________________________
        θ = (2*π*fs[i]*t[i])%(2*π) # grid forming

        V_dq0[:,cntr_cnt_next] = DQ0_transform(V_filt_poc[:,cntr_cnt_next], θ)
        I_dq0[:,cntr_cnt_next] = DQ0_transform(I_filt_inv[:,cntr_cnt_next], θ)

        # Current and Voltage PI Controller
        #_______________________________________________________________________
        Rf_V = [Rf_Va[i]; Rf_Vb[i]; Rf_Vc[i]]
        V_ref_dq0[:,cntr_cnt_next] = DQ0_transform(Rf_V, θ)

        I_ref_dq0[:,cntr_cnt_next], V_err, V_err_t =
        Voltage_Controller(V_dq0[:, cntr_cnt_next],
        V_ref_dq0[:, cntr_cnt_next],
        V_err, V_err_t, μ_cntr)

        #= --- to be removed
        I_ref_dq0[:,cntr_cnt_next] = DQ0_transform(Rf_I, (2*π*fs[i]*t[i])%(2*π))
        # --- =#

        Vd_dq0_new[:,cntr_cnt_next], I_err, I_err_t =
        Current_Controller(I_dq0[:, cntr_cnt_next],
        I_ref_dq0[:, cntr_cnt_next],
        I_err, I_err_t, μ_cntr)

        Vd_abc_new = Inv_DQ0_transform(Vd_dq0_new[:,cntr_cnt_next], θ)

        # Inverter Voltages
        Vd_a[i] = Vd_abc_new[1]
        Vd_b[i] = Vd_abc_new[2]
        Vd_c[i] = Vd_abc_new[3]

        #_______________________________________________________________________
        cntr_cnt = cntr_cnt_next

    elseif i != 1
        Vd_a[i] = Vd_a[i-1]
        Vd_b[i] = Vd_b[i-1]
        Vd_c[i] = Vd_c[i-1]
    else
        Vd_a[i] = 0
        Vd_b[i] = 0
        Vd_c[i] = 0
    end

    # Network Voltages
    Vs_a[i] = 0*Vo_p*sin(2*pi*fs[i]*t[i])
    Vs_b[i] = 0*Vo_p*sin(2*pi*fs[i]*t[i] - 120*π/180)
    Vs_c[i] = 0*Vo_p*sin(2*pi*fs[i]*t[i] + 120*π/180)

end
println("Progress : 100.0 %")

#%% Plots

T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
tsys = 1/fsys

#_______________________________________________________________________________
# Reference DQ0 Current Values
T_plot_start = 0
T_plot_end = 20

N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
N_range = N_plot_start:N_plot_end

N_plot_start_cntr = convert(Int64, round((T_plot_start/fsys + μ_cntr)/μ_cntr))
N_plot_end_cntr = convert(Int64, round((T_plot_end/fsys - μ_cntr)/μ_cntr))
N_range_cntr = N_plot_start_cntr:N_plot_end_cntr

p_I_dq0_d = plot(t_cntr[N_range_cntr], I_ref_dq0[1, N_range_cntr],
    label = "ref d-axis",
    legend = :bottomright,
    xlabel = "Time [s]",
    ylabel = "Current [A]",
    title = "DQ0 Current Transform\nDirect-Axis",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_I_dq0_d = plot!(t_cntr[N_range_cntr], I_dq0[1, N_range_cntr],
    label = "d-axis")

p_I_dq0_q = plot(t_cntr[N_range_cntr], I_ref_dq0[2, N_range_cntr],
        label = "ref q-axis",
        legend = :bottomright,
        xlabel = "Time [s]",
        ylabel = "Current [A]",
        title = "Quadrature-Axis",
        grid = true,
        foreground_color_grid = :black,
        minorgrid = true,
        thickness_scaling = 1.5,
        legendfont = font(5))
p_I_dq0_q = plot!(t_cntr[N_range_cntr], I_dq0[2, N_range_cntr],
        label = "q-axis")

p_I_dq0 = plot(p_I_dq0_d, p_I_dq0_q,# p_I_dq0_0,
    layout = (2, 1),
    legend = true,
    size = (850,900))

#display(p_I_dq0)

#_______________________________________________________________________________
# Reference DQ0 Voltage Values
T_plot_start = 0
T_plot_end = 20

N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
N_range = N_plot_start:N_plot_end

N_plot_start_cntr = convert(Int64, round((T_plot_start/fsys + μ_cntr)/μ_cntr))
N_plot_end_cntr = convert(Int64, round((T_plot_end/fsys - μ_cntr)/μ_cntr))
N_range_cntr = N_plot_start_cntr:N_plot_end_cntr

p_V_dq0_d = plot(t_cntr[N_range_cntr], V_ref_dq0[1, N_range_cntr],
    label = "ref d-axis",
    legend = :bottomright,
    xlabel = "Time [s]",
    ylabel = "Voltage [V]",
    title = "DQ0 Voltage Transform\nDirect-Axis",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_V_dq0_d = plot!(t_cntr[N_range_cntr], V_dq0[1, N_range_cntr],
    label = "d-axis")

p_V_dq0_q = plot(t_cntr[N_range_cntr], V_ref_dq0[2, N_range_cntr],
        label = "ref q-axis",
        legend = :bottomright,
        xlabel = "Time [s]",
        ylabel = "Voltage [V]",
        title = "Quadrature-Axis",
        grid = true,
        foreground_color_grid = :black,
        minorgrid = true,
        thickness_scaling = 1.5,
        legendfont = font(5))
p_V_dq0_q = plot!(t_cntr[N_range_cntr], V_dq0[2, N_range_cntr],
        label = "q-axis")

p_V_dq0 = plot(p_V_dq0_d, p_V_dq0_q,# p_I_dq0_0,
    layout = (2, 1),
    legend = true,
    size = (850,900))

#display(p_V_dq0)

# Output Instantaneous Electrical Potential
#_______________________________________________________________________________
# Phase a Voltage Signals
T_plot_start = 0
T_plot_end = 10

N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
N_range = N_plot_start:N_plot_end

# Phase a Control Signals
p_cntr_a = plot(t[N_range], Rf_Va[N_range],
    label = "Reference_a",
    xlabel = "Time [s]",
    ylabel = "Voltage [V]",
    title = "DC-AC Converter Control Phase A")
p_cntr_a = plot!(t[N_range], x[2,N_range],
    label = "Inverter Phase a")

# Phase b Control Signals
p_cntr_b = plot(t[N_range], Rf_Vb[N_range],
    label = "Reference_b",
    xlabel = "Time [s]",
    ylabel = "Voltage [V]",
    title = "DC-AC Converter Control Phase B")
p_cntr_b = plot!(t[N_range], x[5,N_range],
    label = "Inverter Phase b")

# Phase c Control Signals
p_cntr_c = plot(t[N_range], Rf_Vc[N_range],
    label = "Reference_c",
    xlabel = "Time [s]",
    ylabel = "Voltage [V]",
    title = "DC-AC Converter Control Phase C")
p_cntr_c = plot!(t[N_range], x[8,N_range],
    label = "Inverter Phase c")

p_v_cntr = plot(p_cntr_a, p_cntr_b, p_cntr_c,
    layout = (3, 1),
    legend = true,
    size = (850,900))

#display(p_v_cntr)

# Output Instantaneous Current
#_______________________________________________________________________________
# Phase a Current Signals
T_plot_start = 0 #T_s - 10
T_plot_end = 10 #T_s
N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))

N_range = N_plot_start:N_plot_end

# Phase a Control Signals
p_cntr_a = plot(t[N_range], Rf_Ia[N_range],
    label = "Reference_a",
    xlabel = "Time [s]",
    ylabel = "Current [A]",
    title = "DC-AC Converter Control Phase A")
p_cntr_a = plot!(t[N_range], x[3,N_range],
    label = "Inverter Phase a")

# Phase b Control Signals
p_cntr_b = plot(t[N_range], Rf_Ib[N_range],
    label = "Reference_b",
    xlabel = "Time [s]",
    ylabel = "Current [A]",
    title = "DC-AC Converter Control Phase B")
p_cntr_b = plot!(t[N_range], x[6,N_range],
    label = "Inverter Phase b")

# Phase c Control Signals
p_cntr_c = plot(t[N_range], Rf_Ic[N_range],
    label = "Reference_c",
    xlabel = "Time [s]",
    ylabel = "Current [A]",
    title = "DC-AC Converter Control Phase C")
p_cntr_c = plot!(t[N_range], x[9,N_range],
    label = "Inverter Phase c")

p_i_cntr = plot(p_cntr_a, p_cntr_b, p_cntr_c,
    layout = (3, 1),
    legend = true,
    size = (850,900))

#display(p_i_cntr)

# Phase-Locked-Loop
#_______________________________________________________________________________
T_plot_start = 0
T_plot_end = 20

N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
N_range = N_plot_start:N_plot_end

N_plot_start_cntr = convert(Int64, round((T_plot_start/fsys + μ_cntr)/μ_cntr))
N_plot_end_cntr = convert(Int64, round((T_plot_end/fsys - μ_cntr)/μ_cntr))
N_range_cntr = N_plot_start_cntr:N_plot_end_cntr

p_pll_f = plot(t_cntr[N_range_cntr], fpll[1, N_range_cntr], label = "PLL Frequency",
    legend = :bottomright,
    xlabel = "Time [s]",
    ylabel = "Frequency [Hz]",
    title = "Phase-Locked-Loop",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_pll_f = plot!(t[N_range], fs[N_range], label = "System Frequency")

θs = (2*π*fs.*t .+ Phase).%(2*π)

N_range = i_sp*N_plot_start:i_sp:N_plot_end
θe = sin.(θs[N_range] .- θpll[1, N_range_cntr])
θe = (180/π).*asin.(θe)
p_pll_θ = plot(t[N_range], θe,
    label = "PLL Error",
    legend = :bottomright,
    xlabel = "Time [s]",
    ylabel = "Phase [°]",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
#=p_pll_θ = plot!(t_cntr[N_range_cntr], (180/π).θPLL[N_range_cntr],
    label = "PLL Phase Angle")
p_pll_θ = plot!(t[N_range], (180/π).θs[N_range],
    label = "Source Phase Angle")=#

p_pll = plot(p_pll_f, p_pll_θ,
    layout = (2, 1),
    legend = true,
    size = (900,900))

display(p_pll)

# Instantaneous Powers
#_______________________________________________________________________________
T_plot_start = 0
T_plot_end = 10

N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
N_range = N_plot_start:N_plot_end

Uc = 0.001 # unit conversion
p_p_a = plot(t[N_range], Uc*p_inv[1, N_range], label = "Inverter Phase a",
    legend = :bottomright,
    xlabel = "Time [s]",
    ylabel = "Power [kW]",
    title = "Instantaneous Active Powers",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_p_a = plot!(t[N_range], Uc*p_net[1, N_range], label = "Network Phase a")
p_p_b = plot(t[N_range], Uc*p_inv[2, N_range], label = "Inverter Phase b",
    legend = :bottomright,
    xlabel = "Time [s]",
    ylabel = "Power [kW]",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_p_b = plot!(t[N_range], Uc*p_net[2, N_range],label = "Network Phase b")
p_p_c = plot(t[N_range], Uc*p_inv[3, N_range], label = "Inverter Phase c",
    legend = :bottomright,
    xlabel = "Time [s]",
    ylabel = "Power [kW]",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_p_c = plot!(t[N_range], Uc*p_net[3, N_range], label = "Network Phase c")

p_inv_t = p_inv[1, :] .+ p_inv[2, :] .+ p_inv[3, :]
p_net_t = p_net[1, :] .+ p_net[2, :] .+ p_net[3, :]
p_p_t = plot(t[N_range], Uc*p_inv_t[N_range], label = "Inverter Total",
    legend = :bottomright,
    xlabel = "Time [s]",
    ylabel = "Power [kW]",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_p_t = plot!(t[N_range], Uc*p_net_t[N_range], label = "Network Total")

p_p = plot(p_p_a, p_p_b, p_p_c, p_p_t,
    layout = (4, 1),
    legend = true,
    size = (900,1200))

#display(p_p)

# RMS Voltage Values
#_______________________________________________________________________________
T_plot_start = 0
T_plot_end = 10

N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))

T_sp_rms = 10*fsys #samples in a second, i.e. sampling of microcontroller;
i_sp_rms = convert(Int64, round((1/(μps*T_sp_rms))))
N_range = N_plot_start:i_sp_rms:N_plot_end

N_plot_start_cntr = convert(Int64, round((T_plot_start/fsys + μ_cntr)/μ_cntr))
N_plot_end_cntr = convert(Int64, round((T_plot_end/fsys - μ_cntr)/μ_cntr))
N_range_cntr = N_plot_start_cntr:N_plot_end_cntr

p_v_rms = plot(t[N_range], Vinv_ph[1,2,N_range],
    legend = :bottomright,
    label = "Inverter Phase a",
    xlabel = "Time [s]",
    ylabel = "Electrical Potential [V]",
    title = "RMS Inverter Voltages",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_v_rms = plot!(t[N_range], Vinv_ph[2,2,N_range],
    label = "Inverter Phase b")
p_v_rms = plot!(t[N_range], Vinv_ph[3,2,N_range],
    label = "Inverter Phase c")

p_v_ang = plot(t[N_range], (180/π)*Vinv_ph[1,3,N_range],
    legend = :bottomright,
    label = "Inverter Phase a",
    xlabel = "Time [s]",
    ylabel = "Degrees [°]",
    title = "Inverter Phase Angles",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_v_ang = plot!(t[N_range], (180/π)*Vinv_ph[2,3,N_range],
    label = "Inverter Phase b")
p_v_ang = plot!(t[N_range], (180/π)*Vinv_ph[3,3,N_range],
    label = "Inverter Phase c")

PLL_ph_a = θpll[1, N_range_cntr] .- 2*π*fpll[1, N_range_cntr].*t_cntr[N_range_cntr]
PLL_ph_b = θpll[2, N_range_cntr] .- 2*π*fpll[2, N_range_cntr].*t_cntr[N_range_cntr]
PLL_ph_c = θpll[3, N_range_cntr] .- 2*π*fpll[3, N_range_cntr].*t_cntr[N_range_cntr]

for i in 1:length(PLL_ph_a[N_range_cntr])
    PLL_ph_a[i] = (PLL_ph_a[i] + 2*π*floor(t_cntr[i]/0.02))*180/pi
    PLL_ph_b[i] = (PLL_ph_b[i] + 2*π*floor(t_cntr[i]/0.02))*180/pi
    PLL_ph_c[i] = (PLL_ph_c[i] + 2*π*floor(t_cntr[i]/0.02))*180/pi
    if PLL_ph_a[i] > 180
        PLL_ph_a[i] = PLL_ph_a[i] - 360
    end
    if PLL_ph_b[i] > 180
        PLL_ph_b[i] = PLL_ph_b[i] - 360
    end
    if PLL_ph_c[i] > 180
        PLL_ph_c[i] = PLL_ph_c[i] - 360
    end
    if PLL_ph_a[i] < -180
        PLL_ph_a[i] = PLL_ph_a[i] + 360
    end
    if PLL_ph_b[i] < -180
        PLL_ph_b[i] = PLL_ph_b[i] + 360
    end
    if PLL_ph_c[i] < -180
        PLL_ph_c[i] = PLL_ph_c[i] + 360
    end
end

p_v_ang = plot!(t_cntr[N_range_cntr], PLL_ph_a[N_range_cntr],
    label = "PLL Phase a",
    linestyle = :dash,
    linewidth = 2)
p_v_ang = plot!(t_cntr[N_range_cntr], PLL_ph_b[N_range_cntr],
    label = "PLL Phase b",
    linestyle = :dash,
    linewidth = 2)
p_v_ang = plot!(t_cntr[N_range_cntr], PLL_ph_c[N_range_cntr],
    label = "PLL Phase c",
    linestyle = :dash,
    linewidth = 2)

p_v_rms_ang = plot(p_v_rms, p_v_ang,
    layout = (2, 1),
    legend = :bottomright,
    size = (900,900))

display(p_v_rms_ang)

# RMS Current Values
#_______________________________________________________________________________
T_plot_start = 0
T_plot_end = 10
N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
N_range = N_plot_start:i_sp_rms:N_plot_end

p_i_rms = plot(t[N_range], Iinv_ph[1,2,N_range],
    legend = :bottomright,
    label = "Inverter Phase a",
    xlabel = "Time [s]",
    ylabel = "Current [A]",
    title = "RMS Inverter Current",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_i_rms = plot!(t[N_range], Iinv_ph[2,2,N_range],
    label = "Inverter Phase b")
p_i_rms = plot!(t[N_range], Iinv_ph[3,2,N_range],
    label = "Inverter Phase c")

p_i_ang = plot(t[N_range], (180/π)*Iinv_ph[1,3,N_range],
    legend = :bottomright,
    label = "Inverter Phase a",
    xlabel = "Time [s]",
    ylabel = "Degrees [°]",
    title = "Inverter Phase Angles",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_i_ang = plot!(t[N_range], (180/π)*Iinv_ph[2,3,N_range],
    label = "Inverter Phase b")
p_i_ang = plot!(t[N_range], (180/π)*Iinv_ph[3,3,N_range],
    label = "Inverter Phase c")

p_i_rms_ang = plot(p_i_rms, p_i_ang,
    layout = (2, 1),
    legend = true,
    size = (900,900))

display(p_i_rms_ang)

# Instantaneous Real and Imaginary Powers, and Conventional Active and Reactive
#_______________________________________________________________________________
T_plot_start = 0
T_plot_end = 10

N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
N_range = N_plot_start:N_plot_end

N_plot_start_cntr = convert(Int64, round((T_plot_start/fsys + μ_cntr)/μ_cntr))
N_plot_end_cntr = convert(Int64, round((T_plot_end/fsys - μ_cntr)/μ_cntr))
N_range_cntr = N_plot_start_cntr:N_plot_end_cntr

Uc = 0.001 # unit conversion

P_active_a = (Vinv_ph[1,2,:].*Iinv_ph[1,2,:]).*cos.(Vinv_ph[1,3,:] .- Iinv_ph[1,3,:])
P_active_b = (Vinv_ph[2,2,:].*Iinv_ph[2,2,:]).*cos.(Vinv_ph[2,3,:] .- Iinv_ph[2,3,:])
P_active_c = (Vinv_ph[3,2,:].*Iinv_ph[3,2,:]).*cos.(Vinv_ph[3,3,:] .- Iinv_ph[3,3,:])
P_active = P_active_a + P_active_b + P_active_c

p_p_r_a = plot(t[N_range], Uc*p_q_inst[1, N_range],
    legend = :bottomright,
    label = "Real Power",
    xlabel = "Time [s]",
    ylabel = "Power [kW]",
    title = "Inverter Real and Active Power",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_p_r_a = plot!(t_cntr[N_range_cntr], Uc*p_q_filt[1, N_range_cntr], label = "Filtered Real Power")
p_p_r_a = plot!(t[N_range], Uc*P_active[N_range], label = "Active Power", lw = 2)

Q_reactive_a = (Vinv_ph[1,2,:].*Iinv_ph[1,2,:]).*sin.(Vinv_ph[1,3,:] .- Iinv_ph[1,3,:])
Q_reactive_b = (Vinv_ph[2,2,:].*Iinv_ph[2,2,:]).*sin.(Vinv_ph[2,3,:] .- Iinv_ph[2,3,:])
Q_reactive_c = (Vinv_ph[3,2,:].*Iinv_ph[3,2,:]).*sin.(Vinv_ph[3,3,:] .- Iinv_ph[3,3,:])
Q_reactive = Q_reactive_a + Q_reactive_b + Q_reactive_c
p_p_i_q = plot(t[N_range], Uc*p_q_inst[2 ,N_range], label = "Imaginary Power",
    legend = :bottomright,
    xlabel = "Time [s]",
    ylabel = "Power [kVAi / kVAr]",
    title = "Inverter Imaginary and Reactive Power",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true,
    thickness_scaling = 1.5,
    legendfont = font(5))
p_p_i_q = plot!(t_cntr[N_range_cntr], Uc*p_q_filt[2, N_range_cntr], label = "Filtered Imaginary Power")
p_p_i_q = plot!(t[N_range], Uc*Q_reactive[N_range], label = "Reactive Power", lw = 2)

p_p_real_imag_act_react = plot(p_p_r_a, p_p_i_q,
    layout = (2, 1),
    legend = true,
    size = (900,900))
display(p_p_real_imag_act_react)

# Calculation Checks
Z1 = conj((3/2)*Vinv_ph[1,2,end]^2/(P_active[N] + 1im*Q_reactive[N]))
S1 = conj((3/2)*Vinv_ph[1,2,end]^2/(RL + 1im*XL))
I1 = Vinv_ph[1,2,end]/(abs(RL + 1im*XL)*sqrt(2))

# Fast Fourier Transforms
#_______________________________________________________________________________
T_plot_start = 0
T_plot_end = 1

x_lim = (-1500, +1500)
x_ticks = -1500:250:1500
N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
N_range = N_plot_start:N_plot_end
N_length = length(N_range)
freqs = fftshift(fftfreq(N_length, 1/μps))

Uc = 0.001 # unit conversion

N_plot_start_cntr = convert(Int64, round((T_plot_start/fsys + μ_cntr)/μ_cntr))
N_plot_end_cntr = convert(Int64, round((T_plot_end/fsys - μ_cntr)/μ_cntr))
N_range_cntr = N_plot_start_cntr:N_plot_end_cntr
N_length_cntr = length(N_range_cntr)
freqs_cntr = fftshift(fftfreq(N_length_cntr, 1/μ_cntr))
#=
fft_p_real = (2/N_length_cntr)*fftshift(fft(Uc*p_q_inst[1, N_range]))
FFTs_p_real = plot(freqs_cntr, abs.(fft_p_real),
    seriestype = :scatter, line = :stem,
    title = "Inverter POC Real Power",
    label = "Shifted FFT",
    xlim = (-500, +500),
    xticks = -500:100:500,
    xlabel = "Frequency [Hz]",
    ylabel = "Real Power [W]",
    markershape = :circle, markersize = 2,
    legend = false);
=#

fft_p_real = (2/N_length)*fftshift(fft(Uc*p_q_inst[1, N_range]))
fft_p_real_filt = (2/N_length_cntr)*fftshift(fft(Uc*p_q_filt[1, N_range_cntr]))
FFTs_p_real = plot(freqs, abs.(fft_p_real), label = "Shifted FFT",
    seriestype = :scatter, line = :stem,
    title = "Inverter POC Real Power",
    xlim = x_lim,
    xticks = x_ticks,
    xlabel = "Frequency [Hz]",
    ylabel = "Real Power [W]",
    markershape = :circle, markersize = 2,
    legend = true);
FFTs_p_real = plot!(freqs_cntr, abs.(fft_p_real_filt), label = "Filtered FFT")

fft_p_imag = (2/N_length)*fftshift(fft(Uc*p_q_inst[2, N_range]))
fft_p_imag_filt = (2/N_length_cntr)*fftshift(fft(Uc*p_q_filt[2, N_range_cntr]))
FFTs_p_imag = plot(freqs, abs.(fft_p_imag), label = "Shifted FFT",
    seriestype = :scatter, line = :stem,
    title = "Inverter POC Imaginary Power",
    xlim = x_lim,
    xticks = x_ticks,
    xlabel = "Frequency [Hz]",
    ylabel = "Imaginary Power [VAi]",
    markershape = :circle, markersize = 2,
    legend = false);
FFTs_p_imag = plot!(freqs_cntr, abs.(fft_p_imag_filt), label = "Filtered FFT")

fft_v_poc_a = (2/N_length)*fftshift(fft(x[2,N_range]))
FFTs_v_a = plot(freqs, abs.(fft_v_poc_a),
    seriestype = :scatter, line = :stem,
    title = "Inverter POC Voltage",
    label = "Phase a",
    xlim = x_lim,
    xticks = x_ticks,
    xlabel = "Frequency [Hz]",
    ylabel = "Electrical Potential [V]",
    markershape = :circle, markersize = 2,
    legend = false);

fft_i_poc_a = (2/N_length)*fftshift(fft(x[3,N_range]))
FFTs_i_a = plot(freqs, abs.(fft_i_poc_a),
    seriestype = :scatter, line = :stem,
    title = "Inverter POC Current",
    label = "Phase a",
    xlim = x_lim,
    xticks = x_ticks,
    xlabel = "Frequency [Hz]",
    ylabel = "Current [A]",
    markershape = :circle, markersize = 2,
    legend = false);

p_fft = plot(FFTs_p_real, FFTs_p_imag, FFTs_v_a, FFTs_i_a,
    layout = (2, 2),
    size = (1200,900),
    legend = true,
    margin = 5Plots.mm)
#display(p_fft)

# Output Instantaneous Electrical Potential
#_______________________________________________________________________________
# Phase a Voltage Signals
T_plot_start = 0
T_plot_end = 1

N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
N_range = i_sp*3:i_sp*10

# Phase a Control Signals
p_vsc_a = plot(t[N_range], Vd_a[N_range],
    label = "Reference_a",
    xlabel = "Time [s]",
    ylabel = "Voltage [V]",
    title = "DC-AC Converter Terminal Voltage\nPhase A",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true)

# Phase b Control Signals
p_vsc_b = plot(t[N_range], Vd_b[N_range],
    label = "Reference_b",
    xlabel = "Time [s]",
    ylabel = "Voltage [V]",
    title = " Phase B",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true)

# Phase c Control Signals
p_vsc_c = plot(t[N_range], Vd_c[N_range],
    label = "Reference_c",
    xlabel = "Time [s]",
    ylabel = "Voltage [V]",
    title = "Phase C",
    grid = true,
    foreground_color_grid = :black,
    minorgrid = true)

p_vsc = plot(p_vsc_a, p_vsc_b, p_vsc_c,
    layout = (3, 1),
    legend = false,
    size = (850,900))
#display(p_vsc)

# Save plots
#_______________________________________________________________________________

#=
savefig(p_I_dq0, "p_I_dq0.png")
savefig(p_V_dq0, "p_V_dq0.png")
savefig(p_v_cntr, "p_v_cntr.png")
savefig(p_i_cntr, "p_i_cntr.png")
savefig(p_pll, "p_pll.png")
savefig(p_p, "p_p.png")
savefig(p_v_rms_ang, "p_v_rms_ang.png")
savefig(p_i_rms_ang, "p_i_rms_ang.png")
savefig(p_p_real_imag_act_react, "p_p_real_imag_act_react.png")
savefig(p_fft, "p_fft.png")
savefig(p_vsc, "p_vsc.png")
=#
print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
