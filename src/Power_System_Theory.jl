function RMS(θ, t_signals)

    # Calcutates the DC offset, RMS magnitude, and phase angle relative to the
    # frequency (θ = 2*π*t[:]) for a three phase system
    # note that the reference is a cosine - subtract 90 degrees for sine.

    i_length = size(t_signals, 1)

    rms = Array{Float64, 2}(undef, 3, 3)
    rms = fill!(rms, 0)

    rm = sqrt(2)

    # Least squares method to find RMS of fundamental
    LS = [ones(i_length, 1) rm*cos.(θ) rm*sin.(θ)]

    if !any(isnan, LS)

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

    #v_dq0 = Park_Transform(Clarke_Transform(v_abc), θ)
    v_dq0 = sqrt(2/3)*[cos(θ) cos(θ - 2π/3) cos(θ + 2π/3);
                    -sin(θ) -sin(θ - 2π/3) -sin(θ + 2π/3);
                    sqrt(2)/2 sqrt(2)/2 sqrt(2)/2]*v_abc

    return v_dq0
end

function Inv_DQ0_transform(v_dq0, θ)

    #v_abc  = Inv_Clarke_Transform(Inv_Park_Transform(v_dq0, θ))
    v_abc = sqrt(2/3)*[cos(θ) -sin(θ) sqrt(2)/2;
                        cos(θ - 2π/3) -sin(θ - 2π/3) sqrt(2)/2;
                        cos(θ + 2π/3) -sin(θ + 2π/3) sqrt(2)/2]*v_dq0

    return v_abc
end

function p_q_theory(V_abc, I_abc)

    #= Theory:
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
        powers are constant. The real power is equal to the conventional definition of
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

    #pq0 = [I_αβγ[1] I_αβγ[2] 0; -I_αβγ[2] I_αβγ[1] 0; 0 0 I_αβγ[3]]*V_αβγ # also works

    return pq0 = [V_αβγ[1] V_αβγ[2] 0; V_αβγ[2] -V_αβγ[1] 0; 0 0 V_αβγ[3]]*I_αβγ
end

function Inv_p_q_v(V_αβγ, pq0)

    if V_αβγ[3] != 0

        a_det = (V_αβγ[1]^2 + V_αβγ[2]^2)
        a_adj = [V_αβγ[1] V_αβγ[2] 0; V_αβγ[2] -V_αβγ[1] 0; 0 0 (V_αβγ[1]^2 + V_αβγ[2]^2)/V_αβγ[3]]

        if a_det != 0
            i_αβγ = (1/a_det)*a_adj*pq0
        else
            i_αβγ = [0; 0; 0]
        end

    else

        a_det = (V_αβγ[1]^2 + V_αβγ[2]^2)
        a_adj = [V_αβγ[1] V_αβγ[2]; V_αβγ[2] -V_αβγ[1]]

        if a_det != 0
            i_αβ = (1/a_det)*a_adj*pq0[1:2]
            i_αβγ = [i_αβ ; 0]
        else
            i_αβγ = [0; 0; 0]
        end

    end

    return i_αβγ
end

function Inv_p_q_i(I_αβγ, pq0)

    if I_αβγ[3] != 0

        a_det = (I_αβγ[1]^2 + I_αβγ[2]^2)
        a_adj = [I_αβγ[1] -I_αβγ[2] 0; I_αβγ[2] I_αβγ[1] 0; 0 0 (I_αβγ[1]^2 + I_αβγ[2]^2)/I_αβγ[3]]

        if a_det != 0
            v_αβγ = (1/a_det)*a_adj*pq0
        else
            v_αβγ = [0; 0; 0]
        end
        v_αβγ[3] = 0

    else

        a_det = (I_αβγ[1]^2 + I_αβγ[2]^2)
        a_adj = [I_αβγ[1] -I_αβγ[2]; I_αβγ[2] I_αβγ[1]]

        if a_det != 0
            v_αβ = (1/a_det)*a_adj*pq0[1:2]
            v_αβγ = [v_αβ ; 0]
        else
            v_αβγ = [0; 0; 0]
        end

    end

    return v_αβγ
end

function Series_Load_Impedance(S, pf, vrms; fsys = 50)

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
        X = -imag(Y)^-1

        L_C = -1/(X*ω)

    end

    return R, L_C, X, Z
end

function Parallel_Load_Impedance(S, pf, vrms; fsys = 50)

    ω = 2*π*fsys

    if pf > 0

        P = pf*S #W, Active Power
        Q = sqrt(S^2 - P^2) #VAr, Reactive Power

        Z = conj(3*vrms^2/(P + 1im*Q)) #Load Impedance

        Y = 1/Z # Impedance in parallel

        R = real(Y)^-1
        X = -imag(Y)^-1

        L_C = X/ω

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