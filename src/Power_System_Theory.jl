"""
# Description
Per phase least squares calculation of fundamental rms signal
"""
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

"""
# Description
Three phase DQ0 calculation of rms signals
"""
function DQ_RMS(i_abc)

    return i_rms = sqrt(1/3)*norm(DQ0_transform(i_abc, 0)[1:3])
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

"""
R, L_C, X, Z = Parallel_Load_Impedance(S, pf, vrms; fsys = 50, type_sign = nothing)

# Arguments
- `S::Float`: 3 phase Apparent Power [VA]
- `pf::Float`: power factor. pf > 0 -> inductive load, pf < 0 -> capacitive load
- `Vrms::Float`: Line to Neutral rms voltage of external network [V]

# Keyword Arguments
- `fsys::Float`: system frequency [Hz]
- `type_sign::String`: "L" or "C" -> for purely capacitive or inductive loads

# Return Values
- `R::Float`: resistance [Ω]
- `L_C::Float`: inductance [H] or capacitance [F]
- `X::Float`: reactance [Ω]
- `Z::Complex`: impedance [Ω]

# Theory
Converts a Resistance-Inductance Load, or a Resistance-Capacitance Load
from power to circuit values, where the components are placed in parallel.
"""
function Parallel_Load_Impedance(S, pf, vrms; fsys = 50, type_sign = nothing)

    ω = 2*π*fsys

    if sign(pf) == 1 || type_sign == "L"

        P = pf * S #W, Active Power
        Q = sqrt(S^2 - P^2) #VAr, Reactive Power

        Z = conj(3 * vrms^2 / (P + 1im * Q)) #Load Impedance

        Y = 1/Z # Impedance in parallel

        R = real(Y)^-1
        X = -imag(Y)^-1

        L_C = X / ω

    elseif sign(pf) == -1 || type_sign == "C"
        P = -pf * S #W, Active Power - average instance power
        Q = -sqrt(S^2 - P^2) #VAr, Reactive Power

        Z = conj(3 * vrms^2 / (P + 1im * Q))

        Y = 1/Z # Impedance in parallel

        R = real(Y)^-1
        X = imag(Y)^-1

        L_C = 1 / (X * ω)

    end

    return R, L_C, X, Z
end

"""
    Filter_Design(Sr, fs, fltr; Vrms = 230, Vdc = 800, ΔILf_ILf = 0.15, ΔVCf_VCf = 0.01537)

# Description
Designs the filter inductances and capacitances for the sources

# Arguments
- `Sr::Float`: rated 3 phase Apparent Power [VA]
- `fs::Float`:switching frequency
- `fltr::String`: topology "LCL", "LC", "L"

# Keyword Arguments
- `Vrms::Float`: nominal AC rms voltage
- `Vdc::Float`: nominal DC voltage
- `ΔILf_ILf::Float`: current ripple ratio
- `ΔVCf_VCf::Float`: voltage ripple ratio
"""
function Filter_Design(Sr, fs, fltr; Vrms = 230, Vdc = 800, ΔILf_ILf = 0.15, ΔVCf_VCf = 0.01537)

    #= Theory:
        The filtering capacitors C should be chosen such that the resonant frequency
        1/sqrt(Ls*C) is approximately sqrt(ωn * ωs), where the ωn is the angular
        frequency of the gird voltage, and ωs is the angular switching frequency
        used to turn on/off the switches)

        The attenuation of the LCL-filter is 60db/decade for frequencies above the
        resonant frequency, therefore lower switching frequency for the converter
        can be used. It also provides better decoupling between the filter and the
        grid impedance and lower current ripple accross the grid inductor.

        The LCL filter has good current ripple attenuattion even with small inductance
        values. However it can bring also resonances and unstable states into the system.
        Therefore the filter must be designed precisely accoring to the parameters of the
        specific converter.

        The importatnt parameter of the filter is its cut-off frequency. The cut-off
        frequency of the filter must be minimally one half of the switching frequency of
        the converter, because the filter must have enough attenuation in the range of the
        converter's switching frequency.

        The LCL filter will be vulnerable to oscillation too and it will magnify frequencies
        around its cut-off frequency. Therefore the filter is added with damping. The simplest
        way is to add damping resistor. In general there are four possible places where the
        resistor can be placed series/parallel to the inverter side inductor or series/parallel
        to the filter capacitor.

        The variant with the resistor connected in series with the filter capacitor has been
        chosen. The peak near the resonant frequency can be significantly attenuated. This is
        a simple and reliable solution, but it increases the heat losses in the system and it
        greatly decreases the efficiency of the filter. This problem can be solved by active
        damping. Such a resistor reduces the voltage across the capacitor by a voltage proportional
        to the current that flows through it. This can be also done in the control loop. The
        current through the filter capacitor is measured and differentiatbed by the term
        (s*Cf*R_C). A real resistor is not used and the calculated value is subtracted from the
        demanded current. In this way the filter is actively damped with a virtual resistor
        without any losses. The disadvantage of this method is that an additional current sensor
        is required and the differentiator may bring noise problems because it amplifies high
        frequency signals.
    =#

    #= Practical Example:
            L_filter = 70e-6
            R_filter = 1.1e-3
            R_filter_C = 7e-3
            C_filter = 250e-6
            =#

    #____________________________________________________________
    # ratios

    r_p = 200

    f_p = 5

    Ir = ΔILf_ILf
    Vr = ΔVCf_VCf

    i_lim_r = 1.5
    v_lim_r = 1.5

    #____________________________________________________________
    # Inductor Design
    Vorms = Vrms*1.05
    Vop = Vorms*sqrt(2)

    Zl = 3*Vorms*Vorms/Sr

    Iorms = Vorms/Zl
    Iop = Iorms*sqrt(2)

    ΔIlfmax = Ir*Iop

    Lf_1 = Vdc/(4*fs*ΔIlfmax)

    #____________________________________________________________
    # Capacitor Design
    Vorms = Vrms*0.95
    Vop = Vorms*sqrt(2)

    Zl = 3*Vorms^2/Sr

    Iorms = Vorms/Zl
    Iop = Iorms*sqrt(2)
    Ir = Vdc/(4*fs*Lf_1*Iop)
    ΔIlfmax = Ir*Iop
    ΔVcfmax = Vr*Vop

    Cf = ΔIlfmax/(8*fs*ΔVcfmax)

    fc = 1/(2π*sqrt(Lf_1*Cf))

    #____________________________________________________________
    # Resistor Design

    R_1 = 200 * Lf_1

    ωc = 2π*fc

    R_C = 1/(3*ωc*Cf)

    # v_limit
    v_limit = v_lim_r * Vdc * (1 + ΔVCf_VCf/2)

    if fltr == "LCL"

        i_limit = i_lim_r*Iop

        fc = fs/f_p
        ωc = 2π*fc
        Lf_2 = Lf_1/(ωc^2*Lf_1*Cf - 1)

        R_2 = r_p * Lf_2

        return Lf_1, Lf_2, Cf, fc, R_1, R_2, R_C, i_limit, v_limit

    elseif fltr == "LC"

        i_limit = i_lim_r*Iop*(1 + ΔILf_ILf/2)

        return Lf_1, Cf, fc, R_1, R_C, i_limit, v_limit

    elseif fltr == "L"

        i_limit = i_lim_r*Iop*(1 + ΔILf_ILf/2)

        return Lf_1, R_1, i_limit
    end
end

#############################
# Helper functions for cable_layout

function set_bounds(variable, start_value, low_bound, up_bound)
    if !is_fixed(variable)
        set_lower_bound(variable, low_bound)
        set_upper_bound(variable, up_bound)
        set_start_value(variable, start_value)
    end
end


function get_degree(CM = CM) # how many cables are connected to a node? maybe remove function if not used
    result = zeros(Int, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = count(x -> x != 0, CM[i,:])
    end

    result
end

function get_cable_connections(CM = CM) # which cables are connected to which nodes

    result = Vector{Vector{Int64}}(undef, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = filter(x -> x != 0, abs.(CM[i,:]))
    end

    return result
end

function get_node_connections(CM = CM) # which nodes are connected to each other, including the self-connections

    result = Vector{Vector{Int64}}(undef, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = findall(x -> x != 0, CM[i,:])
        push!(result[i], i)
    end

    return result
end

"""
    p_load_total, q_load_total, s_load_total, s_source_total = CheckPowerBalance(parameters)

# Description
Determines based on the parameters of the grid the total power (active and reactive) drawn from all load
and the total power provided by the sources.
Thereby, steady state is assumed.

# Arguments
- `parameters::Dict`: Completly filled parameter dict which defines the electrical power grid used in the env/nodeconstructor.
- `num_source::Int`: number of sources in the grid (todo: calulate based on parameter dict?)
- `num_load::Int`: number of num_load in the grid (todo: calulate based on parameter dict?)
- `CM::Matrix`: connectivity matrix describing the connections in the grid


# Return Values
- `p_load_total::float`: total active power drawn by all loads (passive components as well as controlled with negative reference value)
- `q_load_total::float`: total reactive power drawn by all loads
-  s_load_total::float`: total apparent power drawn by all loads
- `s_source_total::float`: total apparent power provided by all sources in steady state
"""
function CheckPowerBalance(parameters, num_source, num_load, CM)

    num_nodes = num_source + num_load
    num_cables = Int(maximum(CM))

    p_load_total = 0
    q_load_total = 0
    s_source_total = 0

    for i = 1:num_nodes

        if i <= num_source

            s_source_total = s_source_total + parameters["source"][i]["pwr"]/parameters["grid"]["phase"]

            if parameters["source"][i]["mode"] in ["PQ", 3]

                if parameters["source"][i]["p_set"] < 0
                    p_load_total = p_load_total + parameters["source"][i]["p_set"]/parameters["grid"]["phase"]
                end

                if parameters["source"][i]["q_set"] < 0
                    q_load_total = q_load_total + parameters["source"][i]["q_set"]/parameters["grid"]["phase"]
                end

            elseif parameters["source"][i]["mode"] in ["PV", 4]

                if parameters["source"][i]["p_set"] < 0
                    p_load_total = p_load_total + parameters["source"][i]["p_set"]/parameters["grid"]["phase"]
                end
            end
        else

            p_load_total = p_load_total + parameters["load"][i-num_source]["pwr"]/parameters["grid"]["phase"]
            q_load_total = q_load_total + parameters["load"][i-num_source]["pwr"]/parameters["grid"]["phase"]
        end
    end
    s_load_total = sqrt(p_load_total^2 + q_load_total^2)
    return p_load_total, q_load_total, s_load_total, s_source_total
end

optimizer_status = Dict()
function get_optimizer_status(model)
    
    status = Dict(
                    "termination_status" => termination_status(model),
                    "primal_status"      => (primal_status(model)),
                    "objective_value"    => (objective_value(model))
                )           

    return status
end

function layout_cabels(CM, num_source, num_load, parameters, verbosity = 0)

    model = Model(Ipopt.Optimizer)
    #set_optimizer_attributes(model, "tol" => 1e-1)

    zero_expression = @NLexpression(model, 0.0)

    # Constant values
    omega = 2π*parameters["grid"]["fs"]

    # for every Source: v is fixed 230
    # for one Source: theta is fixed 0

    num_nodes = num_source + num_load
    num_cables = Int(maximum(CM))

    @variable(model, nodes[1 : num_nodes, ["v", "theta", "P", "Q"]])

    # cal total load[pwr]
    total_P_load, total_Q_load, s_load_total, total_S_source = CheckPowerBalance(parameters, num_source, num_load, CM)

    #=
    total_P_load = 0
    total_Q_load = 0
    total_S_source = 0


    for i = 1:num_nodes

        if i <= num_source

            total_S_source = total_S_source + parameters["source"][i]["pwr"]/parameters["grid"]["phase"]

            if parameters["source"][i]["mode"] in ["PQ Control", 3]

                if parameters["source"][i]["p_set"] < 0
                    total_P_load = total_P_load + parameters["source"][i]["p_set"]/parameters["grid"]["phase"]
                end

                if parameters["source"][i]["q_set"] < 0
                    total_Q_load = total_Q_load + parameters["source"][i]["q_set"]/parameters["grid"]["phase"]
                end

            elseif parameters["source"][i]["mode"] in ["PV Control", 4]

                if parameters["source"][i]["p_set"] < 0
                    total_P_load = total_P_load + parameters["source"][i]["p_set"]/parameters["grid"]["phase"]
                end
            end
        else

            total_P_load = total_P_load + parameters["load"][i-num_source]["pwr"]/parameters["grid"]["phase"]
            total_Q_load = total_Q_load + parameters["load"][i-num_source]["pwr"]/parameters["grid"]["phase"]
        end
    end
    =#

    idx_p_mean_cal = []
    idx_q_mean_cal = []

    for i = 1:num_nodes
        if i <= num_source

            if parameters["source"][i]["control_type"] == "classic"

                if parameters["source"][i]["mode"] in ["Swing", "Voltage", 1, 7]

                    fix(nodes[i, "v"], parameters["source"][i]["v_pu_set"] * parameters["grid"]["v_rms"])
                    fix(nodes[i, "theta"], (π/180)*parameters["source"][i]["v_δ_set"])

                    set_bounds(nodes[i, "P"], (total_P_load) / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # come from parameter dict/user?
                    set_bounds(nodes[i, "Q"], (total_Q_load) / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # P and Q are the average from power, excluding cable losses
                    push!(idx_p_mean_cal, i)
                    push!(idx_q_mean_cal, i)

                elseif parameters["source"][i]["mode"] in ["PQ", 2]

                    fix(nodes[i, "P"], parameters["source"][i]["p_set"]/parameters["grid"]["phase"])
                    fix(nodes[i, "Q"], parameters["source"][i]["q_set"]/parameters["grid"]["phase"])

                    set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
                    set_bounds(nodes[i, "v"], parameters["grid"]["v_rms"], 0.95*parameters["grid"]["v_rms"], 1.05*parameters["grid"]["v_rms"])

                elseif parameters["source"][i]["mode"] in ["PV", 9]

                    fix(nodes[i, "P"], parameters["source"][i]["p_set"]/parameters["grid"]["phase"])
                    fix(nodes[i, "v"], parameters["source"][i]["v_pu_set"] * parameters["grid"]["v_rms"])
                    set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
                    set_bounds(nodes[i, "Q"], (total_Q_load) / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # P and Q are the average from power, excluding cable losses
                    push!(idx_q_mean_cal, i)

                elseif parameters["source"][i]["mode"] in ["Semi-Synchronverter", 6]

                    fix(nodes[i, "v"], parameters["source"][i]["v_pu_set"] * parameters["grid"]["v_rms"])
                    set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
                    set_bounds(nodes[i, "P"], total_P_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"])
                    set_bounds(nodes[i, "Q"], total_Q_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # P and Q are the average from power, excluding cable losses
                    push!(idx_p_mean_cal, i)
                    push!(idx_q_mean_cal, i)

                else

                    # all variable -
                    set_bounds(nodes[i, "P"], total_P_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"])
                    set_bounds(nodes[i, "Q"], total_Q_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # P and Q are the average from power, excluding cable losses
                    set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
                    set_bounds(nodes[i, "v"], parameters["grid"]["v_rms"], 0.95*parameters["grid"]["v_rms"], 1.05*parameters["grid"]["v_rms"])
                    push!(idx_p_mean_cal, i)
                    push!(idx_q_mean_cal, i)
                end

            else
                # all variable -
                set_bounds(nodes[i, "P"], total_P_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"])
                set_bounds(nodes[i, "Q"], total_Q_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # P and Q are the average from power, excluding cable losses
                set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
                set_bounds(nodes[i, "v"], parameters["grid"]["v_rms"], 0.95*parameters["grid"]["v_rms"], 1.05*parameters["grid"]["v_rms"])
                push!(idx_p_mean_cal, i)
                push!(idx_q_mean_cal, i)
            end
            #end
        else
            S = parameters["load"][i-num_source]["pwr"]/parameters["grid"]["phase"]
            P = S *parameters["load"][i-num_source]["pf"]
            #println("P_LOAD = $(P)")
            #println("Q_LOAD = $(sqrt(S^2 - P^2))")
            fix(nodes[i, "P"], -P)
            fix(nodes[i, "Q"], -sqrt(S^2 - P^2))

            set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
            set_bounds(nodes[i, "v"], parameters["grid"]["v_rms"], 0.95*parameters["grid"]["v_rms"], 1.05*parameters["grid"]["v_rms"])
        end
    end

    cable_cons = get_cable_connections(CM)
    node_cons = get_node_connections(CM)

    G = Array{NonlinearExpression, 2}(undef, num_nodes, num_nodes) # should be symmetric
    B = Array{NonlinearExpression, 2}(undef, num_nodes, num_nodes) # should be symmetric

    P_node = Array{NonlinearConstraintRef, 1}(undef, num_nodes)
    Q_node = Array{NonlinearConstraintRef, 1}(undef, num_nodes)

    # As radius goes down resistance goes up, inductance goes up, capacitance goes down. Put in formulas for this.
    #@variable(model, cables[1 : num_cables, ["L", "X_R", "C_L"]])
    @variable(model, cables[1 : num_cables, ["radius", "D-to-neutral"]])

    L_cable = Array{NonlinearExpression, 1}(undef, num_cables)
    R_cable = Array{NonlinearExpression, 1}(undef, num_cables)
    C_cable = Array{NonlinearExpression, 1}(undef, num_cables)

    cable_conductance = Array{NonlinearExpression, 1}(undef, num_cables)
    cable_susceptance_0 = Array{NonlinearExpression, 1}(undef, num_cables) # diagonals - where we add capacitances
    cable_susceptance_1 = Array{NonlinearExpression, 1}(undef, num_cables) # off diagonals


    # fixed distance between cable and neutral [2*radius + η - 1.00] m
    D = 0.5 #m

    for i=1:num_cables

        radius_max = 60*(4.1148e-3)/2
        set_bounds(cables[i, "radius"], 25*(3e-3)/2, (2.05232e-3)/2, radius_max) #m
        # set_bounds(cables[i, "radius"], (3e-3)/2, (3e-3)/2, (3e-3)/2) #m
        # assumption: min value of D-to-neutral : 3 * max radius
        set_bounds(cables[i, "D-to-neutral"], 2* 0.7788 * radius_max, 0.7788 * radius_max, 2.00 ) #m
        # assumption to line to line(neutral) --  for low voltages
        #println(parameters["cable"][i]["len"])
        L_cable[i] = @NLexpression(model, parameters["cable"][i]["len"] * 4e-7 * log(cables[i, "D-to-neutral"]/(0.7788 * cables[i, "radius"])))  # m* H/m

        # resistivity remains constant ρ_(T=50) = 1.973e-8
        R_cable[i] = @NLexpression(model, parameters["cable"][i]["len"] * 1.973e-8 / (π * cables[i, "radius"]^2)) # m * Ω/m

        # X_R = 0.38 # ratio of omega*L/R
        # X_R = omega * L(r)/R(r)

        # line to neutral
        C_cable[i] = @NLexpression(model, parameters["cable"][i]["len"] * 2π * 8.854e-12 / (log(cables[i, "D-to-neutral"]/cables[i, "radius"]))) #m * F/m

        cable_conductance[i] = @NLexpression(model, (R_cable[i] / ((R_cable[i])^2 + (omega * L_cable[i])^2)))
        cable_susceptance_1[i] = @NLexpression(model, (-omega * L_cable[i] / (((R_cable[i])^2 + (omega * L_cable[i])^2))))
        cable_susceptance_0[i] = @NLexpression(model, (-omega * L_cable[i] / (((R_cable[i])^2 + (omega * L_cable[i])^2)) + omega*C_cable[i]/2))
    end

    for i in 1:num_nodes

        # diagonal terms
        G[i, i] = @NLexpression(model, sum( cable_conductance[cable_cons[i]][j] for j in 1:length(cable_cons[i])))
        B[i, i] = @NLexpression(model, sum( cable_susceptance_0[cable_cons[i]][j] for j in 1:length(cable_cons[i])))

        # off diagonal terms
        for k in (i+1):num_nodes # this is over the upper triangle

            if CM[i, k] != 0

                cable_num = Int(abs(CM[i, k]))

                G[i, k] = @NLexpression(model, -1*cable_conductance[cable_num])
                B[i, k] = @NLexpression(model, -1*cable_susceptance_1[cable_num])
                G[k, i] = G[i, k]
                B[k, i] = B[i, k]

            else

                G[i, k] = zero_expression # a formula which returns 0.0
                B[i, k] = zero_expression
                G[k, i] = zero_expression
                B[k, i] = zero_expression
            end
        end
    end

    # power flow constraints - this is perfect!!
    for i in 1:num_nodes

        P_node[i] = @NLconstraint(model,

        nodes[i, "P"] == nodes[i,"v"] * sum( nodes[j,"v"] * ((G[i, j] * cos(nodes[i,"theta"] - nodes[j,"theta"]) + B[i, j] * sin(nodes[i,"theta"] - nodes[j,"theta"]))) for j in node_cons[i])

        )

        Q_node[i] = @NLconstraint(model,

        nodes[i, "Q"] == nodes[i,"v"] * sum( nodes[j,"v"] * ((G[i, j] * sin(nodes[i,"theta"] - nodes[j,"theta"]) - B[i, j] * cos(nodes[i,"theta"] - nodes[j,"theta"]))) for j in node_cons[i])

        )

    end

    cable_constraints = Array{NonlinearConstraintRef, 1}(undef, num_cables)

    for i in 1:num_cables

        j, k = Tuple(findfirst(x -> x == i, CM))

        cable_constraints[i] = @NLconstraint(model,
            abs( nodes[j, "v"] * nodes[k, "v"] * (sin(nodes[j, "theta"] - nodes[k, "theta"]))/(omega*L_cable[i]))
            <= 0.93 * nodes[j, "v"] * nodes[k, "v"] * sqrt(C_cable[i]/L_cable[i])
        )

    end

    # non-linear objectives
    @NLexpression(model, P_source_mean, sum(nodes[Int(j),"P"] for j in idx_p_mean_cal) / convert.(Int64,length(idx_p_mean_cal)))
    @NLexpression(model, Q_source_mean, sum(nodes[Int(j),"Q"] for j in idx_q_mean_cal) / convert.(Int64,length(idx_q_mean_cal)))

    @NLexpression(model, v_mean, sum(nodes[j,"v"] for j in 1:num_nodes) / num_nodes)

    # add apparent power
    @NLexpression(model, Power_apparent, sum(sqrt(nodes[i, "P"]^2 + nodes[i, "Q"]^2) for i in 1:num_source))
    # TODO: normalisation, i.e. weighting between minimising P and Q, and minimising cable radius? Maybe use per unit system
    # normalisation : 1. max value
    #                 2. p.u. 0
    # Sbase_1_phase = sum(loads)
    # Vbase_rms = 230
    # Ibase_rms = f(Sbase_1_phase, Vbase_rms)
    # Zbase = f(Vbase_rms, Ibase_rms)

    radius_upper_bound = upper_bound(cables[1, "radius"]);
    radius_lower_bound = lower_bound(cables[1, "radius"]);

    # Lagrangians/weights # TODO make these parameters depending on price ($)?
    λ₁ = 0.99
    λ₂ = 0.00001

    norm_P = length(idx_p_mean_cal)
    norm_Q = length(idx_q_mean_cal)
    @NLobjective(model, Min, λ₁ *  Power_apparent / total_S_source
                            + abs((v_mean - parameters["grid"]["v_rms"]) / parameters["grid"]["v_rms"])
                            #+ abs(sum(nodes[i,"theta"] for i in 2:num_nodes))/π    # maybe helpfull?
                            + sum( ((nodes[i,"P"] - P_source_mean)^2 )/ norm_P for i in idx_p_mean_cal)
                            + sum( ((nodes[i,"Q"] - Q_source_mean)^2) / norm_Q for i in idx_q_mean_cal) # the variance - not exactly right (but good enough)
                            + λ₂ * sum( (cables[i, "radius"]  for i in 1:num_cables)) / (num_cables * (radius_upper_bound - radius_lower_bound) )
                            #D-to-neutral to be minimized
                            )

    set_silent(model)    # turns off output of the solver #TODO put to verbosity?
    optimize!(model)

    # TODO: put these warning/messages with logger // verbosity
    if verbosity > 1
        println("CABLE LAYOUT BASED ON POWER FLOW EQUATIONS:")
        println("""
        termination_status = $(termination_status(model))
        primal_status      = $(primal_status(model))
        objective_value    = $(objective_value(model))
        """)
        global optimizer_status = get_optimizer_status(model)
        @show optimizer_status["termination_status"]


        println()
        println()
        println(value.(nodes))
        println("asdasdasd")
        println()
        println(value.(cables))
    end
    #= OPTIMAL, LOCALLY_SOLVED, INFEASIBLE, DUAL_INFEASIBLE, and TIME_LIMIT.
    if termination_status(model) == OPTIMAL
        println("Solution is optimal")
    elseif termination_status(model) == TIME_LIMIT && has_values(model)
        println("Solution is suboptimal due to a time limit, but a primal solution is available")
    else
        error("The model was not solved correctly.")
    end
    =#
    if termination_status(model) in [OPTIMAL, LOCALLY_SOLVED]
        # if solver was successful, we use the values
        for (index, cable) in enumerate(parameters["cable"])

            cable["L"] = value.(L_cable)[index]
            cable["Lb"] = cable["L"]/cable["len"]

            cable["R"] = value.(R_cable)[index]
            cable["Rb"] = cable["R"]/cable["len"]

            cable["C"] = value.(C_cable)[index]
            cable["Cb"] = cable["C"]/cable["len"]
        end

        for i in 1:num_cables

            j, k = Tuple(findfirst(x -> x == i, CM))

            #= Legacy code
                a =  sqrt(value.(C_cable)[i]/value.(L_cable)[i])
                #println(omega*value.(L_cable)[i])
                println()
                println()
                println(value.(B)[j,j])
                println(value.(G)[j,j])
                println()
                println()

                Y = 1/(value.(R_cable)[i] + omega*value.(L_cable)[i])
                V1 = value.(nodes[k, "v"]) *exp(1im*value.(nodes[k, "theta"]) ) # theta should be replaced with either 0, δ₁ or δ₂
                V2 = value.(nodes[j, "v"]) *exp(1im*value.(nodes[j, "theta"]) ) # theta should be replaced with either 0, δ₁ or δ₂
                println()
                println("Aparent power:")
                println(conj(Y)*conj(V1-V2)*V1)

                println()
                println()

                println(value.(nodes[j, "v"]) * value.(nodes[k, "v"]) * (sin(value.(nodes[j, "theta"]) - value.(nodes[k, "theta"])))/(omega*value.(L_cable)[i]))

                println(value.(nodes[j, "v"]) * value.(nodes[k, "v"]) * (sin(value.(nodes[j, "theta"]) - value.(nodes[k, "theta"]))/(value.(B)[j,k])   +    cos(value.(nodes[j, "theta"]) - value.(nodes[k, "theta"]))/(value.(G)[j,k])))
                println()
                println()

                dn = asin(mod(-a/(omega*value.(L_cable)[i]),2*pi))

                I = min(value.(nodes[j, "v"]), value.(nodes[k, "v"])) * ((omega*value.(L_cable)[i])*sin(dn))
                println(I)
            =#

            Z = value.(R_cable)[i] + 1im*omega*value.(L_cable)[i]
            Zₘ = abs(Z)
            θᵧ = angle(Z)

            Y = 1im*omega*value.(C_cable)[i] # C is the total capacitance of the line, not the halved capacitance
            A = 1 + Y*Z/2
            Aₘ = abs(A)
            θₐ = angle(A)

            vᵣ = value.(nodes[k, "v"]) # magnitude of receiving end voltage
            vₛ = value.(nodes[j, "v"]) # magnitude of sending end voltage

            P = 3.0*vᵣ*vₛ*sqrt(value.(C_cable)[i]/value.(L_cable)[i])
            #P = value.(nodes[j, "P"]) # maybe should be value.(nodes[k, "P"])
            #Q = value.(nodes[j, "Q"])

            # println(Zₘ)
            # println(P)
            # println(value.(C_cable))
            # println(value.(L_cable))
            # println(Aₘ)
            # println(vᵣ)
            # println(θᵧ)
            # println(θₐ)
            # println(vₛ)


            δ = -acos((P*Zₘ + Aₘ*vᵣ*vᵣ*cos(θᵧ - θₐ))/(vᵣ*vₛ)) + θᵧ

            Vr = vᵣ # magnitude of receiving end voltage - set angle to 0.0
            Vs = vₛ*exp(1im*δ) # magnitude and angle of sending end voltage

            Yₗ = 1/Z # for debugging
            Iₗ = (conj(Yₗ)*(Vr - Vs)) # this is almost our answer - the limit through the inductor

            I₂ = ((Vs - A*Vr)/Z)

            #println("Iₗ = ", abs(Iₗ), " This is our answer. The limit through the inductor")

            if !haskey(parameters["cable"][i], "i_limit")
                parameters["cable"][i]["i_limit"] = abs(Iₗ)
            end

            #=
            S = sqrt(P^2 + Q^2)
            println("\nDebugging\n")
            println("1. Cable = ", i)
            println("2. δ = ", δ, " ?= ", value.(nodes[j, "theta"]) - value.(nodes[k, "theta"]))
            println("3. S = ", S, " S₁ ?= ", abs(Vs*conj(Yₗ)*(Vr - Vs)), " ?= ", abs(Vr*conj(Yₗ)*(Vs - Vr)))
            println("4. S = ", S, " S₂ ?= ", abs(Vs*Iₗ), " ?= ", abs(Vs*I₂), " ?= ", abs(Vr*Iₗ), " ?= ", abs(Vr*I₂))
            println("5. |Iₗ| = ", abs(Iₗ), " angle = ", angle(Iₗ)*180/pi, " Iₗ = ", Iₗ, " =? ", abs(conj(Yₗ)*(Vr - Vs)), " =? ", abs(conj(Yₗ)*(Vs - Vr)), " this is probably the right one")
            println("6. |I₂| = ", abs(I₂), " angle = ", angle(I₂)*180/pi, " I₂ = ", I₂, " =? ", abs(conj(Yₗ)*(Vr - Vs)), " =? ", abs(conj(Yₗ)*(Vs - Vr)))
            println("7. vᵣ = ", value.(nodes[j, "v"]))
            println("8. vₛ = ", value.(nodes[k, "v"]))
            println("9. P = ", P, " ?= ", vᵣ*vₛ*cos(θᵧ - δ)/(Zₘ) - Aₘ*vᵣ*vᵣ*cos(θᵧ - θₐ)/(Zₘ), " ?= ", real(Vs*I₂), " ?= ", real(Vr*I₂), " ?= ", real(Vs*Iₗ), " ?= ", real(Vr*Iₗ))
            println("10. Q = ", Q, " ?= ", vᵣ*vₛ*sin(θᵧ - δ)/(Zₘ) - Aₘ*vᵣ*vᵣ*sin(θᵧ - θₐ)/(Zₘ), " ?= ", imag(Vs*I₂), " ?= ", imag(Vr*I₂), " ?= ", imag(Vs*Iₗ), " ?= ", imag(Vr*Iₗ))
            println("11. R = ", value.(R_cable)[i] , " L = ", value.(L_cable)[i], " C = ", value.(C_cable)[i])
            println("12. Vs = ", Vs , " Vr = ", Vr)
            println()
            =#

        end
    else
        if verbosity > 0
            @warn("Power flow equation not solveable! Maybe parameter setting invalid.
                  Default values are used for cable parameters")
        end
    end

    return parameters
end

"""
    R, L = Fault_Level(S, X_R, Vrms; fsys = 50, phase = 3)

# Arguments
- `S::Float`: 3 phase Fault Level [VA]
- `X/R::Float`: the (Short Circuit) ratio of reactance to resistance
- `Vrms::Float`: Line to Neutral rms voltage of external network [V]

# Keyword Arguments
- `fsys::Float`: system frequency [Hz]
- `phase::Int`: single or 3 phase system, i.e. phase = 1, or phase = 3

# Return Values
- `R::Float`: effective resistance [Ω]
- `L::Float`: effective inductance [H]

# Theory
An external network is often characterised by its Fault level and its X/R ratio. In particular the
Fault Level is a measure of the strength of a network. It is the amount of apparent power that the
network can supply when a three-phase bolted to ground fault is applied at the point of connection.
From these values an effective resistance and inductance can be calculated. Typical values of X/R
ratios are in the range of 0.5 to 1.5, for distribution networks. Transmission networks operating
at higher voltages tend to have higher X/R ratios.

"""
function Fault_Level(S, X_R, Vrms; fsys = 50, phase = 3)

    rad = atan(X_R)

    I_fault = S/(phase*Vrms)
    Z = Vrms/I_fault

    R = Z*cos(rad)
    X = Z*sin(rad)

    L = X/(2π*fsys)

    return R, L
end
