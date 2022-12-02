using SpecialFunctions
using LinearAlgebra
using Plots
using StatsBase
using LinearAlgebra
using FFTW
using Expokit #for exponentials of sparse matrices
#using Polynomials
#using CurveFit

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#= The Frequency Dependent Model
    Transmission lines differ from ordinary electric networks in one essential feature.
    Whereas the physical dimensions of electrical networks are very much smaller than
    the operating wavelength, transmission lines are usually a considerable fraction
    of a wavelength and may even be many wavelengths long. The circuit elements in an
    ordinary electric network can be considered discrete and as such may be described
    by lumped parameters. It is assumed that currents flowing in lumped-circuit elements
    do not vary spatially over the elements,and that no standing waves exist. A
    transmission line, on the other hand, is a distributed-parameter network and must
    be described by circuit parameters that are distributed throughout its length.
    Except under matched conditions, standing waves exist in a transmission line.
    A frequency dependent transmission line model is used when the physical geometry
    of the line is available (i.e. conductor radius and positions), and the travelling
    time is greater than the time step.

    The line frequency-dependent surge impedance (or admittance) and line propagation
    matrix are first calculated from the physical line geometry. To obtain the time
    domain response, a convolution must be performed as this is equivalent to a
    multiplication in the frequency domain. It can be achieved efficiently using
    recursive convolutions (which can be shown to be a form of root-matching, even
    though this is not generally recognised). This is performed by fitting a rational
    function in the frequency domain to both the frequency-dependent surge impedance
    and propagation constant.

    As the line parameters are functions of frequency, the relevant equations should
    first be viewed in the frequency domain, making extensive use of curve fitting to
    incorporate the frequency-dependent parameters into the model. Two important
    frequency-dependent parameters influencing wave propaation are the characteristic
    impedance Zc and propagation constant γ. However, quantities such as reactances,
    wavelength, wavenumber, and phase constant would loose their meaning under
    transient conditions. Especially when considering sudden surges, produced by
    breakers opening or closing like when energising the lines, or lightning strikes
    that terminate on the conductors or cause back flash overs. Therefore, rather
    than looking at Zc and γ in the frequency domain and considering each frequency
    independently, they are expressed by a continuous function of frequency that need
    to be approximated by a fitted rational function.

    The frequency dependence of the series impedance is most pronounced in the zero
    sequence mode, thus making frequency-dependent line models more important for
    transients where appreciable zero sequence voltages and zero sequence currents
    exist, such as in single line-to-ground faults.
=#

#%% Cable/trasmission line parameters
#_______________________________________________________________________________
F = 49.0:1:50.0 #don't remove decimal points
l = 1.0 #km, length of the cable/transmission line in km

num_conductors = 3
num_neutrals = 0

transposed = 0 # 1 if transposed

#= Neutrals/Earth/Ground/Shield
    All the neutral conductors are connected in parallel and are grounded at regular intervals
    along the line. Any isolated neutral conductors that carry no current are omitted.The
    phase conductors are insulated from each other and from earth. If the phase currents
    are not balanced, there may be a return current in the grounded neutral wires and in the earth.

    When earth wires are continuous and grounded at each tower then for frequencies
    below 250kHz it is reasonable to assume that the earth wire potential is zero
    along its length.
=#
num_total = num_neutrals + num_conductors

num_bundles = Array{Int8, 1}(undef, num_total*2)
num_bundles = fill!(num_bundles, 1)
gmr_cond = Array{Float64, 1}(undef, num_total*2)
r_cond = Array{Float64, 1}(undef, num_total*2)
earth_resistance = Array{Float64, 1}(undef, length(F))

# Conductors #current data for 19strand 4/0 hard drawn copper
num_bundles[1] = 1
gmr_cond[1] = 5.084064  #mm, geometric mean radius of stranded conductor
r_cond[1] = 6.7056  #mm, outside radius of (stranded) conductor

num_bundles[2] = 1
gmr_cond[2] = 5.084064 #mm, gmr radius of conductor
r_cond[2] = 6.7056 #mm, outside radius of (stranded) conductor

num_bundles[3] = 1
gmr_cond[3] = 5.084064 #mm, gmr radius of conductor
r_cond[3] = 6.7056 #mm, outside radius of (stranded) conductor

#= Neutral #should be 0.375 EBB steel
num_bundles[4] = 1
gmr_cond[4] = 5.084064   #mm, gmr radius of conductor #
r_cond[4] = 6.7056   #mm, outside radius of (stranded) conductor=#

x_pos = Array{Float64, 2}(undef, 2*num_total, maximum(num_bundles))
y_pos = Array{Float64, 3}(undef, 2*num_total, maximum(num_bundles), length(F)) #meters above ground
y_pos_m = Array{Float64, 2}(undef, 2*num_total, maximum(num_bundles)) #meters above ground
x_pos = fill!(x_pos, +100) #+100 is just a placeholder
y_pos = fill!(y_pos, +100)
y_pos_m = fill!(y_pos, +100)

# Conductor Positions # current dat for 69kV circuit
# y_pos = y_pos - 2/3*y_sag
x_pos[1,1] = -3.048 #m
y_pos[1,1,1] = 13.716 #m
#[conductor/neutral number, bundle number, reserved for earth frequency]

x_pos[2,1] = 0 #m
y_pos[2,1,1] = 13.716 #m

x_pos[3,1] = 3.048 #m
y_pos[3,1,1] = 13.716 #m

#= Neutral Positions
x_pos[4,1] = 0 #m
y_pos[4,1,1] = 18.288 #m =#

# System parameters
#_______________________________________________________________________________
f0 = 50
V = 230

#%% Material properties
#_______________________________________________________________________________
μe = 4.263e-3 #m^2/Vs, electron mobilitty
μr = 0.9999994 #relative Permeability of copper
rg = 100 #Ωm, ground resistivity
ρ = 8.96 #g/cm^3, density of the material
x = 1 # number of free electrons per atom, the valence number
M = 63.546 #g/mol, the molar mass
#lp = 3.6147e-10 #lattice parameter of copper in m

#= Ions and Valence Electrons of Metals
    In a metal, the valence electrons are thought of as being shared by all the
    positive ions. Therefore, the electrons are free to move throughout the crystalline
    lattice. The electrons move randomly throughout the crystal, until an electric
    field is applied to the material. Then the electric field forces the electrons to
    move in a direction opposite to the field. The electrons still move somewhat
    randomly, but with a superimposed "drift". This produces current.

    As the temperature of a semiconductor is increased, more valence electrons acquire
    sufficient energy to move to the conduction band, so more current flows.
=#

T_0 = 20 #Temperature in degrees for which μe's value is true
T_op = 32 #Operating temperature in degrees Celcius
T_const = 241.5 #Temperature constant in degrees Celcius
#T_coef = 0.393/100 #Temperature coefficient in percentage per degree celcius (also valid*)

#%% Universal constants
#_______________________________________________________________________________
Na = 6.0221e23 #Avogadro's constant
qe = 1.60217663e-19 #charge of an electron
kB = 1.3806503e-23 #Boltzmann's constant
ϵ0 = 8.85418782e-12 #Permittivity of free space
ϵr = 1.0006 #Relative Permittivity of air
μ0 = 1.25663706e-6 #Permeability of free space

#%% Material equations
#_______________________________________________________________________________
μ = μ0*μr
n = Na*(ρ*1000)*x/(M*0.001) #electrons/m^3, estimating the number of charge carriers in a volume unit
#n = (4*1)/(lp^3) # Copper has a FCC structure with 4 atoms per cell, 8.47e28 m^3 for copper (also valid)
ρe = -n*qe #charge density of the drifting electrons, should be negative

#= Conductivity
    Conductivity is a macroscopic constitutive parameter.

    σ = -ρe*μe + ρh*μh

    The reciprocal of conductivity is called resistivity, in ohm-meters.
    Resistivity depends on the conductor metal. Annealed coper is the international
    standard for measuring resistivity. Hard-drawn aluminum, which has 61% of the
    conductivity of the international standard, has a resistivity at 20 degC of
    2.83e-8 Ωm.

    Conductor resistance depends on the following factors:
    1. Spiraling
    2. temperature
    3. Frequency "Skin Effect"
    4. Current magnitude - magnetic conductors

    If the nuclei are arranged in a perfectly ordered lattice then there should be no
    resistance to oppose the flow of current. Also an increase in temperature causes
    an increase in the thermal population of the higher levels and it should increase
    the eletrical conductivity. But, thermal vibrations of the nuclei increases electrical
    resistance so conductivity actually decreases with temperature.

    Due to impurities and defects in the metal, there will always be a residual resistivity,
    which is not temperature dependent. According to Matthiessen's rule, the resistivity
    arises from independent scattering processes which are additive.

    For dc, the current distribution is uniform throughout the conductor cross-section.
    However, for ac, the current distribution is nonuniform. As frequency increases, the
    current in a solid cylindrical conductor tends to crowd toward the conductor surface,
    with smaller current density at the conductor center. This phenomenon is called
    skin effect. A conductor with a large radius can even have an oscillatory current
    density versus the radial distance from the conductor center.

    With increasing frequency, conductor loss increases which causes the ac resistance
    to increase.

    For magnetic conductors, such as steel conductors used for shield wires, resistance
    depends on current magnitude. The internal flux linkages, and therefore the iron
    or magnetic losses, depend on the current magnitude.

    If we assume that the conductivity of the conductors in a transmission line is high
    enough that the effect of the series resistance on the computation of the propagation
    constant is negligible, the implication being that the waves on the line are
    approximately Transverse Electromagnetic Waves (TEM), where the electric and magnetic
    fields are perpendicular to each other, and both are transverse to the direction of
    propagation.
=#

σT_0 = -ρe*μe # S/m, or A/Vm, conductivity of copper 58.7e6.
rT_0 = 1/σT_0 #Ωm, resistivity of copper is 1.72e-8

#= The effects of heat
    As heat is applied to a crystalline solid, on the atomic level, the kinetic energy
    of the atoms has increased which means the atoms are moving faster. However, in a
    crystalline solid, the atomic movement is limited to vibration around stable lattice
    positions. As the temperature increases, the atoms vibrate at a greater amplitude
    and move farther from their stable lattice positions. This motion has a negative
    effect on the ability of the material to conduct an electrical current, causing
    the resistivity to increase.

    As the temperature increases, the positive ions in the crystal vibrate more, and
    more collisions occur between the valence electrons and the vibrating ions. These
    collisions hinder the "drift" motion of the valence electrons, thus reducing the
    current.

    Resistivity of conductor metals varies linearly over normal opering conditions.
    The resistivity, among other things,depends on the relaxation time. It's the
    average time interval between two successive collisions between an electron and
    the cation in a conductor.

    When the temperature of the metal increases, some electrons gain energy and are
    excited inot empty energy levels in a valence band. This condition creates an
    equal number of empty energy levels, or holes, vacated by the excited electrons.
    Only a small increase in energy is required to cause the excitation of electrons.
    Both excited electrons (free electrons) and the newly created holes can then carry
    an electrical charge. when the electrons are excited into unfilled levels, the
    Fermi energy is unchanged.

    Fermi energy is the energy corresponding to the highest filled state at 0 degrees
    Kelvin. Only electrons above the Fermi energy can be affected by an electric field.

    These are the free electrons.
=#
rT_op = rT_0*((T_op + T_const)/(T_0 + T_const))
#rT_op = rT_0*(1 + T_coef*(T_op - T_0)) # Same as above (also valid*)

σT_op = σT_0*((T_0 + T_const)/(T_op + T_const))

for i in 1:num_total

    gmr_cond[i] = gmr_cond[i]/1000
    gmr_cond[i+num_total] = gmr_cond[i]
    r_cond[i] = r_cond[i]/1000
    r_cond[i+num_total] = r_cond[i]
end

R_dc_c = Array{Float64, 1}(undef, num_total)
S_c = π*(gmr_cond).^2 #m^2, surface area of conductor

# the formula for the resistance of a straight piece of homogeneous material of
# a uniform cross section for steady d.c. conduction current
R_dc_c = 1 ./ (σT_op.*S_c) #Ω/m, conductor resistance per meter for a uniform current
#G_dc = 1/R_dc #Ω/m, Conductance per meter for a uniform current

# if DC resistance per unit length is given at T_0 (usually 20 degrees)
#R_dc_c[1] = 0.00017274119137600003 #Ω/m, for 19 strand 4/0 hard drawn copper
σT_0 = 1/(S_c[1]*R_dc_c[1])
σT_op = σT_0*((T_0 + T_const)/(T_op + T_const))
rT_op = 1/σT_op

#= Diffusion
    The Diffusion constant,
    the Einstein-Smoluchowski relation connects the diffusion constant with the
    electrical mobility. Electrons in a wire are in constant, thermal motion. If we
    imagine putting all the electrons in a small region of a wire, the thermal motion,
    attractive and repulsive coulomb force, quickly spreads them throuhout the whole wire.
    The diffusion constant tells us how quickly this happens. The unit of diffusion
    is area/time. Say that, at some moment in time electrons occupy a certain area.
    The diffusion constant is the velocity of growth in time of this area.
=#

Dif = μe*kB*(T_op + 273.15)/qe*1e6 #mm^2,s

f = 50
α = sqrt(π*μr*μ0*σT_op*f) #Np/m, attenuation constant (approximate)
β = α
#rad/m, phase constant expresses the amount of phase shift that occurs as the wave travels on meter
ηc = (1 + 1im)*α/σT_op
# intrinsic impedance of a good conductor, hence the magnetic field intensity lags
# behind the electric field intensity by 45 degrees
μp = 2*π*f/β #m/s, phase velocity in a good conductor
λ = 2000*π/β #mm, wavelength
δ = 1/α #m, skin depth or depth of penetration

#= Skin Depth
    The amplitude of a wave will be attenuated by 36.8% when a wave travels δ into
    the conductor. Thus at high-frequenies the skin depth or depth of penetration is
    so small that field and current can be considered as, for all practical purposes
    confined in a very thin layer of the conductor surface.
=#

# The internal inductace per unit length of conductor
L_dc = 0.5e-7 #H/m

#%% GMD and GMR
#_______________________________________________________________________________

#= Composite Conductors
    Bundles and stranded subconductors (composite conductors)

    For stranded conductors, alternate layers of strands are spiraled in opposite
    directions to hold the strands together. Spiraling makes the strands 1% or 2%
    longer than the actual conductor length. As a result, the dc resistance of a
    stranded conductor is 1 or 2% lorger.

    It is common practice for EHV lines to use more than one conductor per phase, a
    practice called bundling. Bundling reduces the electric field strength at the
    conductor surfaces which in tur reduces or eliminates corona and its results;
    undesirable power loss, communications interference, and audible noise. Bundling
    also reduces the series reactance of the line by increasing the GMR of the bundle.
    Bundled subconductors are often used to reduce the electric field strength at the
    surface of the conductors, as compared to using one large conductor. This therefore
    reduces the likelihood of corona. The two alternative methods of modelling bundling are:
    1. Replace the bundled subconductors with an equivalent single conductor
    2. Explicitly represent subconductors and use matrix elemination of subconductors.
    The second method is more rigorous, since all subconductors are represented explictly.
    The use of GMR ignores proximity effects and hence is only valid if the subconductor
    spacing is much smaller than the spacing between the phases of the line.

    L_int_d = (μ0*μr)/(8*π) # H/m, the internal inductance per unit length
    The internal inductance arises from the flux linkage internal to the solid
    conductor. It is assumed that the current is uniformly distributed in the inner
    conductor. This assumption does not hold for high-frequency a.c. currents. This
    assumtion introduces a slight error for large conductors, even at power frequencies.
    This complication is due to the skin effect and proximity effect. Briefly stated,
    the skin effect causes the current distribution to become nonuniform, with a
    larger current density appearing on the conductor surface than at the center.
    This reduces the internal flux linkages and lowers the internal inductance as
    compared to the uniform current density (dc). It also increases the resistance.

    There are a number of ways to calculate the electrical parameters from the physical
    geometry of a line, the most common being Carson's series equations.

    Carson showed that the earth can be replaced by a set of "earth return" conductors
    located directly under the overhead the overhead conductors. Each earth return
    conductor carries the negative of its overhead conductor current.

    The GMR of each earth return conductor, is the same as the GMR of its corresponding
    overhead conductor. Also, all the earth return conductors have the same distance
    from their overhead conductors, and the same frequency dependent resistance.

    When shunt admittances (line-to-line capacitances) are concerned we follow a
    similar strategy, by approximating the earth surface as a perfectly conducting
    horizontal plane, even though the earth under the line may have irregular terrain
    and resistivities.

    The effect of the earth plane is accounted for by the method of images. Consider
    a single conductor wiht uniform charge distribution and with height H above a
    perfectly conducting earth plane. When the conductor has a positive charge, an
    equal quantity of negative charge is induced on the earth. The electric field
    lines will terminate at the negative charges on the earth. Also the electric field
    lines are perpendicular to the surfaces of the conductor and earth.

    The earth can be replaced by an image conductor, which has the same radius as the original
    conductor, and lying directly below the original conductor, with a separation of
    2H, and having a equal quantity of negative charge.

    In general, the effect of the earth plane is to slightly increase the capacitance,
    and as the line height increases, the effect of the earth will become negligible.
=#
l = l*1.02*1000 #if stranded

x_pos[(num_total+1):end, :] = x_pos[1:num_total, :]

for i in 1:num_total
    num_bundles[i + num_total] = num_bundles[i] #earth images
end

y_pos_m = y_pos[:,:,1]

for c in 1:num_total
    for b in 1:num_bundles[c]

        y_pos_m[c+num_total, b] = -1*y_pos_m[c, b]

        for f in eachindex(F)
            # all the earth conductors have the same distance from their overhead conductors
            y_pos[c, b, f] = y_pos[c, b, 1]
            y_pos[c+num_total, b, f] = y_pos[c, b, 1] - 658.5*sqrt(rg/F[f])
        end
    end
end

Dl = Array{Float64, 3}(undef, 2*num_total, 2*num_total, length(F))
Dl = fill!(Dl, 1)
Dc = Array{Float64, 2}(undef, 2*num_total, 2*num_total)
Dc = fill!(Dc, 1)

# Calculating the GMD and GMR of the conductors, neutrals, and their earth mirrors.
for f in eachindex(F)
    for c1 in 1:num_total*2
        for c2 in 1:num_total*2
            for b1 in 1:num_bundles[c1]
                for b2 in 1:num_bundles[c2]

                    if c1 == c2 # if we are in the conductors are the same

                        if b1 == b2 # if the bundles are the same

                            m = sqrt(2*π*μ*σT_op*F[f]) #sqrt(2)/δ, where δ is skin depth
                            mr = m*gmr_cond[c1]

                            # Calculating Kelvin function and their derivates for n=0
                            ber_mr = real(besselj(0,mr*exp(3*π*1im/4)))
                            bei_mr = imag(besselj(0,mr*exp(3*π*1im/4)))
                            ber_mr_1 = real(besselj(1,mr*exp(3*π*1im/4)))
                            bei_mr_1 = imag(besselj(1,mr*exp(3*π*1im/4)))
                            ber_d_mr = (ber_mr_1 + bei_mr_1)/sqrt(2)
                            bei_d_mr = (-ber_mr_1 + bei_mr_1)/sqrt(2)

                            αL = (4/mr)*(ber_mr*ber_d_mr + bei_mr*bei_d_mr)/((ber_d_mr)^2 + (bei_d_mr)^2) #H/m

                            Dl[c1, c2, f] = Dl[c1, c2, f]*exp(-αL*L_dc/(2e-7))*gmr_cond[c1]

                            if f == 1 # assuming that capacitance is not frequency dependent
                                Dc[c1, c2] = Dc[c1, c2]*r_cond[c1]
                            end

                        else # if the bundles aren't the same

                            Dl[c1, c2, f] = Dl[c1, c2, f]*sqrt((x_pos[c1, b1] - x_pos[c2, b2])^2 + (y_pos[c1, b1, f] - y_pos[c2, b2, f])^2)

                            if f == 1 # assuming that capacitance is not frequency dependent
                                Dc[c1, c2] = Dc[c1, c2]*sqrt((x_pos[c1, b1] - x_pos[c2, b2])^2 + (y_pos_m[c1, b1] - y_pos_m[c2, b2])^2)
                            end
                        end

                    elseif c1 != c2 # if the conductors aren't the same

                        Dl[c1, c2, f] = Dl[c1, c2, f]*sqrt((x_pos[c1, b1] - x_pos[c2, b2])^2 + (y_pos[c1, b1, f] - y_pos[c2, b2, f])^2)

                        if f == 1 # assuming that capacitance is not frequency dependent
                            Dc[c1, c2] = Dc[c1, c2]*sqrt((x_pos[c1, b1] - x_pos[c2, b2])^2 + (y_pos_m[c1, b1] - y_pos_m[c2, b2])^2)
                        end
                    end
                end
            end
        end
    end

    for c1 in 1:num_total*2
        for c2 in 1:num_total*2
            Dl[c1, c2, f] = (Dl[c1, c2, f])^(1/(num_bundles[c1]*num_bundles[c2]))
            if f == 1 #
                Dc[c1, c2] = (Dc[c1, c2])^(1/(num_bundles[c1]*num_bundles[c2]))
            end
        end
    end
end

#%% Finding the inductance, resistance, capacitance and impedance matrices
#_______________________________________________________________________________

R = Array{Float64, 3}(undef, num_total, num_total, length(F)) #Ω/m
L = Array{Float64, 3}(undef, num_total, num_total, length(F))
Z = Array{ComplexF64, 3}(undef, num_total, num_total, length(F)) #Ω/m
P = Array{Float64, 2}(undef, num_total, num_total) #m/F

for f in eachindex(F)

    earth_resistance[f] = 9.869e-7*F[f]

    for c1 in 1:num_total
        for c2 in 1:num_total
            if c1 == c2
                S = π*(gmr_cond[c1])^2 #m^2, effective surface area of bundle
                # the formula for the resistance of a straight piece of homogeneous material of
                # a uniform cross section for steady d.c. conduction current
                R_dc = 1/(σT_op*S) #Ω/m, bundle resistance per meter for a uniform current
                #G_dc = 1/R_dc #Ω/m, Conductance per meter for a uniform current

                m = sqrt(2*π*μ*σT_op*F[f]) #sqrt(2)/δ, where δ is skin depth
                mr = m*Dl[c1, c1, 1]

                # Calculating Kelvin function and their derivates for n=0
                ber_mr = real(besselj(0,mr*exp(3*π*1im/4)))
                bei_mr = imag(besselj(0,mr*exp(3*π*1im/4)))
                ber_mr_1 = real(besselj(1,mr*exp(3*π*1im/4)))
                bei_mr_1 = imag(besselj(1,mr*exp(3*π*1im/4)))
                ber_d_mr = (ber_mr_1 + bei_mr_1)/sqrt(2)
                bei_d_mr = (-ber_mr_1 + bei_mr_1)/sqrt(2)

                R_cond = R_dc*(mr/2)*(ber_mr*bei_d_mr - bei_mr*ber_d_mr)/((ber_d_mr)^2 + (bei_d_mr)^2) #Ω/m

                R[c1, c2, f] = R_cond + earth_resistance[f]
            else
                R[c1, c2, f] = earth_resistance[f]
            end

            L[c1, c2, f] = 2e-7*log(Dl[c1, c2 + num_total, f]/Dl[c1, c2, f])
            Z[c1, c2, f] = R[c1, c2, f] + 2im*π*F[f]*L[c1, c2, f]

            if f == 1
                P[c1, c2] = (1/(2*π*ϵ0*ϵr))*log(Dc[c1, c2 + num_total]/Dc[c1, c2])
            end
        end
    end
end

# removing the neutral
Zp = Array{ComplexF64, 3}(undef, num_conductors, num_conductors, length(F)) #Ω/m
Cp = Array{ComplexF64, 2}(undef, num_conductors, num_conductors) #m/F

# P is also known as Maxwell's potential coefficient matrix
PA = P[1:num_conductors, 1:num_conductors]
PB = P[1:num_conductors, (num_conductors+1):end]
PC = P[(num_conductors+1):end, 1:num_conductors]
PD = P[(num_conductors+1):end, (num_conductors+1):end]

Cp = inv(PA - PB*inv(PD)*PC) #F/m
#=
    Capacitances in distrubution systems are too small and are normally neglected.

    The shunt admittance of a line consists of the conductance and the capacitive
    susceptance. The conductance is usually ignored because it is very small compared
    to the capacitive susceptance. The capacitance of a line is the result of the
    potential difference between conductors. A charged conductor creates an electric
    field that emanates outward from the center.
=#
for f in eachindex(F)

    ZA = Z[1:num_conductors, 1:num_conductors, f]
    ZB = Z[1:num_conductors, (num_conductors+1):end, f]
    ZC = Z[(num_conductors+1):end, 1:num_conductors, f]
    ZD = Z[(num_conductors+1):end, (num_conductors+1):end, f]

    Zp[:,:,f] = ZA - ZB*inv(ZD)*ZC #Ω/m
end

#= Transpositions and Twists of Line Conductors
    It is apparent that the phase conductors of a three-phase circuit are mutually
    coupled and that currents in any one conductor will produce voltage drops in the
    adjacent conductors. Furthermore, these induced voltage drops may be unequal even
    for balanced currents since the mutual impedances depend entirely on the physical
    arrangement of the wires.

    One means of equalizing the mutual inductances is the contruct transpositions or
    rotations of the overhead line wires. A transposition is a physical rotation of
    the conductors, arranged so that each conductor is moved to occupy the next physical
    position in a regular sequence such as a-b-c, b-c-a, c-a-b, etc.
=#
#_______________________________________________________________________________
if transposed == 1
    for f in 1:length(F)

        off_diag_Z = (1/num_conductors)*sum(tril((Zp[:,:,f]), -1))
        diag_Z = (1/num_conductors)*tr(Zp[:,:,f])

        off_diag_C = (1/num_conductors)*sum(tril((Cp[:,:]), -1))
        diag_C = (1/num_conductors)*tr(Cp[:,:])

        for i in 1:num_conductors
            for j in i:num_conductors
                if i == j
                    Zp[i,j,f] = diag_Z

                    if f == 1
                        Cp[i,j] = diag_C
                    end
                elseif j > i
                    Zp[i,j,f] = off_diag_Z
                    Zp[j,i,f] = off_diag_Z

                    if f == 1
                        Cp[i,j] = off_diag_C
                        Cp[j,i] = off_diag_C
                    end
                end
            end
        end
    end
end

#Zp = Zp*1000*64.3738

#= Surge Impedance Loading
    For waves launched on an infinitely long line at the start, there can only be
    forward travelling waves, and no reflected waves. This is also true for finite
    lines terminated in a characteristic impedance; that is, when the lines are
    matched. From circuit theory we know that a maximum transfer of power from a given
    voltage source to a load occurs under "matched conditions" when the load impedance
    is the complex conjugate of the source (+line) impedance. In transmission line
    terminology, a line is matched when the load impedance is equal to the characteristic
    impedance (not the complex conjugae of the characteristic impedance) of the line.
    When a finite transmission line is terminated with its own characteristic impedance
    (when a finite transmission line is matched), the voltage and current distributions
    on the line are exactly the same as though the line has been extended to infinity.
=#

# The Characteristic Impedance and propation constant
Zc = Array{ComplexF64, 3}(undef, num_conductors, num_conductors, length(F)) #Ω/m
A = Array{ComplexF64, 3}(undef, num_conductors, num_conductors, length(F)) #Ω/m

for f in eachindex(F)
    for c1 in 1:num_conductors
        for c2 in 1:num_conductors
            Zc[c1,c2,f] = sqrt(Zp[c1,c2,f]/(2im*π*F[f]*Cp[c1,c2]))
            A[c1,c2,f] = sqrt(Zp[c1,c2,f]*2im*π*F[f]*Cp[c1,c2])
        end
    end
end

#=
Zd_w = Rd_w + 1im*w*Ld_w
Yd_w = G_w + 1im*w*C

Zc = sqrt((Zd_w)/(Yd_w))


f1 = curve_fit(RationalPoly, F, real.(Zp[1, 1, :]), 2, 3)
Fit_Zr = f1.(F)

p1 = plot(F, real.(Zp[1, 1, :]),
title = "Frequency Dependent Parameters",
xaxis = :log,
label = "Zp[1,1]",
xlabel = "Frequency [Hz]",
ylabel = "Resistance [Ω/m]")
p1 = plot!(F, Fit_Zr,
label = "Fitted Rational Polynomial")
#display(p1)

# The number of poles is normally one more than the zeros, as the attenuation
# funtion magnitude must to go zero when the frequency approaches infinity.
f2 = curve_fit(RationalPoly, F, real.(A[1, 1, :]), 2, 3)
Fit_Ar = f2.(F)

p1 = plot(log.(2*π*F)./log(10), exp.(-1*real.(A[1, 1, :])),
title = "Frequency Dependent Parameters",
#xaxis = :log,
label = "A[1,1]",
xlabel = "Frequency [Hz]",
ylabel = "Magnitude")
p1 = plot!(log.(2*π*F)./log(10), exp.(-1*Fit_Ar),
label = "Fitted Rational Polynomial")
display(p1)

p2 = plot(F, imag.(Zp[1, 1, :]),
title = "Frequency Dependent Parameters",
#xaxis = :log,
label = "Zp[1,1]",
xlabel = "Frequency [Hz]",
ylabel = "Reactance [Ω/m]")
#display(p2)
=#

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
