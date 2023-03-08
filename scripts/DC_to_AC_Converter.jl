using SpecialFunctions
using LinearAlgebra
using Plots
using StatsBase
using LinearAlgebra
using FFTW

fs = 2500 # the frequency of the triangle waveform
fs_p = 1 # the amplitude of the triangle wave
Nps = 20000 # time intervals, samples in a second
T = 1 # periods
Vd = 800 # DC Voltage
fr = 50 # the frequency of the reference/control waveform
ma = 0.8 # fr_p/fs_p the modulation index

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

Ts = 1/fs # period of the triangle waveform
ms = 4*fs_p/Ts # the gradient of the triangle waveform

# During one switching period, each switch switches on once and off once. Thus
# fs is also called the switching frequency

fr_p = ma*fs_p # the amplitude of the reference

mf = fs/fr # the frequency modulation index

# initialise

t = 0:1/Nps:(T/fr - 1/Nps) # time

N = length(t) # number of samples

Rf = Array{Float64, 1}(undef, N)
Tri = Array{Float64, 1}(undef, N)
Control = Array{Float64, 1}(undef, N)
Vo = Array{Float64, 1}(undef, N)

DFTtest = Array{Float64, 1}(undef, N)

Tri = fill!(Tri, 0.0)

Control = fill!(Control, -1)

toggl = 1

for i in 1:N

    global Rf, toggl
    Rf[i] = fr_p*sin(2*pi*fr*t[i])

    if Tri[i] >= fs_p && toggl == 1
        toggl = -1
    elseif Tri[i] <= -fs_p && toggl == -1
        toggl = 1
    end

    if i < N
        Tri[i+1] = toggl*ms/Nps + Tri[i]
    end

    if Rf[i] >= Tri[i]
        Control[i] = 1
    end

end

Vo = Control*Vd/2

p1 = plot(t, Rf)
#p1 = plot(t, DFTtest)
p1 = plot!(t, Tri)
p1 = plot!(t, Control)

display(p1)

W = Array{ComplexF64, 2}(undef, N, N)
W[:,1] = fill!(W[:,1], 1)
W[1,:] = fill!(W[1,:], 1)

for i in 2:N
    for j in i:N

        W[i,j] = exp(-1im*2*pi*(j-1)*(i-1)/N)

        if i != j
            W[j,i] = W[i,j]
        end

    end
end

F = Array{ComplexF64, 1}(undef, N)
F_abs = Array{ComplexF64, 1}(undef, convert(Int, round(N/2 + 0.51)))
F_angle = Array{ComplexF64, 1}(undef, convert(Int, round(N/2 + 0.51)))
F_real = Array{ComplexF64, 1}(undef, convert(Int, round(N/2 + 0.51)))
F_imag = Array{ComplexF64, 1}(undef, convert(Int, round(N/2 + 0.51)))

F = W*Vo/N

for i in 1:convert(Int, round(N/2 + 0.51))

    if i > 1 && i != round(N/2 + 0.51)
        F[i] = 2*F[i]
    end

    F_abs[i] = abs(F[i])
    F_angle[i] = angle(F[i])
    F_real[i] = real(F[i])
    F_imag[i] = imag(F[i])

end

count = 2

freq = Array{Float64, 1}(undef, 126)
harm = Array{Float64, 1}(undef, 126)

harm[1] = ma*Vd/2
freq[1] = 50

# analytic evaluation of harmonics
for m in 1:1:4
    for n in -15:1:15

        global count

        freq[count] = m*fs + n*fr
        harm[count] = real(abs((Vd/(1im*m*π))*besselj(n,m*π*ma/2)*exp(1im*n*π)*(1 - (-1)^(n)*exp(1im*m*π))))

        count = count + 1
    end
end

f_four = 0:fr/T:round(N/2)*fr/T

DFTs_p = plot(f_four, real(F_abs), seriestype = :scatter, line=:stem,
markershape = :circle, markersize = 2, title = "Frequency Spectrum", label="Shifted DFT")

original_k = 1:N
shifted_k = fftshift(fftfreq(N, Nps))/1000
fft_original = 2/N*fft(Vo)
shifted_fft = 2/N*fftshift(fft(Vo))

FFT_p = plot(original_k, abs.(fft_original), seriestype = :scatter, line=:stem,
title = "Original FFT",
markershape = :circle, markersize = 2,
#xticks=original_k[1:500:end],
legend=false);

FFTs_p = plot(shifted_k, abs.(shifted_fft), seriestype = :scatter, line=:stem,
title = "Frequency Spectrum",
label = "Shifted FFT",
markershape = :circle, markersize = 2,
#xticks=shifted_k[1:50:end],
legend=true);

p_out = plot(DFTs_p, FFTs_p, layout = (2,1))
plot!(subplot = 1, freq, harm, lw = 1,
line=:stem, seriestype = :scatter, label = "analytic",
markershape = :star, markersize = 2,
linecolor = :red)

display(p_out)

println("\nThe Fundamental is at a frequency of ", f_four[T + 1],
" and has an amplitude of ", round(real(F_abs[T + 1]), digits = 3))
println("\nThe amplitude of the Fundamental is approximately,
ma x Vd/2 =, ", ma * Vd/2)
println("\nThe harmonics appear in groups around integer
multiples of the switching frequency")
println("\nThe magnitude of the harmonic at the switching frequency, ",
 f_four[Int(mf*T) + 1], " is, ",
round(abs(fft_original[Int(mf*T) + 1]), digits = 3))

println("\nWe usually design the output filter in such a way that
it cuts off between the fundamental and switching frequency")


#=
    It is important to calculate the losses in switching components to
    ensure that the maximum junction temperature is not exceeded
    We differentiate between switching losses and conduction losses.
    The conduction losses of Mosfets are calculated in a different
    way than that of an IGBT.

    The conduction losses in the switch are proportional to the
    electrical potential over the switch when in is "on", and the average
    current is flowing throuhg it.

    Pcond = Von*Id

    The switching losses in the switch.

    The switching times of the switch are very short in comparison with
    the total switching period Ts. We can therefore assume that the
    inductor current is constant during switch-on with a value of iL(min),
    or iL(max).

    We use the minimum value because the switch is switching on, and
    maximum when switching off.
    Before the switch switches on the diode conducts the inductor current.
=#


#=
    To calculate the Fourier series expansion analytically,
    through a simulation one needs to use very small time steps.
    This can become impractical.
    A detailed study of the harmonics is important in order to
    understand the working of some types of converters. Examples
    are power electronic amplifiers and multilevel converters.

    It is difficult to calculate the Fourier coefficients.
    To do this one needs to calculate the cross-over points
    between the triangle wave and sinusoidal reference.

    When the frequency modulation index mf is not an integer
    , Vo, is not necessarily a periodic function. In this case
    the standard Fourier techniques are useless.

    The solution is to approach the problem from another perspective.
=#

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
