using ControlSystems
using JSON
using Plots
using LinearAlgebra
using DrWatson
@quickactivate "MicroGridSimWithRL"

include(srcdir("nodeconstructor.jl"))

ts = 1e-8

n = 3

CM = [0.0   1.0    2.0    3.0    4.0   5.0
-1.0   0.0    6.0    7.0    8.0   9.0
-2.0  -6.0    0.0   10.0   11.0  12.0
-3.0  -7.0  -10.0    0.0   13.0  14.0
-4.0  -8.0  -11.0  -13.0    0.0  15.0
-5.0  -9.0  -12.0  -14.0  -15.0   0.0];

para_1LC_1000m = Dict{Any, Any}("source" => Any[Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LC")], "load" => Any[Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R")], "cable" => Any[Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722)])

para_1LC_10m = Dict{Any, Any}("source" => Any[Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LC")], "load" => Any[Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R")], "cable" => Any[Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722)])

para_1LC_100m = Dict{Any, Any}("source" => Any[Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LC")], "load" => Any[Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R")], "cable" => Any[Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722)])

list = [para_1LC_10m, para_1LC_100m, para_1LC_1000m]
labels = ["para_1LC_10m", "para_1LC_100m", "para_1LC_1000m"]
global p1 = plot()
global p2 = plot()
global p3 = plot()
global p4 = plot()
global p5 = plot()
global p6 = plot()
global p7 = plot()
global p8 = plot()
global p9 = plot()
global p10 = plot()
global p11 = plot()
global p12 = plot()

for (idx, item) in enumerate(list)
    nc = NodeConstructor(num_source=n, num_loads=n, CM=CM, parameters=item)

    if idx==1
        print("States: $(get_states(nc))")
    end

    A, B, C, D = get_sys(nc)
    ns = length(A[1,:])
    na = length(B[1,:])
    ####### plot step response
    # use therefore the last (!) model
    Ad = exp(A*ts)
    Bd = A \ (Ad - C) * B
    sys_d = StateSpace(Ad, Bd, C, D, ts)

    t = collect(0:ts:0.005)

    x0 = [0.0 for i = 1:ns]
    u = [250.0 for i = 1:length(t)]
    uu = [u for i = 1:na ]
    uuu = mapreduce(permutedims, vcat, uu)
    ttt = t
    xout, _, _, _ = lsim(sys_d,uuu,ttt,x0=x0)

    label = labels[idx]
    helper = label[10:end]

    states = get_states(nc)

    
    global p1 = plot!(p1, t, xout[10,:], title = "Ts: $ts", xlabel="time", ylabel=states[10], label = "$helper")
    global p2 = plot!(p2, t, xout[end,:], xlabel="time", ylabel=states[end], label = "$helper")
    global p11 = plot!(p11, t, xout[15,:], xlabel="time", ylabel=states[15], label = "$helper")

    global p5 = plot!(p5, t, xout[10,:], xlims = (0.0005,0.001), ylims = (140,160), title = "Ts: $ts", xlabel="time", ylabel=states[10], label = "$helper")
    global p6 = plot!(p6, t, xout[end,:], xlims = (0.0005,0.001), ylims = (135,155), title = "Ts: $ts", xlabel="time", ylabel=states[end], label = "$helper")
    
    global p9 = plot!(p9, t, xout[end,:], xlims = (0,0.0002), ylims = (0,40), title = "Ts: $ts", xlabel="time", ylabel=states[end], label = "$helper")
end
display(plot(p1, p2, p11, layout= (3,1)))

ts = 1e-4

n = 3

CM = [0.0   1.0    2.0    3.0    4.0   5.0
-1.0   0.0    6.0    7.0    8.0   9.0
-2.0  -6.0    0.0   10.0   11.0  12.0
-3.0  -7.0  -10.0    0.0   13.0  14.0
-4.0  -8.0  -11.0  -13.0    0.0  15.0
-5.0  -9.0  -12.0  -14.0  -15.0   0.0];

para_1LC_1000m = Dict{Any, Any}("source" => Any[Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LC")], "load" => Any[Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R")], "cable" => Any[Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722)])

para_1LC_10m = Dict{Any, Any}("source" => Any[Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LC")], "load" => Any[Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R")], "cable" => Any[Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722), Dict{Any, Any}("C" => 4.0e-9, "L" => 2.64e-6, "R" => 0.00722)])

para_1LC_100m = Dict{Any, Any}("source" => Any[Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LC")], "load" => Any[Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R"), Dict{Any, Any}("R" => 14, "impedance" => "R")], "cable" => Any[Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722), Dict{Any, Any}("C" => 4.0e-8, "L" => 2.6400000000000005e-5, "R" => 0.0722)])

list = [para_1LC_10m, para_1LC_100m, para_1LC_1000m]
labels = ["para_1LC_10m", "para_1LC_100m", "para_1LC_1000m"]
global p = plot()

for (idx, item) in enumerate(list)
    nc = NodeConstructor(num_source=n, num_loads=n, CM=CM, parameters=item)
    A, B, C, D = get_sys(nc)
    ns = length(A[1,:])
    na = length(B[1,:])
    ####### plot step response
    # use therefore the last (!) model
    Ad = exp(A*ts)
    Bd = A \ (Ad - C) * B
    sys_d = StateSpace(Ad, Bd, C, D, ts)

    t = collect(0:ts:0.005)

    x0 = [0.0 for i = 1:ns]
    u = [250.0 for i = 1:length(t)]
    uu = [u for i = 1:na ]
    uuu = mapreduce(permutedims, vcat, uu)
    ttt = t
    xout, _, _, _ = lsim(sys_d,uuu,ttt,x0=x0)

    label = labels[idx]
    helper = label[10:end]

    states = get_states(nc)


    global p3 = plot!(p3, t, xout[10,:], title = "Ts: $ts", xlabel="time", ylabel=states[10], label = "$helper")
    global p4 = plot!(p4, t, xout[end,:], xlabel="time", ylabel=states[end], label = "$helper")
    
    global p7 = plot!(p7, t, xout[10,:], xlims = (0.0005,0.001), ylims = (140,160), title = "Ts: $ts", xlabel="time", ylabel=states[10], label = "$helper")
    global p8 = plot!(p8, t, xout[end,:], xlims = (0.0005,0.001), ylims = (135,155), title = "Ts: $ts", xlabel="time", ylabel=states[end], label = "$helper")

    global p10 = plot!(p10, t, xout[end,:], xlims = (0,0.0002), ylims = (0,40), title = "Ts: $ts", xlabel="time", ylabel=states[end], label = "$helper")
    global p12 = plot!(p12, t, xout[15,:], xlabel="time", ylabel=states[15], label = "$helper")
end

display(plot(p3, p4, p12, layout= (3,1)))

display(plot(p5, p7, layout= (2,1)))
display(plot(p6, p8, layout= (2,1)))
display(plot(p9, p10, layout= (2,1)))
display(plot(p1, p2, p3, p4, layout= (2,2)))