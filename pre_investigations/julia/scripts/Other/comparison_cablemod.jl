using DrWatson
@quickactivate "MicroGridSimWithRL"
include(srcdir("nodeconstructor.jl"))

using ControlSystems
using JSON
using Plots
using LinearAlgebra

discrete = false
julia = true
cut_outliers = false
num_cm = 1
num_mat_start = 1 
num_mat_end = 30

ts=1e-8

# only source to load connections:

CM = [ 0.  0.  2.
       0.  0.  1.
      -2. -1.  0.]
      parameter = Dict{Any, Any}("source" => Any[Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LC")], "load" => Any[Dict{Any, Any}("R" => 14, "impedance" => "R")], "cable" => Any[Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722)])


# fully connected
# CM = [ 0.  1.  3.
#       -1.  0.  2.
#       -3. -2.  0.]
# parameter = Dict{Any, Any}("source" => Any[Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "L2" => 0.0023, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LCL", "R2" => 0.4), Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LC")], "load" => Any[Dict{Any, Any}("R" => 14, "impedance" => "R")], "cable" => Any[Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722), Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722)])


Grid_FC = NodeConstructor(num_source=2, num_loads=1, CM=CM ,parameters=parameter )
A, B, C, D = get_sys(Grid_FC)
Ad = exp(A*ts)
Bd = A \ (Ad - C) * B

println(imag(eigvals(A)))
println(imag(eigvals(A))*(2*pi)^(-1))

global sys_d = StateSpace(Ad, Bd, C, D, ts)

t = collect(0:ts:0.01)

ns = length(A[1,:])
na = length(B[1,:])

#x0 = real(eigenvec[:,2])

#global x0 = [0.0 for i = 1:length(t)]

global x0 = [0.0 for i = 1:ns]
#global u = rand(Float64, ( length(t) )) .*2 .-1
global u = [250.0 for i = 1:length(t)]
global uu = [u for i = 1:na ]
global uuu = mapreduce(permutedims, vcat, uu)
global ttt = t

#xout = lsim(sys_d, uuu, ttt, x0)
xout, _, _, _ = lsim(sys_d,uuu,ttt,x0=x0)

#p4 = plot(t, xout[2,:], xlabel="time", ylabel="v_LC / V")
#display(plot(p4))
p4 = plot(t, xout[3,:], xlims = (0,0.01), ylims = (0,15), title = "cable_mod", xlabel="time", ylabel="i_LCL / A")
##p4 = plot(t, xout[3,:], xlabel="time", ylabel="i_LCL / A")
#display(plot(p4))

#p5 = plot(t, xout[ns,:], xlabel="time", ylabel="v_cable / V") # -> R_loads
#display(plot(p5))
p5 = plot(t, xout[6,:], xlims = (0,0.005), ylims = (0,350), title = "cable_mod", xlabel="time", ylabel="v_LC / V")
##p5 = plot(t, xout[6,:], xlabel="time", ylabel="v_LC / V") # -> RLC_loads
#display(plot(p5))




## without cablemod
## only source to load connections:

R_11= parameter["source"][1]["R1"]
R_22= parameter["source"][1]["R2"]
R_F1= parameter["source"][1]["R_C"]
L_11= parameter["source"][1]["L1"]
L_12= parameter["source"][1]["L2"]
C_F1= parameter["source"][1]["C"]

L_21= parameter["source"][2]["L1"]
R_21= parameter["source"][2]["R1"]
R_F2= parameter["source"][2]["R_C"]
C_F2= parameter["source"][2]["C"]

R_b1= parameter["cable"][1]["R"]
R_b2= parameter["cable"][2]["R"]
#R_b3= parameter["cable"][3]["R"]

R_L= parameter["load"][1]["R"]

## fully connected:

# R_11= parameter["source"][1]["R1"]
# R_22= parameter["source"][1]["R2"]
# R_F1= parameter["source"][1]["R_C"]
# L_11= parameter["source"][1]["L1"]
# L_12= parameter["source"][1]["L2"]
# C_F1= parameter["source"][1]["C"]

# L_21= parameter["source"][2]["L1"]
# R_21= parameter["source"][2]["R1"]
# R_F2= parameter["source"][2]["R_C"]
# C_F2= parameter["source"][2]["C"]

# R_b1= parameter["cable"][1]["R"]*(2/3)
# R_b2= parameter["cable"][2]["R"]*(2/3)
# R_b3= parameter["cable"][3]["R"]*(2/3)


# R_L= parameter["load"][1]["R"] + R_b3

A_no = zeros(5,5)

helpvar= 1 + (R_L+R_b2)*(R_F2)^(-1)
A_no[1,1]= -(R_11+R_F1)*(L_11)^(-1)
A_no[1,2]= - (L_11)^(-1)
A_no[1,3]= (R_F1)*(L_11)^(-1)
A_no[2,1]= (C_F1)^(-1)
A_no[2,3]= -(C_F1)^(-1)
A_no[3,1]= (R_F1)*(L_12)^(-1)
A_no[3,2]= (L_12)^(-1)
A_no[3,3]= -(R_F1)*(L_12)^(-1)-(R_22+R_b1)*(L_12)^(-1)-R_L*(L_12)^(-1)+(R_L*R_L)*(helpvar*R_F2)^(-1)*(L_12)^(-1)
A_no[3,4]= -R_L*(L_12)^(-1) + (R_L)*(R_L+R_b2)*(R_F2*helpvar)^(-1)*(L_12)^(-1)
A_no[3,5]= -(R_L)*(helpvar*R_F2)^(-1)*(L_12)^(-1)
A_no[4,3]= -(R_L)*(helpvar*L_21)^(-1)
A_no[4,4]= (-R_21 -(R_L+R_b2)*(helpvar)^(-1))*(L_21)^(-1)
A_no[4,5]= (-1+ helpvar^(-1))*(L_21)^(-1)
A_no[5,3]= (R_L*(R_F2*helpvar)^(-1))*(C_F2)^(-1)
A_no[5,4]= ((R_L+R_b2)*(R_F2*helpvar)^(-1))*(C_F2)^(-1)
A_no[5,5]= -(R_F2*helpvar*C_F2)^(-1)


println(imag(eigvals(A_no)))
println(imag(eigvals(A_no))*(2*pi)^(-1))

B_no= zeros(5,2)
B_no[1,1]= (L_11)^(-1)
B_no[4,2]= (L_21)^(-1)

C_no = Diagonal(ones(5))
D_no = 0

Ad_no = exp(A_no*ts)
Bd_no = A_no \ (Ad_no - C_no) * B_no

global sys_d_no = StateSpace(Ad_no, Bd_no, C_no, D_no, ts)

t_no = collect(0:ts:0.1)

ns_no = length(A_no[1,:])
na_no = length(B_no[1,:])

#x0 = real(eigenvec[:,2])

#global x0 = [0.0 for i = 1:length(t)]

global x0_no = [0.0 for i = 1:ns_no]
#global u = rand(Float64, ( length(t) )) .*2 .-1
global u_no = [250.0 for i = 1:length(t_no)]
global uu_no = [u_no for i = 1:na_no ]
global uuu_no = mapreduce(permutedims, vcat, uu_no)
global ttt_no = t_no

#xout = lsim(sys_d, uuu, ttt, x0)


xout_no, _, _, _ = lsim(sys_d_no,uuu_no,ttt_no,x0=x0_no)

#p4 = plot(t, xout[2,:], xlabel="time", ylabel="v_LC / V")


p7 = plot(t_no, xout_no[5,:], xlims = (0,0.005), ylims = (0,350), title = "no_cable_mod", xlabel="time", ylabel="v_LC / V")
##p7 = plot(t_no, xout_no[5,:], xlabel="time", ylabel="v_LC / V")
#display(plot(p6))

#p5 = plot(t, xout[ns,:], xlabel="time", ylabel="v_cable / V") # -> R_loads
#display(plot(p5))

p6 = plot(t_no, xout_no[3,:], xlims = (0,0.01), ylims = (0,15), title = "no_cable_mod", xlabel="time", ylabel="i_LCL / A")
##p6 = plot(t_no, xout_no[3,:], xlabel="time", ylabel="i_LCL / A") # 
#display(plot(p7))


p20= plot()
p20= plot!(t, xout[3,:], xlims = (0,0.01), ylims = (0,15), title = "LCL_current", xlabel="time", ylabel="i_LCL / A", label="cable" )
p20= plot!(t_no, xout_no[3,:], xlims = (0,0.01), ylims = (0,15), title = "LCL_current", xlabel="time", ylabel="i_LCL / A", label="no_cable" )
display(p20)

p8= plot()
p8 = plot!(t, xout[3,:], xlims = (0.003,0.007), ylims = (6.5,9), title = "LCL_current", xlabel="time", ylabel="i_LCL / A", label="cable")
p8 = plot!(t_no, xout_no[3,:], xlims = (0.003,0.007), ylims = (6.5,9), title = "LCL_current", xlabel="time", ylabel="i_LCL / A", label="no_cable")
display(p8)

p21= plot()
p21= plot!(t, xout[6,:], xlims = (0,0.005), ylims = (0,350) , title = "LC_voltage",xlabel="time", ylabel="v_LC / V", label="cable" )
p21= plot!(t_no, xout_no[5,:], xlims = (0,0.005), ylims = (0,350), title = "LC_voltage", xlabel="time", ylabel="v_LC / V", label="no_cable")
display(p21)

p12 = plot()
p12 = plot!(t, xout[6,:], xlims = (0,0.0001), ylims = (0,160), title = "LC_voltage", xlabel="time", ylabel="v_LC / V", label="cable" )
p12 = plot!(t_no, xout_no[5,:], xlims = (0,0.0001), ylims = (0,160), title = "LC_voltage", xlabel="time", ylabel="v_LC / V", label="no_cable" )
display(p12)

p10 = plot()
p10 = plot!(t, xout[6,:], xlims = (0.001,0.004), ylims = (230,270), title = "LC_voltage", xlabel="time", ylabel="v_LC / V", label="cable" )
p10 = plot!(t_no, xout_no[5,:], xlims = (0.001,0.004), ylims = (230,270), title = "LC_voltage", xlabel="time", ylabel="v_LC / V", label="no_cable" )
display(p10)

