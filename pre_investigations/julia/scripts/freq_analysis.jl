using DrWatson
@quickactivate "MicroGridSimWithRL"

using Plots

num_cable_parallel = [2]#[300]#[1, 5, 10, 20, 50]
cable_len = [1,2,3,4,5,6,7,8]#[40, 45, 50, 55, 60, ] #[1, 5, 10, 50]

arr = []
ts = 1e-4
t_end = 0.05
t = collect(0:ts:t_end)
t_range = 0:ts:t_end

Rf = 0.4
Lf = 2.3e-3
Cf = 10e-6

f_0_filter_cable_LCL = []
f_0_cable_LCL = []
f_0_cable_LC = []

result = Dict()

f_0_filter_cable_LCL_ndarray = zeros(length(cable_len), length(num_cable_parallel))
f_0_cable_LCL_ndarray = zeros(length(cable_len), length(num_cable_parallel))
f_0_cable_LC_ndarray = zeros(length(cable_len), length(num_cable_parallel))

for l in 1:length(cable_len)

    for n in 1:length(num_cable_parallel)

        #Rb = l * 0.722
        #Cb = l * 8e-9 
        #Lb = l * 0.955e-3

        Rb = cable_len[l] * 0.722
        Cb = cable_len[l] * 8e-9 
        Lb = cable_len[l] * 0.955e-3

        #append!(arr, l)


        
        # Series connection of Cf and n_cable*Cb
        C_LCL_fb = 1/Cf
        Cb_sum = 0
        C_bb_LC_sum = 0
        for i in 1:num_cable_parallel[n]#n_cable
            #C_LCL_fb = C_LCL_fb + 1/Cb
            #C_LCL_bb = C_LCL_bb + 1/Cb
            Cb_sum = Cb_sum + 1/Cb
            C_bb_LC_sum = C_bb_LC_sum + 1/(Cb+Cf)
        end
        C_fb1 = 1/(1/Cf + 1/Cb +1/Cb)
        C_LCL_fb = 1/(1/Cf + Cb_sum)
        C_LCL_fb1 = 1/(1/Cf + 1/Cb_sum)
        

        # Bei LCL Filter
        #C_LCL_bb1 = 1/(1/Cb + 1/Cb)
        C_LCL_bb = 1/Cb_sum

        #Bei LC Filter
        #C_bb_LC1 = 1/(1/(Cb+Cf) + 1/(Cb+Cf)) 
        C_bb_LC = 1/C_bb_LC_sum

        append!(f_0_filter_cable_LCL, 1/sqrt(Lf*C_LCL_fb)/2/pi)
        append!(f_0_cable_LCL, 1/sqrt(Lb*C_LCL_bb)/2/pi)
        append!(f_0_cable_LC, 1/sqrt(Lb*C_bb_LC)/2/pi)

        f_0_filter_cable_LCL_ndarray[l, n] = 1/sqrt(Lf*C_LCL_fb)/2/pi
        f_0_cable_LCL_ndarray[l, n] = 1/sqrt(Lb*C_LCL_bb)/2/pi
        f_0_cable_LC_ndarray[l, n] = 1/sqrt(Lb*C_bb_LC)/2/pi

        println(" ")
        println("Anzahl an parallel geschlossenen Kabeln: $(num_cable_parallel[n])")

        println("Länge Kabel: $(cable_len[l]) km")

        println("Resonanz Freq Filter-Kabel mit LCL Filter: $(f_0_filter_cable_LCL) Hz")
        println("Resonanz Freq Kabel mit LCL Filter: $(f_0_cable_LCL) Hz")
        println("Resonanz Freq Kabel mit LC Filter: $(f_0_cable_LC) Hz")

        println(" ")

    end
end

result["cable_length"] = cable_len
result["num_cable_parallel"] = num_cable_parallel
result["f_0_filter_cable_LCL"] = f_0_filter_cable_LCL_ndarray
result["f_0_cable_LCL"] = f_0_cable_LCL_ndarray
result["f_0_cable_LC"] = f_0_cable_LC_ndarray

#=
for ll in 1:length(cable_len)
    p1 = plot(num_cable_parallel, f_0_filter_cable_LCL_ndarray[ll, :], ylabel="f0_fc_LCL", title="Cable_len: $(cable_len[ll])")
    p2 = plot(num_cable_parallel, f_0_cable_LCL_ndarray[ll, :], ylabel="f0_c_LCL")
    p3 = plot(num_cable_parallel, f_0_cable_LC_ndarray[ll, :], xlabel="num_cable_parallel", ylabel="f0_c_LC", title="Cable_len: $(cable_len[ll]) km", label="f0_c_LC")
    
    #display(p3)
    display(plot(p1, p2, layout=(2,1)))
end
=#
for nn in 1:length(num_cable_parallel)
    p1 = plot(cable_len, f_0_filter_cable_LCL_ndarray[:, nn], ylabel="f0_fc_LCL", title="Num_cable_parallel: $(num_cable_parallel[nn])")
    p2 = plot(cable_len, f_0_cable_LCL_ndarray[:, nn], ylabel="f0_c_LCL")
    p3 = plot(cable_len, f_0_cable_LC_ndarray[:, nn], xlabel="cable_len", ylabel="f0_c_LC")
    #display(p3)
    display(plot(p1, p2, layout=(2,1)))
end
#= 

#Plots.default(overwrite_figure=false)
nn = 1
p1 = plot(cable_len, f_0_filter_cable_LCL_ndarray[:, nn], ylabel="f0_fc_LCL", title="Num_cable_parallel: $(num_cable_parallel[nn])", ylim=(0.5e4, 2e4))
p2 = plot(cable_len, f_0_cable_LCL_ndarray[:, nn], ylabel="f0_c_LCL", ylim=(0.5e4, 2e4))
#p3 = plot(cable_len, f_0_cable_LC_ndarray[:, nn], xlabel="cable_len", ylabel="f0_c_LC")

display(plot(p1, p2))


arr = zeros(length(cable_len), 2)
arr[:,1] = f_0_filter_cable_LCL_ndarray[:, nn]
#arr[:,2] = f_0_cable_LCL_ndarray[:, nn]

p3 = plot(cable_len, arr, xlabel="cable_len", ylabel="Freq / Hz" , title="Num_cable_parallel: $(num_cable_parallel[nn])", label=["f0_fc_LCL" "f0_c_LCL"])
display(p3)
=#