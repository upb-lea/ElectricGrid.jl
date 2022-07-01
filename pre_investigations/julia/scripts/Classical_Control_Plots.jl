function Plot_I_dq0(num_source, T_plot_start, T_plot_end, Source::Source_Controller)

    T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
    tsys = 1/fsys

    N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
    N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
    N_range = N_plot_start:N_plot_end

    p_I_dq0_d = plot(t[N_range], Source.I_ref_dq0[num_source, 1, N_range],
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
    p_I_dq0_d = plot!(t[N_range], Source.I_dq0[num_source, 1, N_range],
        label = "d-axis")

    p_I_dq0_q = plot(t[N_range], Source.I_ref_dq0[num_source, 2, N_range],
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
    p_I_dq0_q = plot!(t[N_range], Source.I_dq0[num_source, 2, N_range],
            label = "q-axis")

    p_I_dq0 = plot(p_I_dq0_d, p_I_dq0_q,# p_I_dq0_0,
        layout = (2, 1),
        legend = true,
        size = (850,900))

    display(p_I_dq0)

    return nothing
end

function Plot_V_dq0(num_source, T_plot_start, T_plot_end, Source::Source_Controller)

    T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
    tsys = 1/fsys

    N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
    N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
    N_range = N_plot_start:N_plot_end


    p_V_dq0_d = plot(t[N_range], Source.V_ref_dq0[num_source, 1, N_range], label = "ref d-axis",
        legend = :bottomright,
        xlabel = "Time [s]",
        ylabel = "Voltage [V]",
        title = "DQ0 Voltage Transform\nDirect-Axis",
        grid = true,
        foreground_color_grid = :black,
        minorgrid = true,
        thickness_scaling = 1.5,
        legendfont = font(5))
    p_V_dq0_d = plot!(t[N_range], Source.V_dq0[num_source, 1, N_range], label = "d-axis")

    p_V_dq0_q = plot(t[N_range], Source.V_ref_dq0[num_source, 2, N_range], label = "ref q-axis",
            legend = :bottomright,
            xlabel = "Time [s]",
            ylabel = "Voltage [V]",
            title = "Quadrature-Axis",
            grid = true,
            foreground_color_grid = :black,
            minorgrid = true,
            thickness_scaling = 1.5,
            legendfont = font(5))
    p_V_dq0_q = plot!(t[N_range], Source.V_dq0[num_source, 2, N_range], label = "q-axis")

    p_V_dq0 = plot(p_V_dq0_d, p_V_dq0_q,# p_I_dq0_0,
        layout = (2, 1),
        legend = true,
        size = (850,900))

    display(p_V_dq0)

    return nothing
end

function Inst_Vout_Vref(num_source, num_node, T_plot_start, T_plot_end, Source::Source_Controller, Env::Environment)

    T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
    tsys = 1/fsys

    N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
    N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
    N_range = N_plot_start:N_plot_end

    # Phase a Control Signals
    p_cntr_a = plot(t[N_range], Source.V_ref[1, 1, N_range],
        label = "Reference_a",
        xlabel = "Time [s]",
        ylabel = "Voltage [V]",
        title = "DC-AC Converter Control Phase A")
    p_cntr_a = plot!(t[N_range], Env.x[2,N_range],
        label = "Inverter Phase a")

    # Phase b Control Signals
    p_cntr_b = plot(t[N_range], Source.V_ref[1, 2, N_range],
        label = "V_ref b",
        xlabel = "Time [s]",
        ylabel = "Voltage [V]",
        title = "DC-AC Converter Control Phase B")
    p_cntr_b = plot!(t[N_range], Env.x[5,N_range],
        label = "Inverter Phase b")

    # Phase c Control Signals
    p_cntr_c = plot(t[N_range], Source.V_ref[1, 3, N_range],
        label = "Reference_c",
        xlabel = "Time [s]",
        ylabel = "Voltage [V]",
        title = "DC-AC Converter Control Phase C")
    p_cntr_c = plot!(t[N_range], Env.x[8,N_range],
        label = "Inverter Phase c")

    p_v_cntr = plot(p_cntr_a, p_cntr_b, p_cntr_c,
        layout = (3, 1),
        legend = true,
        size = (850,900))

    display(p_v_cntr)

    return nothing
end

function Plot_Irms(num_node, T_plot_start, T_plot_end, Env::Environment)
    T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
    tsys = 1/fsys

    N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
    N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
    N_range = N_plot_start:N_plot_end

    p_i_rms = plot(t[N_range], Env.I_ph[num_node,1,2,N_range],
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
    p_i_rms = plot!(t[N_range], Env.I_ph[num_node,2,2,N_range],
        label = "Inverter Phase b")
    p_i_rms = plot!(t[N_range], Env.I_ph[num_node,3,2,N_range],
        label = "Inverter Phase c")

    p_i_ang = plot(t[N_range], (180/π)*Env.I_ph[num_node,1,3,N_range],
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
    p_i_ang = plot!(t[N_range], (180/π)*Env.I_ph[num_node,2,3,N_range],
        label = "Inverter Phase b")
    p_i_ang = plot!(t[N_range], (180/π)*Env.I_ph[num_node,3,3,N_range],
        label = "Inverter Phase c")

    p_i_rms_ang = plot(p_i_rms, p_i_ang,
        layout = (2, 1),
        legend = true,
        size = (900,900))

    display(p_i_rms_ang)

    return nothing
end

function Plot_Vrms(num_node, num_source, T_plot_start, T_plot_end, Env::Environment, Source::Source_Controller)

    T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
    tsys = 1/fsys

    N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
    N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))

    T_sp_rms = 5*fsys #samples in a second
    i_sp_rms = convert(Int64, round((1/(μps*T_sp_rms))))
    N_range = N_plot_start:N_plot_end

    p_v_rms = plot(t[N_range], Env.V_ph[num_node, 1,2,N_range],
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
    p_v_rms = plot!(t[N_range], Env.V_ph[num_node, 2,2,N_range],
        label = "Inverter Phase b")
    p_v_rms = plot!(t[N_range], Env.V_ph[num_node, 3,2,N_range],
        label = "Inverter Phase c")

    p_v_ang = plot(t[N_range], (180/π)*Env.V_ph[num_node, 1,3,N_range],
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
    p_v_ang = plot!(t[N_range], (180/π)*Env.V_ph[num_node, 2,3,N_range],
        label = "Inverter Phase b")
    p_v_ang = plot!(t[N_range], (180/π)*Env.V_ph[num_node, 3,3,N_range],
        label = "Inverter Phase c")

    PLL_ph_a = Source_Control.θpll[num_source, 1, N_range] .- 2*π*Source.fpll[num_source, 1, N_range].*t[N_range]
    PLL_ph_b = Source_Control.θpll[num_source, 2, N_range] .- 2*π*Source.fpll[num_source, 2, N_range].*t[N_range]
    PLL_ph_c = Source_Control.θpll[num_source, 3, N_range] .- 2*π*Source.fpll[num_source, 3, N_range].*t[N_range]

    for i in 1:length(PLL_ph_a[N_range])
        PLL_ph_a[i] = (PLL_ph_a[i] + 2*π*floor(t[i]/0.02))*180/pi
        PLL_ph_b[i] = (PLL_ph_b[i] + 2*π*floor(t[i]/0.02))*180/pi
        PLL_ph_c[i] = (PLL_ph_c[i] + 2*π*floor(t[i]/0.02))*180/pi
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

    p_v_ang = plot!(t[N_range], PLL_ph_a[N_range],
        label = "PLL Phase a",
        linestyle = :dash,
        linewidth = 2)
    p_v_ang = plot!(t[N_range], PLL_ph_b[N_range],
        label = "PLL Phase b",
        linestyle = :dash,
        linewidth = 2)
    p_v_ang = plot!(t[N_range], PLL_ph_c[N_range],
        label = "PLL Phase c",
        linestyle = :dash,
        linewidth = 2)

    p_v_rms_ang = plot(p_v_rms, p_v_ang,
        layout = (2, 1),
        legend = :bottomright,
        size = (900,900))

    display(p_v_rms_ang)

    return nothing
end

function Inst_Iout_Vref(num_source, num_node, T_plot_start, T_plot_end, Source::Source_Controller, Env::Environment)

    T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
    tsys = 1/fsys

    N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
    N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))

    N_range = N_plot_start:N_plot_end

    # Phase a Control Signals
    p_cntr_a = plot(t[N_range], Source.I_ref[num_source, 1, N_range],
        label = "Reference_a",
        xlabel = "Time [s]",
        ylabel = "Current [A]",
        title = "DC-AC Converter Control Phase A")
    p_cntr_a = plot!(t[N_range], Env.x[3, N_range],
        label = "Inverter Phase a")

    # Phase b Control Signals
    p_cntr_b = plot(t[N_range], Source.I_ref[num_source, 2, N_range],
        label = "Reference_b",
        xlabel = "Time [s]",
        ylabel = "Current [A]",
        title = "DC-AC Converter Control Phase B")
    p_cntr_b = plot!(t[N_range], Env.x[6, N_range],
        label = "Inverter Phase b")

    # Phase c Control Signals
    p_cntr_c = plot(t[N_range], Source.I_ref[num_source, 3, N_range],
        label = "Reference_c",
        xlabel = "Time [s]",
        ylabel = "Current [A]",
        title = "DC-AC Converter Control Phase C")
    p_cntr_c = plot!(t[N_range], Env.x[9, N_range],
        label = "Inverter Phase c")

    p_i_cntr = plot(p_cntr_a, p_cntr_b, p_cntr_c,
        layout = (3, 1),
        legend = true,
        size = (850,900))

    display(p_i_cntr)

    return nothing
end

function Plot_PLL(num_source, T_plot_start, T_plot_end, Source::Source_Controller, Env::Environment)

    T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
    tsys = 1/fsys

    N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
    N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
    N_range = N_plot_start:N_plot_end

    p_pll_f = plot(t[N_range], Source.fpll[1,1,N_range], label = "PLL Frequency",
        legend = :bottomright,
        xlabel = "Time [s]",
        ylabel = "Frequency [Hz]",
        title = "Phase-Locked-Loop",
        grid = true,
        foreground_color_grid = :black,
        minorgrid = true,
        thickness_scaling = 1.5,
        legendfont = font(5))
    p_pll_f = plot!(t[N_range], Env.fs[N_range], label = "System Frequency")

    θs = (2*π*Env.fs.*t).%(2*π)

    N_range = N_plot_start:N_plot_end
    θe = sin.(θs[N_range] .- Source_Control.θpll[1, 1, N_range])
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

    return nothing
end

function Plot_P_inst(num_node, T_plot_start, T_plot_end, Env::Environment)

    T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
    tsys = 1/fsys

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

    display(p_p)

    return nothing
end

function Plot_Real_Imag_Active_Reactive(num_node, num_source, T_plot_start, T_plot_end, Env::Environment, Source::Source_Controller)

    T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
    tsys = 1/fsys

    N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
    N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
    N_range = N_plot_start:N_plot_end

    Uc = 0.001 # unit conversion

    p_p_r_a = plot(t[N_range], Uc*Env.p_q_inst[1, 1, N_range],
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
    #p_p_r_a = plot!(t[N_range], Uc*Source.p_q_filt[1, 1, N_range], label = "Filtered Real Power")
    p_p_r_a = plot!(t[N_range], Uc*Env.P[1, 4, N_range], label = "Active Power", lw = 1)

    p_p_i_q = plot(t[N_range], Uc*Env.p_q_inst[1, 2 ,N_range], label = "Imaginary Power",
        legend = :bottomright,
        xlabel = "Time [s]",
        ylabel = "Power [kVAi / kVAr]",
        title = "Inverter Imaginary and Reactive Power",
        grid = true,
        foreground_color_grid = :black,
        minorgrid = true,
        thickness_scaling = 1.5,
        legendfont = font(5))
    #p_p_i_q = plot!(t[N_range], Uc*Source.p_q_filt[1, 2, N_range], label = "Filtered Imaginary Power")
    p_p_i_q = plot!(t[N_range], Uc*Env.Q[1, 4, N_range], label = "Reactive Power", lw = 1)

    p_p_real_imag_act_react = plot(p_p_r_a, p_p_i_q,
        layout = (2, 1),
        legend = true,
        size = (900,900))

    display(p_p_real_imag_act_react)

    return nothing
end

function Plot_fft(num_node, num_source, T_plot_start, T_plot_end, Env::Environment, Source::Source_Controller)

    T_s = t_final*fsys # periods for the reference stage, i.e. simulation time
    tsys = 1/fsys

    x_lim = (-1500, +1500)
    x_ticks = -1500:250:1500

    N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
    N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
    N_range = N_plot_start:N_plot_end
    N_length = length(N_range)
    freqs = fftshift(fftfreq(N_length, 1/μps))

    Uc = 0.001 # unit conversion

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

    fft_p_real = (2/N_length)*fftshift(fft(Uc*Env.p_q_inst[1, 1, N_range]))
    fft_p_real_filt = (2/N_length)*fftshift(fft(Uc*Source.p_q_filt[1, 1, N_range]))
    FFTs_p_real = plot(freqs, abs.(fft_p_real), label = "Shifted FFT",
        seriestype = :scatter, line = :stem,
        title = "Inverter POC Real Power",
        xlim = x_lim,
        xticks = x_ticks,
        xlabel = "Frequency [Hz]",
        ylabel = "Real Power [W]",
        markershape = :circle, markersize = 2,
        legend = true);
    FFTs_p_real = plot!(freqs, abs.(fft_p_real_filt), label = "Filtered FFT")

    fft_p_imag = (2/N_length)*fftshift(fft(Uc*Env.p_q_inst[1, 2, N_range]))
    fft_p_imag_filt = (2/N_length)*fftshift(fft(Uc*Source.p_q_filt[1, 2, N_range]))
    FFTs_p_imag = plot(freqs, abs.(fft_p_imag), label = "Shifted FFT",
        seriestype = :scatter, line = :stem,
        title = "Inverter POC Imaginary Power",
        xlim = x_lim,
        xticks = x_ticks,
        xlabel = "Frequency [Hz]",
        ylabel = "Imaginary Power [VAi]",
        markershape = :circle, markersize = 2,
        legend = false);
    FFTs_p_imag = plot!(freqs, abs.(fft_p_imag_filt), label = "Filtered FFT")

    fft_v_poc_a = (2/N_length)*fftshift(fft(Env.x[2,N_range]))
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

    fft_i_poc_a = (2/N_length)*fftshift(fft(Env.x[3,N_range]))
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

    display(p_fft)

    return nothing
end
