mutable struct NodeConstructor
    num_connections
    num_sources
    num_loads
    num_fltr_LCL
    num_fltr_LC
    num_fltr_L
    num_loads_RLC
    num_loads_LC
    num_loads_RL
    num_loads_RC
    num_loads_L
    num_loads_C
    num_loads_R
    num_impedance
    num_fltr
    num_spp
    cntr
    tot_ele
    CM
    parameters
    S2S_p
    S2L_p
    L2L_p
    ts
    invoke_pfe
    verbosity
end

"""
    NodeConstructor(;
        num_sources,
        num_loads,
        CM = nothing,
        parameters = nothing,
        S2S_p = 0.1,
        S2L_p = 0.8,
        L2L_p = 0.3,
        ts = 10000,
        verbosity = 0
        )

Create a mutable struct NodeConstructor, which serves as a basis for the creation of an
energy grid: `num_sources` corresponse to the amount of sources and `num_loads` is the
amount of loads in the grid. `CM` is the connection matrix which indicates how the elements
in the grid are connected to each other. To specify the elements of the net in more detail,
values for the elements can be passed via `parameters`. If no connection matrix is entered,
it can be generated automatically. `S2S_p` is the probability that a source is connected to
another source and `S2L_p` is the probability that a source is connected to a load.
"""
function NodeConstructor(; num_sources, num_loads, CM=nothing, parameters=nothing,
    S2S_p=0, S2L_p=1, L2L_p=0, ts=10000, invoke_pfe=true, verbosity=0)
    tot_ele = num_sources + num_loads
    cntr = 0
    num_connections = 0

    if CM === nothing
        cntr, CM = GenerateCM(num_sources, num_loads, S2L_p, S2S_p, L2L_p)
        num_connections = cntr
    else
        # if size(CM)[1] != tot_ele
        #     throw("Expect the number of elements in the node to match the specified
        #         structure in the CM, but got $tot_ele and $(size(CM)[1])")
        # end

        @assert(size(CM)[1] == tot_ele,
            "Expect the number of elements in the grid to match the specified
            structure in the CM, but got $tot_ele and $(size(CM)[1]).")

        @assert(sum(CM) == 0,
            "Expect the checksum over the CM matrix to be 0,
            but got $(sum(CM)).")

        @assert(size(CM)[1] == size(CM)[2],
            "Expect the CM matrix to be a square matrix,
            but got $(size(CM)[1]) columns and $(size(CM)[2]) rows.")

        @assert(-transpose(CM .* LowerTriangular(ones(size(CM)))) ==
                CM .* UpperTriangular(ones(size(CM))),
            "Expect the CM matrix to be a antisymetric: transpose(CM) = -CM.")

        num_connections = Int(maximum(CM))
    end

    if parameters === nothing || isa(parameters, Dict)
        # Checks if all entries are given, if not, fills up with random values
        parameters = CheckParameters(
            parameters,
            num_sources,
            num_loads,
            num_connections,
            CM,
            ts,
            invoke_pfe,
            verbosity
        )

        @assert(length(keys(parameters)) == 5,
            "Expect parameters to have the four entries 'cable', 'load', 'grid', 'weather' and
            'source' but got $(keys(parameters))")

        @assert(length(keys(parameters["grid"])) == 9,
            "Expect parameters['grid'] to have the 9 entries but got
            $(length(keys(parameters["grid"])))")

        @assert(length(parameters["source"]) == num_sources,
            "Expect the number of sources to match the number of sources in the parameters,
            but got $num_sources and $(length(parameters["source"]))")

        @assert(length(parameters["load"]) == num_loads,
            "Expect the number of loads to match the number of loads in the parameters,
            but got $num_loads and $(length(parameters["load"]))")

        @assert(length(parameters["cable"]) == num_connections,
            "Expect the number of cables to match the number of cables in the parameters, but
            got $num_connections and $(length(parameters["cable"]))")

        num_fltr_LCL, num_fltr_LC, num_fltr_L = CountFilters(parameters["source"])
        loads = CountLoads(parameters["load"])
        num_loads_RLC = loads[1]
        num_loads_LC = loads[2]
        num_loads_RL = loads[3]
        num_loads_RC = loads[4]
        num_loads_L = loads[5]
        num_loads_C = loads[6]
        num_loads_R = loads[7]

        @assert(num_fltr_LCL + num_fltr_LC + num_fltr_L == num_sources,
            "Expect the number of sources to be identical to the sum of the filter types, but
            the number of sources is $num_sources and the sum of the filters is
            $(num_fltr_LCL + num_fltr_LC + num_fltr_L)")

        @assert((num_loads_RLC + num_loads_LC + num_loads_RL + num_loads_RC + num_loads_R +
                 num_loads_C + num_loads_L == num_loads),
            "Expect the number of loads to be identical to the sum of the loads types, but
            the number of loads is $num_loads and the sum of the loads is
            $(num_loads_RLC + num_loads_RL + num_loads_RC + num_loads_R)")
    else
        throw("Expect parameters to be a dict or nothing, not $(typeof(parameters))")
    end

    num_fltr = 4 * num_fltr_LCL + 3 * num_fltr_LC + 2 * num_fltr_L

    # Equivalent to the number of load states
    num_impedance = (2 * (num_loads_RLC + num_loads_LC + num_loads_RL + num_loads_L) +
                     num_loads_RC + num_loads_C + num_loads_R)

    num_spp = num_fltr_LCL * 4 + num_fltr_LC * 3 + num_fltr_L * 2 + num_connections +
              (num_loads_RLC + num_loads_LC + num_loads_RL + num_loads_L) * 2 +
              (num_loads_RC + num_loads_C + num_loads_R)

    p_load_total, q_load_total, s_load_total, s_source_total = CheckPowerBalance(
        parameters)

    if s_load_total > s_source_total
        if verbosity > 0
            @warn("The apparent power drawn from the loads exceeds the apparent power provided
                by all sources in steady state. Stable grid operation maybe not possible.")
        end
    end

    nc = NodeConstructor(
        num_connections,
        num_sources,
        num_loads,
        num_fltr_LCL,
        num_fltr_LC,
        num_fltr_L,
        num_loads_RLC,
        num_loads_LC,
        num_loads_RL,
        num_loads_RC,
        num_loads_L,
        num_loads_C,
        num_loads_R,
        num_impedance,
        num_fltr,
        num_spp,
        cntr,
        tot_ele,
        CM,
        parameters,
        S2S_p,
        S2L_p,
        L2L_p,
        ts,
        invoke_pfe,
        verbosity
    )

    return nc
end


"""
    GetFltrDistr(num)

Calculates the distribution of filters based on a Dirichlet distribution and the total `num`
of filters.
"""
function GetFltrDistr(num)
    # a = [.49 .02 .49] # probability for LC should be lower
    di_di = Dirichlet(ones(3))
    smpl = rand(di_di, 1) * num
    num_fltr_L = Int(floor(smpl[1]))
    num_fltr_LC = Int(ceil(clamp(smpl[2], 1, num - 1)))
    num_fltr_LCL = num - (num_fltr_LC + num_fltr_L)
    return num_fltr_L, num_fltr_LC, num_fltr_LCL
end


"""
    GetLoadsDistr(num)

Calculates the distribution of loads based on a Dirichlet distribution and the total `num`
of loads.
"""
function GetLoadsDistr(num)
    di_di = Dirichlet(ones(7)) # create dirichlet distribution
    smpl = rand(di_di, 1) * num
    num_loads_R = Int(floor(smpl[1]))
    num_loads_C = Int(floor(smpl[2]))
    num_loads_L = Int(floor(smpl[3]))
    num_loads_RL = Int(floor(smpl[4]))
    num_loads_RC = Int(floor(smpl[5]))
    num_loads_LC = Int(floor(smpl[6]))
    num_loads_RLC = num -
                    (num_loads_R + num_loads_C + num_loads_L + num_loads_RL + num_loads_RC +
                     num_loads_LC)
    return (num_loads_R, num_loads_C, num_loads_L, num_loads_RL, num_loads_RC, num_loads_LC,
        num_loads_RLC)
end


"""
    CheckParameters(parameters, num_sources, num_loads, num_connections, CM, ts, verbosity)

Checks the parameter dict for completeness.

Gets `parameters` and controls the entries based on the given inputs `num_sources`,
`num_loads`, `num_connections`, `CM` and `ts`. Messages can be suppressed by `verbosity`.

# Arguments
- `parameters::Dict`: dict containing all parameters of the node
- `num_sources::Int`: number of sources
- `num_loads::Int`: number of loads
- `num_connections::Int`: number of connections
- `CM::Matrix`: connectivity matrix describing the connections in the grid
"""
function CheckParameters(
    parameters, num_sources, num_loads, num_connections, CM, ts, invoke_pfe, verbosity
)
    # Variable generation of the parameter dicts

    # check if parameters have been specified
    if parameters === nothing
        parameters = Dict()
    end

    # check environment
    if !haskey(parameters, "weather")
        environment_properties = Dict()
        environment_properties["G"] = 1000 # Irradiance in W/m^2
        environment_properties["T"] = 23 # Temperature in °C
        parameters["weather"] = environment_properties
    else
        if !haskey(parameters["weather"], "G")
            parameters["weather"]["G"] = 1000
        end

        if !haskey(parameters["weather"], "T")
            parameters["weather"]["T"] = 23
        end
    end

    # check grid
    if !haskey(parameters, "grid")
        grid_properties = Dict()
        grid_properties["fs"] = 1 / ts
        grid_properties["v_rms"] = 230
        grid_properties["phase"] = 3
        grid_properties["f_grid"] = 50 # Gridfrequency
        grid_properties["Δfmax"] = 0.5 # The drop in frequency
        grid_properties["ΔEmax"] = 5 # The drop in rms voltage
        grid_properties["ramp_end"] = 2 / 50
        grid_properties["process_start"] = 2 / 50
        parameters["grid"] = grid_properties
    else
        if !haskey(parameters["grid"], "fs")
            parameters["grid"]["fs"] = 1 / ts
        end

        if !haskey(parameters["grid"], "v_rms")
            parameters["grid"]["v_rms"] = 230
        end

        if !haskey(parameters["grid"], "phase")
            parameters["grid"]["phase"] = 3
        end

        if !haskey(parameters["grid"], "f_grid")
            parameters["grid"]["f_grid"] = 50.0
        end

        if !haskey(parameters["grid"], "Δfmax")
            parameters["grid"]["Δfmax"] = 0.5
        end

        if !haskey(parameters["grid"], "ΔEmax")
            parameters["grid"]["ΔEmax"] = 5
        end

        if !haskey(parameters["grid"], "ramp_end")
            parameters["grid"]["ramp_end"] = 2 / parameters["grid"]["f_grid"]
        end

        if !haskey(parameters["grid"], "process_start")
            parameters["grid"]["process_start"] = 2 / parameters["grid"]["f_grid"]
        end
    end

    # check sources
    if !haskey(parameters, "source")
        num_fltr_L, num_fltr_LC, num_fltr_LCL = GetFltrDistr(num_sources)
        source_list = []

        for s in 1:num_fltr_LCL
            push!(source_list, SampleFilterLCL(parameters["grid"]))
        end

        for s in 1:num_fltr_LC
            push!(source_list, SampleFilterLC(parameters["grid"]))
        end

        for s in 1:num_fltr_L
            push!(source_list, SampleFilterL(parameters["grid"]))
        end

        parameters["source"] = source_list
    else
        num_def_sources = length(parameters["source"])
        num_undef_sources = num_sources - num_def_sources

        @assert(num_undef_sources >= 0,
            "Expect the number of defined sources within the parameter dict to be less or
            equal to the number of sources in the env, but the entries within the parameter
            dict is $num_def_sources and the number of env sources is $num_sources.")

        if num_undef_sources > 0
            if verbosity > 0
                @warn("The number of defined sources $num_def_sources is smaller than the
                    number specified sources in the environment $num_sources, therefore the
                    remaining $num_undef_sources sources are selected randomly!")
            end
        end

        num_LC_defined = 0
        num_LCL_defined = 0
        source_type_fixed = 0

        for (index, source) in enumerate(parameters["source"])
            if !haskey(source, "pwr")
                source["pwr"] = rand(range(start=5, step=5, stop=50)) * 1e3
            end

            if !haskey(source, "vdc")
                source["vdc"] = 800
            end

            if !haskey(source, "i_rip")
                source["i_rip"] = 0.15
            end

            if !haskey(source, "v_rip")
                source["v_rip"] = 0.01537
            end

            #Inductor design
            if !haskey(source, "L1")
                Vorms = parameters["grid"]["v_rms"] * 1.05
                Vop = Vorms * sqrt(2)
                Zl = 3 * Vorms^2 / source["pwr"]
                Iorms = Vorms / Zl
                Iop = Iorms * sqrt(2)
                ΔIlfmax = source["i_rip"] * Iop
                source["L1"] = (source["vdc"] * (4 * parameters["grid"]["fs"] * ΔIlfmax)^-1)
            end

            if !haskey(source, "fltr")

                if haskey(source, "mode")

                    if (source["mode"] == "Swing")

                        source["fltr"] = "L"
                    else

                        source["fltr"] = "LCL"
                    end
                else
                    source["fltr"] = "LCL"
                end

            elseif !(source["fltr"] in ["L", "LC", "LCL"])
                # TODO: Raise warning: False key
                source["fltr"] = "LCL"
                @warn "filterType not known! set to LCL filter, please choose L, LC, or LCL!"
            end

            if !haskey(source, "i_limit")
                Vorms = parameters["grid"]["v_rms"] * 1.05
                Vop = Vorms * sqrt(2)
                i_lim_r = 1.5
                Zl = 3 * Vorms^2 / source["pwr"]
                Iorms = Vorms / Zl
                Iop = Iorms * sqrt(2)

                if source["fltr"] == "LCL"
                    source["i_limit"] = 1.15 * Iop
                else
                    source["i_limit"] = i_lim_r * Iop * (1 + source["i_rip"] / 2)
                end
            end

            if !haskey(source, "R1")
                source["R1"] = 200 * source["L1"]
            end

            if (source["fltr"] == "LC" || source["fltr"] == "LCL")
                if source["fltr"] == "LC"
                    num_LC_defined += 1
                end

                if source["fltr"] == "LCL"
                    num_LCL_defined += 1
                end

                if !haskey(source, "C")
                    Vorms = parameters["grid"]["v_rms"] * 0.95
                    Vop = Vorms * sqrt(2)
                    Zl = 3 * Vorms * Vorms / source["pwr"] #TODO
                    Iorms = Vorms / Zl
                    Iop = Iorms * sqrt(2)
                    Ir_d = source["vdc"] /
                           (4 * parameters["grid"]["fs"] * source["L1"] * Iop)
                    ΔIlfmax = Ir_d * Iop
                    ΔVcfmax = source["v_rip"] * Vop
                    source["C"] = ΔIlfmax / (8 * parameters["grid"]["fs"] * ΔVcfmax)
                end

                if source["fltr"] == "LC" && !haskey(source, "R_C")
                    fc = parameters["grid"]["fs"] / 5
                    omega_c = 2 * pi * fc
                    source["R_C"] = 1 / (3 * omega_c * source["C"])
                end
            end

            if (source["fltr"] == "LC" &&
                (1 / sqrt(source["L1"] * source["C"]) > parameters["grid"]["fs"] / 2))

                if verbosity > 0
                    @warn ("The LC filter parameters have been poorly chosen.
                        The filtering capacitors should be chosen such that the resonant
                        frequency 1/sqrt(L*C) is approximately sqrt(ωn * ωs), where ωn
                        is the angular frequency of the grid, and ωs is the angular
                        switching frequency.")
                end
            end

            if !haskey(source, "v_limit")
                v_lim_r = 1.5
                source["v_limit"] = v_lim_r * source["vdc"] * (1 + source["v_rip"] / 2)
            end

            if source["fltr"] == "LCL" && !haskey(source, "L2")
                #TODO: add user warnings if the L, C, parameters they choose are stupid
                # (more than 0.5*fs) also calculate the maximal power factor variation.
                # If more than 5% add warning

                fc = parameters["grid"]["fs"] / 5
                omega_c = 2 * pi * fc

                if !haskey(source, "L2")
                    source["L2"] = source["L1"] / (omega_c^2 * source["L1"] * source["C"] - 1)
                end

                if !haskey(source, "R2")
                    source["R2"] = 200 * source["L2"]
                end

                if !haskey(source, "R_C")
                    source["R_C"] = 1 / (3 * omega_c * source["C"]) #*** filter layout
                end
            end

            if source["fltr"] == "LCL"

                fc = (1 / (2 * pi)) * sqrt((source["L1"] + source["L2"]) /
                                           (source["L1"] * source["L2"] * source["C"]))

                if fc > parameters["grid"]["fs"] / 2
                    if verbosity > 0
                        @warn("The LCL filter parameters have been poorly chosen.
                            The cut-off frequency of the filter must be minimally one half
                            of the switching frequency of the converter, because the filter
                            must have enough attenuation in the range of the converter's
                            switching frequency.")
                    end
                end
            end

            if source["source_type"] == "ideal"
                if !haskey(source, "vdc")
                    source["vdc"] = 700
                end
                if !haskey(source, "pwr")
                    source["pwr"] = 5e3
                end
            end

            if source["source_type"] == "pv"
                if !haskey(source, "vdc")
                    source["vdc"] = 700
                end
                if !haskey(source, "pwr")
                    source["pwr"] = 5e3
                end
            end

            if source["source_type"] == "battery"
                if !haskey(source, "module_R")
                    source["module_R"] = [0.010342 0.0012244
                        0.0067316 0.0011396
                        0.0051156 0.0012661
                        0.0043447 0.0012265
                        0.0038826 0.0011163
                        0.0034226 0.0009968
                        0.003346 0.0011458
                        0.0033222 0.001345
                        0.0033201 0.0013091
                        0.0032886 0.0010986
                        0.0028114 0.0010309]
                end
                if !haskey(source, "module_C")
                    source["module_C"] = [2287.7 11897
                        6122 24515
                        18460 42098
                        20975 44453
                        15254 33098
                        10440 24492
                        13903 32975
                        16694 40007
                        15784 35937
                        12165 26430
                        9118 24795]
                end
                if !haskey(source, "n")
                    source["n"] = 1
                end
                if !haskey(source, "R_0")
                    source["R_0"] = 0.016
                end
                if !haskey(source, "V_0")
                    source["V_0"] = [3.5042 3.5136
                        3.5573 3.5646
                        3.6009 3.6153
                        3.6393 3.6565
                        3.6742 3.6889
                        3.7121 3.7214
                        3.7937 3.8078
                        3.8753 3.8945
                        3.97 3.9859
                        4.0764 4.0821
                        4.1924 4.193]
                end
                if !haskey(source, "Q_0")
                    source["Q_0"] = source["n"] * 26 * 3600
                end
                if !haskey(source, "Q")
                    source["Q"] = 0
                end
                if !haskey(source, "SOC_OP")
                    source["SOC_OP"] = vec([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
                end
                if !haskey(source, "T_OP")
                    source["T_OP"] = vec([293.15 313.15])
                end
                if !haskey(source, "i_limit")
                    source["i_limit"] = 500
                end
            end

            if source["source_type"] == "pv"
                if !haskey(source, "module_I_0")
                    source["module_I_0"] = 2.0381e-10
                end
                if !haskey(source, "module_ni")
                    source["module_ni"] = 1.2
                end
                if !haskey(source, "module_N_cell")
                    source["module_N_cell"] = 36
                end
                if !haskey(source, "module_I_ph_ref")
                    source["module_I_ph_ref"] = 3.11
                end
                if !haskey(source, "serial")
                    source["serial"] = 10
                end
                if !haskey(source, "parallel")
                    source["parallel"] = 4
                end
            end

            if !haskey(source, "source_type")
                source["source_type"] = "ideal"
                source_type_fixed += 1
            end

            if !haskey(source, "τv")
                source["τv"] = 0.002
            end

            if !haskey(source, "τf")
                source["τf"] = 0.002
            end

            if !haskey(source, "pf")
                default_pf = 0.8
                if !haskey(source, "p_set") && !haskey(source, "q_set")
                    source["pf"] = default_pf
                elseif haskey(source, "q_set") && !haskey(source, "p_set")
                    p_set = sqrt(source["pwr"]^2 - source["q_set"]^2)
                    source["pf"] = p_set / source["pwr"]
                elseif haskey(source, "p_set") && !haskey(source, "q_set")
                    source["pf"] = source["p_set"] / source["pwr"]
                elseif haskey(source, "p_set") && haskey(source, "q_set")
                    s_set = sqrt(source["p_set"]^2 + source["q_set"]^2) *
                            sign(source["q_set"])
                    if s_set == 0
                        source["pf"] = 1 / sqrt(2)
                    else
                        source["pf"] = source["p_set"] / s_set
                    end
                end
            end

            if !haskey(source, "p_set")
                source["p_set"] = 0#source["pwr"]*source["pf"]
            end

            if !haskey(source, "q_set")
                source["q_set"] = 0#sqrt(source["pwr"]^2 - source["p_set"]^2)
            end

            if !haskey(source, "v_pu_set")
                source["v_pu_set"] = 1.0
            end

            if !haskey(source, "v_δ_set")
                source["v_δ_set"] = 0.0
            end

            if !haskey(source, "control_type")
                source["control_type"] = "classic"
            end

            if source["control_type"] == "classic"
                if !haskey(source, "mode")
                    source["mode"] = "Synchronverter"
                end
            elseif source["control_type"] == "RL"
                if !haskey(source, "mode")
                    source["mode"] = "ElectricGrid_ddpg"
                end
            else
                @assert("Invalid control type, please choose RL or classic")
            end

            if !haskey(source, "γ") # asymptotic mean
                source["γ"] = source["p_set"]
            end

            if !haskey(source, "std_asy") || haskey(source, "κ")# asymptotic standard deviation
                if !haskey(source, "σ")
                    source["std_asy"] = 0.0
                elseif !haskey(source, "κ")
                    source["std_asy"] = source["pwr"] / 4
                else
                    source["std_asy"] = source["σ"] / sqrt(2 * source["κ"])
                end
            end

            if !haskey(source, "κ") # mean reversion parameter
                if source["std_asy"] == 0.0
                    source["κ"] = 0.0
                else
                    source["κ"] = source["σ"]^2 / (2 * source["std_asy"]^2)
                end
            end

            if !haskey(source, "σ")
                source["σ"] = 0.0
            end

            if !haskey(source, "X₀")
                source["X₀"] = source["p_set"]
            end

            if !haskey(source, "Δt")
                steps = 1 # ... cycles for 1 step
                source["Δt"] = round(steps * parameters["grid"]["fs"] /
                                     (parameters["grid"]["f_grid"])) / parameters["grid"]["fs"]

            elseif haskey(source, "Δt")
                if typeof(source["Δt"]) == Int
                    steps = source["Δt"] # ... cycles for 1 step
                    source["Δt"] = round(steps * parameters["grid"]["fs"] / (parameters["grid"]["f_grid"])) / parameters["grid"]["fs"]
                else
                    source["Δt"] = round(source["Δt"] * (parameters["grid"]["fs"])) / parameters["grid"]["fs"]
                end
            end

            if !haskey(source, "k") # degree of polynomial
                if source["σ"] == 0
                    source["k"] = 0
                else
                    source["k"] = 2
                end
            else
                source["k"] = round(source["k"])
            end
        end

        if num_undef_sources > 0
            num_fltr_L_undef, num_fltr_LC_undef, num_fltr_LCL_undef = GetFltrDistr(
                num_undef_sources
            )

            for s in 1:num_fltr_LCL_undef
                push!(parameters["source"], SampleFilterLCL(parameters["grid"]))
            end

            for s in 1:num_fltr_LC_undef
                push!(parameters["source"], SampleFilterLC(parameters["grid"]))
            end

            for s in 1:num_fltr_L_undef
                push!(parameters["source"], SampleFilterL(parameters["grid"]))
            end
        end
        num_fltr_LCL, num_fltr_LC, num_fltr_L = CountFilters(parameters["source"])
    end

    # calculate grid power
    parameters["grid"]["pwr"] = sum(
        [parameters["source"][i]["pwr"] for i in 1:length(parameters["source"])]
    )

    # check loads
    if !haskey(parameters, "load")
        loads_distr = GetLoadsDistr(num_loads)
        num_loads_R = loads_distr[1]
        num_loads_C = loads_distr[2]
        num_loads_L = loads_distr[3]
        num_loads_RL = loads_distr[4]
        num_loads_RC = loads_distr[5]
        num_loads_LC = loads_distr[6]
        num_loads_RLC = loads_distr[7]
        load_list = []

        for l in 1:num_loads_RLC
            push!(load_list, SampleLoadRLC(parameters["grid"], num_loads))  # TODO: Sampling based on Grid_pwr
        end

        for l in 1:num_loads_LC
            push!(load_list, SampleLoadLC(parameters["grid"], num_loads)) # TODO: Sampling based on Grid_pwr
        end

        for l in 1:num_loads_RL
            push!(load_list, SampleLoadRL(parameters["grid"], num_loads)) # TODO: Sampling based on Grid_pwr
        end

        for l in 1:num_loads_L
            push!(load_list, SampleLoadL(parameters["grid"], num_loads)) # TODO: Sampling based on Grid_pwr
        end

        for l in 1:num_loads_RC
            push!(load_list, SampleLoadRC(parameters["grid"], num_loads)) # TODO: Sampling based on Grid_pwr
        end

        for l in 1:num_loads_C
            push!(load_list, SampleLoadC(parameters["grid"], num_loads)) # TODO: Sampling based on Grid_pwr
        end

        for l in 1:num_loads_R
            push!(load_list, SampleLoadR(parameters["grid"], num_loads)) # TODO: Sampling based on Grid_pwr
        end

        parameters["load"] = load_list
    else
        num_def_loads = length(parameters["load"])
        num_undef_loads = num_loads - num_def_loads

        @assert(
            num_undef_loads >= 0,
            "Expect the number of defined loads within the
parameter dict to be less or equal to the number of loads in the env, but the
entries within the parameter dict is $num_def_loads and the number of env loads
is $num_loads."
        )

        if num_undef_loads > 0
            @warn("The number of defined loads $num_def_loads is smaller than the number
                specified loads in the environment $num_loads, therefore the remaining
                $num_undef_loads loads are selected randomly!")
        end

        for (index, load) in enumerate(parameters["load"])
            if !haskey(load, "impedance")
                load["impedance"] = rand(["RLC", "RL", "RC", "LC", "R", "L", "C"])
                @warn "The type of load was not specified and is therefore drawn randomly!"

                # TODO
                # if any(keys(load) in ["R", "L", "C"]
                #   @warn "In the load values were defined which were not used because the type load["impedance"] does not consider the value/s [LIST OF VALUES THAT ARE SPEZIFIED BUT NOT USED]!"

            end

            char = split(load["impedance"], "")
            for (_, value) in enumerate(char)
                if value == "R"
                    if !haskey(load, "R")
                        load["R"] = round(rand(Uniform(10, 1e5)), digits=3)
                    end
                end

                if value == "L"
                    if !haskey(load, "L")
                        load["L"] = rand(Uniform(1e-6, 1e-3))
                    end
                end

                if value == "C"
                    if !haskey(load, "C")
                        load["C"] = rand(Uniform(1e-9, 1e-4))
                    end
                end
            end
        end

        if num_undef_loads > 0

            loads_distr_undef = GetLoadsDistr(num_undef_loads)

            num_loads_R_undef = loads_distr_undef[1]
            num_loads_C_undef = loads_distr_undef[2]
            num_loads_L_undef = loads_distr_undef[3]
            num_loads_RL_undef = loads_distr_undef[4]
            num_loads_RC_undef = loads_distr_undef[5]
            num_loads_LC_undef = loads_distr_undef[6]
            num_loads_RLC_undef = loads_distr_undef[7]

            for l in 1:num_loads_RLC_undef
                push!(parameters["load"], SampleLoadRLC(parameters["grid"], num_loads))
            end

            for l in 1:num_loads_LC_undef
                push!(parameters["load"], SampleLoadLC(parameters["grid"], num_loads))
            end

            for l in 1:num_loads_RL_undef
                push!(parameters["load"], SampleLoadRL(parameters["grid"], num_loads))
            end

            for l in 1:num_loads_L_undef
                push!(parameters["load"], SampleLoadL(parameters["grid"], num_loads))
            end

            for l in 1:num_loads_RC_undef
                push!(parameters["load"], SampleLoadRC(parameters["grid"], num_loads))
            end

            for l in 1:num_loads_C_undef
                push!(parameters["load"], SampleLoadC(parameters["grid"], num_loads))
            end

            for l in 1:num_loads_R_undef
                push!(parameters["load"], SampleLoadR(parameters["grid"], num_loads))
            end
        end
    end

    for (index, load) in enumerate(parameters["load"])
        if load["impedance"] == "R"
            load["Z"] = load["R"]
            load["pf"] = 1
        elseif load["impedance"] == "L"
            load["Z"] = 1im * 2 * pi * parameters["grid"]["f_grid"] * load["L"]
            load["pf"] = 0
        elseif load["impedance"] == "C"
            load["Z"] = 1 / (1im * 2 * pi * parameters["grid"]["f_grid"] * load["C"])
            load["pf"] = 0
        elseif load["impedance"] == "RL"
            load["Z"] = 1im * parameters["grid"]["f_grid"] * 2 * pi * load["R"] *
                        load["L"] / (load["R"] + 1im * parameters["grid"]["f_grid"] *
                                                 2 * pi * load["L"])
            load["pf"] = sign(imag(load["Z"]))abs(cos(atan(imag(load["Z"]) /
                                                           real(load["Z"]))))
        elseif load["impedance"] == "RC"
            load["Z"] = load["R"] / (1 + 1im * parameters["grid"]["f_grid"] *
                                         2 * pi * load["C"] * load["R"])
            load["pf"] = sign(imag(load["Z"]))abs(cos(atan(imag(load["Z"]) /
                                                           real(load["Z"]))))
        elseif load["impedance"] == "LC"
            load["Z"] = 1im * parameters["grid"]["f_grid"] * 2 * pi * load["L"] /
                        (1 - (parameters["grid"]["f_grid"] * 2 * pi)^2 * load["L"] * load["C"])
            load["pf"] = 0
        elseif load["impedance"] == "RLC"
            load["Z"] = 1im * parameters["grid"]["f_grid"] * 2 * pi * load["L"] /
                        (1 + 1im * parameters["grid"]["f_grid"] * 2 * pi * load["L"] / load["R"] -
                         (parameters["grid"]["f_grid"] * 2 * pi)^2 * load["L"] * load["C"])
            load["pf"] = sign(imag(load["Z"]))abs(cos(atan(imag(load["Z"]) /
                                                           real(load["Z"]))))
        end

        if !haskey(load, "pwr")
            load["pwr"] = parameters["grid"]["v_rms"]^2 / abs(load["Z"]) * parameters["grid"]["phase"]
        end

        if !haskey(load, "pwr")
            load["i_limit"] = sqrt(2) * 0.9 * parameters["grid"]["v_rms"] / abs(load["Z"])
        end
    end

    # check cables
    if !haskey(parameters, "cable")
        # no cable params defined -- invoke PFE from here ??
        cable_list = []
        for c in 1:num_connections
            push!(cable_list, SampleCable(parameters))
        end
        parameters["cable"] = cable_list

        if parameters["grid"]["pwr"] > 1e6
            @warn("Power of the network is greater than 1e6, therefore it may happen that
                the automatic cable parameterization cannot be solved.")
        end

        # invoke PFE
        if invoke_pfe
            parameters = LayoutCabels(CM, parameters, verbosity)
        end
    else
        num_def_cables = length(parameters["cable"])
        num_undef_cables = num_connections - num_def_cables

        @assert(
            num_undef_cables >= 0,
            "Expect the number of defined cables within the
parameter dict to be less or equal to the number of sources in the env, but the
entries within the parameter dict is $num_def_cables and the number of env
cables is $num_cables."
        )

        if num_undef_cables > 0
            if verbosity > 0
                @warn("The number of defined cables $num_def_cables is smaller than the
                    number specified cables in the environment $num_connections, therefore
                    the remaining $num_undef_cables cables are selected randomly!")
            end
        end

        cable_from_pfe_idx = []

        for (idx, cable) in enumerate(parameters["cable"])
            if !haskey(cable, "len")
                cable["len"] = rand(Uniform(1e-3, 1e1))
            end

            if !haskey(cable, "R") | !haskey(cable, "L") | !haskey(cable, "C")
                @info("Parameters from cable $(idx) missing. All cable parameters are
                    calculate based on power flow equation. Create a counter - we don't want
                    to see this every time.")
                push!(cable_from_pfe_idx, idx)
            end
        end

        # Solve PFE for all cables which do not have complete parameters
        if !isempty(cable_from_pfe_idx)
            println("START PFE")
            parameters = LayoutCabels(CM, parameters, verbosity)
        end

        if num_undef_cables > 0
            for c in 1:num_undef_cables
                push!(parameters["cable"], SampleCable(parameters))
            end
        end
    end

    return parameters
end


"""
    CountFilters(source_list)

Counts the number of filter types, if `parameters` is passed.
"""
function CountFilters(source_list)
    cntr_LCL = 0
    cntr_LC = 0
    cntr_L = 0

    for (i, source) in enumerate(source_list)
        if source["fltr"] == "LCL"
            cntr_LCL += 1
        elseif source["fltr"] == "LC"
            cntr_LC += 1
        elseif source["fltr"] == "L"
            cntr_L += 1
        end
    end

    return cntr_LCL, cntr_LC, cntr_L
end


"""
    CountLoads(load_list)

Counts the number of load types, if `parameters` is passed.
"""
function CountLoads(load_list)
    cntr_RLC = 0
    cntr_LC = 0
    cntr_RL = 0
    cntr_RC = 0
    cntr_L = 0
    cntr_C = 0
    cntr_R = 0

    for (i, source) in enumerate(load_list)
        if source["impedance"] == "RLC"
            cntr_RLC += 1
        elseif source["impedance"] == "LC"
            cntr_LC += 1
        elseif source["impedance"] == "RL"
            cntr_RL += 1
        elseif source["impedance"] == "RC"
            cntr_RC += 1
        elseif source["impedance"] == "L"
            cntr_L += 1
        elseif source["impedance"] == "C"
            cntr_C += 1
        elseif source["impedance"] == "R"
            cntr_R += 1
        end
    end

    return cntr_RLC, cntr_LC, cntr_RL, cntr_RC, cntr_L, cntr_C, cntr_R
end


"""
    SampleFilterLCL()

Sample parameters for the LCL filter.
"""
function SampleFilterLCL(grid_properties)
    source = SampleSource(grid_properties, "LCL")

    Lf_1, Lf_2, Cf, fc, R_1, R_2, R_C, i_limit, v_limit = FilterDesign(
        source["pwr"],
        grid_properties["fs"],
        source["fltr"];
        Vrms=grid_properties["v_rms"],
        Vdc=source["vdc"],
        ΔILf_ILf=source["i_rip"],
        ΔVCf_VCf=source["v_rip"]
    )

    source["L1"] = Lf_1
    source["L2"] = Lf_2
    source["C"] = Cf
    source["R1"] = R_1
    source["R2"] = R_2
    source["R_C"] = R_C
    source["i_limit"] = i_limit
    source["v_limit"] = v_limit
    source["control_type"] == "classic"
    source["mode"] = "Synchronverter"

    return source
end


"""
    SampleFilterLC()

Sample parameters for the LC filter.
"""
function SampleFilterLC(grid_properties)
    source = SampleSource(grid_properties, "LC")

    Lf_1, Cf, fc, R_1, R_C, i_limit, v_limit = FilterDesign(
        source["pwr"],
        grid_properties["fs"],
        source["fltr"];
        Vrms=grid_properties["v_rms"],
        Vdc=source["vdc"],
        ΔILf_ILf=source["i_rip"],
        ΔVCf_VCf=source["v_rip"]
    )

    source["L1"] = Lf_1
    source["C"] = Cf
    source["R1"] = R_1
    source["R_C"] = R_C
    source["i_limit"] = i_limit
    source["v_limit"] = v_limit
    source["control_type"] == "classic"
    source["mode"] = "Synchronverter"

    return source
end


"""
    SampleFilterL()

Sample parameters for the L filter.
"""
function SampleFilterL(grid_properties)
    source = SampleSource(grid_properties, "L")

    Lf_1, R_1, i_limit = FilterDesign(
        source["pwr"],
        grid_properties["fs"],
        source["fltr"];
        Vrms=grid_properties["v_rms"],
        Vdc=source["vdc"],
        ΔILf_ILf=source["i_rip"],
        ΔVCf_VCf=source["v_rip"])

    source["L1"] = Lf_1
    source["R1"] = R_1
    source["i_limit"] = i_limit
    source["v_limit"] = 1.1 * grid_properties["v_rms"] * sqrt(2)
    source["control_type"] == "classic"
    source["mode"] = "Synchronverter"

    return source
end


"""
    SampleLoadRLC()

Sample parameters for the RLC load.
"""
function SampleLoadRLC(grid_properties, num_loads)
    load = Dict()
    a = 0.3
    b = 1 - a
    S1 = grid_properties["pwr"] / num_loads * 0.7 * a
    S2 = grid_properties["pwr"] / num_loads * 0.7 * b
    pf1 = round(rand(Uniform(0.95, 0.99)), digits=3)
    pf2 = round(rand(Uniform(0.95, 0.99)), digits=3)

    R1, L, _, _ = ParallelLoadImpedance(S1, pf1, grid_properties["v_rms"],
        fsys=grid_properties["f_grid"], type_sign="L")
    R2, C, _, _ = ParallelLoadImpedance(S2, -pf2, grid_properties["v_rms"],
        fsys=grid_properties["f_grid"], type_sign="C")

    load["impedance"] = "RLC"
    load["R"] = R1 + R2
    load["L"] = L
    load["C"] = C

    return load
end


"""
    SampleLoadLC()

Sample parameters for the LC load.
"""
function SampleLoadLC(grid_properties, num_loads)
    load = Dict()
    a = 0.3
    b = 1 - a

    S1 = grid_properties["pwr"] / num_loads * 0.7 * a
    S2 = grid_properties["pwr"] / num_loads * 0.7 * b

    _, L, _, _ = ParallelLoadImpedance(S1, 0, grid_properties["v_rms"],
        fsys=grid_properties["f_grid"], type_sign="L")
    _, C, _, _ = ParallelLoadImpedance(S2, -0, grid_properties["v_rms"],
        fsys=grid_properties["f_grid"], type_sign="C")

    load["impedance"] = "LC"
    load["L"] = L
    load["C"] = C

    return load
end


"""
    SampleLoadRL()

Sample parameters for the RL load.
"""
function SampleLoadRL(grid_properties, num_loads)
    load = Dict()
    pf = round(rand(Uniform(0.95, 0.99)), digits=3)
    S = grid_properties["pwr"] / num_loads * 0.7

    R, L, _, _ = ParallelLoadImpedance(S, pf, grid_properties["v_rms"],
        fsys=grid_properties["f_grid"])

    load["impedance"] = "RL"
    load["R"] = R
    load["L"] = L

    return load
end


"""
    SampleLoadRC()

Sample parameters for the RC load.
"""
function SampleLoadRC(grid_properties, num_loads)
    load = Dict()
    pf = round(rand(Uniform(0.9, 0.98)), digits=3)
    S = grid_properties["pwr"] / num_loads * 0.7

    R, C, _, _ = ParallelLoadImpedance(S, -pf, grid_properties["v_rms"],
        fsys=grid_properties["f_grid"])

    load["impedance"] = "RC"
    load["R"] = R
    load["C"] = C

    return load
end


"""
    SampleLoadL()

Sample parameters for the L load.
"""
function SampleLoadL(grid_properties, num_loads)
    load = Dict()
    load["impedance"] = "L"
    pf = 0
    S = grid_properties["pwr"] / num_loads * 0.7

    _, L, _, _ = ParallelLoadImpedance(S, pf, grid_properties["v_rms"],
        fsys=grid_properties["f_grid"], type_sign=load["impedance"])

    load["L"] = L

    return load
end


"""
    SampleLoadC()

Sample parameters for the C load.
"""
function SampleLoadC(grid_properties, num_loads)
    load = Dict()
    load["impedance"] = "C"
    pf = 0
    S = grid_properties["pwr"] / num_loads * 0.7

    _, C, _, _ = ParallelLoadImpedance(S, -pf, grid_properties["v_rms"],
        fsys=grid_properties["f_grid"], type_sign=load["impedance"])

    load["C"] = C

    return load
end


"""
    SampleLoadR()

Sample parameters for the R load.
"""
function SampleLoadR(grid_properties, num_loads)
    load = Dict()
    pf = 1
    S = grid_properties["pwr"] / num_loads * 0.7

    R, _, _, _ = ParallelLoadImpedance(S, pf, grid_properties["v_rms"],
        fsys=grid_properties["f_grid"])

    load["impedance"] = "R"
    load["R"] = R

    return load
end


"""
    SampleCable()

Sample parameters for the cable.
"""
function SampleCable(parameters)
    cable = Dict()
    cable["len"] = 1.0#rand(Uniform(1e-3, 1e1))

    cable["R"] = 0.208   # Ω, line resistance
    cable["L"] = 0.00025 # H, line inductance
    cable["C"] = 0.4e-3  # F, line capacitance
    cable["i_limit"] = 10e12   # limits will be overwritten if PFE is solved
    cable["v_limit"] = 1.15 * parameters["grid"]["v_rms"] * sqrt(2)

    cable["Rb"] = 0.722 / cable["len"]
    cable["Cb"] = 0.4e-6 / cable["len"]
    cable["Lb"] = 0.264e-3 / cable["len"]

    return cable
end


"""
    SetConnection(cntr, x, p)

Sets x to zero or to the value of the counter as a function of p and increases it as well.
"""
function SetConnection(cntr, x, p)
    if x < p
        cntr += 1
        return cntr, cntr
    else
        x = 0
        return cntr, x
    end
end


"""
    GenerateCM(num_sources, num_loads, S2L_p, S2S_p, L2L_p)

Returns the constructed `CM` and the total number of connections.

Gets `num_sources` and `num_loads` to calculate the total number of elements. Depending on
the probabilities `S2L_p`, `S2S_p` and `L2L_p`, the entries are then set in the CM matrix.
After the entries have been set randomly, it is checked that all elements have at least one
connection and that no subnets have been created.

# Arguments
- `num_sources::Int`: number of sources
- `num_loads::Int`: number of loads
- `S2L_p::Int`: probability that a source is connected to a load
- `S2S_p::Int`: probability that a source is connected to a source
- `L2L_p::Int`: probability that a load is connected to a load

# Return Values
- `cntr::Int`: number of connections
- `CM::Matrix`: connectivity matrix describing the connections in the grid
"""
function GenerateCM(num_sources, num_loads, S2L_p, S2S_p, L2L_p)
    # counting the connections
    cntr = 0

    # get total elements
    tot_ele = num_sources + num_loads

    # get a upper triangular matrix
    mask = UpperTriangular(ones(tot_ele, tot_ele))
    CM = rand(tot_ele, tot_ele) .* mask # fill matrix with random entries between [0,1]
    CM = CM - Diagonal(CM) # delet diagonal bc no connection with itself

    # go through the matrix
    # -1 bc last entry is 0 anyway
    for i in 1:tot_ele-1
        # start at i, bc we need to check only upper triangle
        if i <= num_sources
            for j in i:tot_ele-1
                if j > num_sources - 1  # select propability according to column
                    cntr, x = SetConnection(cntr, CM[i, j+1], S2L_p)
                    CM[i, j+1] = x
                else
                    cntr, x = SetConnection(cntr, CM[i, j+1], S2S_p)
                    CM[i, j+1] = x
                end
            end
        else
            for j in i:tot_ele-1
                cntr, x = SetConnection(cntr, CM[i, j+1], L2L_p)
                CM[i, j+1] = x
            end
        end
    end

    # make sure that no objects disappear or subnets are formed
    if (S2L_p < 1) || (num_loads == 0)
        for i in 1:tot_ele
            # save rows and columns entries
            Col = CM[1:i-1, i]
            Row = CM[i, i+1:tot_ele]

            # get one list in the form of: [column, row]-entries
            entries = vcat(Col, Row)
            non_zero = count(i -> (i != 0), entries) # number of non_zero entries
            zero = count(i -> (i == 0), entries) # number of zero entries
            val_to_set = min(2, zero) # minimum of connections is 2

            if non_zero <= 2 # we need to set values if there are less then 2 entries
                idx_row_entries = findall(x -> x == 0, Col) # Get rows of the entries = 0
                idx_col_entries = findall(x -> x == 0, Row) # Get col of the entries = 0
                idx_list = vcat([(j, i) for j in idx_row_entries],
                    [(i, i + j) for j in idx_col_entries])
                samples = min(val_to_set, length(idx_list))

                # draw samples from the list
                idx_rnd = sample(1:length(idx_list), samples, replace=false)

                for (a, ix) in enumerate(idx_rnd)
                    # Based on the random sample, select an indize
                    # from the list and write into the corresponding CM cell.
                    cntr += 1
                    CM[idx_list[ix][1], idx_list[ix][2]] = cntr
                end
            end
        end
    end

    CM = CM - CM' # copy with negative sign to lower triangle

    return cntr, CM
end


"""
    GetASource(self::NodeConstructor, source_i)

Create the A_src entry for a source in the A matrix.
"""
function GetASource(self::NodeConstructor, source_i)
    parameter_i = self.parameters["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        A_src = zeros(4, 4)
        A_src[1, 1] = -(parameter_i["R1"] + parameter_i["R_C"]) / parameter_i["L1"]
        A_src[1, 2] = -1 / parameter_i["L1"]
        A_src[1, 3] = parameter_i["R_C"] / parameter_i["L1"]
        A_src[2, 1] = 1 / parameter_i["C"]
        A_src[2, 3] = -1 / parameter_i["C"]
        A_src[3, 1] = parameter_i["R_C"] / parameter_i["L2"]
        A_src[3, 2] = 1 / parameter_i["L2"]
        A_src[3, 3] = -(parameter_i["R_C"] + parameter_i["R2"]) / parameter_i["L2"]
        A_src[3, 4] = -1 / parameter_i["L2"]

        C_sum = 0
        CM_row = self.CM[source_i, :]
        indizes = CM_row[CM_row.!=0]
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_src[4, 3] = C_sum^(-1)
    elseif parameter_i["fltr"] == "LC"
        A_src = zeros(3, 3)
        A_src[1, 1] = -(parameter_i["R1"]) / parameter_i["L1"]
        A_src[1, 3] = -1 / parameter_i["L1"]
        A_src[2, 2] = -1 / (parameter_i["C"] * parameter_i["R_C"])
        A_src[2, 3] = 1 / (parameter_i["C"] * parameter_i["R_C"])

        C_sum = 0
        CM_row = self.CM[source_i, :]
        indizes = CM_row[CM_row.!=0]
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_src[3, 1] = C_sum^(-1)
        A_src[3, 2] = 1 / parameter_i["R_C"] * (C_sum)^(-1)
        A_src[3, 3] = -1 / parameter_i["R_C"] * (C_sum)^(-1)

    elseif parameter_i["fltr"] == "L"
        A_src = zeros(2, 2)
        A_src[1, 1] = -parameter_i["R1"] / parameter_i["L1"]
        A_src[1, 2] = -1 / parameter_i["L1"]

        C_sum = 0
        CM_row = self.CM[source_i, :]
        indizes = CM_row[CM_row.!=0]
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_src[2, 1] = C_sum^(-1)
    else
        throw("Expect filter to be \"LCL\", \"LC\" or \"L\", not $(parameter_i["fltr"]).")
    end

    return A_src
end


"""
    GetBSource(self::NodeConstructor, source_i)

Create the B_source entry for a source in the B matrix.
"""
function GetBSource(self::NodeConstructor, source_i)
    parameter_i = self.parameters["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        B_source = zeros(4, 1)
        B_source[1, 1] = 1 / parameter_i["L1"]
    elseif parameter_i["fltr"] == "LC"
        B_source = zeros(3, 1)
        B_source[1, 1] = 1 / parameter_i["L1"]
    elseif parameter_i["fltr"] == "L"
        B_source = zeros(2, 1)
        B_source[1, 1] = 1 / parameter_i["L1"]
    end

    return B_source
end


"""
    GetASourceTrnC(self::NodeConstructor, source_i)

Create the A_src_trn_c entry in the A matrix.
"""
function GetASourceTrnC(self::NodeConstructor, source_i)
    parameter_i = self.parameters["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        A_src_trn_c = zeros(4, self.num_connections)
        CM_row = self.CM[source_i, :]
        indizes = CM_row[CM_row.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        C_sum = 0

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_c[4, idx] = sign * -(C_sum^(-1))
        end
    elseif parameter_i["fltr"] == "LC"
        A_src_trn_c = zeros(3, self.num_connections)
        CM_row = self.CM[source_i, :]
        indizes = CM_row[CM_row.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        C_sum = 0

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_c[3, idx] = sign * -(C_sum^(-1))
        end
    elseif parameter_i["fltr"] == "L"
        A_src_trn_c = zeros(2, self.num_connections)
        CM_row = self.CM[source_i, :]
        indizes = CM_row[CM_row.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        C_sum = 0

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_c[2, idx] = sign * -(C_sum^(-1))
        end
    end

    return A_src_trn_c
end


"""
    GetASourceTrnL(self::NodeConstructor, source_i)

Create the A_src_trn_l entry in the A matrix.
"""
function GetASourceTrnL(self::NodeConstructor, source_i)
    parameter_i = self.parameters["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        A_src_trn_l = zeros(4, self.num_connections)
        CM_col = self.CM[source_i, :]
        indizes = CM_col[CM_col.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_l[4, idx] = sign * 1 / self.parameters["cable"][idx]["L"]
        end

    elseif parameter_i["fltr"] == "LC"
        A_src_trn_l = zeros(3, self.num_connections)
        CM_col = self.CM[source_i, :]
        indizes = CM_col[CM_col.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_l[3, idx] = sign * 1 / self.parameters["cable"][idx]["L"]
        end

    elseif parameter_i["fltr"] == "L"
        A_src_trn_l = zeros(2, self.num_connections)
        CM_col = self.CM[source_i, :]
        indizes = CM_col[CM_col.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_l[2, idx] = sign * 1 / self.parameters["cable"][idx]["L"]
        end
    end

    return A_src_trn_l'
end


"""
    GenerateATrn(self::NodeConstructor)

Create the A_tran entry in the A matrix.
"""
function GenerateATrn(self::NodeConstructor)
    vec = zeros(self.num_connections)
    for (i, ele) in enumerate(self.parameters["cable"])
        vec[i] = -ele["R"] / ele["L"]
    end

    A_trn = Diagonal(vec)
    return A_trn
end


"""
    GenerateATrnLoadC(self::NodeConstructor, load_i)

Create the A_tran_load_c entry in the A matrix.
"""
function GenerateATrnLoadC(self::NodeConstructor, load_i)
    parameter_i = self.parameters["load"][load_i]

    if parameter_i["impedance"] == "RLC" || parameter_i["impedance"] == "LC"
        A_tran_load_c = zeros(2, self.num_connections)
        CM_row = self.CM[self.num_sources+load_i, :]
        indizes = CM_row[CM_row.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        C_sum = parameter_i["C"]

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_tran_load_c[1, idx] = sign * -(C_sum^(-1))
        end

    elseif parameter_i["impedance"] == "RL" || parameter_i["impedance"] == "L"
        A_tran_load_c = zeros(2, self.num_connections)
        CM_row = self.CM[self.num_sources+load_i, :]
        indizes = CM_row[CM_row.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        C_sum = 0

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_tran_load_c[1, idx] = sign * -(C_sum^(-1))
        end

    elseif parameter_i["impedance"] == "RC" || parameter_i["impedance"] == "C"
        A_tran_load_c = zeros(1, self.num_connections)
        CM_row = self.CM[self.num_sources+load_i, :]
        indizes = CM_row[CM_row.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        C_sum = parameter_i["C"]

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_tran_load_c[idx] = sign * -(C_sum^(-1))
        end

    elseif parameter_i["impedance"] == "R"
        A_tran_load_c = zeros(1, self.num_connections)
        CM_row = self.CM[self.num_sources+load_i, :]
        indizes = CM_row[CM_row.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        C_sum = 0

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_tran_load_c[idx] = sign * -(C_sum^(-1))
        end
    end

    return A_tran_load_c
end


"""
    GenerateATrnLoadL(self::NodeConstructor, load_i)

Create the A_tran_load_l entry in the A matrix.
"""
function GenerateATrnLoadL(self::NodeConstructor, load_i)
    parameter_i = self.parameters["load"][load_i]

    if (parameter_i["impedance"] == "RLC" || parameter_i["impedance"] == "LC" ||
        parameter_i["impedance"] == "RL" || parameter_i["impedance"] == "L")

        A_tran_load_l = zeros(self.num_connections, 2)
        CM_col = self.CM[self.num_sources+load_i, :]
        indizes = CM_col[CM_col.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)

            A_tran_load_l[idx, 1] = sign * 1 / self.parameters["cable"][idx]["L"]
        end

    elseif (parameter_i["impedance"] == "RC" || parameter_i["impedance"] == "C" ||
            parameter_i["impedance"] == "R")

        A_tran_load_l = zeros(self.num_connections, 1)
        CM_col = self.CM[self.num_sources+load_i, :]
        indizes = CM_col[CM_col.!=0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)

            A_tran_load_l[idx, 1] = sign * 1 / self.parameters["cable"][idx]["L"]
        end
    end

    return A_tran_load_l
end


"""
    GenerateALoad(self::NodeConstructor, load_i)

Create the A_tran_load_l entry in the A matrix.
"""
function GenerateALoad(self::NodeConstructor, load_i)
    parameter_i = self.parameters["load"][load_i]

    if parameter_i["impedance"] == "RLC"
        A_load = zeros(2, 2)
        A_load[2, 1] = 1 / parameter_i["L"]
        C_sum = parameter_i["C"]
        CM_row = self.CM[self.num_sources+load_i, :]
        indizes = CM_row[CM_row.!=0]
        signs = [sign(x) for x in indizes]
        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1, 1] = -((parameter_i["R"]) * C_sum)^(-1)
        A_load[1, 2] = -(C_sum)^(-1)

    elseif parameter_i["impedance"] == "LC"
        A_load = zeros(2, 2)
        A_load[2, 1] = 1 / parameter_i["L"]
        C_sum = parameter_i["C"]
        CM_row = self.CM[self.num_sources+load_i, :]
        indizes = CM_row[CM_row.!=0]
        signs = [sign(x) for x in indizes]
        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1, 2] = -(C_sum)^(-1)
    elseif parameter_i["impedance"] == "RL"
        A_load = zeros(2, 2)
        A_load[2, 1] = 1 / parameter_i["L"]
        C_sum = 0
        CM_row = self.CM[self.num_sources+load_i, :]
        indizes = CM_row[CM_row.!=0]
        signs = [sign(x) for x in indizes]
        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1, 1] = -((parameter_i["R"]) * C_sum)^(-1)
        A_load[1, 2] = -(C_sum)^(-1)
    elseif parameter_i["impedance"] == "RC"
        A_load = zeros(1, 1)
        C_sum = parameter_i["C"]
        CM_row = self.CM[self.num_sources+load_i, :]
        indizes = CM_row[CM_row.!=0]
        signs = [sign(x) for x in indizes]
        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1, 1] = -((parameter_i["R"]) * C_sum)^(-1)

    elseif parameter_i["impedance"] == "L"
        A_load = zeros(2, 2)
        A_load[2, 1] = 1 / parameter_i["L"]
        C_sum = 0
        CM_row = self.CM[self.num_sources+load_i, :]
        indizes = CM_row[CM_row.!=0]
        signs = [sign(x) for x in indizes]
        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1, 2] = -(C_sum)^(-1)
    elseif parameter_i["impedance"] == "C"
        A_load = zeros(1, 1)
    elseif parameter_i["impedance"] == "R"
        A_load = zeros(1, 1)
        C_sum = 0
        CM_row = self.CM[self.num_sources+load_i, :]
        indizes = CM_row[CM_row.!=0]
        signs = [sign(x) for x in indizes]
        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1, 1] = -((parameter_i["R"]) * C_sum)^(-1)
    else
        throw("Expect Impedance to be \"RLC\", \"LC\", \"RL\", \"RC\", \"L\", \"C\" or
            \"R\", not $(parameter_i["impedance"]).")
    end

    return A_load
end


"""
    GenerateA(self::NodeConstructor)

Returns the A matrix by joining the individual sub-matrices together.
"""
function GenerateA(self::NodeConstructor)
    """
    The previously constructed matrices are now plugged together in the form:
        [[      A_src,   A_src_trn_c,             0],
         [A_src_trn_l,         A_trn, A_tran_load_l],
         [          0, A_tran_load_c,        A_load]]

    with A_src:
        [[LCL, 0 , 0]
         [0,   LC, 0]
         [0,   0,  L]]

    with A_load:
         [[RLC, 0 , 0,  0,  0, 0, 0]
          [0,   LC, 0,  0,  0, 0, 0]
          [0,   0,  RL, 0,  0, 0, 0]
          [0,   0,  0,  L,  0, 0, 0]
          [0,   0 , 0,  0, RC, 0, 0]
          [0,   0 , 0,  0,  0, C, 0]
          [0,   0 , 0,  0,  0, 0, R]]
    """
    # get A_src
    self.num_fltr = 4 * self.num_fltr_LCL + 3 * self.num_fltr_LC + 2 * self.num_fltr_L
    A_src = zeros(self.num_fltr, self.num_fltr) # construct matrix of zeros
    A_src_list = [GetASource(self, i) for i in 1:self.num_sources]

    # TODO: maybe do in one loop like:
    #   for (i, vals) in enumerate(zip([1, 4, 2, 5], 2:12, (:a, :b, :c)))
    cntr = 0
    for (i, ele) in enumerate(A_src_list)
        #if i <= self.num_fltr_LCL

        if self.parameters["source"][i]["fltr"] == "LCL"
            start = 1 + cntr
            stop = 4 + cntr
            cntr += 4
            A_src[start:stop, start:stop] = ele
        elseif self.parameters["source"][i]["fltr"] == "LC"
            start = 1 + cntr
            stop = 3 + cntr
            cntr += 3
            A_src[start:stop, start:stop] = ele
        elseif self.parameters["source"][i]["fltr"] == "L"
            start = 1 + cntr
            stop = 2 + cntr
            cntr += 2
            A_src[start:stop, start:stop] = ele
        end
    end

    # get A_src_trn_c
    A_src_trn_c = zeros(self.num_fltr, self.num_connections)

    # start at 1 bc Source 1
    A_src_trn_c_list = [GetASourceTrnC(self, i) for i in 1:self.num_sources]
    cntr = 0

    for (i, ele) in enumerate(A_src_trn_c_list)
        if self.parameters["source"][i]["fltr"] == "LCL"
            start = 1 + cntr
            stop = 4 + cntr
            cntr += 4
            A_src_trn_c[start:stop, :] = ele
        elseif self.parameters["source"][i]["fltr"] == "LC"
            start = 1 + cntr
            stop = 3 + cntr
            cntr += 3
            A_src_trn_c[start:stop, :] = ele
        elseif self.parameters["source"][i]["fltr"] == "L"
            start = 1 + cntr
            stop = 2 + cntr
            cntr += 2
            A_src_trn_c[start:stop, :] = ele
        end
    end

    # get A_src_trn_l
    A_src_trn_l = zeros(self.num_connections, self.num_fltr)

    # start at 1 bc Source 1
    A_src_trn_l_list = [GetASourceTrnL(self, i) for i in 1:self.num_sources]

    cntr = 0
    for (i, ele) in enumerate(A_src_trn_l_list)
        if self.parameters["source"][i]["fltr"] == "LCL"
            start = 1 + cntr
            stop = 4 + cntr
            cntr += 4
            A_src_trn_l[:, start:stop] = ele
        elseif self.parameters["source"][i]["fltr"] == "LC"
            start = 1 + cntr
            stop = 3 + cntr
            cntr += 3
            A_src_trn_l[:, start:stop] = ele
        elseif self.parameters["source"][i]["fltr"] == "L"
            start = 1 + cntr
            stop = 2 + cntr
            cntr += 2
            A_src_trn_l[:, start:stop] = ele
        end
    end

    A_trn = GenerateATrn(self)

    if self.num_loads > 0
        A_tran_load_l_list = []

        for i in 1:self.num_loads
            push!(A_tran_load_l_list, GenerateATrnLoadL(self, i))
        end

        A_tran_load_l = reduce(hcat, A_tran_load_l_list) # i-> idx // i+1 -> num of load
        A_tran_load_c_list = []

        for i in 1:self.num_loads
            push!(A_tran_load_c_list, GenerateATrnLoadC(self, i))
        end

        A_tran_load_c = reduce(vcat, A_tran_load_c_list)
        A_load_diag = zeros(self.num_impedance, self.num_impedance)
        A_load_list = [GenerateALoad(self, i) for i in 1:self.num_loads]

        cntr = 0
        for (i, ele) in enumerate(A_load_list)
            if (self.parameters["load"][i]["impedance"] == "RLC" ||
                self.parameters["load"][i]["impedance"] == "LC" ||
                self.parameters["load"][i]["impedance"] == "RL" ||
                self.parameters["load"][i]["impedance"] == "L")

                start = 1 + cntr
                stop = 2 + cntr
                cntr += 2
                A_load_diag[start:stop, start:stop] = ele

            elseif (self.parameters["load"][i]["impedance"] == "RC" ||
                    self.parameters["load"][i]["impedance"] == "C" ||
                    self.parameters["load"][i]["impedance"] == "R")

                start = 1 + cntr
                stop = 1 + cntr
                cntr += 1
                start = (i + self.num_loads_RLC + self.num_loads_LC + self.num_loads_RL +
                         self.num_loads_L)
                A_load_diag[start:start, start:start] = ele
            end
        end
        A_load_zeros = zeros(self.num_fltr, self.num_impedance)
        A_load_zeros_t = A_load_zeros'
    end

    if self.num_loads == 0
        A = [A_src A_src_trn_c
            A_src_trn_l A_trn]
    else
        A = [A_src A_src_trn_c A_load_zeros
            A_src_trn_l A_trn A_tran_load_l
            A_load_zeros_t A_tran_load_c A_load_diag]
    end

    if self.parameters["grid"]["phase"] === 1
        return A
    elseif self.parameters["grid"]["phase"] === 3
        z = zeros(size(A))
        A_ = [A z z; z A z; z z A]
        return A_
    end
end


"""
    GenerateB(self::NodeConstructor)

Returns the B matrix by joining the individual sub-matrices together.
"""
function GenerateB(self::NodeConstructor)
    """
    The previously constructed matrices are now plugged together in the form:
        [[B_source_1,          0, ...,           0],
         [         0, B_source_2, ...,           0],
         [         0,          0, ...,           0],
         [         0,          0, ...,  B_source_n]]

    """
    B = zeros(self.num_fltr + self.num_connections + self.num_impedance, self.num_sources)

    # start at 1 bc Source 1
    B_source_list = [GetBSource(self, i) for i in 1:self.num_sources]
    cntr = 0

    for (i, ele) in enumerate(B_source_list)
        if self.parameters["source"][i]["fltr"] == "LCL"
            start = 1 + cntr
            stop = 4 + cntr
            cntr += 4
            B[start:stop, i] = ele
        elseif self.parameters["source"][i]["fltr"] == "LC"
            start = 1 + cntr
            stop = 3 + cntr
            cntr += 3
            B[start:stop, i] = ele
        elseif self.parameters["source"][i]["fltr"] == "L"
            start = 1 + cntr
            stop = 2 + cntr
            cntr += 2
            B[start:stop, i:i] = ele
        end
    end

    if self.parameters["grid"]["phase"] === 1
        return B
    elseif self.parameters["grid"]["phase"] === 3
        z = zeros(size(B))
        B_ = [B z z; z B z; z z B]
        return B_
    end
end


"""
    GenerateC(self::NodeConstructor)

Returns the C matrix.
"""
function GenerateC(self::NodeConstructor)
    C = Diagonal(ones(self.num_fltr + self.num_connections + self.num_impedance))

    if self.parameters["grid"]["phase"] === 1
        return C
    elseif self.parameters["grid"]["phase"] === 3
        z = zeros(size(C))
        C_ = [C z z
            z C z
            z z C]
        return C_
    end
end

"""
    GenerateD(self::NodeConstructor)

Returns the D matrix.
"""
function GenerateD(self::NodeConstructor)
    return 0
end


"""
    GetSystem(self::NodeConstructor)

Returns the system matrices A, B, C and D.
"""
function GetSystem(self::NodeConstructor)
    # Returns state space matrices
    A = GenerateA(self)
    B = GenerateB(self)
    C = GenerateC(self)
    D = GenerateD(self)
    return (A, B, C, D)
end

"""
    GetStateIds(self::NodeConstructor)

Creates the State Vector for an related NodeConstructor and outputs it as a list of strings.
"""
function GetStateIds(self::NodeConstructor)
    states = []
    for s in 1:self.num_sources
        if self.parameters["source"][s]["fltr"] == "LCL"
            push!(states, "source$s" * "_i_L1")
            push!(states, "source$s" * "_v_C_filt")
            push!(states, "source$s" * "_i_L2")
            push!(states, "source$s" * "_v_C_cables")
        elseif self.parameters["source"][s]["fltr"] == "LC"
            push!(states, "source$s" * "_i_L1")
            push!(states, "source$s" * "_v_C_filt")
            push!(states, "source$s" * "_v_C_cables")
        elseif self.parameters["source"][s]["fltr"] == "L"
            push!(states, "source$s" * "_i_L1")
            push!(states, "source$s" * "_v_C_cables")
        end
    end

    for c in 1:self.num_connections
        push!(states, "cable$c" * "_i_L")
    end

    for l in 1:self.num_loads
        if (self.parameters["load"][l]["impedance"] == "RLC" ||
            self.parameters["load"][l]["impedance"] == "LC" ||
            self.parameters["load"][l]["impedance"] == "RL" ||
            self.parameters["load"][l]["impedance"] == "L")
            push!(states, "load$l" * "_v_C_total")
            push!(states, "load$l" * "_i_L")
        elseif (self.parameters["load"][l]["impedance"] == "RC" ||
                self.parameters["load"][l]["impedance"] == "C" ||
                self.parameters["load"][l]["impedance"] == "R")
            push!(states, "load$l" * "_v_C_total")
        end
    end

    if self.parameters["grid"]["phase"] === 3
        A = ["_a", "_b", "_c"]
        states = vcat([broadcast(*, states, A[i]) for i in 1:3]...)
    end

    return states
end


"""
    GetStateParameters(self::NodeConstructor)

Creates a Vector containing the related L or C Parameters for the states of a
NodeConstructor.
"""
function GetStateParameters(self::NodeConstructor)
    state_paras = []
    for s in 1:self.num_sources
        if self.parameters["source"][s]["fltr"] == "LCL"
            push!(state_paras, self.parameters["source"][s]["L1"])
            push!(state_paras, self.parameters["source"][s]["C"])
            push!(state_paras, self.parameters["source"][s]["L2"])
            push!(state_paras, GetCSumCableNode(s, self))
        elseif self.parameters["source"][s]["fltr"] == "LC"
            push!(state_paras, self.parameters["source"][s]["L1"])
            push!(state_paras, self.parameters["source"][s]["C"])
            push!(state_paras, GetCSumCableNode(s, self))
        elseif self.parameters["source"][s]["fltr"] == "L"
            push!(state_paras, self.parameters["source"][s]["L1"])
            push!(state_paras, GetCSumCableNode(s, self))
        end
    end

    for c in 1:self.num_connections
        push!(state_paras, self.parameters["cable"][c]["L"])
    end

    for l in 1:self.num_loads
        if (self.parameters["load"][l]["impedance"] == "RLC" ||
            self.parameters["load"][l]["impedance"] == "LC" ||
            self.parameters["load"][l]["impedance"] == "RL" ||
            self.parameters["load"][l]["impedance"] == "L")
            c = 0

            if haskey(self.parameters["load"][l], "C")
                c = self.parameters["load"][l]["C"]
            end

            push!(state_paras, GetCSumCableNode(self.num_sources + l, self) + c)
            push!(state_paras, self.parameters["load"][l]["L"])
        elseif (self.parameters["load"][l]["impedance"] == "RC" ||
                self.parameters["load"][l]["impedance"] == "C" ||
                self.parameters["load"][l]["impedance"] == "R")
            c = 0

            if haskey(self.parameters["load"][l], "C")
                c = self.parameters["load"][l]["C"]
            end

            push!(state_paras, GetCSumCableNode(self.num_sources + l, self) + c)
        end
    end

    if self.parameters["grid"]["phase"] === 3
        state_paras = vcat([state_paras for i in 1:3]...)
    end

    return state_paras
end

"""
    GetCSumCableNode(node_i, self::NodeConstructor)

Returns the sum of the capacities at a node point.
"""
function GetCSumCableNode(node_i, self::NodeConstructor)
    CM_row = self.CM[node_i, :]
    C_sum = 0
    indizes = CM_row[CM_row.!=0]
    signs = [sign(x) for x in indizes] # get signs
    indizes_ = indizes .* signs # delet signs from indices

    for idx in indizes_
        idx = Int(idx)
        C_sum += self.parameters["cable"][idx]["C"] * 0.5
    end

    return C_sum
end


"""
    GetActionIds(self::NodeConstructor)

Creates the State Vector for an related NodeConstructor and outputs it as a list of strings.
"""
function GetActionIds(self::NodeConstructor)
    actions = []

    for s in 1:self.num_sources
        if self.parameters["source"][s]["fltr"] == "LCL"
            push!(actions, "source$s" * "_u")
        elseif self.parameters["source"][s]["fltr"] == "LC"
            push!(actions, "source$s" * "_u")
        elseif self.parameters["source"][s]["fltr"] == "L"
            push!(actions, "source$s" * "_u")
        end
    end

    if self.parameters["grid"]["phase"] === 3
        A = ["_a", "_b", "_c"]
        actions = vcat([broadcast(*, actions, A[i]) for i in 1:3]...)
    end

    return actions
end


"""
    GetSourceStateIndices(self::NodeConstructor,sources)

Returns all state indices for passed sources.
"""
function GetSourceStateIndices(self::NodeConstructor, sources)
    state_ids = GetStateIds(self)
    action_ids = GetActionIds(self)
    source_indices = Dict()

    for idx in sources
        source = Dict()
        state_indices = findall(x -> occursin("source$idx" * "_", x), state_ids)
        action_indices = findall(x -> occursin("source$idx" * "_", x), action_ids)
        source["state_indices"] = state_indices
        source["action_indices"] = action_indices
        source_indices["source$idx"] = source
    end

    return source_indices
end


"""
    GetCableStateIndices(self::NodeConstructor,cables)

Returns all state indices for passed cables.
"""
function GetCableStateIndices(self::NodeConstructor, cables)
    state_ids = GetStateIds(self)
    cable_indices = Dict()

    for idx in cables
        cable = Dict()
        state_indices = findall(x -> occursin("cable$idx" * "_", x), state_ids)
        cable["state_indices"] = state_indices
        cable_indices["cable$idx"] = cable
    end

    return cable_indices
end


"""
    GetLoadStateIndices(self::NodeConstructor,loads)

Returns all state indices for passed loads.
"""
function GetLoadStateIndices(self::NodeConstructor, loads)
    state_ids = GetStateIds(self)
    load_indices = Dict()

    for idx in loads
        load = Dict()
        state_indices = findall(x -> occursin("load$idx" * "_", x), state_ids)
        load["state_indices"] = state_indices
        load_indices["load$idx"] = load
    end

    return load_indices
end

"""
    DrawGraph(self::NodeConstructor)

Plotting a graphical representation of the grid.
"""
function DrawGraph(self::NodeConstructor)
    CM = self.CM
    parameters = self.parameters
    CMtemp = CM + -2 * LowerTriangular(CM)
    G = SimpleGraph(CMtemp)

    # Position nodes
    pos_x, pos_y = GraphPlot.shell_layout(G)

    # Create plot points
    edge_x = []
    edge_y = []

    for edge in edges(G)
        push!(edge_x, pos_x[src(edge)])
        push!(edge_x, pos_x[dst(edge)])
        push!(edge_x, nothing)
        push!(edge_y, pos_y[src(edge)])
        push!(edge_y, pos_y[dst(edge)])
        push!(edge_y, nothing)
    end

    #  Color nodes
    color_map = []
    node_descriptions = []

    for source in parameters["source"]
        push!(node_descriptions, "Source: " * source["fltr"])

        if source["fltr"] == "LCL"
            push!(color_map, "#FF8800")
        elseif source["fltr"] == "LC"
            push!(color_map, "#FF6600")
        elseif source["fltr"] == "L"
            push!(color_map, "#FF3300")
        end
    end

    for load in parameters["load"]
        push!(node_descriptions, "Load: " * load["impedance"])

        if load["impedance"] == "RLC"
            push!(color_map, "#8F00D1")
        elseif load["impedance"] == "LC"
            push!(color_map, "#4900A8")
        elseif load["impedance"] == "RL"
            push!(color_map, "#3A09C0")
        elseif load["impedance"] == "RC"
            push!(color_map, "#0026FF")
        elseif load["impedance"] == "L"
            push!(color_map, "#0066FF")
        elseif load["impedance"] == "C"
            push!(color_map, "#00CCFF")
        elseif load["impedance"] == "R"
            push!(color_map, "#00F3E7")
        end
    end

    # Create edges
    edges_trace = scatter(
        mode="lines",
        x=edge_x,
        y=edge_y,
        line=attr(
            width=0.8,
            color="#113"
        ),
    )

    # Create nodes
    nodes_trace = scatter(
        x=pos_x,
        y=pos_y,
        mode="markers",
        text=node_descriptions,
        marker=attr(
            color=color_map,
            size=13,
            line=attr(
                color="Black",
                width=1
            )
        )
    )

    # Create Plot
    pl = PlotlyBase.Plot(
        [edges_trace, nodes_trace],
        PlotlyBase.Layout(
            plot_bgcolor="#f1f3f7",
            hovermode="closest",
            showlegend=false,
            showarrow=false,
            dragmode="select",
            xaxis=attr(showgrid=false, zeroline=false, showticklabels=false),
            yaxis=attr(showgrid=false, zeroline=false, showticklabels=false)
        )
    )

    display(pl)
end


"""
    GetYBus(self::NodeConstructor)

Returns the Admittance Matrix of the power grid based on the CM matrix and the parameters.
"""
function GetYBus(self::NodeConstructor)
    Y_bus = zeros(Complex{Float64}, self.tot_ele, self.tot_ele)  # Y -> tot_ele x tot_ele
    omega = 2 * π * self.parameters["grid"]["f_grid"]

    for col in 1:self.tot_ele
        for row in 1:self.tot_ele
            if self.CM[row, col] != 0  # we have a cable connected
                # CM index defines the number of the cable
                cable_idx = abs(Int(self.CM[row, col]))
                G_RL = self.parameters["cable"][cable_idx]["R"] /
                       (self.parameters["cable"][cable_idx]["R"]^2 +
                        omega^2 * self.parameters["cable"][1]["L"]^2)
                B_RL = (omega * self.parameters["cable"][cable_idx]["L"]) /
                       (self.parameters["cable"][cable_idx]["R"]^2 +
                        omega^2 * self.parameters["cable"][1]["L"]^2)
                Y_bus[row, col] = -G_RL - im * B_RL
            elseif row == col  # diagonal elements
                # Go through all col elements of that row to find the connected cable to
                # that bus (non zero elements in CM[:, row])
                cable_idxs = filter(n -> n != 0, self.CM[row, :])
                G = 0
                B = 0

                for idx in cable_idxs
                    # add all RL
                    idx = abs(Int(idx))
                    G += self.parameters["cable"][idx]["R"] /
                         (self.parameters["cable"][idx]["R"]^2 +
                          omega^2 * self.parameters["cable"][1]["L"]^2)
                    B += (omega * self.parameters["cable"][idx]["L"]) /
                         (self.parameters["cable"][idx]["R"]^2 +
                          omega^2 * self.parameters["cable"][1]["L"]^2)
                    # and add all shunt C connected to that bus since diagonal element
                    B += omega * self.parameters["cable"][idx]["C"] / 2
                end

                Y_bus[row, col] = G + im * B
            end
        end
    end

    return Y_bus
end


"""
under development
"""
function SourceSetup(num_sources; random=nothing, awg_pwr=200e3, mode=3)
    #= Modes:
        1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
        2 -> "PQ" - grid following controllable source/load (active and reactive Power)
        3 -> "Droop" - simple grid forming with power balancing
        4 -> "Synchronverter" - enhanced droop control
    =#

    source_list = []
    total_gen = 0.0

    if random != 0 && !isnothing(random)
        Random.seed!(1234)
        pwrs = rand(Uniform(0.5 * awg_pwr, 1.5 * awg_pwr), num_sources)
    end

    for i in 1:num_sources
        source = Dict()

        if random == 0 || isnothing(random)
            source["mode"] = mode
            source["fltr"] = "LCL"  # Filter type
            pwr = awg_pwr
            source["pwr"] = pwr # Rated Apparent Power, VA
            source["p_set"] = 0   # Real Power Set Point, Watt
            source["q_set"] = 0   # Imaginary Power Set Point, VAi
            source["v_pu_set"] = 1.00   # Voltage Set Point, p.u.
            source["v_δ_set"] = 0      # Voltage Angle, degrees
            source["τv"] = 0.002  # Time constant of the voltage loop, seconds
            source["τf"] = 0.002  # Time constant of the frequency loop, seconds
            source["Observer"] = true   # Discrete Luenberger Observer
        else
            source["mode"] = mode
            source["fltr"] = "LCL"  # Filter type
            pwr = pwrs[i]
            source["pwr"] = pwr # Rated Apparent Power, VA
            source["p_set"] = 0 # Real Power Set Point, Watt
            source["q_set"] = 0 # Imaginary Power Set Point, VAi
            source["τv"] = 0.002 # Time constant of the voltage loop, seconds
            source["τf"] = 0.002 # Time constant of the frequency loop, seconds
            source["Observer"] = true # Discrete Luenberger Observer
            source["v_pu_set"] = 1.00 # Voltage Set Point, p.u.
            source["v_δ_set"] = 0  # Voltage Angle, degrees

            if random == 2
                source["std_asy"] = pwr / 8 # Asymptotic Standard Deviation
                source["σ"] = pwr / 8 # Brownian motion scale
                source["Δt"] = 0.1 # Time Step, seconds
                # source["X₀"] = 0 # Initial Process Values, Watt
                source["k"] = 1 # Interpolation degree
                source["γ"] = 0 # asymptotic mean
            end
        end

        total_gen += pwr
        push!(source_list, source)
    end

    return source_list, total_gen
end


"""
    LoadSetup(num_loads, total_gen; gen_load_ratio=6, random=nothing, Vrms=230)

    under development
"""
function LoadSetup(num_loads, total_gen; gen_load_ratio=6, random=nothing, Vrms=230)
    load_list = []
    avg_load = total_gen / (num_loads * gen_load_ratio)

    if random != 0 && !isnothing(random)
        Random.seed!(1234)
        pwrs = rand(Uniform(0.5 * avg_load, 1.5 * avg_load), num_loads)
        Random.seed!(1234)
        pfs = rand(Uniform(0.7, 1.0), num_loads)
    end

    for i in 1:num_loads
        load = Dict()

        if random == 0 || isnothing(random)
            R_load, L_load, _, _ = ParallelLoadImpedance(avg_load, 0.8, Vrms)
            load["impedance"] = "RL"
            load["R"] = R_load
            load["L"] = L_load
            load["S"] = avg_load
        else
            R_load, L_load, _, _ = ParallelLoadImpedance(pwrs[i], pfs[i], Vrms)
            load["impedance"] = "RL"
            load["R"] = R_load
            load["L"] = L_load
            load["S"] = pwrs[i]
        end

        push!(load_list, load)
    end

    return load_list
end


"""
    CableLengthSetup(num_cables; random=0, length_bounds=[0.5; 1.5])

    under development
"""
function CableLengthSetup(num_cables; random=0, length_bounds=[0.5; 1.5])
    cable_list = []

    if random != 0 && !isnothing(random)
        Random.seed!(1234)
        lengths = rand(Uniform(length_bounds[1], length_bounds[2]), num_cables)
    end

    for i in 1:num_cables
        cable = Dict()

        if random == 0
            cable["len"] = sum(length_bounds) / 2   # km
            cable["R"] = 0.208   # Ω, line resistance
            cable["L"] = 0.00025 # H, line inductance
            cable["C"] = 0.4e-3  # F, line capacitance
            cable["i_limit"] = 10e12   # A, line current limit
        else
            cable["len"] = lengths[i]
        end

        push!(cable_list, cable)
    end

    return cable_list
end


"""
    SetupMG(num_sources, num_cables; random=nothing, avg_pwr=200e3, Vrms=230)

under development
"""
function SetupMG(num_sources, num_cables; random=nothing, avg_pwr=200e3, Vrms=230)

    #= Modes:
        1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
        2 -> "PQ" - grid following controllable source/load (active and reactive Power)
        3 -> "Droop" - simple grid forming with power balancing
        4 -> "Synchronverter" - enhanced droop control
    =#

    gen_static_load_ratio = 6
    gen_contr_load_ratio = 3

    source_list = []
    load_list = []
    cable_list = []

    num_contr_loads = num_sources
    num_static_loads = num_sources
    num_loads = num_contr_loads + num_static_loads

    total_gen = 0.0

    if random != 0 && !isnothing(random)

        Random.seed!(1234)
        pwrs = rand(Uniform(0.5 * avg_pwr, 1.5 * avg_pwr), num_sources)
    end

    # Grid Forming sources
    for i in 1:num_sources
        source = Dict()
        if random == 0 || isnothing(random)

            # Grid Forming sources
            source["mode"] = "Synchronverter"
            source["fltr"] = "LCL"  # Filter type
            pwr = avg_pwr
            source["pwr"] = pwr # Rated Apparent Power, VA
            source["p_set"] = 0   # Real Power Set Point, Watt
            source["q_set"] = 0   # Imaginary Power Set Point, VAi
            source["v_pu_set"] = 1.00   # Voltage Set Point, p.u.
            source["Observer"] = true   # Discrete Luenberger Observer
            source["τv"] = 0.002  # Time constant of the voltage loop, seconds
            source["τf"] = 0.002  # Time constant of the frequency loop, seconds
        else
            source["mode"] = "Synchronverter"
            source["fltr"] = "LCL"  # Filter type
            pwr = pwrs[i]
            source["pwr"] = pwr # Rated Apparent Power, VA
            source["p_set"] = 0 # Real Power Set Point, Watt
            source["q_set"] = 0 # Imaginary Power Set Point, VAi
            source["τv"] = 0.002 # Time constant of the voltage loop, seconds
            source["τf"] = 0.002 # Time constant of the frequency loop, seconds
            source["Observer"] = true # Discrete Luenberger Observer
            source["v_pu_set"] = 1.00 # Voltage Set Point, p.u.
        end

        total_gen += pwr
        push!(source_list, source)
    end

    # Static loads
    avg_static_load = total_gen / (num_static_loads * gen_static_load_ratio)

    for i in 1:num_static_loads

        load = Dict()
        R_load, L_load, _, _ = ParallelLoadImpedance(avg_static_load, 1.0, Vrms)
        load["impedance"] = "R"
        load["R"] = R_load
        load["L"] = L_load
        load["S"] = avg_static_load

        push!(load_list, load)
    end

    # Cables
    MG_cables = num_sources
    inter_cables = num_cables - MG_cables

    for i in 1:num_cables
        cable = Dict()

        if i <= inter_cables
            l = 1
        else
            l = 0.1
        end

        cable["len"] = l
        cable["R"] = l * 0.722  # Ω, line resistance
        cable["L"] = l * 0.264e-3# H, line inductance
        cable["C"] = l * 0.4e-6  # F, line capacitance
        cable["i_limit"] = 10e12   # A, line current limit

        #= cable["Rb"] =  0.722
        cable["Cb"] = 0.4e-6
        cable["Lb"] =  0.264e-3 =#

        push!(cable_list, cable)
    end

    # parameters
    parameters = Dict{Any,Any}()
    parameters["source"] = source_list
    parameters["load"] = load_list
    parameters["cable"] = cable_list

    return parameters
end

"""
default source settings
"""
function SampleSource(grid_properties, fltr)
    source = Dict()
    source["source_type"] = "ideal"
    source["fltr"] = fltr

    source["pwr"] = rand(range(start=5, step=5, stop=50)) * 1e3
    source["vdc"] = 800
    source["i_rip"] = 0.15
    source["v_rip"] = 0.01537

    source["τv"] = 0.002 # time constant of the voltage loop # 0.02
    source["τf"] = 0.002 # time constant of the frequency loop # 0.002
    source["pf"] = 0.95 # power factor
    source["p_set"] = 0#source["pwr"] * source["pf"]
    source["q_set"] = 0#sqrt(source["pwr"]^2 - source["p_set"]^2)
    source["v_pu_set"] = 1.0
    source["v_δ_set"] = 0.0
    source["mode"] = "Synchronverter"
    source["control_type"] = "classic"
    source["γ"] = source["p_set"]
    source["std_asy"] = source["pwr"] / 4
    source["σ"] = 0.0
    source["κ"] = source["σ"]^2 / (2 * source["std_asy"]^2)
    source["X₀"] = source["p_set"]
    source["Δt"] = round(grid_properties["fs"] / (grid_properties["f_grid"])) / grid_properties["fs"]
    source["k"] = 0

    return source
end
