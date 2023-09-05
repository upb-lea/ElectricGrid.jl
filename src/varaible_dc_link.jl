"""
Basic example for a diffrent source
"""
Base.@kwdef mutable struct DCLinkVoltages
    sources
end;

Base.@kwdef mutable struct IdealVoltage
    num
    state_idx
    V
    P
end

Base.@kwdef mutable struct PVVoltage
    num
    state_idx
    SolarArray
end

Base.@kwdef mutable struct BatteryVoltage
    num
    state_idx
    BatteryBlock
end

"""
Initializes the VDC link struct by extracting the structs from the parameterdict and
collecting them in a separate struct which will be queried later at runtime.
"""
function DCLinkVoltagesInit(nc)
    source_list = nc.parameters["source"]
    v_dc_vector = []
    vdc_fixed = 0

    states = GetStateIds(nc)

    for (num, source) in enumerate(source_list)

        @assert haskey(source, "source_type") "Source $(num) has no sourceType defined!"

        state = "source$(num)_i_L1"
        state_idx = findall(x -> occursin(state, x), states)

        if source["source_type"] == "ideal"
            push!(v_dc_vector, IdealVoltage(num, state_idx, source["vdc"], source["pwr"]))
        elseif source["source_type"] == "pv"
            #TODO : how to calculate i_dc in 3-phase grid? Which current of env to use?
            ModulePV = solar_module(I_0=source["module_I_0"], ni=source["module_ni"],
                N_cell=source["module_N_cell"], I_ph_ref=source["module_I_ph_ref"])
            SolarArr = solar_array(; solar_module=ModulePV, serial=source["serial"],
                parallel=source["parallel"])

            push!(v_dc_vector, PVVoltage(num, state_idx, SolarArr))
        elseif source["source_type"] == "battery"
            #TODO : how to calculate i_dc in 3-phase grid? Which current of env to use?
            Cell = battery_module(source["module_R"], source["module_C"])
            Battery = battery_block(battery_module=Cell, n=source["n"], R_0=source["R_0"],
                V_0=source["V_0"], Q_0=source["Q_0"], Q=source["Q"],
                SOC_BP=source["SOC_BP"], T_BP=source["T_BP"],
                tau=nc.parameters["grid"]["fs"])

            push!(v_dc_vector, BatteryVoltage(num, state_idx, Battery))
        else
            @warn "sourceType not known! vdc set to fixed value"
            push!(v_dc_vector, IdealVoltage(700, num, 10e3))
            vdc_fixed += 1
        end
    end
    vdc_fixed > 0 && @warn("$vdc_fixed DC-link voltages set to 800 V - please define in
                            nc.parameters -> source -> vdc.")

    DCLinkVoltages(sources=v_dc_vector)
end

function step(DCLinkVoltages, env)

    V = []

    for source in DCLinkVoltages.sources
        if typeof(source) == IdealVoltage

            push!(V, source.V)
        elseif typeof(source) == PVVoltage
            I_ac = env.state[source.state_idx][1]

            println("I: $I_ac")

            # P_ac = I_ac * V_ac
            # I_dc = P_ac / V_dc

            I_dc = I_ac

            println("V: $(get_V(source.SolarArray, I_ac, env.nc.parameters["weather"]["G"], env.nc.parameters["weather"]["T"]))")
            push!(V, get_V(source.SolarArray, I_ac, env.nc.parameters["weather"]["G"], env.nc.parameters["weather"]["T"]))
        elseif typeof(source) == BatteryVoltage
            I = env.state[source.state_idx][1]

            println("I: $I")
            println("V: $(get_V(source.BatteryBlock, I, env.nc.parameters["weather"]["T"]))")
            push!(V, get_V(source.BatteryBlock, I, env.nc.parameters["weather"]["T"]))
        else
            @assert false "Type of source '$typeof(source)' is not valid"
        end

    end

    return V
end
