"""
Basic example for a diffrent source
"""
Base.@kwdef mutable struct DCLinkVoltages
    sources
end;

Base.@kwdef mutable struct IdealVoltage
    num
    state_i_idx
    state_v_idx
    V
    P
end

Base.@kwdef mutable struct PVVoltage
    num
    state_i_idx
    state_v_idx
    SolarArray
end

Base.@kwdef mutable struct BatteryVoltage
    num
    state_i_idx
    state_v_idx
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

        state_i = "source$(num)_i_L1"
        state_i_idx = findall(x -> occursin(state_i, x), states)

        if source["fltr"] == "L"
            state_v = "source$(num)_v_C_cables"
        else
            state_v = "source$(num)_v_RC"
        end

        state_v_idx = findall(x -> occursin(state_v, x), states)


        if source["source_type"] == "ideal"
            push!(v_dc_vector, IdealVoltage(num, state_i_idx, state_v_idx, source["vdc"], source["pwr"]))
        elseif source["source_type"] == "pv"
            #TODO : how to calculate i_dc in 3-phase grid? Which current of env to use?
            ModulePV = solar_module(I_0=source["module_I_0"], ni=source["module_ni"],
                N_cell=source["module_N_cell"], I_ph_ref=source["module_I_ph_ref"])
            SolarArr = solar_array(; solar_module=ModulePV, serial=source["serial"],
                parallel=source["parallel"])

            push!(v_dc_vector, PVVoltage(num, state_i_idx, state_v_idx, SolarArr))
        elseif source["source_type"] == "battery"
            #TODO : how to calculate i_dc in 3-phase grid? Which current of env to use?
            Cell = battery_module(source["module_R"], source["module_C"])
            Battery = battery_block(battery_module=Cell, R_0=source["R_0"],
                V_0=source["V_0"], Q_0=source["Q_0"], Q=source["Q"],
                SOC_BP=source["SOC_BP"], T_BP=source["T_BP"],
                tau=1/nc.parameters["grid"]["fs"], i_bat_limit=source["i_bat_limit"],
                serial = source["serial"], parallel = source["parallel"])

            push!(v_dc_vector, BatteryVoltage(num, state_i_idx, state_v_idx, Battery))
        else
            @warn "sourceType not known! vdc set to fixed value"
            push!(v_dc_vector, IdealVoltage(800, num, 10e3))
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
            if env.nc.parameters["grid"]["phase"] == 3

                V_αβγ = ClarkeTransform(env.x[source.state_v_idx])
                I_αβγ = ClarkeTransform(env.x[source.state_i_idx])

                pq0 = [I_αβγ[1] I_αβγ[2] 0; -I_αβγ[2] I_αβγ[1] 0; 0 0 I_αβγ[3]]*V_αβγ
                # ac case
            else
                throw("unimplemented")
            end

            if source.SolarArray.v_next == 0
                I_dc = 0.
            else
                I_dc = pq0[1]/source.SolarArray.v_next
            end

            push!(V, get_V(source.SolarArray, I_dc, env.nc.parameters["weather"]["G"], env.nc.parameters["weather"]["T"]))
        elseif typeof(source) == BatteryVoltage
            if env.nc.parameters["grid"]["phase"] == 3
                # ac case
                V_αβγ = ClarkeTransform(env.x[source.state_v_idx])
                I_αβγ = ClarkeTransform(env.x[source.state_i_idx])

                pq0 = [I_αβγ[1] I_αβγ[2] 0; -I_αβγ[2] I_αβγ[1] 0; 0 0 I_αβγ[3]]*V_αβγ
            else
                throw("unimplemented")
            end

            if source.BatteryBlock.v_next == 0
                I_dc = 0.
            else
                I_dc = pq0[1]/source.BatteryBlock.v_next
            end

            push!(V, get_V(source.BatteryBlock, I_dc, env.nc.parameters["weather"]["T"]))
        else
            @assert false "Type of source '$typeof(source)' is not valid"
        end

    end

    return V
end
