# MUST disable threading in Qt
#ENV["QSG_RENDER_LOOP"] = "basic"

using QML

# path to the QML file
qml_file = joinpath(dirname(@__FILE__), "qml", "ElectricGridGUI.qml")


sources = []
loads = []
cables = []


function addSource(uid, source_type, control_type, mode, filter, pwr)

    push!(sources, Dict(
        "uid" => uid,
        "source_type" => source_type,
        "control_type" => control_type,
        "mode" => mode,
        "filter" => filter,
        "pwr" => pwr,
    ))

    updated()
end

function addLoad(uid, impedance, pwr)

    push!(loads, Dict(
        "uid" => uid,
        "impedance" => impedance,
        "pwr" => pwr,
    ))

    updated()
end

function addCable(uid, from, to, length, capacity, inductance, resistance)
    push!(cables, Dict(
        "uid" => uid,
        "from" => from,
        "to" => to,
        "length" => length,
        "capacity" => capacity,
        "inductance" => inductance,
        "resistance" => resistance,
    ))

    updated()
end

function updateSource(uid, source_type, control_type, mode, filter, pwr)
    index = findall(x -> x["uid"] == uid, sources)[1]
    
    sources[index] = Dict(
        "uid" => uid,
        "source_type" => source_type,
        "control_type" => control_type,
        "mode" => mode,
        "filter" => filter,
        "pwr" => pwr,
    )

    updated()
end

function updateLoad(uid, impedance, pwr)
    index = findall(x -> x["uid"] == uid, loads)[1]
    
    loads[index] = Dict(
        "uid" => uid,
        "impedance" => impedance,
        "pwr" => pwr,
    )

    updated()
end

function updateCable(uid, from, to, length, capacity, inductance, resistance)
    index = findall(x -> x["uid"] == uid, cables)[1]
    
    cables[index] = Dict(
        "uid" => uid,
        "from" => from,
        "to" => to,
        "length" => length,
        "capacity" => capacity,
        "inductance" => inductance,
        "resistance" => resistance,
    )

    updated()
end

function deleteNode(uid)
    if startswith(uid, "source")
        deleteSource(uid)
    else
        deleteLoad(uid)
    end
end

function deleteSource(uid)
    index = findall(x -> x["uid"] == uid, sources)[1]
    
    deleteat!(sources, index)

    cables_to_delete = []
    for (i, cable) in enumerate(cables)
        if uid == cable["from"] || uid == cable["to"]
            push!(cables_to_delete, i)
        end
    end
    deleteat!(cables, cables_to_delete)

    updated()
end

function deleteLoad(uid)
    index = findall(x -> x["uid"] == uid, loads)[1]
    
    deleteat!(loads, index)

    cables_to_delete = []
    for (i, cable) in enumerate(cables)
        if uid == cable["from"] || uid == cable["to"]
            push!(cables_to_delete, i)
        end
    end
    deleteat!(cables, cables_to_delete)

    updated()
end

function deleteCable(uid)
    index = findall(x -> x["uid"] == uid, cables)[1]
    
    deleteat!(cables, index)

    updated()
end

function updated()
    global CM
    global parameters

    CM = zeros(length(sources) + length(loads), length(sources) + length(loads))

    cableindex = 1
    for cable in cables
        from = startswith(cable["from"], "source") ? findall(x -> x["uid"] == cable["from"], sources)[1] : length(sources) + findall(x -> x["uid"] == cable["from"], loads)[1]
        to = startswith(cable["to"], "source") ? findall(x -> x["uid"] == cable["to"], sources)[1] : length(sources) + findall(x -> x["uid"] == cable["to"], loads)[1]

        CM[from, to] = from < to ? Float64(cableindex) : -1 * Float64(cableindex)
        CM[to, from] = from < to ? -1 * Float64(cableindex) : Float64(cableindex)

        cableindex += 1
    end

    parameters = Dict{Any, Any}(
        "source" => Any[],
        "load"   => Any[],
        "cable"   => Any[]
    )

    for source in sources
        if source["control_type"] == "Classic"
            push!(parameters["source"], Dict{Any, Any}("source_type" => source["source_type"],
                                                        "control_type" => source["control_type"],
                                                        "mode" => source["mode"],
                                                        "fltr" => source["filter"],
                                                        "pwr" => source["pwr"]))
        else
            push!(parameters["source"], Dict{Any, Any}("control_type" => source["control_type"],
                                                        "fltr" => source["filter"],
                                                        "pwr" => source["pwr"]))
        end
    end

    for load in loads
        push!(parameters["load"], Dict{Any, Any}("impedance" => load["impedance"],
                                                    "pwr" => load["pwr"]))
    end

    for cable in cables
        push!(parameters["cable"], Dict{Any, Any}("len" => cable["length"],
                                                    "R" => cable["resistance"] * cable["length"],
                                                    "L" => cable["inductance"] * cable["length"],
                                                    "C" => cable["capacity"] * cable["length"]))
    end

    #println(CM)
    #println(parameters)
end

@qmlfunction addSource addLoad addCable updateSource updateLoad updateCable deleteNode deleteSource deleteLoad deleteCable


loadqml(qml_file)

# Start the GUI
#exec_async()
exec()