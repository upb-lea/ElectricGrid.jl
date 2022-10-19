using Distributions
using LinearAlgebra
using StatsBase
# using Intervals
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
end

"""
    NodeConstructor(;
        num_sources,
        num_loads,
        CM = nothing,
        parameters = nothing,
        S2S_p = 0.1,
        S2L_p = 0.8
        )

Create a mutable struct NodeConstructor, which serves as a basis for the creation of an energy grid: `num_sources` corresponse to the amount of sources and `num_loads` is the amount of loads in the grid. `CM` is the connection matrix which indicates how the elements in the grid are connected to each other. To specify the elements of the net in more detail, values for the elements can be passed via `parameters`. If no connection matrix is entered, it can be generated automatically. `S2S_p` is the probability that a source is connected to another source and `S2L_p` is the probability that a source is connected to a load.
"""
function NodeConstructor(;num_sources, num_loads, CM=nothing, parameters=nothing, S2S_p=0.1, S2L_p=0.8)

    tot_ele = num_sources + num_loads

    cntr = 0
    num_connections = 0

    if CM === nothing
        cntr, CM = CM_generate(tot_ele, num_sources, S2L_p, S2S_p)
        num_connections = cntr
    else
        if size(CM)[1] != tot_ele
            throw("Expect the number of elements in the node to match the specified structure in the CM, but got $tot_ele and $(size(A)[1])")
        end
        num_connections = Int(maximum(CM))
    end

    if parameters === "random"
        
        num_fltr_L, num_fltr_LC, num_fltr_LCL = get_fltr_distr(num_sources)
        
        loads_distr = get_loads_distr(num)
        num_loads_R = loads_distr[1]
        num_loads_C = loads_distr[2]
        num_loads_L = loads_distr[3]
        num_loads_RL = loads_distr[4]
        num_loads_RC = loads_distr[5]
        num_loads_LC = loads_distr[6]
        num_loads_RLC  = loads_distr[7]
        

        parameters = generate_parameters(num_fltr_LCL, num_fltr_LC, num_fltr_L, num_connections, num_loads_RLC, num_loads_LC, num_loads_RL, num_loads_RC,
                                        num_loads_L, num_loads_C, num_loads_R)
    end

    if parameters === nothing || isa(parameters, Dict) 

        parameters = check_parameters(parameters, num_sources, num_loads, num_connections)  # Checks if all entries are given, if not, fills up with random values

        @assert length(keys(parameters)) == 4 "Expect parameters to have the four entries 'cable', 'load', 'grid' and 'source' but got $(keys(parameters))"

        @assert length(keys(parameters["grid"])) == 3 "Expect parameters['grid'] to have the three entries 'fs', 'v_rms' and 'phase' but got $(keys(parameters))"

        @assert length(parameters["source"]) == num_sources "Expect the number of sources to match the number of sources in the parameters, but got $num_sources and $(length(parameters["source"]))"

        @assert length(parameters["load"]) == num_loads "Expect the number of loads to match the number of sources in the parameters, but got $num_loads and $(length(parameters["load"]))"

        @assert length(parameters["cable"]) == num_connections "Expect the number of sources to match the number of sources in the parameters, but got $num_connections and $(length(parameters["cable"]))"

        # println(parameters)

        num_fltr_LCL, num_fltr_LC, num_fltr_L = cntr_fltrs(parameters["source"])
        num_loads_RLC, num_loads_LC, num_loads_RL, num_loads_RC, num_loads_L, num_loads_C, num_loads_R = cntr_loads(parameters["load"])
        

        @assert num_fltr_LCL + num_fltr_LC + num_fltr_L == num_sources "Expect the number of sources to be identical to the sum of the filter types, but the number of sources is $num_sources and the sum of the filters is $(num_fltr_LCL + num_fltr_LC + num_fltr_L)"

        @assert num_loads_RLC + num_loads_RL + num_loads_RC + num_loads_R == num_loads "Expect the number of loads to be identical to the sum of the loads types, but the number of loads is $num_loads and the sum of the loads is $(num_loads_RLC + num_loads_RL + num_loads_RC + num_loads_R)"
        
        

        # valid_realistic_para(parameters) TODO: revise that values are checked independently of network parameters, if necessary only via A-matrix

    else
        throw("Expect parameters to be a dict or nothing, not $(typeof(parameters))")
    end

    num_fltr = 4 * num_fltr_LCL + 3 * num_fltr_LC + 2 * num_fltr_L
    num_impedance = (2 * (num_loads_RLC
                        + num_loads_LC
                        + num_loads_RL
                        + num_loads_L)
                        + num_loads_RC + num_loads_C + num_loads_R)

    num_spp = num_fltr_LCL * 4 + num_fltr_LC * 3 + num_fltr_L * 2 + num_connections + (num_loads_RLC + num_loads_LC + num_loads_RL + num_loads_L) * 2 + (num_loads_RC + num_loads_C + num_loads_R)

    NodeConstructor(num_connections, num_sources, num_loads, num_fltr_LCL, num_fltr_LC, num_fltr_L, num_loads_RLC, num_loads_LC, num_loads_RL, num_loads_RC, num_loads_L, num_loads_C, num_loads_R, num_impedance, num_fltr, num_spp, cntr, tot_ele, CM, parameters, S2S_p, S2L_p)
end

function get_fltr_distr(num)
    
    # a = [.49 .02 .49] # probability for LC should be lower

    di_di= Dirichlet(ones(3))
    smpl = rand(di_di, 1) * num
    num_fltr_L = Int(floor(smpl[1]))
    num_fltr_LC = Int(ceil(clamp(smpl[2], 1, num-1)))
    num_fltr_LCL =  num - (num_fltr_LC + num_fltr_L)

    return num_fltr_L, num_fltr_LC, num_fltr_LCL
end

function get_loads_distr(num)
    di_di= Dirichlet(ones(7)) # create dirichlet distribution
    smpl = rand(di_di, 1) * num
    
    num_loads_R = Int(floor(smpl[1]))
    num_loads_C = Int(floor(smpl[2]))
    num_loads_L = Int(floor(smpl[3]))
    num_loads_RL = Int(floor(smpl[4]))
    num_loads_RC = Int(floor(smpl[5]))
    num_loads_LC = Int(floor(smpl[6]))
    num_loads_RLC =  num - (num_loads_R + num_loads_C + num_loads_L + num_loads_RL + num_loads_RC + num_loads_LC)

    return num_loads_R, num_loads_C, num_loads_L, num_loads_RL, num_loads_RC, num_loads_LC, num_loads_RLC
end

function check_parameters(parameters, num_sources, num_loads, num_connections)

    # Variable generation of the parameter dicts

    # check if grid parameters have been specified

    if parameters === nothing
        parameters = Dict()
    end

    ##############
    # CHECK GRID #
    ##############
    if !haskey(parameters, "grid") 
        grid_properties = Dict()
        grid_properties["fs"] =  10e3
        grid_properties["v_rms"] = 230
        grid_properties["phase"] = 3
        parameters["grid"] = grid_properties

    else
        if !haskey(parameters["grid"], "fs")
            parameters["grid"]["fs"] = 10e3
        end
        if !haskey(parameters["grid"], "v_rms")
            parameters["grid"]["v_rms"] = 230
        end
        if !haskey(parameters["grid"], "phase")
            parameters["grid"]["phase"] = 3
        end
    end

    ################
    # CHECK SOURCE #
    ################

    if !haskey(parameters, "source")

        num_fltr_L, num_fltr_LC, num_fltr_LCL = get_fltr_distr(num_sources)

        source_list = []

        for s in 1:num_fltr_LCL
        push!(source_list, _sample_fltr_LCL(parameters["grid"]))
        end
    
        for s in 1:num_fltr_LC
        push!(source_list, _sample_fltr_LC(parameters["grid"]))
        end
    
        for s in 1:num_fltr_L
        push!(source_list, _sample_fltr_L(parameters["grid"]))
        end

        parameters["source"] = source_list
    
    else
        num_def_sources = length(parameters["source"])

        num_undef_sources = num_sources - num_def_sources

        @assert num_undef_sources >= 0 "Expect the number of defined sources within the parameter dict to be less or equal to the number of sources in the env, but the entries within the parameter dict is $num_def_sources and the number of env sources is $num_sources."

        if num_undef_sources > 0
            println("WARNING: The number of defined sources $num_def_sources is smaller than the number specified sources in the environment $num_sources, therefore the remaining $num_undef_sources sources are selected randomly!")
        end

        num_LC_defined = 0
        source_type_fixed = 0

        for (index, source) in enumerate(parameters["source"])
            if !haskey(source, "pwr")
                source["pwr"] = rand(range(start=5,step=5,stop=50))*1e3
            end
            
            if !haskey(source, "vdc")
                source["vdc"] = rand(range(start=690,step=10,stop=800))
            end

            if !haskey(source, "i_rip")
                source["i_rip"] = rand(Uniform(0.1, 0.15))
            end
                
            if !haskey(source, "v_rip")                    
                source["v_rip"] = rand(Uniform(0.014, 0.016))
            end
        
            #Inductor design
            if !haskey(source, "L1") 
                Vorms = parameters["grid"]["v_rms"]*1.05
                Vop = Vorms*sqrt(2)
            
                Zl = 3*Vorms*Vorms/source["pwr"]
            
                Iorms = Vorms/Zl
                Iop = Iorms*sqrt(2)
            
                ΔIlfmax = source["i_rip"]*Iop
            
                source["L1"] = 0.5*(source["vdc"]*(4*parameters["grid"]["fs"]*ΔIlfmax)^-1)
                
            end

            if !haskey(source, "i_limit")
                Vorms = parameters["grid"]["v_rms"]*1.05
                Vop = Vorms*sqrt(2)
            
                Zl = 3*Vorms*Vorms/source["pwr"]
            
                Iorms = Vorms/Zl
                Iop = Iorms*sqrt(2)
            
                ΔIlfmax = source["i_rip"]*Iop
                source["i_limit"]= Iop+ ΔIlfmax
            end

            if !haskey(source, "R1")
                #TODO: @Septimus: why?
                source["R1"] = 400 * source["L1"]
            end
            
            if !haskey(source, "fltr")
                source["fltr"] = rand(["L", "LC", "LCL"])
            elseif !(source["fltr"] in ["L", "LC", "LCL"])
                # TODO: Raise warning: False key
                source["fltr"] = "L"
                println("WARNING: filterType not known! set to L filter, please choose L, LC, or LCL!")
            end
            

            if (source["fltr"] == "LC" || source["fltr"] == "LCL") && !haskey(source, "C1")
                #Capacitor design
                if source["fltr"] == "LC"
                    num_LC_defined += 1
                end
                Vorms = parameters["grid"]["v_rms"]*0.95
                Vop = Vorms*sqrt(2)

                Zl = 3*Vorms*Vorms/source["pwr"]

                Iorms = Vorms/Zl
                Iop = Iorms*sqrt(2)
                Ir_d = source["vdc"]/(4*parameters["grid"]["fs"]*source["L1"]*Iop)
                ΔIlfmax = Ir_d*Iop
                ΔVcfmax = source["v_rip"]*Vop

                source["C"] = ΔIlfmax/(8*parameters["grid"]["fs"]*ΔVcfmax)
            end

            if !haskey(source, "v_limit")
                Vorms = parameters["grid"]["v_rms"]*0.95
                Vop = Vorms*sqrt(2)

                Zl = 3*Vorms*Vorms/source["pwr"]

                Iorms = Vorms/Zl
                Iop = Iorms*sqrt(2)
                Ir_d = source["vdc"]/(4*parameters["grid"]["fs"]*source["L1"]*Iop)
                ΔIlfmax = Ir_d*Iop
                ΔVcfmax = source["v_rip"]*Vop
                
                source["v_limit"]= Vop+ΔVcfmax
            end

            if  source["fltr"] == "LCL" && !haskey(source, "L2")
                if !haskey(source, "L2")
                    source["L2"] = deepcopy(source["L1"])
                end
                
                if !haskey(source, "R2")
                    source["R2"] = deepcopy(source["R1"]) 
                end  

                if !haskey(source, "R_C")
                    source["R_C"] = 28*source["C"] 
                end  
            end

            if !haskey(source, "source_type")
                source["source_type"] = "ideal"
                source_type_fixed += 1
            end
        end

        if num_undef_sources > 0
        
            num_fltr_L_undef, num_fltr_LC_undef, num_fltr_LCL_undef = get_fltr_distr(num_undef_sources)

            for s in 1:num_fltr_LCL_undef
            push!(parameters["source"], _sample_fltr_LCL(parameters["grid"]))
            end
        
            for s in 1:num_fltr_LC_undef
            push!(parameters["source"], _sample_fltr_LC(parameters["grid"]))
            end
        
            for s in 1:num_fltr_L_undef
            push!(parameters["source"], _sample_fltr_L(parameters["grid"]))
            end
        end

        # Validierung ob LC vorhanden ist?
        if num_LC_defined == 0 &&  num_fltr_LC_undef == 0
            println("WARNING: No LC filter defined/set random, if wanted please set in parameter dict!")
        end

        source_type_fixed > 0 && println("WARNING: $vdc_fixed sourceType not defined! set to ideal!")

        num_fltr_LCL, num_fltr_LC, num_fltr_L = cntr_fltrs(parameters["source"])
    
    end

    ###############
    # CHECK LOADS #
    ###############

    if !haskey(parameters, "load")
        loads_distr = get_loads_distr(num_loads)
        num_loads_R = loads_distr[1]
        num_loads_C = loads_distr[2]
        num_loads_L = loads_distr[3]
        num_loads_RL = loads_distr[4]
        num_loads_RC = loads_distr[5]
        num_loads_LC = loads_distr[6]
        num_loads_RLC  = loads_distr[7]

        load_list = []

        for l in 1:num_loads_RLC
        push!(load_list, _sample_load_RLC())
        end

        for l in 1:num_loads_LC
        push!(load_list, _sample_load_LC())
        end

        for l in 1:num_loads_RL
        push!(load_list, _sample_load_RL())
        end

        for l in 1:num_loads_L
        push!(load_list, _sample_load_L())
        end

        for l in 1:num_loads_RC
        push!(load_list, _sample_load_RC())
        end

        for l in 1:num_loads_C
        push!(load_list, _sample_load_C())
        end

        for l in 1:num_loads_R
        push!(load_list, _sample_load_R())
        end
        parameters["load"] = load_list

    else
        num_def_loads = length(parameters["loads"])

        num_undef_loads = num_loads - num_def_loads

        @assert num_undef_loads >= 0 "Expect the number of defined loads within the parameter dict to be less or equal to the number of loads in the env, but the entries within the parameter dict is $num_def_loads and the number of env loads is $num_loads."

        if num_undef_loads > 0
            println("WARNING: The number of defined loads $num_def_loads is smaller than the number specified loads in the environment $num_loads, therefore the remaining $num_undef_loads loads are selected randomly!")
        end

        for (index, load) in enumerate(parameters["load"])
            if !haskey(load, "impedance")
                load["impedance"] = rand(["RLC","RL","RC","LC","R","L","C"])
                print("WARNING! The type of load was not specified and is therefore drawn randomly!")

                # TODO 
                # if any(keys(load) in ["R", "L", "C"]
                #   print(WARNING! In the load values were defined which were not used because the type load["impedance"] does not consider the value/s [LIST OF VALUES THAT ARE SPEZIFIED BUT NOT USED]!)

            end

            char = [split(load["impedance"], "")]

            for (_, value) in enumerate(char)
                if value == "R"
                    if !haskey(load, "R")
                        load["R"] = round(rand(Uniform(10, 1e5)), digits=3)
                    end
                end
                if !(value == "L")
                    if !haskey(load, "L")
                        load["L"] = round(rand(Uniform(1e-6, 1e-3)), digits=3)
                    end
                end
                if !(value == "C")
                    if !haskey(load, "C")
                        load["C"] = round(rand(Uniform(1e-9, 1e-4)), digits=3)
                    end
                end
            end
        end

        if num_undef_loads > 0

            loads_distr_undef = get_loads_distr(num_undef_loads)

            num_loads_R_undef = loads_distr_undef[1]
            num_loads_C_undef = loads_distr_undef[2]
            num_loads_L_undef = loads_distr_undef[3]
            num_loads_RL_undef = loads_distr_undef[4]
            num_loads_RC_undef = loads_distr_undef[5]
            num_loads_LC_undef = loads_distr_undef[6]
            num_loads_RLC_undef = loads_distr_undef[7]

            for l in 1:num_loads_RLC_undef
            push!(parameters["load"], _sample_load_RLC())
            end

            for l in 1:num_loads_LC_undef
            push!(parameters["load"], _sample_load_LC())
            end

            for l in 1:num_loads_RL_undef
            push!(parameters["load"], _sample_load_RL())
            end

            for l in 1:num_loads_L_undef
            push!(parameters["load"], _sample_load_L())
            end

            for l in 1:num_loads_RC_undef
            push!(parameters["load"], _sample_load_RC())
            end

            for l in 1:num_loads_C_undef
            push!(parameters["load"], _sample_load_C())
            end

            for l in 1:num_loads_R_undef
            push!(parameters["load"], _sample_load_R())
            end
        end

    end


    if !haskey(parameters, "cable")

        cable_list = []
        for c in 1:num_connections
            push!(cable_list, _sample_cable())
        end
        parameters["cable"] = cable_list
    else
        num_def_cables = length(parameters["cable"])

        num_undef_cables = num_connections - num_def_cables

        @assert num_undef_cables >= 0 "Expect the number of defined cables within the parameter dict to be less or equal to the number of sources in the env, but the entries within the parameter dict is $num_def_cables and the number of env cables is $num_cables."

        if num_undef_cables > 0
            println("WARNING: The number of defined cables $num_def_cables is smaller than the number specified cablesin the environment $num_cables, therefore the remaining $num_undef_cables cables are selected randomly!")
        end

        for (idx, cable) in enumerate(parameters["cable"])
            
            if !haskey(cable, "len")
                cable["len"] = rand(Uniform(1e-3, 1e1))
            end
            
            if !haskey(cable, "Rb")
                cable["Rb"] = 0.722 # TODO: Fixed?!
            end

            if !haskey(cable, "Lb")
                cable["Cb"] = 0.4e-6 # TODO: Fixed?!
            end

            if !haskey(cable, "Cb")
                cable["Lb"] = 0.264e-3 # TODO: Fixed?!
            end

            if !haskey(cable, "R")
                cable["R"] = cable["len"] * cable["Rb"]
            end

            if !haskey(cable, "L")
                cable["L"] = cable["len"] * cable["Lb"]
            end

            if !haskey(cable, "C")
                cable["C"] = cable["len"] * cable["Cb"]
            end

        end

        if num_undef_cables > 0
            for c in 1:num_undef_cables
                push!(parameters["cable"], _sample_cable())
            end
        end

    end
    
    # source_list = []
    # cable_list = []
    # load_list = []

    # parameters = Dict()
    # parameters["source"] = source_list
    # parameters["cable"] = cable_list
    # parameters["load"] = load_list
    # # parameters["grid"] = grid_properties

    return parameters
end

"""
    generate_parameters(
        num_fltr_LCL,
        num_fltr_LC,
        num_fltr_L,
        num_connections,
        num_loads_RLC,
        num_loads_LC,
        num_loads_RL,
        num_loads_RC,
        num_loads_L,
        num_loads_C,
        num_loads_R
        )
    
Function that samples the parameters for the individual elements if no parameters are entered. 
"""
function generate_parameters(num_fltr_LCL, num_fltr_LC, num_fltr_L, num_connections, num_loads_RLC, num_loads_LC, num_loads_RL, num_loads_RC,
                            num_loads_L, num_loads_C, num_loads_R)

    source_list = []
    cable_list = []
    load_list = []
    
    grid_properties = Dict()
    grid_properties["fs"] =  10e3
    grid_properties["v_rms"] = 230
    grid_properties["phase"] = 3


    for s in 1:num_fltr_LCL
        push!(source_list, _sample_fltr_LCL(grid_properties))
    end
    
    for s in 1:num_fltr_LC
        push!(source_list, _sample_fltr_LC(grid_properties))
    end
    
    for s in 1:num_fltr_L
        push!(source_list, _sample_fltr_L(grid_properties))
    end

    for c in 1:num_connections
        push!(cable_list, _sample_cable())
    end

    for l in 1:num_loads_RLC
        push!(load_list, _sample_load_RLC())
    end

    for l in 1:num_loads_LC
        push!(load_list, _sample_load_LC())
    end

    for l in 1:num_loads_RL
        push!(load_list, _sample_load_RL())
    end

    for l in 1:num_loads_L
        push!(load_list, _sample_load_L())
    end

    for l in 1:num_loads_RC
        push!(load_list, _sample_load_RC())
    end

    for l in 1:num_loads_C
        push!(load_list, _sample_load_C())
    end

    for l in 1:num_loads_R
        push!(load_list, _sample_load_R())
    end

    parameters = Dict()
    parameters["source"] = source_list
    parameters["cable"] = cable_list
    parameters["load"] = load_list
    parameters["grid"] = grid_properties

    parameters
end


# """
#     valid_realistic_para(para)

# Checks if the passed parameters e.g. for the filters have logical and realistic values.
# """
# function valid_realistic_para(para)

#     para["source"] = source_list
#     para["net"] = net_para

#     for (i, source) in enumerate(source_list)

#         ZL= 3*(net_para["v_rms"])^2 *(source["pwr"])^-1
#         i_peak = sqrt(2)*net_para["v_rms"]*(ZL)^-1
#         ilfmax = source["i_rip"] * i_peak
#         vcfmax = source["v_rip"] * sqrt(2)*net_para["v_rms"]

#         if source["fltr"] == "LCL"
#             In_L = Interval{Closed, Closed}(0.001*(0.5*source["vdc"]*(4*net_para["fs"]*ilfmax)^-1),(0.5*source["vdc"]*(4*net_para["fs"]*ilfmax)^-1)*1000)
#             In_C = Interval{Closed, Closed}(0.001*(ilfmax*(8*net_para["fs"]*vcfmax)^-1),(ilfmax*(8*net_para["fs"]*vcfmax)^-1)*1000)
#             In_R_L = Interval{Closed, Closed}(0.001*(400 * source["L1"]),(400 * source["L1"])*1000)
#             In_R_C = Interval{Closed, Closed}(0.001*(28* source["C"]),(28* source["C"])*1000)

#             if (!(source["L1"] in In_L) || !(source["L2"] in In_L) || !(source["C"] in In_C) || !(source["R1"] in In_R_L) || !(source["R2"] in In_R_L) || !(source["R_C"] in In_R_C))
#                 @warn " Source $i contains filter parameters that are not recommended."
#             end

#         elseif source["fltr"] == "LC"
#             In_L = Interval{Closed, Closed}(0.001*(source["vdc"]*(4*net_para["fs"]*ilfmax)^-1),(source["vdc"]*(4*net_para["fs"]*ilfmax)^-1)*1000)
#             In_C = Interval{Closed, Closed}(0.001*(ilfmax*(8*net_para["fs"]*vcfmax)^-1),(ilfmax*(8*net_para["fs"]*vcfmax)^-1)*1000)
#             In_R_L = Interval{Closed, Closed}(0.001*(400 * source["L1"]),(400 * source["L1"])*1000)
#             In_R_C = Interval{Closed, Closed}(0.001*(28* source["C"]),(28* source["C"])*1000)

#             if (!(source["L1"] in In_L) || !(source["C"] in In_C) || !(source["R1"] in In_R_L) || !(source["R_C"] in In_R_C))
#                 @warn " Source $i contains filter parameters that are not recommended."
#             end

#         elseif source["fltr"] == "L"
#             In_L = Interval{Closed, Closed}(0.001*(source["vdc"]*(4*net_para["fs"]*ilfmax)^-1),(source["vdc"]*(4*net_para["fs"]*ilfmax)^-1)*1000)
#             In_R_L = Interval{Closed, Closed}(0.001*(400 * source["L1"]),(400 * source["L1"])*1000)

#             if (!(source["L1"] in In_L) || !(source["R1"] in In_R_L))
#                 @warn " Source $i contains filter parameters that are not recommended."
#             end

#         end
#     end

# end

"""
    cntr_fltrs(source_list)

Counts the number of filter types, if `parameters` is passed.
"""
function cntr_fltrs(source_list)
    cntr_LCL = 0
    cntr_LC = 0
    cntr_L = 0


    
    for (i, source) in enumerate(source_list)
        # println(source)
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
    cntr_loads(load_list)

Counts the number of load types, if `parameters` is passed.
"""
function cntr_loads(load_list)
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
    _sample_fltr_LCL()

Sample parameters for the LCL filter.
"""
function _sample_fltr_LCL(grid_properties)

    source = Dict()
    source["source_type"] = "ideal"
    source["fltr"] = "LCL"
    
    source["pwr"] = rand(range(start=5,step=5,stop=50))*1e3
    source["vdc"] = rand(range(start=690,step=10,stop=800))
    source["i_rip"] = rand(Uniform(0.1, 0.15))
    source["v_rip"] = rand(Uniform(0.014, 0.016))

   #Inductor design

   Vorms = grid_properties["v_rms"]*1.05
   Vop = Vorms*sqrt(2)

   Zl = 3*Vorms*Vorms/source["pwr"]

   Iorms = Vorms/Zl
   Iop = Iorms*sqrt(2)

   ΔIlfmax = source["i_rip"]*Iop

   source["L1"] = 0.5*(source["vdc"]*(4*grid_properties["fs"]*ΔIlfmax)^-1)
   source["L2"] = deepcopy(source["L1"]) 
   source["i_limit"]= Iop+ ΔIlfmax

   #Capacitor design

   Vorms = grid_properties["v_rms"]*0.95
   Vop = Vorms*sqrt(2)

   Zl = 3*Vorms*Vorms/source["pwr"]

   Iorms = Vorms/Zl
   Iop = Iorms*sqrt(2)
   Ir_d = source["vdc"]/(4*grid_properties["fs"]*source["L1"]*Iop)
   ΔIlfmax = Ir_d*Iop
   ΔVcfmax = source["v_rip"]*Vop

   source["C"] = ΔIlfmax/(8*grid_properties["fs"]*ΔVcfmax)

   source["R1"] = 400 * source["L1"]
   source["R2"] = deepcopy(source["R1"])   
   source["R_C"] = 28* source["C"]
   source["v_limit"]= Vop+ΔVcfmax
   
   

   source
end


"""
    _sample_fltr_LC()

Sample parameters for the LC filter.
"""
function _sample_fltr_LC(grid_properties)  

    source = Dict()
    source["source_type"] = "ideal"
    source["fltr"] = "LC"

    source["pwr"] = rand(range(start=5,step=5,stop=50))*1e3
    source["vdc"] = rand(range(start=690,step=10,stop=800))
    source["i_rip"] = rand(Uniform(0.1, 0.15))
    source["v_rip"] = rand(Uniform(0.014, 0.016))

    #Inductor design

    Vorms = grid_properties["v_rms"]*1.05
    Vop = Vorms*sqrt(2)

    Zl = 3*Vorms*Vorms/source["pwr"]

    Iorms = Vorms/Zl
    Iop = Iorms*sqrt(2)

    ΔIlfmax = source["i_rip"]*Iop

    source["L1"] = 0.5*(source["vdc"]*(4*grid_properties["fs"]*ΔIlfmax)^-1)
    source["i_limit"]= Iop+ ΔIlfmax

    #Capacitor design

    Vorms = grid_properties["v_rms"]*0.95
    Vop = Vorms*sqrt(2)

    Zl = 3*Vorms*Vorms/source["pwr"]

    Iorms = Vorms/Zl
    Iop = Iorms*sqrt(2)
    Ir_d = source["vdc"]/(4*grid_properties["fs"]*source["L1"]*Iop)
    ΔIlfmax = Ir_d*Iop
    ΔVcfmax = source["v_rip"]*Vop

    source["C"] = ΔIlfmax/(8*grid_properties["fs"]*ΔVcfmax)

    source["R1"] = 400 * source["L1"]
    source["R_C"] = 28* source["C"]
    source["v_limit"]= Vop+ΔVcfmax

    source
end


"""
    _sample_fltr_L()

Sample parameters for the L filter.
"""
function _sample_fltr_L(grid_properties)  

    source = Dict()
    source["source_type"] = "ideal"
    source["fltr"] = "L"

    source["pwr"] = rand(range(start=5,step=5,stop=50))*1e3
    source["vdc"] = rand(range(start=690,step=10,stop=800))
    source["i_rip"] = rand(Uniform(0.1, 0.15))
    source["v_rip"] = rand(Uniform(0.014, 0.016))

    #Inductor design

    Vorms = grid_properties["v_rms"]*1.05
    Vop = Vorms*sqrt(2)

    Zl = 3*Vorms*Vorms/source["pwr"]

    Iorms = Vorms/Zl
    Iop = Iorms*sqrt(2)

    ΔIlfmax = source["i_rip"] *Iop

    source["L1"] = 0.5*(source["vdc"]*(4*grid_properties["fs"]*ΔIlfmax)^-1)

    source["R1"] = 400 * source["L1"] 
    source["i_limit"]= Iop+ ΔIlfmax

    source
end


"""
    _sample_load_RLC()

Sample parameters for the RLC load.
"""
function _sample_load_RLC()

    load = Dict()
    load["impedance"] = "RLC"
    load["R"] = round(rand(Uniform(10, 1e5)), digits=3)
    load["L"] = round(rand(Uniform(1e-3, 1e2)), digits=3)
    load["C"] = round(rand(Uniform(1e-3, 1e2)), digits=3)

    load
end


"""
    _sample_load_LC()

Sample parameters for the LC load.
"""
function _sample_load_LC()

    load = Dict()
    load["impedance"] = "LC"
    load["L"] = round(rand(Uniform(1e-3, 1e2)), digits=3)
    load["C"] = round(rand(Uniform(1e-3, 1e2)), digits=3)


    load
end


"""
    _sample_load_RL()

Sample parameters for the RL load.
"""
function _sample_load_RL()

    load = Dict()
    load["impedance"] = "RL"
    load["R"] = round(rand(Uniform(10, 1e5)), digits=3)
    load["L"] = round(rand(Uniform(1e-3, 1e2)), digits=3)

    load
end


"""
    _sample_load_RC()

Sample parameters for the RC load.
"""
function _sample_load_RC()

    load = Dict()
    load["impedance"] = "RC"
    load["R"] = round(rand(Uniform(10, 1e5)), digits=3)
    load["C"] = round(rand(Uniform(1e-3, 1e2)), digits=3)

    load
end


"""
    _sample_load_L()

Sample parameters for the L load.
"""
function _sample_load_L()

    load = Dict()

    load["impedance"] = "L"
    load["L"] = round(rand(Uniform(1e-3, 1e2)), digits=3)

    load
end


"""
    _sample_load_C()

Sample parameters for the C load.
"""
function _sample_load_C()

    load = Dict()

    load["impedance"] = "C"
    load["C"] = round(rand(Uniform(1e-3, 1e2)), digits=3)

    load
end


"""
    _sample_load_R()

Sample parameters for the R load.
"""
function _sample_load_R()

    load = Dict()

    load["impedance"] = "R"
    load["R"] = round(rand(Uniform(10, 1e5)), digits=3)

    load
end


"""
    _sample_cable()

Sample parameters for the cable.
"""
function _sample_cable()
    
    l =  rand(Uniform(1e-3, 1e1))

    Rb = 0.722
    Cb = 0.4e-6
    Lb = 0.264e-3

    cable = Dict()
    cable["R"] = l * Rb
    cable["L"] = l * Lb
    cable["C"] = l * Cb

    cable
end


"""
    tobe_or_n2b(cntr, x, p)

Sets x to zero or to the value of the counter as a function of p and increases it as well.
"""
function tobe_or_n2b(cntr, x, p)
    
    if x < p
        cntr += 1  
        return cntr, cntr
    else
        x = 0
        return cntr, x
    end
end


"""
    CM_generate(tot_ele, num_sources, S2L_p, S2S_p)

Returns the constructed CM and the total number of connections.
"""
function CM_generate(tot_ele, num_sources, S2L_p, S2S_p)

    # counting the connections 
    cntr = 0
    
    # get a upper triangular matrix
    mask = UpperTriangular(ones(tot_ele, tot_ele))
    CM = rand(tot_ele,tot_ele) .* mask # fill matrix with random entries between [0,1]
    CM = CM - Diagonal(CM) # delet diagonal bc no connection with itself

    # go through the matrix
    # -1 bc last entry is 0 anyway
    for i in 1:tot_ele-1
        # start at i, bc we need to check only upper triangle
        for j in i:tot_ele-1
            if j >= num_sources-1 # select propability according to column
                cntr, x = tobe_or_n2b(cntr, CM[i, j+1], S2L_p)
                CM[i, j+1] = x
            else
                cntr, x = tobe_or_n2b(cntr, CM[i, j+1], S2S_p)
                CM[i, j+1] = x
            end
        end
    end

    # make sure that no objects disappear or subnets are formed
    for i in 1:tot_ele
        # save rows and columns entries
        Col = CM[1:i-1, i]
        Row = CM[i, i+1:tot_ele]
        
        # get one list in the form of: [column, row]-entries
        entries = vcat(Col, Row)

        non_zero = count(i->(i != 0), entries) # number of non_zero entries
        zero = count(i->(i == 0), entries) # number of zero entries

        val_to_set = min(2, zero) # minimum of connections is 2
        
        if non_zero <= 2 # we need to set values if there are less then 2 entries
            idx_row_entries = findall(x->x==0, Col) # Get rows of the entries = 0
            idx_col_entries = findall(x->x==0, Row) # Get col of the entries = 0

            idx_list = vcat([(j,i) for j in idx_row_entries], [(i,i+j) for j in idx_col_entries])
            
            samples = min(val_to_set, length(idx_list))
            idx_rnd = sample(1:length(idx_list), samples) # draw samples from the list
            
            for (a, ix) in enumerate(idx_rnd)
                # Based on the random sample, select an indize
                # from the list and write into the corresponding CM cell.
                cntr += 1
                CM[idx_list[ix][1], idx_list[ix][2]] = cntr
            end
        end
    end

    CM = CM - CM' # copy with negative sign to lower triangle

    return cntr, CM
end


"""
    get_A_src(self::NodeConstructor, source_i)

Create the A_src entry for a source in the A matrix.
"""
function get_A_src(self::NodeConstructor, source_i)

    parameter_i = self.parameters["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        A_src = zeros(4,4)
        A_src[1,1] = -(parameter_i["R1"]+parameter_i["R_C"])/parameter_i["L1"]
        A_src[1,2] = -1/parameter_i["L1"]
        A_src[1,3] = parameter_i["R2"]/parameter_i["L1"]
        A_src[2,1] = 1/parameter_i["C"]
        A_src[2,3] = -1/parameter_i["C"]
        A_src[3,1] = parameter_i["R_C"]/parameter_i["L2"]
        A_src[3,2] = 1/parameter_i["L2"]
        A_src[3,3] = -(parameter_i["R_C"]+parameter_i["R2"])/parameter_i["L2"]
        A_src[3,4] = -1/parameter_i["L2"]
        
        C_sum =  0
        
        CM_row = self.CM[source_i,:]
        
        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        
        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end
        
        A_src[4,3] = C_sum^(-1)

    elseif parameter_i["fltr"] == "LC"

        A_src = zeros(3,3)
        A_src[1,1] = -parameter_i["R1"]/parameter_i["L1"]
        A_src[1,3] = -1/parameter_i["L1"]
        A_src[2,2] = -1/(parameter_i["C"]*parameter_i["R_C"])
        A_src[2,3] = 1/(parameter_i["C"]*parameter_i["R_C"])
        
        C_sum =  0
        
        CM_row = self.CM[source_i, :]
        
        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        
        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end
        
        A_src[3,1] = C_sum^(-1)
        A_src[3,2] = (C_sum * parameter_i["R_C"])^(-1)
        A_src[3,3] = -(C_sum * parameter_i["R_C"])^(-1)

    elseif parameter_i["fltr"] == "L"

        A_src = zeros(2,2)
        A_src[1,1] = -parameter_i["R1"]/parameter_i["L1"]
        A_src[1,2] = -1/parameter_i["L1"]
        
        C_sum =  0
        
        CM_row = self.CM[source_i, :]
        
        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        
        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end
        
        A_src[2,1] = C_sum^(-1)
    else
        throw("Expect filter to be \"LCL\", \"LC\" or \"L\", not $(parameter_i["fltr"]).")
    end

    return A_src

end


"""
    get_B_source(self::NodeConstructor, source_i)

Create the B_source entry for a source in the B matrix.
"""
function get_B_source(self::NodeConstructor, source_i)

    parameter_i = self.parameters["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        B_source = zeros(4,1)
        B_source[1,1] =  1/parameter_i["L1"]

    elseif parameter_i["fltr"] == "LC"
        B_source = zeros(3,1)
        B_source[1,1] =  1/parameter_i["L1"]

    elseif parameter_i["fltr"] == "L" 
        B_source = zeros(2,1)
        B_source[1,1] =  1/parameter_i["L1"]
    end

    return B_source
end


"""
    get_A_src_trn_c(self::NodeConstructor, source_i)

Create the A_src_trn_c entry in the A matrix.
"""
function get_A_src_trn_c(self::NodeConstructor, source_i)

    parameter_i = self.parameters["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        A_src_trn_c = zeros(4, self.num_connections)

        CM_row = self.CM[source_i,:]

        indizes = CM_row[CM_row .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        C_sum = 0 

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_c[4,idx] = sign * -(C_sum^(-1))
        end

    elseif parameter_i["fltr"] == "LC"

        A_src_trn_c = zeros(3, self.num_connections)

        CM_row = self.CM[source_i,:]

        indizes = CM_row[CM_row .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        C_sum = 0

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_c[3,idx] = sign * -(C_sum^(-1))
        end

    elseif parameter_i["fltr"] == "L"

        A_src_trn_c = zeros(2, self.num_connections)

        CM_row = self.CM[source_i,:]

        indizes = CM_row[CM_row .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        C_sum = 0

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_c[2,idx] = sign * -(C_sum^(-1))
        end
    end

    return A_src_trn_c
end


"""
    get_A_src_trn_l(self::NodeConstructor, source_i)

Create the A_src_trn_l entry in the A matrix.
"""
function get_A_src_trn_l(self::NodeConstructor, source_i)

    parameter_i = self.parameters["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        A_src_trn_l = zeros(4, self.num_connections)
        CM_col = self.CM[source_i,:]
        
        indizes = CM_col[CM_col .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        
        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_l[4,idx] = sign * 1/self.parameters["cable"][idx]["L"]
        end

    elseif parameter_i["fltr"] == "LC"
        A_src_trn_l = zeros(3, self.num_connections)
        
        CM_col = self.CM[source_i,:]
        
        indizes = CM_col[CM_col .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        
        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_l[3,idx] = sign * 1/self.parameters["cable"][idx]["L"]
        end

    elseif parameter_i["fltr"] == "L"
        A_src_trn_l = zeros(2, self.num_connections)
        
        CM_col = self.CM[source_i,:]
        
        indizes = CM_col[CM_col .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        
        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_src_trn_l[2,idx] = sign * 1/self.parameters["cable"][idx]["L"]
        end
    end

    return A_src_trn_l'
end


"""
    generate_A_trn(self::NodeConstructor)

Create the A_tran entry in the A matrix.
"""
function generate_A_trn(self::NodeConstructor)
    """Create A_trn"""

    vec = zeros(self.num_connections)
    for (i, ele) in enumerate(self.parameters["cable"])
        vec[i] = -ele["R"]/ele["L"]
    end
    A_trn = Diagonal(vec)

    return A_trn
end


"""
    get_A_tran_load_c(self::NodeConstructor, load_i)

Create the A_tran_load_c entry in the A matrix.
"""
function get_A_tran_load_c(self::NodeConstructor, load_i)

    parameter_i = self.parameters["load"][load_i]
    
    if parameter_i["impedance"] == "RLC" || parameter_i["impedance"] == "LC"
        A_tran_load_c = zeros(2, self.num_connections)

        CM_row = self.CM[self.num_sources + load_i, :]

        indizes = CM_row[CM_row .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs

        indizes_ = indizes .* signs # delet signs from indices

        C_sum = parameter_i["C"]

        for ele in self.parameters["cable"]
            C_sum += ele["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_tran_load_c[1,idx] = sign * -(C_sum^(-1))
        end

    elseif parameter_i["impedance"] == "RL" || parameter_i["impedance"] == "L"
        A_tran_load_c = zeros(2, self.num_connections)

        CM_row = self.CM[self.num_sources + load_i, :]

        indizes = CM_row[CM_row .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs

        indizes_ = indizes .* signs # delet signs from indices

        C_sum = 0

        for ele in self.parameters["cable"]
            C_sum += ele["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_tran_load_c[1, idx] = sign * -(C_sum^(-1))
        end

    elseif parameter_i["impedance"] == "RC" || parameter_i["impedance"] == "C"
        A_tran_load_c = zeros(1, self.num_connections)

        CM_row = self.CM[self.num_sources + load_i, :]

        indizes = CM_row[CM_row .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs

        indizes_ = indizes .* signs # delet signs from indices

        C_sum = parameter_i["C"]

        for ele in self.parameters["cable"]
            C_sum += ele["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_tran_load_c[idx] = sign * -(C_sum^(-1))
        end

    elseif parameter_i["impedance"] == "R"
        A_tran_load_c = zeros(1, self.num_connections)

        CM_row = self.CM[self.num_sources + load_i, :]

        indizes = CM_row[CM_row .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs

        indizes_ = indizes .* signs # delet signs from indices

        C_sum = 0

        for ele in self.parameters["cable"]
            C_sum += ele["C"] * 0.5
        end

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_tran_load_c[idx] = sign * -(C_sum^(-1))
        end
    end

    return A_tran_load_c
end


"""
    get_A_tran_load_l(self::NodeConstructor, load_i)

Create the A_tran_load_l entry in the A matrix.
"""
function get_A_tran_load_l(self::NodeConstructor, load_i)
    
    parameter_i = self.parameters["load"][load_i]
    
    if parameter_i["impedance"] == "RLC" || parameter_i["impedance"] == "LC" || parameter_i["impedance"] == "RL" || parameter_i["impedance"] == "L"
        A_tran_load_l = zeros(self.num_connections, 2)

        CM_col = self.CM[self.num_sources + load_i, :]

        indizes = CM_col[CM_col .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)

            A_tran_load_l[idx, 1] = sign * 1/self.parameters["cable"][idx]["L"] 
        end 
        
    elseif parameter_i["impedance"] == "RC" || parameter_i["impedance"] == "C" || parameter_i["impedance"] == "R"
        A_tran_load_l = zeros(self.num_connections, 1)

        CM_col = self.CM[self.num_sources + load_i, :]

        indizes = CM_col[CM_col .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)

            A_tran_load_l[idx, 1] = sign * 1/self.parameters["cable"][idx]["L"] 
        end 
    end

    return A_tran_load_l
end


"""
    get_A_load(self::NodeConstructor, load_i)

Create the A_tran_load_l entry in the A matrix.
"""
function get_A_load(self::NodeConstructor, load_i)
    parameter_i = self.parameters["load"][load_i]

    if parameter_i["impedance"] == "RLC"
        A_load = zeros(2, 2)
        A_load[2,1] = 1 / parameter_i["L"]

        C_sum = parameter_i["C"]

        CM_row = self.CM[self.num_sources + load_i, :]

        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes]

        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1,1] = - ((parameter_i["R"]) * C_sum)^(-1)
        A_load[1,2] = - (C_sum)^(-1)

    elseif parameter_i["impedance"] == "LC"
        A_load = zeros(2, 2)
        A_load[2,1] = 1 / parameter_i["L"]

        C_sum = parameter_i["C"]

        CM_row = self.CM[self.num_sources + load_i, :]

        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes]

        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1,2] = - (C_sum)^(-1)

    elseif parameter_i["impedance"] == "RL"
        A_load = zeros(2, 2)
        A_load[2,1] = 1 / parameter_i["L"]

        C_sum = 0

        CM_row = self.CM[self.num_sources + load_i, :]

        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes]

        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1,1] = - ((parameter_i["R"]) * C_sum)^(-1)
        A_load[1,2] = - (C_sum)^(-1)

    elseif parameter_i["impedance"] == "RC"
        A_load = zeros(1, 1)

        C_sum = parameter_i["C"]

        CM_row = self.CM[self.num_sources + load_i, :]

        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes]

        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1,1] = - ((parameter_i["R"]) * C_sum)^(-1)

    elseif parameter_i["impedance"] == "L"
        A_load = zeros(2, 2)
        A_load[2,1] = 1 / parameter_i["L"]

        C_sum = 0

        CM_row = self.CM[self.num_sources + load_i, :]

        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes]

        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1,2] = - (C_sum)^(-1)

    elseif parameter_i["impedance"] == "C"
        A_load = zeros(1, 1)

    elseif parameter_i["impedance"] == "R"
        A_load = zeros(1, 1)

        C_sum = 0

        CM_row = self.CM[self.num_sources + load_i, :]

        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes]

        indizes_ = indizes .* signs

        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameters["cable"][idx]["C"] * 0.5
        end

        A_load[1,1] = - ((parameter_i["R"]) * C_sum)^(-1)

    else
        throw("Expect Impedance to be \"RLC\", \"LC\", \"RL\", \"RC\", \"L\", \"C\" or \"R\", not $(parameter_i["impedance"]).")
    end

    return A_load
end


"""
    generate_A(self::NodeConstructor)

Generates the A matrix by joining the individual sub-matrices together.
"""
function generate_A(self::NodeConstructor)
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
    A_src_list = [get_A_src(self, i) for i in 1:self.num_sources]
    println(A_src_list)
    for (i, ele) in enumerate(A_src_list)
        if i <= self.num_fltr_LCL
        # IDEE: if self.parameters["source"][i]["ftlr"]  == "LCL"
            start = 4 * i - 3
            stop = 4 * i

            A_src[start:stop, start:stop] = ele

        elseif i <= self.num_fltr_LCL + self.num_fltr_LC
            start = 3 * i + self.num_fltr_LCL - 2
            stop = 3 * i + self.num_fltr_LCL

            A_src[start:stop, start:stop] = ele

        elseif i <= self.num_fltr_LCL + self.num_fltr_LC + self.num_fltr_L
            start = 2 * i + 2 * self.num_fltr_LCL + self.num_fltr_LC - 1
            stop = 2 * i + 2 * self.num_fltr_LCL + self.num_fltr_LC

            A_src[start:stop,start:stop] = ele
        end
    end

    # get A_src_trn_c
    A_src_trn_c = zeros(self.num_fltr, self.num_connections)
    A_src_trn_c_list = [get_A_src_trn_c(self, i) for i in 1:self.num_sources] # start at 1 bc Source 1 ...

    for (i, ele) in enumerate(A_src_trn_c_list)
        if i <= self.num_fltr_LCL
            start = 4 * i - 3
            stop = 4 * i

            A_src_trn_c[start:stop,:] = ele

        elseif i <= self.num_fltr_LCL+self.num_fltr_LC
            start = 3 * i + self.num_fltr_LCL - 2
            stop = 3 * i + self.num_fltr_LCL
            A_src_trn_c[start:stop,:] = ele

        elseif i <= self.num_fltr_LCL+self.num_fltr_LC+self.num_fltr_L
            start = 2 * i + 2 * self.num_fltr_LCL + self.num_fltr_LC - 1
            stop = 2 * i + 2 * self.num_fltr_LCL + self.num_fltr_LC

            A_src_trn_c[start:stop,:] = ele
        end
    end

    # get A_src_trn_l
    A_src_trn_l = zeros(self.num_connections, self.num_fltr)
    A_src_trn_l_list = [get_A_src_trn_l(self, i) for i in 1:self.num_sources] # start at 1 bc Source 1 ...

    for (i, ele) in enumerate(A_src_trn_l_list)
        if i <= self.num_fltr_LCL
            start = 4 * i - 3
            stop = 4 * i

            A_src_trn_l[:,start:stop] = ele

        elseif i <= self.num_fltr_LCL+self.num_fltr_LC
            start = 3 * i + self.num_fltr_LCL - 2
            stop = 3 * i + self.num_fltr_LCL

            A_src_trn_l[:,start:stop] = ele

        elseif i <= self.num_fltr_LCL+self.num_fltr_LC+self.num_fltr_L
            start = 2 * i + 2 * self.num_fltr_LCL + self.num_fltr_LC - 1
            stop = 2 * i + 2 * self.num_fltr_LCL + self.num_fltr_LC

            A_src_trn_l[:,start:stop] = ele
        end
    end

    A_trn = generate_A_trn(self)

    A_tran_load_l_list = []
    for i in 1:self.num_loads
        push!(A_tran_load_l_list, get_A_tran_load_l(self, i))
    end
    A_tran_load_l = reduce(hcat, A_tran_load_l_list) # i-> idx // i+1 -> num of load

    A_tran_load_c_list = []
    for i in 1:self.num_loads
        push!(A_tran_load_c_list, get_A_tran_load_c(self, i))
    end
    A_tran_load_c = reduce(vcat, A_tran_load_c_list)

    # get A_load_diag
    self.num_impedance = (2 * (self.num_loads_RLC
                                + self.num_loads_LC
                                + self.num_loads_RL
                                + self.num_loads_L)
                                + self.num_loads_RC + self.num_loads_C + self.num_loads_R)
    
    A_load_diag = zeros(self.num_impedance, self.num_impedance)
    A_load_list = [get_A_load(self, i) for i in 1:self.num_loads]

    for (i, ele) in enumerate(A_load_list)
        if i <= self.num_loads_RLC + self.num_loads_LC + self.num_loads_RL + self.num_loads_L

            start = 2 * i - 1
            stop = 2 * i
            A_load_diag[start:stop, start:stop] = ele

        elseif i <= self.num_loads_RLC + self.num_loads_LC + self.num_loads_RL + self.num_loads_L + self.num_loads_RC + self.num_loads_C + self.num_loads_R

            start = i + self.num_loads_RLC + self.num_loads_LC + self.num_loads_RL + self.num_loads_L
            A_load_diag[start:start, start:start] = ele

        end
    end

    #         A_transition = np.block([A_trn, A_tran_load_c],
    #                                 [A_tran_load_c, A_load_diag])

    A_load_zeros = zeros(self.num_fltr, self.num_impedance)

    A_load_zeros_t = A_load_zeros'

    A = [A_src       A_src_trn_c        A_load_zeros
        A_src_trn_l           A_trn  A_tran_load_l
        A_load_zeros_t  A_tran_load_c   A_load_diag]

    if self.parameters["grid"]["phase"] === 1
            return A
        elseif self.parameters["grid"]["phase"] === 3
            z = zeros(size(A))
            A_ = [A z z; z A z; z z A]
            return A_
        end
end


"""
    generate_B(self::NodeConstructor)

Generates the B matrix by joining the individual sub-matrices together.
"""
function generate_B(self::NodeConstructor)
    """
    The previously constructed matrices are now plugged together in the form:
        [[B_source_1,          0, ...,           0],
         [         0, B_source_2, ...,           0],
         [         0,          0, ...,           0],
         [         0,          0, ...,  B_source_n]]

    """
    B = zeros(self.num_fltr + self.num_connections + self.num_impedance, self.num_sources)

    B_source_list = [get_B_source(self, i) for i in 1:self.num_sources] # start at 1 bc Source 1 ...

    for (i, ele) in enumerate(B_source_list)
        if i <= self.num_fltr_LCL
            start = 4 * i - 3
            stop = 4 * i

            B[start:stop,i] = ele

        elseif i <= self.num_fltr_LCL + self.num_fltr_LC
            start = 3 * i + self.num_fltr_LCL - 2
            stop = 3 * i + self.num_fltr_LCL

            B[start:stop,i] = ele

        elseif i <= self.num_fltr_LCL + self.num_fltr_LC + self.num_fltr_L
            start = 2 * i + 2 * self.num_fltr_LCL + self.num_fltr_LC - 1
            stop = 2 * i + 2 * self.num_fltr_LCL + self.num_fltr_LC

            B[start:stop,i:i] = ele
        end
    end
    if self.parameters["grid"]["phase"] === 1
        return B
    elseif self.parameters["grid"]["phase"] === 3
        z = zeros(size(B))
        B_ = [B z z ; z B z; z z B]
        return B_
    end
end


"""
    generate_C(self::NodeConstructor)

Generates the C matrix.
"""
function generate_C(self::NodeConstructor)
    """Generate the C matrix

    Retruns:
        C: Identity matrix (2*num_sources+num_connections)
    """
    C =  Diagonal(ones(self.num_fltr + self.num_connections + self.num_impedance))
    if self.parameters["grid"]["phase"] === 1
        return C
    elseif self.parameters["grid"]["phase"] === 3
        z = zeros(size(C))
        C_ = [C z z;
              z C z;
              z z C]
        return C_
    end
end

"""
    generate_D(self::NodeConstructor)

Generates the D matrix.
"""
function generate_D(self::NodeConstructor)
    return 0
end


"""
    get_sys(self::NodeConstructor)

Generates the system matrices A, B, C and D.
"""
function get_sys(self::NodeConstructor)
    """Returns state space matrices"""

    A = generate_A(self)
    B = generate_B(self)
    C = generate_C(self)
    D = generate_D(self)
    return (A, B, C, D)
end

"""
    get_state_ids(self::NodeConstructor)

Creates the State Vector for an related NodeConstructor and outputs it as a list of strings.
"""

#TODO: remove legacy workaround
get_states(self::NodeConstructor) = get_state_ids(self)

function get_state_ids(self::NodeConstructor)
    states = []
    for s in 1:self.num_sources
        if s <= self.num_fltr_LCL
            push!(states, "source$s"*"_i_L1")    # i_f1; dann i_f2....
            push!(states, "source$s"*"_v_C")
            push!(states, "source$s"*"_i_L2")
            push!(states, "source$s"*"_v_C_cables")
        
        elseif s <= self.num_fltr_LCL + self.num_fltr_LC
            push!(states, "source$s"*"_i_L1")
            push!(states, "source$s"*"_v_C")
            push!(states, "source$s"*"_v_C_cables")
        
        elseif s <= self.num_fltr_LCL + self.num_fltr_LC + self.num_fltr_L
            push!(states, "source$s"*"_i_L1")
            push!(states, "source$s"*"_v_C_cables")
        end
    end

    for c in 1:self.num_connections
        push!(states, "cable$c"*"_i_L")
    end

    for l in 1:self.num_loads
        if l <= self.num_loads_RLC + self.num_loads_LC + self.num_loads_RL + self.num_loads_L
            push!(states, "load$l"*"_v_C_total")
            push!(states, "load$l"*"_i_L")
        elseif l <= self.num_loads_RLC + self.num_loads_LC + self.num_loads_RL + self.num_loads_L + self.num_loads_RC + self.num_loads_C + self.num_loads_R
            push!(states, "load$l"*"_v_C_total")
        end
    end

    if self.parameters["grid"]["phase"] === 3
        A = ["_a", "_b", "_c"]
        states = vcat([broadcast(*, states, A[i]) for i in 1:3]...)
    end
    
    return states
end

"""
    get_state_paras(self::NodeConstructor)

Creates a Vector containing the related L or C Parameters for the states of a NodeConstructor.
"""

function get_state_paras(self::NodeConstructor)
    state_paras = []
    for s in 1:self.num_sources
        if s <= self.num_fltr_LCL
            push!(state_paras, self.parameters["source"][s]["L1"]) 
            push!(state_paras, self.parameters["source"][s]["C"]) 
            push!(state_paras, self.parameters["source"][s]["L2"])
            push!(state_paras, get_C_sum_cable_node(s,self))
        
        elseif s <= self.num_fltr_LCL + self.num_fltr_LC
            push!(state_paras, self.parameters["source"][s]["L1"])
            push!(state_paras, self.parameters["source"][s]["C"])
            push!(state_paras, get_C_sum_cable_node(s,self))  
        
        elseif s <= self.num_fltr_LCL + self.num_fltr_LC + self.num_fltr_L
            push!(state_paras, self.parameters["source"][s]["L1"])
            push!(state_paras, get_C_sum_cable_node(s,self))
        end
    end

    for c in 1:self.num_connections
        push!(state_paras, self.parameters["cable"][c]["L"]) 
    end

    for l in 1:self.num_loads
        if l <= self.num_loads_RLC + self.num_loads_LC + self.num_loads_RL + self.num_loads_L
            c=0
                if haskey(self.parameters["load"][l], "C")
                    c= self.parameters["load"][l]["C"]
                end
            push!(state_paras, get_C_sum_cable_node(self.num_sources+l,self) + c) 
            push!(state_paras, self.parameters["load"][l]["L"])  
        elseif l <= self.num_loads_RLC + self.num_loads_LC + self.num_loads_RL + self.num_loads_L + self.num_loads_RC + self.num_loads_C + self.num_loads_R
            c=0
            if haskey(self.parameters["load"][l], "C")
                c= self.parameters["load"][l]["C"]
            end
            push!(state_paras, get_C_sum_cable_node(self.num_sources+l,self) + c)  
        end
    end

    if self.parameters["grid"]["phase"] === 3
        state_paras = vcat([state_paras for i in 1:3]...)
    end

    return state_paras
end

function get_C_sum_cable_node(node_i,self::NodeConstructor)
    CM_row = self.CM[node_i,:]
    C_sum=0
    indizes = CM_row[CM_row .!= 0]
    signs = [sign(x) for x in indizes] # get signs
    indizes_ = indizes .* signs # delet signs from indices
        
    for idx in indizes_
        idx = Int(idx)
        C_sum += self.parameters["cable"][idx]["C"] * 0.5
    end
    return C_sum
end

"""
    get_action_ids(self::NodeConstructor)

Creates the State Vector for an related NodeConstructor and outputs it as a list of strings.
"""

function get_action_ids(self::NodeConstructor)
    actions = []
    
    for s in 1:self.num_sources
        if s <= self.num_fltr_LCL
            push!(actions, "source$s"*"_u")
        
        elseif s <= self.num_fltr_LCL + self.num_fltr_LC
            push!(actions, "source$s"*"_u")
        
        elseif s <= self.num_fltr_LCL + self.num_fltr_LC + self.num_fltr_L
            push!(actions, "source$s"*"_u")
        end
    end

    if self.parameters["grid"]["phase"] === 3
        A = ["_a", "_b", "_c"]
        actions = vcat([broadcast(*, actions, A[i]) for i in 1:3]...)
    end
    

    return actions
end

"""
    get_source_state_indices(self::NodeConstructor,sources)

Returns all state indices for passed sources.
"""

function get_source_state_indices(self::NodeConstructor,sources)
    state_ids=get_state_ids(self)
    action_ids=get_action_ids(self)
    source_indices = Dict()
    for idx in sources
        source = Dict()
        state_indices = findall(x -> occursin("source$idx"*"_", x), state_ids)
        action_indices = findall(x -> occursin("source$idx"*"_", x), action_ids)
        source["state_indices"]= state_indices
        source["action_indices"]= action_indices
        source_indices["source$idx"]= source
    end

    return source_indices
end

"""
    get_cable_state_indices(self::NodeConstructor,cables)

Returns all state indices for passed cables.
"""

function get_cable_state_indices(self::NodeConstructor,cables)
    state_ids=get_state_ids(self)
    cable_indices = Dict()
    for idx in cables
        cable = Dict()
        state_indices = findall(x -> occursin("cable$idx"*"_", x), state_ids)
        cable["state_indices"]= state_indices
        cable_indices["cable$idx"]= cable
    end

    return cable_indices
end

"""
    get_load_state_indices(self::NodeConstructor,loads)

Returns all state indices for passed loads.
"""

function get_load_state_indices(self::NodeConstructor,loads)
    state_ids=get_state_ids(self)
    load_indices = Dict()
    for idx in loads
        load = Dict()
        state_indices = findall(x -> occursin("load$idx"*"_", x), state_ids)
        load["state_indices"]= state_indices
        load_indices["load$idx"]= load
    end

    return load_indices
end

function draw_graph(self::NodeConstructor)
    """Plots a graph according to the CM matrix

    Red nodes corresponse to a source.
    Lightblue nodes corresponse to a load.
    """

    # edges = []
    # color = []
    # for i in range(1, self.num_connections+1):
    #     (row, col) = np.where(self.CM==i)
    #     (row_idx, col_idx) = (row[0]+1, col[0]+1)
    #     edges.append((row_idx, col_idx))
    #     if row_idx <= self.num_sources:
    #         color.append('red')
    #     else:
    #         color.append('blue')
    #     end
    # end

    # G = nx.Graph(edges)

    # color_map = []

    # for node in G:
    #     if node <= self.num_sources:
    #         color_map.append("red")
    #     else:
    #         color_map.append("lightblue")
    #     end
    # end

    # nx.draw(G, node_color=color_map, with_labels = True)
    # plt.show()

    # pass
end
