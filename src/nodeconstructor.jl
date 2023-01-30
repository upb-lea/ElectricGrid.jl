

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
    L2L_p
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
        verbosity = 0
        )

Create a mutable struct NodeConstructor, which serves as a basis for the creation of an
energy grid: `num_sources` corresponse to the amount of sources and `num_loads` is the
amount of loads in the grid. `CM` is the connection matrix which indicates how the elements
in the grid are connected to each other. To specify the elements of the net in more detail,
values for the elements can be passed via `parameters`. If no connection matrix is entered,
it can be generated automatically. `S2S_p` is the probability that a source is connected
to another source and `S2L_p` is the probability that a source is connected to a load.
"""
function NodeConstructor(;num_sources, num_loads, CM=nothing, parameters=nothing,
                        S2S_p=0.1, S2L_p=0.8, L2L_p=0.3, verbosity=0)

    tot_ele = num_sources + num_loads

    cntr = 0
    num_connections = 0

    if CM === nothing
        cntr, CM = CM_generate(num_sources, num_loads, S2L_p, S2S_p, L2L_p)
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

        parameters = check_parameters(parameters, num_sources, num_loads, num_connections, CM, verbosity)  # Checks if all entries are given, if not, fills up with random values

        @assert length(keys(parameters)) == 4 "Expect parameters to have the four entries 'cable', 'load', 'grid' and 'source' but got $(keys(parameters))"

        @assert length(keys(parameters["grid"])) == 8 "Expect parameters['grid'] to have the 8 entries 'fs', 'v_rms', 'phase' and 'f_grid' but got $(keys(parameters["grid"]))"

        @assert length(parameters["source"]) == num_sources "Expect the number of sources to match the number of sources in the parameters, but got $num_sources and $(length(parameters["source"]))"

        @assert length(parameters["load"]) == num_loads "Expect the number of loads to match the number of loads in the parameters, but got $num_loads and $(length(parameters["load"]))"

        @assert length(parameters["cable"]) == num_connections "Expect the number of cables to match the number of cables in the parameters, but got $num_connections and $(length(parameters["cable"]))"

        # println(parameters)

        num_fltr_LCL, num_fltr_LC, num_fltr_L = cntr_fltrs(parameters["source"])
        num_loads_RLC, num_loads_LC, num_loads_RL, num_loads_RC, num_loads_L, num_loads_C, num_loads_R = cntr_loads(parameters["load"])
        

        @assert num_fltr_LCL + num_fltr_LC + num_fltr_L == num_sources "Expect the number of sources to be identical to the sum of the filter types, but the number of sources is $num_sources and the sum of the filters is $(num_fltr_LCL + num_fltr_LC + num_fltr_L)"

        @assert num_loads_RLC + num_loads_LC + num_loads_RL + num_loads_RC + num_loads_R + num_loads_C + num_loads_L == num_loads "Expect the number of loads to be identical to the sum of the loads types, but the number of loads is $num_loads and the sum of the loads is $(num_loads_RLC + num_loads_RL + num_loads_RC + num_loads_R)"
        
        # valid_realistic_para(parameters) TODO: revise that values are checked independently of network parameters, if necessary only via A-matrix

    else
        throw("Expect parameters to be a dict or nothing, not $(typeof(parameters))")
    end

    num_fltr = 4 * num_fltr_LCL + 3 * num_fltr_LC + 2 * num_fltr_L
    # Equivalent to the number of load states
    num_impedance = (2 * (num_loads_RLC
                        + num_loads_LC
                        + num_loads_RL
                        + num_loads_L)
                        + num_loads_RC + num_loads_C + num_loads_R)

    num_spp = num_fltr_LCL * 4 + num_fltr_LC * 3 + num_fltr_L * 2 + num_connections + (num_loads_RLC + num_loads_LC + num_loads_RL + num_loads_L) * 2 + (num_loads_RC + num_loads_C + num_loads_R)

    NodeConstructor(num_connections, num_sources, num_loads, num_fltr_LCL, num_fltr_LC, num_fltr_L, num_loads_RLC, num_loads_LC, num_loads_RL, num_loads_RC, num_loads_L, num_loads_C, num_loads_R, num_impedance, num_fltr, num_spp, cntr, tot_ele, CM, parameters, S2S_p, S2L_p, L2L_p, verbosity)
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

function check_parameters(parameters, num_sources, num_loads, num_connections, CM, verbosity)

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
        grid_properties["fs"] =  1e-4 # TODO: this should be 1/env.ts
        grid_properties["v_rms"] = 230
        grid_properties["phase"] = 3
        grid_properties["f_grid"] = 50
        grid_properties["Δfmax"] = 0.005 # Hz # The drop in frequency, Hz, which will cause a 100% increase in active power
        grid_properties["ΔEmax"] = 5/100 # V # The drop in rms voltage, which will cause a 100% decrease in reactive power
        grid_properties["ramp_end"] = 2/50
        grid_properties["process_start"] = 2/50
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
        if !haskey(parameters["grid"], "f_grid")
            parameters["grid"]["f_grid"] = 50.0
        end
        if !haskey(parameters["grid"], "Δfmax")
            parameters["grid"]["Δfmax"] = 0.005
        end
        if !haskey(parameters["grid"], "ΔEmax")
            parameters["grid"]["ΔEmax"] = 5/100
        end
        if !haskey(parameters["grid"], "ramp_end")
            parameters["grid"]["ramp_end"] = 2/parameters["grid"]["f_grid"]
        end
        if !haskey(parameters["grid"], "process_start")
            parameters["grid"]["process_start"] = 2/parameters["grid"]["f_grid"]
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
            @warn "The number of defined sources $num_def_sources is smaller than the number specified sources in the environment $num_sources, therefore the remaining $num_undef_sources sources are selected randomly!"
        end

        num_LC_defined = 0
        num_LCL_defined = 0
        source_type_fixed = 0

        for (index, source) in enumerate(parameters["source"])
            if !haskey(source, "pwr")
                source["pwr"] = rand(range(start=5,step=5,stop=50))*1e3
            end
            
            if !haskey(source, "vdc")
                source["vdc"] = 800#rand(range(start=690,step=10,stop=800)) # if randomized then the classic controllers go unstable - maybe range is too wide
            end

            if !haskey(source, "i_rip")
                source["i_rip"] = 0.15#rand(Uniform(0.1, 0.15)) # no reason to randomize
            end
                
            if !haskey(source, "v_rip")                    
                source["v_rip"] = 0.01537#rand(Uniform(0.014, 0.016)) # no reason to randomize
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
                source["i_limit"]= 1.15*Iop + ΔIlfmax
            end

            if !haskey(source, "R1")
                
                #= Practical Example:
                L_filter = 70e-6
                R_filter = 1.1e-3
                R_filter_C = 7e-3
                C_filter = 250e-6
                =#

                source["R1"] = 200 * source["L1"] # can be as low as 15
                
            end
            
            if !haskey(source, "fltr")
                source["fltr"] = "LCL"#rand(["LC", "LCL"]) 
            elseif !(source["fltr"] in ["L", "LC", "LCL"])
                # TODO: Raise warning: False key
                source["fltr"] = "L"
                @warn "filterType not known! set to L filter, please choose L, LC, or LCL!"
            end
            
            if (source["fltr"] == "LC" || source["fltr"] == "LCL")
                #Capacitor design
                if source["fltr"] == "LC"
                    num_LC_defined += 1
                end
                if source["fltr"] == "LCL"
                    num_LCL_defined += 1
                end
                if !haskey(source, "C")

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

                if source["fltr"] == "LC" && !haskey(source, "R_C")
                    source["R_C"] = 28*source["C"] # TODO: actually design the damping resistance
                end  
            end

            if source["fltr"] == "LC" && 1/sqrt(source["L1"]*source["C"]) > parameters["grid"]["fs"]/2
                if verbosity > 0
                    @warn ("The LC filter parameters have been poorly chosen.
                    The filtering capacitors should be chosen such that the resonant 
                    frequency 1/sqrt(L*C) is approximately sqrt(ωn * ωs), where ωn 
                    is the angular frequency of the grid, and ωs is the angular 
                    switching frequency.")
                end

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
                
                source["v_limit"]= 1.5*2*(sqrt(2)*parameters["grid"]["v_rms"] + ΔVcfmax)
            end

            if  source["fltr"] == "LCL" && !haskey(source, "L2")

                #= Theory:
                    The attenuation of the LCL-filter is 60db/decade for frequencies above the 
                    resonant frequency, therefore lower switching frequency for the converter 
                    can be used. It also provides better decoupling between the filter and the 
                    grid impedance and lower current ripple accross the grid inductor. 

                    The LCL filter has good current ripple attenuattion even with small inductance 
                    values. However it can bring also resonances and unstable states into the system. 
                    Therefore the filter must be designed precisely accoring to the parameters of the 
                    specific converter. 

                    The importatnt parameter of the filter is its cut-off frequency. The cut-off 
                    frequency of the filter must be minimally one half of the switching frequency of 
                    the converter, because the filter must have enough attenuation in the range of the 
                    converter's switching frequency.

                    The LCL filter will be vulnerable to oscillation too and it will magnify frequencies 
                    around its cut-off frequency. Therefore the filter is added with damping. The simplest 
                    way is to add damping resistor. In general there are four possible places where the 
                    resistor can be placed series/parallel to the inverter side inductor or series/parallel 
                    to the filter capacitor.

                    The variant with the resistor connected in series with the filter capacitor has been 
                    chosen. The peak near the resonant frequency can be significantly attenuated. This is 
                    a simple and reliable solution, but it increases the heat losses in the system and it 
                    greatly decreases the efficiency of the filter. This problem can be solved by active 
                    damping. Such a resistor reduces the voltage across the capacitor by a voltage proportional 
                    to the current that flows through it. This can be also done in the control loop. The 
                    current through the filter capacitor is measured and differentiatbed by the term 
                    (s*Cf*R_C). A real resistor is not used and the calculated value is subtracted from the 
                    demanded current. In this way the filter is actively damped with a virtual resistor 
                    without any losses. The disadvantage of this method is that an additional current sensor 
                    is required and the differentiator may bring noise problems because it amplifies high 
                    frequency signals. 

                =#

                #TODO: add user warnings if the L, C, parameters they choose are stupid  (more than 0.5*fs)  
                # also calculate the maximal power factor variation. If more than 5% add warning             

                fc = parameters["grid"]["fs"]/5
                ωc = 2π*fc

                if !haskey(source, "L2")

                    source["L2"] = source["L1"]/(ωc^2*source["L1"]*source["C"] - 1)

                    if source["L2"] < 0
                        source["L2"] == 1e-6
                        # TODO: add warning - the user choose bad ripple values
                    end
                end
                
                if !haskey(source, "R2")

                    source["R2"] = 200 * source["L2"]
                end  

                if !haskey(source, "R_C")

                    source["R_C"] = 1/(3*ωc*source["C"])
                end 

            end

            if source["fltr"] == "LCL"

                fc = (1/2π)*sqrt((source["L1"] + source["L2"])/(source["L1"]*source["L2"]*source["C"]))

                if fc > parameters["grid"]["fs"]/2
                    if verbosity > 0
                        @warn ("The LCL filter parameters have been poorly chosen.
                        The cut-off frequency of the filter must be minimally one half 
                        of the switching frequency of the converter, because the filter 
                        must have enough attenuation in the range of the converter's 
                        switching frequency.")
                    end
                end
            end

            if !haskey(source, "source_type")
                source["source_type"] = "ideal"
                source_type_fixed += 1
            end

            if !haskey(source, "τv")
                source["τv"] = 0.002 # time constant of the voltage loop # 0.02
            end

            if !haskey(source, "τf")
                source["τf"] = 0.002 # time constant of the frequency loop # 0.002
            end

            if !haskey(source, "pf") # power factor

                default_pf = 0.8

                if !haskey(source, "p_set") && !haskey(source, "q_set")

                    source["pf"] = default_pf # power factor

                elseif haskey(source, "q_set") && !haskey(source, "p_set")

                    p_set = sqrt(source["pwr"]^2 - source["q_set"]^2)
                    source["pf"] = p_set/source["pwr"]

                elseif haskey(source, "p_set") && !haskey(source, "q_set")

                    source["pf"] = source["p_set"]/source["pwr"]

                elseif haskey(source, "p_set") && haskey(source, "q_set")

                    s_set = sqrt(source["p_set"]^2 + source["q_set"]^2)*sign(source["p_set"]*source["q_set"])

                    if s_set == 0
                        source["pf"] = 1/sqrt(2)
                    else
                        source["pf"] = source["p_set"]/s_set
                    end
                end
            end

            if !haskey(source, "p_set")
                source["p_set"] = source["pwr"]*source["pf"]
            end

            if !haskey(source, "q_set")
                source["q_set"] = sqrt(source["pwr"]^2 - source["p_set"]^2)
            end

            if !haskey(source, "v_pu_set")
                source["v_pu_set"] = 1.0
            end

            if !haskey(source, "v_δ_set")
                source["v_δ_set"] = 0.0
            end

            if !haskey(source, "mode")
                source["mode"] = "Droop"
            end

            if !haskey(source, "control_type")
                source["control_type"] = "classic"
            end

            if !haskey(source, "γ") # asymptotic mean
                source["γ"] = source["p_set"]
            end

            if !haskey(source, "std_asy") || haskey(source, "κ")# asymptotic standard deviation

                #std_asy = sqrt(σ^2/(2*κ)) # asymptotic standard deviation
                if !haskey(source, "σ")

                    source["std_asy"] = 0.0
                elseif !haskey(source, "κ")

                    source["std_asy"] = source["pwr"]/4
                else

                    source["std_asy"] = source["σ"]/sqrt(2*source["κ"])
                end
            end

            if !haskey(source, "κ") # mean reversion parameter

                if source["std_asy"] == 0.0

                    source["κ"] = 0.0
                else

                    source["κ"] = source["σ"]^2/(2*source["std_asy"]^2)
                end
            end

            if !haskey(source, "σ") # Brownian motion scale i.e. ∝ diffusion parameter
                source["σ"] = 0.0
            end

            if !haskey(source, "X₀") # initial values
                source["X₀"] = source["p_set"]
            end

            if !haskey(source, "Δt") # time step

                steps = 1 # ... cycles for 1 step
                source["Δt"] = round(steps*parameters["grid"]["fs"]/(parameters["grid"]["f_grid"]))/parameters["grid"]["fs"]

            elseif haskey(source, "Δt")

                if typeof(source["Δt"]) == Int

                    steps = source["Δt"] # ... cycles for 1 step
                    source["Δt"] = round(steps*parameters["grid"]["fs"]/(parameters["grid"]["f_grid"]))/parameters["grid"]["fs"]
                else
                    source["Δt"] = round(source["Δt"]*(parameters["grid"]["fs"]))/parameters["grid"]["fs"]
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
            # Validierung ob LC vorhanden ist?
            if num_LC_defined == 0 && num_fltr_LC_undef == 0 # What is this? What if the user defined an L or LCL filter
                #@warn "Bla Bla Bla Bla Bla Bla .... My name is Plop. No LC filter defined/set random, if wanted please set in parameter dict!"
            end
        else

            if num_LC_defined == 0 
                #@warn "Bla Bla Bla Bla Bla Bla .... My name is Plop. No LC filter defined/set random, if wanted please set in parameter dict!"
            end
        end

        #source_type_fixed > 0 && @warn "Wagga Wagga. Poopy-di scoop. Scoop-diddy-whoop. Whoop-di-scoop-di-poop. $source_type_fixed sourceType not defined! set to ideal! Why do I care??"

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
        num_def_loads = length(parameters["load"])

        num_undef_loads = num_loads - num_def_loads

        @assert num_undef_loads >= 0 "Expect the number of defined loads within the parameter dict to be less or equal to the number of loads in the env, but the entries within the parameter dict is $num_def_loads and the number of env loads is $num_loads."

        if num_undef_loads > 0
            @warn "The number of defined loads $num_def_loads is smaller than the number specified loads in the environment $num_loads, therefore the remaining $num_undef_loads loads are selected randomly!"
        end

        for (index, load) in enumerate(parameters["load"])
            if !haskey(load, "impedance")
                load["impedance"] = rand(["RLC","RL","RC","LC","R","L","C"])
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

        # add Z and pwr to the parameter_load dict; needed for solving the power flow equations
        # assuming all passive loads as parrallel connection of devices
        for (index, load) in enumerate(parameters["load"])
            if load["impedance"] == "R"
                load["Z"] = load["R"]
                load["pf"] = 1
            elseif load["impedance"] == "L"
                load["Z"] = 1im*2*pi*parameters["grid"]["f_grid"]*load["L"]
                load["pf"] = 0
            elseif load["impedance"] == "C"
                load["Z"] = 1/(1im*2*pi*parameters["grid"]["f_grid"]*load["C"])
                load["pf"] = 0
            elseif load["impedance"] == "RL"
                load["Z"] = 1im*parameters["grid"]["f_grid"]*2*pi*load["R"]*load["L"]/(load["R"]+1im*parameters["grid"]["f_grid"]*2*pi*load["L"])
                load["pf"] = cos(atan(load["R"]/(parameters["grid"]["f_grid"]*2*pi*load["L"])))
            elseif load["impedance"] == "RC"
                load["Z"] = load["R"]/(1+1im*parameters["grid"]["f_grid"]*2*pi*load["C"]*load["R"])
                load["pf"] = cos(-atan(load["R"]*parameters["grid"]["f_grid"]*2*pi*load["C"]))
            elseif load["impedance"] == "LC"
                load["Z"] = 1im*parameters["grid"]["f_grid"]*2*pi*load["L"]/(1-(parameters["grid"]["f_grid"]*2*pi)^2*load["L"]*load["C"])
            elseif load["impedance"] == "RLC"
                load["Z"] = 1im*parameters["grid"]["f_grid"]*2*pi*load["L"]/(1+1im*parameters["grid"]["f_grid"]*2*pi*load["L"]/load["R"]-(parameters["grid"]["f_grid"]*2*pi)^2*load["L"]*load["C"])
                # TODO PF RLC parallel             
            end
            
            load["pwr"] = parameters["grid"]["v_rms"]^2 / abs(load["Z"]) * parameters["grid"]["phase"]
            #println(load["pf"])
        end
    end


    ################
    # CHECK CABLES #
    ################
    
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
            @warn "The number of defined cables $num_def_cables is smaller than the number specified cables in the environment $num_connections, therefore the remaining $num_undef_cables cables are selected randomly!"
        end
        cable_from_pfe_idx = []
        for (idx, cable) in enumerate(parameters["cable"])
            # solve powerflow equation (pfe) only if needed - but if not all values are give - take all values from pfe and overwrite the rest
            
            if !haskey(cable, "len")
                cable["len"] = rand(Uniform(1e-3, 1e1))
            end
            #=
            if !haskey(cable, "Rb")
                cable["Rb"] = 0.722 # TODO: Fixed?!
            end

            if !haskey(cable, "Cb")
                cable["Cb"] = 0.4e-6 # TODO: Fixed?!
            end

            if !haskey(cable, "Lb")
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
            =#
            if !haskey(cable, "R") | !haskey(cable, "L") | !haskey(cable, "C")
                @info "Parameters from cable $(idx) missing. All cable parameters are calculate based on power flow equation. Create a counter - we don't want to see this every time."
                push!(cable_from_pfe_idx, idx)
            end

        end

        if ! isempty(cable_from_pfe_idx)
            println("START PFE")
            parameters = layout_cabels(CM, num_sources, num_loads, parameters)
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
function generate_parameters(num_fltr_LCL, num_fltr_LC, num_fltr_L, num_connections,
                             num_loads_RLC, num_loads_LC, num_loads_RL, num_loads_RC,
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
    #TODO: why are these things randomized again?? - maybe I'm not following the code, but surely these have been randomized if the user did not define them
    source["pwr"] = rand(range(start=5,step=5,stop=50))*1e3
    source["vdc"] = 800 #rand(range(start=690,step=10,stop=800))
    source["i_rip"] = 0.15 #rand(Uniform(0.1, 0.15))
    source["v_rip"] = 0.01537 #rand(Uniform(0.014, 0.016))

   #Inductor design

   Vorms = grid_properties["v_rms"]*1.05
   Vop = Vorms*sqrt(2)

   Zl = 3*Vorms*Vorms/source["pwr"]

   Iorms = Vorms/Zl
   Iop = Iorms*sqrt(2)

   ΔIlfmax = source["i_rip"]*Iop

   source["L1"] = 0.5*(source["vdc"]*(4*grid_properties["fs"]*ΔIlfmax)^-1)
   source["L2"] = deepcopy(source["L1"]) 
   source["i_limit"]= Iop + ΔIlfmax

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
   source["v_limit"]= 2*(sqrt(2)*grid_properties["v_rms"] + ΔVcfmax)
   
   

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
    source["vdc"] = 800 #rand(range(start=690,step=10,stop=800))
    source["i_rip"] = 0.15 #rand(Uniform(0.1, 0.15))
    source["v_rip"] = 0.01537 #rand(Uniform(0.014, 0.016))

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
    source["v_limit"]= 2*(sqrt(2)*grid_properties["v_rms"] + ΔVcfmax)

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
    source["vdc"] = 800#rand(range(start=690,step=10,stop=800))
    source["i_rip"] = 0.15#rand(Uniform(0.1, 0.15))
    source["v_rip"] = 0.01537#rand(Uniform(0.014, 0.016))

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
    # TODO: do these values make sence?!, if smaller 1, take out round command!!
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
    
    cable = Dict()
    cable["len"] = 1.0#rand(Uniform(1e-3, 1e1))

    cable["Rb"] =  0.722
    cable["Cb"] = 0.4e-6
    cable["Lb"] =  0.264e-3

    cable["R"] = cable["len"] * cable["Rb"]
    cable["L"] = cable["len"] * cable["Lb"]
    cable["C"] = cable["len"] * cable["Cb"]

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
    CM_generate(num_sources, num_loads, S2L_p, S2S_p)

Returns the constructed CM and the total number of connections.
"""
function CM_generate(num_sources, num_loads,  S2L_p, S2S_p, L2L_p)

    # counting the connections 
    cntr = 0

    # get total elements
    tot_ele = num_sources + num_loads
    
    # get a upper triangular matrix
    mask = UpperTriangular(ones(tot_ele, tot_ele))
    CM = rand(tot_ele,tot_ele) .* mask # fill matrix with random entries between [0,1]
    CM = CM - Diagonal(CM) # delet diagonal bc no connection with itself

    # go through the matrix
    # -1 bc last entry is 0 anyway
    for i in 1:tot_ele-1
        # start at i, bc we need to check only upper triangle
        if i <= num_sources
            for j in i:tot_ele-1
                if j > num_sources-1  # select propability according to column
                    cntr, x = tobe_or_n2b(cntr, CM[i, j+1], S2L_p)
                    CM[i, j+1] = x
                else
                    cntr, x = tobe_or_n2b(cntr, CM[i, j+1], S2S_p)
                    CM[i, j+1] = x
                end
            end
        else
            for j in i:tot_ele-1
                cntr, x = tobe_or_n2b(cntr, CM[i, j+1], L2L_p)
                CM[i, j+1] = x
            end
        end
    end

    # make sure that no objects disappear or subnets are formed
    if S2L_p < 1

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

                idx_rnd = sample(1:length(idx_list), samples, replace=false) # draw samples from the list
                
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
    get_A_src(self::NodeConstructor, source_i) 

Create the A_src entry for a source in the A matrix.
"""
function get_A_src(self::NodeConstructor, source_i)

    parameter_i = self.parameters["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        A_src = zeros(4,4)
        A_src[1,1] = -(parameter_i["R1"]+parameter_i["R_C"])/parameter_i["L1"]
        A_src[1,2] = -1/parameter_i["L1"]
        A_src[1,3] = parameter_i["R_C"]/parameter_i["L1"]
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
        A_src[1,1] = -(parameter_i["R1"])/parameter_i["L1"]
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
        A_src[3,2] = 1/parameter_i["R_C"]*(C_sum)^(-1)
        A_src[3,3] = -1/parameter_i["R_C"]*(C_sum)^(-1)

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
    
    # TODO: maybe do in one loop like: for (i, vals) in enumerate(zip([1, 4, 2, 5], 2:12, (:a, :b, :c)))
    cntr = 0
    for (i, ele) in enumerate(A_src_list)
        #if i <= self.num_fltr_LCL
        
        if self.parameters["source"][i]["fltr"]  == "LCL"
            start = 1 + cntr
            stop = 4 + cntr
            cntr += 4
            A_src[start:stop, start:stop] = ele

        elseif self.parameters["source"][i]["fltr"] == "LC"
            start = 1 + cntr
            stop = 3 + cntr
            cntr += 3
            A_src[start:stop, start:stop] = ele

        elseif self.parameters["source"][i]["fltr"]  == "L"
            start = 1 + cntr
            stop = 2 + cntr
            cntr += 2
            A_src[start:stop,start:stop] = ele
        end
    end
    
    # get A_src_trn_c
    A_src_trn_c = zeros(self.num_fltr, self.num_connections)
    A_src_trn_c_list = [get_A_src_trn_c(self, i) for i in 1:self.num_sources] # start at 1 bc Source 1 ...
    cntr = 0
    for (i, ele) in enumerate(A_src_trn_c_list)
        if self.parameters["source"][i]["fltr"]  == "LCL"
            start = 1 + cntr
            stop = 4 + cntr
            cntr += 4
            A_src_trn_c[start:stop,:] = ele

        elseif self.parameters["source"][i]["fltr"] == "LC"
            start = 1 + cntr
            stop = 3 + cntr
            cntr += 3
            A_src_trn_c[start:stop,:] = ele

        elseif self.parameters["source"][i]["fltr"]  == "L"
            start = 1 + cntr
            stop = 2 + cntr
            cntr += 2
            A_src_trn_c[start:stop,:] = ele
        end
    end

    # get A_src_trn_l
    A_src_trn_l = zeros(self.num_connections, self.num_fltr)
    A_src_trn_l_list = [get_A_src_trn_l(self, i) for i in 1:self.num_sources] # start at 1 bc Source 1 ...

    cntr = 0
    for (i, ele) in enumerate(A_src_trn_l_list)
        if self.parameters["source"][i]["fltr"]  == "LCL"
            start = 1 + cntr
            stop = 4 + cntr
            cntr += 4
            A_src_trn_l[:,start:stop] = ele

        elseif self.parameters["source"][i]["fltr"] == "LC"
            start = 1 + cntr
            stop = 3 + cntr
            cntr += 3
            A_src_trn_l[:,start:stop] = ele

        elseif self.parameters["source"][i]["fltr"]  == "L"
            start = 1 + cntr
            stop = 2 + cntr
            cntr += 2
            A_src_trn_l[:,start:stop] = ele
        end
    end

    A_trn = generate_A_trn(self)

    if self.num_loads > 0
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

        A_load_diag = zeros(self.num_impedance, self.num_impedance)
        A_load_list = [get_A_load(self, i) for i in 1:self.num_loads]

        cntr = 0
        for (i, ele) in enumerate(A_load_list)
            if (self.parameters["load"][i]["impedance"]  == "RLC" ||
                self.parameters["load"][i]["impedance"]  == "LC" ||
                self.parameters["load"][i]["impedance"]  == "RL" ||
                self.parameters["load"][i]["impedance"]  == "L")

                start = 1 + cntr
                stop = 2 + cntr
                cntr += 2
                A_load_diag[start:stop, start:stop] = ele

            elseif (self.parameters["load"][i]["impedance"]  == "RC" ||
                self.parameters["load"][i]["impedance"]  == "C" ||
                self.parameters["load"][i]["impedance"]  == "R")

                start = 1 + cntr
                stop = 1 + cntr
                cntr += 1
                start = i + self.num_loads_RLC + self.num_loads_LC + self.num_loads_RL + self.num_loads_L
                A_load_diag[start:start, start:start] = ele

            end
        end

        A_load_zeros = zeros(self.num_fltr, self.num_impedance)
        A_load_zeros_t = A_load_zeros'
    end

    

    if self.num_loads == 0
        A = [A_src       A_src_trn_c        
        A_src_trn_l           A_trn ]
    else
        A = [A_src       A_src_trn_c        A_load_zeros
        A_src_trn_l           A_trn  A_tran_load_l
        A_load_zeros_t  A_tran_load_c   A_load_diag]
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
    cntr = 0
    for (i, ele) in enumerate(B_source_list)
        if self.parameters["source"][i]["fltr"]  == "LCL"
            start = 1 + cntr
            stop = 4 + cntr
            cntr += 4
            B[start:stop,i] = ele

        elseif self.parameters["source"][i]["fltr"]  == "LC"
            start = 1 + cntr
            stop = 3 + cntr
            cntr += 3
            B[start:stop,i] = ele

        elseif self.parameters["source"][i]["fltr"]  == "L"
            start = 1 + cntr
            stop = 2 + cntr
            cntr += 2
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
        if self.parameters["source"][s]["fltr"]  == "LCL"
            push!(states, "source$s"*"_i_L1")    # i_f1; dann i_f2....
            push!(states, "source$s"*"_v_C_filt")
            push!(states, "source$s"*"_i_L2")
            push!(states, "source$s"*"_v_C_cables")
        
        elseif self.parameters["source"][s]["fltr"]  == "LC"
            push!(states, "source$s"*"_i_L1")
            push!(states, "source$s"*"_v_C_filt")
            push!(states, "source$s"*"_v_C_cables")

        elseif self.parameters["source"][s]["fltr"]  == "L"
            push!(states, "source$s"*"_i_L1")
            push!(states, "source$s"*"_v_C_cables")
        end
    end

    for c in 1:self.num_connections
        # TODO: correct?
        push!(states, "cable$c"*"_i_L")
    end

    for l in 1:self.num_loads
        if (self.parameters["load"][l]["impedance"]  == "RLC" ||
            self.parameters["load"][l]["impedance"]  == "LC" ||
            self.parameters["load"][l]["impedance"]  == "RL" ||
            self.parameters["load"][l]["impedance"]  == "L")
            push!(states, "load$l"*"_v_C_total")
            push!(states, "load$l"*"_i_L")
        elseif (self.parameters["load"][l]["impedance"]  == "RC" ||
            self.parameters["load"][l]["impedance"]  == "C" ||
            self.parameters["load"][l]["impedance"]  == "R")
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
        if self.parameters["source"][s]["fltr"]  == "LCL"
            push!(state_paras, self.parameters["source"][s]["L1"]) 
            push!(state_paras, self.parameters["source"][s]["C"]) 
            push!(state_paras, self.parameters["source"][s]["L2"])
            push!(state_paras, get_C_sum_cable_node(s,self))
        
        elseif self.parameters["source"][s]["fltr"]  == "LC"
            push!(state_paras, self.parameters["source"][s]["L1"])
            push!(state_paras, self.parameters["source"][s]["C"])
            push!(state_paras, get_C_sum_cable_node(s,self))  
        
        elseif self.parameters["source"][s]["fltr"]  == "L"
            push!(state_paras, self.parameters["source"][s]["L1"])
            push!(state_paras, get_C_sum_cable_node(s,self))
        end
    end

    for c in 1:self.num_connections
        push!(state_paras, self.parameters["cable"][c]["L"]) 
    end

    for l in 1:self.num_loads
        if (self.parameters["load"][l]["impedance"]  == "RLC" ||
            self.parameters["load"][l]["impedance"]  == "LC" ||
            self.parameters["load"][l]["impedance"]  == "RL" ||
            self.parameters["load"][l]["impedance"]  == "L")
            c=0
                if haskey(self.parameters["load"][l], "C")
                    c= self.parameters["load"][l]["C"]
                end
            push!(state_paras, get_C_sum_cable_node(self.num_sources+l,self) + c) 
            push!(state_paras, self.parameters["load"][l]["L"])  
        elseif (self.parameters["load"][l]["impedance"]  == "RC" ||
            self.parameters["load"][l]["impedance"]  == "C" ||
            self.parameters["load"][l]["impedance"]  == "R")
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
        if self.parameters["source"][s]["fltr"]  == "LCL"
            push!(actions, "source$s"*"_u")
        
        elseif self.parameters["source"][s]["fltr"]  == "LC"
            push!(actions, "source$s"*"_u")
        
        elseif self.parameters["source"][s]["fltr"]  == "L"
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
        text = node_descriptions,
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

function get_Y_bus(self::NodeConstructor)
    """
    get_Y_bus(self::NodeConstructor)
    Returns the Admittance Matrix of the power grid based on the CM Matrix and the parameter dict 
    """

    Y_bus = zeros(Complex{Float64}, self.tot_ele, self.tot_ele)  # Y -> tot_ele x tot_ele
    omega = 2*π*self.parameters["grid"]["f_grid"]

    for col in 1:self.tot_ele
        for row in 1:self.tot_ele
            if self.CM[row, col] != 0  # we have a cable connected
                # CM index defines the number of the cable
                cable_idx = abs(Int(self.CM[row, col]))
                G_RL = self.parameters["cable"][cable_idx]["R"] / (self.parameters["cable"][cable_idx]["R"]^2 + omega^2 *self.parameters["cable"][1]["L"]^2)
                B_RL = (omega *self.parameters["cable"][cable_idx]["L"])/ (self.parameters["cable"][cable_idx]["R"]^2 + omega^2 *self.parameters["cable"][1]["L"]^2)
                Y_bus[row, col] = -G_RL - im*B_RL
            elseif row == col  # diagonal elements
                # Go through all col elements of that row to find the connected cable to that bus (non zero elements in CM[:, row])
                cable_idxs = filter(n -> n !=0, self.CM[row, :])
                G = 0
                B = 0
                for idx in cable_idxs
                    # add all RL 
                    idx = abs(Int(idx))
                    G += self.parameters["cable"][idx]["R"] / (self.parameters["cable"][idx]["R"]^2 + omega^2 *self.parameters["cable"][1]["L"]^2)
                    B += (omega *self.parameters["cable"][idx]["L"])/ (self.parameters["cable"][idx]["R"]^2 + omega^2 *self.parameters["cable"][1]["L"]^2)
                    # and add all shunt C connected to that bus since diagonal element
                    B += omega * self.parameters["cable"][idx]["C"] / 2 
                end
                Y_bus[row, col] = G + im*B
            end
        end
    end
    return Y_bus
end

function Source_Setup(num_sources; random = nothing, awg_pwr = 200e3, mode = 3)

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
        pwrs = rand(Uniform(0.5*awg_pwr, 1.5*awg_pwr), num_sources)
    end

    for i in 1:num_sources

        source = Dict()

        if random == 0 || isnothing(random)

            pwr = 200e3

            source["mode"]     = mode

            source["fltr"]     = "LCL"  # Filter type

            pwr = awg_pwr
            source["pwr"]      = pwr # Rated Apparent Power, VA
            source["p_set"]    = 0   # Real Power Set Point, Watt
            source["q_set"]    = 0   # Imaginary Power Set Point, VAi

            source["v_pu_set"] = 1.00   # Voltage Set Point, p.u.
            source["v_δ_set"]  = 0      # Voltage Angle, degrees

            source["τv"]       = 0.002  # Time constant of the voltage loop, seconds
            source["τf"]       = 0.002  # Time constant of the frequency loop, seconds

            source["Observer"] = true   # Discrete Luenberger Observer

        else

            source["mode"]     = mode

            source["fltr"]     = "LCL"  # Filter type

            pwr = pwrs[i]
            source["pwr"]      = pwr  # Rated Apparent Power, VA
            source["p_set"]    = 0   # Real Power Set Point, Watt
            source["q_set"]    = 0   # Imaginary Power Set Point, VAi

            source["τv"]       = 0.002  # Time constant of the voltage loop, seconds
            source["τf"]       = 0.002  # Time constant of the frequency loop, seconds

            source["Observer"] = true   # Discrete Luenberger Observer

            source["v_pu_set"] = 1.00   # Voltage Set Point, p.u.
            source["v_δ_set"]  = 0      # Voltage Angle, degrees

            if random == 2

                source["std_asy"]  = pwr/8   # Asymptotic Standard Deviation
                source["σ"]        = pwr/8   # Brownian motion scale i.e. ∝ diffusion, volatility parameter
                source["Δt"]       = 0.1   # Time Step, seconds
                #source["X₀"]       = 0      # Initial Process Values, Watt
                source["k"]        = 1      # Interpolation degree
                source["γ"]        = 0      # asymptotic mean

            end

        end

        total_gen += pwr

        push!(source_list, source)

    end

    return source_list, total_gen
end

function Load_Setup(num_loads, total_gen; gen_load_ratio = 6, random = nothing, Vrms = 230)

    load_list = []

    avg_load = total_gen/(num_loads*gen_load_ratio)

    if random != 0 && !isnothing(random)

        Random.seed!(1234)
        pwrs = rand(Uniform(0.5*avg_load, 1.5*avg_load), num_loads)
        Random.seed!(1234)
        pfs = rand(Uniform(0.7, 1.0), num_loads)
    end

    for i in 1:num_loads

        load = Dict()

        if random == 0 || isnothing(random)

            R_load, L_load, _, _ = Parallel_Load_Impedance(avg_load, 0.6, Vrms)

            load["impedance"] = "RL"
            load["R"] = R_load
            load["L"] = L_load
            load["S"] = avg_load

        else

            R_load, L_load, _, _ = Parallel_Load_Impedance(pwrs[i], pfs[i], Vrms)

            load["impedance"] = "RL"
            load["R"] = R_load
            load["L"] = L_load
            load["S"] = pwrs[i]

        end

        push!(load_list, load)

    end

    return load_list
end

function Cable_Length_Setup(num_cables; random = 0, avg_length = 1.0)

    cable_list = []

    if random != 0 && !isnothing(random)

        Random.seed!(1234)
        lengths = rand(Uniform(0.5*avg_length, 1.5*avg_length), num_cables)
    end

    for i in 1:num_cables

        cable = Dict()
        
        if random == 0

            cable["len"]     = 1   # km
            cable["R"]       = 0.208   # Ω, line resistance
            cable["L"]       = 0.00025 # H, line inductance
            cable["C"]       = 0.4e-3  # F, line capacitance
            cable["i_limit"] = 10e12   # A, line current limit

        else

            cable["len"] = lengths[i]
        end

        push!(cable_list, cable)

    end

    return cable_list

end
