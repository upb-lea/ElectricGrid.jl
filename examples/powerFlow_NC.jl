using JuMP
import Ipopt

#using dare   warum geht das nicht?
using DrWatson
@quickactivate "dare"

using ReinforcementLearning
using PlotlyJS

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("Classical_Control.jl"))
include(srcdir("Power_System_Theory.jl"))
include(srcdir("MultiAgentGridController.jl"))


function set_bounds(variable, start_value, low_bound, up_bound)
    if !is_fixed(variable)
        set_lower_bound(variable, low_bound)
        set_upper_bound(variable, up_bound)
        set_start_value(variable, start_value)
    end
end



function get_degree(CM = CM) # how many cables are connected to a node? maybe remove function if not used
    result = zeros(Int, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = count(x -> x != 0, CM[i,:])
    end

    result
end

function get_cable_connections(CM = CM) # which cables are connected to which nodes

    result = Vector{Vector{Int64}}(undef, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = filter(x -> x != 0, abs.(CM[i,:]))
    end

    return result
end

function get_node_connections(CM = CM) # which nodes are connected to each other, including the self-connections

    result = Vector{Vector{Int64}}(undef, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = findall(x -> x != 0, CM[i,:])
        push!(result[i], i)
    end

    return result
end

function layout_cabels(CM, num_source, num_load, parameters)

    model = Model(Ipopt.Optimizer)
    #set_optimizer_attributes(model, "tol" => 1e-1)

    zero_expression = @NLexpression(model, 0.0)

    # Constant values
    omega = 2π*parameters["grid"]["fs"]

    # for every Source: v is fixed 230  #TODO: change depending on control mode?
    # for one Source: theta is fixed 0
    #TODO: user specified disctances
    #distances = [1.0 2.3 .323]  -> get from parameter dict

    num_nodes = num_source + num_load
    num_cables = maximum(CM)

    @variable(model, nodes[1 : num_nodes, ["v", "theta", "P", "Q"]])
   
    # cal total load[pwr]

    total_P_load = 0 #3000 / 3    #TODO get from parameter dict!
    total_Q_load = 0 #300 / 3
    total_S_source = 0 #3000 / 3

    for i = 1:num_nodes
        if i <= num_source
            total_S_source = total_S_source + parameters["source"][i]["pwr"]/parameters["grid"]["phase"]
            if parameters["source"][i]["mode"] in ["PQ Control", 3]
                if parameters["source"][i]["p_set"] < 0
                    total_P_load = total_P_load + parameters["source"][i]["p_set"]/parameters["grid"]["phase"]
                end
                if parameters["source"][i]["q_set"] < 0
                    total_Q_load = total_Q_load + parameters["source"][i]["q_set"]/parameters["grid"]["phase"]
                end
            elseif parameters["source"][i]["mode"] in ["PV Control", 4]
                if parameters["source"][i]["p_set"] < 0
                    total_P_load = total_P_load + parameters["source"][i]["p_set"]/parameters["grid"]["phase"]
                end
            end
        else
            total_P_load = total_P_load + parameters["load"][i-num_source]["pwr"]/parameters["grid"]["phase"]
            total_Q_load = total_Q_load + parameters["load"][i-num_source]["pwr"]/parameters["grid"]["phase"]
        end

    
    end
    println()
    println(total_P_load)
    println(total_Q_load)
    println(total_S_source)
    println()
    

    idx_p_mean_cal = []
    idx_q_mean_cal = []

    for i = 1:num_nodes
        if i <= num_source
            if parameters["source"][i]["control_type"] == "classic" 
                
                if parameters["source"][i]["mode"] in ["Swing", "Voltage Control", 1, 2] 

                    fix(nodes[i, "v"], parameters["source"][i]["v_pu_set"] * parameters["grid"]["v_rms"])
                    #fix(nodes[1, "theta"], 0.0) 
                    set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # question, does this limit too much, should be more or less.
                    
                    set_bounds(nodes[i, "P"], (total_P_load) / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # come from parameter dict/user?
                    set_bounds(nodes[i, "Q"], (total_Q_load) / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # P and Q are the average from power, excluding cable losses
                    push!(idx_p_mean_cal, i)
                    push!(idx_q_mean_cal, i)

                elseif parameters["source"][i]["mode"] in ["PQ Control", 3] 
                    fix(nodes[i, "P"], parameters["source"][i]["p_set"]/parameters["grid"]["phase"])
                    fix(nodes[i, "Q"], parameters["source"][i]["q_set"]/parameters["grid"]["phase"])
            
                    set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
                    set_bounds(nodes[i, "v"], parameters["grid"]["v_rms"], 0.95*parameters["grid"]["v_rms"], 1.05*parameters["grid"]["v_rms"])

                elseif parameters["source"][i]["mode"] in ["PV Control", 4] 
                    fix(nodes[i, "P"], parameters["source"][i]["p_set"]/parameters["grid"]["phase"])
                    fix(nodes[i, "v"], parameters["source"][i]["v_pu_set"] * parameters["grid"]["v_rms"])
                    set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
                    set_bounds(nodes[i, "Q"], (total_Q_load) / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # P and Q are the average from power, excluding cable losses
                    push!(idx_q_mean_cal, i)


                elseif parameters["source"][i]["mode"] in ["Semi-Synchronverter", 7] 
                    fix(nodes[i, "v"], parameters["source"][i]["v_pu_set"] * parameters["grid"]["v_rms"])
                    set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
                    set_bounds(nodes[i, "P"], total_P_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"])
                    set_bounds(nodes[i, "Q"], total_Q_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # P and Q are the average from power, excluding cable losses
                    push!(idx_p_mean_cal, i)
                    push!(idx_q_mean_cal, i)
                else
                    # all variable - 
                    set_bounds(nodes[i, "P"], total_P_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) 
                    set_bounds(nodes[i, "Q"], total_Q_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # P and Q are the average from power, excluding cable losses
                    set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
                    set_bounds(nodes[i, "v"], parameters["grid"]["v_rms"], 0.95*parameters["grid"]["v_rms"], 1.05*parameters["grid"]["v_rms"])
                    push!(idx_p_mean_cal, i)
                    push!(idx_q_mean_cal, i)
                end

            else
                # all variable - 
                set_bounds(nodes[i, "P"], total_P_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) 
                set_bounds(nodes[i, "Q"], total_Q_load / num_source, -parameters["source"][i]["pwr"]/parameters["grid"]["phase"], parameters["source"][i]["pwr"]/parameters["grid"]["phase"]) # P and Q are the average from power, excluding cable losses
                set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
                set_bounds(nodes[i, "v"], parameters["grid"]["v_rms"], 0.95*parameters["grid"]["v_rms"], 1.05*parameters["grid"]["v_rms"])
                push!(idx_p_mean_cal, i)
                push!(idx_q_mean_cal, i)
            end
            #end
        else
            println(-parameters["load"][i-num_source]["pwr"]/parameters["grid"]["phase"])
            fix(nodes[i, "P"], -parameters["load"][i-num_source]["pwr"]/parameters["grid"]["phase"])
            fix(nodes[i, "Q"], -parameters["load"][i-num_source]["pwr"]/parameters["grid"]["phase"])
    
            set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
            set_bounds(nodes[i, "v"], parameters["grid"]["v_rms"], 0.95*parameters["grid"]["v_rms"], 1.05*parameters["grid"]["v_rms"])
        end
    end

    cable_cons = get_cable_connections(CM)
    node_cons = get_node_connections(CM) 

    G = Array{NonlinearExpression, 2}(undef, num_nodes, num_nodes) # should be symmetric
    B = Array{NonlinearExpression, 2}(undef, num_nodes, num_nodes) # should be symmetric

    P_node = Array{NonlinearConstraintRef, 1}(undef, num_nodes)
    Q_node = Array{NonlinearConstraintRef, 1}(undef, num_nodes)

    # As radius goes down resistance goes up, inductance goes up, capacitance goes down. Put in formulas for this.
    #@variable(model, cables[1 : num_cables, ["L", "X_R", "C_L"]]) 
    @variable(model, cables[1 : num_cables, ["radius"]])

    L_cable = Array{NonlinearExpression, 1}(undef, num_cables)
    R_cable = Array{NonlinearExpression, 1}(undef, num_cables)
    C_cable = Array{NonlinearExpression, 1}(undef, num_cables)

    cable_conductance = Array{NonlinearExpression, 1}(undef, num_cables)
    cable_susceptance_0 = Array{NonlinearExpression, 1}(undef, num_cables) # diagonals - where we add capacitances
    cable_susceptance_1 = Array{NonlinearExpression, 1}(undef, num_cables) # off diagonals


    # fixed distance between cables
    D = 0.5 #m

    for i=1:num_cables
        

        set_bounds(cables[i, "radius"], (3e-3)/2, (2.05232e-3)/2, (4.1148e-3)/2) #m 
        # assumption to line to line
        println(parameters["cable"][i]["len"])
        L_cable[i] = @NLexpression(model, parameters["cable"][i]["len"] * 4e-7 * log(D/(0.7788 * cables[i, "radius"])))  # m* H/m

        # resistivity remains constant ρ_(T=50) = 1.973e-8 
        R_cable[i] = @NLexpression(model, parameters["cable"][i]["len"] * 1.973e-8 / (π * cables[i, "radius"]^2)) # m * Ω/m

        # X_R = 0.38 # ratio of omega*L/R
        # X_R = omega * L(r)/R(r)

        # line to neutral
        C_cable[i] = @NLexpression(model, parameters["cable"][i]["len"] * 2π * 8.854e-12 / (log(D/cables[i, "radius"]))) #m * F/m
    
        cable_conductance[i] = @NLexpression(model, (R_cable[i] / ((R_cable[i])^2 + (omega * L_cable[i])^2)))
        cable_susceptance_1[i] = @NLexpression(model, (-omega * L_cable[i] / (((R_cable[i])^2 + (omega * L_cable[i])^2))))
        cable_susceptance_0[i] = @NLexpression(model, (-omega * L_cable[i] / (((R_cable[i])^2 + (omega * L_cable[i])^2)) + omega*C_cable[i]/2))
    end

    for i in 1:num_nodes

        # diagonal terms
        G[i, i] = @NLexpression(model, sum( cable_conductance[cable_cons[i]][j] for j in 1:length(cable_cons[i])))
        B[i, i] = @NLexpression(model, sum( cable_susceptance_0[cable_cons[i]][j] for j in 1:length(cable_cons[i])))

        # off diagonal terms
        for k in (i+1):num_nodes # this is over the upper triangle

            if CM[i, k] != 0

                cable_num = abs(CM[i, k])

                G[i, k] = @NLexpression(model, -1*cable_conductance[cable_num])
                B[i, k] = @NLexpression(model, -1*cable_susceptance_1[cable_num])
                G[k, i] = G[i, k]
                B[k, i] = B[i, k]
                
            else

                G[i, k] = zero_expression # a formula which returns 0.0
                B[i, k] = zero_expression
                G[k, i] = zero_expression
                B[k, i] = zero_expression
            end
        end
    end

    # power flow constraints - this is perfect!!
    for i in 1:num_nodes

        P_node[i] = @NLconstraint(model,

        nodes[i, "P"] == nodes[i,"v"] * sum( nodes[j,"v"] * ((G[i, j] * cos(nodes[i,"theta"] - nodes[j,"theta"]) + B[i, j] * sin(nodes[i,"theta"] - nodes[j,"theta"]))) for j in node_cons[i])
        
        )

        Q_node[i] = @NLconstraint(model,

        nodes[i, "Q"] == nodes[i,"v"] * sum( nodes[j,"v"] * ((G[i, j] * sin(nodes[i,"theta"] - nodes[j,"theta"]) - B[i, j] * cos(nodes[i,"theta"] - nodes[j,"theta"]))) for j in node_cons[i])

        )

    end

    cable_constraints = Array{NonlinearConstraintRef, 1}(undef, num_cables)
    # maybe remove this? but add as check after optimisation has been completed.
    #TODO: this back in? -> Septimus
    #=
    for i in 1:num_cables

        j, k = Tuple(findfirst(x -> x == i, CM))

        cable_constraints[i] = @NLconstraint(model,
            abs( nodes[j, "v"] * nodes[k, "v"] * (sin(nodes[j, "theta"] - nodes[k, "theta"]))/(omega*cables[i, "L"])) # this formula is not quite correct - missing resistances and capacitances
            <= 0.93 * nodes[j, "v"] * nodes[k, "v"] * sqrt(cables[i, "C_L"]) # check if there should be a 2 in the equation
        )

    end
    =#
    

    # non-linear objectives 
    @NLexpression(model, P_source_mean, sum(nodes[Int(j),"P"] for j in idx_p_mean_cal) / convert.(Int64,length(idx_p_mean_cal)))
    @NLexpression(model, Q_source_mean, sum(nodes[Int(j),"Q"] for j in idx_q_mean_cal) / convert.(Int64,length(idx_q_mean_cal)))

    @NLexpression(model, v_mean, sum(nodes[j,"v"] for j in 1:num_nodes) / num_nodes)

    # add apparent power
    @NLexpression(model, Power_apparent, sum(sqrt(nodes[i, "P"]^2 + nodes[i, "Q"]^2) for i in 1:num_source))
    # TODO: normalisation, i.e. weighting between minimising P and Q, and minimising cable radius? Maybe use per unit system
    # normalisation : 1. max value
    #                 2. p.u. 0
    # Sbase_1_phase = sum(loads)
    # Vbase_rms = 230
    # Ibase_rms = f(Sbase_1_phase, Vbase_rms)
    # Zbase = f(Vbase_rms, Ibase_rms)

    radius_upper_bound = upper_bound(cables[1, "radius"])
    radius_lower_bound = lower_bound(cables[1, "radius"]);

    # Lagrangians
    λ₁ = 0.01
    λ₂ = 0.99

    norm_P = length(idx_p_mean_cal)
    norm_Q = length(idx_q_mean_cal)
    @NLobjective(model, Min, λ₁ *  Power_apparent / total_S_source
                            + abs((v_mean - parameters["grid"]["v_rms"]) / parameters["grid"]["v_rms"])
                            #+ abs(sum(nodes[i,"theta"] for i in 2:num_nodes))/π    # maybe helpfull?
                            + sum( ((nodes[i,"P"] - P_source_mean)^2 )/ norm_P for i in idx_p_mean_cal) 
                            + sum( ((nodes[i,"Q"] - Q_source_mean)^2) / norm_Q for i in idx_q_mean_cal) # the variance - not exactly right (but good enough)
                            + λ₂ * sum( (cables[i, "radius"]  for i in 1:num_cables)) / (num_cables * (radius_upper_bound- radius_lower_bound) )
                            )

    optimize!(model)
    println("""
    termination_status = $(termination_status(model))
    primal_status      = $(primal_status(model))
    objective_value    = $(objective_value(model))
    """)



    println()
    println("Radius : $(value.(cables))")
    println("")
    println.("R :  $(value.(R_cable))")
    println()
    println.("L :  $(value.(L_cable))")
    println()
    println.("C :  $(value.(C_cable))")
    

    println()
    println()
    println(value.(nodes))

    #R = (omega*L)/X_R 
    #C = C_L*L
    
    for (index, cable) in enumerate(parameters["cable"])

        cable["L"] = value.(L_cable)[index]
        cable["Lb"] = cable["L"]/cable["len"]

        cable["R"] = value.(R_cable)[index]
        cable["Rb"] = cable["R"]/cable["len"]

        cable["C"] = value.(C_cable)[index]
        cable["Cb"] = cable["C"]/cable["len"]
    end
    

    return parameters
end


CM = [  0   0   0   1
        0   0   0   2
        0   0   0   3
        -1  -2  -3  0  ] # user specified



num_source = 3 # user
num_load = 1 # user

parameters = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("fltr"=>"LC", "pwr"=>1000, "v_pu_set" => 1.0, "mode" => 1, "control_type" => "classic"),
                    Dict{Any, Any}("fltr"=>"LC", "pwr"=>1000, "v_pu_set" => 1.0, "mode" => 1, "control_type" => "classic",
                                    "p_set"=>500, "q_set"=>500),
                    Dict{Any, Any}("fltr"=>"LC", "pwr"=>1000, "v_pu_set" => 1.0, "mode" => 1, "control_type" => "classic"),
                    #Dict{Any, Any}("fltr"=>"LC", "pwr"=>1000, "v_pu_set" => 1.0, "mode" => "PQ Control”"),
                    #Dict{Any, Any}("fltr"=>"LC", "pwr"=>1000, "v_pu_set" => 1.0, "mode" => "Voltage Control”"),
                    ],
    "load"   => Any[
                    #Dict{Any, Any}("R"=>10, "L" => 0.16, "impedance"=>"RL")#, "pwr"=>1000)
                    Dict{Any, Any}("R"=>58.7, "L"=>0.3863, "impedance"=>"RL")
                    ],
    "grid"   => Dict{Any, Any}("fs"=>50.0, "phase"=>1, "v_rms"=>230, "fg" => 50),
    "cable" => Any[
                    Dict{Any, Any}("len"=>1),
                    Dict{Any, Any}("len"=>1),
                    Dict{Any, Any}("len"=>1)
                    ]
)
#TODO: shift this for every load to the nodeconstructor
for (index, load) in enumerate(parameters["load"])
    # example for RL load
    print(load)
    if !haskey(load, "pwr")
        # parallel R||L
        Z = 1im*parameters["grid"]["fg"]*2*pi*load["R"]*load["L"]/(load["R"]+1im*parameters["grid"]["fg"]*2*pi*load["L"])
        load["pwr"] = parameters["grid"]["v_rms"]^2 / abs(Z) * parameters["grid"]["phase"]
    end
end

#TODO ensure that len is defined in nodeconstructor before this function is called!
#TODO Take len in PFE into acount! (currently assuming 1km)
parameters = layout_cabels(CM, num_source, num_load, parameters)


