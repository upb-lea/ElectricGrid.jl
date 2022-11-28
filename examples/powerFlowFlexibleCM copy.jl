using JuMP
import Ipopt

model = Model(Ipopt.Optimizer)
#set_optimizer_attributes(model, "tol" => 1e-1)

zero_expression = @NLexpression(model, 0.0)

# Constant values
f = 50
omega = 2π*f

# for every Source: v is fixed 230
# for one Source: theta is fixed 0

CM = [  0   0   0   1
        0   0   0   2
        0   0   0   3
        -1  -2  -3  0  ]

num_source = 3
num_load = 1
num_nodes = num_source + num_load
num_cables = maximum(CM)

@variable(model, nodes[1 : num_source + num_load, ["v", "theta", "P", "Q"]])
fix(nodes[1, "theta"], 0.0)
fix(nodes[1, "v"], 230.0)

function set_bounds(variable, start_value, low_bound, up_bound)
    if !is_fixed(variable)
        set_lower_bound(variable, low_bound)
        set_upper_bound(variable, up_bound)
        set_start_value(variable, start_value)
    end
end

for i = 1:num_nodes
    if i <= num_source
        fix(nodes[i, "v"], 230.0)

        set_bounds(nodes[i, "theta"], 0.0, -pi/2, pi/2)
        set_bounds(nodes[i, "P"], (num_load * 1000.0) / num_source, -1000.0, 1000.0)
        set_bounds(nodes[i, "Q"], (num_load * -100.0) / num_source, -1000.0, 1000.0)
    else
        fix(nodes[i, "P"], -1000.0)
        fix(nodes[i, "Q"], -1000.0)

        set_bounds(nodes[i, "theta"], 0.0, -pi/2, pi/2)
        set_bounds(nodes[i, "v"], 230.0, 0.95*230, 1.05*230)
    end
end

function get_degree(CM = CM)
    result = zeros(Int, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = count(x -> x != 0, CM[i,:])
    end

    result
end

function get_cable_connections(CM = CM)

    result = Vector{Vector{Int64}}(undef, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = filter(x -> x != 0, abs.(CM[i,:]))
    end

    return result
end

function get_node_connections(CM = CM)

    result = Vector{Vector{Int64}}(undef, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = findall(x -> x != 0, CM[i,:])
        push!(result[i], i)
    end

    return result
end

cable_cons = get_cable_connections()
node_cons = get_node_connections() 

G = Array{NonlinearExpression, 2}(undef, num_nodes, num_nodes) # should be symmetric
B = Array{NonlinearExpression, 2}(undef, num_nodes, num_nodes) # should be symmetric

P_node = Array{NonlinearConstraintRef, 1}(undef, num_nodes)
Q_node = Array{NonlinearConstraintRef, 1}(undef, num_nodes)

@variable(model, cables[1 : num_cables, ["L", "X_R", "C_L"]])
cable_conductance = Array{NonlinearExpression, 1}(undef, num_cables)
cable_susceptance_0 = Array{NonlinearExpression, 1}(undef, num_cables)
cable_susceptance_1 = Array{NonlinearExpression, 1}(undef, num_cables)

for i=1:num_cables

    set_bounds(cables[i, "L"], 0.00025, 0.00023, 0.00026)
    set_bounds(cables[i, "X_R"], 0.38, 0.37, 0.4)
    set_bounds(cables[i, "C_L"], 0.0016, 0.0015, 0.0017)

    #R = (omega*L)/X_R
    #C = C_L*L
    cable_conductance[i] = @NLexpression(model, ((omega * cables[i, "L"]) / cables[i, "X_R"]) / (((omega*cables[i, "L"]) / cables[i, "X_R"])^2 + omega^2 * cables[i, "L"]^2))
    cable_susceptance_1[i] = @NLexpression(model, (-omega * cables[i, "L"] / (((omega*cables[i, "L"])/cables[i, "X_R"])^2 + omega^2 * cables[i, "L"]^2)))
    cable_susceptance_0[i] = @NLexpression(model, (-omega * cables[i, "L"] / (((omega*cables[i, "L"])/cables[i, "X_R"])^2 + omega^2 * cables[i, "L"]^2)) + omega*cables[i, "C_L"]*cables[i, "L"]/2)
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

            G[i, k] = zero_expression
            B[i, k] = zero_expression
            G[k, i] = zero_expression
            B[k, i] = zero_expression
        end
    end
end

for i in 1:num_nodes

    P_node[i] = @NLconstraint(model,

    nodes[i, "P"] == nodes[i,"v"] * sum( nodes[j,"v"] * ((G[i, j] * cos(nodes[i,"theta"] - nodes[j,"theta"]) + B[i, j] * sin(nodes[i,"theta"] - nodes[j,"theta"]))) for j in node_cons[i])
    
    )

    Q_node[i] = @NLconstraint(model,

    nodes[i, "Q"] == nodes[i,"v"] * sum( nodes[j,"v"] * ((G[i, j] * sin(nodes[i,"theta"] - nodes[j,"theta"]) - B[i, j] * cos(nodes[i,"theta"] - nodes[j,"theta"]))) for j in node_cons[i])

    )

end

cable_constraints = Array{NonlinearConstraintRef, 1}(undef, num_cables)

for i in 1:num_cables

    j, k = Tuple(findfirst(x -> x == i, CM))

    cable_constraints[i] = @NLconstraint(model,
        abs( nodes[j, "v"] * nodes[k, "v"] * (sin(nodes[j, "theta"] - nodes[k, "theta"]))/(omega*cables[i, "L"]))
        <= 0.93 * nodes[j, "v"] * nodes[k, "v"] * sqrt(cables[i, "C_L"]) # check if there should be a 2 in the equation
    )

end
#0.93 * value(nodes[j, "v"] * nodes[k, "v"]) * sqrt(value(cables[i, "C_L"]))
# non-linear objectives
@NLexpression(model, P_source_mean, sum(nodes[j,"P"] for j in 1:num_source) / num_source)
@NLexpression(model, Q_source_mean, sum(nodes[j,"Q"] for j in 1:num_source) / num_source)

@NLobjective(model, Min, abs(sum(nodes[i,"P"] for i in 1:num_source))/1000
                        + abs(sum(nodes[i,"Q"] for i in 1:num_source))/1000
                        + sum(nodes[i,"v"] for i in num_source+1:num_nodes)/230
                        + abs(sum(nodes[i,"theta"] for i in 2:num_nodes))/π
                        + sum(cables[i, "X_R"] for i in 1:num_cables)
                        + sum(1/cables[i, "L"] for i in 1:num_cables)
                        + sum(cables[i, "C_L"] for i in 1:num_cables)
                        + sum( (nodes[i,"P"] - P_source_mean)^2 for i in 1:num_source)/num_source 
                        + sum( (nodes[i,"Q"] - Q_source_mean)^2 for i in 1:num_source)/num_source ) # the variance 

optimize!(model)
println("""
termination_status = $(termination_status(model))
primal_status      = $(primal_status(model))
objective_value    = $(objective_value(model))
""")


println()
println()
println(value.(nodes))

println()
println(value.(cables))