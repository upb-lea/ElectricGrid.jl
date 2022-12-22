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
        -1  -2  -3  0  ] # user specified

#TODO: user specified disctances
#distances = [1.0 2.3 .323]

num_source = 3 # user
num_load = 1 # user


num_nodes = num_source + num_load
num_cables = maximum(CM)

@variable(model, nodes[1 : num_nodes, ["v", "theta", "P", "Q"]])
fix(nodes[1, "theta"], 0.0) # reference
fix(nodes[1, "v"], 230.0) # reference should not be a load bus - ensure, TODO

function set_bounds(variable, start_value, low_bound, up_bound)
    if !is_fixed(variable)
        set_lower_bound(variable, low_bound)
        set_upper_bound(variable, up_bound)
        set_start_value(variable, start_value)
    end
end

for i = 1:num_nodes
    if i <= num_source
        fix(nodes[i, "v"], 230.0) # user may specify 1.05*230

        set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # question, does this limit too much, should be more or less.
        set_bounds(nodes[i, "P"], (num_load * 1000.0) / num_source, -1000.0, 1000.0) # come from parameter dict/user?
        set_bounds(nodes[i, "Q"], (num_load * -100.0) / num_source, -1000.0, 1000.0) # P and Q are the average from power, excluding cable losses
    else
        fix(nodes[i, "P"], -1000.0)
        fix(nodes[i, "Q"], -1000.0)

        set_bounds(nodes[i, "theta"], 0.0, -0.25*pi/2, 0.25*pi/2) # same as above
        set_bounds(nodes[i, "v"], 230.0, 0.95*230, 1.05*230)
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

cable_cons = get_cable_connections()
node_cons = get_node_connections() 

G = Array{NonlinearExpression, 2}(undef, num_nodes, num_nodes) # should be symmetric
B = Array{NonlinearExpression, 2}(undef, num_nodes, num_nodes) # should be symmetric

P_node = Array{NonlinearConstraintRef, 1}(undef, num_nodes)
Q_node = Array{NonlinearConstraintRef, 1}(undef, num_nodes)

# As radius goes down resistance goes up, inductance goes up, capacitance goes down. Put in formulas for this.
#TODO ["radius"] for futute : ["radius", "conductivity"] 
@variable(model, cables[1 : num_cables, ["radius"]]) # this is wrong - 

# defining indivdual expressions for R, L, C as a function of radius with fixed length 
#=
    
    # L = 0.00025
    # assumption for single phase two wire system
    # space b/w conductors
    
    # r' = 1e-(1/4)r

    L =  l * 4e-7 * log(D/0.7788 * r)  # m* H/m

    # resistivity remains constant ρ_(T=50) = 1.973e-8 
    R = l * 1.973e-8 / (π * r^2) # m * Ω/m

    # X_R = 0.38 # ratio of omega*L/R
    # X_R = omega * L(r)/R(r)

    # line to neutral
    C = l * 2π * 8.854e-12 / (log(D/r)) #l * F/m 
=#

# fixed length
l = 1000 #m
# fixed distance b/w cables
D = 0.5 #m

L_cable = Array{NonlinearExpression, 1}(undef, num_cables)
R_cable = Array{NonlinearExpression, 1}(undef, num_cables)
C_cable = Array{NonlinearExpression, 1}(undef, num_cables)
cable_conductance = Array{NonlinearExpression, 1}(undef, num_cables)
cable_susceptance_0 = Array{NonlinearExpression, 1}(undef, num_cables) # diagonals - where we add capacitances
cable_susceptance_1 = Array{NonlinearExpression, 1}(undef, num_cables) # off diagonals


for i=1:num_cables

    # set_bounds(cables[i, "L"], 0.00025, 0.00023, 0.00026)
    # set_bounds(cables[i, "X_R"], 0.38, 0.37, 0.4)
    # set_bounds(cables[i, "C_L"], 0.0016, 0.0015, 0.0017)

    #TODO: AWG 6 - AWG 12 as discussed
    set_bounds(cables[i, "radius"], (2.05232e-3)/2, (2.05232e-3)/2, (4.1148e-3)/2) #m 

    # assumption to line to line
    L_cable[i] = @NLexpression(model, l * 4e-7 * log(D/(0.7788 * cables[i, "radius"])))  # m* H/m

    # resistivity remains constant ρ_(T=50) = 1.973e-8 
    R_cable[i] = @NLexpression(model, l * 1.973e-8 / (π * cables[i, "radius"]^2)) # m * Ω/m

    # X_R = 0.38 # ratio of omega*L/R
    # X_R = omega * L(r)/R(r)

    # line to neutral
    C_cable[i] = @NLexpression(model, l * 2π * 8.854e-12 / (log(D/cables[i, "radius"]))) #m * F/m

    # B1 = -omega * L /(((omega*L)/X_R)^2 + omega^2 * L^2)
    # B1 = -omega * L(r) / (R(r)^2 + (omega * L(r))^2)
    # # G1 = ((omega*L)/X_R) /(((omega*L)/X_R)^2 + omega^2 * L^2)
    # G1 = R(r) / (R(r)^2 + (omega*L(r))^2)
    # B0 = (-omega * L /(((omega*L)/X_R)^2 + omega^2 * L^2)) + omega*C_L*L
    # B0 = (-omega*L(r) / (R(r)^2 + (omega*L(r))^2 )) + omega*C(r)  

    # cable_conductance[i] = @NLexpression(model, ((omega * cables[i, "L"]) / cables[i, "X_R"]) / (((omega*cables[i, "L"]) / cables[i, "X_R"])^2 + omega^2 * cables[i, "L"]^2))
    cable_conductance[i] = @NLexpression(model, (R_cable[i] / ((R_cable[i])^2 + (omega * L_cable[i])^2)))
    # cable_susceptance_1[i] = @NLexpression(model, (-omega * cables[i, "L"] / (((omega*cables[i, "L"])/cables[i, "X_R"])^2 + omega^2 * cables[i, "L"]^2)))
    cable_susceptance_1[i] = @NLexpression(model, (-omega * L_cable[i] / (((R_cable[i])^2 + (omega * L_cable[i])^2))))
    # cable_susceptance_0[i] = @NLexpression(model, (-omega * cables[i, "L"] / (((omega*cables[i, "L"])/cables[i, "X_R"])^2 + omega^2 * cables[i, "L"]^2)) + omega*cables[i, "C_L"]*cables[i, "L"]/2)
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
# for i in 1:num_cables

#     j, k = Tuple(findfirst(x -> x == i, CM))

#     cable_constraints[i] = @NLconstraint(model,
#     # fins suitable expression
#         abs( nodes[j, "v"] * nodes[k, "v"] * (sin(nodes[j, "theta"] - nodes[k, "theta"]))/(omega*L_cable[i])) # this formula is not quite correct - missing resistances and capacitances
#         <= 0.93 * nodes[j, "v"] * nodes[k, "v"] * sqrt(cables[i, "C_L"]) # check if there should be a 2 in the equation
#     )

# end
#0.93 * value(nodes[j, "v"] * nodes[k, "v"]) * sqrt(value(cables[i, "C_L"]))
# non-linear objectives
@NLexpression(model, P_source_mean, sum(nodes[j,"P"] for j in 1:num_source) / num_source)
@NLexpression(model, Q_source_mean, sum(nodes[j,"Q"] for j in 1:num_source) / num_source)
# add apparent power
@NLexpression(model, Power_apparent, sum(sqrt(nodes[i, "P"]^2 + nodes[i, "Q"]^2) for i in 1:num_source))
# TODO: normalisation, i.e. weighting between minimising P and Q, and minimising cable radius? Maybe use per unit system
# normalisation : 1. max value
#                 2. p.u. 

# Sbase_1_phase = sum(loads)
# Vbase_rms = 230
# Ibase_rms = f(Sbase_1_phase, Vbase_rms)
# Zbase = f(Vbase_rms, Ibase_rms)
radius_upper_bound = upper_bound(cables[1, "radius"]);

# Lagrangians
λ₁ = 0.99
λ₂ = 0.01

@NLobjective(model, Min, λ₁ *  Power_apparent / (num_source * sqrt(2) * 1000)
                        # + abs(sum(nodes[i,"P"]/1000 for i in 1:num_source)) # replace sum by mean or divide by 1000*num_source
                        # + abs(sum(nodes[i,"Q"]/1000 for i in 1:num_source)) # take apparent power -- λ_1  
                        #  + sum(nodes[i,"v"]/230 for i in num_source+1:num_nodes) # TODO: maybe wrong- minimise the deviation of voltage, excluding reference node (which does not neeed to be node 1)
                        # + abs(sum(nodes[i,"theta"]/pi for i in 2:num_nodes)) # variances (?)
                        # + sum(cables[i, "X_R"] for i in 1:num_cables) # replaced with radius, how do we normalise radius??? idea: upper bound of radius times by number of cables
                        # + sum(1/cables[i, "L"] for i in 1:num_cables)  # -- λ_2 
                        # + sum(cables[i, "C_L"] for i in 1:num_cables) 
                        + sum( ((nodes[i,"P"] - P_source_mean)^2)/num_source for i in 1:num_source) 
                        + sum( ((nodes[i,"Q"] - Q_source_mean)^2)/num_source for i in 1:num_source) 
                        + λ₂ * sum( (cables[i, "radius"]  for i in 1:num_cables)) / (num_cables * radius_upper_bound )# num_cables
                        ) # the variance - not exactly right (but good enough)

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

println("")
println.("R :  $(value.(R_cable))")


println()
println.("L :  $(value.(L_cable))")

println()
println.("C :  $(value.(C_cable))")