using JuMP
import Ipopt


# Constant values
f = 50
omega = 2*π*f

# Reference Bus 1
v1 = 230
theta1 = 0
PL_1 = 0
QL_1 = 0
Pmax_1 = 50000
Qmax_1 = 50000

# Load Bus 2
PL_2 = 100
PG_2 = 0
QL_2 = 500
QG_2 = 0

# PV Bus 3 
v3 = 235
PL_3 = 100
PG_3 = 100
QL_3 = 100


# Concuctance and Susceptance matrices -> get from nc.get_Ybus() -> real()/imag()
L = 0.0025
R = 0.208
B1 = -omega * L /(R^2 + omega^2 * L^2)
G1 = R /(R^2 + omega^2 * L^2)

G = [[2*G1, -G1, -G1],
     [-G1, 2*G1, -G1],
     [-G1, -G1, 2*G1]]

b_variable = 1e9
b_variable2 = 1e9

#B = [[2*B1, -B1, -B1],
#B = [[B1+b_variable, b_variable, -B1],
 B = [[b_variable2+b_variable, -b_variable, -b_variable2],
      [-b_variable,          2*B1, -B1],
      [-b_variable2,                  -B1, 2*B1]]

model = Model(Ipopt.Optimizer)

# Variables including var_constraints
@variable(model, 99 <= PG_1 <= 101)
@variable(model, 104 <= QG_1 <= 107)
@variable(model, 227 <= v2 <= 231)
@variable(model, -0.1 <= theta2 <= 0.1)
@variable(model, -0.1 <= theta3 <= 0.1)
@variable(model, 490 <= QG_3 <= 500)
@variable(model, -2.1 <= b_var <= 100.0)
@variable(model, -2.1 <= b_var2 <= 100.0)

# non-linear objectives
@NLobjective(model, Min, abs(PG_1)/200 + abs(QG_1)/500 + abs(v2)/240 + abs(theta2)/0.1 + abs(theta3)/0.1 + abs(QG_3)/700 + abs(b_var)/1.5 + abs(b_var2)/1.5)


@NLconstraint(model, P_Bus1,
    v1 * v1 * (G[1][1] * cos(theta1 - theta1) + (b_var2 + b_var) * sin(theta1 - theta1)) + 
    v1 * v2 * (G[1][2] * cos(theta1 - theta2) - b_var * sin(theta1 - theta2)) +
    v1 * v3 * (G[1][3] * cos(theta1 - theta3) + b_var2 * sin(theta1 - theta3)) 
    == PG_1 - PL_1)

@NLconstraint(model, P_Bus2,
    v2 * v1 * (G[2][1] * cos(theta2 - theta1) - b_var * sin(theta2 - theta1)) + 
    v2 * v2 * (G[2][2] * cos(theta2 - theta2) + B[2][2] * sin(theta2 - theta2)) +
    v2 * v3 * (G[2][3] * cos(theta2 - theta3) + B[2][3] * sin(theta2 - theta3))
    == PG_2 - PL_2)

@NLconstraint(model, P_Bus3,
    v3 * v1 * (G[3][1] * cos(theta3 - theta1) - b_var2 * sin(theta3 - theta1)) + 
    v3 * v2 * (G[3][2] * cos(theta3 - theta2) + B[3][2] * sin(theta3 - theta2)) +
    v3 * v3 * (G[3][3] * cos(theta3 - theta3) + B[3][3] * sin(theta3 - theta3)) 
    == PG_3 - PL_3)

@NLconstraint(model, Q_Bus1,
    v1 * v1 * (G[1][1] * sin(theta1 - theta1) - (b_var2 + b_var) * cos(theta1 - theta1)) + 
    v1 * v2 * (G[1][2] * sin(theta1 - theta2) + b_var * cos(theta1 - theta2)) +
    v1 * v3 * (G[1][3] * sin(theta1 - theta3) - B[1][3] * cos(theta1 - theta3))
    == QG_1 - QL_1)

@NLconstraint(model, Q_Bus2,
    v2 * v1 * (G[2][1] * sin(theta2 - theta1) + b_var * cos(theta2 - theta1)) + 
    v2 * v2 * (G[2][2] * sin(theta2 - theta2) - B[2][2] * cos(theta2 - theta2)) +
    v2 * v3 * (G[2][3] * sin(theta2 - theta3) - B[2][3] * cos(theta2 - theta3))
    == QG_2 - QL_2)

@NLconstraint(model, Q_Bus3,
    v3 * v1 * (G[3][1] * sin(theta3 - theta1) + b_var2 * cos(theta3 - theta1)) + 
    v3 * v2 * (G[3][2] * sin(theta3 - theta2) - B[3][2] * cos(theta3 - theta2)) + 
    v3 * v3 * (G[3][3] * sin(theta3 - theta3) - B[3][3] * cos(theta3 - theta3))
    == QG_3 - QL_3)

# Linear constraints:
#@constraint(model, lc1, PG_1 + QG_1 <= Smax_1)
#@constraints(model, Bmatrix,
#    B .== [[B1+b_var, b_var, -B1],
#         [b_var, 2*B1, -B1],
#         [-B1, -B1, 2*B1]])

optimize!(model)
println("""
termination_status = $(termination_status(model))
primal_status      = $(primal_status(model))
objective_value    = $(objective_value(model))
""")
println("PG1, QG1, v2, theta2, theta3, QG_3, b_var, b_var2")
println("$(value(PG_1)), $(value(QG_1)), $(value(v2)), $(value(theta2)), $(value(theta3)), $(value(QG_3)), $(value(b_var)), $(value(b_var2))")

