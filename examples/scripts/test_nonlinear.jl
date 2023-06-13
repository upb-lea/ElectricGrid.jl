using DifferentialEquations
using PlotlyJS

function NonlinearInductance(x,value, drop, length)
    return value*(1-drop/100) + value*(drop/100)*1/(1+exp(-(abs(x)-length/2)/2))
end

p = 80
l(x) = NonlinearInductance(x,0.0001,p,600)

R_S = 1.1e-3
R_K = 1e-3
R_L = 100
L_S = 70e-6
L_K(x) = l(x)
L_L = 1e-2
C_K = 1e-4/4
C_SK = C_K
C_L = 1e-2
C_KL = C_K+C_L
U_i = 250

A(x) = [-R_S/L_S -1/L_S 0 0 0
    1/C_SK 0 -1/C_SK 0 0
    0 1/L_K(x[3]) -R_K/L_K(x[3]) -1/L_K(x[3]) 0
    0 0 1/C_KL -1/(R_L*C_KL) -1/C_KL
    0 0 0 1/L_L 0]
    
B = [1/L_S;0;0;0;0]
u = U_i

function f(dx, x, p, t)
    dx.= A(x)*x+B*u
end

x0 = [0,0,0,0,0]
tspan = (0.0,0.1)
prob = ODEProblem(f, x0, tspan)
alg = Tsit5()

sol = solve(prob, alg,reltol=1e-8, abstol=1e-8);

x = sol.u
x_plot(j) = [x[i][j] for i = 1:length(x)]
t = sol.t
a = ["i_LS","u_CK","i_LK","u_CKL","i_LL"]
trace(i) = scatter(x=t, y=x_plot(i), mode="lines",name = "$(a[i])")

layout = Layout(title="Percentage of nonlinearity $(p)%",xaxis_title="time", yaxis_title="Amplitude")

# plot([trace(1),trace(2),trace(3),trace(4),trace(5)], layout)
plot([trace(2),trace(4)], layout)