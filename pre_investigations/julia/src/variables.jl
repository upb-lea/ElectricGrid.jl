using LinearAlgebra

ts = 1e-4
t_end = 0.1
t = collect(0:ts:t_end+ts)
t = round.(t,digits=4)
#num_samples = len(t)
num_samples = size(t)[1]


R = 0.4
L = 2.3e-3
C = 10e-6
LT = 2.3e-3
RLoad = 14

vi1 = 230
vi2 = 230


R1 = R
L1 = L
C1 = C
R2 = R
L2 = L
C2 = C
LT1 = LT
LT2 = LT
RT1 = R
RT2 = R

i10 = 0
v10 = 0
iT10 = 0
i20 = 0
v20 = 0
iT20 = 0
t0 = 0

# x0 = np.array([i10, v10, iT10, i20, v20, iT20])
x0 = [i10, v10, iT10, i20, v20, iT20]





A1 = [-R1/L1 -1/L1 0
    1/C1 0 -1/C1
    0 1/LT1 -(RLoad/LT1)-(RT1/LT1)]

AT1 = [0 0 0
       0 0 0
       0 0 -(RLoad/LT1)]

A2 = [-R2/L2 -1/L2 0
    1/C2 0 -1/C2
    0 1/LT2 -(RLoad/LT2)-(RT2/LT2)]

AT2 = [0 0 0
       0 0 0
       0 0 -(RLoad/LT2)]


A = [A1 AT1
     AT2 A2]

B0 = zeros(3,1)

B1 = [1/L1
       0
       0]

B2 = [1/L2
       0
       0]

B = [B1 B0
     B0 B2]

C = diagm(ones(6))


Ad = exp(A*ts)

Bd = A \ (Ad - C) * B


f0 = 50
V_eff = 230 * sqrt(2)

v_sin1 = V_eff * sin.(2*pi * f0 * t)
v_sin2 = V_eff * sin.(2*pi * f0 * t .+ 1)

u = [v_sin1, v_sin2]
u = mapreduce(permutedims, vcat, u);