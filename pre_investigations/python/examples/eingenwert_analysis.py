import numpy as np
from pre_investigations.python.dare.utils.nodeconstructor import NodeConstructor
import control
import matplotlib.pyplot as plt

"""
File to excite S2_L1_2C grid with specific real part of the eigenvalue (EV) as initial state to investigate the frequency-
response and compare it to the imaginary part of the eigenvalue (Im(EV) = omega!) 
"""

# define parameters
R = 0.4
L = 2.3e-3
C = 10e-6
LT = 2.3e-3
R_load = 14

parameter = dict()
parameter['R_source'] = R
parameter['L_source'] = L
parameter['C_source'] = C
parameter['L_cable'] = LT
parameter['R_cable'] = R
parameter['R_load'] = R_load

# set up CM
#An entry in the CM reads as follows: Object x (row) is connected to object y (column) via connection line z:
# x ---z---> y
CM = np.array([[0, 0, 1],
               [0, 0, 2],
               [-1, -2, 0]])

# grid with 2 sources, 1 load and a total of 2 connections
Grid_S2_L1_2C = NodeConstructor(2, 1, parameter, CM = CM)
Grid_S2_L1_2C.draw_graph()


A, B, C, D = Grid_S2_L1_2C.get_sys()

sys = control.ss(A, B, C, D)

# define time vector
ts = 1e-6
t_end = 0.005
steps = int(1/ts)
t = np.arange(0, t_end+ts, ts)
num_samples = len(t)

# generate init state
#x0 = np.zeros((A.shape[0],1))
x01 = np.ones((A.shape[0],1)) * 100
#x0 = np.array([-0.007688, -0.70557, -0.007688, -0.70557, -0.019226, -0.019226]) #EW6
x0 = np.array([-0.033502, -0.69945, -0.033502,-0.69945, -0.09821, -0.09821])    #EW1 (REAL)
x0 = np.array([0.5, 0, -0.5, 0, 0.5, -0.5])    #EW2 (REAL)

# simple input signal of constant 230V from all sources
u = np.array([230]).repeat(Grid_S2_L1_2C.num_source)[:,None] * np.zeros((Grid_S2_L1_2C.num_source,len(t)))
#u = np.array([230]).repeat(Grid_S2_L1_2C.num_source)[:,None] * np.ones((Grid_S2_L1_2C.num_source,len(t)))

T, yout, xout = control.forced_response(sys, T=t, U=u, X0=x0, return_x=True, squeeze=True)

plt.plot(t, xout[1], label='$v_1$')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
plt.ylabel('$v_{\mathrm{1}}\,/\,\mathrm{V}$')
plt.title('Plot current $v_1$')
plt.legend()
plt.grid()
plt.show()

plt.plot(t, xout[0], label='$i_1$')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
plt.ylabel('$i_{\mathrm{1}}\,/\,\mathrm{A}$')
plt.legend()
plt.grid()
plt.show()

plt.plot(t, xout[3], label='$v_2$')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
plt.ylabel('$v_{\mathrm{2}}\,/\,\mathrm{V}$')
plt.legend()
plt.grid()
plt.show()

plt.plot(t, xout[2], label='$i_2$')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
plt.ylabel('$i_{\mathrm{2}}\,/\,\mathrm{A}$')
plt.legend()
plt.grid()
plt.show()


plt.plot(t, xout[4], label='$i_T1$')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
plt.ylabel('$i_{\mathrm{T1}}\,/\,\mathrm{A}$')
plt.legend()
plt.grid()
plt.show()
plt.plot(t, xout[5], label='$i_T2$')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
plt.ylabel('$i_{\mathrm{T2}}\,/\,\mathrm{A}$')
plt.legend()
plt.grid()
plt.show()

