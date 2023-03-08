import numpy as np
from dare.utils.nodeconstructorcableloads import NodeConstructorCableLoads

"""
Checks for a given setting the accruing resonance frequencies
"""


CM = np.array([[ 0.,  1.],
            [-1.,  0.]])

parameter = {'source': [{'fltr': 'LCL', 'R': 0.4, 'L1': 2.3e-3, 'L2': 2.3e-3,'C': 10e-6}],
             'cable': [{'R': 0.722, 'L': 0.955e-3, 'C': 8e-09}],
             'load': [{'impedance': 'R', 'R': 14}]}

Grid_FC = NodeConstructorCableLoads(1, 1, CM=CM, parameter=parameter)

Grid_FC.draw_graph()

A,B,C,D = Grid_FC.get_sys()

w, v = np.linalg.eig(A)

sort_w_im = np.sort(np.imag(w)/(2*np.pi))

print(w)
print(" ")
print(" ")
print(np.sort(np.imag(w)/(2*np.pi)))
print(" ")
print(" ")
print(np.count_nonzero(np.real(w)>=0))

print(" ")
print(" ")

