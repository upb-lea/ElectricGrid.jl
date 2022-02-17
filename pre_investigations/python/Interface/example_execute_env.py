import numpy as np
import matplotlib.pyplot as plt

from pre_investigations.python.Interface.env_dare import Env_DARE

"""
Runs the "2 Nodes with 2 Connections" example using a gym env for comparison to the jupyter notebook solutions.
Actions are constant v_i to both sources. 
"""

# define time vector
ts = 1e-4
t_end = 0.01
steps = int(1 / ts)
t = np.arange(0, t_end + ts, ts)
num_samples = len(t)

# network params
R = 0.4
L = 2.3e-3
C = 10e-6
LT = 2.3e-3
R_load = 14

# x0 = np.array([i10, v10, iT10, i20, v20, iT20])
x0 = [0, 0, 0, 0, 0, 0]

vi1 = 230
vi2 = 230

parameter = dict()
parameter['R_source'] = R
parameter['L_source'] = L
parameter['C_source'] = C
parameter['L_cabel'] = LT
parameter['R_cabel'] = R
parameter['R_load'] = R_load
parameter['V_dc'] = 300

# limits for normalization
i_lim = 20
v_lim = 600

limits = dict()
limits['i_lim'] = i_lim
limits['v_lim'] = v_lim

# define CM
CM = np.array([[0, 0, 1],
               [0, 0, 2],
               [-1, -2, 0]])

# simple ref - TEMPORARY IMPLEMENTATION
ref = 200  # W

env = Env_DARE(CM, ts, parameter, x0, limits, ref)

v1_list = []

obs = env.reset()
v1_list.append(obs[1])
for i in range(num_samples - 1):
    action = np.array([[vi1], [vi2]])
    obs, rewards, dones, info = env.step(action)
    # store v1 for comparison
    v1_list.append(obs[1])

plt.figure()
plt.step(t, v1_list, 'r', label='v1_discret')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
plt.ylabel('$v_{\mathrm{1}}\,/\,\mathrm{V}$')
plt.legend()
plt.grid()
plt.ylim([0, 340])
plt.show()
