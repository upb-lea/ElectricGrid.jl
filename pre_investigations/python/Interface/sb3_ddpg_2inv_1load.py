import gym
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from stable_baselines3 import DDPG
from stable_baselines3.common.noise import NormalActionNoise, OrnsteinUhlenbeckActionNoise

from pre_investigations.python.Interface.env_dare import Env_DARE

# define time vector
ts = 1e-4
t_end = 0.01
steps = int(1/ts)
t = np.arange(0, t_end+ts, ts)
num_samples = len(t)

# network params
R = 0.4
L = 2.3e-3
C = 10e-6
LT = 2.3e-3
R_load = 14

# limits for normalization
i_lim = 20
v_lim = 600

limits = dict()
limits['i_lim'] = i_lim
limits['v_lim'] = v_lim

# x0 = np.array([i10, v10, i20, v20, iT10, iT20])
x0 = [0, 0, 0, 0, 0, 0]

vi1 = 200
vi2 = 200

parameter = dict()
parameter['R_source'] = R
parameter['L_source'] = L
parameter['C_source'] = C
parameter['L_cabel'] = LT
parameter['R_cabel'] = R
parameter['R_load'] = R_load

# define CM
CM = np.array([[0, 0, 1],
               [0, 0, 2],
               [-1, -2, 0]])



ref = 200 # W

env = Env_DARE(CM, ts, parameter, x0, limits, ref)


"""
# The noise objects for DDPG
n_actions = env.action_space.shape[-1]
action_noise = NormalActionNoise(mean=np.zeros(n_actions), sigma=0.1 * np.ones(n_actions))

model = DDPG("MlpPolicy", env, action_noise=action_noise, verbose=1)
model.learn(total_timesteps=1000, log_interval=10)
#model.save("ddpg_pendulum")
env = model.get_env()

#del model # remove to demonstrate saving and loading

#model = DDPG.load("ddpg_pendulum")
"""

v1_list = []


obs = env.reset()
obs_matrix = np.array(obs)
v1_list.append(obs[1])
for i in range(num_samples-1):
    #action, _states = model.predict(obs)
    action = np.array([[vi1], [vi2]])
    #action = np.array([action.squeeze(),action.squeeze()])#np.array([[vi1], [vi2]])
    obs, rewards, dones, info = env.step(action)
    obs_matrix = np.vstack((obs_matrix, obs))
    v1_list.append(obs[1])
    #env.render()


# x0 = np.array([i10, v10, iT10, i20, v20, iT20])


# denormalize!
obs_matrix_denorm = obs_matrix*env.norm_array

fig, ax1 = plt.subplots(2, 2)

ax1[0, 0].step(t, obs_matrix_denorm[:, 1], 'r', label='v1')
# plt.plot(t,result[:steps,0], label = 'i1')
ax1[0, 0].set_xlabel(r'$t\,/\,\mathrm{s}$')
ax1[0, 0].set_ylabel('$v_{\mathrm{1}}\,/\,\mathrm{V}$')
# plt.title('{}'.format())
ax1[0, 0].legend()
ax1[0, 0].grid()
#ax1[0].set_ylim([0, 340])

ax1[0, 1].step(t, obs_matrix_denorm[:, 0], 'r', label='i1')
ax1[0, 1].step(t, obs_matrix_denorm[:, 4], '--b', label='iT1')
# plt.plot(t,result[:steps,0], label = 'i1')
ax1[0, 1].set_xlabel(r'$t\,/\,\mathrm{s}$')
ax1[0, 1].set_ylabel('$i_{\mathrm{}}\,/\,\mathrm{A}$')
# plt.title('{}'.format())
ax1[0, 1].legend()
ax1[0, 1].grid()
#ax1[0].set_ylim([0, 340])


ax1[1, 0].step(t, obs_matrix_denorm[:, 3], 'r', label='v2')
# plt.plot(t,result[:steps,0], label = 'i1')
ax1[1, 0].set_xlabel(r'$t\,/\,\mathrm{s}$')
ax1[1, 0].set_ylabel('$v_{\mathrm{2}}\,/\,\mathrm{V}$')
# plt.title('{}'.format())
ax1[1, 0].legend()
ax1[1, 0].grid()
#ax1[0].set_ylim([0, 340])

ax1[1, 1].step(t, obs_matrix_denorm[:, 2], 'r', label='i2')
ax1[1, 1].step(t, obs_matrix_denorm[:, 5], '--b', label='iT2')
# plt.plot(t,result[:steps,0], label = 'i1')
ax1[1, 1].set_xlabel(r'$t\,/\,\mathrm{s}$')
ax1[1, 1].set_ylabel('$i_{\mathrm{}}\,/\,\mathrm{A}$')
# plt.title('{}'.format())
ax1[1, 1].legend()
ax1[1, 1].grid()
#ax1[0].set_ylim([0, 340])
plt.show()

P_ref = [ref] * len(t)
P_load = (obs_matrix_denorm[:, 4]+obs_matrix_denorm[:, 5])*parameter["R_load"]
# plot power via Load
plt.step(t, P_load, 'b', label='P_{load}')
plt.step(t, P_ref, '--', color='gray', label='P_{ref}')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
plt.ylabel('$P_{\mathrm{}}\,/\,\mathrm{W}$')
# plt.title('{}'.format())
plt.legend()
plt.grid()
#ax1[0].set_ylim([0, 340])
plt.show()





