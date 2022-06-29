import matplotlib.pyplot as plt
import numpy as np
from stable_baselines3 import DDPG
from stable_baselines3.common.noise import NormalActionNoise

from pre_investigations.python.Interface.env_dare import Env_DARE

"""
Quick basic example to interact with dare_env V_0 with and without RL agent
"""


def plt_obs(obs_matrix_denorm, title=""):
    fig, ax1 = plt.subplots(2, 2, figsize=(7, 5))

    ax1[0, 0].step(t, obs_matrix_denorm[:, 1], 'r', label='v1')
    # plt.plot(t,result[:steps,0], label = 'i1')
    ax1[0, 0].set_xlabel(r'$t\,/\,\mathrm{s}$')
    ax1[0, 0].set_ylabel('$v_{\mathrm{1}}\,/\,\mathrm{V}$')
    # plt.title('{}'.format())
    ax1[0, 0].legend()
    ax1[0, 0].grid()
    # plt.title(title)
    # ax1[0].set_ylim([0, 340])

    ax1[0, 1].step(t, obs_matrix_denorm[:, 0], 'r', label='i1')
    ax1[0, 1].step(t, obs_matrix_denorm[:, 4], '--b', label='iT1')
    # plt.plot(t,result[:steps,0], label = 'i1')
    ax1[0, 1].set_xlabel(r'$t\,/\,\mathrm{s}$')
    ax1[0, 1].set_ylabel('$i_{\mathrm{}}\,/\,\mathrm{A}$')
    # plt.title('{}'.format())
    ax1[0, 1].legend()
    ax1[0, 1].grid()
    # ax1[0].set_ylim([0, 340])

    ax1[1, 0].step(t, obs_matrix_denorm[:, 3], 'r', label='v2')
    # plt.plot(t,result[:steps,0], label = 'i1')
    ax1[1, 0].set_xlabel(r'$t\,/\,\mathrm{s}$')
    ax1[1, 0].set_ylabel('$v_{\mathrm{2}}\,/\,\mathrm{V}$')
    # plt.title('{}'.format())
    ax1[1, 0].legend()
    ax1[1, 0].grid()
    # ax1[0].set_ylim([0, 340])

    ax1[1, 1].step(t, obs_matrix_denorm[:, 2], 'r', label='i2')
    ax1[1, 1].step(t, obs_matrix_denorm[:, 5], '--b', label='iT2')
    # plt.plot(t,result[:steps,0], label = 'i1')
    ax1[1, 1].set_xlabel(r'$t\,/\,\mathrm{s}$')
    ax1[1, 1].set_ylabel('$i_{\mathrm{}}\,/\,\mathrm{A}$')
    # plt.title('{}'.format())
    ax1[1, 1].legend()
    ax1[1, 1].grid()
    # ax1[0].set_ylim([0, 340])
    fig.suptitle(title, fontsize=14)
    plt.show()

    P_ref = [ref] * len(t)
    P_load = (obs_matrix_denorm[:, 4] + obs_matrix_denorm[:, 5]) * parameter["R_load"]
    # plot power via Load
    plt.step(t, P_load, 'b', label='P_{load}')
    plt.step(t, P_ref, '--', color='gray', label='P_{ref}')
    plt.xlabel(r'$t\,/\,\mathrm{s}$')
    plt.ylabel('$P_{\mathrm{}}\,/\,\mathrm{W}$')
    # plt.title('{}'.format())
    plt.legend()
    plt.grid()
    plt.title(title)

    # ax1[0].set_ylim([0, 340])
    plt.show()


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

# limits for normalization
i_lim = 20
v_lim = 600

limits = dict()
limits['i_lim'] = i_lim
limits['v_lim'] = v_lim

# x0 = np.array([i10, v10, i20, v20, iT10, iT20])
x0 = [0, 0, 0, 0, 0, 0]

vi1 = 0.7  # 200
vi2 = 0.7  # 200

parameter = dict()
parameter['R_source'] = R
parameter['L_source'] = L
parameter['C_source'] = C
parameter['L_cable'] = LT
parameter['R_cable'] = R
parameter['R_load'] = R_load
parameter['V_dc'] = 300

# define CM
CM = np.array([[0, 0, 1],
               [0, 0, 2],
               [-1, -2, 0]])

ref = 200  # W

env = Env_DARE(CM=CM, ts=ts, parameter=parameter, x0=x0, limits=limits, refs=ref)

v1_list = []

obs = env.reset()
obs_matrix = np.array(obs)
v1_list.append(obs[1])
for i in range(num_samples - 1):
    # action, _states = model.predict(obs)
    action = np.array([[vi1], [vi2]])
    # action = np.array([action.squeeze(),action.squeeze()])#np.array([[vi1], [vi2]])
    obs, rewards, dones, info = env.step(action)
    obs_matrix = np.vstack((obs_matrix, obs))
    v1_list.append(obs[1])
    # env.render()

# denormalize!
obs_matrix_denorm = obs_matrix * env.norm_array
plt_obs(obs_matrix_denorm, "Constant action")

"""
Train DDPG agent to supply the load with p_ref by adjusting v_i1,2
"""

# The noise objects for DDPG
n_actions = env.action_space.shape[-1]
action_noise = NormalActionNoise(mean=np.zeros(n_actions), sigma=0.1 * np.ones(n_actions))

# model = DDPG("MlpPolicy", env, action_noise=action_noise, verbose=1)
model = DDPG("MlpPolicy",
             env=env,
             learning_rate=1e-3,
             buffer_size=1_000,  # 1e6
             learning_starts=100,
             batch_size=10,  # 100
             tau=0.005,
             gamma=0.99,
             train_freq=(1, "episode"),
             gradient_steps=-1,
             action_noise=action_noise,
             replay_buffer_class=None,
             replay_buffer_kwargs=None,
             optimize_memory_usage=False,
             tensorboard_log=None,
             create_eval_env=False,
             policy_kwargs=None,
             verbose=0,
             seed=None,
             device="auto",
             _init_setup_model=True,
             )

# test dumm agent
obs = env.reset()
obs_matrix_dumm_ddpg = np.array(obs)
for i in range(num_samples - 1):
    action, _states = model.predict(obs)
    obs, rewards, dones, info = env.step(action)
    obs_matrix_dumm_ddpg = np.vstack((obs_matrix_dumm_ddpg, obs))

obs_matrix_dumm_ddpg_denorm = obs_matrix_dumm_ddpg * env.norm_array
plt_obs(obs_matrix_dumm_ddpg_denorm, "untrained agent")
"""
# From here on training!
# Uncomment to train, not tested yet, just coded

# train
model.learn(total_timesteps=100, log_interval=10)

# test agent
obs = env.reset()
obs_matrix_ddpg = np.array(obs)
for i in range(num_samples - 1):
    action, _states = model.predict(obs)
    obs, rewards, dones, info = env.step(action)
    obs_matrix_ddpg = np.vstack((obs_matrix_ddpg, obs))

obs_matrix_ddpg_denorm = obs_matrix_ddpg * env.norm_array
plt_obs(obs_matrix_ddpg_denorm, "Trained agent 100 steps")

# model.save("ddpg_pendulum")
# env = model.get_env()
# del model # remove to demonstrate saving and loading
# model = DDPG.load("ddpg_pendulum")
"""
