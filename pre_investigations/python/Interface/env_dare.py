import copy

import control
import gym
import numpy as np
import scipy

from pre_investigations.python.dare.utils.nodeconstructor import NodeConstructor


class Reward:
    def __init__(self, params, limits, refs, gamma=0):
        self.params = params
        self.limits = limits
        self.refs = refs
        self.gamma = gamma

    def rew_function(self, obs):
        """

        :param obs: (if normalized=True) normalized observations
        :return: reward to get specific power via load
        """

        # hard coded!!! toDO make flexible depending on e.g. history
        # P_load = ((obs[4] + obs[5])**2)*self.params['R_load']

        # denormalize
        i_T1 = obs[4] * self.limits['i_lim']
        i_T2 = obs[5] * self.limits['i_lim']
        P_load = ((i_T1 + i_T2) ** 2) * self.params['R_load']

        # to normalize, we take the Power by i_max through the load
        P_max = (2 * self.limits['i_lim']) ** 2 * self.params['R_load']

        # if error = P_max -> rew = -1; if error = 0 -> rew = +1
        # gamma normalization for training
        rew = (1 - 2 * np.sqrt(np.abs(P_load - self.refs) / P_max)) * (1 - self.gamma)

        return rew


class Env_DARE(gym.Env):
    TIMEOUT = 300

    def __init__(self, num_sources=2, num_loads=1, CM=None, ts=1e-4, parameter=None, x0=None, limits=None, refs=None, gamma=0, time_start=0):
        """

        :param num_sources:
        :param num_loads:
        :param CM: if set to None, creates random grid
        :param ts:
        :param parameter:
        :param x0: if is None zero array created depending on randomly (?) generated grid
        :param limits:
        :param refs:
        :param gamma:
        :param time_start:
        """


        # toDo shift gamma to env wrapper (or kwargs?)
        super(Env_DARE, self).__init__()

        power_grid = NodeConstructor(num_sources, num_loads, parameter, CM=CM)  # 2 Source with 2 connections
        # power_grid.draw_graph()

        A, B, C, D = power_grid.get_sys()
        # discretize
        A_d = scipy.linalg.expm(A * ts)
        A_inv = scipy.linalg.inv(A)
        B_d = A_inv @ (A_d - np.eye(A.shape[0])) @ B
        C_d = copy.copy(C)

        if x0 is None:
            self.x0 = np.zeros((A_d.shape[0],))
        else:
            self.x0 = x0
        self.time_step_size = ts
        self.sim_time_interval = None
        self.time_start = time_start
        self.sys_d = control.ss(A_d, B_d, C_d, 0, dt=True)
        self.v_dc = parameter['V_dc']

        self.number_of_steps = 0
        self.refs = refs

        if limits is None:
            self.i_lim = None
            self.v_lim = None
        else:
            self.i_lim = limits['i_lim']
            self.v_lim = limits['v_lim']
            # HINT: Currently only LC filters are considered! x is sorted depending on node constructor
            # Therefore, limits are sorted in dependence of power_grid [sources (i, v),..., transitions (i - since
            # only RL connections are considered jet)]
            self.norm_array = np.array([self.i_lim, self.v_lim]*power_grid.num_source +
                                       [self.i_lim]*power_grid.num_connections)

        self.rew = Reward(parameter, limits, self.refs, gamma)

        self.action_space = gym.spaces.Box(
            low=-1,
            high=1,
            shape=(power_grid.num_source,),
            dtype=np.float32
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf,
            high=np.inf,
            shape=(A.shape[0],),
            dtype=np.float32
        )

    def step(self, action):

        self.number_of_steps += 1

        # toDo: check all states depending on the limits (i, v)
        limit_crashed = False

        self.sim_time_interval = np.array([self.current_timestep, self.current_timestep + self.time_step_size])

        self.current_timestep += self.time_step_size

        is_done = False  # bool(
        # limit_crashed and self._state[self.V_INDEX] < 0 and self._h_max > self._r.H0
        # ) or self.number_of_steps >= self.TIMEOUT

        # todo reward design
        if is_done:
            reward = 1
            if self.number_of_steps >= self.TIMEOUT:
                reward = -1
        else:
            reward = 0.0

        obs = self._simulate(action * self.v_dc)

        reward = self.rew.rew_function(obs)

        return obs, reward, is_done, {}

    def _simulate(self, action, normalize=True):
        # actions is "doubled" since in sipy ltisys.py L.3505 last is discarded
        # toDo more performant stepwise interaction possible?
        T, yout, xout = control.forced_response(self.sys_d, T=self.sim_time_interval,
                                                U=np.array([action.squeeze(), action.squeeze()]).T,
                                                X0=self._state,
                                                return_x=True, squeeze=True)

        # get the last solution of the solver (2 are given)
        self._state = xout[:, -1]

        if normalize:
            # toDo
            return (self._state / self.norm_array)
        else:
            return self._state

    def reset(self):
        """
        Resets env to initial state x0 and starttime
        :return: Environment State
        """

        self._state = self.x0
        self.current_timestep = self.time_start

        self.number_of_steps = 0
        return self._state
