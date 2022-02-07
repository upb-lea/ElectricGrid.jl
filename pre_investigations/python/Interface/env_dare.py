import copy

import gym
import numpy as np
import control
import scipy

from pre_investigations.python.Interface.util import NodeConstructor


class Env_DARE(gym.Env):

    TIMEOUT = 300

    def __init__(self, CM, ts, parameter, x0, time_start=0):

        super(Env_DARE, self).__init__()


        power_grid = NodeConstructor(2, 1, parameter, CM=CM)  # 2 Source with 2 connections
        power_grid.draw_graph()

        A, B, C, D = power_grid.get_sys()
        # discretize
        A_d = scipy.linalg.expm(A * ts)
        A_inv = scipy.linalg.inv(A)
        B_d = A_inv @ (A_d - np.eye(A.shape[0])) @ B
        C_d = copy.copy(C)

        self.time_step_size = ts
        self.sim_time_interval = None
        self.time_start = time_start
        self.x0 = x0
        self.sys_d = control.ss(A_d, B_d, C_d, 0, dt=True)

        self.number_of_steps = 0

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

        is_done = False#bool(
            #limit_crashed and self._state[self.V_INDEX] < 0 and self._h_max > self._r.H0
        #) or self.number_of_steps >= self.TIMEOUT

        # todo reward design
        if is_done:
            reward = 1
            if self.number_of_steps >= self.TIMEOUT:
                reward = -1
        else:
            reward = 0.0

        return self._simulate(action), reward, is_done, {}


    def _simulate(self,action, normalize=True):
        # actions is "doubled" since in sipy ltisys.py L.3505 last is discarded
        # toDo more performant stepwise interaction possible?
        T, yout, xout = control.forced_response(self.sys_d, T=self.sim_time_interval,
                                                       U=np.array([action.squeeze(),action.squeeze()]),
                                                       X0=self._state,
                                                       return_x=True, squeeze=True)

        # get the last solution of the solver (2 are given)
        self._state = xout[:, -1]

        if normalize:
            # toDo
            return (self._state/1)
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



