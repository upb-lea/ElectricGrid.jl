
import gym
import numpy as np
import control


class Env_2_Nodes(gym.Env):

    TIMEOUT = 300

    def __init__(self):
        super(Env_2_Nodes, self).__init__()

        self.number_of_steps = 0

        self.U_INDEX = 0
        # todo check spaces
        self.action_space = gym.spaces.Box(
            low=np.array([0.0]),
            high=np.array([1.0]),
            shape=(1,),
            dtype=np.float32
        )

        self.V_INDEX, self.H_INDEX, self.M_INDEX = 0, 1, 2
        self.observation_space = gym.spaces.Box(
            low=np.array([np.finfo(np.float).min, 0.0, self._r.M1]),
            high=np.array([np.finfo(np.float).max, np.finfo(np.float).max, self._r.M0]),
            dtype=np.float32
        )

        self.reset()


    def step(self, action):

        self.number_of_steps += 1

        # toDo: check all states depending on the limits (i, v)
        limit_crashed = False

        is_done = bool(
            limit_crashed and self._state[self.V_INDEX] < 0 and self._h_max > self._r.H0
        ) or self.number_of_steps >= self.TIMEOUT

        # todo reward design
        if is_done:
            reward = 1
            if self.number_of_steps >= self.TIMEOUT:
                reward = -1
        else:
            reward = 0.0

        return self._simulate(action), reward, is_done, {}


    def _simulate(self,action, normalize=True):
        # normalize
        # solver step

        #todo here,
        T, yout, xout = control.forced_response(sys_1, T=t, U=action, X0=x0, return_x=True, squeeze=True)

        self._state = 0
        state = np.array(self._state)
        if normalize:
            # toDo
            return (state/1)
        else:
            return state

    def reset(self):

        # set x0

        #self._state =

        self.number_of_steps = 0
        # todo how to run the first step without action using control.py?
        return self._simulate()



