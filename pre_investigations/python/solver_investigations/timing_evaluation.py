import copy
import json

import time
import timeit
from os import makedirs

import control
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.integrate import ode, odeint, solve_ivp
from stable_baselines3 import DDPG
from stable_baselines3.common.noise import NormalActionNoise

from pre_investigations.python.Interface.env_dare import Env_DARE
from pre_investigations.python.dare.utils.nodeconstructor import NodeConstructor


# time measure based on https://note.nkmk.me/en/python-timeit-measure/
from pre_investigations.python.solver_investigations.visualize import plot_result


def timing_experiment_simulation(repeat: int = 5, loops: int = 10, num_nodes: list[int] = np.arange(2, 12, 2).tolist(),
                                 t_end: list[float] = np.arange(2, 12, 2).tolist(), num_cm: int = 1, ts: float = 1e-4,
                                 methode=None,
                                 methode_args=None,
                                 parameter: dict = None,
                                 limits: dict = None,
                                 save_data: bool = False,
                                 save_folder_name='saves',
                                 debug=False):

    steps = int(1 / ts)
    makedirs(save_folder_name, exist_ok=True)
    # pathlib.Path(self.save_folder.mkdir(exist_ok=True))

    t_result_mean = np.zeros([len(methode), len(num_nodes), len(t_end)])
    t_result_std = np.zeros([len(methode), len(num_nodes), len(t_end)])

    ####################################################################################################################
    # Define different functions for time-meas/debug

    def env_model_ode(t, x):
        #u_ode = np.array([230] * power_grid.num_source)   # use this for debug in case of consant input
        return A_sys @ x + B_sys @ u_ode_constant
        # for random action:
        #return A_sys @ x + B_sys @ u_random_vec[:, int(t / ts)]

    def run_scipy_ode(t_end_ode, ode_solver_):
        ode_solver_.set_initial_value(x0, 0)
        while ode_solver_.successful() and ode_solver_.t <= steps * ts * t_end_ode:
            ode_solver_.integrate(ode_solver_.t + ts)

    def run_scipy_ode_debug(t_end_ode, ode_solver_):
        res_list = []
        ode_solver_.set_initial_value(x0, 0)
        while ode_solver_.successful() and ode_solver_.t <= steps * ts * t_end_ode:
            result = ode_solver_.integrate(ode_solver_.t + ts)
            res_list.append(result[1])
        plt.plot(res_list)
        plt.title('ODE')
        plt.show()
        time.sleep(0.5)

    def run_odeint_debug(env_model_ode, x0, t):
        res_list = odeint(env_model_ode, x0, t, tfirst=True)
        plt.plot(t, res_list[:steps, 1], label='v1')
        plt.title('odeint')
        plt.show()
        time.sleep(0.5)

    def run_solve_ivp_debug(env_model_ode, t_int, x0, t_eval,
                      method):
        result_ivp = solve_ivp(env_model_ode, t_int, x0, t_eval=t_eval, method=method)
        plt.plot(result_ivp.y[1], label='v1')
        plt.title('solve_ivp')
        plt.show()
        time.sleep(0.5)

    def run_control_debug(sys, t, u_random, x0):
        T, yout_d, xout_d = control.forced_response(sys, T=t, U=u_random, X0=x0, return_x=True, squeeze=True)
        plt.step(t, xout_d[1], 'r', label='v1_discret')
        plt.title('Control')
        plt.show()
        time.sleep(0.5)

    def execute_env(env, num_samples, u_env):
        # env.reset()
        for i in range(num_samples - 1):
            _, _, _, _ = env.step(u_env[i])

    def execute_env_debug(env, num_samples, u_env):
        v1_list = []

        obs = env.reset()
        v1_list.append(obs[1])

        for i in range(num_samples - 1):
            # u_random = np.random.uniform(-1, 1, env.action_space.shape[0])
            obs, _, _, _ = env.step(u_env[i])
            v1_list.append(obs[1] * limits['v_lim'])
        plt.step(v1_list, 'r', label='v1_discret')
        plt.title('env_standalone')
        plt.show()
        time.sleep(0.5)

    def execute_agent_in_env(agent, env, num_samples):
        obs = env.reset()
        for i in range(num_samples - 1):
            action, _states = agent.predict(obs)
            obs, _, _, _ = env.step(action)

    def execute_agent_in_env_debug(agent, env, num_samples):
        v1_list = []
        obs = env.reset()
        v1_list.append(obs[1])
        for i in range(num_samples - 1):
            action, _states = agent.predict(obs)
            obs, _, _, _ = env.step(action)
            v1_list.append(obs[1] * limits['v_lim'])
        plt.step(v1_list, 'r', label='v1_discret')
        plt.title('env_agent_interaction')
        plt.show()
        time.sleep(0.5)

    for n in range(len(methode)):

        for k in range(len(num_nodes)):

            # load Cm (num_nodes[k]) [c]
            f = open('CM_matrices/CM_nodes' + str(num_nodes[k]) + '.json', "r")
            CM_list = json.loads(f.read())

            for l in range(len(t_end)):
                # define time vector
                t = np.arange(0, t_end[l] + ts, ts)
                time_list = []

                for c in range(num_cm):

                    CM = np.array(CM_list[c])

                    if methode[n] in ['env_standalone', 'env_agent_interaction']:
                        env = Env_DARE(num_sources=num_nodes[k], num_loads=num_nodes[k], CM=CM, ts=ts,
                                       parameter=parameter, x0=None, limits=limits, refs=200)

                        if methode[n] in ['env_agent_interaction']:
                            n_actions = env.action_space.shape[-1]
                            action_noise = NormalActionNoise(mean=np.zeros(n_actions), sigma=0.1 * np.ones(n_actions))
                            agent = DDPG("MlpPolicy", env, action_noise=action_noise, verbose=1)

                    else:

                        power_grid = NodeConstructor(num_nodes[k], num_nodes[k], parameter, CM=CM)
                        A_sys, B_sys, C_sys, D_sys = power_grid.get_sys()

                        if methode[n] in ['control.py'] and methode_args[n] == 'discrete':
                            A_d = scipy.linalg.expm(A_sys * ts)
                            A_inv = scipy.linalg.inv(A_sys)
                            B_d = A_inv @ (A_d - np.eye(A_sys.shape[0])) @ B_sys
                            C_d = copy.copy(C_sys)

                            sys = control.ss(A_d, B_d, C_d, 0, dt=True)

                        else:
                            sys = control.ss(A_sys, B_sys, C_sys, D_sys)

                        # generate init state
                        x0 = np.zeros((A_sys.shape[0],))

                    ############################################################################
                    # predefine actions to not measure the time for that
                    u_random_vec = np.random.uniform(-1, 1, (power_grid.num_source, len(t) + 1)) * parameter['V_dc']
                    u_ode_constant = np.array([230] * power_grid.num_source)

                    if methode[n] in ['scipy_ode']:
                        if methode_args[n]['integrator'] == 'lsoda':
                            ode_solver = ode(env_model_ode).set_integrator(methode_args[n]['integrator'], atol=1e-10,
                                                                           rtol=1e-6)
                        else:
                            ode_solver = ode(env_model_ode).set_integrator(methode_args[n]['integrator'],
                                                                           methode=methode_args[n]['methode'],
                                                                           atol=1e-10, rtol=1e-6)
                        ode_solver.set_initial_value(x0, 0)

                    if methode[n] in ['control.py']:
                        if debug:
                            u_fix = np.array([230] * power_grid.num_source)[:, None] * np.ones(
                                (power_grid.num_source, len(t)))
                            #u_random = np.random.uniform(-1, 1, (power_grid.num_source, len(t))) * parameter['V_dc']
                            res_list = timeit.repeat(
                                lambda: run_control_debug(sys, t, u_fix, x0), repeat=repeat, number=loops)
                        else:
                            u_fix = np.array([230] * power_grid.num_source)[:, None] * np.ones(
                                (power_grid.num_source, len(t)))
                            #u_random = np.random.uniform(-1, 1, (power_grid.num_source, len(t))) * parameter['V_dc']
                            res_list = timeit.repeat(
                                lambda: control.forced_response(sys, T=t, U=u_fix, X0=x0, return_x=True, squeeze=True)
                                , repeat=repeat, number=loops)

                    if methode[n] in ['scipy_ode']:
                        if debug:
                            res_list = timeit.repeat(lambda: run_scipy_ode_debug(t_end[l], ode_solver), repeat=repeat,
                                                     number=loops)
                        else:
                            res_list = timeit.repeat(lambda: run_scipy_ode(t_end[l], ode_solver), repeat=repeat,
                                                     number=loops)

                    if methode[n] in ['scipy_odeint']:
                        if debug:
                            res_list = timeit.repeat(lambda: run_odeint_debug(env_model_ode, x0, t), repeat=repeat,
                                                     number=loops)
                        else:
                            res_list = timeit.repeat(lambda: odeint(env_model_ode, x0, t, tfirst=True), repeat=repeat,
                                                     number=loops)

                    if methode[n] in ['scipy_solve_ivp']:
                        if debug:
                            res_list = timeit.repeat(lambda: run_solve_ivp_debug(env_model_ode, [0, t_end[l]], x0, t_eval=t,
                                                                           method=methode_args[n]), repeat=repeat,
                                                     number=loops)
                        else:
                            res_list = timeit.repeat(lambda: solve_ivp(env_model_ode, [0, t_end[l]], x0, t_eval=t,
                                                                       method=methode_args[n]), repeat=repeat,
                                                     number=loops)

                    if methode[n] in ['env_standalone']:
                        num_samples = int(t_end[l] / ts)
                        env.reset()
                        if debug:
                            u_vec = np.ones([num_samples, env.action_space.shape[0]]) * 0.76666
                            res_list = timeit.repeat(lambda: execute_env_debug(env, num_samples, u_vec), repeat=repeat,
                                                     number=loops)
                        else:
                            # u_vec = np.random.uniform(-1, 1, (num_samples, env.action_space.shape[0]))
                            u_vec = np.ones([num_samples, env.action_space.shape[0]]) * 0.76666
                            res_list = timeit.repeat(lambda: execute_env(env, num_samples, u_vec), repeat=repeat,
                                                     number=loops)

                    if methode[n] in ['env_agent_interaction']:
                        num_samples = int(t_end[l] / ts)
                        if debug:
                            res_list = timeit.repeat(lambda: execute_agent_in_env_debug(agent, env, num_samples),
                                                     repeat=repeat,
                                                     number=loops)
                        else:
                            res_list = timeit.repeat(lambda: execute_agent_in_env(agent, env, num_samples), repeat=repeat,
                                                     number=loops)
                    # correct for loops-time execution:
                    result_per_loop = ([x / loops for x in res_list])

                    time_list.append(result_per_loop)

                t_result_mean[n][k][l] = np.mean(time_list)
                t_result_std[n][k][l] = np.std(time_list)

                print("done for " + str(n) + ", " + str(k) + " and " + str(l))

                # save time_array to pkl
        if save_data:
            time_result_df = pd.DataFrame(t_result_mean[n][:][:], index=[str(int) for int in num_nodes],
                                          columns=[str(int) for int in t_end])
            if methode[n] in ['scipy_ode']:
                time_result_df.to_pickle(
                    save_folder_name + '/sim_time_mean_' + methode[n] + methode_args[n]['methode'] +
                    '_repreat' + str(repeat) + '_loops' + str(loops) + '.pkl.bz2')
            else:
                time_result_df.to_pickle(save_folder_name + '/sim_time_mean_' + methode[n] + methode_args[n] +
                                         '_repreat' + str(repeat) + '_loops' + str(loops) + '.pkl.bz2')

    result_dict = dict()
    result_dict['methode'] = methode
    result_dict['methode_args'] = methode_args
    result_dict['times_mean'] = t_result_mean.tolist()
    result_dict['times_std'] = t_result_std.tolist()
    result_dict['t_end'] = t_end
    result_dict['num_grid_nodes'] = num_nodes
    # result_dict['t_s'] = ts
    # result_dict['repeats'] = repeat
    # result_dict['loops'] = loops
    result_dict['info'] = 'Logs the mean and std of the execution time for all defined methods to simulate the power ' \
                          'grid for different simulation times (t_end) and grid size (num_grid_nodes). num_grid_nodes ' \
                          'thereby defines the number of sources and the number of loads ' \
                          '(grid size = 2*num_grid_nodes). ' \
                          'Each experiment is executed loops*repeats-times ' \
                          'while the mean and std is calculated based' \
                          'on repeats'

    return result_dict


if __name__ == '__main__':
    """
    Small execution exxample of the two functions
    """
    parameter_used = dict()
    parameter_used['R_source'] = 0.4
    parameter_used['L_source'] = 2.3e-3
    parameter_used['C_source'] = 10e-6
    parameter_used['L_cabel'] = 2.3e-3
    parameter_used['R_cabel'] = 0.4
    parameter_used['R_load'] = 14
    parameter_used['V_dc'] = 300

    repeat_used = 2
    loops_used = 3

    t_s_used = 1e-4
    t_end_used = [0.001, 0.002, 0.004]  # ,0.005, 0.01, 0.05]

    num_nodes_used = [2, 4, ]  # 6, 7, 10]

    """
    time_result = timing_experiment_simulation(repeat=repeat_used, loops=loops_used, num_nodes=num_nodes_used,
                                               t_end=t_end_used, ts=t_s_used,
                                               methode=['control.py', 'control.py', 'scipy_ode', 'scipy_odeint',
                                                        'env_standalone', 'env_agent_interaction'],
                                               methode_args=['continuous', 'discrete', 'LSODA', 'LSODA', '', ''],
                                               parameter=parameter_used)
    """
    time_result = timing_experiment_simulation(repeat=repeat_used, loops=loops_used, num_nodes=num_nodes_used,
                                               t_end=t_end_used, ts=t_s_used,
                                               methode=['env_standalone'],
                                               # ['scipy_solve_ivp'],#['scipy_ode'],  #
                                               methode_args=[''],  # {'integrator': 'lsoda', 'methode': 'lsoda'}],
                                               # methode_args=[{'integrator': 'lsoda', 'methode': 'bdf'}],
                                               parameter=parameter_used,
                                               save_data=True, save_folder_name='saves_python', debug=True)

    used_methode = ['control.py', 'scipy_ode', 'scipy_odeint', 'scipy_solve_ivp', 'env_standalone',
                    'env_agent_interaction']
    used_methode_args = ['discrete', {'integrator': 'lsoda', 'methode': 'bdf'}, 'LSODA', 'LSODA', '', '']

    # possible combinations for scipy.ode:
    # [{'integrator': 'lsoda', 'methode': 1}]
    # [{'integrator': 'vode', 'methode': 'adam'}]
    # [{'integrator': 'vode', 'methode': 'bdf'}]
    # [{'integrator': 'zvode', 'methode': 'bdf'}]

    plot_result(time_result, num_nodes_used, t_end_used)
