import json
from os import makedirs
import numpy as np

a = np.dot()
from pre_investigations.python.solver_investigations.timing_evaluation import timing_experiment_simulation

parameter_used = dict()
parameter_used['R_source'] = 0.4
parameter_used['L_source'] = 2.3e-3
parameter_used['C_source'] = 10e-6
parameter_used['L_cable'] = 2.3e-3
parameter_used['R_cable'] = 0.4
parameter_used['R_load'] = 14
parameter_used['V_dc'] = 300

# define required data for Dare_env
limits = dict()
limits['i_lim'] = 20
limits['v_lim'] = 600

debug = False

# hyperparameters:
repeat_used = 3
loops_used = 3
num_cm = 3
t_s_used = 1e-4
t_end_used = [0.03]  # , 0.1]
num_nodes_used = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
                  29, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]#, 105, 110, 115, 120, 125, 130]
#used_methode = ['control.py', 'env_standalone', 'scipy_ode', 'scipy_odeint', 'scipy_solve_ivp', 'env_agent_interaction']
#used_methode_args = ['discrete', '', {'integrator': 'lsoda', 'methode': 'bdf'}, 'LSODA', 'LSODA', '']
used_methode = ['control.py', 'scipy_solve_ivp', 'scipy_solve_ivp']#, 'scipy_solve_ivp']#, 'env_standalone', 'scipy_solve_ivp']
used_methode_args = ['' , 'RK45', 'RK23']#, 'RK23']#, '', 'RK23']

# store HPs to json
python_simulation_HPs = {
    't_s_used': t_s_used,
    't_end_used': t_end_used,
    'num_nodes_used': num_nodes_used,
    'used_methode': used_methode,
    'used_methode_args': used_methode_args,
    'parameter_used': parameter_used,
    'repeat': repeat_used,
    'loops': loops_used
}
makedirs('saves', exist_ok=True)
file_name = 'saves/sim_setting_python_time_experiment'
with open(file_name + '.json', 'w') as outfile: json.dump(python_simulation_HPs, outfile)

time_result = timing_experiment_simulation(repeat=repeat_used, loops=loops_used, num_nodes=num_nodes_used,
                                           t_end=t_end_used, num_cm= num_cm, ts=t_s_used,
                                           methode=used_methode,
                                           methode_args=used_methode_args,
                                           parameter=parameter_used,
                                           limits=limits,
                                           save_data=False, save_folder_name='saves_python',
                                           debug=debug)

if debug:
    with open('saves/time_result_all_python_debug_nodes'+str(num_nodes_used[0])+'t_end'
              +str(t_end_used[0])+'.json', 'w') as outfile: json.dump(time_result, outfile)
else:
    with open('saves/time_result_all_python_control.json', 'w') as outfile: json.dump(time_result, outfile)

