import json
from os import makedirs

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

# hyperparameters:
repeat_used = 1
loops_used = 1
num_cm = 1
t_s_used = 1e-4
t_end_used = [0.01, 0.1]  # , 0.1]
num_nodes_used = [4, 6, 8]
used_methode = ['control.py', 'env_standalone', 'scipy_ode', 'scipy_odeint', 'scipy_solve_ivp', 'env_agent_interaction']
used_methode_args = ['discrete', '', {'integrator': 'lsoda', 'methode': 'bdf'}, 'LSODA', 'LSODA', '']
#used_methode = ['scipy_ode']
#used_methode_args = [{'integrator': 'lsoda', 'methode': 'bdf'}]

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
                                           debug=True)


with open('saves/time_result_all_python.json', 'w') as outfile: json.dump(time_result, outfile)

