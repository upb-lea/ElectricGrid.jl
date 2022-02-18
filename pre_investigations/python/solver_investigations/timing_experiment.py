import json

import pandas as pd

from pre_investigations.python.solver_investigations.timing_evaluation import plot_result, timing_experiment_simulation

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

file_name = 'sim_time_python_repeats_'+str(repeat_used)+'_loops_'+str(loops_used)

t_s_used = 1e-4
t_end_used = [0.05 , 0.1]

num_nodes_used = [2, 4, 6]

time_result = timing_experiment_simulation(repeat=repeat_used, loops=loops_used, num_nodes=num_nodes_used,
                                           t_end=t_end_used, ts=t_s_used,
                                           methode=['control.py', 'control.py', 'scipy_ode', 'scipy_odeint',
                                                    'env_standalone'],
                                           methode_args=['continuous', 'discrete', 'LSODA', 'LSODA', ''],
                                           parameter=parameter_used)

#time_result_df = pd.DataFrame(time_result['times_mean'])
# here, instead of orint = 'index' better store only result arrays in pd.DF, most of config stuff
#time_result_df.to_pickle(file_name+'_2.pkl.bz2')

#with open(file_name + '.json', 'w') as outfile:
#    json.dump(time_result, outfile)

# todo: store data-arrays with result to pkl.bz2 and rest to json - files connected via good naming
#  (better: postgresql instead of json with unique ID)


plot_result(time_result, num_nodes_used, t_end_used)
