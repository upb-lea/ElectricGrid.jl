import json

import pandas as pd

from pre_investigations.python.solver_investigations.timing_evaluation import plot_result

# JSON file
f = open('saves_python/sim_setting_python_time_experiment.json', "r")

# Reading from file
data = json.loads(f.read())

# Future: better store
time_result = pd.read_pickle('saves_python/time_result_all.pkl.bz2')

time_result_dict = {}

for i in range(time_result.shape[0]):
    time_result_dict[time_result[0][i]] = time_result[1][i]

plot_result(time_result_dict, data['num_nodes_used'], data['t_end_used'])
