using DrWatson
@quickactivate "MicroGridSimWithRL"

using JSON
include(srcdir("timing_evaluation.jl"))

#r = timing_experiment_simulation(5, 5, [2,4,6,8,10], [0.001,0.01,0.03], 3, 1e-4, ["control", "env_without_agent", "env_with_agent", "lsoda"], Dict(), false)
r = timing_experiment_simulation(3, 3, [2,4,6,8,10], [0.01], 1, 1e-4, ["control", "env_without_agent", "env_with_agent", "lsoda"], Dict(), false)

open(datadir(savename("julia_", r, "json", ignores=["info","times_mean","times_std","num_samples"])),"w") do f 
    write(f, JSON.json(r)) 
end