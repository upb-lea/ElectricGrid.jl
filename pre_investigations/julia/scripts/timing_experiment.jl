using DrWatson
@quickactivate "MicroGridSimWithRL"

using JSON
include(srcdir("timing_evaluation.jl"))

#r = timing_experiment_simulation(5, 5, [2,4,6,8,10], [0.001,0.01,0.03], 3, 1e-4, ["control", "env_without_agent", "env_with_agent", "lsoda"], Dict(), false)
#r = timing_experiment_simulation(3, 3, [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30], [0.03], 6, 1e-4, ["control", "env_without_agent"], Dict(), false)
#r = timing_experiment_simulation(3, 3, [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], [0.003], 1, 1e-4, ["env_without_agent"], Dict(), false)
r = timing_experiment_simulation(1, 1, [8], [0.03], 30, 1e-4, ["control"], Dict(), false)

#open(datadir(savename("julia_", r, "json", allowedtypes = (Vector,), ignores=["info","times_mean","times_std","num_samples"])),"w") do f 
open(datadir(r["name"] * ".json"),"w") do f 
    write(f, JSON.json(r)) 
end