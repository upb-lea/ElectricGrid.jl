using DrWatson
@quickactivate "MicroGridSimWithRL"

# using DifferentialEquations
# using Sundials
using Plots
# using LinearAlgebra
# using ControlSystems
# using BenchmarkTools
using ReinforcementLearning
using Flux
using StableRNGs
# using IntervalSets

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent.jl"))

CM = [ 0.  1.
        -1.  0.]

parameters = Dict()
# LC filter
parameters["source"] = [Dict("fltr" => "LC", "R" => 0.4, "L1" => 2.3e-3, "C" => 10e-6)]
parameters["cable"] = [Dict("R" => 0.722, "L" => 0.955e-3, "C" => 8e-09)]
parameters["load"] = [Dict("impedance" => "R", "R" => 14)]

nc = NodeConstructor(num_source=1, num_loads=1, CM=CM, parameters=parameters)

#draw_graph(Grid_FC)   ---   not yet implemented

A, B, C, D = get_sys(nc)

limits = Dict("i_lim" => 20, "v_lim" => 600)

norm_array = vcat([limits[i] for j = 1:nc.num_source for i in ["i_lim", "v_lim"]], [limits["i_lim"] for i = 1:nc.num_connections] )
norm_array = vcat( norm_array, [limits["v_lim"] for i = 1:nc.num_loads] )

ns = length(A[1,:])
na = length(B[1,:])

# time step
ts = 1e-5


V_source = 230
# V_rms
v_dc = V_source / sqrt(2)

global env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, v_dc=v_dc, ts=rationalize(ts))
global agent = create_agent(na, ns)

# ----------------------------------------------------------------------------------------
function execute_env(env::SimEnv, agent::Agent, t_len::Int, debug::Bool)
    if debug
        output = zeros(length(env.Ad[1,:]),t_len+1)
    else
        output = 0.0
    end

    RLBase.reset!(env)
    for i = 1:t_len
        action = agent(env)
        env(action)
        if debug output[:,i+1] = env.state.*env.norm_array end
    end

    return output
end

result = execute_env(env, agent, 300, true)


display(plot(result[4, :]))

