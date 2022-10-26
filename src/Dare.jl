module Dare

using ReinforcementLearning
using PlotlyJS

#export create_setup, Classical_Policy, create_agent_ddpg, Source_Initialiser, MultiAgentGridController, DataHook, plot_hook_results, plot_best_results, NodeConstructor, dare_setup, SimEnv

#code to export all, taken from https://discourse.julialang.org/t/exportall/4970/18
macro make_public(module_name::Symbol)

    eval(Meta.parse("export " * string(module_name)))

    as_module = eval(module_name)
    @assert as_module isa Module

    for name in names(as_module; all = true)
        if (string(name)[1] != '#')
            #println("export " * string(name))
            as_module.eval(Meta.parse("export " * string(name)))
        end
    end

    return nothing
end
export make_public
@make_public Dare

include("./nodeconstructor.jl")
include("./env.jl")
include("./agent_ddpg.jl")
include("./Classical_Control.jl")
include("./Power_System_Theory.jl")
include("./MultiAgentGridController.jl")
include("./plotting.jl")
include("./data_hook.jl")
include("./Dare_Wrapper.jl")

end # module