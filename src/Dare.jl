module Dare

using ReinforcementLearning
using PlotlyJS

export create_setup, Classical_Policy, create_agent_ddpg, Source_Initialiser, MultiAgentGridController, DataHook, plot_hook_results, plot_best_results, NodeConstructor, dare_setup, SimEnv

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
