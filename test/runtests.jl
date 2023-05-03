using ElectricGrid
using Test
using Logging

global_logger(CustomLogger(file_name=split(string(@__FILE__), "\\")[end], add_timestamp=false, log_file=false, log_level=Logging.Error, log_level_file=Logging.Debug))


@testset "NodeConstructor" begin

    include("./nc_test_function.jl")

end

@testset "Env" begin

    include("./env_test_basics.jl")
    include("./env_test_state.jl")

end

@testset "Agent" begin
    #env = ElectricGridEnv()
    #agent = CreateAgentDdpg(na = length(env.action_space), ns = length(env.state_space), use_gpu = false)

    #@test length(agent(env)) == length(env.action_space)
end

@testset "ClassicalController" begin

    include("./classic_control_test.jl")

end

@testset "DataHook" begin

end

@testset "Optimizer" begin
    include("./env_test_optimizer_status.jl")
end

@testset " P Q Test" begin
    include("./power_flow_test_P_Q.jl")
end
