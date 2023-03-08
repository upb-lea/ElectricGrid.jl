using JEG
using Test
using Logging

global_logger(JEGLogger(file_name = split(string(@__FILE__), "\\")[end], add_timestamp = false, log_file = false, log_level = Logging.Error, log_level_file = Logging.Debug))


@testset "NodeConstructor" begin

    include("./nc_test_function.jl")

end

@testset "Env" begin
    #env_test_function.jl

    include("./env_test_basics.jl")
    include("./env_test_state.jl")

    #env = ElectricGridEnv()
    #@test size(env.A)              == (33,33)
    #@test size(env.state_space)    == (33,)
    #@test size(env.action_space)   == (6,)

    #start_state = env.state
    #env(ones(size(env.action_space)))
    #JEG.RLBase.reset!(env)
    #@test env.state                 == start_state
end

@testset "Agent" begin
    #env = ElectricGridEnv()
    #agent = CreateAgentDdpg(na = length(env.action_space), ns = length(env.state_space), use_gpu = false)

    #@test length(agent(env)) == length(env.action_space)
end

@testset "ClassicalController" begin

    include("./classic_control_test.jl")

end

@testset "data_hook" begin

end

# Need to define where this should go
@testset "Optimizer" begin
    include("./env_test_optimizer_status.jl")
end

#@testset "env_run_1" begin
#    @test 1==1
#end
#@test π ≈ 3.14 atol=0.01
