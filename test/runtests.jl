using Dare
using Test
using Logging

global_logger(DareLogger(file_name = split(string(@__FILE__), "\\")[end], add_timestamp = false, log_file = false, log_level = Logging.Error, log_level_file = Logging.Debug))


@testset "NodeConstructor" begin
    # TODO shift to separat file_name
    #include("./nc_test_function.jl")
    # TODO check A, B, C and D matrix with analytical calculation

    #sum of the output of get_fltr_distr should add up to the argument
    for i in 1:4:300
        @test i == sum(get_fltr_distr(i))
    end

    #same for get_loads_distr
    for i in 1:4:300
        @test i == sum(get_loads_distr(i))
    end

    #test lengths of check_parameters output
    #temp = check_parameters(Dict(), 3,2,6) #TODO add CM
    #@test length(temp) == 4
    #@test length(temp["source"])        == 3
    #@test length(temp["source"][1])     == 12
    #@test length(temp["load"])          == 2
    #@test length(temp["load"][1])       == 4
    #@test length(temp["cable"])         == 6
    #@test length(temp["cable"][1])      == 7
    #@test length(temp["grid"])          == 4
end

@testset "Env" begin
    #env_test_function.jl
    
    include("./env_test_state.jl")

    #env = SimEnv()
    #@test size(env.A)              == (33,33)
    #@test size(env.state_space)    == (33,)
    #@test size(env.action_space)   == (6,)

    #start_state = env.state
    #env(ones(size(env.action_space)))
    #Dare.RLBase.reset!(env)
    #@test env.state                 == start_state
end

@testset "Agent" begin
    #env = SimEnv()
    #agent = create_agent_ddpg(na = length(env.action_space), ns = length(env.state_space), use_gpu = false)

    #@test length(agent(env)) == length(env.action_space)
end

@testset "ClassicalController" begin
    
end

@testset "DataHook" begin
    
end


#@testset "env_run_1" begin
#    @test 1==1
#end
#@test π ≈ 3.14 atol=0.01