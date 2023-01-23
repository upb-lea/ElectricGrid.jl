using Test
using MAT
using Dare

@testset "get_fltr_distr" begin
    for i in [1,14,17,42]
        @test i == sum(get_fltr_distr(i))
    end
end

@testset "get_load_distr" begin
    for i in [1,14,17,42]
        @test i == sum(get_load_distr(i))
    end
end

@testset "CM_generate" begin
    for i in [2,14,17,42]
        num_connections, CM = CM_generate(i, i, 1, 1)
        @test CM[1,end] == -CM[end,1]
        @test tr(CM) == 0
    end
    
    for i in [2,14,17,42]
        num_connections, CM = CM_generate(i, 1, 1, 0)
        @test num_connections == i
        num_connections, CM = CM_generate(i, 2, 1, 0)
        @test num_connections == 2*i
end

@testset "matrices" begin

    file = matread("./test/nc_test.mat")

    parameters = Dict()

    grid_properties = Dict()
    grid_properties["fs"] =  10e3
    grid_properties["v_rms"] = 230
    grid_properties["phase"] = 1;
    parameters["grid"] = grid_properties

    source = Dict()
    source["fltr"] = "LC"
    source["R1"] = 1.1e-3
    source["L1"] = 70e-6
    source["C"] = 250e-6
    source["R_C"] = 7e-3
    source_list = []
    push!(source_list, source)
    parameters["source"] = source_list

    cable = Dict()
    cable["R"] = 1e-3
    cable["L"] = 1e-4
    cable["C"] = 1e-4;
    cable_list = []
    push!(cable_list, cable);
    parameters["cable"] = cable_list
    
    load = Dict()
    load["impedance"] = "RLC"
    load["R"] = 100;
    load["L"] = 1e-2
    load["C"] = 1e-2;
    load["pf"] = 0
    load_list = []

    push!(load_list, load);
    parameters["load"] = load_list;

    S1_L1 = NodeConstructor(num_sources=1, num_loads=1, parameters=parameters);
    A, B, C, D = get_sys(S1_L1);

    @test file["A"] == A
    @test file["B"] == B
    
end