using Test
using MAT
using Dare
using LinearAlgebra

@testset "get_fltr_distr" begin
    for i in [1,14,17,42]
        @test i == sum(get_fltr_distr(i))
    end
end

@testset "get_loads_distr" begin
    for i in [1,14,17,42]
        @test i == sum(get_loads_distr(i))
    end
end

@testset "CM_generate" begin
    for i in [2,14,17,42]
        num_connections, CM = CM_generate(i, i, 1, 1, 1)
        @test CM[1,end] == -CM[end,1]
        @test tr(CM) == 0
    end

    for i in [2,14,17,42]
        num_connections, CM = CM_generate(i, 1, 1, 0, 0)
        @test num_connections == i
        num_connections, CM = CM_generate(i, 2, 1, 0, 0)
        @test num_connections == 2*i
    end
end


# 3 Sources (LCl, LC and L) 1 Load (RLC)
@testset "matrices" begin

    file = matread("./test/nc_test_2.mat")

    # Filter
    R_1 = 1.1e-3;
    L_1 = 70e-6;
    R_c = 7e-3;
    C_1 = 250e-6;

    # Cable
    C_b = 1e-4;
    L_b = 1e-4;
    R_b = 1e-3;

    # Load
    R_l = 100;
    C_l = 1e-2;
    L_l = 1e-2;

    parameters = Dict()

    grid_properties = Dict()
    grid_properties["fs"] =  10e3
    grid_properties["v_rms"] = 230
    grid_properties["phase"] = 1;
    parameters["grid"] = grid_properties

    source0 = Dict()
    source1 = Dict()
    source2 = Dict()

    source_list = []

    source0["fltr"] = "L"
    source0["R1"] = R_1
    source0["L1"] = L_1

    source1["fltr"] = "LC"
    source1["R1"] = R_1
    source1["L1"] = L_1
    source1["C"] = C_1
    source1["R_C"] = R_c
    # parameters["source"] = source_list

    source2["fltr"] = "LCL"
    source2["R1"] = R_1
    source2["L1"] = L_1
    source2["C"] = C_1
    source2["R_C"] = R_c
    source2["R2"] = R_1
    source2["L2"] = L_1
    push!(source_list, source2, source1, source0)

    parameters["source"] = source_list

    cable = Dict()
    cable["R"] = R_b
    cable["L"] = L_b
    cable["C"] = C_b
    cable_list = []

    push!(cable_list, cable, cable, cable);
    parameters["cable"] = cable_list

    load = Dict()
    load["impedance"] = "RLC"
    load["R"] = R_l;
    load["L"] = L_l;
    load["C"] = C_l;
    load["pf"] = 0
    load_list = []

    push!(load_list, load);
    parameters["load"] = load_list;

    S3_L1 = NodeConstructor(num_sources=3, num_loads=1, parameters=parameters, S2L_p = 1, S2S_p = 0.0);

    S3_L1.parameters["source"][1]

    A, B, C, D = get_sys(S3_L1);

    @test file["A"] ≈ A atol=1e-9
    @test file["B"] ≈ B atol=1e-9

end
