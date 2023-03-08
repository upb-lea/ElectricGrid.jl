using DrWatson
@quickactivate "MicroGridSimWithRL"

using BenchmarkTools
using CUDA

function test(a::CuArray)
    CUDA.@sync begin

        a = sin.(a)

    end 

    a
end

function test(a::AbstractArray)
    a = sin.(a)

    a
end

# @benchmark test(CuArray(1:0.0001:10))

# @benchmark test(collect(1:0.0001:10))
