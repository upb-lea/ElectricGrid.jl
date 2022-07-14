using DrWatson
@quickactivate "MicroGridSimWithRL"

<<<<<<< HEAD
ENV["JULIA_DEBUG"] = "CUDA"

using CUDA

#CUDA.version()
#CUDA.versioninfo()

#println(has_device())


dev = CuDevice(0)
ctx = CuContext(dev)
synchronize(ctx)

println(has_device())

#buf = Mem.alloc(Device, ...)

event = CuEvent()
record(event)

stream = CuStream()

hallo = CUDA.fill(1.0f0, 3)
println(has_device())
hallo .= 5

a = 12
=======
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
>>>>>>> develop
