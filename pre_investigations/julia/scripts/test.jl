using DrWatson
@quickactivate "MicroGridSimWithRL"

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