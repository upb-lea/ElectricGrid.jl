using DrWatson
@quickactivate "MicroGridSimWithRL"

using PyCall
@pyinclude(srcdir("nodeconstructor.py"))

using ControlSystems
using JSON
using Profile
using PProf
using LinearAlgebra
using BenchmarkTools
using CUDA

ts=1e-4
node = 20
cm = 1
num_threads = 8
useCUDA = true

if !useCUDA
    ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
    BLAS.set_num_threads(num_threads)
end

parameter=Dict()
parameter["R_source"] = 0.4
parameter["L_source"] = 2.3e-3
parameter["C_source"] = 10e-6
parameter["L_cable"] = 2.3e-3
parameter["R_cable"] = 0.4
parameter["R_load"] = 14
parameter["V_dc"] = 300

results = Dict()

CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(node) * ".json"))

function prepareCM(n)
    global CM = reduce(hcat, CM_list[n])'
    global CM = convert(Matrix{Int}, CM)
    global nc = py"NodeConstructor"(node, node, parameter, CM=CM)

    global A, B, C, D = nc.get_sys()
    global ns = length(A[1,:])
    global na = length(B[1,:])

    global t = collect(0:ts:0.03)

    global hallo = CuArray{Float32}(undef, (1,2));
    global hallo .= 5

    global Ad = exp(A*ts)
    global Bd = A \ (Ad - C) * B
    global sys_d = ss(Ad, Bd, C, D, ts)

    global x0 = [0.0 for i = 1:ns]
    global u = [230.0 for i = 1:length(t)]
    global uu = [u for i = 1:na ]
    global uuu = mapreduce(permutedims, vcat, uu)
    global ttt = t

    if useCUDA
        #=
        global Ad = CuArray(Ad)
        global Bd = CuArray(Bd)
        global C = CuArray(C)
        global D = CUDA.zeros(Float64, size(Bd))
        global sys_d = HeteroStateSpace(Ad, Bd, C, D, ts)

        global x0 = CuArray(x0)
        global uuu = CuArray(uuu)
        global ttt = CuArray(ttt)
        =#
        
        global Ad = CuArray(Float32.(Ad))
        global Bd = CuArray(Float32.(Bd))
        global C = CuArray(Float32.(C))
        global D = CUDA.zeros(Float32, size(Bd))
        global sys_d = HeteroStateSpace(Ad, Bd, C, D, ts)

        global x0 = CuArray(Float32.(x0))
        global uuu = CuArray(Float32.(uuu))
        #global ttt = CuArray(Float32.(ttt))
    end
    return nothing
end

@views function ltitr2(A::AbstractMatrix, B::AbstractMatrix, u::AbstractVecOrMat,
    x0::AbstractVecOrMat=zeros(eltype(A), size(A, 1)))

    x = nothing
    for i = 1:1
        T = promote_type(LinearAlgebra.promote_op(LinearAlgebra.matprod, eltype(A), eltype(x0)),
                        LinearAlgebra.promote_op(LinearAlgebra.matprod, eltype(B), eltype(u)))

        n = size(u, 2)

        # Using similar instead of Matrix{T} to allow for CuArrays to be used.
        # This approach is problematic if x0 is sparse for example, but was considered
        # to be good enough for now
        x = similar(x0, T, (length(x0), n))

        x[:,1] .= x0
        mul!(x[:, 2:end], B, u[:, 1:end-1]) # Do all multiplications B*u[:,k] to save view allocations

        for k=1:n-1
            mul!(x[:, k+1], A, x[:,k], true, true)
        end
        
    end
    return x
end

function investigate(sys_d,uuu,ttt,x0)
    CUDA.@sync begin
        result, _, _, _ = lsim(sys_d,uuu,ttt,x0=x0);
    end
    return nothing
end

function investigate2(Ad,Bd,uuu,x0)
    CUDA.@sync begin
        result = ltitr2(Ad,Bd,uuu,x0)
    end
    return nothing
end

#investigate(sys_d,uuu,ttt,x0)
#@profile investigate(sys_d,uuu,ttt,x0)
#pprof(;webport=58699)
#mul!(x0,B,uuu[:,1])

prepareCM(cm)

#@benchmark investigate2(Ad,Bd,uuu,x0)
CUDA.allowscalar(false)
investigate(sys_d,uuu,ttt,x0)
@benchmark investigate(sys_d,uuu,ttt,x0)