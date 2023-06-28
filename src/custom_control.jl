# ?
function CustomLsim(sys::AbstractStateSpace, u::AbstractVecOrMat, t::AbstractVector;
    x0::AbstractVecOrMat=zeros(Bool, nstates(sys)), method::Symbol=:zoh)
    ny, nu = size(sys)
    nx = sys.nx

    if length(x0) != nx
        error("size(x0) must match the number of states of sys")
    end
    if size(u) != (nu, length(t))
        error("u must be of size (nu, length(t))")
    end

    dt = Float64(t[2] - t[1])
    if !all(x -> x ≈ dt, diff(t))
        error("time vector t must be uniformly spaced")
    end

    if iscontinuous(sys)
        if method === :zoh
            dsys = c2d(sys, dt, :zoh)
        elseif method === :foh
            dsys, x0map = c2d_x0map(sys, dt, :foh)
            x0 = x0map * [x0; u[:, 1]]
        else
            error("Unsupported discretization method: $method")
        end
    else
        if !(isapprox(sys.Ts, dt))
            error("Time vector must match sample time of discrete-time system")
        end
        dsys = sys
    end

    x = CustomLtitr(dsys.A, dsys.B, u, x0)
    #y = sys.C * x + sys.D * u
    #return SimResult(y, t, x, u, dsys) # saves the system that actually produced the simulation
    return x
end

function CustomNonlinearsim(AAA,B,u,tspan,x0)
    (rows, columns) = size(AAA)
    AA = Matrix{Any}(undef, (rows, columns)) 
    for row in 1:rows
        for column in 1:columns
            h = AAA[row, column]
            if isa(h, Number)
                AA[row, column] = x -> h
            else
                AA[row, column] = h
            end
        end
    end

    A(x) = (|>).(x, AA)

    function f(dx, x, p, t)
        dx .= A(x) * x + B*p
    end

    prob = ODEProblem(f, x0, tspan, u)
    alg = Tsit5()
    sol = solve(prob, alg, reltol=1e-8, abstol=1e-8)
    xout_d = sol.u[2]
    help = Array{Float64}(undef,0)

    for i in xout_d
        append!(help,i)
    end

    return help
end

@views function CustomLtitr(A::AbstractMatrix, B::AbstractMatrix, u::AbstractVecOrMat,
    x0::AbstractVecOrMat=zeros(eltype(A), size(A, 1)))

    T = promote_type(LinearAlgebra.promote_op(LinearAlgebra.matprod, eltype(A), eltype(x0)),
        LinearAlgebra.promote_op(LinearAlgebra.matprod, eltype(B), eltype(u)))

    n = size(u, 2)

    # Using similar instead of Matrix{T} to allow for CuArrays to be used.
    # This approach is problematic if x0 is sparse for example, but was considered
    # to be good enough for now
    x = similar(x0, T, (length(x0), n))

    x[:, 1] .= x0

    pythoncompare = false

    if pythoncompare
        print("u")
    else
        mul!(x[:, 2:end], B, u[:, 1:end-1]) # Do all multiplications B*u[:,k] to save view allocations
    end

    for k = 1:n-1
        if pythoncompare
            mul!(x[:, k+1], B, u[:, k])
            mul!(x[:, k+1], A, x[:, k], true, true)
        else
            mul!(x[:, k+1], A, x[:, k], true, true)
        end
    end
    return x
end