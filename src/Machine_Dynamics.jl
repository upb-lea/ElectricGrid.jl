
function Shift_Operator(coords, eigenvalues; index_map = nothing, return_eigendecomposition = false)
    """
    Compute a shift operator for the dynamical system, assuming coords are the time series of causal states
    Handles NaN values and split / multi series with index maps
    
    Arguments:
        coords : the coordinates of the causal states
        index_map : see series_Gxy
        return_eigendecomposition : whether to return the eigenvalues and eigenvectors of the shift operator. Default is False
        
    Returns:
        The forward operator for state distributions, represented in the eigenbasis.
        if return_eigendecomposition is True, also return the eigenvalues and left, right eigenvectors
    """
    
    # TODO: generator, meaning the log in eigen decomposition form
    # This would allow to compute powers more efficiently

    # Shift operator, expressing the new coordinates as a combination of the old ones. 
    # This definition is similar to that of the kernel Koopman operator in RKHS, 
    # but done here in diffusion map coordinates space
    
    if index_map === nothing

        valid_prev = range(1, size(coords, 1) - 1)
        indices_no_next = [size(coords, 1)]
    else

        # build reverse operators
        valid_prev = deleteat!(collect(range(1, size(coords, 1) - 1)), findall(diff(index_map) .!= 1))
        # also known as gaplocs, plus the latest state
        indices_no_next = setdiff(range(1, size(coords, 1)), valid_prev)
    end

    # To replace in formulas:
    valid_next = valid_prev .+ 1
    
    # - shift_op is a stochastic transition matrix: rows sum to 1
    shift_op = zeros(size(coords, 2), size(coords, 2))
    shift_op[1,1] = 1.0

    # we want shift_op * transpose(coords[valid_prev,:]) = transpose(coords[valid_next,:])
    # and some generalization outside valid_prev
    # below are raw/naive options, they diverge and do not generalize much
    # even when coords are replaced by coords*eigenvalues
    # This is the fastest
    # shift_op[2:end,:] = transpose(coords[valid_next, 2:end]) * pinv(transpose(coords[valid_prev,:]))

    # The idea now is to express the coordinates with no valid next values
    # as a combination of known coordinates, so that the shift operator
    # expresses coordinates only using the known / good examples

    # Using inverses yields the least accurate results
    
    # Using non-negative least squares instead works better, see the paper

    # weights for each entry with no successor

    X = nonneg_lsq(transpose(transpose(eigenvalues).*coords[valid_prev,:]), 
                    transpose(transpose(eigenvalues).*coords[indices_no_next,:]), 
                    alg = :nnls) #fnnls - for fast Non Neg Least Squares - less accurate

    X = X ./ sum(X, dims = 1) # should already be close to 1
    
    T = diagm(-1 => ones(size(coords, 1)-1))
    
    # put that back into T. Process column by column, due to indexing scheme
    # that would make the matrix read-only

    for (ic, c) in enumerate(indices_no_next)

        T[:, c] .= 0
        T[valid_next, c] = X[:,ic]
    end
    
    shift_op[2:end,:] = transpose(coords \ (transpose(T)*coords[:, 2:end]))
    
    # This may be another option - gives similar results
    #= b = (transpose(T)*coords[:, 2:end])
    ldiv!(factorize(coords), b)
    rows_to_keep = b[:,1:end] .!= 0
    b = b[findall(rows_to_keep[:,1]),:]
    shift_op[2:end,:] = transpose(b) =#
    
    # Now ensure the operator power does not blow up
    # Eigenvalues should have modulus less than 1...
    eigval, right_eigvec, info = eigsolve(shift_op)
    right_eigvec = reduce(hcat, right_eigvec)

    # This is to check that the eigendecomposition isn't bullshit
    #= Λ = Diagonal(eigval)
    println(shift_op * right_eigvec ≈ right_eigvec * Λ) # should return true
    println(left_eigvec * shift_op ≈  Λ * left_eigvec) # should return true =#
    
    if maximum(abs.(eigval)) > 1

        # ... but there may be numerical innacuracies, or irrelevant
        # component in the eigenbasis decomposition (the smaller eigenvalues
        # there have no significance, according to the MMD test, and we use
        # implicitly an inverse here). 
        # This should not happen, but if it does, clip eigvals
        n = size(eigval, 1)

        for i in 1:n

            if abs(eigval[i]) > 1

               eigval[i] /= abs(eigval[i])
            end
        end

        # reconstruct best approximation of the shift operator
        # faster formula using solve adapted from

        Λ = Diagonal(eigval)
        shift_op = real.((right_eigvec * Λ) * inv(right_eigvec))

        if return_eigendecomposition

            eigval, right_eigvec, info = eigsolve(shift_op)
            right_eigvec = reduce(hcat, right_eigvec)
            left_eigvec = inv(right_eigvec) # may cause issues, i.e. slow
        end
    end
    
    if return_eigendecomposition
        left_eigvec = inv(right_eigvec) # may cause issues, i.e. slow
        return shift_op, eigval, left_eigvec, right_eigvec
    end
    
    return shift_op
end

function immediate_future(data, indices)
    """
    Equivalent to lambda d,i: d[i+1,:] . This function is used as a default argument in expectation_operator. See the documentation there.
    """
    return data[indices .+ 1, :]
end

function Expectation_Operator(coords, index_map, targets; func::Function = immediate_future)
    """
    Builds the expectation operator, mapping a state distribution expressed in the eigenbasis, into numerical values, expressed in the original series units.

    Parameters
    ----------
        
    coords : array
        coordinates in the eigenbasis, as returned by the spectral_basis function
        
    index_map : int array
        This is the index_map returned by the series_Gxy function. It indicates the index in the 
        series for the (past, future) pair matching each coordinate entry. That index is that of 
        the present, the last value in the past series. If the series consist of several discontiguous 
        blocks, the index refers to the valid entries in the series that would be made by concatenating 
        these blocks (that series would have invalid entries at each time discontinuity, at nan values, etc). 
        The targets parameter may have its own validity pattern, and one that differs for each heterogenous 
        data source - it may be that a computed state, has no matching target for one data source, but a valid 
        target for another. Expectation operators are computed using all valid data, for each heterogenous data source.
        
    targets : array or list of array
        Target values that we wish to predict from causal states, to build an operator from. Each array must have 
        the same temporal structure (number of samples, discontinuous blocks, etc) as the data series that were 
        used to build the causal states - one target is needed for each original series measurement. Some targets
        may be NaN, in which case they will be ignored for building the operator. These NaN patterns may differ 
        from these implied by the causal state index_map. If a list is provided as the targets argument, 
        these are observations from multiple (possibly heterogenous) data source, and possibly different sources 
        from those used for building the causal states. Then, a list of expectation operators is returned, one of 
        each data source.
        
    function : callable, or list of callable, optional
        This is the function of the targets, which expectation is computed by the operator. Different functions 
        can be provided for each heterogenous data source by providing a list of callables. Each function takes 
        as parameter a data array (which is in fact the stacked targets) and an index array (the index_map), and 
        should return the function computed on the data at each specified index. The default is to use the 
        immediate_future function, which is equivalent to lambda d,i: d[i+1,:]. 
        Due to the way the index_map is built, the function is applied only on target data values that match 
        valid causal states. Thus, it is guaranteed that data in the range i-npast+1:i+nfuture is valid: i is 
        the "present", the last entry of the past sequence in the consecutive (past, future) pair. 
        Assuming nfuture>=1, then the next value at i+1 always exists and it is the first entry in the future sequence.


    Returns
    -------
    expect_ops : array or list of arrays
        This is the expectation operator for the function (or functions) that compute a value from the current state. 
        The expectation is implicitly performed by the reproducing property and using mean maps for the probability 
        distributions. If a list of data sources was provided as the series argument, return one search operator for 
        each data source.
        
    Notes
    ----
    
    - As suggested by the 'targets' parameter name, these need not be the same as the 'series' argument of the series_Gxy 
    that were used to build the causal states. You may very well construct the causal states from some data sources, 
    and try to predict only a subset of these sources, or maybe other observables entirely that you assume depend 
    on the causal states. The only restriction is that the series and the targets variable share the same temporal 
    pattern: the same number of samples in the same number of continuous data blocks. Use NaN if some target values 
    are not available for all sample times.
    
    - Theoretically, any function of the causal state is useable, so long as it can be estimated for each observed 
    (past, future) data pair in the target measurements. The function could thus very well take non-numeric data as 
    argument, but, currently, it should only return numeric values (scalar or vectorial).
    
    - TODO: Given the above, it would be possible to extend other machine learning algorithms taking as input the 
    (past,future) pairs and producing some value of interest. The machine learning trained instance could be fed 
    as the function parameter, then it would benefit from the causal state machinery.
    
    """
    
    if index_map === nothing
        println("You must provide a valid index map")
    end
    
    targets_list = targets

    if (isa(func, Dict) || isa(func, Tuple)) 
        f_list = func 
    else 
        f_list = Array{Function, 1}(undef, length(targets_list))
        f_list = fill!(f_list, func)
    end
    
    eoplist = Vector{Matrix{Float64}}(undef, length(targets_list))
    
    rtol = sqrt(eps(real(float(one(eltype(transpose(coords)))))))
    sci = pinv(transpose(coords), rtol = rtol)

    e_cnt = 1
    
    for (tar, fun) in zip(targets_list, f_list)
        
        if isa(tar, Dict) || isa(tar, Tuple)
            println("This functionality does not exist yet")
            # Something like tar = reduce(hcat, tar)
        end
        
        if ndims(tar) == 1
            tar = reshape(tar, :, 1)
        end    

        fvalues = fun(tar, index_map)
        
        valid_f = (findall(vec(prod(.!isnan.(fvalues), dims = 2))))

        eoplist[e_cnt] = transpose(fvalues[valid_f,:]) * sci[valid_f,:]
        e_cnt += 1
    end
    
    if isa(targets, Dict) || isa(targets, Tuple) 

        return eoplist
    
    else

        return eoplist
    end
end

function Predict(npred, state_dist, shift_op, expect_op; return_dist = 0, bounds = nothing, knn_convexity = nothing, coords = coords, knndim = nothing, extent = nothing)
    """
    Predict values from the current causal states distribution

    Parameters
    ----------
    npred : int
        Number of predictions to generate
        
    state_dist : array of size (num_basis, 1)
        Current state distribution, expressed in the eigen basis
        
    shift_op : array of size (num_basis, num_basis)
        Operator to evolve the distributions in time
        
    expect_op : array of size (data_dim, num_basis) or a list of such arrays
        Operator, or list of operators, that take a distribution of causal states and generate a data value compatible with the original series
        
    return_dist : {0, 1, 2}, optional
        Whether to return the state distribution, or update it:
        - 0 (default): do not return a state vector
        - 1: return the updated state vector
        - 2: return an array of state_dist vectors, one row for each prediction

    Returns
    -------
    predictions: array, or list of arrays, matching the expect_op type
        As many series of npred values as there are expectation operators. 
        Either a list, or a single array, matching the expect_op type
    
    updated_dist: optional, array
        If return_dist > 0, the updated state distribution or all such distributions 
        are returned as a second argument.
    
    Notes
    -----
    
    The causal state space (conditional distributions) is not convex: linear combinations of causal states do not 
    necessarily correspond to a valid value in the conditionned domain. We are dealing here with distributions of 
    causal states, but these distributions are ultimately represented as linear combinations of the data span in 
    a Reproducing Kernel Hilbert Space, or, in this case, as combinations of eigenbasis vectors 
    (themselves represented in RKHS). Therefore, the state distributions are also linear combinations of other 
    causal states. Ultimately, the continuous-time model converges to a limit distribution, which is an average 
    distribution that need not correspond to any single state, so the non-convexity is not an issue for that limit 
    distribution. Making predictions with the expectation operator is still feasible, and we get the expected value 
    from the limit distribution as a result. 

    If, instead, one wants a trajectory, and not the limit average, then a method is required to ensure that each 
    predicted state remain valid as a result of applying a linear shift operator. The knn_convexity argument is an 
    attempt to solve this issue. The API is not definitive and may change in the future. The preimage issue is a 
    recurrent problem in machine learning and no single answer can currently solve all cases.
    
    """
    
    all_state_dists = Array{Float64, 2}(undef, length(state_dist), npred)

    pred = Vector{Matrix{Float64}}(undef, length(expect_op))
    num_basis = size(coords, 2)

    problem = nothing

    if isa(knn_convexity, Int) && knn_convexity > 0

        if isnothing(knndim) || knndim > num_basis || !isa(knndim, Int)
            knndim = num_basis
        end
        if isnothing(extent) extent = 0.0 end

        # Euclidean(3.0), Chebyshev, Minkowski(3.5) and Cityblock
        balltree = BallTree(transpose(coords[:, 1:knndim]); leafsize = 30)

        if knn_convexity > 1

            problem = Model(Ipopt.Optimizer)
            set_silent(problem)

            # Optimisation problem
            @variable(problem, -extent <= w[i = 1:knn_convexity] <= 1 + extent) # weights between nearest neighbours
            @NLparameter(problem, M[i = 1:num_basis, j = 1:knn_convexity] == 0) # nearest neighbours

            @NLexpression(problem, sum_w, sum(w[i] for i in 1:knn_convexity))
            @NLconstraint(problem, sum_w == 1) # weights have to sum to 1

            # x is defined as a combination of neighours
            x = Array{NonlinearExpression, 1}(undef, num_basis) # the end result

            for i in 1:num_basis # maybe parallelize?
                x[i] = @NLexpression(problem, sum(M[i, j]*w[j] for j in 1:knn_convexity)) # matrix multiplication
            end

            # Preserve the state distribution normalization
            # This is always 1, in whatever scaled or coords units, since 
            # the first eigenvatlue is 1.
            @NLconstraint(problem, x[1] == 1.0)

            @NLparameter(problem, target_x[i = 1:num_basis] == 0) # nearest neighbours

            @NLexpression(problem, sum_squares, sum((x[i] - target_x[i])^2 for i in 1:num_basis))
        end
    else

        knn_convexity = nothing
    end

    for p in 1:npred

        if return_dist==2
            all_state_dists[:, p] = state_dist
        end

        # Apply the expectation operator to the current distribution
        # in state space, to make a prediction in data space

        for (eidx, eop) in enumerate(expect_op)

            new_pred = (eop * state_dist)

            if ndims(new_pred) == 1
                new_pred = reshape(new_pred, :, 1)
            end   

            if p == 1
                pred[eidx] = cat(new_pred, dims = 1)
            else
                pred[eidx] = cat(pred[eidx], new_pred, dims = 1)
            end
        end

        # Evolve the distribution
        state_dist = shift_op * state_dist
        # Normalize - not needed anymore by construction of the shift op
        state_dist /= state_dist[1, 1]

        if !isnothing(knn_convexity) 

            idxs, _ = knn(balltree, state_dist[1:knndim], knn_convexity)
            #idxs, dista = knn(balltree, coords[2, :], knn_convexity)

            if !isnothing(problem)

                temp_c = Matrix(transpose(coords[idxs, :]))

                # Setting the parameters to their newest values - maybe parallelize?
                for i in 1:num_basis

                    for j in 1:knn_convexity

                        set_value(M[i, j], temp_c[i, j])
                    end

                    set_value(target_x[i], state_dist[i])
                end

                @NLobjective(problem, Min, sum_squares)
                optimize!(problem)

                state_dist = value.(x)
                # should not be needed mathematically if the 
                # constraint was perfectly respected.
                # There could be numerical inaccuracies in practice
                state_dist /= state_dist[1, 1]

            elseif knn_convexity == 1

                state_dist = coords[idxs[1], :]
            end
        end
    end
        
    if return_dist == 1

        return pred, state_dist
    
    elseif return_dist == 2

        return pred, all_state_dists
    end
    
    return pred
end
