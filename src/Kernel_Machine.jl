#_______________________________________________________________________________
# Kernel Machines

using LinearAlgebra

function series_Gxy(series, scale, npast, nfuture; decay = 1, qcflags = nothing, localdiff = 0, take2skipX = 0)

    #=
        Compute the past and the future Gram matrices, using the kernel set by set_kernel (defaults to Gaussian)
        
        Arguments:
        series: array or list. Each array is a time-contiguous block of data, 
        with as many rows as time samples and each column is a feature. NaN values are allowed and detected. 
        A list may be provided, in which case each list entry is a different data source. 
        Heterogenous data are supported, so each list entry can be of completely different type. 
        However, each list entry must have the same number of values.
        Each first-level list entry can itself be a second-level list instead of a unique array. 
        In that case, these second-level lists gather data acquired in different time-continuous blocks. 
        These blocks must be consistent across heterogenous data sources: at each time, 
        each data source must be acquired, even though some sensors may fail (NaN values).
        So, in order to enter a unique data source, but with multiple time blocks, 
        the syntax is to use a double list : [[data1, data2, ...]] . 
        
        There, the first-level list specifies a unique data source, and the second level list specifies the contiguous data blocks.
        Multiple sources and multiple blocks read: [[src1_block1, src1_block2], [src2_block1, src2_block2]. 
        Same-numbered blocks from different sources must be of the same length, but different blocks may have different lengths.

        scale: float, or list of floats. If a single float is provided it applies to all series. 
        A list of floats specifies a different scale for each data source.
         
        npast: take that consecutive number of samples for building an observed past. 
        If a single integer is provided it applies to all series. Different values can be provided 
        for different data sources (=each series). 
        It is assumed that data samples are measured at matching times for each data source, 
            hence the first causal state can only be computed after the-largest-npast values have been observed.
        
        nfuture: take that consecutive number of samples, just after the past samples, 
        for building an observed future. If a single integer is provided it applies to all series. 
        Similarly as for npast, the largest nfuture value sets the time of the last computable causal state.
        
        decay: ratio for the furthermost weight compared to the immediate past/future. 
        Defaults to 1 = no decay. If a single float is provided it applies to all series
        
        qcflags : quality control flags. <= 0 means invalid data. optional, but if present there must be one qcflag per sample
        
        localdiff: how to compute the differences between the sequences, pasts and futures. 
        The default is 0 = just take the difference. This assumes that what matters are the 
        series absolute values. localdiff = 1 removes the average (weigthed using decay) of 
        each sequence before taking their differences. This may be useful for eliminating a 
        slow drift due to the sensor itself, which should be removed from data. localdiff = 2 
        means that the "present" value is used as reference. This assumes that what matters 
        are relative variations around the "present" value. Warning: this gives equal importance 
        to all such variations. For example, a fluctuation of -3°C in summer (from ~25°C) is not 
        the same as a fluctuation of -3°C in winter (which may very well cross the freezing point).
        
        take2skipX : int, optional (default 0)
        performs a special kind of subsampling: two consecutive samples are retained, X samples are discarded, 
        then this (X+2) sequence repeats. This scheme is designed to preserve consecutive entries for building a shift operator while still allowing subsampling of very large series. The classic subsampling scheme (take one out of X) can also be applied a priori to the original series, but it is equivalent to a (bad) low-pass filtering. Then, the shift operator would be computed on consecutive entries of the subsampled series, hence at a different time scale. The take2skipX allows to still work on the original time scale. Both can be combined (after appropriate low-pass filtering).
        
    Returns:
        dpast, dfuture: two square matrices of size N, between each observed sequences of a past/future sample pairs.
        If lists are passed as arguments, each list entry is considered its own
        series object, with a matching scale. The data are then combined into
        a unique Gram matrix

        idxmap : N array of int
            Returns the indices of the valid (not NaN) and contiguous (past,future) pairs. For each returned row/col index of the dpast, dfuture matrices, the idxmap specifies what data in the original series that (past,future) reference pair refers to.
    
    =#

    series_list = series

    nseries = size(series_list, 1) # the number of sources

    if length(scale) != 1 scales_list = scale else scales_list = scale*ones(nseries) end

    if length(npast) != 1 npasts_list = npast else npasts_list = Int.(npast*ones(nseries)) end

    if length(nfuture) != 1 nfutures_list = nfuture else nfutures_list = Int.(nfuture*ones(nseries)) end

    if length(decay) != 1 decays_list = decay else decays_list = decay*ones(nseries) end

    if length(localdiff) != 1 localdiff_list = localdiff else localdiff_list = localdiff*ones(nseries) end

    if length(localdiff) != 1 localdiff_list = localdiff else localdiff_list = localdiff*ones(nseries) end

    if qcflags === nothing

        qcflags_list = [ nothing for ser in series_list ]
    else

        if isa(qcflags, Dict) || isa(qcflags, Tuple) qcflags_list = qcflags else qcflags_list = [ qcflags ] end
    end

    total_lx, total_ly = nothing, nothing

    index_map = compute_index_map_multiple_sources(series_list, npasts_list, nfutures_list, qcflags_list, take2skipX = take2skipX)

    for (ser, sca, npa, nfu, dec, ldiff) in zip(series_list, scales_list, npasts_list, nfutures_list, decays_list, localdiff_list)
        
        # ser may itself be a list, this is handled by series_xy_logk

        lx, ly = series_xy_logk(ser, sca, npa, nfu, decay = dec, concat_valid_map = index_map, localdiff = ldiff)
            
        if total_lx === nothing

            total_lx, total_ly = lx, ly
        else
            # This modifies the arguments - only adds the lower-triangular part
            total_lx = parallel_add_lowtri(total_lx, lx)
            total_ly = parallel_add_lowtri(total_ly, ly)
        end
    end
   
    # This modifies the arguments - restores the upper part
    parallel_exp_lowtri(total_lx, nseries)
    parallel_exp_lowtri(total_ly, nseries)

    return total_lx, total_ly, index_map
end

function compute_index_map_multiple_sources(series_list, npasts_list, nfutures_list, qcflags_list; take2skipX = 0)
    #=
        Index mapping helper - for internal use
        
        See parameters of compute_index_map_single_source
        
        This version takes lists of multiple sources as argument.
    =#

    # Each source may have its own NaN patterns and the computed 
    # index maps do not match. => Compute these index maps for all 
    # heterogenous sources and then retain only these indices that
    # are common to all sources.
    # Then, pass that global index map to the _logk function
    
    max_npast = maximum(npasts_list)
    max_nfuture = maximum(nfutures_list)
    
    valid_map = nothing

    for (ser, npa, nfu, qc) in zip(series_list, npasts_list, nfutures_list, qcflags_list)

        skip_start = max_npast - npa
        skip_end = max_nfuture - nfu
        
        sidxmap = compute_index_map_single_source(ser, npa, nfu, skip_start = skip_start, skip_end = skip_end, qcflags = qc, take2skipX = take2skipX)
        
        # combine global indices
        # retain only indices that are valid for all sources
        # A NaN in one source prevents the kernel combination
        # TODO: another possibility is to ignore the source, but how?
        # replace by mean? divide by nsource-1 ?

        if valid_map === nothing valid_map = sidxmap else valid_map = intersect(valid_map, sidxmap) end
    end

    return valid_map
end

function compute_index_map_single_source(series, npast, nfuture; skip_start = 0, skip_end = 0, qcflags = nothing, take2skipX = 0)

    #=
        Index mapping helper - for internal use

        Parameters
        ----------
        series : list of arrays
            list of contiguous data blocks
            
        npast : int
            take that consecutive number of samples for building an observed past. This must be an integer, 
                at least 1 value is needed (the present is part of the past)
            
        nfuture : int
            take that consecutive number of samples, just after the past samples, for building an observed future. 
            nfuture can be 0, in which case, the function only looks for valid pasts. combined with npast==1, the function 
            only looks for valid non-NaN values across all data sources
            
        skip_start : int, optional (default 1)
            number of data to mark as invalid at the beginning of each contiguous data block. 
            Useful to align heterogenous data with different npast/nfuture.
            
        skip_end : int, optional (default 1)
            number of data to mark as invalid at the end of each contiguous data block. 
            Useful to align heterogenous data with different npast/nfuture.
            
        qcflags : list of arrays, optional
            quality control flags. <= 0 means invalid data. optional, but if present there must be one qcflag per sample
            
        take2skipX : int, optional (default 0)
            performs a special kind of subsampling: two consecutive samples are retained, X samples are discarded, 
            then this (X+2) sequence repeats. This scheme is designed to preserve consecutive entries for building 
            a shift operator while still allowing subsampling of very large series. The classic subsampling scheme 
            (take one out of X) can also be applied a priori to the original series, but it is equivalent to a (bad) 
            low-pass filtering. Then, the shift operator would be computed on consecutive entries of the subsampled 
            series, hence at a different time scale. The take2skipX allows to still work on the original time scale. 
            Both can be combined (after appropriate low-pass filtering).

        Returns
        -------
        concat_idx_map: (N,) array of int
            Returns the indices of the valid entries, but stacking all series passed as argument in a single series. 
            The indices refer to valid entries in that global stacked array. Hence, discontiguous data blocks generate 
            invalid entries in the global stacked array.
    =#
    
    if npast < 1

        println("npast must be a strictly positive integer")
    end    

    # nfuture can be 0
    if nfuture < 0

        println("nfuture must be a strictly positive integer")
    end

    if !isa(series, Dict) && !isa(series, Tuple)
        
        #println("series must be a list of arrays")
        series = [series]
    end

    if qcflags === nothing

        qcflags_list = [ nothing for ser in series ]

    else 

        if isa(qcflags, Dict) || isa(qcflags, Tuple) 

            qcflags_list = qcflags 
        
        else

            qcflags_list = [ qcflags ]  
        end
    end
        
    concat_idx_map = nothing
    
    nbefore = 0

    for (sidx, s) in enumerate(series)

        # s is an array
        if ndims(s) == 1

            s = reshape(s, :, 1) # turning the vector into a matrix
        end

        n = size(s, 1) # number of rows

        if n < npast + nfuture

            valid_pf = Array{Float64, 2}(undef, 1, n)
            valid_pf = fill!(valid_pf, false)

        else

            # valid instantaneous time points
            valid_t = .!isnan.(s[:,sidx])

            if qcflags !== nothing && qcflags_list[sidx] <= nothing

                valid_t = fill!(valid_t, false) # if quality flag is raised the entire sample is invalidated
            end

            valid_t[1:skip_start] = fill!(valid_t[1:skip_start], false)
            valid_t[n-skip_end+1:end] = fill!(valid_t[n-skip_end+1:end], false)

            # valid past/future combinations
            valid_pf = copy(valid_t)

            for i in 1:nfuture

                valid_pf[1:n-i] = valid_pf[1:n-i] .& valid_t[i+1:n]
            end

            for i in 1:npast-1

                valid_pf[i+1:n] = valid_pf[i+1:n] .& valid_t[1:n-i]
            end
                
            valid_pf[1:npast-1] .= false
            valid_pf[n-nfuture+1:end] .= false

            if take2skipX > 0

                i = npast

                while i < n

                    i += 2
                    if i>=n

                        break;
                    end

                    for k in range(take2skipX)

                        valid_pf[i] = false
                        i += 1

                        if i >= n
                            break;
                        end
                    end
                end
            end

        end
        
        valid_idx = (1:n)[valid_pf]
            
        if concat_idx_map === nothing

            concat_idx_map = valid_idx .+ nbefore
        else

            concat_idx_map = cat(concat_idx_map, valid_idx .+ nbefore, dims = 1)
        end
            
        nbefore += n
    end
      
    return concat_idx_map
end

function set_kernel(;kernel_type = "Gaussian", kernel_params_ = [-0.5; 2])

    #=
    Sets a kernel type, amongst supported kernels
    kernel_type:
        "Gaussian" string: classical Gausian kernel exp(-0.5 * dist**2 / scale**2)
        "Laplacian" string: classical Laplacian kernel exp(- dist / scale )
        callable: should return the log of the kernel. This is applied to callable(dist**2, scale) in the logk function below
    =#

    if typeof(kernel_type) == String

        if kernel_type == "Gaussian"

            kernel_params_[1] = -0.5
            kernel_params_[2] = 2

            return kernel_params_

        elseif kernel_type == "Laplacian"

            kernel_params_[1] = -1
            kernel_params_[2] = 1

            return kernel_params_
        else

            println("Invalid kernel type: only Gaussian and Laplacian are supported, or you should provide your own kernel function taking a squared distance as argument")
        end

    end

    return nothing
end

function series_xy_logk(series, scale, npast, nfuture; decay = 1, concat_valid_map = nothing, localdiff = 0)
    
    if isa(series, Dict) || isa(series, Tuple)
        # Concatenate all the series and use the global validity index below
        println("This functionality does not exist yet")
    end

    if ndims(series) == 1
        series = reshape(series, :, 1) # turning the vector into a matrix
    end

    #
    kernel_params_ = [-0.5, 2]
    #
    return series_xy_logk_indx(series, scale, npast, nfuture, decay,
     concat_valid_map, kernel_params_[1], kernel_params_[2], localdiff)
end

function series_xy_logk_indx(series, scale, npast, nfuture, decay, concat_valid_map, kernel_params_1, kernel_params_2, localdiff)
    
    if npast <= 1 factor_r_past = 1 else factor_r_past = exp(log(decay)/(npast - 1.)) end
    if nfuture <= 1 factor_r_future = 1 else factor_r_future = exp(log(decay)/(nfuture - 1.)) end
    
    sum_r_past = 1
    r = 1

    for t in npast-2:-1:0

        r *= factor_r_past
        sum_r_past += r
    end

    sum_past_factor = kernel_params_1 / (sum_r_past * scale^kernel_params_2)
    sum_r_future = 1
    r = 1

    for t in npast+1:npast + nfuture - 1

        r *= factor_r_future
        sum_r_future += r
    end
    
    sum_future_factor = kernel_params_1 / (sum_r_future * scale^kernel_params_2)    
    
    # The job done for each entry in the matrix
    # computes sum of sq diff for past (sx) and future (sy) sequences
    # n is the number of valid (past, future) pairs
    n = length(concat_valid_map)
    sx = zeros(n, n)
    sy = zeros(n, n)
    
    # Triangular indexing, folded
    # x
    # y y     =>    z z z y y
    # z z z         w w w w x
    # w w w w
    # outer loops can now be parallelized - all have about the same duration
   if n%2 == 1  m = n else m = n-1 end
   
   Threads.@threads for k in (n+1)÷2:n - 1

        for l in 0:k - 1

            i = k + 1
            j = l + 1

            sumx, sumy = sxy_logk(i, j, series, concat_valid_map, npast, nfuture, 
            localdiff, kernel_params_2, factor_r_past, 
            factor_r_future, sum_r_past, sum_past_factor, sum_future_factor)

            sx[i,j] = sumx
            sx[j,i] = sumx
            sy[i,j] = sumy
            sy[j,i] = sumy
        end
            
        for l in k:m - 1

            i = m - k + 1
            j = m - l

            sumx, sumy = sxy_logk(i, j, series, concat_valid_map, npast, nfuture, 
            localdiff, kernel_params_2, factor_r_past, 
            factor_r_future, sum_r_past, sum_past_factor, sum_future_factor)

            sx[i,j] = sumx
            sx[j,i] = sumx
            sy[i,j] = sumy
            sy[j,i] = sumy
        end
    end
    
    return sx, sy
end

function sxy_logk(i, j, series, concat_valid_map, npast, nfuture, localdiff, kernel_params_2, factor_r_past, factor_r_future, sum_r_past, sum_past_factor, sum_future_factor)
            
    # from index of (past,future) pairs to index in data series
    i = concat_valid_map[i]
    j = concat_valid_map[j]

    delta = zeros(size(series, 2))

    if localdiff == 1

        # weighted average over each past series
        # diff of these => weighted avg of diffs
        r = 1

        for t in 0:npast-1 # TODO: Can this loop be parallelized?

            d = series[i-t, :] - series[j-t, :]
            delta += d * r
            r *= factor_r_past
        end

        delta /= sum_r_past
        

    elseif localdiff == 2
        # value of the "present"
        delta = series[i,:] - series[j,:]
    end
    
    r = 1
    sumx = 0

    for t in 0:npast - 1

        d = series[i-t,:] .- series[j-t,:]
        d = d - delta
        ds = sum(d.*d)

        if kernel_params_2 != 2
            ds = ds^(0.5*kernel_params_2)
        end

        sumx += ds * r
        r *= factor_r_past

    end
    
    r = 1
    sumy = 0

    for t in 0:nfuture - 1

        d = series[i+1+t,:] - series[j+1+t,:]
        d = d - delta
        ds = sum(d.*d)

        if kernel_params_2 != 2
            ds = ds^(0.5*kernel_params_2)
        end

        sumy += ds * r
        r *= factor_r_future
    end
    
    return sumx * sum_past_factor, sumy * sum_future_factor
end

function parallel_add_lowtri(total, mat)
    """
    WARNING: ONLY ADDS THE LOWER PART
    """
    N = size(mat, 1)
    if N%2 == 1 M = N else M = N-1 end
        
    # outer loops can be parallelized - all have about the same duration
    Threads.@threads for k in (N+1) ÷ 2:N - 1

        for l in 0:k - 1

            i = k + 1
            j = l + 1
            total[i,j] += mat[i,j]
        end
            
        for l in k:M - 1

            i = M - k + 1
            j = M - l
            total[i,j] += mat[i,j]
        end
    end

    Threads.@threads for d in 1:N

        total[d,d] += mat[d,d]
    end

    return total
end

function parallel_exp_lowtri(mat, scale)

    N = size(mat, 1)
    if N%2 == 1 M = N else M = N-1 end
    invscale = 1/scale
        
    # outer loops can be parallelized - all have about the same duration
    Threads.@threads for k in (N+1) ÷ 2:N - 1

        for l in 0:k - 1

            i = k + 1
            j = l + 1
            e = exp(mat[i,j] * invscale)
            mat[i,j] = e
            mat[j,i] = e
        end
            
        for l in k:M-1

            i = M - k + 1
            j = M - l
            e = exp(mat[i,j] * invscale)
            mat[i,j] = e
            mat[j,i] = e
        end
    end

    Threads.@threads for d in 1:N
        mat[d,d] = exp(mat[d,d] * invscale)
    end

    return mat
end

function Embed_States(Gx, Gy; eps = 1e-8, normalize = true, return_embedder = false)

    """
    Compute a similarity matrix for the embedded causal states, seen as distributions P(Y|X), 
    using the conditional mean embedding.

    Arguments:
        Gx: a similarity matrix of pasts X. Gx[i,j] = kernel_X(x_i, x_j) with kernel_X a reproducing kernel for X.
        Gy: a similarity matrix of futures Y.
        eps: amount of regularization. The theory requires a regularizer and shows that the regularized 
        estimator is consistent in the limit of N -> infinity. In practice, using a too small regularizer 
        may cause divergence, too large causes innaccuracies.
        normalize: whether to renormalize the returned similarity matrix, so that each entry along the 
        diagonal is 1. See below.
        return_embedder: whether to return an embedder object for new, unknown data. Default is False
        
    Returns:

        Gs: A similarity matrix between each causal states. Entries Gs[i,j] can be seen as inner products 
        between states S_i and S_j. 
        With normalized kernels, such as used in series_Gxy, then state vectors should also be normalized. 
        Also, theoretically, that matrix should be positive definite. In practice, numerical innacuracies 
        and estimating from finite samples may destroy both previous properties. Renormalization is performed 
        by default, but if positive definiteness issues happens, try using a larger regularizer.

        embedder: (if return_embedder is True) A matrix for embedding new Kx kernel similarity vectors
    """
    
    #Omega = (Gx + I*eps) \ Gx

    Omega = copy(Gx)
    ldiv!(factorize(Gx + I*eps), Omega)

    Omega = Symmetric(Omega) # should not be needed
    
    embedder = transpose(Omega) * Gy
    
    Gs = embedder * Omega
    
    if normalize

        #= This normalization is especially useful for avoiding numerical
        errors later in the Gs eigendecomposition. But Gs should already
        have nearly 1s on the diagonal to start with. For embedding new
        vectors, this refinement is probably overkill
        TODO: tests if this is really necessary for the embedder,
        knowing the distributions will be normalized later on.
        Tests: actually, this impacts a bit the numerical resolution
        when expressing the embedded Ks as a combination of Gs lines
        option 1 : all the weights on the left side:
        embedder /= Gs.diagonal()[:,None]
        if applied to Gs, this would yield 1s on the diagonal
        but some asymetry
        Tests: this options yields some inaccuracies
        option 2 : use the same normalization as GS:
        but that normalisation is not applicable on the right
        to new kvecs
        Tests: this options yields some inaccuracies, but less than option 1 =#

        d = (1) ./ (sqrt.(diag(Gs)))
        dr = reshape(d, :, 1)
        embedder = embedder .* dr

        # Now Gs

        Gs = Gs .* dr
        Gs = Gs .* reshape(d, 1, :)
        Gs = 0.5 * (Gs + transpose(Gs))

        # just to be extra safe, avoid floating-point residual errors
        Gs[diagind(Gs)] .= 1.0
        
    end
    
    if return_embedder

        return Gs, embedder
    end

    return Gs
end