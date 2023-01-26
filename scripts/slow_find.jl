
function series_xy_logk_indx(series, npast, nfuture, concat_valid_map)
       
    
    # The job done for each entry in the matrix
    # computes sum of sq diff for past (sx) and future (sy) sequences
    # n is the number of valid (past, future) pairs

    sum_past_factor = -0.0025
    sum_future_factor = -0.0025

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

   @time begin

    Threads.@threads for k in (n+1)÷2:n - 1

            for l in 0:k - 1

                i = k + 1
                j = l + 1

                # from index of (past,future) pairs to index in data series
                î = concat_valid_map[i]
                ĵ = concat_valid_map[j]

                sumx, sumy = sxy_logk(î, ĵ, series, npast, nfuture, 0, 
                2.0, 1.0, 1.0, 
                200.0, -0.0025, -0.0025)

                #sumx, sumy = sxy_logk(î, ĵ, series, npast, nfuture, sum_past_factor, sum_future_factor)

                sx[i,j] = sumx
                sx[j,i] = sumx
                sy[i,j] = sumy
                sy[j,i] = sumy
            end
                
            for l in k:m - 1

                i = m - k + 1
                j = m - l

                # from index of (past,future) pairs to index in data series
                î = concat_valid_map[i]
                ĵ = concat_valid_map[j]

                sumx, sumy = sxy_logk(î, ĵ, series, npast, nfuture, 0, 
                2.0, 1.0, 1.0, 
                200.0, -0.0025, -0.0025)

                #sumx, sumy = sxy_logk(î, ĵ, series, npast, nfuture, sum_past_factor, sum_future_factor)

                sx[i,j] = sumx
                sx[j,i] = sumx
                sy[i,j] = sumy
                sy[j,i] = sumy
            end
        end
    end
    
    return sx, sy
end

#= function sxy_logk(i, j, series, npast, nfuture, sum_past_factor, sum_future_factor)
    
    sumx = 0

    for t in 0:npast - 1

        d = series[i-t] .- series[j-t]
        ds = sum(abs2, d)

        sumx += ds

    end
    
    sumy = 0

    for t in 0:nfuture - 1

        d = series[i+1+t] .- series[j+1+t]
        ds = sum(abs2, d)

        sumy += ds
    end
    
    return sumx * sum_past_factor, sumy * sum_future_factor
end =#

function sxy_logk(i, j, series, npast, nfuture, localdiff, 
    kernel_params_2, factor_r_past, factor_r_future, 
    sum_r_past, sum_past_factor, sum_future_factor)

    if localdiff == 1

        # weighted average over each past series
        # diff of these => weighted avg of diffs
        r = 1
        delta = zeros(size(series, 2))

        for t in 0:npast-1

            d = series[i-t, :] .- series[j-t, :]
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

        if localdiff != 0
            d = d .- delta
        end
        
        ds = sum(abs2, d)

        if kernel_params_2 != 2
            ds = ds^(0.5*kernel_params_2)
        end

        sumx += ds * r
        r *= factor_r_past

    end
    
    #= r = 1
    sumy = 0

    for t in 0:nfuture - 1

        d = series[i+1+t,:] .- series[j+1+t,:]

        if localdiff != 0
            d = d .- delta
        end

        ds = sum(abs2, d)

        if kernel_params_2 != 2
            ds = ds^(0.5*kernel_params_2)
        end

        sumy += ds * r
        r *= factor_r_future
    end =#

    #= if localdiff == 1

        # weighted average over each past series
        # diff of these => weighted avg of diffs
        r = 1
        delta = zeros(size(series, 2))

        for t in 0:npast-1

            d = series[i-t, :] .- series[j-t, :]

            delta += d * r
            r *= factor_r_past
        end

        delta /= sum_r_past 

    elseif localdiff == 2
        # value of the "present"
        delta = series[i,:] .- series[j,:]
    end =#
    
    #= r = 1
    sumx = 0

    for t in 0:npast - 1

        d = series[i-t] .- series[j-t]

        if localdiff != 0
            d = d .- delta
        end

        ds = sum(abs2, d)

        if kernel_params_2 != 2
            ds = ds^(0.5*kernel_params_2)
        end

        sumx += ds * r
        r *= factor_r_past
    end =#
    
    r = 1
    sumy = 0

    for t in 0:nfuture - 1

        d = series[i+1+t] .- series[j+1+t]

        if localdiff != 0
            d = d .- delta
        end

        ds = sum(abs2, d)

        if kernel_params_2 != 2
            ds = ds^(0.5*kernel_params_2)
        end

        sumy += ds * r
        r *= factor_r_future
    end
    
    return sumx * sum_past_factor, sumy * sum_future_factor
end

N = 1201
series = rand(N)
npast = 200
nfuture = 200

index_map = collect(npast:1:(N - nfuture))

series = reshape(series, :, 1)

sx, sy = series_xy_logk_indx(series, npast, nfuture, index_map);
