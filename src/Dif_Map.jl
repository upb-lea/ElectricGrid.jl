using KrylovKit

function Spectral_Basis(Gs; num_basis = nothing, scaled = false, alpha = 1)

    """
    Spectral decomposition of the similarity matrix, *after normalization*.
    This uses a variant of the diffusion map algorithm, applied directly to the provided similarity matrix, which replaces the "diffusion" step of the diffusion maps.

    Parameters
    ----------
    Gs : Array of size (N, N)
        A precomputed similarity matrix.
        
    num_basis : number, optional
        Specifies the number of components to retain:
            
        - If num_basis is None (default), then N//2 eigenvalues are retained.
        - If num_basis <= 0, all components are retained
        - If num_basis >=1, this specifies directly the number of components
        - If 0 < num_basis < 1, this specifies a threshold below which eigenvalues are discarded.
        
        The default of N//2 is somewhat arbitrary and you should check what value is adequate for your problem. However, the hope is that the exact value does not matter much: preliminary tests on diverse data sets suggest there is a large plateau with acceptable performances (and then, a decrease as num_basis increases). Note that using all eigenvalues is generally a bad idea: the (necessary) eps parameter used in embed_states introduces a minimal and spurious eigenvalue. Due to the MMD test, the distance between states is also inacurate to some extent (depending on N). And the shift operator definition implicitly uses a matrix inverse, which is badly conditionned when spurious and meaningless eigenvalues are retained. Also, the kernel scale used for building Gx and Gy matters, since it sets the similarity between points, states, hence the sensitivity to small details or to large fluctuations, hence ultimately the eigenspectrum of Gs. One possiblity is to use cross-validation, for example using the accuracy of reconstructing the last nfuture values of the series. 
        
    alpha : float, optional
        Normalization exponent for the diffusion-map like algorithm. 
        
        The default of 1 allows the inference of the geometry, with minimal interference from the density.
        
    overwrite : Bool, optional
        Whether to use the memory in Gs for the computations. This has no effect when using Dask. False by default.

    Returns
    -------
    eigenvalues : array of size M
        The eigenvalues of the spectral decomposition. M is the number of components retained, see the num_basis parameter. The first eigenvalue is always 1 due to the normalization. All other eigenvalues should be greater than 0, unless Gs is badly built (in that case, check the data scale parameter in series_Gxy).

    basis : array of size (N, M)
        The M basis vectors, each of size N, represented as columns. These are also the left eigenvectors of the spectral decomposition. The first vector corresponds to a density and is normalized to sum to 1. The others should all sum to 0. The basis vectors all have the same norm.
        
    coords : array of size (N, M)
        The N coordinates, one for each entry in Gs. Each coordinate is a vector or of dimension M, matching the M basis vectors. The first coordinate is always the constant 1 and can be discarded for purposes such as distance comparisons (knn queries), manifold reconstruction, visualization, etc. The coordinates are the right eigenvectors of the spectral decomposition, scaled to be bi-orthogonal with the basis.

    Notes
    -----

    The pre-normalization step makes this *not* equivalent to a classical eigen or singular decomposition. Distances using the resulting coordinates would correspond to "diffusion distances", assuming Gs is a "diffusion matrix" built like the "diffusion map" algorithm. Although, since we use a completely different building  scheme, this is not a diffusion at all.

    The choice is to provide coordinates scaled with a constant of 1 in the dimension of the eigenvalue 1. This provides consistent values for coordinates along each dimension, with a range close to 1. This also greatly simplifies the computations of the shift operator and avoids division by small values, while being mathematically equivalent in the prediction results. One should however remember that the coordinates should be scaled by the eigenvalues for their relative influence in distance-based algorithms, such as finding nearest neighbors. The scaled distance has some interpretation as explained in the "diffusion map litterature".

    The use of a precision for selecting the number of eigenbasis components is not implemented in Dask mode, since all Dask does is create an operation graph and this would trigger a computation.

    """

    N = size(Gs, 1)

    eigen_cutoff = -1

    if num_basis === nothing

        num_basis = N ÷ 2

    elseif num_basis <= 0

        num_basis = N

    elseif num_basis >= 1

        num_basis = convert(Int64, num_basis)

    else
        eigen_cutoff = num_basis
        num_basis = N
    end

    alpha = clamp(alpha, 0., 1.)

    mat = copy(Gs)

    if alpha > 0

        q = vec(sum(mat, dims = 2))

        if alpha != 1

            q = q.^alpha
        end

        q = (1) ./ q

        mat = mat .* reshape(q, :, 1)
        mat = mat .* reshape(q, 1, :)

    end

    eigval, eigvec, info = geneigsolve((mat, diagm(q)), num_basis, :LM; krylovdim = N, issymmetric = true)
    eigvec = reduce(hcat, eigvec)
    
    if eigen_cutoff >= 0

        num_basis = sum(eigval .>= eigen_cutoff)
    end
    
    if num_basis < size(eigval, 1)
        # copies avoid to keep the full matrix referenced
        eigval = copy(eigval[1:num_basis])
        eigvec = copy(eigvec[:, 1:num_basis])
    end
    
    # Normalization so that eigvec_r[:,1] entries are all 1
    # and that eigvec_l[:,1] matches a density

    eigvec_r = eigvec ./ eigvec[:, 1]
    eigvec_l = eigvec .* (q * eigvec[1,1])

    if scaled
        return eigval, eigvec_l, eigvec_r .* transpose(eigval), info
    end

    return eigval, eigvec_l, eigvec_r, info
end
    
