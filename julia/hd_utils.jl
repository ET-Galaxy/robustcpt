using Random, LinearAlgebra, Statistics, Distributions

function spectral_filter(Y, gamma2)
    n, p = size(Y)
    w = ones(n)
    
    while true
        # Efficiently compute M = Y' * Diag(w) * Y
        # w .* Y performs row-wise scaling
        M = Y' * (w .* Y)
        A = M - n * I(p)

        # Compute largest singular value (equivalent to svd(A)$d[1])
        # Using svdl on a symmetric matrix; opnorm is more direct
        op_norm = opnorm(A, 2)

        if op_norm < 5 * gamma2
            break
        end

        # Get the first left singular vector
        # eigvals/eigvecs is faster for symmetric A, but svd matches R logic
        sv = svd(A)
        v = sv.U[:, 1]

        # tau_i = <v, Y_i>^2 * 1[w_i > 0]
        projections = Y * v
        tau = (projections.^2) .* (w .> 0)

        tau1 = maximum(tau)
        if tau1 == 0
            break
        end

        # Filter step
        w .= (1.0 .- tau ./ tau1) .* w
    end
    return w
end

"""
    rowsum_filter(Y, w, u)

Subroutine to remove top contamination indices based on inner products.
"""
function rowsum_filter(Y, w, u)
    n, p = size(Y)
    sqrt_w = sqrt.(w)

    # S_vec = colSums(Y * sqrt_w)
    S_vec = vec(sum(Y .* sqrt_w, dims=1))

    # tau_i definition
    inner_terms = (Y * S_vec) .* sqrt_w
    tau = abs.(inner_terms .- w .* p) .* (w .> 0)

    # Remove top u*n indices
    num_remove = floor(Int, u * n)
    if num_remove > 0
        # Get indices of largest tau values
        ord = sortperm(tau, rev=true)
        w_new = copy(w)
        w_new[ord[1:num_remove]] .= 0.0
        return w_new
    end
    
    return w
end

"""
    robust_mean_test(Y, kappa0, delta, epsilon; C_gamma=1.0)

Robust Mean Testing routine (Canconne et al. 2023).
"""
function robust_mean_test(Y, kappa0, delta, epsilon; C_gamma=1.0)
    n, p = size(Y)
    
    if n < p
        throw(ArgumentError("Require n >= p"))
    end

    # Contamination fraction u
    u = epsilon + 1/n + sqrt(log(1/delta) / (2*n))
    
    if u > 0.1
        # Following your R logic: strict thresholding
        throw(ArgumentError("n is too small. u ($u) needs to be less than 0.1"))
    end

    # Compute gamma_2
    gamma2 = C_gamma * (
        u * n * p * log(1/u) +
        sqrt(n * p) * log(2 * p / delta) +
        kappa0 * n
    )

    w = spectral_filter(Y, gamma2)
    w_prime = rowsum_filter(Y, w, u)

    sqrt_w = sqrt.(max.(w_prime, 0.0))
    Sum_wS = vec(sum(sqrt_w .* Y, dims=1))

    test_stat = abs(sum(Sum_wS.^2) - p * sum(w_prime))

    return test_stat >= 0.1 * kappa0^2 * n^2
end

# --- Data Generation Functions ---

function rt_hd(n, p=1, df=3.0; sd=nothing, mu=zeros(p))
    # Default scale logic from R code
    actual_sd = isnothing(sd) ? sqrt(df/(df-2)) : sd
    
    # Generate t-dist samples and scale them
    # Note: Distributions.jl TDist is for standard t
    d = TDist(df)
    # Correcting scale to match R logic (dividing by theoretical SD then multiplying by target)
    samples = rand(d, n, p) .* (actual_sd / sqrt(df/(df-2)))
    
    return samples .+ mu'
end

function rlaplace_hd(n, p=1, s=1.0; mu=zeros(p))
    d = Laplace(0, s)
    samples = rand(d, n, p)
    return samples .+ mu'
end