"""
    rume(X; epsilon=0.0, delta=0.05)

Robust Univariate Mean Estimator (RUME), as described by Prasad et al. (2021).

# Arguments
- `X::Vector{Float64}`: Univariate data.
- `epsilon::Float64`: Contamination level (0 < ε < 1).
- `delta::Float64`: Confidence level (0 < δ < 1).

# Returns
Robust estimate of the mean.
"""
function rume(X::Vector{Float64}; epsilon::Float64=0.0, delta::Float64=0.05)
    n_total = length(X)
    if isodd(n_total)
        error("X must be of even length.")
    end
    n = n_total ÷ 2

    # Split sample
    X1 = sort(X[1:n])
    X2 = X[(n+1):end]

    # Compute parameters
    vareps = max(epsilon, log(1 / delta) / n)
    threshold = 2 * vareps + 2 * sqrt(vareps * log(1 / delta) / n) + log(1 / delta) / n

    if threshold >= 0.5
        error("Not enough sample size to achieve required confidence level.")
    end

    datapts = floor(Int, n * (1 - threshold))

    # Find the shortest interval in X1
    candidates = [X1[i + datapts] - X1[i] for i in 1:(n - datapts)]
    min_index = argmin(candidates)
    interval = (X1[min_index], X1[min_index + datapts])

    # Filter X2 within interval
    filtered_X2 = [x for x in X2 if interval[1] ≤ x ≤ interval[2]]

    return mean(filtered_X2)
end

# ============================================
# t-distribution contaminated with Normal samples
# ============================================

"""
    contaminated_sample_t(n; df=3, epsilon=0.0, mu=0.0)

Generate samples from a t-distribution contaminated with a Normal(0,100^2)
distribution with probability `epsilon`.

# Arguments
- `n::Int`: Number of samples to generate.
- `df::Float64`: Degrees of freedom of the t-distribution.
- `epsilon::Float64`: Contamination proportion (0 ≤ epsilon ≤ 1).
- `mu::Float64`: Mean shift for inlier distribution.
"""
function contaminated_sample_t(n::Int; df::Float64=3, epsilon::Float64=0.0, mu::Float64=0.0)
    final_sample = rand(TDist(df), n) .+ mu
    if epsilon > 0
        contam_mask = rand(n) .< epsilon
        final_sample[contam_mask] .= rand(Normal(0, 100), sum(contam_mask))
    end
    return final_sample
end


# ============================================
# Generate samples from a change point model
# ============================================

"""
    change_point_model(n; mechanism=contaminated_sample_t, cpt=nothing, kappa=1.0)

Generate samples from a change point model with user-specified contamination
level and underlying data generating mechanism.
"""
function change_point_model(n::Int; mechanism=contaminated_sample_t, cpt::Union{Int,Nothing}=nothing, kappa::Float64=1.0)
    if isnothing(cpt)
        # No change point
        return mechanism(n)
    else
        # Before change point
        first_segment = mechanism(cpt)
        # After change point, shifted mean
        second_segment = mechanism(n - cpt; mu=kappa)
        return vcat(first_segment, second_segment)
    end
end


# ============================================
# RUMEDIAN changepoint algorithm
# ============================================

"""
    rumedian(online_data; sigma, v=2, epsilon=0.0, alpha=0.1, C1=1.59, C2=2.25)

Runs the RUMEDIAN algorithm on `online_data`.
"""
function rumedian(online_data::Vector{Float64};
                  sigma::Float64,
                  v::Int=2,
                  epsilon::Float64=0.0,
                  alpha::Float64=0.1,
                  C1::Float64=1.59,
                  C2::Float64=2.25)

    n = length(online_data)
    if epsilon > 0.1
        error("epsilon > 0.1 is not supported.")
    end

    for t in 2:n
        delta_t = (8 * alpha) / (3 * t^3 - 3 * t)
        h_t = ceil(Int, 20 * log(1 / delta_t))

        for s in 1:floor(Int, t / 2)
            if (s % 2 == 0) && (s >= h_t)
                # Use RUME
                vareps = max(epsilon, 2 * log(1 / delta_t) / s)
                diff_rume = abs(rume(online_data[(t - s + 1):t], epsilon=epsilon) -
                                rume(online_data[1:s], epsilon=epsilon))
                zeta = 2 * sigma * C1 * max(vareps^(1 - 1 / v), sqrt(2 / s * log(1 / delta_t)))
                if diff_rume > zeta
                    return Dict(:method => "RUME", :subsample => s, :location => t)
                end
            else
                # Use median
                diff_median = abs(median(online_data[(t - s + 1):t]) - median(online_data[1:s]))
                chi = 2 * sigma * C2 * (0.5 * exp(-1) * (delta_t / 2)^(2 / s) - epsilon)^(-1 / v)
                if diff_median > chi
                    return Dict(:method => "median", :subsample => s, :location => t)
                end
            end
        end
    end

    return Dict(:method => "no changepoint", :subsample => nothing, :location => nothing)
end
