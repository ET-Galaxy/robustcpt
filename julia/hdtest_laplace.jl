using DataFrames
include("hd_utils.jl")

# Parameters
n_seq = 500
kappa_seq = 1
threshold= 1
iterations = 100
p, epsilon, delta, C_gamma = 700, 0.04, 0.1, 0.9 # smaller C_gamma needed for smaller threshold

# Grid Setup
heatmap_data = DataFrame(Iterators.product(n_seq, kappa_seq))
rename!(heatmap_data, [:n, :kappa0])
heatmap_data.error_rate = zeros(nrow(heatmap_data))
heatmap_data.error_rate2 = zeros(nrow(heatmap_data))

for idx in 1:nrow(heatmap_data)
    curr_n = heatmap_data.n[idx]
    curr_k0 = heatmap_data.kappa0[idx]
    
    rejections_type1 = 0
    rejections_power = 0

    for i in 1:iterations
        # Draw shared n_corrupt for this iteration
        n_corrupt = rand(Binomial(curr_n, epsilon))
        n_clean = curr_n - n_corrupt
        
        # Common outlier matrix
        corrupt = randn(n_corrupt, p) .- 1.0

        # --- Type I Error (Null: mu = 0) ---
        Y_t1 = [rlaplace_hd(n_clean, p, 1/sqrt(2), mu=zeros(p)); corrupt]
        if robust_mean_test(Y_t1, threshold, delta, epsilon, C_gamma=C_gamma)
            rejections_type1+= 1
        end

        # --- Power (Alternative: mu = kappa0/sqrt(p)) ---
        mu_vec = fill(curr_k0 / sqrt(p), p)
        Y_p = [rlaplace_hd(n_clean, p, 1/sqrt(2), mu=mu_vec); corrupt]
        if robust_mean_test(Y_p, threshold, delta, epsilon, C_gamma=C_gamma)
            rejections_power+= 1
        end
    end
    heatmap_data.error_rate[idx] = rejections_type1[] / iterations
    heatmap_data.error_rate2[idx] = 1 - (rejections_power[] / iterations)
end
println(heatmap_data)
#CSV.write("laplacetest_scrtp.csv", heatmap_data)
println("Simulation complete. Results in laplacetest_scrtp.csv")