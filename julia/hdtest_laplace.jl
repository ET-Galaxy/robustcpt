using Distributed, DataFrames, CSV

# Parameters
n_seq = 500:100:3000
kappa_seq = 0.25:0.01:0.43
iterations = 1000
p, epsilon, delta, C_gamma = 10, 0.05, 0.1, 0.1

# Grid Setup
heatmap_data = DataFrame(Iterators.product(n_seq, kappa_seq))
rename!(heatmap_data, [:n, :kappa0])
heatmap_data.error_rate = zeros(nrow(heatmap_data))
heatmap_data.error_rate2 = zeros(nrow(heatmap_data))

println("Simulation started. Workers: $(Threads.nthreads())")

for idx in 1:nrow(heatmap_data)
    curr_n = heatmap_data.n[idx]
    curr_k0 = heatmap_data.kappa0[idx]
    
    # Atomic counters for thread-safe increments
    rejections_type1 = Threads.Atomic{Int}(0)
    rejections_power = Threads.Atomic{Int}(0)

    Threads.@threads for i in 1:iterations
        # Draw shared n_corrupt for this iteration
        n_corrupt = rand(Binomial(curr_n, epsilon))
        n_clean = curr_n - n_corrupt
        
        # Common outlier matrix
        corrupt = randn(n_corrupt, p) .- 1.0

        # --- Type I Error (Null: mu = 0) ---
        Y_t1 = [rlaplace_hd(n_clean, p, 1.0, mu=zeros(p)); corrupt]
        if robust_mean_test(Y_t1, curr_k0, delta, epsilon, C_gamma=C_gamma)
            Threads.atomic_add!(rejections_type1, 1)
        end

        # --- Power (Alternative: mu = kappa0/sqrt(p)) ---
        mu_vec = fill(curr_k0 / sqrt(p), p)
        Y_p = [rlaplace_hd(n_clean, p, 1.0, mu=mu_vec); corrupt]
        if robust_mean_test(Y_p, curr_k0, delta, epsilon, C_gamma=C_gamma)
            Threads.atomic_add!(rejections_power, 1)
        end
    end
    
    heatmap_data.error_rate[idx] = rejections_type1[] / iterations
    heatmap_data.error_rate2[idx] = 1 - (rejections_power[] / iterations)
end

CSV.write("laplacetest_scrtp.csv", heatmap_data)
println("Simulation complete. Results in laplacetest_scrtp.csv")