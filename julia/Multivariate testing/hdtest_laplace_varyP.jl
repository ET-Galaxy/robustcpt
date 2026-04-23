using Distributed, DataFrames, CSV

@everywhere begin
    using Random, Statistics, LinearAlgebra, Distributions
    include("hd_utils.jl") 

    # Define constants for the workers
    const kappa_ = 0.5
    const s_=1/sqrt(2)
    const epsilon_ = 0.04
    const delta_ = 0.1  
end

# Define your grid
n_seq = 260:20:1400
p_seq = 25:25:600
reps = 1000

# 1. Create a single master DataFrame
results_df = DataFrame(Iterators.product(n_seq, p_seq))
rename!(results_df, [:n, :p])
results_df.type1_error = zeros(nrow(results_df))
results_df.type2_error = zeros(nrow(results_df))

for idx in 1:nrow(results_df)
    curr_n = results_df.n[idx]
    curr_p = results_df.p[idx]
    
    # pmap returns a vector of tuples [(t1, pwr), (t1, pwr), ...]
    sim_results = pmap(1:reps) do _
        for C_gamma_ in [0.1,0.05]
            # --- Type I Error ---
            n_corrupt = rand(Binomial(curr_n, epsilon_))
            n_clean = curr_n - n_corrupt
            corrupt = randn(n_corrupt, curr_p) .- 1/sqrt(curr_p)
            
            Y_t1 = [rlaplace_hd(n_clean, curr_p, s_, mu=zeros(curr_p)); corrupt]
            t1_rej = robust_mean_test(Y_t1, kappa_, delta_, epsilon_, C_gamma=C_gamma_)
 
            # --- Power ---
            n_corrupt = rand(Binomial(curr_n, epsilon_))
            n_clean = curr_n - n_corrupt
            corrupt = randn(n_corrupt, curr_p) .- 1/sqrt(curr_p)

            mu_vec = fill(kappa_ / sqrt(curr_p), curr_p)
            Y_p = [rlaplace_hd(n_clean, curr_p, s_, mu=mu_vec); corrupt]
            pwr_rej = robust_mean_test(Y_p, kappa_, delta_, epsilon_, C_gamma=C_gamma_)
            if (t1_rej <= 0.1 && pwr_rej <= 0.1) || C_gamma_ == 0.05
                return (t1_rej, pwr_rej)
            end
        end
    end

    # 2. Extract and assign directly to the DataFrame row
    results_df.type1_error[idx] = sum(r[1] for r in sim_results) / reps
    results_df.type2_error[idx] = 1 - (sum(r[2] for r in sim_results) / reps)
end

# 3. Save as a single combined file
CSV.write("samp_comp/th1varyp_small.csv", results_df)
println("All results saved to th1varyp_small.csv")

