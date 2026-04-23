using Distributed, DataFrames, CSV

# Define your grid
n_seq = 400:20:1200
p_seq = 50:25:750
reps = 1000
epsilon =0.01
delta=0.1

@everywhere begin
    using Random, Statistics, LinearAlgebra, Distributions
    include("hd_utils.jl") 

    # Define constants for the workers
    const kappa_ = 0.5
    const epsilon_ = 0.01
    const delta_ = 0.1
end

# 1. Create a single master DataFrame
results_df = DataFrame(Iterators.product(n_seq, p_seq))
rename!(results_df, [:n, :p])
results_df.type1_error = fill(NaN, nrow(results_df))
results_df.type2_error = fill(NaN, nrow(results_df))

for idx in 1:nrow(results_df)
    curr_n = results_df.n[idx]
    curr_p = results_df.p[idx]
    if epsilon + min(curr_p/(20*curr_n),1/20) + sqrt(log(1/delta) / (2*curr_n))>0.1
        results_df.type1_error[idx] = NaN
        results_df.type2_error[idx] = NaN
    else
    # pmap returns a vector of tuples [(t1, pwr), (t1, pwr), ...]
    sim_results = pmap(1:reps) do _
        for C_gamma_ in [0.1,0.06]
            # --- Type I Error ---
            n_corrupt = rand(Binomial(curr_n, epsilon_))
            n_clean = curr_n - n_corrupt
            corrupt = randn(n_corrupt, curr_p) .- 1/sqrt(curr_p)
            
            Y_t1 = [rt_hd(n_clean, curr_p, 4.1, sd=1, mu=zeros(curr_p)); corrupt]
            t1_rej = robust_mean_test(Y_t1, kappa_, delta_, epsilon_, C_gamma=C_gamma_)
 
            # --- Power ---
            n_corrupt = rand(Binomial(curr_n, epsilon_))
            n_clean = curr_n - n_corrupt
            corrupt = randn(n_corrupt, curr_p) .- 1/sqrt(curr_p)

            mu_vec = fill(kappa_ / sqrt(curr_p), curr_p)
            Y_p = [rt_hd(n_clean, curr_p, 4.1, sd=1, mu=mu_vec); corrupt]
            pwr_rej = robust_mean_test(Y_p, kappa_, delta_, epsilon_, C_gamma=C_gamma_)
            if (t1_rej <= 0.1 && pwr_rej <= 0.1) || C_gamma_ == 0.06
                return (t1_rej, pwr_rej)
            end
         end
    end
    # 2. Extract and assign directly to the DataFrame row
    results_df.type1_error[idx] = sum(r[1] for r in sim_results) / reps
    results_df.type2_error[idx] = 1 - (sum(r[2] for r in sim_results) / reps)
    end
end

# 3. Save as a single combined file
CSV.write("samp_comp/v4varyp.csv", results_df)
println("All results saved to v4varyp.csv")

