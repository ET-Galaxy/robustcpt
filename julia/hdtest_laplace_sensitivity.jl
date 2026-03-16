using DataFrames
include("hd_utils.jl")

# Parameters
n_seq = 500
p_seq = 600
C_gamma_seq = 0.1
threshold_seq= 1
kappa_seq = 1
iterations = 100
epsilon, delta = 0.04, 0.1

heatmap_data = DataFrame(
    (n=n, p=p, t=t, C=C) 
    for n in n_seq, p in p_seq, t in threshold_seq, C in C_gamma_seq
)

rename!(heatmap_data, [:n, :p, :kappa0, :C_gamma])
heatmap_data.error_rate = zeros(nrow(heatmap_data))
heatmap_data.power = zeros(nrow(heatmap_data))

for idx in 1:nrow(heatmap_data)
    curr_n = heatmap_data.n[idx]
    curr_p = heatmap_data.p[idx]
    curr_k0 = heatmap_data.kappa0[idx]
    curr_Cgam = heatmap_data.C_gamma[idx]
    
    rejections_type1 = 0
    rejections_power = 0

    for i in 1:iterations
        n_corrupt = rand(Binomial(curr_n, epsilon))
        n_clean = curr_n - n_corrupt
        corrupt = randn(n_corrupt, curr_p) .- 1.0

        # --- Type I Error (Null: mu = 0) ---
        Y_t1 = [rlaplace_hd(n_clean, curr_p, 1/sqrt(2), mu=zeros(curr_p)); corrupt]
        if robust_mean_test(Y_t1, curr_k0, delta, epsilon, C_gamma=curr_Cgam)
            rejections_type1+= 1
        end

        # --- Power (Alternative: mu = kappa0/sqrt(p)) ---
        n_corrupt = rand(Binomial(curr_n, epsilon))
        n_clean = curr_n - n_corrupt
        corrupt = randn(n_corrupt, curr_p) .- 1.0
        mu_vec = fill(curr_k0 / sqrt(curr_p), curr_p)
        Y_p = [rlaplace_hd(n_clean, curr_p, 1/sqrt(2), mu=mu_vec); corrupt]
        if robust_mean_test(Y_p, curr_k0, delta, epsilon, C_gamma=curr_Cgam)
            rejections_power+= 1
        end
    end
    heatmap_data.error_rate[idx] = rejections_type1[] / iterations
    heatmap_data.power[idx] = rejections_power[] / iterations
end
println(heatmap_data)
#CSV.write("laplacetest_scrtp.csv", heatmap_data)
println("Simulation complete. Results in laplacetest_scrtp.csv")