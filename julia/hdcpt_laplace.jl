using Distributed, DataFrames
include("hd_utils.jl")

# Parameters
p_seq = 5
kappa_seq = 1
iterations = 1
n=3200
cpt=800
alpha=0.2
epsilon, alpha, C_gamma = 0.02, 0.2, 0.1

# Grid Setup
heatmap_data = DataFrame(Iterators.product(p_seq, kappa_seq))
rename!(heatmap_data, [:p, :kappa])
heatmap_data.detect_time = zeros(nrow(heatmap_data))

for idx in 1:nrow(heatmap_data)
    curr_p = heatmap_data.p[idx]
    curr_k = heatmap_data.kappa[idx]
    for i in 1:iterations
        Y=rlaplace_hd_cpt(n, p=curr_p, epsilon=epsilon; cpt=cpt, mu_norm=curr_k)
        for t in 2:n
            current_Y=Y[1:t,:]
            delta_t=4*alpha/(t^2*(t+1))
            K=sqrt(log(1/delta_t)/2)
            threshold_s=((K+sqrt(K^2+4(0.09-2*epsilon)))/(2*(0.09-2*epsilon)))^2
            limit = floor(Int, t/2)
            if threshold_s<=limit
                X_t = zeros(limit, current_p)
                for s in 1:limit
                    X_t[s, :] .= (current_Y[t - s + 1, :] .- current_Y[s, :]) ./ (sqrt(2))
                end
            end
            for s in threshold_s:limit
                if robust_mean_test(X_t[1:s,:], 1, delta, epsilon, C_gamma=C_gamma)
                    heatmap_data.detect_time[i]=t
                    break
                end
            end
        end
        heatmap_data.detect_time[i]=n+1
    end
end

#CSV.write("laplace_cpt.csv", heatmap_data)