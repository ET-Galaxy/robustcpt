using Distributed, DataFrames
include("hd_utils.jl")

# Parameters
p_seq = 10
kappa_seq = 0
kappa_guess = 1
iterations = 10
n=4000
cpt=1500
epsilon, alpha, C_gamma = 0.01, 0.1, 0.1

# Grid Setup
#heatmap_data = DataFrame(Iterators.product(p_seq, kappa_seq))
#rename!(heatmap_data, [:p, :kappa])
#heatmap_data.avg_detect_time = zeros(nrow(heatmap_data))

# for idx in 1:nrow(heatmap_data)
# end
curr_p = p_seq #heatmap_data.p[idx]
curr_k = kappa_seq  #heatmap_data.kappa[idx]

detection_result = fill(n+1,iterations)

for i in 1:iterations
    Y=rlaplace_hd_cpt(n, curr_p, epsilon; cpt=cpt, mu_norm=curr_k)
    r=2
    for t in 2:n
        current_Y=Y[1:t,:]
        delta_t=4*alpha/(t*r*(r+1))
        K=sqrt(log(1/delta_t)/2)
        threshold_s=ceil(Int,((K+sqrt(K^2+4(0.1-2*epsilon)))/(2*(0.1-2*epsilon)))^2)
        limit = floor(Int, t/2)
        if threshold_s<=limit
            X_t = zeros(limit, curr_p)
            for s in 1:limit
                X_t[s, :] .= (current_Y[t - s + 1, :] .- current_Y[s, :]) ./ (sqrt(2))
            end
            for s in threshold_s:limit
                if robust_mean_test(X_t[1:s,:], kappa_guess/sqrt(2), delta_t, 2*epsilon, C_gamma=C_gamma)
                    detection_result[i]=t
                    @goto next_iteration
                end
            end
            r=r+1
        end
    end
    @label next_iteration
end

print(detection_result)

#CSV.write("laplace_cpt.csv", heatmap_data)