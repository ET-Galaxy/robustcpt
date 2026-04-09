using Distributed, DataFrames
include("hd_utils.jl")

# Parameters
curr_p = 10
kappa_seq = 10000:1000:20000
kappa_guess = 0.5
Tu=2.9
iterations = 1
n=2000
cpt=500
epsilon, alpha, C_gamma = 0.01, 0.1, 0.05

# Grid Setup
df = DataFrame(kappa = kappa_seq)
for i in 1:iterations
    df[!, Symbol("sim_$i")] = fill(n+1, length(kappa_seq))
end

for idx in 1:nrow(df)
    curr_k = df.kappa[idx]
    for i in 1:iterations
        Y=rlaplace_hd_cpt(n, curr_p, epsilon; cpt=cpt, mu_norm=curr_k)
        r=2
        for t in 2:n
            current_Y=Y[1:t,:]
            delta_t=4*alpha/(t*r*(r+1))
            limit = floor(Int, t/2)
            L=log(1/delta_t)
            A=1+2*L/3
            D=0.1-2*epsilon
            threshold_s=ceil(Int,(A^2-2*L)/(A*D+L*2*epsilon-sqrt(L*(4*L*epsilon^2+4*A*epsilon*D+2*D^2))))
            if threshold_s<=limit
                X_t = zeros(limit, curr_p)
                for s in 1:limit
                    X_t[s, :] .= (current_Y[t - s + 1, :] .- current_Y[s, :]) ./ (sqrt(2))
                end
                for s in threshold_s:limit
                    if robust_mean_test(X_t[1:s,:], kappa_guess/sqrt(2), delta_t, 2*epsilon, C_gamma=C_gamma, Tu=Tu)
                        df[idx, Symbol("sim_$i")]=t
                        @goto next_iteration
                    end
                end
                r=r+1
            end
        end
        @label next_iteration
    end
end

print(df)
#CSV.write("laplace_cpt.csv", heatmap_data)