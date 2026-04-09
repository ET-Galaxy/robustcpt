using Distributed, DataFrames

# Define your grid
kappa_seq = 0
reps = 1

@everywhere begin
    using Random, Statistics, LinearAlgebra, Distributions
    include("hd_utils.jl") 

    # Define constants for the workers
    const p_ = 10
    const kappa_guess_ = 0.5
    const Tu_= 2
    const epsilon_ = 0
    const alpha_ = 0.1   
    const C_gamma_ = 0.05
    const n_=70000
    const cpt_=35000
    
    function run_one_sim(idx, curr_k, i)
        Y = rt_hd_cpt(n_, p_, epsilon_; cpt=cpt_, mu_norm=curr_k)
        r = 2
        low_limit = 17384
        for t in 34768:n_
            current_Y = Y[1:t, :]
            delta_t = 4*alpha_/(t*r*(r+1))
            limit = fld(t, 2)
            k = ceil(Int, 8*log(1/delta_t))
            for s in low_limit:limit
                samplesize=fld(s,2)
                block_size = div(samplesize, k)
                q = 2*epsilon_ + min(p_/(20*block_size),1/20)
                u= q + sqrt(2*q*log(16k)/block_size)+2*log(16k)/(3*block_size)
                if u <= 0.1 
                    println(s, t)
                    X_t = zeros(samplesize, p_)

                    for c in 1:samplesize
                        @inbounds X_t[c, :] .= (current_Y[t-c+1, :] .- current_Y[c, :]) ./ sqrt(2)
                    end

                    if robust_mean_test_mom(X_t, kappa_guess_/sqrt(2),
                                            delta_t, 2*epsilon_,
                                            C_gamma=C_gamma_, Tu=Tu_)
                        return (idx, i, t)
                    end
                    r += 1
                else 
                   low_limit=s+1
                end
            end
        end
        return (idx, i, n_+1)  # no detection case
    end
end

# 1. Create a single master DataFrame
df = DataFrame(kappa = kappa_seq)
for i in 1:reps
    df[!, Symbol("sim_$i")] = zeros(length(kappa_seq))
end

for idx in 1:nrow(df)
    curr_k = df.kappa[idx]

    sim_results = pmap(i -> run_one_sim(idx, curr_k, i), 1:reps)

    for (idx_res, i_res, t_res) in sim_results
        df[idx_res, Symbol("sim_$i_res")] = t_res
    end
end

# 3. Save as a single combined file
print(df)
#CSV.write("cpt/hdcpt_t_p10.csv", df)
#println("All results saved to cpt/hdcpt_t_p10.csv")

