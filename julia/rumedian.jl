using Random, Distributions, Statistics, StatsBase
include("change_point_utils.jl")
Random.seed!(12)

epsilon = 0.0
df = 2.1
sigma = sqrt(df / (df - 2))
alpha = 0.2
v = 2
n = 2000

mechanism_df = n -> contaminated_sample_t(n; df=df, epsilon=epsilon)

detected = String[]
for i in 1:1
    online_data = change_point_model(n; mechanism=mechanism_df, cpt=nothing)
    result = rumedian_v(online_data, sigma; v=v, epsilon=epsilon, alpha=alpha, C1=0.8, C2=0.04)
    push!(detected, result["method"])
end

# Summarize results
println(countmap(detected))