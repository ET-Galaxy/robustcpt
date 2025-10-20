using Statistics, Random
include("change_point_utils.jl")
Random.seed!(1)

mu=0.0
n=100
df=3.0
epsilon=0.1
t_samples = rand(TDist(df), n) .+ mu
contam_mask = rand(Binomial(1, epsilon), n) .== 1
contam_idx = findall(contam_mask)
if !isempty(contam_idx)
    t_samples[contam_idx] = rand(Normal(0, 100), length(contam_idx))
end
t_samples