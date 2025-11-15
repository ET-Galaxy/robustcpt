using Statistics, Random
include("change_point_utils.jl")
Random.seed!(1)

n=100
println(contaminated_sample_t(100, 0.0; df=3.0, epsilon=0.05))