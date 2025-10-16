using Statistics, Random
include("change_point_utils.jl")
Random.seed!(1)

X = randn(200)
result = rume(X; epsilon=0, delta=0.05)
println(result)