using Statistics, Random
include("change_point_utils.jl")
Random.seed!(1)

X = randn(2000)
result = rume(X; epsilon=0.0, delta=0.1)
println(result)