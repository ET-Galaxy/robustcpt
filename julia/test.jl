using Random, Statistics, LinearAlgebra, Distributions
include("hd_utils.jl") 
n_=10000
p_=5
Y = rt_hd_cpt(n_, p_, 0.01; cpt=1, mu_norm=2)
robust_mean_test_mom(Y, 2, 0.1, 0.01, C_gamma=0.1, Tu=6)