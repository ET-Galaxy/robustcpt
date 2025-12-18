using Random, Statistics, StatsBase, DelimitedFiles
include("change_point_utils.jl")

# ==== Weibull ====
epsilon = 0.1
sigma = 1
theta = 1
alpha = 0.2
n = 2400
C1=0.53 # RUME constant
C2=0.075 # median constant
#eps=0.1: C1=0.21, C2=0.013 for theta=1
mechanism_df = (n, mu=0.0) -> contaminated_laplace(n, mu; epsilon=epsilon)

detected = String[]
for i in 1:100
    online_data = change_point_model(n; mechanism=mechanism_df, cpt=nothing)
    result = rumedian_theta(online_data, sigma; theta=theta, epsilon=epsilon, alpha=alpha, C1=C1, C2=C2)   
    push!(detected, result["method"])
end
println(countmap(detected))

# reps = 200
# kappa_sizes = 0 # sigma*(0)
# locations = zeros(length(kappa_sizes), reps)
# for k in 1:length(kappa_sizes)
#     for i in 1:reps
#         online_data = change_point_model(n; mechanism=mechanism_df, cpt=600, kappa=kappa_sizes[k])
#         result = rumedian_v(online_data, sigma; v=v, epsilon=epsilon, alpha=alpha, C1=C1, C2=C2)
#         locations[k,i]=result["location"]
#     end
# end
# println(locations)
# writedlm("locations_v2e10new3.csv", locations, ',')

# c2=0.021
# done=false
# while !done
#     detected = String[]
#     for i in 1:100
#         online_data = change_point_model(n; mechanism=mechanism_df, cpt=nothing)
#         result = rumedian_v(online_data, sigma; v=v, epsilon=epsilon, alpha=alpha, C1=0.24, C2=c2)   
#         #eps=0: C1=0.24, C2=0.04 for df=2.1, C1=0.195, C2=0.11 for df=4.1; C1=0.55, C2=0.49 for df=6.1
#         #eps=0.05: 
#         push!(detected, result["method"])
#     end
#     nocpt=countmap(detected)["no changepoint"]
#     if nocpt>=95
#         global c2=c2/2
#     elseif nocpt<=89
#         global c2=1.5*c2
#     else
#         global done=true
#         println(countmap(detected))
#         writedlm("t1e_v2e1.csv", countmap(detected), ',')
#     end
# end