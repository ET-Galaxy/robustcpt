using Random, Statistics, StatsBase, DelimitedFiles
include("change_point_utils.jl")

epsilon = 0.0
df = 4.1
sigma = sqrt(df / (df - 2))
alpha = 0.2
v = 4
n = 2000

mechanism_df = (n, mu=0.0) -> contaminated_sample_t(n, mu; df=df, epsilon=epsilon)

detected = String[]
for i in 1:100
    online_data = change_point_model(n; mechanism=mechanism_df, cpt=nothing)
    result = rumedian_v(online_data, sigma; v=v, epsilon=epsilon, alpha=alpha, C1=0.55, C2=0.3)   
    #C1=0.24, C2=0.04 for df=2.1, C1=0.5, C2=0.3 for df=4.1
    push!(detected, result["method"])
end

# Summarize results
println(countmap(detected))
writedlm("type_1_error_v4.csv", countmap(detected), ',')

# reps = 100
# kappa_sizes = [5*sigma, 2*sigma, sigma, 0.2*sigma, 0.15*sigma, 0.1*sigma]
# locations = zeros(length(kappa_sizes), reps)
# for k in 1:length(kappa_sizes)
#     for i in 1:reps
#         online_data = change_point_model(n; mechanism=mechanism_df, cpt=600, kappa=kappa_sizes[k])
#         result = rumedian_v(online_data, sigma; v=v, epsilon=epsilon, alpha=alpha, C1=0.24, C2=0.04)
#         locations[k,i]=result["location"]
#     end
# end
# println(locations)
# writedlm("locations_v4.csv", locations, ',')