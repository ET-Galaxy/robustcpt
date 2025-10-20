using Random, Statistics, StatsBase, DelimitedFiles
include("change_point_utils.jl")

function sigma_t(v)
    if v==2
        return (df / (df - 2))^(1/2)
    elseif v==4
        return (3*df^2 / ((df - 2)*(df-4)))^(1/4)
    elseif v==6
        return (15*df^3 / ((df - 2)*(df-4)*(df-6)))^(1/6)
    end
end



epsilon = 0.0
df = 2.1
v = df-0.1
sigma = sigma_t(v)
alpha = 0.2
n = 2000

mechanism_df = (n, mu=0.0) -> contaminated_sample_t(n, mu; df=df, epsilon=epsilon)

# detected = String[]
# for i in 1:100
#     online_data = change_point_model(n; mechanism=mechanism_df, cpt=nothing)
#     result = rumedian_v(online_data, sigma; v=v, epsilon=epsilon, alpha=alpha, C1=0.19, C2=0.11)   
#     #C1=0.24, C2=0.04 for df=2.1, C1=0.195, C2=0.11 for df=4.1; C1=0.55, C2=0.49 for df=6.1
#     push!(detected, result["method"])
# end

# # Summarize results
# println(countmap(detected))
# writedlm("type_1_error_.csv", countmap(detected), ',')

reps = 100
kappa_sizes = [5*sigma, 2*sigma, sigma, 0.2*sigma, 0.15*sigma, 0.1*sigma]
locations = zeros(length(kappa_sizes), reps)
for k in 1:length(kappa_sizes)
    for i in 1:reps
        online_data = change_point_model(n; mechanism=mechanism_df, cpt=600, kappa=kappa_sizes[k])
        result = rumedian_v(online_data, sigma; v=v, epsilon=epsilon, alpha=alpha, C1=0.195, C2=0.11)
        locations[k,i]=result["location"]
    end
end
println(locations)
# writedlm("locations_v4.csv", locations, ',')