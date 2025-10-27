using Random, Statistics, StatsBase, DelimitedFiles
include("change_point_utils.jl")

function sigma_t(df,v)
    if v==2
        return (df / (df - 2))^(1/2)
    elseif v==4
        return (3*df^2 / ((df - 2)*(df-4)))^(1/4)
    elseif v==6
        return (15*df^3 / ((df - 2)*(df-4)*(df-6)))^(1/6)
    end
end

epsilon = 0.1
df = 2.1
v = 2
sigma = sigma_t(df,v)
alpha = 0.2
n = 2000
C1=0.22 # RUME constant
C2=0.05 # median constant
#eps=0: C1=0.22, C2=0.04 for df=2.1; C1=0.195, C2=0.11 for df=4.1; C1=0.17, C2=0.14 for df=6.1
#eps=0.05: C1=0.22, C2=0.024 for df=2.1; C1=0.2, C2=0.05 for df=4.1; C1=0.19, C2=0.077 for df=6.1
#eps=0.1: C1=0.22, C2=0.05 for df=6.1
mechanism_df = (n, mu=0.0) -> contaminated_sample_t(n, mu; df=df, epsilon=epsilon)

# detected = String[]
# for i in 1:200
#     online_data = change_point_model(n; mechanism=mechanism_df, cpt=nothing)
#     result = rumedian_v(online_data, sigma; v=v, epsilon=epsilon, alpha=alpha, C1=C1, C2=C2)   
#     push!(detected, result["method"])
# end
# println(countmap(detected))
# writedlm("t1e_v6e10.csv", countmap(detected), ',')

reps = 100
kappa_sizes = [4*sigma, 2*sigma, sigma, 0.5*sigma, 0.25*sigma, sigma/8, sigma/16]
locations = zeros(length(kappa_sizes), reps)
for k in 1:length(kappa_sizes)
    for i in 1:reps
        online_data = change_point_model(n; mechanism=mechanism_df, cpt=600, kappa=kappa_sizes[k])
        result = rumedian_v(online_data, sigma; v=v, epsilon=epsilon, alpha=alpha, C1=C1, C2=C2)
        locations[k,i]=result["location"]
    end
end
println(locations)
writedlm("locations_v6e10.csv", locations, ',')

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