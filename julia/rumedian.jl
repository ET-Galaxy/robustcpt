using Distributed
addprocs()  # add worker processes

@everywhere begin
    using Random
    # define contaminated_sample_t, change_point_model, and rumedian in Julia
    
end

n = 500
alpha = 0.25
df = 2.1
epsilon = 0.0

@sync @distributed (vcat) for i in 1:100
    online_data = change_point_model(n; df=df, epsilon=epsilon)
    rumedian(online_data, alpha=alpha, C1=1.5, C2=2).method
end