# library(purrr)
set.seed(12)
# Example 1: t-distribution with deg of freedom df, no contamination
# Underlying distribution
epsilon=0
df=3
sigma=sqrt(df/(df-1))
alpha=0.25
v=2
# mechanism<-partial(contaminated_sample_t, df=v,epsilon=epsilon)
n=2000

# Type I error
detected<-c()
for (i in 1:100){
  online_data<-change_point_model(n, cpt=NA, kappa=0)
  detected[i]<-rumedian(online_data, alpha=alpha)$method
}
summary(factor(detected))
