library(purrr)
# Example 1: t-distribution with deg of freedom df, no contamination
# Underlying distribution
epsilon=0
df=2.1
sigma=sqrt(df/(df-2))
alpha=0.2
v=2
mechanism_df<-partial(contaminated_sample_t, df=df,epsilon=epsilon)
n=1200

# Type I error
detected<-c()
for (i in 1:40) {
  online_data <- change_point_model(n, mechanism=mechanism_df, cpt=NA)
  detected[i]<- rumedian_v(online_data, v=2, sigma=sigma, alpha=alpha, C1=1, C2=0.035)$method
}
summary(factor(detected))

# Power to detect changepoint
# Large signal
n<-1500
detected<-c()
locations<-c()
for (i in 1:10){
  online_data<-change_point_model(n, cpt=500, kappa=1)
  test<-rumedian(online_data, alpha=alpha)
  detected[i]<-test$method
  locations[i]<-test$location
}
summary(factor(detected))
hist(locations)

