library(robustcpt)
library(purrr)
# Example 1: t-distribution with deg of freedom df, no contamination
# Underlying distribution
epsilon=0
df=2.1
sigma=sqrt(df/(df-2))
alpha=0.2
v=2
mechanism_df<-partial(contaminated_sample_t, df=df,epsilon=epsilon)
n<-2000

# Type I error
detected<-c()
for (i in 1:100) {
  online_data <- change_point_model(n, mechanism=mechanism_df, cpt=NA)
  detected[i]<- rumedian_v(online_data, v=2, sigma=sigma, alpha=alpha, C1=1, C2=0.04)$method
}
summary(factor(detected))

# Power to detect changepoint
# Large signal
n<-2000
locations<-matrix(nrow=6,ncol=50)
size_kappa<-c(4*sigma, 2*sigma, sigma, sigma/2, sigma/4, sigma/8, sigma/16)
j<-0
for (kap in size_kappa){
  j<-j+1
  for (i in 1:100){
    online_data<-change_point_model(n, cpt=500, kappa=kap)
    test<-rumedian_v(online_data, v=2, sigma=sigma, alpha=alpha, C1=1, C2=0.04)
    locations[j,i]<-test$location
  }
}
