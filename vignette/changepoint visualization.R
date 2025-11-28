contaminated_sample_normal<-function(n,epsilon=0,mu=0){
  all_indices<-1:n
  final_sample<-rnorm(n)+mu
  # final_sample<-rweibull(n,1)
  # flipsign=rbinom(n,size=1,prob=0.5)
  # final_sample<-(-1)^flipsign*final_sample+mu
  contaminated_indices<-all_indices[rbinom(n,1,epsilon)==1]
  for (i in contaminated_indices){
    if (mu>0){
      final_sample[i]<- -rexp(1,0.5)
    }
    final_sample[i]<- rexp(1,0.4)
  }
  return(final_sample)
}

contaminated_sample_t2<-function(n,df=3,epsilon=0,mu=0){
  all_indices<-1:n
  final_sample<-rt(n,df)+mu
  contaminated_indices<-all_indices[rbinom(n,1,epsilon)==1]
  for (i in contaminated_indices){
    final_sample[i]<-rnorm(1,0,100)
  }
  return(final_sample)
}

change_point_model2<-function(n, mechanism=contaminated_sample_t,cpt=NA,kappa=1,
                              epsilon=0){
  if (is.na(cpt)){
    # No change point: generate all n points using mechanism
    final_sample<-mechanism(n=n, epsilon=epsilon)
  } else {
    # Generate first segment (before change point)
    final_sample<-mechanism(n=cpt,epsilon=epsilon)

    # Generate second segment (after change point) with shift mu = kappa
    final_sample[(cpt+1):n]<-mechanism(n=n-cpt, mu=kappa, epsilon=epsilon)
  }
  return(final_sample)
}

data<-change_point_model2(n=120,mechanism=contaminated_sample_normal,
                          cpt=100,kappa=3, epsilon=0.2)
plot(data)
# draw piecewise horizontal lines
segments(x0 = 1, x1 = 100, y0 = 0, y1 = 0, col = "red", lwd = 2)
segments(x0 = 101, x1 = 200, y0 = 3, y1 = 3, col = "blue", lwd = 2)
