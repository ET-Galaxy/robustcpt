# Cutoff (Easy case - no contamination)
# RUME Constant determined as 1.59
s=280
v=2
delta=0.000005
epsilon=0
mean_dist<-c()
const<-c()
for (trial1 in 1:1){
  for (trial in 1:200000){
    online_data1=contaminated_sample(s,epsilon=epsilon)
    online_data2=contaminated_sample(s,epsilon=epsilon)
    vareps<-max(epsilon, 2*log(1/delta)/s)
    mean_dist[trial]<-abs(rume(online_data1,epsilon=0,delta=delta)-rume(online_data2,epsilon=0,delta=delta))
    }
  theoretical_bound<-2*sqrt(3)*max(vareps^(1-1/v),sqrt(2*log(1/delta)/s))
  const[trial1]<-quantile(mean_dist,1-delta)/theoretical_bound
}

# RUME Constant determined as 1.2
s=400
v=4
delta=0.000001
epsilon=0
mean_dist<-c()
const<-c()
for (trial1 in 1:1){
  for (trial in 1:500000){
    online_data1=contaminated_sample_t(s,epsilon=epsilon, df=v+1)
    online_data2=contaminated_sample_t(s,epsilon=epsilon, df=v+1)
    vareps<-max(epsilon, 2*log(1/delta)/s)
    mean_dist[trial]<-abs(rume(online_data1,epsilon=0,delta=delta)-rume(online_data2,epsilon=0,delta=delta))
  }
  theoretical_bound<-2*sqrt(3)*max(vareps^(1-1/v),sqrt(2*log(1/delta)/s))
  const[trial1]<-quantile(mean_dist,1-delta)/theoretical_bound
}
const

# Median constant for v=2
# Determined around 0.13
s=30
v=2
delta=0.0001
epsilon=0
mean_dist<-c()
const<-c()
for (trial1 in 1:3){
  for (trial in 1:40000){
    online_data1=contaminated_sample_t(s,epsilon=epsilon)
    online_data2=contaminated_sample_t(s,epsilon=epsilon)
    mean_dist[trial]<-abs(median(online_data1)-median(online_data2))
  }
  theoretical_bound<-2*sqrt(3)*(2*exp(1)*(1-epsilon)/(delta^(2/s)-2*exp(1)*epsilon))^(1/v)
  const[trial1]<-quantile(mean_dist,1-delta)/theoretical_bound
}

# Median constant for v=4
# Determined around 0.3
s=30
v=4
delta=0.00001
epsilon=0
mean_dist<-c()
const<-c()
for (trial1 in 1:3){
  for (trial in 1:50000){
    online_data1=contaminated_sample_t(s,epsilon=epsilon, df=v+1)
    online_data2=contaminated_sample_t(s,epsilon=epsilon, df=v+1)
    mean_dist[trial]<-abs(median(online_data1)-median(online_data2))
  }
  theoretical_bound<-2*sqrt(3)*(2*exp(1)*(1-epsilon)/(delta^(2/s)-2*exp(1)*epsilon))^(1/v)
  const[trial1]<-quantile(mean_dist,1-delta)/theoretical_bound
}
const

