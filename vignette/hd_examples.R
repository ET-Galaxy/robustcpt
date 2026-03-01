library(MASS)

# ==== Laplace distribution ====
n<-1000
p<-4
delta=0.1
epsilon=0.05

iterations <- 100  # Number of trials per kappa0
kappa_seq <- seq(0.1, 2.0, by = 0.1)

# Initialize a vector to store the error rates
type_1_errors <- numeric(length(kappa_seq))

# --- Nested Simulation Loop ---
for (j in seq_along(kappa_seq)) {
  current_kappa <- kappa_seq[j]
  rejections <- 0

  for (i in 1:iterations) {
    # Generate Corrupted Data
    corrupt_index <- rbinom(1, n, prob = epsilon)
    corrupt_data <- mvrnorm(n = corrupt_index, mu = rep(1, p), Sigma = diag(p))

    # Generate Clean Data (Null Hypothesis is True: mu = 0)
    clean <- rlaplace_hd(n = n - corrupt_index, p = p, mu = rep(0, p))

    # Combine
    Y <- rbind(clean, corrupt_data)

    # Run the test
    # (Assuming RobustMeanTest returns 1/TRUE for rejection and 0/FALSE for fail-to-reject)
    res <- RobustMeanTest(Y, kappa0 = current_kappa, delta = delta, epsilon = epsilon, C_gamma = 1)

    if (res == 1 || res == TRUE) {
      rejections <- rejections + 1
    }
  }

  # Calculate proportion of rejections (Type I Error Rate)
  type_1_errors[j] <- rejections / iterations
}

for (i in 1:100){
  corrupt_index<-rbinom(1,n,prob=epsilon)
  corrupt_data<-mvrnorm(n=corrupt_index, mu=rep(1,p), Sigma=diag(p))
  clean<-rlaplace_hd(n=n-corrupt_index, p=p, mu=rep(0,p))
  Y<-rbind(clean,corrupt_data)
  res[i]<-RobustMeanTest(Y,kappa0=0.5, delta=delta,epsilon=epsilon,C_gamma = 1)
}
summary(factor(res))


# ==== t-distribution ====
n<-1000
p<-4
delta=0.1
epsilon=0.05
res<-c()
for (i in 1:100){
  corrupt_index<-rbinom(1,n,prob=epsilon)
  corrupt_data<-contaminated_sample_t_hd(n=corrupt_index, p=p, mu=rep(0,p))
  clean<-contaminated_sample_t_hd(n=n-corrupt_index, p=p, mu=rep(0,p))
  Y<-rbind(clean,corrupt_data)
  res[i]<-RobustMeanTest(Y,kappa0=0.3, delta=delta,epsilon=epsilon,C_gamma = 1)
}
summary(factor(res))
