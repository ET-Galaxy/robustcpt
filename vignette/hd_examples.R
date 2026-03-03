library(MASS)
library(ggplot2)
library(dplyr)

# ==== 1. Laplace distribution ====
p<-5
delta=0.1
epsilon=0.05
C_gamma<-0.1
iterations<-100
## ===== 1.1 Illustrate effect of misspecification ====
mu1<-1
n<-1000
kappa_seq <- seq(0.1, 3.0, by = 0.1)
error_rates<-numeric(length(kappa_seq))
powers <- numeric(length(kappa_seq))
for (j in seq_along(kappa_seq)) {
  current_kappa <- kappa_seq[j]
  rejections <- 0

  for (i in 1:iterations) {
    # Generate Corrupted Data
    corrupt_index <- rbinom(1, n, prob = epsilon)
    corrupt_data <- mvrnorm(n = corrupt_index, mu = rep(-1, p), Sigma = diag(p))

    # Generate Clean Data
    clean <- rlaplace_hd(n = n - corrupt_index, s=1, p = p, mu = rep(0, p))

    # Combine
    Y <- rbind(clean, corrupt_data)

    # Run the test: RobustMeanTest returns TRUE for rejection
    res <- RobustMeanTest(Y, kappa0 = current_kappa, delta = delta,
                          epsilon = epsilon, C_gamma = C_gamma)

    if (res) {
      rejections <- rejections + 1
    }
  }
  # Calculate proportion of rejections
  error_rates[j] <- rejections / iterations
}

for (j in seq_along(kappa_seq)) {
  current_kappa <- kappa_seq[j]
  rejections <- 0

  for (i in 1:iterations) {
    # Generate Corrupted Data
    corrupt_index <- rbinom(1, n, prob = epsilon)
    corrupt_data <- mvrnorm(n = corrupt_index, mu = rep(-1, p), Sigma = diag(p))

    # Generate Clean Data
    clean <- rlaplace_hd(n = n - corrupt_index, s=1, p = p, mu = rep(mu1/sqrt(p), p))

    # Combine
    Y <- rbind(clean, corrupt_data)

    # Run the test: RobustMeanTest returns TRUE for rejection
    res <- RobustMeanTest(Y, kappa0 = current_kappa, delta = delta,
                          epsilon = epsilon, C_gamma = C_gamma)

    if (res) {
      rejections <- rejections + 1
    }
  }
  # Calculate proportion of rejections
  powers[j] <- rejections / iterations
}


plot_df<- data.frame(kappa0 = kappa_seq, error_rate=error_rates, power= powers)

# For fixed n, fixed mu, vary kappa0 to see the effect of misspecification
# REMOVE "COLOUR" FROM LEGEND
ggplot() +
  geom_line(data = plot_df, aes(x = kappa0, y = error_rate, color = "Type I Error"), size = 1) +
  geom_line(data = plot_df, aes(x = kappa0, y = power, color = "Power"), size = 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", size = 0.8) +
  geom_vline(xintercept = mu1, linetype = "dotted", color = "black", size = 1) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "darkblue") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title = "Type I Error Rate/Power vs. Kappa0 (100 iterations each)",
    subtitle = paste("n =", n, ", p =", p, ", eps =", epsilon, ", C_gamma=", C_gamma),
    x = expression(kappa[0]),
    y = "Type I Error Rate/Power"
  ) +
  theme_minimal()

## === 1.2 Sample complexity vs signal size ====
# Define the grid ranges
n_seq <- seq(500, 3000, by = 100)      # Example values for n
kappa_seq <- seq(0.25, 0.43, by = 0.005)
iterations <- 1000

# Create an empty data frame to store results
heatmap_data <- expand.grid(n = n_seq, kappa0 = kappa_seq)
heatmap_data$error_rate <- 0
heatmap_data$error_rate2 <- 0

# --- Type I error ---
# Nested Loop
for (idx in 1:nrow(heatmap_data)) {
  current_n <- heatmap_data$n[idx]
  current_k <- heatmap_data$kappa0[idx]

  rejections <- 0
  for (i in 1:iterations) {
    # Generate Corrupted Data
    corrupt_index <- rbinom(1, current_n, prob = epsilon)
    corrupt_data <- mvrnorm(n = corrupt_index, mu = rep(-1, p), Sigma = diag(p))
    # Generate Clean Data (Null Hypothesis is True: mu = 0)
    clean <- rlaplace_hd(n = current_n - corrupt_index, p = p, mu = rep(0, p))

    # Combine
    Y <- rbind(clean, corrupt_data)

    # Run the test: RobustMeanTest returns TRUE for rejection
    res <- RobustMeanTest(Y, kappa0 = current_k, delta = delta,
                          epsilon = epsilon, C_gamma = C_gamma)
    if (res == 1) rejections <- rejections + 1
  }
  heatmap_data$error_rate[idx] <- rejections / iterations
}

# --- Power ---
# Nested Loop
for (idx in 1:nrow(heatmap_data)) {
  current_n <- heatmap_data$n[idx]
  current_k <- heatmap_data$kappa0[idx]

  rejections <- 0
  for (i in 1:iterations) {
    # Generate Corrupted Data
    corrupt_index <- rbinom(1, current_n, prob = epsilon)
    corrupt_data <- mvrnorm(n = corrupt_index, mu = rep(-1, p), Sigma = diag(p))
    # Generate Clean Data (Null Hypothesis is True: mu = 0)
    clean <- rlaplace_hd(n = current_n - corrupt_index, p = p, mu = rep(current_k/sqrt(p), p))

    # Combine
    Y <- rbind(clean, corrupt_data)

    # Run the test: RobustMeanTest returns TRUE for rejection
    res <- RobustMeanTest(Y, kappa0 = current_k, delta = delta,
                          epsilon = epsilon, C_gamma = C_gamma)
    if (res == 1) rejections <- rejections + 1
  }
  heatmap_data$error_rate2[idx] <- 1-rejections / iterations
}

# --- Plotting ---
# # heatmap for type I error
# ggplot(heatmap_data, aes(x = factor(n), y = kappa0, fill = error_rate)) +
#   geom_tile() +
#   # Add text labels inside tiles if the grid is small enough
#   geom_text(aes(label = round(error_rate2, 2)), size = 3, color = "white") +
#   # Color scale: high error (Yellow/Purple) vs low error (Dark)
#   scale_fill_viridis_c(option = "magma", labels = scales::percent) +
#   labs(
#     title = "Type I Error Rate Heatmap",
#     subtitle = "Variation across sample size (n) and Kappa0",
#     x = "Sample Size (n)",
#     y = expression(kappa[0]),
#     fill = "Error Rate"
#   ) +
#   theme_minimal()
#
# # heatmap for type II error
# ggplot(heatmap_data, aes(x = factor(n), y = kappa0, fill = error_rate2)) +
#   geom_tile() +
#   # Add text labels inside tiles if the grid is small enough
#   geom_text(aes(label = round(error_rate2, 2)), size = 3, color = "white") +
#   # Color scale: high error (Yellow/Purple) vs low error (Dark)
#   scale_fill_viridis_c(option = "magma", labels = scales::percent) +
#   labs(
#     title = "Type II Error Rate Heatmap",
#     subtitle = "Variation across sample size (n) and Kappa0",
#     x = "Sample Size (n)",
#     y = expression(kappa[0]),
#     fill = "Error Rate"
#   ) +
#   theme_minimal()

# Heatmap for when both error guarantees are satisfied.
heatmap_data <- heatmap_data %>%
  mutate(region = case_when(
    error_rate <= 0.1  & error_rate2 <= 0.1  ~ "Errors controlled",
    error_rate > 0.1 | error_rate2 > 0.1  ~ "Errors not controlled"
  ))

region_colors <- c(
  "Errors controlled" = "#3498db", # Blue
  "Errors not controlled"   = "#e74c3c"  # Alizarin Red
)

ggplot(heatmap_data, aes(x = factor(n), y = kappa0, fill = region)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_manual(values = region_colors) +
  labs(
    title = "Error Rate Thresholds",
    x = "Sample Size (n)",
    y = expression(kappa[0]),
    fill = "Regions"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Sample complexity plot of algorithm
# 1. Calculate the minimum kappa0 for the 'green' condition for each n
boundary_df <- heatmap_data %>%
  # Filter for the 'Green Region' (both error rates below 10%)
  filter(error_rate < 0.1 & error_rate2 < 0.1) %>%
  # Group by sample size
  group_by(n) %>%
  # Find the smallest kappa0 that satisfies the condition
  summarize(min_kappa0 = min(kappa0), .groups = 'drop')

mod<-lm(n~0+I(1/min_kappa0^2), data=boundary_df)
beta_hat <- coef(mod)

# 2. Create a sequence for a smooth curve
n_range <- seq(min(boundary_df$n), max(boundary_df$n), length.out = 200)
curve_df <- data.frame(
  n = n_range,
  # Calculate kappa0 based on the fitted beta: kappa0 = sqrt(beta / n)
  fitted_kappa0 = sqrt(beta_hat / n_range)
)

ggplot(boundary_df, aes(x = n, y = min_kappa0)) +
  geom_point(color = "#27ae60", size = 3) +
  # The Fitted Model Line (The Curve)
  geom_line(data = curve_df, aes(x = n, y = fitted_kappa0),
            color = "firebrick", size = 1.2, linetype = "solid") +
  # Add a shaded area below to show the 'Unsafe Zone'
  geom_ribbon(aes(ymin = 0, ymax = min_kappa0), fill = "#e74c3c", alpha = 0.1) +
  annotate("text", x = mean(boundary_df$n), y = min(boundary_df$min_kappa0)/2,
           label = "Insufficient Robustness", color = "#e74c3c", alpha = 0.6) +
  labs(
    title = "Minimum $\\kappa_0$ Boundary for Valid Inference",
    subtitle = "Smallest $\\kappa_0$ required to maintain Type I error < 10% for both tests",
    x = "Sample Size (n)",
    y = expression(paste("Minimum ", kappa[0]))
  ) +
  theme_minimal()

# ==== 2. t-distribution ====
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

mu1<-1
n<-1000
kappa_seq <- seq(0.1, 3.0, by = 0.1)
error_rates<-numeric(length(kappa_seq))
powers <- numeric(length(kappa_seq))
for (j in seq_along(kappa_seq)) {
  current_kappa <- kappa_seq[j]
  rejections <- 0

  for (i in 1:iterations) {
    # Generate Corrupted Data
    corrupt_index <- rbinom(1, n, prob = epsilon)
    corrupt_data <- mvrnorm(n = corrupt_index, mu = rep(-1, p), Sigma = diag(p))

    # Generate Clean Data
    clean <- rlaplace_hd(n = n - corrupt_index, s=1, p = p, mu = rep(0, p))

    # Combine
    Y <- rbind(clean, corrupt_data)

    # Run the test: RobustMeanTest returns TRUE for rejection
    res <- RobustMeanTest(Y, kappa0 = current_kappa, delta = delta,
                          epsilon = epsilon, C_gamma = C_gamma)

    if (res) {
      rejections <- rejections + 1
    }
  }
  # Calculate proportion of rejections
  error_rates[j] <- rejections / iterations
}

for (j in seq_along(kappa_seq)) {
  current_kappa <- kappa_seq[j]
  rejections <- 0

  for (i in 1:iterations) {
    # Generate Corrupted Data
    corrupt_index <- rbinom(1, n, prob = epsilon)
    corrupt_data <- mvrnorm(n = corrupt_index, mu = rep(-1, p), Sigma = diag(p))

    # Generate Clean Data
    clean <- rlaplace_hd(n = n - corrupt_index, s=1, p = p, mu = rep(mu1/sqrt(p), p))

    # Combine
    Y <- rbind(clean, corrupt_data)

    # Run the test: RobustMeanTest returns TRUE for rejection
    res <- RobustMeanTest(Y, kappa0 = current_kappa, delta = delta,
                          epsilon = epsilon, C_gamma = C_gamma)

    if (res) {
      rejections <- rejections + 1
    }
  }
  # Calculate proportion of rejections
  powers[j] <- rejections / iterations
}


plot_df<- data.frame(kappa0 = kappa_seq, error_rate=error_rates, power= powers)

# For fixed n, fixed mu, vary kappa0 to see the effect of misspecification
# REMOVE "COLOUR" FROM LEGEND
ggplot() +
  geom_line(data = plot_df, aes(x = kappa0, y = error_rate, color = "Type I Error"), size = 1) +
  geom_line(data = plot_df, aes(x = kappa0, y = power, color = "Power"), size = 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", size = 0.8) +
  geom_vline(xintercept = mu1, linetype = "dotted", color = "black", size = 1) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "darkblue") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title = "Type I Error Rate/Power vs. Kappa0 (100 iterations each)",
    subtitle = paste("n =", n, ", p =", p, ", eps =", epsilon, ", C_gamma=", C_gamma),
    x = expression(kappa[0]),
    y = "Type I Error Rate/Power"
  ) +
  theme_minimal()

