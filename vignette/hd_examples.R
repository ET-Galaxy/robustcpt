library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(robustcpt)

# ==== 1. Laplace distribution ====
p<-5
delta=4*0.2/(1500^2*1501)
epsilon=0
iterations<-100
## ===== 1.1 Illustrate effect of misspecification ====
mu1<-1
n<-1500
C_gamma_seq<-0.1
kappa_seq <- seq(0.1, 3.5, by = 0.1)
results <- expand.grid(C_gamma = C_gamma_seq, kappa0 = kappa_seq) %>%
  mutate(error_rate = 0, power = 0)

for (i in 1:nrow(results)) {
  curr_C <- results$C_gamma[i]
  curr_K <- results$kappa0[i]

  rej_null <- 0
  rej_alt  <- 0

  for (m in 1:iterations) {
    # 1. Null Hypothesis (mu = 0)
    if (epsilon>0){
      n_corrupt <- rbinom(1, n, prob = epsilon)
      corrupt_data <- mvrnorm(n = n_corrupt, mu = rep(-1, p), Sigma = diag(p))
      clean <- rlaplace_hd(n = n - n_corrupt, s = 1/sqrt(2), p = p, mu = rep(0, p))
      Y <- rbind(clean, corrupt_data)
    } else {
      Y <- rlaplace_hd(n = n, s = 1/sqrt(2), p = p, mu = rep(0, p))
    }

    if(RobustMeanTest(Y, kappa0 = curr_K, delta = delta, epsilon = epsilon, C_gamma = curr_C)) {
      rej_null <- rej_null + 1
    }

    # 2. Alternative Hypothesis (mu = mu1/sqrt(p))
    if (epsilon>0){
      n_corrupt <- rbinom(1, n, prob = epsilon)
      corrupt_data <- mvrnorm(n = n_corrupt, mu = rep(-1, p), Sigma = diag(p))
      clean <- rlaplace_hd(n = n - n_corrupt, s = 1/sqrt(2), p = p, mu = rep(mu1/sqrt(p), p))
      Y <- rbind(clean, corrupt_data)
    } else {
      Y <- rlaplace_hd(n = n, s = 1/sqrt(2), p = p, mu = rep(mu1/sqrt(p), p))
    }

    if(RobustMeanTest(Y, kappa0 = curr_K, delta = delta, epsilon = epsilon, C_gamma = curr_C)) {
      rej_alt <- rej_alt + 1
    }
  }

  results$error_rate[i] <- rej_null / iterations
  results$power[i]      <- rej_alt / iterations
}
#plot_df<-read.csv("data/Latest_format/hd_misspec_p100.csv")
#results<-rbind(plot_df,results)
#write.csv(results,"data/Latest_format/hd_misspec_p100.csv", row.names = FALSE)

# For fixed n, fixed mu, vary kappa0 to see the effect of misspecification
# Professional Color Palette

# 1. Transform data for plotting
# This creates a 'Metric' column (Type I Error vs Power)
# and a 'Value' column (the actual rates)
plot_df_long <- results %>%
  pivot_longer(
    cols = c(error_rate, power),
    names_to = "Metric",
    values_to = "Rate"
  ) %>%
  mutate(
    # Clean up names for the legend
    Metric = ifelse(Metric == "error_rate", "Type I Error", "Power"),
    C_gamma_fact = factor(C_gamma)
  )

# 2. Create the plot
ggplot(plot_df_long, aes(x = kappa0, y = Rate, color = C_gamma_fact, linetype = Metric)) +
  geom_line(size = 1.1) +
  # Add a horizontal line for the nominal alpha/delta level
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "blue", size = 0.8) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", size = 0.8) +
  geom_vline(xintercept = mu1, color = "black", size = 1) +
  scale_color_viridis_d(option = "plasma", end = 0.8) + # Professional color palette
  labs(
    x = expression(Signal~size~input~(kappa[0])),
    y = "Empirical Probability",
    color = expression(C[gamma]),
    linetype = "Metric"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

ggplot(plot_df_long, aes(x = kappa0, y = Rate, color = Metric)) +
  geom_line(size = 1) +
  facet_wrap(~C_gamma_fact, labeller = label_both) +
  theme_bw()

## === 1.2 Sample complexity vs signal size ====
# Define the grid ranges
C_gamma=0.1
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
heatmap_data<-read.csv("data/Latest_format/rawdata/hd_test/th1e5p100.csv")

# Heatmap for when both error guarantees are satisfied.
heatmap_data <- heatmap_data %>%
  mutate(region = case_when(
    type1_error <= 0.1  & type2_error <= 0.1  ~ "Reliably Detectable",
    type1_error > 0.1 | type2_error > 0.1  ~ "Not Reliably Detectable"
  ))

region_colors <- c(
  "Reliably Detectable" = "#3498db", # Blue
  "Not Reliably Detectable"   = "#e74c3c"  # Alizarin Red
)

ggplot(heatmap_data, aes(x = factor(n), y = kappa0, fill = region)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_manual(values = region_colors) +
  scale_x_discrete(breaks = seq(500, 3000, by = 500)) +
  labs(
    x = "n",
    y = expression(kappa[0]),
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top"
  )

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
df<-2.1
p<-100
delta<-0.1
epsilon<-0.05
C_gamma<-0.1
iterations<-1000
## ===== 2.1 Illustrate effect of misspecification ====
mu1<-1
n<-1000
kappa_seq <- seq(0.1, 3.5, by = 0.1)
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
    clean <- rt_hd(n = n - corrupt_index, df=df, sd=1, p = p, mu = rep(0, p))

    # Combine
    Y <- rbind(clean, corrupt_data)

    # Run the test: RobustMeanTest returns TRUE for rejection
    res <- RobustMeanTest(Y, kappa0 = current_kappa, delta = delta,
                          epsilon = epsilon, C_gamma = C_gamma, finite_moment = TRUE)

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
    clean <- rt_hd(n = n - corrupt_index, df=2.1, sd=1, p = p, mu = rep(mu1/sqrt(p), p))

    # Combine
    Y <- rbind(clean, corrupt_data)

    # Run the test: RobustMeanTest returns TRUE for rejection
    res <- RobustMeanTest(Y, kappa0 = current_kappa, delta = delta,
                          epsilon = epsilon, C_gamma = C_gamma, finite_moment = TRUE)

    if (res) {
      rejections <- rejections + 1
    }
  }
  # Calculate proportion of rejections
  powers[j] <- rejections / iterations
}

plot_df<- data.frame(kappa0 = kappa_seq, error_rate=error_rates, power= powers)
#plot_df<-read.csv("data/Latest_format/hd_misspec_p100_t.csv")

write.csv(plot_df,"data/Latest_format/hd_misspec_p100_t.csv", row.names = FALSE)

# For fixed n, fixed mu, vary kappa0 to see the effect of misspecification
# Professional Color Palette

ggplot(plot_df, aes(x = kappa0)) +
  geom_line(aes(y = error_rate, color = "Type I Error"), linewidth = 1) +
  geom_line(aes(y = power, color = "Power"), linewidth = 1) +

  geom_vline(xintercept = mu1, color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_hline(yintercept = delta, linetype = "dotted", color = "blue", linewidth = 1) +
  geom_hline(yintercept = 1 - delta, linetype = "dotted", color = "red", linewidth = 1)+

  # Scales
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_x_continuous(limits = c(0.1, 3.5)) +
  # Labels
  labs(
    #title = bquote(n == .(n) ~ "," ~ p == .(p) ~ "," ~ C[gamma] == .(C_gamma)),
    x = expression(kappa[0]),
    y = "Empirical Probability",
    color = NULL
  ) +
  # Theme Customization
  theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size=13),
    plot.title = element_text(hjust = 0.5)
  )
