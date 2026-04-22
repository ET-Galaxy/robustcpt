library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(robustcpt)

# ==== 1. Laplace distribution ====
p<-5
delta=0.1 #4*0.2/(1500^2*1501)
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

#results<-rbind(plot_df,results)
#results<-results[order(results$C_gamma), ]
#write.csv(results,"data/Latest_format/hd_misspec_p100_t.csv", row.names = FALSE)
#plot_df<-read.csv("data/Latest_format/hd_misspec_p100_tmore.csv")
results<-read.csv("data/Latest_format/rawdata/hd_test/v4_sensitive_p10.csv")

# For fixed n, fixed mu, vary kappa0 to see the effect of misspecification
# Professional Color Palette

# 1. Transform data for plotting
# This creates a 'Metric' column (Type I Error vs Power)
# and a 'Value' column (the actual rates)
plot_df_long <- results %>%
#filter(C_gamma %in% c(0.01,0.03,0.05)) %>%
  filter(kappa0<=3) %>%
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
  geom_line(size = 1.1, alpha = 0.7) +
  # Add a horizontal line for the nominal alpha/delta level
  geom_hline(yintercept = 0.9, linetype = "solid", color = "blue", linewidth = 0.7) +
  geom_hline(yintercept = 0.1, linetype = "solid", color = "red", linewidth = 0.7) +
  geom_vline(xintercept = mu1, color = "black", linewidth = 1.5) +
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
  facet_wrap(
    ~C_gamma_fact,
    ncol = 4,
    labeller = labeller(
      C_gamma_fact = function(x) {
        paste0("C[gamma] == ", x)
      },
      .default = label_parsed
    )
  ) +
  geom_vline(xintercept = 1, color = "black", linewidth = 1)+
  geom_hline(yintercept = 0.9, linetype = "solid", color = "blue", linewidth = 0.7) +
  geom_hline(yintercept = 0.1, linetype = "solid", color = "red", linewidth = 0.7) +
  labs(
    x = expression(Signal~size~input~(kappa[0])),
    y = "Empirical Probability"
  )+
  theme_bw()+
  theme(
  legend.position = "top"
  )

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

# --- Plotting (Vary kappa, fix p) ---
heatmap_data<-read.csv("data/Latest_format/rawdata/hd_test/v4e4p100.csv")
heatmap_data2<-read.csv("data/Latest_format/rawdata/hd_test/v4e4p100more.csv")
heatmap_data<-rbind(heatmap_data,heatmap_data2)
#write.csv(heatmap_data,"data/Latest_format/v2e4p100.csv", row.names = FALSE)

# Heatmap for when both error guarantees are satisfied.
heatmap_data <- heatmap_data %>%
  filter(kappa0>=0.25)%>%
  mutate(region = case_when(
    type1_error <= 0.1  & type2_error <= 0.1  ~ "Reliably Detectable",
    type1_error > 0.1 | type2_error > 0.1  ~ "Not Reliably Detectable"
  ))

region_colors <- c(
  "Reliably Detectable" = "#3498db", # Blue
  "Not Reliably Detectable"   = "#e74c3c"  # Alizarin Red
)

ggplot(heatmap_data, aes(x = kappa0, y = factor(n), fill = region)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_manual(values = region_colors) +
  scale_x_continuous(breaks = seq(0.2, 1.2, by = 0.05)) +
  scale_y_discrete(breaks = seq(500, 3000, by = 500)) +
  labs(
    x = expression(kappa),
    y = "n",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top"
  )

# --- Plotting (Vary p, fix kappa) ---
heatmap_data<-read.csv("data/Latest_format/rawdata/hd_test/v4varyp.csv")
#heatmap_data2<-read.csv("data/Latest_format/rawdata/hd_test/th1varyp2.csv")
heatmap_data<-rbind(heatmap_data,new_df1)
#write.csv(heatmap_data,"data/Latest_format/th1varyp.csv", row.names = FALSE)

new_df <- expand.grid(
  n = seq(2050, 2200, by = 50),
  p = seq(50, 1000, by = 50)
)
new_df$type1_error <- 0
new_df$type2_error <- 0

new_df1 <- expand.grid(
  n = seq(500, 1450, by = 50),
  p = seq(1050, 1250, by = 50)
)
new_df1$type1_error <- 1
new_df1$type2_error <- 1


# Heatmap for when both error guarantees are satisfied.
heatmap_data <- heatmap_data %>%
  filter(p>=100) %>%
  filter(n<=1000) %>%
  mutate(region = case_when(
    type1_error <= 0.1  & type2_error <= 0.1  ~ "Reliably Detectable",
    type1_error > 0.1 | type2_error > 0.1  ~ "Not Reliably Detectable",
    is.na(type1_error) | is.na(type2_error) ~ "Not Reliably Detectable"
  ))

region_colors <- c(
  "Reliably Detectable" = "#3498db", # Blue
  "Not Reliably Detectable"   = "#e74c3c"  # Alizarin Red
)

ggplot(heatmap_data, aes(x = p, y = factor(n), fill = region)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_manual(values = region_colors) +
  scale_x_continuous(breaks = seq(100, 1200, by = 100)) +
  scale_y_discrete(breaks = seq(200, 3000, by = 200)) +
  labs(
    x = "p",
    y = "n",
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
  filter(type1_error < 0.1 & type2_error < 0.1) %>%
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

ggplot(boundary_df, aes(y = n, x = min_kappa0)) +
  geom_point(color = "#27ae60", size = 2) +
  # The Fitted Model Line (The Curve)
  geom_line(data = curve_df, aes(y = n, x = fitted_kappa0),
            color = "firebrick", size = 1.2, linetype = "solid") +
  # Add a shaded area below to show the 'Unsafe Zone'
  geom_ribbon(aes(xmin = 0, xmax = min_kappa0), fill = "#e74c3c", alpha = 0.1) +
  annotate("text", y = mean(boundary_df$n), x = min(boundary_df$min_kappa0)/2,
           label = "Insufficient Robustness", color = "#e74c3c", alpha = 0.6) +
  labs(
    title = "Minimum $\\kappa_0$ Boundary for Valid Inference",
    subtitle = "Smallest $\\kappa_0$ required to maintain Type I error < 10% for both tests",
    y = "Sample Size (n)",
    x = expression(paste("Minimum ", kappa[0]))
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



# ===== 3. HD CPT =====
cpt=3000
n=9000
get_proportions <- function(x) {
  total <- length(x)
  c(
    not_detected = sum(x >= n+1) / total,
    false_alarm = sum(x >= 0 & x <= cpt) / total,
    detected = sum(x > cpt & x <= n) / total
  )
}

# Detectability
rawdata<-read.csv("data/Latest_format/rawdata/hd_test/hdcpt_t_p10.csv")
rawdata2<-read.csv("data/Latest_format/rawdata/hd_test/hdcpt_t_p10_test.csv")
rawdata3<-read.csv("data/Latest_format/rawdata/hd_test/hdcpt_t_p10more1.csv")

fulldata<-rbind(rawdata,rawdata3,rawdata2)
#write.csv(fulldata,"data/Latest_format/hdcpt_p10t_full.csv", row.names = FALSE)

data1 <- fulldata %>%
  pivot_longer(-kappa, values_to = "stoppingT")

props <- data1 %>%
  filter(kappa<=0.7) %>%
  group_by(kappa) %>%
  reframe(
    tibble::as_tibble_row(get_proportions(stoppingT))
  )
props$false_alarm[1]<-props$false_alarm[1]+props$detected[1]
props$detected[1]<-0


props_longer<- props %>%
  pivot_longer(
    cols = c(not_detected, false_alarm, detected),
    names_to = "category",
    values_to = "proportion") %>%
  mutate(category = factor(category,
                           levels = c("not_detected","detected","false_alarm")))

snr_step <- min(diff(sort(unique(props_longer$kappa))))
ggplot(props_longer,
      aes(x = kappa,
          y = proportion,
          fill = category)) +
  geom_col(width = snr_step * 0.98,
           colour = "white",
           linewidth = 0.2) +
  geom_hline(yintercept = 0.1,
             linetype = "dashed",
             linewidth = 0.5) +
  scale_x_continuous(
    breaks = seq(0, 1.2, by = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = expression(kappa),
    y = "Empirical probability"
  ) +
  scale_fill_manual(
    values = c(
      not_detected = "#999999",
      detected     = "#0072B2",
      false_alarm  = "#D55E00"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
#    axis.text.x = element_text(angle=45, vjust=0.5)
  )


# Detection delay
detected_long <- data1 %>%
  filter(kappa >=0.5)%>%
  filter(kappa <=2)%>%
  filter(stoppingT > cpt)

result <- detected_long %>%
  group_by(kappa) %>%
  summarise(meanT = mean(stoppingT), sdT=sd(stoppingT))

result$meanT<-result$meanT-cpt

ggplot(result, aes(x = kappa, y = meanT)) +
  geom_point(size = 1.6, colour = "black") +
  # geom_line(data = grid,
  #           aes(y = fit_inv_sq,
  #               colour = "Regime 2"),
  #           linewidth = 0.9) +
  # geom_line(data = grid2,
  #           aes(y = fit_log_inv,
  #               colour = "Regime 3"),
  #           linewidth = 0.9) +
  scale_x_continuous(
    breaks = seq(0, 5, by = 0.2))+
  scale_y_continuous(
    limits=c(0,3000),
    breaks = seq(0, 3000, by = 500))+
  labs(
    x = expression(kappa),
    y = "Mean detection delay",
    colour = NULL
  ) +
  # scale_colour_manual(
  #   breaks = c("Regime 2", "Regime 3"),
  #   values = c(
  #     "Regime 2" = "pink",
  #     "Regime 3" = "#2C7BB6"
  #   )
  # ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.6, "cm")
  )
