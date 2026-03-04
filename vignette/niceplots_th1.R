library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# Function to compute proportions for each row
get_proportions <- function(x) {
  total <- length(x)
  c(
    not_detected = sum(x >= 2400) / total,
    false_alarm = sum(x >= 0 & x <= 600) / total,
    detected = sum(x > 600 & x < 2400) / total
  )
}

# Read CSV with no header
data1 <- read.csv("data/Latest_format/locations_th1e10all.csv")

# ==== Data treatment ====
# data2 <- read.csv("data/Latest_format/rawdata/locations_th1e10_R2R3.csv", header = FALSE)
data_bars <- read.csv("data/Latest_format/rawdata/locations_th1e10_R3R4.csv", header = FALSE)
data_bars[data_bars == -1] <- 2401

exponent_seq <- seq(0, 20, by = 0.25)
powers_of_two <- 2.0^exponent_seq

detected_long <- data_bars %>%
  mutate(snr = powers_of_two) %>%  #seq(from = 0.085, to = 0.495, by = 0.01)
  pivot_longer(-snr, values_to = "stoppingT")

colnames(detected_long)[2]<-"trial"

#total<-as.data.frame(detected_long)
total<-rbind(data1,as.data.frame(detected_long))
# total<-total[order(total$snr),]

write.csv(total,"data/Latest_format/locations_th1e10all.csv", row.names = FALSE)
#write.csv(detected_long,"data/Latest_format/locations_v2e10all.csv", row.names = FALSE)

# ==== Proportion plot ====
# Group by snr and compute proportions
props <- data1 %>%
  filter(snr<=0.23) %>%
  group_by(snr) %>%
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

snr_step <- min(diff(sort(unique(props_longer$snr))))

ggplot(props_longer,
       aes(x = snr,
           y = proportion,
           fill = category)) +
  geom_col(width = snr_step * 0.98,
           colour = "white",
           linewidth = 0.2) +
  geom_hline(yintercept = 0.2,
             linetype = "dashed",
             linewidth = 0.4) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = expression(kappa/phi),
    y = "Empirical probability",
    title = "\u03B8 = 1, \u03b5 = 0.1"
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
  )

ggsave(
  filename = "data/presentation illustrations/th1e10_R1R2plot.jpg",
  plot = last_plot(),         # Saves the last ggplot you displayed
  device = agg_png,           # This is the magic part
  width = 5.43,
  height = 3.92,
  units = "in",
  res = 300                   # High resolution for publication
)

# ==== Boxplot ====
detected_long <- data1 %>%
  filter(stoppingT > 600, snr>=0.25) %>%
  mutate(snr = factor(snr))

ggplot(detected_long, aes(x = snr, y = stoppingT-600)) +
  geom_boxplot(fill = "steelblue", outlier.color = "red") +
  labs(x = "Signal size", y = "Detection delay", title = "Detection delay against signal size (v=2, eps=0.05)") +
  theme_minimal()

# ====== Mean plot =====
# log-log
long_data <- data1 %>%
  filter(snr >=0.23)%>%
  filter(snr <=1)%>%
  filter(stoppingT > 600)

result <- long_data %>%
  group_by(snr) %>%
  summarise(meanT = mean(stoppingT), sdT=sd(stoppingT))

result$meanT<-result$meanT-600

mod1<-lm(meanT~I(snr^(-2.63))+0, data=result[result$snr<=0.44,], weights = 1/sdT)
mod2<-lm(meanT~I(snr^(-1))+0, data=result[result$snr>=0.51,], weights = 1/sdT)

# prediction grid
grid <- data.frame(snr = seq(0.23,0.44,length.out = 1000))
grid2 <- data.frame(snr = seq(0.51,1,length.out = 1000))

grid$fit_inv <- predict(mod1, newdata = grid)
grid2$fit_inv_sq <- predict(mod2, newdata = grid2)

ggplot(result, aes(x = snr, y = meanT)) +
  theme_bw(base_size = 14) +
  annotation_logticks() +
  geom_point(alpha = 0.7, size = 1.8) +
  # Using the new labels here
  geom_line(data = grid, aes(y = fit_inv, color = "Inverse-square (slope = -2.63)"),
            linewidth = 1, linetype = "solid") +
  geom_line(data = grid2, aes(y = fit_inv_sq, color = "Inversely proportional"),
            linewidth = 1, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(kappa/phi),
    y = "Mean detection delay",
    title = "\u03B8 = 1, \u03b5 = 0.1",
    color = NULL
  ) +
  scale_color_manual(
    values = c(
      "Inverse-square (slope = -2.63)" = "pink",
      "Inversely proportional" = "#2C7BB6"
    )
  ) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    plot.title = element_text(hjust = 0.5)
  )

# R3-R4 plot
# log-log
long_data <- data1 %>%
  filter(snr >=0.51)%>%
  filter(snr <=1e5)%>%
  filter(stoppingT > 600)

result <- long_data %>%
  group_by(snr) %>%
  summarise(meanT = mean(stoppingT), sdT=sd(stoppingT))

result$meanT<-result$meanT-600

#mod1<-lm(meanT~I(snr^(-2.63))+0, data=result[result$snr<=0.46,], weights = 1/sdT)
mod2<-lm(meanT~I(snr^(-1))+0, data=result[result$snr<=1.2,], weights = 1/sdT)
mod3<-lm(meanT~1, data=result[result$snr>=2,])

# prediction grid
grid <- data.frame(snr = seq(0.5,1,length.out = 10000))
grid2 <- data.frame(snr = seq(2,100000,length.out = 1000))

grid$fit_inv <- predict(mod2, newdata = grid)
grid2$fit_inv_sq <- predict(mod3, newdata = grid2)

ggplot(result, aes(x = snr, y = meanT)) +
  theme_bw(base_size = 14) +
  annotation_logticks() +
  geom_point(alpha = 0.7, size = 1.8) +
  # Using the new labels here
  geom_line(data = grid, aes(y = fit_inv, color = "Inversely proportional"),
            linewidth = 1, linetype = "solid") +
  geom_line(data = grid2, aes(y = fit_inv_sq, color = "Constant"),
            linewidth = 1, linetype = "dashed") +
  scale_x_log10(breaks = 10^(0:5)) +
  scale_y_log10(limits = c(20, 75)) +
  labs(
    x = expression(kappa/phi),
    y = "Mean detection delay",
    title = "\u03B8 = 1, \u03b5 = 0.1",
    color = NULL
  ) +
  scale_color_manual(
    # BREAKS controls the order in the legend
    breaks = c("Inversely proportional", "Constant"),
    values = c(
      "Inversely proportional" = "#2C7BB6",
      "Constant" = "#D55E00"
    )
  ) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    plot.title = element_text(hjust = 0.5)
  )
