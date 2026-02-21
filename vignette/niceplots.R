library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# Read CSV with no header
data1 <- read.csv("data/Latest_format/locations_v2e10R2R3.csv")

# ==== Data treatment ====
data2 <- read.csv("data/Latest_format/rawdata/locations_th1e10_R1R2_tunedup.csv", header = FALSE)
data_bars<-data2
data_bars[data_bars == -1] <- 2401

detected_long <- data_bars %>%
  mutate(snr = seq(from = 0, to = 0.1, by = 0.01)) %>%
  pivot_longer(-snr, values_to = "stoppingT")
colnames(detected_long)[2]<-"trial"

total<-as.data.frame(detected_long)
# total<-rbind(data1,as.data.frame(detected_long))
# total<-total[order(total$snr),]

write.csv(total,"data/Latest_format/locations_th1e10.csv", row.names = FALSE)
# write.csv(detected_long,
#         "data/Latest_format/locations_th1e10.csv", row.names = FALSE)

# ==== Proportion plot ====
# Function to compute proportions for each row
get_proportions <- function(x) {
  total <- length(x)
  c(
    not_detected = sum(x >= 2400) / total,
    false_alarm = sum(x >= 0 & x <= 600) / total,
    detected = sum(x > 600 & x < 2400) / total
  )
}

# Apply row-wise
# Group by snr and compute proportions
props <- data1 %>%
  filter(snr<=0.1) %>%
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
    y = "Proportion",
    title = "Detection performance (v = 2, ε = 0.1)"
  ) +
  scale_fill_manual(
    values = c(
      not_detected = "#999999",
      detected     = "#0072B2",
      false_alarm  = "#D55E00"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )
# ==== Boxplot ====
detected_long <- data1 %>%
  filter(stoppingT > 600, snr>=0.08) %>%
  mutate(snr = factor(snr))

ggplot(detected_long, aes(x = snr, y = stoppingT-600)) +
  geom_boxplot(fill = "steelblue", outlier.color = "red") +
  labs(x = "Signal size", y = "Detection delay", title = "Detection delay against signal size (v=2, eps=0.05)") +
  theme_minimal()

# ====== Mean plot =====
detected_long <- data1 %>%
  filter(snr >=0.08)%>%
  filter(stoppingT > 600)

result <- detected_long %>%
  group_by(snr) %>%
  summarise(meanT = mean(stoppingT), sdT=sd(stoppingT))

result$meanT<-result$meanT-600

mod<-lm(meanT~I((snr)^(-2))+0, data=result, weights = 1/sdT)
summary(mod)

mod2<-lm(meanT~I((log(snr*4.582576*4))^(-1))+0, data=result, weights = 1/sdT)
#mod2<-lm(meanT~I((log(snr*6))^(-1))+0, data=result[result$snr>=0.2,])
summary(mod2)

# dense grid for smooth curves
grid <- data.frame(
  snr = seq(min(result$snr),
            max(result$snr),
            length.out = 400)
)

grid$fit_inv_sq  <- predict(mod,  newdata = grid)
grid$fit_log_inv <- predict(mod2, newdata = grid)

ggplot(result, aes(x = snr, y = meanT)) +
  geom_point(aes(colour = "Observed"),
             size = 2.2) +
  geom_line(data = grid,
            aes(y = fit_inv_sq,
                colour = "Inverse-square"),
            linewidth = 1) +
  geom_line(data = grid,
            aes(y = fit_log_inv,
                colour = "Inverse-log"),
            linewidth = 1,
            linetype = "22") +
  labs(
    x = expression(kappa/phi),
    y = "Mean detection delay",
    colour = NULL
  ) +
  scale_colour_manual(
    values = c(
      "Observed"       = "black",
      "Inverse-square" = "#1f78b4",
      "Inverse-log"    = "#b15928"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(colour = "grey93", linewidth = 0.2),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    plot.title = element_text(face = "bold")
  )


# prediction grid
grid <- data.frame(
  snr = seq(min(result$snr[result$snr > 0]),
            max(result$snr),
            length.out = 400)
)

grid$fit_inv_sq  <- predict(mod,  newdata = grid)
grid$fit_log_inv <- predict(mod2, newdata = grid)

ggplot(result, aes(x = snr, y = meanT)) +
  geom_point(aes(colour = "Observed"),
             size = 2.5) +
  geom_line(data = grid,
            aes(y = fit_inv_sq,
                colour = "Inverse-square"),
            linewidth = 1) +
  geom_line(data = grid,
            aes(y = fit_log_inv,
                colour = "Inverse-log"),
            linewidth = 1,
            linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(kappa/phi),
    y = "Mean detection delay",
    colour = NULL,
    title = "Log–log comparison of scaling models"
  ) +
  scale_colour_manual(
    values = c(
      "Observed"       = "black",
      "Inverse-square" = "#0072B2",
      "Inverse-log"    = "#D55E00"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )
