library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# Read CSV with no header
data1 <- read.csv("data/Latest_format/locations_v2e10all.csv")

# ==== Data treatment ====
data2 <- read.csv("data/Latest_format/locations_v2e10new2.csv", header = FALSE)
data_bars<-data2
data_bars[data_bars == -1] <- 2401

detected_long <- data_bars %>%
  mutate(snr = c(2:50)/5) %>%
  pivot_longer(-snr, values_to = "stoppingT")
colnames(detected_long)[2]<-"trial"

write.csv(rbind(data1,as.data.frame(detected_long)),
          "data/Latest_format/locations_v2e10all.csv", row.names = FALSE)

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
props <- detected_long %>%
  filter(snr<=0.08) %>%
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
                           levels = c("detected","not_detected","false_alarm")))

ggplot(props_longer, aes(x = factor(snr), y = proportion, fill = category)) +
  geom_bar(stat = "identity") +
  labs(
    x = "SNR",
    y = "Proportion",
    fill = "Category",
    title = "Proportion of Outcomes by SNR"
  ) +
  scale_fill_manual(
    values = c(
      false_alarm  = "firebrick2",
      not_detected = "gold",
      detected     = "forestgreen"
    )
  ) +
  theme_minimal(base_size = 14)

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
  summarise(meanT = mean(stoppingT))

mod<-lm(meanT~I((snr)^(-2)), data=result[result$snr<=0.3,])
summary(mod)

mod2<-lm(meanT~I((log(snr*4.582576))^(-1)), data=result[result$snr>=0.31,])
summary(mod2)

plot(result$snr, result$meanT,
     xlab = "Signal-to-noise ratio",
     ylab = "Detection time",
     pch = 19, xlim=c(0.1,0.6), ylim = c(630,850))

plot(result$snr, result$meanT,
     xlab = "Signal-to-noise ratio",
     ylab = "Detection time",
     pch = 19, xlim=c(1,10), ylim = c(624,630))

# add fitted regression line
curve(predict(mod, newdata = data.frame(snr = x)), from = 0,
      to = 0.3,add = TRUE, lwd = 2, col = "blue")

curve(predict(mod2, newdata = data.frame(snr = x)), from = 0.3,
      to = 10,add = TRUE, lwd = 2, col = "green")
