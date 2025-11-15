library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# Read CSV with no header
data <- read.csv("data/Latest_format/locations_v2new2.csv", header = FALSE)

# Function to compute proportions for each row
get_proportions <- function(x) {
  total <- length(x)
  c(
    not_detected = sum(x == -1) / total,
    false_alarm = sum(x >= 0 & x <= 600) / total,
    detected = sum(x > 600) / total
  )
}

# Apply row-wise
props <- t(apply(data, 1, get_proportions))
props <- as.data.frame(props)

desired_order <-2:50/5
rownames(props) <- desired_order
props[1,2]<-props[1,2]+props[1,3]
props[1,3]<-0

# Prepare proportions data for ggplot
props_long <- props %>%
  rownames_to_column("row_label") %>%
  pivot_longer(cols = -row_label, names_to = "category", values_to = "proportion") %>%
  mutate(
    row_label = factor(row_label, levels = desired_order),
    # define stacking order: zero_600 at bottom, above600 above it, neg1 on top
    category = factor(category, levels = c("detected","not_detected", "false_alarm"))
  )


ggplot(props_long, aes(x = row_label, y = proportion, fill = category)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.2, linetype = "dotted", color = "yellow", size = 1)+
  labs(x = "Signal size", y = "Proportion", title = "Proportion of Successful Detection against signal size (v=2, eps=0.05)") +
  scale_fill_manual(values = c("not_detected" = "orange", "false_alarm" = "red", "detected" = "steelblue")) +
  theme_minimal()

# Boxplot
detected_long <- data[-1,] %>%
  mutate(row_label = rownames(props)[-1]) %>%
  pivot_longer(-row_label, values_to = "value") %>%
  filter(value > 600) %>%
  mutate(row_label = factor(row_label, levels = desired_order[-1]))

ggplot(detected_long, aes(x = row_label, y = value-600)) +
  geom_boxplot(fill = "steelblue", outlier.color = "red") +
  labs(x = "Signal size", y = "Detection delay", title = "Detection delay against signal size (v=2, eps=0.05)") +
  theme_minimal()
