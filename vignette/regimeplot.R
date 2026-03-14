library(ggplot2)
library(dplyr)
library(tidyr)

# Parameters
nu <- 2
Delta <- 100
eps <- seq(0.001, 0.1, length.out = 999) # Higher resolution for smoother curves

# Define Boundaries
df <- data.frame(eps = eps) %>%
  mutate(
    b1 = 5 * pmax(eps * log(1/eps), Delta^(-1/2)),
    b2 = 2.5,
    b3 = 1.5 * log(1/eps),
    xmax = 10.4 # Plot limit
  )

# Plotting
ggplot(df, aes(y = eps)) +
  # Use geom_ribbon with explicit boundaries to avoid gaps
  geom_ribbon(aes(xmin = 0, xmax = b1, fill = "Regime 1"), alpha = 0.15) +
  geom_ribbon(aes(xmin = b1, xmax = b2, fill = "Regime 2"), alpha = 0.15) +
  geom_ribbon(aes(xmin = b2, xmax = b3, fill = "Regime 3"), alpha = 0.15) +
  geom_ribbon(aes(xmin = b3, xmax = xmax, fill = "Regime 4"), alpha = 0.15) +

  # Boundary lines - using a consistent color or subtle grey helps the labels pop
  geom_line(aes(x = b1), color = "grey30", size = 0.8, linetype = "solid") +
  geom_line(aes(x = b2), color = "grey30", size = 0.8, linetype = "dashed") +
  geom_line(aes(x = b3), color = "grey30", size = 0.8, linetype = "dotdash") +

  # Direct Labels (instead of legend for a cleaner look)
  annotate("text", x = c(0.4, 1.6, 3.5, 7.5), y = 0.05,
           label = paste0("R", 1:4), size = 6, fontface = "bold") +

  # Styling
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(
    y = expression(textstyle(Contamination~Level)~~(epsilon)),
    x = expression(textstyle(Signal~Strength)~~(kappa))
  ) +
  theme_classic(base_size = 18) + # Classic is often preferred over BW for theory plots
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
