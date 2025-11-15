library(ggplot2)

# === Without shading ====
nu <- 2
Delta <- 100

eps <- seq(0, 0.1, length.out = 500)

boundary1 <- 5*pmax(eps^(1 - 1/nu),rep(Delta^(-1/2), length(eps)))
boundary2 <- rep(2.5,length(eps))
boundary3 <- eps^(-1/nu)

df <- data.frame(
  eps = eps,
  b1 = boundary1,
  b2 = boundary2,
  b3 = boundary3
)

ggplot(df, aes(y = eps)) +
  geom_line(aes(x = b1), color = "blue", size = 1) +
  geom_line(aes(x = b2), color = "red", size = 1) +
  geom_line(aes(x = b3), color = "purple", size = 1) +
  labs(
    y = expression(epsilon),
    x = expression(kappa)
  ) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 0.1), expand = FALSE)+
  theme_bw()+
  theme(text=element_text(size=20), axis.text.x = element_blank(), axis.text.y = element_blank())

# ==== With shading =====
eps <- seq(0, 0.1, length.out = 500)
Delta <- 100

# for v
nu <- 2
boundary1 <- 5*pmax(eps^(1 - 1/nu),rep(Delta^(-1/2), length(eps)))
boundary2 <- rep(2.5,length(eps))
boundary3 <- eps^(-1/nu)

# For theta
boundary1 <- 5*pmax(eps*log(1/eps),rep(Delta^(-1/2), length(eps)))
boundary2 <- rep(2.5,length(eps))
boundary3 <- 1.5*log(1/eps)

# Construct shaded regions (note: x = kappa, y = eps)
region1 <- df %>% mutate(xmin = 0,  xmax = b1, region = "1")
region2 <- df %>% mutate(xmin = b1, xmax = b2, region = "2")
region3 <- df %>% mutate(xmin = b2, xmax = b3, region = "3")
region4 <- df %>% mutate(xmin = b3, xmax = 10, region = "4")

regions <- bind_rows(region1, region2, region3, region4)

# Plot
ggplot(df, aes(y = eps)) +
  
  # Shading for each region
  geom_ribbon(data = regions,
              aes(xmin = xmin, xmax = xmax, fill = region),
              alpha = 0.25, show.legend = FALSE) +
  
  # Boundary lines
  geom_line(aes(x = b1), color = "red", size = 1) +
  geom_line(aes(x = b2), color = "blue", size = 1) +
  geom_line(aes(x = b3), color = "purple", size = 1) +
  
  # Region labels 
  annotate("text", x = 0.8,  y = 0.07, label = "R1", size = 7) +
  annotate("text", x = 2, y = 0.05, label = "R2", size = 7) +
  annotate("text", x = 4,  y = 0.03, label = "R3", size = 7) +
  annotate("text", x = 7,  y = 0.05,  label = "R4", size = 7) +
  
  labs(
    y = expression(epsilon),
    x = expression(kappa)
  ) +
  
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 0.1), expand = FALSE) +
  
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )