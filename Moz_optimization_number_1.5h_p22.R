# Pakete laden
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

number <- c(20000,	30000,	40000, 50000, 60000, 70000, 80000)
avgFI <- c(31654,
           37123.5,
           46193.5,
           51732.5,
           62661.25,
           68960.5,
           77281.75
)
SD <- c(3937.554656,
        9027.627688,
        4053.707316,
        5517.501639,
        2158.158687,
        4550.270798,
        5044.867978
)

# Data Frame erstellen
df <- data.frame(
  number = number,
  avgFI = avgFI,
  SD = SD
)

# Plot
ggplot(df, aes(x = number, y = avgFI)) +
  geom_point(color = "darkred", size = 4) +
  geom_errorbar(aes(ymin = avgFI - SD, ymax = avgFI + SD),
                width = 2000, color = "black") +
  labs(
    x = "Cell Number per Well",
    y = "Moz-FI (p22) after 90 minutes"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15)
  )