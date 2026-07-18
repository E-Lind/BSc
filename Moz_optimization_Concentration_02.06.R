# Pakete laden
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Daten
number <- c(20, 25, 30)
avgFI <- c(62277.5,	75378.75,	82994.75)
SD <- c(16698.08857,	3561.333833,	5189.567387)

# Data Frame erstellen
df <- data.frame(
  number = factor(number, levels = c(20, 25, 30)),
  avgFI = avgFI,
  SD = SD
)

# Plot
ggplot(df, aes(x = number, y = avgFI)) +
  geom_point(color = "blue", size = 4) +
  geom_errorbar(
    aes(ymin = avgFI - SD, ymax = avgFI + SD),
    width = 0.2,
    color = "black"
  ) +
  labs(
    x = "Resazurin concentration [µM]",
    y = "Moz-FI after 3 hours"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15)
  )