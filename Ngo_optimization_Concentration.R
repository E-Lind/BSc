# Pakete laden
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Daten
number <- c(20, 25, 30)
avgFI <- c(65005.25, 85954.83333, 86524)
SD <- c(17608.18935, 3244.05076, 304.7630227)

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
    y = "Ngousso-FI after 3 hours"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15)
  )