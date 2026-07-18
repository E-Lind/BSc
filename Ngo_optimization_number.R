# Pakete laden
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

number <- c(15000,	20000,	25000, 30000, 35000)
avgFI <- c(40262.25, 51567, 56186, 65005.25, 63224.25)
SD <- c(7216.08003,  15593.91538,  14736.68914,  17608.18935, 15748.57269)

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
    y = "Ngousso-FI after 3 hours"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15)
  )