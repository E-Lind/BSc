# Pakete laden
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

time <- c(60, 120, 180)
avgFI <- c(27765.33333, 58044, 86063)
SD <- c(4902.236938, 5446.809555, 5316.241279)


# Data Frame erstellen
df <- data.frame(
  time = time,
  avgFI = avgFI,
  SD = SD
)

# Plot
ggplot(df, aes(x = time, y = avgFI)) +
  geom_point(color = "darkred", size = 4) +
  geom_errorbar(aes(ymin = avgFI - SD, ymax = avgFI + SD),
                width = 60, color = "black") +
  labs(
    x = "incubation time [min]",
    y = "Gaoua-FI with 100.000 cells per well"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15)
  )