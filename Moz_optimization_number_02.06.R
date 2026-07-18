# Pakete laden
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

number <- c(15000,	20000,	25000, 30000)
avgFI <- c(26978.75, 40327.5,  47564,  62277.5)
SD <- c(12716.44173,  10863.00373,  14325.4516, 16698.08857)

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
    y = "Moz-FI after 3 hours"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15)
  )