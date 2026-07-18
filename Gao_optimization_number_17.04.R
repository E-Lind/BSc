# Pakete laden
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

number <- c(30000,	50000,	80000, 100000, 130000, 150000)
avgFI <- c(41615.25,	55188.25,	73358.5, 86063,	77911,	80468)
SD <- c(11110.44195,	15841.3638,	18481.00612, 5316.241279,	1762.513073,	6519.298301)

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
    y = "Gaoua-FI after 3 hours"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15)
  )