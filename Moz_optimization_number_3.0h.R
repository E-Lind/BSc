# Pakete laden
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

number <- c(20000,	30000,	50000, 80000, 100000, 120000, 150000)
avgFI <- c(18690.125,
           21646.25,
           36128.25,
           43698.875,
           47359.75,
           53083.25,
           65200.875
           
)
SD <- c(1626.086186,
        2874.370611,
        4811.85944,
        2905.881265,
        2503.465366,
        5664.708754,
        5348.48544
        
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
    y = "Moz-FI (p52) after 180 minutes"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15)
  )