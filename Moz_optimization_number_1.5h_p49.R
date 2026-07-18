# Pakete laden
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

number <- c(20000,	30000,	40000, 50000, 60000, 70000, 80000)
avgFI <- c(32384.5,
           38974,
           44748.75,
           52182.75,
           61056.5,
           64102.75,
           74398.75
           
)
SD <- c(574.6088525,
        3946.160877,
        854.2549873,
        3929.537494,
        8094.55006,
        3176.178981,
        10623.13571
        
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
    y = "Moz-FI (p49) after 90 minutes"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15)
  )