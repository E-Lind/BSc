# Pakete laden
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

number <- c(20000,	30000,	40000)
avgFI <- c(63399.25,
           68796.75,
           78209.25
           
           
)
SD <- c(4609.457126,
        10071.77985,
        7846.977359
        
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
    y = "Moz-FI (p49) after 150 minutes"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15)
  )