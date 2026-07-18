library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(purrr)
library(tibble)
library(tidyverse)

Ngo <- list(
  c01 = c(109.1121891, 37.30775749, 74.38778969,  68.17247807),
  c02 = c(54.28229462, 59.67355845, 57.97674024),
  c03 = c(126.5930178, 69.62780155, 97.7153328,  92.45752659),
  c04 = c(50.28492141, 47.25352709, 53.15108249, 51.18734907),
  c05 = c(71.22717451,  76.4616786,  73.63788501),
  c06 = c(63.25361183, 67.8928526, 64.37211371),
  c07 = c(91.11871372,  80.08621785, 75.05719612),
  c08 = c(84.7360505, 86.12358598,      90.17391857),
  c09 = c(64.36152184,     57.2162437, 63.10956234,      67.33360166),
  c10 = c(82.20035589,    79.86802525,    73.10829132),
  c11 = c(83.58789137,     80.53955006,   77.1670974),
  c12 = c(68.48176079,     72.76087785,   62.25797568),
  c13 = c(67.06668644,  61.64152862,  62.81510825), 
  c14 = c(33.33792315, 35.73804177,  35.21268483,  35.79523789),
  ctrl = c(128.5821718,  88.01953142, 90.09977545, 93.29852137)
)

Tia <- list(
  c01 = c(20.48455768,     30.34251618,   26.65427643),
  c02 = c(43.17824948,      53.05886048,       45.14901745),
  c03 = c(103.1950329,      110.3655805,   98.70159958,      90.49315534),
  c04 = c(76.3621106,                 81.06353513,                 74.88969774,                 70.67015378),
  c05 = c(78.79416595,                 79.45932588,                 75.27890897),
  c06 = c(82.30942293,                 75.17800236,                 74.09480074),
  c07 = c(84.15251314,                92.80782953,              81.39714475,             75.94200959),
  c08 = c(29.22430614,                 50.72153378,                 73.95476707,                 76.06968734         ),
  c09 = c(48.9484604,                 56.99421846,                 59.43245177         ),
  c10 = c(61.52060091,                 62.52554842,                 57.35459923         ),
  c11 = c(81.14590788,                 70.58572171,                 68.46462348,                 78.10017556         ),
  c12 = c(72.9436416,                 70.0214684,                 70.38184917         ),
  c13 = c(47.53782711,                 50.57326284,                 50.96865202,                 52.6428781         ), 
  c14 = c(71.64421151,                 67.62236214,                 66.97367676,                 60.42092474         ),
  ctrl = c(95.50759631,                  100.4437832,                  94.35025922         )
)

Tie <- list(
  c01 = c(6.178199781,                 5.227269492,                 2.684812163,                 3.137500534         ),
  c02 = c(0.808575455,                 0.614972881,                 1.010719319,                 0.931000612         ),
  c03 = c(55.61518641,                 55.23082836,                 62.03823651,                 61.71651458         ),
  c04 = c(13.10803308,                 19.95245349,                 27.83891127,                 14.30666078         ),
  c05 = c(90.07075035,                 83.10675189,                 85.31040471,                 83.96372799         ),
  c06 = c(57.14692442,                 57.8501573,                 65.2099022,                 71.66427036         ),
  c07 = c(84.54453571,                 90.50350905,                 91.60248836         ),
  c08 = c(45.00405711,                 40.90139081,                 40.09850954,                 51.76021752         ),
  c09 = c(25.56977522,                 24.85230686,                 37.92048059,                 35.83071163         ),
  c10 = c(97.4760488,                 92.06087093,                 85.1509673,                 96.07527724         ),
  c11 = c(81.01413583,                 74.90711347,                 84.21711959         ),
  c12 = c(65.63411961,                 66.99218472,                 79.00123849,                 82.76794739         ),
  c13 = c(92.32849801,                 89.20238587,                 99.46901647,                 101.6726693         ), 
  c14 = c(79.99487523,                 81.57216678,                 91.89289222,                 90.3696955         ),
  ctrl = c(100.7900693,                  112.2837986,                  113.7842185         )
)



# ---------------------------
# 1. Alle Datensätze bündeln
# ---------------------------
all_data <- list(
  Ngo = Ngo,
  Tia = Tia,
  Tie = Tie
)

# ---------------------------
# 2. DataFrame erstellen
# ---------------------------
df <- imap_dfr(all_data, function(dataset, dataset_name) {
  imap_dfr(dataset, function(values, compound) {
    tibble(
      Dataset = dataset_name,
      Compound = compound,
      Wert = values
    )
  })
})

# ---------------------------
# 3. Reihenfolge nach Mittelwert
# ---------------------------
compound_order <- df %>%
  group_by(Compound) %>%
  summarise(mean_value = mean(Wert, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_value) %>%
  pull(Compound)

# Kontrolle immer links
compound_order <- c("ctrl", setdiff(compound_order, "ctrl"))

# Faktor-Reihenfolge setzen
df$Compound <- factor(df$Compound, levels = compound_order)

# ---------------------------
# 4. Hintergrundbereiche
# ---------------------------
bg <- data.frame(
  xmin = seq(0.5, length(compound_order) - 0.5, by = 2),
  xmax = seq(1.5, length(compound_order) + 0.5, by = 2),
  ymin = -Inf,
  ymax = Inf
)

# ---------------------------
# 5. Plot
# ---------------------------
ggplot(df, aes(x = Compound, y = Wert)) +
  
  # Hintergrund
  geom_rect(
    data = bg,
    inherit.aes = FALSE,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "grey95"
  ) +
  
  # Boxplots
  geom_boxplot(
    aes(fill = Dataset, colour = Dataset),
    position = position_dodge(width = 0.8),
    width = 0.65,
    alpha = 0.5,
    outlier.shape = NA
  ) +
  
  # Einzelwerte
  geom_jitter(
    aes(colour = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 2.5,
    alpha = 0.9
  ) +
  
  scale_fill_manual(
    name = NULL,
    values = c(
      Ngo = "orange",
      Tia = "blue",
      Tie = "forestgreen"
    ),
    labels = c(
      Ngo = "N´Gousso",
      Tia = "Tiassale",
      Tie = "Tiefora"
    )
  ) +
  
  scale_colour_manual(
    values = c(
      Ngo = "orange",
      Tia = "blue",
      Tie = "forestgreen"
    ),
    guide = "none"
  ) +
  
  labs(
    x = "Potential P450 Compounds",
    y = "Normalised FI"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "top"
  )