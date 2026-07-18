library(tidyverse)
library(survival)
library(survminer)
library(zoo)
library(cowplot)


## Gruppengrößen
n <- c(Ctrl = 85,
       c11  = 78,
       c25  = 88,
       c28  = 60)

day <- c(0, 1, 2, 3, 4, 5,6 ,7 ,8 ,9 ,10, 11, 12, 13, 14, 15, 16, 17)
avg_Control <- c(0, 9.333333333,	14,	14,	19.55555556,	23.33333333,	37.33333333,	51.11111111,	64.44444444,	72.66666667,	80.66666667,	89.55555556,	94.44444444,	94.44444444,	94.44444444,	100,	100,	100)
avg_C11 <- c(0, 6.45410628,	11.57487923,	12.68599034,	13.79710145,	22.24154589,	33.37198068,	42.83091787,	54.05797101,	64.4057971,	79.64251208,	97.55555556,	100,	100,	NA,	100,	100,	NA)
avg_C25 <- c(0, 8.968253968,	17.93650794,	25.87301587,	32.53968254,	35.87301587,	40.3968254,	43.80952381,	70.47619048,	72.46031746,	88.88888889,	94.44444444,	98.88888889,	100,	100,	NA,	100,	100)
avg_C28 <- c(0, 11.66666667,	20,	30,	43.33333333,	58.33333333,	61.66666667,	68.33333333,	81.66666667,	85,	85,	91.66666667,	98.33333333,	98.33333333,	100,	NA,	100,	100)
sd_Control <- c(0, 7.23161782,	9.863888498,	9.863888498,	11.46976702,	12.58305739,	19.23730943,	26.87419249,	32.48931448,	36.50824361,	41.35975858,	45.32475409,	47.87135539,	54.2968521,	65.99663291,	70.71067812,	70.71067812,	70.71067812)
sd_C11 <- c(0, 2.355031586,	6.866649741,	7.531008967,	8.584213468,	9.630926404,	1.391706866,	5.670822697,	5.259491974,	4.061282293,	3.856960781,	2.143033502,	0,	0,	NA,	NA,	NA,	NA)
sd_C25 <- c(0, 6.844289099,	11.48135957,	13.04391095,	15.89601508,	18.98999219,	19.59096079,	19.4539983,	3.433858358,	10.30799321,	11.70628195,	5.091750772,	1.924500897,	0,	NA,	NA,	NA,	NA)
sd_C28 <- c(0, 7.071067812,	9.428090416,	9.428090416,	4.714045208,	7.071067812,	7.071067812,	2.357022604,	2.357022604,	2.357022604,	11.78511302,	2.357022604,	2.357022604,	2.357022604,	NA,	NA,	NA,	NA)

#Number at risk
CTRL <- c(85,	77,	73,	73,	68,	65,	53,	41,	30,	23,	16,	9,	5)	
C11 <- c(78,	80,	75,	72,	66,	58,	52,	47,	35,	27,	18,	5,	0)	
C25 <- c(88,	78,	73,	67,	59,	50,	48,	44,	36,	34,	34,	30,	26,	26)
C28 <- c(60,	79,	76,	74,	73,	72,	67,	63,	46,	45,	29,	28,	27,	27)



mortality <- list(
  Control = avg_Control,
  C11 = avg_C11,
  C25 = avg_C25,
  C28 = avg_C28
)



# ==========================
# Funktion: Mortalität -> Survivaldaten
# ==========================

mortality <- data.frame(
  Day = day,
  Control = avg_Control,
  C11 = avg_C11,
  C25 = avg_C25,
  C28 = avg_C28
)


# NA ersetzen (letzter bekannter Wert)
mortality <- mortality %>%
  fill(
    Control,
    C11,
    C25,
    C28,
    .direction = "downup"
  )


# Überleben (%)
survival_data <- mortality %>%
  mutate(
    Control = 100 - Control,
    C11 = 100 - C11,
    C25 = 100 - C25,
    C28 = 100 - C28
  ) %>%
  pivot_longer(
    cols = -Day,
    names_to = "Group",
    values_to = "Survival"
  )


# ==========================================
# Kaplan-Meier-ähnliche Treppenkurve
# ==========================================

km_plot <- ggplot(
  survival_data,
  aes(
    x = Day,
    y = Survival/100,
    colour = Group
  )
)+
  
  geom_step(
    linewidth = 1.2
  )+
  
  scale_colour_manual(
    values = c(
      Control = "black",
      C11 = "#D55E00",
      C25 = "#0072B2",
      C28 = "#009E73"
    ),
    labels = c(
      Control = "Control",
      C11 = "C11",
      C25 = "C25",
      C28 = "C28"
    )
  )+
  
  scale_y_continuous(
    limits = c(0,1),
    breaks = seq(0,1,0.1)
  )+
  
  scale_x_continuous(
    breaks = 0:17
  )+
  
  labs(
    x = "Time [d]",
    y = "Survival probability",
    colour = NULL
  )+
  
  theme_classic(
    base_size = 14
  )+
  
  theme(
    legend.position = "top"
  )



# ==========================================
# Number at risk Tabelle
# ==========================================

# Zeitpunkte der Risk Table
risk_table <- bind_rows(
  
  data.frame(
    Day = 0:(length(CTRL)-1),
    Group = "Control",
    At_risk = CTRL
  ),
  
  data.frame(
    Day = 0:(length(C11)-1),
    Group = "C11",
    At_risk = C11
  ),
  
  data.frame(
    Day = 0:(length(C25)-1),
    Group = "C25",
    At_risk = C25
  ),
  
  data.frame(
    Day = 0:(length(C28)-1),
    Group = "C28",
    At_risk = C28
  )
)


# Risk Table Plot

risk_plot <- ggplot(
  risk_table,
  aes(
    x = Day,
    y = Group,
    label = At_risk
  )
)+
  
  geom_text(
    size = 4
  )+
  
  scale_x_continuous(
    breaks = 0:13
  )+
  
  labs(
    x = "Time [d]",
    y = NULL
  )+
  
  theme_classic(
    base_size = 14
  )+
  
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 12)
  )

# ==========================================
# Plot zusammenfügen
# ==========================================

final_plot <- plot_grid(
  km_plot,
  risk_plot,
  ncol = 1,
  rel_heights = c(3,1)
)

final_plot











# ==========================================
# ANOVA: Mortalität Tag 17
# Control vs C11, C25, C28
# ==========================================

library(tidyverse)
library(multcomp)


# ==========================================
# Daten zusammenstellen
# ==========================================

anova_data <- data.frame(
  
  Day = rep(
    0:17,
    4
  ),
  
  Group = rep(
    c(
      "Control",
      "C11",
      "C25",
      "C28"
    ),
    each = 18
  ),
  
  Mortality = c(
    avg_Control,
    avg_C11,
    avg_C25,
    avg_C28
  )
)



# ==========================================
# NA entfernen
# ==========================================

anova_data <- na.omit(anova_data)



# ==========================================
# Faktor definieren
# ==========================================

anova_data$Group <- factor(
  anova_data$Group,
  levels = c(
    "Control",
    "C11",
    "C25",
    "C28"
  )
)



# ==========================================
# One-way ANOVA
# ==========================================

anova_model <- aov(
  Mortality ~ Group,
  data = anova_data
)


anova_result <- summary(anova_model)


print(anova_result)



# ==========================================
# Dunnett Test
# Vergleich jede Gruppe gegen Control
# ==========================================

dunnett <- glht(
  anova_model,
  linfct = mcp(
    Group = "Dunnett"
  )
)


dunnett_summary <- summary(dunnett)


print(dunnett_summary)



# ==========================================
# Dunnett Ergebnisse + Signifikanzsterne
# ==========================================

dunnett_results <- data.frame(
  
  Comparison = names(
    coef(dunnett)
  ),
  
  Estimate = dunnett_summary$test$coefficients,
  
  p_value = dunnett_summary$test$pvalues
  
)



# Funktion für Signifikanzsterne

get_stars <- function(p){
  
  if(p < 0.001){
    
    return("***")
    
  } else if(p < 0.01){
    
    return("**")
    
  } else if(p < 0.05){
    
    return("*")
    
  } else {
    
    return("ns")
    
  }
  
}



# Sterne ergänzen

dunnett_results$Significance <- sapply(
  
  dunnett_results$p_value,
  
  get_stars
  
)



# ==========================================
# Ergebnis anzeigen
# ==========================================

print(dunnett_results)