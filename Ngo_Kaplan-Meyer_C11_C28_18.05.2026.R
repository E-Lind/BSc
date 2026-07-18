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

day <- c(0, 1, 2, 3, 4, 5,6 ,7 ,8 ,9 ,10, 11, 12, 13, 14, 15)
avg_Control <- c(0, 25.71428571,	34.28571429,	48.57142857,	51.42857143,	54.28571429,	54.28571429,	57.14285714,	62.85714286,	65.71428571,	68.57142857,	77.14285714, NA, NA,	94.28571429,	100)
avg_C11 <- c(0, 8.823529412,	14.70588235,	20.58823529,	38.23529412,	38.23529412,	38.23529412,	64.70588235,	76.47058824,	85.29411765,	85.29411765,	88.23529412, NA, NA,	94.11764706,	100)
avg_C28 <- c(0, 35.29411765,	35.29411765,	47.05882353,	47.05882353,	50,	55.88235294,	79.41176471,	79.41176471,	82.35294118,	82.35294118,	82.35294118, NA, NA,	100,	100)


#Number at risk
CTRL <- c(35, 26, 22, 18, 17, 16, 16, 15, 13, 12,11, 8, NA, NA, 2, 0)	
C11 <- c(34, 31, 29, 28, 21, 21, 21, 12, 8, 5, 5, 4, NA, NA, 2, 0)	
C28 <- c(34, 22, 22, 18, 18, 17, 15, 7, 7, 6, 6, 6, NA, NA, 0, 0)



mortality <- list(
  Control = avg_Control,
  C11 = avg_C11,
  C28 = avg_C28
)



# ==========================
# Funktion: Mortalität -> Survivaldaten
# ==========================

mortality <- data.frame(
  Day = day,
  Control = avg_Control,
  C11 = avg_C11,
  C28 = avg_C28
)


# NA ersetzen (letzter bekannter Wert)
mortality <- mortality %>%
  fill(
    Control,
    C11,
    C28,
    .direction = "downup"
  )


# Überleben (%)
survival_data <- mortality %>%
  mutate(
    Control = 100 - Control,
    C11 = 100 - C11,
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
      C28 = "#009E73"
    ),
    labels = c(
      Control = "Control",
      C11 = "C11",
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
    breaks = 0:15
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
      "C28"
    ),
    each = 18
  ),
  
  Mortality = c(
    avg_Control,
    avg_C11,
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