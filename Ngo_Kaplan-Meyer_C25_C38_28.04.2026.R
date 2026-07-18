library(tidyverse)
library(survival)
library(survminer)
library(zoo)
library(cowplot)


## Gruppengrößen
n <- c(Ctrl = 80,
       c25  = 74,
       c38  = 31)

day <- c(0, 1, 2, 3, 4)
avg_Control <- c(0, 11.25, NA, 16.25, 26.25)
avg_C25 <- c(0, 13.333333333, NA, 16, NA)  
avg_C38 <- c(0, 21.875, NA, 40.3225806, NA)


#Number at risk
CTRL <- c(80, 71, NA, 67, 67)	
C25 <- c(74, 64, NA, 62, NA)	
C38 <- c(32, 25, NA, 14, NA)



mortality <- list(
  Control = avg_Control,
  C25 = avg_C25,
  C38 = avg_C38
)



# ==========================
# Funktion: Mortalität -> Survivaldaten
# ==========================

mortality <- data.frame(
  Day = day,
  Control = avg_Control,
  C25 = avg_C25,
  C38 = avg_C38
)


# NA ersetzen (letzter bekannter Wert)
mortality <- mortality %>%
  fill(
    Control,
    C25,
    C38,
    .direction = "downup"
  )


# Überleben (%)
survival_data <- mortality %>%
  mutate(
    Control = 100 - Control,
    C25 = 100 - C25,
    C38 = 100 - C38
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
      C25 = "#D55E00",
      C38 = "#009E73"
    ),
    labels = c(
      Control = "Control",
      C25 = "C25",
      C38 = "C38"
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
    Day = 0:(length(C25)-1),
    Group = "C25",
    At_risk = C25
  ),
  
  
  data.frame(
    Day = 0:(length(C38)-1),
    Group = "C38",
    At_risk = C38
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
    breaks = 0:4
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
      "C25",
      "C38"
    ),
    each = 18
  ),
  
  Mortality = c(
    avg_Control,
    avg_C25,
    avg_C38
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
    "C25",
    "C38"
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