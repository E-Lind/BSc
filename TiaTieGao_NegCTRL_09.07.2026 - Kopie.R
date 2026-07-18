# ==================================================
# Libraries
# ==================================================

library(dplyr)
library(ggplot2)
library(emmeans)
library(purrr)


# ==================================================
# Daten
# ==================================================

datasets <- list(
  
  Tiassale = list(
    vec_list = list(
      c95013967 = c(76.91441335, 98.47938494, 110.2783483, 75.27517769),
      c48344124 = c(88.97333245, 93.81784798, 98.51345681),
      c24530475 = c(104.8647063, 102.3850311, 112.9044063),
      ctrl = c(98.65731583, 101.3426842)
    ),
    neg_ctrl = c(-4.263400815, -1.942727706, -0.670711124, -0.309801658)
  ),
  
  
  Tiefora = list(
    vec_list = list(
      c95013967 = c(96.70916447, 91.8775875, 79.58772425, 81.84591017),
      c48344124 = c(88.9733095, 85.50291682, 79.6237925),
      c24530475 = c(86.32464559, 82.87150295, 87.67014804),
      ctrl = c(101.8418329, 106.1135679, 102.3687429, 89.67585623)
    ),
    neg_ctrl = c(0.642171622, -0.112125204, -0.907194831, -0.093306988)
  ),
  
  
  "N´Gousso" = list(
    vec_list = list(
      c95013967 = c(88.12204322, 102.6195854, 101.3302988, 88.13013014),
      c48344124 = c(69.74508864, 73.75273656, 65.21410128),
      c24530475 = c(101.3996153, 92.64840948),
      ctrl = c(103.0724531, 96.92754695)
    ),
    neg_ctrl = c(0.485793009, 0.074515218, -0.351781145, -0.823133221)
  )
  
)


# ==================================================
# Alle Daten zusammenführen
# ==================================================

all_data <- imap_dfr(
  datasets,
  function(dat, assay_name){
    
    bind_rows(
      lapply(names(dat$vec_list), function(comp){
        
        data.frame(
          Assay = assay_name,
          Compound = comp,
          Wert = dat$vec_list[[comp]]
        )
        
      })
    )
    
  }
)


# ==================================================
# Reihenfolge der Compounds
# ==================================================

all_data$Compound <- factor(
  all_data$Compound,
  levels = c(
    "ctrl",
    "c95013967",
    "c48344124",
    "c24530475"
  )
)


# ==================================================
# Mittelwerte und SD
# ==================================================

summary_df <- all_data %>%
  group_by(
    Assay,
    Compound
  ) %>%
  summarise(
    mean = mean(Wert),
    sd = sd(Wert),
    .groups = "drop"
  )


# ==================================================
# Kontrolllinien + Z'
# ==================================================

z_df <- imap_dfr(
  datasets,
  function(dat, assay_name){
    
    pos_ctrl <- dat$vec_list$ctrl
    neg_ctrl <- dat$neg_ctrl
    
    data.frame(
      Assay = assay_name,
      ctrl_mean = mean(pos_ctrl),
      z_prime =
        1 -
        (3 * (sd(pos_ctrl) + sd(neg_ctrl))) /
        abs(mean(pos_ctrl) - mean(neg_ctrl))
    )
    
  }
)


# ==================================================
# One-way ANOVA + Dunnett
# ==================================================

stats_df <- imap_dfr(
  datasets,
  function(dat, assay_name){
    
    df <- bind_rows(
      lapply(names(dat$vec_list), function(comp){
        
        data.frame(
          Compound = comp,
          Wert = dat$vec_list[[comp]]
        )
        
      })
    )
    
    
    anova_model <- aov(
      Wert ~ Compound,
      data = df
    )
    
    
    emm <- emmeans(
      anova_model,
      "Compound"
    )
    
    
    dunnett <- contrast(
      emm,
      method = "trt.vs.ctrl",
      ref = "ctrl",
      adjust = "dunnett"
    )
    
    
    result <- summary(dunnett) %>%
      as.data.frame() %>%
      mutate(
        
        Compound = sub(
          " - ctrl",
          "",
          contrast
        ),
        
        sig = case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01 ~ "**",
          p.value < 0.05 ~ "*",
          TRUE ~ NA_character_
        ),
        
        Assay = assay_name
        
      ) %>%
      filter(!is.na(sig))
    
    
    if(nrow(result) > 0){
      
      result <- result %>%
        left_join(
          
          all_data %>%
            group_by(
              Assay,
              Compound
            ) %>%
            summarise(
              y_pos = max(Wert) + 5,
              .groups = "drop"
            ),
          
          by = c(
            "Assay",
            "Compound"
          )
          
        )
      
    }
    
    
    return(result)
    
  }
)


# ==================================================
# Plot
# ==================================================

final_plot <- ggplot(
  all_data,
  aes(
    x = Compound,
    y = Wert
  )
) +
  
  
  geom_jitter(
    width = 0.15,
    size = 3
  ) +
  
  
  geom_errorbar(
    data = summary_df,
    
    aes(
      x = Compound,
      ymin = mean - sd,
      ymax = mean + sd
    ),
    
    inherit.aes = FALSE,
    
    width = 0.2,
    linewidth = 1,
    colour = "red"
  ) +
  
  
  geom_point(
    data = summary_df,
    
    aes(
      x = Compound,
      y = mean
    ),
    
    inherit.aes = FALSE,
    
    colour = "red",
    size = 3
  ) +
  
  
  geom_hline(
    data = z_df,
    
    aes(
      yintercept = ctrl_mean
    ),
    
    inherit.aes = FALSE,
    
    colour = "blue",
    linetype = "dashed"
  ) +
  
  
  geom_text(
    data = stats_df,
    
    aes(
      x = Compound,
      y = y_pos,
      label = sig
    ),
    
    inherit.aes = FALSE,
    
    size = 6
  ) +
  
  
  facet_wrap(
    ~Assay,
    nrow = 1
  ) +
  
  
  labs(
    x = "Negative Control Compounds [5 µM]",
    y = "Normalised-FI after 48h"
  ) +
  
  
  # Nur die Anzeige der x-Achse ändern
  scale_x_discrete(
    labels = c(
      "ctrl" = "ctrl",
      "c95013967" = "N01",
      "c48344124" = "N03",
      "c24530475" = "N02"
    )
  ) +
  
  
  theme_bw() +
  
  theme(
    
    text = element_text(size = 14),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    strip.text = element_text(
      size = 14,
      face = "bold"
    )
    
  )


# ==================================================
# Plot anzeigen
# ==================================================

print(final_plot)