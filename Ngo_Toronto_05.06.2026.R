
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
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
      
       neg_ctrl <- c(-2.921238826,
                     -1.595136211,
                     -0.025420497,
                     0.872770411
       )
       pos_ctrl <- c(128.5821718,
                     88.01953142,
                     90.09977545,
                     93.29852137
       ) 
       
       # ---------------------------
       # 2. DataFrame
       # ---------------------------
       df <- bind_rows(lapply(names(vec_list), function(name) {
         data.frame(
           Compound = name,
           Wert = vec_list[[name]]
         )
       }))
       
       # Mittelwerte berechnen (ohne Kontrolle)
       mittelwerte <- df %>%
         filter(Compound != "ctrl") %>%
         group_by(Compound) %>%
         summarise(mean = mean(Wert), .groups = "drop") %>%
         arrange(mean)
       
       # Kontrolle immer links, Rest nach Mittelwert
       levels_order <- c("ctrl", mittelwerte$Compound)
       
       df$Compound <- factor(df$Compound, levels = levels_order)
       
       ctrl_mean <- mean(vec_list$ctrl)
       # ---------------------------
       # One-way ANOVA
       # ---------------------------
       anova_model <- aov(Wert ~ Compound, data = df)
       
       # ANOVA-Tabelle
       summary(anova_model)
       
       # ---------------------------
       # Post-hoc: Dunnett-Test
       # ---------------------------
       emm <- emmeans(anova_model, "Compound")
       
       dunnett <- contrast(
         emm,
         method = "trt.vs.ctrl",
         ref = "ctrl",
         adjust = "dunnett"
       )
       
       stats <- summary(dunnett) %>%
         as.data.frame() %>%
         mutate(
           Compound = sub(" - ctrl", "", contrast),
           sig = case_when(
             p.value < 0.001 ~ "***",
             p.value < 0.01  ~ "**",
             p.value < 0.05  ~ "*",
             TRUE ~ NA_character_
           )
         ) %>%
         filter(!is.na(sig))
       
       stats$Compound <- factor(stats$Compound,
                                levels = levels(df$Compound))
       
       stats$y_pos <- sapply(as.character(stats$Compound), function(x) {
         max(df$Wert[df$Compound == x]) + 5
       })
       
       # ---------------------------
       # Z'-Faktor
       # ---------------------------
       z_prime <- 1 - (3 * (sd(pos_ctrl) + sd(neg_ctrl))) /
         abs(mean(pos_ctrl) - mean(neg_ctrl))
       
       # ---------------------------
       # Plot
       # ---------------------------
       ggplot(df, aes(x = Compound, y = Wert)) +
         geom_jitter(width = 0.15, size = 3) +
         stat_summary(fun = mean,
                      geom = "crossbar",
                      width = 0.5,
                      colour = "red") +
         stat_summary(fun.data = mean_sdl,
                      fun.args = list(mult = 1),
                      geom = "errorbar",
                      width = 0.25,
                      colour = "red") +
         geom_text(data = stats,
                   aes(x = Compound,
                       y = y_pos,
                       label = sig),
                   size = 5) +
         geom_hline(yintercept = ctrl_mean,
                    linetype = "dashed",
                    colour = "blue") +
         labs(x = "Potential P450 Compounds [5µM]", y = "Ngousso Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))