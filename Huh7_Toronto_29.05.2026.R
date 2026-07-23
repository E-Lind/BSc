
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c01 = c(99.20510494,
                 103.1746886,
                 108.0805363,
                 115.635219
         ),
         c02 = c(80.20099839,
                 95.81078073,
                 123.7798361,
                 110.6163739
         ),
         c03 = c(43.69697133,
                 54.61222573,
                 50.13195258,
                 51.85332841
         ),
         c04 = c(45.45943712,
                 46.79632832,
                 53.64807917,
                 57.26986565
         ),
         c05 = c(35.88107392,
                 33.34817137,
                 45.96865912,
                 33.36431386
         ),
         c06 = c(36.90098543,
                 66.76751643,
                 45.61792696,
                 53.73466158
         ),
         c07 = c(74.98697595,
                 72.09600376,
                 68.08239515
         ),
         c08 = c(78.42385761,
                 74.18425423,
                 65.06228309,
                 71.66309168
         ),
         c09 = c(98.60049553,
                 101.0556206,
                 87.17308412,
                 96.96129962
         ),
         c10 = c(48.4663417,
                 48.57640409,
                 58.27950467,
                 64.14509648
         ),
         c11 = c(100.5669436,
                 102.6904141,
                 85.48839576
         ),
         c12 = c(73.11004527,
                 81.47038466,
                 96.20260285,
                 85.00852372
         ),
         c13 = c(118.0404492,
                 123.5655813,
                 
                 123.0710343
         ), 
         c14 = c(78.31526272,
                 83.44270274,
                 79.78276129,
                 75.71485524
         ),
         ctrl = c(82.2452239,
                  99.18749496,
                  118.5672811
         )
         )
      
       pos_ctrl <- c(-11.48207573,
                     1.895641285,
                     4.295001455,
                     2.972785239
       )
       neg_ctrl <- c(82.2452239,
                     99.18749496,
                     118.5672811
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
         labs(x = "Potential P450 Compounds [5µM]", y = "Huh7 Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
