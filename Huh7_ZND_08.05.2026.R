
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c31 = c(20.54679872,
                 21.51923958,
                 15.15851982,
                 14.88751171
         ),
         c32 = c(33.25229659,
                 47.1853018,
                 33.88199191,
                 50.99535699
         ),
         c38 = c(14.81577427,
                 21.68662695,
                 30.66177789,
                 25.74377777
         ),
         c74 = c(-7.725723851,
                 -9.288005898,
                 8.518821115,
                 6.581910209
         ),
         ctrl = c(76.03969472,
                  100.1673874,
                  90.74587012
         )
         )
       pos_ctrl <- c(-11.30263237,
                     -10.06715422,
                     -2.861526812,
                     -1.355040552
       )
       neg_ctrl <- c(76.03969472,
                     100.1673874,
                     133.0470478,
                     90.74587012
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
         labs(x = "ZND Compounds [5µM]", y = "Huh7 Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))
