
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c01 = c(2.083962801,
                 1.761489504,
                 2.488307558,
                 2.269426719
         ),
         c02 = c(45.75152569,
                 1.200085213,
                 4.11236997,
                 42.21267365
         ),
         c03 = c(55.845107,
                 62.46499926,
                 59.52096842,
                 58.10409093
         ),
         c04 = c(
           34.85927995,
           44.41484881,
           42.63706856,
           36.42820444
         ),
         c05 = c(41.72645743,
                 50.97626145,
                 47.38227149,
                 41.10657353
         ),
         c06 = c(12.30285756,
                 17.16000718,
                 23.68131952,
                 29.14832797
         ),
         c07 = c(53.64794437,
                 60.08404356,
                 57.73149225,
                 49.66631811
         ),
         c08 = c(35.17841155,
                 51.51260305,
                 48.09238108,
                 41.97040096
         ),
         c09 = c(52.07400734,
                 57.62288725,
                 49.4557915,
                 51.48921121
         ),
         c10 = c(46.52345659,
                 42.73397763,
                 39.46246096
         ),
         c11 = c(52.51343987,
                 51.05312037,
                 52.97459339
         ),
         c12 = c(49.3288072,
                 47.94367578,
                 46.98461012
         ),
         c13 = c(32.94783201,
                 30.55685123,
                 39.36555189
         ), 
         c14 = c(31.65459712,
                 42.72061086,
                 40.25110032,
                 34.92778464
         ),
         ctrl = c(135.0666041,
                  130.9262475,
                  117.3890523
         )
         )
      
       pos_ctrl <- c(53.38896323,
                     45.62454139,
                     43.55603391,
                     63.01303678
       )
       neg_ctrl <- c(34.00714844,
                     135.0666041,
                     130.9262475,
                     117.3890523
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
         labs(x = "Potential P450 Compounds  [5µM]", y = "Tiassale Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
