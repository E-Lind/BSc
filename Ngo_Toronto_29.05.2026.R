
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c01 = c(93.8896475,
                 95.78017256,
                 97.88403418
         ),
         c02 = c(107.4138238,
                 97.94077263,
                 85.08157097,
                 78.10501152
         ),
         c03 = c(105.359892,
                 108.0651811,
                 107.250417
         ),
         c04 = c(85.42200166,
                 67.71733662,
                 76.0669665,
                 80.19979499
         ),
         c05 = c(81.9019484,
                 83.65857072,
                 80.10220486
         ),
         c06 = c(77.6125218,
                 73.20054015,
                 82.39670766,
                 77.61479133
         ),
         c07 = c(92.69360104,
                 104.7698121,
                 105.4166304,
                 94.03716747
         ),
         c08 = c(112.7358901,
                 116.7189291,
                 114.0794565
         ),
         c09 = c(119.0996743,
                 85.13603988,
                 50.70488064,
                 86.71790779
         ),
         c10 = c(81.73173306,
                 80.96235971,
                 74.51914166,
                 85.98257751
         ),
         c11 = c(93.93276872,
                 85.89633507,
                 86.61123951,
                 101.0727349
         ),
         c12 = c(87.93664963,
                 84.96128547,
                 83.20012407
         ),
         c13 = c(100.811738,
                 100.7050698,
                 98.36063713
         ), 
         c14 = c(39.03491684,
                 50.07394911,
                 59.76260634
         ),
         ctrl = c(71.34178862,
                  111.2175692,
                  117.4406421,
                  147.8025199
         )
         )
      
       pos_ctrl <- c(1.224415689,
                     2.908412799,
                     2.095918236,
                     0.908949923
       )
       neg_ctrl <- c(71.34178862,
                     111.2175692,
                     117.4406421,
                     147.8025199
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
