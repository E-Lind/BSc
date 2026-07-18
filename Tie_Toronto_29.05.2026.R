
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c01 = c(0.134216671,
                 -0.961579712,
                 -1.079249122,
                 -0.395295675
         ),
         c02 = c(59.49843414,
                 62.32249998,
                 0.152602516,
                 -3.009762884
         ),
         c03 = c(71.60367472,
                 74.95725291,
                 73.09660536
         ),
         c04 = c(2.502313552,
                 -1.899257825,
                 -1.380776986,
                 -0.512965085
         ),
         c05 = c(100.7305309,
                 74.98299309,
                 93.5526969,
                 90.2358904
         ),
         c06 = c(29.88251445,
                 11.92322071,
                 10.05154165
         ),
         c07 = c(124.374728,
                 140.4917601,
                 104.1135265,
                 92.78416856
         ),
         c08 = c(4.969693998,
                 14.57813678,
                 15.39814548,
                 21.78171099
         ),
         c09 = c(27.26804724,
                 25.01761977,
                 38.50915309),
         c10 = c(27.81226826,
                 32.90147025,
                 36.69263157,
                 42.80040939
         ),
         c11 = c(106.3161507,
                 108.2834362,
                 109.3902641,
                 117.958068
         ),
         c12 = c(71.53380851,
                 73.44225925,
                 60.89575839,
                 83.95896279
         ),
         c13 = c(89.04081045,
                 89.94539404,
                 63.12412284,
                 120.9549608
         ), 
         c14 = c(38.14143618,
                 56.37284043,
                 54.43497233),
         ctrl = c(60.67880541,
                  126.0368085,
                  113.2843861,
                  175.8698037
         )
         )
      
       neg_ctrl <- c(3.912507891,
                     8.46116603,
                     5.250997432,
                     6.585809805
       )
       pos_ctrl <- c(60.67880541,
                     126.0368085,
                     113.2843861,
                     175.8698037
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
         labs(x = "Potential P450 Compounds [5µM]", y = "Tiefora Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))