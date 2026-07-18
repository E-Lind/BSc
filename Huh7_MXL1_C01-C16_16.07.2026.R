
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c01 = c(112.7475871,                 90.43930969,                 100.2236663         ),
         c02 = c(105.9570241,                 113.3338427,                 107.3532017         ),
         c03 = c(37.63740367,                 35.7438813,                 37.31788045         ),
         c04 = c(70.37047328,                 105.1207068,                 98.93862723         ),
         c05 = c(27.83220863,                 39.34615825,                 41.37582529,                 53.76221221         ),
         c06 = c(71.69580001,                 98.55103168,                 58.68981589         ),
         c07 = c(41.48557457,                 59.16215456,                 57.82154628         ),
         c08 = c(32.78342925,                 30.7134745,                 32.54587068         ),
         c09 = c(52.98702111,                 83.83073605,                 91.9007811         ),
         c10 = c(37.4262405,                 39.81432923,                 38.38758861         ),
         c11 = c(81.20786722),
         c12 = c(14.37472432,                 101.2364159,                 105.4791285         ),
         c13 = c(112.680904,                 89.30430766,                 115.8650223         ), 
         c14 = c(67.24609192,                 53.50103671,                 57.85072014,                 75.29946619         ),
         c15 = c(113.4116397,                 89.80026326,                 115.6024576         ),   
         c16 = c(100.822425,                 97.26599266,                 108.6618576         ), 
         ctrl = c(84.82542571,                  115.4524206,                  99.72215373         )
         )
      
       neg_ctrl <- c(-4.513960039,                     3.59915118,                     3.607486568,                     4.146508341       )
       pos_ctrl <- c(84.82542571,                     115.4524206,                     99.72215373       ) 
       
       
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
         labs(x = "MXL1 Compounds C01 - C16 [5µM]", y = "Huh7 Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14),axis.text.x = element_text(angle = 45, hjust = 1))