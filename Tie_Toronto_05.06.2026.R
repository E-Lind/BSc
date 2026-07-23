
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c01 = c(6.178199781,                 5.227269492,                 2.684812163,                 3.137500534         ),
         c02 = c(0.808575455,                 0.614972881,                 1.010719319,                 0.931000612         ),
         c03 = c(55.61518641,                 55.23082836,                 62.03823651,                 61.71651458         ),
         c04 = c(13.10803308,                 19.95245349,                 27.83891127,                 14.30666078         ),
         c05 = c(90.07075035,                 83.10675189,                 85.31040471,                 83.96372799         ),
         c06 = c(57.14692442,                 57.8501573,                 65.2099022,                 71.66427036         ),
         c07 = c(84.54453571,                 90.50350905,                 91.60248836         ),
         c08 = c(45.00405711,                 40.90139081,                 40.09850954,                 51.76021752         ),
         c09 = c(25.56977522,                 24.85230686,                 37.92048059,                 35.83071163         ),
         c10 = c(97.4760488,                 92.06087093,                 85.1509673,                 96.07527724         ),
         c11 = c(81.01413583,                 74.90711347,                 84.21711959         ),
         c12 = c(65.63411961,                 66.99218472,                 79.00123849,                 82.76794739         ),
         c13 = c(92.32849801,                 89.20238587,                 99.46901647,                 101.6726693         ), 
         c14 = c(79.99487523,                 81.57216678,                 91.89289222,                 90.3696955         ),
         ctrl = c(100.7900693,                  112.2837986,                  113.7842185         )
         )
       
       pos_ctrl <- c(-3.623642291,
                     2.463450396,
                     0.077583384,
                     3.590900679
       )
       neg_ctrl <- c(73.14191353,
                     100.7900693,
                     112.2837986,
                     113.7842185
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
