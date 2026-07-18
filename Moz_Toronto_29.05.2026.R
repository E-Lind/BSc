
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c01 = c(90.70508949,
                 79.63681769,
                 84.51095995
         ),
         c02 = c(112.0500867,
                 99.54588275,
                 113.6354093,
                 106.5560855
         ),
         c03 = c(98.17461217,
                 100.2058764,
                 95.71747367
         ),
         c04 = c(63.29305613,
                 47.62489618,
                 40.43631674,
                 68.23631935
         ),
         c05 = c(83.40725436,
                 84.51987878,
                 89.7351664
         ),
         c06 = c(58.03763376,
                 65.46925139,
                 56.86926663
         ),
         c07 = c(87.94917008,
                 79.80850522,
                 84.89224006
         ),
         c08 = c(90.70508949,
                 79.63681769,
                 84.51095995
         ),
         c09 = c(112.0500867,
                 99.54588275,
                 113.6354093,
                 106.5560855
         ),
         c10 = c(98.17461217,
                 100.2058764,
                 95.71747367
         ),
         c11 = c(63.29305613,
                 47.62489618,
                 40.43631674,
                 68.23631935
         ),
         c12 = c(83.40725436,
                 84.51987878,
                 89.7351664
         ),
         c13 = c(58.03763376,
                 65.46925139,
                 56.86926663
         ), 
         c14 = c(87.94917008,
                 79.80850522,
                 84.89224006
         ),
         ctrl = c(95.70409542,
                  103.9228001,
                  100.3731045
         )
         )
      
       neg_ctrl <- c(-0.316061146,
                     1.122100682,
                     0.734131445,
                     0.424201997
       )
       pos_ctrl <- c(95.70409542,
                     103.9228001,
                     100.3731045,
                     119.1160322
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
         labs(x = "Potential P450 Compounds [5µM]", y = "Moz Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))