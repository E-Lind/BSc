
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c03 = c(36.25693993,
                 29.98677565,
                 33.90619834
         ),
         c08 = c(2.012722733,
                 3.513002041,
                 4.773875076,
                 4.605150655
         ),
         c14 = c(81.31548047,
                 110.2266379,
                 95.88734225,
                 100.0986126
         ),
         c26 = c(80.07056785,
                 105.0554625,
                 73.07534457,
                 100.0689718
         ),
         c38 = c(51.47177856,
                 79.98620564,
                 57.16280767,
                 72.17016086
         ),
         c62 = c(73.14146631,
                 75.40328557,
                 75.85701745
         ),
         c70 = c(81.79657307,
                 91.44122576,
                 110.9562571,
                 100.2080555
         ),
         c74 = c(59.88291893,
                 84.50528404,
                 74.03524972,
                 71.88743345
         ),
         c76 = c(104.8160562,
                 102.7799628,
                 118.3436506
         ),
         c77 = c(99.35987323,
                 108.8403616,
                 113.1633549
         ),
         c78 = c(100.2764572,
                 102.8369643,
                 103.892632
         ),
         ctrl = c(101.1542802,
                  101.52821,
                  100.8897933,
                  96.42771641
         )
         )
      
       neg_ctrl <- c(-0.576855114,
                     0.339728901,
                     0.768380132,
                     1.618842414
       )
       pos_ctrl <- c(101.1542802,
                     101.52821,
                     100.8897933,
                     96.42771641
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
         labs(x = "VEC Compounds [5µM]", y = "Moz Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))