
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c03 = c(-2.417476894,
                 10.06601849,
                 24.96939143,
                 15.0018005
         ),
         c08 = c(42.38866883,
                 55.36190133,
                 61.42119794,
                 70.72620334
         ),
         c14 = c(53.71023887,
                 58.56919938,
                 51.40559357
         ),
         c26 = c(94.85776017,
                 120.9194574,
                 108.9833153
         ),
         c38 = c(50.71419998,
                 75.68119073,
                 101.1091106
         ),
         c57 = c(52.66354579,
                 74.09674709,
                 90.04681311
         ), 
         c62 = c(47.35325891,
                 75.18185092,
                 125.1542432,
                 103.2505101
         ),
         c66 = c(29.19457448,
                 56.96555035,
                 119.7575321,
                 94.19517465
         ),   
         c70 = c(74.83615412,
                 72.93482175,
                 84.5060617
         ), 
         c74 = c(46.40259273,
                 53.89268995,
                 63.25531149
         ),
         c76 = c(109.2137799,
                 98.97731365,
                 102.2710359
         ),
         c77 = c(104.335614,
                 83.3057256,
                 90.0564158
         ),
         c78 = c(120.2760773,
                 104.5180651,
                 83.26731485
         ),
         ctrl = c(105.5839635,
                  109.4154363,
                  110.9614692
         )
         )
      
       pos_ctrl <- c(-0.468131077,
                     5.495138639,
                     -0.929060137,
                     0.761013084
       )
       neg_ctrl <- c(74.03913096,
                     105.5839635,
                     109.4154363,
                     110.9614692
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
         labs(x = "VEC Compounds [5µM]", y = "Huh7 Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))
