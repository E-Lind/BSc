
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c49 = c(49.17942907,                 62.82612669,                 44.95851457         ),
         c50 = c(107.8964874),
         c51 = c(85.83135349,                 53.64742863,                 80.83117073,                 41.51248218         ),
         c52 = c(69.58368361,                 70.97116123,                 71.7372711         ),
         c53 = c(84.73628422,                 126.303593,                 111.8893234         ),
         c54 = c(34.31631273,                 37.41876531,                 37.3763661,                 38.70243795         ),
         c55 = c(42.82978179,                 88.29343178,                 42.36631456,                 91.92514346         ),
         c56 = c(19.37278409,                 30.15680398,                 45.71731423         ),
         c57 = c(35.20377207,                 40.29898754,                 39.7477978         ),
         c58 = c(109.8585475,                 112.925911,                 102.348039         ),
         c59 = c(121.505172,                 120.3574692,                 112.575021         ),
         c60 = c(33.64231149,                 35.93040681,                 29.41554881,                 24.87298512         ),
         c61 = c(88.23933623,                 87.94692788,                 88.32121057         ), 
         c62 = c(70.03253043,                 68.19181988,                 59.94444241,                 56.86830659         ),
         c63 = c(105.0235754,                 57.64172667,                 76.44943163,                 83.77279871         ),   
         c64 = c(110.9565408,                 103.8861069,                 95.74253445         ),
         ctrl = c(119.5358017,                  117.0810337         )
         )
      
       pos_ctrl <- c(-2.946014109,                     0.425454147,                     -0.124273548,                     -0.956175299       )
       neg_ctrl <- c(63.38316459,                     119.5358017,                     117.0810337       ) 
       
       
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
         labs(x = "MXL1 Compounds C49 - C64 [5µM]", y = "Huh7 Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14),axis.text.x = element_text(angle = 45, hjust = 1))
