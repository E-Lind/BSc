
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c17 = c(16.86113472,                 45.7255419,                 34.1940762         ),
         c18 = c(109.9992851,                 108.5632134,                 116.4679624         ),
         c19 = c(32.22900897,                 28.33858887,                 29.73018054),
         c20 = c(88.1071335,                 87.05390829,                 95.05714898         ),
         c21 = c(6.691077769,                 8.217698314,                 8.505230383),
         c22 = c(112.2582387,                 105.1080628,                 106.7125235         ),
         c23 = c(106.456763,                 99.06671221,                 99.96743421         ),
         c24 = c(19.12961977,                 62.82178572,                 88.42167134,                 58.24192408         ),
         c25 = c(52.69938602, 132.3267063,                 121.7563285         ),
         c26 = c(39.63017975,                 138.4157142, 126.650728         ),
         c27 = c(127.6308787,                 121.0176411,                 134.7556375         ),
         c28 = c(127.9152336,                 105.9404761,                 138.0026847         ),
         c29 = c(74.59153766,                 118.3663095         ), 
         c30 = c(6.60211757,                 6.451202948,                 36.11942907,                 61.23321075         ),
         c31 = c(113.3877155,                 118.0533602,                 124.0645279         ),   
         c32 = c(92.34704009,                 105.0206912,                 117.2352441         ), 
         ctrl = c(109.6418558,                  92.20724549,                  94.80456556,                  103.3463332         )
         )
      
       neg_ctrl <- c(0.787138897,                     3.775248413,                     1.271654263,                     5.031811214       )
       pos_ctrl <- c(109.6418558,                     92.20724549,                     94.80456556,                     103.3463332       ) 
       
       
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
         labs(x = "MXL1 Compounds C17 - 32 [5µM]", y = "Huh7 Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14),axis.text.x = element_text(angle = 45, hjust = 1))