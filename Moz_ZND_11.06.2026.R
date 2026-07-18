
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
    library(emmeans)
     
       vec_list <- list(
         c05 = c(14.87681977,
                 19.57602289,
                 13.72082949,
                 12.07690641
         ),
         c11 = c(97.96903679,
                 108.8335214,
                 105.3267896
         ),
         c12 = c(71.79167094,
                 93.31771492,
                 57.09440587,
                 111.8751211
         ),
         c18 = c(90.90997184,
                 91.81743562,
                 86.73062234
         ),
         c25 = c(92.19820559,
                 92.84574256,
                 87.30291734
         ),
         c28 = c(90.8170532,
                 105.263845,
                 97.28321504
         ),
         c31 = c(98.30631187,
                 98.46535476,
                 103.8894296
         ),
         c32 = c(90.54406912,
                 89.3999098,
                 93.27628362
         ),
         c34 = c(107.3124125,
                 98.28732167,
                 104.3974173
         ),
         c38 = c(59.85828566,
                 91.5814086,
                 103.8015999
         ),
         c40 = c(79.19267928,
                 85.65409357,
                 85.53777862
         ),
         c53 = c(80.60982268,
                 81.70888029,
                 81.8940347
         ), 
         c66 = c(71.22629193,
                 78.12685451,
                 80.77836067
         ),
         c67 = c(110.6143329,
                 103.2413891,
                 103.8514492
         ),
         c69 = c(101.0361526,
                 91.62888409,
                 96.78472239
         ),   
         c70 = c(101.316258,
                 98.86414888,
                 94.67918437
         ), 
         c74 = c(-2.427184466,
                 0.713319249,
                 0.829634201,
                 0.777411161
         ),
         c79 = c(90.3043179,
                 86.25228476,
                 85.86773328
         ),
         c80 = c(64.13582738,
                 62.93707124,
                 74.88902604
         ),
         ctrl = c(101.1542802,
                  101.52821,
                  100.8897933,
                  96.42771641,
                  99.29380207,
                  100.3904859,
                  99.39350061,
                  100.9222114)
         )
      
       neg_ctrl <- c(-0.576855114,
                     0.339728901,
                     0.768380132,
                     1.618842414,
                     -1.541173119,
                     -0.07180668,
                     0.307997246,
                     1.221900444
                     
       )
       pos_ctrl <- c(101.1542802,
                     101.52821,
                     100.8897933,
                     96.42771641,
                     99.29380207,
                     100.3904859,
                     99.39350061,
                     100.9222114
                     
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
         labs(x = "ZND Compounds [5µM]", y = "Moz Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14),axis.text.x = element_text(angle = 45, hjust = 1))