
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
    library(emmeans)
     
      
         vec_list <- list(
           c03 = c(16.27260827,
                   29.67326548,
                   24.5758863
           ),
           c08 = c(-14.87804226,
                   -23.41641955,
                   16.55045283,
                   10.29895007
           ),
           c14 = c(105.5248324,
                   153.9552777,
                   175.3920547
           ),
           c23 = c(20.13037322,
                   12.36141166,
                   28.14512035
           ),
           c26 = c(-1.060618204,
                   -2.107878496,
                   -1.552189362
           ),
           c31 = c(51.64435895,
                   89.9334776,
                   95.75752718,
                   133.9077235
           ),
           c35 = c(31.87464935,
                   125.7006225,
                   35.01643023,
                   191.1223318
           ),
           c38 = c(35.40113809,
                   78.50979135,
                   127.0257273,
                   92.58368732
           ),
           c41 = c(7.969330234,
                   15.0543667,
                   18.4953648
           ),
           c46 = c(-19.59071358,
                   0.82017579,
                   3.128422965
           ),
           c52 = c(41.97323074,
                   108.2391601,
                   34.89888061,
                   159.6297187
           ),
           c55 = c(-1.285031124,
                   27.71766717,
                   144.9253293,
                   187.8416286
           ),
           c57 = c(126.0960167,
                   165.2721007,
                   174.526462
           ), 
           c62 = c(64.89540755,
                   85.07119767,
                   133.7260559,
                   177.5506933
           ),
           c66 = c(26.0612861,
                   32.04563063,
                   31.21209693
           ),   
           c70 = c(148.0350512,
                   120.464321,
                   188.3225134
           ), 
           ctrl = c(85.70169111,
                    107.0209185,
                    101.4853998,
                    105.7919906,
                    83.5045499,
                    110.3677521,
                    79.72926224),
         c74 = c(6.02241107,
                 21.96585696,
                 22.01398812,
                 33.59554787
         ),
         c75 = c(125.1079191,
                 152.5186132,
                 146.5503497
         ),
         c76 = c(135.7388885,
                 122.1839513,
                 137.7874709
         ),
         c77 = c(123.2247875,
                 143.2292998,
                 113.4180642
         ),
         c78 = c(146.5954727,
                 130.0082725,
                 120.1955328
         )
         )
      
         pos_ctrl <- c(-11.23561706,
                       -2.821689103,
                       -0.646762428,
                       2.367451305,
                       -32.79634528,
                       8.517004622,
                       22.53746894,
                       38.89823943
                       
         )
         neg_ctrl <- c(83.5045499,
                       126.3984357,
                       110.3677521,
                       79.72926224,
                       85.70169111,
                       107.0209185,
                       101.4853998,
                       105.7919906)         
         
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
           labs(x = "VEC Compounds [5 µM]",
                y = "Gaoua Normalised-FI after 48 h",
                title = paste0("Z' = ", round(z_prime, 3))) +
           theme_bw() +
           theme(axis.text.x = element_text(angle = 45, hjust = 1))
