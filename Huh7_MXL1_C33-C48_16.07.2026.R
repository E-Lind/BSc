
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c33 = c(39.00422675,                 53.81443726,                 48.10417703,                 46.44869882         ),
         c34 = c(73.570363,                 68.63086358,                 82.97281618         ),
         c35 = c(23.82521134,                 36.99237527,                 65.01947621         ),
         c36 = c(113.7825294,                 128.8683076,                 114.5988729         ),
         c37 = c(102.4428145,                 110.4943643,                 121.409332         ),
         c38 = c(24.22095143,                 25.38952428,                 22.838969,                 28.37933035         ),
         c39 = c(37.85430134,                 30.71026024,                 30.54450522,                 33.75808056         ),
         c40 = c(87.00273496,                 104.4484502,                 107.514918         ),
         c41 = c(37.86258909,                 60.80101111,                 56.12257583,                 40.73636665         ),
         c42 = c(95.12265871,                 108.7290734,                 130.7869219         ),
         c43 = c(96.20213824, 118.2330515,                 99.09663517,                 127.0325709         ),
         c44 = c(113.2604011,                 123.137328,                 105.9671805         ),
         c45 = c(63.26247306,                 60.59796121,                 59.51640975,                 63.47381071         ), 
         c46 = c(64.84128957,                 65.16865573,                 68.87120835         ),
         c47 = c(114.2445715,                 115.8627548,                 123.804492         ),   
         c48 = c(114.4123985,                 114.1057517,                 96.85065473),
         ctrl = c(80.03273662,                  120.9908006,                  96.38239682,                  102.594066         )
         )
      
       pos_ctrl <- c(-3.176280457,                     1.877175535,                     2.020139234,                     4.371788497       )
       neg_ctrl <- c(80.03273662,                     120.9908006,                     96.38239682,                     102.594066       ) 
       
       
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
         labs(x = "MXL1 Compounds C33 - C48 [5µM]", y = "Huh7 Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14),axis.text.x = element_text(angle = 45, hjust = 1))
