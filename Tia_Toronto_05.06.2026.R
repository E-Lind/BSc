
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c01 = c(20.48455768,     30.34251618,   26.65427643),
         c02 = c(43.17824948,      53.05886048,       45.14901745),
         c03 = c(103.1950329,      110.3655805,   98.70159958,      90.49315534),
         c04 = c(76.3621106,                 81.06353513,                 74.88969774,                 70.67015378),
         c05 = c(78.79416595,                 79.45932588,                 75.27890897),
         c06 = c(82.30942293,                 75.17800236,                 74.09480074),
         c07 = c(84.15251314,                92.80782953,              81.39714475,             75.94200959),
         c08 = c(29.22430614,                 50.72153378,                 73.95476707,                 76.06968734         ),
         c09 = c(48.9484604,                 56.99421846,                 59.43245177         ),
         c10 = c(61.52060091,                 62.52554842,                 57.35459923         ),
         c11 = c(81.14590788,                 70.58572171,                 68.46462348,                 78.10017556         ),
         c12 = c(72.9436416,                 70.0214684,                 70.38184917         ),
         c13 = c(47.53782711,                 50.57326284,                 50.96865202,                 52.6428781         ), 
         c14 = c(71.64421151,                 67.62236214,                 66.97367676,                 60.42092474         ),
         ctrl = c(95.50759631,                  100.4437832,                  94.35025922         )
         )
      
       neg_ctrl <- c(-0.427308625,
                     1.034807634,
                     1.137773568,
                     1.561993215
       )
       pos_ctrl <- c(95.50759631,
                     109.6983613,
                     100.4437832,
                     94.35025922
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
         labs(x = "Potential P450 Compounds [5µM]", y = "Tiassale Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))