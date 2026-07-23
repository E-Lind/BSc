
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c65 = c(81.05074775,                 102.6641299,                 88.41027262,                 94.85973658),
         c66 = c(109.3100587,                 136.0315247,                 118.9684471,                 121.5066815         ),
         c67 = c(108.1019133,                 119.0610339,                 122.6696626,                 115.8927573         ),
         c68 = c(91.84953509,                 83.91416522,                 74.77290819,                 90.60074182         ),
         c69 = c(112.8238422,                 106.4805144,                 100.277196         ),
         c70 = c(97.81799911,                 115.5404757,                 123.1958268         ),
         c71 = c(97.28957664,                 109.6600821,                 111.3627767         ),
         c72 = c(10.29633435,                 67.00238806,                 86.63531511,                 33.17657553         ),
         c73 = c(68.03439262,                 83.22540944,                 83.91190701         ),
         c74 = c(128.8391078,                 132.097713,                 141.1463832         ),
         ctrl = c(85.44975188,                  99.30616324,                  96.4337129         )
         )
      
       pos_ctrl <- c(-1.630431714,                     2.350802513,                     1.85173685,                     0.733920093       )
       neg_ctrl <- c(85.44975188,                     99.30616324,                     96.4337129,                     118.810372       ) 
       
       
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
         labs(x = "MXL1 Compounds C65 - 74 [5µM]", y = "Huh7 Normalised-FI after 48h",
              title = paste0("Z' = ", round(z_prime, 3))) +
         theme_bw() +
         theme(text = element_text(size = 14),axis.text.x = element_text(angle = 45, hjust = 1))
