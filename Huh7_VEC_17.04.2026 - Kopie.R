
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
    library(purrr)
    library(tibble)
    library(tidyverse)
     
       vec_list <- list(
         c03 = c(67.24688433, 118.4859409, 134.4917087),
         c08 = c(10.18436502, 98.51684005, 54.75126172,  40.14007622),
         c14 = c(0.144196107,  1.050571635,  9.036976002, 5.401174168),
         c23 = c(9.354207436,   10.04634875,  30.10608714, 42.01668555),
         c26 = c(17.94417551,  28.78566279, 77.44772891 ),
         c31 = c(77.83087857,  47.48171799, 134.0961994 ),
         c35 = c(21.36780307, 22.40189515,  22.30301782,  51.58718715),
         c38 = c(16.34977856,  22.6037697, 67.33546194 ),
         c41 = c(42.96425996, 40.18539499, 71.98269647,  58.26758678),
         c46 = c(14.64620455, 6.284890308, 39.9443815,  72.04037491),
         c52 = c(19.1554228,  77.36739108,  60.19157483,  57.39005047),
         c55 = c(58.39118344, 77.25615408, 91.63662581),
         c57 = c(35.40632403,  27.40756,  32.94263055,   66.22721187), 
         c62 = c(51.89823875,  93.84900608,  90.4294984),
         c66 = c(19.47265424,   14.0405809, 26.10979504, 48.29951591),   
         c70 = c(29.86095375,   23.63992172,   68.77742301,  78.19136883 ), 
         c74 = c(39.60426719,   28.71497499,  48.44702168),
         c75 = c(50.23829549,   66.4040252,  97.23488705),
         c76 = c(70.63579892,  130.9666091, 82.39121774 ),
         c77 = c(50.37699379,  88.54115946,  83.42850398),
         c78 = c(61.73845042,  47.40088236,  92.28011155),
         ctrl = c(38.00481018,  44.35099523,    89.60500494,  60.46554743, 84.09928932, 100,  95.20032959)
         )
      
       # ---------------------------
       # 2. DataFrame
       # ---------------------------
       df <- do.call(rbind, lapply(names(vec_list), function(name){
         data.frame(Compound = as.character(name), Wert = vec_list[[name]])
       }))
       
       df$Compound <- factor(df$Compound, levels = c("ctrl", setdiff(names(vec_list), "ctrl")))
       
       ctrl_mean <- mean(vec_list$ctrl)
       
       # ---------------------------
       # t-Test
       # ---------------------------
       compounds <- setdiff(names(vec_list), "ctrl")
       
       stats <- data.frame(Compound = compounds,
                           p_value = sapply(compounds, function(comp){
                             t.test(vec_list[[comp]], mu = ctrl_mean, alternative = "less")$p.value
                           }))
       
       stats$sig <- ifelse(stats$p_value < 0.001, "***",
                           ifelse(stats$p_value < 0.01, "**",
                                  ifelse(stats$p_value < 0.05, "*", NA)))
       
       stats <- stats %>% filter(!is.na(sig))
       
       # gleiche Reihenfolge
       stats$Compound <- factor(stats$Compound, levels = levels(df$Compound))
       
       # 🔥 Position für Sterne: max + Offset
       stats$y_pos <- sapply(stats$Compound, function(c){
         max(vec_list[[as.character(c)]]) + 5
       })
       
       # ---------------------------
       # Plot
       # ---------------------------
       ggplot(df, aes(x = Compound, y = Wert)) +
         geom_jitter(width = 0.15, size = 3) +
         stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "red") +
         stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                      geom = "errorbar", width = 0.25, colour = "red") +
         # ⭐ Sterne direkt über Punkten
         geom_text(data = stats,
                   aes(x = Compound, y = y_pos, label = sig),
                   size = 5) +
         geom_hline(yintercept = ctrl_mean, linetype = "dashed", color = "blue") +
         labs(x = "VEC Compounds [5µM]", y = "Huh7 Normalised-FI after 48h") +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))