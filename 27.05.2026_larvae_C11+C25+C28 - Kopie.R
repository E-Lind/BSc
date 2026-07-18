# Pakete laden
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)


H16 = list(
  C11_20uM = list(
    dead = c(1.224489796),
    puppae = c(14.28571429),
    mosquito = c(0),
    larvae = c(84.48979592)
  ),
  C11_0.2uM = list(
    dead = c(0.865800866),
    puppae = c(0),
    mosquito = c(0),
    larvae = c(99.13419913)
  ),
  C25_20uM = list(
    dead = c(1.754385965),
    puppae = c(7.30994152),
    mosquito = c(0),
    larvae = c(90.93567251)
  ),
  C25_0.2uM = list(
    dead = c(1.632653061),
    puppae = c(5.714285714),
    mosquito = c(0),
    larvae = c(92.65306122)
  ),
  C28_20uM = list(
    dead = c(2.554744526),
    puppae = c(3.649635036),
    mosquito = c(0),
    larvae = c(93.79562044)
  ),
  C28_0.2uM = list(
    dead = c(0.896860987),
    puppae = c(0),
    mosquito = c(0),
    larvae = c(99.10313901)
  ),
  ctrl56 = list(
    dead = c(0.512820513),
    puppae = c(0.512820513),
    mosquito = c(0),
    larvae = c(98.97435897)
  ),
  ctrl61 = list(
    dead = c(0.387596899),
    puppae = c(0),
    mosquito = c(0),
    larvae = c(99.6124031)
  )
)

H24 = list(
  C11_20uM = list(
    dead = c(1.224489796),
    puppae = c(62.04081633),
    mosquito = c(0),
    larvae = c(36.73469388)
  ),
  C11_0.2uM = list(
    dead = c(0.865800866),
    puppae = c(40.25974026),
    mosquito = c(0),
    larvae = c(58.87445887)
  ),
  C25_20uM = list(
    dead = c(1.754385965),
    puppae = c(64.61988304),
    mosquito = c(0),
    larvae = c(33.62573099)
  ),
  C25_0.2uM = list(
    dead = c(1.632653061),
    puppae = c(32.65306122),
    mosquito = c(0),
    larvae = c(65.71428571)
  ),
  C28_20uM = list(
    dead = c(2.554744526),
    puppae = c(56.20437956),
    mosquito = c(0),
    larvae = c(41.24087591)
  ),
  C28_0.2uM = list(
    dead = c(1.34529148),
    puppae = c(33.1838565),
    mosquito = c(0),
    larvae = c(65.47085202)
  ),
  ctrl56 = list(
    dead = c(1.025641026),
    puppae = c(44.61538462),
    mosquito = c(0),
    larvae = c(54.35897436)
  ),
  ctrl61 = list(
    dead = c(0.387596899),
    puppae = c(18.21705426),
    mosquito = c(0),
    larvae = c(81.39534884)
  )
)

H40 = list(
  C11_20uM = list(
    dead = c(5.306122449),
    puppae = c(75.51020408),
    mosquito = c(13.06122449),
    larvae = c(6.12244898)
  ),
  C11_0.2uM = list(
    dead = c(3.03030303),
    puppae = c(82.68398268),
    mosquito = c(0),
    larvae = c(14.28571429)
  ),
  C25_20uM = list(
    dead = c(1.754385965),
    puppae = c(87.13450292),
    mosquito = c(10.52631579),
    larvae = c(0.584795322)
  ),
  C25_0.2uM = list(
    dead = c(3.265306122),
    puppae = c(65.71428571),
    mosquito = c(5.306122449),
    larvae = c(25.71428571)
  ),
  C28_20uM = list(
    dead = c(8.759124088),
    puppae = c(77.73722628),
    mosquito = c(3.649635036),
    larvae = c(9.854014599)
  ),
  C28_0.2uM = list(
    dead = c(3.139013453),
    puppae = c(78.02690583),
    mosquito = c(0),
    larvae = c(18.83408072)
  ),
  ctrl56 = list(
    dead = c(2.564102564),
    puppae = c(77.43589744),
    mosquito = c(0.512820513),
    larvae = c(19.48717949)
  ),
  ctrl61 = list(
    dead = c(1.937984496),
    puppae = c(65.89147287),
    mosquito = c(0),
    larvae = c(32.17054264)
  )
)


# =========================================================
# 1. LIST -> LONG FORMAT (FOR PLOT)
# =========================================================

convert_to_df <- function(lst){
  
  do.call(rbind,
          lapply(names(lst), function(n){
            
            values <- unlist(lst[[n]])
            
            data.frame(
              sample = n,
              category = names(values),
              value = as.numeric(values)
            )
          }))
}

# =========================================================
# 2. CREATE DATASETS
# =========================================================

df_H16 <- convert_to_df(H16)
df_H24 <- convert_to_df(H24)
df_H40 <- convert_to_df(H40)

# =========================================================
# 3. NORMALIZE TO 100%
# =========================================================

normalize <- function(df){
  df %>%
    group_by(sample) %>%
    mutate(value = value / sum(value) * 100) %>%
    ungroup()
}

df_H16 <- normalize(df_H16)
df_H24 <- normalize(df_H24)
df_H40 <- normalize(df_H40)

# =========================================================
# 4. ORDERING
# =========================================================

sample_order <- c(
  "ctrl56","ctrl61",
  "C11_20uM","C11_0.2uM",
  "C25_20uM","C25_0.2uM",
  "C28_20uM","C28_0.2uM"
)

cat_order <- c("larvae","puppae","mosquito","dead")

df_H16$sample <- factor(df_H16$sample, levels = sample_order)
df_H24$sample <- factor(df_H24$sample, levels = sample_order)
df_H40$sample <- factor(df_H40$sample, levels = sample_order)

df_H16$category <- factor(df_H16$category, levels = cat_order)
df_H24$category <- factor(df_H24$category, levels = cat_order)
df_H40$category <- factor(df_H40$category, levels = cat_order)

# =========================================================
# 5. WIDE FORMAT (FOR GLM)  <<< FIX FOR YOUR ERROR
# =========================================================

df_glm_H16 <- df_H16 %>%
  pivot_wider(names_from = category, values_from = value)

df_glm_H24 <- df_H24 %>%
  pivot_wider(names_from = category, values_from = value)

df_glm_H40 <- df_H40 %>%
  pivot_wider(names_from = category, values_from = value)

# =========================================================
# 6. COLORS
# =========================================================

cols <- c(
  larvae   = "gray60",
  puppae   = "#fc8d62",
  mosquito = "#8da0cb",
  dead     = "#e78ac3"
)

# =========================================================
# 7. STATISTICS (GLM + emmeans)
# =========================================================

make_stats <- function(df_glm){
  
  model <- glm(
    cbind(larvae, puppae + mosquito + dead) ~ sample,
    family = binomial,
    data = df_glm
  )
  
  emm <- emmeans(model, ~ sample)
  
  sig <- contrast(
    emm,
    method = "trt.vs.ctrl",
    ref = "ctrl56"
  ) %>%
    summary(infer = TRUE, adjust = "holm") %>%
    as.data.frame()
  
  # FIX: emmeans uses "contrast" column
  sig$sample <- sig$contrast
  
  sig$label <- ifelse(sig$p.value < 0.001, "***",
                      ifelse(sig$p.value < 0.01, "**",
                             ifelse(sig$p.value < 0.05, "*", "ns")))
  
  sig$y <- 108
  
  sig
}

sig_H16 <- make_stats(df_glm_H16)
sig_H24 <- make_stats(df_glm_H24)
sig_H40 <- make_stats(df_glm_H40)

# =========================================================
# 8. LABELS
# =========================================================

sample_labels <- c(
  ctrl56 = "Control 56 µl DMSO",
  ctrl61 = "Control 61 µl DMSO",
  C11_20uM = "C11 20 µM",
  C11_0.2uM = "C11 0.2 µM",
  C25_20uM = "C25 20 µM",
  C25_0.2uM = "C25 0.2 µM",
  C28_20uM = "C28 20 µM",
  C28_0.2uM = "C28 0.2 µM"
)

# =========================================================
# 9. PLOT FUNCTION
# =========================================================

plot_fun <- function(df_plot, sig, title){
  
  ggplot(df_plot,
         aes(x = sample,
             y = value,
             fill = category)) +
    
    geom_bar(stat = "identity", width = 0.8) +
    
    scale_fill_manual(
      values = cols,
      labels = c(
        larvae = "no change",
        puppae = "puppae",
        mosquito = "mosquito",
        dead = "dead"
      ),
      name = NULL
    ) +
    
    scale_x_discrete(labels = sample_labels) +
    
    coord_cartesian(ylim = c(0,110)) +
    
    labs(
      x = title,
      y = "Normalized Percentage"
    ) +
    
    geom_text(
      data = sig,
      aes(x = sample, y = y, label = label),
      inherit.aes = FALSE,
      size = 5
    ) +
    
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# =========================================================
# 10. FINAL PLOTS
# =========================================================

plot_H16 <- plot_fun(df_H16, sig_H16,
                     "Larvae after 16 hours of exposure")

plot_H24 <- plot_fun(df_H24, sig_H24,
                     "Larvae after 24 hours of exposure")

plot_H40 <- plot_fun(df_H40, sig_H40,
                     "Larvae after 40 hours of exposure")

# =========================================================
# 11. OUTPUT
# =========================================================

plot_H16
plot_H24
plot_H40