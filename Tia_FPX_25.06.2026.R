# Pakete laden
library(ggplot2)
library(scales)

# Beispiel: Vektoren (durch deine Daten ersetzen)
# Beispiel: Vektoren (durch deine Daten ersetzen)
q <- c(50, 16.665, 5.555, 1.85, 0.615, 0.205)
w <- c(53.16086279,	92.86639019,	115.8304259,	120.5832401,	123.7925965,	132.5995655)
e <- c(3.413147756,	2.515210719,	14.03208252,	6.530576733,	6.889731915,	19.38446796)

neg_ctrl <- c(79.6712966,
              93.54968951,
              91.87150286,
              97.48008689
)
pos_ctrl <- c(86.15890988,
              98.35033902,
              109.2392425,
              106.2515086
)

data <- data.frame(
  Konzentration = q,
  Wirkung = w,
  SD = e
)

# Achsengrenzen
y_min <- floor(min(w - e)/20)*20
y_max <- ceiling(max(w + e)/20)*20

# 4-Parameter-logistische Funktion
sigmoid <- function(conc, bottom, top, logLC50, slope) {
  bottom + (top - bottom) /
    (1 + exp(slope * (log(conc) - logLC50)))
}

# Fit
fit <- nls(
  Wirkung ~ sigmoid(Konzentration, bottom, top, logLC50, slope),
  data = data,
  start = list(
    bottom = min(w),
    top = max(w),
    logLC50 = log(median(q)),
    slope = 1
  ),
  algorithm = "port",
  lower = c(0, 0, log(min(q)/100), 0.001),
  upper = c(200, 200, log(max(q)*100), 20),
  control = nls.control(maxiter = 5000, warnOnly = TRUE)
)

# Parameter
coef_fit <- coef(fit)

print(coef_fit)
print(summary(fit))

# ---------------------------
# LC50
# ---------------------------
target_y <- coef_fit["bottom"] +
  (coef_fit["top"] - coef_fit["bottom"]) / 2

LC50_numeric <- tryCatch({
  uniroot(
    function(x)
      sigmoid(x,
              coef_fit["bottom"],
              coef_fit["top"],
              coef_fit["logLC50"],
              coef_fit["slope"]) - target_y,
    lower = min(q) / 10,
    upper = max(q) * 10
  )$root
}, error = function(e) NA)


# Vorhersagekurve
q_pred <- exp(seq(
  log(min(q)),
  log(max(q)),
  length.out = 400
))

pred <- sigmoid(
  q_pred,
  coef_fit["bottom"],
  coef_fit["top"],
  coef_fit["logLC50"],
  coef_fit["slope"]
)

curve_data <- data.frame(
  Konzentration = q_pred,
  Wirkung = pred
)

# ---------------------------
# Z'-Faktor
# ---------------------------
z_prime <- 1 - (3 * (sd(pos_ctrl) + sd(neg_ctrl))) /
  abs(mean(pos_ctrl) - mean(neg_ctrl))

# Plot
p <- ggplot(data, aes(x = Konzentration, y = Wirkung)) +
  
  geom_point(size = 2, colour = "black") +
  
  geom_errorbar(
    aes(
      ymin = Wirkung - SD,
      ymax = Wirkung + SD
    ),
    width = 0.05,
    colour = "black"
  ) +
  
  geom_line(
    data = curve_data,
    aes(x = Konzentration, y = Wirkung),
    linewidth = 1,
    colour = "black"
  ) +
  
  geom_vline(
    xintercept = LC50,
    colour = "red",
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = target_y,
    colour = "red",
    linetype = "dashed"
  ) +
  
  scale_x_log10(
    expand = c(0, 0),
    breaks = sort(unique(c(
      log_breaks()(range(q)),
      LC50
    ))),
    labels = function(x) round(x, 4)
  ) +
  
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(y_min, y_max, by = 20),
    expand = c(0, 0)
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = expression("Fenpyroximate ["*mu*"M]"),
    y = "Tiassale Normalised-FI after 48h",
    title = paste0(
      "IC50 = ",
      ifelse(is.na(LC50_numeric), "nicht bestimmbar", round(LC50_numeric, 4)),
      " µM | Z' = ",
      round(z_prime, 3)
    )
  )
print (p)