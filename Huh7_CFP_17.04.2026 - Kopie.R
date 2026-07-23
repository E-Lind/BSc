# Pakete laden
library(ggplot2)

# Beispiel: Vektoren (durch deine Daten ersetzen)
q <- c(50, 16.665,	5.555, 	1.85, 	0.615, 	0.205, 	0.065, 	0.0225)
w <- c(1.085635539,	24.99403199,	81.88437249,	77.26079125,	98.24160681,	117.12331,	92.55571953,	67.99084182)
e <- c(0.615048665,	8.055855293,	17.60495298,	26.32530594,	21.5551648,	45.99434934,	29.16201603,	15.52860549)

pos_ctrl <- c(-1.76056338,
              2.085006185,
              1.800711821,
              2.284663296
              
)
neg_ctrl <- c(55.89748041,
              65.23145034,
              147.0800148,
              131.7910545
              
)

data <- data.frame(Konzentration = q, Wirkung = w, SD = e)

# ---------------------------
# y-Achse
# ---------------------------
y_min <- floor(min(w - e) / 20) * 20
y_max <- ceiling(max(w + e) / 20) * 20

# ---------------------------
# Sigmoid-Funktion
# ---------------------------
sigmoid <- function(conc, bottom, top, logLC50, slope) {
  bottom + (top - bottom) / (1 + exp(slope * (log(conc) - logLC50)))
}

# ---------------------------
# Fit
# ---------------------------
fit <- nls(
  Wirkung ~ sigmoid(Konzentration, bottom, top, logLC50, slope),
  data = data,
  start = list(
    bottom = min(w),
    top = max(w),
    logLC50 = log(median(q)),
    slope = 1.5
  ),
  control = nls.control(maxiter = 5000, warnOnly = TRUE)
)

coef_fit <- coef(fit)

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

# ---------------------------
# Z'-Faktor
# ---------------------------
z_prime <- 1 - (3 * (sd(pos_ctrl) + sd(neg_ctrl))) /
  abs(mean(pos_ctrl) - mean(neg_ctrl))

# ---------------------------
# Kurve
# ---------------------------
q_pred <- exp(seq(log(min(q)), log(max(q)), length.out = 400))

curve_data <- data.frame(
  Konzentration = q_pred,
  Wirkung = sigmoid(
    q_pred,
    coef_fit["bottom"],
    coef_fit["top"],
    coef_fit["logLC50"],
    coef_fit["slope"]
  )
)

# ---------------------------
# Plot
# ---------------------------
ggplot(data, aes(x = Konzentration, y = Wirkung)) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = Wirkung - SD, ymax = Wirkung + SD),
                width = 0.05, color = "black") +
  geom_line(data = curve_data,
            aes(x = Konzentration, y = Wirkung),
            color = "black", linewidth = 1) +
  geom_vline(xintercept = LC50_numeric, color = "red", linetype = "dashed") +
  geom_hline(yintercept = target_y, color = "red", linetype = "dashed") +
  scale_x_log10(
    expand = c(0, 0),
    breaks = scales::log_breaks()(q),
    labels = function(x) round(x, 3)
  ) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(y_min, y_max, by = 20),
    expand = c(0, 0)
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  theme_classic(base_size = 14) +
  labs(
    x = "Chlorfenapyr [µM]",
    y = "Huh7 Normalised-FI after 48h",
    title = paste0(
      "IC50 = ",
      ifelse(is.na(LC50_numeric), "nicht bestimmbar", round(LC50_numeric, 4)),
      " µM | Z' = ",
      round(z_prime, 3))
  )
