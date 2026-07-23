# Pakete laden
library(ggplot2)

# Beispiel: Vektoren (durch deine Daten ersetzen)
q <- c(10,	5,	2.5, 1.25,	0.625, 0.3125, 0.15625)
w <- c(64.58262792,	69.98824865,	88.51952915,	93.97295198,	104.3758838,	92.15448045,	97.96442727)
e <- c(11.57072426,	6.789899204,	11.99276494,	10.09630493,	21.88139715,	6.265139074,	6.611168057)
pos_ctrl <- c(-1.001852331,  -4.33206525,  -4.296213675,  -2.053498516)
neg_ctrl <- c(84.07492979,  107.0836736, 110.903858, 97.93753859)

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
    x = "ZND C11 [µM]",
    y = "Gaoua Normalised-FI after 48h",
    title = paste0(
      "IC50 = ",
      ifelse(is.na(LC50_numeric), "nicht bestimmbar", round(LC50_numeric, 4)),
      " µM | Z' = ",
      round(z_prime, 3)
    )
  )
