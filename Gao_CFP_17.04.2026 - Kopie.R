# Pakete laden
library(ggplot2)

# Beispiel: Vektoren (durch deine Daten ersetzen)
q <- c(50, 16.665,	5.555, 	1.85, 	0.615, 	0.205, 	0.065, 	0.0225)
w <- c(9.585620817,	9.571331879,	6.599984959,	27.78145446,	65.15078589,	96.73760999,	111.0039859,	74.89283297
	)
e <- c(7.747074991,	3.144989752,	5.284419017,	8.634559907,	18.6633,	13.396889,	21.7981314,	39.18011548
)
pos_ctrl <- c(-11.23561706,
              -2.821689103,
              -0.646762428,
              2.367451305
)
neg_ctrl <- c(83.5045499,
              126.3984357,
              110.3677521,
              79.72926224
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
    y = "Gaoua Normalised-FI after 48h",
    title = paste0(
      "IC50 = ",
      ifelse(is.na(LC50_numeric), "nicht bestimmbar", round(LC50_numeric, 4)),
      " µM | Z' = ",
      round(z_prime, 3)
    )
  )
