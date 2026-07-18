# Pakete laden
library(ggplot2)

# Beispiel: Vektoren (durch deine Daten ersetzen)
q <- c(50, 16.665, 5.555, 1.85, 0.615, 0.205, 0.065, 0.0225)
w <- c(18.34735071,	52.60618167,	80.59083263,	104.3424096,	130.6507569,	146.7672414,	141.0697014,	140.966148)
e <- c(1.143869063,	7.569317971,	15.13482405,	20.05936188,	18.37468019,	12.18458786,	9.068865199,	11.09572717)

neg_ctrl <- c(95.91673675,
              100.7884777,
              101.1248949,
              86.4529016)
pos_ctrl <- c(81.16852397,
              106.2316022,
              107.3060345,
              105.2938394)

data <- data.frame(Konzentration = q, Wirkung = w, SD = e)
y_min <- floor(min(w - e)/20)*20
y_max <- ceiling(max(w + e)/20)*20


# --- Sigmoid-Funktion ---
sigmoid <- function(conc, bottom, top, logLC50, slope){
  bottom + (top - bottom) / (1 + exp(slope * (log(conc) - logLC50)))
}

# --- Fit ---
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
LC50 <- coef_fit["LC50"]


# --- Numerische Berechnung des LC50 (y = 50) ---
target_y <- 50
LC50_numeric <- uniroot(
  function(x) sigmoid(x, coef_fit["bottom"], coef_fit["top"], coef_fit["logLC50"], coef_fit["slope"]) - target_y,
  lower = min(q)/10,
  upper = max(q)*10
)$root

# ---------------------------
# Z'-Faktor
# ---------------------------
z_prime <- 1 - (3 * (sd(pos_ctrl) + sd(neg_ctrl))) /
  abs(mean(pos_ctrl) - mean(neg_ctrl))


# --- Glatte Kurve ---
q_pred <- exp(seq(log(min(q)), log(max(q)), length.out = 400))
pred <- sigmoid(q_pred, coef_fit["bottom"], coef_fit["top"], coef_fit["logLC50"], coef_fit["slope"])
curve_data <- data.frame(Konzentration = q_pred, Wirkung = pred)

# --- Plot ---
p <- ggplot(data, aes(x = Konzentration, y = Wirkung)) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = pmax(Wirkung - SD, 0), ymax = Wirkung + SD),
                width = 0.05, color = "black") +
  geom_line(data = curve_data, aes(x = Konzentration, y = Wirkung),
            color = "black", linewidth = 1) +
  geom_vline(xintercept = LC50_numeric, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 50, color = "red", linetype = "dashed") +
  scale_x_log10(
    expand = c(0, 0),
    breaks = sort(c(scales::log_breaks()(range(q)), LC50_numeric)),  # LC50 als zusätzlicher Tick
    labels = function(x) round(x, 3)
  ) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(y_min, y_max, by = 20),
    expand = c(0, 0)) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  theme_classic(base_size = 14) +
  labs(
    x = "Fenpyroximate [µm]",
    y = "Tiefora Normalised-FI after 48h",
    title = paste0(
      "IC50 = ",
      ifelse(is.na(LC50_numeric), "nicht bestimmbar", round(LC50_numeric, 4)),
      " µM | Z' = ",
      round(z_prime, 3)
    )
  )
print (p)