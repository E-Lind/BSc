# Pakete laden
library(ggplot2)

# Beispiel: Vektoren (durch deine Daten ersetzen)
q <- c(50, 16.665, 5.555, 1.85, 0.615, 0.205, 0.065, 0.0225)
w <- c(18.18870911,	128.6015187,	141.7478333,	160.4631195,	167.9385012,	167.4394269,	158.3373581,	150.2354504)
e <- c(4.578889883,	23.09864793,	30.20422246,	32.46454878,	35.78093502,	33.54561318,	19.6465983,	38.00302886)

neg_ctrl <- c(-733.75,
              -1029.75,
              -636.75,
              410.25
)
pos_ctrl <- c(37924.25,
              34146.25,
              42450.25,
              34555.25
)

data <- data.frame(
  Konzentration = q,
  Wirkung = w,
  SD = e
)

y_min <- floor(min(w - e) / 20) * 20
y_max <- ceiling(max(w + e) / 20) * 20

# ---------------------------
# Sigmoid-Funktion
# ---------------------------
sigmoid <- function(conc, bottom, top, logLC50, slope){
  bottom + (top - bottom) /
    (1 + exp(slope * (log(conc) - logLC50)))
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
  control = nls.control(
    maxiter = 5000,
    warnOnly = TRUE
  )
)

coef_fit <- coef(fit)

print(coef_fit)

# ---------------------------
# Konzentration bei y = 50 bestimmen
# ---------------------------
target_y <- 50

f <- function(x){
  sigmoid(
    x,
    coef_fit["bottom"],
    coef_fit["top"],
    coef_fit["logLC50"],
    coef_fit["slope"]
  ) - target_y
}

lower <- min(q) / 10
upper <- max(q) * 10

if (f(lower) * f(upper) < 0) {
  LC50_numeric <- uniroot(
    f,
    lower = lower,
    upper = upper
  )$root
} else {
  LC50_numeric <- NA
  message("y = 50 wird von der Fit-Kurve nicht erreicht.")
}

# ---------------------------
# Z'-Faktor
# ---------------------------
z_prime <- 1 -
  (3 * (sd(pos_ctrl) + sd(neg_ctrl))) /
  abs(mean(pos_ctrl) - mean(neg_ctrl))

# ---------------------------
# Glatte Kurve
# ---------------------------
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
# x-Achsen-Breaks
# ---------------------------
x_breaks <- log_breaks()(range(q))

if (!is.na(LC50_numeric)) {
  x_breaks <- sort(unique(c(x_breaks, LC50_numeric)))
}

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(data,
            aes(x = Konzentration,
                y = Wirkung)) +
  
  geom_point(size = 2) +
  
  geom_errorbar(
    aes(
      ymin = pmax(Wirkung - SD, 0),
      ymax = Wirkung + SD
    ),
    width = 0.05
  ) +
  
  geom_line(
    data = curve_data,
    aes(x = Konzentration,
        y = Wirkung),
    linewidth = 1
  ) +
  
  geom_hline(
    yintercept = 50,
    colour = "red",
    linetype = "dashed"
  ) +
  
  {
    if (!is.na(LC50_numeric))
      geom_vline(
        xintercept = LC50_numeric,
        colour = "red",
        linetype = "dashed"
      )
  } +
  
  scale_x_log10(
    breaks = x_breaks,
    labels = function(x) round(x, 3),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(y_min, y_max, 20),
    expand = c(0, 0)
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Fenpyroximate [µm]",
    y = "Tiefora Normalised-FI after 48h",
    title = paste0(
      "IC50 = ",
      ifelse(is.na(LC50_numeric), "not determinable", round(LC50_numeric, 4)),
      " | Z' = ",
      round(z_prime, 3)
    )
  )
print (p)