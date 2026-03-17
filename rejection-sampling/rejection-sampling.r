library(latex2exp)
set.seed(42)

# --- Example 1: Beta(2.7, 6.3) via Uniform(0,1) proposal ---
# Target: f(x) = Beta(2.7, 6.3)
# Proposal: g(x) = Uniform(0,1), so g(x) = 1
# Constant: c = max f(x)/g(x) = max f(x), attained at the mode

a <- 2.7; b <- 6.3
mode_beta <- (a - 1) / (a + b - 2)  # mode of Beta
c_beta <- dbeta(mode_beta, a, b)     # c = f(mode) since g = 1

N <- 10000
x_prop <- runif(N)                         # proposals from g
u <- runif(N)                              # uniform for accept/reject
accept_prob <- dbeta(x_prop, a, b) / c_beta
accepted <- u <= accept_prob

cat(sprintf("Beta example: accepted %d / %d = %.1f%% (theory: 1/c = %.1f%%)\n",
            sum(accepted), N, 100 * mean(accepted), 100 / c_beta))

# --- Example 2: N(0,1) via t_2 proposal ---
# Target: f(x) = dnorm(x)
# Proposal: g(x) = dt(x, df=2)
# Constant: c = max f(x)/g(x), found numerically

ratio_fn <- function(x) dnorm(x) / dt(x, df = 2)
opt <- optimize(ratio_fn, interval = c(-5, 5), maximum = TRUE)
c_norm <- opt$objective

M <- 10000
t_prop <- rt(M, df = 2)
u2 <- runif(M)
accept_prob2 <- dnorm(t_prop) / (c_norm * dt(t_prop, df = 2))
accepted2 <- u2 <= accept_prob2

cat(sprintf("Normal example: accepted %d / %d = %.1f%% (theory: 1/c = %.1f%%)\n",
            sum(accepted2), M, 100 * mean(accepted2), 100 / c_norm))

xgrid  <- seq(0, 1, length.out = 500)
xgrid2 <- seq(-5, 5, length.out = 500)

# ============================================================
# Figure 1: 2x2 panel — both examples, envelope + histogram
# ============================================================

pdf("fig_envelopes.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Beta envelope with accept/reject points
plot(NULL, xlim = c(0, 1), ylim = c(0, c_beta + 0.3),
     xlab = "x", ylab = "Density",
     main = TeX("(a) Envelope: $c \\cdot g(x) \\geq f(x)$"))
polygon(c(0, 0, 1, 1), c(0, c_beta, c_beta, 0),
        col = rgb(0.85, 0.85, 0.85, 0.4), border = NA)
lines(xgrid, dbeta(xgrid, a, b), col = "firebrick", lwd = 2)
abline(h = c_beta, col = "steelblue", lwd = 2, lty = 2)

n_show <- 2000
points(x_prop[1:n_show], u[1:n_show] * c_beta,
       pch = 20, cex = 0.35,
       col = ifelse(accepted[1:n_show],
                    rgb(0.2, 0.6, 0.2, 0.5),
                    rgb(0.8, 0.2, 0.2, 0.3)))
legend("topright", cex = 0.7,
       legend = c(TeX("$f(x) = \\mathrm{Beta}(2.7, 6.3)$"),
                  TeX("$c \\cdot g(x)$"), "Accepted", "Rejected"),
       col = c("firebrick", "steelblue", rgb(0.2, 0.6, 0.2), rgb(0.8, 0.2, 0.2)),
       lwd = c(2, 2, NA, NA), lty = c(1, 2, NA, NA),
       pch = c(NA, NA, 20, 20), bg = "white")

# Panel B: Beta histogram vs true density
hist(x_prop[accepted], breaks = 50, freq = FALSE,
     col = rgb(0.2, 0.6, 0.2, 0.3), border = "white",
     main = "(b) Accepted Samples vs True Density",
     xlab = "x", ylab = "Density", xlim = c(0, 1))
lines(xgrid, dbeta(xgrid, a, b), col = "firebrick", lwd = 2)
legend("topright", cex = 0.7,
       legend = c("Accepted samples", TeX("$\\mathrm{Beta}(2.7, 6.3)$")),
       fill = c(rgb(0.2, 0.6, 0.2, 0.3), NA),
       border = c("white", NA),
       col = c(NA, "firebrick"), lwd = c(NA, 2), lty = c(NA, 1),
       bg = "white")

# Panel C: Normal envelope
plot(xgrid2, c_norm * dt(xgrid2, df = 2), type = "l",
     col = "steelblue", lwd = 2, lty = 2,
     ylim = c(0, 0.55), xlab = "x", ylab = "Density",
     main = TeX("(c) Envelope: $c \\cdot t_2(x) \\geq \\phi(x)$"))
lines(xgrid2, dnorm(xgrid2), col = "firebrick", lwd = 2)
polygon(c(xgrid2, rev(xgrid2)),
        c(c_norm * dt(xgrid2, df = 2), rev(dnorm(xgrid2))),
        col = rgb(0.8, 0.2, 0.2, 0.15), border = NA)
legend("topright", cex = 0.7,
       legend = c(TeX("$f(x) = \\phi(x)$"), TeX("$c \\cdot t_2(x)$"), "Wasted area"),
       col = c("firebrick", "steelblue", rgb(0.8, 0.2, 0.2, 0.4)),
       lwd = c(2, 2, NA), lty = c(1, 2, NA),
       pch = c(NA, NA, 15), bg = "white")

# Panel D: Normal histogram
hist(t_prop[accepted2], breaks = 60, freq = FALSE,
     col = rgb(0.2, 0.6, 0.2, 0.3), border = "white",
     main = "(d) Accepted Samples vs True Density",
     xlab = "x", ylab = "Density", xlim = c(-4, 4))
lines(xgrid2, dnorm(xgrid2), col = "firebrick", lwd = 2)
legend("topright", cex = 0.7,
       legend = c("Accepted samples", TeX("$\\phi(x)$")),
       fill = c(rgb(0.2, 0.6, 0.2, 0.3), NA),
       border = c("white", NA),
       col = c(NA, "firebrick"), lwd = c(NA, 2), lty = c(NA, 1),
       bg = "white")

dev.off()


# ============================================================
# Figure 2: Acceptance rate vs c
# ============================================================

pdf("fig_efficiency.pdf", width = 6, height = 4)
par(mar = c(4, 4, 2.5, 1), family = "serif")

c_vals <- seq(c_beta, c_beta * 3, length.out = 30)
eff <- 1 / c_vals

plot(c_vals, eff * 100, type = "b", pch = 19, cex = 0.7,
     col = "steelblue", lwd = 1.5,
     xlab = TeX("Envelope constant $c$"),
     ylab = "Acceptance Rate (%)",
     main = TeX("Efficiency: acceptance rate $= 1/c$"))
abline(v = c_beta, lty = 2, col = "firebrick")
text(c_beta + 0.1, max(eff * 100) - 2,
     TeX(sprintf("$c^* = %.2f$", c_beta)),
     col = "firebrick", adj = 0, cex = 0.8)

dev.off()


# ============================================================
# Figure 3: Sequential sampling (step-by-step)
# ============================================================

pdf("fig_steps.pdf", width = 10, height = 5.5)
par(mfrow = c(2, 3), mar = c(3.5, 3.5, 2, 1), family = "serif")

set.seed(7)
steps <- 6
x_seq <- runif(steps)
u_seq <- runif(steps)
acc_seq <- u_seq <= dbeta(x_seq, a, b) / c_beta

for (k in 1:steps) {
  plot(NULL, xlim = c(0, 1), ylim = c(0, c_beta + 0.2),
       xlab = "x", ylab = "",
       main = paste0("Step ", k,
                      ifelse(acc_seq[k], " (Accept)", " (Reject)")))
  polygon(c(0, 0, 1, 1), c(0, c_beta, c_beta, 0),
          col = rgb(0.9, 0.9, 0.9, 0.4), border = NA)
  lines(xgrid, dbeta(xgrid, a, b), col = "firebrick", lwd = 2)
  abline(h = c_beta, col = "steelblue", lwd = 1.5, lty = 2)

  if (k > 1) {
    prev_cols <- ifelse(acc_seq[1:(k-1)],
                        rgb(0.2, 0.6, 0.2, 0.4),
                        rgb(0.8, 0.2, 0.2, 0.3))
    points(x_seq[1:(k-1)], u_seq[1:(k-1)] * c_beta,
           pch = 20, cex = 1.2, col = prev_cols)
  }

  pt_col <- ifelse(acc_seq[k], "green4", "red3")
  points(x_seq[k], u_seq[k] * c_beta,
         pch = 8, cex = 2, col = pt_col, lwd = 2)
  arrows(x_seq[k], 0, x_seq[k], u_seq[k] * c_beta,
         length = 0.08, col = pt_col, lty = 3)
}

dev.off()

cat("Figures saved: fig_envelopes.pdf, fig_efficiency.pdf, fig_steps.pdf\n")


# ============================================================
# Figure 4: 2D Contour — rejection sampling in two dimensions
# Target: bivariate normal N([0,0], I)
# Proposal: bivariate uniform on [-3,3]^2 (bounding box)
# ============================================================

pdf("fig_contour_2d.pdf", width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), family = "serif")

set.seed(42)
n2d <- 5000
x1_prop <- runif(n2d, -3, 3)
x2_prop <- runif(n2d, -3, 3)

# Target: N([0,0], I), evaluated on proposals
f_2d <- dnorm(x1_prop) * dnorm(x2_prop)
g_2d <- 1 / 36  # uniform on [-3,3]^2
c_2d <- dnorm(0)^2 / g_2d  # max of f / g
u_2d <- runif(n2d) 
acc_2d <- u_2d <= f_2d / (c_2d * g_2d) # acceptance condition

cat(sprintf("2D rejection sampling: accepted %d / %d = %.1f%%  (theory: %.1f%%)\n",
            sum(acc_2d), n2d, 100 * mean(acc_2d), 100 / c_2d))

# Panel A: Contour of target + accepted/rejected points
g1 <- seq(-3, 3, length.out = 100)
g2 <- seq(-3, 3, length.out = 100)
z_grid <- outer(g1, g2, function(a, b) dnorm(a) * dnorm(b))

plot(x1_prop, x2_prop, pch = 20, cex = 0.3,
     col = ifelse(acc_2d, rgb(0.2, 0.6, 0.2, 0.4), rgb(0.8, 0.2, 0.2, 0.15)),
     xlab = TeX("$x_1$"), ylab = TeX("$x_2$"),
     main = TeX("(a) 2D Rejection: $N(0, I_2)$ via Uniform"))
contour(g1, g2, z_grid, add = TRUE, col = "firebrick",
        levels = c(0.01, 0.03, 0.05, 0.1, 0.12, 0.15),
        lwd = 1.5, labcex = 0.6)
legend("topright", cex = 0.6, pch = 20,
       legend = c("Accepted", "Rejected"),
       col = c(rgb(0.2, 0.6, 0.2), rgb(0.8, 0.2, 0.2)), bg = "white")

# Panel B: Acceptance rate vs dimension d
dims <- 1:15
acc_rates <- numeric(length(dims))
for (di in seq_along(dims)) {
  d <- dims[di]
  # For N(0,I_d) with Uniform[-3,3]^d proposal:
  # c = (2*pi)^{-d/2} / (1/6)^d = (6/sqrt(2*pi))^d
  # acceptance = 1/c = (sqrt(2*pi)/6)^d
  acc_rates[di] <- (sqrt(2 * pi) / 6)^d
}

plot(dims, acc_rates * 100, type = "b", pch = 19, cex = 0.8,
     col = "steelblue", lwd = 1.5, log = "y",
     xlab = "Dimension d", ylab = "Acceptance Rate (%, log scale)",
     main = "(b) Curse of Dimensionality",
     xaxt = "n")
axis(1, at = dims)
abline(h = 1, col = "firebrick", lty = 2)
text(10, 2, "1% threshold", col = "firebrick", cex = 0.7)

dev.off()


# ============================================================
# Figure 5: Where Rejection Sampling Fails
# Target: Cauchy(0,1) — heavy tails
# Proposal: N(0,1) — light tails
# The envelope condition c*g(x) >= f(x) is IMPOSSIBLE
# because g(x)/f(x) -> 0 as |x| -> infinity
# ============================================================

pdf("fig_rs_failure.pdf", width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), family = "serif")

xf <- seq(-8, 8, length.out = 500)

# Panel A: No finite c exists
plot(xf, dcauchy(xf), type = "l", col = "firebrick", lwd = 2,
     xlab = "x", ylab = "Density",
     main = "(a) Cannot Envelope Cauchy with Normal",
     ylim = c(0, 0.45))
# Try c = 2 (arbitrary)
c_try <- pi / 2  # f(0)/g(0) = (1/pi)/(1/sqrt(2pi)) = sqrt(2pi)/pi
lines(xf, c_try * dnorm(xf), col = "steelblue", lwd = 2, lty = 2)

# Show where envelope fails
fail_region <- abs(xf) > 2.5
polygon(c(xf[fail_region & xf > 0], rev(xf[fail_region & xf > 0])),
        c(dcauchy(xf[fail_region & xf > 0]),
          pmin(c_try * dnorm(xf[fail_region & xf > 0]),
               dcauchy(xf[fail_region & xf > 0]))),
        col = rgb(0.8, 0.2, 0.2, 0.3), border = NA)
polygon(c(xf[fail_region & xf < 0], rev(xf[fail_region & xf < 0])),
        c(dcauchy(xf[fail_region & xf < 0]),
          pmin(c_try * dnorm(xf[fail_region & xf < 0]),
               dcauchy(xf[fail_region & xf < 0]))),
        col = rgb(0.8, 0.2, 0.2, 0.3), border = NA)

legend("topright", cex = 0.65,
       legend = c(TeX("Target: Cauchy(0,1)"),
                  TeX(sprintf("$c \\cdot g(x)$, $c=%.2f$, $g=N(0,1)$", c_try)),
                  "Envelope violation"),
       col = c("firebrick", "steelblue", rgb(0.8, 0.2, 0.2, 0.4)),
       lwd = c(2, 2, NA), lty = c(1, 2, NA), pch = c(NA, NA, 15), bg = "white")

# Panel B: Ratio f(x)/g(x) diverges
ratio_cauchy <- dcauchy(xf) / dnorm(xf)
plot(xf, ratio_cauchy, type = "l", col = "darkorange", lwd = 2,
     xlab = "x", ylab = TeX("$f(x)/g(x)$"),
     main = TeX("(b) Ratio $f(x)/g(x) \\rightarrow \\infty$"),
     ylim = c(0, 15))
abline(h = c_try, col = "steelblue", lty = 2, lwd = 1.5)
text(0, c_try + 0.8, sprintf("c = %.2f", c_try),
     col = "steelblue", cex = 0.7)
text(5, 10, TeX("No finite $c$ works!"), col = "firebrick", cex = 0.85)

dev.off()


cat("Figures saved: fig_contour_2d.pdf, fig_rs_failure.pdf\n")
