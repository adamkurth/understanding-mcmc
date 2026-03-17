library(latex2exp)
set.seed(42)

# ============================================================
# IMPORTANCE SAMPLING: THEORY AND DEMONSTRATIONS
# ============================================================

# ============================================================
# Example 1: Estimating E_p[h(X)] where p is hard to sample
# Target: p(x) = Beta(5, 2) — pretend we can evaluate but not sample
# Proposal: q(x) = Uniform(0,1)
# Goal: estimate E_p[X] = a/(a+b) = 5/7 ≈ 0.7143
# ============================================================

a <- 5; b <- 2
true_mean <- a / (a + b)  # 5/7 ≈ 0.7143

N <- 10000
x_q <- runif(N)  # draw from proposal q(x) = 1
w_raw <- dbeta(x_q, a, b) / dunif(x_q)  # importance weights w(x) = p(x)/q(x)
w_norm <- w_raw / sum(w_raw)             # self-normalized weights

# unnormalized estimator (when normalizing constant is known)
IS_unnorm <- mean(w_raw * x_q)
# self-normalized estimator
IS_selfnorm <- sum(w_norm * x_q)

cat(sprintf("Example 1: E_p[X] where p = Beta(%g, %g)\n", a, b))
cat(sprintf("  True mean:              %.4f\n", true_mean))
cat(sprintf("  Unnormalized IS est:    %.4f\n", IS_unnorm))
cat(sprintf("  Self-normalized IS est: %.4f\n", IS_selfnorm))
cat(sprintf("  Naive MC (from q):      %.4f  (biased — wrong distribution)\n\n",
            mean(x_q)))

# Effective sample size
ESS <- (sum(w_raw))^2 / sum(w_raw^2)
cat(sprintf("  ESS = %.0f / %d (%.1f%%)\n\n", ESS, N, 100 * ESS / N))

# ============================================================
# Example 2: Estimating a tail probability P(X > 4) for X ~ N(0,1)
# Direct MC is very inefficient (rare event)
# Proposal: shifted normal N(4, 1) puts mass in the tail
# ============================================================

threshold <- 4
true_prob <- pnorm(threshold, lower.tail = FALSE)  # ≈ 3.167e-5

M <- 50000

# --- Direct MC ---
x_direct <- rnorm(M)
direct_est <- mean(x_direct > threshold)

# --- IS with shifted Normal proposal ---
mu_q <- 4
x_is <- rnorm(M, mean = mu_q, sd = 1)
w_tail <- dnorm(x_is) / dnorm(x_is, mean = mu_q, sd = 1)
IS_tail <- mean(w_tail * (x_is > threshold))

cat(sprintf("Example 2: P(X > %g) for X ~ N(0,1)\n", threshold))
cat(sprintf("  True probability:  %.6e\n", true_prob))
cat(sprintf("  Direct MC est:     %.6e  (often 0 with small M)\n", direct_est))
cat(sprintf("  IS estimate:       %.6e\n\n", IS_tail))

# ============================================================
# Example 3: Proposal choice matters — good vs bad proposal
# Target: p(x) = N(0,1), estimating E[X^2] = 1
# Good proposal: t_3 (heavier tails, covers p well)
# Bad proposal: N(3,0.5^2) (shifted away, poor coverage)
# ============================================================

K <- 5000
n_rep <- 200  # replications to see variance

est_good <- est_bad <- numeric(n_rep)
ess_good <- ess_bad <- numeric(n_rep)

for (r in 1:n_rep) {
  # --- Good proposal: t_3 ---
  x_g <- rt(K, df = 3)
  w_g <- dnorm(x_g) / dt(x_g, df = 3)
  w_g_norm <- w_g / sum(w_g)
  est_good[r] <- sum(w_g_norm * x_g^2)
  ess_good[r] <- (sum(w_g))^2 / sum(w_g^2)

  # --- Bad proposal: N(3, 0.5^2) ---
  x_b <- rnorm(K, mean = 3, sd = 0.5)
  w_b <- dnorm(x_b) / dnorm(x_b, mean = 3, sd = 0.5)
  w_b_norm <- w_b / sum(w_b)
  est_bad[r] <- sum(w_b_norm * x_b^2)
  ess_bad[r] <- (sum(w_b))^2 / sum(w_b^2)
}

cat("Example 3: E[X^2] where X ~ N(0,1), true value = 1\n")
cat(sprintf("  Good proposal (t_3):    mean=%.3f, sd=%.4f, avg ESS=%.0f\n",
            mean(est_good), sd(est_good), mean(ess_good)))
cat(sprintf("  Bad proposal (N(3,.5)): mean=%.3f, sd=%.4f, avg ESS=%.0f\n\n",
            mean(est_bad), sd(est_bad), mean(ess_bad)))


# ============================================================
# FIGURES
# ============================================================

xgrid <- seq(0, 1, length.out = 500)
xgrid2 <- seq(-6, 8, length.out = 500)

# ============================================================
# Figure 1: Importance Sampling Mechanics (2x2 panel)
# ============================================================

pdf("fig_is_mechanics.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Target p(x) vs Proposal q(x)
plot(xgrid, dbeta(xgrid, a, b), type = "l", col = "firebrick", lwd = 2,
     xlab = "x", ylab = "Density",
     main = TeX("(a) Target $p(x)$ vs Proposal $q(x)$"),
     ylim = c(0, 3.5))
lines(xgrid, dunif(xgrid), col = "steelblue", lwd = 2, lty = 2)
legend("topleft", cex = 0.7,
       legend = c(TeX("$p(x) = \\mathrm{Beta}(5,2)$"),
                  TeX("$q(x) = \\mathrm{Uniform}(0,1)$")),
       col = c("firebrick", "steelblue"), lwd = 2, lty = c(1, 2), bg = "white")

# Panel B: Importance weights w(x) = p(x)/q(x)
plot(xgrid, dbeta(xgrid, a, b) / 1, type = "l", col = "darkorange", lwd = 2,
     xlab = "x", ylab = "Weight w(x)",
     main = TeX("(b) Importance Weights $w(x) = p(x)/q(x)$"))
abline(h = 1, col = "gray50", lty = 3)
legend("topleft", cex = 0.7,
       legend = c(TeX("$w(x) = \\mathrm{Beta}(5,2)/1$"), "w = 1 (no reweighting)"),
       col = c("darkorange", "gray50"), lwd = c(2, 1), lty = c(1, 3), bg = "white")

# Panel C: Weighted samples scatter
n_show <- 500
plot(x_q[1:n_show], w_raw[1:n_show], pch = 20, cex = 0.6,
     col = rgb(0.2, 0.4, 0.8, 0.5),
     xlab = "x", ylab = "Weight",
     main = "(c) Samples Colored by Weight")
abline(v = true_mean, col = "firebrick", lwd = 2, lty = 2)
text(true_mean - 0.04, max(w_raw[1:n_show]) * 0.95,
     TeX(sprintf("$E[X] = %.3f$", true_mean)),
     col = "firebrick", adj = 1, cex = 0.8)

# Panel D: Self-normalized weights as discrete distribution
ord <- order(x_q[1:n_show])
sub_x <- x_q[1:n_show][ord]
sub_w <- w_norm[1:n_show][ord]
plot(sub_x, sub_w, type = "h", col = rgb(0.2, 0.6, 0.2, 0.6),
     lwd = 0.8,
     xlab = "x", ylab = "Normalized weight",
     main = TeX("(d) Self-Normalized Weights $\\tilde{w}_i$"))
lines(xgrid, dbeta(xgrid, a, b) / N, col = "firebrick", lwd = 2)
legend("topleft", cex = 0.7,
       legend = c(TeX("$\\tilde{w}_i$"), TeX("$p(x)/N$ (scaled)")),
       col = c(rgb(0.2, 0.6, 0.2), "firebrick"), lwd = c(1, 2),
       lty = c(1, 1), bg = "white")

dev.off()


# ============================================================
# Figure 2: Tail probability — direct MC vs IS
# ============================================================

pdf("fig_is_tail.pdf", width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: The rare-event problem
xg <- seq(-4, 7, length.out = 500)
plot(xg, dnorm(xg), type = "l", col = "firebrick", lwd = 2,
     xlab = "x", ylab = "Density",
     main = TeX("(a) Rare Event: $P(X > 4)$"),
     ylim = c(0, 0.5))
lines(xg, dnorm(xg, mean = mu_q), col = "steelblue", lwd = 2, lty = 2)
polygon(c(threshold, xg[xg >= threshold], max(xg)),
        c(0, dnorm(xg[xg >= threshold]), 0),
        col = rgb(0.8, 0.2, 0.2, 0.2), border = NA)
abline(v = threshold, col = "gray40", lty = 3)
legend("topright", cex = 0.7,
       legend = c(TeX("Target $p(x) = N(0,1)$"),
                  TeX(sprintf("Proposal $q(x) = N(%d,1)$", mu_q)),
                  TeX("$P(X > 4)$ region")),
       col = c("firebrick", "steelblue", rgb(0.8, 0.2, 0.2, 0.4)),
       lwd = c(2, 2, NA), lty = c(1, 2, NA),
       pch = c(NA, NA, 15), bg = "white")

# Panel B: Convergence comparison
set.seed(99)
ns <- seq(100, 50000, by = 200)
direct_conv <- is_conv <- numeric(length(ns))
for (j in seq_along(ns)) {
  nn <- ns[j]
  direct_conv[j] <- mean(rnorm(nn) > threshold)
  x_tmp <- rnorm(nn, mean = mu_q, sd = 1)
  w_tmp <- dnorm(x_tmp) / dnorm(x_tmp, mean = mu_q, sd = 1)
  is_conv[j] <- mean(w_tmp * (x_tmp > threshold))
}

plot(ns, is_conv, type = "l", col = "steelblue", lwd = 1.5,
     xlab = "Number of samples", ylab = TeX("$\\hat{P}(X > 4)$"),
     main = "(b) Convergence: Direct MC vs IS",
     ylim = c(0, max(is_conv) * 2.5))
lines(ns, direct_conv, col = "firebrick", lwd = 1.5)
abline(h = true_prob, col = "black", lty = 2, lwd = 1.5)
legend("topright", cex = 0.7,
       legend = c("IS (N(4,1) proposal)", "Direct MC", "True value"),
       col = c("steelblue", "firebrick", "black"),
       lwd = c(1.5, 1.5, 1.5), lty = c(1, 1, 2), bg = "white")

dev.off()


# ============================================================
# Figure 3: Good vs Bad Proposal — weight distributions & ESS
# ============================================================

pdf("fig_is_proposal_choice.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Densities
xg3 <- seq(-5, 7, length.out = 500)
plot(xg3, dnorm(xg3), type = "l", col = "firebrick", lwd = 2,
     xlab = "x", ylab = "Density",
     main = "(a) Target and Two Proposals",
     ylim = c(0, 0.85))
lines(xg3, dt(xg3, df = 3), col = "steelblue", lwd = 2, lty = 2)
lines(xg3, dnorm(xg3, mean = 3, sd = 0.5), col = "darkorange", lwd = 2, lty = 4)
legend("topright", cex = 0.65,
       legend = c(TeX("Target: $p(x) = N(0,1)$"),
                  TeX("Good: $q(x) = t_3$"),
                  TeX("Bad: $q(x) = N(3, 0.25)$")),
       col = c("firebrick", "steelblue", "darkorange"),
       lwd = 2, lty = c(1, 2, 4), bg = "white")

# Panel B: Histogram of weight distribution — good proposal
set.seed(42)
x_g_one <- rt(K, df = 3)
w_g_one <- dnorm(x_g_one) / dt(x_g_one, df = 3)
x_b_one <- rnorm(K, mean = 3, sd = 0.5)
w_b_one <- dnorm(x_b_one) / dnorm(x_b_one, mean = 3, sd = 0.5)

hist(w_g_one, breaks = 60, freq = FALSE, col = rgb(0.27, 0.51, 0.71, 0.4),
     border = "white", xlim = c(0, quantile(w_g_one, 0.99)),
     main = TeX("(b) Weight Distribution: Good ($t_3$)"),
     xlab = "Weight w(x)", ylab = "Density")
abline(v = mean(w_g_one), col = "steelblue", lwd = 2, lty = 2)

# Panel C: Weight distribution — bad proposal (log scale)
hist(log10(w_b_one + 1e-300), breaks = 60, freq = FALSE,
     col = rgb(0.93, 0.49, 0.19, 0.4), border = "white",
     main = TeX("(c) Log-Weight Distribution: Bad ($N(3, 0.25)$)"),
     xlab = TeX("$\\log_{10}(w)$"), ylab = "Density")

# Panel D: Repeated estimation — boxplot comparison
boxplot(list(Good = est_good, Bad = est_bad),
        col = c(rgb(0.27, 0.51, 0.71, 0.4), rgb(0.93, 0.49, 0.19, 0.4)),
        main = TeX("(d) IS Estimates of $E[X^2]$ (200 reps)"),
        ylab = TeX("$\\hat{\\mu}$"), outline = FALSE)
abline(h = 1, col = "firebrick", lwd = 2, lty = 2)
text(1.5, 1.05, "True value = 1", col = "firebrick", cex = 0.75)

dev.off()


# ============================================================
# Figure 4: Effective Sample Size vs proposal mismatch
# ============================================================

pdf("fig_is_ess.pdf", width = 7, height = 4.5)
par(mar = c(4, 4, 2.5, 1), family = "serif")

mu_vals <- seq(0, 6, length.out = 50)
ess_vals <- numeric(length(mu_vals))

set.seed(42)
n_ess <- 5000
for (j in seq_along(mu_vals)) {
  x_tmp <- rnorm(n_ess, mean = mu_vals[j], sd = 1)
  w_tmp <- dnorm(x_tmp) / dnorm(x_tmp, mean = mu_vals[j], sd = 1)
  ess_vals[j] <- (sum(w_tmp))^2 / sum(w_tmp^2)
}

plot(mu_vals, ess_vals / n_ess * 100, type = "b", pch = 19, cex = 0.6,
     col = "steelblue", lwd = 1.5,
     xlab = TeX("Proposal mean $\\mu_q$"),
     ylab = "ESS / N (%)",
     main = TeX("ESS vs Proposal Mismatch: target $N(0,1)$, proposal $N(\\mu_q, 1)$"))
abline(v = 0, col = "firebrick", lty = 2, lwd = 1.5)
text(0.3, max(ess_vals / n_ess * 100) - 3,
     TeX("$\\mu_q = 0$ (perfect)"), col = "firebrick", adj = 0, cex = 0.75)

dev.off()


# ============================================================
# Figure 5: Step-by-step IS visualization (6 panels)
# ============================================================

pdf("fig_is_steps.pdf", width = 10, height = 5.5)
par(mfrow = c(2, 3), mar = c(3.5, 3.5, 2.5, 1), family = "serif")

set.seed(7)
steps <- 6
x_step <- runif(steps)
w_step <- dbeta(x_step, a, b)  # since q = Uniform, w = p(x)/1

# Running weighted mean
running_est <- numeric(steps)
for (k in 1:steps) {
  running_est[k] <- sum(w_step[1:k] * x_step[1:k]) / sum(w_step[1:k])
}

for (k in 1:steps) {
  plot(NULL, xlim = c(0, 1), ylim = c(0, 3.5),
       xlab = "x", ylab = "Density",
       main = bquote("Step" ~ .(k) ~ ": " ~ hat(mu) == .(sprintf("%.3f", running_est[k]))))
  lines(xgrid, dbeta(xgrid, a, b), col = "firebrick", lwd = 2)
  abline(h = 1, col = "steelblue", lwd = 1.5, lty = 2)

  if (k > 1) {
    for (j in 1:(k - 1)) {
      pt_size <- 0.5 + 2 * w_step[j] / max(dbeta(xgrid, a, b))
      points(x_step[j], dbeta(x_step[j], a, b),
             pch = 20, cex = pt_size,
             col = rgb(0.4, 0.4, 0.4, 0.4))
    }
  }

  # current sample: star, size proportional to weight
  pt_size_k <- 0.5 + 2 * w_step[k] / max(dbeta(xgrid, a, b))
  points(x_step[k], dbeta(x_step[k], a, b),
         pch = 8, cex = 1.5 + pt_size_k * 0.5, col = "darkorange", lwd = 2)
  arrows(x_step[k], 0, x_step[k], dbeta(x_step[k], a, b),
         length = 0.08, col = "darkorange", lty = 3)

  # show weight text
  text(x_step[k], dbeta(x_step[k], a, b) + 0.25,
       sprintf("w=%.2f", w_step[k]),
       col = "darkorange", cex = 0.7)

  # running mean line vs true mean
  abline(v = running_est[k], col = "green4", lwd = 1.5, lty = 2)
  abline(v = true_mean, col = "firebrick", lwd = 1, lty = 3)
}

dev.off()


cat("Figures saved:\n")
cat("  fig_is_mechanics.pdf\n")
cat("  fig_is_tail.pdf\n")
cat("  fig_is_proposal_choice.pdf\n")
cat("  fig_is_ess.pdf\n")
cat("  fig_is_steps.pdf\n")


# ============================================================
# Figure 6: 2D Contour — IS weight landscape
# Target: bivariate N([1,1], [[1,.7],[.7,1]])
# Proposal: bivariate N([0,0], 2*I)
# Contour shows target; points sized/colored by weight
# ============================================================

pdf("fig_is_contour_2d.pdf", width = 9, height = 4.5)
# Use layout() so filled.contour and regular plot share one page
# Columns: panel A | panel B plot | panel B color bar
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(4, 3.5, 1.5))
par(family = "serif")

set.seed(42)
n2d <- 2000

# Target: bivariate normal with correlation
mu_t <- c(1, 1)
rho <- 0.7
Sigma_t <- matrix(c(1, rho, rho, 1), 2, 2)
Sigma_t_inv <- solve(Sigma_t)
det_Sigma_t <- det(Sigma_t)

# Evaluate bivariate normal density manually (avoid extra packages)
dbvnorm <- function(x1, x2, mu, Sinv, detS) {
  dx <- cbind(x1 - mu[1], x2 - mu[2])
  quad <- rowSums((dx %*% Sinv) * dx)
  exp(-0.5 * quad) / (2 * pi * sqrt(detS))
}

# Proposal: N([0,0], 2*I)
sd_q <- sqrt(2)
x1_q <- rnorm(n2d, 0, sd_q)
x2_q <- rnorm(n2d, 0, sd_q)

# Weights
f_vals <- dbvnorm(x1_q, x2_q, mu_t, Sigma_t_inv, det_Sigma_t)
q_vals <- dnorm(x1_q, 0, sd_q) * dnorm(x2_q, 0, sd_q)
w_2d <- f_vals / q_vals
w_2d_norm <- w_2d / sum(w_2d)

ESS_2d <- (sum(w_2d))^2 / sum(w_2d^2)
cat(sprintf("\n2D IS contour: ESS = %.0f / %d (%.1f%%)\n", ESS_2d, n2d, 100 * ESS_2d / n2d))

# Panel A: Contour of target + weighted samples
g1 <- seq(-3, 5, length.out = 100)
g2 <- seq(-3, 5, length.out = 100)
z_target <- outer(g1, g2, function(a, b) dbvnorm(a, b, mu_t, Sigma_t_inv, det_Sigma_t))

# Color by weight quantile
w_q <- cut(w_2d, quantile(w_2d, probs = c(0, 0.5, 0.8, 0.95, 1)),
           include.lowest = TRUE, labels = FALSE)
pt_cols <- c(rgb(0.6, 0.6, 0.6, 0.3), rgb(0.3, 0.5, 0.8, 0.5),
             rgb(0.1, 0.3, 0.8, 0.7), rgb(0.8, 0.2, 0.1, 0.9))
pt_sizes <- c(0.3, 0.5, 0.9, 1.5)

par(mar = c(4, 4, 2.5, 1))
plot(x1_q, x2_q, pch = 20, cex = pt_sizes[w_q], col = pt_cols[w_q],
     xlab = TeX("$x_1$"), ylab = TeX("$x_2$"),
     main = "(a) IS Weights on 2D Contour",
     xlim = c(-3, 5), ylim = c(-3, 5))
contour(g1, g2, z_target, add = TRUE, col = "firebrick",
        levels = c(0.01, 0.03, 0.06, 0.1, 0.14),
        lwd = 1.5, labcex = 0.5)
# proposal contour
z_prop <- outer(g1, g2, function(a, b) dnorm(a, 0, sd_q) * dnorm(b, 0, sd_q))
contour(g1, g2, z_prop, add = TRUE, col = "steelblue", lty = 2,
        levels = c(0.01, 0.03, 0.05), lwd = 1, labcex = 0.5)
legend("topright", cex = 0.55,
       legend = c("Target contour", "Proposal contour",
                  "Low weight", "Mid weight", "High weight", "Very high weight"),
       col = c("firebrick", "steelblue", pt_cols),
       lwd = c(1.5, 1, NA, NA, NA, NA), lty = c(1, 2, NA, NA, NA, NA),
       pch = c(NA, NA, 20, 20, 20, 20),
       pt.cex = c(NA, NA, pt_sizes), bg = "white")

# Panel B: Weight heatmap using image() + manual color bar
z_weight <- outer(g1, g2, function(a, b) {
  f_v <- dbvnorm(a, b, mu_t, Sigma_t_inv, det_Sigma_t)
  q_v <- dnorm(a, 0, sd_q) * dnorm(b, 0, sd_q)
  ifelse(q_v > 1e-10, f_v / q_v, 0)
})

n_cols <- 20
heat_pal <- hcl.colors(n_cols, "YlOrRd", rev = TRUE)

par(mar = c(4, 4, 2.5, 1))
image(g1, g2, z_weight, col = heat_pal,
      xlab = TeX("$x_1$"), ylab = TeX("$x_2$"),
      main = "(b) Weight Surface w(x) = p(x)/q(x)")
contour(g1, g2, z_target, add = TRUE, col = "black",
        levels = c(0.03, 0.1, 0.14), lwd = 1, labcex = 0.5)

# Color bar
par(mar = c(4, 0.5, 2.5, 2.5))
z_range <- range(z_weight)
bar_seq <- seq(z_range[1], z_range[2], length.out = n_cols + 1)
image(1, bar_seq[-length(bar_seq)],
      matrix(bar_seq[-length(bar_seq)], nrow = 1),
      col = heat_pal, xaxt = "n", yaxt = "n",
      xlab = "", ylab = "")
axis(4, las = 1, cex.axis = 0.7)

dev.off()


# ============================================================
# Figure 7: Where Importance Sampling Fails — infinite variance
# Target: Cauchy(0,1), proposal: N(0,1)
# q has lighter tails than p => w(x) = p(x)/q(x) has infinite variance
# Self-normalized estimates jump wildly
# ============================================================

pdf("fig_is_failure.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Tail behavior
xf <- seq(-10, 10, length.out = 500)
plot(xf, dcauchy(xf), type = "l", col = "firebrick", lwd = 2,
     xlab = "x", ylab = "Density",
     main = "(a) Cauchy Target vs Normal Proposal",
     ylim = c(0, 0.42))
lines(xf, dnorm(xf), col = "steelblue", lwd = 2, lty = 2)
legend("topright", cex = 0.65,
       legend = c("Target: Cauchy(0,1)", "Proposal: N(0,1)"),
       col = c("firebrick", "steelblue"), lwd = 2, lty = c(1, 2), bg = "white")

# Panel B: Weight ratio diverges in tails
w_ratio <- dcauchy(xf) / dnorm(xf)
plot(xf, w_ratio, type = "l", col = "darkorange", lwd = 2,
     xlab = "x", ylab = "w(x) = p(x)/q(x)",
     main = TeX("(b) Weights Diverge: $w(x) \\rightarrow \\infty$"),
     ylim = c(0, 30))
abline(h = 1, col = "gray50", lty = 3)
text(6, 20, TeX("$\\mathrm{Var}(w) = \\infty$"), col = "firebrick", cex = 0.85)

# Panel C: Single run — a few huge weights dominate
set.seed(123)
n_fail <- 10000
x_fail <- rnorm(n_fail)
w_fail <- dcauchy(x_fail) / dnorm(x_fail)
w_fail_norm <- w_fail / sum(w_fail)

plot(x_fail, w_fail_norm, pch = 20, cex = 0.3,
     col = rgb(0.2, 0.4, 0.8, 0.4),
     xlab = "x", ylab = "Normalized weight",
     main = "(c) Weight Concentration",
     ylim = c(0, max(w_fail_norm) * 1.1))
# Annotate the biggest weights
big_idx <- order(w_fail_norm, decreasing = TRUE)[1:3]
points(x_fail[big_idx], w_fail_norm[big_idx],
       pch = 1, cex = 2, col = "firebrick", lwd = 2)
ess_fail <- (sum(w_fail))^2 / sum(w_fail^2)
text(mean(range(x_fail)), max(w_fail_norm) * 0.9,
     sprintf("ESS = %.0f / %d", ess_fail, n_fail),
     col = "firebrick", cex = 0.8)

# Panel D: Repeated runs — estimates of E[X^2] jump wildly
# True E_Cauchy[X^2] is undefined (infinite), but the self-normalized
# estimator will give finite but wildly variable values
n_reps_fail <- 300
est_fail <- numeric(n_reps_fail)
ess_fail_all <- numeric(n_reps_fail)
for (r in 1:n_reps_fail) {
  xx <- rnorm(2000)
  ww <- dcauchy(xx) / dnorm(xx)
  ww_n <- ww / sum(ww)
  est_fail[r] <- sum(ww_n * xx^2)
  ess_fail_all[r] <- (sum(ww))^2 / sum(ww^2)
}

plot(1:n_reps_fail, est_fail, pch = 20, cex = 0.4,
     col = rgb(0.2, 0.4, 0.8, 0.5),
     xlab = "Replication", ylab = TeX("IS est. of $E[X^2]$"),
     main = TeX("(d) Unstable Estimates ($E_{Cauchy}[X^2] = \\infty$)"))
abline(h = median(est_fail), col = "darkorange", lty = 2, lwd = 1.5)
text(n_reps_fail * 0.7, max(est_fail) * 0.85,
     sprintf("Avg ESS = %.0f / 2000", mean(ess_fail_all)),
     col = "firebrick", cex = 0.7)

dev.off()


# ============================================================
# SAMPLING IMPORTANCE RESAMPLING (SIR)
# ============================================================

cat("\n=== Sampling Importance Resampling (SIR) ===\n")

# ============================================================
# SIR Example 1: Beta(5,2) posterior from Uniform proposal
# Draw N proposal samples, weight them, resample M << N
# to get approximately i.i.d. draws from the target
# ============================================================

set.seed(42)
N_sir <- 10000
M_sir <- 1000  # resample this many

x_sir <- runif(N_sir)
w_sir <- dbeta(x_sir, a, b)  # q = Uniform => w = p(x)
w_sir_norm <- w_sir / sum(w_sir)

# Resample with replacement, proportional to weights
idx_sir <- sample(1:N_sir, size = M_sir, replace = TRUE, prob = w_sir_norm)
x_resampled <- x_sir[idx_sir]

cat(sprintf("SIR Example 1: Beta(5,2), N=%d proposals, M=%d resamples\n", N_sir, M_sir))
cat(sprintf("  Weighted IS mean:    %.4f\n", sum(w_sir_norm * x_sir)))
cat(sprintf("  Resampled mean:      %.4f\n", mean(x_resampled)))
cat(sprintf("  True mean:           %.4f\n\n", true_mean))


# ============================================================
# SIR Example 2: Bimodal target from broad proposal
# Target: 0.3*N(-2,0.5^2) + 0.7*N(2,0.8^2)
# Proposal: N(0, 3^2) — broad, covers both modes
# ============================================================

target_bimodal <- function(x) {
  0.3 * dnorm(x, -2, 0.5) + 0.7 * dnorm(x, 2, 0.8)
}

set.seed(42)
N_sir2 <- 20000
M_sir2 <- 2000
sd_prop <- 3

x_sir2 <- rnorm(N_sir2, 0, sd_prop)
w_sir2 <- target_bimodal(x_sir2) / dnorm(x_sir2, 0, sd_prop)
w_sir2_norm <- w_sir2 / sum(w_sir2)

idx_sir2 <- sample(1:N_sir2, size = M_sir2, replace = TRUE, prob = w_sir2_norm)
x_resamp2 <- x_sir2[idx_sir2]

ESS_sir2 <- (sum(w_sir2))^2 / sum(w_sir2^2)
cat(sprintf("SIR Example 2: Bimodal mixture, N=%d, M=%d, ESS=%.0f\n",
            N_sir2, M_sir2, ESS_sir2))


# ============================================================
# Figure 8: SIR Demonstration (2x2)
# ============================================================

pdf("fig_sir.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Beta(5,2) — weighted vs resampled histogram
hist(x_resampled, breaks = 40, freq = FALSE,
     col = rgb(0.27, 0.51, 0.71, 0.4), border = "white",
     main = "(a) SIR: Resampled Draws vs Target",
     xlab = "x", ylab = "Density", xlim = c(0, 1))
lines(xgrid, dbeta(xgrid, a, b), col = "firebrick", lwd = 2)
legend("topleft", cex = 0.65,
       legend = c(sprintf("SIR resamples (M=%d)", M_sir),
                  TeX("$p(x) = \\mathrm{Beta}(5,2)$")),
       fill = c(rgb(0.27, 0.51, 0.71, 0.4), NA),
       border = c("white", NA),
       col = c(NA, "firebrick"), lwd = c(NA, 2), bg = "white")

# Panel B: Bimodal — proposal vs target vs resampled
xg_sir <- seq(-6, 6, length.out = 500)
plot(xg_sir, target_bimodal(xg_sir), type = "l", col = "firebrick", lwd = 2,
     xlab = "x", ylab = "Density",
     main = "(b) SIR: Bimodal Mixture Target",
     ylim = c(0, 0.4))
lines(xg_sir, dnorm(xg_sir, 0, sd_prop), col = "steelblue", lwd = 2, lty = 2)
hist(x_resamp2, breaks = 50, freq = FALSE, add = TRUE,
     col = rgb(0.27, 0.51, 0.71, 0.2), border = rgb(0.27, 0.51, 0.71, 0.5))
legend("topright", cex = 0.55,
       legend = c("Target (mixture)", TeX("Proposal $N(0, 9)$"), "SIR resamples"),
       col = c("firebrick", "steelblue", rgb(0.27, 0.51, 0.71, 0.5)),
       lwd = c(2, 2, 1), lty = c(1, 2, 1), bg = "white")

# Panel C: Weights before resampling — bimodal case
ord2 <- order(x_sir2[1:2000])
plot(x_sir2[1:2000][ord2], w_sir2_norm[1:2000][ord2], type = "h",
     col = rgb(0.2, 0.6, 0.2, 0.5), lwd = 0.5,
     xlab = "x", ylab = "Normalized weight",
     main = TeX("(c) Weights Before Resampling"))
lines(xg_sir, target_bimodal(xg_sir) / N_sir2, col = "firebrick", lwd = 2)

# Panel D: IS weighted ECDF vs resampled ECDF vs true CDF
# Use Beta(5,2) example for clarity
plot(ecdf(x_resampled), col = "steelblue", lwd = 1.5,
     main = "(d) CDF: Resampled vs True",
     xlab = "x", ylab = "CDF", do.points = FALSE)
lines(xgrid, pbeta(xgrid, a, b), col = "firebrick", lwd = 2, lty = 2)
# weighted ECDF from original IS
x_sorted <- sort(x_sir)
w_sorted <- w_sir_norm[order(x_sir)]
lines(x_sorted, cumsum(w_sorted), col = "darkorange", lwd = 1.5, lty = 3)
legend("bottomright", cex = 0.6,
       legend = c("SIR resampled ECDF", "True CDF", "IS weighted ECDF"),
       col = c("steelblue", "firebrick", "darkorange"),
       lwd = c(1.5, 2, 1.5), lty = c(1, 2, 3), bg = "white")

dev.off()


# ============================================================
# Figure 9: SIR — effect of resample size M
# ============================================================

pdf("fig_sir_convergence.pdf", width = 7, height = 4.5)
par(mar = c(4, 4, 2.5, 1), family = "serif")

set.seed(42)
M_vals <- c(50, 100, 200, 500, 1000, 2000, 5000)
n_reps_sir <- 200
sir_sd <- matrix(NA, length(M_vals), n_reps_sir)

for (j in seq_along(M_vals)) {
  Mj <- M_vals[j]
  for (r in 1:n_reps_sir) {
    x_prop_sir <- runif(N_sir)
    w_prop_sir <- dbeta(x_prop_sir, a, b)
    w_prop_norm <- w_prop_sir / sum(w_prop_sir)
    idx_tmp <- sample(1:N_sir, size = Mj, replace = TRUE, prob = w_prop_norm)
    sir_sd[j, r] <- mean(x_prop_sir[idx_tmp])
  }
}

sir_means <- rowMeans(sir_sd)
sir_sds <- apply(sir_sd, 1, sd)

plot(M_vals, sir_sds, type = "b", pch = 19, cex = 0.8,
     col = "steelblue", lwd = 1.5, log = "x",
     xlab = "Resample size M (log scale)", ylab = TeX("SD of $\\hat{\\mu}$ (200 reps)"),
     main = "SIR Precision vs Resample Size M")
abline(h = 0, col = "gray70")

dev.off()


cat("New figures saved:\n")
cat("  fig_is_contour_2d.pdf\n")
cat("  fig_is_failure.pdf\n")
cat("  fig_sir.pdf\n")
cat("  fig_sir_convergence.pdf\n")
