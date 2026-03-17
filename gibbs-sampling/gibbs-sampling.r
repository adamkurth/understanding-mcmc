library(latex2exp)
set.seed(42)

# ============================================================
# GIBBS SAMPLING: THEORY AND DEMONSTRATIONS
# ============================================================
# Helper: Gibbs sampler for bivariate normal
# Target: (theta1, theta2) ~ N(mu, Sigma) with correlation rho
# Conditionals:
#   theta1 | theta2 ~ N(mu1 + rho*(sigma1/sigma2)*(theta2-mu2), sigma1^2*(1-rho^2))
#   theta2 | theta1 ~ N(mu2 + rho*(sigma2/sigma1)*(theta1-mu1), sigma2^2*(1-rho^2))
# ============================================================

gibbs_bivariate <- function(n_iter, rho, mu1 = 0, mu2 = 0,
                            sigma1 = 1, sigma2 = 1, init = c(-2, -2)) {
  theta <- matrix(NA, n_iter, 2)
  theta[1, ] <- init
  cond_sd1 <- sigma1 * sqrt(1 - rho^2)
  cond_sd2 <- sigma2 * sqrt(1 - rho^2)
  for (t in 2:n_iter) {
    # draw theta1 | theta2
    cond_mean1 <- mu1 + rho * (sigma1 / sigma2) * (theta[t - 1, 2] - mu2)
    theta[t, 1] <- rnorm(1, cond_mean1, cond_sd1)
    # draw theta2 | theta1
    cond_mean2 <- mu2 + rho * (sigma2 / sigma1) * (theta[t, 1] - mu1)
    theta[t, 2] <- rnorm(1, cond_mean2, cond_sd2)
  }
  theta
}

# ============================================================
# Figure 1: Bivariate Normal — Gibbs Mechanics (2x2)
# ============================================================

cat("Generating Figure 1: Bivariate normal mechanics...\n")

rho1 <- 0.7
n1 <- 5000
samp1 <- gibbs_bivariate(n1, rho = rho1)

# Compute contour grid for bivariate normal
bvn_density <- function(x1, x2, rho) {
  z <- x1^2 - 2 * rho * x1 * x2 + x2^2
  exp(-z / (2 * (1 - rho^2))) / (2 * pi * sqrt(1 - rho^2))
}

g <- seq(-3.5, 3.5, length.out = 150)
z_bvn <- outer(g, g, function(a, b) bvn_density(a, b, rho1))

pdf("plots/fig_gibbs_bivariate.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Staircase trace on joint contour
plot(NULL, xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5),
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"),
     main = TeX("(a) Gibbs Staircase on Joint ($\\rho = 0.7$)"))
contour(g, g, z_bvn, col = "gray70", lwd = 0.8, drawlabels = FALSE, add = TRUE,
        levels = c(0.01, 0.03, 0.06, 0.1, 0.14))
# Draw staircase for first 30 iterations
n_show <- 30
for (t in 2:n_show) {
  # horizontal move: theta1 changes
  segments(samp1[t - 1, 1], samp1[t - 1, 2], samp1[t, 1], samp1[t - 1, 2],
           col = rgb(0.13, 0.55, 0.13, 0.6), lwd = 1.2)
  # vertical move: theta2 changes
  segments(samp1[t, 1], samp1[t - 1, 2], samp1[t, 1], samp1[t, 2],
           col = rgb(0.27, 0.51, 0.71, 0.6), lwd = 1.2)
}
points(samp1[1:n_show, 1], samp1[1:n_show, 2], pch = 20, cex = 0.6,
       col = rgb(0.2, 0.2, 0.2, 0.7))
points(samp1[1, 1], samp1[1, 2], pch = 17, cex = 1.5, col = "firebrick")
legend("topleft", cex = 0.6,
       legend = c(TeX("Update $\\theta_1 \mid \\theta_2$"),
                  TeX("Update $\\theta_2 \mid \\theta_1$"),
                  "Start"),
       col = c(rgb(0.13, 0.55, 0.13), rgb(0.27, 0.51, 0.71), "firebrick"),
       lwd = c(1.5, 1.5, NA), lty = c(1, 1, NA),
       pch = c(NA, NA, 17), bg = "white")

# Panel B: Trace plot for theta1
plot(1:n1, samp1[, 1], type = "l", col = rgb(0.13, 0.55, 0.13, 0.4),
     lwd = 0.5, xlab = "Iteration", ylab = TeX("$\\theta_1$"),
     main = TeX("(b) Trace Plot: $\\theta_1$"))
abline(h = 0, col = "firebrick", lty = 2, lwd = 1.5)

# Panel C: Trace plot for theta2
plot(1:n1, samp1[, 2], type = "l", col = rgb(0.27, 0.51, 0.71, 0.4),
     lwd = 0.5, xlab = "Iteration", ylab = TeX("$\\theta_2$"),
     main = TeX("(c) Trace Plot: $\\theta_2$"))
abline(h = 0, col = "firebrick", lty = 2, lwd = 1.5)

# Panel D: Marginal histograms vs true density
burn <- 500
xg <- seq(-3.5, 3.5, length.out = 300)
hist(samp1[burn:n1, 1], breaks = 50, freq = FALSE,
     col = rgb(0.13, 0.55, 0.13, 0.3), border = "white",
     main = "(d) Marginal Posteriors vs Truth",
     xlab = TeX("$\\theta$"), ylab = "Density",
     xlim = c(-3.5, 3.5), ylim = c(0, 0.5))
hist(samp1[burn:n1, 2], breaks = 50, freq = FALSE,
     col = rgb(0.27, 0.51, 0.71, 0.3), border = "white", add = TRUE)
lines(xg, dnorm(xg), col = "firebrick", lwd = 2, lty = 2)
legend("topright", cex = 0.6,
       legend = c(TeX("$\\theta_1$ samples"), TeX("$\\theta_2$ samples"),
                  TeX("True $N(0,1)$")),
       fill = c(rgb(0.13, 0.55, 0.13, 0.3), rgb(0.27, 0.51, 0.71, 0.3), NA),
       border = c("white", "white", NA),
       col = c(NA, NA, "firebrick"), lwd = c(NA, NA, 2), lty = c(NA, NA, 2),
       bg = "white")

dev.off()


# ============================================================
# Figure 2: Step-by-step Gibbs (6 panels) — bivariate normal
# ============================================================

cat("Generating Figure 2: Step-by-step visualization...\n")

set.seed(7)
samp_step <- gibbs_bivariate(8, rho = rho1, init = c(-2.5, 2))

pdf("plots/fig_gibbs_steps.pdf", width = 10, height = 5.5)
par(mfrow = c(2, 3), mar = c(3.5, 3.5, 2.5, 1), family = "serif")

for (k in 1:6) {
  t <- k + 1  # iteration index (step 1 = iteration 2)
  plot(NULL, xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5),
       xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"),
       main = bquote("Step" ~ .(k)))
  contour(g, g, z_bvn, col = "gray75", lwd = 0.6, drawlabels = FALSE, add = TRUE,
          levels = c(0.01, 0.03, 0.06, 0.1, 0.14))

  # Draw all previous staircase segments
  if (t > 2) {
    for (s in 2:(t - 1)) {
      segments(samp_step[s - 1, 1], samp_step[s - 1, 2],
               samp_step[s, 1], samp_step[s - 1, 2],
               col = rgb(0.5, 0.5, 0.5, 0.3), lwd = 0.8)
      segments(samp_step[s, 1], samp_step[s - 1, 2],
               samp_step[s, 1], samp_step[s, 2],
               col = rgb(0.5, 0.5, 0.5, 0.3), lwd = 0.8)
    }
    points(samp_step[2:(t - 1), 1], samp_step[2:(t - 1), 2],
           pch = 20, cex = 0.7, col = "gray50")
  }

  # Current step: horizontal then vertical
  # Horizontal: theta1 update (green)
  segments(samp_step[t - 1, 1], samp_step[t - 1, 2],
           samp_step[t, 1], samp_step[t - 1, 2],
           col = "green4", lwd = 2, lty = 1)
  # Intermediate point
  points(samp_step[t, 1], samp_step[t - 1, 2],
         pch = 4, cex = 1.2, col = "green4", lwd = 2)

  # Vertical: theta2 update (blue)
  segments(samp_step[t, 1], samp_step[t - 1, 2],
           samp_step[t, 1], samp_step[t, 2],
           col = "steelblue", lwd = 2, lty = 1)

  # New sample
  points(samp_step[t, 1], samp_step[t, 2],
         pch = 8, cex = 1.5, col = "darkorange", lwd = 2)

  # Start point
  points(samp_step[1, 1], samp_step[1, 2],
         pch = 17, cex = 1.2, col = "firebrick")

  # Conditional distribution annotation (on first panel)
  if (k == 1) {
    legend("bottomright", cex = 0.5,
           legend = c(TeX("Draw $\\theta_1 \\mid \\theta_2$"),
                      TeX("Draw $\\theta_2 \\mid \\theta_1$"),
                      "New sample", "Start"),
           col = c("green4", "steelblue", "darkorange", "firebrick"),
           lwd = c(2, 2, NA, NA), lty = c(1, 1, NA, NA),
           pch = c(NA, NA, 8, 17), bg = "white")
  }
}

dev.off()


# ============================================================
# Figure 3: Effect of correlation on mixing (2x2)
# ============================================================

cat("Generating Figure 3: Correlation effect on mixing...\n")

set.seed(42)
n_corr <- 2000
rho_vals <- c(0.2, 0.7, 0.95)
samp_low <- gibbs_bivariate(n_corr, rho = rho_vals[1], init = c(-2, -2))
samp_med <- gibbs_bivariate(n_corr, rho = rho_vals[2], init = c(-2, -2))
samp_high <- gibbs_bivariate(n_corr, rho = rho_vals[3], init = c(-2, -2))

pdf("plots/fig_gibbs_correlation.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Low correlation staircase
g_corr <- seq(-4, 4, length.out = 100)
z_low <- outer(g_corr, g_corr, function(a, b) bvn_density(a, b, rho_vals[1]))
plot(NULL, xlim = c(-4, 4), ylim = c(-4, 4),
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"),
     main = TeX(sprintf("(a) $\\rho = %.1f$: Fast Mixing", rho_vals[1])))
contour(g_corr, g_corr, z_low, col = "gray70", lwd = 0.8,
        drawlabels = FALSE, add = TRUE)
n_draw <- 50
for (t in 2:n_draw) {
  segments(samp_low[t - 1, 1], samp_low[t - 1, 2],
           samp_low[t, 1], samp_low[t - 1, 2],
           col = rgb(0.13, 0.55, 0.13, 0.4), lwd = 0.8)
  segments(samp_low[t, 1], samp_low[t - 1, 2],
           samp_low[t, 1], samp_low[t, 2],
           col = rgb(0.27, 0.51, 0.71, 0.4), lwd = 0.8)
}
points(samp_low[1:n_draw, 1], samp_low[1:n_draw, 2],
       pch = 20, cex = 0.4, col = "gray30")

# Panel B: Moderate correlation
z_med <- outer(g_corr, g_corr, function(a, b) bvn_density(a, b, rho_vals[2]))
plot(NULL, xlim = c(-4, 4), ylim = c(-4, 4),
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"),
     main = TeX(sprintf("(b) $\\rho = %.1f$: Moderate Mixing", rho_vals[2])))
contour(g_corr, g_corr, z_med, col = "gray70", lwd = 0.8,
        drawlabels = FALSE, add = TRUE)
for (t in 2:n_draw) {
  segments(samp_med[t - 1, 1], samp_med[t - 1, 2],
           samp_med[t, 1], samp_med[t - 1, 2],
           col = rgb(0.13, 0.55, 0.13, 0.4), lwd = 0.8)
  segments(samp_med[t, 1], samp_med[t - 1, 2],
           samp_med[t, 1], samp_med[t, 2],
           col = rgb(0.27, 0.51, 0.71, 0.4), lwd = 0.8)
}
points(samp_med[1:n_draw, 1], samp_med[1:n_draw, 2],
       pch = 20, cex = 0.4, col = "gray30")

# Panel C: High correlation
z_high <- outer(g_corr, g_corr, function(a, b) bvn_density(a, b, rho_vals[3]))
plot(NULL, xlim = c(-4, 4), ylim = c(-4, 4),
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"),
     main = TeX(sprintf("(c) $\\rho = %.2f$: Slow Mixing", rho_vals[3])))
contour(g_corr, g_corr, z_high, col = "gray70", lwd = 0.8,
        drawlabels = FALSE, add = TRUE)
for (t in 2:n_draw) {
  segments(samp_high[t - 1, 1], samp_high[t - 1, 2],
           samp_high[t, 1], samp_high[t - 1, 2],
           col = rgb(0.13, 0.55, 0.13, 0.4), lwd = 0.8)
  segments(samp_high[t, 1], samp_high[t - 1, 2],
           samp_high[t, 1], samp_high[t, 2],
           col = rgb(0.27, 0.51, 0.71, 0.4), lwd = 0.8)
}
points(samp_high[1:n_draw, 1], samp_high[1:n_draw, 2],
       pch = 20, cex = 0.4, col = "gray30")

# Panel D: Autocorrelation comparison
max_lag <- 50
acf_low <- acf(samp_low[, 1], lag.max = max_lag, plot = FALSE)
acf_med <- acf(samp_med[, 1], lag.max = max_lag, plot = FALSE)
acf_high <- acf(samp_high[, 1], lag.max = max_lag, plot = FALSE)

plot(0:max_lag, acf_low$acf, type = "l", col = "steelblue", lwd = 2,
     xlab = "Lag", ylab = TeX("ACF of $\\theta_1$"),
     main = TeX("(d) Autocorrelation vs $\\rho$"),
     ylim = c(-0.1, 1))
lines(0:max_lag, acf_med$acf, col = "darkorange", lwd = 2, lty = 2)
lines(0:max_lag, acf_high$acf, col = "firebrick", lwd = 2, lty = 4)
abline(h = 0, col = "gray50", lty = 3)
legend("topright", cex = 0.65,
       legend = c(TeX(sprintf("$\\rho = %.1f$", rho_vals[1])),
                  TeX(sprintf("$\\rho = %.1f$", rho_vals[2])),
                  TeX(sprintf("$\\rho = %.2f$", rho_vals[3]))),
       col = c("steelblue", "darkorange", "firebrick"),
       lwd = 2, lty = c(1, 2, 4), bg = "white")

dev.off()


# ============================================================
# Figure 4: Bayesian Normal Model — Gibbs Sampler (2x2)
# Target: y_i ~ N(mu, sigma^2), i=1,...,n
# Prior: mu ~ N(mu0, tau0^2), sigma^2 ~ Inv-Chi^2(nu0, sigma0^2)
# ============================================================

cat("Generating Figure 4: Bayesian normal model...\n")

set.seed(42)
# Generate data
n_data <- 30
true_mu <- 5
true_sigma2 <- 4
y_data <- rnorm(n_data, mean = true_mu, sd = sqrt(true_sigma2))
ybar <- mean(y_data)
s2 <- var(y_data)

# Priors
mu0 <- 0; tau0_sq <- 100  # vague prior on mu
nu0 <- 1; sigma0_sq <- 1  # weakly informative on sigma^2

# Gibbs sampler
n_gibbs <- 10000
mu_samp <- numeric(n_gibbs)
sigma2_samp <- numeric(n_gibbs)
mu_samp[1] <- 0
sigma2_samp[1] <- 1

for (t in 2:n_gibbs) {
  # Full conditional: mu | sigma^2, y
  #   precision = n/sigma^2 + 1/tau0^2
  #   mean = (n*ybar/sigma^2 + mu0/tau0^2) / precision
  prec <- n_data / sigma2_samp[t - 1] + 1 / tau0_sq
  cond_mu_mean <- (n_data * ybar / sigma2_samp[t - 1] + mu0 / tau0_sq) / prec
  mu_samp[t] <- rnorm(1, cond_mu_mean, sqrt(1 / prec))

  # Full conditional: sigma^2 | mu, y ~ Inv-Chi^2(nu_n, sigma_n^2)
  #   nu_n = nu0 + n
  #   nu_n * sigma_n^2 = nu0*sigma0^2 + sum((y_i - mu)^2)
  nu_n <- nu0 + n_data
  ss <- sum((y_data - mu_samp[t])^2)
  scale_n <- (nu0 * sigma0_sq + ss) / nu_n
  # draw from Inv-Chi^2(nu_n, scale_n) = Inv-Gamma(nu_n/2, nu_n*scale_n/2)
  sigma2_samp[t] <- 1 / rgamma(1, shape = nu_n / 2, rate = nu_n * scale_n / 2)
}

cat(sprintf("Normal model: n=%d, ybar=%.2f, s2=%.2f\n", n_data, ybar, s2))
cat(sprintf("  Posterior mu:     mean=%.3f, sd=%.3f\n",
            mean(mu_samp[burn:n_gibbs]), sd(mu_samp[burn:n_gibbs])))
cat(sprintf("  Posterior sigma2: mean=%.3f, sd=%.3f\n",
            mean(sigma2_samp[burn:n_gibbs]), sd(sigma2_samp[burn:n_gibbs])))

pdf("plots/fig_gibbs_normal_model.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Trace plot of mu
plot(1:n_gibbs, mu_samp, type = "l", col = rgb(0.13, 0.55, 0.13, 0.3),
     lwd = 0.5, xlab = "Iteration", ylab = TeX("$\\mu$"),
     main = TeX("(a) Trace Plot: $\\mu$"))
abline(h = true_mu, col = "firebrick", lwd = 2, lty = 2)
abline(v = burn, col = "gray40", lty = 3)
text(burn + 200, max(mu_samp) * 0.95, "burn-in", col = "gray40", adj = 0, cex = 0.7)

# Panel B: Trace plot of sigma^2
plot(1:n_gibbs, sigma2_samp, type = "l", col = rgb(0.27, 0.51, 0.71, 0.3),
     lwd = 0.5, xlab = "Iteration", ylab = TeX("$\\sigma^2$"),
     main = TeX("(b) Trace Plot: $\\sigma^2$"))
abline(h = true_sigma2, col = "firebrick", lwd = 2, lty = 2)
abline(v = burn, col = "gray40", lty = 3)

# Panel C: Joint posterior with Gibbs path
plot(mu_samp[burn:n_gibbs], sigma2_samp[burn:n_gibbs],
     pch = 20, cex = 0.15, col = rgb(0.2, 0.4, 0.8, 0.15),
     xlab = TeX("$\\mu$"), ylab = TeX("$\\sigma^2$"),
     main = "(c) Joint Posterior Samples")
# Overlay staircase for a few iterations
n_path <- 40
start_path <- burn + 1
for (t in start_path:(start_path + n_path)) {
  segments(mu_samp[t - 1], sigma2_samp[t - 1],
           mu_samp[t], sigma2_samp[t - 1],
           col = rgb(0.13, 0.55, 0.13, 0.5), lwd = 0.8)
  segments(mu_samp[t], sigma2_samp[t - 1],
           mu_samp[t], sigma2_samp[t],
           col = rgb(0.27, 0.51, 0.71, 0.5), lwd = 0.8)
}
points(true_mu, true_sigma2, pch = 4, cex = 2, col = "firebrick", lwd = 2)

# Panel D: Marginal posterior histograms
xg_mu <- seq(min(mu_samp[burn:n_gibbs]), max(mu_samp[burn:n_gibbs]),
             length.out = 200)
hist(mu_samp[burn:n_gibbs], breaks = 50, freq = FALSE,
     col = rgb(0.13, 0.55, 0.13, 0.3), border = "white",
     main = TeX("(d) Marginal Posterior of $\\mu$"),
     xlab = TeX("$\\mu$"), ylab = "Density")
# Analytic posterior of mu (approximately, with large tau0)
# p(mu|y) ~ N(ybar, sigma2/n) approximately
approx_sd <- sqrt(mean(sigma2_samp[burn:n_gibbs]) / n_data)
lines(xg_mu, dnorm(xg_mu, ybar, approx_sd), col = "firebrick", lwd = 2, lty = 2)
abline(v = true_mu, col = "darkorange", lwd = 1.5, lty = 4)
legend("topright", cex = 0.6,
       legend = c("Gibbs histogram", "Normal approx", TeX("True $\\mu$")),
       fill = c(rgb(0.13, 0.55, 0.13, 0.3), NA, NA),
       border = c("white", NA, NA),
       col = c(NA, "firebrick", "darkorange"),
       lwd = c(NA, 2, 1.5), lty = c(NA, 2, 4), bg = "white")

dev.off()


# ============================================================
# Figure 5: Hierarchical Normal Model — Gibbs  (2x2)
# Eight-schools style: y_j | theta_j ~ N(theta_j, sigma_j^2)
#                      theta_j | mu, tau ~ N(mu, tau^2)
#                      mu ~ N(0, 100), tau^2 ~ Inv-Chi^2(1, 1)
# ============================================================

cat("Generating Figure 5: Hierarchical normal model...\n")

set.seed(42)
# Eight schools data (from BDA3)
J <- 8
y_schools <- c(28, 8, -3, 7, -1, 1, 18, 12)
se_schools <- c(15, 10, 16, 11, 9, 11, 10, 18)
sigma2_schools <- se_schools^2

# Gibbs sampler for hierarchical model
n_hier <- 15000
burn_hier <- 2000
theta_hier <- matrix(NA, n_hier, J)
mu_hier <- numeric(n_hier)
tau2_hier <- numeric(n_hier)

# Initialize
theta_hier[1, ] <- y_schools
mu_hier[1] <- mean(y_schools)
tau2_hier[1] <- var(y_schools)

# Prior hyperparameters
mu_prior_mean <- 0; mu_prior_var <- 10000
tau2_prior_nu <- 1; tau2_prior_s2 <- 1

for (t in 2:n_hier) {
  tau2_curr <- tau2_hier[t - 1]

  # 1. Update each theta_j | mu, tau^2, y_j
  for (j in 1:J) {
    prec_j <- 1 / sigma2_schools[j] + 1 / tau2_curr
    mean_j <- (y_schools[j] / sigma2_schools[j] + mu_hier[t - 1] / tau2_curr) / prec_j
    theta_hier[t, j] <- rnorm(1, mean_j, sqrt(1 / prec_j))
  }

  # 2. Update mu | theta, tau^2
  prec_mu <- J / tau2_curr + 1 / mu_prior_var
  mean_mu <- (sum(theta_hier[t, ]) / tau2_curr + mu_prior_mean / mu_prior_var) / prec_mu
  mu_hier[t] <- rnorm(1, mean_mu, sqrt(1 / prec_mu))

  # 3. Update tau^2 | theta, mu  ~ Inv-Chi^2(nu_n, s_n^2)
  nu_n_tau <- tau2_prior_nu + J
  ss_tau <- tau2_prior_nu * tau2_prior_s2 + sum((theta_hier[t, ] - mu_hier[t])^2)
  scale_tau <- ss_tau / nu_n_tau
  tau2_hier[t] <- 1 / rgamma(1, shape = nu_n_tau / 2, rate = nu_n_tau * scale_tau / 2)
}

tau_hier <- sqrt(tau2_hier)

cat(sprintf("Hierarchical model: posterior mu=%.2f (sd=%.2f), posterior tau=%.2f (sd=%.2f)\n",
            mean(mu_hier[burn_hier:n_hier]), sd(mu_hier[burn_hier:n_hier]),
            mean(tau_hier[burn_hier:n_hier]), sd(tau_hier[burn_hier:n_hier])))

pdf("plots/fig_gibbs_hierarchical.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Trace plots for selected theta_j's
theta_cols <- c("steelblue", "darkorange", "green4", "firebrick")
plot(NULL, xlim = c(1, n_hier), ylim = range(theta_hier[burn_hier:n_hier, ]),
     xlab = "Iteration", ylab = TeX("$\\theta_j$"),
     main = TeX("(a) Trace Plots: Selected $\\theta_j$"))
for (jj in c(1, 3, 5, 7)) {
  idx <- which(c(1, 3, 5, 7) == jj)
  lines(1:n_hier, theta_hier[, jj], col = adjustcolor(theta_cols[idx], 0.3), lwd = 0.4)
}
legend("topright", cex = 0.55,
       legend = paste0("School ", c(1, 3, 5, 7)),
       col = theta_cols, lwd = 1.5, bg = "white")

# Panel B: Posterior of mu (population mean)
xg_mu_h <- seq(quantile(mu_hier[burn_hier:n_hier], 0.001),
               quantile(mu_hier[burn_hier:n_hier], 0.999), length.out = 200)
hist(mu_hier[burn_hier:n_hier], breaks = 60, freq = FALSE,
     col = rgb(0.27, 0.51, 0.71, 0.4), border = "white",
     main = TeX("(b) Posterior of $\\mu$ (Population Mean)"),
     xlab = TeX("$\\mu$"), ylab = "Density")
abline(v = mean(y_schools), col = "darkorange", lwd = 2, lty = 2)
abline(v = mean(mu_hier[burn_hier:n_hier]), col = "firebrick", lwd = 2, lty = 4)
legend("topright", cex = 0.6,
       legend = c(TeX("Posterior of $\\mu$"),
                  TeX("$\\bar{y}$ (raw mean)"),
                  "Posterior mean"),
       fill = c(rgb(0.27, 0.51, 0.71, 0.4), NA, NA),
       border = c("white", NA, NA),
       col = c(NA, "darkorange", "firebrick"),
       lwd = c(NA, 2, 2), lty = c(NA, 2, 4), bg = "white")

# Panel C: Posterior of tau
hist(tau_hier[burn_hier:n_hier], breaks = 60, freq = FALSE,
     col = rgb(0.93, 0.49, 0.19, 0.4), border = "white",
     main = TeX("(c) Posterior of $\\tau$ (Between-School SD)"),
     xlab = TeX("$\\tau$"), ylab = "Density",
     xlim = c(0, quantile(tau_hier[burn_hier:n_hier], 0.99)))
abline(v = median(tau_hier[burn_hier:n_hier]), col = "firebrick", lwd = 2, lty = 2)

# Panel D: Shrinkage plot
theta_means <- colMeans(theta_hier[burn_hier:n_hier, ])
theta_lo <- apply(theta_hier[burn_hier:n_hier, ], 2, quantile, 0.025)
theta_hi <- apply(theta_hier[burn_hier:n_hier, ], 2, quantile, 0.975)

plot(y_schools, 1:J, pch = 4, cex = 1.2, col = "firebrick", lwd = 2,
     xlim = range(c(theta_lo, theta_hi, y_schools)) + c(-2, 2),
     ylim = c(0.5, J + 0.5), yaxt = "n",
     xlab = "Effect estimate", ylab = "",
     main = "(d) Shrinkage: Raw vs Posterior")
axis(2, at = 1:J, labels = paste("School", 1:J), las = 1, cex.axis = 0.7)
points(theta_means, 1:J, pch = 19, cex = 1, col = "steelblue")
segments(theta_lo, 1:J, theta_hi, 1:J, col = "steelblue", lwd = 1.5)
# Connect raw to posterior
segments(y_schools, 1:J, theta_means, 1:J, col = "gray60", lty = 3, lwd = 0.8)
# Grand mean
abline(v = mean(mu_hier[burn_hier:n_hier]), col = "darkorange", lty = 2, lwd = 1.5)
legend("bottomright", cex = 0.55,
       legend = c(TeX("Raw estimate $y_j$"),
                  TeX("Posterior mean $E[\\theta_j \\mid y]$"),
                  "95% CI",
                  TeX("Grand mean $\\hat{\\mu}$")),
       pch = c(4, 19, NA, NA),
       col = c("firebrick", "steelblue", "steelblue", "darkorange"),
       lwd = c(2, NA, 1.5, 1.5), lty = c(NA, NA, 1, 2), bg = "white")

dev.off()


# ============================================================
# Figure 6: Convergence Diagnostics (2x2)
# Run multiple chains from dispersed starting points
# ============================================================

cat("Generating Figure 6: Convergence diagnostics...\n")

set.seed(42)
n_diag <- 3000
n_chains <- 4
starts <- list(c(-4, 4), c(4, -4), c(-4, -4), c(4, 4))
chain_cols <- c("steelblue", "firebrick", "green4", "darkorange")

chains <- list()
for (ch in 1:n_chains) {
  chains[[ch]] <- gibbs_bivariate(n_diag, rho = rho1, init = starts[[ch]])
}

pdf("plots/fig_gibbs_diagnostics.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Multiple chains overlaid — trace of theta1
plot(NULL, xlim = c(1, n_diag), ylim = c(-4.5, 4.5),
     xlab = "Iteration", ylab = TeX("$\\theta_1$"),
     main = TeX("(a) Multiple Chains: $\\theta_1$"))
for (ch in 1:n_chains) {
  lines(1:n_diag, chains[[ch]][, 1],
        col = adjustcolor(chain_cols[ch], 0.4), lwd = 0.5)
}
abline(h = 0, col = "gray40", lty = 3)
legend("topright", cex = 0.55,
       legend = paste("Chain", 1:n_chains),
       col = chain_cols, lwd = 1.5, bg = "white")

# Panel B: Running mean plot
plot(NULL, xlim = c(1, n_diag), ylim = c(-3, 3),
     xlab = "Iteration", ylab = TeX("Running mean of $\\theta_1$"),
     main = TeX("(b) Running Mean: $\\theta_1$"))
for (ch in 1:n_chains) {
  rmean <- cumsum(chains[[ch]][, 1]) / (1:n_diag)
  lines(1:n_diag, rmean, col = chain_cols[ch], lwd = 1.2)
}
abline(h = 0, col = "gray40", lty = 2, lwd = 1.5)

# Panel C: R-hat over iterations
# Compute R-hat at each iteration t using all 4 chains
rhat_seq <- numeric(n_diag)
eval_pts <- seq(50, n_diag, by = 10)
rhat_eval <- numeric(length(eval_pts))

for (idx in seq_along(eval_pts)) {
  t_cur <- eval_pts[idx]
  chain_means <- sapply(chains, function(ch) mean(ch[1:t_cur, 1]))
  chain_vars <- sapply(chains, function(ch) var(ch[1:t_cur, 1]))
  W <- mean(chain_vars)
  B <- t_cur * var(chain_means)
  var_plus <- ((t_cur - 1) / t_cur) * W + (1 / t_cur) * B
  rhat_eval[idx] <- sqrt(var_plus / W)
}

plot(eval_pts, rhat_eval, type = "l", col = "steelblue", lwd = 2,
     xlab = "Iteration", ylab = TeX("$\\hat{R}$"),
     main = TeX("(c) $\\hat{R}$ Convergence Diagnostic"),
     ylim = c(0.95, max(rhat_eval) * 1.02))
abline(h = 1, col = "firebrick", lty = 2, lwd = 1.5)
abline(h = 1.1, col = "darkorange", lty = 3, lwd = 1)
text(n_diag * 0.7, 1.12, TeX("$\\hat{R} = 1.1$ threshold"),
     col = "darkorange", cex = 0.7)

# Panel D: Autocorrelation plot (post-burn-in, chain 1)
burn_diag <- 500
acf_chain <- acf(chains[[1]][burn_diag:n_diag, 1], lag.max = 40, plot = FALSE)
plot(0:40, acf_chain$acf, type = "h", col = "steelblue", lwd = 3,
     xlab = "Lag", ylab = "ACF",
     main = TeX("(d) Autocorrelation ($\\rho = 0.7$, Chain 1)"))
abline(h = 0, col = "gray40")
# add approximate 95% CI band
n_eff_approx <- n_diag - burn_diag
abline(h = c(-1, 1) * 1.96 / sqrt(n_eff_approx),
       col = "firebrick", lty = 2, lwd = 1)

dev.off()


# ============================================================
# Figure 7: Where Gibbs Fails — Multimodal + Strong Correlation (2x2)
# ============================================================

cat("Generating Figure 7: Failure modes...\n")

set.seed(42)

# --- Failure 1: Multimodal target ---
# Bimodal: mixture of two bivariate normals
# Gibbs can get trapped in one mode
n_fail <- 5000
theta_bimodal <- matrix(NA, n_fail, 2)
theta_bimodal[1, ] <- c(3, 3)

# Target: 0.5*N([3,3], 0.5*I) + 0.5*N([-3,-3], 0.5*I)
# Component-wise full conditionals are mixtures — we approximate by
# treating it like a single-component sampler to show the failure
for (t in 2:n_fail) {
  x2_prev <- theta_bimodal[t - 1, 2]
  # Conditional p(x1 | x2) is a mixture of two normals
  # Component means for x1: 3 and -3
  # Component weights depend on x2
  w1_cond <- 0.5 * dnorm(x2_prev, 3, sqrt(0.5))
  w2_cond <- 0.5 * dnorm(x2_prev, -3, sqrt(0.5))
  p1 <- w1_cond / (w1_cond + w2_cond)
  if (runif(1) < p1) {
    theta_bimodal[t, 1] <- rnorm(1, 3, sqrt(0.5))
  } else {
    theta_bimodal[t, 1] <- rnorm(1, -3, sqrt(0.5))
  }

  x1_cur <- theta_bimodal[t, 1]
  w1_cond2 <- 0.5 * dnorm(x1_cur, 3, sqrt(0.5))
  w2_cond2 <- 0.5 * dnorm(x1_cur, -3, sqrt(0.5))
  p1_2 <- w1_cond2 / (w1_cond2 + w2_cond2)
  if (runif(1) < p1_2) {
    theta_bimodal[t, 2] <- rnorm(1, 3, sqrt(0.5))
  } else {
    theta_bimodal[t, 2] <- rnorm(1, -3, sqrt(0.5))
  }
}

# --- Failure 2: Well separated modes with low bridge probability ---
# Use widely separated modes where conditional overlap is negligible
n_fail2 <- 5000
theta_stuck <- matrix(NA, n_fail2, 2)
theta_stuck[1, ] <- c(5, 5)

for (t in 2:n_fail2) {
  x2_prev <- theta_stuck[t - 1, 2]
  # modes at (5,5) and (-5,-5), sigma=0.5 — conditional overlap ~ 0
  w1 <- 0.5 * dnorm(x2_prev, 5, 0.5)
  w2 <- 0.5 * dnorm(x2_prev, -5, 0.5)
  p1 <- w1 / (w1 + w2)
  if (runif(1) < p1) {
    theta_stuck[t, 1] <- rnorm(1, 5, 0.5)
  } else {
    theta_stuck[t, 1] <- rnorm(1, -5, 0.5)
  }

  x1_cur <- theta_stuck[t, 1]
  w1_2 <- 0.5 * dnorm(x1_cur, 5, 0.5)
  w2_2 <- 0.5 * dnorm(x1_cur, -5, 0.5)
  p1_2 <- w1_2 / (w1_2 + w2_2)
  if (runif(1) < p1_2) {
    theta_stuck[t, 2] <- rnorm(1, 5, 0.5)
  } else {
    theta_stuck[t, 2] <- rnorm(1, -5, 0.5)
  }
}


pdf("plots/fig_gibbs_failure.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), family = "serif")

# Panel A: Bimodal — close modes, Gibbs can mix (but slowly)
g_bm <- seq(-6, 6, length.out = 100)
z_bm <- outer(g_bm, g_bm, function(a, b) {
  0.5 * exp(-((a - 3)^2 + (b - 3)^2) / 1) / (2 * pi * 0.5) +
  0.5 * exp(-((a + 3)^2 + (b + 3)^2) / 1) / (2 * pi * 0.5)
})
plot(NULL, xlim = c(-6, 6), ylim = c(-6, 6),
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"),
     main = "(a) Bimodal: Close Modes (Mixing)")
contour(g_bm, g_bm, z_bm, col = "gray65", lwd = 0.8,
        drawlabels = FALSE, add = TRUE)
points(theta_bimodal[, 1], theta_bimodal[, 2],
       pch = 20, cex = 0.15, col = rgb(0.2, 0.4, 0.8, 0.15))
points(c(3, -3), c(3, -3), pch = 4, cex = 2, col = "firebrick", lwd = 2)

# Panel B: Well-separated modes — Gibbs gets stuck
g_sep <- seq(-8, 8, length.out = 100)
z_sep <- outer(g_sep, g_sep, function(a, b) {
  0.5 * exp(-((a - 5)^2 + (b - 5)^2) / 0.5) / (2 * pi * 0.25) +
  0.5 * exp(-((a + 5)^2 + (b + 5)^2) / 0.5) / (2 * pi * 0.25)
})
plot(NULL, xlim = c(-8, 8), ylim = c(-8, 8),
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"),
     main = "(b) Separated Modes: Gibbs Stuck")
contour(g_sep, g_sep, z_sep, col = "gray65", lwd = 0.8,
        drawlabels = FALSE, add = TRUE)
points(theta_stuck[, 1], theta_stuck[, 2],
       pch = 20, cex = 0.15, col = rgb(0.8, 0.2, 0.2, 0.2))
points(c(5, -5), c(5, -5), pch = 4, cex = 2, col = "firebrick", lwd = 2)
text(0, 0, "Unreachable\nmode", col = "firebrick", cex = 0.8, font = 2)

# Panel C: Trace plot showing stuck behavior
plot(1:n_fail2, theta_stuck[, 1], type = "l",
     col = rgb(0.8, 0.2, 0.2, 0.5), lwd = 0.5,
     xlab = "Iteration", ylab = TeX("$\\theta_1$"),
     main = TeX("(c) Trace: $\\theta_1$ (Stuck in One Mode)"))
abline(h = c(5, -5), col = "gray50", lty = 2, lwd = 1)
text(n_fail2 * 0.8, -5, "Missed mode", cex = 0.7, col = "gray40")

# Panel D: Strong correlation — compare ESS
# Compute effective number of independent samples via ACF
compute_ess <- function(x, max_lag = 500) {
  n <- length(x)
  acf_vals <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1]
  # Use initial monotone sequence estimator (truncate at first negative pair)
  sum_rho <- 0
  for (k in seq(1, length(acf_vals) - 1, by = 2)) {
    pair_sum <- acf_vals[k] + acf_vals[k + 1]
    if (pair_sum < 0) break
    sum_rho <- sum_rho + pair_sum
  }
  n / (1 + 2 * sum_rho)
}

rho_sweep <- seq(0, 0.98, by = 0.02)
ess_sweep <- numeric(length(rho_sweep))
n_ess_test <- 5000

set.seed(42)
for (ir in seq_along(rho_sweep)) {
  samp_tmp <- gibbs_bivariate(n_ess_test, rho = rho_sweep[ir])
  ess_sweep[ir] <- compute_ess(samp_tmp[500:n_ess_test, 1])
}

plot(rho_sweep, ess_sweep / (n_ess_test - 500) * 100, type = "b",
     pch = 19, cex = 0.5, col = "steelblue", lwd = 1.5,
     xlab = TeX("Correlation $\\rho$"),
     ylab = "ESS / N (%)",
     main = "(d) Effective Sample Size vs Correlation")
abline(v = 0.95, col = "firebrick", lty = 2, lwd = 1)
text(0.88, max(ess_sweep / (n_ess_test - 500) * 100) * 0.3,
     TeX("$\\rho = 0.95$"), col = "firebrick", cex = 0.7)

dev.off()


cat("\nAll figures saved to plots/:\n")
cat("  fig_gibbs_bivariate.pdf\n")
cat("  fig_gibbs_steps.pdf\n")
cat("  fig_gibbs_correlation.pdf\n")
cat("  fig_gibbs_normal_model.pdf\n")
cat("  fig_gibbs_hierarchical.pdf\n")
cat("  fig_gibbs_diagnostics.pdf\n")
cat("  fig_gibbs_failure.pdf\n")
