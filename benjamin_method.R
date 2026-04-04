# ============================================================
# The Benjamin Method (Math 310 Multiple Comparisons Project)
# Weighted Step-Down Pairwise Comparisons with fixed weights
# Author: Wyatt Benjamin
# ============================================================

# ------------------------------------------------------------
# Helper: build all pair indices (i < j)
# ------------------------------------------------------------
.all_pairs <- function(I) {
  pairs <- t(combn(1:I, 2))
  colnames(pairs) <- c("group1", "group2")
  pairs
}

# ------------------------------------------------------------
# Helper: fixed weights based only on index distance |i - j|
# gamma controls how concentrated the weights are:
#   gamma = 0   => uniform weights
#   gamma > 0   => more weight on far-apart comparisons
# ------------------------------------------------------------
.benjamin_weights <- function(I, gamma = 0.8) {
  pairs <- .all_pairs(I)
  d <- abs(pairs[, "group1"] - pairs[, "group2"])   # 1, 2, ..., I-1
  w_raw <- d^gamma
  w <- w_raw / sum(w_raw)
  list(pairs = pairs, d = d, w = w)
}

# ------------------------------------------------------------
# Core step-down rule (weighted)
# Input:
#   p: vector of p-values length K
#   w: vector of weights length K, sum(w) = 1 (fixed)
# Output:
#   reject: TRUE/FALSE decisions length K
#   step_order: the order induced by sorting p-values (1 = smallest p)
#   alpha_step: per-hypothesis step threshold used in its position
# ------------------------------------------------------------
.benjamin_stepdown <- function(p, w, alpha = 0.05) {
  K <- length(p)
  ord <- order(p)               # smallest p first
  p_s <- p[ord]
  w_s <- w[ord]

  # cumulative weight spent before each step (sorted space)
  cum_prev <- c(0, cumsum(w_s)[1:(K-1)])

  # step thresholds in sorted space
  alpha_s <- alpha * w_s / (1 - cum_prev)

  # decide: reject sequentially until first failure
  reject_s <- rep(FALSE, K)
  for (m in 1:K) {
    if (p_s[m] <= alpha_s[m]) {
      reject_s[m] <- TRUE
    } else {
      break
    }
  }

  # map back to original hypothesis order
  reject <- rep(FALSE, K); reject[ord] <- reject_s
  step_order <- rep(NA_integer_, K); step_order[ord] <- 1:K
  alpha_step <- rep(NA_real_, K); alpha_step[ord] <- alpha_s

  list(reject = reject, step_order = step_order, alpha_step = alpha_step, ord = ord)
}

# ------------------------------------------------------------
# REQUIRED FUNCTION: my_comparison(means, J, MSE, alpha)
# Implements the Benjamin method and returns a data frame.
# ------------------------------------------------------------
my_comparison <- function(means, J, MSE, alpha = 0.05, gamma = 0.8) {
  means <- as.numeric(means)
  I <- length(means)
  if (I < 2) stop("Need at least 2 groups.")
  if (J <= 1) stop("Need J > 1.")
  if (MSE <= 0) stop("Need MSE > 0.")

  df_error <- I * (J - 1)
  se <- sqrt(2 * MSE / J)

  W <- .benjamin_weights(I, gamma = gamma)
  pairs <- W$pairs
  w <- W$w
  d <- W$d

  # compute t and p for each pair
  diff <- means[pairs[, "group1"]] - means[pairs[, "group2"]]
  tstat <- abs(diff) / se
  pval <- 2 * pt(-abs(tstat), df = df_error)

  # step-down decisions
  SD <- .benjamin_stepdown(pval, w, alpha = alpha)

  out <- data.frame(
    group1 = pairs[, "group1"],
    group2 = pairs[, "group2"],
    index_dist = d,
    diff = diff,
    t = tstat,
    p = pval,
    w = w,
    step_order = SD$step_order,
    alpha_step = SD$alpha_step,
    significant = SD$reject
  )

  # nice display: sorted by step_order (i.e., by p-value)
  out <- out[order(out$step_order), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------
# Convenience: compute MSE from a data matrix (I x J)
# (equal sample sizes)
# ------------------------------------------------------------
mse_from_matrix <- function(X) {
  # X is I x J, each row a group
  vars <- apply(X, 1, var)
  mean(vars)  # equals pooled MSE when J is constant
}

# ------------------------------------------------------------
# Simulation: estimate FWER under the global null (all means equal)
# Returns: estimated FWER = P(at least one rejection)
# ------------------------------------------------------------
simulate_fwer <- function(I, J, n_sim = 10000, sigma = 1, alpha = 0.05, gamma = 0.8, seed = 1) {
  set.seed(seed)
  any_reject <- 0

  for (s in 1:n_sim) {
    X <- matrix(rnorm(I * J, mean = 0, sd = sigma), nrow = I, ncol = J)
    means <- rowMeans(X)
    MSE <- mse_from_matrix(X)

    res <- my_comparison(means, J, MSE, alpha = alpha, gamma = gamma)
    if (any(res$significant)) any_reject <- any_reject + 1
  }

  any_reject / n_sim
}

# ------------------------------------------------------------
# Power analysis:
# I = 5, J = 10, sigma = 1
# true means: (0, 0, 0, delta, delta)
# reports power for detecting mu4 != mu1
# Tukey power computed via Tukey HSD critical difference
# ------------------------------------------------------------
simulate_power <- function(deltas = c(0.5, 1.0, 1.5, 2.0),
                           n_sim = 1000, alpha = 0.05, gamma = 0.8, seed = 1) {
  set.seed(seed)
  I <- 5; J <- 10; sigma <- 1
  df_error <- I * (J - 1)

  # Tukey critical factor for equal-n HSD
  qcrit <- qtukey(1 - alpha, nmeans = I, df = df_error)

  out <- data.frame(delta = deltas, power_benjamin = NA_real_, power_tukey = NA_real_)

  for (k in seq_along(deltas)) {
    delta <- deltas[k]
    mu <- c(0, 0, 0, delta, delta)

    ben_count <- 0
    tuk_count <- 0

    for (s in 1:n_sim) {
      X <- matrix(rnorm(I * J, mean = rep(mu, each = J), sd = sigma), nrow = I, byrow = TRUE)
      means <- rowMeans(X)
      MSE <- mse_from_matrix(X)

      # Benjamin decision for pair (1,4)
      res <- my_comparison(means, J, MSE, alpha = alpha, gamma = gamma)
      ben_sig <- any(res$significant & res$group1 == 1 & res$group2 == 4)

      # Tukey decision for pair (1,4): |diff| > qcrit * sqrt(MSE / J)
      hsd <- qcrit * sqrt(MSE / J)
      tuk_sig <- abs(means[1] - means[4]) > hsd

      if (ben_sig) ben_count <- ben_count + 1
      if (tuk_sig) tuk_count <- tuk_count + 1
    }

    out$power_benjamin[k] <- ben_count / n_sim
    out$power_tukey[k] <- tuk_count / n_sim
  }

  out
}

# ------------------------------------------------------------
# DEMO: Example 11.5 (oil filters) from the handout
# ------------------------------------------------------------
demo_oil_filters <- function() {
  means <- c(14.5, 13.8, 13.3, 14.3, 13.1)
  J <- 9
  MSE <- 0.088
  my_comparison(means, J, MSE, alpha = 0.05, gamma = 0.8)
}

# ------------------------------------------------------------
# Run everything (optional)
# ------------------------------------------------------------
run_all <- function() {
  cat("=== Benjamin method demo (oil filters) ===\\n")
  print(demo_oil_filters())

  cat("\\n=== FWER simulations (n = 10,000) ===\\n")
  scenarios <- data.frame(
    scenario = c("A","B","C","D","E"),
    I = c(3,5,7,5,5),
    J = c(10,10,10,5,20)
  )
  scenarios$fwer <- NA_real_

  for (r in 1:nrow(scenarios)) {
    scenarios$fwer[r] <- simulate_fwer(
      I = scenarios$I[r],
      J = scenarios$J[r],
      n_sim = 10000,
      sigma = 1,
      alpha = 0.05,
      gamma = 0.8,
      seed = 100 + r
    )
  }
  print(scenarios)

  cat("\\n=== Power vs Tukey (n = 1,000 each delta) ===\\n")
  pow <- simulate_power(deltas = c(0.5,1.0,1.5,2.0), n_sim = 1000, alpha = 0.05, gamma = 0.8, seed = 999)
  print(pow)

  # quick plot
  plot(pow$delta, pow$power_benjamin, type = "b", xlab = "Effect size delta",
       ylab = "Power: P(reject mu4 = mu1)", ylim = c(0,1))
  lines(pow$delta, pow$power_tukey, type = "b")
  legend("topleft", legend = c("Benjamin", "Tukey HSD"), lty = 1, pch = 1)
}
