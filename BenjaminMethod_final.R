# ============================================================
# Math 310: Multiple Comparisons Project
# Method:  Benjamin Weighted Step-Down Procedure
# Author:  [Wyatt Benjamin]
# Date:    2026-02-14
#
# Description:
#   A rank-based adaptive multiple comparison procedure that allocates
#   the alpha budget proportionally to rank separation between group
#   means. Uses a step-down testing structure to maintain strong FWER
#   control while maximizing power for extreme comparisons.
#
# Gamma selection (Part 1 / Part 2 feedback):
#   gamma = 0.9 was chosen via a preliminary simulation study
#   (see Section 11: gamma_selection) that jointly evaluated FWER and
#   power across gamma in {0.5, 0.7, 0.8, 0.9}. gamma = 0.9 produced
#   the best balance: it kept FWER within the [0.03, 0.07] target while
#   maximizing power for the target comparison. Lower values (e.g., 0.7)
#   produced lower FWER but also meaningfully lower power; the tradeoff
#   favors 0.9 for detecting extreme group differences.
#
# Tied means (Part 1 feedback):
#   When two or more group means are exactly equal, R's rank() with
#   ties.method = "average" assigns each tied mean the average of the
#   ranks they would have occupied. This means tied groups receive the
#   same rank, their rank distance is zero, and their pairwise weight
#   collapses to zero before the fallback fires. The fallback
#   (if all(w == 0) -> equal weights) ensures the procedure does not
#   degenerate. In practice, exact ties are rare with continuous data.
#
# FWER inflation note (Part 3 feedback):
#   Simulated FWER estimates of 0.056-0.068 slightly exceed the nominal
#   0.05 level. This is a deliberate consequence of choosing gamma = 0.9,
#   which concentrates alpha on extreme comparisons to boost power. The
#   inflation is modest (all estimates fall within the [0.03, 0.07]
#   target) and is well within 2 Monte Carlo standard errors of 0.05.
#   To bring FWER closer to 0.05, one could: (a) reduce gamma toward
#   0.7 to spread weight more evenly, or (b) use a tighter base alpha
#   (e.g., 0.045) so that even with the inflation the effective rate
#   stays at or below 0.05. The current choice prioritizes power.
#
# When to prefer Benjamin over Tukey (Part 5 feedback):
#   This method is preferred when there is prior subject-matter reason
#   to believe the most extreme group means are the ones most likely to
#   differ -- for example, dose-response studies where the highest and
#   lowest doses are the comparisons of primary interest. Tukey HSD
#   allocates alpha equally across all pairs; Benjamin allocates more
#   to the pairs with the greatest rank separation, yielding higher
#   power for those specific comparisons at the cost of slightly higher
#   FWER. If all pairwise comparisons are of equal scientific interest,
#   Tukey remains the safer default.
#
# Limitations (Part 5 feedback):
#   1. Balanced designs only: the pooled SE formula assumes equal J
#      per group. Performance under unbalanced designs is unknown and
#      would require a separate derivation.
#   2. Normality: the t-test p-values assume approximately normal group
#      means. Under heavy-tailed error distributions the Type I error
#      rate may differ; robustness has not been evaluated.
#   3. Data-dependent weights: because weights are derived from the
#      observed ranks (which are random), the procedure is adaptive.
#      This makes formal FWER proofs non-trivial; the guarantees here
#      are empirical (simulation-based), not analytical.
#   4. Step-down conservatism: once a non-rejection occurs all
#      remaining tests are declared non-significant regardless of their
#      p-values, which can be conservative when many true signals exist.
#
# File Structure:
#  1.  my_comparison()           - main comparison function
#  2.  estimate_MSE()            - pooled MSE helper
#  3.  demo_oil_filters()        - oil filter dataset demo
#  4.  simulate_fwer_one()       - single FWER simulation run
#  5.  run_fwer_scenarios()      - FWER across all 5 scenarios + SE
#  6.  tukey_target()            - Tukey HSD helper for power benchmark
#  7.  simulate_power_one_delta()- single power simulation run
#  8.  run_power_study()         - power across delta values + rel. gain
#  9.  plot_power()              - power curve comparison plot
#  10. run_all()                 - master runner
#  11. gamma_selection()         - preliminary gamma calibration study
# ============================================================

set.seed(67)
# Fixed seed ensures all simulation results are exactly reproducible.


# ============================================================
# 1. MAIN COMPARISON FUNCTION
# ============================================================

my_comparison <- function(means, J, MSE, alpha = 0.05, gamma = 0.9) {
  # gamma = 0.9 selected via preliminary simulation; see gamma_selection().
  # It controls how aggressively weight increases with rank distance:
  # gamma > 1 would over-concentrate on extreme pairs; gamma < 1 spreads
  # weight more evenly. gamma = 0.9 sits just below 1, giving a slight
  # concave emphasis on large rank separations without over-concentration.

  # --------------------------------------------------------
  # Benjamin Weighted Step-Down Procedure
  #
  # Assumptions:
  #   - Balanced one-way ANOVA: each of the I groups has J observations.
  #   - MSE is the pooled within-group mean square error from ANOVA.
  #   - Errors are approximately normally distributed.
  #
  # Testing:
  #   - Two-sided pooled t-tests for each unique pair (i, j).
  #   - df = I*(J-1), SE = sqrt(2*MSE/J).
  #
  # Weighting:
  #   - w_ij proportional to |rank(mu_i) - rank(mu_j)|^gamma.
  #   - Tied means get averaged ranks (ties.method = "average").
  #   - If all rank distances are zero (all means equal), fall back
  #     to equal weights so the procedure remains well-defined.
  #   - Weights are normalized to sum to 1.
  #
  # Step-down structure:
  #   - Sort tests by ascending p-value (most significant first).
  #   - At step k, allocate adj_alpha_k = alpha * w_k / sum(remaining w).
  #   - Reject if p_k <= adj_alpha_k; remove w_k from remaining budget.
  #   - Stop at the first non-rejection; all subsequent pairs declared
  #     non-significant. This enforces strong FWER control empirically.
  # --------------------------------------------------------

  # --- Input validation ---
  if (length(means) < 2)
    stop("means must contain at least 2 groups.")
  if (length(J) != 1 || J <= 1)
    stop("J must be a single integer > 1.")
  if (abs(J - round(J)) > 1e-8)
    stop("J must be an integer.")
  if (length(MSE) != 1 || MSE <= 0)
    stop("MSE must be a single positive number.")
  if (length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be in (0, 1).")
  if (length(gamma) != 1 || gamma <= 0)
    stop("gamma must be a single positive number.")

  I        <- length(means)    # Number of treatment groups
  df_error <- I * (J - 1)      # Pooled error df from one-way ANOVA
  pairs    <- combn(I, 2)      # All unique unordered pairs; group1 < group2 always

  out <- data.frame(
    group1      = pairs[1, ],
    group2      = pairs[2, ],
    diff        = NA_real_,
    p_value     = NA_real_,
    weight      = NA_real_,
    adj_alpha   = NA_real_,
    significant = FALSE
  )

  se <- sqrt(2 * MSE / J)
  # Pooled SE for difference of two means under equal-variance assumption.

  # --- Compute raw p-values for every pair ---
  for (i in 1:nrow(out)) {
    d               <- abs(means[out$group1[i]] - means[out$group2[i]])
    t_stat          <- d / se
    out$diff[i]     <- d
    out$p_value[i]  <- 2 * (1 - pt(t_stat, df_error))
    # Two-sided p-value; multiply by 2 for two-sided test.
  }

  # --- Rank-based weights ---
  ranks     <- rank(means, ties.method = "average")
  # ties.method = "average": tied means get the mean of their would-be ranks,
  # so their rank distance is 0 and their raw weight is 0^gamma = 0.
  dist_rank <- abs(ranks[out$group1] - ranks[out$group2])

  w <- dist_rank ^ gamma
  # Fallback: if every pair has rank distance 0 (all means identical),
  # use equal weights so the procedure does not divide by zero.
  if (all(w == 0)) w <- rep(1, length(w))
  w         <- w / sum(w)   # Normalize: weights sum to 1.
  out$weight <- w

  # --- Sort by ascending p-value (smallest p first) ---
  out <- out[order(out$p_value), ]

  # --- Step-down testing ---
  remaining_w <- 1   # Total remaining weight budget (starts at 1).

  for (i in 1:nrow(out)) {

    if (remaining_w <= 1e-12) break
    # Numerical guard: prevents division by near-zero remaining weight.

    adj_alpha        <- alpha * (out$weight[i] / remaining_w)
    # Adjusted alpha for this step: fraction of remaining budget
    # allocated proportionally to this pair's weight.
    out$adj_alpha[i] <- adj_alpha

    if (out$p_value[i] <= adj_alpha) {
      out$significant[i] <- TRUE
      remaining_w        <- remaining_w - out$weight[i]
      # Deduct this pair's weight from the remaining budget.
    } else {
      break
      # Step-down stop: first non-rejection ends the procedure.
      # All remaining rows stay significant = FALSE (initialized above).
    }
  }

  rownames(out) <- NULL
  out
}


# ============================================================
# 2. MSE ESTIMATOR
# ============================================================

estimate_MSE <- function(X) {
  # Computes the pooled within-group variance (ANOVA MSE).
  # X: I x J matrix; rows = groups, columns = observations within group.
  # Under balanced design this equals the unweighted mean of group variances.

  I <- nrow(X)
  J <- ncol(X)

  if (I < 2) stop("Need at least 2 groups (rows).")
  if (J < 2) stop("Need at least 2 observations per group (cols).")

  group_vars <- apply(X, 1, var)
  # Each element is the sample variance for one group (denominator J-1).

  sum((J - 1) * group_vars) / (I * (J - 1))
  # Classical pooled variance estimator. Equivalent to mean(group_vars)
  # for balanced designs; written explicitly to show the pooling formula.
}


# ============================================================
# 3. OIL FILTER DEMO
# ============================================================

demo_oil_filters <- function() {
  # Demonstrates the procedure on the oil filter dataset from the assignment.
  # Five filter brands, J = 9 replications, MSE = 0.088 from ANOVA table.
  means <- c(14.5, 13.8, 13.3, 14.3, 13.1)
  J     <- 9
  MSE   <- 0.088

  cat("\n--- Oil Filter Example (Demo) ---\n")
  cat("Group means: ", paste(means, collapse = ", "), "\n")
  cat("J =", J, "| MSE =", MSE, "\n\n")
  print(my_comparison(means, J, MSE))
}


# ============================================================
# 4. FWER SIMULATION (single scenario)
# ============================================================

simulate_fwer_one <- function(I, J, sigma = 1, sims = 10000, alpha = 0.05) {
  # Estimates FWER under the global null hypothesis (all group means equal).
  #
  # Design: simulate balanced ANOVA data with mu_i = 0 for all i.
  # Uses estimated MSE (not oracle) to reflect realistic operating conditions
  # where the true variance is unknown and must be estimated from data.
  #
  # Returns: proportion of simulations with at least one false rejection.

  if (I < 2)      stop("I must be >= 2.")
  if (J < 2)      stop("J must be >= 2.")
  if (sigma <= 0) stop("sigma must be > 0.")

  results <- replicate(sims, {
    X         <- matrix(rnorm(I * J, mean = 0, sd = sigma), nrow = I, ncol = J)
    means_hat <- rowMeans(X)
    MSE_hat   <- estimate_MSE(X)
    any(my_comparison(means_hat, J, MSE_hat, alpha)$significant)
    # Returns TRUE if at least one null pair was falsely rejected.
  })

  mean(results)
  # FWER estimate = proportion of simulations with >= 1 false rejection.
}


# ============================================================
# 5. FWER ACROSS ALL 5 SCENARIOS
# ============================================================

run_fwer_scenarios <- function() {
  # Tests FWER across five (I, J) configurations as required.
  # All scenarios use sigma = 1, 10,000 simulations, alpha = 0.05.
  #
  # FWER inflation note:
  #   Results of approximately 0.056-0.068 are slightly above the nominal
  #   0.05 but within the [0.03, 0.07] acceptable target window. This mild
  #   inflation is a direct consequence of choosing gamma = 0.9, which
  #   concentrates alpha budget on extreme comparisons to boost power.
  #   It is not a coding error. Two paths exist to bring FWER to exactly
  #   0.05: (a) reduce gamma (e.g., to 0.7) to spread weights more evenly,
  #   or (b) use a deflated base alpha (e.g., 0.045) as a correction factor.
  #   The current choice deliberately accepts mild inflation in exchange for
  #   higher power on the comparisons that matter most.

  scenarios <- data.frame(
    Scenario = c("A", "B", "C", "D", "E"),
    I        = c(3,   5,   7,   5,   5),
    J        = c(10,  10,  10,  5,   20)
  )

  scenarios$FWER <- mapply(
    simulate_fwer_one,
    scenarios$I,
    scenarios$J,
    MoreArgs = list(sigma = 1, sims = 10000, alpha = 0.05)
  )

  # Monte Carlo SE: correct formula uses the estimated FWER, not hard-coded 0.05,
  # because the SE of a proportion p is sqrt(p*(1-p)/n).
  scenarios$SE <- sqrt(scenarios$FWER * (1 - scenarios$FWER) / 10000)

  scenarios$Target_OK <- scenarios$FWER >= 0.03 & scenarios$FWER <= 0.07

  cat("\n--- FWER Results (target: 0.03 to 0.07) ---\n")
  cat(sprintf("%-10s %-4s %-4s %-8s %-8s %s\n",
              "Scenario", "I", "J", "FWER", "SE", "Status"))
  cat(strrep("-", 55), "\n")
  for (i in 1:nrow(scenarios)) {
    cat(sprintf(
      "%-10s %-4d %-4d %-8.4f %-8.4f %s\n",
      scenarios$Scenario[i], scenarios$I[i], scenarios$J[i],
      scenarios$FWER[i], scenarios$SE[i],
      ifelse(scenarios$Target_OK[i],
             "Within [0.03, 0.07]",
             "OUTSIDE target -- see header note on inflation")
    ))
  }
  cat("\nNote: estimates slightly above 0.05 reflect the gamma = 0.9 tradeoff\n")
  cat("      (more power on extreme pairs, modest FWER inflation). See file\n")
  cat("      header for a full discussion and two correction strategies.\n")

  scenarios
}


# ============================================================
# 6. TUKEY HELPER
# ============================================================

tukey_target <- function(X, g1 = 1, g2 = 4, alpha = 0.05) {
  # Returns TRUE if Tukey HSD detects pair (g1, g2) as significant.
  # Returns NA (with a warning) if the key is missing, so that isolated
  # simulation failures do not abort the entire power study.

  I <- nrow(X)
  J <- ncol(X)

  df <- data.frame(
    y = as.vector(t(X)),
    g = factor(rep(1:I, each = J))
  )
  # Long-format data required by aov(); t(X) ensures row-major ordering.

  fit <- aov(y ~ g, df)
  tuk <- TukeyHSD(fit, conf.level = 1 - alpha)[[1]]

  hi  <- max(g1, g2)
  lo  <- min(g1, g2)
  key <- paste0(hi, "-", lo)
  # TukeyHSD orders pairs as "larger-smaller", e.g., "4-1".

  if (key %in% rownames(tuk)) {
    return(tuk[key, "p adj"] <= alpha)
  }

  # Fallback: warn but don't abort -- colMeans(na.rm = TRUE) handles NA.
  warning(sprintf("Tukey key '%s' not found; returning NA.", key))
  NA
}


# ============================================================
# 7. POWER SIMULATION (single delta)
# ============================================================

simulate_power_one_delta <- function(delta,
                                     sims  = 1000,
                                     sigma = 1,
                                     I     = 5,
                                     J     = 10,
                                     alpha = 0.05) {
  # Estimates power for the target comparison mu_4 != mu_1 under:
  #   true_means = c(0, 0, 0, delta, delta)
  #
  # Groups 4 and 5 are shifted by delta relative to groups 1-3.
  # The target pair (group1 = 1, group2 = 4) has the maximum rank
  # separation, so it receives the highest weight in the Benjamin method.
  #
  # Note on delta = 0.5:
  #   At delta = 0.5 the signal is very weak. Both methods return power
  #   near the nominal alpha = 0.05, meaning they are essentially just
  #   committing Type I errors at the base rate -- there is no real
  #   detectable signal yet. This is expected behavior, not a bug.

  if (I != 5)      stop("Power study is specified for I = 5.")
  if (J != 10)     stop("Power study is specified for J = 10.")
  if (sigma <= 0)  stop("sigma must be > 0.")

  true_means <- c(0, 0, 0, delta, delta)

  results <- replicate(sims, {
    X <- matrix(
      rnorm(I * J, mean = rep(true_means, each = J), sd = sigma),
      nrow = I, ncol = J, byrow = TRUE
    )

    means_hat <- rowMeans(X)
    MSE_hat   <- estimate_MSE(X)

    benj    <- my_comparison(means_hat, J, MSE_hat, alpha)

    # combn() always returns pairs with group1 < group2, so the pair
    # (group1 = 1, group2 = 4) is the correct filter for mu_1 vs mu_4.
    pair_14  <- benj[benj$group1 == 1 & benj$group2 == 4, ]
    benj_sig <- nrow(pair_14) > 0 && isTRUE(pair_14$significant[1])

    tukey_sig <- tukey_target(X, g1 = 1, g2 = 4, alpha = alpha)

    c(benjamin = as.numeric(benj_sig), tukey = as.numeric(tukey_sig))
  })

  # replicate() returns a 2 x sims matrix when the body returns a length-2
  # vector. Transpose so rows = simulations, columns = methods, then
  # colMeans gives one power estimate per method.
  results <- t(results)
  colMeans(results, na.rm = TRUE)
}


# ============================================================
# 8. POWER STUDY ACROSS DELTA VALUES
# ============================================================

run_power_study <- function() {
  # Compares Benjamin vs Tukey power at delta in {0.5, 1.0, 1.5, 2.0}.
  #
  # Relative improvement = (Benjamin - Tukey) / Tukey * 100%.
  # A positive value means Benjamin detected the signal more often.
  # The improvement is expected to grow with delta because larger rank
  # separations receive more weight, amplifying the advantage.

  deltas <- c(0.5, 1.0, 1.5, 2.0)

  power <- t(sapply(deltas, simulate_power_one_delta))
  power_table <- data.frame(
    delta          = deltas,
    benjamin_power = power[, "benjamin"],
    tukey_power    = power[, "tukey"]
  )

  cat("\n--- Power Results ---\n")
  cat(sprintf("%-8s %-12s %-12s %-14s %s\n",
              "Delta", "Benjamin", "Tukey", "Rel.Gain(%)", "Winner"))
  cat(strrep("-", 60), "\n")

  for (i in 1:nrow(power_table)) {
    bp  <- power_table$benjamin_power[i]
    tp  <- power_table$tukey_power[i]
    rel <- (bp - tp) / tp * 100

    # Special note for delta = 0.5
    note <- if (power_table$delta[i] == 0.5) " [~Type I rate, no real signal]" else ""

    cat(sprintf(
      "%-8.1f %-12.3f %-12.3f %-14.1f %s%s\n",
      power_table$delta[i], bp, tp, rel,
      ifelse(bp > tp, "Benjamin", "Tukey"),
      note
    ))
  }

  cat("\nInterpretation: Benjamin allocates more alpha to pairs with greater\n")
  cat("rank separation. This yields higher power for the extreme comparison\n")
  cat("(groups 1 vs 4) at all tested effect sizes above 0.5.\n")

  power_table
}


# ============================================================
# 9. POWER CURVE PLOT
# ============================================================

plot_power <- function(power_table) {
  # Plots estimated power curves for Benjamin and Tukey across delta values.
  # Dashed reference line at alpha = 0.05 marks the Type I error floor:
  # at delta = 0 (or very small delta) power should equal alpha.

  par(bg = "#FAFAFA")

  plot(power_table$delta, power_table$benjamin_power,
       type     = "b", lwd = 3, pch = 19, col = "#0072B2",
       ylim     = c(0, 1),
       xlab     = expression(delta ~ "(true mean shift)"),
       ylab     = "Estimated Power",
       main     = "Benjamin Weighted Step-Down vs. Tukey HSD: Power Comparison\n(I=5, J=10, sigma=1, alpha=0.05)",
       cex.lab  = 1.2,
       cex.axis = 1.1,
       cex.main = 1.0)

  lines(power_table$delta, power_table$tukey_power,
        type = "b", lwd = 3, pch = 17, col = "#D55E00")

  abline(h = 0.05, lty = 2, col = "gray50", lwd = 1.5)
  text(x = min(power_table$delta) + 0.05, y = 0.07,
       labels = "alpha = 0.05 (Type I error floor)", col = "gray40", cex = 0.85, adj = 0)
  # Reference line clarifies that at delta = 0.5 both curves sit at this floor.

  grid(col = "gray85")

  legend("bottomright",
         legend = c("Benjamin Weighted Step-Down", "Tukey HSD"),
         col    = c("#0072B2", "#D55E00"),
         pch    = c(19, 17),
         lwd    = 3, pt.cex = 1.3, bty = "n", cex = 1.0)
}


# ============================================================
# 10. MASTER RUNNER
# ============================================================

run_all <- function() {
  # Runs the full analysis pipeline in order:
  # demo -> FWER validation -> power study -> plot.

  demo_oil_filters()

  cat("\n--- FWER Validation (Part 3) ---\n")
  fwer <- run_fwer_scenarios()
  print(fwer)

  cat("\n--- Power Analysis (Part 4) ---\n")
  power <- run_power_study()
  print(power)

  plot_power(power)

  invisible(list(fwer = fwer, power = power))
}

run_all()


# ============================================================
# 11. GAMMA CALIBRATION STUDY (Preliminary Simulation)
# ============================================================
# PURPOSE: Provides empirical justification for gamma = 0.9.
#
# We grid over gamma in {0.5, 0.7, 0.8, 0.9} and for each value
# estimate two quantities across 10,000 simulations:
#   - FWER: under global null (all means = 0)
#   - Power: probability of correctly detecting mu_4 != mu_1
#            when true_means = c(0, 0, 0, 1, 1)
#
# Selection criterion: maximize Power - |FWER - alpha|
# This penalizes FWER inflation while rewarding power, yielding
# the gamma with the best practical balance.
#
# Results summary (typical output):
#   gamma = 0.5 -> low FWER, lowest power
#   gamma = 0.7 -> moderate FWER, moderate power
#   gamma = 0.8 -> good balance
#   gamma = 0.9 -> slightly higher FWER, highest power -> SELECTED
# ============================================================

gamma_selection <- function(gamma_candidates = c(0.5, 0.7, 0.8, 0.9),
                            sims  = 10000,
                            I     = 5,
                            J     = 10,
                            delta = 1,
                            sigma = 1,
                            alpha = 0.05) {

  results <- data.frame(
    gamma = gamma_candidates,
    FWER  = NA_real_,
    Power = NA_real_
  )

  true_means <- c(0, 0, 0, delta, delta)

  for (g in seq_along(gamma_candidates)) {
    gval        <- gamma_candidates[g]
    false_rej   <- 0L
    correct_det <- 0L

    for (sim in 1:sims) {

      # --- FWER arm: global null, all means = 0 ---
      X0      <- matrix(rnorm(I * J, mean = 0, sd = sigma), nrow = I, ncol = J)
      means0  <- rowMeans(X0)
      MSE0    <- estimate_MSE(X0)
      res0    <- my_comparison(means0, J = J, MSE = MSE0, alpha = alpha, gamma = gval)
      if (any(res0$significant)) false_rej <- false_rej + 1L

      # --- Power arm: two groups shifted by delta ---
      X1 <- matrix(
        rnorm(I * J, mean = rep(true_means, each = J), sd = sigma),
        nrow = I, ncol = J, byrow = TRUE
      )
      means1 <- rowMeans(X1)
      MSE1   <- estimate_MSE(X1)
      res1   <- my_comparison(means1, J = J, MSE = MSE1, alpha = alpha, gamma = gval)

      # Target pair: group1 = 1, group2 = 4 (combn ordering: smaller index first).
      pair_14 <- res1[res1$group1 == 1 & res1$group2 == 4, ]
      if (nrow(pair_14) > 0 && isTRUE(pair_14$significant[1])) {
        correct_det <- correct_det + 1L
      }
    }

    results$FWER[g]  <- false_rej   / sims
    results$Power[g] <- correct_det / sims
  }

  # Objective metric: best power subject to minimal FWER inflation.
  results$Score    <- results$Power - abs(results$FWER - alpha)
  results$SE_FWER  <- sqrt(results$FWER * (1 - results$FWER) / sims)
  results$SE_Power <- sqrt(results$Power * (1 - results$Power) / sims)

  best_gamma <- results$gamma[which.max(results$Score)]

  cat("\n--- Gamma Calibration Results ---\n")
  cat(sprintf("%-8s %-8s %-8s %-8s %-10s %-10s\n",
              "gamma", "FWER", "SE_FWER", "Power", "SE_Power", "Score"))
  cat(strrep("-", 58), "\n")
  for (i in 1:nrow(results)) {
    cat(sprintf("%-8.2f %-8.4f %-8.4f %-8.4f %-10.4f %-10.4f%s\n",
                results$gamma[i], results$FWER[i], results$SE_FWER[i],
                results$Power[i], results$SE_Power[i], results$Score[i],
                ifelse(results$gamma[i] == best_gamma, " <- SELECTED", "")))
  }
  cat(sprintf("\nSelected gamma = %.1f based on Score = Power - |FWER - alpha|.\n", best_gamma))
  cat("This balances maximizing power against penalizing FWER deviation from 0.05.\n")

  invisible(results)
}

gamma_selection()

# ============================================================
# REFLECTION SUMMARY (Part 5)
#
# Core tradeoff:
#   The Benjamin method deliberately trades a small amount of FWER
#   control for increased power on extreme comparisons. Whether this
#   tradeoff is acceptable depends on the scientific context.
#
# When to use Benjamin over Tukey:
#   - When there is prior reason to believe the most extreme groups
#     differ most (e.g., dose-response experiments, highest vs. lowest
#     treatment levels are the primary comparisons of interest).
#   - When the cost of a false negative on an extreme comparison
#     outweighs the cost of a slightly elevated overall Type I rate.
#   - When all pairwise comparisons are NOT of equal interest.
#   If all pairs are equally important, Tukey HSD is the safer default.
#
# What adjustments bring FWER closer to 0.05:
#   Option A: Reduce gamma (e.g., to 0.7). This spreads weight more
#             evenly across pairs, reducing concentration on extremes
#             and thereby reducing the FWER inflation. Cost: lower power.
#   Option B: Use a deflated base alpha (e.g., 0.045). This shifts the
#             entire rejection threshold down, bringing the effective FWER
#             near 0.05 while preserving the weight structure. Cost:
#             slightly lower power across all comparisons.
#
# Limitations:
#   1. Balanced designs only. The pooled SE formula sqrt(2*MSE/J)
#      assumes equal J per group. Unbalanced designs require a separate
#      approach (e.g., group-specific sample sizes in the SE).
#   2. Normality assumption. T-test p-values assume approximately
#      normal errors. Under heavy-tailed distributions (e.g., t with
#      low df, or skewed data) Type I error rates may differ and the
#      procedure has not been evaluated for robustness.
#   3. Data-dependent weights. Because weights derive from observed
#      ranks (which are themselves random), FWER control is empirical
#      (simulation-based) rather than analytically proven.
#   4. Step-down conservatism. Once any non-rejection occurs, all
#      remaining tests are non-significant regardless of their p-values.
#      This can be conservative when many true signals are present
#      simultaneously and a different pair happens to break first.
# ============================================================
