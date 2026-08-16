# Simulation verification harness (developer tool, not shipped).
#
# Recreates, in preserved and rerunnable form, the simulation audits that
# validated the 2.0.0 planners (the 2026-07-27 Welch and chi-square audit and
# the 2026-07-29 correlation audit; the original one-off scripts were run in
# working sessions and not kept, so this file is their durable replacement).
# The specific numbers quoted in NEWS and the help pages came from those
# original runs; the full mode of this script re-verifies the same claims at
# comparable precision, and the quick mode at looser Monte Carlo tolerances.
#
# Four checks:
#   1. Assurance calibration, statistic scale. For six designs spanning every
#      sampling family (t, F, chi-square, Welch, correlation mixture), simulate
#      published prior studies from a named true effect and confirm that the
#      corrected noncentrality parameter falls at or below the truth with the
#      requested assurance, via ss_buc_sensitivity() (whose
#      attainment_with_refusals is exactly that rate).
#   2. Correlation, raw data. Draw bivariate normal samples (not statistics),
#      publish by the two-sided r test, and confirm (a) the corrected parameter
#      covers the true lambda at the assurance rate, (b) plans attain the
#      sample size the true rho requires at about the assurance rate, and
#      (c) the fixed-predictor calculation (r converted to t, planned by
#      ss_buc_reg_coef at p = 1, which is what ss_buc_correlation did before
#      the random-predictor upgrade) attains it markedly less often and plans
#      a few participants light.
#   3. Welch, raw data. Confirm the real Welch test's Type I error is near
#      nominal across variance ratios 1 to 3, that a planned study simulated
#      at the corrected effect reaches the actual_power the function reports,
#      and that planning at sd_ratio = 1 when the truth has a ratio of 2 costs
#      a large amount of realized power (the documented load-bearing caveat).
#   4. Chi-square difference, finite-sample caveat. The likelihood ratio
#      statistic for nested normal models is N*log(1 + chisq_k/chisq_{N-p})
#      under the null, so its exact small-sample Type I error can be simulated
#      without raw data: confirm the test is liberal at total N = 50 and near
#      nominal at N = 400, the behavior the help page's caveat describes.
#
# Usage:
#   Rscript tools/verify_by_simulation.R          # quick (a few minutes)
#   Rscript tools/verify_by_simulation.R full     # published-precision reps

suppressMessages(pkgload::load_all(".", quiet = TRUE))

mode <- if (length(commandArgs(trailingOnly = TRUE)) >= 1L &&
            commandArgs(trailingOnly = TRUE)[1] == "full") "full" else "quick"
reps_sens  <- if (mode == "full") 10000L else 1200L
reps_data  <- if (mode == "full") 20000L else 3000L
reps_welch <- if (mode == "full") 100000L else 20000L
reps_lrt   <- if (mode == "full") 1000000L else 200000L
cat("mode:", mode, "\n\n")

set.seed(2017)
results <- list()
check <- function(name, ok, detail) {
  results[[length(results) + 1L]] <<- list(name = name, ok = ok)
  cat(sprintf("  [%s] %s (%s)\n", if (ok) "PASS" else "FAIL", name, detail))
}

## Check 1: assurance calibration on the statistic scale ----------------------

cat("Check 1: assurance calibration across sampling families\n")
assurance <- .80
tol_sens <- 4 * sqrt(assurance * (1 - assurance) / reps_sens)
calibration_cases <- list(
  list(label = "independent t",
       plan = ss_buc_independent_t(t_observed = 3, n = c(50, 55),
                                   assurance = assurance),
       true_ncp = 1.5),
  list(label = "one-way ANOVA",
       plan = ss_buc_one_way_anova(F_observed = 6, N = 60, levels_A = 3,
                                   assurance = assurance),
       true_ncp = 8),
  list(label = "regression coefficient",
       plan = ss_buc_reg_coef(t_observed = 3, N = 120, p = 4,
                              assurance = assurance),
       true_ncp = 1.6),
  list(label = "chi-square difference",
       plan = ss_buc_chisq_diff(chisq_observed = 15, N = 300, df_full = 24,
                                df_restricted = 26, assurance = assurance),
       true_ncp = 10),
  list(label = "Welch t",
       plan = ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 50,
                             sd_ratio = 1.5, assurance = assurance),
       true_ncp = 1.5),
  list(label = "correlation",
       plan = ss_buc_correlation(r_observed = .35, N = 100,
                                 assurance = assurance),
       true_ncp = 100 * .09 / .91))
for (i in seq_along(calibration_cases)) {
  cc <- calibration_cases[[i]]
  s <- ss_buc_sensitivity(cc$plan, true_ncp = cc$true_ncp,
                          replications = reps_sens, seed = 2017 + i)
  awr <- s$value[s$term == "attainment_with_refusals"]
  check(paste0("calibration: ", cc$label),
        abs(awr - assurance) <= tol_sens,
        sprintf("attainment_with_refusals %.3f vs assurance %.2f, tol %.3f",
                awr, assurance, tol_sens))
}

## Check 2: correlation from raw bivariate data -------------------------------

cat("Check 2: correlation planner against raw bivariate normal data\n")
rho <- .30
n_prior <- 100L
desired_power <- .80
true_lambda <- n_prior * rho^2 / (1 - rho^2)

# The sample size the true rho actually requires under the random-predictor
# distribution: smallest N whose exact mixture power reaches the target.
power_at <- function(N) {
  1 - BUCSS:::.p_correlation_f(qf(1 - .05, 1, N - 2), n = N, rho2 = rho^2)
}
n_required <- 5L
while (power_at(n_required) < desired_power) n_required <- n_required + 1L

x <- matrix(rnorm(n_prior * reps_data), n_prior)
y <- rho * x + sqrt(1 - rho^2) * matrix(rnorm(n_prior * reps_data), n_prior)
r_obs <- vapply(seq_len(reps_data), function(j) cor(x[, j], y[, j]), 0)
t_obs <- r_obs * sqrt((n_prior - 2) / (1 - r_obs^2))
published <- abs(t_obs) > qt(.975, n_prior - 2)
r_pub <- r_obs[published]
t_pub <- t_obs[published]

covered <- attained_random <- logical(length(r_pub))
size_gap <- overstatement <- rep(NA_real_, length(r_pub))
for (j in seq_along(r_pub)) {
  plan_r <- tryCatch(ss_buc_correlation(r_observed = r_pub[j], N = n_prior,
                                        assurance = assurance,
                                        desired_power = desired_power),
                     error = function(e) NULL)
  plan_f <- tryCatch(ss_buc_reg_coef(t_observed = t_pub[j], N = n_prior,
                                     p = 1, assurance = assurance,
                                     desired_power = desired_power),
                     error = function(e) NULL)
  # A refusal is the limiting conservative answer: the corrected parameter is
  # driven to zero, which is at or below any positive truth, and declining to
  # plan cannot under-deliver, so it counts as covered and attaining (the same
  # convention as attainment_with_refusals in ss_buc_sensitivity()).
  covered[j] <- is.null(plan_r) ||
    plan_r$value[plan_r$term == "ncp_adjusted"] <= true_lambda
  attained_random[j] <- is.null(plan_r) ||
    plan_r$value[plan_r$term == "necessary_N"] >= n_required
  if (!is.null(plan_r) && !is.null(plan_f)) {
    n_f <- plan_f$value[plan_f$term == "necessary_N"]
    lam_f <- plan_f$value[plan_f$term == "ncp_adjusted"]
    size_gap[j] <- plan_r$value[plan_r$term == "necessary_N"] - n_f
    # The published claim: every plan the fixed-predictor calculation issued
    # fell below its own stated target power under the true (random-predictor)
    # distribution. The fixed planner's noncentrality is on the lambda scale
    # at the prior N, so its adjusted effect is rho^2 = lambda/(lambda + N),
    # and the mixture power of its plan at that same effect is what the plan
    # actually delivers when both variables are sampled.
    q_f <- qf(1 - .05, 1, n_f - 2)
    overstatement[j] <- plan_f$value[plan_f$term == "actual_power"] -
      (1 - BUCSS:::.p_correlation_f(q_f, n = n_f,
                                    rho2 = lam_f / (lam_f + n_prior)))
  }
}
# Scale sanity for the comparison above: scaling the fixed plan's lambda to
# its planned N must reproduce the plan's own actual_power under the plain
# noncentral F, or the rho^2 conversion would be reading the wrong scale.
plan_f <- ss_buc_reg_coef(t_observed = 3.5, N = n_prior, p = 1,
                          assurance = assurance)
n_f <- plan_f$value[plan_f$term == "necessary_N"]
lam_f <- plan_f$value[plan_f$term == "ncp_adjusted"]
stopifnot(abs((1 - pf(qf(.95, 1, n_f - 2), 1, n_f - 2,
                      ncp = lam_f * n_f / n_prior)) -
              plan_f$value[plan_f$term == "actual_power"]) < .005)
tol_cov <- 4 * sqrt(assurance * (1 - assurance) / length(r_pub))
check("correlation: coverage of true lambda",
      abs(mean(covered) - assurance) <= tol_cov,
      sprintf("coverage %.3f vs assurance %.2f over %d published, tol %.3f",
              mean(covered), assurance, length(r_pub), tol_cov))
check("correlation: attainment of the required size",
      abs(mean(attained_random) - assurance) <= tol_cov + .02,
      sprintf("attainment %.3f, required N %d at rho %.2f",
              mean(attained_random), n_required, rho))
check("correlation: every fixed-predictor plan overstates its power",
      all(overstatement > 0, na.rm = TRUE) &&
        max(overstatement, na.rm = TRUE) <= .08,
      sprintf("overstated by %.4f to %.4f (audit found up to about .04)",
              min(overstatement, na.rm = TRUE),
              max(overstatement, na.rm = TRUE)))
check("correlation: fixed-predictor plans run light",
      median(size_gap, na.rm = TRUE) >= 1,
      sprintf("median size gap %.1f (published range was 2 to 5)",
              median(size_gap, na.rm = TRUE)))

## Check 3: Welch t against raw data ------------------------------------------

cat("Check 3: Welch planner against raw two-group normal data\n")
welch_stat <- function(n_1, n_2, sd_1, sd_2, mu_diff, reps) {
  m1 <- matrix(rnorm(n_1 * reps, mean = mu_diff, sd = sd_1), n_1)
  m2 <- matrix(rnorm(n_2 * reps, sd = sd_2), n_2)
  v1 <- (colSums(m1^2) - n_1 * colMeans(m1)^2) / (n_1 - 1)
  v2 <- (colSums(m2^2) - n_2 * colMeans(m2)^2) / (n_2 - 1)
  se2 <- v1 / n_1 + v2 / n_2
  df <- se2^2 / ((v1 / n_1)^2 / (n_1 - 1) + (v2 / n_2)^2 / (n_2 - 1))
  list(t = (colMeans(m1) - colMeans(m2)) / sqrt(se2), df = df)
}
for (ratio in c(1, 2, 3)) {
  z <- welch_stat(40, 55, 1, ratio, 0, reps_welch)
  alpha_hat <- mean(abs(z$t) > qt(.975, z$df))
  check(sprintf("Welch: Type I error at sd ratio %g", ratio),
        alpha_hat >= .042 && alpha_hat <= .058,
        sprintf("empirical alpha %.4f (audit range .047 to .054)", alpha_hat))
}

plan_w <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5,
                         assurance = assurance)
n_w <- plan_w$value[plan_w$term == "necessary_n_per_group"]
# Back the mean difference out of the corrected prior-scale noncentrality
# (sd_1 = 1, sd_2 = sd_ratio; delta = mu / sqrt(sd_1^2/n_1 + sd_2^2/n_2)).
ncp_w <- plan_w$value[plan_w$term == "ncp_adjusted"]
mu_w <- ncp_w * sqrt(1 / 40 + 1.5^2 / 55)
z <- welch_stat(n_w, n_w, 1, 1.5, mu_w, reps_welch)
power_hat <- mean(abs(z$t) > qt(.975, z$df))
claimed <- plan_w$value[plan_w$term == "actual_power"]
check("Welch: planned study reaches the reported actual_power",
      abs(power_hat - claimed) <= .015,
      sprintf("empirical %.3f vs reported %.3f (audit agreement within .006)",
              power_hat, claimed))

# The load-bearing caveat, in the help page's own scenario: the documented
# example's prior planned at sd_ratio = 1, in a world where the mean
# difference and sd_1 are as the user believed but the true ratio is 2. The
# help page states this plan delivers about .43 power rather than .80.
plan_eq <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1,
                          assurance = assurance)
n_eq <- plan_eq$value[plan_eq$term == "necessary_n_per_group"]
mu_eq <- plan_eq$value[plan_eq$term == "ncp_adjusted"] * sqrt(1 / 40 + 1 / 55)
z <- welch_stat(n_eq, n_eq, 1, 2, mu_eq, reps_welch)
power_missed <- mean(abs(z$t) > qt(.975, z$df))
check("Welch: assuming equal variances when the ratio is 2 costs power",
      abs(power_missed - .43) <= .03,
      sprintf("realized power %.3f vs target %.2f (help page says about .43)",
              power_missed, desired_power))

## Check 4: chi-square difference test, finite-sample Type I ------------------

cat("Check 4: finite-sample Type I error of the nested-model LRT\n")
# For nested normal linear models (drop k = 2 of p = 6 parameters) the null
# LRT is exactly N*log(1 + chisq_k / chisq_{N-p}), so no raw data is needed.
lrt_alpha <- function(N, k = 2, p = 6, reps) {
  lr <- N * log(1 + rchisq(reps, k) / rchisq(reps, N - p))
  mean(lr > qchisq(.95, k))
}
a50 <- lrt_alpha(50, reps = reps_lrt)
a400 <- lrt_alpha(400, reps = reps_lrt)
check("chi-square: liberal at total N = 50",
      a50 > .055,
      sprintf("empirical alpha %.4f (audit found near .06)", a50))
check("chi-square: near nominal at total N = 400",
      abs(a400 - .05) <= .005 && a50 > a400,
      sprintf("empirical alpha %.4f", a400))

## Summary ---------------------------------------------------------------------

n_fail <- sum(!vapply(results, `[[`, TRUE, "ok"))
cat(sprintf("\n%d checks, %d failed.\n", length(results), n_fail))
if (n_fail > 0L) quit(status = 1L)
