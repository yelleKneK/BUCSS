# Simulation verification harness (ships with the package, in inst/tools/).
#
# Recreates, in preserved and rerunnable form, the simulation audits that
# validated the 2.0.0 planners (the 2026-07-27 Welch and chi square audit and
# the 2026-07-29 correlation audit; the original one-off scripts were run in
# working sessions and not kept, so this file is their durable replacement).
# Some numbers quoted in NEWS and the help pages came from those original
# runs, and others from the deep quality control audit and fix wave of
# 2026-09-28 and 2026-09-29 (the Welch accuracy figures, the correlation's 1/n
# law, the chi square Type I figures). This script re-verifies the claims its
# checks name, the full mode at published precision and the quick mode at
# looser Monte Carlo tolerances. Checks 3 (Type I error) and 4 verify the
# statistical facts the help pages state, not BUCSS code, so they do not call
# the package.
#
# Four checks:
#   1. Assurance calibration, statistic scale. For six designs spanning every
#      sampling family (t, F, chi square, Welch, correlation mixture), simulate
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
#      the random-predictor upgrade) overstates the power of every plan it
#      issues, by roughly 1/n for a plan of n, and plans a few participants
#      light.
#   3. Welch, raw data. Confirm the real Welch test's Type I error is near
#      nominal across variance ratios 1 to 3, that a planned study simulated
#      at the corrected effect reaches the actual_power the function reports,
#      and that planning at sd_ratio = 1 when the truth has a ratio of 2 costs
#      a large amount of realized power (the documented load-bearing caveat).
#   4. Chi-square difference, finite-sample caveat. The likelihood ratio
#      statistic for nested normal models is N*log(1 + chisq_k/chisq_{N-p})
#      under the null, so its small-sample Type I error can be simulated
#      without raw data, and it also has a closed form: confirm the two agree
#      at total N = 50, 200, and 400 and equal the figures the help page
#      quotes (.07, .055, and .052).
#
# This file ships with the package (inst/tools/), so anyone who installs
# BUCSS can rerun the verification without the source repository.
#
# Usage, from an installed package (quick mode runs in a few seconds; the
# full mode, with published-precision replications, in about half a minute):
#   source(system.file("tools", "verify_by_simulation.R", package = "BUCSS"))
#   Sys.setenv(BUCSS_VERIFY = "full")   # before source(), for the full mode
# or from the root of the BUCSS source repository, which loads that source
# with pkgload rather than the installed package:
#   Rscript inst/tools/verify_by_simulation.R          # quick
#   Rscript inst/tools/verify_by_simulation.R full     # full
#
# Sourcing this file leaves the session as it found it, apart from attaching
# BUCSS if it was not already attached. Everything runs inside
# local(), so nothing is created in, or overwritten in, the global environment;
# the seed is set with withr::local_seed(), so the caller's random number state
# and generator kind are restored when the run ends, even on an error. The run
# is pinned to R's default generators (Mersenne-Twister, Inversion, Rejection),
# so its numbers do not depend on an RNGkind() the caller may have set. A failed
# check ends the run with an error, which is a nonzero exit status under
# Rscript, and never quits the R session it was sourced into.

invisible(local({
  # The BUCSS source tree is recognized by its own DESCRIPTION, not by the
  # presence of any package's: a user sourcing this from inside another
  # package's directory gets the installed BUCSS.
  desc_package <- if (file.exists("DESCRIPTION")) tryCatch(
    unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
    error = function(e) NA_character_) else NA_character_
  in_bucss_source <- dir.exists("R") && identical(desc_package, "BUCSS")
  if (in_bucss_source && requireNamespace("pkgload", quietly = TRUE)) {
    suppressMessages(pkgload::load_all(".", quiet = TRUE))
  } else {
    library(BUCSS)
  }

  # Say what is being verified, so a run can never be mistaken for a check of
  # a different build.
  cat("BUCSS ", as.character(utils::packageVersion("BUCSS")), " from ",
      getNamespaceInfo("BUCSS", "path"), "\n", R.version.string, "\n", sep = "")

  cmd_args <- commandArgs(trailingOnly = TRUE)
  mode <- if ((length(cmd_args) >= 1L && cmd_args[1] == "full") ||
              identical(Sys.getenv("BUCSS_VERIFY"), "full")) "full" else "quick"
  reps_sens  <- if (mode == "full") 10000L else 1200L
  reps_data  <- if (mode == "full") 20000L else 3000L
  reps_welch <- if (mode == "full") 100000L else 20000L
  reps_lrt   <- if (mode == "full") 1000000L else 200000L
  cat("mode:", mode, "\n\n")

  withr::local_seed(2017, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion",
                    .rng_sample_kind = "Rejection")
  # Check results accumulate in a small environment, not by superassignment.
  state <- new.env(parent = emptyenv())
  state$results <- list()
  check <- function(name, ok, detail) {
    state$results[[length(state$results) + 1L]] <- list(name = name, ok = ok)
    cat(sprintf("  [%s] %s (%s)\n", if (ok) "PASS" else "FAIL", name, detail))
  }

  ## Check 1: assurance calibration on the statistic scale ---------------------

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
    list(label = "chi square difference",
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

  ## Check 2: correlation from raw bivariate data ------------------------------

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
  size_gap <- overstatement <- size_fixed <- rep(NA_real_, length(r_pub))
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
      size_fixed[j] <- n_f
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
  # The published claim: the overstatement is roughly 1/n for a plan of n.
  # Across 1,546 plans on a grid of r from .10 to .90 and N from 10 to 2,000,
  # overstatement times n lay between .86 and 1.05 (2026-09-28).
  scaled <- overstatement * size_fixed
  check("correlation: every fixed-predictor plan overstates its power, by about 1/n",
        all(overstatement > 0, na.rm = TRUE) &&
          all(scaled >= .8 & scaled <= 1.1, na.rm = TRUE),
        sprintf("overstated by %.4f to %.4f; times the plan's n, %.3f to %.3f (on a grid of 1,546 plans, .86 to 1.05)",
                min(overstatement, na.rm = TRUE),
                max(overstatement, na.rm = TRUE),
                min(scaled, na.rm = TRUE), max(scaled, na.rm = TRUE)))
  check("correlation: fixed-predictor plans run light",
        median(size_gap, na.rm = TRUE) >= 1,
        sprintf("median size gap %.1f (typically 2 to 5)",
                median(size_gap, na.rm = TRUE)))

  ## Check 3: Welch t against raw data -----------------------------------------

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

  # The floor: this prior would otherwise be planned at 6 per group, where the
  # real Welch test has power of about .79, short of the .80 target the
  # approximation reports. The planner returns at least 10 per group, and the
  # real test at the returned size must reach the target (it does not at 6).
  plan_s <- ss_buc_welch_t(t_observed = 4.5, n_1 = 12, n_2 = 12, assurance = .5)
  n_s <- plan_s$value[plan_s$term == "necessary_n_per_group"]
  mu_s <- plan_s$value[plan_s$term == "ncp_adjusted"] * sqrt(1 / 12 + 1 / 12)
  z <- welch_stat(n_s, n_s, 1, 1, mu_s, reps_welch)
  power_s <- mean(abs(z$t) > qt(.975, z$df))
  mcse_s <- sqrt(power_s * (1 - power_s) / reps_welch)
  check("Welch: a prior that would plan 6 per group is planned at 10, which reaches the target",
        n_s >= 10 && power_s >= desired_power - 3 * mcse_s,
        sprintf("planned %d per group; real Welch power %.3f (target %.2f)",
                as.integer(n_s), power_s, desired_power))

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

  ## Check 4: chi square difference test, finite-sample Type I -----------------

  cat("Check 4: finite-sample Type I error of the nested-model LRT\n")
  # For nested normal linear models (drop k = 2 of p = 6 parameters) the null
  # LRT is exactly N*log(1 + chisq_k / chisq_{N-p}), so no raw data is needed.
  lrt_alpha <- function(N, k = 2, p = 6, reps) {
    lr <- N * log(1 + rchisq(reps, k) / rchisq(reps, N - p))
    mean(lr > qchisq(.95, k))
  }
  # The same Type I error in closed form: the LRT rejects exactly when the
  # F statistic for the k constraints exceeds (N - p)/k * (exp(crit/N) - 1).
  lrt_alpha_exact <- function(N, k = 2, p = 6) {
    pf((N - p) / k * expm1(qchisq(.95, k) / N), k, N - p, lower.tail = FALSE)
  }
  for (N_lrt in c(50, 200, 400)) {
    a_sim <- lrt_alpha(N_lrt, reps = reps_lrt)
    a_exact <- lrt_alpha_exact(N_lrt)
    check(sprintf("chi square: Type I error at total N = %d", N_lrt),
          abs(a_sim - a_exact) <= 5 * sqrt(a_exact * (1 - a_exact) / reps_lrt),
          sprintf("simulated %.4f, exact %.4f (help page: about .07 at 50, .055 at 200, .052 at 400)",
                  a_sim, a_exact))
  }
  check("chi square: the figures the help page quotes (.07, .055, .052) are the exact values",
        abs(lrt_alpha_exact(50) - .07) < .005 &&
          abs(lrt_alpha_exact(200) - .055) < .0025 &&
          abs(lrt_alpha_exact(400) - .052) < .0025,
        sprintf("exact %.4f, %.4f, %.4f", lrt_alpha_exact(50),
                lrt_alpha_exact(200), lrt_alpha_exact(400)))

  ## Summary -------------------------------------------------------------------

  n_fail <- sum(!vapply(state$results, `[[`, TRUE, "ok"))
  cat(sprintf("\n%d checks, %d failed.\n", length(state$results), n_fail))
  if (n_fail > 0L)
    stop(n_fail, " of ", length(state$results), " checks failed.",
         call. = FALSE)
}))
