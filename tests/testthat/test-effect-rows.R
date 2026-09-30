# The planning effect size: every plan carries one row, directly after
# ncp_adjusted, restating that parameter in the design's usual metric (a
# standardized mean difference, Cohen's f, Cohen's f squared, a squared
# multiple correlation, a correlation, a squared Mahalanobis distance, or the
# root mean square error of approximation of a difference test). It is the
# effect size the plan is built on, so the ordinary power of a study of the
# returned size at that effect, computed below from the effect alone with a
# textbook power function and never from ncp_adjusted, is exactly
# actual_power. That one statement pins the metric's definition, the branch an
# unbalanced prior's row comes from (the reported parameter's), and the
# planned study's scaling. It holds in the frame the plan was made in: with
# random_predictors = TRUE the power is the exact random-predictor
# (negative binomial mixture) power, and with epsilon below 1 it is the
# Muller and Barton approximation, corrected degrees of freedom and
# noncentrality, with the effect on the sphericity scale.

val <- function(x, term) x$value[x$term == term]
size <- function(x) x$value[x$term %in% BUCSS:::.BUCSS_SIZE_TERMS]
effect_term <- function(x) x$term[x$term %in% BUCSS:::.BUCSS_EFFECT_TERMS]
effect_value <- function(x) x$value[x$term %in% BUCSS:::.BUCSS_EFFECT_TERMS]

# Two-sided t power, F power under the Muller and Barton approximation (the
# plain noncentral F at epsilon = 1), and the distribution of the F statistic
# for q coefficients when the predictors are sampled along with the outcome:
# a chi square (m df) mixture of noncentral F's, whose Poisson weights mixed
# over the gamma become negative binomial weights with size m/2 and
# probability 1/(1 + f^2). Written out here from the formula, summing the
# whole series, independently of the package's windowed implementation.
t_power <- function(delta, df, alpha) {
  crit <- qt(1 - alpha / 2, df)
  (1 - pt(crit, df, ncp = delta)) + pt(-crit, df, ncp = delta)
}
f_power <- function(lambda, df1, df2, alpha, epsilon = 1) {
  d1 <- epsilon * df1
  d2 <- epsilon * df2
  1 - pf(qf(1 - alpha, d1, d2), d1, d2, ncp = epsilon * lambda)
}
mixture_cdf <- function(x, df1, df2, m, prob) {
  k <- 0:qnbinom(1 - 1e-15, size = m / 2, prob = prob)
  sum(dnbinom(k, size = m / 2, prob = prob) *
        pbeta(df1 * x / (df1 * x + df2), df1 / 2 + k, df2 / 2))
}
random_power <- function(f2, n, p, q, alpha) {
  df2 <- n - p - 1
  1 - mixture_cdf(qf(1 - alpha, q, df2), q, df2, n - 1 - (p - q), 1 / (1 + f2))
}
welch_df <- function(a, b, ratio) {
  v1 <- 1 / a
  v2 <- ratio^2 / b
  (v1 + v2)^2 / (v1^2 / (a - 1) + v2^2 / (b - 1))
}

# The ordinary power of the returned size at the planning effect size, in
# the plan's own frame, from the effect alone.
power_from_effect <- function(r) {
  e <- effect_value(r)
  a <- val(r, "alpha_planned")
  n <- size(r)
  eff <- attr(r, "effect")
  eps <- if (length(val(r, "epsilon"))) val(r, "epsilon") else 1
  random <- length(val(r, "random_predictors")) && val(r, "random_predictors") == 1
  switch(attr(r, "design"),
    "Independent t test" = t_power(e * sqrt(n / 2), 2 * n - 2, a),
    "Dependent (paired) t test" = ,
    "One-sample t test" = t_power(e * sqrt(n), n - 1, a),
    "Welch (unequal variance) t test" = {
      ratio <- if (length(val(r, "sd_ratio"))) val(r, "sd_ratio")
               else if (length(val(r, "sd_1"))) val(r, "sd_2") / val(r, "sd_1")
               else sqrt(val(r, "var_2") / val(r, "var_1"))
      t_power(e * sqrt(n / 2), welch_df(n, n, ratio), a)
    },
    "One-way between-subjects ANOVA" = {
      g <- val(r, "levels_A")
      f_power(n * g * e^2, g - 1, n * g - g, a)
    },
    "Two-way between-subjects ANOVA" = {
      A <- val(r, "levels_A"); B <- val(r, "levels_B"); cells <- A * B
      df1 <- switch(eff, factor_A = A - 1, factor_B = B - 1,
                    interaction = (A - 1) * (B - 1))
      f_power(n * cells * e^2, df1, n * cells - cells, a)
    },
    "Between-subjects ANOVA (general)" = {
      cells <- val(r, "cells")
      nuisance <- val(r, "N") - cells - val(r, "df_denominator")
      f_power(n * cells * e^2, val(r, "df_numerator"), n * cells - cells - nuisance, a)
    },
    "Analysis of covariance" = {
      cells <- val(r, "cells")
      f_power(n * cells * e^2, val(r, "df_numerator"),
              n * cells - cells - val(r, "n_covariates"), a)
    },
    "One- or two-way within-subjects ANOVA" = {
      A <- val(r, "levels_A"); B <- val(r, "levels_B")
      df1 <- if (!length(B)) A - 1 else
        switch(eff, factor_A = A - 1, factor_B = B - 1,
               interaction = (A - 1) * (B - 1))
      # k = df1 error degrees of freedom per subject
      f_power(n * df1 * e^2, df1, df1 * (n - 1), a, eps)
    },
    "Within-subjects ANOVA (any number of factors)" = {
      df1 <- val(r, "df_numerator")
      f_power(n * df1 * e^2, df1, df1 * (n - 1), a, eps)
    },
    "Two-factor split-plot (mixed) ANOVA" = {
      g <- val(r, "levels_between"); b <- val(r, "levels_within")
      df1 <- switch(eff, between = g - 1, within = b - 1,
                    interaction = (g - 1) * (b - 1))
      k <- if (eff == "between") 1 else b - 1
      f_power(n * g * k * e^2, df1, (n * g - g) * k, a, eps)
    },
    "Split-plot (mixed) ANOVA (any number of factors)" = {
      g <- val(r, "num_groups"); df1 <- val(r, "df_numerator")
      k <- switch(eff, between_only = 1, within_only = df1,
                  between_within = val(r, "df_num_within"))
      f_power(n * g * k * e^2, df1, (n * g - g) * k, a, eps)
    },
    "Multiple regression: single coefficient" = {
      p <- val(r, "p")
      if (random) random_power(e, n, p, 1, a) else f_power(n * e, 1, n - p - 1, a)
    },
    "Multiple regression: joint test of predictors" = {
      p <- val(r, "p"); q <- val(r, "p_joint")
      if (random) random_power(e, n, p, q, a) else f_power(n * e, q, n - p - 1, a)
    },
    "Multiple regression: omnibus test of model R2" = {
      p <- val(r, "p")
      f2 <- e / (1 - e)
      if (random) random_power(f2, n, p, p, a) else f_power(n * f2, p, n - p - 1, a)
    },
    "Pearson correlation" = random_power(e^2 / (1 - e^2), n, 1, 1, a),
    "Two-group multivariate comparison (Hotelling's T squared)" = {
      pv <- val(r, "p_variables")
      f_power(n / 2 * e, pv, 2 * n - pv - 1, a)
    },
    "Nested model chi square difference test" = {
      d <- val(r, "df_restricted") - val(r, "df_full")
      1 - pchisq(qchisq(1 - a, d), d, ncp = d * n * e^2)
    },
    stop("no power formula for ", attr(r, "design")))
}

# One or more plans per design, echoed form, and frame, with unbalanced prior
# studies, sphericity corrections, sampled predictors, negative statistics,
# and other alpha levels, assurances, and target powers.
effect_plans <- list(
  ss_buc_independent_t(t_observed = 3, n = 20),
  ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .9),
  ss_buc_independent_t(t_observed = -2.8, N = 41, alpha_prior = .1),
  ss_buc_independent_t(t_observed = 4, n = c(14, 21), desired_power = .9,
                       alpha_planned = .01),
  ss_buc_smd(t_observed = 3, n = 20, alpha_prior = 1, assurance = .5),
  ss_buc_paired_t(t_observed = 3, N = 40),
  ss_buc_smd_paired(t_observed = -3.5, N = 25, assurance = .7),
  ss_buc_one_sample_t(t_observed = 3, N = 40),
  ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5),
  ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 30, sd_1 = 1, sd_2 = 2),
  ss_buc_welch_t(t_observed = -3.2, n_1 = 30, n_2 = 40, var_1 = 1, var_2 = 3),
  ss_buc_welch_t(t_observed = 6, n_1 = 20, n_2 = 20, sd_ratio = 1),  # at the floor
  ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 4),
  ss_buc_one_way_anova(F_observed = 5, N = 121, levels_A = 3),
  ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                         effect = "factor_B"),
  ss_buc_factorial_anova(F_observed = 6, N = 125, levels_A = 2, levels_B = 3,
                         effect = "interaction", assurance = .6),
  ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6, df_numerator = 2,
                                 df_denominator = 104),
  ss_buc_factorial_anova_general(F_observed = 6, N = 125, cells = 6, df_numerator = 2,
                                 df_denominator = 119),
  ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2),
  ss_buc_ancova(F_observed = 9, N = 121, cells = 4, n_covariates = 3, df_numerator = 1),
  ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3),
  ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3, effect = "factor_B"),
  ss_buc_rm_anova(F_observed = 6, N = 40, levels_A = 2, levels_B = 3,
                  effect = "interaction", epsilon = .7),
  ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = 1.73 / 2),
  ss_buc_rm_anova_general(F_observed = 6.5, N = 80, df_numerator = 1, assurance = .5),
  ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2, epsilon = .7),
  ss_buc_mixed_anova(F_observed = 9, N = 60, levels_between = 2, levels_within = 3,
                     effect = "between"),
  ss_buc_mixed_anova(F_observed = 5, N = 61, levels_between = 2, levels_within = 3,
                     effect = "within"),
  ss_buc_mixed_anova(F_observed = 8, N = 91, levels_between = 2, levels_within = 3,
                     effect = "interaction", epsilon = .8),
  ss_buc_mixed_anova_general(F_observed = 5, N = 91, df_numerator = 2, num_groups = 3,
                             effect = "between_only"),
  ss_buc_mixed_anova_general(F_observed = 8, N = 90, df_numerator = 2, num_groups = 3,
                             effect = "within_only", epsilon = .75),
  ss_buc_mixed_anova_general(F_observed = 4, N = 61, df_numerator = 6, num_groups = 3,
                             effect = "between_within", df_num_within = 3,
                             epsilon = .6),
  ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3),
  ss_buc_reg_coef(t_observed = -3, N = 150, p = 3, random_predictors = TRUE),
  ss_buc_R2(F_observed = 8, N = 60, p = 3),
  ss_buc_R2(F_observed = 5, N = 150, p = 4, random_predictors = TRUE),
  ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2),
  ss_buc_reg_joint(F_observed = 6, N = 60, p = 6, p_joint = 3, random_predictors = TRUE),
  ss_buc_correlation(r_observed = .35, N = 100),
  ss_buc_correlation(r_observed = -.3, N = 151, alpha_prior = 1, assurance = .9),
  ss_buc_correlation(t_observed = 3.7, N = 60),
  ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3),
  ss_buc_manova(T2_observed = 16, N = 61, p_variables = 3),
  ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42, df_restricted = 43),
  ss_buc_chisq_diff(chisq_observed = 15, N = 200, df_full = 42, df_restricted = 44,
                    desired_power = .9)
)

test_that("every plan carries exactly one planning effect size row, its design's own, after ncp_adjusted", {
  covered <- unique(vapply(effect_plans, attr, character(1), which = "design"))
  expect_setequal(covered, names(BUCSS:::.BUCSS_DESIGN_SPECS))
  specs <- BUCSS:::.BUCSS_DESIGN_SPECS
  for (nm in names(specs))
    expect_true(specs[[nm]]$planning_term %in% BUCSS:::.BUCSS_EFFECT_TERMS, info = nm)
  expect_setequal(names(BUCSS:::.BUCSS_EFFECT_LABELS), BUCSS:::.BUCSS_EFFECT_TERMS)
  expect_true(all(BUCSS:::.BUCSS_EFFECT_TERMS %in% BUCSS:::.BUCSS_COMPUTED_TERMS))
  for (r in effect_plans) {
    design <- attr(r, "design")
    expect_identical(sum(r$term %in% BUCSS:::.BUCSS_EFFECT_TERMS), 1L, info = design)
    expect_identical(effect_term(r), specs[[design]]$planning_term, info = design)
    expect_identical(which(r$term %in% BUCSS:::.BUCSS_EFFECT_TERMS),
                     which(r$term == "ncp_adjusted") + 1L, info = design)
    expect_true(is.finite(effect_value(r)) && effect_value(r) > 0, info = design)
    expect_true(BUCSS:::.is_plan_shape(r), info = design)
  }
})

test_that("actual_power is the ordinary power at the planning effect size, in the plan's own frame", {
  for (r in effect_plans) {
    info <- paste(attr(r, "design"), attr(r, "effect"),
                  paste(r$term, signif(r$value, 6), collapse = " "))
    expect_equal(power_from_effect(r), val(r, "actual_power"), tolerance = 1e-12,
                 info = info)
  }
})

test_that("the identity is frame-aware: the other frame's power is not actual_power", {
  # With sampled predictors the fixed-frame noncentral F overstates the power
  # at the same f squared; under a sphericity correction the uncorrected test
  # overstates it. The planning effect size keeps its meaning; the power
  # statement is the frame's.
  rand <- ss_buc_R2(F_observed = 5, N = 150, p = 4, random_predictors = TRUE)
  f2 <- effect_value(rand) / (1 - effect_value(rand))
  n <- size(rand)
  expect_gt(f_power(n * f2, 4, n - 5, .05), val(rand, "actual_power") + .003)
  eps <- ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2, epsilon = .7)
  n <- size(eps)
  expect_gt(f_power(n * 2 * effect_value(eps)^2, 2, 2 * (n - 1), .05),
            val(eps, "actual_power") + .003)
})

test_that("every planner's row is the documented conversion of ncp_adjusted at the balanced prior design", {
  ncp <- function(x) val(x, "ncp_adjusted")
  r <- ss_buc_independent_t(t_observed = 3, n = 20)
  expect_equal(val(r, "smd_planning"), ncp(r) * sqrt(2 / 20), tolerance = 1e-14)
  r <- ss_buc_paired_t(t_observed = 3, N = 40)
  expect_equal(val(r, "dz_planning"), ncp(r) / sqrt(40), tolerance = 1e-14)
  r <- ss_buc_one_sample_t(t_observed = 3, N = 40)
  expect_equal(val(r, "std_mean_planning"), ncp(r) / sqrt(40), tolerance = 1e-14)
  r <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 30, sd_ratio = 1.5)
  expect_equal(val(r, "smd_planning"),
               ncp(r) * sqrt(1 / 40 + 1.5^2 / 30) / sqrt((1 + 1.5^2) / 2), tolerance = 1e-14)
  r <- ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 4)
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / 120), tolerance = 1e-14)
  # f^2/(1 + f^2) is the population partial eta squared, lambda/(lambda + N)
  expect_equal(val(r, "f_planning")^2 / (1 + val(r, "f_planning")^2),
               ncp(r) / (ncp(r) + 120), tolerance = 1e-14)
  r <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                              effect = "interaction")
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / 120), tolerance = 1e-14)
  r <- ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                      df_numerator = 2, df_denominator = 104)
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / 120), tolerance = 1e-14)
  r <- ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2)
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / 120), tolerance = 1e-14)
  r <- ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3)
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / (2 * 30)), tolerance = 1e-14)
  # on the sphericity scale for every epsilon: the same k and N
  r <- ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2, epsilon = .7)
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / (2 * 30)), tolerance = 1e-14)
  r <- ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2, levels_within = 3,
                          effect = "within")
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / (2 * 60)), tolerance = 1e-14)
  r <- ss_buc_mixed_anova(F_observed = 9, N = 60, levels_between = 2, levels_within = 3,
                          effect = "between")
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / 60), tolerance = 1e-14)
  r <- ss_buc_mixed_anova(F_observed = 8, N = 60, levels_between = 2, levels_within = 3,
                          effect = "interaction", epsilon = .8)
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / (2 * 60)), tolerance = 1e-14)
  r <- ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 6, num_groups = 3,
                                  effect = "between_within", df_num_within = 3)
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / (3 * 90)), tolerance = 1e-14)
  r <- ss_buc_mixed_anova_general(F_observed = 8, N = 90, df_numerator = 2, num_groups = 3,
                                  effect = "within_only")
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / (2 * 90)), tolerance = 1e-14)
  r <- ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2, num_groups = 3,
                                  effect = "between_only")
  expect_equal(val(r, "f_planning"), sqrt(ncp(r) / 90), tolerance = 1e-14)
  r <- ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3)
  expect_equal(val(r, "f2_planning"), ncp(r) / 60, tolerance = 1e-14)
  r <- ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                        random_predictors = TRUE)
  expect_equal(val(r, "f2_planning"), ncp(r) / 150, tolerance = 1e-14)
  r <- ss_buc_R2(F_observed = 8, N = 60, p = 3)
  expect_equal(val(r, "R2_planning"), ncp(r) / (ncp(r) + 60), tolerance = 1e-14)
  r <- ss_buc_correlation(r_observed = .35, N = 100)
  expect_equal(val(r, "rho_planning")^2, ncp(r) / (ncp(r) + 100), tolerance = 1e-14)
  r <- ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3)
  expect_equal(val(r, "D2_planning"), 2 * ncp(r) / 40, tolerance = 1e-14)
  r <- ss_buc_chisq_diff(chisq_observed = 15, N = 200, df_full = 42, df_restricted = 44)
  expect_equal(val(r, "rmsea_d_planning"), sqrt(ncp(r) / (2 * 200)), tolerance = 1e-14)
})

test_that("an unbalanced prior's row is that of the balanced design whose parameter is reported", {
  # Option A: the row restates ncp_adjusted, so it comes from the same
  # rounding branch. Each branch is itself the plan of a balanced prior study.
  check <- function(unbalanced, round_down, round_up, term) {
    use_rd <- val(round_down, "ncp_adjusted") <= val(round_up, "ncp_adjusted")
    branch <- if (use_rd) round_down else round_up
    expect_identical(val(unbalanced, "ncp_adjusted"), val(branch, "ncp_adjusted"))
    expect_identical(val(unbalanced, term), val(branch, term))
    # the returned size is the larger of the two branches' requirements, so an
    # ordinary power analysis at the row can ask for fewer, never more
    expect_identical(size(unbalanced), max(size(round_down), size(round_up)))
    expect_lte(size(branch), size(unbalanced))
  }
  # the harmonic mean of 50 and 55 is 52.38: branches of 52 and 53 per group
  check(ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .9),
        ss_buc_independent_t(t_observed = 3, n = 52, assurance = .9),
        ss_buc_independent_t(t_observed = 3, n = 53, assurance = .9), "smd_planning")
  check(ss_buc_one_way_anova(F_observed = 5, N = 121, levels_A = 3),
        ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 3),
        ss_buc_one_way_anova(F_observed = 5, N = 123, levels_A = 3), "f_planning")
  check(ss_buc_mixed_anova(F_observed = 5, N = 61, levels_between = 2, levels_within = 3,
                           effect = "within"),
        ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2, levels_within = 3,
                           effect = "within"),
        ss_buc_mixed_anova(F_observed = 5, N = 62, levels_between = 2, levels_within = 3,
                           effect = "within"), "f_planning")
  check(ss_buc_manova(F_observed = 4.5, N = 81, p_variables = 3),
        ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3),
        ss_buc_manova(F_observed = 4.5, N = 82, p_variables = 3), "D2_planning")

  # For the flagship the reported parameter is the round-down branch's (52
  # per group), and an ordinary power analysis at the row asks for that
  # branch's own 1,475 per group; the plan is the round-up branch's 1,485.
  flagship <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .9)
  d <- val(flagship, "smd_planning")
  ordinary <- function(n) t_power(d * sqrt(n / 2), 2 * n - 2, .05)
  n_ordinary <- BUCSS:::.smallest_n_for_power(ordinary, .8)
  expect_identical(n_ordinary, 1475)
  expect_identical(n_ordinary, size(ss_buc_independent_t(t_observed = 3, n = 52,
                                                         assurance = .9)))
  expect_identical(size(flagship), 1485)
})

test_that("identities among the planners hold for the planning effect size", {
  # One outcome variable: Hotelling's T squared is the squared two-group t,
  # and the squared Mahalanobis distance is the squared standardized mean
  # difference (the two roots are solved on different scales, so equal to the
  # root finder's tolerance).
  for (N in c(40, 41, 100)) {
    uni <- ss_buc_manova(F_observed = 9, N = N, p_variables = 1)
    tt <- ss_buc_independent_t(t_observed = 3, N = N)
    expect_equal(val(uni, "D2_planning"), val(tt, "smd_planning")^2, tolerance = 1e-6,
                 info = N)
    expect_identical(size(uni), size(tt))
  }
  # An ANCOVA with no covariates is the one-way ANOVA.
  for (N in c(60, 61)) {
    expect_identical(
      val(ss_buc_ancova(F_observed = 5, N = N, cells = 3, n_covariates = 0), "f_planning"),
      val(ss_buc_one_way_anova(F_observed = 5, N = N, levels_A = 3), "f_planning"))
  }
  # The general between-subjects planner with no nuisance parameters is the
  # two-way planner.
  expect_identical(
    val(ss_buc_factorial_anova_general(F_observed = 5, N = 125, cells = 6, df_numerator = 2,
                                       df_denominator = 119), "f_planning"),
    val(ss_buc_factorial_anova(F_observed = 5, N = 125, levels_A = 2, levels_B = 3,
                               effect = "factor_B"), "f_planning"))
  # The two within-subjects planners agree, and so do the two split-plot
  # planners for all three effects.
  expect_identical(
    val(ss_buc_rm_anova(F_observed = 7, N = 60, levels_A = 2, levels_B = 3,
                        effect = "interaction", epsilon = .8), "f_planning"),
    val(ss_buc_rm_anova_general(F_observed = 7, N = 60, df_numerator = 2, epsilon = .8),
        "f_planning"))
  expect_identical(
    val(ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3), "f_planning"),
    val(ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2), "f_planning"))
  pairs <- list(
    list(ss_buc_mixed_anova(F_observed = 9, N = 91, levels_between = 3, levels_within = 4,
                            effect = "between"),
         ss_buc_mixed_anova_general(F_observed = 9, N = 91, df_numerator = 2, num_groups = 3,
                                    effect = "between_only")),
    list(ss_buc_mixed_anova(F_observed = 5, N = 91, levels_between = 3, levels_within = 4,
                            effect = "within", epsilon = .7),
         ss_buc_mixed_anova_general(F_observed = 5, N = 91, df_numerator = 3, num_groups = 3,
                                    effect = "within_only", epsilon = .7)),
    list(ss_buc_mixed_anova(F_observed = 4, N = 91, levels_between = 3, levels_within = 4,
                            effect = "interaction"),
         ss_buc_mixed_anova_general(F_observed = 4, N = 91, df_numerator = 6, num_groups = 3,
                                    effect = "between_within", df_num_within = 3)))
  for (pr in pairs) expect_identical(val(pr[[1]], "f_planning"), val(pr[[2]], "f_planning"))
  # The one-sample planner is the paired planner under another name.
  one <- ss_buc_one_sample_t(t_observed = 3, N = 40)
  paired <- ss_buc_paired_t(t_observed = 3, N = 40)
  expect_identical(val(one, "std_mean_planning"), val(paired, "dz_planning"))
  expect_false("dz_planning" %in% one$term)
  # A coefficient and the omnibus R squared with one predictor share the
  # parameter, so R squared is f squared / (1 + f squared), in either frame.
  for (random in c(FALSE, TRUE)) {
    coef <- ss_buc_reg_coef(t_observed = 3, N = 60, p = 1, random_predictors = random)
    omni <- ss_buc_R2(F_observed = 9, N = 60, p = 1, random_predictors = random)
    f2 <- val(coef, "f2_planning")
    expect_equal(val(omni, "R2_planning"), f2 / (1 + f2), tolerance = 1e-10, info = random)
  }
  # With one sampled predictor the coefficient planner is the correlation
  # planner, and its f squared is rho^2/(1 - rho^2) of the correlation's rho.
  coef <- ss_buc_reg_coef(t_observed = 3, N = 60, p = 1, random_predictors = TRUE)
  corr <- ss_buc_correlation(t_observed = 3, N = 60)
  rho <- val(corr, "rho_planning")
  expect_equal(val(coef, "f2_planning"), rho^2 / (1 - rho^2), tolerance = 1e-12)
  # The Welch planner with a ratio of 1 and equal prior groups is the
  # equal-variance planner, above its floor of 10 per group.
  welch <- ss_buc_welch_t(t_observed = 3, n_1 = 20, n_2 = 20, sd_ratio = 1)
  indep <- ss_buc_independent_t(t_observed = 3, n = 20)
  expect_equal(val(welch, "smd_planning"), val(indep, "smd_planning"), tolerance = 1e-12)
  expect_identical(size(welch), size(indep))
})

test_that("the t-based and correlation rows are magnitudes", {
  flips <- list(
    list(ss_buc_independent_t(t_observed = 3, n = c(50, 55)),
         ss_buc_independent_t(t_observed = -3, n = c(50, 55))),
    list(ss_buc_paired_t(t_observed = 3, N = 40), ss_buc_paired_t(t_observed = -3, N = 40)),
    list(ss_buc_one_sample_t(t_observed = 3, N = 40),
         ss_buc_one_sample_t(t_observed = -3, N = 40)),
    list(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 30, sd_ratio = 1.5),
         ss_buc_welch_t(t_observed = -3, n_1 = 40, n_2 = 30, sd_ratio = 1.5)),
    list(ss_buc_reg_coef(t_observed = 3, N = 150, p = 3),
         ss_buc_reg_coef(t_observed = -3, N = 150, p = 3)),
    list(ss_buc_correlation(r_observed = .35, N = 100),
         ss_buc_correlation(r_observed = -.35, N = 100)))
  for (pr in flips) {
    expect_identical(effect_value(pr[[2]]), effect_value(pr[[1]]))
    expect_gt(effect_value(pr[[2]]), 0)
  }
})

test_that("the documented within-subjects hand-off reproduces the plan", {
  # Power functions that standardize by the raw-score standard deviation and
  # take the correlation among the measures separately (lambda = N m f^2 /
  # (1 - rho) for the within-subjects tests and N m f^2 / (1 + (m - 1) rho)
  # for the between-subjects test of a split-plot design, with m measurements
  # per subject) reproduce the plan, at a correlation of 0 and the same
  # epsilon, from f_planning times sqrt(k/m), where k is the error degrees of
  # freedom per subject: (a - 1)/a for a one-way within-subjects design.
  rm <- ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = 1.73 / 2)
  a <- 3
  e <- 1.73 / 2
  n <- size(rm)
  f_raw <- val(rm, "f_planning") * sqrt((a - 1) / a)
  expect_equal(f_power(n * a * f_raw^2, a - 1, (a - 1) * (n - 1), .05, e),
               val(rm, "actual_power"), tolerance = 1e-12)
  for (effect in c("between", "within", "interaction")) {
    b <- 3
    g <- 2
    sp <- ss_buc_mixed_anova(F_observed = 8, N = 60, levels_between = g, levels_within = b,
                             effect = effect, epsilon = if (effect == "between") 1 else .8)
    e <- if (effect == "between") 1 else .8
    k <- if (effect == "between") 1 else b - 1
    df1 <- switch(effect, between = g - 1, within = b - 1, interaction = (g - 1) * (b - 1))
    n_total <- size(sp) * g
    f_raw <- val(sp, "f_planning") * sqrt(k / b)
    expect_equal(f_power(n_total * b * f_raw^2, df1, (n_total - g) * k, .05, e),
                 val(sp, "actual_power"), tolerance = 1e-12, info = effect)
  }
})

test_that("the help-page examples' planning effect sizes are pinned", {
  pin <- function(x, term, value) {
    expect_identical(effect_term(x), term)
    expect_equal(effect_value(x), value, tolerance = 1e-9)
  }
  pin(ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = .05, alpha_planned = .05,
                           assurance = .80, desired_power = .80),
      "smd_planning", 0.349438801615)
  pin(ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .90),
      "smd_planning", 0.103204674334)
  pin(ss_buc_paired_t(t_observed = 3, N = 40), "dz_planning", 0.176343409855)
  pin(ss_buc_one_sample_t(t_observed = 3, N = 40), "std_mean_planning", 0.176343409855)
  pin(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5),
      "smd_planning", 0.266784209225)
  pin(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1),
      "smd_planning", 0.272381681217)
  pin(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 2),
      "smd_planning", 0.25901053189)
  pin(ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 4),
      "f_planning", 0.176461168656)
  pin(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                             effect = "factor_B"),
      "f_planning", 0.0494152017458)
  pin(ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6, df_numerator = 2,
                                     df_denominator = 114),
      "f_planning", 0.0494152017458)
  pin(ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2),
      "f_planning", 0.145221463116)
  pin(ss_buc_ancova(F_observed = 9, N = 120, cells = 3, n_covariates = 2, df_numerator = 1),
      "f_planning", 0.123452947972)
  pin(ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3, effect = "factor_B"),
      "f_planning", 0.0503531394155)
  pin(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = 1.73 / 2),
      "f_planning", 0.1826927721)
  pin(ss_buc_rm_anova_general(F_observed = 6.5, N = 80, df_numerator = 1, assurance = .50),
      "f_planning", 0.175849813952)
  pin(ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2, epsilon = .7),
      "f_planning", 0.240593933138)
  pin(ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2, levels_within = 3,
                         effect = "within"),
      "f_planning", 0.0498935525155)
  pin(ss_buc_mixed_anova(F_observed = 7, N = 60, levels_between = 2, levels_within = 3,
                         effect = "within", epsilon = .8),
      "f_planning", 0.1626423004)
  pin(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2, num_groups = 3,
                                 effect = "between_only"),
      "f_planning", 0.0464675947857)
  pin(ss_buc_mixed_anova_general(F_observed = 8, N = 90, df_numerator = 2, num_groups = 3,
                                 effect = "within_only", epsilon = .75),
      "f_planning", 0.170141671047)
  pin(ss_buc_reg_coef(t_observed = 3, N = 150, p = 3), "f2_planning", 0.0126209268614)
  pin(ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, random_predictors = TRUE),
      "f2_planning", 0.012741986309)
  pin(ss_buc_R2(F_observed = 5, N = 150, p = 4), "R2_planning", 0.0495284462904)
  pin(ss_buc_R2(F_observed = 5, N = 150, p = 4, random_predictors = TRUE),
      "R2_planning", 0.0491587351307)
  pin(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2),
      "f2_planning", 0.00243520057088)
  pin(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2, random_predictors = TRUE),
      "f2_planning", 0.00247675238128)
  pin(ss_buc_correlation(r_observed = .35, N = 100), "rho_planning", 0.250821974826)
  pin(ss_buc_correlation(t_observed = .35 * sqrt(98 / (1 - .35^2)), N = 100),
      "rho_planning", 0.250821974826)
  pin(ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3), "D2_planning", 0.0708792315057)
  pin(ss_buc_manova(T2_observed = 4.5 * 3 * 78 / 76, N = 80, p_variables = 3),
      "D2_planning", 0.0708792315057)
  pin(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42, df_restricted = 43),
      "rmsea_d_planning", 0.104019240822)
})

test_that("the row is a computed row: required by the shape check, compared and never replayed", {
  shape <- "no longer has the shape of a single plan"
  keep_attrs <- function(y, x) {
    for (nm in setdiff(names(attributes(x)), c("names", "row.names")))
      attr(y, nm) <- attr(x, nm)
    y
  }
  plan <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .9)
  # removed
  trimmed <- keep_attrs(as.data.frame(plan)[plan$term != "smd_planning", ], plan)
  expect_false(BUCSS:::.is_plan_shape(trimmed))
  expect_error(planning_sentence(trimmed), shape, fixed = TRUE)
  expect_error(glance(trimmed), shape, fixed = TRUE)
  # renamed to another design's metric
  renamed <- plan
  renamed$term[renamed$term == "smd_planning"] <- "dz_planning"
  renamed <- keep_attrs(as.data.frame(renamed), plan)
  expect_false(BUCSS:::.is_plan_shape(renamed))
  # doubled
  doubled <- keep_attrs(rbind(as.data.frame(plan),
                              data.frame(term = "f_planning", value = .1)), plan)
  expect_false(BUCSS:::.is_plan_shape(doubled))
  # an object from no registry design still needs exactly one such row
  bare <- structure(as.data.frame(plan), class = c("bucss_plan", "data.frame"))
  expect_true(BUCSS:::.is_plan_shape(bare))
  expect_false(BUCSS:::.is_plan_shape(
    structure(as.data.frame(plan)[plan$term != "smd_planning", ],
              class = c("bucss_plan", "data.frame"))))

  for (r in effect_plans) {
    inputs <- BUCSS:::.bucss_design_inputs(r)
    expect_false(any(names(inputs) %in% BUCSS:::.BUCSS_EFFECT_TERMS), info = attr(r, "design"))
    spec <- BUCSS:::.BUCSS_DESIGN_SPECS[[attr(r, "design")]]
    replay <- do.call(get(spec$planner, envir = asNamespace("BUCSS")),
                      BUCSS:::.bucss_design_args(inputs, spec, attr(r, "effect")))
    expect_true(BUCSS:::.replay_matches(replay, r), info = attr(r, "design"))
  }
  # a changed row no longer replays
  tampered <- plan
  tampered$value[tampered$term == "smd_planning"] <- .2
  spec <- BUCSS:::.BUCSS_DESIGN_SPECS[[attr(plan, "design")]]
  replay <- do.call(ss_buc_independent_t,
                    BUCSS:::.bucss_design_args(BUCSS:::.bucss_design_inputs(plan), spec, NULL))
  expect_true(BUCSS:::.replay_matches(replay, plan))
  expect_false(BUCSS:::.replay_matches(replay, tampered))
})

test_that("the output views show the row, and the legacy positions are unchanged", {
  plan <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .9)
  expect_output(print(plan), "smd_planning +0.103")
  expect_identical(glance(plan)$smd_planning, val(plan, "smd_planning"))
  expect_identical(names(tidy(plan)), c("term", "estimate", "power"))
  # glance() keeps every row as a column, the planning effect size among them
  expect_identical(which(names(glance(plan)) == "smd_planning"),
                   which(names(glance(plan)) == "ncp_adjusted") + 1L)
  # [[1]] and [[2]] still read the sample size and ncp_adjusted
  expect_identical(plan[[1]], 1485)
  expect_identical(plan[[2]], val(plan, "ncp_adjusted"))
  # the deprecated 1.x wrappers return the same object, row included
  suppressWarnings(old <- ss.power.it(t.observed = 3, n = c(50, 55), assurance = .9))
  expect_identical(val(old, "smd_planning"), val(plan, "smd_planning"))
  skip_if_not_installed("knitr")
  expect_match(as.character(knitr::knit_print(plan)), "smd_planning", fixed = TRUE)
})

test_that("planning_sentence() names the planning effect size in a second sentence", {
  flagship <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .9)
  expect_identical(
    planning_sentence(flagship),
    paste0("A planned study with 1,485 participants per group (total N = 2,970) ",
           "attains 80% power with 90% assurance, after correcting the prior ",
           "study's result for publication bias (assuming its literature ",
           "published results significant at the .05 level) and for uncertainty. ",
           "The plan is built on a planning effect size of 0.103, expressed as a ",
           "standardized mean difference."))
  labels <- c(
    "Dependent (paired) t test" = "a standardized mean of the paired differences",
    "One-sample t test" = "a standardized mean",
    "One-way between-subjects ANOVA" = "Cohen's f",
    "Multiple regression: single coefficient" = "Cohen's f squared",
    "Multiple regression: omnibus test of model R2" = "a squared multiple correlation",
    "Pearson correlation" = "a correlation",
    "Two-group multivariate comparison (Hotelling's T squared)" =
      "a squared Mahalanobis distance",
    "Nested model chi square difference test" =
      "the root mean square error of approximation of the difference test")
  for (r in effect_plans) {
    s <- planning_sentence(r)
    design <- attr(r, "design")
    expect_match(s, "The plan is built on a planning effect size of ", fixed = TRUE)
    if (design %in% names(labels))
      expect_match(s, paste0(", expressed as ", labels[[design]], "( on|\\.$)"), info = design)
    # the value is never called corrected, and never printed in scientific
    # notation
    expect_no_match(s, "orrected in this way", fixed = TRUE)
    expect_no_match(sub(".*The plan is built on", "", s), "correct", fixed = TRUE)
    expect_no_match(s, "e-", fixed = TRUE)
  }
  # A metric that cannot exceed 1 is written without the leading zero.
  expect_match(planning_sentence(ss_buc_R2(F_observed = 8, N = 60, p = 3)),
               "planning effect size of .175, expressed as a squared multiple correlation.",
               fixed = TRUE)
  expect_match(planning_sentence(ss_buc_correlation(r_observed = .35, N = 100)),
               "planning effect size of .251, expressed as a correlation.", fixed = TRUE)
  # Under a sphericity correction the row is on the sphericity scale.
  expect_match(planning_sentence(ss_buc_rm_anova_general(F_observed = 8, N = 30,
                                                         df_numerator = 2, epsilon = .7)),
               "0.241, expressed as Cohen's f on the sphericity scale.", fixed = TRUE)
  expect_match(planning_sentence(ss_buc_rm_anova_general(F_observed = 8, N = 30,
                                                         df_numerator = 2)),
               ", expressed as Cohen's f.", fixed = TRUE)
  # Its digits follow the display option; the stored row is untouched.
  old <- options(bucss.digits = 5)
  on.exit(options(old), add = TRUE)
  expect_match(planning_sentence(flagship), "planning effect size of 0.1032, expressed",
               fixed = TRUE)
  expect_identical(val(flagship, "smd_planning"),
                   val(ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .9),
                       "smd_planning"))
})
