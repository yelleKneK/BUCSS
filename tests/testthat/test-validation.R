# Validation tests: every guard that should stop() or warning() before any
# sample size planning happens. Messages are matched as fixed substrings.

test_that("out-of-range alpha_prior is rejected", {
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = 1.5),
               "alpha_prior", fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = 0),
               "alpha_prior", fixed = TRUE)
})

test_that("out-of-range alpha_planned is rejected", {
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, alpha_planned = 1),
               "alpha_planned", fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, alpha_planned = 0),
               "alpha_planned", fixed = TRUE)
})

test_that("out-of-range assurance is rejected", {
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, assurance = -.1),
               "assurance", fixed = TRUE)
})

test_that("assurance below .5 warns but still computes", {
  expect_warning(ss_buc_independent_t(t_observed = 3, n = 20, assurance = .4),
                 "< .5", fixed = TRUE)
})

test_that("out-of-range power is rejected", {
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, power = -.1),
               "power", fixed = TRUE)
})

test_that("a nonsignificant prior t is rejected", {
  expect_error(ss_buc_independent_t(t_observed = 1, n = 20),
               "nonsignificant", fixed = TRUE)
})

test_that("a nonsignificant prior F is rejected", {
  expect_error(
    ss_buc_factorial_anova(F_observed = 1, N = 120, levels_A = 2, levels_B = 3,
                           effect = "factor_B"),
    "nonsignificant", fixed = TRUE)
})

test_that("between-subjects ANOVA requires N", {
  expect_error(
    ss_buc_factorial_anova(F_observed = 5, levels_A = 2, levels_B = 3,
                           effect = "factor_B"),
    "must specify 'N'", fixed = TRUE)
})

test_that("the two-way function requires a second factor and points to the one-way function", {
  expect_error(
    ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
                           effect = "interaction"),
    "levels_B", fixed = TRUE)
  expect_error(
    ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
                           effect = "interaction"),
    "ss_buc_one_way_anova", fixed = TRUE)
})

test_that("the general between-subjects function requires df_denominator", {
  expect_error(
    ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                   df_numerator = 2),
    "must specify 'df_denominator'", fixed = TRUE)
})

test_that("df_denominator cannot exceed N - cells", {
  # The 1.2.1 documented example used df_denominator = 117 with N = 120 and
  # cells = 6, which is impossible (the maximum residual df is 120 - 6 = 114).
  # Now that df_denominator functions, that value is correctly rejected.
  expect_error(
    ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                   df_numerator = 2, df_denominator = 117),
    "cannot exceed", fixed = TRUE)
})

test_that("df_denominator must be positive", {
  expect_error(
    ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                   df_numerator = 2, df_denominator = 0),
    "must be a positive number", fixed = TRUE)
})

test_that("a corrected noncentrality parameter of zero is rejected", {
  # t = 2.03 with n = 20 (df = 38) just clears the alpha = .05 critical value
  # (2.0244), so the prior effect is significant but so weak that at assurance
  # .95 the bias/uncertainty correction drives the ncp to zero.
  expect_error(ss_buc_independent_t(t_observed = 2.03, n = 20, assurance = .95),
               "noncentrality parameter is zero", fixed = TRUE)
})

test_that("the zero-ncp error reports the closed-form assurance ceiling (independent t)", {
  # The largest assurance a prior result can support is 1 - p/alpha_prior, where
  # p is the prior study's two-sided p value. For t = 3 with n = 20 (df = 38)
  # and the default alpha_prior = .05 that ceiling is about .90, so the
  # requested assurance of .95 is out of reach and the corrected ncp is zero.
  p_prior <- 2 * (1 - pt(3, df = 38))
  ceiling_assurance <- floor((1 - p_prior / .05) * 100 + 1e-9) / 100
  fragment <- paste0("usable plan is about ",
                     sub("^0", "", sprintf("%.2f", ceiling_assurance)))
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, assurance = .95),
               fragment, fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, assurance = .95),
               "noncentrality parameter is zero", fixed = TRUE)
})

test_that("the zero-ncp error reports the closed-form assurance ceiling (joint F test)", {
  # The same closed form holds for an F test, where p is the upper-tail F p
  # value. F = 5 on 2 and 145 degrees of freedom gives a ceiling near .84, below
  # the requested assurance of .95, confirming the ceiling on the single-branch
  # F path as well as the two-branch t path above.
  p_prior <- 1 - pf(5, df1 = 2, df2 = 145)
  ceiling_assurance <- floor((1 - p_prior / .05) * 100 + 1e-9) / 100
  fragment <- paste0("usable plan is about ",
                     sub("^0", "", sprintf("%.2f", ceiling_assurance)))
  expect_error(
    ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2, assurance = .95),
    fragment, fixed = TRUE)
})
