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

test_that("endpoint assurance and power are rejected, not silently planned", {
  # assurance = 0 has no root (TM is strictly positive): the old behavior
  # returned the uniroot bracket endpoint as the adjusted NCP.
  expect_error(ss_buc_independent_t(t_observed = 3, n = 50, assurance = 0),
               "assurance", fixed = TRUE)
  # desired_power = 100 coerces to a target of exactly 1.0, which the search can only
  # "reach" by floating-point error in pt(); desired_power = 0 trivially returned n = 2.
  expect_error(ss_buc_paired_t(t_observed = 3, N = 25, desired_power = 100),
               "power", fixed = TRUE)
  expect_error(ss_buc_paired_t(t_observed = 3, N = 25, desired_power = 0),
               "power", fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, assurance = 100),
               "assurance", fixed = TRUE)
  # The documented percentage rule is untouched: exactly 1 still means 1
  # percent (coercion runs before the range check), so desired_power = 1 plans at a
  # .01 target without error, and assurance = 1 reaches the < .5 warning
  # rather than the range stop.
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, desired_power = 1), NA)
  expect_warning(ss_buc_independent_t(t_observed = 3, n = 20, assurance = 1),
                 "< .5", fixed = TRUE)
})

test_that("assurance below .5 warns but still computes", {
  expect_warning(ss_buc_independent_t(t_observed = 3, n = 20, assurance = .4),
                 "< .5", fixed = TRUE)
})

test_that("out-of-range power is rejected", {
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, desired_power = -.1),
               "power", fixed = TRUE)
})

test_that("a nonsignificant prior t is rejected", {
  expect_error(ss_buc_independent_t(t_observed = 1, n = 20),
               "nonsignificant", fixed = TRUE)
})

test_that("a negative prior t plans from its magnitude (two-sided rule)", {
  # The publication rule is two-sided and TM is symmetric in the sign of
  # t, so t = -3 is the same evidence as t = 3; the sign is the analyst's
  # coding of the comparison, and the planned sample size must not depend
  # on it. The signed value is still echoed in the stored inputs.
  neg <- ss_buc_independent_t(t_observed = -3, n = 20)
  pos <- ss_buc_independent_t(t_observed = 3, n = 20)
  expect_equal(neg$value[neg$term != "t_observed"],
               pos$value[pos$term != "t_observed"])
  expect_identical(neg$value[neg$term == "t_observed"][1], -3)

  neg_p <- ss_buc_paired_t(t_observed = -3, N = 40)
  pos_p <- ss_buc_paired_t(t_observed = 3, N = 40)
  expect_equal(neg_p$value[neg_p$term != "t_observed"],
               pos_p$value[pos_p$term != "t_observed"])
  expect_identical(neg_p$value[neg_p$term == "t_observed"][1], -3)

  neg_r <- ss_buc_reg_coef(t_observed = -3, N = 150, p = 3)
  pos_r <- ss_buc_reg_coef(t_observed = 3, N = 150, p = 3)
  expect_equal(neg_r$value[neg_r$term != "t_observed"],
               pos_r$value[pos_r$term != "t_observed"])
  expect_identical(neg_r$value[neg_r$term == "t_observed"][1], -3)

  # A magnitude below the two-tailed cutoff is still rejected, whatever
  # its sign, and the alpha_prior remedy in the message is now reachable.
  expect_error(ss_buc_independent_t(t_observed = -1, n = 20),
               "nonsignificant", fixed = TRUE)
  expect_error(ss_buc_paired_t(t_observed = -1, N = 40),
               "nonsignificant", fixed = TRUE)
  expect_error(ss_buc_reg_coef(t_observed = -1, N = 150, p = 3),
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

test_that("within-subjects ANOVA rejects a non-default effect when 'levels_B' is absent", {
  # Previously this silently computed the single-factor result; now it stops.
  expect_error(ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 3, effect = "factor_B"),
               "levels_B", fixed = TRUE)
  expect_error(ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 3, effect = "interaction"),
               "levels_B", fixed = TRUE)
})

test_that("the general split-plot function requires 'df_num_within' for a between-within effect", {
  expect_error(
    ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                               num_groups = 3, effect = "between_within"),
    "df_num_within", fixed = TRUE)
})

test_that("the general planners require their structural arguments", {
  expect_error(ss_buc_factorial_anova_general(F_observed = 5, N = 120,
                                              df_numerator = 2, df_denominator = 114),
               "must specify 'cells'", fixed = TRUE)
  expect_error(ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                              df_denominator = 114),
               "must specify 'df_numerator'", fixed = TRUE)
  expect_error(ss_buc_rm_anova_general(F_observed = 6.5, N = 80),
               "must specify 'df_numerator'", fixed = TRUE)
  expect_error(ss_buc_mixed_anova_general(F_observed = 5, N = 90, num_groups = 3),
               "must specify 'df_numerator'", fixed = TRUE)
  expect_error(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2),
               "must specify 'num_groups'", fixed = TRUE)
})

test_that("the regression planners require 'p' (and 'p_joint')", {
  expect_error(ss_buc_reg_coef(t_observed = 3, N = 150), "specify 'p'", fixed = TRUE)
  expect_error(ss_buc_R2(F_observed = 5, N = 150), "specify 'p'", fixed = TRUE)
  expect_error(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4),
               "specify 'p_joint'", fixed = TRUE)
})

test_that("non-scalar, non-finite, and non-whole inputs get friendly errors", {
  # A fractional predictor count used to flow through the planned-n search and
  # return a non-integer "necessary sample size" (624.5).
  expect_error(ss_buc_reg_coef(t_observed = 3, N = 150, p = 2.5),
               "whole number", fixed = TRUE)
  expect_error(ss_buc_one_way_anova(F_observed = 5, N = 120.5, levels_A = 3),
               "whole number", fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20.5),
               "whole number", fixed = TRUE)
  # NA, Inf, and vector inputs used to die inside base R ("missing value where
  # TRUE/FALSE needed", "the condition has length > 1").
  expect_error(ss_buc_paired_t(t_observed = NA, N = 40),
               "single finite numeric value", fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = Inf, n = 20),
               "single finite numeric value", fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, assurance = NA),
               "single finite numeric value", fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20,
                                    alpha_prior = c(.05, .10)),
               "single finite numeric value", fixed = TRUE)
  expect_error(ss_buc_paired_t(t_observed = 3, N = NA),
               "whole number", fixed = TRUE)
})

test_that("degenerate sample sizes are rejected before any qt/qf call", {
  # rm designs: N = 1 used to reach qf() with 0 denominator df (NaN) and die
  # with "missing value where TRUE/FALSE needed".
  expect_error(ss_buc_rm_anova(F_observed = 200, N = 1, levels_A = 3,
                               alpha_prior = 1, assurance = .5),
               "at least 2", fixed = TRUE)
  expect_error(ss_buc_rm_anova_general(F_observed = 200, N = 1,
                                       df_numerator = 2),
               "at least 2", fixed = TRUE)
  # independent t: N = 3 used to reach qt() with 0 df on the round-down branch.
  expect_error(ss_buc_independent_t(t_observed = 3, N = 3),
               "at least 4", fixed = TRUE)
})

test_that("a prior 'N' too small for the design is rejected with a clear message", {
  expect_error(ss_buc_one_way_anova(F_observed = 5, N = 5, levels_A = 4),
               "too small", fixed = TRUE)
  expect_error(ss_buc_mixed_anova(F_observed = 5, N = 3, levels_between = 2,
                                  levels_within = 3, effect = "within"),
               "too small", fixed = TRUE)
  expect_error(ss_buc_factorial_anova_general(F_observed = 5, N = 11, cells = 6,
                                              df_numerator = 2, df_denominator = 5),
               "too small", fixed = TRUE)
})

test_that("every remaining documented guard fires with its message", {
  # Data-driven sweep of the stop() branches not covered by the targeted tests
  # above: one row per (call, fixed substring the error must contain). The
  # nonsignificant-statistic rows are the load-bearing ones; if one of those
  # guards silently vanished, planning would proceed from a nonsignificant
  # prior.
  guard_cases <- list(
    # validator: percentages above 100 are out of range after the coercion
    list(quote(ss_buc_independent_t(t_observed = 3, n = 20, assurance = 150)),
         "assurance"),
    list(quote(ss_buc_independent_t(t_observed = 3, n = 20, desired_power = 150)),
         "power"),
    # independent t structural guards
    list(quote(ss_buc_independent_t(t_observed = 3, n = 20, N = 40)),
         "should not specify 'n'"),
    list(quote(ss_buc_independent_t(t_observed = 3, n = c(10, 10, 10))),
         "vector of length two"),
    # paired t
    list(quote(ss_buc_paired_t(t_observed = 3)), "number of pairs"),
    list(quote(ss_buc_paired_t(t_observed = 1, N = 40)), "nonsignificant"),
    # one-way ANOVA
    list(quote(ss_buc_one_way_anova(F_observed = 5, levels_A = 4)),
         "must specify 'N'"),
    list(quote(ss_buc_one_way_anova(F_observed = 5, N = 120)), "levels_A"),
    list(quote(ss_buc_one_way_anova(F_observed = 1, N = 120, levels_A = 4)),
         "nonsignificant"),
    # factorial ANOVA
    list(quote(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_B = 3)),
         "levels_A"),
    list(quote(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2)),
         "levels_B"),
    list(quote(ss_buc_factorial_anova(F_observed = 5, N = 10, levels_A = 2,
                                      levels_B = 3)),
         "too small"),
    # factorial general
    list(quote(ss_buc_factorial_anova_general(F_observed = 5, cells = 6,
                                              df_numerator = 2,
                                              df_denominator = 114)),
         "must specify 'N'"),
    list(quote(ss_buc_factorial_anova_general(F_observed = 1, N = 120,
                                              cells = 6, df_numerator = 2,
                                              df_denominator = 114)),
         "nonsignificant"),
    # within-subjects designs
    list(quote(ss_buc_rm_anova(F_observed = 5, levels_A = 3)),
         "must specify 'N'"),
    list(quote(ss_buc_rm_anova(F_observed = 5, N = 60)), "levels_A"),
    list(quote(ss_buc_rm_anova(F_observed = 1, N = 60, levels_A = 3)),
         "nonsignificant"),
    list(quote(ss_buc_rm_anova_general(F_observed = 5, df_numerator = 2)),
         "must specify 'N'"),
    list(quote(ss_buc_rm_anova_general(F_observed = 1, N = 60,
                                       df_numerator = 2)),
         "nonsignificant"),
    # split-plot designs
    list(quote(ss_buc_mixed_anova(F_observed = 5, levels_between = 2,
                                  levels_within = 3)),
         "must specify 'N'"),
    list(quote(ss_buc_mixed_anova(F_observed = 5, N = 60,
                                  levels_between = 2)),
         "within-subjects factor"),
    list(quote(ss_buc_mixed_anova(F_observed = 5, N = 60,
                                  levels_within = 3)),
         "between-subjects factor"),
    list(quote(ss_buc_mixed_anova(F_observed = 1, N = 60, levels_between = 2,
                                  levels_within = 3)),
         "nonsignificant"),
    list(quote(ss_buc_mixed_anova_general(F_observed = 5, df_numerator = 1,
                                          num_groups = 2)),
         "must specify 'N'"),
    list(quote(ss_buc_mixed_anova_general(F_observed = 5, N = 60,
                                          df_numerator = 1, num_groups = 2,
                                          effect = "between_within")),
         "df_num_within"),
    list(quote(ss_buc_mixed_anova_general(F_observed = 5, N = 3,
                                          df_numerator = 1, num_groups = 2)),
         "too small"),
    list(quote(ss_buc_mixed_anova_general(F_observed = 1, N = 60,
                                          df_numerator = 1, num_groups = 2)),
         "nonsignificant"),
    # regression designs
    list(quote(ss_buc_reg_coef(t_observed = 3, p = 3)), "total sample size"),
    list(quote(ss_buc_reg_coef(t_observed = 1, N = 150, p = 3)),
         "nonsignificant"),
    list(quote(ss_buc_reg_coef(t_observed = 3, N = 10, p = 9)),
         "degrees of freedom"),
    list(quote(ss_buc_R2(F_observed = 5, p = 3)), "total sample size"),
    list(quote(ss_buc_R2(F_observed = 1, N = 150, p = 3)), "nonsignificant"),
    list(quote(ss_buc_reg_joint(F_observed = 5, N = 150, p = 3, p_joint = 4)),
         "cannot exceed"),
    list(quote(ss_buc_reg_joint(F_observed = 1, N = 150, p = 4, p_joint = 2)),
         "nonsignificant"),
    # an unrecognized effect value falls through .match_effect to match.arg
    list(quote(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
                                      levels_B = 3, effect = "bogus")),
         "should be one of")
  )
  for (case in guard_cases) {
    expect_error(suppressWarnings(eval(case[[1]])), case[[2]], fixed = TRUE,
                 label = paste(deparse(case[[1]]), collapse = " "))
  }
})

test_that("the Mode 2 zero-NCP message says alpha_prior is the only lever", {
  # A barely significant prior has a ceiling below the .5 assurance floor, so
  # the message must switch from "re-run at or below the ceiling" (Mode 1) to
  # naming alpha_prior as the only remaining remedy.
  expect_error(ss_buc_independent_t(t_observed = 2.03, n = 20,
                                    assurance = .95),
               "only remaining lever", fixed = TRUE)
})

test_that("the noncentrality safety cap warns and names the remedy", {
  # A finite but absurd statistic drives the root past the 1e7 bracket cap;
  # the planner warns (rather than looping forever) and asks the user to
  # verify the inputs. The paired design has a single branch, so exactly one
  # warning fires.
  expect_warning(ss_buc_paired_t(t_observed = 5e7, N = 25),
                 "verify", fixed = TRUE)
})
