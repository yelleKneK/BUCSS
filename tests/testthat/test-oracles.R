# Regression oracles. These values pin the BUCSS 2.0.0 output for each
# documented example so it cannot drift unnoticed. The bias and uncertainty
# adjusted noncentrality parameter is found by direct root finding (uniroot)
# rather than the fixed 0 to 100 grid used in 1.x, so the adjusted parameter is
# exact rather than snapped to the grid resolution. The rounded (three-decimal)
# noncentrality parameters are therefore identical to the 1.2.1 grid values; a
# handful of sample sizes move by a few units because the slightly different
# parameter propagates through the planned-n search. Cases whose sample size
# differs from the 1.2.1 grid value are flagged inline. See NEWS for the full
# accounting. The one exception is 'alpha_prior = 1', which 1.x computed as
# .999 and 2.0.0 computes as exactly 1 (no publication filter), so a plan at
# 'alpha_prior = 1' can differ from 1.2.1 in its parameter as well (flagged
# inline).

ss  <- function(x) x$value[x$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1]
ncp <- function(x) x$value[x$term == "ncp_adjusted"]

expect_oracle <- function(result, sample_size, ncp_adjusted) {
  expect_equal(ss(result), sample_size)
  expect_equal(round(ncp(result), 3), ncp_adjusted)
}

test_that("independent t test (scalar n)", {
  expect_oracle(
    ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = .05, alpha_planned = .05,
                   assurance = .80, desired_power = .80),
    130, 1.105)
})

test_that("independent t test (unequal n vector)", {
  # 1.2.1 grid gave 1482; uniroot gives 1485 (same rounded ncp).
  expect_oracle(
    ss_buc_independent_t(t_observed = 3, n = c(50, 55), alpha_prior = .05,
                   alpha_planned = .05, desired_power = .80, assurance = .90),
    1485, 0.526)
})

test_that("independent t test (scalar n, assurance .90)", {
  # The vignette "Assurance and the Publication Bias Correction" tabulates it.
  # 1.2.1 grid gave 5069; uniroot gives 5042 (same rounded ncp).
  expect_oracle(
    ss_buc_independent_t(t_observed = 3, n = 20, assurance = .90, desired_power = .80),
    5042, 0.176)
})

test_that("paired t test at alpha_prior = 1 (the safeguard power table)", {
  # The vignette "Assurance and the Publication Bias Correction" tabulates it:
  # a prior t with a two-sided p value of .15 from 30 pairs, planned at
  # alpha_prior = 1 and assurance .80. With no publication filter the adjusted
  # parameter solves P(|T| <= t | df 29, ncp) = .80: .5093695352, and 910
  # pairs (independent check 2026-10-13, base R only, exactly 1). The 1.2.1
  # grid gave 915, and the builds that still computed 1 as .999 gave 913 with
  # .5084551131.
  plan <- ss_buc_paired_t(t_observed = qt(1 - .15 / 2, df = 29), N = 30, alpha_prior = 1,
                          assurance = .80)
  expect_oracle(plan, 910, 0.509)
  expect_equal(ncp(plan), 0.5093695352, tolerance = 1e-7)
  # The row above it in the same table, a two-sided p value of .10: the
  # parameter is .7936462474 and the plan 376 pairs, whose power is .80026
  # against .79921 for 375 (independent check 2026-10-13, base R only, exactly
  # 1). The 1.2.1 grid gave 377 with .793, and so did the builds that computed
  # 1 as .999 (.7930526); NEWS lists both rows among the documented sizes
  # that moved.
  plan <- ss_buc_paired_t(t_observed = qt(1 - .10 / 2, df = 29), N = 30, alpha_prior = 1,
                          assurance = .80)
  expect_oracle(plan, 376, 0.794)
  expect_equal(ncp(plan), 0.7936462474, tolerance = 1e-7)
})

test_that("dependent t test", {
  expect_oracle(
    ss_buc_paired_t(t_observed = 3, N = 40, alpha_prior = .05, alpha_planned = .05,
                    assurance = .80, desired_power = .80),
    255, 1.115)
})

test_that("one-way between-subjects ANOVA", {
  expect_oracle(
    ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 4, alpha_prior = .05,
                         alpha_planned = .05, assurance = .80, desired_power = .80),
    89, 3.737)
})

test_that("two-way between-subjects ANOVA", {
  expect_oracle(
    ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                           effect = "factor_B", alpha_prior = .05,
                           alpha_planned = .05, assurance = .80, desired_power = .80),
    659, 0.293)
})

test_that("between-subjects ANOVA (general) matches the two-way result at df_denominator = N - cells", {
  expect_oracle(
    ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                   df_numerator = 2, df_denominator = 114,
                                   alpha_prior = .05, alpha_planned = .05,
                                   assurance = .80, desired_power = .80),
    659, 0.293)
})

test_that("within-subjects ANOVA", {
  # 1.2.1 grid gave 1904; uniroot gives 1902 (same rounded ncp).
  expect_oracle(
    ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3,
                    effect = "factor_B", alpha_prior = .05, alpha_planned = .05,
                    assurance = .80, desired_power = .80),
    1902, 0.304)
})

test_that("within-subjects ANOVA (general)", {
  expect_oracle(
    ss_buc_rm_anova_general(F_observed = 6.5, N = 80, df_numerator = 1,
                            alpha_prior = .05, alpha_planned = .05,
                            assurance = .50, desired_power = .80),
    256, 2.474)
})

test_that("split-plot ANOVA", {
  # 1.2.1 grid gave 968; uniroot gives 969 (same rounded ncp).
  expect_oracle(
    ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
                       levels_within = 3, effect = "within", alpha_prior = .05,
                       alpha_planned = .05, assurance = .80, desired_power = .80),
    969, 0.299)
})

test_that("split-plot ANOVA (general)", {
  # 1.2.1 grid gave 1491; uniroot gives 1489 (same rounded ncp).
  expect_oracle(
    ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                               num_groups = 3, effect = "between_only",
                               alpha_prior = .05,
                               alpha_planned = .05, assurance = .80, desired_power = .80),
    1489, 0.194)
})

test_that("regression single coefficient", {
  expect_oracle(
    ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, alpha_prior = .05,
                    alpha_planned = .05, assurance = .80, desired_power = .80),
    624, 1.893)
})

test_that("regression model R2", {
  expect_oracle(
    ss_buc_R2(F_observed = 5, N = 150, p = 4, alpha_prior = .05,
              alpha_planned = .05, assurance = .80, desired_power = .80),
    234, 7.816)
})

test_that("regression joint test", {
  # 1.2.1 grid gave 3963; uniroot gives 3960 (same rounded ncp).
  expect_oracle(
    ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                     alpha_prior = .05, alpha_planned = .05, assurance = .80,
                     desired_power = .80),
    3960, 0.365)
})

test_that("effect size aliases are identical to their test-named functions", {
  expect_identical(
    ss_buc_smd(t_observed = 3, n = 20, assurance = .80, desired_power = .80),
    ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, desired_power = .80))
  expect_identical(
    ss_buc_smd_paired(t_observed = 3, N = 40, assurance = .80, desired_power = .80),
    ss_buc_paired_t(t_observed = 3, N = 40, assurance = .80, desired_power = .80))
})

test_that("df_denominator functions: more nuisance parameters call for a larger sample", {
  base <- ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                         df_numerator = 2, df_denominator = 114,
                                         assurance = .80, desired_power = .80)
  with_covariates <- ss_buc_factorial_anova_general(F_observed = 5, N = 120,
                                                    cells = 6, df_numerator = 2,
                                                    df_denominator = 104,
                                                    assurance = .80, desired_power = .80)
  expect_gt(ss(with_covariates), ss(base))
  # df_denominator = 104 with N - cells = 114 carries 10 nuisance parameters
  # (e.g., covariates) into the planned study. These pin the 2.0.0 output so the
  # functional df_denominator path cannot drift unnoticed.
  expect_oracle(with_covariates, 737, 0.262)
})

test_that("the adjusted noncentrality parameter is no longer capped at 100", {
  # The flagship engine fix: the 1.x grid stopped at 100, so this prior (whose
  # implied parameter is 148.48) returned the capped, wrong pair 5 / 100. The
  # uniroot bracket auto-expands, verified against an independent root solve.
  expect_oracle(
    ss_buc_one_way_anova(F_observed = 60, N = 120, levels_A = 4,
                         assurance = .80, desired_power = .80),
    4, 148.482)
})

test_that("the round-down branch binds when the prior N does not divide evenly", {
  # With N = 121 in 3 groups the conservative two-sided rounding matters: the
  # round-up branch alone would plan 1277; the round-down branch requires 1280,
  # and the max-n / min-NCP rule must select it. Guards the conservatism rule,
  # which no evenly divisible case can (the branches coincide there).
  expect_oracle(
    ss_buc_one_way_anova(F_observed = 5, N = 121, levels_A = 3,
                         assurance = .80, desired_power = .80),
    1280, 0.302)
})

test_that("previously unpinned planner arms match their verified values", {
  # Each of these three arms had no oracle: the within_only effect of the
  # general split-plot planner, the one-way within-subjects success path, and
  # the independent t total-N path with odd N (its own two-branch rounding).
  # All three values are verified matches of the CRAN 1.2.1 grid at step .001.
  expect_oracle(
    ss_buc_mixed_anova_general(F_observed = 5, N = 90, num_groups = 3,
                               df_numerator = 2, effect = "within_only",
                               assurance = .80, desired_power = .80),
    704, 0.411)
  expect_oracle(
    ss_buc_rm_anova(F_observed = 6.5, N = 80, levels_A = 3,
                    assurance = .80, desired_power = .80),
    200, 3.897)
  expect_oracle(
    ss_buc_independent_t(t_observed = 3, N = 41, assurance = .80, desired_power = .80),
    132, 1.105)
})

test_that("the sphericity correction (epsilon) pins its plans", {
  # epsilon below 1 applies the Muller and Barton (1989) approximation: both
  # degrees of freedom and the noncentrality parameter are multiplied by
  # epsilon, for the prior study and the planned one, with ncp_adjusted on the
  # sphericity scale. Pinned to full precision (the size exactly, the
  # parameter to 1e-10); the spec's independent prototype gave the same
  # values to ten decimals. The first is the help page's and the vignette's
  # example (F(1.73, 50.2) with 30 subjects), which plans 92 under
  # sphericity; the fourth plans 60 per group under sphericity.
  pin <- function(result, sample_size, ncp_adjusted) {
    expect_equal(ss(result), sample_size)
    expect_equal(ncp(result), ncp_adjusted, tolerance = 1e-10)
  }
  pin(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = 1.73 / 2),
      162, 2.0025989386642826)
  pin(ss_buc_rm_anova(F_observed = 9, N = 40, levels_A = 2, levels_B = 3,
                      effect = "interaction", epsilon = .7),
      82, 6.1816562006893037)
  pin(ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2, epsilon = .7),
      109, 3.4731264397829271)
  pin(ss_buc_mixed_anova(F_observed = 7, N = 60, levels_between = 2, levels_within = 3,
                         effect = "within", epsilon = .8),
      108, 3.1743021455404503)
  pin(ss_buc_mixed_anova(F_observed = 8, N = 61, levels_between = 2, levels_within = 3,
                         effect = "interaction", epsilon = .8),
      64, 5.5142729566147777)
  pin(ss_buc_mixed_anova_general(F_observed = 8, N = 91, df_numerator = 2, num_groups = 3,
                                 effect = "within_only", epsilon = .75),
      71, 5.2106738807869997)
  pin(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 6, num_groups = 3,
                                 effect = "between_within", df_num_within = 3,
                                 epsilon = .6),
      56, 10.512636680923496)
  expect_oracle(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3), 92, 3.213)
  expect_oracle(ss_buc_mixed_anova(F_observed = 7, N = 60, levels_between = 2,
                                   levels_within = 3, effect = "within"), 60, 4.937)
})
