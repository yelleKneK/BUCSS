# Regression oracles. These values pin the BUCSS 2.0.0 output for each
# documented example so it cannot drift unnoticed. The bias and uncertainty
# adjusted noncentrality parameter is found by direct root finding (uniroot)
# rather than the fixed 0 to 100 grid used in 1.x, so the adjusted parameter is
# exact rather than snapped to the grid resolution. The rounded (three-decimal)
# noncentrality parameters are therefore identical to the 1.2.1 grid values; a
# handful of sample sizes move by a few units because the slightly different
# parameter propagates through the planned-n search. Cases whose sample size
# differs from the 1.2.1 grid value are flagged inline. See NEWS for the full
# accounting.

ss  <- function(x) x$value[x$term == "necessary_sample_size"]
ncp <- function(x) x$value[x$term == "ncp_adjusted"]

expect_oracle <- function(result, sample_size, ncp_adjusted) {
  expect_equal(ss(result), sample_size)
  expect_equal(round(ncp(result), 3), ncp_adjusted)
}

test_that("independent t test (scalar n)", {
  expect_oracle(
    ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = .05, alpha_planned = .05,
                   assurance = .80, power = .80),
    130, 1.105)
})

test_that("independent t test (unequal n vector)", {
  # 1.2.1 grid gave 1482; uniroot gives 1485 (same rounded ncp).
  expect_oracle(
    ss_buc_independent_t(t_observed = 3, n = c(50, 55), alpha_prior = .05,
                   alpha_planned = .05, power = .80, assurance = .90),
    1485, 0.526)
})

test_that("dependent t test", {
  expect_oracle(
    ss_buc_paired_t(t_observed = 3, N = 40, alpha_prior = .05, alpha_planned = .05,
                    assurance = .80, power = .80),
    255, 1.115)
})

test_that("one-way between-subjects ANOVA", {
  expect_oracle(
    ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 4, alpha_prior = .05,
                         alpha_planned = .05, assurance = .80, power = .80),
    89, 3.737)
})

test_that("two-way between-subjects ANOVA", {
  expect_oracle(
    ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                           effect = "factor_B", alpha_prior = .05,
                           alpha_planned = .05, assurance = .80, power = .80),
    659, 0.293)
})

test_that("between-subjects ANOVA (general) matches the two-way result at df_denominator = N - cells", {
  expect_oracle(
    ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                   df_numerator = 2, df_denominator = 114,
                                   alpha_prior = .05, alpha_planned = .05,
                                   assurance = .80, power = .80),
    659, 0.293)
})

test_that("within-subjects ANOVA", {
  # 1.2.1 grid gave 1904; uniroot gives 1902 (same rounded ncp).
  expect_oracle(
    ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3,
                    effect = "factor_B", alpha_prior = .05, alpha_planned = .05,
                    assurance = .80, power = .80),
    1902, 0.304)
})

test_that("within-subjects ANOVA (general)", {
  expect_oracle(
    ss_buc_rm_anova_general(F_observed = 6.5, N = 80, df_numerator = 1,
                            alpha_prior = .05, alpha_planned = .05,
                            assurance = .50, power = .80),
    256, 2.474)
})

test_that("split-plot ANOVA", {
  # 1.2.1 grid gave 968; uniroot gives 969 (same rounded ncp).
  expect_oracle(
    ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
                       levels_within = 3, effect = "within", alpha_prior = .05,
                       alpha_planned = .05, assurance = .80, power = .80),
    969, 0.299)
})

test_that("split-plot ANOVA (general)", {
  # 1.2.1 grid gave 1491; uniroot gives 1489 (same rounded ncp).
  expect_oracle(
    ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                               num_groups = 3, effect = "between_only",
                               df_num_within = 3, alpha_prior = .05,
                               alpha_planned = .05, assurance = .80, power = .80),
    1489, 0.194)
})

test_that("regression single coefficient", {
  expect_oracle(
    ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, alpha_prior = .05,
                    alpha_planned = .05, assurance = .80, power = .80),
    624, 1.893)
})

test_that("regression model R2", {
  expect_oracle(
    ss_buc_R2(F_observed = 5, N = 150, p = 4, alpha_prior = .05,
              alpha_planned = .05, assurance = .80, power = .80),
    234, 7.816)
})

test_that("regression joint test", {
  # 1.2.1 grid gave 3963; uniroot gives 3960 (same rounded ncp).
  expect_oracle(
    ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                     alpha_prior = .05, alpha_planned = .05, assurance = .80,
                     power = .80),
    3960, 0.365)
})

test_that("effect-size aliases are identical to their test-named functions", {
  expect_identical(
    ss_buc_smd(t_observed = 3, n = 20, assurance = .80, power = .80),
    ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, power = .80))
  expect_identical(
    ss_buc_smd_paired(t_observed = 3, N = 40, assurance = .80, power = .80),
    ss_buc_paired_t(t_observed = 3, N = 40, assurance = .80, power = .80))
})

test_that("df_denominator functions: more nuisance parameters call for a larger sample", {
  base <- ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                         df_numerator = 2, df_denominator = 114,
                                         assurance = .80, power = .80)
  with_covariates <- ss_buc_factorial_anova_general(F_observed = 5, N = 120,
                                                    cells = 6, df_numerator = 2,
                                                    df_denominator = 104,
                                                    assurance = .80, power = .80)
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
                         assurance = .80, power = .80),
    4, 148.482)
})

test_that("the round-down branch binds when the prior N does not divide evenly", {
  # With N = 121 in 3 groups the conservative two-sided rounding matters: the
  # round-up branch alone would plan 1277; the round-down branch requires 1280,
  # and the max-n / min-NCP rule must select it. Guards the conservatism rule,
  # which no evenly divisible case can (the branches coincide there).
  expect_oracle(
    ss_buc_one_way_anova(F_observed = 5, N = 121, levels_A = 3,
                         assurance = .80, power = .80),
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
                               assurance = .80, power = .80),
    704, 0.411)
  expect_oracle(
    ss_buc_rm_anova(F_observed = 6.5, N = 80, levels_A = 3,
                    assurance = .80, power = .80),
    200, 3.897)
  expect_oracle(
    ss_buc_independent_t(t_observed = 3, N = 41, assurance = .80, power = .80),
    132, 1.105)
})
