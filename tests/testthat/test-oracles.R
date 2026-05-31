# Regression oracles: these values reproduce BUCSS 1.2.1 exactly and must not
# change. The ss_buc_* 2.0.0 rewrite is verified against them. The one
# documented exception is the between-subjects ANOVA (general) example, whose
# df_denominator now functions: at df_denominator = N - cells (no nuisance
# parameters) it reproduces the 1.2.1 omnibus value, which is what is pinned
# here.

ss  <- function(x) x$value[x$term == "necessary_sample_size"]
ncp <- function(x) x$value[x$term == "ncp_adjusted"]

expect_oracle <- function(result, sample_size, ncp_adjusted) {
  expect_equal(ss(result), sample_size)
  expect_equal(round(ncp(result), 3), ncp_adjusted)
}

test_that("independent t test reproduces 1.2.1 (scalar n)", {
  expect_oracle(
    ss_buc_indep_t(t_observed = 3, n = 20, alpha_prior = .05, alpha_planned = .05,
                   assurance = .80, power = .80),
    130, 1.105)
})

test_that("independent t test reproduces 1.2.1 (unequal n vector)", {
  expect_oracle(
    ss_buc_indep_t(t_observed = 3, n = c(50, 55), alpha_prior = .05,
                   alpha_planned = .05, power = .80, assurance = .90),
    1482, 0.526)
})

test_that("dependent t test reproduces 1.2.1", {
  expect_oracle(
    ss_buc_paired_t(t_observed = 3, N = 40, alpha_prior = .05, alpha_planned = .05,
                    assurance = .80, power = .80),
    255, 1.115)
})

test_that("one-way between-subjects ANOVA matches the omnibus computation", {
  expect_oracle(
    ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 4, alpha_prior = .05,
                         alpha_planned = .05, assurance = .80, power = .80),
    89, 3.737)
})

test_that("two-way between-subjects ANOVA reproduces 1.2.1", {
  expect_oracle(
    ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                           effect = "factor_B", alpha_prior = .05,
                           alpha_planned = .05, assurance = .80, power = .80),
    659, 0.293)
})

test_that("between-subjects ANOVA (general) reproduces the two-way result at df_denominator = N - cells", {
  expect_oracle(
    ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                   df_numerator = 2, df_denominator = 114,
                                   alpha_prior = .05, alpha_planned = .05,
                                   assurance = .80, power = .80),
    659, 0.293)
})

test_that("within-subjects ANOVA reproduces 1.2.1", {
  expect_oracle(
    ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3,
                    effect = "factor_B", alpha_prior = .05, alpha_planned = .05,
                    assurance = .80, power = .80),
    1904, 0.304)
})

test_that("within-subjects ANOVA (general) reproduces 1.2.1", {
  expect_oracle(
    ss_buc_rm_anova_general(F_observed = 6.5, N = 80, df_numerator = 1,
                            alpha_prior = .05, alpha_planned = .05,
                            assurance = .50, power = .80),
    256, 2.474)
})

test_that("split-plot ANOVA reproduces 1.2.1", {
  expect_oracle(
    ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
                       levels_within = 3, effect = "within", alpha_prior = .05,
                       alpha_planned = .05, assurance = .80, power = .80),
    968, 0.299)
})

test_that("split-plot ANOVA (general) reproduces 1.2.1", {
  expect_oracle(
    ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                               num_groups = 3, effect = "between_only",
                               df_num_within = 3, alpha_prior = .05,
                               alpha_planned = .05, assurance = .80, power = .80),
    1491, 0.194)
})

test_that("regression single coefficient reproduces 1.2.1", {
  expect_oracle(
    ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, alpha_prior = .05,
                    alpha_planned = .05, assurance = .80, power = .80),
    624, 1.893)
})

test_that("regression model R2 reproduces 1.2.1", {
  expect_oracle(
    ss_buc_R2(F_observed = 5, N = 150, p = 4, alpha_prior = .05,
              alpha_planned = .05, assurance = .80, power = .80),
    234, 7.816)
})

test_that("regression joint test reproduces 1.2.1", {
  expect_oracle(
    ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                     alpha_prior = .05, alpha_planned = .05, assurance = .80,
                     power = .80),
    3963, 0.365)
})

test_that("effect-size aliases are identical to their test-named functions", {
  expect_identical(
    ss_buc_smd(t_observed = 3, n = 20, assurance = .80, power = .80),
    ss_buc_indep_t(t_observed = 3, n = 20, assurance = .80, power = .80))
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
})
