# Guard coverage added after the 2026-09-28 deep QC audit. The audit
# instrumented every stop(), warning(), zero-NCP refusal, and argument check
# and found 105 of 218 sites that no test ever made fire (inventory:
# Quality Control/BUCSS/BUCSS_QC_deep_2026-09-28_evidence/guard_inventory/).
# Each row below makes one of those sites, or one of the guards the fix wave
# added, fire with its own message. Messages are matched as fixed substrings;
# the refusals are also matched by condition class.
#
# The round-up nonsignificance checks in ss_buc_mixed_anova() and
# ss_buc_mixed_anova_general() come after the round-down check, but they are
# reachable. The critical F at 1 - alpha_prior falls as the error degrees of
# freedom grow only when it is an upper quantile; with a large 'alpha_prior'
# (1 is read as .999) it is a lower quantile, and with several numerator
# degrees of freedom a lower quantile rises. Their own test is at the end of
# this file.

expect_guard <- function(call, fragment, class = NULL) {
  label <- paste(deparse(call, width.cutoff = 500L), collapse = " ")
  if (is.null(class)) {
    expect_error(suppressWarnings(eval(call, parent.frame())), fragment,
                 fixed = TRUE, label = label)
  } else {
    expect_error(suppressWarnings(eval(call, parent.frame())), fragment,
                 fixed = TRUE, class = class, label = label)
  }
}

test_that("argument checks the audit found never fire now fire", {
  cases <- list(
    # the shared validator
    list(quote(ss_buc_paired_t(t_observed = 3, N = 40, alpha_planned = NA)),
         "'alpha_planned' must be a single finite numeric value"),
    list(quote(ss_buc_paired_t(t_observed = 3, N = 40, desired_power = c(.8, .9))),
         "'desired_power' must be a single finite numeric value"),
    list(quote(print(ss_buc_paired_t(t_observed = 3, N = 40), digits = 0)),
         "'digits' must be a single positive number"),
    # the one-sample t names its own sample size, not the paired t's pairs
    list(quote(ss_buc_one_sample_t(t_observed = 3)), "specify 'N'"),
    # one-way and two-way between subjects
    list(quote(ss_buc_one_way_anova(F_observed = NA, N = 120, levels_A = 4)),
         "'F_observed' must be a single finite numeric value"),
    list(quote(ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 2.5)),
         "'levels_A' must be a single whole number"),
    list(quote(ss_buc_factorial_anova(F_observed = Inf, N = 120, levels_A = 2, levels_B = 3)),
         "'F_observed' must be a single finite numeric value"),
    list(quote(ss_buc_factorial_anova(F_observed = 5, N = 120.5, levels_A = 2, levels_B = 3)),
         "'N' must be a single whole number"),
    list(quote(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 1, levels_B = 3)),
         "'levels_A' must be at least 2"),
    list(quote(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 1)),
         "'levels_B' must be at least 2"),
    # between-subjects general
    list(quote(ss_buc_factorial_anova_general(F_observed = NA, N = 120, cells = 6,
                                              df_numerator = 2, df_denominator = 114)),
         "'F_observed' must be a single finite numeric value"),
    list(quote(ss_buc_factorial_anova_general(F_observed = 5, N = 1, cells = 6,
                                              df_numerator = 2, df_denominator = 114)),
         "'N' must be at least 2"),
    list(quote(ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 1,
                                              df_numerator = 2, df_denominator = 114)),
         "'cells' must be at least 2"),
    list(quote(ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                              df_numerator = 0, df_denominator = 114)),
         "'df_numerator' must be at least 1"),
    # ANCOVA
    list(quote(ss_buc_ancova(F_observed = 6, cells = 3, n_covariates = 2)),
         "You must specify 'N'"),
    list(quote(ss_buc_ancova(F_observed = 6, N = 120, n_covariates = 2)),
         "You must specify 'cells'"),
    list(quote(ss_buc_ancova(F_observed = 6, N = 120, cells = 3)),
         "You must specify 'n_covariates'"),
    list(quote(ss_buc_ancova(F_observed = 6, N = 120, cells = 2.5, n_covariates = 2)),
         "'cells' must be a single whole number"),
    list(quote(ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = -1)),
         "'n_covariates' must be at least 0"),
    list(quote(ss_buc_ancova(F_observed = 6, N = 1, cells = 3, n_covariates = 2)),
         "'N' must be at least 2"),
    list(quote(ss_buc_ancova(F_observed = 6, N = 6, cells = 3, n_covariates = 3)),
         "leaves no error degrees of freedom in the prior study"),
    list(quote(ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2,
                             df_numerator = 0)),
         "'df_numerator' must be at least 1"),
    # within subjects
    list(quote(ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 3, effect = "bogus")),
         "should be one of"),
    list(quote(ss_buc_rm_anova(F_observed = NA, N = 60, levels_A = 3)),
         "'F_observed' must be a single finite numeric value"),
    list(quote(ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 1)),
         "'levels_A' must be at least 2"),
    list(quote(ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 1,
                               effect = "factor_B")),
         "'levels_B' must be at least 2"),
    list(quote(ss_buc_rm_anova_general(F_observed = NA, N = 80, df_numerator = 1)),
         "'F_observed' must be a single finite numeric value"),
    list(quote(ss_buc_rm_anova_general(F_observed = 6.5, N = 80, df_numerator = 0)),
         "'df_numerator' must be at least 1"),
    # split plot
    list(quote(ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
                                  levels_within = 3, effect = "bogus")),
         "should be one of"),
    list(quote(ss_buc_mixed_anova(F_observed = NA, N = 60, levels_between = 2,
                                  levels_within = 3)),
         "'F_observed' must be a single finite numeric value"),
    list(quote(ss_buc_mixed_anova(F_observed = 5, N = 60.5, levels_between = 2,
                                  levels_within = 3)),
         "'N' must be a single whole number"),
    list(quote(ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 1,
                                  levels_within = 3)),
         "'levels_between' must be at least 2"),
    list(quote(ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
                                  levels_within = 1)),
         "'levels_within' must be at least 2"),
    list(quote(ss_buc_mixed_anova_general(F_observed = NA, N = 90, df_numerator = 2,
                                          num_groups = 3)),
         "'F_observed' must be a single finite numeric value"),
    list(quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90.5, df_numerator = 2,
                                          num_groups = 3)),
         "'N' must be a single whole number"),
    list(quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                                          num_groups = 1)),
         "'num_groups' must be at least 2"),
    list(quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 0,
                                          num_groups = 3)),
         "'df_numerator' must be at least 1"),
    list(quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                                          num_groups = 3, effect = "between_within",
                                          df_num_within = 0)),
         "'df_num_within' must be at least 1"),
    # regression
    list(quote(ss_buc_reg_coef(t_observed = NA, N = 150, p = 3)),
         "'t_observed' must be a single finite numeric value"),
    list(quote(ss_buc_reg_coef(t_observed = 3, N = 150.5, p = 3)),
         "'N' must be a single whole number"),
    list(quote(ss_buc_R2(F_observed = NA, N = 150, p = 4)),
         "'F_observed' must be a single finite numeric value"),
    list(quote(ss_buc_R2(F_observed = 5, N = 1, p = 4)),
         "'N' must be at least 2"),
    list(quote(ss_buc_R2(F_observed = 5, N = 150, p = 0)),
         "'p' must be at least 1"),
    list(quote(ss_buc_R2(F_observed = 5, N = 5, p = 4)),
         "0 or negative degrees of freedom"),
    list(quote(ss_buc_reg_joint(F_observed = 5, p = 4, p_joint = 2)),
         "specify 'N'"),
    list(quote(ss_buc_reg_joint(F_observed = 5, N = 150, p_joint = 2)),
         "specify 'p'"),
    list(quote(ss_buc_reg_joint(F_observed = NA, N = 150, p = 4, p_joint = 2)),
         "'F_observed' must be a single finite numeric value"),
    list(quote(ss_buc_reg_joint(F_observed = 5, N = 150.5, p = 4, p_joint = 2)),
         "'N' must be a single whole number"),
    list(quote(ss_buc_reg_joint(F_observed = 5, N = 150, p = 0, p_joint = 2)),
         "'p' must be at least 1"),
    list(quote(ss_buc_reg_joint(F_observed = 5, N = 5, p = 4, p_joint = 2)),
         "0 or negative degrees of freedom"),
    list(quote(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 0)),
         "'p_joint' must be at least 1"),
    # correlation
    list(quote(ss_buc_correlation(r_observed = .35)), "specify 'N'"),
    list(quote(ss_buc_correlation(r_observed = .35, N = 3)), "'N' must be at least 4"),
    list(quote(ss_buc_correlation(N = 100)), "either 'r_observed'"),
    list(quote(ss_buc_correlation(r_observed = .35, t_observed = 3.7, N = 100)),
         "not both"),
    list(quote(ss_buc_correlation(r_observed = NA, N = 100)),
         "'r_observed' must be a single finite numeric value"),
    list(quote(ss_buc_correlation(r_observed = 1, N = 100)),
         "strictly between -1 and 1"),
    list(quote(ss_buc_correlation(t_observed = Inf, N = 100)),
         "'t_observed' must be a single finite numeric value"),
    # MANOVA
    list(quote(ss_buc_manova(F_observed = 4.5, p_variables = 3)), "You must specify 'N'"),
    list(quote(ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 1.5)),
         "'p_variables' must be a single whole number"),
    list(quote(ss_buc_manova(F_observed = 4.5, N = 3, p_variables = 1)),
         "'N' must be at least 4"),
    list(quote(ss_buc_manova(N = 80, p_variables = 3)), "either 'F_observed' or 'T2_observed'"),
    list(quote(ss_buc_manova(F_observed = 4.5, T2_observed = 14, N = 80, p_variables = 3)),
         "not both"),
    list(quote(ss_buc_manova(T2_observed = NA, N = 80, p_variables = 3)),
         "'T2_observed' must be a single finite numeric value"),
    list(quote(ss_buc_manova(F_observed = c(4.5, 5), N = 80, p_variables = 3)),
         "'F_observed' must be a single finite numeric value"),
    list(quote(ss_buc_manova(F_observed = 4.5, N = 5, p_variables = 4)),
         "leaves no error degrees of freedom"),
    # an odd N passes the check above (5 - 3 - 1 = 1) and reaches the one
    # inside the rounding branches, where 2 * floor(5 / 2) - 3 - 1 = 0
    list(quote(ss_buc_manova(F_observed = 40, N = 5, p_variables = 3)),
         "leaves no error degrees of freedom"),
    # chi square difference
    list(quote(ss_buc_chisq_diff(chisq_observed = 9.5, df_full = 42, df_restricted = 43)),
         "You must specify 'N'"),
    list(quote(ss_buc_chisq_diff(chisq_observed = NA, N = 250, df_full = 42,
                                 df_restricted = 43)),
         "'chisq_observed' must be a single finite numeric value"),
    list(quote(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250.5, df_full = 42,
                                 df_restricted = 43)),
         "'N' must be a single whole number"),
    list(quote(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = -1,
                                 df_restricted = 43)),
         "'df_full' must be at least 0"),
    list(quote(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 1.5,
                                 df_restricted = 43)),
         "'df_full' must be a single whole number"),
    list(quote(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42,
                                 df_restricted = 0)),
         "'df_restricted' must be at least 1"),
    # Welch
    list(quote(ss_buc_welch_t(t_observed = NA, n_1 = 40, n_2 = 55)),
         "'t_observed' must be a single finite numeric value"),
    list(quote(ss_buc_welch_t(t_observed = 3, n_1 = 1, n_2 = 55)),
         "'n_1' must be at least 2"),
    list(quote(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55.5)),
         "'n_2' must be a single whole number"),
    list(quote(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_1 = NA, sd_2 = 12.6)),
         "'sd_1' must be a single finite numeric value"),
    list(quote(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_1 = 8.4, sd_2 = Inf)),
         "'sd_2' must be a single finite numeric value"),
    list(quote(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, var_1 = NA,
                              var_2 = 158.76)),
         "'var_1' must be a single finite numeric value"),
    list(quote(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, var_1 = 70.56,
                              var_2 = c(1, 2))),
         "'var_2' must be a single finite numeric value"),
    list(quote(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = NA)),
         "'sd_ratio' must be a single finite numeric value")
  )
  for (case in cases) expect_guard(case[[1]], case[[2]])
})

test_that("the design guards added by the fix wave fire", {
  cases <- list(
    # df_numerator cannot exceed the degrees of freedom among the cell means
    list(quote(ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                              df_numerator = 6, df_denominator = 114)),
         "'df_numerator' cannot exceed 'cells' - 1"),
    list(quote(ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2,
                             df_numerator = 3)),
         "the omnibus test uses 'cells' - 1 and a contrast uses 1"),
    # ANCOVA checks its own error df, in its own terms (it has no
    # 'df_denominator' argument to tell the user about): N - cells -
    # n_covariates = 1 passes the first check, but at the rounded-down cell
    # size (floor(7 / 3) - 1) * 3 - 3 = 0
    list(quote(ss_buc_ancova(F_observed = 6, N = 7, cells = 3, n_covariates = 3)),
         "no error degrees of freedom remain once the covariates are estimated"),
    # split plot general: impossible degrees of freedom
    list(quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 3,
                                          num_groups = 3, effect = "between_only")),
         "cannot exceed 'num_groups' - 1"),
    list(quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 3,
                                          num_groups = 3, effect = "between_within",
                                          df_num_within = 2)),
         "must be a multiple of 'df_num_within'"),
    list(quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 6,
                                          num_groups = 3, effect = "between_within",
                                          df_num_within = 2)),
         "must be a multiple of 'df_num_within'")
  )
  for (case in cases) expect_guard(case[[1]], case[[2]])

  expect_warning(
    ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                               num_groups = 3, effect = "between_only",
                               df_num_within = 2),
    "'df_num_within' is used only when", fixed = TRUE)
  # ...and the value it ignores does not change the plan.
  expect_identical(
    suppressWarnings(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                                                num_groups = 3, effect = "between_only",
                                                df_num_within = 2)),
    ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                               num_groups = 3, effect = "between_only"))
})

test_that("every planner's zero-NCP refusal fires, with its class", {
  # The help pages' try() calls, plus unbalanced (odd N) cases for the
  # planners that round the prior design two ways, so that the refusal is
  # reached with the two branches' ceilings different.
  refusals <- list(
    quote(ss_buc_independent_t(t_observed = 3, n = 20, assurance = .95)),
    quote(ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .95)),
    quote(ss_buc_independent_t(t_observed = 3, N = 41, assurance = .95)),
    quote(ss_buc_paired_t(t_observed = 3, N = 40, assurance = .95)),
    quote(ss_buc_one_sample_t(t_observed = 3, N = 40, assurance = .95)),
    quote(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5,
                         assurance = .95)),
    quote(ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 4, assurance = .95)),
    quote(ss_buc_one_way_anova(F_observed = 6, N = 13, levels_A = 4, assurance = .95)),
    quote(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                                 effect = "factor_B", assurance = .95)),
    quote(ss_buc_factorial_anova(F_observed = 5, N = 121, levels_A = 2, levels_B = 3,
                                 effect = "factor_B", assurance = .95)),
    quote(ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                                         df_numerator = 2, df_denominator = 114,
                                         assurance = .95)),
    quote(ss_buc_factorial_anova_general(F_observed = 6, N = 13, cells = 4,
                                         df_numerator = 3, df_denominator = 9,
                                         assurance = .95)),
    quote(ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2,
                        assurance = .95)),
    quote(ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3,
                          effect = "factor_B", assurance = .95)),
    quote(ss_buc_rm_anova_general(F_observed = 6.5, N = 80, df_numerator = 1,
                                  assurance = .95)),
    quote(ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
                             levels_within = 3, effect = "within", assurance = .95)),
    quote(ss_buc_mixed_anova(F_observed = 4, N = 61, levels_between = 2,
                             levels_within = 3, effect = "within", assurance = .6)),
    quote(ss_buc_mixed_anova(F_observed = 7, N = 61, levels_between = 2,
                             levels_within = 3, effect = "between", assurance = .8)),
    quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                                     num_groups = 3, effect = "between_only",
                                     assurance = .95)),
    quote(ss_buc_mixed_anova_general(F_observed = 4, N = 61, df_numerator = 2,
                                     num_groups = 2, effect = "within_only",
                                     assurance = .6)),
    quote(ss_buc_mixed_anova_general(F_observed = 6, N = 91, df_numerator = 2,
                                     num_groups = 3, effect = "between_within",
                                     df_num_within = 2, assurance = .95)),
    quote(ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, assurance = .95)),
    quote(ss_buc_R2(F_observed = 5, N = 150, p = 4, assurance = .99)),
    quote(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2, assurance = .95)),
    quote(ss_buc_correlation(r_observed = .27, N = 100, assurance = .95)),
    quote(ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3, assurance = .95)),
    quote(ss_buc_manova(F_observed = 8, N = 81, p_variables = 3, assurance = .999)),
    quote(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42,
                            df_restricted = 43, assurance = .99))
  )
  for (call in refusals)
    expect_guard(call, "noncentrality parameter is zero", class = "bucss_zero_ncp")
})

test_that("every nonsignificance refusal fires, with its class", {
  # 'aware' rows are the rounding band: significant at the prior study's own
  # degrees of freedom, not at the rounded-down design's (see the checks at
  # the end of this test). The rest are plainly nonsignificant.
  plain <- "nonsignificant based on your specified 'alpha_prior'"
  aware <- "significant at the prior study's own degrees of freedom"
  cases <- list(
    list(quote(ss_buc_independent_t(t_observed = 1, n = 20)), plain),
    list(quote(ss_buc_independent_t(t_observed = 2.11, n = c(5, 200))), aware),
    list(quote(ss_buc_paired_t(t_observed = 1, N = 40)), plain),
    list(quote(ss_buc_one_sample_t(t_observed = 1, N = 40)), plain),
    list(quote(ss_buc_welch_t(t_observed = 1, n_1 = 40, n_2 = 55, sd_ratio = 1.5)), plain),
    list(quote(ss_buc_one_way_anova(F_observed = 1, N = 120, levels_A = 4)), plain),
    list(quote(ss_buc_one_way_anova(F_observed = 3.7, N = 13, levels_A = 4)), plain),
    list(quote(ss_buc_one_way_anova(F_observed = 3.95, N = 13, levels_A = 4)), aware),
    list(quote(ss_buc_factorial_anova(F_observed = 1, N = 120, levels_A = 2, levels_B = 3)),
         plain),
    list(quote(ss_buc_factorial_anova(F_observed = 5, N = 13, levels_A = 2, levels_B = 2)),
         plain),
    list(quote(ss_buc_factorial_anova(F_observed = 5.2, N = 13, levels_A = 2, levels_B = 2)),
         aware),
    list(quote(ss_buc_factorial_anova_general(F_observed = 1, N = 120, cells = 6,
                                              df_numerator = 2, df_denominator = 114)),
         plain),
    list(quote(ss_buc_factorial_anova_general(F_observed = 3.7, N = 13, cells = 4,
                                              df_numerator = 3, df_denominator = 9)),
         plain),
    list(quote(ss_buc_factorial_anova_general(F_observed = 3.95, N = 13, cells = 4,
                                              df_numerator = 3, df_denominator = 9)),
         aware),
    list(quote(ss_buc_ancova(F_observed = 1, N = 120, cells = 3, n_covariates = 2)), plain),
    list(quote(ss_buc_rm_anova(F_observed = 1, N = 60, levels_A = 3)), plain),
    list(quote(ss_buc_rm_anova_general(F_observed = 1, N = 60, df_numerator = 2)), plain),
    list(quote(ss_buc_mixed_anova(F_observed = 1, N = 60, levels_between = 2,
                                  levels_within = 3)), plain),
    list(quote(ss_buc_mixed_anova(F_observed = 4.005, N = 61, levels_between = 2,
                                  levels_within = 3, effect = "between")), aware),
    list(quote(ss_buc_mixed_anova_general(F_observed = 1, N = 60, df_numerator = 1,
                                          num_groups = 2)), plain),
    list(quote(ss_buc_mixed_anova_general(F_observed = 4.005, N = 61, df_numerator = 1,
                                          num_groups = 2, effect = "between_only")), aware),
    list(quote(ss_buc_reg_coef(t_observed = 1, N = 150, p = 3)), plain),
    list(quote(ss_buc_R2(F_observed = 1, N = 150, p = 3)), plain),
    list(quote(ss_buc_reg_joint(F_observed = 1, N = 150, p = 4, p_joint = 2)), plain),
    list(quote(ss_buc_correlation(r_observed = .1, N = 50)), plain),
    list(quote(ss_buc_manova(F_observed = 1, N = 80, p_variables = 3)), plain),
    list(quote(ss_buc_chisq_diff(chisq_observed = 2, N = 250, df_full = 42,
                                 df_restricted = 43)), plain)
  )
  for (case in cases) expect_guard(case[[1]], case[[2]], class = "bucss_nonsignificant")

  # The band rows really are in the band (independent check 2026-09-28, from
  # the prior studies' own and rounded-down degrees of freedom):
  expect_lt(2 * pt(2.11, 203, lower.tail = FALSE), .05)      # own df 203
  expect_gt(2 * pt(2.11, 16, lower.tail = FALSE), .05)       # rounded down: 2 * 9 - 2
  expect_lt(2 * pt(2.11, 18, lower.tail = FALSE), .05)      # round-up df 18 passes
  expect_lt(pf(3.95, 3, 9, lower.tail = FALSE), .05)         # one-way own df (3, 9)
  expect_gt(pf(3.95, 3, 8, lower.tail = FALSE), .05)         # rounded down (3, 8)
  expect_gt(pf(3.7, 3, 9, lower.tail = FALSE), .05)          # plainly nonsignificant
  expect_lt(pf(3.7, 3, 12, lower.tail = FALSE), .05)         # but passes the round-up branch
  expect_lt(pf(5.2, 1, 9, lower.tail = FALSE), .05)          # 2 x 2, N = 13: own (1, 9)
  expect_gt(pf(5.2, 1, 8, lower.tail = FALSE), .05)
  expect_gt(pf(5, 1, 9, lower.tail = FALSE), .05)            # plainly nonsignificant
  expect_lt(pf(5, 1, 12, lower.tail = FALSE), .05)           # but passes the round-up branch
  expect_lt(pf(4.005, 1, 59, lower.tail = FALSE), .05)       # split plot, N = 61: own (1, 59)
  expect_gt(pf(4.005, 1, 58, lower.tail = FALSE), .05)       # rounded down (1, 58)
})

test_that("the plot method's guards fire", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 40)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  # attr<- keeps the class, so these objects still dispatch to the method.
  unknown <- plan
  attr(unknown, "design") <- "Not a design"
  expect_error(plot(unknown), "No plot is defined for the design", fixed = TRUE)

  no_ceiling <- plan
  attr(no_ceiling, "assurance_ceiling") <- NULL
  expect_error(plot(no_ceiling), "does not carry an assurance ceiling", fixed = TRUE)

  tampered <- plan
  class(tampered) <- "data.frame"
  tampered$value[tampered$term == "ncp_adjusted"] <- 99
  class(tampered) <- c("bucss_plan", "data.frame")
  expect_error(plot(tampered), "cannot be plotted", fixed = TRUE)

  expect_error(plot(plan, lower = NA), "'lower' must be a single finite numeric value",
               fixed = TRUE)

  # Every requested assurance is above this plan's ceiling (about .90).
  expect_warning(curve_values <- plot(plan, assurance = c(.96, .97)),
                 "nothing was drawn", fixed = TRUE)
  expect_true(all(is.na(curve_values$ncp_adjusted)))
})

test_that("the sensitivity audit's guards fire", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)

  expect_error(BUCSS:::.draw_prior_statistic(list(family = "nonesuch"), 10, 1, 5),
               "Unknown sampling family", fixed = TRUE)
  expect_error(ss_buc_sensitivity(plan, true_ncp = 2, replications = 10, seed = 1.5),
               "'seed' must be a single whole number", fixed = TRUE)

  unknown <- plan
  attr(unknown, "design") <- "Not a design"
  expect_error(ss_buc_sensitivity(unknown, true_ncp = 2, replications = 10),
               "No simulation is defined for the design", fixed = TRUE)

  tampered <- plan
  class(tampered) <- "data.frame"
  tampered$value[tampered$term == "ncp_adjusted"] <- 99
  class(tampered) <- c("bucss_plan", "data.frame")
  expect_error(ss_buc_sensitivity(tampered, true_ncp = 2, replications = 10),
               "cannot be simulated", fixed = TRUE)

  # A true parameter beyond anything the planner can return (its cap is 1e7).
  expect_error(
    suppressWarnings(ss_buc_sensitivity(plan, true_ncp = 1e8, replications = 2, seed = 2017)),
    "Could not determine the sample size that 'true_ncp' requires", fixed = TRUE)

  # A true parameter so small that the planner refuses at it.
  expect_error(ss_buc_sensitivity(plan, true_ncp = 1e-6, replications = 2, seed = 2017),
               "At 'true_ncp' = 1e-06 the planner refuses to plan", fixed = TRUE)

  # A literature that publishes only at p < .0001 publishes almost nothing
  # from a small true effect (rate about .0003), so the draw budget runs out.
  strict <- ss_buc_paired_t(t_observed = 6, N = 30, alpha_prior = .0001)
  expect_error(ss_buc_sensitivity(strict, true_ncp = .5, replications = 2, seed = 2017),
               "Fewer than 'replications' simulated prior studies reached significance",
               fixed = TRUE)

  # Every simulated prior study refused: a weak truth and a demanding
  # assurance. The rates that are undefined come back NA, with a warning.
  demanding <- ss_buc_paired_t(t_observed = 5, N = 30, assurance = .95)
  expect_warning(
    s <- ss_buc_sensitivity(demanding, true_ncp = .3, replications = 2, seed = 2017),
    "Every simulated prior study was refused", fixed = TRUE)
  expect_identical(s$value[s$term == "refusal_rate"], 1)
  expect_true(is.na(s$value[s$term == "attainment"]))
  expect_true(is.na(s$value[s$term == "attainment_mcse"]))
  expect_identical(s$value[s$term == "attainment_with_refusals"], 1)
})

test_that("the correlation series refuses a window it cannot sum", {
  # rho^2 = 1 - 5e-7 with n = 5 puts the negative binomial's upper quantile
  # near 7e7 terms, beyond the 5e7 the series will evaluate. (Reaching this
  # through ss_buc_correlation() takes about 25 s of near-limit sums first,
  # so the helper is called directly.)
  expect_error(BUCSS:::.p_correlation_f(1, 5, 1 - 5e-7),
               "imply a noncentrality parameter too large to evaluate", fixed = TRUE)
})

test_that("an unexpected planner error stops the audit rather than counting as a refusal", {
  # Only the classed refusals count as refusals (N-07). Any other error during
  # the replications is reported with its message. The mock lets the first
  # call (the replay that checks the stored inputs) through and fails after.
  # local_mocked_bindings() needs testthat 3.1.7, the floor in Suggests.
  real_planner <- ss_buc_paired_t
  state <- new.env(parent = emptyenv())
  state$calls <- 0
  local_mocked_bindings(
    ss_buc_paired_t = function(...) {
      state$calls <- state$calls + 1
      if (state$calls > 1) stop("an unexpected failure", call. = FALSE)
      real_planner(...)
    },
    .package = "BUCSS")
  plan <- real_planner(t_observed = 3, N = 30)
  expect_error(ss_buc_sensitivity(plan, true_ncp = 2, replications = 5, seed = 2017),
               "While replaying the planner: an unexpected failure", fixed = TRUE)
})

test_that("the split-plot planners' round-up nonsignificance checks fire", {
  # Twenty groups, N = 41: 2 and 3 per group, so the between effect's df are
  # (19, 20) rounded down and (19, 40) rounded up. With 'alpha_prior' of .8 or
  # 1 (read as .999) the critical F is a lower quantile, which rises with the
  # error df, so each statistic below passes the round-down check and falls
  # short at the round-up one (independent check 2026-09-29, from qf() alone).
  expect_gt(.24, qf(1 - .999, 19, 20))
  expect_lte(.24, qf(1 - .999, 19, 40))
  expect_gt(.69, qf(1 - .8, 19, 20))
  expect_lte(.69, qf(1 - .8, 19, 40))
  cases <- list(
    quote(ss_buc_mixed_anova(F_observed = .24, N = 41, levels_between = 20,
                             levels_within = 2, effect = "between", alpha_prior = 1)),
    quote(ss_buc_mixed_anova(F_observed = .69, N = 41, levels_between = 20,
                             levels_within = 2, effect = "between", alpha_prior = .8)),
    quote(ss_buc_mixed_anova_general(F_observed = .24, N = 41, df_numerator = 19,
                                     num_groups = 20, effect = "between_only",
                                     alpha_prior = 1)))
  for (call in cases)
    expect_error(eval(call), class = "bucss_nonsignificant",
                 label = paste(deparse(call, width.cutoff = 500L), collapse = " "))
})
