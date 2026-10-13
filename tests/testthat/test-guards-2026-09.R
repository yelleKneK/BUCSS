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
# below 1 (.8, .999) it is a lower quantile, and with several numerator
# degrees of freedom a lower quantile rises. (At exactly 1, no publication
# filter, the critical value is 0 at every degrees of freedom, and only a
# statistic of 0 is refused.) Their own test follows the audit's blocks.

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
         "'effect' should be one of"),
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
         "'effect' should be one of"),
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
    # an odd N passes the check above (5 - 3 - 1 = 1) and reaches the one for
    # the rounded-down balanced design, where 2 * floor(5 / 2) - 3 - 1 = 0
    list(quote(ss_buc_manova(F_observed = 40, N = 5, p_variables = 3)),
         "leaves no error degrees of freedom once an odd 'N' is rounded down"),
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
  # .999 the critical F is a lower quantile, which rises with the error df, so
  # each statistic below passes the round-down check and falls short at the
  # round-up one (independent check 2026-09-29, from qf() alone). The .999
  # cases were written as 'alpha_prior = 1' while 1 was read as .999; at
  # exactly 1 the critical value is 0 in both designs, the checks cannot
  # differ, and the same statistic reaches the zero noncentrality refusal
  # instead (its largest supportable assurance is 1 - p, far below .80).
  expect_gt(.24, qf(1 - .999, 19, 20))
  expect_lte(.24, qf(1 - .999, 19, 40))
  expect_gt(.69, qf(1 - .8, 19, 20))
  expect_lte(.69, qf(1 - .8, 19, 40))
  cases <- list(
    quote(ss_buc_mixed_anova(F_observed = .24, N = 41, levels_between = 20,
                             levels_within = 2, effect = "between", alpha_prior = .999)),
    quote(ss_buc_mixed_anova(F_observed = .69, N = 41, levels_between = 20,
                             levels_within = 2, effect = "between", alpha_prior = .8)),
    quote(ss_buc_mixed_anova_general(F_observed = .24, N = 41, df_numerator = 19,
                                     num_groups = 20, effect = "between_only",
                                     alpha_prior = .999)))
  for (call in cases)
    expect_error(eval(call), class = "bucss_nonsignificant",
                 label = paste(deparse(call, width.cutoff = 500L), collapse = " "))
  stopifnot(qf(1 - 1, 19, 20) == 0, qf(1 - 1, 19, 40) == 0)
  for (call in cases[c(1, 3)]) {
    call$alpha_prior <- 1
    expect_error(eval(call), class = "bucss_zero_ncp",
                 label = paste(deparse(call, width.cutoff = 500L), collapse = " "))
  }
})

# Guards added 2026-10-12, from the triage of the unfinished 2026-09-29
# review. Each block names its finding and fails without the guard it covers.

# Evaluate 'code', returning the error it raises (NULL if none) and the text
# of every warning raised on the way.
conditions_of <- function(code) {
  record <- new.env(parent = emptyenv())
  record$warnings <- character(0)
  error <- tryCatch(
    withCallingHandlers({
      code
      NULL
    }, warning = function(w) {
      record$warnings <- c(record$warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) e)
  list(error = error, warnings = record$warnings)
}

test_that("GF-1: a MANOVA prior with no error df once balanced stops cleanly, with no stray warning", {
  # N = 5 with three outcomes has 5 - 3 - 1 = 1 error df, but the rounded-down
  # balanced design, 2 * floor(5 / 2) cases, has 4 - 3 - 1 = 0. A statistic
  # short of the round-up branch's critical value used to be called
  # nonsignificant, with R's "NaNs produced" beside the error (the message
  # evaluated pf() at 0 df); a larger one reached the df message alone.
  stopifnot(2 * floor(5 / 2) - 3 - 1 == 0)
  no_df <- "leaves no error degrees of freedom once an odd 'N' is rounded down to two equal groups"
  cases <- list(
    quote(ss_buc_manova(F_observed = 5, N = 5, p_variables = 3)),
    quote(ss_buc_manova(F_observed = 5, N = 7, p_variables = 5)),
    quote(ss_buc_manova(F_observed = 5, N = 9, p_variables = 7)),
    quote(ss_buc_manova(T2_observed = 50, N = 5, p_variables = 3)),
    quote(ss_buc_manova(F_observed = 40, N = 5, p_variables = 3)))
  for (call in cases) {
    label <- paste(deparse(call, width.cutoff = 500L), collapse = " ")
    seen <- conditions_of(eval(call))
    expect_s3_class(seen$error, "error")
    expect_false(inherits(seen$error, "bucss_nonsignificant"), label = label)
    expect_match(conditionMessage(seen$error), no_df, fixed = TRUE, label = label)
    expect_identical(seen$warnings, character(0), label = label)
  }
  # One error df in the rounded-down design is enough (N = 7 with four
  # outcomes: 6 - 4 - 1 = 1), so the guard does not reach a design that plans.
  expect_s3_class(ss_buc_manova(F_observed = 500, N = 7, p_variables = 4, alpha_prior = 1,
                                assurance = .5), "bucss_plan")
})

test_that("GF-3: a malformed 'effect' gets one message naming the argument and its choices", {
  # NULL used to end in "missing value where TRUE/FALSE needed" (or, in the
  # general split-plot planner, to plan the first choice); a number or NA in
  # match.arg()'s "'arg' must be NULL or a character vector"; and a vector of
  # two effects planned the first with no notice.
  planners <- list(
    list(choices = c("factor_A", "factor_B", "interaction"),
         run = function(e) ss_buc_factorial_anova(F_observed = 12, N = 120, levels_A = 2,
                                                  levels_B = 3, effect = e)),
    list(choices = c("factor_A", "factor_B", "interaction"),
         run = function(e) ss_buc_rm_anova(F_observed = 12, N = 60, levels_A = 2,
                                           levels_B = 3, effect = e)),
    list(choices = c("between", "within", "interaction"),
         run = function(e) ss_buc_mixed_anova(F_observed = 12, N = 60, levels_between = 2,
                                              levels_within = 3, effect = e)),
    list(choices = c("between_only", "within_only", "between_within"),
         run = function(e) ss_buc_mixed_anova_general(F_observed = 12, N = 90,
                                                      df_numerator = 2, num_groups = 3,
                                                      effect = e)))
  for (p in planners) {
    msg <- paste0("'effect' should be one of \"", p$choices[1], "\", \"", p$choices[2],
                  "\", or \"", p$choices[3], "\", given as a single character string.")
    malformed <- list(NULL, 2, NA, NA_character_, TRUE, character(0), "bogus", "",
                      p$choices[2:1], rev(p$choices), list("bogus"))
    for (bad in malformed)
      expect_error(p$run(bad), msg, fixed = TRUE,
                   label = paste(p$choices[1], paste(deparse(bad), collapse = " ")))
    # The untouched default, the whole vector of choices, still means the
    # first; a single value, an unambiguous prefix, and a one-element factor
    # (a row of expand.grid(), say) are read as before.
    expect_identical(p$run(p$choices), p$run(p$choices[1]))
    expect_identical(attr(p$run(p$choices), "effect"), p$choices[1])
    expect_identical(attr(p$run(p$choices[2]), "effect"), p$choices[2])
    expect_identical(p$run(factor(p$choices[2])), p$run(p$choices[2]))
  }
  expect_identical(
    attr(ss_buc_factorial_anova(F_observed = 12, N = 120, levels_A = 2, levels_B = 3,
                                effect = "inter"), "effect"), "interaction")
})

test_that("GF-5: a spread too extreme to compute with is refused in the terms the user gave", {
  # A ratio near 1e77 or more overflows the squares in the Welch-Satterthwaite
  # degrees of freedom and used to end in "missing value where TRUE/FALSE
  # needed"; a pair of standard deviations or variances whose ratio overflows
  # or underflows was reported as a bad 'sd_ratio', which the user had not
  # given.
  w <- function(...) ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, ...)
  too_far <- paste0("The two groups' spread differs by too many orders of magnitude to plan ",
                    "from: the Welch-Satterthwaite degrees of freedom cannot be computed. ",
                    "Please check ")
  for (r in c(1e78, 1e100, 1e150, 1e200))
    expect_error(w(sd_ratio = r), paste0(too_far, "'sd_ratio'."), fixed = TRUE)
  pairs <- list(list(sd_1 = 1e-200, sd_2 = 1e200),      # the ratio is Inf
                list(sd_1 = 1e200, sd_2 = 1e-200),      # the ratio is 0
                list(sd_1 = 1e-100, sd_2 = 1))          # a finite ratio of 1e100
  for (p in pairs) {
    msg <- tryCatch(do.call(w, p), error = conditionMessage)
    expect_identical(msg, paste0(too_far, "'sd_1' and 'sd_2'."))
  }
  msg <- tryCatch(w(var_1 = 1e-320, var_2 = 1), error = conditionMessage)
  expect_identical(msg, paste0(too_far, "'var_1' and 'var_2'."))

  # Extreme ratios that can be computed still plan, at the limit they
  # approach: with one group's variance negligible the test is a one-sample t
  # on the other group, so the plan is the paired t planner's from that
  # group's size (df n - 1, and the noncentrality scaling as n).
  one_sample <- function(N) ss_buc_paired_t(t_observed = 3, N = N)
  computed <- c("actual_power", "ncp_adjusted")
  for (case in list(list(ratio = 1e-150, N = 40), list(ratio = 1e-200, N = 40),
                    list(ratio = 1e77, N = 55))) {
    plan <- w(sd_ratio = case$ratio)
    limit <- one_sample(case$N)
    expect_identical(plan$value[plan$term == "necessary_n_per_group"],
                     limit$value[limit$term == "necessary_N"])
    expect_equal(plan$value[plan$term %in% computed], limit$value[limit$term %in% computed],
                 tolerance = 1e-12)
  }
})

test_that("GF-6: 'n = NULL' with 'N' is 'N' alone, and two group sizes get their own message", {
  # 'N = NULL' with 'n' always planned; 'n = NULL' with 'N' stopped with
  # "Because you specified 'N' you should not specify 'n'."
  expect_identical(ss_buc_independent_t(t_observed = 3, n = NULL, N = 40),
                   ss_buc_independent_t(t_observed = 3, N = 40))
  expect_identical(ss_buc_independent_t(t_observed = 3, N = NULL, n = 20),
                   ss_buc_independent_t(t_observed = 3, n = 20))
  expect_identical(ss_buc_independent_t(t_observed = 3, n = NULL, N = 41),
                   ss_buc_independent_t(t_observed = 3, N = 41))
  # Neither size, however it is spelled, and both sizes are still refused.
  expect_error(ss_buc_independent_t(t_observed = 3, n = NULL), "either 'n'", fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = 3, n = NULL, N = NULL), "either 'n'",
               fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, N = 40),
               "should not specify 'n'", fixed = TRUE)
  # Two sizes were told "'n' must be a single whole number.", which 'n' need
  # not be.
  two <- "Each of the two group sizes in 'n' must be a whole number of at least 2."
  for (bad in list(c(20, NA), c(20, 2.5), c(1, 20), c(20, Inf), list(20, "a")))
    expect_error(ss_buc_independent_t(t_observed = 3, n = bad), two, fixed = TRUE,
                 label = paste(deparse(bad), collapse = " "))
  # A single size keeps the count check's own messages.
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20.5),
               "'n' must be a single whole number.", fixed = TRUE)
  expect_error(ss_buc_independent_t(t_observed = 3, n = 1), "'n' must be at least 2.",
               fixed = TRUE)
})

test_that("GF-7: an alpha_planned too small to compute with is refused by the validator", {
  # At 1e-17 the confidence level 1 - alpha_planned is exactly 1, the critical
  # value is infinite, and the planned power is 0 or not a number at every
  # sample size: the random-predictor planners ended in "missing value where
  # TRUE/FALSE needed", and the others in the 2^52 stop, which blamed the
  # assurance.
  stopifnot(1 - 1e-17 == 1)
  small <- paste0("'alpha_planned' is too small to plan with: 1 - 'alpha_planned' cannot be ",
                  "told apart from 1 in double precision arithmetic")
  cases <- list(
    quote(ss_buc_correlation(r_observed = .35, N = 100, alpha_planned = 1e-17)),
    quote(ss_buc_R2(F_observed = 5, N = 150, p = 4, alpha_planned = 1e-17,
                    random_predictors = TRUE)),
    quote(ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, alpha_planned = 1e-17,
                          random_predictors = TRUE)),
    quote(ss_buc_R2(F_observed = 5, N = 150, p = 4, alpha_planned = 1e-17)),
    quote(ss_buc_paired_t(t_observed = 3, N = 40, alpha_planned = 1e-17)),
    quote(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42, df_restricted = 43,
                            alpha_planned = 1e-300)))
  for (call in cases) expect_guard(call, small)
  # The smallest values for which 1 - alpha_planned is below 1 are untouched:
  # the correlation planner still plans there.
  stopifnot(1 - 1e-16 < 1)
  expect_s3_class(ss_buc_correlation(r_observed = .35, N = 100, alpha_planned = 1e-16),
                  "bucss_plan")
})

test_that("ENG-3 and CL-33: the 2^52 stop names every cause, and a power that is not a number has its own message", {
  # A prior N of 4e15 at the default assurance of .80, nowhere near this
  # prior's ceiling (about .95): the stop used to say only that the assurance
  # was essentially at its ceiling and should be lowered. The same prior at
  # N = 1e15 plans, so the sample size is the cause.
  msg <- tryCatch(ss_buc_paired_t(t_observed = 3, N = 4e15), error = conditionMessage)
  expect_match(msg, "The planned sample size would exceed about 4.5e15 (2^52). Either",
               fixed = TRUE)
  expect_match(msg, paste0("the requested 'assurance' is essentially at the largest value ",
                           "the prior result can support (lower 'assurance' a little)"),
               fixed = TRUE)
  expect_match(msg, "or the prior study's sample size is implausibly large", fixed = TRUE)
  expect_match(msg, "or 'alpha_planned' is too small to compute with", fixed = TRUE)
  expect_s3_class(ss_buc_paired_t(t_observed = 3, N = 1e15), "bucss_plan")
  # The third cause: for a t planner the two-sided critical value is already
  # infinite a little above what the validator refuses.
  stopifnot(1 - 1e-16 < 1, 1 - 1e-16 / 2 == 1)
  expect_error(ss_buc_paired_t(t_observed = 3, N = 40, alpha_planned = 1e-16),
               "or 'alpha_planned' is too small to compute with", fixed = TRUE)

  # A power that is not a finite number used to end the search with base R's
  # "missing value where TRUE/FALSE needed", wherever it was met: at the
  # start, while the bracket doubles, or while it is bisected.
  search <- BUCSS:::.smallest_n_for_power
  not_a_number <- "The planned study's power could not be computed (it is not a finite number)"
  expect_error(search(function(n) NaN, .5), not_a_number, fixed = TRUE)
  expect_error(search(function(n) if (n > 8) NA_real_ else 0, .5), not_a_number, fixed = TRUE)
  expect_error(search(function(n) if (n == 6) NaN else as.numeric(n >= 8), .5), not_a_number,
               fixed = TRUE)
  # Reached through a planner (the nuisance degrees of freedom of this
  # impossible prior study overflow the planned study's).
  expect_guard(quote(ss_buc_factorial_anova_general(F_observed = 12, N = 1e17, cells = 2,
                                                    df_numerator = 1, df_denominator = 1000)),
               not_a_number)
})

test_that("ENG-9: a user already at alpha_prior = 1 is not told to increase it", {
  # 'alpha_prior = 1' is computed as exactly 1: the critical value is 0, so
  # only a statistic of exactly 0 is refused as nonsignificant there, and
  # "Please increase 'alpha_prior'" could not be followed. One case per
  # planner; the premise is from qt(), qf(), and qchisq() alone. (While 1 was
  # read as .999 the cases were statistics with p above .999, such as a t of
  # .0005; those now pass this check and are refused for a zero noncentrality
  # parameter instead, as the second half of this block pins.)
  stopifnot(qt(1 - 1 / 2, 19) == 0, qf(1 - 1, 1, 1e6) == 0, qf(1 - 1, 2, 57) == 0,
            qchisq(1 - 1, 1) == 0)
  too_small <- paste0(" An 'alpha_prior' of 1 already assumes no publication filter, so it ",
                      "cannot be increased: ")
  ending <- " is too small to plan from even with no publication filter."
  planners <- list(
    list(quote(ss_buc_independent_t(n = 20, alpha_prior = 1)), "t_observed", "'t_observed'"),
    list(quote(ss_buc_paired_t(N = 30, alpha_prior = 1)), "t_observed", "'t_observed'"),
    list(quote(ss_buc_one_sample_t(N = 30, alpha_prior = 1)), "t_observed", "'t_observed'"),
    list(quote(ss_buc_welch_t(n_1 = 20, n_2 = 30, alpha_prior = 1)), "t_observed",
         "'t_observed'"),
    list(quote(ss_buc_one_way_anova(N = 41, levels_A = 4, alpha_prior = 1)), "F_observed",
         "'F_observed'"),
    list(quote(ss_buc_factorial_anova(N = 120, levels_A = 2, levels_B = 3, alpha_prior = 1)),
         "F_observed", "'F_observed'"),
    list(quote(ss_buc_factorial_anova_general(N = 120, cells = 6, df_numerator = 2,
                                              df_denominator = 114, alpha_prior = 1)),
         "F_observed", "'F_observed'"),
    list(quote(ss_buc_ancova(N = 120, cells = 3, n_covariates = 2, alpha_prior = 1)),
         "F_observed", "'F_observed'"),
    list(quote(ss_buc_rm_anova(N = 60, levels_A = 3, alpha_prior = 1)), "F_observed",
         "'F_observed'"),
    list(quote(ss_buc_rm_anova_general(N = 60, df_numerator = 2, alpha_prior = 1)),
         "F_observed", "'F_observed'"),
    list(quote(ss_buc_mixed_anova(N = 60, levels_between = 2, levels_within = 3,
                                  alpha_prior = 1)), "F_observed", "'F_observed'"),
    list(quote(ss_buc_mixed_anova_general(N = 60, df_numerator = 1, num_groups = 2,
                                          alpha_prior = 1)), "F_observed", "'F_observed'"),
    list(quote(ss_buc_reg_coef(N = 150, p = 3, alpha_prior = 1)), "t_observed",
         "'t_observed'"),
    list(quote(ss_buc_R2(N = 150, p = 3, alpha_prior = 1)), "F_observed", "'F_observed'"),
    list(quote(ss_buc_reg_joint(N = 150, p = 4, p_joint = 2, alpha_prior = 1)), "F_observed",
         "'F_observed'"),
    list(quote(ss_buc_correlation(N = 100, alpha_prior = 1)), "r_observed",
         "the prior result"),
    list(quote(ss_buc_manova(N = 80, p_variables = 3, alpha_prior = 1)), "F_observed",
         "the prior result"),
    list(quote(ss_buc_chisq_diff(N = 250, df_full = 42, df_restricted = 43, alpha_prior = 1)),
         "chisq_observed", "the prior result"))
  expect_length(planners, 18L)
  # The statistics the .999 reading refused as nonsignificant, by argument.
  tiny <- c(t_observed = .0005, F_observed = 1e-9, r_observed = 1e-6, chisq_observed = 1e-9)
  for (case in planners) {
    # A statistic of exactly 0: nonsignificant even with no publication filter.
    call <- case[[1]]
    call[[case[[2]]]] <- 0
    label <- paste(deparse(call, width.cutoff = 500L), collapse = " ")
    err <- tryCatch(eval(call), error = function(e) e)
    expect_s3_class(err, "bucss_nonsignificant")
    expect_match(conditionMessage(err), paste0(too_small, case[[3]], ending), fixed = TRUE,
                 label = label)
    expect_match(conditionMessage(err), "is nonsignificant based on your specified 'alpha_prior'",
                 fixed = TRUE, label = label)
    expect_no_match(conditionMessage(err), "Please increase", fixed = TRUE, label = label)
    # The string spelling of 1 is refused the same way.
    call$alpha_prior <- "No Publication Filter"
    expect_identical(conditionMessage(tryCatch(eval(call), error = function(e) e)),
                     conditionMessage(err), label = label)
    # A statistic barely above 0 is not nonsignificant at 1: it is too weak to
    # support the assurance (its ceiling is 1 - p, below .001 here).
    call <- case[[1]]
    call[[case[[2]]]] <- tiny[[case[[2]]]]
    label <- paste(deparse(call, width.cutoff = 500L), collapse = " ")
    err <- tryCatch(eval(call), error = function(e) e)
    expect_s3_class(err, "bucss_zero_ncp")
    expect_match(conditionMessage(err), "the prior result is too weak to plan from",
                 fixed = TRUE, label = label)
  }
  # Below 1 the advice can be followed, and is still given: at .999, the
  # largest value the old reading used, as at any smaller value.
  expect_error(ss_buc_paired_t(t_observed = .0005, N = 30, alpha_prior = .999),
               "Please increase 'alpha_prior' so 't_observed' exceeds the critical value.",
               fixed = TRUE, class = "bucss_nonsignificant")
  expect_error(ss_buc_paired_t(t_observed = 1, N = 40),
               "Please increase 'alpha_prior' so 't_observed' exceeds the critical value.",
               fixed = TRUE, class = "bucss_nonsignificant")
  stopifnot(pf(.5, 3, 116, lower.tail = FALSE) > .5)
  expect_error(ss_buc_one_way_anova(F_observed = .5, N = 120, levels_A = 4, alpha_prior = .5),
               "Please increase 'alpha_prior' so 'F_observed' exceeds the critical value.",
               fixed = TRUE, class = "bucss_nonsignificant")
})
