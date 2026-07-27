# The designs added after the original twelve. Each is either an exact
# relabeling of an existing planner (verified by identity) or new plumbing
# (verified against a hand computation), plus the contract every planner
# shares.

test_that("the one-sample t planner is the paired t planner relabeled", {
  one <- ss_buc_one_sample_t(t_observed = 3, N = 40)
  pair <- ss_buc_paired_t(t_observed = 3, N = 40)
  expect_equal(one$value[one$term == "necessary_N"],
               pair$value[pair$term == "necessary_N"])
  expect_equal(one$value[one$term == "ncp_adjusted"],
               pair$value[pair$term == "ncp_adjusted"])
  expect_equal(one$value[one$term == "actual_power"],
               pair$value[pair$term == "actual_power"])
  expect_identical(attr(one, "design"), "One-sample t test")
  expect_identical(attr(one, "sample_size_unit"), "total")
  # pinned value (equal to the paired planner's documented example)
  expect_equal(one$value[one$term == "necessary_N"], 255)
})

test_that("the correlation planner matches its own t form and the regression planner", {
  r_form <- ss_buc_correlation(r_observed = .35, N = 100)
  t_form <- ss_buc_correlation(t_observed = .35 * sqrt(98 / (1 - .35^2)), N = 100)
  expect_equal(r_form$value[r_form$term == "necessary_N"],
               t_form$value[t_form$term == "necessary_N"])
  expect_equal(r_form$value[r_form$term == "ncp_adjusted"],
               t_form$value[t_form$term == "ncp_adjusted"])
  # a correlation is the slope in a simple regression, so the single-predictor
  # regression planner must agree
  reg <- ss_buc_reg_coef(t_observed = .35 * sqrt(98 / (1 - .35^2)), N = 100, p = 1)
  expect_equal(r_form$value[r_form$term == "necessary_N"],
               reg$value[reg$term == "necessary_N"])
  expect_equal(r_form$value[r_form$term == "ncp_adjusted"],
               reg$value[reg$term == "ncp_adjusted"])
  expect_equal(r_form$value[r_form$term == "necessary_N"], 119)
  # a correlation of 1 or more is not a correlation
  expect_error(ss_buc_correlation(r_observed = 1, N = 100), "between -1 and 1",
               fixed = TRUE)
  expect_error(ss_buc_correlation(N = 100), "either 'r_observed'", fixed = TRUE)
  expect_error(ss_buc_correlation(r_observed = .3, t_observed = 3, N = 100),
               "not both", fixed = TRUE)
})

test_that("the ANCOVA planner reduces to the ANOVA planner with no covariates", {
  none <- ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 0)
  anova <- ss_buc_one_way_anova(F_observed = 6, N = 120, levels_A = 3)
  expect_equal(none$value[none$term == "necessary_n_per_cell"],
               anova$value[anova$term == "necessary_n_per_group"])
  expect_equal(none$value[none$term == "ncp_adjusted"],
               anova$value[anova$term == "ncp_adjusted"])

  # covariates cost error degrees of freedom, so they call for a larger sample
  two <- ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2)
  expect_gte(two$value[two$term == "necessary_n_per_cell"],
             none$value[none$term == "necessary_n_per_cell"])
  # the planned study's error df carry the covariates through
  expect_equal(attr(two, "df_error"),
               two$value[two$term == "total_N"] - 3 - 2)
  expect_equal(two$value[two$term == "necessary_n_per_cell"], 154)
  # a single-degree-of-freedom contrast among the adjusted means
  contrast <- ss_buc_ancova(F_observed = 9, N = 120, cells = 3,
                            n_covariates = 2, df_numerator = 1)
  expect_equal(attr(contrast, "df_effect"), 1)
  expect_error(ss_buc_ancova(F_observed = 6, N = 120, cells = 3),
               "n_covariates", fixed = TRUE)
})

test_that("the multivariate planner agrees with itself through T squared", {
  f_form <- ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3)
  t2_form <- ss_buc_manova(T2_observed = 4.5 * 3 * 78 / 76, N = 80, p_variables = 3)
  # the design results agree; only the echoed input row differs (F versus T2)
  design_rows <- c("necessary_n_per_group", "total_N", "actual_power",
                   "ncp_adjusted")
  expect_equal(f_form$value[f_form$term %in% design_rows],
               t2_form$value[t2_form$term %in% design_rows])
  # one outcome variable is the equal-variance two-group t test: the same plan,
  # with the noncentrality on the F scale (lambda) rather than the t scale
  # (delta), so lambda equals delta squared
  uni <- ss_buc_manova(F_observed = 9, N = 80, p_variables = 1)
  tt <- ss_buc_independent_t(t_observed = 3, N = 80)
  expect_equal(uni$value[uni$term == "necessary_n_per_group"],
               tt$value[tt$term == "necessary_n_per_group"])
  expect_equal(uni$value[uni$term == "ncp_adjusted"],
               tt$value[tt$term == "ncp_adjusted"]^2)
  expect_equal(attr(f_form, "df_effect"), 3)
  expect_equal(attr(f_form, "df_error"),
               2 * f_form$value[f_form$term == "necessary_n_per_group"] - 3 - 1)
  expect_error(ss_buc_manova(F_observed = 4.5, N = 80), "p_variables",
               fixed = TRUE)
  expect_error(ss_buc_manova(F_observed = 4.5, N = 6, p_variables = 6),
               "degrees of freedom", fixed = TRUE)
})

test_that("the chi-square difference planner matches a hand computation", {
  res <- ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_difference = 1)
  ncp <- res$value[res$term == "ncp_adjusted"]
  size <- res$value[res$term == "necessary_N"]
  # the noncentrality is linear in N, and the test has no denominator df
  crit <- qchisq(.95, df = 1)
  expect_equal(res$value[res$term == "actual_power"],
               1 - pchisq(crit, df = 1, ncp = (size / 250) * ncp))
  expect_lt(1 - pchisq(crit, df = 1, ncp = ((size - 1) / 250) * ncp), .80)
  expect_null(attr(res, "df_error"))
  expect_equal(attr(res, "df_effect"), 1)
  # the ceiling is the same closed form as everywhere else
  p_prior <- 1 - pchisq(9.5, df = 1)
  expect_equal(attr(res, "assurance_ceiling"), 1 - p_prior / .05)
  expect_error(ss_buc_chisq_diff(chisq_observed = 1, N = 250, df_difference = 1),
               "nonsignificant", fixed = TRUE)
})

test_that("the Welch planner reduces sensibly and responds to the variance ratio", {
  equal_sd <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1)
  pooled <- ss_buc_independent_t(t_observed = 3, n = c(40, 55))
  # Welch degrees of freedom are smaller than pooled, so the Welch plan is a
  # little more conservative, but the two must be close
  expect_gte(equal_sd$value[equal_sd$term == "necessary_n_per_group"],
             pooled$value[pooled$term == "necessary_n_per_group"] - 5)
  expect_lte(equal_sd$value[equal_sd$term == "necessary_n_per_group"],
             pooled$value[pooled$term == "necessary_n_per_group"] + 15)
  # a more discrepant variance ratio needs more participants
  ratio_2 <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 2)
  expect_gt(ratio_2$value[ratio_2$term == "necessary_n_per_group"],
            equal_sd$value[equal_sd$term == "necessary_n_per_group"])
  expect_error(ss_buc_welch_t(t_observed = 3, n_1 = 40), "n_1", fixed = TRUE)
  expect_error(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 0),
               "positive", fixed = TRUE)
})

test_that("every new planner honors the shared result contract", {
  results <- list(
    ss_buc_correlation(r_observed = .35, N = 100),
    ss_buc_one_sample_t(t_observed = 3, N = 40),
    ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2),
    ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3),
    ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_difference = 1),
    ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5)
  )
  for (res in results) {
    expect_s3_class(res, "bucss_power")
    expect_type(res$value, "double")
    expect_true(res$term[1] %in% BUCSS:::.BUCSS_SIZE_TERMS)
    expect_true(all(c("actual_power", "ncp_adjusted", "alpha_prior",
                      "alpha_planned", "assurance", "desired_power") %in% res$term))
    expect_null(attr(res, "inputs"))
    expect_false(is.null(attr(res, "design")))
    expect_false(is.null(attr(res, "sample_size_unit")))
    expect_false(is.null(attr(res, "assurance_ceiling")))
    # the power of the reported plan is never below what was asked for
    expect_gte(res$value[res$term == "actual_power"],
               res$value[res$term == "desired_power"])
    # the shared verbs and display all work
    expect_identical(nrow(tidy(res)), 1L)
    expect_identical(nrow(glance(res)), 1L)
    expect_type(planning_sentence(res), "character")
    expect_output(print(res), "Design: ", fixed = TRUE)
    expect_identical(res[[1]], res$value[res$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1])
  }
})

test_that("the new planners reject an assurance above their ceiling", {
  # this prior supports assurance up to about .993
  expect_error(ss_buc_correlation(r_observed = .35, N = 100, assurance = .999),
               "noncentrality parameter is zero", fixed = TRUE)
  expect_error(ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3,
                             assurance = .95),
               "noncentrality parameter is zero", fixed = TRUE)
  expect_error(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250,
                                 df_difference = 1, assurance = .99),
               "noncentrality parameter is zero", fixed = TRUE)
  expect_error(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55,
                              sd_ratio = 1.5, assurance = .95),
               "noncentrality parameter is zero", fixed = TRUE)
})
