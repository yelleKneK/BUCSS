# ss_buc_sensitivity(): the Monte Carlo companion to the planners.

# One plan per design, used by the table checks below. Every design in
# .BUCSS_DESIGN_SPECS must appear here.
sensitivity_plans <- list(
  ss_buc_independent_t(t_observed = 3, n = c(50, 50), assurance = .80),
  ss_buc_independent_t(t_observed = 3, N = 100, assurance = .80),
  ss_buc_paired_t(t_observed = 3, N = 30, assurance = .80),
  ss_buc_one_sample_t(t_observed = 3, N = 30, assurance = .80),
  ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 40, sd_ratio = 1.5),
  ss_buc_one_way_anova(F_observed = 8, N = 60, levels_A = 3),
  ss_buc_factorial_anova(F_observed = 12, N = 60, levels_A = 3, levels_B = 2,
                         effect = "factor_A"),
  ss_buc_factorial_anova(F_observed = 12, N = 60, levels_A = 3, levels_B = 2,
                         effect = "interaction"),
  ss_buc_factorial_anova_general(F_observed = 12, N = 60, cells = 6,
                                 df_numerator = 2, df_denominator = 54),
  ss_buc_ancova(F_observed = 12, N = 60, cells = 6, df_numerator = 2,
                n_covariates = 2),
  ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3),
  ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2),
  ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                     levels_within = 3, effect = "between"),
  ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                     levels_within = 3, effect = "within"),
  ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                     levels_within = 3, effect = "interaction"),
  ss_buc_mixed_anova_general(F_observed = 20, N = 60, num_groups = 2,
                             df_numerator = 1, effect = "between_only"),
  ss_buc_mixed_anova_general(F_observed = 20, N = 60, num_groups = 2,
                             df_numerator = 2, effect = "within_only"),
  ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3),
  ss_buc_R2(F_observed = 8, N = 60, p = 3),
  ss_buc_reg_joint(F_observed = 12, N = 60, p = 5, p_joint = 2),
  ss_buc_manova(F_observed = 8, N = 40, p_variables = 3),
  ss_buc_chisq_diff(chisq_observed = 15, df_difference = 2, N = 200),
  ss_buc_correlation(r_observed = .35, N = 100),
  ss_buc_correlation(t_observed = 3.7, N = 100)
)

test_that("every design has a simulation spec, and no spec is orphaned", {
  covered <- unique(vapply(sensitivity_plans, attr, character(1), "design"))
  expect_setequal(covered, names(BUCSS:::.BUCSS_DESIGN_SPECS))
  # each spec names a real exported planner
  for (spec in BUCSS:::.BUCSS_DESIGN_SPECS)
    expect_true(is.function(get(spec$planner, envir = asNamespace("BUCSS"))))
})

test_that("each spec's prior distribution reproduces the plan's own ceiling", {
  # The assurance ceiling a planner stores is 1 - p/alpha_prior, where p is the
  # prior study's p value. Computing that p independently from the spec's
  # family and degrees of freedom pins every entry of the design table: a wrong
  # family, or degrees of freedom off by one, moves it immediately.
  for (o in sensitivity_plans) {
    spec <- BUCSS:::.BUCSS_DESIGN_SPECS[[attr(o, "design")]]
    inputs <- BUCSS:::.bucss_design_inputs(o)
    effect <- attr(o, "effect")
    df <- if (length(formals(spec$df)) > 1L) spec$df(inputs, effect) else
      spec$df(inputs)
    args <- BUCSS:::.bucss_design_args(inputs, spec, effect)
    stat <- args[[spec$statistic]]
    alpha_prior <- if (inputs$alpha_prior == 1) .999 else inputs$alpha_prior
    p_spec <- switch(
      spec$family,
      t = 2 * (1 - pt(abs(stat), df)),
      f = 1 - pf(stat, df[1], df[2]),
      f_from_t = 1 - pf(stat^2, 1, df[2]),
      correlation = 1 - pf(stat^2, 1, df - 2),
      chisq = 1 - pchisq(stat, df))
    expect_equal(p_spec, alpha_prior * (1 - attr(o, "assurance_ceiling")),
                 tolerance = 1e-9,
                 info = paste("design:", attr(o, "design"),
                              "effect:", if (is.null(effect)) "-" else effect))
  }
})

test_that("the rebuilt planner call reproduces every plan exactly", {
  computed <- c(BUCSS:::.BUCSS_SIZE_TERMS, "total_N", "actual_power",
                "ncp_adjusted")
  for (o in sensitivity_plans) {
    spec <- BUCSS:::.BUCSS_DESIGN_SPECS[[attr(o, "design")]]
    args <- BUCSS:::.bucss_design_args(BUCSS:::.bucss_design_inputs(o), spec,
                                      attr(o, "effect"))
    replay <- do.call(get(spec$planner, envir = asNamespace("BUCSS")), args)
    expect_equal(replay$value[replay$term %in% computed],
                 o$value[o$term %in% computed],
                 info = paste("design:", attr(o, "design")))
  }
})

test_that("ss_buc_sensitivity returns the documented shape", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30, assurance = .80)
  s <- ss_buc_sensitivity(plan, true_ncp = 2, replications = 200, seed = 2017)
  expect_s3_class(s, "bucss_sensitivity")
  expect_s3_class(s, "data.frame")
  expect_identical(names(s), c("term", "value"))
  expect_true(is.numeric(s$value))
  expect_identical(
    s$term,
    c("true_ncp", "ncp_adjusted", "replications", "publication_rate",
      "refusal_rate", "attainment", "attainment_mcse",
      "attainment_with_refusals", "size_at_true_effect", "size_q25",
      "size_median", "size_q75"))
  expect_equal(s$value[s$term == "true_ncp"], 2)
  expect_equal(s$value[s$term == "replications"], 200)
  expect_equal(s$value[s$term == "ncp_adjusted"],
               plan$value[plan$term == "ncp_adjusted"])
  expect_identical(attr(s, "design"), attr(plan, "design"))
  expect_identical(attr(s, "sample_size_unit"), attr(plan, "sample_size_unit"))
  expect_equal(attr(s, "assurance"), .80)
  # rates are proportions; the size rows are whole numbers
  for (t in c("publication_rate", "refusal_rate", "attainment",
              "attainment_with_refusals")) {
    expect_gte(s$value[s$term == t], 0)
    expect_lte(s$value[s$term == t], 1)
  }
  for (t in c("size_at_true_effect", "size_median")) {
    v <- s$value[s$term == t]
    expect_equal(v, floor(v))
  }
  expect_output(print(s), "Requested assurance")
})

test_that("ss_buc_sensitivity validates its inputs", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  expect_error(ss_buc_sensitivity(plan), "must specify 'true_ncp'", fixed = TRUE)
  expect_error(ss_buc_sensitivity(plan, true_ncp = -1), "positive number",
               fixed = TRUE)
  expect_error(ss_buc_sensitivity(plan, true_ncp = c(1, 2)),
               "single finite numeric", fixed = TRUE)
  expect_error(ss_buc_sensitivity(plan, true_ncp = 2, replications = 10.5),
               "whole number", fixed = TRUE)
  expect_error(ss_buc_sensitivity(data.frame(term = "x", value = 1), true_ncp = 2),
               "returned by one of the ss_buc_", fixed = TRUE)
})

test_that("the seed is honored and the caller's random stream is restored", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  a <- ss_buc_sensitivity(plan, true_ncp = 2, replications = 100, seed = 2017)
  b <- ss_buc_sensitivity(plan, true_ncp = 2, replications = 100, seed = 2017)
  expect_equal(a$value, b$value)

  set.seed(99)
  before <- get(".Random.seed", envir = .GlobalEnv)
  ss_buc_sensitivity(plan, true_ncp = 2, replications = 100, seed = 2017)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
})

test_that("size_at_true_effect is the size the true effect requires", {
  # By construction it is the plan the method issues when its corrected
  # parameter lands exactly on true_ncp, so a larger true effect must never
  # require a larger sample.
  plan <- ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3)
  ncp <- plan$value[plan$term == "ncp_adjusted"]
  sizes <- vapply(c(.6, 1, 1.6), function(m) {
    s <- ss_buc_sensitivity(plan, true_ncp = ncp * m, replications = 50,
                            seed = 2017)
    s$value[s$term == "size_at_true_effect"]
  }, numeric(1))
  expect_true(all(diff(sizes) < 0))
  # and the plan the object itself reports is the one its own parameter needs
  s <- ss_buc_sensitivity(plan, true_ncp = ncp, replications = 50, seed = 2017)
  expect_equal(s$value[s$term == "size_at_true_effect"],
               plan$value[plan$term == "necessary_N"])
})

test_that("the assurance guarantee holds when refusals are counted", {
  skip_on_cran()
  # The method's claim is P(corrected parameter <= true parameter) = assurance
  # over published prior studies, a refusal counting as attaining because
  # declining to plan is the limiting case of an arbitrarily conservative
  # plan. Checked here on designs whose required sample size is large enough
  # that the integer sample size does not blunt the comparison.
  reps <- 4000L
  se <- sqrt(.8 * .2 / reps)
  cases <- list(
    ss_buc_paired_t(t_observed = 3, N = 30, assurance = .80),
    ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3),
    ss_buc_chisq_diff(chisq_observed = 15, df_difference = 2, N = 200),
    ss_buc_correlation(r_observed = .35, N = 100)
  )
  for (plan in cases) {
    ncp <- plan$value[plan$term == "ncp_adjusted"]
    s <- ss_buc_sensitivity(plan, true_ncp = ncp * 1.4, replications = reps,
                            seed = 2017)
    expect_lt(abs(s$value[s$term == "attainment_with_refusals"] - .80), 4 * se)
  }
})
