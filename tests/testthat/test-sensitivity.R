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
  # The within-subjects designs under a sphericity correction: the registry
  # scales the df by the echoed epsilon, so the ceiling pin below checks the
  # corrected df, and the replay carries epsilon.
  ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3, epsilon = .7),
  ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2, epsilon = .7),
  ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                     levels_within = 3, effect = "within", epsilon = .7),
  ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                     levels_within = 3, effect = "interaction", epsilon = .7),
  ss_buc_mixed_anova_general(F_observed = 20, N = 60, num_groups = 2,
                             df_numerator = 2, effect = "within_only", epsilon = .7),
  ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 6, num_groups = 3,
                             effect = "between_within", df_num_within = 3,
                             epsilon = .6),
  ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3),
  ss_buc_R2(F_observed = 8, N = 60, p = 3),
  ss_buc_reg_joint(F_observed = 12, N = 60, p = 5, p_joint = 2),
  ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3, random_predictors = TRUE),
  ss_buc_R2(F_observed = 8, N = 60, p = 3, random_predictors = TRUE),
  ss_buc_reg_joint(F_observed = 12, N = 60, p = 5, p_joint = 2,
                   random_predictors = TRUE),
  ss_buc_manova(F_observed = 8, N = 40, p_variables = 3),
  ss_buc_chisq_diff(chisq_observed = 15, df_full = 42, df_restricted = 44, N = 200),
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
    inputs <- BUCSS:::.bucss_design_inputs(o)
    # The regression entries name a fixed and a random family; the plan's
    # random_predictors row picks one.
    spec <- BUCSS:::.bucss_resolve_spec(
      BUCSS:::.BUCSS_DESIGN_SPECS[[attr(o, "design")]], inputs)
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
      # At a zero effect the random-predictor mixture is the central F.
      f_random = 1 - pf(stat, df[1], df[2]),
      f_from_t_random = 1 - pf(stat^2, 1, df[2]),
      correlation = 1 - pf(stat^2, 1, df - 2),
      chisq = 1 - pchisq(stat, df))
    expect_equal(p_spec, alpha_prior * (1 - attr(o, "assurance_ceiling")),
                 tolerance = 1e-9,
                 info = paste("design:", attr(o, "design"),
                              "effect:", if (is.null(effect)) "-" else effect))
  }
})

test_that("the rebuilt planner call reproduces every plan exactly", {
  # every computed row, the planning effect size included
  computed <- BUCSS:::.BUCSS_COMPUTED_TERMS
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
      "attainment_with_refusals", "attainment_with_refusals_mcse",
      "size_at_true_effect", "size_q25",
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
  # Called without a plan, the same message, not R's 'argument "object" is
  # missing' (GF-12 of the 2026-10-12 triage).
  no_plan <- "'object' must be a plan returned by one of the ss_buc_* functions."
  expect_error(ss_buc_sensitivity(true_ncp = 2), no_plan, fixed = TRUE)
  expect_error(ss_buc_sensitivity(), no_plan, fixed = TRUE)
})

test_that("print.bucss_sensitivity validates 'digits', as print.bucss_plan does", {
  # GF-9 and CODE-1 of the 2026-10-12 triage: 'digits = NA' printed NA for
  # every value that is not a whole number, 0 and -1 printed one significant
  # figure, and a string, a vector, or NULL reached signif() or vapply() and
  # failed with base R's message.
  s <- ss_buc_sensitivity(ss_buc_paired_t(t_observed = 3, N = 30), true_ncp = 2,
                          replications = 30, seed = 2017)
  for (bad in list(NA, 0, -1, "a", c(2, 3), NULL, NA_real_))
    expect_error(capture.output(print(s, digits = bad)),
                 "'digits' must be a single positive number.", fixed = TRUE,
                 label = paste(deparse(bad), collapse = " "))
  # A valid value prints at that many significant figures.
  rate <- s$value[s$term == "publication_rate"]
  shown <- capture.output(print(s, digits = 5))
  expect_match(grep("publication_rate", shown, value = TRUE),
               format(signif(rate, 5), scientific = FALSE), fixed = TRUE)
})

test_that("the audit's footers never show an assurance or a power below 1 as 1.00", {
  # CL-37 of the 2026-10-12 triage: the two footers were formatted to two
  # decimals, so a plan made at .996 or .999 was described as "Requested
  # assurance: 1.00" and "Target power: 1.00", values no planner accepts.
  footer <- BUCSS:::.format_proportion
  expect_identical(footer(.999), ".999")
  expect_identical(footer(.996), ".996")
  expect_identical(footer(.9999), ".9999")
  expect_identical(footer(1 - 1e-12), ".999999999999")
  # The last double below 1 stops at the first decimal count whose text reads
  # below 1, sixteen. The rule reads the text and does not parse it: parsed,
  # "0.9999999999999999" is 1 on a build of R without extended precision and
  # below 1 on a build with it, so the footer had 17 decimals on the first
  # and 16 on the second.
  expect_identical(footer(1 - .Machine$double.eps / 2), ".9999999999999999")
  expect_identical(footer(1 - .Machine$double.eps), ".9999999999999998")
  # Every other value prints at two decimals, as before.
  expect_identical(vapply(c(.8, .5, .95, .825, .05, .99, .994), footer, character(1)),
                   c(".80", ".50", ".95", sub("^0", "", sprintf("%.2f", .825)), ".05", ".99",
                     ".99"))
  for (level in c(.996, .999)) {
    plan <- ss_buc_paired_t(t_observed = 6, N = 30, assurance = level, desired_power = level)
    s <- ss_buc_sensitivity(plan, true_ncp = 5, replications = 20, seed = 2017)
    shown <- capture.output(print(s))
    stated <- sub("^0", "", format(level))
    expect_true(paste0("Requested assurance: ", stated) %in% shown)
    expect_true(paste0("Target power: ", stated) %in% shown)
    footers <- grep("^(Requested assurance|Target power): ", shown, value = TRUE)
    expect_length(footers, 2L)
    expect_false(any(grepl("1.00", footers, fixed = TRUE)))
  }
})

test_that("the seed is honored and the caller's random stream is restored", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  a <- ss_buc_sensitivity(plan, true_ncp = 2, replications = 100, seed = 2017)
  b <- ss_buc_sensitivity(plan, true_ncp = 2, replications = 100, seed = 2017)
  expect_equal(a$value, b$value)

  # The caller's stream: seed 99 in the global state, restored when this test
  # ends so no seed is left behind.
  withr::local_seed(99)
  before <- get(".Random.seed", envir = .GlobalEnv)
  ss_buc_sensitivity(plan, true_ncp = 2, replications = 100, seed = 2017)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
})

test_that("a seeded call restores the caller's generator kind", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  withr::local_seed(99, .rng_kind = "L'Ecuyer-CMRG")
  kind <- RNGkind()
  expected <- withr::with_preserve_seed(runif(3))
  ss_buc_sensitivity(plan, true_ncp = 2, replications = 100, seed = 2017)
  expect_identical(RNGkind(), kind)
  expect_identical(runif(3), expected)
})

test_that("a seeded call leaves no seed behind when the caller had none", {
  # Without this, the caller's next unseeded draw would be a fixed
  # continuation of seed 2017 in every session.
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  withr::local_preserve_seed()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    rm(".Random.seed", envir = globalenv())
  ss_buc_sensitivity(plan, true_ncp = 2, replications = 100, seed = 2017)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("loading BUCSS enforces the withr version its seed handling needs", {
  # CP-1 of the 2026-10-12 triage. In withr 3.0.2 and earlier, local_seed()
  # called inside a function restores nothing, so the three tests above fail
  # there, and DESCRIPTION asks for withr (>= 3.0.3). R applies a version
  # requirement in Imports only to the namespaces NAMESPACE imports from: a
  # 'withr::' call alone imports nothing, and with no import an older withr
  # already loaded, or earlier on the library path, was used without
  # complaint. With the import, loading BUCSS stops with "namespace 'withr'
  # 3.0.2 is being loaded, but >= 3.0.3 is required".
  imports <- getNamespaceImports(asNamespace("BUCSS"))
  expect_true("local_seed" %in% unlist(imports[names(imports) == "withr"], use.names = FALSE))
  expect_match(utils::packageDescription("BUCSS")$Imports, "withr (>= 3.0.3)", fixed = TRUE)
})

test_that("a seed too large for set.seed() gets the package's own error", {
  # 2^31 passed the whole-number check and then failed inside set.seed() with
  # base R's error and a coercion warning (QC 2026-09-30, finding 11).
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  expect_error(ss_buc_sensitivity(plan, true_ncp = 2, replications = 10, seed = 2^31),
               "'seed' must be at most 2147483647.", fixed = TRUE)
  expect_s3_class(ss_buc_sensitivity(plan, true_ncp = 2, replications = 10,
                                     seed = .Machine$integer.max), "bucss_sensitivity")
})

test_that("the default seed = NULL draws from the caller's own stream", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  run <- function(s) withr::with_seed(s, ss_buc_sensitivity(
    plan, true_ncp = 2, replications = 100))
  expect_identical(run(5), run(5))
  expect_false(identical(run(5)$value, run(6)$value))
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

test_that("a simulated prior study in the rounding band counts as a refusal", {
  # TQ-2 of the 2026-10-12 triage. A one-way prior of 13 in 4 groups is
  # planned at 3 and 4 per group; a simulated study whose F is significant at
  # its own df (3, 9) and not at the rounded-down design's (3, 8) meets the
  # planner's nonsignificance refusal, which the audit counts with the
  # refusals, as it counts the zero noncentrality parameter refusal. With the
  # handler for that class removed the audit stops ("While replaying the
  # planner: Your observed F statistic is significant at the prior study's own
  # degrees of freedom ...") and the whole suite still passed.
  plan <- ss_buc_one_way_anova(F_observed = 6, N = 13, levels_A = 4, assurance = .61)
  s <- ss_buc_sensitivity(plan, true_ncp = 3, replications = 50, seed = 2017)
  expect_s3_class(s, "bucss_sensitivity")
  expect_gt(s$value[s$term == "refusal_rate"], 0)
  # The band is real for this design: such statistics exist, and the planner
  # refuses one with that class.
  f_band <- (qf(.95, 3, 9) + qf(.95, 3, 8)) / 2
  stopifnot(pf(f_band, 3, 9, lower.tail = FALSE) < .05,
            pf(f_band, 3, 8, lower.tail = FALSE) > .05)
  expect_error(ss_buc_one_way_anova(F_observed = f_band, N = 13, levels_A = 4, assurance = .61),
               "significant at the prior study's own degrees of freedom", fixed = TRUE,
               class = "bucss_nonsignificant")
})

test_that("the shipped verification script is one invisible(local({ ... })) expression", {
  # TQ-4 of the 2026-10-12 triage. inst/tools/verify_by_simulation.R runs
  # inside local() so that sourcing it leaves nothing in the user's workspace
  # (about 90 objects otherwise). Nothing checked that: with local() replaced
  # by identity() the suite passed, the script printed the same "24 checks, 0
  # failed.", and the release gate's static detector found nothing.
  f <- system.file("tools", "verify_by_simulation.R", package = "BUCSS")
  skip_if(!nzchar(f), "the shipped verification script was not found")
  e <- parse(f, keep.source = FALSE)
  expect_length(e, 1L)
  expect_identical(e[[1]][[1]], as.name("invisible"))
  expect_identical(e[[1]][[2]][[1]], as.name("local"))
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
    ss_buc_chisq_diff(chisq_observed = 15, df_full = 42, df_restricted = 44, N = 200),
    ss_buc_correlation(r_observed = .35, N = 100)
  )
  for (plan in cases) {
    ncp <- plan$value[plan$term == "ncp_adjusted"]
    s <- ss_buc_sensitivity(plan, true_ncp = ncp * 1.4, replications = reps,
                            seed = 2017)
    expect_lt(abs(s$value[s$term == "attainment_with_refusals"] - .80), 4 * se)
  }
})
