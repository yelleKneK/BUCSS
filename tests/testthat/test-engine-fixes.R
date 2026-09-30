# The numerical fixes made after the 2026-09-28 deep QC audit. Every expected
# value below is what the package returns now, and was also confirmed without
# BUCSS:
# either it is an auditor's independent value (finding id named), or it was
# computed from the method's definition in a separate script that never loads
# BUCSS (bisection for the root, a linear walk for the planned n, the two-sided
# rounding rule re-derived from its description, and quadrature where the
# comment says so), marked "independent check 2026-09-28".

size <- function(x) x$value[x$term %in% BUCSS:::.BUCSS_SIZE_TERMS]
val <- function(x, term) x$value[x$term == term]

# with_time_limit() is in helper-time-limit.R.

# Power of a two-sided noncentral t test, from first principles.
power_t2 <- function(df, ncp, alpha = .05) {
  crit <- qt(1 - alpha / 2, df)
  pt(crit, df, ncp, lower.tail = FALSE) + pt(-crit, df, ncp)
}

test_that("the Welch planner never returns fewer than 10 per group", {
  expect_identical(BUCSS:::.welch_min_n, 10)
  # Strong priors whose target power is reached below 10 per group. The
  # unfloored answers are 2, 5, and 8 per group (independent check
  # 2026-09-28); each plan is raised to 10, and its actual_power is the power
  # at 10, above the target.
  strong <- ss_buc_welch_t(t_observed = 12, n_1 = 8, n_2 = 8, assurance = .5)
  expect_identical(size(strong), 10)
  expect_equal(val(strong, "ncp_adjusted"), 11.7255892480, tolerance = 1e-9)
  # With sd_ratio = 1 the Welch-Satterthwaite df at 10 per group is 18, and
  # the ncp rescales by the ratio of the two standard errors.
  expect_equal(val(strong, "actual_power"),
               power_t2(18, sqrt((1 / 8 + 1 / 8) / (2 / 10)) * val(strong, "ncp_adjusted")))

  p8 <- ss_buc_welch_t(t_observed = 8, n_1 = 20, n_2 = 20)
  expect_identical(size(p8), 10)
  expect_equal(val(p8, "actual_power"), 0.9950870091, tolerance = 1e-8)
  p6 <- ss_buc_welch_t(t_observed = 6, n_1 = 20, n_2 = 20)
  expect_identical(size(p6), 10)
  expect_equal(val(p6, "actual_power"), 0.9097893690, tolerance = 1e-8)
  expect_gt(val(p6, "actual_power"), .8)

  # A plan above the floor is unchanged by it.
  p5 <- ss_buc_welch_t(t_observed = 5, n_1 = 20, n_2 = 20)
  expect_identical(size(p5), 11)
  expect_equal(val(p5, "ncp_adjusted"), 3.9738703315, tolerance = 1e-8)
  expect_equal(val(p5, "actual_power"), 0.8004339265, tolerance = 1e-8)
})

test_that("the Welch help page's five example plans are unchanged by the floor", {
  # Sizes as reported in the audit (CC-03); sizes and ncp_adjusted confirmed by
  # independent check 2026-09-28.
  ratio_15 <- list(
    ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5,
                   alpha_prior = .05, alpha_planned = .05, assurance = .80,
                   desired_power = .80),
    ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_1 = 8.4, sd_2 = 12.6),
    ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, var_1 = 70.56, var_2 = 158.76))
  for (r in ratio_15) {
    expect_identical(size(r), 222)
    expect_equal(val(r, "ncp_adjusted"), 1.3246896477, tolerance = 1e-8)
  }
  r1 <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1)
  expect_identical(size(r1), 213)
  expect_equal(val(r1, "ncp_adjusted"), 1.3107728398, tolerance = 1e-8)
  r2 <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 2)
  expect_identical(size(r2), 236)
  expect_equal(val(r2, "ncp_adjusted"), 1.3100248639, tolerance = 1e-8)
})

test_that("two equal group sizes plan exactly like the common size", {
  # The harmonic mean of 49 and 49 computed as 2 / (1/49 + 1/49) is
  # 49.000000000000007, which used to open a spurious second rounding branch
  # at 50 (plan 222 instead of 219). The identity is checked for every k in
  # 2:300, 27 of which have an inexact floating-point harmonic mean; a prior
  # that the planner refuses must be refused the same way both times.
  computed <- function(n) {
    tryCatch({
      r <- ss_buc_independent_t(t_observed = 3, n = n)
      list(terms = r$term[!r$term %in% c("n", "n_1", "n_2")],
           values = r$value[!r$term %in% c("n", "n_1", "n_2")],
           ceiling = attr(r, "assurance_ceiling"), df_error = attr(r, "df_error"))
    }, error = function(e) class(e)[1])
  }
  pair <- ss_buc_independent_t(t_observed = 3, n = c(49, 49))
  single <- ss_buc_independent_t(t_observed = 3, n = 49)
  expect_identical(computed(c(49, 49)), computed(49))
  expect_identical(size(pair), 219)                      # independent check 2026-09-28
  expect_equal(val(pair, "ncp_adjusted"), 1.3300239601, tolerance = 1e-8)
  expect_identical(pair$value[pair$term %in% c("n_1", "n_2")], c(49, 49))
  expect_identical(single$value[single$term == "n"], 49)

  mismatched <- Filter(function(k) !identical(computed(c(k, k)), computed(k)), 2:300)
  expect_identical(mismatched, integer(0))
})

test_that("an exact integer harmonic mean gives a single rounding branch", {
  # 2 * 10 * 15 / (10 + 15) = 12 exactly; the floating-point quotient was
  # 11.999999999999998, which planned from 11 and 12 and returned 177 with ncp
  # .7023. From the harmonic mean itself (independent check 2026-09-28):
  # 154 per group, ncp .7858559562, actual power .8013395467.
  expect_identical((2 * 10 * 15) %/% (10 + 15), 12)
  expect_identical((2 * 10 * 15) %% (10 + 15), 0)
  r <- ss_buc_independent_t(t_observed = 3, n = c(10, 15))
  expect_identical(size(r), 154)
  expect_equal(val(r, "ncp_adjusted"), 0.7858559562, tolerance = 1e-8)
  expect_equal(val(r, "actual_power"), 0.8013395467, tolerance = 1e-8)
  twelve <- ss_buc_independent_t(t_observed = 3, n = 12)
  keep <- !r$term %in% c("n", "n_1", "n_2")
  expect_identical(r$value[keep], twelve$value[twelve$term != "n"])
})

test_that("the planned sample size search is exact up to 2^53 and stops only above 2^52", {
  search <- BUCSS:::.smallest_n_for_power
  # Doubles hold every whole number exactly up to 2^53; beyond it adjacent
  # doubles are 2 apart, and the old bisection cycled forever at hi - lo == 2.
  # The bracket now stops doubling once it has passed 2^52 without reaching
  # the target, so it never extends past 2^53, and the search stops with a
  # message instead. That happens only when the answer exceeds 2^52; where
  # exactly depends on 'start', at the first start * 2^k above 2^52.
  expect_identical((2^52 + 1) - 2^52, 1)                # still exact above 2^52
  expect_identical((2^53 + 1) - 2^53, 0)                # not above 2^53
  # From the default start of 2 the last bracket is 2^53: an answer above it
  # stops, and an answer between 2^52 and 2^53 is returned exactly.
  expect_error(
    with_time_limit(search(function(n) as.numeric(n >= 2^53 + 12345), target = .5)),
    "would exceed about 4.5e15 (2^52)", fixed = TRUE)
  expect_identical(with_time_limit(search(function(n) as.numeric(n >= 2^52 + 3), .5)),
                   2^52 + 3)
  expect_identical(with_time_limit(search(function(n) as.numeric(n >= 2^53 - 1), .5)),
                   2^53 - 1)
  # From a start of 3 the last bracket is 3 * 2^51, between 2^52 and 2^53: an
  # answer just above it stops, and one just below it is returned exactly.
  expect_error(
    with_time_limit(search(function(n) as.numeric(n >= 3 * 2^51 + 1), .5, start = 3)),
    "2^52", fixed = TRUE)
  expect_identical(
    with_time_limit(search(function(n) as.numeric(n >= 3 * 2^51 - 1), .5, start = 3)),
    3 * 2^51 - 1)
  # Below 2^52 every answer is found exactly.
  expect_identical(with_time_limit(search(function(n) as.numeric(n >= 2^52 - 345), .5)),
                   2^52 - 345)
  expect_identical(search(function(n) as.numeric(n >= 12345), .5), 12345)
  expect_identical(search(function(n) 1, .5, start = 7), 7)
})

test_that("the correlation series sums over the mass it needs", {
  # At rho = .99 and N = 50000 the negative binomial mass sits near k = 1.23e6,
  # above the old fixed cap of 1e6 terms, which returned 2.3e-217 here (R-01).
  # Reference: nested finite-limit quadrature of the exact representation
  # F = (Z + sqrt(c W))^2 / (V / (N - 2)) over W ~ chi square(N - 1) and
  # V ~ chi square(N - 2), with no series and no pf() (independent check
  # 2026-09-28: 0.499116838735014; the audit's refuters found .4991 by a
  # windowed series and by Monte Carlo).
  N <- 50000
  q <- (N - 2) * .99^2 / (1 - .99^2)
  expect_equal(BUCSS:::.p_correlation_f(q, N, .99^2), 0.499116838735014,
               tolerance = 1e-9)
  # Below a million terms the window starts at 0, as the full sum does; the
  # same nested quadrature gives 0.0883937021306566 at rho = .3, N = 100,
  # q = 3 (independent check 2026-09-28).
  expect_equal(BUCSS:::.p_correlation_f(3, 100, .09), 0.0883937021306566,
               tolerance = 1e-12)
})

test_that("strong correlations from very large priors plan correctly", {
  skip_on_cran()  # about 1.5 s
  # The audit's independent values (R-01): n = 9 at r = .8, N = 2e6 (the
  # truncated series gave 13) and n = 29 at r = .5, N = 1e7 (gave 44), with
  # corrected parameters 3.55027e6 and 3.32979e6.
  r8 <- ss_buc_correlation(r_observed = .8, N = 2e6)
  expect_identical(size(r8), 9)
  expect_equal(val(r8, "ncp_adjusted"), 3.55027e6, tolerance = 1e-5)
  r5 <- ss_buc_correlation(r_observed = .5, N = 1e7)
  expect_identical(size(r5), 29)
  expect_equal(val(r5, "ncp_adjusted"), 3.32979e6, tolerance = 1e-5)
})

test_that("R's pnt precision warning is muffled, and only that warning", {
  expect_no_warning(r <- ss_buc_independent_t(t_observed = 15, n = 50))
  expect_no_warning(ss_buc_paired_t(t_observed = 12, N = 250, desired_power = .9))
  expect_no_warning(ss_buc_one_sample_t(t_observed = 15, N = 50))
  expect_no_warning(ss_buc_welch_t(t_observed = 15, n_1 = 50, n_2 = 50))
  # The muffled evaluations are accurate. Reference: noncentral t tail areas
  # by quadrature over the chi square mixture, no pt() (independent check
  # 2026-09-28): ncp 13.7270669032588; planned 4 per group, power
  # .895703411828067.
  expect_equal(val(r, "ncp_adjusted"), 13.7270669032588, tolerance = 1e-9)
  expect_identical(size(r), 4)
  expect_equal(val(r, "actual_power"), 0.895703411828067, tolerance = 1e-9)

  expect_warning(BUCSS:::.quiet_pnt(warning("something else")), "something else",
                 fixed = TRUE)
})

test_that(".quiet_pnt() matches the warning by its call and routine, in any language", {
  # The warning as R raises it: its call is pt(), and its message, translated
  # or not, names the routine 'pnt{final}'. Here in German, as R's own
  # translation words it.
  from_pt <- function(msg) warningCondition(msg, call = quote(pt(15, 98, ncp = 1)))
  expect_no_warning(BUCSS:::.quiet_pnt(warning(from_pt(
    "full precision may not have been achieved in 'pnt{final}'"))))
  expect_no_warning(BUCSS:::.quiet_pnt(warning(from_pt(
    "eventuell wurde in 'pnt{final}' nicht die volle Genauigkeit erreicht"))))
  # Any other warning reaches the user: another warning from pt(), a pnt
  # message from any other call, and one with no call at all.
  expect_warning(BUCSS:::.quiet_pnt(warning(from_pt("NaNs produced"))),
                 "NaNs produced", fixed = TRUE)
  expect_warning(BUCSS:::.quiet_pnt(warning(warningCondition(
    "eventuell wurde in 'pnt{final}' nicht die volle Genauigkeit erreicht",
    call = quote(qt(.5, 3))))), "pnt{final}", fixed = TRUE)
  expect_warning(BUCSS:::.quiet_pnt(warning(warningCondition(
    "full precision may not have been achieved in 'pnt{final}'", call = NULL))),
    "pnt{final}", fixed = TRUE)
})

test_that("the pnt warning is muffled in a German session", {
  # testthat runs each test in English, so a test that only calls the planners
  # cannot see a translated message; this one switches the language itself.
  # Translations need NLS, and an LC_ALL of C overrides the language.
  skip_if_not(capabilities("NLS"), "R was built without NLS")
  skip_if(Sys.getenv("LC_ALL") %in% c("C", "C.UTF-8", "POSIX"),
          "LC_ALL is set to a C locale")
  msgid <- "full precision may not have been achieved in '%s'\n"
  translated <- withr::with_language("de", gettext(msgid, domain = "R", trim = FALSE))
  skip_if(identical(translated, msgid), "R's German translations are not installed")
  withr::local_language("de")
  # The raw call does warn, in German, so the planners below meet the warning.
  w <- tryCatch(pt(15, 98, ncp = 1), warning = identity)
  skip_if_not(inherits(w, "warning"), "pt() raises no precision warning in this R")
  expect_match(conditionMessage(w), "pnt{final}", fixed = TRUE)
  expect_no_match(conditionMessage(w), "full precision", fixed = TRUE)
  expect_no_warning(ss_buc_independent_t(t_observed = 15, n = 50))
  expect_no_warning(ss_buc_paired_t(t_observed = 12, N = 250, desired_power = .9))
  expect_no_warning(ss_buc_one_sample_t(t_observed = 15, N = 50))
  expect_no_warning(ss_buc_welch_t(t_observed = 15, n_1 = 50, n_2 = 50))
})

test_that("the two refusals carry their own condition classes", {
  expect_error(ss_buc_independent_t(t_observed = 3, n = 20, assurance = .95),
               class = "bucss_zero_ncp")
  expect_error(ss_buc_independent_t(t_observed = 1, n = 20),
               class = "bucss_nonsignificant")
  zero <- tryCatch(ss_buc_paired_t(t_observed = 3, N = 40, assurance = .95),
                   error = identity)
  expect_s3_class(zero, c("bucss_zero_ncp", "error", "condition"), exact = TRUE)
  expect_match(conditionMessage(zero), "noncentrality parameter is zero", fixed = TRUE)
  expect_null(conditionCall(zero))
  nonsig <- tryCatch(ss_buc_paired_t(t_observed = 1, N = 40), error = identity)
  expect_s3_class(nonsig, c("bucss_nonsignificant", "error", "condition"), exact = TRUE)
  expect_match(conditionMessage(nonsig),
               "is nonsignificant based on your specified 'alpha_prior'", fixed = TRUE)
  expect_null(conditionCall(nonsig))
  # An input error is neither.
  bad <- tryCatch(ss_buc_paired_t(t_observed = 3, N = NA), error = identity)
  expect_false(inherits(bad, "bucss_zero_ncp"))
  expect_false(inherits(bad, "bucss_nonsignificant"))
})

test_that("the zero-NCP message names the smaller branch ceiling, a feasible one", {
  # One-way, F = 6, N = 13, four groups: the rounded-down design (3 per group,
  # df 3 and 8) supports .6171712195 and the rounded-up one (4 per group, df 3
  # and 12) .805383381 (independent check 2026-09-28; the audit, B-01, found
  # .6172 and .8054). The message used to name .80, which then failed too.
  ceiling_rd <- 1 - pf(6, 3, 8, lower.tail = FALSE) / .05
  ceiling_ru <- 1 - pf(6, 3, 12, lower.tail = FALSE) / .05
  expect_equal(ceiling_rd, 0.6171712195, tolerance = 1e-9)
  expect_equal(ceiling_ru, 0.805383381, tolerance = 1e-8)

  msg <- tryCatch(ss_buc_one_way_anova(F_observed = 6, N = 13, levels_A = 4,
                                       assurance = .95),
                  bucss_zero_ncp = conditionMessage)
  expect_match(msg, "usable plan is about .61;", fixed = TRUE)
  expect_no_match(msg, "about .80", fixed = TRUE)

  plan <- ss_buc_one_way_anova(F_observed = 6, N = 13, levels_A = 4, assurance = .61)
  expect_identical(size(plan), 213)                     # independent check 2026-09-28
  expect_equal(val(plan, "ncp_adjusted"), 0.1547029449, tolerance = 1e-8)
  expect_equal(attr(plan, "assurance_ceiling"), ceiling_rd, tolerance = 1e-10)
  expect_error(ss_buc_one_way_anova(F_observed = 6, N = 13, levels_A = 4, assurance = .62),
               class = "bucss_zero_ncp")

  # The general planner takes the same path.
  expect_error(ss_buc_factorial_anova_general(F_observed = 6, N = 13, cells = 4,
                                              df_numerator = 3, df_denominator = 9,
                                              assurance = .95),
               "usable plan is about .61;", fixed = TRUE)
})

test_that("every two-branch planner names the smaller branch ceiling, and it plans", {
  # For each planner that rounds an unbalanced prior study two ways, a case
  # whose two branch ceilings floor to different two-decimal values, so a
  # planner that named one branch's ceiling alone (as the round-up branch's
  # used to be named) would name the larger one. Each ceiling is
  # 1 - p/alpha_prior at that branch's degrees of freedom, computed here from
  # pt() and pf() alone (independent check 2026-09-29); in every case the
  # round-down branch, with fewer error degrees of freedom, has the smaller.
  # The message must name the smaller, floored; a re-run at that value must
  # return a plan, and a re-run at the larger branch's floored value refuses.
  fl <- function(x) floor(x * 100 - 1e-9) / 100
  ceiling_f <- function(F, df1, df2) 1 - pf(F, df1, df2, lower.tail = FALSE) / .05
  ceiling_t <- function(t, df) 1 - 2 * pt(t, df, lower.tail = FALSE) / .05
  cases <- list(
    # n = c(5, 200): harmonic mean 9.76, so 9 and 10 per group, df 16 and 18
    # (.6132 and .6380).
    list(label = "independent t", rd = ceiling_t(2.6, 16), ru = ceiling_t(2.6, 18),
         named = ".61",
         run = function(a) ss_buc_independent_t(t_observed = 2.6, n = c(5, 200),
                                                assurance = a)),
    # 2 x 2, N = 13: 3 and 4 per cell, df (1, 8) and (1, 12) (.6586 and .7787).
    list(label = "two-way", rd = ceiling_f(9, 1, 8), ru = ceiling_f(9, 1, 12),
         named = ".65",
         run = function(a) ss_buc_factorial_anova(F_observed = 9, N = 13, levels_A = 2,
                                                  levels_B = 2, effect = "factor_A",
                                                  assurance = a)),
    # Four cells, N = 13: df (3, 8) and (3, 12) (.6172 and .8054).
    list(label = "between-subjects general", rd = ceiling_f(6, 3, 8),
         ru = ceiling_f(6, 3, 12), named = ".61",
         run = function(a) ss_buc_factorial_anova_general(F_observed = 6, N = 13, cells = 4,
                                                          df_numerator = 3,
                                                          df_denominator = 9,
                                                          assurance = a)),
    # Four cells and one covariate, N = 13: df (3, 7) and (3, 11) (.6736 and
    # .8663).
    list(label = "ANCOVA", rd = ceiling_f(7, 3, 7), ru = ceiling_f(7, 3, 11),
         named = ".67",
         run = function(a) ss_buc_ancova(F_observed = 7, N = 13, cells = 4,
                                         n_covariates = 1, assurance = a)),
    # Four groups, N = 13, between effect: df (3, 8) and (3, 12).
    list(label = "split plot", rd = ceiling_f(6, 3, 8), ru = ceiling_f(6, 3, 12),
         named = ".61",
         run = function(a) ss_buc_mixed_anova(F_observed = 6, N = 13, levels_between = 4,
                                              levels_within = 3, effect = "between",
                                              assurance = a)),
    list(label = "split plot general", rd = ceiling_f(6, 3, 8),
         ru = ceiling_f(6, 3, 12), named = ".61",
         run = function(a) ss_buc_mixed_anova_general(F_observed = 6, N = 13,
                                                      df_numerator = 3, num_groups = 4,
                                                      effect = "between_only",
                                                      assurance = a)),
    # Three outcomes, N = 13: 6 and 7 per group, df (3, 8) and (3, 10) (.7485
    # and .8384).
    list(label = "MANOVA", rd = ceiling_f(7, 3, 8), ru = ceiling_f(7, 3, 10),
         named = ".74",
         run = function(a) ss_buc_manova(F_observed = 7, N = 13, p_variables = 3,
                                         assurance = a))
  )
  for (case in cases) {
    info <- case$label
    expect_lt(case$rd, case$ru)
    expect_false(fl(case$rd) == fl(case$ru), info = info)
    expect_identical(sub("^0", "", sprintf("%.2f", fl(case$rd))), case$named, info = info)
    msg <- tryCatch(case$run(.95), bucss_zero_ncp = conditionMessage)
    expect_match(msg, paste0("usable plan is about ", case$named, ";"), fixed = TRUE,
                 info = info)
    plan <- case$run(fl(case$rd))
    expect_s3_class(plan, "bucss_plan")
    expect_equal(attr(plan, "assurance_ceiling"), case$rd, tolerance = 1e-10, info = info)
    expect_error(case$run(fl(case$ru)), class = "bucss_zero_ncp", info = info)
  }
})

test_that("a prior significant at its own df is not called nonsignificant", {
  # n = c(5, 200): the prior study's own df are 203, where t = 2.05 has
  # p = .04165; the harmonic mean 9.76 rounds to 9 and 10 per group (df 16 and
  # 18), where it is not significant. The message says so, reports the own-df
  # p, and the ceiling that p supports, 1 - .04165/.05 = .167, floored to .16
  # (independent check 2026-09-28; T-04).
  p_own <- 2 * pt(2.05, 203, lower.tail = FALSE)
  expect_equal(p_own, 0.0416498694, tolerance = 1e-8)
  expect_gt(2 * pt(2.05, 18, lower.tail = FALSE), .05)
  msg <- tryCatch(ss_buc_independent_t(t_observed = 2.05, n = c(5, 200)),
                  bucss_nonsignificant = conditionMessage)
  expect_match(msg, "is significant at the prior study's own degrees of freedom (p = .0416)",
               fixed = TRUE)
  expect_match(msg, "supports an assurance of at most about .16, below the recommended floor",
               fixed = TRUE)
  expect_no_match(msg, "is nonsignificant", fixed = TRUE)

  # A prior that is nonsignificant at its own df still gets the old message.
  expect_gt(2 * pt(1.9, 203, lower.tail = FALSE), .05)
  msg <- tryCatch(ss_buc_independent_t(t_observed = 1.9, n = c(5, 200)),
                  bucss_nonsignificant = conditionMessage)
  expect_match(msg, "is nonsignificant based on your specified 'alpha_prior'", fixed = TRUE)
  expect_no_match(msg, "own degrees of freedom", fixed = TRUE)
})

test_that("the band message reports the balanced design's p value and an alpha_prior that plans", {
  # The message reports the p value at the rounded-down balanced design (the
  # one with the fewest degrees of freedom, where the statistic fell short) and
  # an 'alpha_prior' at which planning at the requested assurance works: the
  # threshold p/(1 - assurance), rounded up to three significant digits,
  # strictly above it (at the threshold itself the ceiling equals the
  # assurance and the planner refuses). It mentions the recommended floor of
  # .5 only when the own-df ceiling is below it. The p values are from pt()
  # alone (independent check 2026-09-29); the printed values are pinned as
  # text and re-run (QC 2026-09-30, finding 2).
  fmt3 <- function(p) sub("^0", "", format(signif(p, 3), scientific = FALSE))
  printed_alpha <- function(msg)
    as.numeric(sub(".*an 'alpha_prior' of ([.0-9]+) or more.*", "\\1", msg))
  opening <- paste0("Your observed t statistic is significant at the prior study's own ",
                    "degrees of freedom (p = ")
  advice <- paste0("Near that value the plan is very large, because the result then ",
                   "barely supports the requested assurance; a larger 'alpha_prior' ",
                   "gives a practical plan, and so does a lower 'assurance', which ",
                   "lowers the 'alpha_prior' needed.")

  # Below .5: t = 2.05, n = c(5, 200). The harmonic mean 9.76 rounds down to 9
  # per group (df 16), where p = .0571; the own df are 203 (p = .0416, ceiling
  # 1 - .0416/.05, floored to .16). The threshold is .0571/.2 = .28561.
  p_bal <- 2 * pt(2.05, 16, lower.tail = FALSE)
  expect_equal(p_bal / .2, 0.2856145, tolerance = 1e-6)
  low <- tryCatch(ss_buc_independent_t(t_observed = 2.05, n = c(5, 200)),
                  bucss_nonsignificant = conditionMessage)
  expect_match(low, paste0(opening, ".0416)"), fixed = TRUE)
  expect_match(low, paste0("and there its p value is ", fmt3(p_bal), ", above 'alpha_prior'."),
               fixed = TRUE)
  expect_match(low, paste0("Planning from it at an 'assurance' of .80 would need an ",
                           "'alpha_prior' of .286 or more, which assumes a more lenient ",
                           "publication threshold. "), fixed = TRUE)
  expect_match(low, advice, fixed = TRUE)
  expect_match(low, paste0("Even at its own degrees of freedom the result supports an ",
                           "assurance of at most about .16, below the recommended floor ",
                           "of .5."), fixed = TRUE)
  # The advice is true: at the printed value the planner plans (59,048 per
  # group, very large), and at the threshold itself it still refuses.
  expect_identical(printed_alpha(low), .286)
  expect_identical(size(ss_buc_independent_t(t_observed = 2.05, n = c(5, 200),
                                             alpha_prior = .286)), 59048)
  expect_error(ss_buc_independent_t(t_observed = 2.05, n = c(5, 200),
                                    alpha_prior = p_bal / .2 * (1 - 1e-9)),
               class = "bucss_zero_ncp")

  # .5 or more: t = 2.77, n = c(2, 1000). The own df are 1000 (p = .00571,
  # ceiling .88); the harmonic mean 3.99 rounds down to 3 per group (df 4),
  # where p = .0503.
  p_bal <- 2 * pt(2.77, 4, lower.tail = FALSE)
  high <- tryCatch(ss_buc_independent_t(t_observed = 2.77, n = c(2, 1000)),
                   bucss_nonsignificant = conditionMessage)
  expect_match(high, paste0(opening, ".00571)"), fixed = TRUE)
  expect_match(high, paste0("its p value is ", fmt3(p_bal), ","), fixed = TRUE)
  expect_match(high, "an 'alpha_prior' of .252 or more", fixed = TRUE)
  expect_gt(.252, p_bal / .2)
  expect_s3_class(ss_buc_independent_t(t_observed = 2.77, n = c(2, 1000), alpha_prior = .252),
                  "bucss_plan")
  expect_no_match(high, "below the recommended floor", fixed = TRUE)
  expect_no_match(high, "at least about", fixed = TRUE)

  # The audit's case: N = 9 rounds to 4 per group (df 6), where p = .0529, so
  # the threshold is .26439. Rounded to the nearest three digits it was .264,
  # at which the planner refuses; the message now names .265, which plans.
  t9 <- (qt(.975, 7) + qt(.975, 6)) / 2
  n9 <- tryCatch(ss_buc_independent_t(t9, N = 9, assurance = .8),
                 bucss_nonsignificant = conditionMessage)
  expect_match(n9, "an 'alpha_prior' of .265 or more", fixed = TRUE)
  expect_identical(size(ss_buc_independent_t(t9, N = 9, assurance = .8, alpha_prior = .265)),
                   11246)
  expect_error(ss_buc_independent_t(t9, N = 9, assurance = .8, alpha_prior = .264),
               class = "bucss_zero_ncp")

  # An assurance so high that no practical alpha_prior will do: the threshold
  # is 1.26 here, and at .999 or more even 'alpha_prior = 1' (read as .999)
  # refuses. The assurance prints as entered, never as 1.00.
  none <- tryCatch(ss_buc_independent_t(t_observed = 2.77, n = c(2, 1000), assurance = .96),
                   bucss_nonsignificant = conditionMessage)
  expect_match(none, paste0("Planning from it at an 'assurance' of .96 would need an ",
                            "'alpha_prior' above .999, and 'alpha_prior = 1' (no publication ",
                            "bias) is treated as .999, so no practical 'alpha_prior' allows ",
                            "it; a lower 'assurance' lowers the 'alpha_prior' needed."),
               fixed = TRUE)
  p9 <- 2 * pt(t9, 6, lower.tail = FALSE)
  edge <- 1 - p9 / .9995                         # a threshold of .9995
  msg <- tryCatch(ss_buc_independent_t(t9, N = 9, assurance = edge),
                  bucss_nonsignificant = conditionMessage)
  expect_match(msg, "would need an 'alpha_prior' above .999", fixed = TRUE)
  expect_no_match(msg, "or more", fixed = TRUE)
  expect_error(ss_buc_independent_t(t9, N = 9, assurance = edge, alpha_prior = 1),
               class = "bucss_zero_ncp")
  top <- tryCatch(ss_buc_independent_t(t_observed = 2.77, n = c(2, 1000), assurance = .999),
                  bucss_nonsignificant = conditionMessage)
  expect_match(top, "at an 'assurance' of .999 would need", fixed = TRUE)
  expect_no_match(top, "1.00", fixed = TRUE)
})

test_that("every rounding planner's band message names an alpha_prior that plans", {
  # One unbalanced prior study per planner that reduces it to a balanced
  # design, with a statistic halfway between the critical values at the prior
  # study's own error df and at the rounded-down design's, at five assurances.
  # Rounding the threshold to the nearest three digits named a refused value
  # in 30 of these 40 cases (QC 2026-09-30, finding 2); every printed value
  # must now be above the threshold and plan.
  mid <- function(a, b) (a + b) / 2
  f_mid <- mid(qf(.95, 2, 10), qf(.95, 2, 9))
  f_mid2 <- mid(qf(.95, 2, 19), qf(.95, 2, 18))
  designs <- list(
    function(...) ss_buc_independent_t(t_observed = mid(qt(.975, 7), qt(.975, 6)), N = 9, ...),
    function(...) ss_buc_independent_t(t_observed = 2.05, n = c(5, 200), ...),
    function(...) ss_buc_one_way_anova(F_observed = f_mid, N = 13, levels_A = 3, ...),
    function(...) ss_buc_factorial_anova(F_observed = f_mid2, N = 25, levels_A = 2, levels_B = 3,
                                         effect = "factor_B", ...),
    function(...) ss_buc_factorial_anova_general(F_observed = f_mid2, N = 25, cells = 6,
                                                 df_numerator = 2, df_denominator = 19, ...),
    function(...) ss_buc_mixed_anova(F_observed = f_mid, N = 13, levels_between = 3,
                                     levels_within = 3, effect = "between", ...),
    function(...) ss_buc_mixed_anova_general(F_observed = f_mid, N = 13, df_numerator = 2,
                                             num_groups = 3, effect = "between_only", ...),
    function(...) ss_buc_manova(F_observed = f_mid, N = 13, p_variables = 2, ...))
  named <- 0L
  for (d in designs) for (a in c(.5, .6, .7, .8, .9)) {
    msg <- tryCatch(d(assurance = a), bucss_nonsignificant = conditionMessage)
    expect_type(msg, "character")
    expect_match(msg, "or more, which assumes", fixed = TRUE)
    alpha <- as.numeric(sub(".*an 'alpha_prior' of ([.0-9]+) or more.*", "\\1", msg))
    expect_s3_class(d(assurance = a, alpha_prior = alpha), "bucss_plan")
    named <- named + 1L
  }
  expect_identical(named, 40L)
})

test_that("a stated ceiling is strictly below the true one, so it plans", {
  # A ceiling a hair below .87 (the audit's construction: N = 61, whose
  # rounded-down design has 58 error df) and one a hair above it (N = 60, df
  # 58, stored as 0.87000000000000066). Both used to be stated as .87, where
  # the first refuses again (a plan needs the ceiling above the assurance)
  # and the second plans more than 2^52 per group. Both are now stated as .86,
  # in the zero-NCP message and in the print and knit_print footers, and .86
  # plans 4,700 per group (QC 2026-09-30, finding 8).
  below <- qt(.05 * (1 - (.87 - 5e-12)) / 2, 58, lower.tail = FALSE)
  above <- qt(.05 * .13 / 2, 58, lower.tail = FALSE)
  for (case in list(list(t = below, N = 61), list(t = above, N = 60))) {
    msg <- tryCatch(ss_buc_independent_t(t_observed = case$t, N = case$N, assurance = .95),
                    bucss_zero_ncp = conditionMessage)
    expect_match(msg, "usable plan is about .86;", fixed = TRUE)
    plan <- ss_buc_independent_t(t_observed = case$t, N = case$N, assurance = .86)
    expect_identical(size(plan), 4700)
    expect_equal(attr(plan, "assurance_ceiling"), .87, tolerance = 1e-10)
    expect_true(any(grepl("Largest supportable assurance: .86", capture.output(print(plan)),
                          fixed = TRUE)))
    if (requireNamespace("knitr", quietly = TRUE))
      expect_match(as.character(knitr::knit_print(plan)), "Largest supportable assurance: .86",
                   fixed = TRUE)
  }
  expect_error(ss_buc_independent_t(t_observed = below, N = 61, assurance = .87),
               class = "bucss_zero_ncp")
  # A ceiling of 1 in double precision (t = 30 on 38 df, p = 4.6e-28) is stated as
  # .99, an assurance a planner accepts, not as 1.00, which none does.
  strong <- ss_buc_independent_t(t_observed = 30, n = 20)
  expect_identical(attr(strong, "assurance_ceiling"), 1)
  expect_true(any(grepl("Largest supportable assurance: .99", capture.output(print(strong)),
                        fixed = TRUE)))
})
