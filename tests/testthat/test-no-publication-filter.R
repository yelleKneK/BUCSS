# 'alpha_prior = 1' means no publication filter, and it is computed as exactly
# 1 (Ken Kelley's decision, 2026-10-13; BUCSS 1.x and the earlier builds of
# 2.0.0 computed with .999 in its place). The string "No Publication Filter"
# is accepted as a synonym for 1, so a call can say what the number means.
#
# With no filter the publication critical value is 0, the power term of the
# truncated likelihood is 1, and the truncated distribution is the statistic's
# ordinary distribution. The adjusted noncentrality parameter is therefore the
# value at which the ordinary distribution function of the statistic, at its
# observed value, equals the assurance, and the largest supportable assurance
# is exactly 1 - p. The blocks below check every planner against that
# definition, computed here from pt(), pf(), pchisq(), and uniroot() alone
# (for sampled predictors, from a quadrature over the chi square quantiles
# rather than the negative binomial series the package sums), and against
# sizes and parameters computed on 2026-10-13 by a separate script that never
# loads BUCSS (the same definition, with the planned study's ordinary power
# analysis written out for each design). Each case also carries the parameter
# the .999 reading returned ('old'), so the tests show what moved and by how
# much, and fail if a substitute for 1 ever comes back.

npf_size <- function(x) x$value[x$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1]
npf_val <- function(x, term) x$value[x$term == term]
npf_with <- function(call, ...) {
  extra <- list(...)
  for (nm in names(extra)) call[[nm]] <- extra[[nm]]
  eval(call)
}

# The root of G(lambda) = assurance, for G decreasing in lambda.
npf_root <- function(G, assurance) {
  upper <- 1
  while (G(upper) > assurance) upper <- 2 * upper
  uniroot(function(lambda) G(lambda) - assurance, c(0, upper), tol = 1e-12)$root
}
# The ordinary distribution function of each kind of statistic at its observed
# value, as a function of the noncentrality parameter: P(|T| <= |t|) for the t
# planners (the package's two-sided construction); the F with both degrees of
# freedom and the parameter scaled by epsilon for a within-subjects effect;
# the chi square; and, for sampled predictors, the noncentral F averaged over
# the chi square (m df) that scales its parameter, lambda/N times W.
npf_t <- function(t, df) function(lambda) pt(abs(t), df, lambda) - pt(-abs(t), df, lambda)
npf_f <- function(F, df1, df2, epsilon = 1)
  function(lambda) pf(F, epsilon * df1, epsilon * df2, epsilon * lambda)
npf_chisq <- function(x, df) function(lambda) pchisq(x, df, lambda)
npf_mixture <- function(F, df1, df2, m, N) function(lambda)
  integrate(function(u) pf(F, df1, df2, (lambda / N) * qchisq(u, m)), 0, 1,
            rel.tol = 1e-11, subdivisions = 2000L)$value
# Draw to a null device: nothing is written and no display is needed.
npf_null_device <- function(code) {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  force(code)
}
npf_welch_df <- function(a, b, ratio) {
  v1 <- 1 / a
  v2 <- ratio^2 / b
  (v1 + v2)^2 / (v1^2 / (a - 1) + v2^2 / (b - 1))
}

# One prior study per planner (both predictor frames for the regression
# planners), at an assurance of .80. 'G' holds the distribution function at
# each balanced design the planner reads the prior study as (one, or the
# rounded-down and rounded-up pair); the plan reports the smaller root and the
# smaller ceiling. 'size' and 'ncp' are the separate script's values; 'old' is
# the parameter the .999 reading returned.
npf_cases <- list(
  # 50 and 55 per group: harmonic mean 52.4, read as 52 and 53 (df 102, 104).
  list(plan = quote(ss_buc_independent_t(t_observed = 3, n = c(50, 55))),
       G = list(npf_t(3, 102), npf_t(3, 104)), size = 93, ncp = 2.1326836375,
       old = 2.1326085374),
  list(plan = quote(ss_buc_paired_t(t_observed = 2.4, N = 30)),
       G = list(npf_t(2.4, 29)), size = 108, ncp = 1.4967522784, old = 1.4965077115),
  list(plan = quote(ss_buc_one_sample_t(t_observed = 2.6, N = 25)),
       G = list(npf_t(2.6, 24)), size = 72, ncp = 1.6742676240, old = 1.6740800370),
  list(plan = quote(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 2)),
       G = list(npf_t(3, npf_welch_df(40, 55, 2))), size = 91, ncp = 2.1271596016,
       old = 2.1270832724),
  # N = 46 in four groups: 11 and 12 per group (error df 40 and 44).
  list(plan = quote(ss_buc_one_way_anova(F_observed = 3.2, N = 46, levels_A = 4)),
       G = list(npf_f(3.2, 3, 40), npf_f(3.2, 3, 44)), size = 44, ncp = 3.0258326969,
       old = 3.0251610187),
  # N = 77 in six cells: 12 and 13 per cell (error df 66 and 72).
  list(plan = quote(ss_buc_factorial_anova(F_observed = 4.5, N = 77, levels_A = 2,
                                           levels_B = 3, effect = "factor_B")),
       G = list(npf_f(4.5, 2, 66), npf_f(4.5, 2, 72)), size = 36, ncp = 3.6154613159,
       old = 3.6149714506),
  # Nine nuisance parameters (125 - 6 - 110); 20 and 21 per cell.
  list(plan = quote(ss_buc_factorial_anova_general(F_observed = 7, N = 125, cells = 6,
                                                   df_numerator = 1,
                                                   df_denominator = 110)),
       G = list(npf_f(7, 1, 120 - 6 - 9), npf_f(7, 1, 126 - 6 - 9)), size = 53,
       ncp = 3.1823527995, old = 3.1818249244),
  # Three covariates; 15 and 16 per cell.
  list(plan = quote(ss_buc_ancova(F_observed = 4, N = 62, cells = 4, n_covariates = 3,
                                  df_numerator = 1)),
       G = list(npf_f(4, 1, 60 - 4 - 3), npf_f(4, 1, 64 - 4 - 3)), size = 99,
       ncp = 1.2753065955, old = 1.2744232925),
  list(plan = quote(ss_buc_rm_anova(F_observed = 5, N = 30, levels_A = 4, epsilon = .7)),
       G = list(npf_f(5, 3, 3 * 29, .7)), size = 66, ncp = 6.5014891498,
       old = 6.5009963384),
  list(plan = quote(ss_buc_rm_anova_general(F_observed = 3.5, N = 40, df_numerator = 3,
                                            epsilon = .8)),
       G = list(npf_f(3.5, 3, 3 * 39, .8)), size = 142, ncp = 3.6323786945,
       old = 3.6315716387),
  # N = 47 in three groups: 15 and 16 per group.
  list(plan = quote(ss_buc_mixed_anova(F_observed = 3.5, N = 47, levels_between = 3,
                                       levels_within = 3, effect = "interaction",
                                       epsilon = .8)),
       G = list(npf_f(3.5, 4, (45 - 3) * 2, .8), npf_f(3.5, 4, (48 - 3) * 2, .8)),
       size = 45, ncp = 5.1179841364, old = 5.1174327376),
  # N = 64 in three groups: 21 and 22 per group.
  list(plan = quote(ss_buc_mixed_anova_general(F_observed = 3, N = 64, df_numerator = 4,
                                               num_groups = 3, effect = "between_within",
                                               df_num_within = 2, epsilon = .8)),
       G = list(npf_f(3, 4, (63 - 3) * 2, .8), npf_f(3, 4, (66 - 3) * 2, .8)),
       size = 82, ncp = 3.7966853630, old = 3.7958693861),
  list(plan = quote(ss_buc_reg_coef(t_observed = 2.5, N = 100, p = 3)),
       G = list(npf_f(2.5^2, 1, 96)), size = 295, ncp = 2.6837560462, old = 2.6831340531),
  # One of five sampled predictors: m = N - p = 35.
  list(plan = quote(ss_buc_reg_coef(t_observed = -2.2, N = 40, p = 5,
                                    random_predictors = TRUE)),
       G = list(npf_mixture(2.2^2, 1, 34, 35, 40)), size = 169, ncp = 1.9637691712,
       old = 1.9628023681),
  list(plan = quote(ss_buc_R2(F_observed = 5, N = 150, p = 4)),
       G = list(npf_f(5, 4, 145)), size = 181, ncp = 10.1851391443, old = 10.1851057392),
  # All three sampled predictors: m = N - 1 = 44.
  list(plan = quote(ss_buc_R2(F_observed = 3.1, N = 45, p = 3, random_predictors = TRUE)),
       G = list(npf_mixture(3.1, 3, 41, 44, 45)), size = 177, ncp = 2.8856949924,
       old = 2.8849224015),
  list(plan = quote(ss_buc_reg_joint(F_observed = 4.5, N = 150, p = 4, p_joint = 2)),
       G = list(npf_f(4.5, 2, 145)), size = 393, ncp = 3.7099604690, old = 3.7094968134),
  # Three of six sampled predictors: m = N - 1 - (p - p_joint) = 56. The one
  # case here whose size moves: the .999 reading planned 196.
  list(plan = quote(ss_buc_reg_joint(F_observed = 3.3, N = 60, p = 6, p_joint = 3,
                                     random_predictors = TRUE)),
       G = list(npf_mixture(3.3, 3, 53, 56, 60)), size = 195, ncp = 3.5301619090,
       old = 3.5294889530),
  # A correlation samples both variables: m = N - 1 = 99, F = t^2.
  list(plan = quote(ss_buc_correlation(r_observed = .35, N = 100)),
       G = list(npf_mixture(.35^2 * 98 / (1 - .35^2), 1, 98, 99, 100)), size = 103,
       ncp = 7.9755147145, old = 7.9754198536),
  # N = 81 in two groups: 40 and 41 per group (error df 74 and 76).
  list(plan = quote(ss_buc_manova(F_observed = 3.2, N = 81, p_variables = 5)),
       G = list(npf_f(3.2, 5, 80 - 5 - 1), npf_f(3.2, 5, 82 - 5 - 1)), size = 91,
       ncp = 6.0005545722, old = 6.0003231618),
  list(plan = quote(ss_buc_chisq_diff(chisq_observed = 9, N = 250, df_full = 42,
                                      df_restricted = 43)),
       G = list(npf_chisq(9, 1)), size = 422, ncp = 4.6585969789, old = 4.6582966763))

test_that("the cases cover every planner", {
  planners <- unique(vapply(npf_cases, function(case) as.character(case$plan[[1]]),
                            character(1)))
  expect_setequal(planners, vapply(BUCSS:::.BUCSS_DESIGN_SPECS, `[[`, character(1),
                                   "planner"))
  expect_length(planners, 18L)
})

test_that("at alpha_prior = 1 the parameter solves the statistic's ordinary distribution", {
  for (case in npf_cases) {
    label <- paste(deparse(case$plan, width.cutoff = 500L), collapse = " ")
    plan <- npf_with(case$plan, alpha_prior = 1)
    # The root at each balanced design; the plan reports the smaller.
    roots <- vapply(case$G, npf_root, numeric(1), assurance = .80)
    expect_equal(npf_val(plan, "ncp_adjusted"), min(roots), tolerance = 1e-7, label = label)
    # The separate script's values, and the value .999 gave, which differs.
    expect_equal(npf_val(plan, "ncp_adjusted"), case$ncp, tolerance = 1e-7, label = label)
    expect_identical(npf_size(plan), case$size, label = label)
    expect_gt(npf_val(plan, "ncp_adjusted"), case$old + 2e-6, label = label)
    # The echoed value is the number 1, the value computed with.
    expect_identical(npf_val(plan, "alpha_prior"), 1, label = label)
  }
})

test_that("at alpha_prior = 1 the largest supportable assurance is exactly 1 - p", {
  # At a zero noncentrality parameter each distribution function above is 1
  # minus the prior study's p value at that design; the plan stores the
  # smaller of the two designs' values.
  for (case in npf_cases) {
    label <- paste(deparse(case$plan, width.cutoff = 500L), collapse = " ")
    plan <- npf_with(case$plan, alpha_prior = 1)
    one_minus_p <- min(vapply(case$G, function(G) G(0), numeric(1)))
    expect_equal(attr(plan, "assurance_ceiling"), one_minus_p, tolerance = 1e-12,
                 label = label)
    expect_lt(attr(plan, "assurance_ceiling"), 1)
  }
  # Written from the p value itself, for one planner of each family.
  ceiling_of <- function(x) attr(x, "assurance_ceiling")
  expect_equal(ceiling_of(ss_buc_paired_t(t_observed = 2.4, N = 30, alpha_prior = 1)),
               1 - 2 * pt(2.4, 29, lower.tail = FALSE), tolerance = 1e-12)
  expect_equal(ceiling_of(ss_buc_rm_anova_general(F_observed = 3.5, N = 40, df_numerator = 3,
                                                  epsilon = .8, alpha_prior = 1)),
               1 - pf(3.5, .8 * 3, .8 * 117, lower.tail = FALSE), tolerance = 1e-12)
  expect_equal(ceiling_of(ss_buc_chisq_diff(chisq_observed = 9, N = 250, df_full = 42,
                                            df_restricted = 43, alpha_prior = 1)),
               1 - pchisq(9, 1, lower.tail = FALSE), tolerance = 1e-12)
  # The same in both predictor frames: at a zero effect the mixture is the
  # central F.
  expect_equal(ceiling_of(ss_buc_R2(F_observed = 3.1, N = 45, p = 3, alpha_prior = 1,
                                    random_predictors = TRUE)),
               1 - pf(3.1, 3, 41, lower.tail = FALSE), tolerance = 1e-12)
})

test_that("the endpoint is monotone: 1 is at least .9999, which is at least .999", {
  # A more lenient publication threshold leaves less bias to correct, so the
  # adjusted parameter and the ceiling rise with 'alpha_prior' and the plan
  # does not grow. Exactly 1 is the end of that sequence, not a value apart
  # from it: every prior study above moves by more than the root finder's
  # tolerance at each step (the smallest relative step from .9999 to 1 is
  # 3e-7, for the R squared prior with fixed predictors).
  for (case in npf_cases) {
    label <- paste(deparse(case$plan, width.cutoff = 500L), collapse = " ")
    plans <- lapply(c(.999, .9999, 1), function(a) npf_with(case$plan, alpha_prior = a))
    ncp <- vapply(plans, npf_val, numeric(1), term = "ncp_adjusted")
    expect_gt(ncp[3], ncp[2], label = label)
    expect_gt(ncp[2], ncp[1], label = label)
    # .999 is the value the old reading returned for 1.
    expect_equal(ncp[1], case$old, tolerance = 1e-7, label = label)
    ceilings <- vapply(plans, attr, numeric(1), "assurance_ceiling")
    expect_gt(ceilings[3], ceilings[2], label = label)
    expect_gt(ceilings[2], ceilings[1], label = label)
    sizes <- vapply(plans, npf_size, numeric(1))
    expect_lte(sizes[3], sizes[2], label = label)
    expect_lte(sizes[2], sizes[1], label = label)
  }
})

test_that("a weak prior plans exactly up to 1 - p, and a statistic of 0 never plans", {
  # Weak prior studies (the pilot-study use of 'alpha_prior = 1'), where the
  # truncation point matters most. Each plans at an assurance below 1 - p,
  # with the parameter the ordinary distribution gives (separate script:
  # sizes 176, 365, 1,887, and 9,780), and is refused, for a zero
  # noncentrality parameter, at an assurance at or just above 1 - p.
  weak <- list(
    list(plan = quote(ss_buc_paired_t(t_observed = .9, N = 12)), G = npf_t(.9, 11),
         assurance = .5, size = 176, ncp = 0.7357433086,
         one_minus_p = 1 - 2 * pt(.9, 11, lower.tail = FALSE)),
    list(plan = quote(ss_buc_one_way_anova(F_observed = .9, N = 46, levels_A = 4)),
         G = npf_f(.9, 3, 40), assurance = .5, size = 365, ncp = 0.3552637421,
         one_minus_p = 1 - pf(.9, 3, 40, lower.tail = FALSE)),
    list(plan = quote(ss_buc_correlation(r_observed = .1, N = 60)),
         G = npf_mixture(.1^2 * 58 / (1 - .1^2), 1, 58, 59, 60), assurance = .5,
         size = 1887, ncp = 0.2502304093,
         one_minus_p = 1 - pf(.1^2 * 58 / (1 - .1^2), 1, 58, lower.tail = FALSE)),
    list(plan = quote(ss_buc_chisq_diff(chisq_observed = 1.2, N = 400, df_full = 20,
                                        df_restricted = 24)),
         G = npf_chisq(1.2, 4), assurance = .1, size = 9780, ncp = 0.4881543651,
         one_minus_p = 1 - pchisq(1.2, 4, lower.tail = FALSE)))
  for (case in weak) {
    label <- paste(deparse(case$plan, width.cutoff = 500L), collapse = " ")
    stopifnot(case$assurance < case$one_minus_p, case$one_minus_p < .62)
    # An assurance below .5 warns, as it always has; that is not under test.
    plan <- suppressWarnings(npf_with(case$plan, alpha_prior = 1,
                                      assurance = case$assurance))
    expect_equal(npf_val(plan, "ncp_adjusted"), npf_root(case$G, case$assurance),
                 tolerance = 1e-6, label = label)
    expect_equal(npf_val(plan, "ncp_adjusted"), case$ncp, tolerance = 1e-6, label = label)
    expect_identical(npf_size(plan), case$size, label = label)
    expect_equal(attr(plan, "assurance_ceiling"), case$one_minus_p, tolerance = 1e-12,
                 label = label)
    # Just below 1 - p a plan is still issued; at 1 - p, and above, none is.
    below <- suppressWarnings(npf_with(case$plan, alpha_prior = 1,
                                       assurance = case$one_minus_p - 1e-4))
    expect_gt(npf_val(below, "ncp_adjusted"), 0, label = label)
    for (a in c(case$one_minus_p + 1e-9, case$one_minus_p + 1e-4, .80))
      expect_error(suppressWarnings(npf_with(case$plan, alpha_prior = 1, assurance = a)),
                   class = "bucss_zero_ncp", label = label)
  }
  # While 1 was read as .999, the ceiling was 1 - p/.999, so the last of the
  # assurances that now plan was refused: the paired t above, at 1 - p - 1e-4.
  p <- 2 * pt(.9, 11, lower.tail = FALSE)
  stopifnot(1 - p - 1e-4 > 1 - p / .999)
  expect_error(ss_buc_paired_t(t_observed = .9, N = 12, alpha_prior = .999,
                               assurance = 1 - p - 1e-4), class = "bucss_zero_ncp")

  # A statistic of exactly 0 has p = 1 and supports no assurance at all: with
  # no publication filter it is the one statistic refused as nonsignificant
  # (every planner is checked in test-guards-2026-09.R), at any assurance.
  for (a in c(.01, .5, .8)) {
    expect_error(suppressWarnings(ss_buc_paired_t(t_observed = 0, N = 30, alpha_prior = 1,
                                                  assurance = a)),
                 "is too small to plan from even with no publication filter.",
                 fixed = TRUE, class = "bucss_nonsignificant")
    expect_error(suppressWarnings(ss_buc_one_way_anova(F_observed = 0, N = 46, levels_A = 4,
                                                       alpha_prior = 1, assurance = a)),
                 class = "bucss_nonsignificant")
  }
})

test_that("\"No Publication Filter\" is 1, for every planner", {
  for (case in npf_cases) {
    label <- paste(deparse(case$plan, width.cutoff = 500L), collapse = " ")
    one <- npf_with(case$plan, alpha_prior = 1)
    named <- npf_with(case$plan, alpha_prior = "No Publication Filter")
    expect_identical(named, one, label = label)
    # The echoed row is the number 1: the value column is numeric throughout.
    expect_true(is.numeric(named$value), label = label)
    expect_identical(npf_val(named, "alpha_prior"), 1, label = label)
  }
})

test_that("the accepted spellings of \"No Publication Filter\", and the refused strings", {
  one <- ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = 1)
  # Case does not matter, surrounding whitespace is trimmed, and a run of
  # spaces, underscores, or hyphens between the three words is read alike.
  accepted <- c("No Publication Filter", "no publication filter", "NO PUBLICATION FILTER",
                "nO pUBLICATION fILTER", "  No Publication Filter  ",
                "\tNo Publication Filter\n", "No_Publication_Filter",
                "no-publication-filter", "No  Publication   Filter",
                "no__publication--filter", "No - Publication _ Filter",
                "no_publication filter")
  for (spelling in accepted)
    expect_identical(ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = spelling), one,
                     label = deparse(spelling))

  # Anything else is refused, with one message that names the numeric range
  # and the string.
  refusal <- paste0("'alpha_prior' must be a number above 0 and at most 1, or the string ",
                    "\"No Publication Filter\", which means 1 (any result would have been ",
                    "published, whatever its p value).")
  refused <- list("none", "No Filter", "no publication bias", "NoPublicationFilter",
                  "No Publication Filters", "No.Publication.Filter", "_No Publication Filter",
                  "No Publication Filter.", "1", ".05", "", " ", NA_character_,
                  character(0), c("No Publication Filter", "No Publication Filter"))
  for (bad in refused) {
    err <- tryCatch(ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = bad),
                    error = function(e) e)
    expect_s3_class(err, "error")
    expect_identical(conditionMessage(err), refusal, label = paste(deparse(bad), collapse = ""))
    # call. = FALSE: the internal helper is not named.
    expect_null(conditionCall(err))
  }
  # The numeric range is unchanged: above 0 and at most 1.
  for (bad in list(0, -.05, 1.0001, 5))
    expect_error(ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = bad),
                 "There is a problem with 'alpha_prior'", fixed = TRUE)
  # Only 'alpha_prior' takes the string, and nothing that is neither a number
  # nor a string is read as one.
  expect_error(ss_buc_paired_t(t_observed = 3, N = 30,
                               alpha_planned = "No Publication Filter"),
               "'alpha_planned' must be a single finite numeric value.", fixed = TRUE)
  for (bad in list(TRUE, NULL, factor("No Publication Filter"), list(1)))
    expect_error(ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = bad),
                 "'alpha_prior' must be a single finite numeric value.", fixed = TRUE)
})

test_that("a plan made with \"No Publication Filter\" works everywhere a plan does", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = "No Publication Filter")
  one <- ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = 1)
  # Printing, the verbs, and the sentence read the numeric row.
  expect_identical(capture.output(print(plan)), capture.output(print(one)))
  expect_match(capture.output(print(plan)), "^ alpha_prior +1 *$", all = FALSE)
  expect_identical(glance(plan)$alpha_prior, 1)
  expect_identical(tidy(plan), tidy(one))
  expect_match(planning_sentence(plan), "with no publication bias assumed", fixed = TRUE)
  # The replays rebuild the planner call from the numeric row.
  audit <- ss_buc_sensitivity(plan, true_ncp = 2, replications = 50, seed = 2017)
  expect_identical(audit, ss_buc_sensitivity(one, true_ncp = 2, replications = 50,
                                             seed = 2017))
  drawn <- npf_null_device(plot(plan, assurance = c(.6, .8, .9)))
  expect_identical(drawn$sample_size[2], npf_size(plan))
  expect_identical(drawn$ncp_adjusted[2], npf_val(plan, "ncp_adjusted"))
  expect_equal(attr(drawn, "assurance_ceiling"), 1 - 2 * pt(3, 29, lower.tail = FALSE),
               tolerance = 1e-12)
  # The whole ceiling figure: plans up to the ceiling 1 - p, none beyond it.
  curve_values <- npf_null_device(plot(plan, which = "size", n_points = 12))
  expect_false(anyNA(curve_values$sample_size))
  expect_true(all(curve_values$assurance < 1 - 2 * pt(3, 29, lower.tail = FALSE)))
})

test_that("with no publication filter every simulated prior study is published", {
  # ss_buc_sensitivity() at 'alpha_prior = 1': the publication rate is exactly
  # 1 and no draw is discarded, in every sampling family.
  for (family in c("t", "f", "f_from_t", "correlation", "f_random", "f_from_t_random",
                   "chisq"))
    expect_identical(BUCSS:::.prior_publication(list(family = family), c(2, 30, 31, 34), 3, 1),
                     list(crit = -Inf, rate = 1))
  families <- unlist(lapply(BUCSS:::.BUCSS_DESIGN_SPECS, `[[`, "family"), use.names = FALSE)
  expect_setequal(families, c("t", "f", "f_from_t", "correlation", "f_random",
                              "f_from_t_random", "chisq"))
  # Below 1 there is a threshold, and .999 is not "no filter".
  below <- BUCSS:::.prior_publication(list(family = "t"), 29, 2, .999)
  expect_equal(below$crit, qt(1 - .999 / 2, 29))
  expect_lt(below$rate, 1)

  # The audit itself, reproduced from the definition. With seed 2017 the
  # simulated studies are the first 300 draws of rt(), none discarded; a study
  # is refused exactly when its largest supportable assurance, 1 - p, is at or
  # below the requested .80, that is, when its two-sided p is .20 or more.
  plan <- ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = 1)
  audit <- ss_buc_sensitivity(plan, true_ncp = 2, replications = 300, seed = 2017)
  row <- function(term) audit$value[audit$term == term]
  expect_identical(row("publication_rate"), 1)
  draws <- withr::with_seed(2017, rt(1000, df = 29, ncp = 2))[1:300]
  p <- 2 * pt(abs(draws), 29, lower.tail = FALSE)
  stopifnot(sum(p >= .20) > 30, sum(p >= .20) < 150)
  expect_identical(row("refusal_rate"), mean(p >= .20))
  # Each issued plan is the planner's own plan from that statistic.
  issued <- vapply(draws[p < .20], function(t)
    npf_size(ss_buc_paired_t(t_observed = t, N = 30, alpha_prior = 1)), numeric(1))
  expect_identical(c(row("size_q25"), row("size_median"), row("size_q75")),
                   quantile(issued, c(.25, .5, .75), names = FALSE))
  # With .999 in place of 1, about one study in a thousand was discarded
  # before planning; under this seed and noncentrality parameter none of the
  # first 300 draws is that small, so the audits differ only where the plans
  # differ.
  stopifnot(all(abs(draws) > qt(1 - .999 / 2, 29)))
})
