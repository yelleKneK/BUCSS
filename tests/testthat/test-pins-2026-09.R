# Pins the 2026-09-28 deep QC audit asked for. Each block names its finding,
# and where the finding reported one-line mutants that passed the whole suite,
# the block names each mutant it kills and says in words what the mutant
# changes. The audit ran more than one numbered series, so a label carries its
# series: H-M10 and the like are auditor H's (the evidence folder's
# mutants/mutants_run.txt), B-M04 and the like are auditor B's (finding B-03 in
# the findings file), and C02 is from mutants/mutants.txt. The M-04 and CC-03
# blocks pin values no mutant was reported against.
#
# Every expected value is what the package returns now, and was also
# confirmed without BUCSS: either it is an auditor's independent value
# (finding id named) or it was computed in a separate script that never loads
# BUCSS, from the method's definition (bisection for the corrected parameter,
# a linear walk for the planned n, the two-sided rounding rule re-derived from
# its description), marked "independent check 2026-09-28". The two agree to
# about 1e-9 when both use stats::pf(); a root from an exact series for the
# noncentral F (Poisson-weighted central beta tails) differs by up to about
# 7e-8 (relative) in ncp_adjusted, which is the accuracy of R's noncentral
# pf(), and gives the same sizes, with powers within about 2e-8. The tolerances
# leave room only for platform differences in pt(), pf(), and pchisq().

size <- function(x) x$value[x$term %in% BUCSS:::.BUCSS_SIZE_TERMS]
val <- function(x, term) x$value[x$term == term]

expect_plan <- function(x, n, ncp, power = NULL, tol = 1e-8) {
  expect_identical(size(x), n)
  expect_equal(val(x, "ncp_adjusted"), ncp, tolerance = tol)
  if (!is.null(power)) expect_equal(val(x, "actual_power"), power, tolerance = tol)
}

# First-principles power of the planned study for the reported pair: the
# returned n at the returned ncp, rescaled from the prior per-group size
# 'n_prior' the ncp came from.
power_f <- function(n, n_prior, ncp, df1, df2, alpha = .05) {
  pf(qf(1 - alpha, df1, df2), df1, df2, ncp = (n / n_prior) * ncp, lower.tail = FALSE)
}
power_t2 <- function(df, ncp, alpha = .05) {
  crit <- qt(1 - alpha / 2, df)
  pt(crit, df, ncp, lower.tail = FALSE) + pt(-crit, df, ncp)
}

test_that("H-03: actual_power follows the displayed ncp where branch orderings disagree", {
  # In each case the round-down branch has the smaller ncp (the one shown)
  # while the round-up branch has the smaller ncp per prior participant, so a
  # planner that picked the branch by ncp / prior n (mutant H-M10) would report
  # the other branch's power (shown in the comments; independent check
  # 2026-09-28). The first-principles line recomputes each pinned power from
  # the displayed ncp and the round-down prior n.

  # One-way, F = 4, N = 61, four groups: prior n 15 and 16 (H-03 found .82699;
  # the other branch gives .81021).
  r <- ss_buc_one_way_anova(F_observed = 4, N = 61, levels_A = 4, assurance = .5)
  expect_plan(r, 38, 4.7080325303, 0.8269883262)
  expect_equal(val(r, "actual_power"),
               power_f(38, 15, val(r, "ncp_adjusted"), 3, 4 * 38 - 4))

  # Independent t, n = c(50, 55): harmonic mean 52.4, prior n 52 and 53
  # (other branch .8008896678). This is also the plan in the
  # ss_buc_sensitivity() help page (CC-03: 233, ncp 1.338).
  r <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .8)
  expect_plan(r, 233, 1.3380186786, 0.8068829772)
  expect_equal(val(r, "actual_power"),
               power_t2(2 * 233 - 2, sqrt(233 / 52) * val(r, "ncp_adjusted")))

  # Two-way, 3 x 4, N = 41: prior n 3 and 4, at 'alpha_prior = 1', computed
  # as exactly 1 since 2026-10-13: n 9, ncp 11.9158701383, power .9838166011;
  # other branch .9255350521, which is also what mutant B-M07, taking the
  # two-way planner's actual_power from the other branch, reports (independent
  # check 2026-10-13: the root of pf(2.9, 6, 24, ncp) = .5 from pf() alone).
  # While 1 was computed as .999, B-03's verifier found n 9, ncp 11.9158296797,
  # power .9838162280, and other branch .9255340798.
  r <- ss_buc_factorial_anova(F_observed = 2.9, N = 41, levels_A = 3, levels_B = 4,
                              effect = "interaction", alpha_prior = 1,
                              alpha_planned = .01, assurance = .5, desired_power = .9)
  expect_plan(r, 9, 11.9158701383, 0.9838166011)
  expect_equal(pf(2.9, 6, 3 * 12 - 12, ncp = val(r, "ncp_adjusted")), .5, tolerance = 1e-8)
  expect_equal(val(r, "actual_power"),
               power_f(9, 3, val(r, "ncp_adjusted"), 6, 12 * 9 - 12, alpha = .01))

  # General between subjects with two nuisance parameters, and the ANCOVA
  # that is the same computation: N = 121, three cells, prior n 40 and 41
  # (B-03's verifier found n 157, ncp 2.530713, power .810023; other branch
  # .8019732672).
  r <- ss_buc_factorial_anova_general(F_observed = 6, N = 121, cells = 3,
                                      df_numerator = 2, df_denominator = 116)
  expect_plan(r, 157, 2.5307128019, 0.8100230523)
  expect_equal(val(r, "actual_power"),
               power_f(157, 40, val(r, "ncp_adjusted"), 2, 3 * 157 - 3 - 2))
  a <- ss_buc_ancova(F_observed = 6, N = 121, cells = 3, n_covariates = 2)
  expect_plan(a, 157, 2.5307128019, 0.8100230523)

  # Split plot, N = 61, two groups: prior n 30 and 31 (other branch
  # .8009631178).
  r <- ss_buc_mixed_anova(F_observed = 4, N = 61, levels_between = 2,
                          levels_within = 3, effect = "within", assurance = .5)
  expect_plan(r, 412, 0.7096438949, 0.8040385549)
  expect_equal(val(r, "actual_power"),
               power_f(412, 30, val(r, "ncp_adjusted"), 2, 2 * (412 - 1) * 2))

  # Split plot general, between x within, N = 91, three groups: prior n 30
  # and 31 (W-03 found n 111, ncp 2.710585; other branch .8029680924).
  r <- ss_buc_mixed_anova_general(F_observed = 6, N = 91, df_numerator = 2,
                                  num_groups = 3, effect = "between_within",
                                  df_num_within = 2)
  expect_plan(r, 111, 2.7105848403, 0.8147553949)
  expect_equal(val(r, "actual_power"),
               power_f(111, 30, val(r, "ncp_adjusted"), 2, (3 * 111 - 3) * 2))

  # Two-group multivariate, N = 81, three outcomes: prior n 40 and 41 (M-04's
  # verifier found n 51, ncp 9.1425607265, power .8116580900; other branch
  # .8025675559).
  r <- ss_buc_manova(F_observed = 8, N = 81, p_variables = 3, assurance = .9)
  expect_plan(r, 51, 9.1425607265, 0.8116580900)
  expect_equal(val(r, "actual_power"),
               power_f(51, 40, val(r, "ncp_adjusted"), 3, 2 * 51 - 3 - 1))
})

test_that("H-04: the independent t planned study has 2n - 2 degrees of freedom", {
  # Mutant H-M17 (planned df 2n - 1) moves these to 7, 2, and 4 (H-04;
  # independent check 2026-09-28 gives 8, 3, and 5, and 7 and 2 under the
  # mutant).
  r <- ss_buc_independent_t(t_observed = 6, n = 20, assurance = .5, desired_power = .9)
  expect_plan(r, 8, 5.9562942527, 0.9380150514)
  expect_equal(val(r, "actual_power"),
               power_t2(2 * 8 - 2, sqrt(8 / 20) * val(r, "ncp_adjusted")))
  expect_lt(power_t2(2 * 7 - 2, sqrt(7 / 20) * val(r, "ncp_adjusted")), .9)

  expect_plan(ss_buc_independent_t(t_observed = 10, n = 10, assurance = .5),
              3, 9.8270793368, 0.9745124116)
  expect_plan(ss_buc_independent_t(t_observed = 10, n = 20, assurance = .8,
                                   desired_power = .9),
              5, 8.6521407970, 0.9645947641)

  # The paired t planned study has n - 1.
  p <- ss_buc_paired_t(t_observed = 6, N = 20, assurance = .5, desired_power = .9)
  expect_plan(p, 9, 5.9082490495, 0.9327564026)
  expect_equal(val(p, "actual_power"),
               power_t2(9 - 1, sqrt(9 / 20) * val(p, "ncp_adjusted")))
  expect_lt(power_t2(8 - 1, sqrt(8 / 20) * val(p, "ncp_adjusted")), .9)
})

test_that("H-05 and W-03: split-plot planners with an unbalanced prior N", {
  # N = 61 and 91 with two groups do not divide evenly, so the two rounding
  # branches differ. Mutant H-M11 (scale the ncp by the unrounded N instead of
  # the branch's n) moves 412 to 415 (H-05; independent check 2026-09-28).
  within_61 <- ss_buc_mixed_anova(F_observed = 4, N = 61, levels_between = 2,
                                  levels_within = 3, effect = "within", assurance = .5)
  expect_plan(within_61, 412, 0.7096438949, 0.8040385549)
  # The stored ceiling is the smaller one, the round-down branch's.
  expect_equal(attr(within_61, "assurance_ceiling"),
               1 - pf(4, 2, (30 * 2 - 2) * 2, lower.tail = FALSE) / .05, tolerance = 1e-10)
  # With two groups and three occasions the interaction has the same df as
  # the within effect, so the same plan.
  expect_plan(ss_buc_mixed_anova(F_observed = 4, N = 61, levels_between = 2,
                                 levels_within = 3, effect = "interaction",
                                 assurance = .5),
              412, 0.7096438949, 0.8040385549)
  for (eff in c("within", "interaction"))
    expect_plan(ss_buc_mixed_anova(F_observed = 4, N = 91, levels_between = 2,
                                   levels_within = 3, effect = eff, assurance = .5),
                489, 0.8998609509, 0.8055550558)

  # Between effect (W-03: n 411, branches 408 / 411, ncp .579611, ceiling
  # .790536).
  between <- ss_buc_mixed_anova(F_observed = 7, N = 61, levels_between = 2,
                                levels_within = 3, effect = "between", assurance = .7)
  expect_plan(between, 411, 0.5796109655, 0.8036294659)
  expect_equal(attr(between, "assurance_ceiling"), 0.7905360266, tolerance = 1e-9)
  # Three groups, four occasions, N = 47 (W-03: n 26, ncp 8.565264).
  expect_plan(ss_buc_mixed_anova(F_observed = 4, N = 47, levels_between = 3,
                                 levels_within = 4, effect = "interaction"),
              26, 8.5652642384, 0.8264203714)

  # The general split-plot planner, all three effects.
  expect_plan(ss_buc_mixed_anova_general(F_observed = 4, N = 61, df_numerator = 2,
                                         num_groups = 2, effect = "within_only",
                                         assurance = .5),
              412, 0.7096438949, 0.8040385549)
  expect_plan(ss_buc_mixed_anova_general(F_observed = 7, N = 61, df_numerator = 1,
                                         num_groups = 2, effect = "between_only",
                                         assurance = .7),
              411, 0.5796109655, 0.8036294659)
  # between_within (W-03: n 111, ncp 2.710585, pinned above; and n 27, ncp
  # 15.588020 for a balanced six-df interaction).
  expect_plan(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 6,
                                         num_groups = 3, effect = "between_within",
                                         df_num_within = 3),
              27, 15.5880196631, 0.8011680112)
})

test_that("B-03: unbalanced, nuisance-df, and two-alpha paths of the between-subjects planners", {
  # Each case kills one of auditor B's seven mutants that survived (B-03; the
  # verifier's own implementation and independent check 2026-09-28 agree).
  # The seventh, B-M07, is killed by the two-way case under H-03 above.

  # Nuisance df carried into the planned study at small n: general, N = 28,
  # six cells, three nuisance parameters. Mutant B-M04 (planned df ignores
  # nuisance_df) gives n 2 with power .9408.
  r <- ss_buc_factorial_anova_general(F_observed = 37.82, N = 28, cells = 6,
                                      df_numerator = 1, df_denominator = 19,
                                      assurance = .5)
  expect_plan(r, 3, 36.3405647838, 0.9959600640)
  expect_equal(val(r, "actual_power"),
               power_f(3, 4, val(r, "ncp_adjusted"), 1, 6 * 3 - 6 - 3))

  # The round-down prior df also carries the nuisance parameters (mutant
  # B-M16, whose round-down prior df ignore them, gives n 10), and
  # actual_power is the displayed branch's (mutant B-M13, which takes the
  # general planner's actual_power from the other branch, gives .9609).
  expect_plan(ss_buc_factorial_anova_general(F_observed = 7.006, N = 28, cells = 6,
                                             df_numerator = 3, df_denominator = 16),
              15, 3.1688137973, 0.8154749525)

  # A case where the round-up branch sets n (mutant B-M11, which takes n from
  # the round-down branch alone, gives 42), and one where the round-down
  # branch alone would plan 2.
  expect_plan(ss_buc_one_way_anova(F_observed = 4.2, N = 151, levels_A = 5,
                                   assurance = .5, desired_power = .9),
              43, 11.3652625080, 0.9104546150)
  expect_plan(ss_buc_one_way_anova(F_observed = 30.357, N = 22, levels_A = 4,
                                   assurance = .5),
              3, 85.5459642713, 0.9981570015)

  # alpha_prior different from alpha_planned: mutant B-M12 (the round-down
  # critical value uses alpha_planned) gives 1405.
  expect_plan(ss_buc_one_way_anova(F_observed = 4.2, N = 62, levels_A = 5,
                                   alpha_prior = .1, alpha_planned = .01,
                                   assurance = .5, desired_power = .9),
              24, 12.1309954513, 0.9335839098)

  # The stored ceiling is the smaller branch's (mutant B-M09, which stores the
  # larger, gives .7685386387 here and .8054 for F = 6, N = 13).
  expect_equal(attr(ss_buc_one_way_anova(F_observed = 4, N = 61, levels_A = 4,
                                         assurance = .5), "assurance_ceiling"),
               0.7619448215, tolerance = 1e-9)
})

test_that("H-10: alpha_prior = 1 is exactly 1, pinned where the difference shows", {
  # A weak prior (the pilot-study use of alpha_prior = 1), where the
  # truncation point matters. With no publication filter the adjusted
  # parameter solves P(|T| <= 1.5 | df 18, ncp) = .8: independent check
  # 2026-10-13 (base R only) gives 304 per group and .5094066426, power
  # .8007379476, and a ceiling of 1 - p = .849049547816. Until 2026-10-13 a 1
  # was computed as .999, which gives 305 and .5084746316 (independent check
  # 2026-09-28), and mutant H-M13 (.99 instead of .999) gave 315 and
  # .4999740138; any substitute for 1 moves this plan.
  one <- ss_buc_independent_t(t_observed = 1.5, n = 10, alpha_prior = 1, assurance = .8)
  expect_plan(one, 304, 0.5094066426, 0.8007379476, tol = 1e-7)
  expect_equal(attr(one, "assurance_ceiling"), 1 - 2 * pt(1.5, 18, lower.tail = FALSE),
               tolerance = 1e-12)
  by_hand <- uniroot(function(ncp) pt(1.5, 18, ncp) - pt(-1.5, 18, ncp) - .8, c(0, 5),
                     tol = 1e-12)$root
  expect_equal(val(one, "ncp_adjusted"), by_hand, tolerance = 1e-7)
  # .999 is a different, slightly stricter filter, with its own plan (the one
  # 'alpha_prior = 1' used to return), and .9999 lies between the two.
  point999 <- ss_buc_independent_t(t_observed = 1.5, n = 10, alpha_prior = .999,
                                   assurance = .8)
  expect_plan(point999, 305, 0.5084746316)
  point9999 <- ss_buc_independent_t(t_observed = 1.5, n = 10, alpha_prior = .9999,
                                    assurance = .8)
  expect_plan(point9999, 304, 0.5093135513, tol = 1e-7)
  expect_gt(val(one, "ncp_adjusted"), val(point9999, "ncp_adjusted"))
  expect_gt(val(point9999, "ncp_adjusted"), val(point999, "ncp_adjusted"))
  expect_gt(attr(one, "assurance_ceiling"), attr(point9999, "assurance_ceiling"))
  expect_gt(attr(point9999, "assurance_ceiling"), attr(point999, "assurance_ceiling"))
  # ...and alpha_prior is echoed as entered.
  expect_identical(val(one, "alpha_prior"), 1)
})

test_that("M-04: MANOVA and chi square difference against independent references", {
  # M-04's verifier reproduced these to 10 decimals with its own
  # implementation (uniroot at tol 1e-13 plus a linear n search); independent
  # check 2026-09-28 agrees.
  m <- ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3)
  expect_plan(m, 310, 1.4175846301, 0.8005703109)
  x <- ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42, df_restricted = 43)
  expect_plan(x, 726, 2.7050006153, 0.8003226100)
  expect_equal(val(x, "actual_power"),
               pchisq(qchisq(.95, 1), 1, ncp = (726 / 250) * val(x, "ncp_adjusted"),
                      lower.tail = FALSE))
  expect_lt(pchisq(qchisq(.95, 1), 1, ncp = (725 / 250) * val(x, "ncp_adjusted"),
                   lower.tail = FALSE), .8)
  # The odd-N case is pinned under H-03 above (n 51, branches 51 and 50).
})

test_that("CC-03: every help-page example of the new planners is pinned", {
  # MANOVA: the F form and the equivalent T squared form.
  expect_plan(ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3,
                            alpha_prior = .05, alpha_planned = .05, assurance = .80,
                            desired_power = .80),
              310, 1.4175846301, 0.8005703109)
  expect_plan(ss_buc_manova(T2_observed = 4.5 * 3 * 78 / 76, N = 80, p_variables = 3),
              310, 1.4175846301, 0.8005703109)
  # ANCOVA: the omnibus example and the contrast example (independent check
  # 2026-09-28).
  expect_plan(ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2,
                            alpha_prior = .05, alpha_planned = .05, assurance = .80,
                            desired_power = .80),
              154, 2.5307128019, 0.8019478588)
  expect_plan(ss_buc_ancova(F_observed = 9, N = 120, cells = 3, n_covariates = 2,
                            df_numerator = 1),
              173, 1.8288756436, 0.8015673583)
  # The plan in the ss_buc_sensitivity() examples.
  expect_plan(ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .80),
              233, 1.3380186786, 0.8068829772)
  # The correlation and chi square examples' ncp_adjusted, unpinned before.
  # The correlation reference computes the mixture by quadrature over the
  # predictor's chi square, not by the package's series (independent check
  # 2026-09-28: 122, 6.71352534, .8016708332).
  expect_plan(ss_buc_correlation(r_observed = .35, N = 100, alpha_prior = .05,
                                 alpha_planned = .05, assurance = .80,
                                 desired_power = .80),
              122, 6.7135253503, 0.8016708335)
  expect_plan(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42,
                                df_restricted = 43, alpha_prior = .05,
                                alpha_planned = .05, assurance = .80,
                                desired_power = .80),
              726, 2.7050006153, 0.8003226100)
  # The five Welch examples are pinned in test-engine-fixes.R.
})

test_that("CC-01: the Welch planned df depend on the sd ratio", {
  # Mutant C02 computes the planned study's df as if the ratio were 1, which
  # moves 236 to 235 at sd_ratio = 2 and 81 to 80 below (CC-01; independent
  # check 2026-09-28 gives 236 and 81, and 235 and 80 under the mutant).
  ws_df <- function(a, b, r) {
    v1 <- 1 / a
    v2 <- r^2 / b
    (v1 + v2)^2 / (v1^2 / (a - 1) + v2^2 / (b - 1))
  }
  r <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 2)
  expect_plan(r, 236, 1.3100248639, 0.8011641987)
  scale <- sqrt((1 / 40 + 4 / 55) / ((1 + 4) / 236))
  expect_equal(val(r, "actual_power"),
               power_t2(ws_df(236, 236, 2), scale * val(r, "ncp_adjusted")))
  expect_equal(attr(r, "df_error"), ws_df(236, 236, 2))
  expect_plan(ss_buc_welch_t(t_observed = 3.2, n_1 = 20, n_2 = 30, sd_ratio = 2),
              81, 1.6493043975, 0.8048639379)
})

# Pins added 2026-10-12, from the triage of the unfinished 2026-09-29 review.
# The same rule holds: every expected value was first computed without BUCSS
# (bisection for the corrected parameter and a linear walk for the planned n,
# from the method's definition; scratch script recompute_pins.R of that
# date), and agrees with what the package returns to the digits shown.

test_that("WEL-4: the Welch planned study's critical value uses alpha_planned", {
  # The surviving mutant builds the planned study's critical t from
  # 'alpha_prior' in place of 'alpha_planned'. No Welch test set the two
  # apart, so it passed the whole suite while returning 222 for the first
  # call below and 107 for the second (independent check 2026-10-12 gives 331
  # and 136, and 222 and 107 under the mutant).
  ws_df <- function(a, b, r) {
    v1 <- 1 / a
    v2 <- r^2 / b
    (v1 + v2)^2 / (v1^2 / (a - 1) + v2^2 / (b - 1))
  }
  scale <- function(n) sqrt((1 / 40 + 1.5^2 / 55) / ((1 + 1.5^2) / n))
  strict <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5,
                           alpha_planned = .01)
  expect_plan(strict, 331, 1.3246896477, 0.8013231680)
  expect_equal(val(strict, "actual_power"),
               power_t2(ws_df(331, 331, 1.5), scale(331) * val(strict, "ncp_adjusted"), .01))
  expect_lt(power_t2(ws_df(330, 330, 1.5), scale(330) * val(strict, "ncp_adjusted"), .01), .8)
  # alpha_planned leaves the correction alone: the parameter is the default
  # call's (pinned in test-engine-fixes.R).
  expect_identical(val(strict, "ncp_adjusted"),
                   val(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5),
                       "ncp_adjusted"))
  # The other direction: a lenient alpha_prior with the planned study at .05.
  lenient <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5,
                            alpha_prior = .10)
  expect_plan(lenient, 136, 1.6995658132, 0.8026224932)
  expect_equal(val(lenient, "actual_power"),
               power_t2(ws_df(136, 136, 1.5), scale(136) * val(lenient, "ncp_adjusted")))
  expect_lt(power_t2(ws_df(135, 135, 1.5), scale(135) * val(lenient, "ncp_adjusted")), .8)
})

test_that("VIG-5: the Welch page's '.76 rather than .80' and '.43' are what the plan delivers", {
  # The help page and the worked-examples vignette say that the plan built
  # assuming equal variances (213 per group) delivers about .76 power if the
  # ratio is in fact 2 in both studies, and about .43 if the prior study's
  # groups were equally variable but the planned study's ratio is 2, with the
  # mean difference fixed in units of the first group's standard deviation.
  # Neither figure was computed or pinned anywhere. Recomputed 2026-10-12
  # without BUCSS: .759520 and .425674 by the noncentral t on the population
  # Welch-Satterthwaite df, .759517 and .425673 for the real Welch test (its
  # df estimated from the two sample variances, by quadrature over them), and
  # .7606 (SE .0010) and .4271 (SE .0011) in 200,000 real tests on raw data.
  ws_df <- function(a, b, r) {
    v1 <- 1 / a
    v2 <- r^2 / b
    (v1 + v2)^2 / (v1^2 / (a - 1) + v2^2 / (b - 1))
  }
  # Power of the real Welch test with n per group and ratio r, at the
  # noncentrality 'ncp' of the mean difference over its true standard error:
  # a midpoint rule over the two sample variances, each a scaled chi square.
  real_welch_power <- function(n, r, ncp, K = 300L) {
    u <- (seq_len(K) - .5) / K
    g <- expand.grid(a = qchisq(u, n - 1) / (n - 1), b = r^2 * qchisq(u, n - 1) / (n - 1))
    df_hat <- (g$a + g$b)^2 * (n - 1) / (g$a^2 + g$b^2)
    k <- qt(.975, df_hat) * sqrt((g$a + g$b) / (1 + r^2))
    mean(pnorm(-k + ncp) + pnorm(-k - ncp))
  }
  equal <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1)
  unequal <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 2)
  n <- size(equal)
  expect_identical(n, 213)

  # The ratio is 2 in both studies: the parameter is the ratio-2 plan's,
  # rescaled to 213 per group at a ratio of 2.
  ncp_both <- val(unequal, "ncp_adjusted") * sqrt((1 / 40 + 4 / 55) / ((1 + 4) / n))
  expect_equal(power_t2(ws_df(n, n, 2), ncp_both), 0.7595196888, tolerance = 1e-8)
  expect_lt(abs(real_welch_power(n, 2, ncp_both) - .759517), 2e-5)
  expect_equal(round(real_welch_power(n, 2, ncp_both), 2), .76)

  # Equal variances in the prior study, a ratio of 2 in the planned one: the
  # mean difference in units of the first group's standard deviation is the
  # ratio-1 plan's parameter times the prior standard error.
  difference <- val(equal, "ncp_adjusted") * sqrt(1 / 40 + 1 / 55)
  ncp_planned <- difference / sqrt((1 + 4) / n)
  expect_equal(power_t2(ws_df(n, n, 2), ncp_planned), 0.4256736461, tolerance = 1e-8)
  expect_lt(abs(real_welch_power(n, 2, ncp_planned) - .425673), 2e-5)
  expect_equal(round(real_welch_power(n, 2, ncp_planned), 2), .43)
})

test_that("TQ-5: every planner page states the scale of ncp_adjusted", {
  # The returned value's scale differs by planner (t, F, or chi square), and
  # only the help page says which. With the statement on the regression
  # coefficient page changed to "the t statistic" (in the roxygen source and
  # the Rd together, so the release gate's documentation check saw nothing)
  # the whole suite passed. The statement is pinned here, with the fact it
  # states for that page: the coefficient planner's parameter is on the F
  # scale, equal to the R2 planner's for one predictor.
  man <- file.path("..", "..", "man")
  db <- if (dir.exists(man)) tools::Rd_db(dir = file.path("..", "..")) else tools::Rd_db("BUCSS")
  skip_if(length(db) == 0L, "no Rd sources available")
  names(db) <- basename(names(db))
  value_text <- function(rd) {
    v <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\value"), rd)
    gsub("\\s+", " ", paste(unlist(v), collapse = ""))
  }
  t_pages <- c("ss_buc_independent_t.Rd", "ss_buc_paired_t.Rd", "ss_buc_one_sample_t.Rd",
               "ss_buc_welch_t.Rd")
  pages <- setdiff(grep("^ss_buc_.*[.]Rd$", names(db), value = TRUE), "ss_buc_sensitivity.Rd")
  expect_length(pages, 18L)
  for (nm in pages) {
    scale <- if (nm %in% t_pages) "on the scale of the t statistic" else
      if (nm == "ss_buc_chisq_diff.Rd") "on the scale of the chi square statistic" else
        "on the scale of the F statistic"
    expect_match(value_text(db[[nm]]), paste0("noncentrality parameter, ", scale),
                 fixed = TRUE, info = nm)
  }
  expect_match(value_text(db[["ss_buc_reg_coef.Rd"]]), "(the square of the prior t)",
               fixed = TRUE)
  expect_equal(val(ss_buc_reg_coef(t_observed = 3, N = 100, p = 1), "ncp_adjusted"),
               val(ss_buc_R2(F_observed = 9, N = 100, p = 1), "ncp_adjusted"),
               tolerance = 1e-10)
})
