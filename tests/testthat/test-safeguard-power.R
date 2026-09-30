# Safeguard power (Perugini, Gallucci, & Costantini, 2014, Perspectives on
# Psychological Science, 9, 319-332) plans from the lower limit of a two-sided
# 60% confidence interval for the prior effect. It is the special case
# alpha_prior = 1, assurance = .80: with no publication threshold the
# truncation disappears, TM is the ordinary distribution function of the
# statistic, and its .80 root is the lower limit of a one-sided 80% (two-sided
# 60%) interval for the noncentrality parameter; a two-sided level C is
# assurance (1 + C)/2. The limits below are computed from scratch, by
# inverting pt(), pf(), and pchisq() with uniroot(), with no call into the
# package and no MBESS; the plans are the ordinary power analysis at the limit,
# found by a linear walk.
#
# Two departures are deliberate and one-directional, each making the package's
# parameter smaller (so its plan is never smaller): alpha_prior = 1 is read as
# .999, and the effects on one degree of freedom bound the magnitude of the
# effect from |t|, which adds P(T < -t) to the signed interval's probability.
# Every comparison therefore asserts the direction and a relative difference of
# at most 1e-3. The largest relative difference among the priors pinned here is
# 2.3e-4 (the regression coefficient, t = 2.5), then 1.4e-4 (confidence level
# .80); every F and chi square prior on more than one degree of freedom is
# within 9.0e-5 (measured 2026-09-30, and against MBESS 5.0.1's
# conf.limits.nct, conf.limits.ncf, and conf.limits.nc.chisq, which agree with
# these inversions to 1.3e-8). The weak-prior region, where the |t|
# construction declines a prior the signed interval would still plan from, is
# pinned in the last block.

sg_size <- function(x) x$value[x$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1]
sg_ncp <- function(x) x$value[x$term == "ncp_adjusted"]
sg_rel <- function(a, b) abs(a - b) / b

# The signed lower limit of a one-sided interval at 'level'.
lower_t <- function(t, df, level = .80)
  uniroot(function(x) pt(t, df, x) - level, c(0, t + 10), tol = 1e-12)$root
lower_f <- function(F, df1, df2, level = .80)
  uniroot(function(x) pf(F, df1, df2, x) - level, c(0, 10 * F * df1 + 100),
          tol = 1e-12)$root
lower_chisq <- function(x2, df, level = .80)
  uniroot(function(x) pchisq(x2, df, x) - level, c(0, 10 * x2 + 100),
          tol = 1e-12)$root

# The ordinary power plan: the smallest size whose power reaches .80.
first_n <- function(power_at, start) {
  n <- start
  while (power_at(n) < .80) n <- n + 1
  n
}
power_t <- function(ncp, df) {
  cr <- qt(.975, df)
  (1 - pt(cr, df, ncp)) + pt(-cr, df, ncp)
}
power_f <- function(lambda, df1, df2)
  1 - pf(qf(.95, df1, df2), df1, df2, lambda)

test_that("alpha_prior = 1 and assurance = .80 is safeguard power: t tests", {
  # Paired t = 3 with N = 30 pairs; the parameter scales as sqrt(n / N).
  plan <- ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = 1,
                          assurance = .80)
  lim <- lower_t(3, 29)
  expect_lte(sg_rel(sg_ncp(plan), lim), 1e-3)
  expect_lte(sg_ncp(plan), lim)
  expect_identical(sg_size(plan),
                   first_n(function(n) power_t(lim * sqrt(n / 30), n - 1), 2))
  expect_identical(sg_size(plan), 57)

  # Independent t = 3 with 20 per group; the parameter scales as sqrt(n / 20).
  plan <- ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = 1,
                               assurance = .80)
  lim <- lower_t(3, 38)
  expect_lte(sg_rel(sg_ncp(plan), lim), 1e-3)
  expect_lte(sg_ncp(plan), lim)
  expect_identical(sg_size(plan),
                   first_n(function(n) power_t(lim * sqrt(n / 20), 2 * n - 2), 2))
  expect_identical(sg_size(plan), 37)
})

test_that("an unequal-n prior is safeguard power at each rounded design", {
  # Independent t = 3 with 50 and 55: the planner reads the prior at the
  # harmonic mean, 52.38, rounded down and up (52 and 53 per group, df 102 and
  # 104), reports the smaller parameter, and returns the larger plan.
  plan <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), alpha_prior = 1,
                               assurance = .80)
  lim_rd <- lower_t(3, 102)
  lim_ru <- lower_t(3, 104)
  expect_lte(lim_rd, lim_ru)
  expect_lte(sg_rel(sg_ncp(plan), lim_rd), 1e-3)
  expect_lte(sg_ncp(plan), lim_rd)
  plan_rd <- first_n(function(n) power_t(lim_rd * sqrt(n / 52), 2 * n - 2), 2)
  plan_ru <- first_n(function(n) power_t(lim_ru * sqrt(n / 53), 2 * n - 2), 2)
  expect_identical(c(plan_rd, plan_ru), c(91, 93))
  expect_identical(sg_size(plan), max(plan_rd, plan_ru))
  # Safeguard power at the actual group sizes (df 103, d = delta times
  # sqrt(1/50 + 1/55)) plans 92: the rounding costs one participant per group,
  # in the conservative direction.
  lim_actual <- lower_t(3, 103)
  d <- lim_actual * sqrt(1 / 50 + 1 / 55)
  plan_actual <- first_n(function(n) power_t(d * sqrt(n / 2), 2 * n - 2), 2)
  expect_identical(plan_actual, 92)
  expect_identical(sg_size(plan) - plan_actual, 1)
})

test_that("a very unequal prior is read on the harmonic mean design's df", {
  # Independent t = 3 with 4 and 100: the harmonic mean, 7.69, is read as 7 and
  # 8 per group on df 12 and 14, not on the prior study's own df 102. The fewer
  # degrees of freedom lower the limit, so the plan is 18 per group where
  # safeguard power at the actual group sizes plans 15, well beyond the rounding
  # term (8 / 7.69 - 1, about 4 percent). The package page and the assurance
  # vignette quote these numbers.
  plan <- ss_buc_independent_t(t_observed = 3, n = c(4, 100), alpha_prior = 1,
                               assurance = .80)
  lim_rd <- lower_t(3, 12)
  lim_ru <- lower_t(3, 14)
  expect_lte(lim_rd, lim_ru)
  expect_lte(sg_rel(sg_ncp(plan), lim_rd), 1e-3)
  expect_lte(sg_ncp(plan), lim_rd)
  plan_rd <- first_n(function(n) power_t(lim_rd * sqrt(n / 7), 2 * n - 2), 2)
  plan_ru <- first_n(function(n) power_t(lim_ru * sqrt(n / 8), 2 * n - 2), 2)
  expect_identical(c(plan_rd, plan_ru), c(16, 18))
  expect_identical(sg_size(plan), max(plan_rd, plan_ru))
  expect_identical(sg_size(plan), 18)
  lim_actual <- lower_t(3, 102)
  d <- lim_actual * sqrt(1 / 4 + 1 / 100)
  plan_actual <- first_n(function(n) power_t(d * sqrt(n / 2), 2 * n - 2), 2)
  expect_identical(plan_actual, 15)
})

test_that("alpha_prior = 1 and assurance = .80 is safeguard power: F and chi square", {
  # One-way ANOVA, F = 5 with N = 60 in three groups; lambda scales with N.
  plan <- ss_buc_one_way_anova(F_observed = 5, N = 60, levels_A = 3,
                               alpha_prior = 1, assurance = .80)
  lim <- lower_f(5, 2, 57)
  expect_lte(sg_rel(sg_ncp(plan), lim), 1e-3)
  expect_lte(sg_ncp(plan), lim)
  expect_identical(sg_size(plan),
                   first_n(function(n) power_f(lim * n / 20, 2, 3 * n - 3), 2))
  expect_identical(sg_size(plan), 47)

  # The same F with N = 61, an unbalanced prior read as 20 and as 21 per group
  # (df 57 and 60); the reported parameter is the smaller, the plan the larger.
  plan <- ss_buc_one_way_anova(F_observed = 5, N = 61, levels_A = 3,
                               alpha_prior = 1, assurance = .80)
  lim_rd <- lower_f(5, 2, 57)
  lim_ru <- lower_f(5, 2, 60)
  expect_lte(sg_rel(sg_ncp(plan), min(lim_rd, lim_ru)), 1e-3)
  expect_lte(sg_ncp(plan), min(lim_rd, lim_ru))
  expect_identical(sg_size(plan), max(
    first_n(function(n) power_f(lim_rd * n / 20, 2, 3 * n - 3), 2),
    first_n(function(n) power_f(lim_ru * n / 21, 2, 3 * n - 3), 2)))
  expect_identical(sg_size(plan), 49)

  # Omnibus R squared, F = 5 with N = 150 and four fixed predictors.
  plan <- ss_buc_R2(F_observed = 5, N = 150, p = 4, alpha_prior = 1,
                    assurance = .80, random_predictors = FALSE)
  lim <- lower_f(5, 4, 145)
  expect_lte(sg_rel(sg_ncp(plan), lim), 1e-3)
  expect_lte(sg_ncp(plan), lim)
  expect_identical(sg_size(plan),
                   first_n(function(n) power_f(lim * n / 150, 4, n - 5), 7))
  expect_identical(sg_size(plan), 181)

  # A single regression coefficient, t = 2.5 with N = 100 and three fixed
  # predictors, is planned on the F = t^2 scale with one numerator df: its
  # interval is the F interval, which is the |t| construction.
  plan <- ss_buc_reg_coef(t_observed = 2.5, N = 100, p = 3, alpha_prior = 1,
                          assurance = .80, random_predictors = FALSE)
  lim <- lower_f(2.5^2, 1, 96)
  expect_lte(sg_rel(sg_ncp(plan), lim), 1e-3)
  expect_lte(sg_ncp(plan), lim)
  expect_identical(sg_size(plan),
                   first_n(function(n) power_f(lim * n / 100, 1, n - 4), 5))
  expect_identical(sg_size(plan), 295)

  # Nested model chi square difference of 12 on 2 df with N = 300.
  plan <- ss_buc_chisq_diff(chisq_observed = 12, N = 300, df_full = 10,
                            df_restricted = 12, alpha_prior = 1,
                            assurance = .80)
  lim <- lower_chisq(12, 2)
  expect_lte(sg_rel(sg_ncp(plan), lim), 1e-3)
  expect_lte(sg_ncp(plan), lim)
  expect_identical(sg_size(plan), first_n(
    function(n) 1 - pchisq(qchisq(.95, 2), 2, lim * n / 300), 2))
  expect_identical(sg_size(plan), 483)
})

test_that("with sampled predictors the equivalence holds in that frame", {
  # The F statistic for R squared with four multivariate normal predictors is
  # a chi square (N - 1 df) mixture of noncentral F's: negative binomial
  # weights with size (N - 1)/2 and probability N/(N + lambda). Its .80 limit
  # and the power analysis in the same frame reproduce the plan.
  mixture_cdf <- function(x, df1, df2, m, prob) {
    k <- 0:qnbinom(1 - 1e-15, size = m / 2, prob = prob)
    sum(dnbinom(k, size = m / 2, prob = prob) *
          pbeta(df1 * x / (df1 * x + df2), df1 / 2 + k, df2 / 2))
  }
  lim <- uniroot(function(l) mixture_cdf(5, 4, 145, 149, 150 / (150 + l)) - .80,
                 c(0, 100), tol = 1e-12)$root
  plan <- ss_buc_R2(F_observed = 5, N = 150, p = 4, alpha_prior = 1,
                    assurance = .80, random_predictors = TRUE)
  expect_lte(sg_rel(sg_ncp(plan), lim), 1e-3)
  expect_lte(sg_ncp(plan), lim)
  prob <- 150 / (150 + lim)
  expect_identical(sg_size(plan), first_n(function(n) {
    1 - mixture_cdf(qf(.95, 4, n - 5), 4, n - 5, n - 1, prob)
  }, 7))
  expect_identical(sg_size(plan), 184)
  # The fixed-frame limit is smaller; the two frames are not interchangeable.
  expect_gt(sg_ncp(plan), lower_f(5, 4, 145))
})

test_that("a two-sided confidence level C corresponds to assurance (1 + C)/2", {
  # C = .50 is assurance .75, and C = .80 is assurance .90.
  for (C in c(.50, .80)) {
    plan <- ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = 1,
                            assurance = (1 + C) / 2)
    lim <- lower_t(3, 29, level = (1 + C) / 2)
    expect_lte(sg_rel(sg_ncp(plan), lim), 1e-3)
    expect_lte(sg_ncp(plan), lim)
  }
  # And assurance .5 is the median unbiased estimate, the zero-width limit.
  plan <- ss_buc_paired_t(t_observed = 3, N = 30, alpha_prior = 1,
                          assurance = .50)
  lim <- lower_t(3, 29, level = .50)
  expect_lte(sg_rel(sg_ncp(plan), lim), 1e-3)
})

test_that("a weak t prior is declined where the signed limit is still positive", {
  # At alpha_prior = 1 the largest supportable assurance is 1 - p/.999, so an
  # assurance of .80 is declined once the two-sided p reaches .1998; the
  # signed one-sided 80% limit stays above zero until the one-sided p reaches
  # .20, a two-sided p of .40.
  t_at <- function(p) qt(1 - p / 2, 29)
  expect_error(ss_buc_paired_t(t_observed = t_at(.21), N = 30, alpha_prior = 1,
                               assurance = .80),
               class = "bucss_zero_ncp")
  expect_gt(lower_t(t_at(.21), 29), 0)
  expect_gt(pt(t_at(.39), 29), .80)
  expect_lt(pt(t_at(.41), 29), .80)
  # Just inside the boundary a plan is issued, and it is larger than the
  # safeguard plan (the |t| construction matters most for a weak prior).
  plan <- ss_buc_paired_t(t_observed = t_at(.19), N = 30, alpha_prior = 1,
                          assurance = .80)
  lim <- lower_t(t_at(.19), 29)
  expect_lt(sg_ncp(plan), lim)
  expect_gt(sg_size(plan),
            first_n(function(n) power_t(lim * sqrt(n / 30), n - 1), 2))
})
