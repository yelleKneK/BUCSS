# The documented input coercions are load-bearing and must be preserved: a
# percentage (>= 1) for 'assurance' or 'power' is divided by 100, and
# 'alpha_prior = 1' models "no publication filter", computed as exactly 1.
# These pin that behavior so a refactor cannot quietly drop it.

test_that("assurance and power entered as percentages equal the proportion form", {
  pct  <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = 80, desired_power = 80)
  prop <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, desired_power = .80)

  expect_equal(pct$value, prop$value)
  # the echoed planning inputs are the coerced proportions
  expect_equal(pct$value[pct$term == "assurance"], .80)
  expect_equal(pct$value[pct$term == "desired_power"], .80)
})

test_that("alpha_prior = 1 models no publication filter and is echoed as entered", {
  no_bias <- ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = 1,
                                  assurance = .80, desired_power = .80)
  default <- ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = .05,
                                  assurance = .80, desired_power = .80)

  expect_s3_class(no_bias, "bucss_plan")
  # the value entered is echoed, and it is the value computed with
  expect_identical(no_bias$value[no_bias$term == "alpha_prior"], 1)
  # assuming no filter needs fewer participants than assuming the default one
  expect_lt(no_bias$value[no_bias$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1],
            default$value[default$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1])
  # With no publication filter the adjusted parameter solves the ordinary
  # distribution of |t| at the observed value: P(|T| <= 3 | df 38, ncp) = .80.
  # Independent check 2026-10-13 (base R, exactly 1): 2.0903759148 and 37 per
  # group. Until then 1 was computed as .999, which gave 2.0902910290 (and 37),
  # so the plan cannot drift back while the class and echo assertions pass.
  by_hand <- uniroot(function(ncp) suppressWarnings(pt(3, 38, ncp) - pt(-3, 38, ncp)) - .80,
                     c(0, 10), tol = 1e-12)$root
  stopifnot(abs(by_hand - 2.0903759148) < 1e-8)
  expect_equal(no_bias$value[no_bias$term == "ncp_adjusted"], by_hand, tolerance = 1e-7)
  expect_gt(no_bias$value[no_bias$term == "ncp_adjusted"], 2.0903)
  expect_equal(no_bias$value[no_bias$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1], 37)
  # and the largest supportable assurance is exactly 1 - p
  expect_equal(attr(no_bias, "assurance_ceiling"), 1 - 2 * pt(3, 38, lower.tail = FALSE),
               tolerance = 1e-14)
})

test_that("effect-argument shorthand is accepted and normalized to the canonical value", {
  base <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                                 effect = "factor_B", assurance = .80, desired_power = .80)
  for (syn in c("B", "b")) {
    alt <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                                  effect = syn, assurance = .80, desired_power = .80)
    expect_equal(alt$value, base$value)
    expect_identical(attr(alt, "effect"), "factor_B")
  }
  expect_identical(
    attr(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                                effect = "AxB", assurance = .80, desired_power = .80),
         "effect"),
    "interaction")
  expect_identical(
    attr(ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
                            levels_within = 3, effect = "ws", assurance = .80,
                            desired_power = .80),
         "effect"),
    "within")
})
