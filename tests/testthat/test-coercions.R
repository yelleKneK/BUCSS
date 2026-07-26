# The documented input coercions are load-bearing and must be preserved: a
# percentage (>= 1) for 'assurance' or 'power' is divided by 100, and
# 'alpha_prior == 1' models "no publication bias". These pin that behavior so a
# refactor cannot quietly drop it.

test_that("assurance and power entered as percentages equal the proportion form", {
  pct  <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = 80, power = 80)
  prop <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, power = .80)

  expect_equal(pct$value, prop$value)
  # the echoed planning inputs are the coerced proportions
  expect_equal(attr(pct, "inputs")$assurance, .80)
  expect_equal(attr(pct, "inputs")$power, .80)
})

test_that("alpha_prior = 1 models no publication bias and is echoed as entered", {
  no_bias <- ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = 1,
                                  assurance = .80, power = .80)
  default <- ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = .05,
                                  assurance = .80, power = .80)

  expect_s3_class(no_bias, "bucss_power")
  # the user's value (1) is echoed, not the internal .999 used to avoid a
  # degenerate truncation
  expect_equal(attr(no_bias, "inputs")$alpha_prior, 1)
  # assuming no bias needs fewer participants than assuming the strongest bias
  expect_lt(no_bias$value[no_bias$term == "necessary_sample_size"],
            default$value[default$term == "necessary_sample_size"])
  # pin the .999-substitution path's output (a verified 1.2.1 match), so the
  # no-publication-bias computation cannot drift while the class and echo
  # assertions above keep passing
  expect_equal(no_bias$value[no_bias$term == "necessary_sample_size"], 37)
  expect_equal(round(no_bias$value[no_bias$term == "ncp_adjusted"], 3), 2.090)
})

test_that("effect-argument shorthand is accepted and normalized to the canonical value", {
  base <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                                 effect = "factor_B", assurance = .80, power = .80)
  for (syn in c("B", "b")) {
    alt <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                                  effect = syn, assurance = .80, power = .80)
    expect_equal(alt$value, base$value)
    expect_identical(attr(alt, "effect"), "factor_B")
  }
  expect_identical(
    attr(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                                effect = "AxB", assurance = .80, power = .80),
         "effect"),
    "interaction")
  expect_identical(
    attr(ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
                            levels_within = 3, effect = "ws", assurance = .80,
                            power = .80),
         "effect"),
    "within")
})
