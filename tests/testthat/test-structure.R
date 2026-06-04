# Smoke tests for the tidy bucss_power object that every ss_buc_* function
# returns: its class, columns, attribute payload, and print method.

test_that("ss_buc_* returns a tidy bucss_power data.frame", {
  result <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, power = .80)

  expect_s3_class(result, "bucss_power")
  expect_s3_class(result, "data.frame")
  expect_identical(result$term, c("necessary_sample_size", "ncp_adjusted"))
  expect_type(result$value, "double")
  expect_identical(nrow(result), 2L)
})

test_that("design, unit, and inputs travel on attributes (not in value)", {
  result <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, power = .80)

  expect_identical(attr(result, "design"), "Independent t test")
  expect_identical(attr(result, "sample_size_unit"), "per group")
  expect_null(attr(result, "effect"))
  expect_type(attr(result, "inputs"), "list")
  expect_identical(attr(result, "inputs")$t_observed, 3)
})

test_that("two-way ANOVA results also carry the effect tested", {
  result <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
                                   levels_B = 3, effect = "factor_B",
                                   assurance = .80, power = .80)

  expect_identical(attr(result, "effect"), "factor_B")
  expect_identical(attr(result, "sample_size_unit"), "per cell")
})

test_that("the two quantities stay reachable as a plain data.frame", {
  result <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, power = .80)

  expect_equal(result$value[result$term == "necessary_sample_size"], 130)
  expect_equal(round(result$value[result$term == "ncp_adjusted"], 3), 1.105)
})

test_that("print.bucss_power shows the human-readable summary", {
  result <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
                                   levels_B = 3, effect = "factor_B",
                                   assurance = .80, power = .80)

  expect_output(print(result),
                "Bias and uncertainty corrected sample size", fixed = TRUE)
  expect_output(print(result),
                "Design: Two-way between-subjects ANOVA", fixed = TRUE)
  expect_output(print(result), "Effect of interest: factor_B", fixed = TRUE)
  expect_output(print(result), "Necessary sample size", fixed = TRUE)
  expect_output(print(result), "Planning inputs", fixed = TRUE)
})

test_that("print returns its argument invisibly", {
  result <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, power = .80)

  out <- capture.output(vis <- withVisible(print(result)))
  expect_false(vis$visible)
  expect_identical(vis$value, result)
})

test_that("results carry the assurance ceiling and, for cell designs, the implied total N", {
  res <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, power = .80)

  # The ceiling is the closed form 1 - p/alpha_prior for the prior result.
  p_prior <- 2 * (1 - pt(3, df = 38))
  expect_equal(attr(res, "assurance_ceiling"), 1 - p_prior / .05)
  expect_equal(attr(res, "total_n"),
               2 * res$value[res$term == "necessary_sample_size"])

  # "total"-unit designs report the ceiling but no separate total_n.
  rj <- ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                         assurance = .80, power = .80)
  expect_null(attr(rj, "total_n"))
  expect_false(is.null(attr(rj, "assurance_ceiling")))
})

test_that("print shows the implied total N and the assurance ceiling", {
  res <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                                effect = "factor_B", assurance = .80, power = .80)

  expect_output(print(res), "total N = ", fixed = TRUE)
  expect_output(print(res), "Assurance this prior can support", fixed = TRUE)
})
