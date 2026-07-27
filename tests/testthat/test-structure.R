# Smoke tests for the tidy bucss_power object that every ss_buc_* function
# returns: its class, rows (design results followed by echoed planning
# inputs, in DMAR's style), attribute payload, and print method.

test_that("ss_buc_* returns a tidy bucss_power data.frame", {
  result <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, desired_power = .80)

  expect_s3_class(result, "bucss_power")
  expect_s3_class(result, "data.frame")
  expect_type(result$value, "double")
  # design results first (the size row is named for its unit), then the rows
  # echoing the planning inputs
  expect_identical(result$term,
                   c("necessary_n_per_group", "total_N", "actual_power",
                     "ncp_adjusted", "t_observed", "n", "alpha_prior",
                     "alpha_planned", "assurance", "desired_power"))

  total_unit <- ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                                 assurance = .80, desired_power = .80)
  expect_identical(total_unit$term,
                   c("necessary_N", "actual_power", "ncp_adjusted",
                     "F_observed", "N", "p", "p_joint", "alpha_prior",
                     "alpha_planned", "assurance", "desired_power"))

  # a length-2 prior n echoes as n_1/n_2 rows
  unequal <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .90)
  expect_true(all(c("n_1", "n_2") %in% unequal$term))
  expect_equal(unequal$value[unequal$term == "n_1"], 50)
  expect_equal(unequal$value[unequal$term == "n_2"], 55)
})

test_that("design, unit, and effect travel on attributes (not in value)", {
  result <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, desired_power = .80)

  expect_identical(attr(result, "design"), "Independent t test")
  expect_identical(attr(result, "sample_size_unit"), "per group")
  expect_null(attr(result, "effect"))
  # the planning inputs are rows, not an attribute
  expect_null(attr(result, "inputs"))
  expect_identical(result$value[result$term == "t_observed"], 3)
})

test_that("two-way ANOVA results also carry the effect tested", {
  result <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
                                   levels_B = 3, effect = "factor_B",
                                   assurance = .80, desired_power = .80)

  expect_identical(attr(result, "effect"), "factor_B")
  expect_identical(attr(result, "sample_size_unit"), "per cell")
  expect_identical(result$term[1], "necessary_n_per_cell")
})

test_that("the headline quantities stay reachable as a plain data.frame", {
  result <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, desired_power = .80)

  expect_equal(result$value[result$term == "necessary_n_per_group"], 130)
  expect_equal(round(result$value[result$term == "ncp_adjusted"], 3), 1.105)
})

test_that("actual_power meets the desired power and is conservative", {
  # single-branch design: the achieved power at the returned N is exact
  res <- ss_buc_paired_t(t_observed = 3, N = 40, assurance = .80, desired_power = .80)
  expect_gte(res$value[res$term == "actual_power"], .80)
  expect_lt(res$value[res$term == "actual_power"], .85)

  # two-branch design: evaluated under the conservative branch, still >= target
  res2 <- ss_buc_one_way_anova(F_observed = 5, N = 121, levels_A = 3,
                               assurance = .80, desired_power = .80)
  expect_gte(res2$value[res2$term == "actual_power"], .80)
})

test_that("print.bucss_power shows the aligned table and factual footers", {
  result <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
                                   levels_B = 3, effect = "factor_B",
                                   assurance = .80, desired_power = .80)

  expect_output(print(result), "term", fixed = TRUE)
  expect_output(print(result), "necessary_n_per_cell", fixed = TRUE)
  expect_output(print(result),
                "Design: Two-way between-subjects ANOVA (factor_B)", fixed = TRUE)
  expect_output(print(result), "Sample size unit: per cell", fixed = TRUE)
  expect_output(print(result), "Largest supportable assurance:", fixed = TRUE)
})

test_that("print returns its argument invisibly", {
  result <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, desired_power = .80)

  out <- capture.output(vis <- withVisible(print(result)))
  expect_false(vis$visible)
  expect_identical(vis$value, result)
})

test_that("results carry the assurance ceiling and, for cell designs, a total_N row", {
  res <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, desired_power = .80)

  # The ceiling is the closed form 1 - p/alpha_prior for the prior result.
  p_prior <- 2 * (1 - pt(3, df = 38))
  expect_equal(attr(res, "assurance_ceiling"), 1 - p_prior / .05)
  expect_equal(res$value[res$term == "total_N"],
               2 * res$value[res$term == "necessary_n_per_group"])

  # "total"-unit designs report the ceiling but carry no total_N row.
  rj <- ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                         assurance = .80, desired_power = .80)
  expect_false("total_N" %in% rj$term)
  expect_false(is.null(attr(rj, "assurance_ceiling")))
})

test_that("tidy() gives the compact estimate view and glance() the wide view", {
  res <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                                effect = "factor_B", assurance = .80, desired_power = .80)

  # tidy(): term/estimate/power, as in DMAR's sample size planners
  td <- tidy(res)
  expect_s3_class(td, "data.frame")
  expect_identical(nrow(td), 1L)
  expect_identical(td$term, "sample_size")
  expect_equal(td$estimate, res$value[res$term == "necessary_n_per_cell"])
  expect_equal(td$power, res$value[res$term == "actual_power"])

  # glance(): the estimate plus every echoed input and the design metadata
  gl <- glance(res)
  expect_identical(nrow(gl), 1L)
  expect_identical(gl$design, "Two-way between-subjects ANOVA")
  expect_identical(gl$effect, "factor_B")
  expect_equal(gl$necessary_n_per_cell, res$value[res$term == "necessary_n_per_cell"])
  expect_equal(gl$total_N, res$value[res$term == "total_N"])
  expect_equal(gl$F_observed, 5)
  expect_equal(gl$desired_power, .80)
  expect_equal(gl$df_effect, 2)     # levels_B - 1
  expect_equal(gl$df_error, 3948)   # 659 per cell * 6 cells - 6
  expect_equal(gl$assurance_ceiling, attr(res, "assurance_ceiling"))

  # a total-unit design: no total_N column anywhere, df still present
  rj <- ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                         assurance = .80, desired_power = .80)
  expect_false("total_N" %in% names(glance(rj)))
  expect_equal(glance(rj)$df_effect, 2)   # p_joint
})

test_that("planning_sentence() writes the manuscript sentence", {
  res <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .90)
  s <- planning_sentence(res)
  expect_match(s, "1485 participants per group", fixed = TRUE)
  expect_match(s, "total N = 2970", fixed = TRUE)
  expect_match(s, "80% power", fixed = TRUE)
  expect_match(s, "90% assurance", fixed = TRUE)
  expect_error(planning_sentence(1), "bucss_power", fixed = TRUE)
})

test_that("row and column subsetting is not hijacked by the legacy [[ method", {
  # [.data.frame extracts columns with [[, which the legacy positional [[
  # contract would otherwise intercept for columns 1 and 2 (term and value),
  # silently replacing them with the sample size and the NCP.
  res <- ss_buc_paired_t(t_observed = 3, N = 40)
  expect_identical(head(res, 3)$term, res$term[1:3])
  expect_identical(head(res, 3)$value, res$value[1:3])
  expect_equal(res[res$term == "ncp_adjusted", ]$value,
               res$value[res$term == "ncp_adjusted"])
  expect_identical(res[, "term"], res$term)
  expect_identical(nrow(res[res$term %in% c("necessary_N", "ncp_adjusted"), ]), 2L)
  # a subset is a plain data.frame, and the legacy contract still holds on the
  # unsubsetted object
  expect_false(inherits(head(res, 2), "bucss_power"))
  expect_identical(res[[1]], res$value[res$term == "necessary_N"])
  expect_identical(res[[2]], res$value[res$term == "ncp_adjusted"])
})
