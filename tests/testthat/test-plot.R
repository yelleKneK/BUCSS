# plot.bucss_power(): the assurance ceiling, drawn.
#
# Every test draws to a null device, so nothing is written and no display is
# needed; what is checked is the returned curve, which is what the figure is
# made of.

with_null_device <- function(code) {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  force(code)
}

test_that("plot returns the plotted curve invisibly, in the documented shape", {
  plan <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80)
  with_null_device({
    vis <- withVisible(plot(plan))
    expect_false(vis$visible)
    curve_values <- vis$value
  })
  expect_s3_class(curve_values, "data.frame")
  expect_identical(names(curve_values),
                   c("assurance", "ncp_adjusted", "sample_size", "actual_power"))
  expect_equal(attr(curve_values, "assurance_ceiling"),
               attr(plan, "assurance_ceiling"))
  expect_identical(attr(curve_values, "design"), attr(plan, "design"))
  expect_identical(attr(curve_values, "sample_size_unit"),
                   attr(plan, "sample_size_unit"))
  expect_equal(nrow(curve_values), 60)
})

test_that("the curve has the shape the figure claims", {
  plan <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80)
  with_null_device(cv <- plot(plan, n_points = 40))
  ok <- !is.na(cv$ncp_adjusted)
  expect_true(all(ok))
  # every evaluated assurance is below the ceiling, and none is a percentage
  expect_true(all(cv$assurance < attr(plan, "assurance_ceiling")))
  expect_true(all(cv$assurance > 0 & cv$assurance < 1))
  # the correction falls, the sample size climbs, and the target is always met
  expect_true(all(diff(cv$ncp_adjusted) < 0))
  expect_true(all(diff(cv$sample_size) >= 0))
  expect_true(all(cv$actual_power >= .80))
  # the grid packs points against the wall: over a third of them fall in the
  # last tenth of the assurance range, which is what makes the wall visible
  span <- range(cv$assurance)
  in_last_tenth <- mean(cv$assurance > span[2] - .1 * diff(span))
  expect_gt(in_last_tenth, 1 / 3)
  # the plan's own values lie on its own curve
  with_null_device(at_plan <- plot(plan, assurance = c(.5, .80)))
  expect_equal(at_plan$ncp_adjusted[at_plan$assurance == .80],
               plan$value[plan$term == "ncp_adjusted"])
  expect_equal(at_plan$sample_size[at_plan$assurance == .80],
               plan$value[plan$term == "necessary_n_per_group"])
})

test_that("plot works for every design in the registry", {
  plans <- list(
    ss_buc_independent_t(t_observed = 3, n = c(50, 50), assurance = .80),
    ss_buc_paired_t(t_observed = 3, N = 30),
    ss_buc_one_sample_t(t_observed = 3, N = 30),
    ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 40, sd_ratio = 1.5),
    ss_buc_one_way_anova(F_observed = 8, N = 60, levels_A = 3),
    ss_buc_factorial_anova(F_observed = 12, N = 60, levels_A = 3, levels_B = 2,
                           effect = "factor_A"),
    ss_buc_factorial_anova_general(F_observed = 12, N = 60, cells = 6,
                                   df_numerator = 2, df_denominator = 54),
    ss_buc_ancova(F_observed = 12, N = 60, cells = 6, df_numerator = 2,
                  n_covariates = 2),
    ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3),
    ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2),
    ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                       levels_within = 3, effect = "between"),
    ss_buc_mixed_anova_general(F_observed = 20, N = 60, num_groups = 2,
                               df_numerator = 1, effect = "between_only"),
    ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3),
    ss_buc_R2(F_observed = 8, N = 60, p = 3),
    ss_buc_reg_joint(F_observed = 12, N = 60, p = 5, p_joint = 2),
    ss_buc_manova(F_observed = 8, N = 40, p_variables = 3),
    ss_buc_chisq_diff(chisq_observed = 15, df_full = 42, df_restricted = 44,
                      N = 200),
    ss_buc_correlation(r_observed = .35, N = 100)
  )
  expect_setequal(unique(vapply(plans, attr, character(1), "design")),
                  names(BUCSS:::.BUCSS_DESIGN_SPECS))
  for (plan in plans) {
    with_null_device(cv <- plot(plan, n_points = 12))
    info <- attr(plan, "design")
    expect_equal(nrow(cv), 12, info = info)
    expect_false(anyNA(cv$sample_size), info = info)
    expect_true(all(diff(cv$ncp_adjusted) < 0), info = info)
    expect_true(all(diff(cv$sample_size) >= 0), info = info)
  }
})

test_that("plot accepts its arguments and restores the graphics state", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  with_null_device({
    before <- par("mfrow")
    invisible(plot(plan, which = "both", n_points = 8))
    expect_identical(par("mfrow"), before)
    invisible(plot(plan, which = "ncp", n_points = 8))
    expect_identical(par("mfrow"), before)
    invisible(plot(plan, which = "size", n_points = 8))
    expect_identical(par("mfrow"), before)
    # arguments in ... override the defaults rather than colliding with them
    expect_silent(invisible(plot(plan, n_points = 8, main = "custom",
                                 ylab = "custom", col = "navy")))
    # an explicit grid is used as given, sorted
    cv <- plot(plan, assurance = c(.8, .5, .7))
    expect_equal(cv$assurance, c(.5, .7, .8))
    # a lower bound below the ceiling but above the default is honored
    cv <- plot(plan, lower = .75, n_points = 10)
    expect_equal(min(cv$assurance), .75)
  })
})

test_that("plot validates its inputs", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  ceiling_value <- attr(plan, "assurance_ceiling")
  with_null_device({
    expect_error(plot(plan, which = "elsewhere"), "should be one of")
    expect_error(plot(plan, lower = ceiling_value + .01),
                 "below the plan's assurance ceiling", fixed = TRUE)
    expect_error(plot(plan, lower = 0), "below the plan's assurance ceiling",
                 fixed = TRUE)
    expect_error(plot(plan, n_points = 2.5), "whole number", fixed = TRUE)
    expect_error(plot(plan, assurance = .8),
                 "at least two finite values", fixed = TRUE)
    # an assurance of exactly 1 would be read by the planners as 1 percent
    expect_error(plot(plan, assurance = c(.5, 1)),
                 "strictly between 0 and 1", fixed = TRUE)
    expect_error(plot(plan, assurance = c(.5, NA)),
                 "at least two finite values", fixed = TRUE)
  })
})
