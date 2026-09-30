# The worked examples of Anderson, Kelley, and Maxwell (2017, Psychological
# Science, 28, 1547-1562) and of its supplemental material ("Conceptual Logic
# of Taylor and Muller's Method"), checked against the article's printed
# values on 2026-09-30.
#
# The article's three researchers plan from effect sizes it describes as
# simulation medians (its notes 18 to 20), printed rounded: d = 0.68,
# f^2 = 0.09, d_z = 0.50. A sample size far from the assurance ceiling does not
# notice the rounding; one near it does (TM 5th percentile for the second
# researcher moves from 1,418 to 1,400 between f^2 = .08996 and .09). Each
# printed size is pinned at an unrounded value that rounds to the printed
# effect size, and the rounded input is pinned beside it, so a reader can see
# that the difference is the input's rounding and not the method. BUCSS 1.2.1
# gives the same sizes at the rounded inputs (330, 64, 50, and 80; 1,401 where
# 2.0.0 gives 1,400 is the 1.x grid).

plan_size <- function(x) x$value[x$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1]
plan_ncp  <- function(x) x$value[x$term == "ncp_adjusted"]

test_that("supplement Example 1: d = 0.68 with 64 per group (Table S1)", {
  t_obs <- 0.68 * sqrt(64 / 2)
  plans <- lapply(c(.50, .80, .95), function(a)
    ss_buc_independent_t(t_observed = t_obs, n = 64, assurance = a))
  # "50% assurance requires only 36 participants per group", 80% gives 65,
  # "95% assurance requires 169 participants per group".
  expect_identical(vapply(plans, plan_size, numeric(1)), c(36, 65, 169))
  # Table S1's delta at each assurance (the t-scale parameter times sqrt(2/n)).
  delta <- vapply(plans, plan_ncp, numeric(1)) * sqrt(2 / 64)
  expect_identical(round(delta, 2), c(0.67, 0.50, 0.31))
  # Table S1's lambda (the F scale, the square of the t-scale parameter).
  lambda <- vapply(plans, plan_ncp, numeric(1))^2
  expect_identical(round(lambda[c(1, 3)], 1), c(14.4, 3.0))
})

test_that("supplement Example 2: d = 0.68 with 25 per group (Table S2)", {
  t_obs <- 0.68 * sqrt(25 / 2)
  plan <- ss_buc_independent_t(t_observed = t_obs, n = 25, assurance = .5)
  # "a delta of 0.22 corresponds to 50% assurance"; lambda .60.
  expect_identical(round(plan_ncp(plan) * sqrt(2 / 25), 2), 0.22)
  expect_identical(round(plan_ncp(plan)^2, 2), 0.60)
  # "the smallest value of the conditional probability ... is .4024", which
  # prevents any assurance above about 60%.
  expect_identical(round(1 - attr(plan, "assurance_ceiling"), 4), 0.4024)
  expect_error(ss_buc_independent_t(t_observed = t_obs, n = 25, assurance = .8),
               class = "bucss_zero_ncp")
})

test_that("the article's first researcher (underpowered prior study)", {
  # The TM 20th percentile "returns a lambda of zero in this case".
  expect_error(ss_buc_independent_t(t_observed = 0.68 * sqrt(12.5), n = 25,
                                    assurance = .8),
               class = "bucss_zero_ncp")
  # The 50th percentile "suggests 332 participants per group" with power .73
  # at the true delta of 0.2.
  plan <- ss_buc_independent_t(t_observed = 0.6799 * sqrt(12.5), n = 25,
                               assurance = .5)
  expect_identical(plan_size(plan), 332)
  expect_identical(plan_size(ss_buc_independent_t(
    t_observed = 0.68 * sqrt(12.5), n = 25, assurance = .5)), 330)
  df <- 2 * 332 - 2
  crit <- qt(.975, df)
  true_ncp <- 0.2 * sqrt(332 / 2)
  power <- pt(crit, df, ncp = true_ncp, lower.tail = FALSE) +
    pt(-crit, df, ncp = true_ncp)
  expect_identical(round(power, 2), 0.73)
})

test_that("the article's second researcher (3 x 2 between-subjects main effect)", {
  # Sample f^2 = SS_A / SS_error with 25 per cell (N = 150), so
  # F = f^2 * 144 / 2. The article prints 24, 64, and 1,418 per group.
  # F = 6.4773 is f^2 = .08996; the rounded f^2 = .09 is F = 6.48.
  plans_at <- function(F) vapply(c(.50, .80, .95), function(a)
    plan_size(ss_buc_factorial_anova(F_observed = F, N = 150,
                                     levels_A = 3, levels_B = 2,
                                     effect = "factor_A", assurance = a)),
    numeric(1))
  expect_identical(plans_at(6.4773), c(24, 64, 1418))
  expect_identical(plans_at(0.09 * 144 / 2), c(24, 64, 1400))
  # The 50th percentile plan "will result in a power of .76" at f^2 = .0625.
  N <- 24 * 6
  power <- pf(qf(.95, 2, N - 6), 2, N - 6, ncp = 0.0625 * N, lower.tail = FALSE)
  expect_identical(round(power, 2), 0.76)
})

test_that("the article's third researcher (paired t, N = 100)", {
  # The article prints 34, 49, and 79 with power .81, .93, and .99 at d_z = .5.
  plans_at <- function(dz) vapply(c(.50, .80, .95), function(a)
    plan_size(ss_buc_paired_t(t_observed = dz * sqrt(100), N = 100,
                              assurance = a)),
    numeric(1))
  expect_identical(plans_at(0.502), c(34, 49, 79))
  expect_identical(plans_at(0.50), c(34, 50, 80))
  power <- vapply(c(34, 49, 79), function(N) {
    crit <- qt(.975, N - 1)
    pt(crit, N - 1, ncp = 0.5 * sqrt(N), lower.tail = FALSE) +
      pt(-crit, N - 1, ncp = 0.5 * sqrt(N))
  }, numeric(1))
  expect_identical(round(power, 2), c(0.81, 0.93, 0.99))
})
