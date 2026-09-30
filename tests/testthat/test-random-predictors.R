# The random_predictors switch of the three regression planners. Sampled
# (multivariate normal) predictors make the statistic a chi square mixture of
# noncentral F distributions, in closed form a negative binomial mixture of
# betas (.p_random_f() in R/bucss-internal.R); ss_buc_correlation() is its
# one-predictor case. The default, FALSE, is the fixed-predictor frame of
# BUCSS 1.x, and its plans are pinned in test-oracles.R. The flag's validation
# errors are in test-validation.R.

row <- function(x, term) x$value[x$term == term]
computed <- function(x) c(x$value[x$term %in% BUCSS:::.BUCSS_SIZE_TERMS],
                          row(x, "actual_power"), row(x, "ncp_adjusted"),
                          attr(x, "assurance_ceiling"))

test_that("one sampled predictor is exactly the correlation planner", {
  # From N = 4, the correlation planner's minimum, as the help pages say; at
  # N = 3 the regression planner still plans and the correlation one refuses.
  for (t in c(2.5, 3.7, 8)) for (N in c(4, 6, 30, 100, 1000)) for (a in c(.5, .8)) {
    reg <- tryCatch(ss_buc_reg_coef(t_observed = t, N = N, p = 1, assurance = a,
                                    random_predictors = TRUE),
                    bucss_zero_ncp = function(e) "refused",
                    bucss_nonsignificant = function(e) "nonsignificant")
    cor <- tryCatch(ss_buc_correlation(t_observed = t, N = N, assurance = a),
                    bucss_zero_ncp = function(e) "refused",
                    bucss_nonsignificant = function(e) "nonsignificant")
    if (is.character(reg)) {
      expect_identical(reg, cor, info = paste(t, N, a))
      next
    }
    expect_identical(computed(reg), computed(cor), info = paste(t, N, a))
    expect_identical(attr(reg, "df_error"), attr(cor, "df_error"))
  }
  # and so is the omnibus test with one sampled predictor
  expect_identical(computed(ss_buc_R2(F_observed = 3.7^2, N = 100, p = 1,
                                      random_predictors = TRUE)),
                   computed(ss_buc_correlation(t_observed = 3.7, N = 100)))
  # while the fixed frame is the regression slope, a few participants lighter
  expect_lt(row(ss_buc_reg_coef(t_observed = 3.7, N = 100, p = 1), "necessary_N"),
            row(ss_buc_correlation(t_observed = 3.7, N = 100), "necessary_N"))
})

test_that("the three regression planners agree where their tests coincide", {
  for (rp in c(FALSE, TRUE)) {
    expect_identical(
      computed(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 4,
                                random_predictors = rp)),
      computed(ss_buc_R2(F_observed = 5, N = 150, p = 4, random_predictors = rp)))
    expect_identical(
      computed(ss_buc_reg_joint(F_observed = 9, N = 150, p = 4, p_joint = 1,
                                random_predictors = rp)),
      computed(ss_buc_reg_coef(t_observed = 3, N = 150, p = 4,
                               random_predictors = rp)))
  }
})

test_that("the mixture matches an independent quadrature and reduces to the central F", {
  # The chi square mixture written out over equal-probability strata of the
  # tested predictors' spread, m = N - 1 - (p - q) df, with R's own noncentral
  # F for the conditional part. 50000 strata keep this under a quarter of a
  # second; the quadrature's own error is then of order 1e-7.
  by_quadrature <- function(x, q, df2, m, c, M = 50000L) {
    w <- qchisq((seq_len(M) - 0.5) / M, df = m)
    mean(pf(x, q, df2, ncp = c * w))
  }
  # N, p, q, c, x: a joint test, the omnibus test, and a large prior study
  for (cfg in list(c(40, 5, 2, .5, 1), c(30, 4, 4, .03, 5), c(300, 6, 3, .5, 5),
                   c(300, 6, 3, .03, 1))) {
    N <- cfg[1]; p <- cfg[2]; q <- cfg[3]; c <- cfg[4]; x <- cfg[5]
    df2 <- N - p - 1; m <- N - 1 - (p - q)
    expect_lt(abs(BUCSS:::.p_random_f(x, q, df2, m, 1 / (1 + c)) -
                    by_quadrature(x, q, df2, m, c)), 1e-6)
  }
  for (q in c(1, 3))
    expect_equal(BUCSS:::.p_random_f(2.5, q, 30, 35, 1), pf(2.5, q, 30))
})

test_that("the frame leaves the assurance ceiling where it was", {
  # At a zero effect the mixture is the central F, so every example's stated
  # ceiling holds in both frames (to rounding in the last few bits).
  calls <- list(
    quote(ss_buc_reg_coef(t_observed = 3, N = 150, p = 3)),
    quote(ss_buc_R2(F_observed = 5, N = 150, p = 4)),
    quote(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2)))
  for (cl in calls) {
    fixed <- cl; fixed$random_predictors <- FALSE
    random <- cl; random$random_predictors <- TRUE
    expect_equal(attr(eval(random), "assurance_ceiling"),
                 attr(eval(fixed), "assurance_ceiling"), tolerance = 1e-12)
  }
})

test_that("the random-predictor plans are pinned, and actual_power is their exact power", {
  # The documented examples in both frames. The fixed plans are the defaults
  # pinned in test-oracles.R; they are repeated here with the flag written out.
  rc <- ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, random_predictors = TRUE)
  r2 <- ss_buc_R2(F_observed = 5, N = 150, p = 4, random_predictors = TRUE)
  rj <- ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                         random_predictors = TRUE)
  expect_equal(row(rc, "necessary_N"), 623)
  expect_equal(round(row(rc, "ncp_adjusted"), 3), 1.911)
  expect_equal(row(r2, "necessary_N"), 239)
  expect_equal(round(row(r2, "ncp_adjusted"), 3), 7.755)
  expect_equal(row(rj, "necessary_N"), 3898)
  expect_equal(round(row(rj, "ncp_adjusted"), 3), 0.372)
  fc <- ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, random_predictors = FALSE)
  f2 <- ss_buc_R2(F_observed = 5, N = 150, p = 4, random_predictors = FALSE)
  fj <- ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                         random_predictors = FALSE)
  expect_equal(c(row(fc, "necessary_N"), row(f2, "necessary_N"),
                 row(fj, "necessary_N")), c(624, 234, 3960))
  expect_equal(round(c(row(fc, "ncp_adjusted"), row(f2, "ncp_adjusted"),
                       row(fj, "ncp_adjusted")), 3), c(1.893, 7.816, 0.365))

  # actual_power is the mixture power of the reported plan: q numerator df,
  # n - p - 1 residual df, and m = n - 1 - (p - q) chi square df, at the
  # effect ncp_adjusted/N of the prior study of N = 150.
  mixture_power <- function(plan, q, p) {
    n <- row(plan, "necessary_N"); ncp <- row(plan, "ncp_adjusted")
    1 - BUCSS:::.p_random_f(qf(.95, q, n - p - 1), q, n - p - 1,
                            n - 1 - (p - q), 1 - ncp / (ncp + 150))
  }
  expect_equal(row(rc, "actual_power"), mixture_power(rc, 1, 3))
  expect_equal(row(r2, "actual_power"), mixture_power(r2, 4, 4))
  expect_equal(row(rj, "actual_power"), mixture_power(rj, 2, 4))
  for (plan in list(rc, r2, rj)) {
    expect_gte(row(plan, "actual_power"), .80)
    # and one participant fewer would fall short
    shorter <- plan
    shorter$value[shorter$term == "necessary_N"] <- row(plan, "necessary_N") - 1
    q <- attr(plan, "df_effect"); p <- row(plan, "p")
    expect_lt(mixture_power(shorter, q, p), .80)
  }
})

test_that("the frame is echoed as a numeric row and named in the footers", {
  for (rp in c(FALSE, TRUE)) {
    r <- ss_buc_R2(F_observed = 5, N = 150, p = 4, random_predictors = rp)
    expect_type(r$value, "double")
    expect_identical(row(r, "random_predictors"), as.numeric(rp))
    expect_identical(glance(r)$random_predictors, as.numeric(rp))
    footer <- if (rp) "Predictors: random (sampled along with the outcome)" else
      "Predictors: fixed by design"
    expect_output(print(r), footer, fixed = TRUE)
    if (requireNamespace("knitr", quietly = TRUE))
      expect_match(as.character(knitr::knit_print(r)), footer, fixed = TRUE)
  }
  # the row sits among the echoed inputs, after p (and p_joint)
  rj <- ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2)
  expect_identical(rj$term[which(rj$term == "p_joint") + 1L], "random_predictors")
  rc <- ss_buc_reg_coef(t_observed = 3, N = 150, p = 3)
  expect_identical(rc$term[which(rc$term == "p") + 1L], "random_predictors")
  # the 1.x positional contract is unaffected
  expect_identical(rc[[1]], row(rc, "necessary_N"))
  expect_identical(rc[[2]], row(rc, "ncp_adjusted"))
  # plans from the other planners have no such row and no such footer
  cor <- ss_buc_correlation(r_observed = .35, N = 100)
  expect_false("random_predictors" %in% cor$term)
  expect_false(any(grepl("Predictors:", capture.output(print(cor)), fixed = TRUE)))
})

test_that("the audit's publication rate, replay, and plot follow the plan's frame", {
  # The ceiling pin in test-sensitivity.R cannot tell the two families apart
  # (both are the central F at a zero effect), and neither can the replay, so
  # pin the publication rate, which the audit computes analytically at the
  # true parameter, against each family written out with its own degrees of
  # freedom: m = N - p for a coefficient, N - 1 for the omnibus test, and
  # N - 1 - (p - p_joint) for a joint test. The draws themselves are pinned
  # in the next block.
  true_ncp <- 12
  cases <- list(
    list(plan = function(rp) ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3,
                                             random_predictors = rp),
         q = 1, df2 = 56, m = 57),
    list(plan = function(rp) ss_buc_R2(F_observed = 8, N = 60, p = 3,
                                       random_predictors = rp),
         q = 3, df2 = 56, m = 59),
    list(plan = function(rp) ss_buc_reg_joint(F_observed = 12, N = 60, p = 5,
                                              p_joint = 2, random_predictors = rp),
         q = 2, df2 = 54, m = 56))
  for (cs in cases) for (rp in c(FALSE, TRUE)) {
    s <- ss_buc_sensitivity(cs$plan(rp), true_ncp = true_ncp, replications = 20,
                            seed = 2017)
    crit <- qf(.95, cs$q, cs$df2)
    expected <- if (rp)
      1 - BUCSS:::.p_random_f(crit, cs$q, cs$df2, cs$m,
                              1 - true_ncp / (true_ncp + 60)) else
      1 - pf(crit, cs$q, cs$df2, ncp = true_ncp)
    expect_equal(row(s, "publication_rate"), expected,
                 info = paste(attr(s, "design"), rp))
  }
  # A random plan's rebuilt call carries the flag as TRUE and replays.
  plan <- ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3, random_predictors = TRUE)
  spec <- BUCSS:::.BUCSS_DESIGN_SPECS[[attr(plan, "design")]]
  args <- BUCSS:::.bucss_design_args(BUCSS:::.bucss_design_inputs(plan), spec, NULL)
  expect_identical(args$random_predictors, TRUE)
  expect_true(BUCSS:::.replay_matches(do.call(ss_buc_reg_coef, args), plan))
  # and plot() redraws it in its own frame: at the plan's own assurance the
  # curve passes through the plan, not through its fixed-predictor twin.
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  drawn <- plot(plan, assurance = c(.6, .8))
  expect_identical(drawn$sample_size[2], row(plan, "necessary_N"))
  expect_identical(drawn$ncp_adjusted[2], row(plan, "ncp_adjusted"))
  twin <- ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3)
  expect_false(identical(row(twin, "ncp_adjusted"), row(plan, "ncp_adjusted")))
})

test_that("the simulated prior studies are drawn from the plan's mixture", {
  # The publication rate above is computed analytically, so it cannot see the
  # draws. Pin them two ways. First, under one seed each random family's draw
  # is exactly the chi square mixture written out by hand, with the degrees
  # of freedom the audit reads from the registry: each study's noncentrality
  # parameter is the partial f squared (ncp/N) times a chi square on m df,
  # drawn first, then the F statistics (a t is the root of an F on 1 df). The
  # correlation, whose registry df is N alone, shares the draw.
  M <- 50L
  ncp <- 12
  draw <- function(plan) {
    inputs <- BUCSS:::.bucss_design_inputs(plan)
    spec <- BUCSS:::.bucss_resolve_spec(
      BUCSS:::.BUCSS_DESIGN_SPECS[[attr(plan, "design")]], inputs)
    list(family = spec$family,
         x = withr::with_seed(2017, BUCSS:::.draw_prior_statistic(
           spec, spec$df(inputs), ncp, M)))
  }
  by_hand <- function(q, df2, m, N) withr::with_seed(2017,
    rf(M, df1 = q, df2 = df2, ncp = (ncp / N) * rchisq(M, df = m)))
  coef <- draw(ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3,
                               random_predictors = TRUE))
  expect_identical(coef$family, "f_from_t_random")
  expect_identical(coef$x, sqrt(by_hand(1, 56, 57, 60)))
  r2 <- draw(ss_buc_R2(F_observed = 8, N = 60, p = 3, random_predictors = TRUE))
  expect_identical(r2$family, "f_random")
  expect_identical(r2$x, by_hand(3, 56, 59, 60))
  joint <- draw(ss_buc_reg_joint(F_observed = 12, N = 60, p = 5, p_joint = 2,
                                 random_predictors = TRUE))
  expect_identical(joint$family, "f_random")
  expect_identical(joint$x, by_hand(2, 54, 56, 60))
  cor <- draw(ss_buc_correlation(r_observed = .35, N = 100))
  expect_identical(cor$family, "correlation")
  expect_identical(cor$x, sqrt(by_hand(1, 98, 99, 100)))

  # Second, the draws have the distribution the planner computes with. At a
  # prior N of 30 with five predictors, two of them tested, and a parameter
  # of 15, the empirical CDF of 20,000 draws at F = 5 lies within 4 standard
  # errors of the mixture's exact CDF, while the fixed frame's noncentral F
  # is about 30 standard errors away, so a draw that lost its chi square
  # mixing, or mixed on the wrong df, fails here.
  plan <- ss_buc_reg_joint(F_observed = 6, N = 30, p = 5, p_joint = 2,
                           random_predictors = TRUE)
  inputs <- BUCSS:::.bucss_design_inputs(plan)
  spec <- BUCSS:::.bucss_resolve_spec(
    BUCSS:::.BUCSS_DESIGN_SPECS[[attr(plan, "design")]], inputs)
  expect_identical(spec$df(inputs), c(2, 24, 26, 30))
  x <- withr::with_seed(2017, BUCSS:::.draw_prior_statistic(
    spec, spec$df(inputs), 15, 20000L))
  exact <- BUCSS:::.p_random_f(5, 2, 24, 26, 1 - 15 / (15 + 30))
  se <- sqrt(exact * (1 - exact) / 20000)
  expect_lt(abs(mean(x <= 5) - exact), 4 * se)
  expect_gt(abs(mean(x <= 5) - pf(5, 2, 24, ncp = 15)), 20 * se)
})

test_that("a plan without its random_predictors row is refused, not simulated in another frame", {
  # The row is required by the shape check (the registry's 'echo'), so a
  # trimmed random plan cannot be audited with prior studies drawn from the
  # fixed family while its replay runs in the planner's default frame.
  drop_row <- function(x, term) {
    y <- as.data.frame(x)[x$term != term, ]
    for (nm in setdiff(names(attributes(x)), c("names", "row.names")))
      attr(y, nm) <- attr(x, nm)
    y
  }
  shape <- "no longer has the shape of a single plan"
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  for (plan in list(
    ss_buc_reg_coef(t_observed = 3.5, N = 60, p = 3, random_predictors = TRUE),
    ss_buc_R2(F_observed = 8, N = 60, p = 3, random_predictors = TRUE),
    ss_buc_reg_joint(F_observed = 12, N = 60, p = 5, p_joint = 2,
                     random_predictors = TRUE))) {
    trimmed <- drop_row(plan, "random_predictors")
    expect_false(BUCSS:::.is_plan_shape(trimmed))
    expect_error(ss_buc_sensitivity(trimmed, true_ncp = 10, replications = 10,
                                    seed = 2017), shape, fixed = TRUE)
    expect_error(plot(trimmed), shape, fixed = TRUE)
    expect_error(planning_sentence(trimmed), shape, fixed = TRUE)
  }
})

test_that("the assurance guarantee holds under random predictors", {
  skip_on_cran()
  # Plans of several hundred participants, so integer granularity does not
  # push attainment above nominal; 4000 replications, 4 standard errors.
  reps <- 4000L
  se <- sqrt(.8 * .2 / reps)
  cases <- list(
    ss_buc_reg_coef(t_observed = 3.2, N = 150, p = 4, random_predictors = TRUE),
    ss_buc_R2(F_observed = 4, N = 200, p = 4, random_predictors = TRUE),
    ss_buc_reg_joint(F_observed = 6, N = 150, p = 5, p_joint = 2,
                     random_predictors = TRUE))
  for (plan in cases) {
    ncp <- row(plan, "ncp_adjusted")
    s <- ss_buc_sensitivity(plan, true_ncp = ncp * 1.4, replications = reps,
                            seed = 2017)
    expect_lt(abs(row(s, "attainment_with_refusals") - .80), 4 * se,
              label = attr(plan, "design"))
  }
})
