# The sphericity parameter 'epsilon' of the within-subjects and split-plot
# planners: the Muller and Barton (1989) approximation, in which both degrees
# of freedom and the noncentrality parameter of a within-subjects F are
# multiplied by epsilon, for the prior study and the planned one alike, with
# 'ncp_adjusted' kept on the sphericity scale. epsilon is echoed as a row for
# every effect it applies to (1 when sphericity is assumed) and is required
# there by the shape check.

pick <- function(x, term) x$value[x$term == term]
size <- function(x) x$value[x$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1]

eps_calls <- list(
  quote(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3)),
  quote(ss_buc_rm_anova(F_observed = 9, N = 40, levels_A = 2, levels_B = 3,
                        effect = "interaction")),
  quote(ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2)),
  quote(ss_buc_mixed_anova(F_observed = 8, N = 60, levels_between = 2,
                           levels_within = 3, effect = "within")),
  quote(ss_buc_mixed_anova(F_observed = 8, N = 61, levels_between = 2,
                           levels_within = 3, effect = "interaction")),
  quote(ss_buc_mixed_anova_general(F_observed = 8, N = 91, df_numerator = 2,
                                   num_groups = 3, effect = "within_only")),
  quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 6,
                                   num_groups = 3, effect = "between_within",
                                   df_num_within = 3))
)
with_eps <- function(cl, e) { cl$epsilon <- e; cl }

test_that("epsilon = 1, given or not, is the sphericity-assumed plan, echoed as a row", {
  for (cl in eps_calls) {
    plan <- eval(cl)
    expect_identical(eval(with_eps(cl, 1)), plan, info = deparse(cl))
    # The row sits after the design counts, just before alpha_prior.
    expect_identical(pick(plan, "epsilon"), 1, info = deparse(cl))
    expect_identical(which(plan$term == "epsilon"),
                     which(plan$term == "alpha_prior") - 1L, info = deparse(cl))
    # The sentence reads as it does for any plan that assumes sphericity.
    expect_false(grepl("sphericity", planning_sentence(plan), fixed = TRUE))
  }
  # Every within-subjects effect echoes it, including a one-df effect, where
  # it can only be 1.
  expect_identical(pick(ss_buc_rm_anova(F_observed = 12, N = 60, levels_A = 2), "epsilon"), 1)
  expect_identical(pick(ss_buc_rm_anova_general(F_observed = 6.5, N = 80, df_numerator = 1,
                                                assurance = .5), "epsilon"), 1)
  # The between-subjects effects involve no sphericity and echo no row.
  expect_false("epsilon" %in% ss_buc_mixed_anova(F_observed = 9, N = 60, levels_between = 2,
                                                 levels_within = 3, effect = "between")$term)
  expect_false("epsilon" %in% ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                                                         num_groups = 3,
                                                         effect = "between_only")$term)
})

test_that("the correction is the sphericity-assumed one on corrected df at epsilon * ncp", {
  # TM depends on the noncentrality parameter only through epsilon * ncp, and
  # the critical value is taken at the corrected degrees of freedom. The root
  # is found here independently of the package's closure.
  for (e in c(.9, .7, .55)) {
    r <- ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2,
                                 epsilon = e)
    crit <- qf(.95, e * 2, e * 58)
    tm <- function(x) {
      pw <- pf(crit, e * 2, e * 58, ncp = x, lower.tail = FALSE)
      (pw - pf(8, e * 2, e * 58, ncp = x, lower.tail = FALSE)) / pw
    }
    root <- uniroot(function(x) tm(x) - .8, c(0, 200), tol = 1e-12)$root
    expect_equal(e * pick(r, "ncp_adjusted"), root, tolerance = 1e-6)
    # the ceiling is 1 - p/alpha_prior with p at the corrected df
    expect_equal(attr(r, "assurance_ceiling"),
                 1 - pf(8, e * 2, e * 58, lower.tail = FALSE) / .05,
                 tolerance = 1e-12)
    # df attributes are the planned test's corrected df
    n <- size(r)
    expect_equal(attr(r, "df_effect"), e * 2)
    expect_equal(attr(r, "df_error"), e * 2 * (n - 1))
    # actual_power is the Muller-Barton power of the returned plan, and one
    # subject fewer falls short
    mb_power <- function(n) {
      d1 <- e * 2; d2 <- e * 2 * (n - 1)
      pf(qf(.95, d1, d2), d1, d2, ncp = e * (n / 30) * pick(r, "ncp_adjusted"),
         lower.tail = FALSE)
    }
    expect_equal(pick(r, "actual_power"), mb_power(n), tolerance = 1e-12)
    expect_lt(mb_power(n - 1), .8)
    expect_equal(pick(r, "epsilon"), e)
  }

  # The same for a balanced split-plot interaction, whose corrected df are
  # epsilon times (g - 1)(w - 1) and (N - g)(w - 1), and whose planned power
  # is per group.
  e <- .75
  r <- ss_buc_mixed_anova(F_observed = 6, N = 60, levels_between = 3,
                          levels_within = 4, effect = "interaction", epsilon = e)
  d1 <- e * 6; d2 <- e * (60 - 3) * 3
  crit <- qf(.95, d1, d2)
  tm <- function(x) {
    pw <- pf(crit, d1, d2, ncp = x, lower.tail = FALSE)
    (pw - pf(6, d1, d2, ncp = x, lower.tail = FALSE)) / pw
  }
  root <- uniroot(function(x) tm(x) - .8, c(0, 500), tol = 1e-12)$root
  expect_equal(e * pick(r, "ncp_adjusted"), root, tolerance = 1e-6)
  n <- size(r)
  p2 <- e * (3 * n - 3) * 3
  expect_equal(attr(r, "df_error"), p2)
  expect_equal(pick(r, "actual_power"),
               pf(qf(.95, d1, p2), d1, p2, ncp = e * (n / 20) * pick(r, "ncp_adjusted"),
                  lower.tail = FALSE), tolerance = 1e-12)
})

test_that("a smaller epsilon never plans a smaller study", {
  for (cl in eps_calls) {
    sizes <- vapply(c(1, .95, .85, .75, .65), function(e) {
      r <- tryCatch(eval(with_eps(cl, e)), bucss_zero_ncp = function(err) NULL)
      if (is.null(r)) Inf else size(r)
    }, numeric(1))
    expect_true(all(diff(sizes) >= 0), info = deparse(cl))
    expect_gt(sizes[5], sizes[1])
  }
})

test_that("an epsilon plan echoes its value, replays from its rows, and says so", {
  r <- ss_buc_mixed_anova(F_observed = 8, N = 61, levels_between = 2,
                          levels_within = 3, effect = "interaction", epsilon = .8)
  expect_equal(pick(r, "epsilon"), .8)
  expect_identical(which(r$term == "epsilon"), which(r$term == "alpha_prior") - 1L)
  spec <- BUCSS:::.BUCSS_DESIGN_SPECS[[attr(r, "design")]]
  args <- BUCSS:::.bucss_design_args(BUCSS:::.bucss_design_inputs(r), spec,
                                    attr(r, "effect"))
  expect_equal(args$epsilon, .8)
  expect_true(BUCSS:::.replay_matches(do.call(ss_buc_mixed_anova, args), r))
  expect_match(planning_sentence(r),
               "and for uncertainty, allowing for a departure from sphericity (epsilon = .8).",
               fixed = TRUE)
  flagship <- ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = 1.73 / 2)
  expect_identical(
    planning_sentence(flagship),
    paste0("A planned study with a total sample size of 162 attains 80% power with ",
           "80% assurance, after correcting the prior study's result for publication ",
           "bias (assuming its literature published results significant at the .05 ",
           "level) and for uncertainty, allowing for a departure from sphericity ",
           "(epsilon = .865). The plan is built on a planning effect size of ",
           "0.183, expressed as Cohen's f on the sphericity scale."))
  expect_equal(glance(flagship)$epsilon, 1.73 / 2)
  expect_output(print(flagship), "epsilon", fixed = TRUE)
})

test_that("epsilon is validated against the effect's own lower bound", {
  expect_error(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = .49),
               "at least 1/2", fixed = TRUE)
  expect_error(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = 1.01),
               "at most 1", fixed = TRUE)
  expect_error(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = NA),
               "single finite", fixed = TRUE)
  expect_error(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = c(.8, .9)),
               "single finite", fixed = TRUE)
  expect_error(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = "0.8"),
               "single finite", fixed = TRUE)
  expect_error(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = 0),
               "at least 1/2", fixed = TRUE)
  # two levels: always spherical
  expect_error(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 2, epsilon = .9),
               "always satisfies sphericity", fixed = TRUE)
  expect_error(ss_buc_rm_anova_general(F_observed = 6.5, N = 30, df_numerator = 1,
                                       epsilon = .9),
               "always satisfies sphericity", fixed = TRUE)
  expect_error(ss_buc_mixed_anova(F_observed = 7, N = 60, levels_between = 2,
                                  levels_within = 2, effect = "within", epsilon = .9),
               "always satisfies sphericity", fixed = TRUE)
  expect_error(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                                          num_groups = 3, effect = "between_within",
                                          df_num_within = 1, epsilon = .9),
               "always satisfies sphericity", fixed = TRUE)
  # the bound of a two-way within-subjects effect is its own df: the
  # interaction of a 3 x 4 design has 6 df, bound 1/6
  expect_s3_class(ss_buc_rm_anova(F_observed = 8, N = 40, levels_A = 3, levels_B = 4,
                                  effect = "interaction", epsilon = 1 / 6), "bucss_plan")
  # the interaction of a split-plot is bounded by the within factor's df, not
  # by its own numerator df: with 4 groups and 3 occasions df1 = 6, bound 1/2
  expect_error(ss_buc_mixed_anova(F_observed = 7, N = 60, levels_between = 4,
                                  levels_within = 3, effect = "interaction",
                                  epsilon = .3), "at least 1/2", fixed = TRUE)
  expect_s3_class(ss_buc_mixed_anova(F_observed = 7, N = 60, levels_between = 4,
                                     levels_within = 3, effect = "interaction",
                                     epsilon = .5), "bucss_plan")
  expect_error(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 6,
                                          num_groups = 3, effect = "between_within",
                                          df_num_within = 3, epsilon = .3),
               "at least 1/3", fixed = TRUE)
  expect_error(ss_buc_mixed_anova_general(F_observed = 8, N = 90, df_numerator = 2,
                                          num_groups = 3, effect = "within_only",
                                          epsilon = .4),
               "at least 1/2", fixed = TRUE)
  # the bound itself, written as the fraction, and a ratio of df a rounding
  # error beyond either end, are accepted as the value they meant
  expect_s3_class(ss_buc_rm_anova(F_observed = 12, N = 30, levels_A = 4, epsilon = 1/3),
                  "bucss_plan")
  expect_identical(ss_buc_rm_anova(F_observed = 12, N = 30, levels_A = 4,
                                   epsilon = 1/3 - 1e-12),
                   ss_buc_rm_anova(F_observed = 12, N = 30, levels_A = 4, epsilon = 1/3))
  expect_identical(ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3,
                                   epsilon = 1 + 1e-12),
                   ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3))
  # a one-df effect given 1 plus a rounding error is the default plan
  expect_identical(ss_buc_rm_anova(F_observed = 12, N = 30, levels_A = 2,
                                   epsilon = 1 - 1e-12),
                   ss_buc_rm_anova(F_observed = 12, N = 30, levels_A = 2))
  # the one-df message fits every one-df effect, not only a two-level factor
  expect_error(ss_buc_rm_anova_general(F_observed = 9, N = 30, df_numerator = 1,
                                       epsilon = .9),
               "(as for a two-level factor or a single contrast)", fixed = TRUE)
})

test_that("an epsilon a rounding error below 1 is read as exactly 1", {
  # the departure test is epsilon < 1, so a value a rounding error below 1
  # must become exactly 1 (the plan that assumes sphericity, bit for bit, with
  # no departure claimed), just as a value a rounding error above 1 does
  base <- ss_buc_rm_anova(F_observed = 9, N = 30, levels_A = 4)
  near <- ss_buc_rm_anova(F_observed = 9, N = 30, levels_A = 4, epsilon = 1 - 1e-10)
  expect_identical(near, base)
  expect_identical(pick(near, "epsilon"), 1)
  expect_identical(planning_sentence(near), planning_sentence(base))
  expect_false(grepl("sphericity", planning_sentence(near), fixed = TRUE))
  # a value genuinely below 1 is a departure, and the sentence never shows it
  # as 1: it gets the figures it needs
  r <- ss_buc_rm_anova(F_observed = 9, N = 30, levels_A = 4, epsilon = .9996)
  expect_equal(pick(r, "epsilon"), .9996)
  expect_match(planning_sentence(r),
               "allowing for a departure from sphericity (epsilon = .9996).",
               fixed = TRUE)
  expect_false(grepl("(epsilon = 1)", planning_sentence(r), fixed = TRUE))
  # the usual three figures are unchanged
  expect_match(planning_sentence(ss_buc_rm_anova(F_observed = 9, N = 30, levels_A = 4,
                                                 epsilon = 1.73 / 2)),
               "(epsilon = .865)", fixed = TRUE)
})

test_that("an epsilon reported at the lower bound to three or four decimals plans", {
  # the bound of a 4-level factor is 1/3; the error names it as about 0.3333,
  # so that value, and a published .333, must be read as the bound itself
  at_bound <- ss_buc_rm_anova(F_observed = 30, N = 30, levels_A = 4, epsilon = 1/3)
  expect_identical(ss_buc_rm_anova(F_observed = 30, N = 30, levels_A = 4,
                                   epsilon = .3333), at_bound)
  expect_identical(ss_buc_rm_anova(F_observed = 30, N = 30, levels_A = 4,
                                   epsilon = .333), at_bound)
  expect_identical(pick(at_bound, "epsilon"), 1/3)
  # the same holds for the split-plot interaction's own bound (1/(levels_within - 1))
  sp <- ss_buc_mixed_anova(F_observed = 12, N = 60, levels_between = 2,
                           levels_within = 4, effect = "interaction", epsilon = 1/3)
  expect_identical(ss_buc_mixed_anova(F_observed = 12, N = 60, levels_between = 2,
                                      levels_within = 4, effect = "interaction",
                                      epsilon = .333), sp)
  # beyond half a unit in the third decimal place it is still refused, and
  # the value the message names is one the check accepts
  err <- tryCatch(ss_buc_rm_anova(F_observed = 30, N = 30, levels_A = 4, epsilon = .3328),
                  error = function(e) conditionMessage(e))
  expect_match(err, "at least 1/3 (about 0.3333)", fixed = TRUE)
  named <- as.numeric(sub(".*\\(about ([0-9.]+)\\).*", "\\1", err))
  expect_identical(ss_buc_rm_anova(F_observed = 30, N = 30, levels_A = 4,
                                   epsilon = named), at_bound)
  for (k in 2:12) {
    msg <- tryCatch(BUCSS:::.check_epsilon(1 / k - 1e-3, k), error = conditionMessage)
    shown <- as.numeric(sub(".*\\(about ([0-9.]+)\\).*", "\\1", msg))
    expect_identical(BUCSS:::.check_epsilon(shown, k), max(shown, 1 / k))
  }
  # the allowance never exceeds half the bound, so zero and negative values are
  # refused however many degrees of freedom the effect has (with a fixed 5e-4
  # they passed from 2,000 df on and were clamped to the bound)
  for (k in c(1000, 2000, 3000)) {
    expect_error(BUCSS:::.check_epsilon(0, k), paste0("at least 1/", k), fixed = TRUE)
    expect_error(BUCSS:::.check_epsilon(-1e-4, k), paste0("at least 1/", k), fixed = TRUE)
    expect_identical(BUCSS:::.check_epsilon(signif(1 / k, 4), k), 1 / k)
  }
  expect_error(ss_buc_rm_anova_general(F_observed = 30, N = 50, df_numerator = 3000,
                                       epsilon = 0),
               "at least 1/3000", fixed = TRUE)
})

test_that("the between-subjects effects ignore epsilon, with a warning", {
  base <- ss_buc_mixed_anova(F_observed = 9, N = 60, levels_between = 2,
                             levels_within = 3, effect = "between")
  expect_warning(r <- ss_buc_mixed_anova(F_observed = 9, N = 60, levels_between = 2,
                                         levels_within = 3, effect = "between",
                                         epsilon = .7),
                 "subject means, which involves no sphericity, so 'epsilon' was ignored",
                 fixed = TRUE)
  expect_identical(r, base)
  expect_no_warning(r1 <- ss_buc_mixed_anova(F_observed = 9, N = 60, levels_between = 2,
                                             levels_within = 3, effect = "between",
                                             epsilon = 1))
  expect_identical(r1, base)
  expect_error(ss_buc_mixed_anova(F_observed = 9, N = 60, levels_between = 2,
                                  levels_within = 3, effect = "between", epsilon = NA),
               "single finite", fixed = TRUE)
  base_g <- ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                                       num_groups = 3, effect = "between_only")
  expect_warning(rg <- ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                                                  num_groups = 3, effect = "between_only",
                                                  epsilon = .7),
                 "'effect = \"between_only\"' is tested on subject means", fixed = TRUE)
  expect_identical(rg, base_g)
})

test_that("a result significant under sphericity can be nonsignificant after correction", {
  expect_s3_class(ss_buc_rm_anova(F_observed = 3.3, N = 30, levels_A = 4,
                                  assurance = .5), "bucss_plan")
  expect_error(ss_buc_rm_anova(F_observed = 3.3, N = 30, levels_A = 4, epsilon = .4),
               class = "bucss_nonsignificant")
})

test_that("an unbalanced split-plot prior is judged at its own corrected degrees of freedom", {
  # 2 groups x 3 occasions with N = 7, which the planner reduces to 6: at
  # epsilon = .8 the prior study's own corrected df are 1.6 and 8, the
  # balanced design's 1.6 and 6.4. A statistic between the two critical
  # values is significant at its own df, and the refusal says so and names an
  # alpha_prior at which the planner plans.
  e <- .8
  f_mid <- (qf(.95, e * 2, e * 10) + qf(.95, e * 2, e * 8)) / 2
  for (a in c(.5, .8)) {
    err <- tryCatch(ss_buc_mixed_anova(F_observed = f_mid, N = 7, levels_between = 2,
                                       levels_within = 3, effect = "within",
                                       epsilon = e, assurance = a),
                    bucss_nonsignificant = function(err) err)
    expect_s3_class(err, "bucss_nonsignificant")
    msg <- conditionMessage(err)
    expect_match(msg, "significant at the prior study's own degrees of freedom", fixed = TRUE)
    p_own <- pf(f_mid, e * 2, e * 10, lower.tail = FALSE)
    expect_match(msg, paste0("(p = ", sub("^0", "", format(signif(p_own, 3))), ")"),
                 fixed = TRUE)
    named <- as.numeric(sub(".*an 'alpha_prior' of ([.0-9]+) or more.*", "\\1", msg))
    expect_s3_class(ss_buc_mixed_anova(F_observed = f_mid, N = 7, levels_between = 2,
                                       levels_within = 3, effect = "within",
                                       epsilon = e, assurance = a, alpha_prior = named),
                    "bucss_plan")
  }
})

test_that("the registry scales df by epsilon and reproduces the ceiling", {
  plans <- list(
    ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3, epsilon = .7),
    ss_buc_rm_anova(F_observed = 9, N = 40, levels_A = 2, levels_B = 3,
                    effect = "interaction", epsilon = .7),
    ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2, epsilon = .7),
    ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                       levels_within = 3, effect = "within", epsilon = .7),
    ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                       levels_within = 3, effect = "interaction", epsilon = .7),
    ss_buc_mixed_anova_general(F_observed = 20, N = 60, num_groups = 2,
                               df_numerator = 2, effect = "within_only", epsilon = .7),
    ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 6, num_groups = 3,
                               effect = "between_within", df_num_within = 3,
                               epsilon = .6))
  for (o in plans) {
    spec <- BUCSS:::.BUCSS_DESIGN_SPECS[[attr(o, "design")]]
    inputs <- BUCSS:::.bucss_design_inputs(o)
    effect <- attr(o, "effect")
    e <- inputs$epsilon
    df <- if (length(formals(spec$df)) > 1L) spec$df(inputs, effect) else spec$df(inputs)
    # the prior study's own df times epsilon (the priors are balanced, so the
    # planner's ceiling is at those df too)
    expect_equal(1 - pf(inputs$F_observed, df[1], df[2]),
                 .05 * (1 - attr(o, "assurance_ceiling")), tolerance = 1e-9,
                 info = attr(o, "design"))
    expect_identical(spec$ncp_scale(inputs), e)
    args <- BUCSS:::.bucss_design_args(inputs, spec, effect)
    replay <- do.call(get(spec$planner, envir = asNamespace("BUCSS")), args)
    expect_true(BUCSS:::.replay_matches(replay, o), info = attr(o, "design"))
  }
  # a between-subjects plan has no epsilon row and a scale of 1
  b <- ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                          levels_within = 3, effect = "between")
  spec <- BUCSS:::.BUCSS_DESIGN_SPECS[[attr(b, "design")]]
  expect_identical(spec$ncp_scale(BUCSS:::.bucss_design_inputs(b)), 1)
})

test_that("a plan with its epsilon row removed is refused, not replayed under sphericity", {
  shape <- "no longer has the shape of a single plan"
  drop_eps <- function(p) {
    out <- as.data.frame(p)[p$term != "epsilon", , drop = FALSE]
    rownames(out) <- NULL
    for (a in c("design", "sample_size_unit", "effect", "assurance_ceiling",
                "df_effect", "df_error"))
      attr(out, a) <- attr(p, a)
    class(out) <- class(p)
    out
  }
  plans <- list(
    ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3, epsilon = .7),
    ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3),
    ss_buc_rm_anova_general(F_observed = 8, N = 30, df_numerator = 2, epsilon = .7),
    ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                       levels_within = 3, effect = "within", epsilon = .7),
    ss_buc_mixed_anova(F_observed = 20, N = 60, levels_between = 2,
                       levels_within = 3, effect = "interaction"),
    ss_buc_mixed_anova_general(F_observed = 20, N = 60, num_groups = 2,
                               df_numerator = 2, effect = "within_only", epsilon = .7),
    ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 6, num_groups = 3,
                               effect = "between_within", df_num_within = 3,
                               epsilon = .6))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  for (p in plans) {
    trimmed <- drop_eps(p)
    info <- paste(attr(p, "design"), attr(p, "effect"))
    expect_true(BUCSS:::.is_plan_shape(p), info = info)
    expect_false(BUCSS:::.is_plan_shape(trimmed), info = info)
    expect_error(planning_sentence(trimmed), shape, fixed = TRUE)
    expect_error(glance(trimmed), shape, fixed = TRUE)
    expect_error(ss_buc_sensitivity(trimmed, true_ncp = 5, replications = 10, seed = 2017),
                 shape, fixed = TRUE)
    expect_error(plot(trimmed), shape, fixed = TRUE)
  }
})

test_that("ss_buc_sensitivity draws prior studies under the same approximation", {
  plan <- ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3, epsilon = .7)
  s <- ss_buc_sensitivity(plan, true_ncp = 6, replications = 200, seed = 2017)
  crit <- qf(.95, .7 * 2, .7 * 58)
  expect_equal(pick(s, "publication_rate"),
               pf(crit, .7 * 2, .7 * 58, ncp = .7 * 6, lower.tail = FALSE))
  # the split-plot interaction, at the prior study's own N
  plan2 <- ss_buc_mixed_anova(F_observed = 20, N = 61, levels_between = 2,
                              levels_within = 4, effect = "interaction", epsilon = .6)
  s2 <- ss_buc_sensitivity(plan2, true_ncp = 8, replications = 50, seed = 2017)
  d1 <- .6 * 3; d2 <- .6 * (61 - 2) * 3
  expect_equal(pick(s2, "publication_rate"),
               pf(qf(.95, d1, d2), d1, d2, ncp = .6 * 8, lower.tail = FALSE))
})

test_that("the assurance guarantee holds under epsilon within the approximation", {
  # A check of internal consistency only: the prior studies are drawn under
  # the same approximation the planner corrects under (see the help pages for
  # its accuracy against raw data).
  skip_on_cran()
  reps <- 4000L
  se <- sqrt(.8 * .2 / reps)
  cases <- list(
    ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3, epsilon = .7),
    ss_buc_mixed_anova(F_observed = 8, N = 60, levels_between = 2,
                       levels_within = 4, effect = "interaction", epsilon = .6))
  for (plan in cases) {
    ncp <- pick(plan, "ncp_adjusted")
    s <- ss_buc_sensitivity(plan, true_ncp = ncp * 1.4, replications = reps,
                            seed = 2017)
    expect_lt(abs(pick(s, "attainment_with_refusals") - .80), 4 * se)
  }
})

test_that("plot() replays an epsilon plan", {
  plan <- ss_buc_rm_anova(F_observed = 8, N = 30, levels_A = 3, epsilon = .7)
  pdf(NULL)
  on.exit(dev.off())
  v <- plot(plan, n_points = 8)
  expect_equal(attr(v, "assurance_ceiling"), attr(plan, "assurance_ceiling"))
  expect_true(all(diff(v$ncp_adjusted[!is.na(v$ncp_adjusted)]) <= 0))
  # the replay at the plan's own assurance is the plan
  at_plan <- plot(plan, assurance = c(.7, .8), which = "size")
  expect_equal(at_plan$sample_size[at_plan$assurance == .8], size(plan))
})

test_that("the deprecated 1.x wrappers have no epsilon and plan under sphericity", {
  for (f in c("ss.power.wa", "ss.power.wa.general", "ss.power.spa",
              "ss.power.spa.general"))
    expect_false("epsilon" %in% names(formals(get(f, envir = asNamespace("BUCSS")))))
  old <- suppressWarnings(ss.power.wa(F.observed = 5, N = 60, levels.A = 2, levels.B = 3,
                                      effect = "factor.B"))
  expect_identical(pick(old, "epsilon"), 1)
  # epsilon is the last formal of each planner, so no 2.0.0 positional call
  # can reach it by accident
  for (f in c("ss_buc_rm_anova", "ss_buc_rm_anova_general", "ss_buc_mixed_anova",
              "ss_buc_mixed_anova_general")) {
    fm <- formals(get(f, envir = asNamespace("BUCSS")))
    expect_identical(names(fm)[length(fm)], "epsilon", info = f)
    expect_identical(fm$epsilon, 1, info = f)
  }
})
