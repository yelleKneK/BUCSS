# The dot-named 1.x functions are deprecated (not defunct) as of 2.0.0: each
# warns once per session, then forwards to its ss_buc_* replacement and returns
# that function's bucss_plan object. The bucss_plan object keeps 1.x positional
# extraction working (result[[1]] = sample size, result[[2]] = adjusted NCP).

# Map each deprecated function to the ss_buc_* call it should reproduce. Each
# entry is list(old = <call>, new = <call>) as quoted expressions evaluated below.
legacy_pairs <- list(
  it = list(
    old = quote(ss.power.it(t.observed = 3, n = 20, alpha.prior = .05, assurance = .80, power = .80)),
    new = quote(ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = .05, assurance = .80, desired_power = .80))),
  dt = list(
    old = quote(ss.power.dt(t.observed = 3, N = 40)),
    new = quote(ss_buc_paired_t(t_observed = 3, N = 40))),
  ba_oneway = list(
    old = quote(ss.power.ba(F.observed = 5, N = 120, levels.A = 4)),
    new = quote(ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 4))),
  ba_factorial = list(
    old = quote(ss.power.ba(F.observed = 5, N = 120, levels.A = 2, levels.B = 3, effect = "factor.B")),
    new = quote(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3, effect = "factor_B"))),
  ba_general = list(
    old = quote(ss.power.ba.general(F.observed = 5, N = 120, cells = 6, df.numerator = 2, df.denominator = 114)),
    new = quote(ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6, df_numerator = 2, df_denominator = 114))),
  wa = list(
    old = quote(ss.power.wa(F.observed = 5, N = 60, levels.A = 2, levels.B = 3, effect = "factor.B")),
    new = quote(ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3, effect = "factor_B"))),
  wa_general = list(
    old = quote(ss.power.wa.general(F.observed = 6.5, N = 80, df.numerator = 1, assurance = .50)),
    new = quote(ss_buc_rm_anova_general(F_observed = 6.5, N = 80, df_numerator = 1, assurance = .50))),
  spa = list(
    old = quote(ss.power.spa(F.observed = 5, N = 60, levels.between = 2, levels.within = 3, effect = "within")),
    new = quote(ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2, levels_within = 3, effect = "within"))),
  spa_general = list(
    old = quote(ss.power.spa.general(F.observed = 5, N = 90, df.numerator = 2, num.groups = 3, effect = "between.only", df.num.within = 3)),
    new = quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2, num_groups = 3, effect = "between_only"))),
  reg1 = list(
    old = quote(ss.power.reg1(t.observed = 3, N = 150, p = 3)),
    new = quote(ss_buc_reg_coef(t_observed = 3, N = 150, p = 3,
                                random_predictors = FALSE))),
  reg_all = list(
    old = quote(ss.power.reg.all(F.observed = 5, N = 150, p = 4)),
    new = quote(ss_buc_R2(F_observed = 5, N = 150, p = 4,
                          random_predictors = FALSE))),
  reg_joint = list(
    old = quote(ss.power.reg.joint(F.observed = 5, N = 150, p = 4, p.joint = 2)),
    new = quote(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                                 random_predictors = FALSE)))
)

test_that("each deprecated function forwards to its ss_buc_* replacement", {
  for (nm in names(legacy_pairs)) {
    old_res <- suppressWarnings(eval(legacy_pairs[[nm]]$old))
    new_res <- eval(legacy_pairs[[nm]]$new)
    expect_s3_class(old_res, "bucss_plan")
    expect_equal(old_res$value, new_res$value, info = nm)
    expect_identical(attr(old_res, "design"), attr(new_res, "design"), info = nm)
    expect_identical(attr(old_res, "effect"), attr(new_res, "effect"), info = nm)
  }
})

test_that("the regression wrappers keep the fixed-predictor frame of 1.x", {
  # BUCSS 1.x computed every regression plan with fixed predictors. The
  # wrappers pass random_predictors = FALSE explicitly, so a later change of
  # the planners' default cannot move them.
  for (nm in c("reg1", "reg_all", "reg_joint")) {
    old_res <- suppressWarnings(eval(legacy_pairs[[nm]]$old))
    expect_identical(old_res$value[old_res$term == "random_predictors"], 0,
                     info = nm)
  }
  for (f in c("ss.power.reg1", "ss.power.reg.all", "ss.power.reg.joint"))
    expect_true(any(grepl("random_predictors = FALSE",
                          deparse(body(get(f, envir = asNamespace("BUCSS")))),
                          fixed = TRUE)), info = f)
})

test_that("the wrappers forward \"No Publication Filter\" and the refusal of a low target power", {
  # The wrappers pass 'alpha.prior' and 'power' to the planners untouched, so
  # the rules of 2026-10-13 reach every 1.x call: 'alpha.prior' may be the
  # string "No Publication Filter", which is 1 (computed as exactly 1), and a
  # 'power' at or below 'alpha.planned' is not allowable (1.x returned its
  # smallest size for it; note that a 1.x 'power = 1' is read as 1 percent).
  for (nm in names(legacy_pairs)) {
    with_old <- function(...) {
      call <- legacy_pairs[[nm]]$old
      extra <- list(...)
      for (arg in names(extra)) call[[arg]] <- extra[[arg]]
      suppressWarnings(eval(call))
    }
    one <- with_old(alpha.prior = 1)
    named <- with_old(alpha.prior = "No Publication Filter")
    expect_s3_class(named, "bucss_plan")
    expect_identical(named, one, info = nm)
    expect_identical(named$value[named$term == "alpha_prior"], 1, info = nm)
    # The 1.x positional reading still works on it.
    expect_identical(named[[1]], one[[1]], info = nm)
    expect_identical(named[[2]], one[[2]], info = nm)
    # ...and it is the replacement's plan at 'alpha_prior = 1'.
    new_call <- legacy_pairs[[nm]]$new
    new_call$alpha_prior <- 1
    expect_identical(named, eval(new_call), info = nm)
    expect_error(with_old(alpha.prior = "none"),
                 "or the string \"No Publication Filter\", which means 1", fixed = TRUE,
                 info = nm)

    expect_error(with_old(power = .03),
                 "A 'desired_power' at or below 'alpha_planned' is not allowable:",
                 fixed = TRUE, info = nm)
    expect_error(with_old(power = 1),
                 "it was read as a percentage (1 percent, that is, .01).", fixed = TRUE,
                 info = nm)
    expect_error(with_old(power = .09, alpha.planned = .10),
                 "'desired_power' is .09 and 'alpha_planned' is .10.", fixed = TRUE,
                 info = nm)
  }
})

test_that("deprecation warning fires once per session and names the replacement", {
  # Clear the per-session warning record so the first call below warns.
  rm(list = ls(envir = BUCSS:::.bucss_deprecated_warned),
     envir = BUCSS:::.bucss_deprecated_warned)
  w <- testthat::capture_warnings(ss.power.it(t.observed = 3, n = 20))
  expect_match(w, "deprecated", all = FALSE)
  expect_match(w, "ss_buc_independent_t", all = FALSE)
  # A subsequent call in the same session is silent.
  expect_no_warning(ss.power.it(t.observed = 3, n = 20))
})

test_that("ss.power.spa.general takes df.num.within with any effect, warning only once", {
  # 1.x accepted 'df.num.within' with every effect, and the web apps pass it
  # on every call. The wrapper forwards it only for the between-within
  # interaction, so a 1.x call with another effect raises the deprecation
  # warning and nothing else (not the 2.x planner's "'df_num_within' is used
  # only when" warning).
  if (exists("ss.power.spa.general", envir = BUCSS:::.bucss_deprecated_warned,
             inherits = FALSE))
    rm("ss.power.spa.general", envir = BUCSS:::.bucss_deprecated_warned)
  record <- new.env(parent = emptyenv())
  record$classes <- list()
  collect <- function(code) withCallingHandlers(code, warning = function(w) {
    record$classes <- c(record$classes, list(class(w)))
    invokeRestart("muffleWarning")
  })
  between <- collect(ss.power.spa.general(F.observed = 5, N = 90, df.numerator = 2,
                                          num.groups = 3, effect = "between.only",
                                          df.num.within = 3))
  expect_length(record$classes, 1L)
  expect_true("bucss_deprecatedWarning" %in% record$classes[[1]])
  # After the once-per-session warning, the same calls are silent.
  expect_no_warning(ss.power.spa.general(F.observed = 5, N = 90, df.numerator = 2,
                                         num.groups = 3, effect = "between.only",
                                         df.num.within = 3))
  expect_no_warning(within <- ss.power.spa.general(F.observed = 5, N = 90,
                                                   df.numerator = 2, num.groups = 3,
                                                   effect = "within.only",
                                                   df.num.within = 3))
  # ...and the value it ignores does not change the plan.
  expect_identical(between, ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                                                       num_groups = 3, effect = "between_only"))
  expect_identical(within, ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                                                      num_groups = 3, effect = "within_only"))
})

test_that("ss.power.ba dispatches on levels.B", {
  one_way <- suppressWarnings(ss.power.ba(F.observed = 5, N = 120, levels.A = 4))
  two_way <- suppressWarnings(ss.power.ba(F.observed = 5, N = 120, levels.A = 2,
                                          levels.B = 3, effect = "factor.B"))
  expect_identical(attr(one_way, "design"), "One-way between-subjects ANOVA")
  expect_identical(attr(two_way, "design"), "Two-way between-subjects ANOVA")
})

test_that("1.x dotted effect values are translated to 2.x snake_case", {
  r <- suppressWarnings(ss.power.spa.general(F.observed = 5, N = 90,
                                             df.numerator = 2, num.groups = 3,
                                             effect = "between.within",
                                             df.num.within = 2))
  expect_identical(attr(r, "effect"), "between_within")
})

test_that("the legacy 'step' argument is accepted and ignored", {
  with_step <- suppressWarnings(ss.power.dt(t.observed = 3, N = 40, step = .01))
  without    <- suppressWarnings(ss.power.dt(t.observed = 3, N = 40))
  expect_equal(with_step$value, without$value)
})

test_that("bucss_plan keeps 1.x positional extraction working", {
  r <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, desired_power = .80)
  expect_identical(r[[1]], r$value[r$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1])
  expect_identical(r[[2]], r$value[r$term == "ncp_adjusted"])
  # Named/column access is unaffected by the [[ method.
  expect_identical(r[["value"]], r$value)
  expect_identical(r[["term"]], r$term)
})

test_that("the legacy positional [[ applies only to an object with one plan's shape", {
  # CL15-1 and PC-3 of the 2026-10-12 triage. Tools built on vctrs
  # (dplyr::bind_rows(), dplyr::filter(), purrr::map_dfr()) keep the class on
  # a stack of plans or a trimmed plan. The legacy reading then returned the
  # first plan's size for [[1]] (222 for this stack, the second plan's 31
  # lost), a vector of two parameters for [[2]], and NA or an empty vector for
  # a plan with its row removed, with no condition raised. Such an object is
  # now read as the plain data frame print() already shows it to be.
  a <- ss_buc_independent_t(t_observed = 3, n = c(50, 50))
  b <- ss_buc_one_way_anova(F_observed = 8, N = 60, levels_A = 3)
  expect_identical(a[[1]], 222)
  expect_identical(b[[1]], 31)
  keep_class <- function(df, like) {
    for (nm in setdiff(names(attributes(like)), c("names", "row.names")))
      attr(df, nm) <- attr(like, nm)
    df
  }
  stack <- keep_class(rbind(as.data.frame(a), as.data.frame(b)), a)
  expect_s3_class(stack, "bucss_plan")
  expect_identical(nrow(stack), nrow(a) + nrow(b))
  expect_identical(stack[[1]], c(a$term, b$term))
  expect_identical(stack[[2]], c(a$value, b$value))
  expect_length(stack[[1]], nrow(a) + nrow(b))

  no_ncp <- keep_class(as.data.frame(a)[a$term != "ncp_adjusted", ], a)
  expect_s3_class(no_ncp, "bucss_plan")
  expect_identical(no_ncp[[2]], a$value[a$term != "ncp_adjusted"])
  no_size <- keep_class(as.data.frame(a)[a$term %in% c("ncp_adjusted", "assurance"), ], a)
  expect_identical(no_size[[1]], c("ncp_adjusted", "assurance"))

  # A plan with a column added is not one plan's shape either.
  wide <- keep_class(cbind(as.data.frame(a), extra = 1), a)
  expect_identical(wide[[1]], a$term)

  # One plan is read as 1.x read it, by a double 1 or 2 only, and a plan with
  # its rows reordered (still one plan's shape) is too.
  expect_identical(a[[1]], a$value[a$term == "necessary_n_per_group"])
  expect_identical(a[[2]], a$value[a$term == "ncp_adjusted"])
  expect_identical(b[[1]], b$value[b$term == "necessary_n_per_group"])
  shuffled <- keep_class(as.data.frame(a)[rev(seq_len(nrow(a))), ], a)
  expect_identical(shuffled[[1]], a[[1]])
  expect_identical(shuffled[[2]], a[[2]])
  expect_identical(a[[1L]], a$term)
  # The relabeled and wrapper results are plans like any other.
  one <- ss_buc_one_sample_t(t_observed = 3, N = 40)
  expect_identical(one[[1]], one$value[one$term == "necessary_N"])
  old <- suppressWarnings(ss.power.dt(t.observed = 3, N = 40))
  expect_identical(old[[1]], 255)
  expect_identical(old[[2]], old$value[old$term == "ncp_adjusted"])
})

test_that("a malformed 'effect' in a 1.x wrapper gets one message, in the dotted spelling", {
  # GF-3 of the 2026-10-12 triage: these reached match.arg() and failed with
  # "'arg' must be of length 1" or "'arg' must be NULL or a character vector".
  # The message names the argument and the 1.x choices the user typed.
  dotted <- function(choices)
    paste0("'effect' should be one of \"", choices[1], "\", \"", choices[2], "\", or \"",
           choices[3], "\", given as a single character string.")
  ab <- dotted(c("factor.A", "factor.B", "interaction"))
  for (bad in list(c("factor.B", "factor.A"), 2, NA, NULL, "bogus", "factor_B")) {
    label <- paste(deparse(bad), collapse = " ")
    expect_error(suppressWarnings(ss.power.ba(F.observed = 12, N = 120, levels.A = 2,
                                              levels.B = 3, effect = bad)),
                 ab, fixed = TRUE, label = label)
    expect_error(suppressWarnings(ss.power.wa(F.observed = 12, N = 60, levels.A = 2,
                                              levels.B = 3, effect = bad)),
                 ab, fixed = TRUE, label = label)
    expect_error(suppressWarnings(ss.power.spa(F.observed = 12, N = 60, levels.between = 2,
                                               levels.within = 3, effect = bad)),
                 dotted(c("between", "within", "interaction")), fixed = TRUE, label = label)
    expect_error(suppressWarnings(ss.power.spa.general(F.observed = 12, N = 90,
                                                       df.numerator = 2, num.groups = 3,
                                                       effect = bad)),
                 dotted(c("between.only", "within.only", "between.within")), fixed = TRUE,
                 label = label)
  }
  # The untouched default still means the first choice, and a dotted value is
  # translated as before.
  expect_identical(
    attr(suppressWarnings(ss.power.wa(F.observed = 12, N = 60, levels.A = 2, levels.B = 3)),
         "effect"), "factor_A")
  expect_identical(
    attr(suppressWarnings(ss.power.spa.general(F.observed = 12, N = 90, df.numerator = 2,
                                               num.groups = 3)), "effect"), "between_only")
  expect_identical(
    attr(suppressWarnings(ss.power.ba(F.observed = 12, N = 120, levels.A = 2, levels.B = 3,
                                      effect = "factor.B")), "effect"), "factor_B")
})

test_that("ss.power.ba keeps the 1.x refusal to plan a two-factor effect alone", {
  # 1.2.1 stopped with "You cannot select 'effect=interaction' if you do not
  # specify 'levels.B'". Dispatching to the one-way planner instead would
  # silently return the Factor A plan to a user who asked for something else.
  expect_error(
    suppressWarnings(ss.power.ba(F.observed = 5, N = 120, levels.A = 4,
                                 effect = "interaction")),
    "levels.B", fixed = TRUE)
  expect_error(
    suppressWarnings(ss.power.ba(F.observed = 5, N = 120, levels.A = 4,
                                 effect = "factor.B")),
    "levels.B", fixed = TRUE)
  # the one-way dispatch itself still works
  expect_s3_class(
    suppressWarnings(ss.power.ba(F.observed = 5, N = 120, levels.A = 4)),
    "bucss_plan")
})
