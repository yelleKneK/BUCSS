# The dot-named 1.x functions are deprecated (not defunct) as of 2.0.0: each
# warns once per session, then forwards to its ss_buc_* replacement and returns
# that function's bucss_power object. The bucss_power object keeps 1.x positional
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
    new = quote(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2, num_groups = 3, effect = "between_only", df_num_within = 3))),
  reg1 = list(
    old = quote(ss.power.reg1(t.observed = 3, N = 150, p = 3)),
    new = quote(ss_buc_reg_coef(t_observed = 3, N = 150, p = 3))),
  reg_all = list(
    old = quote(ss.power.reg.all(F.observed = 5, N = 150, p = 4)),
    new = quote(ss_buc_R2(F_observed = 5, N = 150, p = 4))),
  reg_joint = list(
    old = quote(ss.power.reg.joint(F.observed = 5, N = 150, p = 4, p.joint = 2)),
    new = quote(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2)))
)

test_that("each deprecated function forwards to its ss_buc_* replacement", {
  for (nm in names(legacy_pairs)) {
    old_res <- suppressWarnings(eval(legacy_pairs[[nm]]$old))
    new_res <- eval(legacy_pairs[[nm]]$new)
    expect_s3_class(old_res, "bucss_power")
    expect_equal(old_res$value, new_res$value, info = nm)
    expect_identical(attr(old_res, "design"), attr(new_res, "design"), info = nm)
    expect_identical(attr(old_res, "effect"), attr(new_res, "effect"), info = nm)
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

test_that("bucss_power keeps 1.x positional extraction working", {
  r <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80, desired_power = .80)
  expect_identical(r[[1]], r$value[r$term %in% BUCSS:::.BUCSS_SIZE_TERMS][1])
  expect_identical(r[[2]], r$value[r$term == "ncp_adjusted"])
  # Named/column access is unaffected by the [[ method.
  expect_identical(r[["value"]], r$value)
  expect_identical(r[["term"]], r$term)
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
    "bucss_power")
})
