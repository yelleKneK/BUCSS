# The output layer after the 2026-09-28 deep QC audit: planning_sentence()
# names only the corrections a plan made (N-02), counts carry thousands
# separators (N-13), a proportion below 1 never prints as 100% (N-14/R-09),
# the broom verbs and the sentence refuse an object that is no longer one plan
# (N-03), renaming drops the class (N-12/L-03), glance() begins like tidy()
# (A-02), ss_buc_sensitivity() takes a T squared MANOVA plan (M-01) and
# reports its rates and standard errors as documented (N-05, N-06, N-07), and
# plot() draws a plan whose ceiling is below .5 (N-04). The follow-ups from
# the review of that fix wave (2026-09-29) are here too: a plan with an added
# column or reordered rows, stacks refused by the audit and the figure,
# stacked audits, an alpha of .10, and the warning summary's counts.

val <- function(x, term) x$value[x$term == term]
size <- function(x) x$value[x$term %in% BUCSS:::.BUCSS_SIZE_TERMS]

# Evaluate 'code', muffling and recording every warning it raises.
collect_warnings <- function(code) {
  record <- new.env(parent = emptyenv())
  record$seen <- character(0)
  value <- withCallingHandlers(code, warning = function(w) {
    record$seen <- c(record$seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, seen = record$seen)
}

test_that("planning_sentence() names the corrections the plan made, and only those", {
  both <- planning_sentence(ss_buc_paired_t(t_observed = 3, N = 40))
  expect_match(both, paste0("after correcting the prior study's result for publication bias ",
                            "(assuming its literature published results significant at the ",
                            ".05 level) and for uncertainty."), fixed = TRUE)
  expect_no_match(both, "no publication bias", fixed = TRUE)
  expect_no_match(both, "not for uncertainty", fixed = TRUE)

  # alpha_prior = 1 assumes no publication bias, so none is corrected.
  no_bias <- planning_sentence(ss_buc_paired_t(t_observed = 3, N = 40, alpha_prior = 1))
  expect_match(no_bias, "after correcting the prior study's result for uncertainty, with no publication bias assumed.",
               fixed = TRUE)
  expect_no_match(no_bias, "for publication bias", fixed = TRUE)
  expect_no_match(no_bias, "significant at the", fixed = TRUE)

  # assurance = .5 takes the median, which corrects for bias only.
  no_uncertainty <- planning_sentence(ss_buc_paired_t(t_observed = 3, N = 40, assurance = .5))
  expect_match(no_uncertainty, "for publication bias (assuming its literature published results significant at the .05 level) but not for uncertainty.",
               fixed = TRUE)
  expect_no_match(no_uncertainty, "and for uncertainty", fixed = TRUE)

  neither <- planning_sentence(ss_buc_paired_t(t_observed = 3, N = 40, alpha_prior = 1,
                                               assurance = .5))
  expect_match(neither, "using the prior study's result with no correction for publication bias or for uncertainty.",
               fixed = TRUE)
  expect_no_match(neither, "after correcting", fixed = TRUE)

  # The threshold is the plan's own alpha_prior.
  expect_match(planning_sentence(ss_buc_paired_t(t_observed = 3, N = 40, alpha_prior = .025,
                                                 assurance = .6)),
               "significant at the .025 level) and for uncertainty.", fixed = TRUE)
})

test_that("planning_sentence() writes counts with thousands separators", {
  indt <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .9)
  expect_identical(size(indt), 1485)
  s <- planning_sentence(indt)
  expect_match(s, "A planned study with 1,485 participants per group (total N = 2,970) attains 80% power with 90% assurance",
               fixed = TRUE)
  pairs <- ss_buc_paired_t(t_observed = 3, N = 40, assurance = .88)
  expect_gt(size(pairs), 999)
  expect_match(planning_sentence(pairs),
               paste0("A planned study with ", formatC(size(pairs), big.mark = ",", format = "d"),
                      " pairs attains"), fixed = TRUE)
  total <- ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3,
                           effect = "factor_B")
  expect_match(planning_sentence(total), "with a total sample size of 1,902 attains",
               fixed = TRUE)
  for (s in list(s, planning_sentence(pairs), planning_sentence(total)))
    expect_no_match(s, "e+", fixed = TRUE)
})

test_that("a proportion below 1 never prints as 100%", {
  # It gets as many decimals as it needs to stay below 100.
  high_power <- planning_sentence(ss_buc_paired_t(t_observed = 3, N = 40, desired_power = .9996))
  expect_match(high_power, "attains 99.96% power with 80% assurance", fixed = TRUE)
  expect_no_match(high_power, "100%", fixed = TRUE)
  high_assurance <- planning_sentence(ss_buc_paired_t(t_observed = 10, N = 40, assurance = .9996))
  expect_match(high_assurance, "attains 80% power with 99.96% assurance", fixed = TRUE)
  expect_no_match(high_assurance, "100%", fixed = TRUE)
  higher <- planning_sentence(ss_buc_paired_t(t_observed = 10, N = 40, assurance = .9999))
  expect_match(higher, "attains 80% power with 99.99% assurance", fixed = TRUE)
  expect_no_match(higher, "100%", fixed = TRUE)
  # However close to 1 (QC 2026-09-30, finding 10: 1 - 1e-11 printed as 100%).
  near_one <- planning_sentence(ss_buc_independent_t(t_observed = 12, n = c(50, 55),
                                                     desired_power = 1 - 1e-11))
  expect_match(near_one, "attains 99.999999999% power", fixed = TRUE)
  expect_no_match(near_one, "100%", fixed = TRUE)
  # A proportion that needs no decimals still gets none.
  expect_match(planning_sentence(ss_buc_paired_t(t_observed = 3, N = 40, assurance = .825)),
               "attains 80% power with 82.5% assurance", fixed = TRUE)
})

test_that("a value below 1 never prints as 1 in the table, the kable, or the sentence", {
  # PC-6 of the 2026-10-12 triage. Whole numbers print clean, so a printed 1
  # reads as exactly 1; three significant figures showed this plan's power of
  # .99957 as "1". A value below 1 in magnitude now gets as many more figures
  # as it takes to stay below 1.
  fmt <- BUCSS:::.format_bucss_value
  expect_identical(fmt(c(.99957412839, .9996, .99999, 1 - 1e-12), 3),
                   c("0.9996", "0.9996", "0.99999", "0.999999999999"))
  expect_identical(fmt(c(-.99957, -.99999), 3), c("-0.9996", "-0.99999"))
  # The last double below 1 prints as stored, and a 'digits' above the
  # session's own print precision is covered too.
  expect_identical(fmt(1 - .Machine$double.eps / 2, 3), "0.99999999999999989")
  expect_identical(fmt(.999999996, 10), "0.999999996")
  for (x in c(.9995, .99949999, .99957, 1 - 1e-7, 1 - 1e-10, 1 - 1e-15, -.9999))
    for (d in c(1, 2, 3, 3.5, 6, 10, 15, 20))
      expect_lt(abs(as.numeric(fmt(x, d))), 1)
  # Every other value prints exactly as before: whole numbers clean, the rest
  # at 'digits' significant figures.
  expect_identical(fmt(c(1, 0, 2, 255, .9994, .5, .05, 12.3456, 123456.7, .0005, 1.0004, NA), 3),
                   c("1", "0", "2", "255", "0.999", "0.5", "0.05", "12.3", "123000", "0.0005",
                     "1", NA))

  plan <- ss_buc_welch_t(t_observed = 8, n_1 = 20, n_2 = 20, alpha_prior = 1, assurance = .5)
  expect_equal(val(plan, "actual_power"), 0.99957412839, tolerance = 1e-9)
  shown <- capture.output(print(plan))
  expect_match(grep("actual_power", shown, value = TRUE), "^ actual_power +0\\.9996 *$")
  if (requireNamespace("knitr", quietly = TRUE))
    expect_match(as.character(knitr::knit_print(plan)), "|actual_power          |0.9996 |",
                 fixed = TRUE)
  # The sentence's planning effect size too: a planning correlation short of
  # 1 is not written as 1.
  expect_identical(BUCSS:::.format_effect_size(.99996, "rho_planning"), ".99996")
  expect_identical(BUCSS:::.format_effect_size(.99996, "f_planning"), "0.99996")
  expect_identical(BUCSS:::.format_effect_size(.176343, "dz_planning"), "0.176")
  expect_identical(BUCSS:::.format_effect_size(.25, "R2_planning"), ".25")
})

test_that("the below-1 rule reads the text, so display options cannot turn it into an error", {
  # Round 2 of the 2026-10-12 fixes. The rule first decided whether a value
  # printed as 1 by parsing its own text back to a number. Under
  # options(OutDec = ",") that text is not a number R can parse, so print(),
  # knit_print(), and planning_sentence() stopped on every plan with "missing
  # value where TRUE/FALSE needed"; and a 'bucss.digits' option that is no
  # single positive number, which print() refuses by name but which reaches
  # knit_print() and planning_sentence() unvalidated, became that error or
  # "invalid value 0 for 'digits' argument" where those two had returned
  # output before the rule existed.
  fmt <- BUCSS:::.format_bucss_value
  sig <- BUCSS:::.format_signif
  plan <- ss_buc_correlation(r_observed = .99, N = 20)
  rho <- val(plan, "rho_planning")
  expect_equal(rho, 0.984296678, tolerance = 1e-6)
  second <- function(x) sub(".*The plan is built on a planning effect size of ", "",
                            planning_sentence(x))
  has_knitr <- requireNamespace("knitr", quietly = TRUE)
  kable_row <- function(x, term) {
    lines <- strsplit(as.character(knitr::knit_print(x)), "\n")[[1]]
    gsub(" +", "", grep(term, lines, value = TRUE, fixed = TRUE))
  }

  # A decimal comma: every row prints, and a value that would round to 1 is
  # still given the figures that keep it below 1.
  withr::with_options(list(OutDec = ","), {
    expect_identical(fmt(c(.801, .99957, -.99957, .5, 12.3456), 3),
                     c("0,801", "0,9996", "-0,9996", "0,5", "12,3"))
    shown <- capture.output(print(plan))
    expect_match(grep("rho_planning", shown, value = TRUE), "^ rho_planning +0,984 *$")
    expect_match(grep("actual_power", shown, value = TRUE), "^ actual_power +0,883 *$")
    # (format() warns that the sentence's thousands separator and the decimal
    # mark are both commas; that warning is base R's and is not the subject.)
    expect_identical(suppressWarnings(second(plan)), "0,984, expressed as a correlation.")
    if (has_knitr) expect_identical(kable_row(plan, "rho_planning"), "|rho_planning|0,984|")
  })

  # 'bucss.digits' below 1 is read as one significant figure, and the
  # planning correlation of .984, which one figure rounds to 1, is not
  # written as 1.
  for (bad in list(0, .4, -1)) withr::with_options(list(bucss.digits = bad), {
    expect_identical(second(plan), ".98, expressed as a correlation.", label = deparse(bad))
    if (has_knitr) {
      expect_identical(kable_row(plan, "rho_planning"), "|rho_planning|0.98|",
                       label = deparse(bad))
      expect_identical(kable_row(plan, "actual_power"), "|actual_power|0.9|",
                       label = deparse(bad))
    }
  })
  expect_identical(vapply(c(.98, -.98, .96, .5, .04), sig, character(1), digits = 0),
                   c("0.98", "-0.98", "0.96", "0.5", "0.04"))

  # 'bucss.digits' missing or of length 2: what format(signif()) itself
  # gives, as before the rule existed, never a condition from inside the rule.
  withr::with_options(list(bucss.digits = NA), {
    expect_identical(second(plan), "NA, expressed as a correlation.")
    if (has_knitr) expect_identical(kable_row(plan, "rho_planning"), "|rho_planning|NA|")
  })
  expect_identical(sig(.98, NA), "NA")
  expect_identical(sig(.98, NaN), "NaN")
  expect_identical(sig(rho, c(2, 3)), c("0.980", "0.984"))
  withr::with_options(list(bucss.digits = c(2, 3)),
                      expect_length(planning_sentence(plan), 2L))

  # print() still refuses each of these by name, and a valid option still
  # prints as it did.
  for (bad in list(0, NA, c(2, 3), -1)) withr::with_options(list(bucss.digits = bad),
    expect_error(capture.output(print(plan)), "'digits' must be a single positive number.",
                 fixed = TRUE, label = paste(deparse(bad), collapse = " ")))
  expect_identical(second(plan), ".984, expressed as a correlation.")
})

test_that("planning_sentence() writes an alpha of .1 as APA does, .10", {
  s <- planning_sentence(ss_buc_paired_t(t_observed = 3, N = 40, alpha_prior = .1))
  expect_match(s, "significant at the .10 level)", fixed = TRUE)
  expect_no_match(s, ".1 level", fixed = TRUE)
})

test_that("a stack of plans is not tidied, glanced, or narrated as one plan", {
  a <- ss_buc_paired_t(t_observed = 3, N = 40)
  b <- ss_buc_paired_t(t_observed = 4, N = 40)
  shape <- "no longer has the shape of a single plan"

  # What vctrs-based tools do: keep the class on the stacked rows.
  stack <- structure(rbind(as.data.frame(a), as.data.frame(b)),
                     class = c("bucss_plan", "data.frame"))
  expect_error(tidy(stack), shape, fixed = TRUE)
  expect_error(glance(stack), shape, fixed = TRUE)
  expect_error(planning_sentence(stack), shape, fixed = TRUE)
  shown <- capture.output(print(stack))
  expect_length(grep("necessary_N", shown, fixed = TRUE), 2L)
  expect_false(any(grepl("Design:", shown, fixed = TRUE)))
  expect_false(any(grepl("Sample size unit:", shown, fixed = TRUE)))
  expect_false(any(grepl("Largest supportable assurance", shown, fixed = TRUE)))
  capture.output(vis <- withVisible(print(stack)))
  expect_false(vis$visible)
  expect_identical(vis$value, stack)

  # A plan with a row removed is not one plan either.
  trimmed <- structure(as.data.frame(a)[-2, ], class = c("bucss_plan", "data.frame"))
  expect_error(tidy(trimmed), shape, fixed = TRUE)
  expect_error(planning_sentence(trimmed), shape, fixed = TRUE)

  # The single plan is untouched: its print keeps the footers.
  expect_true(any(grepl("Design: Dependent (paired) t test", capture.output(print(a)),
                        fixed = TRUE)))

  # The audit and the figure refuse a stack too, before replaying anything.
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_error(ss_buc_sensitivity(stack, true_ncp = 2, replications = 10, seed = 2017),
               shape, fixed = TRUE)
  expect_error(plot(stack), shape, fixed = TRUE)
})

# Copy a plan's attributes (all but the data frame's own) onto 'y'.
restore_plan_attributes <- function(y, x) {
  for (nm in setdiff(names(attributes(x)), c("names", "row.names")))
    attr(y, nm) <- attr(x, nm)
  y
}

test_that("a plan with an added column is not printed or tidied as one plan", {
  a <- ss_buc_paired_t(t_observed = 3, N = 40)
  # What dplyr::bind_cols() would leave: the rows intact, a third column, and
  # the class and attributes kept.
  wide <- restore_plan_attributes(cbind(as.data.frame(a), extra = 1), a)
  expect_s3_class(wide, "bucss_plan")
  expect_identical(attr(wide, "design"), attr(a, "design"))
  shown <- capture.output(print(wide))
  expect_true(any(grepl("extra", shown, fixed = TRUE)))
  expect_length(grep("necessary_N", shown, fixed = TRUE), 1L)
  expect_false(any(grepl("Design:", shown, fixed = TRUE)))
  expect_false(any(grepl("Sample size unit:", shown, fixed = TRUE)))
  expect_false(any(grepl("Largest supportable assurance", shown, fixed = TRUE)))
  expect_error(tidy(wide), "no longer has the shape of a single plan", fixed = TRUE)
  # The refusal names this cause as well: the object is neither a stack nor
  # missing rows (PC-8 of the 2026-10-12 triage; the message used to offer
  # only those two).
  cause <- paste0("it may be several plans stacked together, a plan with rows removed, or a ",
                  "plan with columns added or renamed.")
  expect_error(tidy(wide), cause, fixed = TRUE)
  expect_error(glance(wide), cause, fixed = TRUE)
  expect_error(planning_sentence(wide), cause, fixed = TRUE)
  # A renamed column, with the class kept (what dplyr::rename() would leave).
  renamed <- as.data.frame(a)
  names(renamed) <- c("term", "estimate")
  renamed <- restore_plan_attributes(renamed, a)
  expect_s3_class(renamed, "bucss_plan")
  expect_error(tidy(renamed), cause, fixed = TRUE)
})

test_that("a plan missing any row its planner echoed is not one plan", {
  # dplyr::filter() and the like keep the class and the attributes, so rows
  # are removed here the way they remove them. Without its 'alpha_prior' row
  # a plan corrected at .05 was narrated "with no publication bias assumed",
  # ss_buc_sensitivity() failed with "argument is of length zero", and with
  # 'alpha_prior = .10' the audit and the figure blamed the package (QC
  # 2026-09-30, finding 3). Every row of every echoed form is now required.
  drop_row <- function(x, term)
    restore_plan_attributes(as.data.frame(x)[x$term != term, ], x)
  shape <- "no longer has the shape of a single plan"
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  x <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .9,
                            alpha_prior = .10)
  trimmed <- drop_row(x, "alpha_prior")
  expect_s3_class(trimmed, "bucss_plan")
  expect_identical(attr(trimmed, "design"), "Independent t test")
  expect_error(planning_sentence(trimmed), shape, fixed = TRUE)
  expect_error(tidy(trimmed), shape, fixed = TRUE)
  expect_error(glance(trimmed), shape, fixed = TRUE)
  expect_error(ss_buc_sensitivity(trimmed, true_ncp = 2, replications = 10, seed = 2017),
               shape, fixed = TRUE)
  expect_error(plot(trimmed), shape, fixed = TRUE)
  shown <- capture.output(print(trimmed))
  expect_false(any(grepl("Design:", shown, fixed = TRUE)))

  # One plan in each echoed form of every design. Each plan as returned has
  # the shape; with any one row removed it does not.
  plans <- list(
    x,
    ss_buc_independent_t(t_observed = 3, n = 20),
    ss_buc_independent_t(t_observed = 3, N = 41),
    ss_buc_paired_t(t_observed = 3, N = 30),
    ss_buc_one_sample_t(t_observed = 3, N = 30),
    ss_buc_welch_t(t_observed = 3, n_1 = 30, n_2 = 40, sd_ratio = 1.5),
    ss_buc_welch_t(t_observed = 3, n_1 = 30, n_2 = 40, sd_1 = 1, sd_2 = 2),
    ss_buc_welch_t(t_observed = 3, n_1 = 30, n_2 = 40, var_1 = 1, var_2 = 3),
    ss_buc_one_way_anova(F_observed = 5, N = 60, levels_A = 3),
    ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                           effect = "factor_B"),
    ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6, df_numerator = 2,
                                   df_denominator = 114),
    ss_buc_ancova(F_observed = 5, N = 120, cells = 6, n_covariates = 2, df_numerator = 2),
    ss_buc_rm_anova(F_observed = 5, N = 30, levels_A = 3),
    ss_buc_rm_anova(F_observed = 10, N = 30, levels_A = 3, levels_B = 2, effect = "factor_B"),
    ss_buc_rm_anova(F_observed = 10, N = 30, levels_A = 3, levels_B = 2, effect = "interaction"),
    ss_buc_rm_anova_general(F_observed = 5, N = 60, df_numerator = 2),
    ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2, levels_within = 3,
                       effect = "within"),
    ss_buc_mixed_anova(F_observed = 8, N = 60, levels_between = 2, levels_within = 3,
                       effect = "interaction", epsilon = .8),
    ss_buc_mixed_anova(F_observed = 9, N = 60, levels_between = 2, levels_within = 3,
                       effect = "between"),
    ss_buc_mixed_anova_general(F_observed = 5, N = 60, df_numerator = 2, num_groups = 3,
                               effect = "between_only"),
    ss_buc_mixed_anova_general(F_observed = 8, N = 90, df_numerator = 2, num_groups = 3,
                               effect = "within_only"),
    ss_buc_mixed_anova_general(F_observed = 5, N = 60, df_numerator = 4, num_groups = 3,
                               effect = "between_within", df_num_within = 2),
    ss_buc_reg_coef(t_observed = 3, N = 100, p = 3),
    ss_buc_R2(F_observed = 5, N = 100, p = 3),
    ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2),
    ss_buc_correlation(r_observed = .35, N = 100),
    ss_buc_correlation(t_observed = 3.7, N = 100),
    ss_buc_manova(F_observed = 5, N = 60, p_variables = 3),
    ss_buc_manova(T2_observed = 16, N = 60, p_variables = 3),
    ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42, df_restricted = 43))
  covered <- unique(vapply(plans, attr, character(1), which = "design"))
  expect_setequal(covered, names(BUCSS:::.BUCSS_DESIGN_SPECS))
  for (p in plans) {
    expect_true(BUCSS:::.is_plan_shape(p), info = attr(p, "design"))
    for (term in p$term)
      expect_false(BUCSS:::.is_plan_shape(drop_row(p, term)),
                   info = paste(attr(p, "design"), "without", term))
  }

  # The one row that may be absent: levels_B of a two-way within-subjects
  # plan for Factor A, which is the one-way plan, so the replay still
  # reproduces it (and a one-way plan has no levels_B row to begin with).
  factor_a <- ss_buc_rm_anova(F_observed = 5, N = 30, levels_A = 3, levels_B = 2)
  expect_true(BUCSS:::.is_plan_shape(drop_row(factor_a, "levels_B")))
})

test_that("every registry entry states its echoed rows in the shape check's form", {
  # 'echo' is the design's requirements and 'echo_effect' a list named by
  # effects its planner accepts, each holding that effect's further
  # requirements, so a planning input added to a design is required by the
  # shape check once it is listed there (QC 2026-09-30, finding 3). A name
  # that is not one of the planner's effects would require nothing, silently.
  specs <- BUCSS:::.BUCSS_DESIGN_SPECS
  for (nm in names(specs)) {
    spec <- specs[[nm]]
    expect_type(spec$echo, "character")
    expect_gte(length(spec$echo), 2L)
    if (is.null(spec$echo_effect)) next
    expect_type(spec$echo_effect, "list")
    effects <- eval(formals(get(spec$planner, envir = asNamespace("BUCSS")))$effect)
    expect_true(length(effects) > 1L && all(names(spec$echo_effect) %in% effects),
                info = nm)
    expect_true(all(vapply(spec$echo_effect, is.character, logical(1))), info = nm)
    expect_identical(BUCSS:::.bucss_echo_requirements(spec, names(spec$echo_effect)[1]),
                     c(spec$echo, spec$echo_effect[[1]]))
  }
  rm <- specs[["One- or two-way within-subjects ANOVA"]]
  expect_identical(BUCSS:::.bucss_echo_requirements(rm, "factor_A"), rm$echo)
  expect_identical(BUCSS:::.bucss_echo_requirements(rm, NULL), rm$echo)
})

test_that("a plan whose rows were reordered still audits and plots as itself", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  shuffled <- restore_plan_attributes(as.data.frame(plan)[rev(seq_len(nrow(plan))), ],
                                      plan)
  expect_s3_class(shuffled, "bucss_plan")
  expect_false(identical(shuffled$term, plan$term))
  expect_identical(
    ss_buc_sensitivity(shuffled, true_ncp = 2, replications = 50, seed = 2017),
    ss_buc_sensitivity(plan, true_ncp = 2, replications = 50, seed = 2017))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_identical(plot(shuffled), plot(plan))
})

test_that("two stacked sensitivity audits print as a plain data frame", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  a <- ss_buc_sensitivity(plan, true_ncp = 2, replications = 50, seed = 2017)
  b <- ss_buc_sensitivity(plan, true_ncp = 3, replications = 50, seed = 2017)
  both <- rbind(a, b)
  expect_s3_class(both, "bucss_sensitivity")          # rbind() keeps the class
  expect_identical(attr(both, "design"), attr(a, "design"))
  shown <- capture.output(vis <- withVisible(print(both)))
  expect_length(grep("true_ncp", shown, fixed = TRUE), 2L)
  for (footer in c("Design:", "Sample size unit:", "Requested assurance:", "Target power:"))
    expect_false(any(grepl(footer, shown, fixed = TRUE)), info = footer)
  expect_false(vis$visible)
  expect_identical(vis$value, both)
  # One audit keeps its footers.
  expect_true(any(grepl("Requested assurance: .80", capture.output(print(a)), fixed = TRUE)))
})

test_that("knit_print() renders a stack as a plain data frame, without footers", {
  skip_if_not_installed("knitr")
  a <- ss_buc_paired_t(t_observed = 3, N = 40)
  b <- ss_buc_paired_t(t_observed = 4, N = 40)
  stack <- structure(rbind(as.data.frame(a), as.data.frame(b)),
                     class = c("bucss_plan", "data.frame"))
  # knitr's default method prints the plain data frame, as a chunk would.
  printed <- capture.output(out <- knitr::knit_print(stack))
  expect_false(inherits(out, "knit_asis"))
  expect_length(grep("necessary_N", printed, fixed = TRUE), 2L)
  expect_false(any(grepl("Largest supportable assurance", printed, fixed = TRUE)))
  expect_match(as.character(knitr::knit_print(a)), "Largest supportable assurance: .90",
               fixed = TRUE)
})

test_that("renaming a plan's columns drops the class", {
  a <- ss_buc_paired_t(t_observed = 3, N = 40)
  x <- a
  names(x) <- c("term", "value")
  expect_identical(class(x), "data.frame")
  expect_identical(class(setNames(a, c("term", "value"))), "data.frame")
  y <- a
  colnames(y) <- c("what", "number")
  expect_identical(class(y), "data.frame")
  expect_identical(y$number, a$value)
  # The renamed object prints its rows, not an empty table under footers.
  shown <- capture.output(print(y))
  expect_true(any(grepl("necessary_N", shown, fixed = TRUE)))
  expect_false(any(grepl("<0 rows>", shown, fixed = TRUE)))
})

test_that("glance() begins with tidy()'s three columns, then keeps every row", {
  plans <- list(
    ss_buc_paired_t(t_observed = 3, N = 40),
    ss_buc_independent_t(t_observed = 3, n = c(50, 55)),
    ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                           effect = "factor_B"),
    ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42, df_restricted = 43))
  for (p in plans) {
    g <- glance(p)
    td <- tidy(p)
    expect_identical(names(td), c("term", "estimate", "power"))
    expect_identical(g[1:3], td)
    expect_identical(nrow(g), 1L)
    expect_identical(unlist(g[p$term], use.names = FALSE), p$value)
    expect_identical(g$design, attr(p, "design"))
    expect_identical(g$assurance_ceiling, attr(p, "assurance_ceiling"))
  }
})

test_that("ss_buc_sensitivity() audits a MANOVA plan given as T squared", {
  # It used to stop, because the replay could not rebuild the call. The T
  # squared plan replays through its exact F, so the audit is the F plan's.
  t2 <- ss_buc_manova(T2_observed = 4.5 * 3 * 78 / 76, N = 80, p_variables = 3)
  f <- ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3)
  s_t2 <- ss_buc_sensitivity(t2, true_ncp = 1.4, replications = 50, seed = 2017)
  s_f <- ss_buc_sensitivity(f, true_ncp = 1.4, replications = 50, seed = 2017)
  expect_s3_class(s_t2, "bucss_sensitivity")
  expect_identical(s_t2$term, s_f$term)
  expect_identical(s_t2$value, s_f$value)
  expect_false(anyNA(s_t2$value))
})

test_that("the sensitivity rates and standard errors are the documented ones", {
  s <- ss_buc_sensitivity(ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 4),
                          true_ncp = 3.737, replications = 400, seed = 2017)
  R <- val(s, "replications")
  refusal <- val(s, "refusal_rate")
  att <- val(s, "attainment")
  att_r <- val(s, "attainment_with_refusals")
  expect_identical(R, 400)
  expect_gt(refusal, 0)                  # the identities below need both kinds
  expect_lt(refusal, 1)
  # A refusal counts as attaining in attainment_with_refusals.
  expect_equal(att_r, refusal + (1 - refusal) * att)
  # Each rate's Monte Carlo standard error uses its own denominator: the
  # issued plans for attainment, all replications for the other.
  issued <- R * (1 - refusal)
  expect_equal(issued, round(issued))
  expect_equal(val(s, "attainment_mcse"), sqrt(att * (1 - att) / issued))
  expect_equal(val(s, "attainment_with_refusals_mcse"), sqrt(att_r * (1 - att_r) / R))
  expect_false(isTRUE(all.equal(val(s, "attainment_mcse"), sqrt(att * (1 - att) / R))))
})

test_that("a plan of exactly the size the truth requires counts as attaining", {
  # With the Welch floor of 10 per group, a strong truth requires 10 per group
  # and no simulated plan can be smaller, so every issued plan attains: the
  # comparison is 'at least', not 'more than' (mutant C05).
  plan <- ss_buc_welch_t(t_observed = 8, n_1 = 20, n_2 = 20)
  s <- ss_buc_sensitivity(plan, true_ncp = 6.8, replications = 50, seed = 2017)
  expect_identical(val(s, "size_at_true_effect"), 10)
  expect_identical(val(s, "size_q25"), 10)
  expect_lt(val(s, "refusal_rate"), 1)
  expect_identical(val(s, "attainment"), 1)
  expect_identical(val(s, "attainment_with_refusals"), 1)
})

test_that("the planner's warnings reach the user once, summarized and counted", {
  # A planner that warns on every call: every simulated prior study warns once
  # (a refused one too, when the warning comes before the refusal), and so
  # does every call of the search for the sample size the true effect
  # requires. (This test used the planners' own warning for an assurance below
  # .5 as its source; since 2026-10-12 that one is left out of the summary,
  # see the next test but one, so the warning is injected.)
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  planner <- ss_buc_paired_t
  local_mocked_bindings(ss_buc_paired_t = function(...) {
    warning("a planner warning")
    planner(...)
  }, .package = "BUCSS")
  reps <- 20
  seen <- collect_warnings(
    ss_buc_sensitivity(plan, true_ncp = 1.5, replications = reps, seed = 2017))$seen
  # The replay of the plan itself, which checks the stored inputs, is silent,
  # so the summary is the only warning.
  expect_length(seen, 1L)
  pattern <- paste0("^The planner raised warnings on ([0-9]+) of the ([0-9]+) simulated ",
                    "prior studies and ([0-9]+) times while finding the sample size the ",
                    "true effect requires; the first was: a planner warning$")
  expect_match(seen, pattern)
  counts <- as.numeric(regmatches(seen, regexec(pattern, seen))[[1]][-1])
  expect_identical(counts[1], reps)                     # every simulated study
  expect_identical(counts[2], reps)
  expect_gt(counts[3], 0)                               # and the search, counted apart
})

test_that("a planner run that warns twice is counted once", {
  # Every run warns twice here, both warnings injected into the planner, and
  # the summary counts runs, not warnings (QC 2026-09-30, finding 9: it
  # reported warnings on 40 of the 20 studies).
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  planner <- ss_buc_paired_t
  local_mocked_bindings(ss_buc_paired_t = function(...) {
    warning("a first warning")
    out <- planner(...)
    warning("a second warning")
    out
  }, .package = "BUCSS")
  reps <- 20
  seen <- collect_warnings(
    ss_buc_sensitivity(plan, true_ncp = 1.5, replications = reps, seed = 2017))$seen
  expect_length(seen, 1L)
  pattern <- paste0("^The planner raised warnings on ([0-9]+) of the ([0-9]+) simulated ",
                    "prior studies and ([0-9]+) times while finding")
  counts <- as.numeric(regmatches(seen, regexec(pattern, seen))[[1]][-1])
  expect_identical(counts[1], reps)
  expect_identical(counts[2], reps)
  expect_match(seen, "the first was: a first warning", fixed = TRUE)
})

test_that("the plan's own assurance-below-.5 warning is not re-reported by the audit", {
  # ENG-2 of the 2026-10-12 triage. Every simulated study is planned at the
  # plan's own assurance, so a plan made at an assurance below .5 raised the
  # planners' warning about it on every run, and the audit reported "warnings
  # on 50 of the 50 simulated prior studies ...; the first was: The
  # 'assurance' you have entered is < .5", the warning the user had already
  # seen when the plan was built. It is now muffled without being counted, as
  # plot() muffles it.
  expect_warning(plan <- ss_buc_paired_t(t_observed = 3, N = 30, assurance = .4),
                 "The 'assurance' you have entered is < .5", fixed = TRUE)
  audit <- collect_warnings(
    ss_buc_sensitivity(plan, true_ncp = 2, replications = 50, seed = 2017))
  expect_length(audit$seen, 0L)
  expect_s3_class(audit$value, "bucss_sensitivity")
  expect_gt(audit$value$value[audit$value$term == "refusal_rate"], 0)

  # Any other warning is still summarized, counted by run, and named: here
  # the planner raises its own assurance warning first and an injected one
  # after it, and the summary's first warning is the injected one.
  planner <- ss_buc_paired_t
  local_mocked_bindings(ss_buc_paired_t = function(...) {
    out <- planner(...)
    warning("a later warning")
    out
  }, .package = "BUCSS")
  reps <- 20
  seen <- collect_warnings(
    ss_buc_sensitivity(plan, true_ncp = 3, replications = reps, seed = 2017))$seen
  expect_length(seen, 1L)
  expect_match(seen, "the first was: a later warning", fixed = TRUE)
  expect_no_match(seen, "you have entered is < .5", fixed = TRUE)
  pattern <- "^The planner raised warnings on ([0-9]+) of the ([0-9]+) simulated prior studies"
  counts <- as.numeric(regmatches(seen, regexec(pattern, seen))[[1]][-1])
  expect_lte(counts[1], reps)
  expect_gt(counts[1], 0)
  expect_identical(counts[2], reps)
})

test_that("an audit of a plan whose planner never warns raises no summary", {
  plan <- ss_buc_paired_t(t_observed = 3, N = 30)
  seen <- collect_warnings(
    ss_buc_sensitivity(plan, true_ncp = 2, replications = 20, seed = 2017))$seen
  expect_length(seen, 0L)
})

test_that("plot() draws a plan whose ceiling is below .5", {
  # The planner's own warning, at building the plan:
  expect_warning(plan <- ss_buc_paired_t(t_observed = 2.3, N = 30, assurance = .3),
                 "< .5", fixed = TRUE)
  ceiling_value <- attr(plan, "assurance_ceiling")
  expect_lt(ceiling_value, .5)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  # Plotting it draws: the grid lies below .5 and every point of it is a
  # plan, so there is no "nothing was drawn" warning, and the grid's own
  # below-.5 warnings are muffled, as is the replay of the plan itself that
  # checks the stored inputs. So plotting raises no warning at all.
  drawn <- collect_warnings(plot(plan))
  curve_values <- drawn$value
  expect_length(drawn$seen, 0L)
  expect_identical(nrow(curve_values), 60L)
  expect_true(all(curve_values$assurance < ceiling_value))
  expect_false(anyNA(curve_values$ncp_adjusted))
  expect_false(anyNA(curve_values$sample_size))
  expect_true(all(curve_values$ncp_adjusted > 0))
})
