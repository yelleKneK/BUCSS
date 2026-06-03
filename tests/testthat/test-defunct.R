# Each dotted-name function from BUCSS 1.x is defunct as of 2.0.0: calling it
# must raise an error that says so and names its ss_buc_* replacement.

defunct_map <- list(
  ss.power.it          = "ss_buc_independent_t",
  ss.power.dt          = "ss_buc_paired_t",
  ss.power.ba.general  = "ss_buc_factorial_anova_general",
  ss.power.wa          = "ss_buc_rm_anova",
  ss.power.wa.general  = "ss_buc_rm_anova_general",
  ss.power.spa         = "ss_buc_mixed_anova",
  ss.power.spa.general = "ss_buc_mixed_anova_general",
  ss.power.reg1        = "ss_buc_reg_coef",
  ss.power.reg.all     = "ss_buc_R2",
  ss.power.reg.joint   = "ss_buc_reg_joint"
)

test_that("each dotted 1.x function is defunct and names its replacement", {
  for (old in names(defunct_map)) {
    fn <- get(old)
    expect_error(fn(), "defunct", fixed = TRUE)
    expect_error(fn(), defunct_map[[old]], fixed = TRUE)
  }
})

test_that("ss.power.ba is defunct and points to both between-subjects replacements", {
  expect_error(ss.power.ba(), "defunct", fixed = TRUE)
  expect_error(ss.power.ba(), "ss_buc_one_way_anova", fixed = TRUE)
  expect_error(ss.power.ba(), "ss_buc_factorial_anova", fixed = TRUE)
})

test_that("the defunct error points to the CRAN archive for the old 1.x behavior", {
  expect_error(ss.power.it(),
               "cran.r-project.org/src/contrib/Archive/BUCSS", fixed = TRUE)
})
