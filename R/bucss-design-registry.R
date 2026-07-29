# The design registry: the one place that records, per design, how to rebuild
# the planner call that produced a result and what distribution the prior
# study's statistic follows. Two exported features read it, ss_buc_sensitivity()
# and plot.bucss_power(), so it lives here rather than inside either.
#
# The degrees of freedom recorded here are the prior study's OWN, which is what
# a simulated prior study must be drawn from. They can differ from the degrees
# of freedom a planner uses internally: the between-subjects and split-plot
# planners round an unbalanced prior N to a balanced design in two directions
# and correct under each, which is a deliberate conservatism of the planner
# rather than a fact about the prior study.
#
# Every entry is pinned in tests/testthat/test-sensitivity.R against an
# independently written computation of the prior study's p value, and both
# features additionally check at run time that rebuilding the planner call from
# the result object reproduces that object.

.BUCSS_DESIGN_SPECS <- list(
  "Independent t test" = list(
    planner = "ss_buc_independent_t", statistic = "t_observed", family = "t",
    df = function(i) {
      if (!is.null(i$n_1)) i$n_1 + i$n_2 - 2
      else if (!is.null(i$n)) 2 * i$n - 2
      else i$N - 2
    }),
  "Dependent (paired) t test" = list(
    planner = "ss_buc_paired_t", statistic = "t_observed", family = "t",
    df = function(i) i$N - 1),
  "One-sample t test" = list(
    planner = "ss_buc_one_sample_t", statistic = "t_observed", family = "t",
    df = function(i) i$N - 1),
  "Welch (unequal variance) t test" = list(
    planner = "ss_buc_welch_t", statistic = "t_observed", family = "t",
    df = function(i) {
      v1 <- 1 / i$n_1
      v2 <- i$sd_ratio^2 / i$n_2
      (v1 + v2)^2 / (v1^2 / (i$n_1 - 1) + v2^2 / (i$n_2 - 1))
    }),
  "One-way between-subjects ANOVA" = list(
    planner = "ss_buc_one_way_anova", statistic = "F_observed", family = "f",
    df = function(i) c(i$levels_A - 1, i$N - i$levels_A)),
  "Two-way between-subjects ANOVA" = list(
    planner = "ss_buc_factorial_anova", statistic = "F_observed", family = "f",
    df = function(i, effect) {
      df1 <- switch(effect, factor_A = i$levels_A - 1, factor_B = i$levels_B - 1,
                    interaction = (i$levels_A - 1) * (i$levels_B - 1))
      c(df1, i$N - i$levels_A * i$levels_B)
    }),
  "Between-subjects ANOVA (general)" = list(
    planner = "ss_buc_factorial_anova_general", statistic = "F_observed",
    family = "f", df = function(i) c(i$df_numerator, i$df_denominator)),
  "Analysis of covariance" = list(
    planner = "ss_buc_ancova", statistic = "F_observed", family = "f",
    df = function(i) c(i$df_numerator, i$N - i$cells - i$n_covariates)),
  "One or two-way within-subjects ANOVA" = list(
    planner = "ss_buc_rm_anova", statistic = "F_observed", family = "f",
    df = function(i, effect) {
      df1 <- if (is.null(i$levels_B)) i$levels_A - 1 else
        switch(effect, factor_A = i$levels_A - 1, factor_B = i$levels_B - 1,
               interaction = (i$levels_A - 1) * (i$levels_B - 1))
      c(df1, df1 * (i$N - 1))
    }),
  "Within-subjects ANOVA (any number of factors)" = list(
    planner = "ss_buc_rm_anova_general", statistic = "F_observed", family = "f",
    df = function(i) c(i$df_numerator, i$df_numerator * (i$N - 1))),
  "Two-factor split-plot (mixed) ANOVA" = list(
    planner = "ss_buc_mixed_anova", statistic = "F_observed", family = "f",
    df = function(i, effect) {
      df1 <- switch(effect, between = i$levels_between - 1,
                    within = i$levels_within - 1,
                    interaction = (i$levels_between - 1) * (i$levels_within - 1))
      df2 <- if (effect == "between") i$N - i$levels_between else
        (i$N - i$levels_between) * (i$levels_within - 1)
      c(df1, df2)
    }),
  "Split-plot (mixed) ANOVA (any number of factors)" = list(
    planner = "ss_buc_mixed_anova_general", statistic = "F_observed",
    family = "f",
    df = function(i, effect) {
      mult <- switch(effect, between_only = 1, within_only = i$df_numerator,
                     between_within = i$df_num_within)
      c(i$df_numerator, (i$N - i$num_groups) * mult)
    }),
  "Multiple regression: single coefficient" = list(
    planner = "ss_buc_reg_coef", statistic = "t_observed", family = "f_from_t",
    df = function(i) c(1, i$N - i$p - 1)),
  "Multiple regression: omnibus test of model R2" = list(
    planner = "ss_buc_R2", statistic = "F_observed", family = "f",
    df = function(i) c(i$p, i$N - i$p - 1)),
  "Multiple regression: joint test of predictors" = list(
    planner = "ss_buc_reg_joint", statistic = "F_observed", family = "f",
    df = function(i) c(i$p_joint, i$N - i$p - 1)),
  "Two-group multivariate comparison (Hotelling's T squared)" = list(
    planner = "ss_buc_manova", statistic = "F_observed", family = "f",
    df = function(i) c(i$p_variables, i$N - i$p_variables - 1)),
  "Nested model chi-square difference test" = list(
    planner = "ss_buc_chisq_diff", statistic = "chisq_observed",
    family = "chisq", df = function(i) i$df_difference),
  "Pearson correlation" = list(
    planner = "ss_buc_correlation", statistic = "t_observed",
    family = "correlation", df = function(i) i$N)
)

# Pull the planning inputs back off a result as a named list, one element per
# echoed row. The design table above reads these row names, so they stay
# exactly as the constructor wrote them; rebuilding the planner's own argument
# names is .bucss_design_args()'s job.
.bucss_design_inputs <- function(object) {
  drop <- c(.BUCSS_SIZE_TERMS, "total_N", "actual_power", "ncp_adjusted")
  keep <- !(object$term %in% drop)
  out <- as.list(object$value[keep])
  names(out) <- object$term[keep]
  out
}

# Turn those rows into the planner's arguments. The two places the row names
# and the formals differ: a length-2 prior 'n' is echoed as n_1/n_2 by the
# independent t (the Welch planner's formals really are n_1 and n_2), and the
# effect travels as an attribute rather than a row.
.bucss_design_args <- function(inputs, spec, effect) {
  args <- inputs
  if (spec$planner == "ss_buc_independent_t" && !is.null(args$n_1)) {
    args$n <- c(args$n_1, args$n_2)
    args$n_1 <- NULL
    args$n_2 <- NULL
  }
  if (!is.null(effect)) args$effect <- effect
  # The correlation planner accepts either the correlation or its t; simulated
  # statistics arrive as a t, so a call stated in r is restated in t.
  if (spec$family == "correlation" && !is.null(args$r_observed)) {
    r <- abs(args$r_observed)
    args$t_observed <- r * sqrt((args$N - 2) / (1 - r^2))
    args$r_observed <- NULL
  }
  args
}

