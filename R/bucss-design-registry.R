# The design registry: the one place that records, per design, how to rebuild
# the planner call that produced a result and what distribution the prior
# study's statistic follows. Two exported features read it, ss_buc_sensitivity()
# and plot.bucss_plan(), so it lives here rather than inside either.
#
# The degrees of freedom recorded here are the prior study's own, which is what
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
#
# 'echo' lists the design's own echoed planning inputs, the rows a replay
# rebuilds the call from, beyond the four every planner echoes (alpha_prior,
# alpha_planned, assurance, desired_power). Each element is one requirement;
# " | " separates the forms a planner may have echoed (the independent t's
# prior size is n, or n_1 and n_2, or N) and " + " joins the rows of one form.
# 'echo_effect' is a list named by the effects that require more: each element
# holds the further requirements, in the same grammar, of a plan for that
# effect (the attribute 'effect' of the plan). .is_plan_shape() reads both
# through .bucss_echo_requirements(), so a plan with one of these rows removed
# is refused with the shape message before a replay can fail on it or
# silently fall back to the planner's default. Every planning input a planner
# echoes belongs here, for every plan that echoes it, unless a replay without
# it reproduces the plan (levels_B of a two-way within-subjects Factor A plan,
# below); tests/testthat/test-output-fixes.R removes each row of one plan per
# echoed form and checks that the result is refused.
#
# 'ncp_scale', where present, maps the echoed inputs to the factor by which
# the plan's 'ncp_adjusted' is multiplied to give the noncentrality parameter
# of the prior study's sampling distribution; absent, it is 1.
#
# 'planning_term' names the planning effect size row the design's planner
# writes (one of .BUCSS_EFFECT_TERMS in R/bucss-internal.R): the effect size
# the plan is built on, ncp_adjusted restated in the design's usual metric.
# .is_plan_shape() requires exactly that row, so a plan with it removed or
# renamed is refused. It is a computed row, like the sample size, so the
# replays never pass it back to the planner; they compare it.

# The sphericity parameter of a within-subjects plan: its echoed 'epsilon'
# row, or 1 for the between-subjects effect of a split-plot design, which
# involves no sphericity and echoes no such row. Every within-subjects and
# interaction plan echoes the row (1 when it assumed sphericity), and it is
# required in 'echo' or 'echo_effect' below, so a plan with the row removed
# fails the shape check rather than reaching this default. Under the Muller
# and Barton (1989) approximation a within-subjects statistic is noncentral F
# on epsilon times its degrees of freedom with noncentrality epsilon times
# the sphericity-scale parameter, so the entries below scale their df by it
# and name it as their 'ncp_scale', and a simulated prior study follows the
# same approximation the planner corrects under.
.bucss_epsilon <- function(i) if (is.null(i$epsilon)) 1 else i$epsilon

.BUCSS_DESIGN_SPECS <- list(
  "Independent t test" = list(
    planner = "ss_buc_independent_t", statistic = "t_observed", family = "t",
    planning_term = "smd_planning",
    echo = c("t_observed", "n | n_1 + n_2 | N"),
    df = function(i) {
      if (!is.null(i$n_1)) i$n_1 + i$n_2 - 2
      else if (!is.null(i$n)) 2 * i$n - 2
      else i$N - 2
    }),
  "Dependent (paired) t test" = list(
    planner = "ss_buc_paired_t", statistic = "t_observed", family = "t",
    planning_term = "dz_planning",
    echo = c("t_observed", "N"),
    df = function(i) i$N - 1),
  "One-sample t test" = list(
    planner = "ss_buc_one_sample_t", statistic = "t_observed", family = "t",
    planning_term = "std_mean_planning",
    echo = c("t_observed", "N"),
    df = function(i) i$N - 1),
  "Welch (unequal variance) t test" = list(
    planner = "ss_buc_welch_t", statistic = "t_observed", family = "t",
    planning_term = "smd_planning",
    echo = c("t_observed", "n_1", "n_2",
             "sd_ratio | sd_1 + sd_2 | var_1 + var_2"),
    df = function(i) {
      # The spread may be echoed as a ratio, as two standard deviations, or as
      # two variances; only the ratio enters the degrees of freedom.
      ratio <- if (!is.null(i$sd_ratio)) i$sd_ratio
               else if (!is.null(i$sd_1)) i$sd_2 / i$sd_1
               else sqrt(i$var_2 / i$var_1)
      v1 <- 1 / i$n_1
      v2 <- ratio^2 / i$n_2
      (v1 + v2)^2 / (v1^2 / (i$n_1 - 1) + v2^2 / (i$n_2 - 1))
    }),
  "One-way between-subjects ANOVA" = list(
    planner = "ss_buc_one_way_anova", statistic = "F_observed", family = "f",
    planning_term = "f_planning",
    echo = c("F_observed", "N", "levels_A"),
    df = function(i) c(i$levels_A - 1, i$N - i$levels_A)),
  "Two-way between-subjects ANOVA" = list(
    planner = "ss_buc_factorial_anova", statistic = "F_observed", family = "f",
    planning_term = "f_planning",
    echo = c("F_observed", "N", "levels_A", "levels_B"),
    df = function(i, effect) {
      df1 <- switch(effect, factor_A = i$levels_A - 1, factor_B = i$levels_B - 1,
                    interaction = (i$levels_A - 1) * (i$levels_B - 1))
      c(df1, i$N - i$levels_A * i$levels_B)
    }),
  "Between-subjects ANOVA (general)" = list(
    planner = "ss_buc_factorial_anova_general", statistic = "F_observed",
    family = "f",
    planning_term = "f_planning",
    echo = c("F_observed", "N", "cells", "df_numerator", "df_denominator"),
    df = function(i) c(i$df_numerator, i$df_denominator)),
  "Analysis of covariance" = list(
    planner = "ss_buc_ancova", statistic = "F_observed", family = "f",
    planning_term = "f_planning",
    echo = c("F_observed", "N", "cells", "n_covariates", "df_numerator"),
    df = function(i) c(i$df_numerator, i$N - i$cells - i$n_covariates)),
  "One- or two-way within-subjects ANOVA" = list(
    planner = "ss_buc_rm_anova", statistic = "F_observed", family = "f",
    planning_term = "f_planning",
    # levels_B is echoed only for a two-way design. It is required for the
    # effects that exist only there; a two-way Factor A plan is the one-way
    # plan, so without it the replay still reproduces the plan.
    echo = c("F_observed", "N", "levels_A", "epsilon"),
    echo_effect = list(factor_B = "levels_B", interaction = "levels_B"),
    df = function(i, effect) {
      df1 <- if (is.null(i$levels_B)) i$levels_A - 1 else
        switch(effect, factor_A = i$levels_A - 1, factor_B = i$levels_B - 1,
               interaction = (i$levels_A - 1) * (i$levels_B - 1))
      .bucss_epsilon(i) * c(df1, df1 * (i$N - 1))
    },
    ncp_scale = .bucss_epsilon),
  "Within-subjects ANOVA (any number of factors)" = list(
    planner = "ss_buc_rm_anova_general", statistic = "F_observed", family = "f",
    planning_term = "f_planning",
    echo = c("F_observed", "N", "df_numerator", "epsilon"),
    df = function(i) .bucss_epsilon(i) * c(i$df_numerator, i$df_numerator * (i$N - 1)),
    ncp_scale = .bucss_epsilon),
  # The between-subjects effect involves no sphericity and echoes no
  # 'epsilon', so only the within-subjects and interaction effects require it.
  "Two-factor split-plot (mixed) ANOVA" = list(
    planner = "ss_buc_mixed_anova", statistic = "F_observed", family = "f",
    planning_term = "f_planning",
    echo = c("F_observed", "N", "levels_between", "levels_within"),
    echo_effect = list(within = "epsilon", interaction = "epsilon"),
    df = function(i, effect) {
      df1 <- switch(effect, between = i$levels_between - 1,
                    within = i$levels_within - 1,
                    interaction = (i$levels_between - 1) * (i$levels_within - 1))
      df2 <- if (effect == "between") i$N - i$levels_between else
        (i$N - i$levels_between) * (i$levels_within - 1)
      .bucss_epsilon(i) * c(df1, df2)
    },
    ncp_scale = .bucss_epsilon),
  "Split-plot (mixed) ANOVA (any number of factors)" = list(
    planner = "ss_buc_mixed_anova_general", statistic = "F_observed",
    family = "f",
    planning_term = "f_planning",
    echo = c("F_observed", "N", "df_numerator", "num_groups"),
    echo_effect = list(within_only = "epsilon",
                       between_within = c("df_num_within", "epsilon")),
    df = function(i, effect) {
      mult <- switch(effect, between_only = 1, within_only = i$df_numerator,
                     between_within = i$df_num_within)
      .bucss_epsilon(i) * c(i$df_numerator, (i$N - i$num_groups) * mult)
    },
    ncp_scale = .bucss_epsilon),
  # The regression planners take a random_predictors switch, echoed as a 1/0
  # row (required in 'echo', so a plan without it is refused rather than
  # simulated in a frame it was not computed in), and the sampling family
  # follows it; see .bucss_resolve_spec(). The df carry, after the numerator
  # and residual df, m = N - 1 - (p - q), the chi square df of the tested
  # predictors' spread once the other predictors are accounted for, and N,
  # which scales the noncentrality parameter to the partial f squared; the
  # fixed families read only the first two.
  "Multiple regression: single coefficient" = list(
    planner = "ss_buc_reg_coef", statistic = "t_observed",
    planning_term = "f2_planning",
    family = c(fixed = "f_from_t", random = "f_from_t_random"),
    echo = c("t_observed", "N", "p", "random_predictors"),
    df = function(i) c(1, i$N - i$p - 1, i$N - i$p, i$N)),
  "Multiple regression: omnibus test of model R2" = list(
    planner = "ss_buc_R2", statistic = "F_observed",
    planning_term = "R2_planning",
    family = c(fixed = "f", random = "f_random"),
    echo = c("F_observed", "N", "p", "random_predictors"),
    df = function(i) c(i$p, i$N - i$p - 1, i$N - 1, i$N)),
  "Multiple regression: joint test of predictors" = list(
    planner = "ss_buc_reg_joint", statistic = "F_observed",
    planning_term = "f2_planning",
    family = c(fixed = "f", random = "f_random"),
    echo = c("F_observed", "N", "p", "p_joint", "random_predictors"),
    df = function(i) c(i$p_joint, i$N - i$p - 1,
                       i$N - 1 - (i$p - i$p_joint), i$N)),
  "Two-group multivariate comparison (Hotelling's T squared)" = list(
    planner = "ss_buc_manova", statistic = "F_observed", family = "f",
    planning_term = "D2_planning",
    echo = c("F_observed | T2_observed", "N", "p_variables"),
    df = function(i) c(i$p_variables, i$N - i$p_variables - 1)),
  "Nested model chi square difference test" = list(
    planner = "ss_buc_chisq_diff", statistic = "chisq_observed",
    family = "chisq",
    planning_term = "rmsea_d_planning",
    echo = c("chisq_observed", "N", "df_full", "df_restricted"),
    df = function(i) i$df_restricted - i$df_full),
  "Pearson correlation" = list(
    planner = "ss_buc_correlation", statistic = "t_observed",
    family = "correlation",
    planning_term = "rho_planning",
    echo = c("r_observed | t_observed", "N"),
    df = function(i) i$N)
)

# Resolve an entry whose sampling family depends on the plan: the regression
# entries name a fixed and a random family, and the echoed random_predictors
# row picks one. Every caller reads the spec through this (directly, or through
# .bucss_design_args()), since an unresolved entry carries two families. The
# row is required in 'echo', so a regression plan without it fails the shape
# check that ss_buc_sensitivity() and plot() apply first; were one to arrive
# here anyway, it would resolve to the fixed family, the default.
.bucss_resolve_spec <- function(spec, inputs) {
  if (length(spec$family) > 1L) {
    random <- isTRUE(as.logical(inputs$random_predictors))
    spec$family <- unname(spec$family[if (random) "random" else "fixed"])
  }
  spec
}

# The echoed rows a plan of this design must carry for the effect it tests:
# the design's 'echo' plus that effect's entry in 'echo_effect', if any. Each
# element is one requirement in the registry's grammar (see the top of this
# file); .is_plan_shape() checks that one of its forms is present, each row
# of that form exactly once.
.bucss_echo_requirements <- function(spec, effect = NULL) {
  extra <- if (is.character(effect) && length(effect) == 1L && !is.na(effect))
    spec$echo_effect[[effect]]
  c(spec$echo, extra)
}

# Pull the planning inputs back off a result as a named list, one element per
# echoed row. The design table above reads these row names, so they stay
# exactly as the constructor wrote them; rebuilding the planner's own argument
# names is .bucss_design_args()'s job.
.bucss_design_inputs <- function(object) {
  keep <- !(object$term %in% .BUCSS_COMPUTED_TERMS)
  out <- as.list(object$value[keep])
  names(out) <- object$term[keep]
  out
}

# Turn those rows into the planner's arguments. The two places the row names
# and the formals differ: a length-2 prior 'n' is echoed as n_1/n_2 by the
# independent t (the Welch planner's formals really are n_1 and n_2), and the
# effect travels as an attribute rather than a row.
.bucss_design_args <- function(inputs, spec, effect) {
  spec <- .bucss_resolve_spec(spec, inputs)
  args <- inputs
  if (spec$planner == "ss_buc_independent_t" && !is.null(args$n_1)) {
    args$n <- c(args$n_1, args$n_2)
    args$n_1 <- NULL
    args$n_2 <- NULL
  }
  if (!is.null(effect)) args$effect <- effect
  # The predictor frame is echoed as a numeric 1/0 row (the value column is
  # numeric throughout); the planners take TRUE or FALSE.
  if (!is.null(args$random_predictors))
    args$random_predictors <- as.logical(args$random_predictors)
  # The correlation planner accepts either the correlation or its t; simulated
  # statistics arrive as a t, so a call stated in r is restated in t.
  if (spec$family == "correlation" && !is.null(args$r_observed)) {
    r <- abs(args$r_observed)
    args$t_observed <- r * sqrt((args$N - 2) / (1 - r^2))
    args$r_observed <- NULL
  }
  # Likewise the MANOVA planner accepts Hotelling's T squared or its exact F;
  # the simulation draws F, so a call stated in T squared is restated in F.
  if (spec$planner == "ss_buc_manova" && !is.null(args$T2_observed)) {
    args$F_observed <- args$T2_observed * (args$N - args$p_variables - 1) /
      (args$p_variables * (args$N - 2))
    args$T2_observed <- NULL
  }
  args
}

# Does a replayed planner call reproduce the plan? Only the computed rows are
# compared (a correlation stated in r is replayed in t, so the echoed input
# rows legitimately carry a different label), the planning effect size among
# them, and they are matched by term, so a plan whose rows were reordered
# still replays.
.replay_matches <- function(replay, object) {
  computed <- .BUCSS_COMPUTED_TERMS
  terms_object <- object$term[object$term %in% computed]
  terms_replay <- replay$term[replay$term %in% computed]
  if (!setequal(terms_object, terms_replay) ||
      length(terms_object) != length(terms_replay)) return(FALSE)
  got <- replay$value[match(terms_object, replay$term)]
  want <- object$value[match(terms_object, object$term)]
  isTRUE(all.equal(got, want, tolerance = 1e-10))
}
