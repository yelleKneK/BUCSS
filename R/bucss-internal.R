# Internal constructor and print method shared by every ss_buc_* function.

# The size-row vocabulary, mirroring DMAR's .SS_POWER_SIZE_TERMS: the size row
# is named for its unit rather than carrying one generic name, so a BUCSS
# table reads like a DMAR table and DMAR's tidiers could absorb it. Every
# consumer that needs "the" sample size (the legacy `[[` method, the broom
# verbs, the characterization harness, the oracles' helper) looks the row up
# through this vector rather than hard-coding a name.
.BUCSS_SIZE_TERMS <- c("necessary_n_per_group", "necessary_n_per_cell",
                       "necessary_N")

# The planning effect size vocabulary: every plan carries exactly one of these
# rows, directly after ncp_adjusted. It is the effect size the plan is built
# on, ncp_adjusted restated in the design's usual metric at the balanced prior
# design that parameter belongs to (for the planners that round an unbalanced
# prior study both ways, the branch whose parameter is reported), so the
# ordinary power of a study of the returned size at this effect is exactly
# actual_power. The names say what the value is: never "corrected" (at an
# assurance above .5 it is a lower confidence limit, not an estimate) and not
# "_adjusted" (R2_adjusted names the Ezekiel/Wherry adjusted R squared).
# Each design's own term is recorded as 'planning_term' in
# R/bucss-design-registry.R. This one table serves the constructor (which
# checks the term), the shape check (which requires exactly one such row, and
# the design's own), the replays (which treat the row as computed, never as
# an input), and planning_sentence() (the value is the metric's name as the
# sentence writes it); print(), knit_print(), tidy(), and glance() reach it
# through the shape check, and glance() gives the row its own column.
.BUCSS_EFFECT_LABELS <- c(
  smd_planning = "a standardized mean difference",
  dz_planning = "a standardized mean of the paired differences",
  std_mean_planning = "a standardized mean",
  f_planning = "Cohen's f",
  f2_planning = "Cohen's f squared",
  R2_planning = "a squared multiple correlation",
  rho_planning = "a correlation",
  D2_planning = "a squared Mahalanobis distance",
  rmsea_d_planning = paste("the root mean square error of approximation of",
                           "the difference test"))
.BUCSS_EFFECT_TERMS <- names(.BUCSS_EFFECT_LABELS)
# The planning effect sizes that cannot exceed 1, which the sentence writes
# without a leading zero.
.BUCSS_EFFECT_BOUNDED <- c("R2_planning", "rho_planning")

# Every row a planner computes rather than echoes. The replays rebuild a
# planner call from the other rows and compare these.
.BUCSS_COMPUTED_TERMS <- c(.BUCSS_SIZE_TERMS, "total_N", "actual_power",
                           "ncp_adjusted", .BUCSS_EFFECT_TERMS)

# The conversions from a noncentrality parameter to each design's usual effect
# size, each at the balanced prior design the parameter belongs to. Every
# planner scales its parameter in proportion to the sample size (n_rep/n_prior
# on the F and chi square scales, its square root on the t scale), so each
# quantity is a population value that does not depend on the sample size. The
# t-scale ones are magnitudes, since the planners correct |t|.
#   .es_smd: Cohen's delta of two groups of n each, lambda = delta sqrt(n/2).
#   .es_std_mean: the standardized mean (of the paired differences, for a
#     paired design), lambda = delta sqrt(N).
#   .es_smd_welch: the mean difference over the root mean square of the two
#     standard deviations (Bonett, 2008), with the first group's standard
#     deviation set to 1 and 'ratio' the second's; with n per group in the
#     planned study its noncentrality parameter is this value times sqrt(n/2).
#   .es_f: Cohen's f of an effect relative to its own error term, with
#     lambda = k N f^2, where N counts participants and k is the error degrees
#     of freedom per participant: 1 for a between-subjects effect (so lambda =
#     N f^2, DMAR's convention), the effect's numerator df for a
#     within-subjects ANOVA, levels_within - 1 for a split-plot within or
#     interaction effect. f^2/(1 + f^2) is the population partial eta squared.
#   .es_f2: Cohen's f squared of a coefficient or set of coefficients,
#     lambda = N f^2, which is also the population value when the predictors
#     are random (the mixture of .p_random_f() has c = f^2).
#   .es_R2, .es_rho: the population squared multiple correlation and the
#     correlation, lambda = N rho^2/(1 - rho^2).
#   .es_D2: the squared Mahalanobis distance of two groups of n each,
#     lambda = (n/2) D^2.
#   .es_rmsea_d: the root mean square error of approximation of a nested model
#     difference test with 'df' constraints, lambda = df N RMSEA^2 (N, not
#     N - 1, the scaling ss_buc_chisq_diff() plans with).
.es_smd <- function(ncp, n_per_group) ncp * sqrt(2 / n_per_group)
.es_std_mean <- function(ncp, N) ncp / sqrt(N)
.es_smd_welch <- function(ncp, n_1, n_2, ratio)
  ncp * sqrt((1 / n_1 + ratio^2 / n_2) / ((1 + ratio^2) / 2))
.es_f <- function(ncp, N, k = 1) sqrt(ncp / (k * N))
.es_f2 <- function(ncp, N) ncp / N
.es_R2 <- function(ncp, N) ncp / (ncp + N)
.es_rho <- function(ncp, N) sqrt(ncp / (ncp + N))
.es_D2 <- function(ncp, n_per_group) 2 * ncp / n_per_group
.es_rmsea_d <- function(ncp, N, df) sqrt(ncp / (df * N))

# Build the tidy result object returned by every ss_buc_* function.
#
# Everything numeric lives as rows of the `value` column, in DMAR's echoed-
# inputs style: the design results first (the size row named by `size_term`,
# the implied `total_N` for per-group and per-cell designs, the conservative
# `actual_power`, the bias and uncertainty adjusted prior-study
# noncentrality parameter `ncp_adjusted`, and the planning effect size, the
# row named by `effect_size_term`), followed by rows echoing the
# user's planning inputs, so the assumptions travel with the result through
# subsetting and CSV export. A length-2 prior `n` echoes as `n_1`/`n_2` rows.
# Only non-numeric metadata (the human design label, the unit, the effect
# tested) plus the assurance ceiling and the planned test's df travel on
# attributes.
.bucss_plan_result <- function(sample_size, size_term, ncp,
                                effect_size_term, effect_size,
                                actual_power = NULL, design, sample_size_unit,
                                effect = NULL, assurance_ceiling = NULL,
                                total_n = NULL, df_effect = NULL, df_error = NULL,
                                inputs = list()) {
  stopifnot(size_term %in% .BUCSS_SIZE_TERMS,
            length(effect_size_term) == 1L,
            effect_size_term %in% .BUCSS_EFFECT_TERMS,
            is.numeric(effect_size), length(effect_size) == 1L,
            is.finite(effect_size))
  terms <- size_term
  values <- sample_size
  if (!is.null(total_n)) {
    terms <- c(terms, "total_N")
    values <- c(values, total_n)
  }
  if (!is.null(actual_power)) {
    terms <- c(terms, "actual_power")
    values <- c(values, actual_power)
  }
  terms <- c(terms, "ncp_adjusted", effect_size_term)
  values <- c(values, ncp, effect_size)
  inputs <- inputs[!vapply(inputs, is.null, logical(1))]
  for (nm in names(inputs)) {
    v <- inputs[[nm]]
    if (nm == "n" && length(v) == 2L) {
      terms <- c(terms, "n_1", "n_2")
      values <- c(values, v[1], v[2])
    } else {
      terms <- c(terms, nm)
      values <- c(values, v)
    }
  }
  out <- data.frame(term = terms, value = values, stringsAsFactors = FALSE)
  attr(out, "design") <- design
  attr(out, "sample_size_unit") <- sample_size_unit
  attr(out, "effect") <- effect
  attr(out, "assurance_ceiling") <- assurance_ceiling
  attr(out, "df_effect") <- df_effect
  attr(out, "df_error") <- df_error
  class(out) <- c("bucss_plan", "data.frame")
  out
}

# Check that a value is a single finite number, so a typo (NA, Inf, a vector)
# gets a friendly error naming the argument rather than a cryptic base-R
# condition failure ("missing value where TRUE/FALSE needed", "the condition
# has length > 1"). Used for the observed statistic and, via
# .validate_planning_inputs(), the four planning inputs.
.check_scalar_finite <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop("'", name, "' must be a single finite numeric value.", call. = FALSE)
  invisible(x)
}

# Check that a design count (a sample size, number of levels, cells, groups,
# predictors, or degrees of freedom) is a single finite whole number of at
# least 'min'. Guards against inputs like p = 2.5, which would otherwise flow
# through the planned-n search and return a non-integer "sample size".
.check_count <- function(x, name, min = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x))
    stop("'", name, "' must be a single whole number.", call. = FALSE)
  if (x < min)
    stop("'", name, "' must be at least ", min, ".", call. = FALSE)
  invisible(x)
}

# Validate the sphericity parameter 'epsilon' of a within-subjects effect whose
# contrasts span 'dim_within' dimensions, and return the value to use. The
# dimension is the effect's numerator df for a purely within-subjects effect,
# and the within-subjects factor's df for a between x within interaction
# (whose numerator df is larger, but whose covariance is that of the
# within-subjects contrasts, so a split-plot interaction with 4 groups and 3
# occasions has df 6 and yet cannot go below 1/2). Box's (1954) lower bound
# is 1/dim_within, and a one-dimensional effect (a two-level factor or a
# single contrast) is always spherical. A ratio of reported degrees of freedom
# can miss 1 by a rounding error; within the tolerance it is read as exactly 1,
# so 1 + 1e-12 and 1 - 1e-12 both give exactly the plan that assumes
# sphericity. Below the lower bound the allowance is wider, half a unit in the
# third decimal place, because the bound is usually a repeating decimal: an
# estimate reported at the bound as .333 (or .3333) for 1/3 is read as the
# bound, the value the error message below names. The allowance never exceeds
# half the bound, so zero and negative values are refused even for an effect
# of 1,000 or more within-subjects degrees of freedom, whose bound is .001 or
# less.
.check_epsilon <- function(epsilon, dim_within) {
  .check_scalar_finite(epsilon, "epsilon")
  tol <- sqrt(.Machine$double.eps)
  if (epsilon > 1 + tol)
    stop("'epsilon' must be at most 1 (1 assumes sphericity; a Huynh-Feldt ",
         "estimate above 1 is truncated to 1).", call. = FALSE)
  if (epsilon >= 1 - tol) return(1)
  if (dim_within == 1)
    stop("This effect's within-subjects contrasts have a single degree of ",
         "freedom (as for a two-level factor or a single contrast), and such ",
         "an effect always satisfies sphericity, so 'epsilon' must be 1.",
         call. = FALSE)
  lower <- 1 / dim_within
  if (epsilon < lower - min(5e-4, lower / 2))
    stop("'epsilon' must be at least 1/", dim_within, " (about ",
         format(signif(lower, 4)), "), the lower bound for an effect whose ",
         "within-subjects contrasts have ", dim_within, " degrees of freedom. ",
         "If a paper reports corrected degrees of freedom, 'epsilon' is the ",
         "reported numerator degrees of freedom divided by the uncorrected ",
         "numerator degrees of freedom.", call. = FALSE)
  min(max(epsilon, lower), 1)
}

# The between-subjects effect of a split-plot design is tested on subject
# means, which involves no sphericity, so 'epsilon' does not apply to it. A
# supplied value is checked to be a number, a value other than 1 is ignored
# with a warning, and the plan is exactly the one without it (it carries no
# 'epsilon' row). 'supplied' is !missing(epsilon) in the planner.
.ignore_epsilon <- function(epsilon, supplied, effect_label) {
  if (supplied) {
    .check_scalar_finite(epsilon, "epsilon")
    if (epsilon != 1)
      warning("'epsilon' applies only to effects with a within-subjects ",
              "component; ", effect_label, " is tested on subject means, ",
              "which involves no sphericity, so 'epsilon' was ignored.",
              call. = FALSE)
  }
  1
}

# Seed discipline for every function that takes a 'seed' argument. A NULL
# seed is a no-op: the draws come from the caller's current random number
# stream, exactly as if the function had no 'seed' argument. A supplied seed
# (a whole number of at least 1) is set for the duration of the calling
# function only, and the caller's generator state, including the generator
# kind, is restored when that function exits, or removed again if the caller
# had none, so a reproducible call does not perturb the script around it. The
# saving and restoring is delegated to withr::local_seed(), which is why
# withr (>= 3.0.3) is in Imports: in withr 3.0.2 and earlier, local_seed()
# called from inside a function did not restore anything. BUCSS itself never
# reads or writes the random number state. Call this from the body of the
# function whose exit should restore the state (the default 'envir').
.bucss_local_seed <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) return(invisible(NULL))
  .check_count(seed, "seed")
  # set.seed() takes an integer, so a larger whole number would reach it as
  # NA, with base R's coercion warning and error rather than this one.
  if (seed > .Machine$integer.max)
    stop("'seed' must be at most ", .Machine$integer.max, ".", call. = FALSE)
  withr::local_seed(seed, .local_envir = envir)
  invisible(seed)
}

# Validate the planning inputs shared by every ss_buc_* function and apply the
# documented coercions, so the checks stay identical across all planners. The
# three coercions are load-bearing: 'alpha_prior == 1' becomes .999 to model
# "no publication bias" without a degenerate truncation, and 'assurance' or
# 'power' entered as a percentage (>= 1) is divided by 100 so that, for example,
# 80 means .80. Returns the (possibly coerced) values; 'alpha_prior_input' is
# the value the user supplied, echoed back in the result's planning inputs.
# call. = FALSE keeps this internal helper out of the error the user sees.
.validate_planning_inputs <- function(alpha_prior, alpha_planned, assurance,
                                      desired_power) {
  .check_scalar_finite(alpha_prior, "alpha_prior")
  .check_scalar_finite(alpha_planned, "alpha_planned")
  .check_scalar_finite(assurance, "assurance")
  .check_scalar_finite(desired_power, "desired_power")
  if (alpha_prior > 1 | alpha_prior <= 0) stop("There is a problem with 'alpha_prior' of the prior study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).", call. = FALSE)
  alpha_prior_input <- alpha_prior
  if (alpha_prior == 1) alpha_prior <- .999

  if (alpha_planned >= 1 | alpha_planned <= 0) stop("There is a problem with 'alpha_planned' of the planned study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).", call. = FALSE)

  if (assurance >= 1) assurance <- assurance / 100
  if (assurance <= 0 | assurance >= 1) stop("There is a problem with 'assurance' (i.e., the proportion of times statistical power is at or above the desired value), please specify as a value strictly between 0 and 1 (the default is .80). An 'assurance' of exactly 0 or 1 does not define a plannable study; note that 100 entered as a percentage means 1.", call. = FALSE)
  if (assurance < .5) warning("The 'assurance' you have entered is < .5, which implies you will have under a 50% chance at achieving your desired level of power.", call. = FALSE)

  if (desired_power >= 1) desired_power <- desired_power / 100
  if (desired_power <= 0 | desired_power >= 1) stop("There is a problem with 'desired_power' (i.e., desired statistical power), please specify as a value strictly between 0 and 1 (the default is .80). A planned power of exactly 0 or 1 is not attainable; note that 100 entered as a percentage means 1.", call. = FALSE)

  list(alpha_prior = alpha_prior, alpha_prior_input = alpha_prior_input,
       assurance = assurance, desired_power = desired_power)
}

# Map common shorthand for the `effect` argument to its canonical value, so a
# user may type "A", "b", or "AxB" instead of the documented
# "factor_A"/"factor_B"/"interaction" (or "between"/"within"/"interaction").
# Case- and separator-insensitive: "A x B", "axb", and "AB" all reach
# "interaction". Anything not in the shorthand table falls through to
# match.arg(), which accepts the canonical values (and their unambiguous
# prefixes) and otherwise errors helpfully. The canonical forms stay the ones
# the documentation suggests; the shorthand is accepted, not advertised.
.match_effect <- function(effect, choices) {
  raw <- effect[1]
  key <- tolower(gsub("[^[:alnum:]]", "", raw))
  syn <- switch(
    paste(choices, collapse = "|"),
    "factor_A|factor_B|interaction" = c(
      a = "factor_A", factora = "factor_A",
      b = "factor_B", factorb = "factor_B",
      axb = "interaction", ab = "interaction", int = "interaction"),
    "between|within|interaction" = c(
      bs = "between", betweensubjects = "between",
      ws = "within", withinsubjects = "within",
      bxw = "interaction", wxb = "interaction", int = "interaction"),
    NULL)
  if (!is.null(syn) && key %in% names(syn)) return(syn[[key]])
  match.arg(raw, choices)
}

# Rebuild a planner's result under a different design label, size-row name,
# unit, and echoed inputs, without recomputing anything. Used by the planners
# that are exactly a relabeling of another design (the one-sample t is the
# paired t; the ANCOVA is the general between-subjects planner with the
# covariates carried as nuisance parameters), so the numerical work lives in
# one place and only the presentation differs. The planning effect size row
# carries over with its value unchanged; 'effect_size_term' renames it where
# the relabeled design names the same quantity differently (the paired
# planner's standardized mean of the differences is the one-sample
# standardized mean). Left NULL, the term carries over too (the ANCOVA's
# f_planning is the general planner's, relative to the covariate-adjusted
# error variance).
.relabel_bucss_result <- function(res, design, sample_size_unit, size_term,
                                  inputs, effect = NULL,
                                  effect_size_term = NULL) {
  size <- res$value[res$term %in% .BUCSS_SIZE_TERMS][1]
  total_n <- res$value[res$term == "total_N"]
  if (length(total_n) == 0L) total_n <- NULL
  actual_power <- res$value[res$term == "actual_power"]
  if (length(actual_power) == 0L) actual_power <- NULL
  is_effect_size <- res$term %in% .BUCSS_EFFECT_TERMS
  if (is.null(effect_size_term)) effect_size_term <- res$term[is_effect_size]
  .bucss_plan_result(
    sample_size = size,
    size_term = size_term,
    effect_size_term = effect_size_term,
    effect_size = res$value[is_effect_size],
    actual_power = actual_power,
    ncp = res$value[res$term == "ncp_adjusted"],
    design = design,
    sample_size_unit = sample_size_unit,
    effect = effect,
    assurance_ceiling = attr(res, "assurance_ceiling"),
    total_n = total_n,
    df_effect = attr(res, "df_effect"),
    df_error = attr(res, "df_error"),
    inputs = inputs
  )
}

# Upper safety bound on the root-finding bracket in .solve_ncp_assurance(). It
# is far beyond any noncentrality parameter a real prior study implies; reaching
# it signals a likely input error (or an essentially deterministic prior) rather
# than a genuine solution, and the solver warns and returns it instead of
# expanding without bound.
.bucss_ncp_upper_cap <- 1e7

# Build the truncated-likelihood mean TM as a function of the noncentrality
# parameter for the distributional shapes the planners use. TM(ncp) is the
# truncated cumulative density evaluated at the observed statistic: the area of
# the (publication-)truncated noncentral distribution between the critical value
# and the observed statistic, divided by the power (the untruncated area beyond
# the critical value). It is monotone decreasing in 'ncp', from its ceiling
# TM(0) = 1 - p_prior / alpha_prior down to 0. Every returned closure must be
# vectorized over 'ncp', because .solve_ncp_assurance() calls it both scalar
# (uniroot) and vectorized (the grid fallback).
#
# .tm_f covers the one-sided noncentral F designs (the ANOVA and regression
# planners; the single regression coefficient passes the squared t as 'stat'
# against an F critical value). .tm_t covers the two-sided noncentral t designs
# (the independent and paired t tests), summing both tails exactly as the 1.x
# code did. .tm_chisq is the analogue for a nested model comparison, and
# .tm_random_f the analogue for sampled (random) predictors, whose statistic is
# not noncentral F at all (see below): the Pearson correlation through
# .tm_correlation, and the regression planners with random_predictors = TRUE
# through .tm_regression.
#
# 'epsilon' (default 1) is the sphericity parameter of a within-subjects
# effect. Under the Muller and Barton (1989) approximation, the convention of
# DMAR's within-subjects power functions, the statistic is noncentral F on
# epsilon * df1 and epsilon * df2 degrees of freedom with noncentrality
# epsilon * ncp, where ncp stays on the sphericity scale (the value it would
# have if the contrast covariance were spherical with the same average
# variance), so 'ncp_adjusted' keeps its meaning for every epsilon. 'crit'
# must already be the critical value at the corrected degrees of freedom. At
# epsilon = 1 every product below is exact (1 * x is x in floating point), so
# the closure is the sphericity-assumed one bit for bit.
.tm_f <- function(stat, crit, df1, df2, epsilon = 1) {
  df1 <- epsilon * df1
  df2 <- epsilon * df2
  function(ncp) {
    ncp <- epsilon * ncp
    power <- 1 - pf(crit, df1 = df1, df2 = df2, ncp = ncp)
    area_above <- 1 - pf(stat, df1 = df1, df2 = df2, ncp = ncp)
    (power - area_above) / power
  }
}

.tm_t <- function(stat, crit, df) {
  # The two tail masses below are taken beyond +stat and -stat, which is
  # only the two-sided tail area when stat is nonnegative; a raw negative
  # statistic would produce TM values outside [0, 1]. The publication rule
  # is symmetric in the sign of t, so the magnitude is the statistic.
  stat <- abs(stat)
  function(ncp) {
    .quiet_pnt({
      power <- (1 - pt(crit, df = df, ncp = ncp)) + pt(-crit, df = df, ncp = ncp)
      area_above <- (1 - pt(stat, df = df, ncp = ncp)) + pt(-stat, df = df, ncp = ncp)
    })
    (power - area_above) / power
  }
}

# Evaluate 'expr' with R's "full precision may not have been achieved in
# 'pnt{final}'" warning muffled. stats::pt() raises it while the root finder
# probes a small noncentrality parameter against a large observed t (t of 12
# or more with a few dozen per group), a region whose tail mass is far below
# anything TM depends on; the solved parameter agrees with independent
# quadrature to about 1e-11 in every case examined (2026-09-28). Only that
# message is muffled, so any other warning still reaches the user.
.quiet_pnt <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    # Matched on the call (a warning from pt()) and on the routine's name,
    # which R's translations keep, so it works in any language.
    cl <- conditionCall(w)
    if (is.call(cl) && identical(cl[[1L]], as.name("pt")) &&
        grepl("pnt", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}

.tm_chisq <- function(stat, crit, df) {
  function(ncp) {
    power <- 1 - pchisq(crit, df = df, ncp = ncp)
    area_above <- 1 - pchisq(stat, df = df, ncp = ncp)
    (power - area_above) / power
  }
}

# Exact distribution of the F statistic for testing q regression coefficients
# (a single coefficient, a subset, or all p of them) when the predictors are
# sampled from a multivariate normal distribution along with the outcome: the
# random-predictor frame of an observational study, in closed form. A Pearson
# correlation is the case q = p = 1.
#
# Conditional on the sampled predictors, the statistic is noncentral
# F(q, N - p - 1) with noncentrality lambda_X = c * W, where c is the
# population partial f squared of the tested set (for a correlation,
# c = rho^2/(1 - rho^2)) and W = beta_2' X_2' (I - H_1) X_2 beta_2 divided by
# its population scale is chi square on m = N - 1 - (p - q) df: the tested
# predictors' spread after removing the intercept and the p - q other
# predictors. So the marginal distribution is a chi square mixture of
# noncentral F's, not the noncentral F itself. Treating it as noncentral F
# (which is exact only when the predictors are fixed by design) overstates the
# power, because it drops the extra dispersion the sampled predictors
# contribute.
#
# The mixture has a closed form. The noncentral F CDF is a Poisson(lambda/2)
# mixture of central betas; mixing a Poisson rate over a gamma gives a negative
# binomial, and lambda/2 = c * W/2 with W/2 ~ Gamma(m/2, scale = c). The
# Poisson weights therefore become negative binomial weights with size m/2 and
# probability 1/(1 + c):
#
#   P(F <= x) = sum_k dnbinom(k, size = m/2, prob = 1/(1 + c)) *
#               pbeta(q x/(q x + df2), q/2 + k, df2/2)
#
# which is Fisher's exact distribution of the squared multiple correlation
# (see Gatsonis & Sampson, 1989) at a sample size reduced by the p - q
# untested predictors. This is exact rather than approximate, needs nothing
# beyond stats, and converges geometrically at rate c/(1 + c). Numerical
# quadrature was the obvious alternative and is a trap: integrate() over
# (0, Inf) silently returns 0 once N is a few hundred, because the adaptive
# rule samples where a chi square on m df has no mass.
#
# The caller passes 'prob' = 1/(1 + c) rather than c, computed as
# 1 - lambda/(lambda + N) from the noncentrality parameter lambda = c * N
# every planner reports; that is the expression the correlation planner has
# always used, so its results are unchanged to the last bit.
#
# The series is truncated where the remaining negative binomial mass is below
# 1e-14; the omitted terms carry pbeta values smaller than the last retained
# one, so the truncation error is far below that mass. The number of terms is
# governed by the noncentrality parameter, lambda/2 on average. A planned study
# near its target power carries lambda near 8, but the prior study's side
# carries lambda = c * N, which reaches the millions for a strong effect in a
# very large sample (r = .9 with N near half a million). There the mass sits
# far from k = 0, so the sum runs over a window around the mode, from the
# 1e-15 quantile to the 1 - 1e-14 quantile, whose width grows only with the
# square root of lambda. An earlier fixed cap of 1e6 terms instead dropped
# nearly all of the mass there and returned a wrong plan with no warning. Below
# a million terms the window starts at 0, exactly as the full sum does.
#
# 'what' and 'arg' name the prior statistic in the one error this can raise.
.p_random_f <- function(x, df1, df2, m, prob, what = "result",
                        arg = "the observed statistic") {
  if (prob <= 0) return(0)
  size <- m / 2
  k_hi <- qnbinom(1 - 1e-14, size = size, prob = prob)
  k_lo <- if (k_hi > 1e6) qnbinom(1e-15, size = size, prob = prob) else 0
  if (k_hi - k_lo > 5e7)
    stop("The prior study's ", what, " and sample size imply a noncentrality ",
         "parameter too large to evaluate. Please check ", arg, " and 'N'.",
         call. = FALSE)
  k <- k_lo:k_hi
  sum(dnbinom(k, size = size, prob = prob) *
        pbeta(df1 * x / (df1 * x + df2), df1 / 2 + k, df2 / 2))
}

# The correlation case, kept under its own name (the shipped verification
# script and the tests call it): q = 1, df2 = n - 2, m = n - 1.
.p_correlation_f <- function(q, n, rho2) {
  .p_random_f(q, 1, n - 2, n - 1, 1 - rho2, what = "correlation",
              arg = "'r_observed' (or 't_observed')")
}

# TM under random predictors, parameterized (like every other planner) by the
# noncentrality parameter the adjusted effect implies at the prior study's
# sample size, lambda = c * N. Solving on that scale rather than on the c or
# rho^2 scale matters: uniroot()'s tolerance is absolute, and a root of order
# 1e-6, which is what a plan near the assurance ceiling implies, would be
# resolved to only a percent or two. On the lambda scale the root is of order 1
# to 100 and the shared solver needs no change. At lambda = 0 the mixture is
# the central F, so the ceiling TM(0) equals .tm_f()'s.
.tm_random_f <- function(stat, crit, df1, df2, m, N, what = "result",
                         arg = "the observed statistic") {
  function(ncp) {
    vapply(ncp, function(lambda) {
      prob <- 1 - lambda / (lambda + N)
      power <- 1 - .p_random_f(crit, df1, df2, m, prob, what, arg)
      area_above <- 1 - .p_random_f(stat, df1, df2, m, prob, what, arg)
      (power - area_above) / power
    }, numeric(1))
  }
}

.tm_correlation <- function(stat, crit, N) {
  .tm_random_f(stat, crit, 1, N - 2, N - 1, N, what = "correlation",
               arg = "'r_observed' (or 't_observed')")
}

# The three regression planners share one frame switch. With random
# predictors the prior study's TM and the planned study's power both use the
# mixture above, with m = N - 1 - (p - q) at the prior N and at the planned
# N alike; with fixed predictors they use the noncentral F, the planned
# study's parameter scaled by n_rep/N. 'q' is the numerator df (1 for a
# coefficient, p for R2, p_joint for a joint test).
.tm_regression <- function(stat, crit, q, N, p, random_predictors, what, arg) {
  if (random_predictors)
    .tm_random_f(stat, crit, q, N - p - 1, N - 1 - (p - q), N, what, arg)
  else
    .tm_f(stat, crit, q, N - p - 1)
}

.power_regression <- function(ncp, N, p, q, alpha_planned, random_predictors,
                              what, arg) {
  prob <- 1 - ncp / (ncp + N)
  function(n_rep) {
    denom_df <- n_rep - p - 1
    critical_F <- qf(1 - alpha_planned, df1 = q, df2 = denom_df)
    if (random_predictors)
      1 - .p_random_f(critical_F, q, denom_df, n_rep - 1 - (p - q), prob,
                      what, arg)
    else
      1 - pf(critical_F, df1 = q, df2 = denom_df, ncp = (n_rep / N) * ncp)
  }
}

# A single TRUE or FALSE, for a logical switch such as 'random_predictors'.
.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  invisible(x)
}

# Solve TM(ncp) = assurance for the bias and uncertainty adjusted noncentrality
# parameter. Because TM is monotone decreasing from its ceiling TM(0) to 0, the
# equation has a unique root whenever 'assurance' is below the ceiling. This root
# is found directly with uniroot() rather than by selecting the closest value on
# the fixed 0..100 grid the 1.x code used, so the
# adjusted parameter is no longer capped at 100 and no longer snapped to the grid
# resolution. The upper bracket starts at 100 (the historical grid maximum) and
# doubles until TM falls below assurance, with a generous safety cap so a
# pathological input cannot loop forever. Returns the root and the ceiling
# TM(0); the caller routes a zero (at-or-above-ceiling) result through
# .stop_zero_ncp(), exactly as the grid code routed a selected NCP of 0. A
# uniroot() failure (it should not occur once bracketed) falls back to a fine
# grid search over the established bracket. The ceiling it returns equals the
# 1.x max(TM), so the zero-NCP message is unchanged.
.solve_ncp_assurance <- function(tm_fun, assurance) {
  ceiling_tm <- tm_fun(0)
  if (!is.finite(ceiling_tm)) ceiling_tm <- 1
  if (ceiling_tm <= assurance) return(list(ncp = 0, ceiling = ceiling_tm))

  f <- function(x) tm_fun(x) - assurance
  upper <- 100
  while (upper < .bucss_ncp_upper_cap && f(upper) > 0) upper <- upper * 2
  if (f(upper) > 0) {
    warning("The bias and uncertainty corrected noncentrality parameter exceeds ",
            .bucss_ncp_upper_cap, "; returning ", .bucss_ncp_upper_cap, ". The ",
            "prior result implies an extremely large effect, so please verify ",
            "that the observed statistic and 'N' of the prior study are correct.",
            call. = FALSE)
    return(list(ncp = .bucss_ncp_upper_cap, ceiling = ceiling_tm))
  }
  root <- tryCatch(
    uniroot(f, lower = 0, upper = upper, tol = .Machine$double.eps^0.5)$root,
    error = function(e) {
      grid <- seq(0, upper, length.out = 200001L)
      tm <- tm_fun(grid)
      min(grid[which(abs(tm - assurance) == min(abs(tm - assurance)))])
    }
  )
  list(ncp = root, ceiling = ceiling_tm)
}

# Smallest planned-study sample size (in the design's own unit) whose power
# reaches 'target'. 'power_at' returns the planned study's power for a candidate
# size; power is monotone increasing in the size, so the first size at or above
# 'target' is the answer. 'start' is the smallest admissible size (2 for most
# designs; larger when nuisance-parameter or predictor degrees of freedom must
# be absorbed before the error df is positive). Every planner routes its
# planned-n search through this helper so the search, and its off-by-one, lives
# in one place rather than being re-derived per design; each planner supplies a
# 'power_at' closure carrying its own degrees-of-freedom and noncentrality
# scaling. The bracket doubles until the target is enclosed and then bisects,
# so the cost is logarithmic in the answer where the 1.x incremental search was
# linear (a plan near the assurance ceiling can need a sample size in the
# millions); monotonicity guarantees the same integer the incremental search
# returned. Arithmetic is in doubles so a huge answer cannot overflow the
# integer counter. Doubles hold every whole number exactly only up to 2^53;
# beyond it adjacent doubles are 2 apart, so a midpoint can equal an end point
# and the bisection would never finish. The bracket therefore stops doubling
# once it has passed 2^52 without reaching the target, so it never extends past
# 2^53, and the search stops with a message instead: the answer then exceeds
# 2^52 (about 4.5e15). Answers below that are always found exactly. Such a plan
# only arises with an assurance a hair's breadth below the ceiling.
.smallest_n_for_power <- function(power_at, target, start = 2L) {
  lo <- as.numeric(start)
  if (power_at(lo) >= target) return(lo)
  hi <- lo * 2                       # invariant: power_at(lo) < target
  while (power_at(hi) < target) {
    if (hi > 2^52)
      stop("The planned sample size would exceed about 4.5e15 (2^52), because ",
           "the requested 'assurance' is essentially at the largest value the ",
           "prior result can support. Please lower 'assurance' a little.",
           call. = FALSE)
    lo <- hi
    hi <- hi * 2
  }
  while (hi - lo > 1) {              # invariant: power_at(hi) >= target
    mid <- floor((lo + hi) / 2)
    if (mid <= lo || mid >= hi) break
    if (power_at(mid) < target) lo <- mid else hi <- mid
  }
  hi
}

# Stop because the prior study's statistic is not significant at 'alpha_prior',
# so the truncated likelihood the correction is built on does not include it.
# 'what' names the statistic ("t statistic", "F statistic", ...) and 'arg' the
# argument that carries it. 'p_own' is the prior study's p value at its own
# degrees of freedom and 'p_balanced' its p value at the balanced design the
# planner substitutes (the rounded-down one, which has the fewest degrees of
# freedom), both supplied only by the planners that reduce an unbalanced prior
# study to a balanced one: a statistic can be significant at the prior study's
# own degrees of freedom and not at the balanced design's, and the message then
# says so rather than calling the prior study nonsignificant, which it is not.
# Planning at the requested 'assurance' needs the balanced design's ceiling,
# 1 - p_balanced/alpha_prior, to exceed the assurance, so 'alpha_prior' must be
# strictly above p_balanced/(1 - assurance); at that threshold itself the
# planner refuses. The message names a value strictly above it (see
# .format_alpha_above()) and says that plans near it are very large; when the
# threshold is .999 or more it says instead that no practical 'alpha_prior'
# will do ('alpha_prior = 1' is read as .999). It mentions the recommended
# floor of .5 only when the own-df ceiling is below it. The condition carries
# the class 'bucss_nonsignificant' so code that loops over plans can catch it
# alone.
.stop_nonsignificant <- function(what, arg, p_own = NULL, alpha_prior = NULL,
                                 p_balanced = NULL, assurance = NULL) {
  fmt2 <- function(p) sub("^0", "", sprintf("%.2f", p))
  fmtp <- function(p) sub("^0", "", format(signif(p, 3), scientific = FALSE))
  # The rounding explanation applies only when the prior study is significant
  # at its own degrees of freedom and less so at the balanced design.
  if (!is.null(p_own) && !is.null(alpha_prior) && !is.null(p_balanced) &&
      p_own < alpha_prior && p_balanced > p_own) {
    ceiling_own <- .floor_assurance(1 - p_own / alpha_prior)
    needed <- if (!is.null(assurance)) p_balanced / (1 - assurance)
    msg <- paste0(
      "Your observed ", what, " is significant at the prior study's own degrees ",
      "of freedom (p = ", fmtp(p_own), "), but this planner reduces a prior ",
      "study with unequal group sizes to a balanced design with fewer error ",
      "degrees of freedom, and there its p value is ", fmtp(p_balanced),
      ", above 'alpha_prior'. ",
      if (!is.null(needed) && needed < .999)
        paste0("Planning from it at an 'assurance' of ", .format_alpha(assurance),
               " would need an 'alpha_prior' of ", .format_alpha_above(needed),
               " or more, which assumes a more lenient publication threshold. ",
               "Near that value the plan is very large, because the result ",
               "then barely supports the requested assurance; a larger ",
               "'alpha_prior' gives a practical plan, and so does a lower ",
               "'assurance', which lowers the 'alpha_prior' needed. ")
      else if (!is.null(needed))
        paste0("Planning from it at an 'assurance' of ", .format_alpha(assurance),
               " would need an 'alpha_prior' above .999, and 'alpha_prior = 1' ",
               "(no publication bias) is treated as .999, so no practical ",
               "'alpha_prior' allows it; a lower 'assurance' lowers the ",
               "'alpha_prior' needed. "),
      if (ceiling_own < .5)
        paste0("Even at its own degrees of freedom the result supports an ",
               "assurance of at most about ", fmt2(ceiling_own), ", below the ",
               "recommended floor of .5."))
  } else {
    msg <- paste0(
      "Your observed ", what, " is nonsignificant based on your specified ",
      "'alpha_prior' of the prior study. Please increase 'alpha_prior' so ",
      arg, " exceeds the critical value.")
  }
  stop(errorCondition(msg, class = "bucss_nonsignificant", call = NULL))
}

# The 'alpha_prior' .stop_nonsignificant() names, as text: the threshold 'x',
# p_balanced/(1 - assurance), plus a margin of one part in a million, rounded
# up to three significant digits. The value named is therefore strictly above
# the threshold, where the planner plans; rounding to the nearest three
# digits named a value at or below it about half the time, and the planner
# refused there. Never above .999, the largest value used as given; the
# caller passes only thresholds below .999, so .999 is itself above them.
# Display only: the value is text in a message, never stored.
.format_alpha_above <- function(x) {
  y <- x * (1 + 1e-6)
  k <- 2 - floor(log10(y))
  .format_alpha(min(ceiling(y * 10^k) / 10^k, .999))
}

# The assurance ceiling as the messages and footers state it: floored to two
# decimals, less a small tolerance, so the value stated is strictly below the
# ceiling and is itself a feasible assurance. A planner refuses when the
# ceiling is at or below the requested assurance, so a ceiling of .87, or one
# a hair either side of .87 that floating-point arithmetic would carry to 87
# hundredths, is stated as .86. Display only.
.floor_assurance <- function(x) floor(x * 100 - 1e-9) / 100

# Stop with a single, informative error when the bias and uncertainty
# correction drives the prior-study noncentrality parameter to zero. At the
# requested assurance the prior effect cannot be distinguished from zero, so no
# planned sample size can attain the target power. Every ss_buc_* function
# routes its zero-NCP guard through this helper so the wording stays identical
# across all call sites. Each passes 'max_assurance', the largest assurance the
# prior result can support: the truncated likelihood TM evaluated at a zero
# noncentrality parameter, max(TM) = 1 - p_prior / alpha_prior, where p_prior is
# the prior study's p-value. The message uses it to tell the user whether
# lowering 'assurance' can help (and how far) or whether 'alpha_prior' is the
# only remaining lever. The two-sided rounding planners pass the smaller of
# their two branches' ceilings, so the value the message names is itself a
# feasible assurance. The condition carries the class 'bucss_zero_ncp', so
# ss_buc_sensitivity() and plot() can catch this refusal and nothing else;
# call = NULL keeps the (uninformative) calling expression out of the message.
.stop_zero_ncp <- function(max_assurance) {
  # Floor to two decimals, strictly below the ceiling, so the reported value
  # is itself feasible (see .floor_assurance()); rounding up could name an
  # assurance that still fails.
  ceiling_assurance <- .floor_assurance(max_assurance)
  fmt <- function(p) sub("^0", "", sprintf("%.2f", p))

  if (ceiling_assurance >= 0.5) {
    lever <- paste0(
      "With this prior result and 'alpha_prior', the largest 'assurance' that ",
      "still yields a usable plan is about ", fmt(ceiling_assurance), "; re-run ",
      "with 'assurance' at or below that value, or raise 'alpha_prior' to lift ",
      "the ceiling. Re-running at a lower 'assurance' is not free: the plan it ",
      "returns carries that lower assurance, not the one you first asked for, ",
      "so the long run proportion of replications reaching your target power ",
      "falls accordingly. "
    )
  } else {
    lever <- paste0(
      "No 'assurance' down to its recommended floor of .5 yields a usable plan ",
      "with this prior result and 'alpha_prior', so lowering 'assurance' will ",
      "not help; raising 'alpha_prior' is the only remaining lever. If even ",
      "'alpha_prior = 1' leaves the corrected parameter at zero, the prior ",
      "result is too weak to plan from. "
    )
  }

  msg <- paste0(
    "Sample size planning is not possible with these settings. After correcting ",
    "the prior study's effect for publication bias and uncertainty at the ",
    "requested level of assurance, the corrected noncentrality parameter is ",
    "zero, which means there is not enough information in the prior result to ",
    "plan the sample size as desired. ",
    lever,
    "'alpha_prior' does not change the prior study, which is what it is; it is ",
    "your assumption about the publication threshold the study's literature ",
    "faced, so raising it toward 1 assumes less publication bias to undo (for ",
    "example, 'alpha_prior = .10' assumes the literature's journals would have ",
    "published a result significant at the .10 level, a more lenient policy than ",
    ".05, and 'alpha_prior = 1' assumes no publication bias at all). See ",
    "Anderson, Kelley, and Maxwell (2017) and Anderson and Kelley (2024)."
  )
  stop(errorCondition(msg, class = "bucss_zero_ncp", call = NULL))
}

# Check that an object still has the shape of one planner's answer before a
# method reads "the" sample size, power, or assurance from it, or replays the
# planner from its echoed inputs: exactly one size row; the implied total for
# a per-group or per-cell size, and none for a total; exactly one planning
# effect size row (one of .BUCSS_EFFECT_TERMS, and for a design in the
# registry the one its planner writes, its 'planning_term'); one row each of
# the results and of the four planning inputs every planner echoes; and, for a
# design in the registry, one row each of the inputs that design echoes (the
# 'echo' and 'echo_effect' fields in R/bucss-design-registry.R, read through
# .bucss_echo_requirements(), so a planning input added to a design is
# required by listing it there). Base R's
# rbind() and subsetting already drop the class from a stacked or modified
# result, but vctrs-based tools (dplyr::bind_rows(), dplyr::filter(), and the
# like) keep it, so a stack of two plans would otherwise be tidied and narrated
# as though every row belonged to the first, and a plan without its
# 'alpha_prior' row narrated as though no publication bias had been assumed.
# Stopping names the problem instead.
.BUCSS_PLAN_ROWS <- c("actual_power", "ncp_adjusted", "alpha_prior",
                      "alpha_planned", "assurance", "desired_power")

.is_plan_shape <- function(x) {
  # Exactly the two columns a plan has; an added column (dplyr::bind_cols(),
  # say) would otherwise be hidden by the print method.
  if (!is.data.frame(x) || !identical(names(x), c("term", "value"))) return(FALSE)
  terms <- x[["term"]]
  if (is.null(terms) || !is.character(terms)) return(FALSE)
  once <- function(nm) sum(terms == nm) == 1L
  size <- terms[terms %in% .BUCSS_SIZE_TERMS]
  if (length(size) != 1L || !all(vapply(.BUCSS_PLAN_ROWS, once, logical(1))))
    return(FALSE)
  if (sum(terms == "total_N") != as.integer(size != "necessary_N")) return(FALSE)
  effect_size <- terms[terms %in% .BUCSS_EFFECT_TERMS]
  if (length(effect_size) != 1L) return(FALSE)
  design <- attr(x, "design")
  spec <- if (is.character(design) && length(design) == 1L)
    .BUCSS_DESIGN_SPECS[[design]]
  if (is.null(spec)) return(TRUE)
  if (!identical(effect_size, spec$planning_term)) return(FALSE)
  requirements <- .bucss_echo_requirements(spec, attr(x, "effect"))
  all(vapply(requirements, function(req) {
    forms <- strsplit(strsplit(req, " | ", fixed = TRUE)[[1L]], " + ", fixed = TRUE)
    any(vapply(forms, function(rows) all(vapply(rows, once, logical(1))), logical(1)))
  }, logical(1)))
}

.check_plan_shape <- function(x) {
  if (!.is_plan_shape(x))
    stop("This object no longer has the shape of a single plan (one sample ",
         "size row, and one row each of the results and of the planning ",
         "inputs its planner echoed, such as 'alpha_prior' and 'assurance'); ",
         "it may be several plans stacked together or a plan with rows ",
         "removed. Apply this function to each plan separately, as its ",
         "planner returned it.", call. = FALSE)
  invisible(x)
}

# The footer line naming the predictor frame, for the plans that echo one (the
# regression planners' random_predictors row, 1 or 0). NULL for every other
# plan, so their footers are unchanged.
.predictor_footer <- function(x) {
  rp <- x$value[x$term == "random_predictors"]
  if (length(rp) != 1L) return(NULL)
  paste0("Predictors: ", if (rp == 1) "random (sampled along with the outcome)"
         else "fixed by design")
}

# Format a value column for display the way DMAR's dmar_tbl print does: whole
# numbers print clean (no decimal part), everything else at 'digits'
# significant figures, never scientific notation. Display only; the stored
# values keep full precision.
.format_bucss_value <- function(v, digits) {
  vapply(v, function(x) {
    if (is.na(x)) return(NA_character_)
    if (is.finite(x) && x == floor(x) && abs(x) < 1e15) {
      return(format(x, scientific = FALSE))
    }
    format(signif(x, digits), scientific = FALSE)
  }, character(1))
}

#' Print a Bias and Uncertainty Corrected Sample Size Result
#'
#' Prints the aligned \code{term}/\code{value} table (the planned-study
#' results followed by the rows echoing the planning inputs), then factual
#' footer lines naming the design, for the regression planners the predictor
#' frame (\code{random_predictors}), the unit the sample size is counted in,
#' and the largest assurance the prior result can support. Whole numbers
#' print clean; other values print at \code{getOption("bucss.digits", 3)}
#' significant figures. Only the display is rounded; the stored values keep
#' full precision.
#'
#' @param x An object of class \code{bucss_plan}.
#' @param digits Significant figures for non-integer values; defaults to
#'   \code{getOption("bucss.digits", 3)}.
#' @param ... Further arguments, ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' result <- ss_buc_independent_t(t_observed = 3, n = c(50, 55),
#'   assurance = .90)
#' result
#'
#' # Display only: the stored values keep full precision.
#' print(result, digits = 6)
#' result$value[result$term == "ncp_adjusted"]
#'
#' @aliases bucss_plan
#' @export
#' @keywords internal
print.bucss_plan <- function(x, digits = getOption("bucss.digits", 3L), ...) {
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 1)
    stop("'digits' must be a single positive number.", call. = FALSE)
  # An object that no longer has one plan's shape (a stack of plans kept in
  # the class by a vctrs-based tool, or a plan with columns renamed) prints as
  # the plain data frame it is, without footers that would describe only one
  # of the plans in it.
  if (!.is_plan_shape(x) || !is.numeric(x[["value"]])) {
    print(structure(x, class = "data.frame"))
    return(invisible(x))
  }
  shown <- data.frame(term = x$term,
                      value = .format_bucss_value(x$value, digits),
                      stringsAsFactors = FALSE)
  print.data.frame(shown, row.names = FALSE, right = FALSE)

  design <- attr(x, "design")
  effect <- attr(x, "effect")
  unit <- attr(x, "sample_size_unit")
  assurance_ceiling <- attr(x, "assurance_ceiling")
  cat("\n")
  if (!is.null(design)) {
    cat("Design: ", design,
        if (!is.null(effect)) paste0(" (", effect, ")"), "\n", sep = "")
  }
  predictors <- .predictor_footer(x)
  if (!is.null(predictors)) cat(predictors, "\n", sep = "")
  if (!is.null(unit)) cat("Sample size unit: ", unit, "\n", sep = "")
  if (!is.null(assurance_ceiling)) {
    cval <- .floor_assurance(assurance_ceiling)
    cat("Largest supportable assurance: ",
        sub("^0", "", sprintf("%.2f", cval)), "\n", sep = "")
  }
  invisible(x)
}
