# Internal constructor and print method shared by every ss_buc_* function.

# The size-row vocabulary, mirroring DMAR's .SS_POWER_SIZE_TERMS: the size row
# is named for its unit rather than carrying one generic name, so a BUCSS
# table reads like a DMAR table and DMAR's tidiers could absorb it. Every
# consumer that needs "the" sample size (the legacy `[[` method, the broom
# verbs, the characterization harness, the oracles' helper) looks the row up
# through this vector rather than hard-coding a name.
.BUCSS_SIZE_TERMS <- c("necessary_n_per_group", "necessary_n_per_cell",
                       "necessary_N")

# Build the tidy result object returned by every ss_buc_* function.
#
# Everything numeric lives as rows of the `value` column, in DMAR's echoed-
# inputs style: the design results first (the size row named by `size_term`,
# the implied `total_N` for per-group and per-cell designs, the conservative
# `actual_power`, and the bias and uncertainty adjusted prior-study
# noncentrality parameter `ncp_adjusted`), followed by rows echoing the
# user's planning inputs, so the assumptions travel with the result through
# subsetting and CSV export. A length-2 prior `n` echoes as `n_1`/`n_2` rows.
# Only non-numeric metadata (the human design label, the unit, the effect
# tested) plus the assurance ceiling and the planned test's df travel on
# attributes.
.bucss_power_result <- function(sample_size, size_term, ncp,
                                actual_power = NULL, design, sample_size_unit,
                                effect = NULL, assurance_ceiling = NULL,
                                total_n = NULL, df_effect = NULL, df_error = NULL,
                                inputs = list()) {
  stopifnot(size_term %in% .BUCSS_SIZE_TERMS)
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
  terms <- c(terms, "ncp_adjusted")
  values <- c(values, ncp)
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
  class(out) <- c("bucss_power", "data.frame")
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
# one place and only the presentation differs.
.relabel_bucss_result <- function(res, design, sample_size_unit, size_term,
                                  inputs, effect = NULL) {
  size <- res$value[res$term %in% .BUCSS_SIZE_TERMS][1]
  total_n <- res$value[res$term == "total_N"]
  if (length(total_n) == 0L) total_n <- NULL
  actual_power <- res$value[res$term == "actual_power"]
  if (length(actual_power) == 0L) actual_power <- NULL
  .bucss_power_result(
    sample_size = size,
    size_term = size_term,
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
# .tm_correlation the analogue for a Pearson correlation, whose statistic is
# not noncentral F at all once both variables are sampled (see below).
.tm_f <- function(stat, crit, df1, df2) {
  function(ncp) {
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
    power <- (1 - pt(crit, df = df, ncp = ncp)) + pt(-crit, df = df, ncp = ncp)
    area_above <- (1 - pt(stat, df = df, ncp = ncp)) + pt(-stat, df = df, ncp = ncp)
    (power - area_above) / power
  }
}

.tm_chisq <- function(stat, crit, df) {
  function(ncp) {
    power <- 1 - pchisq(crit, df = df, ncp = ncp)
    area_above <- 1 - pchisq(stat, df = df, ncp = ncp)
    (power - area_above) / power
  }
}

# Exact distribution of the test statistic for a Pearson correlation when both
# variables are sampled (the random-predictor frame a correlation study is
# actually run in), in closed form.
#
# Conditional on the sampled X, the statistic F = t^2 for testing rho = 0 is
# noncentral F(1, n - 2) with noncentrality lambda = c * SS_x/sigma_x^2, where
# c = rho^2/(1 - rho^2); and SS_x/sigma_x^2 is chi-square on n - 1 df. So the
# marginal distribution is a chi-square mixture of noncentral F's, NOT the
# noncentral F itself. Treating it as noncentral F (which is exact only when
# the predictor is fixed by design) overstates the power, because it drops the
# extra dispersion the sampled predictor contributes.
#
# The mixture has a closed form. The noncentral F CDF is a Poisson(lambda/2)
# mixture of central betas; mixing a Poisson rate over a gamma gives a negative
# binomial, and lambda/2 = c * W/2 with W/2 ~ Gamma((n-1)/2, scale = c). The
# Poisson weights therefore become negative binomial weights with size
# (n - 1)/2 and probability 1/(1 + c) = 1 - rho^2:
#
#   P(F <= q) = sum_k dnbinom(k, size = (n-1)/2, prob = 1 - rho^2) *
#               pbeta(q/(q + n - 2), 1/2 + k, (n-2)/2)
#
# This is exact rather than approximate, needs nothing beyond stats, and
# converges geometrically at rate rho^2. Numerical quadrature was the obvious
# alternative and is a trap: integrate() over (0, Inf) silently returns 0 once
# n is a few hundred, because the adaptive rule samples where a chi-square on
# n - 1 df has no mass.
#
# The series is truncated where the remaining negative binomial mass is below
# 1e-14; the omitted terms carry pbeta values smaller than the last retained
# one, so the truncation error is far below that mass. The term cap is pure
# defense: the count is governed by the noncentrality parameter, not by n, and
# a planned study near the target power carries lambda near 8 whatever its
# sample size.
.p_correlation_f <- function(q, n, rho2) {
  prob <- 1 - rho2
  if (prob <= 0) return(0)
  df2 <- n - 2
  size <- (n - 1) / 2
  k <- 0:min(qnbinom(1 - 1e-14, size = size, prob = prob), 1e6)
  sum(dnbinom(k, size = size, prob = prob) *
        pbeta(q / (q + df2), 0.5 + k, df2 / 2))
}

# TM for a Pearson correlation, parameterized (like every other planner) by the
# noncentrality parameter the adjusted effect implies at the prior study's
# sample size, lambda = c * N. Solving on that scale rather than on the rho^2
# scale matters: uniroot()'s tolerance is absolute, and a rho^2 root of order
# 1e-6, which is what a plan near the assurance ceiling implies, would be
# resolved to only a percent or two. On the lambda scale the root is of order 1
# to 100 and the shared solver needs no change.
.tm_correlation <- function(stat, crit, N) {
  function(ncp) {
    vapply(ncp, function(lambda) {
      rho2 <- lambda / (lambda + N)
      power <- 1 - .p_correlation_f(crit, N, rho2)
      area_above <- 1 - .p_correlation_f(stat, N, rho2)
      (power - area_above) / power
    }, numeric(1))
  }
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
# integer counter.
.smallest_n_for_power <- function(power_at, target, start = 2L) {
  lo <- as.numeric(start)
  if (power_at(lo) >= target) return(lo)
  hi <- lo * 2                       # invariant: power_at(lo) < target
  while (power_at(hi) < target) {
    lo <- hi
    hi <- hi * 2
  }
  while (hi - lo > 1) {              # invariant: power_at(hi) >= target
    mid <- floor((lo + hi) / 2)
    if (power_at(mid) < target) lo <- mid else hi <- mid
  }
  hi
}

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
# only remaining lever. call. = FALSE keeps the (uninformative) calling
# expression out of the message.
.stop_zero_ncp <- function(max_assurance) {
  # Floor to two decimals (with a tiny tolerance against floating-point
  # underflow) so the reported ceiling is itself a feasible value; rounding up
  # could name an assurance that still fails.
  ceiling_assurance <- floor(max_assurance * 100 + 1e-9) / 100
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

  stop(
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
    "Anderson, Kelley, and Maxwell (2017) and Anderson and Kelley (2024).",
    call. = FALSE
  )
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

#' Print a bias and uncertainty corrected sample size result
#'
#' Prints the aligned \code{term}/\code{value} table (the planned-study
#' results followed by the rows echoing the planning inputs), then factual
#' footer lines naming the design, the unit the sample size is counted in,
#' and the largest assurance the prior result can support. Whole numbers
#' print clean; other values print at \code{getOption("bucss.digits", 3)}
#' significant figures. Only the display is rounded; the stored values keep
#' full precision.
#'
#' @param x An object of class \code{bucss_power}.
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
#' @aliases bucss_power
#' @export
#' @keywords internal
print.bucss_power <- function(x, digits = getOption("bucss.digits", 3L), ...) {
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
  if (!is.null(unit)) cat("Sample size unit: ", unit, "\n", sep = "")
  if (!is.null(assurance_ceiling)) {
    cval <- floor(assurance_ceiling * 100 + 1e-9) / 100
    cat("Largest supportable assurance: ",
        sub("^0", "", sprintf("%.2f", cval)), "\n", sep = "")
  }
  invisible(x)
}
