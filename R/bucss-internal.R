# Internal constructor and print method shared by every ss_buc_* function.

# Build the tidy result object returned by every ss_buc_* function.
#
# The planned-study quantities (the necessary sample size, the implied total
# sample size for per-group and per-cell designs, and the bias and uncertainty
# adjusted prior-study noncentrality parameter) live as rows of a numeric
# `value` column so the object behaves like an ordinary data.frame for
# downstream arithmetic; tidy() pivots them to a one-row wide view. The
# `necessary_sample_size` and `ncp_adjusted` term names are fixed (the
# regression oracles read them). The `total_N` row is present only when a
# `total_n` is supplied (the per-group and per-cell designs). Everything else
# (the human design label, the unit, the effect tested, the assurance ceiling
# the prior supports, and the planning inputs) travels on attributes.
.bucss_power_result <- function(sample_size, ncp, design, sample_size_unit,
                                effect = NULL, assurance_ceiling = NULL,
                                total_n = NULL, df_effect = NULL, df_error = NULL,
                                inputs = list()) {
  if (is.null(total_n)) {
    out <- data.frame(term = c("necessary_sample_size", "ncp_adjusted"),
                      value = c(sample_size, ncp), stringsAsFactors = FALSE)
  } else {
    out <- data.frame(term = c("necessary_sample_size", "total_N", "ncp_adjusted"),
                      value = c(sample_size, total_n, ncp), stringsAsFactors = FALSE)
  }
  inputs <- inputs[!vapply(inputs, is.null, logical(1))]
  attr(out, "design") <- design
  attr(out, "sample_size_unit") <- sample_size_unit
  attr(out, "effect") <- effect
  attr(out, "assurance_ceiling") <- assurance_ceiling
  attr(out, "df_effect") <- df_effect
  attr(out, "df_error") <- df_error
  attr(out, "inputs") <- inputs
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
    stop("'", name, "' must be a single finite number.", call. = FALSE)
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
                                      power) {
  .check_scalar_finite(alpha_prior, "alpha_prior")
  .check_scalar_finite(alpha_planned, "alpha_planned")
  .check_scalar_finite(assurance, "assurance")
  .check_scalar_finite(power, "power")
  if (alpha_prior > 1 | alpha_prior <= 0) stop("There is a problem with 'alpha_prior' of the prior study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).", call. = FALSE)
  alpha_prior_input <- alpha_prior
  if (alpha_prior == 1) alpha_prior <- .999

  if (alpha_planned >= 1 | alpha_planned <= 0) stop("There is a problem with 'alpha_planned' of the planned study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).", call. = FALSE)

  if (assurance >= 1) assurance <- assurance / 100
  if (assurance <= 0 | assurance >= 1) stop("There is a problem with 'assurance' (i.e., the proportion of times statistical power is at or above the desired value), please specify as a value strictly between 0 and 1 (the default is .80). An 'assurance' of exactly 0 or 1 does not define a plannable study; note that 100 entered as a percentage means 1.", call. = FALSE)
  if (assurance < .5) warning("The 'assurance' you have entered is < .5, which implies you will have under a 50% chance at achieving your desired level of power.", call. = FALSE)

  if (power >= 1) power <- power / 100
  if (power <= 0 | power >= 1) stop("There is a problem with 'power' (i.e., desired statistical power), please specify as a value strictly between 0 and 1 (the default is .80). A planned power of exactly 0 or 1 is not attainable; note that 100 entered as a percentage means 1.", call. = FALSE)

  list(alpha_prior = alpha_prior, alpha_prior_input = alpha_prior_input,
       assurance = assurance, power = power)
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

# Build the truncated-likelihood mean TM as a function of the noncentrality
# parameter for the two distributional shapes the planners use. TM(ncp) is the
# truncated cumulative density evaluated at the observed statistic: the area of
# the (publication-)truncated noncentral distribution between the critical value
# and the observed statistic, divided by the power (the untruncated area beyond
# the critical value). It is monotone decreasing in 'ncp', from its ceiling
# TM(0) = 1 - p_prior / alpha_prior down to 0. The returned closure is
# vectorized over 'ncp' (pf()/pt() are), so it serves both uniroot() (scalar
# calls) and the grid fallback (vector call) in .solve_ncp_assurance().
#
# Upper safety bound on the root-finding bracket in .solve_ncp_assurance(). It
# is far beyond any noncentrality parameter a real prior study implies; reaching
# it signals a likely input error (or an essentially deterministic prior) rather
# than a genuine solution, and the solver warns and returns it instead of
# expanding without bound.
.bucss_ncp_upper_cap <- 1e7

# .tm_f covers the one-sided noncentral F designs (the ANOVA and regression
# planners; the single regression coefficient passes the squared t as 'stat'
# against an F critical value). .tm_t covers the two-sided noncentral t designs
# (the independent and paired t tests), summing both tails exactly as the 1.x
# code did.
.tm_f <- function(stat, crit, df1, df2) {
  function(ncp) {
    power <- 1 - pf(crit, df1 = df1, df2 = df2, ncp = ncp)
    area_above <- 1 - pf(stat, df1 = df1, df2 = df2, ncp = ncp)
    (power - area_above) / power
  }
}

.tm_t <- function(stat, crit, df) {
  function(ncp) {
    power <- (1 - pt(crit, df = df, ncp = ncp)) + pt(-crit, df = df, ncp = ncp)
    area_above <- (1 - pt(stat, df = df, ncp = ncp)) + pt(-stat, df = df, ncp = ncp)
    (power - area_above) / power
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
# planned-n search through this helper so the incremental search, and its
# off-by-one, lives in one place rather than being re-derived per design; each
# planner supplies a 'power_at' closure carrying its own degrees-of-freedom and
# noncentrality scaling.
.smallest_n_for_power <- function(power_at, target, start = 2L) {
  n_rep <- start
  while (power_at(n_rep) < target) n_rep <- n_rep + 1L
  n_rep
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
      "the ceiling. "
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

#' Print a bias and uncertainty corrected sample size result
#'
#' Formats the object returned by the \code{ss_buc_*} functions for humans:
#' the design, the necessary planned-study sample size with the unit it is
#' counted in, the bias and uncertainty adjusted prior-study noncentrality
#' parameter, and the planning inputs. The object itself is an ordinary
#' \code{data.frame} with \code{term} and numeric \code{value} columns, so the
#' two quantities remain available as \code{x$value[x$term == "..."]}.
#'
#' @param x An object of class \code{bucss_power}.
#' @param ... Further arguments, ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @aliases bucss_power
#' @export
#' @keywords internal
print.bucss_power <- function(x, ...) {
  sample_size <- x$value[x$term == "necessary_sample_size"]
  ncp <- x$value[x$term == "ncp_adjusted"]
  design <- attr(x, "design")
  unit <- attr(x, "sample_size_unit")
  effect <- attr(x, "effect")
  total_n <- x$value[x$term == "total_N"]
  if (length(total_n) == 0L) total_n <- NULL
  assurance_ceiling <- attr(x, "assurance_ceiling")
  inputs <- attr(x, "inputs")

  cat("Bias and uncertainty corrected sample size\n")
  if (!is.null(design)) cat("Design:", design, "\n")
  if (!is.null(effect)) cat("Effect of interest:", effect, "\n")
  cat("\n")
  if (!is.null(total_n)) {
    cat("Necessary sample size (", unit, "): ", sample_size,
        "  (total N = ", total_n, ")\n", sep = "")
  } else {
    cat("Necessary sample size (", unit, "): ", sample_size, "\n", sep = "")
  }
  cat("Adjusted noncentrality parameter: ",
      format(signif(ncp, getOption("bucss.digits", 3L))), "\n", sep = "")
  if (!is.null(assurance_ceiling)) {
    cval <- floor(assurance_ceiling * 100 + 1e-9) / 100
    cat("Assurance this prior can support (ceiling): about ",
        sub("^0", "", sprintf("%.2f", cval)), "\n", sep = "")
  }
  if (length(inputs) > 0L) {
    cat("\nPlanning inputs:\n")
    nms <- names(inputs)
    for (i in seq_along(inputs)) {
      cat("  ", nms[i], " = ",
          paste(format(inputs[[i]]), collapse = ", "), "\n", sep = "")
    }
  }
  cat("\nUse tidy() for these quantities as a one-row data frame, ",
      "or glance() for a summary.\n", sep = "")
  invisible(x)
}
