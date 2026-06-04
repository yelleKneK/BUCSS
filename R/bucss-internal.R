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
                                total_n = NULL, inputs = list()) {
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
  attr(out, "inputs") <- inputs
  class(out) <- c("bucss_power", "data.frame")
  out
}

# Validate the planning inputs shared by every ss_buc_* function and apply the
# documented coercions, so the checks stay identical across all planners. The
# three coercions are load-bearing: 'alpha_prior == 1' becomes .999 to model
# "no publication bias" without a degenerate truncation, and 'assurance' or
# 'power' entered as a percentage (>= 1) is divided by 100 so that, for example,
# 80 means .80. Returns the (possibly coerced) values; 'alpha_prior_input' is
# the value the user supplied, echoed back in the result's planning inputs. The
# grid resolution 'step' is not a user argument; it is read here from
# getOption("bucss.step", .001) so the planner signatures stay focused on the
# study design, and is returned for the planner's NCP grid.
# call. = FALSE keeps this internal helper out of the error the user sees.
.validate_planning_inputs <- function(alpha_prior, alpha_planned, assurance,
                                      power) {
  if (alpha_prior > 1 | alpha_prior <= 0) stop("There is a problem with 'alpha_prior' of the prior study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).", call. = FALSE)
  alpha_prior_input <- alpha_prior
  if (alpha_prior == 1) alpha_prior <- .999

  if (alpha_planned >= 1 | alpha_planned <= 0) stop("There is a problem with 'alpha_planned' of the planned study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).", call. = FALSE)

  if (assurance >= 1) assurance <- assurance / 100
  if (assurance < 0 | assurance > 1) stop("There is a problem with 'assurance' (i.e., the proportion of times statistical power is at or above the desired value), please specify as a value between 0 and 1 (the default is .80).", call. = FALSE)
  if (assurance < .5) warning("The 'assurance' you have entered is < .5, which implies you will have under a 50% chance at achieving your desired level of power.", call. = FALSE)

  if (power >= 1) power <- power / 100
  if (power < 0 | power > 1) stop("There is a problem with 'power' (i.e., desired statistical power), please specify as a value between 0 and 1 (the default is .80).", call. = FALSE)

  step <- getOption("bucss.step", 0.001)
  if (length(step) != 1L || !is.finite(step) || step <= 0 || step >= 1) stop("The 'bucss.step' option (set with options(bucss.step = ...)) must be a single number greater than 0 and less than 1; the default is .001.", call. = FALSE)

  list(alpha_prior = alpha_prior, alpha_prior_input = alpha_prior_input,
       assurance = assurance, power = power, step = step)
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
