# Internal constructor and print method shared by every ss_buc_* function.

# Build the tidy result object returned by every ss_buc_* function.
#
# The two computed quantities (the necessary planned-study sample size and the
# bias and uncertainty adjusted prior-study noncentrality parameter) live in a
# numeric `value` column so the object behaves like an ordinary data.frame for
# downstream arithmetic. Everything non-numeric (the human design label, the
# unit the sample size is counted in, the effect tested, and the planning
# inputs) travels on attributes, never as a string in `value`.
.bucss_power_result <- function(sample_size, ncp, design, sample_size_unit,
                                effect = NULL, inputs = list()) {
  out <- data.frame(
    term = c("necessary_sample_size", "ncp_adjusted"),
    value = c(sample_size, ncp),
    stringsAsFactors = FALSE
  )
  inputs <- inputs[!vapply(inputs, is.null, logical(1))]
  attr(out, "design") <- design
  attr(out, "sample_size_unit") <- sample_size_unit
  attr(out, "effect") <- effect
  attr(out, "inputs") <- inputs
  class(out) <- c("bucss_power", "data.frame")
  out
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
  inputs <- attr(x, "inputs")

  cat("Bias and uncertainty corrected sample size\n")
  if (!is.null(design)) cat("Design:", design, "\n")
  if (!is.null(effect)) cat("Effect of interest:", effect, "\n")
  cat("\n")
  cat("Necessary sample size (", unit, "): ", sample_size, "\n", sep = "")
  cat("Adjusted noncentrality parameter: ", format(ncp), "\n", sep = "")
  if (length(inputs) > 0L) {
    cat("\nPlanning inputs:\n")
    nms <- names(inputs)
    for (i in seq_along(inputs)) {
      cat("  ", nms[i], " = ",
          paste(format(inputs[[i]]), collapse = ", "), "\n", sep = "")
    }
  }
  invisible(x)
}
