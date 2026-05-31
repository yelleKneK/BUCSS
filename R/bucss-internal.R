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
