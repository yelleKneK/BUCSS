# Output helpers that mirror DMAR's display conventions: a knit_print method
# so a bucss_power result renders as a kable in R Markdown with the same
# footer lines the console print shows, and planning_sentence(), the sentence
# an author pastes into a manuscript.

#' @exportS3Method knitr::knit_print
knit_print.bucss_power <- function(x, ...) {
  digits <- getOption("bucss.digits", 3L)
  shown <- data.frame(term = x$term,
                      value = .format_bucss_value(x$value, digits),
                      stringsAsFactors = FALSE)
  footers <- character(0)
  design <- attr(x, "design")
  effect <- attr(x, "effect")
  unit <- attr(x, "sample_size_unit")
  assurance_ceiling <- attr(x, "assurance_ceiling")
  if (!is.null(design)) {
    footers <- c(footers, paste0("Design: ", design,
                                 if (!is.null(effect)) paste0(" (", effect, ")")))
  }
  if (!is.null(unit)) footers <- c(footers, paste0("Sample size unit: ", unit))
  if (!is.null(assurance_ceiling)) {
    cval <- floor(assurance_ceiling * 100 + 1e-9) / 100
    footers <- c(footers, paste0("Largest supportable assurance: ",
                                 sub("^0", "", sprintf("%.2f", cval))))
  }
  out <- c(as.character(knitr::kable(shown, row.names = FALSE)), "",
           paste(footers, collapse = "; "), "")
  knitr::asis_output(paste(out, collapse = "\n"))
}

# Format a proportion as a percentage with no trailing ".0" (0.8 -> "80",
# 0.825 -> "82.5"), for the sentence helper.
.format_pct <- function(p) {
  sub("\\.?0+$", "", sprintf("%.1f", 100 * p))
}

#' The planning sentence for a manuscript
#'
#' Turns a \code{bucss_power} result into the sentence an author writes in a
#' manuscript's planning or method section: the necessary sample size in its
#' design's unit, the desired power, the assurance, and a clause noting the
#' publication bias and uncertainty correction. This mirrors the
#' \code{results_sentence()} helper in the sibling package \code{DMAR}.
#'
#' @param x An object of class \code{bucss_power}, as returned by any
#'   \code{ss_buc_*} function.
#'
#' @return A single character string.
#'
#' @export
#'
#' @examples
#' result <- ss_buc_independent_t(t_observed = 3, n = c(50, 55),
#'   assurance = .90)
#' planning_sentence(result)
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu})
planning_sentence <- function(x) {
  if (!inherits(x, "bucss_power")) {
    stop("'x' must be a bucss_power result from an ss_buc_* function.",
         call. = FALSE)
  }
  size <- x$value[x$term %in% .BUCSS_SIZE_TERMS][1]
  total_n <- x$value[x$term == "total_N"]
  desired_power <- x$value[x$term == "desired_power"]
  assurance <- x$value[x$term == "assurance"]
  unit <- attr(x, "sample_size_unit")

  size_phrase <- if (identical(unit, "total")) {
    paste0("a total sample size of ", size)
  } else if (identical(unit, "number of pairs")) {
    paste0(size, " pairs")
  } else {
    paste0(size, " participants ", unit,
           if (length(total_n) == 1L) paste0(" (total N = ", total_n, ")"))
  }
  paste0("A planned study with ", size_phrase, " attains ",
         .format_pct(desired_power), "% power with ",
         .format_pct(assurance),
         "% assurance, after correcting the prior study's result for ",
         "publication bias and uncertainty.")
}
