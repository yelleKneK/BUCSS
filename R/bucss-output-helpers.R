# Output helpers that mirror DMAR's display conventions: a knit_print method
# so a bucss_plan result renders as a kable in R Markdown with the same
# footer lines the console print shows, and planning_sentence(), the sentence
# an author pastes into a manuscript.

#' @exportS3Method knitr::knit_print
knit_print.bucss_plan <- function(x, ...) {
  # As in print.bucss_plan(): an object that is no longer one plan renders as
  # the plain data frame it is, without one plan's footers.
  if (!.is_plan_shape(x) || !is.numeric(x[["value"]]))
    return(knitr::knit_print(structure(x, class = "data.frame"), ...))
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
  footers <- c(footers, .predictor_footer(x))
  if (!is.null(unit)) footers <- c(footers, paste0("Sample size unit: ", unit))
  if (!is.null(assurance_ceiling)) {
    cval <- .floor_assurance(assurance_ceiling)
    footers <- c(footers, paste0("Largest supportable assurance: ",
                                 sub("^0", "", sprintf("%.2f", cval))))
  }
  out <- c(as.character(knitr::kable(shown, row.names = FALSE)), "",
           paste(footers, collapse = "; "), "")
  knitr::asis_output(paste(out, collapse = "\n"))
}

# Format a proportion as a percentage with no trailing ".0" (0.8 -> "80",
# 0.825 -> "82.5"), for the sentence helper. A proportion below 1 never prints
# as 100: it gets as many decimals as it needs (0.9999 -> "99.99").
.format_pct <- function(p) {
  d <- 1L
  while (p < 1 && round(100 * p, d) >= 100 && d < 15L) d <- d + 1L
  sub("\\.?0+$", "", sprintf(paste0("%.", d, "f"), round(100 * p, d)))
}

# The sphericity parameter as the sentence writes it: three significant
# figures and no leading zero (.865), with as many more figures as it takes
# for a value below 1 never to print as 1 (.9996 -> ".9996"), the rule
# .format_pct() applies to proportions. Display only.
.format_epsilon <- function(e) {
  d <- 3L
  while (e < 1 && signif(e, d) >= 1 && d < 15L) d <- d + 1L
  sub("^0", "", format(signif(e, d), digits = d, scientific = FALSE))
}

# An alpha level as APA writes it: no leading zero, and two decimals when two
# suffice (.1 -> ".10", .05 -> ".05", .025 -> ".025").
.format_alpha <- function(a) {
  s <- if (abs(a * 100 - round(a * 100)) < 1e-9) sprintf("%.2f", a) else
    format(a, scientific = FALSE, trim = TRUE)
  sub("^0", "", s)
}

# A count for prose: whole, with thousands separators, never in scientific
# notation (1485 -> "1,485").
.format_count <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

#' The Planning Sentence for a Manuscript
#'
#' Turns a \code{bucss_plan} result into the sentence an author writes in a
#' manuscript's planning or method section: the necessary sample size in its
#' design's unit, the desired power, the assurance, and a clause stating which
#' corrections the plan made. The clause is built from the plan's own settings:
#' the publication bias correction is named, with its assumed threshold, only
#' when \code{alpha_prior} is below 1, and the uncertainty correction only when
#' \code{assurance} is above .5 (at .5 the plan uses the median estimate). A
#' within-subjects or split-plot plan made with \code{epsilon} below 1 adds the
#' departure from sphericity it allowed for. A second sentence names the
#' planning effect size, the effect size the plan is built on (the result's
#' row whose name ends in \code{_planning}), in the design's usual metric: a
#' standardized mean difference, Cohen's \emph{f}, a correlation, and so on,
#' as the planner's help page describes. It is shown at
#' \code{getOption("bucss.digits", 3)} significant figures, and without a
#' leading zero when the metric cannot exceed 1. This mirrors the
#' \code{results_sentence()} helper in the sibling package \code{DMAR}.
#'
#' @param x An object of class \code{bucss_plan}, as returned by any
#'   \code{ss_buc_*} function.
#'
#' @return A single character string holding the two sentences.
#'
#' @export
#'
#' @examples
#' result <- ss_buc_independent_t(t_observed = 3, n = c(50, 55),
#'   assurance = .90)
#' planning_sentence(result)
#'
#' # With no publication bias assumed (a pilot study, say), the sentence
#' # names the uncertainty correction alone.
#' planning_sentence(ss_buc_independent_t(t_observed = 3, n = 20,
#'   alpha_prior = 1))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu})
planning_sentence <- function(x) {
  if (!inherits(x, "bucss_plan")) {
    stop("'x' must be a bucss_plan result from an ss_buc_* function.",
         call. = FALSE)
  }
  .check_plan_shape(x)
  size <- x$value[x$term %in% .BUCSS_SIZE_TERMS]
  total_n <- x$value[x$term == "total_N"]
  desired_power <- x$value[x$term == "desired_power"]
  assurance <- x$value[x$term == "assurance"]
  alpha_prior <- x$value[x$term == "alpha_prior"]
  unit <- attr(x, "sample_size_unit")

  size_phrase <- if (identical(unit, "total")) {
    paste0("a total sample size of ", .format_count(size))
  } else if (identical(unit, "number of pairs")) {
    paste0(.format_count(size), " pairs")
  } else {
    paste0(.format_count(size), " participants ", unit,
           if (length(total_n) == 1L)
             paste0(" (total N = ", .format_count(total_n), ")"))
  }

  # Which corrections the plan made follows from its own settings:
  # alpha_prior = 1 assumes no publication bias, and an assurance of .5 or
  # less takes the median (or a less conservative) estimate, which is no
  # correction for uncertainty.
  bias <- length(alpha_prior) == 1L && alpha_prior < 1
  uncertainty <- assurance > .5
  threshold <- if (bias) .format_alpha(alpha_prior)
  clause <- if (bias && uncertainty) {
    paste0("after correcting the prior study's result for publication bias ",
           "(assuming its literature published results significant at the ",
           threshold, " level) and for uncertainty")
  } else if (bias) {
    paste0("after correcting the prior study's result for publication bias ",
           "(assuming its literature published results significant at the ",
           threshold, " level) but not for uncertainty")
  } else if (uncertainty) {
    paste0("after correcting the prior study's result for uncertainty, with no ",
           "publication bias assumed")
  } else {
    paste0("using the prior study's result with no correction for publication ",
           "bias or for uncertainty")
  }
  # A within-subjects plan made under a sphericity correction says so; one
  # that assumed sphericity (its epsilon row is 1) reads as it always has.
  epsilon <- x$value[x$term == "epsilon"]
  departure <- length(epsilon) == 1L && epsilon < 1
  sphericity <- if (departure)
    paste0(", allowing for a departure from sphericity (epsilon = ",
           .format_epsilon(epsilon), ")")
  # The second sentence names the planning effect size, the effect size the
  # plan is built on, in the design's usual metric. Under a sphericity
  # correction it stays on the sphericity scale, as ncp_adjusted does.
  effect_size_term <- x$term[x$term %in% .BUCSS_EFFECT_TERMS]
  effect_size <- .format_effect_size(x$value[x$term == effect_size_term],
                                     effect_size_term)
  paste0("A planned study with ", size_phrase, " attains ",
         .format_pct(desired_power), "% power with ",
         .format_pct(assurance), "% assurance, ", clause, sphericity, ". ",
         "The plan is built on a planning effect size of ", effect_size,
         ", expressed as ", .BUCSS_EFFECT_LABELS[[effect_size_term]],
         if (departure) " on the sphericity scale", ".")
}

# The planning effect size as the sentence writes it: at
# getOption("bucss.digits", 3) significant figures, never in scientific
# notation, and, for the metrics that cannot exceed 1 (.BUCSS_EFFECT_BOUNDED),
# without the leading zero, as APA style writes such statistics. Display
# only; the stored row keeps full precision.
.format_effect_size <- function(value, term) {
  shown <- format(signif(value, getOption("bucss.digits", 3L)),
                  scientific = FALSE)
  if (term %in% .BUCSS_EFFECT_BOUNDED) shown <- sub("^0\\.", ".", shown)
  shown
}
