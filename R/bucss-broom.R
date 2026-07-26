# tidy() and glance() methods for bucss_power results, registered with the
# generics package and shaped exactly like the DMAR package's dmar_ss_power
# tidiers: tidy() is the compact estimate view (term/estimate/power), and
# glance() is the one-row wide view carrying the estimate together with every
# echoed planning input and the design metadata as columns.

#' Tidy and summarize a bias and uncertainty corrected sample size result
#'
#' \code{tidy()} returns the compact estimate view, mirroring the sibling
#' package \code{DMAR}'s sample size planners: one row with \code{term}
#' (\code{"sample_size"}), \code{estimate} (the necessary sample size in the
#' design's unit), and \code{power} (the conservative achieved power at that
#' sample size; see the planner's help page). \code{glance()} returns the
#' one-row wide view: the estimate plus one column per row of the result
#' (\code{total_N} where applicable, \code{ncp_adjusted}, and every echoed
#' planning input) together with the design, the effect tested, the unit, the
#' planned test's degrees of freedom, and the assurance ceiling.
#'
#' These methods come from the \pkg{generics} package and need no extra setup;
#' the verbs are re-exported by \pkg{BUCSS}, so \code{tidy(result)} and
#' \code{glance(result)} work as soon as the package is loaded. The result is
#' also an ordinary \code{data.frame}, so the long view and the stored full
#' precision values remain available directly, for example
#' \code{result$value[result$term == "ncp_adjusted"]}.
#'
#' @param x An object of class \code{bucss_power} returned by an \code{ss_buc_*}
#'   function.
#' @param ... Ignored.
#'
#' @return A one-row \code{data.frame}. For \code{tidy()}, the columns
#'   \code{term}, \code{estimate}, and \code{power}. For \code{glance()}, the
#'   wide summary described above.
#'
#' @examples
#' result <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
#'   levels_B = 3, effect = "factor_B")
#'
#' # Compact estimate view, as in DMAR: the sample size and its power.
#' tidy(result)
#' tidy(result)$estimate
#'
#' # One-row wide view: every quantity and echoed input as a column.
#' glance(result)
#' glance(result)$total_N
#'
#' @name bucss_tidiers
#' @author Ken Kelley (\email{kkelley@@nd.edu})
NULL

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @rdname bucss_tidiers
#' @exportS3Method generics::tidy
tidy.bucss_power <- function(x, ...) {
  power <- x$value[x$term == "actual_power"]
  if (length(power) == 0L) power <- NA_real_
  data.frame(
    term     = "sample_size",
    estimate = x$value[x$term %in% .BUCSS_SIZE_TERMS][1],
    power    = power,
    stringsAsFactors = FALSE
  )
}

#' @rdname bucss_tidiers
#' @exportS3Method generics::glance
glance.bucss_power <- function(x, ...) {
  wide <- as.data.frame(as.list(stats::setNames(x$value, x$term)),
                        stringsAsFactors = FALSE)
  effect <- attr(x, "effect")
  wide$design <- attr(x, "design")
  if (!is.null(effect)) wide$effect <- effect
  wide$sample_size_unit <- attr(x, "sample_size_unit")
  df_effect <- attr(x, "df_effect")
  df_error <- attr(x, "df_error")
  if (!is.null(df_effect)) wide$df_effect <- df_effect
  if (!is.null(df_error)) wide$df_error <- df_error
  ceiling <- attr(x, "assurance_ceiling")
  if (!is.null(ceiling)) wide$assurance_ceiling <- ceiling
  wide
}
