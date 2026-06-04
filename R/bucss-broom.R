# tidy() and glance() methods for bucss_power results, registered with the
# generics package (the same mechanism the DMAR package uses). They give the
# programmer-friendly views of the tidy result: tidy() pivots the value rows to
# a one-row wide data.frame (so a single quantity is pulled out by name), and
# glance() returns a one-row model-level summary.

#' Tidy and summarize a bias and uncertainty corrected sample size result
#'
#' \code{tidy()} returns the planned-study quantities as a one-row
#' \code{data.frame} with one column per quantity (the "wide" view), which is
#' the convenient way to pull a single value out by name or to plot or join the
#' result. \code{glance()} returns a one-row model-level summary (the design,
#' the effect tested, the headline sample sizes, and the assurance ceiling the
#' prior result supports).
#'
#' These methods come from the \pkg{generics} package and need no extra setup;
#' the verbs are re-exported by \pkg{BUCSS}, so \code{tidy(result)} and
#' \code{glance(result)} work as soon as the package is loaded. The result is
#' also an ordinary \code{data.frame}, so the long view and the stored full
#' precision values remain available directly, for example
#' \code{result$value[result$term == "necessary_sample_size"]}.
#'
#' @param x An object of class \code{bucss_power} returned by an \code{ss_buc_*}
#'   function.
#' @param ... Ignored.
#'
#' @return A one-row \code{data.frame}. For \code{tidy()}, one column per
#'   planned-study quantity (\code{necessary_sample_size}, \code{total_N} where
#'   applicable, \code{ncp_adjusted}, and \code{assurance_ceiling}). For
#'   \code{glance()}, a one-row summary with the design, effect, headline sample
#'   sizes, unit, and assurance ceiling.
#'
#' @examples
#' result <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
#'   levels_B = 3, effect = "factor_B")
#'
#' # One-row wide view: pull a quantity out by name, no row indexing needed.
#' tidy(result)
#' tidy(result)$necessary_sample_size
#' tidy(result)$total_N
#'
#' # One-row model-level summary.
#' glance(result)
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
  wide <- as.data.frame(as.list(stats::setNames(x$value, x$term)),
                        stringsAsFactors = FALSE)
  ceiling <- attr(x, "assurance_ceiling")
  if (!is.null(ceiling)) wide$assurance_ceiling <- ceiling
  wide
}

#' @rdname bucss_tidiers
#' @exportS3Method generics::glance
glance.bucss_power <- function(x, ...) {
  effect <- attr(x, "effect")
  if (is.null(effect)) effect <- NA_character_
  ceiling <- attr(x, "assurance_ceiling")
  if (is.null(ceiling)) ceiling <- NA_real_
  total_N <- if ("total_N" %in% x$term) x$value[x$term == "total_N"] else NA_real_
  data.frame(
    design                = attr(x, "design"),
    effect                = effect,
    necessary_sample_size = x$value[x$term == "necessary_sample_size"],
    total_N               = total_N,
    sample_size_unit      = attr(x, "sample_size_unit"),
    assurance_ceiling     = ceiling,
    stringsAsFactors      = FALSE
  )
}
