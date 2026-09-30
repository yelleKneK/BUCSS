# tidy() and glance() methods for bucss_plan results, registered with the
# generics package and shaped like the DMAR package's planner tidiers: tidy()
# is the compact estimate view (term/estimate/power), and glance() is the
# one-row wide view, which begins with the same three columns and then carries
# every echoed planning input and the design metadata as columns.

#' Tidy and Summarize a Bias and Uncertainty Corrected Sample Size Result
#'
#' \code{tidy()} returns the compact estimate view, as the sibling package
#' \code{DMAR}'s sample size planners do: one row with \code{term}
#' (\code{"sample_size"}), \code{estimate} (the necessary sample size in the
#' design's unit, which the result's \code{sample_size_unit} attribute names),
#' and \code{power} (the power of the plan the result reports; see
#' \code{actual_power} on the planner's help page). \code{glance()} returns
#' the one-row wide view: the same three columns, then one column per row of
#' the result, namely the sample size row, \code{total_N} where applicable,
#' \code{actual_power}, \code{ncp_adjusted}, the planning effect size (the
#' row whose name ends in \code{_planning}, such as \code{smd_planning}), and
#' every echoed planning input, together with the design, the effect tested,
#' the unit, the planned test's degrees of freedom where the design has them,
#' and the assurance ceiling.
#'
#' These methods come from the \pkg{generics} package and need no extra setup;
#' the verbs are re-exported by \pkg{BUCSS}, so \code{tidy(result)} and
#' \code{glance(result)} work as soon as the package is loaded. The result is
#' also an ordinary \code{data.frame}, so the long view and the stored full
#' precision values remain available directly, for example
#' \code{result$value[result$term == "ncp_adjusted"]}.
#'
#' @param x An object of class \code{bucss_plan} returned by an \code{ss_buc_*}
#'   function.
#' @param ... Ignored.
#'
#' @return A one-row \code{data.frame}. For \code{tidy()}, the columns
#'   \code{term}, \code{estimate}, and \code{power}. For \code{glance()}, the
#'   wide summary described above. Both stop with an explanation if the object
#'   no longer has the shape of a single plan (for example, several plans
#'   stacked with \code{dplyr::bind_rows()}, or a plan with a row removed by
#'   \code{dplyr::filter()}; both keep the class).
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

#' @title Objects Exported From Other Packages
#' @importFrom generics tidy
#' @export
generics::tidy

#' @title Objects Exported From Other Packages
#' @importFrom generics glance
#' @export
generics::glance

#' @rdname bucss_tidiers
#' @exportS3Method generics::tidy
tidy.bucss_plan <- function(x, ...) {
  .check_plan_shape(x)
  data.frame(
    term     = "sample_size",
    estimate = x$value[x$term %in% .BUCSS_SIZE_TERMS],
    power    = x$value[x$term == "actual_power"],
    stringsAsFactors = FALSE
  )
}

#' @rdname bucss_tidiers
#' @exportS3Method generics::glance
glance.bucss_plan <- function(x, ...) {
  .check_plan_shape(x)
  wide <- cbind(tidy.bucss_plan(x),
                as.data.frame(as.list(stats::setNames(x$value, x$term)),
                              stringsAsFactors = FALSE))
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
