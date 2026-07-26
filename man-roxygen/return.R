#' @return An object of class \code{bucss_power}: a \code{data.frame} with a
#'   character \code{term} column and a numeric \code{value} column. The rows
#'   are \code{necessary_sample_size} (the suggested <%= size_phrase %> for the
#'   planned study), then, for the designs whose sample size is counted per
#'   group or per cell, \code{total_N} (the implied total sample size), and
#'   \code{ncp_adjusted} (the publication bias and uncertainty adjusted prior
#'   study noncentrality parameter). The design, the unit the sample size is
#'   counted in, the largest assurance this prior can support (the assurance
#'   ceiling), the planned test's degrees of freedom (\code{df_effect} and
#'   \code{df_error}), and the planning inputs travel on attributes; they are
#'   shown by \code{print.bucss_power}, and \code{tidy()} and \code{glance()}
#'   return them as columns.
