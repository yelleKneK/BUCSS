#' @return An object of class \code{bucss_power}: a \code{data.frame} with a
#'   character \code{term} column and a numeric \code{value} column whose two
#'   rows are \code{necessary_sample_size} (the suggested <%= size_phrase %> for
#'   the planned study) and \code{ncp_adjusted} (the publication bias and
#'   uncertainty adjusted prior study noncentrality parameter). The design, the
#'   effect tested, the unit the sample size is counted in, and the planning
#'   inputs travel on attributes and are shown by \code{print.bucss_power}.
