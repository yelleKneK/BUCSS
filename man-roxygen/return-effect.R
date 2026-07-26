#' @return An object of class \code{bucss_power}: a \code{data.frame} with a
#'   character \code{term} column and a numeric \code{value} column, in the
#'   same shape as the sibling package \code{DMAR}'s planners. The design
#'   results come first: the necessary <%= size_phrase %> for the planned
#'   study, named for its unit (\code{necessary_n_per_group},
#'   \code{necessary_n_per_cell}, or \code{necessary_N}); for the per-group
#'   and per-cell designs, \code{total_N}, the implied total sample size;
#'   \code{actual_power}, the conservative achieved power at the returned
#'   sample size (see Details); and \code{ncp_adjusted}, the publication bias
#'   and uncertainty adjusted prior study noncentrality parameter. Rows
#'   echoing the planning inputs follow, so the assumptions travel with the
#'   result. The design, the effect tested, the unit the sample size is
#'   counted in, the largest assurance this prior can support (the assurance
#'   ceiling), and the planned test's degrees of freedom (\code{df_effect}
#'   and \code{df_error}) travel on attributes; they are shown by
#'   \code{print.bucss_power}, and \code{tidy()} and \code{glance()} give the
#'   compact and wide views.
