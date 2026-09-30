#' @return An object of class \code{bucss_plan}: a \code{data.frame} with a
#'   character \code{term} column and a numeric \code{value} column, in the
#'   same shape as the sibling package \code{DMAR}'s planners. The design
#'   results come first: the necessary <%= size_phrase %> for the planned
#'   study, named for its unit (\code{necessary_n_per_group},
#'   \code{necessary_n_per_cell}, or \code{necessary_N}); for the per-group
#'   and per-cell designs, \code{total_N}, the implied total sample size;
#'   \code{actual_power}, the power of the plan returned, that is, of the
#'   returned sample size at the returned noncentrality parameter (see
#'   Details); \code{ncp_adjusted}, the publication bias and uncertainty
#'   adjusted prior study noncentrality parameter, <%= ncp_scale %>; and the
#'   planning effect size, the effect size the plan is built on, which is
#'   \code{ncp_adjusted} restated in the design's usual metric:
#'   <%= effect_phrase %>. Because it restates \code{ncp_adjusted},
#'   \code{actual_power} is the ordinary power of a study of the returned
#'   sample size when the population effect is the planning effect size. Rows
#'   echoing the planning inputs follow, so the assumptions travel with the
#'   result. The design, the effect tested, the unit the sample size is
#'   counted in, the largest assurance this prior can support (the assurance
#'   ceiling), and the planned test's degrees of freedom (\code{df_effect}
#'   and \code{df_error}) travel on attributes. \code{print.bucss_plan} shows
#'   the design, the effect, the unit, and the ceiling beneath the table;
#'   \code{tidy()} and \code{glance()} give the compact and wide views, the
#'   latter with every attribute as a column.
