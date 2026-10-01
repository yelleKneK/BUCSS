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
#'   planning effect size, <%= effect_phrase %>.
#'
#'   A planning effect size (the row whose name ends in \code{_planning}) is
#'   an adjusted effect size: the value, in the design's usual effect size
#'   metric, used to plan the sample size given the situation at hand (the
#'   prior study's result, the publication bias assumed through
#'   \code{alpha_prior}, and the uncertainty allowed for through
#'   \code{assurance}). It is \code{ncp_adjusted} restated in that metric
#'   and is what the plan is built on, not an estimate of the population
#'   effect to report: above an assurance of .5 it is deliberately
#'   conservative (a lower confidence limit), and at .5 it is the median
#'   unbiased estimate. Because it restates \code{ncp_adjusted},
#'   \code{actual_power} is the ordinary power of a study of the returned
#'   sample size when the population effect is the planning effect size.
#'
#'   Rows echoing the planning inputs follow the design results, so the
#'   assumptions travel with the result. The design, the effect tested, the unit the sample size is
#'   counted in, the largest assurance this prior can support (the assurance
#'   ceiling), and the planned test's degrees of freedom (\code{df_effect}
#'   and \code{df_error}) travel on attributes. \code{print.bucss_plan} shows
#'   the design, the effect, the unit, and the ceiling beneath the table;
#'   \code{tidy()} and \code{glance()} give the compact and wide views, the
#'   latter with every attribute as a column.
