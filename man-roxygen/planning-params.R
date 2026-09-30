#' @param alpha_prior Alpha level \eqn{\alpha} for the previous study or the
#'   assumed statistical significance necessary for publishing in the field; to
#'   assume no publication bias, a value of 1 can be entered.
#' @param alpha_planned Alpha level \eqn{\alpha} assumed for the planned study.
#' @param assurance Desired level of assurance: the long run proportion of
#'   planned studies whose statistical power at the true effect reaches or
#'   surpasses \code{desired_power}, when each is planned from a prior study
#'   published under the rule \code{alpha_prior} describes. The adjusted
#'   noncentrality parameter is the lower limit of a one-sided confidence
#'   interval at this level, computed from the publication-truncated
#'   distribution of the prior statistic, so \code{assurance = .5} plans from
#'   a median unbiased estimate (correcting for publication bias only), a
#'   larger value adds the correction for uncertainty, and a value below .5 is
#'   not recommended. The term follows Kelley and Maxwell (2003), who used
#'   assurance in sample size planning for the probability that a confidence
#'   interval will be no wider than desired. With \code{alpha_prior = 1} and
#'   \code{assurance = .80} the plan corresponds to safeguard power, a special
#'   case (see \code{\link{BUCSS-package}} for how, and for where the two
#'   differ). Enter it as a proportion in (0, 1) or as a percentage greater
#'   than 1 (for example, 80 is read as .80); a value of exactly 1 is read as
#'   1 percent. A percentage is echoed in the result's planning-input rows as
#'   the coerced proportion (80 is echoed as .80); the same applies to
#'   \code{desired_power}.
#' @param desired_power Desired level of statistical power for the planned
#'   study. Enter it as a proportion in (0, 1) or as a percentage greater than
#'   1; a value of exactly 1 is read as 1 percent.
