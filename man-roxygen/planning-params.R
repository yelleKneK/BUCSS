#' @param alpha_prior The significance level \eqn{\alpha} a prior result had
#'   to reach to be published (more generally, to reach you): the correction
#'   assumes you see only results with \emph{p} below \code{alpha_prior}. A
#'   value of 1 means no publication filter: every \emph{p} value is at most
#'   1, so any result would have been published whatever its \emph{p} value,
#'   there is no publication bias to correct, and only the correction for
#'   uncertainty (set by \code{assurance}) remains, as for a pilot study. The
#'   planners compute with .999 in place of 1, as BUCSS 1.x did, so plans
#'   match earlier versions; that filter is so lenient (it drops only results
#'   with \emph{p} above .999) that it changes the plan negligibly, and the
#'   result echoes the 1 you entered.
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
#'   \code{assurance = .80} the plan corresponds to safeguard power, which
#'   plans for the lower limit of a two-sided 60 percent confidence interval
#'   around the prior effect size rather than for the effect size itself (see
#'   \code{\link{BUCSS-package}} for how, and for where the two differ).
#'   Enter it as a proportion in (0, 1) or as a percentage greater
#'   than 1 (for example, 80 is read as .80); a value of exactly 1 is read as
#'   1 percent. A percentage is echoed in the result's planning-input rows as
#'   the coerced proportion (80 is echoed as .80); the same applies to
#'   \code{desired_power}.
#' @param desired_power Desired level of statistical power for the planned
#'   study. Enter it as a proportion in (0, 1) or as a percentage greater than
#'   1; a value of exactly 1 is read as 1 percent.
