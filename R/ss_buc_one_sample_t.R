#' Necessary Sample Size to Reach Desired Power for a One-Sample t Test Using a
#' Publication Bias and Uncertainty Correction Procedure
#'
#' @description \code{ss_buc_one_sample_t} returns the necessary total sample
#'   size to achieve a desired level of statistical power for a planned study
#'   using a one-sample \emph{t} test, based on information obtained from a
#'   previous study. The effect from the previous study can be corrected for
#'   publication bias and/or uncertainty to provide a sample size that will
#'   achieve more accurate statistical power for a planned study, when compared
#'   to approaches that use a sample effect size at face value or rely on sample
#'   size only. The bias and uncertainty adjusted previous study noncentrality
#'   parameter is also returned, together with the planning effect size, the
#'   effect size the plan is built on, as a standardized mean.
#'
#' @details The one-sample \emph{t} test and the dependent (paired) \emph{t}
#'   test are the same test applied to a single column of scores: the statistic
#'   has \emph{N} - 1 degrees of freedom and a noncentrality parameter of
#'   \eqn{d\sqrt{N}}{d*sqrt(N)}. \code{ss_buc_one_sample_t} therefore performs
#'   exactly the computation \code{\link{ss_buc_paired_t}} performs, and
#'   reports it with a
#'   one-sample design label and a total sample size rather than a number of
#'   pairs. See \code{\link{ss_buc_paired_t}} for the full description of the
#'   method, the meaning of \code{assurance}, the role of \code{alpha_prior},
#'   and what happens when the corrected noncentrality parameter is zero.
#'
#'   The returned \code{actual_power} is the statistical power of the plan
#'   this function reports: the power of a study of the returned sample size
#'   when the effect is the returned bias and uncertainty adjusted
#'   noncentrality parameter. Because a sample size must be a whole number,
#'   \code{actual_power} is never below the \code{desired_power} that was
#'   requested and is usually a little above it.
#'
#'   The observed \emph{t} may be entered with either sign; the publication
#'   rule the correction assumes is two-sided, so only the magnitude enters
#'   the computation.
#'
#'   If you are working from a standardized mean (Cohen's \eqn{d}) rather than
#'   a \emph{t} statistic, convert it before calling:
#'   \eqn{t = d\sqrt{N}}{t = d*sqrt(N)}.
#'
#' @param t_observed Observed \emph{t} value from a previous study used to plan
#'   sample size for a planned study. Either sign is accepted: the publication
#'   rule is two-sided, so only the magnitude enters the correction.
#' @param N Total sample size of the previous study.
#' @template planning-params
#'
#' @templateVar size_phrase total sample size
#' @templateVar ncp_scale on the scale of the \emph{t} statistic (the noncentrality parameter of the prior study's noncentral \emph{t})
#' @templateVar effect_phrase \code{std_mean_planning}, the standardized mean (the distance of the mean from its null value in standard deviations, a magnitude), \eqn{\lambda/\sqrt{N}}{lambda/sqrt(N)}, where \eqn{\lambda}{lambda} is \code{ncp_adjusted} and \emph{N} the prior sample size; it is the \code{dz_planning} of \code{\link{ss_buc_paired_t}} under the one-sample name
#' @template return
#'
#' @export
#'
#' @examples
#' result <- ss_buc_one_sample_t(t_observed = 3, N = 40, alpha_prior = .05,
#'   alpha_planned = .05, assurance = .80, desired_power = .80)
#' result
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .90):
#' try(ss_buc_one_sample_t(t_observed = 3, N = 40, assurance = .95))
#'
#' @seealso \code{\link{ss_buc_paired_t}}, which performs the identical
#'   computation for a paired design.
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
ss_buc_one_sample_t <- function(t_observed, N, alpha_prior = .05,
                                alpha_planned = .05, assurance = .80,
                                desired_power = .80) {
  # The mathematics is the paired t planner's, so it lives in one place; only
  # the design label and the unit differ. A missing N is caught here so the
  # message speaks of participants, not the paired planner's pairs.
  if (missing(N)) stop("You need to specify 'N', which is the total sample size (the number of participants) of the previous study.", call. = FALSE)
  res <- ss_buc_paired_t(t_observed = t_observed, N = N,
                         alpha_prior = alpha_prior,
                         alpha_planned = alpha_planned,
                         assurance = assurance, desired_power = desired_power)
  .relabel_bucss_result(
    res,
    design = "One-sample t test",
    sample_size_unit = "total",
    size_term = "necessary_N",
    # The paired planner's standardized mean of the differences is, for one
    # column of scores, the standardized mean.
    effect_size_term = "std_mean_planning",
    inputs = list(t_observed = t_observed, N = N,
                  alpha_prior = res$value[res$term == "alpha_prior"],
                  alpha_planned = res$value[res$term == "alpha_planned"],
                  assurance = res$value[res$term == "assurance"],
                  desired_power = res$value[res$term == "desired_power"])
  )
}
