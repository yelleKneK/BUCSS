#' Necessary Sample Size to Reach Desired Power for a Pearson Correlation Using
#' a Publication Bias and Uncertainty Correction Procedure
#'
#' @description \code{ss_buc_correlation} returns the necessary total sample
#'   size to achieve a desired level of statistical power for a test of a
#'   Pearson correlation in a planned study, based on information obtained from
#'   a previous study. The effect from the previous study can be corrected for
#'   publication bias and/or uncertainty to provide a sample size that will
#'   achieve more accurate statistical power for a planned study, when compared
#'   to approaches that use a sample effect size at face value or rely on sample
#'   size only. The bias and uncertainty adjusted previous study noncentrality
#'   parameter is also returned, together with the planning effect size, the
#'   effect size the plan is built on, as a correlation.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss_buc_correlation} uses the observed
#'   correlation and sample size from a previous study to correct the
#'   noncentrality parameter associated with the effect of interest for
#'   publication bias and/or uncertainty. This new estimated noncentrality
#'   parameter is then used to calculate the necessary total sample size to
#'   achieve the desired level of power in the planned study.
#'
#'   The observed correlation \emph{r} with sample size \emph{N} gives
#'   \eqn{t = r\sqrt{(N - 2)/(1 - r^2)}}{t = r*sqrt((N - 2)/(1 - r^2))} on
#'   \emph{N} - 2 degrees of freedom,
#'   and the planner works from that \emph{t}. Supplying \code{t_observed}
#'   directly instead of \code{r_observed} is equivalent.
#'
#'   The approach uses a likelihood function of a truncated distribution, where
#'   the truncation occurs due to small effect sizes being unobserved due to
#'   publication bias. The numerator of the likelihood function is the density
#'   of the test statistic; the denominator is the power of the test, which
#'   serves to truncate the distribution. (See Taylor & Muller, 1996, Equation
#'   2.1; Anderson, Kelley, & Maxwell, 2017, Equation 3; and Anderson &
#'   Maxwell, 2017, for more details.) Which density belongs in that numerator
#'   is described next; it is the density the regression planners use with
#'   \code{random_predictors = TRUE}.
#'
#'   \strong{Both variables are sampled, and that changes the distribution.}
#'   The test of a correlation is arithmetically the test of the slope in a
#'   simple regression, so it is tempting to plan it as one. That is exact only
#'   when the predictor is fixed by design, as it is in an experiment where the
#'   researcher sets the values of \emph{X}: the estimates, and the
#'   distribution of the statistic when there is no effect, are the same
#'   whether the predictor is fixed or sampled, but the distribution when there
#'   is an effect, and so the power, is not (Sampson, 1974; Gatsonis &
#'   Sampson, 1989; Kelley & Maxwell, 2003). A correlation study samples
#'   \emph{both} variables. Conditional on the sampled \emph{X}, the statistic
#'   is noncentral \emph{F}, but its noncentrality parameter is itself random:
#'   it is proportional to the sampled spread of \emph{X}, which follows a
#'   chi square distribution on \emph{N} - 1 degrees of freedom. The marginal
#'   distribution is therefore a chi square mixture of noncentral \emph{F}
#'   distributions, which is more dispersed than any single noncentral
#'   \emph{F}, and treating the statistic as noncentral \emph{F} overstates the
#'   power of the planned study in every case, by roughly 1/\emph{n} for a plan
#'   of \emph{n} participants: about .01 at \emph{n} = 100, .02 at 50, .05 at
#'   20, .10 at 10, and up to about .2 for the smallest plans.
#'
#'   \code{ss_buc_correlation} uses the mixture, which has an exact closed form
#'   (the Poisson weights of the noncentral \emph{F} become negative binomial
#'   weights), for both the correction and the planned study's power. The
#'   returned sample size is therefore exact for bivariate normal data rather
#'   than a few participants light, and it is typically 2 to 5 participants
#'   larger than the fixed-predictor calculation (1 to 3 for the smallest plans,
#'   up to about 8 for plans in the hundreds). Near the assurance ceiling the
#'   difference reverses and the mixture calls for a sample a few percent
#'   smaller, which can be hundreds of participants because plans there are
#'   very large (7,346 rather than the 7,722 of treating the predictor as
#'   fixed, for \emph{r} = .50 from a prior study of 20 at an
#'   \code{assurance} of .5). If your predictor really is fixed by design, the
#'   analysis is a regression slope rather than a correlation, and
#'   \code{\link{ss_buc_reg_coef}} with \code{p = 1} and
#'   \code{random_predictors = FALSE} is the planner for it. With one sampled
#'   predictor the regression planners reduce to this one:
#'   \code{ss_buc_reg_coef(p = 1, random_predictors = TRUE)} gives this
#'   planner's own result for the same prior \emph{t} and \emph{N}, exactly,
#'   as do \code{\link{ss_buc_R2}} with \code{p = 1} and
#'   \code{\link{ss_buc_reg_joint}} with \code{p = 1} and
#'   \code{p_joint = 1} given \code{F_observed} = \eqn{t^2}{t^2}; each
#'   states the planning effect size in its own metric.
#'
#'   The largest assurance the prior result can support is the same in either
#'   frame, since at a zero effect the mixture is the central \emph{F}
#'   distribution.
#'
#'   Assurance is the proportion of times that power will be at or above the
#'   desired level, if the experiment were to be reproduced many times. For
#'   example, \code{assurance = .5} means that power will be above the desired
#'   level half of the time, but below the desired level the other half of the
#'   time. Selecting \code{assurance = .5} (selecting the noncentrality
#'   parameter at the 50th percentile of the likelihood distribution) results in
#'   a median unbiased estimate of the population noncentrality parameter and
#'   does not correct for uncertainty. In order to correct for uncertainty,
#'   \code{assurance > .5} can be selected, which corresponds to selecting the
#'   noncentrality parameter associated with the (1 - assurance) quantile of the
#'   likelihood distribution.
#'
#'   \code{alpha_prior} is the significance level a prior result had to reach
#'   to be published: the correction assumes that only results with \emph{p}
#'   values below \code{alpha_prior} reach the literature, and so reach you.
#'   Different \eqn{\alpha} levels represent differing amounts of publication
#'   bias. For example, setting \code{alpha_prior = .20} assumes that studies
#'   with \emph{p} values less than .20 are published, whereas those with
#'   larger \emph{p} values are not, which reflects less severe publication
#'   bias than the default of .05. If the previous study of interest has not
#'   been subjected to publication bias (e.g., a pilot study), set
#'   \code{alpha_prior = 1}. Because no \emph{p} value exceeds 1, this means
#'   no publication filter at all: any result would have been published
#'   whatever its \emph{p} value, so there is no publication bias to correct,
#'   and only the correction for uncertainty, set by \code{assurance},
#'   remains.
#'
#'   In some cases, the corrected noncentrality parameter for a given level of
#'   assurance will be estimated to be zero. This is an indication that, at the
#'   desired level of assurance, the previous study's effect cannot be
#'   accurately estimated due to high levels of uncertainty and bias. When this
#'   happens, subsequent sample size planning is not possible with the chosen
#'   specifications, and the function stops with an informative error rather
#'   than returning a result. The error reports the largest assurance the prior
#'   result can support, which is \eqn{1 - p/\alpha}{1 - p/alpha}, where
#'   \emph{p} is the prior study's \emph{p} value and \eqn{\alpha} is
#'   \code{alpha_prior}; at any higher assurance the corrected noncentrality
#'   parameter is driven to zero. Two remedies are available. First, users can
#'   select a lower value of assurance (e.g., .8 instead of .95), at or below
#'   the reported ceiling. Second, users can reduce the influence of publication
#'   bias by setting \code{alpha_prior} at a value greater than .05, which raises
#'   the ceiling. Note that \code{alpha_prior} is not a property of the prior
#'   study, which is fixed, but the analyst's assumption about the publication
#'   threshold its literature faced. When the ceiling is at or below the
#'   recommended floor of .5, lowering assurance cannot help and
#'   \code{alpha_prior} is the only remaining lever. It is possible to correct
#'   for uncertainty only by setting \code{alpha_prior = 1} and choosing the
#'   desired level of assurance. We encourage users to make the adjustments as
#'   minimal as possible.
#'
#'   The returned \code{actual_power} is the statistical power of the plan
#'   this function reports: the power of a study of the returned sample size
#'   when the effect is the returned bias and uncertainty adjusted
#'   noncentrality parameter. Because a sample size must be a whole number,
#'   \code{actual_power} is never below the \code{desired_power} that was
#'   requested and is usually a little above it.
#'
#'   The observed correlation may be entered with either sign, since the
#'   publication rule the correction assumes is two-sided and only the
#'   magnitude enters the computation.
#'
#'   The returned \code{ncp_adjusted} is the noncentrality parameter the
#'   planning effect size implies at the prior study's sample size,
#'   \eqn{N\rho^2/(1 - \rho^2)}{N rho^2/(1 - rho^2)}, on the scale of the
#'   \emph{F} statistic, as in the regression planners (the \emph{t}-based
#'   planners report theirs on the scale of \emph{t} instead). The planning
#'   effect size itself, the correlation
#'   \code{sqrt(ncp_adjusted / (ncp_adjusted + N))}, is returned as the
#'   \code{rho_planning} row.
#'
#' @param r_observed Observed Pearson correlation from a previous study used to
#'   plan sample size for a planned study. Either sign is accepted; only the
#'   magnitude enters the correction. Supply either \code{r_observed} or
#'   \code{t_observed}, not both.
#' @param N Total sample size of the previous study.
#' @param t_observed Observed \emph{t} statistic for the correlation from the
#'   previous study, an alternative to \code{r_observed}.
#' @template planning-params
#'
#' @templateVar size_phrase total sample size
#' @templateVar ncp_scale on the scale of the \emph{F} statistic, \eqn{N\rho^2/(1 - \rho^2)}{N rho^2/(1 - rho^2)} at the prior study's sample size
#' @templateVar effect_phrase \code{rho_planning}, the population correlation (a magnitude), \eqn{\sqrt{\lambda/(\lambda + N)}}{sqrt(lambda/(lambda + N))}, where \eqn{\lambda}{lambda} is \code{ncp_adjusted} and \emph{N} the prior sample size; the power statement that follows is for sampled variables, the exact distribution this planner uses
#' @template return
#'
#' @export
#'
#' @examples
#' result <- ss_buc_correlation(r_observed = .35, N = 100, alpha_prior = .05,
#'   alpha_planned = .05, assurance = .80, desired_power = .80)
#' result
#'
#' # The equivalent call from the observed t statistic:
#' ss_buc_correlation(t_observed = .35 * sqrt(98 / (1 - .35^2)), N = 100)
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance. The correlation
#' # above supports .99, so a weaker prior is needed to reach the ceiling
#' # (here near .86):
#' try(ss_buc_correlation(r_observed = .27, N = 100, assurance = .95))
#'
#' @seealso \code{\link{ss_buc_reg_coef}} for a coefficient in a model with
#'   more than one predictor (with \code{random_predictors = TRUE} for
#'   sampled predictors).
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
#' @references Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation:
#'   Exact power and sample size calculations. \emph{Psychological Bulletin,
#'   106,} 516--524. \doi{10.1037/0033-2909.106.3.516}
#'
#'   Sampson, A. R. (1974). A tale of two regressions. \emph{Journal of the
#'   American Statistical Association, 69,} 682--689.
#'   \doi{10.1080/01621459.1974.10480189}
ss_buc_correlation <- function(r_observed, N, t_observed, alpha_prior = .05,
                               alpha_planned = .05, assurance = .80,
                               desired_power = .80) {
  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, desired_power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  desired_power <- v$desired_power

  if (missing(N)) stop("You need to specify 'N', which is the total sample size of the previous study.", call. = FALSE)
  .check_count(N, "N", min = 4)
  if (missing(r_observed) && missing(t_observed)) {
    stop("You must specify the prior study's result: either 'r_observed' (the observed correlation) or 't_observed' (its t statistic).", call. = FALSE)
  }
  if (!missing(r_observed) && !missing(t_observed)) {
    stop("Specify either 'r_observed' or 't_observed', not both.", call. = FALSE)
  }
  if (!missing(r_observed)) {
    .check_scalar_finite(r_observed, "r_observed")
    if (abs(r_observed) >= 1) stop("'r_observed' must be a correlation strictly between -1 and 1.", call. = FALSE)
    r_input <- r_observed
    # The correlation's t is the simple-regression slope t on N - 2 df.
    t_stat <- abs(r_observed) * sqrt((N - 2) / (1 - r_observed^2))
  } else {
    .check_scalar_finite(t_observed, "t_observed")
    r_input <- NULL
    t_stat <- abs(t_observed)
  }

  df_numerator <- 1
  df_denominator <- N - 2

  value_critical <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator)

  if (t_stat^2 <= value_critical) .stop_nonsignificant("correlation", "the prior result", alpha_prior = alpha_prior)

  # Both the correction and the planned-study power use the exact
  # random-predictor distribution (.p_correlation_f), because a correlation
  # study samples both variables. The adjusted parameter is reported on the
  # usual scale, the noncentrality parameter the adjusted correlation implies
  # at the prior study's N, so it stays comparable with the other planners;
  # the correlation itself, sqrt(ncp / (ncp + N)), is the planning effect size
  # row.
  ncp_solution <- .solve_ncp_assurance(
    .tm_correlation(t_stat^2, value_critical, N), assurance)
  ncp <- ncp_solution$ncp

  if (ncp == 0) .stop_zero_ncp(ncp_solution$ceiling)

  rho2_adjusted <- ncp / (ncp + N)
  power_at <- function(n_rep) {
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = n_rep - 2)
    1 - .p_correlation_f(critical_F, n_rep, rho2_adjusted)
  }
  output_n <- .smallest_n_for_power(power_at, desired_power, start = 4)

  actual_power <- power_at(output_n)
  .bucss_plan_result(
    sample_size = output_n,
    size_term = "necessary_N",
    effect_size_term = "rho_planning",
    effect_size = .es_rho(ncp, N),
    actual_power = actual_power,
    df_effect = df_numerator,
    df_error = output_n - 2,
    ncp = ncp,
    design = "Pearson correlation",
    sample_size_unit = "total",
    assurance_ceiling = ncp_solution$ceiling,
    inputs = list(r_observed = r_input,
                  t_observed = if (is.null(r_input)) t_observed else NULL,
                  N = N, alpha_prior = alpha_prior_input,
                  alpha_planned = alpha_planned, assurance = assurance,
                  desired_power = desired_power)
  )
}
