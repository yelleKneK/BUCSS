#' Necessary sample size to reach desired power for a Pearson correlation using
#' a publication bias and uncertainty correction procedure
#'
#' @description \code{ss_buc_correlation} returns the necessary total sample
#'   size to achieve a desired level of statistical power for a test of a
#'   Pearson correlation in a planned study, based on information obtained from
#'   a previous study. The effect from the previous study can be corrected for
#'   publication bias and/or uncertainty to provide a sample size that will
#'   achieve more accurate statistical power for a planned study, when compared
#'   to approaches that use a sample effect size at face value or rely on sample
#'   size only. The bias and uncertainty adjusted previous study noncentrality
#'   parameter is also returned, which can be transformed to various effect size
#'   metrics.
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
#'   The test of a Pearson correlation is the test of the slope in a simple
#'   regression: the observed correlation \emph{r} with sample size \emph{N}
#'   gives \eqn{t = r\sqrt{(N - 2)/(1 - r^2)}} on \emph{N} - 2 degrees of
#'   freedom, and the planner works from that \emph{t}. Supplying
#'   \code{t_observed} directly instead of \code{r_observed} is equivalent.
#'
#'   The approach uses a likelihood function of a truncated noncentral
#'   \emph{F} distribution, where the truncation occurs due to small effect
#'   sizes being unobserved due to publication bias. The numerator of the
#'   likelihood function is the density of a noncentral \emph{F} distribution.
#'   The denominator is the power of the test, which serves to truncate the
#'   distribution. In the single predictor case, this formula reduces to the
#'   density of a truncated noncentral \emph{t} distribution. (See Taylor &
#'   Muller, 1996, Equation 2.1, and Anderson & Maxwell, 2017, for more
#'   details.)
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
#'   If the previous study of interest has not been subjected to publication
#'   bias (e.g., a pilot study), \code{alpha_prior} can be set to 1 to indicate
#'   no publication bias. Alternative \eqn{\alpha} levels can also be
#'   accommodated to represent differing amounts of publication bias. For
#'   example, setting \code{alpha_prior = .20} would reflect less severe
#'   publication bias than the default of .05. In essence, setting
#'   \code{alpha_prior} at .20 assumes that studies with \emph{p} values less
#'   than .20 are published, whereas those with larger \emph{p} values are not.
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
#'   \strong{One approximation to be aware of.} Like the other planners here,
#'   \code{ss_buc_correlation} plans in the fixed-predictor (regression)
#'   frame, in which the predictor values are treated as fixed by design. A
#'   correlation study usually samples both variables, the random-predictor
#'   frame, in which the exact power is slightly lower. The bias and
#'   uncertainty correction itself is almost unaffected by the distinction
#'   (the corrected correlation moves by less than .005 over a wide range),
#'   but the planned sample size runs a few participants light: across the
#'   range checked, the exact random-predictor plan needed 2 to 5 more
#'   participants. Treat the returned sample size as a close lower bound and
#'   consider adding a small margin.
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
#'   more than one predictor.
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
ss_buc_correlation <- function(r_observed, N, t_observed, alpha_prior = .05,
                               alpha_planned = .05, assurance = .80,
                               desired_power = .80) {
  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, desired_power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  desired_power <- v$desired_power

  if (missing(N)) stop("You need to specify 'N', which is the total sample size of the original study.", call. = FALSE)
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

  if (t_stat^2 <= value_critical) stop("Your observed correlation is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so the prior result exceeds the critical value.", call. = FALSE)

  ncp_solution <- .solve_ncp_assurance(
    .tm_f(t_stat^2, value_critical, df_numerator, df_denominator), assurance)
  ncp <- ncp_solution$ncp

  if (ncp == 0) .stop_zero_ncp(ncp_solution$ceiling)

  power_at <- function(n_rep) {
    denom_df <- n_rep - 2
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = denom_df)
    1 - pf(critical_F, df1 = df_numerator, df2 = denom_df, ncp = (n_rep / N) * ncp)
  }
  output_n <- .smallest_n_for_power(power_at, desired_power, start = 4)

  actual_power <- power_at(output_n)
  .bucss_power_result(
    sample_size = output_n,
    size_term = "necessary_N",
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
