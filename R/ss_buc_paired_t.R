#' Necessary sample size to reach desired power for a dependent (paired) t test
#' (or standardized mean difference) using a publication bias and uncertainty
#' correction procedure
#'
#' @description \code{ss_buc_paired_t} returns the necessary sample size (the
#'   number of pairs) to achieve a desired level of statistical power for a
#'   planned study using a dependent \emph{t} test, based on information obtained
#'   from a previous study. The effect from the previous study can be corrected
#'   for publication bias and/or uncertainty to provide a sample size that will
#'   achieve more accurate statistical power for a planned study, when compared
#'   to approaches that use a sample effect size at face value or rely on sample
#'   size only. The bias and uncertainty adjusted previous study noncentrality
#'   parameter is also returned, which can be transformed to various effect size
#'   metrics.
#'
#'   \code{ss_buc_smd_paired} is an alias for \code{ss_buc_paired_t}: the two are
#'   the same function. The paired \emph{t} test and the standardized mean
#'   difference of the difference scores describe the same comparison, so the
#'   same planner serves both framings. Prefer the test-specific name
#'   (\code{ss_buc_paired_t}) when the prior \emph{t} was computed from raw
#'   (unstandardized) data, since the planner works directly from the observed
#'   \emph{t}.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss_buc_paired_t} uses the observed
#'   \emph{t} value and sample size from a previous study to correct the
#'   noncentrality parameter associated with the effect of interest for
#'   publication bias and/or uncertainty. This new estimated noncentrality
#'   parameter is then used to calculate the necessary sample size to achieve
#'   the desired level of power in the planned study.
#'
#'   The approach uses a likelihood function of a truncated noncentral
#'   \emph{F} distribution, where the truncation occurs due to small effect
#'   sizes being unobserved due to publication bias. The numerator of the
#'   likelihood function is the density of a noncentral \emph{F} distribution.
#'   The denominator is the power of the test, which serves to truncate the
#'   distribution. In the two-group case, this formula reduces to the density of
#'   a truncated noncentral \emph{t} distribution. (See Taylor & Muller, 1996,
#'   Equation 2.1, and Anderson & Maxwell, 2017, for more details.)
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
#'   The observed \emph{t} may be entered with either sign. The publication
#'   rule the correction assumes is two-sided and the truncated likelihood
#'   is symmetric in the sign of \emph{t}, so only the magnitude enters the
#'   computation: \code{t_observed = -3} plans exactly the sample size that
#'   \code{t_observed = 3} plans. The sign of a paired \emph{t} records the
#'   direction the within-pair difference was taken, an arbitrary coding
#'   choice rather than evidence. When planning a replication, confirm that
#'   the magnitude belongs to an effect in the direction the planned study
#'   is designed to detect (Anderson & Kelley, 2024): the correction
#'   adjusts the size of the prior effect, not its direction.
#'
#'   If you are working from a standardized mean difference (Cohen's \eqn{d_z})
#'   of the difference scores rather than a \emph{t} statistic, convert it before
#'   calling: for a paired design \eqn{t = d_z\sqrt{N}}, where \emph{N} is the
#'   number of pairs.
#'
#' @param t_observed Observed \emph{t} value from a previous study used to plan
#'   sample size for a planned study. Either sign is accepted: the
#'   publication rule is two-sided, so only the magnitude enters the
#'   correction (see Details).
#' @param N Total sample size (the number of pairs) of the previous study.
#' @template planning-params
#'
#' @templateVar size_phrase number of pairs
#' @template return
#'
#' @export
#'
#' @examples
#' result <- ss_buc_paired_t(t_observed = 3, N = 40, alpha_prior = .05,
#'   alpha_planned = .05, assurance = .80, desired_power = .80)
#' result
#'
#' # ss_buc_smd_paired is the same function under an effect size name
#' ss_buc_smd_paired(t_observed = 3, N = 40)
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .90):
#' try(ss_buc_paired_t(t_observed = 3, N = 40, assurance = .95))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
ss_buc_paired_t <- function(t_observed, N, alpha_prior = .05, alpha_planned = .05,
                            assurance = .80, desired_power = .80) {
  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, desired_power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  desired_power <- v$desired_power

  if (missing(N)) stop("You need to specify a sample size (i.e., the number of pairs) used in the previous study.", call. = FALSE)
  .check_scalar_finite(t_observed, "t_observed")

  # The publication rule is two-sided and TM is symmetric in the sign of t,
  # so only the magnitude enters the correction; the sign records the
  # direction the paired difference was taken, a coding choice, not
  # evidence. The user's signed value is echoed unchanged in the stored
  # inputs.
  t_stat <- abs(t_observed)
  .check_count(N, "N", min = 2)

  DF <- N - 1

  value_critical <- qt(1 - alpha_prior / 2, df = DF)
  if (t_stat <= value_critical) stop("Your observed t statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 't_observed' exceeds the critical value.", call. = FALSE)

  ncp_solution <- .solve_ncp_assurance(
    .tm_t(t_stat, value_critical, DF), assurance)
  ncp <- ncp_solution$ncp

  if (ncp == 0) .stop_zero_ncp(ncp_solution$ceiling)

  power_at <- function(n_rep) {
    denom_df <- n_rep - 1
    critical_t <- qt(1 - alpha_planned / 2, df = denom_df)
    scaled <- sqrt(n_rep / N) * ncp
    (1 - pt(critical_t, df = denom_df, ncp = scaled)) +
      pt(-1 * critical_t, df = denom_df, ncp = scaled)
  }
  output_n <- .smallest_n_for_power(power_at, desired_power)

  df_error <- output_n - 1
  actual_power <- power_at(output_n)
  .bucss_power_result(
    sample_size = output_n,
    size_term = "necessary_N",
    actual_power = actual_power,
    df_effect = 1,
    df_error = df_error,
    ncp = ncp,
    design = "Dependent (paired) t test",
    sample_size_unit = "number of pairs",
    assurance_ceiling = ncp_solution$ceiling,
    inputs = list(t_observed = t_observed, N = N, alpha_prior = alpha_prior_input,
                  alpha_planned = alpha_planned, assurance = assurance,
                  desired_power = desired_power)
  )
}

#' @rdname ss_buc_paired_t
#' @export
ss_buc_smd_paired <- ss_buc_paired_t
