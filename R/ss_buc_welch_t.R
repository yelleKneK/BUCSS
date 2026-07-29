#' Necessary sample size to reach desired power for a Welch (unequal variance) t
#' test using a publication bias and uncertainty correction procedure
#'
#' @description \code{ss_buc_welch_t} returns the necessary per-group sample
#'   size to achieve a desired level of statistical power for a planned study
#'   using a Welch \emph{t} test, the two-group comparison that does not assume
#'   equal variances, based on information obtained from a previous study. The
#'   effect from the previous study can be corrected for publication bias and/or
#'   uncertainty to provide a sample size that will achieve more accurate
#'   statistical power for a planned study, when compared to approaches that use
#'   a sample effect size at face value or rely on sample size only. The bias
#'   and uncertainty adjusted previous study noncentrality parameter is also
#'   returned, which can be transformed to various effect size metrics.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   \strong{This planner is an approximation, and needs one assumption the
#'   others do not.} A Welch \emph{t} statistic is not exactly distributed as a
#'   noncentral \emph{t}: its denominator degrees of freedom are estimated from
#'   the data by the Welch-Satterthwaite formula, and they depend on the ratio
#'   of the two population standard deviations. \code{ss_buc_welch_t} treats the
#'   statistic as noncentral \emph{t} on the Welch-Satterthwaite degrees of
#'   freedom implied by the prior study's group sizes and the standard deviation
#'   ratio you supply, which is the usual approximation in power analysis for
#'   this test. In simulation that treatment is accurate: across ratios from 1
#'   to 3 the Type I error of the real Welch test stays between .047 and .054
#'   against a nominal .05, and the planned study's power is within about .006
#'   of what this function predicts, with the agreement no worse at the
#'   extreme ratios.
#'
#'   \strong{The assumption that matters is the ratio itself.} Planning assumes
#'   the ratio you supply is the ratio the planned study will have, and getting
#'   that wrong costs far more than the distributional approximation does. For
#'   the example below, a plan built assuming equal variances delivers about
#'   .43 power rather than .80 if the true ratio turns out to be 2. Re-run
#'   across a range of plausible ratios and plan for the least favorable one
#'   you find credible. The required sample size does not always rise with the
#'   ratio: it falls when the group assumed to be more variable is the smaller
#'   one.
#'
#'   When the ratio is 1 and the prior groups were the same size, this planner
#'   returns exactly what \code{\link{ss_buc_independent_t}} returns. With
#'   unequal prior group sizes the two differ, and should: the equal-variance
#'   planner reduces the prior study to an equivalent equal-\emph{n} design at
#'   the harmonic mean, while this one uses the actual Welch-Satterthwaite
#'   degrees of freedom, which are smaller when the groups are unbalanced. The
#'   same prior \emph{t} is then less compelling and the plan is larger.
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
#'   The observed \emph{t} may be entered with either sign; the publication
#'   rule the correction assumes is two-sided, so only the magnitude enters
#'   the computation. The planned study is assumed to have equal group sizes,
#'   which is the allocation that maximizes power when the variances are equal
#'   and is a reasonable default otherwise. It is not the optimal allocation
#'   when they are not: sampling in proportion to the standard deviations would
#'   need about 11 percent fewer participants in total at a ratio of 2, and
#'   about 25 percent fewer at a ratio of 3.
#'
#' @param t_observed Observed Welch \emph{t} value from a previous study used to
#'   plan sample size for a planned study. Either sign is accepted.
#' @param n_1,n_2 Group sample sizes of the previous study.
#' @param sd_ratio Ratio of the second group's standard deviation to the first
#'   group's, in the previous study, assumed to hold in the planned study as
#'   well. A value of 1 means equal variances. Give this, or the two group
#'   standard deviations, or the two group variances, but not more than one of
#'   the three.
#' @param sd_1,sd_2 Group standard deviations of the previous study, an
#'   alternative to \code{sd_ratio}. Only their ratio enters the computation,
#'   so the unit does not matter.
#' @param var_1,var_2 Group variances of the previous study, a second
#'   alternative to \code{sd_ratio}.
#' @template planning-params
#'
#' @templateVar size_phrase per-group sample size
#' @template return
#'
#' @export
#'
#' @examples
#' # A prior study with 40 and 55 participants whose second group was half
#' # again as variable as the first.
#' result <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55,
#'   sd_ratio = 1.5, alpha_prior = .05, alpha_planned = .05,
#'   assurance = .80, desired_power = .80)
#' result
#'
#' # The same study described by the group standard deviations, or by the
#' # group variances, rather than by their ratio.
#' ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_1 = 8.4, sd_2 = 12.6)
#' ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, var_1 = 70.56,
#'   var_2 = 158.76)
#'
#' # How sensitive is the plan to the assumed spread?
#' ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1)
#' ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 2)
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .93):
#' try(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5,
#'   assurance = .95))
#'
#' @seealso \code{\link{ss_buc_independent_t}} for the equal-variance test.
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu})
#'
#' @template references
ss_buc_welch_t <- function(t_observed, n_1, n_2, sd_ratio = 1,
                           sd_1, sd_2, var_1, var_2,
                           alpha_prior = .05, alpha_planned = .05,
                           assurance = .80, desired_power = .80) {
  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, desired_power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  desired_power <- v$desired_power

  if (missing(n_1) || missing(n_2)) {
    stop("You must specify both group sample sizes of the previous study, 'n_1' and 'n_2'.", call. = FALSE)
  }
  .check_scalar_finite(t_observed, "t_observed")
  .check_count(n_1, "n_1", min = 2)
  .check_count(n_2, "n_2", min = 2)

  # The spread may be described three ways. Only the ratio enters the
  # computation, so the two pairs are folded into 'sd_ratio' here and the form
  # the user chose is echoed back unchanged in the returned result.
  gave_sd <- !missing(sd_1) || !missing(sd_2)
  gave_var <- !missing(var_1) || !missing(var_2)
  if (sum(c(!missing(sd_ratio), gave_sd, gave_var)) > 1L) {
    stop("Describe the groups' spread one way only: 'sd_ratio', or 'sd_1' and 'sd_2', or 'var_1' and 'var_2'.", call. = FALSE)
  }
  if (gave_sd) {
    if (missing(sd_1) || missing(sd_2)) stop("You must specify both 'sd_1' and 'sd_2', the group standard deviations of the previous study.", call. = FALSE)
    .check_scalar_finite(sd_1, "sd_1")
    .check_scalar_finite(sd_2, "sd_2")
    if (sd_1 <= 0 || sd_2 <= 0) stop("'sd_1' and 'sd_2' must be positive numbers.", call. = FALSE)
    sd_ratio <- sd_2 / sd_1
  } else if (gave_var) {
    if (missing(var_1) || missing(var_2)) stop("You must specify both 'var_1' and 'var_2', the group variances of the previous study.", call. = FALSE)
    .check_scalar_finite(var_1, "var_1")
    .check_scalar_finite(var_2, "var_2")
    if (var_1 <= 0 || var_2 <= 0) stop("'var_1' and 'var_2' must be positive numbers.", call. = FALSE)
    sd_ratio <- sqrt(var_2 / var_1)
  }
  .check_scalar_finite(sd_ratio, "sd_ratio")
  if (sd_ratio <= 0) stop("'sd_ratio' must be a positive number.", call. = FALSE)

  t_stat <- abs(t_observed)

  # Welch-Satterthwaite degrees of freedom, with the first group's standard
  # deviation set to 1 without loss of generality.
  welch_df <- function(a, b, ratio) {
    v1 <- 1 / a
    v2 <- ratio^2 / b
    (v1 + v2)^2 / (v1^2 / (a - 1) + v2^2 / (b - 1))
  }
  df_prior <- welch_df(n_1, n_2, sd_ratio)

  value_critical <- qt(1 - alpha_prior / 2, df = df_prior)

  if (t_stat <= value_critical) stop("Your observed t statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 't_observed' exceeds the critical value.", call. = FALSE)

  ncp_solution <- .solve_ncp_assurance(
    .tm_t(t_stat, value_critical, df_prior), assurance)
  ncp <- ncp_solution$ncp

  if (ncp == 0) .stop_zero_ncp(ncp_solution$ceiling)

  # The noncentrality parameter is the mean difference divided by the standard
  # error of that difference, so it rescales by the ratio of the prior standard
  # error to the planned one. With equal planned group sizes the planned
  # variance term is (1 + sd_ratio^2) / n_rep.
  se2_prior <- 1 / n_1 + sd_ratio^2 / n_2
  power_at <- function(n_rep) {
    denom_df <- welch_df(n_rep, n_rep, sd_ratio)
    critical_t <- qt(1 - alpha_planned / 2, df = denom_df)
    se2_planned <- (1 + sd_ratio^2) / n_rep
    scaled <- sqrt(se2_prior / se2_planned) * ncp
    (1 - pt(critical_t, df = denom_df, ncp = scaled)) +
      pt(-1 * critical_t, df = denom_df, ncp = scaled)
  }
  output_n <- .smallest_n_for_power(power_at, desired_power, start = 2)

  actual_power <- power_at(output_n)
  .bucss_power_result(
    sample_size = output_n,
    size_term = "necessary_n_per_group",
    actual_power = actual_power,
    df_effect = 1,
    df_error = welch_df(output_n, output_n, sd_ratio),
    ncp = ncp,
    design = "Welch (unequal variance) t test",
    sample_size_unit = "per group",
    assurance_ceiling = ncp_solution$ceiling,
    total_n = 2 * output_n,
    inputs = c(list(t_observed = t_observed, n_1 = n_1, n_2 = n_2),
               # Echo the spread in the form it was given, and only that form,
               # so the echoed rows remain a valid description of the call.
               if (gave_sd) list(sd_1 = sd_1, sd_2 = sd_2)
               else if (gave_var) list(var_1 = var_1, var_2 = var_2)
               else list(sd_ratio = sd_ratio),
               list(alpha_prior = alpha_prior_input,
                    alpha_planned = alpha_planned, assurance = assurance,
                    desired_power = desired_power))
  )
}
