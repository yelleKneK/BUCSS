#' Necessary sample size to reach desired power for a two-group multivariate
#' comparison using a publication bias and uncertainty correction procedure
#'
#' @description \code{ss_buc_manova} returns the necessary per-group sample size
#'   to achieve a desired level of statistical power for a planned study
#'   comparing two groups on several outcome variables at once (Hotelling's
#'   \eqn{T^2}, equivalently a one-way multivariate analysis of variance with
#'   two groups), based on information obtained from a previous study. The
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
#'   \strong{Scope.} The correction this package implements works on a single
#'   test statistic with an exact noncentral \emph{F} distribution. A
#'   multivariate hypothesis has one only when its rank is one, which is the
#'   case for a comparison of two groups: Hotelling's \eqn{T^2} on \emph{p}
#'   outcome variables with total sample size \emph{N} maps exactly to
#'   \eqn{F = T^2 (N - p - 1) / (p (N - 2))} on \emph{p} and \emph{N} - \emph{p}
#'   - 1 degrees of freedom. Supply either the \emph{F} or the \eqn{T^2}. For a
#'   hypothesis of rank greater than one (three or more groups, say) the
#'   multivariate noncentrality is a set of eigenvalues rather than a single
#'   number, the usual test statistics have only approximate \emph{F} forms, and
#'   this correction does not apply; \code{ss_buc_manova} will not accept such a
#'   design.
#'
#'   The approach uses a likelihood function of a truncated noncentral
#'   \emph{F} distribution, where the truncation occurs due to small effect
#'   sizes being unobserved due to publication bias. The numerator of the
#'   likelihood function is the density of a noncentral \emph{F} distribution.
#'   The denominator is the power of the test, which serves to truncate the
#'   distribution. Thus, the ratio of the numerator and the denominator is a
#'   truncated noncentral \emph{F} distribution. (See Taylor & Muller, 1996,
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
#'   \code{ss_buc_manova} assumes that the planned study will have equal
#'   \emph{n} per group. If the prior \emph{N} is odd, the per-group size is
#'   rounded both up and down, the planned sample size is computed under each,
#'   and the larger of the two is returned, with the smaller of the two adjusted
#'   noncentrality parameters, to be conservative.
#'
#' @param F_observed Observed \emph{F} value for the multivariate comparison
#'   from a previous study. Supply either \code{F_observed} or
#'   \code{T2_observed}, not both.
#' @param N Total sample size of the previous study (both groups combined).
#' @param p_variables Number of outcome variables compared.
#' @param T2_observed Observed Hotelling's \eqn{T^2} from the previous study,
#'   an alternative to \code{F_observed}.
#' @template planning-params
#'
#' @templateVar size_phrase per-group sample size
#' @template return
#'
#' @export
#'
#' @examples
#' result <- ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3,
#'   alpha_prior = .05, alpha_planned = .05, assurance = .80,
#'   desired_power = .80)
#' result
#'
#' # The equivalent call from Hotelling's T squared:
#' ss_buc_manova(T2_observed = 4.5 * 3 * 78 / 76, N = 80, p_variables = 3)
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .87):
#' try(ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3,
#'   assurance = .95))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
ss_buc_manova <- function(F_observed, N, p_variables, T2_observed,
                          alpha_prior = .05, alpha_planned = .05,
                          assurance = .80, desired_power = .80) {
  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, desired_power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  desired_power <- v$desired_power

  if (missing(N)) stop("You must specify 'N', which is the total sample size of the previous study.", call. = FALSE)
  if (missing(p_variables)) stop("You must specify 'p_variables', the number of outcome variables compared.", call. = FALSE)
  .check_count(p_variables, "p_variables", min = 1)
  .check_count(N, "N", min = 4)
  if (missing(F_observed) && missing(T2_observed)) {
    stop("You must specify the prior study's result: either 'F_observed' or 'T2_observed'.", call. = FALSE)
  }
  if (!missing(F_observed) && !missing(T2_observed)) {
    stop("Specify either 'F_observed' or 'T2_observed', not both.", call. = FALSE)
  }
  if (N - p_variables - 1 < 1) {
    stop("The combination of 'N' and 'p_variables' leaves no error degrees of freedom. Reduce 'p_variables' or increase 'N'.", call. = FALSE)
  }
  if (!missing(T2_observed)) {
    .check_scalar_finite(T2_observed, "T2_observed")
    # Hotelling's T squared maps exactly to F for a two-group comparison.
    F_stat <- T2_observed * (N - p_variables - 1) / (p_variables * (N - 2))
    T2_input <- T2_observed
    F_input <- NULL
  } else {
    .check_scalar_finite(F_observed, "F_observed")
    F_stat <- F_observed
    T2_input <- NULL
    F_input <- F_observed
  }

  df_numerator <- p_variables

  # Two-sided rounding of the prior per-group size, as in the other
  # between-subjects designs.
  n_ru <- ceiling(N / 2)
  n_rd <- floor(N / 2)

  branch <- function(n_prior) {
    N_branch <- 2 * n_prior
    df_denominator <- N_branch - p_variables - 1
    if (df_denominator < 1) {
      stop("The combination of 'N' and 'p_variables' leaves no error degrees of freedom. Reduce 'p_variables' or increase 'N'.", call. = FALSE)
    }
    crit <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator)
    if (F_stat <= crit) stop("Your observed multivariate test is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so the prior result exceeds the critical value.", call. = FALSE)
    solution <- .solve_ncp_assurance(
      .tm_f(F_stat, crit, df_numerator, df_denominator), assurance)
    if (solution$ncp == 0) .stop_zero_ncp(solution$ceiling)
    solution
  }
  solution_ru <- branch(n_ru)
  solution_rd <- branch(n_rd)

  # The multivariate noncentrality is proportional to the per-group sample
  # size, as in the univariate two-group case.
  power_at <- function(n_rep, n_prior, ncp_b) {
    denom_df <- (2 * n_rep) - p_variables - 1
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = denom_df)
    1 - pf(critical_F, df1 = df_numerator, df2 = denom_df,
           ncp = (n_rep / n_prior) * ncp_b)
  }
  n_start <- max(2, ceiling((p_variables + 2) / 2))
  repn_ru <- .smallest_n_for_power(function(k) power_at(k, n_ru, solution_ru$ncp),
                                   desired_power, start = n_start)
  repn_rd <- .smallest_n_for_power(function(k) power_at(k, n_rd, solution_rd$ncp),
                                   desired_power, start = n_start)

  output_n <- max(repn_ru, repn_rd)
  actual_power <- if (solution_rd$ncp <= solution_ru$ncp)
    power_at(output_n, n_rd, solution_rd$ncp) else
      power_at(output_n, n_ru, solution_ru$ncp)

  .bucss_power_result(
    sample_size = output_n,
    size_term = "necessary_n_per_group",
    actual_power = actual_power,
    df_effect = df_numerator,
    df_error = (2 * output_n) - p_variables - 1,
    ncp = min(solution_ru$ncp, solution_rd$ncp),
    design = "Two-group multivariate comparison (Hotelling's T squared)",
    sample_size_unit = "per group",
    assurance_ceiling = min(solution_ru$ceiling, solution_rd$ceiling),
    total_n = 2 * output_n,
    inputs = list(F_observed = F_input, T2_observed = T2_input, N = N,
                  p_variables = p_variables, alpha_prior = alpha_prior_input,
                  alpha_planned = alpha_planned, assurance = assurance,
                  desired_power = desired_power)
  )
}
