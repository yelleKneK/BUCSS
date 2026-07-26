#' Necessary sample size to reach desired power for a between-subjects ANOVA
#' with any number of factors using a publication bias and uncertainty
#' correction procedure
#'
#' @description \code{ss_buc_factorial_anova_general} returns the necessary
#'   per-cell sample size to achieve a desired level of statistical power for a
#'   planned study testing any type of effect (omnibus or contrast) using a
#'   fully between-subjects ANOVA with any number of factors, based on
#'   information obtained from a previous study. The effect from the previous
#'   study can be corrected for publication bias and/or uncertainty to provide a
#'   sample size that will achieve more accurate statistical power for a planned
#'   study, when compared to approaches that use a sample effect size at face
#'   value or rely on sample size only. The bias and uncertainty adjusted
#'   previous study noncentrality parameter is also returned, which can be
#'   transformed to various effect size metrics.
#'
#'   For the common one-way and two-way omnibus designs,
#'   \code{\link{ss_buc_one_way_anova}} and \code{\link{ss_buc_factorial_anova}}
#'   take the factor levels directly; this general form instead takes the
#'   degrees of freedom, so it covers more factors and single-degree-of-freedom
#'   contrasts.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss_buc_factorial_anova_general} uses the
#'   observed \emph{F} value and sample size from a previous study to correct the
#'   noncentrality parameter associated with the effect of interest for
#'   publication bias and/or uncertainty. This new estimated noncentrality
#'   parameter is then used to calculate the necessary per-cell sample size to
#'   achieve the desired level of power in the planned study.
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
#'   The denominator degrees of freedom you supply describe the prior study's
#'   \emph{F} test for the effect of interest. In a pure cell-means model the
#'   residual degrees of freedom equal \emph{N} minus the number of cells;
#'   supplying a smaller value accommodates additional estimated parameters,
#'   such as covariates or blocking factors. The difference,
#'   \code{(N - cells) - df_denominator}, is treated as the number of such
#'   nuisance parameters and is carried through to the planned study, whose
#'   denominator degrees of freedom are reduced by the same amount. Supplying
#'   \code{df_denominator = N - cells} reproduces the case with no such extra
#'   parameters, and a smaller \code{df_denominator} calls for a larger sample
#'   size. \code{df_denominator} cannot exceed \code{N - cells}.
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
#'   \code{ss_buc_factorial_anova_general} assumes that the planned study will
#'   have equal \emph{n}. Unequal \emph{n} in the previous study is handled in
#'   the following way for between-subjects ANOVA designs. If the user enters an
#'   \emph{N} not equally divisible by the number of cells, the function
#'   calculates \emph{n} by dividing \emph{N} by the number of cells and both
#'   rounding up and rounding down the result, effectively assuming equal
#'   \emph{n}. The suggested sample size for the planned study is calculated
#'   using both of these values of \emph{n}, and the function returns the larger
#'   of these two suggestions, to be conservative. The adjusted noncentrality
#'   parameter that is output is the lower of the two possibilities, again, to
#'   be conservative. Although equal-n previous studies are preferable, this
#'   approach will work well as long as the cell sizes are only slightly
#'   discrepant.
#'
#' @param F_observed Observed \emph{F} value from a previous study used to plan
#'   sample size for a planned study.
#' @param N Total sample size of the previous study.
#' @param cells Number of cells for the design (the product of the number of
#'   levels of each factor).
#' @param df_numerator Numerator degrees of freedom for the effect of interest.
#' @param df_denominator Denominator (error) degrees of freedom for the effect
#'   of interest in the prior study. This equals \code{N - cells} when the cell
#'   means are the only estimated parameters; supply a smaller value to account
#'   for covariates or other nuisance parameters. It cannot exceed
#'   \code{N - cells}.
#' @template planning-params
#'
#' @templateVar size_phrase per-cell sample size
#' @template return
#'
#' @export
#'
#' @examples
#' # With df_denominator = N - cells the general form reproduces the two-way
#' # omnibus result from ss_buc_factorial_anova (here 120 - 6 = 114).
#' result <- ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
#'   df_numerator = 2, df_denominator = 114, alpha_prior = .05,
#'   alpha_planned = .05, assurance = .80, power = .80)
#' result
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .83):
#' try(ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
#'   df_numerator = 2, df_denominator = 114, assurance = .95))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
ss_buc_factorial_anova_general <- function(F_observed, N, cells, df_numerator,
                                           df_denominator, alpha_prior = .05,
                                           alpha_planned = .05, assurance = .80,
                                           power = .80) {
  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  power <- v$power

  if (missing(N)) stop("You must specify 'N', which is the total sample size.")
  if (missing(cells)) stop("You must specify 'cells', the number of cells (groups) in the between-subjects design.")
  if (missing(df_numerator)) stop("You must specify 'df_numerator', the numerator degrees of freedom for the effect of interest.")

  if (missing(df_denominator)) stop("You must specify 'df_denominator', the denominator (error) degrees of freedom for the effect of interest in the prior study.")
  .check_scalar_finite(F_observed, "F_observed")
  .check_count(N, "N", min = 2)
  .check_count(cells, "cells", min = 2)
  .check_count(df_numerator, "df_numerator", min = 1)
  if (df_denominator <= 0) stop("'df_denominator' must be a positive number.")
  .check_count(df_denominator, "df_denominator", min = 1)

  nuisance_df <- (N - cells) - df_denominator
  if (nuisance_df < 0) stop("'df_denominator' cannot exceed 'N' - 'cells', which is the residual degrees of freedom when the cell means are the only estimated parameters. Check 'N', 'cells', and 'df_denominator'.")
  if ((floor(N / cells) - 1) * cells - nuisance_df < 1) stop("Your prior study 'N' is too small for this design: after rounding the per-cell sample size down, no error degrees of freedom remain. Increase 'N', or reduce the nuisance parameters implied by 'df_denominator'.")

  # Smallest planned per-cell n whose denominator df is still positive once the
  # nuisance parameters are carried through. Reduces to 2 when nuisance_df = 0.
  n_rep_start <- max(2, ceiling((cells + nuisance_df + 1) / cells))

  ## Rounding up
  n_ru <- ceiling(N / cells)
  N_ru <- n_ru * cells
  df_denominator_ru <- (n_ru * cells) - cells - nuisance_df

  crit_F_ru <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator_ru)
  if (F_observed <= crit_F_ru) stop("Your observed F statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 'F_observed' exceeds the critical value.")

  solution_ru <- .solve_ncp_assurance(
    .tm_f(F_observed, crit_F_ru, df_numerator, df_denominator_ru), assurance)
  ncp_ru <- solution_ru$ncp

  if (ncp_ru == 0) .stop_zero_ncp(solution_ru$ceiling)

  power_at <- function(n_rep, n_prior, ncp_b) {
    denom_df <- (n_rep * cells) - cells - nuisance_df
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = denom_df)
    1 - pf(critical_F, df1 = df_numerator, df2 = denom_df, ncp = (n_rep / n_prior) * ncp_b)
  }
  repn_ru <- .smallest_n_for_power(function(k) power_at(k, n_ru, ncp_ru), power,
                                   start = n_rep_start)

  ## Rounding down
  n_rd <- floor(N / cells)
  N_rd <- n_rd * cells
  df_denominator_rd <- (n_rd * cells) - cells - nuisance_df

  crit_F_rd <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator_rd)
  if (F_observed <= crit_F_rd) stop("Your observed F statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 'F_observed' exceeds the critical value.")

  solution_rd <- .solve_ncp_assurance(
    .tm_f(F_observed, crit_F_rd, df_numerator, df_denominator_rd), assurance)
  ncp_rd <- solution_rd$ncp

  if (ncp_rd == 0) .stop_zero_ncp(solution_rd$ceiling)

  repn_rd <- .smallest_n_for_power(function(k) power_at(k, n_rd, ncp_rd), power,
                                   start = n_rep_start)

  output_n <- max(repn_ru, repn_rd)
  df_error <- output_n * cells - cells - nuisance_df
  .bucss_power_result(
    sample_size = output_n,
    df_effect = df_numerator,
    df_error = df_error,
    ncp = min(ncp_rd, ncp_ru),
    design = "Between-subjects ANOVA (general)",
    sample_size_unit = "per cell",
    assurance_ceiling = min(solution_ru$ceiling, solution_rd$ceiling),
    total_n = output_n * cells,
    inputs = list(F_observed = F_observed, N = N, cells = cells,
                  df_numerator = df_numerator, df_denominator = df_denominator,
                  alpha_prior = alpha_prior_input, alpha_planned = alpha_planned,
                  assurance = assurance, power = power)
  )
}
