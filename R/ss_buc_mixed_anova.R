#' Necessary sample size to reach desired power for a two-factor split-plot
#' (mixed) ANOVA using a publication bias and uncertainty correction procedure
#'
#' @description \code{ss_buc_mixed_anova} returns the necessary per between-subjects
#'   cell sample size to achieve a desired level of statistical power for a
#'   planned study testing an omnibus effect using a two-factor split-plot
#'   (mixed) ANOVA, based on information obtained from a previous study. The
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
#'   The approach implemented in \code{ss_buc_mixed_anova} uses the observed
#'   \emph{F} value and sample size from a previous study to correct the
#'   noncentrality parameter associated with the effect of interest for
#'   publication bias and/or uncertainty. This new estimated noncentrality
#'   parameter is then used to calculate the necessary per between-subjects cell
#'   sample size to achieve the desired level of power in the planned study.
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
#'   \code{ss_buc_mixed_anova} assumes that the planned study will have equal
#'   \emph{n}. Unequal \emph{n} in the previous study is handled in the
#'   following way for split-plot designs. If the user enters an \emph{N} not
#'   equally divisible by the number of between-subjects cells, the function
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
#'   \code{ss_buc_mixed_anova} assumes sphericity for the within-subjects effects.
#'
#' @param F_observed Observed \emph{F} value from a previous study used to plan
#'   sample size for a planned study.
#' @param N Total sample size of the previous study.
#' @param levels_between Number of levels for the between-subjects factor.
#' @param levels_within Number of levels for the within-subjects factor.
#' @param effect Effect most of interest to the planned study: between-subjects
#'   main effect (\code{between}), within-subjects main effect (\code{within}),
#'   or interaction (\code{interaction}).
#' @template planning-params
#'
#' @templateVar size_phrase per between-subjects cell sample size
#' @template return-effect
#'
#' @export
#'
#' @examples
#' result <- ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
#'   levels_within = 3, effect = "within", alpha_prior = .05,
#'   alpha_planned = .05, assurance = .80, power = .80)
#' result
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .83):
#' try(ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
#'   levels_within = 3, effect = "within", assurance = .95))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
ss_buc_mixed_anova <- function(F_observed, N, levels_between, levels_within,
                         effect = c("between", "within", "interaction"),
                         alpha_prior = .05, alpha_planned = .05,
                         assurance = .80, power = .80, step = .001) {
  effect <- match.arg(effect)

  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, power, step)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  power <- v$power

  if (missing(N)) stop("You must specify 'N', which is the total sample size.")
  if (missing(levels_within)) stop("You must specify the number of levels of the within-subjects factor. If there is no within-subjects factor, use the between-subjects approach.")
  if (missing(levels_between)) stop("You must specify the number of levels of the between-subjects factor. If there is no between-subjects factor, use the within-subjects approach.")
  if (N < 2 * levels_between) stop("Your prior study 'N' is too small for this design: at least two observations per between-subjects cell are required, so 'N' must be at least 2 * 'levels_between'.")

  NCP <- seq(from = 0, to = 100, by = step)

  ## Rounding down
  n_rd <- floor(N / levels_between)
  N_rd <- n_rd * levels_between
  if (effect == "between") {
    df_numerator <- levels_between - 1
    df_denominator_rd <- N_rd - levels_between
  } else if (effect == "within") {
    df_numerator <- levels_within - 1
    df_denominator_rd <- (N_rd - levels_between) * (levels_within - 1)
  } else {
    df_numerator <- (levels_between - 1) * (levels_within - 1)
    df_denominator_rd <- (N_rd - levels_between) * (levels_within - 1)
  }

  crit_F_rd <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator_rd)
  if (F_observed <= crit_F_rd) stop("Your observed F statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 'F_observed' exceeds the critical value.")

  power_values_rd <- 1 - pf(crit_F_rd, df1 = df_numerator, df2 = df_denominator_rd, ncp = NCP)
  area_above_F_rd <- 1 - pf(F_observed, df1 = df_numerator, df2 = df_denominator_rd, ncp = NCP)
  area_between_rd <- power_values_rd - area_above_F_rd

  TM_rd <- area_between_rd / power_values_rd
  ncp_rd <- min(NCP[which(abs(TM_rd - assurance) == min(abs(TM_rd - assurance)))])

  if (ncp_rd == 0) .stop_zero_ncp(max(TM_rd))

  n_rep <- 2
  denom_df <- if (effect == "between") n_rep * levels_between - levels_between else levels_between * (n_rep - 1) * (levels_within - 1)
  diff_rd <- -1
  while (diff_rd < 0) {
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = denom_df)
    powers_rd <- 1 - pf(critical_F, df1 = df_numerator, df2 = denom_df, ncp = (n_rep / n_rd) * ncp_rd)
    diff_rd <- powers_rd - power
    n_rep <- n_rep + 1
    denom_df <- if (effect == "between") n_rep * levels_between - levels_between else levels_between * (n_rep - 1) * (levels_within - 1)
  }
  repn_rd <- n_rep - 1

  ## Rounding up
  n_ru <- ceiling(N / levels_between)
  N_ru <- n_ru * levels_between
  if (effect == "between") {
    df_numerator <- levels_between - 1
    df_denominator_ru <- N_ru - levels_between
  } else if (effect == "within") {
    df_numerator <- levels_within - 1
    df_denominator_ru <- (N_ru - levels_between) * (levels_within - 1)
  } else {
    df_numerator <- (levels_between - 1) * (levels_within - 1)
    df_denominator_ru <- (N_ru - levels_between) * (levels_within - 1)
  }

  crit_F_ru <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator_ru)
  if (F_observed <= crit_F_ru) stop("Your observed F statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 'F_observed' exceeds the critical value.")

  power_values_ru <- 1 - pf(crit_F_ru, df1 = df_numerator, df2 = df_denominator_ru, ncp = NCP)
  area_above_F_ru <- 1 - pf(F_observed, df1 = df_numerator, df2 = df_denominator_ru, ncp = NCP)
  area_between_ru <- power_values_ru - area_above_F_ru

  TM_ru <- area_between_ru / power_values_ru
  ncp_ru <- min(NCP[which(abs(TM_ru - assurance) == min(abs(TM_ru - assurance)))])

  if (ncp_ru == 0) .stop_zero_ncp(max(TM_ru))

  n_rep <- 2
  denom_df <- if (effect == "between") n_rep * levels_between - levels_between else levels_between * (n_rep - 1) * (levels_within - 1)
  diff_ru <- -1
  while (diff_ru < 0) {
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = denom_df)
    powers_ru <- 1 - pf(critical_F, df1 = df_numerator, df2 = denom_df, ncp = (n_rep / n_ru) * ncp_ru)
    diff_ru <- powers_ru - power
    n_rep <- n_rep + 1
    denom_df <- if (effect == "between") n_rep * levels_between - levels_between else levels_between * (n_rep - 1) * (levels_within - 1)
  }
  repn_ru <- n_rep - 1

  output_n <- max(repn_rd, repn_ru)
  .bucss_power_result(
    sample_size = output_n,
    ncp = min(ncp_rd, ncp_ru),
    design = "Two-factor split-plot (mixed) ANOVA",
    sample_size_unit = "per between-subjects cell",
    assurance_ceiling = min(max(TM_ru), max(TM_rd)),
    total_n = output_n * levels_between,
    effect = effect,
    inputs = list(F_observed = F_observed, N = N, levels_between = levels_between,
                  levels_within = levels_within, alpha_prior = alpha_prior_input,
                  alpha_planned = alpha_planned, assurance = assurance,
                  power = power)
  )
}
