#' Necessary sample size to reach desired power for a split-plot (mixed) ANOVA
#' with any number of factors using a publication bias and uncertainty
#' correction procedure
#'
#' @description \code{ss_buc_mixed_anova_general} returns the necessary per-group
#'   sample size to achieve a desired level of statistical power for a planned
#'   study testing any type of effect (omnibus, contrast) using a split-plot
#'   (mixed) ANOVA with any number of factors, based on information obtained
#'   from a previous study. The effect from the previous study can be corrected
#'   for publication bias and/or uncertainty to provide a sample size that will
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
#'   The approach implemented in \code{ss_buc_mixed_anova_general} uses the observed
#'   \emph{F} value and sample size from a previous study to correct the
#'   noncentrality parameter associated with the effect of interest for
#'   publication bias and/or uncertainty. This new estimated noncentrality
#'   parameter is then used to calculate the necessary per-group sample size to
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
#'   \code{ss_buc_mixed_anova_general} assumes that the planned study will have equal
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
#'   \code{ss_buc_mixed_anova_general} assumes sphericity for the within-subjects
#'   effects.
#'
#' @param F_observed Observed \emph{F} value from a previous study used to plan
#'   sample size for a planned study.
#' @param N Total sample size of the previous study.
#' @param df_numerator Numerator degrees of freedom for the effect of interest.
#' @param num_groups Number of distinct groups (product of the number of levels
#'   of the between-subjects factors).
#' @param effect Effect of interest: involves only between-subjects effects
#'   (\code{between_only}), involves only within-subjects effects
#'   (\code{within_only}), or involves both between and within effects
#'   (\code{between_within}).
#' @param df_num_within Numerator degrees of freedom for only the
#'   within-subjects components of the effect of interest. Only needed when
#'   \code{effect = "between_within"}.
#' @template planning-params
#'
#' @return An object of class \code{bucss_power}: a \code{data.frame} with a
#'   character \code{term} column and a numeric \code{value} column whose two
#'   rows are \code{necessary_sample_size} (the suggested per-group sample size
#'   for the planned study) and \code{ncp_adjusted} (the publication bias and
#'   uncertainty adjusted prior study noncentrality parameter). The design, the
#'   effect tested, the unit the sample size is counted in, and the planning
#'   inputs travel on attributes and are shown by \code{print.bucss_power}.
#'
#' @export
#'
#' @examples
#' result <- ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
#'   num_groups = 3, effect = "between_only", df_num_within = 3,
#'   alpha_prior = .05, alpha_planned = .05, assurance = .80, power = .80)
#' result
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .82):
#' try(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
#'   num_groups = 3, effect = "between_only", df_num_within = 3,
#'   assurance = .95))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
ss_buc_mixed_anova_general <- function(F_observed, N, df_numerator, num_groups,
                                 effect = c("between_only", "within_only", "between_within"),
                                 df_num_within, alpha_prior = .05,
                                 alpha_planned = .05, assurance = .80,
                                 power = .80, step = .001) {
  effect <- match.arg(effect)

  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, power, step)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  power <- v$power

  if (missing(N)) stop("You must specify 'N', which is the total sample size.")
  if (missing(df_numerator)) stop("You must specify 'df_numerator', the numerator degrees of freedom for the effect of interest.")
  if (missing(num_groups)) stop("You must specify 'num_groups', the number of between-subjects groups.")
  if (effect == "between_within" && missing(df_num_within)) stop("For 'effect = \"between_within\"' you must specify 'df_num_within', the numerator degrees of freedom of the within-subjects component.")
  if (N < 2 * num_groups) stop("Your prior study 'N' is too small for this design: at least two observations per between-subjects group are required, so 'N' must be at least 2 * 'num_groups'.")

  NCP <- seq(from = 0, to = 100, by = step)

  denom_mult <- switch(effect,
                       between_only = 1,
                       within_only = df_numerator,
                       between_within = df_num_within)

  ## Rounding down
  n_rd <- floor(N / num_groups)
  N_rd <- n_rd * num_groups
  df_denominator_rd <- (N_rd - num_groups) * denom_mult

  crit_F_rd <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator_rd)
  if (F_observed <= crit_F_rd) stop("Your observed F statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 'F_observed' exceeds the critical value.")

  power_values_rd <- 1 - pf(crit_F_rd, df1 = df_numerator, df2 = df_denominator_rd, ncp = NCP)
  area_above_F_rd <- 1 - pf(F_observed, df1 = df_numerator, df2 = df_denominator_rd, ncp = NCP)
  area_between_rd <- power_values_rd - area_above_F_rd

  TM_rd <- area_between_rd / power_values_rd
  ncp_rd <- min(NCP[which(abs(TM_rd - assurance) == min(abs(TM_rd - assurance)))])

  if (ncp_rd == 0) .stop_zero_ncp(max(TM_rd))

  n_rep <- 2
  denom_df <- ((n_rep * num_groups) - num_groups) * denom_mult
  diff_rd <- -1
  while (diff_rd < 0) {
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = denom_df)
    powers_rd <- 1 - pf(critical_F, df1 = df_numerator, df2 = denom_df, ncp = (n_rep / n_rd) * ncp_rd)
    diff_rd <- powers_rd - power
    n_rep <- n_rep + 1
    denom_df <- ((n_rep * num_groups) - num_groups) * denom_mult
  }
  repn_rd <- n_rep - 1

  ## Rounding up
  n_ru <- ceiling(N / num_groups)
  N_ru <- n_ru * num_groups
  df_denominator_ru <- (N_ru - num_groups) * denom_mult

  crit_F_ru <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator_ru)
  if (F_observed <= crit_F_ru) stop("Your observed F statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 'F_observed' exceeds the critical value.")

  power_values_ru <- 1 - pf(crit_F_ru, df1 = df_numerator, df2 = df_denominator_ru, ncp = NCP)
  area_above_F_ru <- 1 - pf(F_observed, df1 = df_numerator, df2 = df_denominator_ru, ncp = NCP)
  area_between_ru <- power_values_ru - area_above_F_ru

  TM_ru <- area_between_ru / power_values_ru
  ncp_ru <- min(NCP[which(abs(TM_ru - assurance) == min(abs(TM_ru - assurance)))])

  if (ncp_ru == 0) .stop_zero_ncp(max(TM_ru))

  n_rep <- 2
  denom_df <- ((n_rep * num_groups) - num_groups) * denom_mult
  diff_ru <- -1
  while (diff_ru < 0) {
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = denom_df)
    powers_ru <- 1 - pf(critical_F, df1 = df_numerator, df2 = denom_df, ncp = (n_rep / n_ru) * ncp_ru)
    diff_ru <- powers_ru - power
    n_rep <- n_rep + 1
    denom_df <- ((n_rep * num_groups) - num_groups) * denom_mult
  }
  repn_ru <- n_rep - 1

  inputs <- list(F_observed = F_observed, N = N, df_numerator = df_numerator,
                 num_groups = num_groups, alpha_prior = alpha_prior_input,
                 alpha_planned = alpha_planned, assurance = assurance,
                 power = power)
  if (effect == "between_within") inputs$df_num_within <- df_num_within

  output_n <- max(repn_rd, repn_ru)
  .bucss_power_result(
    sample_size = output_n,
    ncp = min(ncp_rd, ncp_ru),
    design = "Split-plot (mixed) ANOVA (any number of factors)",
    sample_size_unit = "per group",
    assurance_ceiling = min(max(TM_ru), max(TM_rd)),
    total_n = output_n * num_groups,
    effect = effect,
    inputs = inputs
  )
}
