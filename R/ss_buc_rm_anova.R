#' Necessary sample size to reach desired power for a one or two-way
#' within-subjects ANOVA using a publication bias and uncertainty correction
#' procedure
#'
#' @description \code{ss_buc_rm_anova} returns the necessary total sample size to
#'   achieve a desired level of statistical power for a planned study testing an
#'   omnibus effect using a one or two-way fully within-subjects ANOVA, based on
#'   information obtained from a previous study. The effect from the previous
#'   study can be corrected for publication bias and/or uncertainty to provide a
#'   sample size that will achieve more accurate statistical power for a planned
#'   study, when compared to approaches that use a sample effect size at face
#'   value or rely on sample size only. The bias and uncertainty adjusted
#'   previous study noncentrality parameter is also returned, which can be
#'   transformed to various effect size metrics.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss_buc_rm_anova} uses the observed
#'   \emph{F} value and sample size from a previous study to correct the
#'   noncentrality parameter associated with the effect of interest for
#'   publication bias and/or uncertainty. This new estimated noncentrality
#'   parameter is then used to calculate the necessary total sample size to
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
#'   \code{ss_buc_rm_anova} assumes sphericity for the within-subjects effects.
#'
#' @param F_observed Observed \emph{F} value from a previous study used to plan
#'   sample size for a planned study.
#' @param N Total sample size of the previous study.
#' @param levels_A Number of levels for Factor A.
#' @param levels_B Number of levels for Factor B, which is NULL if a single
#'   factor design.
#' @param effect Effect most of interest to the planned study: main effect of A
#'   (\code{factor_A}), main effect of B (\code{factor_B}), or interaction
#'   (\code{interaction}).
#' @template planning-params
#'
#' @templateVar size_phrase total sample size
#' @template return-effect
#'
#' @export
#'
#' @examples
#' result <- ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3,
#'   effect = "factor_B", alpha_prior = .05, alpha_planned = .05,
#'   assurance = .80, power = .80)
#' result
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .83):
#' try(ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3,
#'   effect = "factor_B", assurance = .95))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
ss_buc_rm_anova <- function(F_observed, N, levels_A, levels_B = NULL,
                        effect = c("factor_A", "factor_B", "interaction"),
                        alpha_prior = .05, alpha_planned = .05,
                        assurance = .80, power = .80, step = .001) {
  effect <- match.arg(effect)

  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, power, step)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  power <- v$power

  if (missing(N)) stop("You must specify 'N', which is the total sample size.")

  if (is.null(levels_B) && effect != "factor_A") stop("For a one-way within-subjects design ('levels_B' not supplied), the only testable effect is the single factor; keep the default 'effect = \"factor_A\"'. To test 'factor_B' or the interaction, supply 'levels_B'.")

  n <- N
  NCP <- seq(from = 0, to = 100, by = step)

  if (is.null(levels_B)) {
    df_numerator <- levels_A - 1
  } else {
    df_numerator <- switch(effect,
                           factor_A = levels_A - 1,
                           factor_B = levels_B - 1,
                           interaction = (levels_A - 1) * (levels_B - 1))
  }
  df_denominator <- df_numerator * (n - 1)

  crit_F <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator)
  if (F_observed <= crit_F) stop("Your observed F statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 'F_observed' exceeds the critical value.")

  power_values <- 1 - pf(crit_F, df1 = df_numerator, df2 = df_denominator, ncp = NCP)
  area_above_F <- 1 - pf(F_observed, df1 = df_numerator, df2 = df_denominator, ncp = NCP)
  area_between <- power_values - area_above_F

  TM <- area_between / power_values
  ncp <- min(NCP[which(abs(TM - assurance) == min(abs(TM - assurance)))])

  if (ncp == 0) .stop_zero_ncp(max(TM))

  n_rep <- 2
  denom_df <- df_numerator * (n_rep - 1)
  diff <- -1
  while (diff < 0) {
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = denom_df)
    powers <- 1 - pf(critical_F, df1 = df_numerator, df2 = denom_df, ncp = (n_rep / n) * ncp)
    diff <- powers - power
    n_rep <- n_rep + 1
    denom_df <- df_numerator * (n_rep - 1)
  }
  output_n <- n_rep - 1

  .bucss_power_result(
    sample_size = output_n,
    ncp = ncp,
    design = "One or two-way within-subjects ANOVA",
    sample_size_unit = "total",
    assurance_ceiling = max(TM),
    effect = effect,
    inputs = list(F_observed = F_observed, N = N, levels_A = levels_A,
                  levels_B = levels_B, alpha_prior = alpha_prior_input,
                  alpha_planned = alpha_planned, assurance = assurance,
                  power = power)
  )
}
