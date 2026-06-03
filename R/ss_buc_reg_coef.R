#' Necessary sample size to reach desired power for a single coefficient in a
#' multiple regression using a publication bias and uncertainty correction
#' procedure
#'
#' @description \code{ss_buc_reg_coef} returns the necessary total sample size to
#'   achieve a desired level of statistical power for a single regression
#'   coefficient in a planned study using multiple regression, based on
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
#'   The approach implemented in \code{ss_buc_reg_coef} uses the observed
#'   \emph{t} value and sample size from a previous study to correct the
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
#' @param t_observed Observed \emph{t} value from a previous study used to plan
#'   sample size for a planned study.
#' @param N Total sample size of the previous study.
#' @param p Number of predictors; be sure to include any product terms or
#'   polynomials that are in the model.
#' @param alpha_prior Alpha level \eqn{\alpha} for the previous study or the
#'   assumed statistical significance necessary for publishing in the field; to
#'   assume no publication bias, a value of 1 can be entered.
#' @param alpha_planned Alpha level \eqn{\alpha} assumed for the planned study.
#' @param assurance Desired level of assurance, or the long run proportion of
#'   times that the planned study power will reach or surpass the desired level
#'   (assurance > .5 corrects for uncertainty; assurance < .5 is not
#'   recommended).
#' @param power Desired level of statistical power for the planned study.
#' @param step Value used in the iterative scheme to determine the noncentrality
#'   parameter necessary for sample size planning (0 < step < .01). Users should
#'   not generally need to change this value; smaller values lead to more
#'   accurate sample size planning results, but unnecessarily small values will
#'   add unnecessary computational time.
#'
#' @return An object of class \code{bucss_power}: a \code{data.frame} with a
#'   character \code{term} column and a numeric \code{value} column whose two
#'   rows are \code{necessary_sample_size} (the suggested total sample size for
#'   the planned study) and \code{ncp_adjusted} (the publication bias and
#'   uncertainty adjusted prior study noncentrality parameter). The design, the
#'   unit the sample size is counted in, and the planning inputs travel on
#'   attributes and are shown by \code{print.bucss_power}.
#'
#' @export
#'
#' @examples
#' result <- ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, alpha_prior = .05,
#'   alpha_planned = .05, assurance = .80, power = .80)
#' result
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .93):
#' try(ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, assurance = .95))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @references Anderson, S. F., & Kelley, K. (2024). Sample size planning for
#'   replication studies: The devil is in the design. \emph{Psychological
#'   Methods, 29,} 844--867. \doi{10.1037/met0000520}
#'
#'   Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size planning
#'   for more accurate statistical power: A method adjusting sample effect sizes
#'   for publication bias and uncertainty. \emph{Psychological Science, 28,}
#'   1547--1562. \doi{10.1177/0956797617723724}
#'
#'   Anderson, S. F., & Maxwell, S. E. (2017). Addressing the 'replication
#'   crisis': Using original studies to design replication studies with
#'   appropriate statistical power. \emph{Multivariate Behavioral Research, 52,}
#'   305--324.
#'
#'   Taylor, D. J., & Muller, K. E. (1996). Bias in linear model power and
#'   sample size calculation due to estimating noncentrality. \emph{Communications
#'   in Statistics: Theory and Methods, 25,} 1595--1610.
ss_buc_reg_coef <- function(t_observed, N, p, alpha_prior = .05,
                          alpha_planned = .05, assurance = .80, power = .80,
                          step = .001) {
  if (alpha_prior > 1 | alpha_prior <= 0) stop("There is a problem with 'alpha_prior' of the prior study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).")
  alpha_prior_input <- alpha_prior
  if (alpha_prior == 1) alpha_prior <- .999

  if (alpha_planned >= 1 | alpha_planned <= 0) stop("There is a problem with 'alpha_planned' of the planned study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).")

  if (assurance >= 1) assurance <- assurance / 100
  if (assurance < 0 | assurance > 1) stop("There is a problem with 'assurance' (i.e., the proportion of times statistical power is at or above the desired value), please specify as a value between 0 and 1 (the default is .80).")
  if (assurance < .5) warning("The 'assurance' you have entered is < .5, which implies you will have under a 50% chance at achieving your desired level of power.")

  if (power >= 1) power <- power / 100
  if (power < 0 | power > 1) stop("There is a problem with 'power' (i.e., desired statistical power), please specify as a value between 0 and 1 (the default is .80).")

  if (missing(N)) stop("You need to specify 'N', which is the total sample size of the original study.")
  if (N <= 1) stop("Your total sample size is too small.")
  if (p < 1) stop("Your number of predictors is too small.")
  if (N - p - 1 < 1) stop("The combination of your sample size and number of predictors leads to 0 or negative degrees of freedom.")

  df_numerator <- 1
  df_denominator <- N - p - 1

  NCP <- seq(from = 0, to = 100, by = step)
  value_critical <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator)

  if (t_observed^2 <= value_critical) stop("Your observed t statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 't_observed' exceeds the critical value.")

  power_values <- 1 - pf(value_critical, df1 = df_numerator, df2 = df_denominator, ncp = NCP)
  area_above_F <- 1 - pf(t_observed^2, df1 = df_numerator, df2 = df_denominator, ncp = NCP)
  area_between <- power_values - area_above_F
  TM <- area_between / power_values

  ncp <- min(NCP[which(abs(TM - assurance) == min(abs(TM - assurance)))])

  if (ncp == 0) .stop_zero_ncp(max(TM))

  n_rep <- 2 + p + 1
  denom_df <- n_rep - p - 1
  diff <- -1
  while (diff < 0) {
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = denom_df)
    powers <- 1 - pf(critical_F, df1 = df_numerator, df2 = denom_df, ncp = (n_rep / N) * ncp)
    diff <- powers - power
    n_rep <- n_rep + 1
    denom_df <- n_rep - p - 1
  }
  output_n <- n_rep - 1

  .bucss_power_result(
    sample_size = output_n,
    ncp = ncp,
    design = "Multiple regression: single coefficient",
    sample_size_unit = "total",
    inputs = list(t_observed = t_observed, N = N, p = p,
                  alpha_prior = alpha_prior_input, alpha_planned = alpha_planned,
                  assurance = assurance, power = power)
  )
}
