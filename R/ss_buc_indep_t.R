#' Necessary sample size to reach desired power for an independent t test (or
#' standardized mean difference) using an uncertainty and publication bias
#' correction procedure
#'
#' @description \code{ss_buc_indep_t} returns the necessary per-group sample size
#'   to achieve a desired level of statistical power for a planned study using
#'   an independent \emph{t} test, based on information obtained from a previous
#'   study. The effect from the previous study can be corrected for publication
#'   bias and/or uncertainty to provide a sample size that will achieve more
#'   accurate statistical power for a planned study, when compared to approaches
#'   that use a sample effect size at face value or rely on sample size only.
#'   The bias and uncertainty adjusted previous study noncentrality parameter is
#'   also returned, which can be transformed to various effect size metrics.
#'
#'   \code{ss_buc_smd} is an alias for \code{ss_buc_indep_t}: the two are the
#'   same function. The two-group \emph{t} test and the standardized mean
#'   difference (Cohen's \eqn{d}) describe the same comparison, so the same
#'   planner serves both framings. Prefer the test-specific name
#'   (\code{ss_buc_indep_t}) when the prior \emph{t} was computed from raw
#'   (unstandardized) data, since the planner works directly from the observed
#'   \emph{t}.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss_buc_indep_t} uses the observed
#'   \emph{t} value and sample size from a previous study to correct the
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
#'   specifications. Two alternatives are recommended. First, users can select a
#'   lower value of assurance (e.g., .8 instead of .95). Second, users can
#'   reduce the influence of publication bias by setting \code{alpha_prior} at a
#'   value greater than .05. It is possible to correct for uncertainty only by
#'   setting \code{alpha_prior = 1} and choosing the desired level of assurance.
#'   We encourage users to make the adjustments as minimal as possible.
#'
#'   If you are working from a standardized mean difference (Cohen's \eqn{d})
#'   rather than a \emph{t} statistic, convert it before calling: for two
#'   independent groups \eqn{t = d / \sqrt{1/n_1 + 1/n_2}}, which equals
#'   \eqn{d\sqrt{n/2}} when the per-group sample sizes are equal.
#'
#'   \code{ss_buc_indep_t} assumes that the planned study will have equal
#'   \emph{n}. Unequal \emph{n} in the previous study is handled in the
#'   following way for the independent \emph{t}. If the user enters an odd value
#'   for \emph{N}, no information is available on the exact group sizes. The
#'   function calculates \emph{n} by dividing \emph{N} by 2 and both rounding up
#'   and rounding down the result, thus assuming equal \emph{n}. The suggested
#'   sample size for the planned study is calculated using both of these values
#'   of \emph{n}, and the function returns the larger of these two suggestions,
#'   to be conservative. If the user enters a vector for \emph{n} with two
#'   different values, specific information is available on the exact group
#'   sizes. \emph{n} is calculated as the harmonic mean of these two values (a
#'   measure of effective sample size). Again, this is rounded both up and down.
#'   The suggested sample size for the planned study is calculated using both of
#'   these values of \emph{n}, and the function returns the larger of these two
#'   suggestions, to be conservative. The adjusted noncentrality parameter that
#'   is output is the lower of the two possibilities, again, to be conservative.
#'   When the individual group sizes of an unequal-\emph{n} previous study are
#'   known, we highly encourage entering these explicitly, especially if the
#'   sample sizes are quite discrepant, as this allows for the greatest
#'   precision in estimating an appropriate planned study \emph{n}.
#'
#' @param t_observed Observed \emph{t} value from a previous study used to plan
#'   sample size for a planned study.
#' @param n Per group sample size (if equal) or the two group sample sizes of
#'   the previous study (enter either a single number or a vector of length 2).
#' @param N Total sample size of the previous study, assumed equal across groups
#'   if specified.
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
#'   rows are \code{necessary_sample_size} (the suggested per-group sample size
#'   for the planned study) and \code{ncp_adjusted} (the publication bias and
#'   uncertainty adjusted prior study noncentrality parameter). The design, the
#'   unit the sample size is counted in, and the planning inputs travel on
#'   attributes and are shown by \code{print.bucss_power}.
#'
#' @export
#'
#' @examples
#' result <- ss_buc_indep_t(t_observed = 3, n = 20, alpha_prior = .05,
#'   alpha_planned = .05, assurance = .80, power = .80)
#' result
#' result$value[result$term == "necessary_sample_size"]
#'
#' # ss_buc_smd is the same function under an effect-size name
#' ss_buc_smd(t_observed = 3, n = 20)
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @references Anderson, S. F., & Maxwell, S. E. (2017). Addressing the
#'   'replication crisis': Using original studies to design replication studies
#'   with appropriate statistical power. \emph{Multivariate Behavioral Research,
#'   52,} 305--322.
#'
#'   Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample size planning
#'   for more accurate statistical power: A method correcting sample effect
#'   sizes for uncertainty and publication bias. \emph{Psychological Science,
#'   28,} 1547--1562. \doi{10.1177/0956797617723724}
#'
#'   Taylor, D. J., & Muller, K. E. (1996). Bias in linear model power and
#'   sample size calculation due to estimating noncentrality. \emph{Communications
#'   in Statistics: Theory and Methods, 25,} 1595--1610.
ss_buc_indep_t <- function(t_observed, n, N, alpha_prior = .05, alpha_planned = .05,
                           assurance = .80, power = .80, step = .001) {
  if (alpha_prior > 1 | alpha_prior <= 0) stop("There is a problem with 'alpha_prior' of the prior study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).")
  alpha_prior_input <- alpha_prior
  if (alpha_prior == 1) alpha_prior <- .999

  if (alpha_planned >= 1 | alpha_planned <= 0) stop("There is a problem with 'alpha_planned' of the planned study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).")

  if (assurance >= 1) assurance <- assurance / 100
  if (assurance < 0 | assurance > 1) stop("There is a problem with 'assurance' (i.e., the proportion of times statistical power is at or above the desired value), please specify as a value between 0 and 1 (the default is .80).")
  if (assurance < .5) warning("The 'assurance' you have entered is < .5, which implies you will have under a 50% chance at achieving your desired level of power.")

  if (power >= 1) power <- power / 100
  if (power < 0 | power > 1) stop("There is a problem with 'power' (i.e., desired statistical power), please specify as a value between 0 and 1 (the default is .80).")

  if (!missing(N)) {
    if (!is.null(N)) {
      if (N <= 2) stop("Your total sample size is too small.")
      if (!missing(n)) stop("Because you specified 'N' you should not specify 'n'.")
      if (missing(n)) {
        n_ru <- ceiling(N / 2)
        N_ru <- 2 * n_ru
        DF_ru <- N_ru - 2

        n_rd <- floor(N / 2)
        N_rd <- 2 * n_rd
        DF_rd <- N_rd - 2
      }
    }
  }

  if (missing(N)) {
    if (!(length(n) %in% c(1, 2))) stop("The value of 'n' should be a vector of length two or a single value (for equal group sample sizes).")
    if (length(n) == 2) {
      n_1 <- n[1]
      n_2 <- n[2]
      n_ru <- ceiling(2 / ((1 / n_1) + (1 / n_2)))
      N_ru <- 2 * n_ru
      DF_ru <- N_ru - 2
      n_rd <- floor(2 / ((1 / n_1) + (1 / n_2)))
      N_rd <- 2 * n_rd
      DF_rd <- N_rd - 2
    }
    if (length(n) == 1) {
      n_ru <- n
      N_ru <- 2 * n_ru
      DF_ru <- N_ru - 2
      n_rd <- n
      N_rd <- 2 * n_rd
      DF_rd <- N_rd - 2
    }
  }

  NCP <- seq(from = 0, to = 100, by = step)

  ## Rounding up
  value_critical_ru <- qt(1 - alpha_prior / 2, df = DF_ru)
  if (t_observed <= value_critical_ru) stop("Your observed t statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 't_observed' exceeds the critical value.")

  area_above_critical_value_ru <- 1 - pt(value_critical_ru, df = DF_ru, ncp = NCP)
  area_other_tail_ru <- pt(-1 * value_critical_ru, df = DF_ru, ncp = NCP)
  power_values_ru <- area_above_critical_value_ru + area_other_tail_ru
  area_above_t_ru <- 1 - pt(t_observed, df = DF_ru, ncp = NCP)
  area_above_t_opp_ru <- pt(-1 * t_observed, df = DF_ru, ncp = NCP)
  area_between_ru <- (area_above_critical_value_ru - area_above_t_ru) + (area_other_tail_ru - area_above_t_opp_ru)

  TM_ru <- area_between_ru / power_values_ru
  ncp_ru <- min(NCP[which(abs(TM_ru - assurance) == min(abs(TM_ru - assurance)))])

  if (ncp_ru == 0) stop("The corrected noncentrality parameter is zero. Please either choose a lower value of assurance and/or a higher value of 'alpha_prior' for the prior study (e.g., accounting for less publication bias).")

  n_rep <- 2
  denom_df <- (2 * n_rep) - 2
  diff_ru <- -1
  while (diff_ru < 0) {
    critical_t <- qt(1 - alpha_planned / 2, df = denom_df)
    powers1_ru <- 1 - pt(critical_t, df = denom_df, ncp = sqrt(n_rep / n_ru) * ncp_ru)
    powers2_ru <- pt(-1 * critical_t, df = denom_df, ncp = sqrt(n_rep / n_ru) * ncp_ru)
    powers_ru <- powers1_ru + powers2_ru
    diff_ru <- powers_ru - power
    n_rep <- n_rep + 1
    denom_df <- (2 * n_rep) - 2
  }
  repn_ru <- n_rep - 1

  ## Rounding down
  value_critical_rd <- qt(1 - alpha_prior / 2, df = DF_rd)
  if (t_observed <= value_critical_rd) stop("Your observed t statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 't_observed' exceeds the critical value.")

  area_above_critical_value_rd <- 1 - pt(value_critical_rd, df = DF_rd, ncp = NCP)
  area_other_tail_rd <- pt(-1 * value_critical_rd, df = DF_rd, ncp = NCP)
  power_values_rd <- area_above_critical_value_rd + area_other_tail_rd
  area_above_t_rd <- 1 - pt(t_observed, df = DF_rd, ncp = NCP)
  area_above_t_opp_rd <- pt(-1 * t_observed, df = DF_rd, ncp = NCP)
  area_between_rd <- (area_above_critical_value_rd - area_above_t_rd) + (area_other_tail_rd - area_above_t_opp_rd)

  TM_rd <- area_between_rd / power_values_rd
  ncp_rd <- min(NCP[which(abs(TM_rd - assurance) == min(abs(TM_rd - assurance)))])

  if (ncp_rd == 0) stop("The corrected noncentrality parameter is zero. Please either choose a lower value of assurance and/or a higher value of 'alpha_prior' for the prior study (e.g., accounting for less publication bias).")

  n_rep <- 2
  denom_df <- (2 * n_rep) - 2
  diff_rd <- -1
  while (diff_rd < 0) {
    critical_t <- qt(1 - alpha_planned / 2, df = denom_df)
    powers1_rd <- 1 - pt(critical_t, df = denom_df, ncp = sqrt(n_rep / n_rd) * ncp_rd)
    powers2_rd <- pt(-1 * critical_t, df = denom_df, ncp = sqrt(n_rep / n_rd) * ncp_rd)
    powers_rd <- powers1_rd + powers2_rd
    diff_rd <- powers_rd - power
    n_rep <- n_rep + 1
    denom_df <- (2 * n_rep) - 2
  }
  repn_rd <- n_rep - 1

  output_n <- max(repn_ru, repn_rd)
  .bucss_power_result(
    sample_size = output_n,
    ncp = min(ncp_ru, ncp_rd),
    design = "Independent t test",
    sample_size_unit = "per group",
    inputs = list(t_observed = t_observed, alpha_prior = alpha_prior_input,
                  alpha_planned = alpha_planned, assurance = assurance,
                  power = power)
  )
}

#' @rdname ss_buc_indep_t
#' @export
ss_buc_smd <- ss_buc_indep_t
