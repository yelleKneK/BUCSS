#' Necessary Sample Size to Reach Desired Power for a Joint Test of Multiple
#' Predictors in a Multiple Regression Using a Publication Bias and Uncertainty
#' Correction Procedure
#'
#' @description \code{ss_buc_reg_joint} returns the necessary total sample
#'   size to achieve a desired level of statistical power for a joint test of
#'   multiple predictors in a planned study using multiple regression, based on
#'   information obtained from a previous study. The effect from the previous
#'   study can be corrected for publication bias and/or uncertainty to provide a
#'   sample size that will achieve more accurate statistical power for a planned
#'   study, when compared to approaches that use a sample effect size at face
#'   value or rely on sample size only. The bias and uncertainty adjusted
#'   previous study noncentrality parameter is also returned, together with the
#'   planning effect size, the effect size the plan is built on, as Cohen's
#'   \eqn{f^2}{f^2}.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss_buc_reg_joint} uses the observed
#'   \emph{F} value and sample size from a previous study to correct the
#'   noncentrality parameter associated with the effect of interest for
#'   publication bias and/or uncertainty. This new estimated noncentrality
#'   parameter is then used to calculate the necessary total sample size to
#'   achieve the desired level of power in the planned study.
#'
#'   \strong{Fixed or random predictors.} Whether the predictors' values are
#'   fixed by design, as in an experiment in which the researcher sets them,
#'   or sampled along with the outcome, as in most observational research,
#'   changes neither the estimates nor the distribution of the test statistic
#'   when the tested coefficients are zero, but it does change the statistic's
#'   distribution, and so the power, when they are not all zero (Sampson,
#'   1974; Gatsonis & Sampson, 1989; Kelley & Maxwell, 2003). With
#'   \code{random_predictors = FALSE}, the default, the prior study's statistic
#'   and the planned study's power come from the noncentral \emph{F}
#'   distribution, which is exact when the predictor values are fixed; this is
#'   the frame of BUCSS 1.x and of Anderson (2021). With
#'   \code{random_predictors = TRUE} the noncentrality parameter is itself
#'   random: it is proportional to the sampled spread of the tested predictors
#'   once the other predictors are accounted for, a chi square variable on
#'   \emph{N} - 1 - (\emph{p} - \code{p_joint}) degrees of freedom, so the
#'   statistic follows a chi square mixture of noncentral \emph{F}
#'   distributions, which is more dispersed. The mixture has an exact closed
#'   form (the Poisson weights of the noncentral \emph{F} become negative
#'   binomial weights), used for both the correction and the planned study's
#'   power. It is exact when the tested predictors are multivariate normal
#'   around their linear regression on the other predictors, which may
#'   themselves be fixed by design (dummy codes for assigned groups, say). In
#'   observational research in the behavioral sciences the predictors are
#'   typically random (Kelley & Maxwell, 2003), and
#'   \code{random_predictors = TRUE} is then the frame to plan in.
#'
#'   Treating sampled predictors as fixed overstates the planned study's
#'   power, by about 1/\emph{n} to 1.4/\emph{n} for a plan of \emph{n}
#'   participants when every predictor is tested, and more the more
#'   predictors are left untested: 2/\emph{n} for two of four and
#'   3.8/\emph{n} for two of eight (at \emph{n} = 100), so the assurance the
#'   plans deliver can fall short of the requested \code{assurance}. The
#'   delivered assurance is what \code{assurance} promises: the proportion of
#'   plans, each made from its own published prior study, whose power at the
#'   true effect reaches \code{desired_power}, counting a refusal to plan as
#'   meeting the target. In simulations with multivariate normal predictors,
#'   an \code{assurance} of .80, and prior studies of 40 to 100 participants
#'   with power .3 to .8 at the true effect, testing two of four or five
#'   predictors or three of six, the fixed frame delivered about .76 to .81,
#'   the lower values from smaller, better powered prior studies;
#'   \code{random_predictors = TRUE} delivered .79 to .81 in the same
#'   configurations.
#'
#'   Which frame plans the larger study depends on the untested predictors
#'   and on how far the plan outgrows the prior study, because the fixed
#'   frame scales the prior study's noncentrality parameter by \emph{N} while
#'   the tested predictors' sampled spread has \emph{N} - 1 - (\emph{p} -
#'   \code{p_joint}) degrees of freedom. With every predictor tested the
#'   random-predictor plan is usually a few participants larger (and at most
#'   a few percent smaller). With predictors left untested the two plans are
#'   within a few participants, the random-predictor plan usually the larger,
#'   when the plan is no more than about twice the prior \emph{N}; beyond
#'   that the random-predictor plan is usually smaller, and the more the plan
#'   outgrows the prior study, the closer the ratio of the two plans comes to
#'   (\emph{N} - 1 - (\emph{p} - \code{p_joint}))/\emph{N}, so the
#'   random-predictor plan is smaller by at most about (1 + \emph{p} -
#'   \code{p_joint})/\emph{N}: 518 rather than 600 participants for \emph{F} =
#'   6 from a prior study of 30 testing two of six predictors, where the limit
#'   is about 17 percent, and, with two of eight predictors tested, a limit of
#'   35 percent from a prior study of 20 and 7 percent from 100. The
#'   largest assurance a prior result can support is the same in both frames,
#'   because at a zero effect both are the central \emph{F} distribution. With
#'   \code{p_joint = p} this function returns what \code{\link{ss_buc_R2}}
#'   returns, and with \code{p_joint = 1} what \code{\link{ss_buc_reg_coef}}
#'   returns for the square root of \code{F_observed}, in either frame.
#'
#'   The approach uses a likelihood function of a truncated noncentral
#'   \emph{F} distribution, where the truncation occurs due to small effect
#'   sizes being unobserved due to publication bias. The numerator of the
#'   likelihood function is the density of a noncentral \emph{F} distribution
#'   (with \code{random_predictors = TRUE}, of the mixture described above).
#'   The denominator is the power of the test, which serves to truncate the
#'   distribution. Thus, the ratio of the numerator and the denominator is a
#'   truncated noncentral \emph{F} distribution. (See Taylor & Muller, 1996,
#'   Equation 2.1; Anderson, Kelley, & Maxwell, 2017, Equation 3; and Anderson
#'   & Maxwell, 2017, for more details.)
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
#'   The search for the planned sample size starts at \emph{p} + 3, the
#'   smallest total that leaves the planned study two error degrees of
#'   freedom (as in BUCSS 1.x), so no plan is smaller than that, even when the
#'   prior effect is so large that a smaller study would already reach the
#'   desired power. A plan at this floor signals a very large prior effect,
#'   not a recommended design.
#'
#' @param F_observed Observed \emph{F} value from a previous study used to plan
#'   sample size for a planned study.
#' @param N Total sample size of the previous study.
#' @param p Number of predictors; be sure to include any product terms or
#'   polynomials that are in the model.
#' @param p_joint Number of predictors tested jointly for significance.
#' @template planning-params
#' @param random_predictors \code{TRUE} if the predictors are sampled along
#'   with the outcome, as in most observational research, and \code{FALSE} if
#'   their values are fixed by design, as in an experiment that assigns them.
#'   The default is \code{FALSE}, the frame of BUCSS 1.x; see Details. The
#'   choice is echoed in the result as a \code{random_predictors} row valued 1
#'   (\code{TRUE}) or 0 (\code{FALSE}) and named in the printed footer.
#'
#' @templateVar size_phrase total sample size
#' @templateVar ncp_scale on the scale of the \emph{F} statistic at the prior study's sample size
#' @templateVar effect_phrase \code{f2_planning}, Cohen's \eqn{f^2}{f^2} for the set of predictors tested, \eqn{\lambda/N}{lambda/N}, where \eqn{\lambda}{lambda} is \code{ncp_adjusted} and \emph{N} the prior sample size; the planned study's power is computed from it at the planned sample size, in the frame \code{random_predictors} names (with sampled predictors it is the population \eqn{f^2}{f^2} of that frame, and the power statement that follows is in that frame)
#' @template return
#'
#' @export
#'
#' @examples
#' result <- ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
#'   alpha_prior = .05, alpha_planned = .05, assurance = .80, desired_power = .80)
#' result
#'
#' # Predictors sampled along with the outcome, as in most observational
#' # research:
#' ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
#'   random_predictors = TRUE)
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .84):
#' try(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
#'   assurance = .95))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
#' @references Anderson, S. F. (2021). Using prior information to plan
#'   appropriately powered regression studies: A tutorial using BUCSS.
#'   \emph{Psychological Methods, 26,} 513--526. \doi{10.1037/met0000366}
#'
#'   Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact power
#'   and sample size calculations. \emph{Psychological Bulletin, 106,}
#'   516--524. \doi{10.1037/0033-2909.106.3.516}
#'
#'   Sampson, A. R. (1974). A tale of two regressions. \emph{Journal of the
#'   American Statistical Association, 69,} 682--689.
#'   \doi{10.1080/01621459.1974.10480189}
ss_buc_reg_joint <- function(F_observed, N, p, p_joint, alpha_prior = .05,
                               alpha_planned = .05, assurance = .80,
                               desired_power = .80, random_predictors = FALSE) {
  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, desired_power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  desired_power <- v$desired_power

  if (missing(N)) stop("You need to specify 'N', which is the total sample size of the previous study.", call. = FALSE)
  if (missing(p)) stop("You need to specify 'p', the number of predictors in the model.", call. = FALSE)
  .check_scalar_finite(F_observed, "F_observed")
  .check_count(N, "N", min = 2)
  .check_count(p, "p", min = 1)
  .check_flag(random_predictors, "random_predictors")
  if (N - p - 1 < 1) stop("The combination of your sample size and number of predictors leads to 0 or negative degrees of freedom.", call. = FALSE)
  if (missing(p_joint)) stop("You need to specify 'p_joint', the number of predictors tested jointly.", call. = FALSE)
  .check_count(p_joint, "p_joint", min = 1)
  if (p_joint > p) stop("The number of tested predictors cannot exceed the number of total predictors.", call. = FALSE)

  df_numerator <- p_joint
  df_denominator <- N - p - 1

  value_critical <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator)

  if (F_observed <= value_critical) .stop_nonsignificant("F statistic", "'F_observed'")

  ncp_solution <- .solve_ncp_assurance(
    .tm_regression(F_observed, value_critical, df_numerator, N, p,
                   random_predictors, "F statistic", "'F_observed'"), assurance)
  ncp <- ncp_solution$ncp

  if (ncp == 0) .stop_zero_ncp(ncp_solution$ceiling)

  power_at <- .power_regression(ncp, N, p, df_numerator, alpha_planned,
                                random_predictors, "F statistic", "'F_observed'")
  output_n <- .smallest_n_for_power(power_at, desired_power, start = 2 + p + 1)

  df_error <- output_n - p - 1
  actual_power <- power_at(output_n)
  .bucss_plan_result(
    sample_size = output_n,
    size_term = "necessary_N",
    effect_size_term = "f2_planning",
    effect_size = .es_f2(ncp, N),
    actual_power = actual_power,
    df_effect = df_numerator,
    df_error = df_error,
    ncp = ncp,
    design = "Multiple regression: joint test of predictors",
    sample_size_unit = "total",
    assurance_ceiling = ncp_solution$ceiling,
    inputs = list(F_observed = F_observed, N = N, p = p, p_joint = p_joint,
                  random_predictors = as.numeric(random_predictors),
                  alpha_prior = alpha_prior_input, alpha_planned = alpha_planned,
                  assurance = assurance, desired_power = desired_power)
  )
}
