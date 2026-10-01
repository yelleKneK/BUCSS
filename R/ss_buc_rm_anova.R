#' Necessary Sample Size to Reach Desired Power for a One- or Two-Way
#' Within-Subjects ANOVA Using a Publication Bias and Uncertainty Correction
#' Procedure
#'
#' @description \code{ss_buc_rm_anova} returns the necessary total sample size to
#'   achieve a desired level of statistical power for a planned study testing an
#'   omnibus effect using a one- or two-way fully within-subjects ANOVA, based on
#'   information obtained from a previous study. The effect from the previous
#'   study can be corrected for publication bias and/or uncertainty to provide a
#'   sample size that will achieve more accurate statistical power for a planned
#'   study, when compared to approaches that use a sample effect size at face
#'   value or rely on sample size only. The bias and uncertainty adjusted
#'   previous study noncentrality parameter is also returned, together with the
#'   planning effect size, the effect size the plan is built on, as Cohen's
#'   \emph{f}.
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
#'   \strong{Sphericity.} The default \code{epsilon = 1} assumes sphericity, as
#'   the uncorrected within-subjects test does. When a prior study reports
#'   Greenhouse-Geisser (1959) or Huynh-Feldt (1976) corrected degrees of
#'   freedom, supply \code{epsilon}, the reported numerator degrees of freedom
#'   divided by the uncorrected ones. For example, \emph{F}(1.73, 50.2) for a
#'   three-level within-subjects factor gives \code{epsilon = 1.73/2}, about
#'   .865, and the denominator agrees: 50.2/58 is also about .865. The ratio of
#'   the two reported degrees of freedom, 50.2/1.73 = 29.0, is one less than
#'   the number of subjects the analysis used, here the study's 30; when it
#'   points to a different number (after dropping subjects with missing data,
#'   say), that number is the \emph{N} to supply. An implied \code{epsilon}
#'   above 1 or below its lower bound (see Arguments) means the reported
#'   degrees of freedom are not a Greenhouse-Geisser or Huynh-Feldt correction
#'   of this effect; a Huynh-Feldt estimate above 1 is reported as the
#'   uncorrected degrees of freedom, which is \code{epsilon = 1}.
#'   \code{ss_buc_rm_anova} then applies the approximation of Muller and Barton
#'   (1989) to both studies: the prior \emph{F} is treated as noncentral
#'   \emph{F} with both degrees of freedom and the noncentrality parameter
#'   multiplied by \code{epsilon}, the publication threshold is the critical
#'   value at the corrected degrees of freedom, and the planned study's power
#'   is computed the same way. Each within-subjects effect has its own
#'   \eqn{\epsilon}, so supply the one reported for the effect you plan from.
#'   The returned \code{ncp_adjusted} stays on the sphericity scale, so it is
#'   comparable across values of \code{epsilon}, and the planned test's
#'   corrected degrees of freedom are the \code{df_effect} and \code{df_error}
#'   attributes. One \code{epsilon} describes both studies, because it is a
#'   property of the covariance matrix of the repeated measures, which the
#'   planned study shares with the prior one; the reported value is the prior
#'   study's estimate of it, and re-running across the values you find
#'   plausible shows how much the plan depends on it. The value used is echoed
#'   as the \code{epsilon} row of the result, 1 when sphericity is assumed.
#'
#'   \strong{The correction under \code{epsilon} is an approximation.} Below
#'   sphericity the distribution of a within-subjects \emph{F} depends on the
#'   whole covariance matrix of the repeated measures, and on how the effect
#'   lines up with it, not on \eqn{\epsilon} alone. In simulations of raw data
#'   (one-way designs of 4 and 6 levels with \eqn{\epsilon} from .45 to .86,
#'   prior studies of 12 to 40 subjects published when their Greenhouse-Geisser
#'   test reached \emph{p} < .05, and a 2 x 4 split-plot interaction), planning
#'   from each published study with its own reported Greenhouse-Geisser
#'   estimate at an \code{assurance} of .80 gave a plan reaching .80 power for
#'   the Greenhouse-Geisser test in .70 to .80 of replications when the effect
#'   followed the dominant pattern of variation (a linear trend over
#'   autocorrelated occasions), and in .92 to .99 when it lay in the least
#'   variable contrasts. Ignoring the departure (\code{epsilon = 1}) reached it
#'   in .26 to .64 of replications for the trend-like effects and in .70 to .87
#'   for the effects in the least variable contrasts. The shortfall does not
#'   shrink as the prior study grows, and is not larger near the lower bound of
#'   \eqn{\epsilon}. The planned study's power is approximated the same way.
#'   Separate simulations of planned one-way studies of 10 to 30 subjects, with
#'   4, 6, or 8 levels, measured the real power of the Greenhouse-Geisser test,
#'   whose estimate of \eqn{\epsilon} is biased downward. For a trend-like
#'   effect (a linear trend over occasions with first-order autoregressive
#'   correlations of .3 to .9, \eqn{\epsilon} from .50 to .94) it fell below
#'   \code{actual_power} by about .04 to .08 with four levels and .06 to .13
#'   with six, the larger shortfalls for the smaller plans; a planned study
#'   analyzed with the Huynh-Feldt correction tracks \code{actual_power} more
#'   closely. Even when sphericity holds and \code{epsilon = 1}, a planned
#'   study analyzed with the Greenhouse-Geisser correction anyway, as some
#'   software does by default, has less power than \code{actual_power}, the
#'   power of the uncorrected test: about .02 to .06 less with four levels, .03
#'   to .11 with six, and .04 to .15 with eight, for plans of 30 down to 10
#'   subjects (the Huynh-Feldt test is within about .02). When a paper reports
#'   both corrections, supply the Greenhouse-Geisser one: the Huynh-Feldt
#'   estimate is larger, and planning from it attained about .02 to .09 less,
#'   the larger differences from the smaller prior studies. When sphericity in
#'   fact holds but the prior study reported a Greenhouse-Geisser correction
#'   anyway, supplying it is conservative, with a median plan 7 to 34 percent
#'   larger.
#'
#'   \strong{The planning effect size and other power software.} Here
#'   \code{f_planning} is Cohen's \emph{f} relative to the effect's own error
#'   term, whose variance already reflects the correlation among the repeated
#'   measures, so \eqn{f^2/(1 + f^2)}{f^2/(1 + f^2)} is the population value of
#'   the partial eta squared a paper reports for the effect. It is not the
#'   \emph{f} that power functions for repeated measures designs usually take,
#'   such as \code{ss_power_rm_anova()} in the sibling package \code{DMAR}:
#'   those standardize the effect by the standard deviation of the raw scores
#'   and take the correlation among the measures as a separate argument, so
#'   giving them \code{f_planning}, with any correlation, plans a different
#'   study. For a one-way design with \emph{a} levels, the value that
#'   reproduces this plan in \code{ss_power_rm_anova()}, with the correlation
#'   set to 0 and the same \code{epsilon}, is \code{f_planning} times
#'   \eqn{\sqrt{(a - 1)/a}}{sqrt((a - 1)/a)}.
#'
#' @param F_observed Observed \emph{F} value from a previous study used to plan
#'   sample size for a planned study.
#' @param N Total sample size of the previous study.
#' @param levels_A Number of levels for Factor A.
#' @param levels_B Number of levels for Factor B, which is \code{NULL} for a
#'   single-factor design.
#' @param effect Effect of most interest to the planned study: main effect of A
#'   (\code{factor_A}), main effect of B (\code{factor_B}), or interaction
#'   (\code{interaction}). Case-insensitive shorthand is also accepted but not
#'   recommended in scripts: \code{"A"} (or \code{"a"}) for \code{factor_A},
#'   \code{"B"} for \code{factor_B}, and \code{"AxB"} (or \code{"AB"}) for
#'   \code{interaction}.
#' @template planning-params
#' @param epsilon The sphericity parameter \eqn{\epsilon} of the effect tested,
#'   from 1 divided by the effect's numerator degrees of freedom (its lower
#'   bound) to 1 (the default, which assumes sphericity). For a prior study
#'   reported with Greenhouse-Geisser or Huynh-Feldt corrected degrees of
#'   freedom, it is the reported numerator degrees of freedom divided by the
#'   uncorrected ones. An effect with one numerator degree of freedom (a
#'   two-level factor) is always spherical and takes only 1. The same value
#'   describes the planned study, and it is echoed as the \code{epsilon} row of
#'   the result. See Details.
#'
#' @templateVar size_phrase total sample size
#' @templateVar ncp_scale on the scale of the \emph{F} statistic, at the prior study's sample size, and on the sphericity scale for every \code{epsilon} (see Details)
#' @templateVar effect_phrase \code{f_planning}, Cohen's \emph{f} for the effect tested, relative to its own error term, \eqn{\sqrt{\lambda/(kN)}}{sqrt(lambda/(k N))}, where \eqn{\lambda}{lambda} is \code{ncp_adjusted}, \emph{N} the prior number of subjects, and \emph{k} the effect's numerator degrees of freedom (the error degrees of freedom per subject), so \eqn{f^2/(1 + f^2)}{f^2/(1 + f^2)} is the population partial eta squared for the effect; like \code{ncp_adjusted} it is on the sphericity scale, and for \code{epsilon} below 1 the power statement that follows holds under the same Muller and Barton (1989) approximation, with the corrected degrees of freedom and noncentrality parameter; it is not the \emph{f} that power functions for repeated measures take (see Details)
#' @template return-effect
#'
#' @export
#'
#' @examples
#' result <- ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3,
#'   effect = "factor_B", alpha_prior = .05, alpha_planned = .05,
#'   assurance = .80, desired_power = .80)
#' result
#'
#' # A prior study reported with Greenhouse-Geisser corrected degrees of
#' # freedom, 1.73 and 50.2, for a three-level factor with 30 subjects: the
#' # correction is 1.73/2, and the plan is larger than under sphericity (92
#' # subjects).
#' ss_buc_rm_anova(F_observed = 6.5, N = 30, levels_A = 3, epsilon = 1.73 / 2)
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
#' @references Greenhouse, S. W., & Geisser, S. (1959). On methods in the
#'   analysis of profile data. \emph{Psychometrika, 24,} 95--112.
#'   \doi{10.1007/BF02289823}
#'
#'   Huynh, H., & Feldt, L. S. (1976). Estimation of the Box correction for
#'   degrees of freedom from sample data in randomized block and split-plot
#'   designs. \emph{Journal of Educational Statistics, 1,} 69--82.
#'   \doi{10.3102/10769986001001069}
#'
#'   Muller, K. E., & Barton, C. N. (1989). Approximate power for
#'   repeated-measures ANOVA lacking sphericity. \emph{Journal of the American
#'   Statistical Association, 84,} 549--555.
#'   \doi{10.1080/01621459.1989.10478802}
ss_buc_rm_anova <- function(F_observed, N, levels_A, levels_B = NULL,
                        effect = c("factor_A", "factor_B", "interaction"),
                        alpha_prior = .05, alpha_planned = .05,
                        assurance = .80, desired_power = .80, epsilon = 1) {
  effect <- .match_effect(effect, c("factor_A", "factor_B", "interaction"))

  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, desired_power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  desired_power <- v$desired_power

  if (missing(N)) stop("You must specify 'N', which is the total sample size.", call. = FALSE)

  if (is.null(levels_B) && effect != "factor_A") stop("For a one-way within-subjects design ('levels_B' not supplied), the only testable effect is the single factor; keep the default 'effect = \"factor_A\"'. To test 'factor_B' or the interaction, supply 'levels_B'.", call. = FALSE)
  if (missing(levels_A)) stop("You must specify 'levels_A', the number of levels of the within-subjects factor.", call. = FALSE)
  .check_scalar_finite(F_observed, "F_observed")
  .check_count(N, "N", min = 2)
  .check_count(levels_A, "levels_A", min = 2)
  if (!is.null(levels_B)) .check_count(levels_B, "levels_B", min = 2)

  n <- N

  if (is.null(levels_B)) {
    df_numerator <- levels_A - 1
  } else {
    df_numerator <- switch(effect,
                           factor_A = levels_A - 1,
                           factor_B = levels_B - 1,
                           interaction = (levels_A - 1) * (levels_B - 1))
  }
  # Every effect of a fully within-subjects design is a within-subjects
  # effect, whose contrasts span its own numerator df.
  epsilon <- .check_epsilon(epsilon, df_numerator)
  df_denominator <- df_numerator * (n - 1)

  # Under a sphericity correction (Muller & Barton, 1989) both degrees of
  # freedom and the noncentrality parameter are multiplied by epsilon, for
  # the prior study and the planned one alike; see .tm_f(). At epsilon = 1
  # every product is exact.
  crit_F <- qf(1 - alpha_prior, df1 = epsilon * df_numerator,
               df2 = epsilon * df_denominator)
  if (F_observed <= crit_F) .stop_nonsignificant("F statistic", "'F_observed'")

  ncp_solution <- .solve_ncp_assurance(
    .tm_f(F_observed, crit_F, df_numerator, df_denominator, epsilon), assurance)
  ncp <- ncp_solution$ncp

  if (ncp == 0) .stop_zero_ncp(ncp_solution$ceiling)

  power_at <- function(n_rep) {
    df1 <- epsilon * df_numerator
    denom_df <- epsilon * (df_numerator * (n_rep - 1))
    critical_F <- qf(1 - alpha_planned, df1 = df1, df2 = denom_df)
    1 - pf(critical_F, df1 = df1, df2 = denom_df,
           ncp = epsilon * ((n_rep / n) * ncp))
  }
  output_n <- .smallest_n_for_power(power_at, desired_power)

  df_error <- epsilon * (df_numerator * (output_n - 1))
  actual_power <- power_at(output_n)
  .bucss_plan_result(
    sample_size = output_n,
    size_term = "necessary_N",
    effect_size_term = "f_planning",
    effect_size = .es_f(ncp, N, k = df_numerator),
    actual_power = actual_power,
    df_effect = epsilon * df_numerator,
    df_error = df_error,
    ncp = ncp,
    design = "One- or two-way within-subjects ANOVA",
    sample_size_unit = "total",
    assurance_ceiling = ncp_solution$ceiling,
    effect = effect,
    inputs = list(F_observed = F_observed, N = N, levels_A = levels_A,
                  levels_B = levels_B, epsilon = epsilon,
                  alpha_prior = alpha_prior_input,
                  alpha_planned = alpha_planned, assurance = assurance,
                  desired_power = desired_power)
  )
}
