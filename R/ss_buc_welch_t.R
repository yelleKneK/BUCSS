#' Necessary Sample Size to Reach Desired Power for a Welch (Unequal Variance) t
#' Test Using a Publication Bias and Uncertainty Correction Procedure
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
#'   returned, together with the planning effect size, the effect size the plan
#'   is built on, as a standardized mean difference.
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
#'   this test, and computes the planned study's power the same way.
#'   \strong{A minimum planned sample size of 10 per group is imposed, because
#'   of estimation accuracy issues at smaller sample sizes.} With fewer per
#'   group the approximation is not accurate enough: the real test estimates
#'   its two variances, and so its degrees of freedom, from too few
#'   observations for the noncentral \emph{t} to describe it, and the
#'   approximation overstates its power, badly so for the smallest (a plan of
#'   2 per group that it puts at .82 has real power near .50). From 10 per
#'   group, at an \code{alpha_planned} of .05, the approximation is accurate:
#'   across ratios from 1 to 5, the power of the real Welch test is within
#'   .006 of what this function reports for every plan it returns, and within
#'   .003 from about 12 per group at a target power of .80, and its Type I
#'   error stays between .047 and .054 against a nominal .05 across ratios
#'   from 1 to 3. At a stricter \code{alpha_planned} the approximation
#'   overstates the power of the smallest plans by more, at 10 per group by
#'   up to about .012 at .01, .016 at .005, and .025 at .001, and it is
#'   within .006 again from 14, 16, and 20 per group, respectively. When the
#'   target power is reached with fewer than 10 per group, the planner returns
#'   10, whose power exceeds the target.
#'
#'   The correction of the prior study's statistic uses the same approximation
#'   and holds up well unless a prior group is very small. Simulating prior
#'   studies as raw data, published by the real Welch test, the corrected
#'   parameter fell at or below the truth in about .79 to .81 of them at an
#'   \code{assurance} of .80 when each prior group had 10 to 20 participants
#'   (a refusal counted as at or below). With 5 to 8 in the smaller prior
#'   group it ranged from about .78 to .86 when the true ratio was supplied,
#'   and fell as low as about .73 when the ratio was taken from the prior
#'   study's own standard deviations.
#'   One consequence of reading the prior study on the Welch-Satterthwaite
#'   degrees of freedom the supplied ratio implies, rather than on those its
#'   data produced: when the prior groups are small and unequal in size or in
#'   spread, including at the default ratio of 1, a prior result its paper
#'   reported as significant can fall short of the critical value here and be
#'   refused as nonsignificant. Such a result is typically close to the
#'   publication threshold, where it could support little assurance in any
#'   case.
#'
#'   \strong{The assumption that matters is the ratio itself.} One
#'   \code{sd_ratio} describes both studies: the prior study, whose statistic
#'   is read under it, and the planned study. Getting it wrong costs far more
#'   than the distributional approximation does. For the example below, the
#'   plan built assuming equal variances (213 per group) delivers about .76
#'   power rather than .80 if the ratio is in fact 2 in both studies, and
#'   about .43 if the prior study's groups were equally variable but the
#'   planned study's ratio is 2, with the mean difference held fixed in units
#'   of the first group's standard deviation. Re-running across a range of
#'   plausible ratios, and planning for the least favorable one you find
#'   credible, guards against the first kind of error. The second, a planned
#'   study more unequal than the prior one, is not something the prior result
#'   can speak to, and needs its own allowance. The required sample size does
#'   not always rise as the two standard deviations move apart: when the group
#'   assumed to be more variable is the smaller one it often falls, though
#'   not always (for the prior study of the example below it rises).
#'
#'   When the ratio is 1 and the prior groups were the same size, this planner
#'   returns what \code{\link{ss_buc_independent_t}} returns, except for the
#'   minimum of 10 per group imposed here. With
#'   unequal prior group sizes the two differ, and should: the equal-variance
#'   planner reduces the prior study to an equivalent equal-\emph{n} design at
#'   the harmonic mean, while this one uses the actual Welch-Satterthwaite
#'   degrees of freedom, which are smaller when the groups are unbalanced. The
#'   same prior \emph{t} is then less compelling, and the plan is usually
#'   larger when one prior group is more than about twice the size of the
#'   other; when the two sizes are close it is more often the same, and it
#'   can be smaller, because the equal-variance planner rounds the harmonic
#'   mean down and up to whole group sizes and returns the larger of the two
#'   plans.
#'
#'   The approach uses a likelihood function of a truncated noncentral
#'   \emph{F} distribution, where the truncation occurs due to small effect
#'   sizes being unobserved due to publication bias. The numerator of the
#'   likelihood function is the density of a noncentral \emph{F} distribution.
#'   The denominator is the power of the test, which serves to truncate the
#'   distribution. In the two-group case, this formula reduces to the density of
#'   a truncated noncentral \emph{t} distribution. (See Taylor & Muller, 1996,
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
#'   \code{alpha_prior} is the significance level a prior result had to reach
#'   to be published: the correction assumes that only results with \emph{p}
#'   values below \code{alpha_prior} reach the literature, and so reach you.
#'   Different \eqn{\alpha} levels represent differing amounts of publication
#'   bias. For example, setting \code{alpha_prior = .20} assumes that studies
#'   with \emph{p} values less than .20 are published, whereas those with
#'   larger \emph{p} values are not, which reflects less severe publication
#'   bias than the default of .05. If the previous study of interest has not
#'   been subjected to publication bias (e.g., a pilot study), set
#'   \code{alpha_prior = 1}. Because no \emph{p} value exceeds 1, this means
#'   no publication filter at all: any result would have been published
#'   whatever its \emph{p} value, so there is no publication bias to correct,
#'   and only the correction for uncertainty, set by \code{assurance},
#'   remains.
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
#'   need about 10 percent fewer participants in total at a ratio of 2, and
#'   about 20 percent fewer at a ratio of 3.
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
#' @templateVar ncp_scale on the scale of the \emph{t} statistic (the noncentrality parameter of the prior study's noncentral \emph{t})
#' @templateVar effect_phrase \code{smd_planning}, a standardized mean difference (a magnitude) whose standardizer is the root mean square of the two groups' standard deviations (Bonett, 2008), \eqn{\lambda\sqrt{(1/n_1 + r^2/n_2)/((1 + r^2)/2)}}{lambda * sqrt((1/n_1 + r^2/n_2)/((1 + r^2)/2))}, where \eqn{\lambda}{lambda} is \code{ncp_adjusted}, \eqn{n_1}{n_1} and \eqn{n_2}{n_2} are the prior group sizes, and \emph{r} is the ratio of standard deviations; with \emph{n} per group in the planned study its noncentrality parameter is this value times \eqn{\sqrt{n/2}}{sqrt(n/2)}, as for the equal-variance test, and with a ratio of 1 it is Cohen's \eqn{\delta}{delta}
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
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
#' @references Bonett, D. G. (2008). Confidence intervals for standardized
#'   linear contrasts of means. \emph{Psychological Methods, 13,} 99--109.
#'   \doi{10.1037/1082-989X.13.2.99}
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
  # A spread too extreme to compute with is reported in terms of what the
  # user supplied: a pair whose ratio overflows or underflows (it used to be
  # reported as a bad 'sd_ratio', an argument the user had not given), and,
  # below, a ratio whose fourth power overflows in the degrees of freedom (it
  # used to end in base R's "missing value where TRUE/FALSE needed").
  spread_given <- if (gave_sd) "'sd_1' and 'sd_2'" else if (gave_var) "'var_1' and 'var_2'" else "'sd_ratio'"
  stop_spread <- function()
    stop("The two groups' spread differs by too many orders of magnitude to plan from: the Welch-Satterthwaite degrees of freedom cannot be computed. Please check ", spread_given, ".", call. = FALSE)
  if (gave_sd) {
    if (missing(sd_1) || missing(sd_2)) stop("You must specify both 'sd_1' and 'sd_2', the group standard deviations of the previous study.", call. = FALSE)
    .check_scalar_finite(sd_1, "sd_1")
    .check_scalar_finite(sd_2, "sd_2")
    if (sd_1 <= 0 || sd_2 <= 0) stop("'sd_1' and 'sd_2' must be positive numbers.", call. = FALSE)
    sd_ratio <- sd_2 / sd_1
    if (!is.finite(sd_ratio) || sd_ratio <= 0) stop_spread()
  } else if (gave_var) {
    if (missing(var_1) || missing(var_2)) stop("You must specify both 'var_1' and 'var_2', the group variances of the previous study.", call. = FALSE)
    .check_scalar_finite(var_1, "var_1")
    .check_scalar_finite(var_2, "var_2")
    if (var_1 <= 0 || var_2 <= 0) stop("'var_1' and 'var_2' must be positive numbers.", call. = FALSE)
    sd_ratio <- sqrt(var_2 / var_1)
    if (!is.finite(sd_ratio) || sd_ratio <= 0) stop_spread()
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
  # The squares in the degrees of freedom overflow for a ratio near 1e77 or
  # more (where exactly depends on the group sizes), in the prior study or at
  # the smallest planned study; larger planned studies are further from it. A
  # very small ratio underflows harmlessly to the limiting degrees of freedom.
  if (!is.finite(df_prior) ||
      !is.finite(welch_df(.welch_min_n, .welch_min_n, sd_ratio))) stop_spread()

  value_critical <- qt(1 - alpha_prior / 2, df = df_prior)

  if (t_stat <= value_critical) .stop_nonsignificant("t statistic", "'t_observed'", alpha_prior = alpha_prior)

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
  # The search starts at 10 per group, a minimum imposed because of estimation
  # accuracy issues at smaller sample sizes. Below it the noncentral t on the
  # population Welch-Satterthwaite df overstates the power of the real Welch
  # test, whose degrees of freedom are estimated from the planned study's data:
  # by up to about .33 at 2 per group, and by more than .006 below 10 per group
  # (exact power of the real test, checked against 400,000-replicate
  # simulation, 2026-09-28). From 10 per group the two agree within .006 at an
  # alpha_planned of .05, so a plan that reaches its target sooner is raised
  # to 10. The figure is for .05 only: at a stricter alpha_planned the
  # overstatement at 10 per group reaches about .012 (.01), .016 (.005), and
  # .025 (.001), and is within .006 from 14, 16, and 20 per group (exact
  # power of the real test over ratios 1 to 5, recomputed 2026-10-12). The
  # minimum itself does not depend on alpha_planned.
  output_n <- .smallest_n_for_power(power_at, desired_power, start = .welch_min_n)

  actual_power <- power_at(output_n)
  .bucss_plan_result(
    sample_size = output_n,
    size_term = "necessary_n_per_group",
    effect_size_term = "smd_planning",
    effect_size = .es_smd_welch(ncp, n_1, n_2, sd_ratio),
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

# The smallest per-group size ss_buc_welch_t() returns (see the comment at its
# planned-n search).
.welch_min_n <- 10
