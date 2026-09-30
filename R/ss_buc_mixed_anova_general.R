#' Necessary Sample Size to Reach Desired Power for a Split-Plot (Mixed) ANOVA
#' With Any Number of Factors Using a Publication Bias and Uncertainty
#' Correction Procedure
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
#'   parameter is also returned, together with the planning effect size, the
#'   effect size the plan is built on, as Cohen's \emph{f}.
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
#'   \code{ss_buc_mixed_anova_general} assumes that the planned study will have
#'   equal \emph{n}. Unequal \emph{n} in the previous study is handled in the
#'   following way for split-plot designs. If the user enters an \emph{N} not
#'   equally divisible by the number of between-subjects cells, the function
#'   calculates \emph{n} by dividing \emph{N} by the number of cells and both
#'   rounding up and rounding down the result, effectively assuming equal
#'   \emph{n}. The suggested sample size for the planned study is calculated
#'   using both of these values of \emph{n}, and the function returns the
#'   larger of these two suggestions, to be conservative. The adjusted
#'   noncentrality parameter that is output is the lower of the two
#'   possibilities, again, to be conservative. The planning effect size belongs
#'   to the same balanced design as the reported noncentrality parameter, so an
#'   ordinary power analysis at it can return a sample size a little smaller
#'   than the one reported here (never a larger one), because the reported size
#'   is the larger of the two rounded designs' requirements. Although
#'   equal-\emph{n} previous studies are preferable, this approach will work
#'   well as long as the cell sizes are only slightly discrepant.
#'
#'   \strong{Sphericity.} The default \code{epsilon = 1} assumes sphericity, as
#'   the uncorrected within-subjects test does. When a prior study reports
#'   Greenhouse-Geisser (1959) or Huynh-Feldt (1976) corrected degrees of
#'   freedom, supply \code{epsilon}, the reported numerator degrees of freedom
#'   divided by the uncorrected ones. For example, \emph{F}(1.5, 130.5) for a
#'   \code{within_only} effect with 2 numerator degrees of freedom in a study
#'   of 90 subjects in 3 groups gives \code{epsilon = 1.5/2}, and the
#'   denominator agrees: 130.5/174 is also .75. The ratio of the two reported
#'   degrees of freedom, 130.5/1.5 = 87, is the number of subjects the analysis
#'   used less the number of groups (for \code{between_within}, that difference
#'   divided by the between-subjects degrees of freedom,
#'   \code{df_numerator}/\code{df_num_within}); when it points to a different
#'   number of subjects than the study reports, that number is the \emph{N} to
#'   supply. An implied \code{epsilon} above 1 or below its lower bound (see
#'   Arguments) means the reported degrees of freedom are not a
#'   Greenhouse-Geisser or Huynh-Feldt correction of this effect; a Huynh-Feldt
#'   estimate above 1 is reported as the uncorrected degrees of freedom, which
#'   is \code{epsilon = 1}. \code{ss_buc_mixed_anova_general} then applies the
#'   approximation of Muller and Barton (1989) to both studies: the prior
#'   \emph{F} is treated as noncentral \emph{F} with both degrees of freedom
#'   and the noncentrality parameter multiplied by \code{epsilon}, the
#'   publication threshold is the critical value at the corrected degrees of
#'   freedom, and the planned study's power is computed the same way. It
#'   applies to the effects with a within-subjects component
#'   (\code{within_only} and \code{between_within}); a \code{between_only}
#'   effect is tested on subject means and does not involve sphericity. For
#'   \code{between_within} it is the \eqn{\epsilon} of the within-subjects
#'   component's contrasts. The returned \code{ncp_adjusted} stays on the
#'   sphericity scale, so it is comparable across values of \code{epsilon}, and
#'   the planned test's corrected degrees of freedom are the \code{df_effect}
#'   and \code{df_error} attributes. One \code{epsilon} describes both studies,
#'   because it is a property of the covariance matrix of the repeated
#'   measures, which the planned study shares with the prior one; the reported
#'   value is the prior study's estimate of it, and re-running across the
#'   values you find plausible shows how much the plan depends on it. The value
#'   used is echoed as the \code{epsilon} row of the result, 1 when sphericity
#'   is assumed.
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
#'   in .26 to .64. The shortfall does not shrink as the prior study grows, and
#'   is not larger near the lower bound of \eqn{\epsilon}. The planned study's
#'   power is approximated the same way, and for the Greenhouse-Geisser test,
#'   whose estimate of \eqn{\epsilon} is biased downward, the real power of a
#'   plan of 10 to 30 subjects can be up to about .08 below \code{actual_power}
#'   for a trend-like effect; a planned study analyzed with the Huynh-Feldt
#'   correction tracks \code{actual_power} more closely. Even when sphericity
#'   holds and \code{epsilon = 1}, a planned study analyzed with the
#'   Greenhouse-Geisser correction anyway, as some software does by default,
#'   has up to about .08 less power than \code{actual_power}, the power of the
#'   uncorrected test, for plans of 10 to 30 subjects (the Huynh-Feldt test is
#'   within about .02). When a paper reports both corrections, supply the
#'   Greenhouse-Geisser one: the Huynh-Feldt estimate is larger, and planning
#'   from it attained .03 to .08 less. When sphericity in fact holds but the
#'   prior study reported a Greenhouse-Geisser correction anyway, supplying it
#'   is conservative, with a median plan 7 to 34 percent larger.
#'
#'   \strong{The planning effect size and other power software.} Here
#'   \code{f_planning} is Cohen's \emph{f} relative to the effect's own error
#'   term, whose variance already reflects the correlation among the repeated
#'   measures, so \eqn{f^2/(1 + f^2)}{f^2/(1 + f^2)} is the population value of
#'   the partial eta squared a paper reports for the effect. It is not the
#'   \emph{f} that power functions for repeated measures designs usually take,
#'   such as \code{ss_power_split_plot_anova()} in the sibling package
#'   \code{DMAR}: those standardize the effect by the standard deviation of the
#'   raw scores and take the correlation among the measures as a separate
#'   argument, so giving them \code{f_planning}, with any correlation, plans a
#'   different study. For a design with one between-subjects factor and one
#'   within-subjects factor of \emph{b} levels, the value that reproduces this
#'   plan in \code{ss_power_split_plot_anova()}, with the correlation set to 0
#'   and the same \code{epsilon}, is \code{f_planning} times
#'   \eqn{\sqrt{k/b}}{sqrt(k/b)}, where \emph{k} is the error degrees of
#'   freedom per subject (1 for \code{between_only}, \emph{b} - 1 for the other
#'   two effects). The \code{between_only} effect is included: its
#'   \code{f_planning} is relative to the variance of subject means, not of raw
#'   scores.
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
#' @param epsilon The sphericity parameter \eqn{\epsilon} of the
#'   within-subjects contrasts, from 1/\code{df_numerator} (for
#'   \code{within_only}) or 1/\code{df_num_within} (for \code{between_within}),
#'   its lower bound, to 1 (the default, which assumes sphericity); when that
#'   bound is 1 it takes only 1. For a prior study reported with
#'   Greenhouse-Geisser or Huynh-Feldt corrected degrees of freedom, it is the
#'   reported numerator degrees of freedom divided by the uncorrected ones. The
#'   same value describes the planned study, and it is echoed as the
#'   \code{epsilon} row of the result. It does not apply to
#'   \code{between_only}, whose test is on subject means: there a value other
#'   than 1 is ignored with a warning, and the result has no \code{epsilon}
#'   row. See Details.
#'
#' @templateVar size_phrase per-group sample size
#' @templateVar ncp_scale on the scale of the \emph{F} statistic, at the prior study's sample size (for a prior \emph{N} that does not divide evenly, at the balanced sample size the planner substitutes for it), and on the sphericity scale for every \code{epsilon} (see Details)
#' @templateVar effect_phrase \code{f_planning}, Cohen's \emph{f} for the effect tested, relative to its own error term, \eqn{\sqrt{\lambda/(kN)}}{sqrt(lambda/(k N))}, where \eqn{\lambda}{lambda} is \code{ncp_adjusted}, \emph{N} the total of the balanced prior design it belongs to, and \emph{k} the error degrees of freedom per subject, 1 for \code{between_only}, \code{df_numerator} for \code{within_only}, and \code{df_num_within} for \code{between_within}, so \eqn{f^2/(1 + f^2)}{f^2/(1 + f^2)} is the population partial eta squared for the effect; like \code{ncp_adjusted} it is on the sphericity scale, and for \code{epsilon} below 1 the power statement that follows holds under the same Muller and Barton (1989) approximation, with the corrected degrees of freedom and noncentrality parameter; it is not the \emph{f} that power functions for repeated measures take (see Details)
#' @template return-effect
#'
#' @export
#'
#' @examples
#' result <- ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
#'   num_groups = 3, effect = "between_only",
#'   alpha_prior = .05, alpha_planned = .05, assurance = .80, desired_power = .80)
#' result
#'
#' # A within-subjects effect with 2 numerator degrees of freedom, from a
#' # prior study of 90 subjects in 3 groups reported with Greenhouse-Geisser
#' # corrected degrees of freedom of 1.5 and 130.5: the correction is 1.5/2
#' # (under sphericity, the plan is 41 per group).
#' ss_buc_mixed_anova_general(F_observed = 8, N = 90, df_numerator = 2,
#'   num_groups = 3, effect = "within_only", epsilon = .75)
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .82):
#' try(ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
#'   num_groups = 3, effect = "between_only",
#'   assurance = .95))
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
ss_buc_mixed_anova_general <- function(F_observed, N, df_numerator, num_groups,
                                 effect = c("between_only", "within_only", "between_within"),
                                 df_num_within, alpha_prior = .05,
                                 alpha_planned = .05, assurance = .80,
                                 desired_power = .80, epsilon = 1) {
  effect <- match.arg(effect)

  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, desired_power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  desired_power <- v$desired_power

  if (missing(N)) stop("You must specify 'N', which is the total sample size.", call. = FALSE)
  if (missing(df_numerator)) stop("You must specify 'df_numerator', the numerator degrees of freedom for the effect of interest.", call. = FALSE)
  if (missing(num_groups)) stop("You must specify 'num_groups', the number of between-subjects groups.", call. = FALSE)
  if (effect == "between_within" && missing(df_num_within)) stop("For 'effect = \"between_within\"' you must specify 'df_num_within', the numerator degrees of freedom of the within-subjects component.", call. = FALSE)
  .check_scalar_finite(F_observed, "F_observed")
  .check_count(N, "N", min = 2)
  .check_count(num_groups, "num_groups", min = 2)
  .check_count(df_numerator, "df_numerator", min = 1)
  if (effect == "between_within") .check_count(df_num_within, "df_num_within", min = 1)
  if (N < 2 * num_groups) stop("Your prior study 'N' is too small for this design: at least two observations per between-subjects group are required, so 'N' must be at least 2 * 'num_groups'.", call. = FALSE)
  if (effect == "between_only" && df_numerator > num_groups - 1) stop("For 'effect = \"between_only\"', 'df_numerator' cannot exceed 'num_groups' - 1, the degrees of freedom available among the between-subjects groups.", call. = FALSE)
  if (effect == "between_within" &&
      (df_numerator %% df_num_within != 0 || df_numerator / df_num_within > num_groups - 1))
    stop("For 'effect = \"between_within\"', 'df_numerator' must be the between-subjects degrees of freedom (at most 'num_groups' - 1) times 'df_num_within', so it must be a multiple of 'df_num_within' no larger than ('num_groups' - 1) * 'df_num_within'.", call. = FALSE)
  if (effect != "between_within" && !missing(df_num_within))
    warning("'df_num_within' is used only when 'effect = \"between_within\"', so it was ignored.", call. = FALSE)
  # Sphericity concerns the within-subjects contrasts only. For within_only
  # they span the effect's own df_numerator dimensions; for between_within
  # they are the within-subjects component's df_num_within. A between_only
  # effect is tested on subject means and has no epsilon.
  epsilon <- switch(effect,
    between_only = .ignore_epsilon(epsilon, !missing(epsilon),
                                   "'effect = \"between_only\"'"),
    within_only = .check_epsilon(epsilon, df_numerator),
    between_within = .check_epsilon(epsilon, df_num_within))

  denom_mult <- switch(effect,
                       between_only = 1,
                       within_only = df_numerator,
                       between_within = df_num_within)

  ## Rounding down
  n_rd <- floor(N / num_groups)
  N_rd <- n_rd * num_groups
  df_denominator_rd <- (N_rd - num_groups) * denom_mult

  # The prior study's own p value, at its own degrees of freedom, so that a
  # statistic significant there but not at the rounded-down design is not
  # called nonsignificant (see .stop_nonsignificant()).
  # Under a sphericity correction (Muller & Barton, 1989) both degrees of
  # freedom and the noncentrality parameter are multiplied by epsilon, for
  # the prior study and the planned one alike; see .tm_f(). At epsilon = 1
  # (always, for between_only) every product is exact.
  p_own <- pf(F_observed, df1 = epsilon * df_numerator,
              df2 = epsilon * ((N - num_groups) * denom_mult), lower.tail = FALSE)

  crit_F_rd <- qf(1 - alpha_prior, df1 = epsilon * df_numerator,
                  df2 = epsilon * df_denominator_rd)
  if (F_observed <= crit_F_rd)
    .stop_nonsignificant("F statistic", "'F_observed'", p_own, alpha_prior,
                         pf(F_observed, epsilon * df_numerator,
                            epsilon * df_denominator_rd, lower.tail = FALSE),
                         assurance)

  solution_rd <- .solve_ncp_assurance(
    .tm_f(F_observed, crit_F_rd, df_numerator, df_denominator_rd, epsilon),
    assurance)
  ncp_rd <- solution_rd$ncp

  ## Rounding up
  n_ru <- ceiling(N / num_groups)
  N_ru <- n_ru * num_groups
  df_denominator_ru <- (N_ru - num_groups) * denom_mult

  crit_F_ru <- qf(1 - alpha_prior, df1 = epsilon * df_numerator,
                  df2 = epsilon * df_denominator_ru)
  if (F_observed <= crit_F_ru)
    .stop_nonsignificant("F statistic", "'F_observed'", p_own, alpha_prior,
                         pf(F_observed, epsilon * df_numerator,
                            epsilon * df_denominator_rd, lower.tail = FALSE),
                         assurance)

  solution_ru <- .solve_ncp_assurance(
    .tm_f(F_observed, crit_F_ru, df_numerator, df_denominator_ru, epsilon),
    assurance)
  ncp_ru <- solution_ru$ncp

  # Both branches are solved before either can refuse, so the refusal names
  # the smaller of the two ceilings, which is a feasible assurance for both.
  if (ncp_ru == 0 || ncp_rd == 0)
    .stop_zero_ncp(min(solution_ru$ceiling, solution_rd$ceiling))

  power_at <- function(n_rep, n_prior, ncp_b) {
    df1 <- epsilon * df_numerator
    denom_df <- epsilon * (((n_rep * num_groups) - num_groups) * denom_mult)
    critical_F <- qf(1 - alpha_planned, df1 = df1, df2 = denom_df)
    1 - pf(critical_F, df1 = df1, df2 = denom_df,
           ncp = epsilon * ((n_rep / n_prior) * ncp_b))
  }
  repn_rd <- .smallest_n_for_power(function(k) power_at(k, n_rd, ncp_rd), desired_power)
  repn_ru <- .smallest_n_for_power(function(k) power_at(k, n_ru, ncp_ru), desired_power)

  inputs <- list(F_observed = F_observed, N = N, df_numerator = df_numerator,
                 df_num_within = if (effect == "between_within") df_num_within,
                 num_groups = num_groups,
                 # echoed for the effects it applies to (1 under sphericity),
                 # not for between_only, which it does not
                 epsilon = if (effect != "between_only") epsilon,
                 alpha_prior = alpha_prior_input,
                 alpha_planned = alpha_planned, assurance = assurance,
                 desired_power = desired_power)

  output_n <- max(repn_rd, repn_ru)
  df_error <- epsilon * (((output_n * num_groups) - num_groups) * denom_mult)
  # The reported noncentrality parameter is the smaller of the two branches',
  # and n_b the prior per-group size it belongs to. The power of the pair this
  # function reports: the returned sample size evaluated at the returned
  # (adjusted) noncentrality parameter, using the prior n that parameter came
  # from. It is at least 'desired_power', because the returned size is the
  # larger of the two roundings' requirements, and usually a little higher,
  # because the size is a whole number. The planning effect size restates the
  # same parameter as Cohen's f at n_b, relative to the effect's own error
  # term, whose degrees of freedom grow by denom_mult per participant (1 for
  # between_only). The ordinary power of the returned size at that effect
  # (under the same Muller-Barton approximation when epsilon is below 1, the
  # effect staying on the sphericity scale) is exactly actual_power.
  use_rd <- ncp_rd <= ncp_ru
  n_b <- if (use_rd) n_rd else n_ru
  ncp_b <- if (use_rd) ncp_rd else ncp_ru
  actual_power <- power_at(output_n, n_b, ncp_b)
  .bucss_plan_result(
    sample_size = output_n,
    size_term = "necessary_n_per_group",
    effect_size_term = "f_planning",
    effect_size = .es_f(ncp_b, n_b * num_groups, k = denom_mult),
    actual_power = actual_power,
    df_effect = epsilon * df_numerator,
    df_error = df_error,
    ncp = ncp_b,
    design = "Split-plot (mixed) ANOVA (any number of factors)",
    sample_size_unit = "per group",
    assurance_ceiling = min(solution_ru$ceiling, solution_rd$ceiling),
    total_n = output_n * num_groups,
    effect = effect,
    inputs = inputs
  )
}
