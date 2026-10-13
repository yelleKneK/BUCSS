#' Necessary Sample Size to Reach Desired Power for a Nested Model Chi Square
#' Difference Test Using a Publication Bias and Uncertainty Correction Procedure
#'
#' @description \code{ss_buc_chisq_diff} returns the necessary total sample size
#'   to achieve a desired level of statistical power for a planned study whose
#'   effect of interest is tested by the chi square difference between two
#'   nested models (the likelihood ratio test used, for example, to test a
#'   constrained path in a structural equation model), based on information
#'   obtained from a previous study. The effect from the previous study can be
#'   corrected for publication bias and/or uncertainty to provide a sample size
#'   that will achieve more accurate statistical power for a planned study, when
#'   compared to approaches that use a sample effect size at face value or rely
#'   on sample size only. The bias and uncertainty adjusted previous study
#'   noncentrality parameter is also returned, together with the planning
#'   effect size, the effect size the plan is built on, as the root mean square
#'   error of approximation of the difference test,
#'   \eqn{\mathrm{RMSEA}_D}{RMSEA_D} (D for difference).
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   \strong{Scope, and one thing this function is not for.} The correction
#'   applies to a test whose statistic is large when the effect is present, so
#'   that publication selects large values. A chi square difference test
#'   between nested models is such a test: under the alternative it is
#'   approximately noncentral chi square with degrees of freedom equal to the
#'   number of constraints released and a noncentrality parameter proportional
#'   to sample size, and a paper reports the constrained path as supported when
#'   the difference test is significant. Steiger, Shapiro, and Browne (1985)
#'   give the asymptotic distribution of sequential chi square statistics for
#'   nested models and of the difference tests between them, and Satorra and
#'   Saris (1985) compute the power of the likelihood ratio test in covariance
#'   structure analysis from the noncentral chi square distribution, as is
#'   done here. The correction is therefore built for the difference test.
#'
#'   It is \emph{not} appropriate for an omnibus model fit chi square, where a
#'   publishable result is a \emph{small} statistic (a model that fits). The
#'   selection region is inverted there, so the truncated likelihood this
#'   function builds does not describe it.
#'
#'   \strong{The function refuses that case rather than only warning about it.}
#'   This is why both models' degrees of freedom are required rather than only
#'   their difference. An omnibus fit chi square is the comparison of a model
#'   against the \emph{saturated} model, and a saturated model has zero degrees
#'   of freedom, so \code{df_full = 0} is its signature. Given only the
#'   difference, nothing distinguishes a fit test from a genuine comparison of
#'   two substantive models, and the mistake would pass silently. Supplying
#'   both is also a more faithful description of what was tested, and the two
#'   values are what a paper reports for each model. A just-identified full
#'   model (every parameter free, as in a path model compared with the same
#'   model with one path fixed at zero) also has zero degrees of freedom, and
#'   there the difference test is a genuine one. The plan depends only on the
#'   number of constraints, so for that case the call with \code{df_full = 1}
#'   and \code{df_restricted} one more than the number of constraints gives
#'   the plan.
#'
#'   The correction uses a likelihood function of a truncated noncentral
#'   chi square distribution, the same construction the rest of the package uses
#'   for the noncentral \emph{F} (see Taylor & Muller, 1996, Equation 2.1;
#'   Anderson, Kelley, & Maxwell, 2017, Equation 3; and Anderson & Maxwell,
#'   2017): the truncated area between the critical value and the observed
#'   statistic, divided by the power of the test. Because the noncentrality
#'   parameter of a likelihood ratio test is proportional to sample size, the
#'   planned study's noncentrality parameter is the corrected one scaled by the
#'   ratio of the planned to the prior sample size.
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
#'   \code{alpha_prior = 1}, which may also be written
#'   \code{alpha_prior = "No Publication Filter"}; Anderson, Kelley, and
#'   Maxwell (2017, p. 1552) note that it can be set equal to 1 when
#'   publication bias is thought not to apply. Because no \emph{p} value
#'   exceeds 1, this means no publication filter at all: any result would
#'   have been published whatever its \emph{p} value, so there is no
#'   publication bias to correct, and only the correction for uncertainty,
#'   set by \code{assurance}, remains.
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
#'   The chi square difference test is asymptotic, so both the correction and
#'   the planned-study power inherit that approximation. The likelihood ratio
#'   test is liberal in small samples, and more so the more parameters a model
#'   estimates. For a test that drops 2 of 6 parameters of a normal linear
#'   model, its actual Type I error is about .07 at a total sample size of 50,
#'   .055 at 200, and .052 at 400; dropping 5 of 21 parameters, it is about
#'   .08 at 200. At such sample sizes the real test's power can also sit above
#'   or below \code{actual_power}. Treat a planned sample size that is small
#'   relative to the number of parameters estimated with caution, and read a
#'   very small one as a signal that the prior effect was large rather than as
#'   a serious recommendation. The search for the planned sample size starts
#'   at the number of constraints plus 2 (\code{df_restricted - df_full + 2}),
#'   so no plan is smaller than that, even when a smaller study would already
#'   reach the desired power.
#'
#' @param chisq_observed Observed chi square difference between the two nested
#'   models in the previous study.
#' @param N Total sample size of the previous study.
#' @param df_full Degrees of freedom of the full (less constrained) model. A
#'   saturated full model has zero, the signature of an omnibus test of model
#'   fit rather than a difference test between two substantive models; the
#'   function refuses that case (see Details, which also covers a
#'   just-identified full model).
#' @param df_restricted Degrees of freedom of the restricted (more
#'   constrained) model. It exceeds \code{df_full} by the number of
#'   constraints imposed, which is the difference test's degrees of freedom.
#' @template planning-params
#'
#' @templateVar size_phrase total sample size
#' @templateVar ncp_scale on the scale of the chi square statistic, at the prior study's sample size
#' @templateVar effect_phrase \code{rmsea_d_planning}, \eqn{\mathrm{RMSEA}_D}{RMSEA_D}, where the D stands for difference: the root mean square error of approximation of the chi square difference test between the two nested models (Savalei, Brace, & Fouladi, 2024), as distinct from the RMSEA of either model's fit; here \eqn{\sqrt{\lambda/(\mathrm{df}_D N)}}{sqrt(lambda/(df_D N))}, where \eqn{\lambda}{lambda} is \code{ncp_adjusted}, \eqn{\mathrm{df}_D}{df_D} is the number of constraints (\code{df_restricted - df_full}, the difference test's degrees of freedom), and \emph{N} the prior sample size, the scaling this planner plans with; the usual sample formula divides by \eqn{N - 1}{N - 1} instead, which gives a value larger by the factor \eqn{\sqrt{N/(N - 1)}}{sqrt(N/(N - 1))}
#' @template return
#'
#' @export
#'
#' @examples
#' # A path constrained to zero: the full model has 42 degrees of freedom and
#' # the restricted model 43, so the difference test has 1.
#' result <- ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42,
#'   df_restricted = 43, alpha_prior = .05, alpha_planned = .05,
#'   assurance = .80, desired_power = .80)
#' result
#'
#' # A saturated full model means the comparison is an omnibus test of model
#' # fit, which this method does not describe, and the function says so:
#' try(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 0,
#'   df_restricted = 43))
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .95):
#' try(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42,
#'   df_restricted = 43, assurance = .99))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
#' @references Satorra, A., & Saris, W. E. (1985). Power of the likelihood
#'   ratio test in covariance structure analysis. \emph{Psychometrika, 50,}
#'   83--90. \doi{10.1007/BF02294150}
#'
#'   Savalei, V., Brace, J. C., & Fouladi, R. T. (2024). We need to
#'   change how we compute RMSEA for nested model comparisons in structural
#'   equation modeling. \emph{Psychological Methods, 29,} 480--493.
#'   \doi{10.1037/met0000537}
#'
#'   Steiger, J. H., Shapiro, A., & Browne, M. W. (1985). On the multivariate
#'   asymptotic distribution of sequential chi-square statistics.
#'   \emph{Psychometrika, 50,} 253--263. \doi{10.1007/BF02294104}
ss_buc_chisq_diff <- function(chisq_observed, N, df_full, df_restricted,
                              alpha_prior = .05, alpha_planned = .05,
                              assurance = .80, desired_power = .80) {
  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, desired_power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  desired_power <- v$desired_power

  if (missing(N)) stop("You must specify 'N', which is the total sample size of the previous study.", call. = FALSE)
  if (missing(df_full)) stop("You must specify 'df_full', the degrees of freedom of the full (less constrained) model.", call. = FALSE)
  if (missing(df_restricted)) stop("You must specify 'df_restricted', the degrees of freedom of the restricted (more constrained) model.", call. = FALSE)
  .check_scalar_finite(chisq_observed, "chisq_observed")
  .check_count(N, "N", min = 2)
  .check_count(df_full, "df_full", min = 0)
  .check_count(df_restricted, "df_restricted", min = 1)

  # The guard against the one misuse this design invites. A saturated model has
  # zero degrees of freedom, so 'df_full = 0' means the comparison is the
  # restricted model's own test of fit against the saturated model, not a
  # difference test between two substantive models. The correction does not
  # describe that test: a publishable fit chi square is a small one, so the
  # selection region is inverted and the truncated likelihood this method
  # builds is the wrong one. Requiring both models' degrees of freedom, rather
  # than only their difference, is what makes this detectable at all.
  if (df_full == 0) stop("'df_full' is 0, so the full model is saturated and this comparison is the restricted model's omnibus test of model fit, not a difference test between two substantive models. This function does not apply to an omnibus fit chi square: the correction assumes a literature publishes significant results, whereas a publishable fit chi square is a small (nonsignificant) one, which inverts the selection region the method is built on. A just-identified full model also has 0 degrees of freedom; if that is your case and the comparison is a genuine difference test, the plan depends only on the number of constraints, so call with 'df_full = 1' and 'df_restricted' equal to 1 plus the number of constraints.", call. = FALSE)

  if (df_restricted <= df_full) stop("'df_restricted' must be greater than 'df_full': the restricted model imposes constraints, so it has more degrees of freedom than the full model. If the two are reversed, swap them.", call. = FALSE)

  df_difference <- df_restricted - df_full

  value_critical <- qchisq(1 - alpha_prior, df = df_difference)

  if (chisq_observed <= value_critical) .stop_nonsignificant("chi square difference", "the prior result", alpha_prior = alpha_prior)

  ncp_solution <- .solve_ncp_assurance(
    .tm_chisq(chisq_observed, value_critical, df_difference), assurance)
  ncp <- ncp_solution$ncp

  if (ncp == 0) .stop_zero_ncp(ncp_solution$ceiling)

  # The noncentrality parameter of a likelihood ratio test is proportional to
  # sample size, so it scales by the ratio of planned to prior N. The degrees
  # of freedom do not depend on the sample size.
  critical_chisq <- qchisq(1 - alpha_planned, df = df_difference)
  power_at <- function(n_rep) {
    1 - pchisq(critical_chisq, df = df_difference, ncp = (n_rep / N) * ncp)
  }
  # A planned study cannot have fewer observations than the constraints the
  # difference test releases, so the search starts there rather than at 2. A
  # returned size near this floor means the prior effect was very large, not
  # that such a study would be sensible; the Details section says so.
  output_n <- .smallest_n_for_power(power_at, desired_power,
                                    start = df_difference + 2)

  actual_power <- power_at(output_n)
  .bucss_plan_result(
    sample_size = output_n,
    size_term = "necessary_N",
    effect_size_term = "rmsea_d_planning",
    effect_size = .es_rmsea_d(ncp, N, df_difference),
    actual_power = actual_power,
    df_effect = df_difference,
    df_error = NULL,
    ncp = ncp,
    design = "Nested model chi square difference test",
    sample_size_unit = "total",
    assurance_ceiling = ncp_solution$ceiling,
    inputs = list(chisq_observed = chisq_observed, N = N,
                  df_full = df_full, df_restricted = df_restricted,
                  alpha_prior = alpha_prior_input,
                  alpha_planned = alpha_planned, assurance = assurance,
                  desired_power = desired_power)
  )
}
