#' Necessary sample size to reach desired power for an analysis of covariance
#' using a publication bias and uncertainty correction procedure
#'
#' @description \code{ss_buc_ancova} returns the necessary per-cell sample size
#'   to achieve a desired level of statistical power for a planned study testing
#'   a group effect (omnibus or contrast) while adjusting for one or more
#'   covariates, based on information obtained from a previous study. The effect
#'   from the previous study can be corrected for publication bias and/or
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
#'   An analysis of covariance is a between-subjects analysis with additional
#'   estimated parameters. \code{ss_buc_ancova} takes the number of covariates
#'   directly and carries them through to the planned study: the prior study's
#'   error degrees of freedom are \emph{N} minus the number of cells minus the
#'   number of covariates, and the planned study's error degrees of freedom are
#'   reduced by the same number of covariates. Supplying
#'   \code{n_covariates = 0} reproduces the corresponding
#'   \code{\link{ss_buc_one_way_anova}} or
#'   \code{\link{ss_buc_factorial_anova_general}} result exactly.
#'
#'   By default the omnibus group effect is planned
#'   (\code{df_numerator = cells - 1}). For a single-degree-of-freedom contrast
#'   among the adjusted means, supply the contrast's \emph{F} (the square of its
#'   \emph{t}) and set \code{df_numerator = 1}.
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
#'   The returned \code{actual_power} is the statistical power of the plan
#'   this function reports: the power of a study of the returned sample size
#'   when the effect is the returned bias and uncertainty adjusted
#'   noncentrality parameter. Because a sample size must be a whole number,
#'   \code{actual_power} is never below the \code{desired_power} that was
#'   requested and is usually a little above it.
#'
#'   The covariates are treated as fixed, which is the usual convention in
#'   sample size planning for an analysis of covariance. When the covariates
#'   are themselves sampled, the planned study's power is very slightly lower
#'   than reported.
#'
#'   \code{ss_buc_ancova} assumes that the planned study will have equal
#'   \emph{n}. If the prior \emph{N} is not equally divisible by the number of
#'   cells, the per-cell size is rounded both up and down, the planned sample
#'   size is computed under each, and the larger of the two is returned, with
#'   the smaller of the two adjusted noncentrality parameters, to be
#'   conservative.
#'
#' @param F_observed Observed \emph{F} value for the effect of interest from a
#'   previous study used to plan sample size for a planned study. For a
#'   single-degree-of-freedom contrast, this is the square of the contrast's
#'   \emph{t}.
#' @param N Total sample size of the previous study.
#' @param cells Number of cells (groups) in the between-subjects design.
#' @param n_covariates Number of covariates adjusted for in the analysis.
#' @param df_numerator Numerator degrees of freedom for the effect of interest.
#'   Defaults to \code{cells - 1}, the omnibus group effect; use 1 for a
#'   single-degree-of-freedom contrast.
#' @template planning-params
#'
#' @templateVar size_phrase per-cell sample size
#' @template return
#'
#' @export
#'
#' @examples
#' # Three groups, two covariates, prior N = 120: the omnibus group effect.
#' result <- ss_buc_ancova(F_observed = 6, N = 120, cells = 3,
#'   n_covariates = 2, alpha_prior = .05, alpha_planned = .05,
#'   assurance = .80, desired_power = .80)
#' result
#'
#' # A single-degree-of-freedom contrast among the adjusted means.
#' ss_buc_ancova(F_observed = 9, N = 120, cells = 3, n_covariates = 2,
#'   df_numerator = 1)
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .93):
#' try(ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2,
#'   assurance = .95))
#'
#' @seealso \code{\link{ss_buc_factorial_anova_general}}, which this function
#'   calls, and which takes the error degrees of freedom directly.
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu})
#'
#' @template references
ss_buc_ancova <- function(F_observed, N, cells, n_covariates,
                          df_numerator = cells - 1, alpha_prior = .05,
                          alpha_planned = .05, assurance = .80,
                          desired_power = .80) {
  if (missing(N)) stop("You must specify 'N', which is the total sample size.", call. = FALSE)
  if (missing(cells)) stop("You must specify 'cells', the number of cells (groups) in the between-subjects design.", call. = FALSE)
  if (missing(n_covariates)) stop("You must specify 'n_covariates', the number of covariates adjusted for in the analysis.", call. = FALSE)
  .check_count(cells, "cells", min = 2)
  .check_count(n_covariates, "n_covariates", min = 0)
  .check_count(N, "N", min = 2)
  if (N - cells - n_covariates < 1) {
    stop("The combination of 'N', 'cells', and 'n_covariates' leaves no error degrees of freedom in the prior study. Check the design.", call. = FALSE)
  }

  # The covariates are the nuisance parameters the general between-subjects
  # planner already carries through to the planned study.
  res <- ss_buc_factorial_anova_general(
    F_observed = F_observed, N = N, cells = cells,
    df_numerator = df_numerator,
    df_denominator = N - cells - n_covariates,
    alpha_prior = alpha_prior, alpha_planned = alpha_planned,
    assurance = assurance, desired_power = desired_power)

  .relabel_bucss_result(
    res,
    design = "Analysis of covariance",
    sample_size_unit = "per cell",
    size_term = "necessary_n_per_cell",
    inputs = list(F_observed = F_observed, N = N, cells = cells,
                  n_covariates = n_covariates, df_numerator = df_numerator,
                  alpha_prior = res$value[res$term == "alpha_prior"],
                  alpha_planned = res$value[res$term == "alpha_planned"],
                  assurance = res$value[res$term == "assurance"],
                  desired_power = res$value[res$term == "desired_power"])
  )
}
