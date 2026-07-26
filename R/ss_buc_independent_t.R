#' Necessary sample size to reach desired power for an independent t test (or
#' standardized mean difference) using a publication bias and uncertainty
#' correction procedure
#'
#' @description \code{ss_buc_independent_t} returns the necessary per-group sample size
#'   to achieve a desired level of statistical power for a planned study using
#'   an independent \emph{t} test, based on information obtained from a previous
#'   study. The effect from the previous study can be corrected for publication
#'   bias and/or uncertainty to provide a sample size that will achieve more
#'   accurate statistical power for a planned study, when compared to approaches
#'   that use a sample effect size at face value or rely on sample size only.
#'   The bias and uncertainty adjusted previous study noncentrality parameter is
#'   also returned, which can be transformed to various effect size metrics.
#'
#'   \code{ss_buc_smd} is an alias for \code{ss_buc_independent_t}: the two are the
#'   same function. The two-group \emph{t} test and the standardized mean
#'   difference (Cohen's \eqn{d}) describe the same comparison, so the same
#'   planner serves both framings. Prefer the test-specific name
#'   (\code{ss_buc_independent_t}) when the prior \emph{t} was computed from raw
#'   (unstandardized) data, since the planner works directly from the observed
#'   \emph{t}.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss_buc_independent_t} uses the observed
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
#'   If you are working from a standardized mean difference (Cohen's \eqn{d})
#'   rather than a \emph{t} statistic, convert it before calling: for two
#'   independent groups \eqn{t = d / \sqrt{1/n_1 + 1/n_2}}, which equals
#'   \eqn{d\sqrt{n/2}} when the per-group sample sizes are equal.
#'
#'   \code{ss_buc_independent_t} assumes that the planned study will have equal
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
#' @template planning-params
#'
#' @templateVar size_phrase per-group sample size
#' @template return
#'
#' @export
#'
#' @examples
#' result <- ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = .05,
#'   alpha_planned = .05, assurance = .80, power = .80)
#' result
#' result$value[result$term == "necessary_sample_size"]
#'
#' # ss_buc_smd is the same function under an effect-size name
#' ss_buc_smd(t_observed = 3, n = 20)
#'
#' # Asking for more assurance than the prior result can support stops with an
#' # informative error that reports the largest workable assurance (here near
#' # .90, the value 1 - p/alpha_prior for this prior result):
#' try(ss_buc_independent_t(t_observed = 3, n = 20, assurance = .95))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
ss_buc_independent_t <- function(t_observed, n, N, alpha_prior = .05, alpha_planned = .05,
                           assurance = .80, power = .80) {
  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  power <- v$power

  .check_scalar_finite(t_observed, "t_observed")

  if (!missing(N)) {
    if (!is.null(N)) {
      # The round-down branch has df = 2 * floor(N / 2) - 2, so N = 3 would
      # already drive the error df to 0; the smallest workable total is 4.
      .check_count(N, "N", min = 4)
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
    for (n_i in n) .check_count(n_i, "n", min = 2)
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

  ## Rounding up
  value_critical_ru <- qt(1 - alpha_prior / 2, df = DF_ru)
  if (t_observed <= value_critical_ru) stop("Your observed t statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 't_observed' exceeds the critical value.")

  solution_ru <- .solve_ncp_assurance(
    .tm_t(t_observed, value_critical_ru, DF_ru), assurance)
  ncp_ru <- solution_ru$ncp

  if (ncp_ru == 0) .stop_zero_ncp(solution_ru$ceiling)

  power_at <- function(n_rep, n_prior, ncp_b) {
    denom_df <- (2 * n_rep) - 2
    critical_t <- qt(1 - alpha_planned / 2, df = denom_df)
    scaled <- sqrt(n_rep / n_prior) * ncp_b
    (1 - pt(critical_t, df = denom_df, ncp = scaled)) +
      pt(-1 * critical_t, df = denom_df, ncp = scaled)
  }
  repn_ru <- .smallest_n_for_power(function(k) power_at(k, n_ru, ncp_ru), power)

  ## Rounding down
  value_critical_rd <- qt(1 - alpha_prior / 2, df = DF_rd)
  if (t_observed <= value_critical_rd) stop("Your observed t statistic is nonsignificant based on your specified 'alpha_prior' of the prior study. Please increase 'alpha_prior' so 't_observed' exceeds the critical value.")

  solution_rd <- .solve_ncp_assurance(
    .tm_t(t_observed, value_critical_rd, DF_rd), assurance)
  ncp_rd <- solution_rd$ncp

  if (ncp_rd == 0) .stop_zero_ncp(solution_rd$ceiling)

  repn_rd <- .smallest_n_for_power(function(k) power_at(k, n_rd, ncp_rd), power)

  output_n <- max(repn_ru, repn_rd)
  df_error <- 2 * output_n - 2
  .bucss_power_result(
    sample_size = output_n,
    df_effect = 1,
    df_error = df_error,
    ncp = min(ncp_ru, ncp_rd),
    design = "Independent t test",
    sample_size_unit = "per group",
    assurance_ceiling = min(solution_ru$ceiling, solution_rd$ceiling),
    total_n = 2 * output_n,
    inputs = list(t_observed = t_observed, alpha_prior = alpha_prior_input,
                  alpha_planned = alpha_planned, assurance = assurance,
                  power = power)
  )
}

#' @rdname ss_buc_independent_t
#' @export
ss_buc_smd <- ss_buc_independent_t
