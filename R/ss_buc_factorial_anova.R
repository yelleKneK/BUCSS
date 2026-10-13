#' Necessary Sample Size to Reach Desired Power for a Two-Way Between-Subjects
#' ANOVA Using a Publication Bias and Uncertainty Correction Procedure
#'
#' @description \code{ss_buc_factorial_anova} returns the necessary per-cell
#'   sample size to achieve a desired level of statistical power for a planned
#'   study testing an omnibus main effect or interaction in a two-way fully
#'   between-subjects ANOVA, based on information obtained from a previous study.
#'   The effect from the previous study can be corrected for publication bias
#'   and/or uncertainty to provide a sample size that will achieve more accurate
#'   statistical power for a planned study, when compared to approaches that use
#'   a sample effect size at face value or rely on sample size only. The bias
#'   and uncertainty adjusted previous study noncentrality parameter is also
#'   returned, together with the planning effect size, the effect size the plan
#'   is built on, as Cohen's \emph{f}.
#'
#'   For a single between-subjects factor, use
#'   \code{\link{ss_buc_one_way_anova}}; for a design with more than two factors,
#'   or for a single-degree-of-freedom contrast, use
#'   \code{\link{ss_buc_factorial_anova_general}}.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss_buc_factorial_anova} uses the observed
#'   \emph{F} value and sample size from a previous study to correct the
#'   noncentrality parameter associated with the effect of interest for
#'   publication bias and/or uncertainty. This new estimated noncentrality
#'   parameter is then used to calculate the necessary per-cell sample size to
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
#'   \code{ss_buc_factorial_anova} assumes that the planned study will have
#'   equal \emph{n}. Unequal \emph{n} in the previous study is handled in the
#'   following way for between-subjects ANOVA designs. If the user enters an
#'   \emph{N} not equally divisible by the number of cells, the function
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
#'   is the larger of the two rounded designs' requirements. For a prior study
#'   with unequal group sizes, or a total sample size that does not divide
#'   evenly, the largest supportable assurance is computed at the balanced
#'   design the planner plans from, so it can be below
#'   \eqn{1 - p/\alpha}{1 - p/alpha} at the study's own \emph{p} value.
#'   Although equal-\emph{n} previous studies are preferable, this approach
#'   will work well as long as the cell sizes are only slightly discrepant.
#'
#' @param F_observed Observed \emph{F} value from a previous study used to plan
#'   sample size for a planned study.
#' @param N Total sample size of the previous study.
#' @param levels_A Number of levels for Factor A.
#' @param levels_B Number of levels for Factor B.
#' @param effect Effect of most interest to the planned study: main effect of A
#'   (\code{factor_A}), main effect of B (\code{factor_B}), or interaction
#'   (\code{interaction}). Case-insensitive shorthand is also accepted but not
#'   recommended in scripts: \code{"A"} (or \code{"a"}) for \code{factor_A},
#'   \code{"B"} for \code{factor_B}, and \code{"AxB"} (or \code{"AB"}) for
#'   \code{interaction}.
#' @template planning-params
#'
#' @templateVar size_phrase per-cell sample size
#' @templateVar ncp_scale on the scale of the \emph{F} statistic, at the prior study's sample size (for a prior \emph{N} that does not divide evenly, at the balanced sample size the planner substitutes for it; see Details)
#' @templateVar effect_phrase \code{f_planning}, Cohen's \emph{f} for the effect tested, relative to the error term of its test, \eqn{\sqrt{\lambda/N}}{sqrt(lambda/N)}, where \eqn{\lambda}{lambda} is \code{ncp_adjusted} and \emph{N} the total of the balanced prior design it belongs to (see Details), so \eqn{\lambda = N f^2}{lambda = N f^2}, the convention of power analysis for between-subjects designs, and \eqn{f^2/(1 + f^2)}{f^2/(1 + f^2)} is the population partial eta squared for the effect
#' @template return-effect
#'
#' @export
#'
#' @examples
#' result <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
#'   levels_B = 3, effect = "factor_B", alpha_prior = .05, alpha_planned = .05,
#'   assurance = .80, desired_power = .80)
#' result
#'
#' # Pull the planned per-cell sample size, the implied total sample size, and
#' # the adjusted noncentrality parameter out of the result:
#' result$value[result$term == "necessary_n_per_cell"]   # per cell
#' result$value[result$term == "total_N"]                 # total across all cells
#' result$value[result$term == "ncp_adjusted"]
#'
#' # Requesting more assurance than the prior result can support stops with an
#' # informative error naming the largest workable assurance (here near .83):
#' try(ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
#'   levels_B = 3, effect = "factor_B", assurance = .95))
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
ss_buc_factorial_anova <- function(F_observed, N, levels_A, levels_B,
                                   effect = c("factor_A", "factor_B", "interaction"),
                                   alpha_prior = .05, alpha_planned = .05,
                                   assurance = .80, desired_power = .80) {
  effect <- .match_effect(effect, c("factor_A", "factor_B", "interaction"))

  v <- .validate_planning_inputs(alpha_prior, alpha_planned, assurance, desired_power)
  alpha_prior <- v$alpha_prior
  alpha_prior_input <- v$alpha_prior_input
  assurance <- v$assurance
  desired_power <- v$desired_power

  if (missing(N)) stop("You must specify 'N', which is the total sample size.", call. = FALSE)

  if (missing(levels_A)) stop("You must specify 'levels_A', the number of levels of the first factor.", call. = FALSE)
  if (missing(levels_B) || is.null(levels_B)) stop("You must specify 'levels_B' for a two-way design. For a single between-subjects factor, use 'ss_buc_one_way_anova'.", call. = FALSE)
  .check_scalar_finite(F_observed, "F_observed")
  .check_count(N, "N", min = 2)
  .check_count(levels_A, "levels_A", min = 2)
  .check_count(levels_B, "levels_B", min = 2)

  cells <- levels_A * levels_B
  if (N < 2 * cells) stop("Your prior study 'N' is too small for this design: at least two observations per cell are required, so 'N' must be at least 2 * 'levels_A' * 'levels_B'.", call. = FALSE)
  df_numerator <- switch(effect,
                         factor_A = levels_A - 1,
                         factor_B = levels_B - 1,
                         interaction = (levels_A - 1) * (levels_B - 1))

  # The prior study's own p value, at its own degrees of freedom, so that a
  # statistic significant there but not at the rounded-down design is not
  # called nonsignificant (see .stop_nonsignificant()).
  p_own <- pf(F_observed, df1 = df_numerator, df2 = N - cells, lower.tail = FALSE)

  ## Rounding up
  n_ru <- ceiling(N / cells)
  N_ru <- n_ru * cells
  df_denominator_ru <- N_ru - cells

  crit_F_ru <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator_ru)
  if (F_observed <= crit_F_ru)
    .stop_nonsignificant("F statistic", "'F_observed'", p_own, alpha_prior,
                         pf(F_observed, df_numerator, floor(N / cells) * cells - cells,
                            lower.tail = FALSE), assurance)

  solution_ru <- .solve_ncp_assurance(
    .tm_f(F_observed, crit_F_ru, df_numerator, df_denominator_ru), assurance)
  ncp_ru <- solution_ru$ncp

  ## Rounding down
  n_rd <- floor(N / cells)
  N_rd <- n_rd * cells
  df_denominator_rd <- N_rd - cells

  crit_F_rd <- qf(1 - alpha_prior, df1 = df_numerator, df2 = df_denominator_rd)
  if (F_observed <= crit_F_rd)
    .stop_nonsignificant("F statistic", "'F_observed'", p_own, alpha_prior,
                         pf(F_observed, df_numerator, df_denominator_rd, lower.tail = FALSE),
                         assurance)

  solution_rd <- .solve_ncp_assurance(
    .tm_f(F_observed, crit_F_rd, df_numerator, df_denominator_rd), assurance)
  ncp_rd <- solution_rd$ncp

  # Both branches are solved before either can refuse, so the refusal names
  # the smaller of the two ceilings, which is a feasible assurance for both.
  if (ncp_ru == 0 || ncp_rd == 0)
    .stop_zero_ncp(min(solution_ru$ceiling, solution_rd$ceiling))

  power_at <- function(n_rep, n_prior, ncp_b) {
    denom_df <- n_rep * cells - cells
    critical_F <- qf(1 - alpha_planned, df1 = df_numerator, df2 = denom_df)
    1 - pf(critical_F, df1 = df_numerator, df2 = denom_df, ncp = (n_rep / n_prior) * ncp_b)
  }
  repn_ru <- .smallest_n_for_power(function(k) power_at(k, n_ru, ncp_ru), desired_power)
  repn_rd <- .smallest_n_for_power(function(k) power_at(k, n_rd, ncp_rd), desired_power)

  output_n <- max(repn_ru, repn_rd)
  df_error <- output_n * cells - cells
  # The reported noncentrality parameter is the smaller of the two branches',
  # and n_b the prior per-group size it belongs to. The power of the pair this
  # function reports: the returned sample size evaluated at the returned
  # (adjusted) noncentrality parameter, using the prior n that parameter came
  # from. It is at least 'desired_power', because the returned size is the
  # larger of the two roundings' requirements, and usually a little higher,
  # because the size is a whole number. The planning effect size restates the
  # same parameter as Cohen's f at n_b, so the ordinary power of the returned
  # size at that effect is exactly actual_power.
  use_rd <- ncp_rd <= ncp_ru
  n_b <- if (use_rd) n_rd else n_ru
  ncp_b <- if (use_rd) ncp_rd else ncp_ru
  actual_power <- power_at(output_n, n_b, ncp_b)
  .bucss_plan_result(
    sample_size = output_n,
    size_term = "necessary_n_per_cell",
    effect_size_term = "f_planning",
    effect_size = .es_f(ncp_b, n_b * cells),
    actual_power = actual_power,
    df_effect = df_numerator,
    df_error = df_error,
    ncp = ncp_b,
    design = "Two-way between-subjects ANOVA",
    sample_size_unit = "per cell",
    assurance_ceiling = min(solution_ru$ceiling, solution_rd$ceiling),
    total_n = output_n * cells,
    effect = effect,
    inputs = list(F_observed = F_observed, N = N, levels_A = levels_A,
                  levels_B = levels_B, alpha_prior = alpha_prior_input,
                  alpha_planned = alpha_planned, assurance = assurance,
                  desired_power = desired_power)
  )
}
