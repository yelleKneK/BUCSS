# Deprecated 1.x API, retained for backward compatibility. Each dot-named
# function from BUCSS 1.x is kept here as a thin wrapper that translates its
# arguments to the 2.x convention and calls the corresponding ss_buc_* planner,
# returning that planner's bucss_power object. The wrappers are not actively
# developed; the statistical work lives in the ss_buc_* functions. A one-time
# (per session, per function) deprecation warning steers users to the new name;
# the ss_buc_* functions themselves are silent. The 1.x grid-resolution 'step'
# argument is accepted for call compatibility and ignored (the 2.0.0 engine
# finds the noncentrality parameter by root finding, not a grid). The 1.x effect
# values (e.g. "factor.A", "between.only") are translated to their snake_case
# 2.x forms with a single dot-to-underscore rule.

# Track which deprecated functions have already warned this session, so the
# warning fires once per function per session rather than on every call.
.bucss_deprecated_warned <- new.env(parent = emptyenv())

.deprecate_once <- function(old, replacement) {
  if (is.null(.bucss_deprecated_warned[[old]])) {
    assign(old, TRUE, envir = .bucss_deprecated_warned)
    msg <- paste0(
      "'", old, "' is deprecated as of BUCSS 2.0.0 and is no longer actively ",
      "developed. It now calls ", replacement, " and returns a 'bucss_power' ",
      "object (its first two elements remain the necessary sample size and the ",
      "adjusted noncentrality parameter, so 'result[[1]]' and 'result[[2]]' ",
      "still work). Please use the new function directly. This warning is shown ",
      "once per session."
    )
    # A classed condition (mirroring R's own .Deprecated(), which signals a
    # 'deprecatedWarning') so callers can catch or suppress it by class; call =
    # NULL keeps the uninformative wrapper call out of the message.
    warning(warningCondition(
      msg, class = c("bucss_deprecatedWarning", "deprecatedWarning"),
      call = NULL))
  }
  invisible(NULL)
}

# Translate a 1.x dotted 'effect' value to its 2.x snake_case form. The single
# dot-to-underscore rule covers every case: "factor.A" -> "factor_A",
# "between.only" -> "between_only", and the dotless values ("between",
# "interaction") are unchanged.
.translate_effect <- function(effect, choices_1x) {
  gsub(".", "_", match.arg(effect, choices_1x), fixed = TRUE)
}

#' Deprecated 1.x functions in BUCSS
#'
#' @description The dot-named functions from BUCSS 1.x are \strong{deprecated} as
#'   of BUCSS 2.0.0 and are no longer actively developed. They are retained so
#'   that scripts written for BUCSS 1.x continue to run unchanged. Each one
#'   translates its arguments to the 2.x convention and calls the corresponding
#'   \code{ss_buc_*} function, returning that function's
#'   \code{\link{bucss_power}} object.
#'
#'   For backward compatibility the returned object still supports the 1.x
#'   positional extraction: \code{result[[1]]} is the necessary sample size and
#'   \code{result[[2]]} is the bias and uncertainty adjusted noncentrality
#'   parameter, the same two values the 1.x functions returned as an unnamed
#'   list (this is provided by a \code{[[} method for \code{bucss_power}). The
#'   richer named view
#'   (\code{result$value}, \code{tidy(result)}, \code{glance(result)}) is also
#'   available.
#'
#'   Calling one of these functions issues a deprecation warning that names its
#'   replacement, once per function per session; the \code{ss_buc_*} functions
#'   themselves never warn. New code should call the \code{ss_buc_*} functions
#'   directly.
#'
#' @details The replacements are:
#'
#' \tabular{ll}{
#'   \strong{Deprecated} \tab \strong{Replacement} \cr
#'   \code{ss.power.it} \tab \code{\link{ss_buc_independent_t}} (alias \code{\link{ss_buc_smd}}) \cr
#'   \code{ss.power.dt} \tab \code{\link{ss_buc_paired_t}} (alias \code{\link{ss_buc_smd_paired}}) \cr
#'   \code{ss.power.ba} \tab \code{\link{ss_buc_one_way_anova}} (when \code{levels.B} is \code{NULL}) or \code{\link{ss_buc_factorial_anova}} \cr
#'   \code{ss.power.ba.general} \tab \code{\link{ss_buc_factorial_anova_general}} \cr
#'   \code{ss.power.wa} \tab \code{\link{ss_buc_rm_anova}} \cr
#'   \code{ss.power.wa.general} \tab \code{\link{ss_buc_rm_anova_general}} \cr
#'   \code{ss.power.spa} \tab \code{\link{ss_buc_mixed_anova}} \cr
#'   \code{ss.power.spa.general} \tab \code{\link{ss_buc_mixed_anova_general}} \cr
#'   \code{ss.power.reg1} \tab \code{\link{ss_buc_reg_coef}} \cr
#'   \code{ss.power.reg.all} \tab \code{\link{ss_buc_R2}} \cr
#'   \code{ss.power.reg.joint} \tab \code{\link{ss_buc_reg_joint}} \cr
#' }
#'
#'   Because BUCSS 2.0.0 finds the noncentrality parameter by direct root
#'   finding rather than the fixed 0 to 100 grid used in 1.x, the values these
#'   wrappers return can differ slightly from 1.x and are no longer capped when
#'   the implied parameter exceeds 100. The \code{step} argument is accepted for
#'   call compatibility and ignored.
#'
#' @param t.observed,F.observed The observed \emph{t} or \emph{F} statistic from
#'   the prior study (the 1.x spelling of \code{t_observed} / \code{F_observed}).
#' @param N Total sample size of the prior study (number of pairs for
#'   \code{ss.power.dt}).
#' @param n Per-group sample size(s) of the prior study, for \code{ss.power.it}.
#' @param p,p.joint Number of predictors, and the number tested jointly, for the
#'   regression planners.
#' @param levels.A,levels.B Number of levels of the first and (optional) second
#'   factor, for \code{ss.power.ba} and \code{ss.power.wa}.
#' @param levels.between,levels.within Number of levels of the between- and
#'   within-subjects factors, for \code{ss.power.spa}.
#' @param cells Number of cells (groups) in the between-subjects design, for
#'   \code{ss.power.ba.general}.
#' @param df.numerator,df.denominator Numerator and denominator degrees of
#'   freedom for the effect of interest in the prior study, for the
#'   \code{.general} planners.
#' @param df.num.within Numerator degrees of freedom of the within-subjects
#'   component, for \code{ss.power.spa.general} with
#'   \code{effect = "between.within"}.
#' @param num.groups Number of between-subjects groups, for
#'   \code{ss.power.spa.general}.
#' @param effect The effect to plan for, in the 1.x dotted spelling (translated
#'   internally to the 2.x form).
#' @param alpha.prior,alpha.planned The prior-study and planned-study Type I
#'   error rates (the 1.x spelling of \code{alpha_prior} / \code{alpha_planned}).
#' @param assurance The desired assurance.
#' @param power The desired statistical power of the planned study.
#' @param step Ignored. Accepted only so that 1.x calls that set it do not error.
#'
#' @return A \code{\link{bucss_power}} object, as returned by the corresponding
#'   \code{ss_buc_*} function.
#'
#' @seealso The \code{ss_buc_*} replacements listed above.
#'
#' @name bucss-deprecated
#' @keywords internal
NULL

#' @rdname bucss-deprecated
#' @export
ss.power.it <- function(t.observed, n, N, alpha.prior = .05, alpha.planned = .05,
                        assurance = .80, power = .80, step = .001) {
  .deprecate_once("ss.power.it", "'ss_buc_independent_t'")
  ss_buc_independent_t(t_observed = t.observed, n = n, N = N,
                       alpha_prior = alpha.prior, alpha_planned = alpha.planned,
                       assurance = assurance, desired_power = power)
}

#' @rdname bucss-deprecated
#' @export
ss.power.dt <- function(t.observed, N, alpha.prior = .05, alpha.planned = .05,
                        assurance = .80, power = .80, step = .001) {
  .deprecate_once("ss.power.dt", "'ss_buc_paired_t'")
  ss_buc_paired_t(t_observed = t.observed, N = N,
                  alpha_prior = alpha.prior, alpha_planned = alpha.planned,
                  assurance = assurance, desired_power = power)
}

#' @rdname bucss-deprecated
#' @export
ss.power.ba <- function(F.observed, N, levels.A, levels.B = NULL,
                        effect = c("factor.A", "factor.B", "interaction"),
                        alpha.prior = .05, alpha.planned = .05,
                        assurance = .80, power = .80, step = .001) {
  .deprecate_once("ss.power.ba",
                  "'ss_buc_one_way_anova' (one factor) or 'ss_buc_factorial_anova' (two factors)")
  if (is.null(levels.B)) {
    ss_buc_one_way_anova(F_observed = F.observed, N = N, levels_A = levels.A,
                         alpha_prior = alpha.prior, alpha_planned = alpha.planned,
                         assurance = assurance, desired_power = power)
  } else {
    effect <- .translate_effect(effect, c("factor.A", "factor.B", "interaction"))
    ss_buc_factorial_anova(F_observed = F.observed, N = N, levels_A = levels.A,
                           levels_B = levels.B, effect = effect,
                           alpha_prior = alpha.prior, alpha_planned = alpha.planned,
                           assurance = assurance, desired_power = power)
  }
}

#' @rdname bucss-deprecated
#' @export
ss.power.ba.general <- function(F.observed, N, cells, df.numerator, df.denominator,
                                alpha.prior = .05, alpha.planned = .05,
                                assurance = .80, power = .80, step = .001) {
  .deprecate_once("ss.power.ba.general", "'ss_buc_factorial_anova_general'")
  ss_buc_factorial_anova_general(F_observed = F.observed, N = N, cells = cells,
                                 df_numerator = df.numerator,
                                 df_denominator = df.denominator,
                                 alpha_prior = alpha.prior,
                                 alpha_planned = alpha.planned,
                                 assurance = assurance, desired_power = power)
}

#' @rdname bucss-deprecated
#' @export
ss.power.wa <- function(F.observed, N, levels.A, levels.B = NULL,
                        effect = c("factor.A", "factor.B", "interaction"),
                        alpha.prior = .05, alpha.planned = .05,
                        assurance = .80, power = .80, step = .001) {
  .deprecate_once("ss.power.wa", "'ss_buc_rm_anova'")
  effect <- .translate_effect(effect, c("factor.A", "factor.B", "interaction"))
  ss_buc_rm_anova(F_observed = F.observed, N = N, levels_A = levels.A,
                  levels_B = levels.B, effect = effect,
                  alpha_prior = alpha.prior, alpha_planned = alpha.planned,
                  assurance = assurance, desired_power = power)
}

#' @rdname bucss-deprecated
#' @export
ss.power.wa.general <- function(F.observed, N, df.numerator,
                                alpha.prior = .05, alpha.planned = .05,
                                assurance = .80, power = .80, step = .001) {
  .deprecate_once("ss.power.wa.general", "'ss_buc_rm_anova_general'")
  ss_buc_rm_anova_general(F_observed = F.observed, N = N,
                          df_numerator = df.numerator,
                          alpha_prior = alpha.prior, alpha_planned = alpha.planned,
                          assurance = assurance, desired_power = power)
}

#' @rdname bucss-deprecated
#' @export
ss.power.spa <- function(F.observed, N, levels.between, levels.within,
                         effect = c("between", "within", "interaction"),
                         alpha.prior = .05, alpha.planned = .05,
                         assurance = .80, power = .80, step = .001) {
  .deprecate_once("ss.power.spa", "'ss_buc_mixed_anova'")
  effect <- .translate_effect(effect, c("between", "within", "interaction"))
  ss_buc_mixed_anova(F_observed = F.observed, N = N,
                     levels_between = levels.between, levels_within = levels.within,
                     effect = effect, alpha_prior = alpha.prior,
                     alpha_planned = alpha.planned, assurance = assurance,
                     desired_power = power)
}

#' @rdname bucss-deprecated
#' @export
ss.power.spa.general <- function(F.observed, N, df.numerator, num.groups,
                                 effect = c("between.only", "within.only", "between.within"),
                                 df.num.within, alpha.prior = .05,
                                 alpha.planned = .05, assurance = .80,
                                 power = .80, step = .001) {
  .deprecate_once("ss.power.spa.general", "'ss_buc_mixed_anova_general'")
  effect <- .translate_effect(effect,
                              c("between.only", "within.only", "between.within"))
  args <- list(F_observed = F.observed, N = N, df_numerator = df.numerator,
               num_groups = num.groups, effect = effect,
               alpha_prior = alpha.prior, alpha_planned = alpha.planned,
               assurance = assurance, desired_power = power)
  if (!missing(df.num.within)) args$df_num_within <- df.num.within
  do.call(ss_buc_mixed_anova_general, args)
}

#' @rdname bucss-deprecated
#' @export
ss.power.reg1 <- function(t.observed, N, p, alpha.prior = .05, alpha.planned = .05,
                          assurance = .80, power = .80, step = .001) {
  .deprecate_once("ss.power.reg1", "'ss_buc_reg_coef'")
  ss_buc_reg_coef(t_observed = t.observed, N = N, p = p,
                  alpha_prior = alpha.prior, alpha_planned = alpha.planned,
                  assurance = assurance, desired_power = power)
}

#' @rdname bucss-deprecated
#' @export
ss.power.reg.all <- function(F.observed, N, p, alpha.prior = .05, alpha.planned = .05,
                             assurance = .80, power = .80, step = .001) {
  .deprecate_once("ss.power.reg.all", "'ss_buc_R2'")
  ss_buc_R2(F_observed = F.observed, N = N, p = p,
            alpha_prior = alpha.prior, alpha_planned = alpha.planned,
            assurance = assurance, desired_power = power)
}

#' @rdname bucss-deprecated
#' @export
ss.power.reg.joint <- function(F.observed, N, p, p.joint, alpha.prior = .05,
                               alpha.planned = .05, assurance = .80, power = .80,
                               step = .001) {
  .deprecate_once("ss.power.reg.joint", "'ss_buc_reg_joint'")
  ss_buc_reg_joint(F_observed = F.observed, N = N, p = p, p_joint = p.joint,
                   alpha_prior = alpha.prior, alpha_planned = alpha.planned,
                   assurance = assurance, desired_power = power)
}

#' Positional extraction from a bias and uncertainty corrected sample size result
#'
#' @description Backward-compatible positional extraction for
#'   \code{\link{bucss_power}} objects. The deprecated 1.x functions (see
#'   \code{\link{bucss-deprecated}}) returned an unnamed two-element list,
#'   \code{list(sample_size, ncp)}, so code written for BUCSS 1.x reads
#'   \code{result[[1]]} and \code{result[[2]]}. The \code{ss_buc_*} functions
#'   return a richer object, but \code{result[[1]]} still yields the necessary
#'   sample size and \code{result[[2]]} the bias and uncertainty adjusted
#'   noncentrality parameter, so that 1.x extraction code keeps working.
#'
#'   Any other subscript (a column name, or a length other than one) falls
#'   through to the ordinary \code{data.frame} method, so \code{result[["value"]]}
#'   and \code{result$value} are unaffected.
#'
#' @param x A \code{bucss_power} object.
#' @param i A subscript. A single \code{1} or \code{2} selects the legacy scalar;
#'   anything else is passed to the \code{data.frame} method.
#' @param ... Passed to the \code{data.frame} method.
#'
#' @return For \code{i} equal to \code{1} or \code{2}, a single numeric value;
#'   otherwise whatever the \code{data.frame} \code{[[} method returns.
#'
#' @keywords internal
#' @export
`[[.bucss_power` <- function(x, i, ...) {
  if (!missing(i) && length(i) == 1L && is.numeric(i) && i %in% c(1, 2)) {
    term <- .subset2(x, "term")
    value <- .subset2(x, "value")
    if (i == 1) return(value[term %in% .BUCSS_SIZE_TERMS][1])
    return(value[term == "ncp_adjusted"])
  }
  NextMethod()
}

#' Subset a bias and uncertainty corrected sample size result
#'
#' Row and column subsetting returns a plain \code{data.frame}. The class is
#' dropped first so that \code{[.data.frame}'s internal column extraction
#' (which uses \code{[[} on columns 1 and 2) cannot collide with the legacy
#' positional \code{[[} contract above; without this, \code{head(result)} and
#' \code{result[i, ]} would silently return the legacy scalars in place of the
#' \code{term} and \code{value} columns. A subset is no longer the full
#' planning result (its attributes and print contract no longer apply), so the
#' honest return is an ordinary \code{data.frame}.
#'
#' @param x A \code{bucss_power} object.
#' @param ... Passed to the \code{data.frame} method.
#'
#' @return A \code{data.frame} (or vector, under the usual \code{drop} rules).
#'
#' @keywords internal
#' @export
`[.bucss_power` <- function(x, ...) {
  class(x) <- setdiff(class(x), "bucss_power")
  x[...]
}
