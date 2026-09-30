# Monte Carlo companion to the planners: simulate the literature a plan came
# from and see how often the plan actually works. The per-design knowledge it
# needs lives in .BUCSS_DESIGN_SPECS, in R/bucss-design-registry.R.

# Draw m prior statistics from the design's sampling distribution at a true
# noncentrality parameter, on the scale the planner's statistic argument
# expects. The publication filter is applied by the caller.
.draw_prior_statistic <- function(spec, df, ncp, m) {
  switch(
    spec$family,
    t = rt(m, df = df, ncp = ncp),
    f = rf(m, df1 = df[1], df2 = df[2], ncp = ncp),
    # The single regression coefficient is corrected as an F on the squared t,
    # so the parameter is on the F scale and the statistic is returned as a t.
    f_from_t = sqrt(rf(m, df1 = 1, df2 = df[2], ncp = ncp)),
    # A correlation samples both variables: the noncentrality parameter is
    # itself random, proportional to the sampled spread of the predictor. See
    # .p_random_f() in R/bucss-internal.R. df is the prior N.
    correlation = sqrt(.draw_random_f(m, c(1, df - 2, df - 1, df), ncp)),
    # Regression with sampled predictors (random_predictors = TRUE): the same
    # mixture with q numerator df, residual df df[2], and m = df[3] chi square
    # df, at the prior N = df[4].
    f_random = .draw_random_f(m, df, ncp),
    f_from_t_random = sqrt(.draw_random_f(m, df, ncp)),
    chisq = rchisq(m, df = df, ncp = ncp),
    stop("Unknown sampling family.", call. = FALSE)
  )
}

# M draws of the random-predictor F statistic, one per study: the
# noncentrality parameter is itself random, the partial f squared (ncp/N)
# times a chi square on m df. df = c(q, df2, m, N). The chi square draws come
# first, then the F draws, as the correlation's draw always did.
.draw_random_f <- function(M, df, ncp) {
  lambda <- (ncp / df[4]) * rchisq(M, df = df[3])
  rf(M, df1 = df[1], df2 = df[2], ncp = lambda)
}

# The publication threshold on the same scale, and the proportion of studies
# that clear it (the prior study's own power at the true effect).
.prior_publication <- function(spec, df, ncp, alpha_prior) {
  switch(
    spec$family,
    t = list(crit = qt(1 - alpha_prior / 2, df = df),
             rate = (1 - pt(qt(1 - alpha_prior / 2, df), df, ncp)) +
               pt(-qt(1 - alpha_prior / 2, df), df, ncp)),
    f = list(crit = qf(1 - alpha_prior, df[1], df[2]),
             rate = 1 - pf(qf(1 - alpha_prior, df[1], df[2]), df[1], df[2], ncp)),
    f_from_t = list(crit = sqrt(qf(1 - alpha_prior, 1, df[2])),
                    rate = 1 - pf(qf(1 - alpha_prior, 1, df[2]), 1, df[2], ncp)),
    correlation = {
      crit_f <- qf(1 - alpha_prior, 1, df - 2)
      list(crit = sqrt(crit_f),
           rate = 1 - .p_correlation_f(crit_f, df, ncp / (ncp + df)))
    },
    f_random = {
      crit_f <- qf(1 - alpha_prior, df[1], df[2])
      list(crit = crit_f,
           rate = 1 - .p_random_f(crit_f, df[1], df[2], df[3],
                                  1 - ncp / (ncp + df[4])))
    },
    f_from_t_random = {
      crit_f <- qf(1 - alpha_prior, 1, df[2])
      list(crit = sqrt(crit_f),
           rate = 1 - .p_random_f(crit_f, 1, df[2], df[3],
                                  1 - ncp / (ncp + df[4])))
    },
    chisq = list(crit = qchisq(1 - alpha_prior, df),
                 rate = 1 - pchisq(qchisq(1 - alpha_prior, df), df, ncp))
  )
}

# Smallest statistic whose corrected parameter reaches 'target_ncp'. The map
# from statistic to corrected parameter is monotone increasing (a refusal is
# read as a corrected parameter of zero, its limiting value), so the bracket
# is found by walking out from the observed statistic and the root by
# bisection. Used to obtain the sample size the true effect itself requires,
# which is what "did this plan attain the target" is measured against.
.stat_for_ncp <- function(run_ncp, stat_observed, target_ncp) {
  lo <- hi <- abs(stat_observed)
  if (run_ncp(lo) > target_ncp) {
    for (i in 1:200) { lo <- lo / 1.2; if (run_ncp(lo) <= target_ncp) break }
    if (run_ncp(lo) > target_ncp) return(NA_real_)
  } else {
    for (i in 1:200) { hi <- hi * 1.2; if (run_ncp(hi) > target_ncp) break }
    if (run_ncp(hi) <= target_ncp) return(NA_real_)
  }
  uniroot(function(s) run_ncp(s) - target_ncp, lower = lo, upper = hi,
          tol = .Machine$double.eps^0.5)$root
}

#' Simulate the Literature a Plan Came From and Audit How Often It Works
#'
#' @description \code{ss_buc_sensitivity} takes a plan returned by any
#'   \code{ss_buc_*} function and asks the question the method's assurance
#'   claim is really about: if the true effect were a particular value, and a
#'   literature published only the prior studies that reached significance,
#'   how often would a plan built this way actually reach the target power?
#'
#' @details The bias and uncertainty correction rests on a claim about the long
#'   run. At \code{assurance = .80}, eight replications in ten should reach the
#'   target power. That claim is a statement about a whole literature, and a
#'   user planning one study never sees it exercised.
#'
#'   This function exercises it on the user's own numbers. It reads the design
#'   and the planning inputs off the plan you pass it, simulates prior studies
#'   of that design from a true effect \emph{you} name, discards the ones a
#'   literature would not have published (those failing to reach significance
#'   at \code{alpha_prior}), runs the same planner on each survivor, and
#'   reports what happened.
#'
#'   \strong{You must name the true effect.} There is no default, deliberately.
#'   Defaulting to the corrected estimate would be convenient and circular: it
#'   would ask the method to grade its own answer. \code{true_ncp} is on the
#'   same scale as the \code{ncp_adjusted} the plan reports, so the plan's own
#'   value is the natural anchor. Run the function at that value, at a fraction
#'   of it, and at a multiple of it; the interesting reading is how fast the
#'   attainment rate falls when the truth is smaller than the corrected
#'   estimate.
#'
#'   \strong{Reading the output.} A plan attains the target when its sample
#'   size is at least the sample size the true effect actually requires, which
#'   is reported as \code{size_at_true_effect}. Two attainment rates are given,
#'   and the difference between them is not a rounding detail:
#'
#'   \itemize{
#'     \item \code{attainment} counts only the prior studies for which the
#'       method issued a plan. This is what a user experiences.
#'     \item \code{attainment_with_refusals} counts a refusal (the planner's
#'       stop when the corrected noncentrality parameter is zero, or, for a
#'       design the planner reduces to a balanced one, its stop when a
#'       published statistic is not significant at the balanced design's
#'       degrees of freedom) as attaining, because declining to plan
#'       is the limiting case of an arbitrarily conservative plan. This is the
#'       quantity the assurance guarantee is about, and it should land on
#'       \code{assurance}.
#'   }
#'
#'   When the refusal rate is high the first number sits well below
#'   \code{assurance}, and that is not a defect: it is the visible price of the
#'   method refusing rather than guessing. Its expected value is
#'   \code{(assurance - refusal_rate)/(1 - refusal_rate)}.
#'
#'   \strong{What this can and cannot detect.} The prior studies are drawn from
#'   the sampling distribution the method assumes for that design, so this is a
#'   check that the correction and the planned sample size are internally
#'   right, and a demonstration of what publication bias does. It cannot detect
#'   a violated distributional assumption, since the simulation shares that
#'   assumption. For \code{\link{ss_buc_welch_t}}, in particular, it holds the
#'   variance ratio fixed at the value you supplied. For the regression
#'   planners the simulated prior studies follow the plan's own
#'   \code{random_predictors} setting, so audit a plan in the frame your
#'   research actually uses. For a within-subjects or split-plot plan made
#'   with \code{epsilon} below 1, the prior studies are drawn under the same
#'   approximation to a departure from sphericity that the planner corrects
#'   under, so the audit checks the method's internal consistency, not the
#'   approximation itself; the planner's help page reports how well the
#'   approximation holds against simulated raw data.
#'
#'   \strong{The package-wide verification ships with the package.} This
#'   function audits one plan; the harness that verifies the planners
#'   themselves, including checks against raw simulated data that do not share
#'   the method's distributional assumptions (bivariate normal samples for the
#'   correlation planner, multivariate normal samples for a regression
#'   coefficient with \code{random_predictors = TRUE}, two-group normal
#'   samples for the Welch planner), is
#'   installed with \code{BUCSS} and can be rerun in a few seconds with
#'   \code{source(system.file("tools", "verify_by_simulation.R", package =
#'   "BUCSS"))}.
#'
#' @param object A plan returned by any \code{ss_buc_*} function.
#' @param true_ncp The assumed true noncentrality parameter of the prior
#'   study's design, on the same scale as the \code{ncp_adjusted} row of
#'   \code{object}. Required; there is deliberately no default.
#' @param replications Number of published prior studies to simulate. At an
#'   assurance of .80 the default of 1000 gives
#'   \code{attainment_with_refusals} a Monte Carlo standard error near .013,
#'   and 10000 about .004; both are reported. The standard error of
#'   \code{attainment}, which counts only the issued plans, is larger when many
#'   prior studies are refused.
#' @param seed Optional random number seed. \code{NULL} (the default) draws
#'   from the caller's own random number stream, advancing it as any
#'   simulation does. When a seed is supplied, the caller's random number state
#'   and generator are restored on exit (and no state is left behind if the
#'   caller had none).
#'
#' @return An object of class \code{bucss_sensitivity}, a \code{data.frame}
#'   with columns \code{term} and \code{value} carrying the assumed true
#'   parameter and the plan's own \code{ncp_adjusted} for comparison, the
#'   number of replications, the publication and refusal rates, the two
#'   attainment rates each with its Monte Carlo standard error, the sample size
#'   the true effect requires, and the quartiles of the sample sizes the
#'   simulated plans called for. The design, the effect tested (for designs
#'   with one), the sample size unit, the requested \code{assurance}, and the
#'   requested \code{desired_power} travel as attributes. If every simulated
#'   prior study is refused, \code{attainment}, its standard error, and the
#'   quartiles are \code{NA}, with a warning.
#'
#' @export
#'
#' @examples
#' plan <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .80)
#' ncp <- plan$value[plan$term == "ncp_adjusted"]
#' ncp
#'
#' # Take the corrected estimate at face value as the truth. (200 replications
#' # keep the example quick: the rate counting refusals is good to about .03,
#' # the rate over issued plans to about .06. A reported audit deserves 10000.)
#' ss_buc_sensitivity(plan, true_ncp = ncp, replications = 200, seed = 113)
#'
#' # And ask what happens if the truth is smaller than the correction assumed.
#' ss_buc_sensitivity(plan, true_ncp = .6 * ncp, replications = 200, seed = 113)
#'
#' @seealso \code{\link{ss_buc_independent_t}} and the other \code{ss_buc_*}
#'   planners, whose results this function takes as input.
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu})
#'
#' @template references
ss_buc_sensitivity <- function(object, true_ncp, replications = 1000L,
                               seed = NULL) {
  if (!inherits(object, "bucss_plan"))
    stop("'object' must be a plan returned by one of the ss_buc_* functions.", call. = FALSE)
  if (missing(true_ncp))
    stop("You must specify 'true_ncp', the assumed true noncentrality parameter of the prior study's design. It is on the same scale as the 'ncp_adjusted' row of 'object'; there is deliberately no default, because defaulting to the corrected estimate would ask the method to grade its own answer.", call. = FALSE)
  .check_scalar_finite(true_ncp, "true_ncp")
  if (true_ncp <= 0) stop("'true_ncp' must be a positive number.", call. = FALSE)
  .check_count(replications, "replications", min = 2)

  .check_plan_shape(object)
  .bucss_local_seed(seed)

  design <- attr(object, "design")
  spec <- .BUCSS_DESIGN_SPECS[[design]]
  if (is.null(spec))
    stop("No simulation is defined for the design '", design, "'.", call. = FALSE)

  inputs <- .bucss_design_inputs(object)
  # The regression entries name a fixed and a random sampling family; the
  # plan's own random_predictors row picks the one it was computed in.
  spec <- .bucss_resolve_spec(spec, inputs)
  effect <- attr(object, "effect")
  df <- if (length(formals(spec$df)) > 1L) spec$df(inputs, effect) else spec$df(inputs)

  planner <- get(spec$planner, envir = asNamespace("BUCSS"))
  args <- .bucss_design_args(inputs, spec, effect)

  # Runtime check that the call was rebuilt correctly: replaying it with the
  # observed statistic must reproduce the plan exactly. This is what protects
  # against a silently wrong reconstruction.
  # Only the computed rows are compared: a correlation stated in r is replayed
  # in t, so the echoed input rows legitimately carry a different label.
  # The replay repeats any warning the plan raised when it was built (an
  # assurance below .5, say), which the user has already seen.
  replay <- withCallingHandlers(do.call(planner, args),
                                warning = function(w) invokeRestart("muffleWarning"))
  if (!.replay_matches(replay, object))
    stop("The planning inputs stored in 'object' do not reproduce it when replayed, so this plan cannot be simulated. Please report this, naming the design '", design, "'.", call. = FALSE)

  alpha_prior <- inputs$alpha_prior
  if (alpha_prior == 1) alpha_prior <- .999
  # A within-subjects plan made under a sphericity correction draws its prior
  # studies under the same approximation: the registry's df are already
  # scaled by epsilon, and the noncentrality parameter is scaled here, so
  # 'true_ncp' stays on the scale of the plan's 'ncp_adjusted'. The factor is
  # 1 for every other design, and for a within-subjects plan that assumed
  # sphericity, so their draws are unchanged.
  ncp_scale <- if (is.null(spec$ncp_scale)) 1 else spec$ncp_scale(inputs)
  pub <- .prior_publication(spec, df, true_ncp * ncp_scale, alpha_prior)

  # One planner run on a candidate statistic: the corrected parameter, the
  # sample size, or NA when the correction refuses. Only the planner's own
  # refusals (the zero noncentrality parameter stop, and the nonsignificance
  # stop a simulated study can meet when the planner rounds an unbalanced
  # design) count as refusals; any other error stops the audit with its
  # message rather than being counted as a refusal, which would count it as
  # attaining. Planner warnings are collected and reported once.
  # Warnings are counted separately for the search for the sample size the
  # true effect requires and for the simulated prior studies themselves, and
  # a planner run that raises several warnings is counted once, so the count
  # is of runs (studies), never more than there were.
  warned <- new.env(parent = emptyenv())
  warned$phase <- "search"
  warned$n <- c(search = 0L, studies = 0L)
  warned$first <- NULL
  run <- function(stat) {
    a <- args
    a[[spec$statistic]] <- stat
    warned$this_run <- FALSE
    out <- withCallingHandlers(
      tryCatch(do.call(planner, a),
               bucss_zero_ncp = function(e) NULL,
               bucss_nonsignificant = function(e) NULL,
               error = function(e) stop("While replaying the planner: ",
                                        conditionMessage(e), call. = FALSE)),
      warning = function(w) {
        if (!warned$this_run) {
          warned$this_run <- TRUE
          warned$n[[warned$phase]] <- warned$n[[warned$phase]] + 1L
        }
        if (is.null(warned$first)) warned$first <- conditionMessage(w)
        invokeRestart("muffleWarning")
      })
    if (is.null(out)) return(c(NA_real_, NA_real_))
    c(out$value[out$term %in% .BUCSS_SIZE_TERMS][1],
      out$value[out$term == "ncp_adjusted"])
  }
  run_ncp <- function(stat) { v <- run(stat)[2]; if (is.na(v)) 0 else v }

  # The sample size the true effect itself requires: run the planner on the
  # statistic whose corrected parameter is exactly 'true_ncp'.
  stat_true <- .stat_for_ncp(run_ncp, args[[spec$statistic]], true_ncp)
  if (is.na(stat_true))
    stop("Could not determine the sample size that 'true_ncp' requires for this design. Try a 'true_ncp' closer to the plan's own 'ncp_adjusted'.", call. = FALSE)
  size_at_true <- run(stat_true)[1]
  if (is.na(size_at_true))
    stop("At 'true_ncp' = ", signif(true_ncp, 3), " the planner refuses to plan (the corrected parameter is zero), so the sample size the true effect requires is undefined. Choose a larger 'true_ncp'.", call. = FALSE)

  # Simulate published prior studies. Oversample by the publication rate, in
  # blocks, so the requested number of published studies is actually reached.
  drawn <- numeric(0)
  budget <- 0L
  block <- max(1000L, ceiling(replications / max(pub$rate, .01)))
  while (length(drawn) < replications && budget < 500L * replications) {
    s <- .draw_prior_statistic(spec, df, true_ncp * ncp_scale, block)
    budget <- budget + block
    drawn <- c(drawn, s[abs(s) > pub$crit])
  }
  if (length(drawn) < replications)
    stop("Fewer than 'replications' simulated prior studies reached significance at 'alpha_prior'. The publication rate at this 'true_ncp' is about ", signif(pub$rate, 2), "; raise 'true_ncp' or lower 'replications'.", call. = FALSE)
  drawn <- drawn[seq_len(replications)]

  warned$phase <- "studies"
  res <- vapply(drawn, run, numeric(2))
  sizes <- res[1, ]
  refused <- is.na(sizes)
  issued <- sizes[!refused]

  attained <- issued >= size_at_true
  refusal_rate <- mean(refused)
  if (length(issued) == 0L) {
    warning("Every simulated prior study was refused, so 'attainment', its standard error, and the size quartiles are undefined and reported as NA.", call. = FALSE)
    attainment <- NA_real_
    quart <- rep(NA_real_, 3)
  } else {
    attainment <- mean(attained)
    quart <- quantile(issued, c(.25, .5, .75), names = FALSE)
  }
  if (sum(warned$n) > 0L)
    warning("The planner raised warnings on ", warned$n[["studies"]], " of the ",
            replications, " simulated prior studies and ", warned$n[["search"]],
            " times while finding the sample size the true effect requires; ",
            "the first was: ", warned$first, call. = FALSE)
  # A refusal counts as attaining: declining to plan is the limiting case of an
  # arbitrarily conservative plan, and it is the assurance guarantee's own
  # convention.
  with_refusals <- refused
  with_refusals[!refused] <- attained
  attainment_with_refusals <- mean(with_refusals)

  out <- data.frame(
    term = c("true_ncp", "ncp_adjusted", "replications", "publication_rate",
             "refusal_rate", "attainment", "attainment_mcse",
             "attainment_with_refusals", "attainment_with_refusals_mcse",
             "size_at_true_effect", "size_q25", "size_median", "size_q75"),
    value = c(true_ncp, object$value[object$term == "ncp_adjusted"],
              replications, pub$rate, refusal_rate, attainment,
              if (length(issued) == 0L) NA_real_ else
                sqrt(attainment * (1 - attainment) / length(issued)),
              attainment_with_refusals,
              sqrt(attainment_with_refusals * (1 - attainment_with_refusals) /
                     replications),
              size_at_true, quart[1], quart[2], quart[3]),
    stringsAsFactors = FALSE)

  attr(out, "design") <- design
  attr(out, "effect") <- effect
  attr(out, "sample_size_unit") <- attr(object, "sample_size_unit")
  attr(out, "assurance") <- inputs$assurance
  attr(out, "desired_power") <- inputs$desired_power
  class(out) <- c("bucss_sensitivity", "data.frame")
  out
}

#' Print a Simulated Sensitivity Audit of a Plan
#'
#' Prints the aligned \code{term}/\code{value} table, then factual footer
#' lines naming the design, the sample size unit, and the requested
#' \code{assurance} and \code{desired_power} the simulated rates should be
#' read against.
#'
#' @param x An object of class \code{bucss_sensitivity}.
#' @param digits Significant figures for non-integer values; defaults to
#'   \code{getOption("bucss.digits", 3)}.
#' @param ... Further arguments, ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' plan <- ss_buc_paired_t(t_observed = 3, N = 30, assurance = .80)
#' ss_buc_sensitivity(plan, true_ncp = 2, replications = 200, seed = 113)
#'
#' @aliases bucss_sensitivity
#' @export
#' @keywords internal
print.bucss_sensitivity <- function(x, digits = getOption("bucss.digits", 3L), ...) {
  # An object that is no longer one audit (a subset, or two audits stacked,
  # which keep the class) prints as the plain data frame it is, without footers
  # describing only one of them.
  one_audit <- is.data.frame(x) && identical(names(x), c("term", "value")) &&
    is.numeric(x[["value"]]) && !is.null(attr(x, "design")) &&
    all(vapply(c("true_ncp", "attainment", "attainment_with_refusals"),
               function(nm) sum(x[["term"]] == nm) == 1L, logical(1)))
  if (!one_audit) {
    print(structure(x, class = "data.frame"))
    return(invisible(x))
  }
  shown <- data.frame(term = x$term,
                      value = .format_bucss_value(x$value, digits),
                      stringsAsFactors = FALSE)
  print.data.frame(shown, row.names = FALSE, right = FALSE)

  design <- attr(x, "design")
  effect <- attr(x, "effect")
  cat("\n")
  cat("Design: ", design,
      if (!is.null(effect)) paste0(" (", effect, ")"), "\n", sep = "")
  cat("Sample size unit: ", attr(x, "sample_size_unit"), "\n", sep = "")
  cat("Requested assurance: ",
      sub("^0", "", sprintf("%.2f", attr(x, "assurance"))), "\n", sep = "")
  cat("Target power: ",
      sub("^0", "", sprintf("%.2f", attr(x, "desired_power"))), "\n", sep = "")
  invisible(x)
}
