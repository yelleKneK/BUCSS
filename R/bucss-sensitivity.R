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
    # .p_correlation_f() in R/bucss-internal.R.
    correlation = {
      N <- df
      lambda <- (ncp / N) * rchisq(m, df = N - 1)
      sqrt(rf(m, df1 = 1, df2 = N - 2, ncp = lambda))
    },
    chisq = rchisq(m, df = df, ncp = ncp),
    stop("Unknown sampling family.", call. = FALSE)
  )
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

#' Simulate the literature a plan came from and audit how often it works
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
#'     \item \code{attainment_with_refusals} counts a refusal (the zero
#'       noncentrality parameter stop) as attaining, because declining to plan
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
#'   variance ratio fixed at the value you supplied.
#'
#' @param object A plan returned by any \code{ss_buc_*} function.
#' @param true_ncp The assumed true noncentrality parameter of the prior
#'   study's design, on the same scale as the \code{ncp_adjusted} row of
#'   \code{object}. Required; there is deliberately no default.
#' @param replications Number of published prior studies to simulate. The
#'   default of 1000 gives a Monte Carlo standard error near .013 on the
#'   attainment rate; 10000 gives about .004.
#' @param seed Optional random number seed. \code{NULL} (the default) leaves
#'   the random number stream alone. When supplied, the caller's random number
#'   state is restored on exit.
#'
#' @return An object of class \code{bucss_sensitivity}, a \code{data.frame}
#'   with columns \code{term} and \code{value} carrying the assumed true
#'   parameter and the plan's own \code{ncp_adjusted} for comparison, the
#'   number of replications, the publication and refusal rates, the two
#'   attainment rates and the Monte Carlo standard error of the first, the
#'   sample size the true effect requires, and the quartiles of the sample
#'   sizes the simulated plans called for. The design, the sample size unit,
#'   the requested \code{assurance}, and the requested \code{desired_power}
#'   travel as attributes.
#'
#' @export
#'
#' @examples
#' plan <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .80)
#' plan$value[plan$term == "ncp_adjusted"]
#'
#' # Take the corrected estimate at face value as the truth.
#' ss_buc_sensitivity(plan, true_ncp = 1.9, replications = 200, seed = 2017)
#'
#' # And ask what happens if the truth is smaller than the correction assumed.
#' ss_buc_sensitivity(plan, true_ncp = 1.2, replications = 200, seed = 2017)
#'
#' @seealso \code{\link{ss_buc_independent_t}} and the other \code{ss_buc_*}
#'   planners, whose results this function takes as input.
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu})
#'
#' @template references
ss_buc_sensitivity <- function(object, true_ncp, replications = 1000L,
                               seed = NULL) {
  if (!inherits(object, "bucss_power"))
    stop("'object' must be a plan returned by one of the ss_buc_* functions.", call. = FALSE)
  if (missing(true_ncp))
    stop("You must specify 'true_ncp', the assumed true noncentrality parameter of the prior study's design. It is on the same scale as the 'ncp_adjusted' row of 'object'; there is deliberately no default, because defaulting to the corrected estimate would ask the method to grade its own answer.", call. = FALSE)
  .check_scalar_finite(true_ncp, "true_ncp")
  if (true_ncp <= 0) stop("'true_ncp' must be a positive number.", call. = FALSE)
  .check_count(replications, "replications", min = 2)

  if (!is.null(seed)) {
    .check_count(seed, "seed")
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  design <- attr(object, "design")
  spec <- .BUCSS_DESIGN_SPECS[[design]]
  if (is.null(spec))
    stop("No simulation is defined for the design '", design, "'.", call. = FALSE)

  inputs <- .bucss_design_inputs(object)
  effect <- attr(object, "effect")
  df <- if (length(formals(spec$df)) > 1L) spec$df(inputs, effect) else spec$df(inputs)

  planner <- get(spec$planner, envir = asNamespace("BUCSS"))
  args <- .bucss_design_args(inputs, spec, effect)

  # Runtime check that the call was rebuilt correctly: replaying it with the
  # observed statistic must reproduce the plan exactly. This is what protects
  # against a silently wrong reconstruction.
  # Only the computed rows are compared: a correlation stated in r is replayed
  # in t, so the echoed input rows legitimately carry a different label.
  computed <- c(.BUCSS_SIZE_TERMS, "total_N", "actual_power", "ncp_adjusted")
  replay <- do.call(planner, args)
  got <- replay$value[replay$term %in% computed]
  want <- object$value[object$term %in% computed]
  if (length(got) != length(want) ||
      !isTRUE(all.equal(got, want, tolerance = 1e-10)))
    stop("The planning inputs stored in 'object' do not reproduce it when replayed, so this plan cannot be simulated. Please report this, naming the design '", design, "'.", call. = FALSE)

  alpha_prior <- inputs$alpha_prior
  if (alpha_prior == 1) alpha_prior <- .999
  pub <- .prior_publication(spec, df, true_ncp, alpha_prior)

  # One planner run on a candidate statistic: the corrected parameter, the
  # sample size, or NA when the correction refuses.
  run <- function(stat) {
    a <- args
    a[[spec$statistic]] <- stat
    out <- tryCatch(do.call(planner, a), error = function(e) NULL)
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

  # Simulate published prior studies. Oversample by the publication rate, in
  # blocks, so the requested number of published studies is actually reached.
  drawn <- numeric(0)
  budget <- 0L
  block <- max(1000L, ceiling(replications / max(pub$rate, .01)))
  while (length(drawn) < replications && budget < 500L * replications) {
    s <- .draw_prior_statistic(spec, df, true_ncp, block)
    budget <- budget + block
    drawn <- c(drawn, s[abs(s) > pub$crit])
  }
  if (length(drawn) < replications)
    stop("Fewer than 'replications' simulated prior studies reached significance at 'alpha_prior'. The publication rate at this 'true_ncp' is about ", signif(pub$rate, 2), "; raise 'true_ncp' or lower 'replications'.", call. = FALSE)
  drawn <- drawn[seq_len(replications)]

  res <- vapply(drawn, run, numeric(2))
  sizes <- res[1, ]
  refused <- is.na(sizes)
  issued <- sizes[!refused]

  attained <- issued >= size_at_true
  attainment <- mean(attained)
  refusal_rate <- mean(refused)
  quart <- quantile(issued, c(.25, .5, .75), names = FALSE)
  # A refusal counts as attaining: declining to plan is the limiting case of an
  # arbitrarily conservative plan, and it is the assurance guarantee's own
  # convention.
  with_refusals <- refused
  with_refusals[!refused] <- attained

  out <- data.frame(
    term = c("true_ncp", "ncp_adjusted", "replications", "publication_rate",
             "refusal_rate", "attainment", "attainment_mcse",
             "attainment_with_refusals", "size_at_true_effect",
             "size_q25", "size_median", "size_q75"),
    value = c(true_ncp, object$value[object$term == "ncp_adjusted"],
              replications, pub$rate, refusal_rate, attainment,
              sqrt(attainment * (1 - attainment) / length(issued)),
              mean(with_refusals),
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

#' Print a simulated sensitivity audit of a plan
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
#' ss_buc_sensitivity(plan, true_ncp = 2, replications = 200, seed = 2017)
#'
#' @aliases bucss_sensitivity
#' @export
#' @keywords internal
print.bucss_sensitivity <- function(x, digits = getOption("bucss.digits", 3L), ...) {
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
