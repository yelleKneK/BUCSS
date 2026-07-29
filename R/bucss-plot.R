# The assurance ceiling, drawn.
#
# The hardest thing to convey about this method is why a weak prior result
# cannot support a high assurance, and why the planners refuse rather than
# returning a large number. Explained in prose it sounds like a limitation.
# Drawn, it is obviously a wall: the corrected noncentrality parameter falls to
# zero at 1 - p/alpha_prior, and the sample size it implies grows without
# bound there. This method re-runs a plan's own planner across a range of
# assurance values and draws both curves.
#
# The per-design knowledge needed to re-run the planner (which function, and
# how to rebuild its arguments from the result's echoed rows) lives in
# .BUCSS_DESIGN_SPECS, in R/bucss-design-registry.R, shared with
# ss_buc_sensitivity().

# Assurance values to evaluate. The interesting behavior is all in the last
# fraction of a percent below the ceiling, so the grid is spaced evenly in the
# logarithm of the headroom (ceiling - assurance) rather than in assurance
# itself: half the points fall in the final few hundredths, which is what makes
# the wall visible instead of a single vertical jump between two grid points.
.assurance_grid <- function(ceiling_value, lower, n_points) {
  gap_max <- ceiling_value - lower
  gap_min <- min(1e-3, gap_max / 1000)
  gaps <- exp(seq(log(gap_max), log(gap_min), length.out = n_points))
  sort(ceiling_value - gaps)
}

#' Plot the assurance ceiling of a plan
#'
#' @description Re-runs the planner that produced \code{x} across a range of
#'   \code{assurance} values and draws what happens as the requested assurance
#'   approaches the largest value the prior result can support: the bias and
#'   uncertainty adjusted noncentrality parameter falling to zero, and the
#'   sample size it implies growing asymptotically.
#'
#' @details Every prior result has a largest supportable assurance, the closed
#'   form \eqn{1 - p/\alpha}{1 - p/alpha} where \emph{p} is the prior study's
#'   \emph{p} value and \eqn{\alpha} is \code{alpha_prior}. Request more than
#'   that and the corrected noncentrality parameter is driven to zero, no
#'   sample size attains the target power, and the planner stops. The ceiling
#'   is reported in the printed footer of every plan and named in that error,
#'   but its shape is the point: assurance does not degrade gracefully as it
#'   rises, it hits a wall.
#'
#'   Two panels are drawn, sharing an assurance axis.
#'
#'   \itemize{
#'     \item \strong{The corrected noncentrality parameter}, which falls to
#'       zero at the ceiling. This is the quantity the correction actually
#'       solves for, and the panel shows how little of it survives near the
#'       ceiling.
#'     \item \strong{The necessary sample size}, on a logarithmic scale,
#'       because it spans orders of magnitude over the plotted range. This is
#'       the same story in the units a researcher budgets in.
#'   }
#'
#'   The ceiling is drawn as a dashed vertical line in both panels and the
#'   plan's own assurance is marked with a filled point, so a reader can see at
#'   a glance how much headroom their prior result has left.
#'
#'   The grid is spaced evenly in the logarithm of the remaining headroom
#'   rather than in assurance, so that roughly half the plotted points fall in
#'   the last few hundredths below the ceiling. Spacing evenly in assurance
#'   would draw the wall as a single jump between two grid points and lose the
#'   shape entirely.
#'
#'   The plotted values are returned invisibly as an ordinary
#'   \code{data.frame}, one row per assurance, so the figure can be redrawn
#'   with any other graphics system. \code{BUCSS} draws it in base graphics and
#'   takes no plotting dependency of its own.
#'
#' @param x A plan returned by any \code{ss_buc_*} function.
#' @param which Which panel to draw: \code{"both"} (the default, side by
#'   side), \code{"ncp"} for the corrected noncentrality parameter alone, or
#'   \code{"size"} for the sample size alone.
#' @param assurance Assurance values to evaluate. \code{NULL} (the default)
#'   builds a grid running from \code{lower} up to just short of the ceiling.
#' @param lower Smallest assurance to evaluate, used only when
#'   \code{assurance} is \code{NULL}. Defaults to .5, the recommended floor,
#'   or to half the ceiling when the ceiling itself is below .5.
#' @param n_points Number of assurance values in the automatic grid.
#' @param col Color for the curve and the marked plan.
#' @param ... Further arguments passed to the underlying \code{plot} calls.
#'
#' @return Invisibly, a \code{data.frame} with columns \code{assurance},
#'   \code{ncp_adjusted}, \code{sample_size}, and \code{actual_power}, one row
#'   per evaluated assurance. The largest supportable assurance travels on the
#'   \code{assurance_ceiling} attribute, and the design and sample size unit on
#'   \code{design} and \code{sample_size_unit}.
#'
#' @examples
#' plan <- ss_buc_independent_t(t_observed = 3, n = 20, assurance = .80)
#'
#' # The ceiling, drawn: the corrected parameter falls to zero and the
#' # sample size grows asymptotically.
#' plot(plan)
#'
#' # The plotted values, for redrawing elsewhere.
#' curve_values <- plot(plan, which = "size", n_points = 20)
#' head(curve_values)
#' attr(curve_values, "assurance_ceiling")
#'
#' @seealso \code{\link{ss_buc_sensitivity}} to audit a plan by simulation, and
#'   \code{vignette("understanding-assurance", package = "BUCSS")} for what the
#'   ceiling means.
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}) and
#'   Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu})
#'
#' @template references
#'
#' @export
plot.bucss_power <- function(x, which = c("both", "ncp", "size"),
                             assurance = NULL, lower = NULL, n_points = 60L,
                             col = "#B03A2E", ...) {
  which <- match.arg(which)
  design <- attr(x, "design")
  spec <- .BUCSS_DESIGN_SPECS[[design]]
  if (is.null(spec))
    stop("No plot is defined for the design '", design, "'.", call. = FALSE)

  ceiling_value <- attr(x, "assurance_ceiling")
  if (is.null(ceiling_value) || !is.finite(ceiling_value))
    stop("This plan does not carry an assurance ceiling, so there is nothing to plot.", call. = FALSE)

  inputs <- .bucss_design_inputs(x)
  args <- .bucss_design_args(inputs, spec, attr(x, "effect"))
  planner <- get(spec$planner, envir = asNamespace("BUCSS"))

  # Replaying the rebuilt call must reproduce the plan, exactly as in
  # ss_buc_sensitivity(): every point drawn comes from that call with only
  # 'assurance' changed, so a wrong reconstruction would draw a wrong curve.
  computed <- c(.BUCSS_SIZE_TERMS, "total_N", "actual_power", "ncp_adjusted")
  replay <- do.call(planner, args)
  if (!isTRUE(all.equal(replay$value[replay$term %in% computed],
                        x$value[x$term %in% computed], tolerance = 1e-10)))
    stop("The planning inputs stored in 'x' do not reproduce it when replayed, so this plan cannot be plotted. Please report this, naming the design '", design, "'.", call. = FALSE)

  if (is.null(assurance)) {
    .check_count(n_points, "n_points", min = 3)
    if (is.null(lower)) lower <- if (ceiling_value > .5) .5 else ceiling_value / 2
    .check_scalar_finite(lower, "lower")
    if (lower <= 0 || lower >= ceiling_value)
      stop("'lower' must be greater than 0 and below the plan's assurance ceiling of ",
           signif(ceiling_value, 3), ".", call. = FALSE)
    assurance <- .assurance_grid(ceiling_value, lower, n_points)
  } else {
    if (!is.numeric(assurance) || length(assurance) < 2L || any(!is.finite(assurance)))
      stop("'assurance' must be a numeric vector of at least two finite values.", call. = FALSE)
    # Strictly inside (0, 1): the planners read an assurance of 1 or more as a
    # percentage, so an unguarded 1 here would quietly become .01.
    if (any(assurance <= 0 | assurance >= 1))
      stop("Every value of 'assurance' must be strictly between 0 and 1.", call. = FALSE)
    assurance <- sort(assurance)
  }

  # An assurance at or above the ceiling is a refusal, not a failure: it is the
  # wall the figure exists to show, and it is recorded as NA.
  at <- function(a) {
    out <- tryCatch(do.call(planner, c(args[names(args) != "assurance"],
                                       list(assurance = a))),
                    error = function(e) NULL, warning = function(w) NULL)
    if (is.null(out)) return(c(NA_real_, NA_real_, NA_real_))
    c(out$value[out$term %in% .BUCSS_SIZE_TERMS][1],
      out$value[out$term == "ncp_adjusted"],
      out$value[out$term == "actual_power"])
  }
  vals <- vapply(assurance, at, numeric(3))

  curve_values <- data.frame(assurance = assurance,
                             ncp_adjusted = vals[2, ],
                             sample_size = vals[1, ],
                             actual_power = vals[3, ])
  attr(curve_values, "assurance_ceiling") <- ceiling_value
  attr(curve_values, "design") <- design
  attr(curve_values, "sample_size_unit") <- attr(x, "sample_size_unit")

  plan_assurance <- inputs$assurance
  ok <- !is.na(curve_values$ncp_adjusted)
  if (!any(ok)) {
    warning("No assurance value in the plotted range yields a usable plan, so nothing was drawn.",
            call. = FALSE)
    return(invisible(curve_values))
  }

  if (which == "both") {
    old_par <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3.2, 1.1))
    on.exit(par(old_par), add = TRUE)
  }

  # The axis runs a little past the ceiling so that the region where no plan
  # exists is visible as a region rather than as the absence of one. That band
  # is the whole point of the figure.
  x_lo <- min(curve_values$assurance)
  x_hi <- ceiling_value + .10 * (ceiling_value - x_lo)

  # Arguments in '...' override the defaults rather than colliding with them,
  # so a caller can set their own main, ylab, or xlim.
  dots <- list(...)
  panel <- function(y, ylab, main, log = "") {
    pargs <- list(x = curve_values$assurance[ok], y = y[ok], type = "n",
                  log = log, xlab = "Assurance", ylab = ylab, main = main,
                  font.main = 1L, xlim = c(x_lo, x_hi))
    pargs[names(dots)] <- dots
    do.call(plot, pargs)
    usr <- par("usr")
    y_lo <- if (log == "y") 10^usr[3] else usr[3]
    y_hi <- if (log == "y") 10^usr[4] else usr[4]
    rect(ceiling_value, y_lo, usr[2], y_hi, col = "grey92", border = NA)
    text(mean(c(ceiling_value, usr[2])),
         if (log == "y") 10^mean(usr[3:4]) else mean(usr[3:4]),
         "no plan exists", srt = 90, cex = .7, col = "grey45")
    lines(curve_values$assurance[ok], y[ok], lwd = 2, col = col)
    abline(v = ceiling_value, lty = 2, col = "grey40")
    mtext(paste0("ceiling ", sub("^0", "", sprintf("%.3f", ceiling_value))),
          side = 3, at = ceiling_value, adj = 1.05, cex = .75, col = "grey40")
    if (!is.na(plan_assurance) && plan_assurance >= x_lo &&
        plan_assurance <= x_hi) {
      y_plan <- if (identical(ylab, "Adjusted noncentrality parameter"))
        x$value[x$term == "ncp_adjusted"] else
          x$value[x$term %in% .BUCSS_SIZE_TERMS][1]
      points(plan_assurance, y_plan, pch = 19, col = col, cex = 1.2)
    }
    box()
  }

  if (which %in% c("both", "ncp"))
    panel(curve_values$ncp_adjusted, "Adjusted noncentrality parameter",
          "The Correction Falls to Zero")
  if (which %in% c("both", "size"))
    panel(curve_values$sample_size,
          paste0("Necessary sample size (", attr(x, "sample_size_unit"), ")"),
          "The Sample Size Grows Asymptotically", log = "y")

  invisible(curve_values)
}
