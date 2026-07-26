#' @param alpha_prior Alpha level \eqn{\alpha} for the previous study or the
#'   assumed statistical significance necessary for publishing in the field; to
#'   assume no publication bias, a value of 1 can be entered.
#' @param alpha_planned Alpha level \eqn{\alpha} assumed for the planned study.
#' @param assurance Desired level of assurance, or the long run proportion of
#'   times that the planned study power will reach or surpass the desired level
#'   (assurance > .5 corrects for uncertainty; assurance < .5 is not
#'   recommended). Enter it as a proportion in (0, 1) or as a percentage greater
#'   than 1 (for example, 80 is read as .80); a value of exactly 1 is read as 1
#'   percent. A percentage is echoed in the result's planning-input rows as
#'   the coerced proportion (80 is echoed as .80); the same applies to
#'   \code{desired_power}.
#' @param desired_power Desired level of statistical power for the planned
#'   study. Enter it as a proportion in (0, 1) or as a percentage greater than
#'   1; a value of exactly 1 is read as 1 percent.
