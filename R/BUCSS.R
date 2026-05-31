#' BUCSS: Bias and Uncertainty Corrected Sample Size
#'
#' Implements sample size planning for a future study from the summary
#' statistics of a prior study, correcting the prior effect for publication
#' bias and uncertainty so the planned study attains its target statistical
#' power more reliably than methods that take a sample effect size at face
#' value. The correction builds on the truncated noncentral distribution of
#' Taylor and Muller (1996) and is described in Anderson, Kelley, and Maxwell
#' (2017).
#'
#' The functions are also available, through a web interface, as Shiny apps at
#' \url{https://designingexperiments.com}.
#'
#' @references
#' Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample size planning
#' for more accurate statistical power: A method correcting sample effect sizes
#' for uncertainty and publication bias. \emph{Psychological Science, 28,}
#' 1547--1562. \doi{10.1177/0956797617723724}
#'
#' @author Samantha F. Anderson (\email{samantha.f.anderson@@asu.edu}),
#'   Ken Kelley (\email{kkelley@@nd.edu}), and Scott E. Maxwell
#'
#' @import stats
#'
"_PACKAGE"
