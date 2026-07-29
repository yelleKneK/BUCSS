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
#' BUCSS follows the design and naming conventions of the \pkg{DMAR} package
#' (Design, Measurement, and Analysis of Human-Centered Research): its
#' \code{ss_buc_*} planners parallel the \code{ss_aipe_*} (accuracy in parameter
#' estimation) and \code{ss_power_*} (power analysis) sample size planning
#' families in \pkg{DMAR}, returning a tidy result object in the same house
#' style. \pkg{DMAR} is a natural companion for effect size estimation,
#' confidence intervals, and sample size planning across additional frameworks.
#'
#' The bias and uncertainty adjusted noncentrality parameter is found by solving
#' the truncated-likelihood equation directly with \code{\link[stats]{uniroot}},
#' so it is no longer capped at 100 (the 1.x grid bound) and needs no
#' resolution argument.
#'
#' @references
#' Anderson, S. F., & Kelley, K. (2024). Sample size planning for replication
#' studies: The devil is in the design. \emph{Psychological Methods, 29,}
#' 844--867. \doi{10.1037/met0000520}
#'
#' Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size planning
#' for more accurate statistical power: A method adjusting sample effect sizes
#' for publication bias and uncertainty. \emph{Psychological Science, 28,}
#' 1547--1562. \doi{10.1177/0956797617723724}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#' experiments and analyzing data: A model comparison perspective} (4th ed.).
#' Routledge.
#'
#' @seealso \pkg{DMAR}
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}), Samantha F. Anderson
#'   (\email{samantha.f.anderson@@asu.edu}), and Scott E. Maxwell
#'
#' @importFrom stats dnbinom pbeta pchisq pf pt qchisq qf qnbinom qt uniroot
#'
"_PACKAGE"
