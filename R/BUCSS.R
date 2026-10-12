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
#' Each \code{ss_buc_*} function plans one analysis. The \emph{t} tests:
#' \code{\link{ss_buc_one_sample_t}}, \code{\link{ss_buc_paired_t}},
#' \code{\link{ss_buc_independent_t}}, and \code{\link{ss_buc_welch_t}}.
#' Between-subjects analysis of variance and covariance:
#' \code{\link{ss_buc_one_way_anova}}, \code{\link{ss_buc_factorial_anova}},
#' \code{\link{ss_buc_factorial_anova_general}}, and
#' \code{\link{ss_buc_ancova}}. Within-subjects and split-plot analysis of
#' variance: \code{\link{ss_buc_rm_anova}},
#' \code{\link{ss_buc_rm_anova_general}}, \code{\link{ss_buc_mixed_anova}},
#' and \code{\link{ss_buc_mixed_anova_general}}. Regression and correlation:
#' \code{\link{ss_buc_reg_coef}}, \code{\link{ss_buc_R2}},
#' \code{\link{ss_buc_reg_joint}}, and \code{\link{ss_buc_correlation}}. And
#' \code{\link{ss_buc_manova}} for a two-group multivariate comparison and
#' \code{\link{ss_buc_chisq_diff}} for the chi square difference between
#' nested models. \code{\link{ss_buc_sensitivity}} simulates the literature a
#' plan came from and reports how often such a plan attains its target power,
#' and \code{\link[=plot.bucss_plan]{plot}} draws the largest assurance a prior
#' result can support.
#'
#' The planners carried over from BUCSS 1.x are also available, through a web
#' interface, as Shiny apps at \url{https://designingexperiments.com}; the
#' designs added in 2.0.0 are available in R only.
#'
#' BUCSS follows the design and naming conventions of the \pkg{DMAR} package
#' (Design, Measurement, and Analysis in R): its \code{ss_buc_*} planners
#' parallel the \code{ss_aipe_*} (accuracy in parameter estimation) and
#' \code{ss_power_*} (power analysis) sample size planning families in
#' \pkg{DMAR}, returning a tidy result object in the same house style.
#' \pkg{DMAR} is a natural companion for effect size estimation, confidence
#' intervals, and sample size planning across additional frameworks.
#'
#' The bias and uncertainty adjusted noncentrality parameter is found by solving
#' the truncated-likelihood equation directly with \code{\link[stats]{uniroot}},
#' so it is no longer capped at 100 (the 1.x grid bound) and needs no
#' resolution argument.
#'
#' \strong{How assurance is used.} Every planner takes an \code{assurance}: the
#' long run proportion of planned studies whose power at the true effect
#' reaches or surpasses the desired power, when each is planned from a prior
#' study published under the rule \code{alpha_prior} describes. A planned study
#' reaches its target power whenever the noncentrality parameter it was planned
#' from is at or below the true one, because power rises with the effect (and,
#' apart from rounding the sample size up to a whole number and from the
#' smallest sample size a planner returns, the design's own minimum or the
#' Welch planner's 10 per group, only then). The planners therefore plan from
#' the lower limit of a one-sided confidence interval for the prior study's
#' noncentrality parameter, at confidence level \code{assurance}, computed from
#' the distribution of the prior statistic truncated at the publication
#' threshold: the value at which a published prior study would exceed the
#' observed statistic with probability 1 - \code{assurance} (what Anderson,
#' Kelley, and Maxwell, 2017, call the (1 - \code{assurance}) quantile of the
#' likelihood distribution). At \code{assurance = .5} this is the median
#' unbiased estimate, which corrects for publication bias only; above .5 it
#' adds the correction for uncertainty. The term follows Kelley and Maxwell
#' (2003), who used assurance in sample size planning for the probability that
#' a confidence interval will be no wider than desired; here it is the
#' probability that the planned study's power will be at least the desired
#' power. \code{\link{ss_buc_sensitivity}} checks the proportion by simulation
#' for any plan, counting a refusal to plan as meeting the target, as the
#' guarantee does, and also reporting the rate among the plans issued.
#'
#' \strong{The planning effect size.} Every plan also reports a planning
#' effect size, in a row whose name ends in \code{_planning}
#' (\code{smd_planning}, \code{f_planning}, \code{rho_planning}, and so on).
#' It is an adjusted effect size: the value, in the design's usual effect size
#' metric, used to plan the sample size given the situation at hand (the prior
#' study's result, the publication bias assumed through \code{alpha_prior},
#' and the uncertainty allowed for through \code{assurance}). It is
#' \code{ncp_adjusted} restated in that metric and is what the plan is built
#' on, not an estimate of the population effect to report: above an assurance
#' of .5 it is deliberately conservative (a lower confidence limit), and at .5
#' it is the median unbiased estimate. Each planner's Value section names its
#' metric.
#'
#' \strong{Safeguard power is a special case.} Safeguard power (Perugini,
#' Gallucci, & Costantini, 2014) is a way of planning a new study from a prior
#' study's effect size that guards against the prior estimate being too large.
#' Instead of planning for the effect size the prior study observed, the
#' researcher plans for the lower limit of a two-sided 60 percent confidence
#' interval around it (a limit that, in repeated sampling, falls below the
#' true effect 80 percent of the time), so the new study keeps its intended
#' power even if the true effect is somewhat smaller than the one observed. It
#' allows for uncertainty in the prior estimate but not for publication bias.
#'
#' It is the special case \code{alpha_prior = 1} and \code{assurance = .80}
#' of the method here, that is, BUCSS without the publication bias correction,
#' at an assurance of .80. \code{alpha_prior} is the significance level a
#' prior result had to reach to be published, and no \emph{p} value exceeds 1,
#' so \code{alpha_prior = 1} means no publication filter: any result would
#' have been published, the truncation is removed, and the prior statistic has
#' its ordinary distribution. At \code{assurance = .80} the adjusted
#' noncentrality parameter is then the lower limit of a one-sided 80 percent
#' confidence interval, which is the lower limit of the two-sided 60 percent
#' interval (either way, 20 percent lies below it), and the plan is the
#' smallest sample size whose ordinary power at that value reaches
#' \code{desired_power}. In general, the lower limit of a two-sided interval
#' at confidence level \emph{C} is \code{assurance = (1 + C)/2}, so a
#' two-sided 80 percent interval is \code{assurance = .90}.
#'
#' A few differences remain, and each can only make the plan here the same as
#' the safeguard plan or larger. First, the planners compute with .999 in
#' place of an \code{alpha_prior} of exactly 1, as BUCSS 1.x did, so plans
#' match earlier versions; that filter is so lenient (it drops only results
#' with \emph{p} above .999) that at an assurance of .80 it lowers the
#' adjusted noncentrality parameter by less than 0.1 percent for a prior
#' significant at .05. Second, for the \emph{t} tests and every other effect
#' on one degree of freedom the interval bounds the size of the effect
#' whatever its sign, because the publication rule is two-sided. Third, an
#' unbalanced prior design is read as a balanced one, and the Welch planner
#' returns at least 10 per group. All of the differences are negligible for a
#' convincing prior result and grow as the prior result weakens;
#' \code{vignette("understanding-assurance", package = "BUCSS")} measures
#' them. For the regression planners the equivalence holds in the frame
#' \code{random_predictors} names. With \code{alpha_prior} below 1 the same
#' construction is applied to the publication-truncated distribution, which is
#' what the bias correction adds to safeguard power.
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
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple regression:
#' Obtaining regression coefficients that are accurate, not simply
#' significant. \emph{Psychological Methods, 8,} 305--321.
#' \doi{10.1037/1082-989X.8.3.305}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#' experiments and analyzing data: A model comparison perspective} (4th ed.).
#' Routledge.
#'
#' Perugini, M., Gallucci, M., & Costantini, G. (2014). Safeguard power as a
#' protection against imprecise power estimates. \emph{Perspectives on
#' Psychological Science, 9,} 319--332. \doi{10.1177/1745691614528519}
#'
#' Taylor, D. J., & Muller, K. E. (1996). Bias in linear model power and
#' sample size calculation due to estimating noncentrality. \emph{Communications
#' in Statistics: Theory and Methods, 25,} 1595--1610.
#' \doi{10.1080/03610929608831787}
#'
#' @seealso \pkg{DMAR}
#'
#' @author Ken Kelley (\email{kkelley@@nd.edu}), Samantha F. Anderson
#'   (\email{samantha.f.anderson@@asu.edu}), and Scott E. Maxwell
#'
#' @importFrom graphics abline box lines mtext par points rect text
#' @importFrom stats dnbinom pbeta pchisq pf pt qchisq qf qnbinom qt quantile
#' @importFrom stats rchisq rf rt uniroot
#' @importFrom withr local_seed
#'
"_PACKAGE"
