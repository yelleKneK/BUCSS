# Within subjects ANOVA (general case)

#' Necessary sample size to reach desired power for a within-subjects ANOVA with
#' any number of factors using an uncertainty and publication bias correction
#' procedure
#'
#' @description \code{ss.power.wa.general} returns the necessary per-group
#'   sample size to achieve a desired level of statistical power for a planned
#'   study testing any type of effect (omnibus, contrast) using a fully
#'   within-subjects ANOVA with any number of factors, based on information
#'   obtained from a previous study. The effect from the previous study can be
#'   corrected for publication bias and/or uncertainty to provide a sample size
#'   that will achieve more accurate statistical power for a planned study, when
#'   compared to approaches that use a sample effect size at face value or rely
#'   on sample size only. The bias and uncertainty adjusted previous study
#'   noncentrality parameter is also returned, which can be transformed to
#'   various effect size metrics.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss.power.wa.general} uses the observed
#'   \eqn{F}-value and sample size from a previous study to correct the
#'   noncentrality parameter associated with the effect of interest for
#'   publication bias and/or uncertainty. This new estimated noncentrality
#'   parameter is then used to calculate the necessary per-group sample size to
#'   achieve the desired level of power in the planned study.
#'
#'   The approach uses a likelihood function of a truncated non-central F
#'   distribution, where the truncation occurs due to small effect sizes being
#'   unobserved due to publication bias. The numerator of the likelihood
#'   function is simply the density of a noncentral F distribution. The
#'   denominator is the power of the test, which serves to truncate the
#'   distribution. Thus, the ratio of the numerator and the denominator is a
#'   truncated noncentral F distribution. (See Taylor & Muller, 1996, Equation
#'   2.1. and Anderson & Maxwell, 2017 for more details.)
#'
#'   Assurance is the proportion of times that power will be at or above the
#'   desired level, if the experiment were to be reproduced many times. For
#'   example, assurance = .5 means that power will be above the desired level
#'   half of the time, but below the desired level the other half of the time.
#'   Selecting assurance = .5 (selecting the noncentrality parameter at the 50th
#'   percentile of the likelihood distribution) results in a median-unbiased
#'   estimate of the population noncentrality parameter and does not correct for
#'   uncertainty. In order to correct for uncertainty, assurance > .5
#'   can be selected, which corresponds to selecting the noncentrality parameter
#'   associated with the (1 - assurance) quantile of the likelihood
#'   distribution.
#'
#'   If the previous study of interest has not been subjected to publication
#'   bias (e.g., a pilot study), \code{alpha.prior} can be set to 1 to indicate
#'   no publication bias. Alternative \eqn{\alpha}-levels can also be
#'   accommodated to represent differing amounts of publication bias. For
#'   example, setting \code{alpha.prior}=.20 would reflect less severe
#'   publication bias than the default of .05. In essence, setting
#'   \code{alpha.prior} at .20 assumes that studies with \eqn{p}-values less
#'   than .20 are published, whereas those with larger \eqn{p}-values are not.
#'
#'   In some cases, the corrected noncentrality parameter for a given level of
#'   assurance will be estimated to be zero. This is an indication that, at the
#'   desired level of assurance, the previous study's effect cannot be
#'   accurately estimated due to high levels of uncertainty and bias. When this
#'   happens, subsequent sample size planning is not possible with the chosen
#'   specifications. Two alternatives are recommended. First, users can select a
#'   lower value of assurance (e.g. .8 instead of .95). Second, users can reduce
#'   the influence of publciation bias by setting \code{alpha.prior} at a value
#'   greater than .05. It is possible to correct for uncertainty only by setting
#'   \code{alpha.prior}=1 and choosing the desired level of assurance. We
#'   encourage users to make the adjustments as minimal as possible.
#'
#'   \code{ss.power.wa.general} assumes sphericity for the within-subjects
#'   effects.
#'
#' @param F.observed Observed F-value from a previous study used to plan sample
#'   size for a planned study
#' @param N Total sample size of the previous study
#' @param df.numerator Numerator degrees of freedom for the effect of interest
#' @param alpha.prior Alpha-level \eqn{\alpha} for the previous study or the
#'   assumed statistical significance necessary for publishing in the field; to
#'   assume no publication bias, a value of 1 can be entered 
#' @param alpha.planned Alpha-level (\eqn{\alpha}) assumed for the planned study
#' @param assurance Desired level of assurance, or the long run proportion of
#'   times that the planned study power will reach or surpass desired level
#'   (assurance > .5 corrects for uncertainty; assurance < .5 not recommended)
#' @param power Desired level of statistical power for the planned study
#' @param step Value used in the iterative scheme to determine the noncentrality
#'   parameter necessary for sample size planning (0 < step < .01) (users should
#'   not generally need to change this value; smaller values lead to more
#'   accurate sample size planning results, but unnecessarily small values will
#'   add unnecessary computational time)
#'
#' @return Suggested per-group sample size for planned study
#' Publication bias and uncertainty- adjusted prior study noncentrality parameter
#'
#' @export
#' @import stats
#'
#' @examples ss.power.wa.general(F.observed=6.5, N=80,  df.numerator=1,
#' alpha.prior=.05, alpha.planned=.05, assurance=.50, power=.80, step=.001)
#'
#' @author Samantha F. Anderson \email{samantha.f.anderson@asu.edu},
#' Ken Kelley \email{kkelley@@nd.edu}
#'
#' @references Anderson, S. F., & Maxwell, S. E. (2017).
#'   Addressing the 'replication crisis': Using original studies to design
#'   replication studies with appropriate statistical power. \emph{Multivariate
#'   Behavioral Research, 52,} 305-322.
#'
#'   Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample size
#'   planning for more accurate statistical power: A method correcting sample
#'   effect sizes for uncertainty and publication bias. \emph{Psychological
#'   Science, 28,} 1547-1562.
#'
#'   Taylor, D. J., & Muller, K. E. (1996). Bias in linear model power and
#'   sample size calculation due to estimating noncentrality.
#'   \emph{Communications in Statistics: Theory and Methods, 25,} 1595-1610.

ss.power.wa.general <- function(F.observed, N, df.numerator, alpha.prior=.05, alpha.planned=.05, assurance=.80, power=.80, step=.001)
{
  if(alpha.prior > 1 | alpha.prior <= 0) stop("There is a problem with 'alpha' of the prior study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).")
  if(alpha.prior == 1) {alpha.prior <- .999 }
  if(alpha.planned >= 1 | alpha.planned <= 0) stop("There is a problem with 'alpha' of the planned study (i.e., the Type I error rate), please specify as a value between 0 and 1 (the default is .05).")

  if(assurance >= 1)
  {
    assurance <- assurance/100
  }

  if(assurance<0 | assurance>1)
  {
    stop("There is a problem with 'assurance' (i.e., the proportion of times statistical power is at or above the desired value), please specify as a value between 0 and 1 (the default is .80).")
  }
  
  if(assurance <.5)
  {
    warning( "THe assurance you have entered is < .5, which implies you will have under a 50% chance at achieving your desired level of power" )
  }

  if(power >= 1) power <- power/100

  if(power<0 | power>1) stop("There is a problem with 'power' (i.e., desired statistical power), please specify as a value between 0 and 1 (the default is .80).")

  if(missing(N)) stop("You must specify 'N', which is the total sample size.")

  NCP <- seq(from=0, to=100, by=step) # sequence of possible values for the noncentral parameter.

    n <- N
    df.denominator <- df.numerator*(n-1)

  f.density <- df(F.observed, df1=df.numerator, df2=df.denominator, ncp=NCP) # density of F using F observed
  critF <- qf(1-alpha.prior, df1=df.numerator, df2=df.denominator)

  if(F.observed <= critF) stop("Your observed F statistic is nonsignificant based on your specfied alpha of the prior study. Please increase 'alpha.prior' so 'F.observed' exceeds the critical value")

  power.values <- 1 - pf(critF, df1=df.numerator, df2=df.denominator, ncp = NCP) # area above critical F
  area.above.F <- 1 - pf(F.observed, df1=df.numerator, df2=df.denominator, ncp = NCP) # area above observed F
  area.area.between <- power.values - area.above.F

  TM <- area.area.between/power.values
  TM.Percentile <- min(NCP[which(abs(TM-assurance)==min(abs(TM-assurance)))])

  if(TM.Percentile==0) stop("The corrected noncentrality parameter is zero. Please either choose a lower value of assurance and/or a higher value of alpha for the prior study (e.g. accounting for less publication bias)")

  if (TM.Percentile > 0)
  {
    nrep <- 2
    denom.df <- df.numerator*(nrep-1)
    diff <- -1
    while (diff < 0 )
    {
      criticalF <- qf(1-alpha.planned, df1 = df.numerator, df2 = denom.df)
      powers <- 1 - pf(criticalF, df1 = df.numerator, df2 = denom.df, ncp = (nrep/n)*TM.Percentile)
      diff <- powers - power
      nrep <- nrep + 1
      denom.df <- df.numerator*(nrep-1)
    }
  }
  repn <- nrep-1
  return(list(repn, TM.Percentile)) # This is the same as total N needed
}
