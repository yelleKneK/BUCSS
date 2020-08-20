#' Necessary sample size to reach desired power for an independent t-test using
#' an uncertainty and publication bias correction procedure
#'
#' @description \code{ss.power.it} returns the necessary per-group sample size
#'   to achieve a desired level of statistical power for a planned study using
#'   an independent t-test, based on information obtained from a previous study.
#'   The effect from the previous study can be corrected for publication bias
#'   and/or uncertainty to provide a sample size that will achieve more accurate
#'   statistical power for a planned study, when compared to approaches that use
#'   a sample effect size at face value or rely on sample size only. The bias
#'   and uncertainty adjusted previous study noncentrality parameter is also
#'   returned, which can be transformed to various effect size metrics.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss.power.it} uses the observed
#'   \eqn{t}-value and sample size from a previous study to correct the
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
#'   distribution. In the two-group case, this formula reduces to the density of
#'   a truncated noncentral \eqn{t}-distribution.(See Taylor & Muller, 1996,
#'   Equation 2.1. and Anderson & Maxwell, 2017, for more details.)
#'
#'   Assurance is the proportion of times that power will be at or above the
#'   desired level, if the experiment were to be reproduced many times. For
#'   example, assurance = .5 means that power will be above the desired level
#'   half of the time, but below the desired level the other half of the time.
#'   Selecting assurance = .5 (selecting the noncentrality parameter at the 50th
#'   percentile of the likelihood distribution) results in a median-unbiased
#'   estimate of the population noncentrality parameter and does not corrects for
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
#'   \code{ss.power.it} assumes that the planned study will have equal n.
#'   Unequal n in the previous study is handled in the following way for the
#'   independent-t. If the user enters an odd value for N, no information is
#'   available on the exact group sizes. The function calculates n by dividing N
#'   by 2 and both rounding up and rounding down the result, thus assuming equal
#'   n. The suggested sample size for the planned study is calculated using both
#'   of these values of n, and the function returns the larger of these two
#'   suggestions, to be conservative. If the user enters a vector for n with two
#'   different values, specific information is available on the exact group
#'   sizes. n is calcualted as the harmonic mean of these two values (a measure
#'   of effective sample size). Again, this is rounded both up and down. The
#'   suggested sample size for the planned study is calculated using both of
#'   these values of n, and the function returns the larger of these two
#'   suggestions, to be conservative. The adjusted noncentrality parameter
#'   that is output is the lower of the two possibilities, again, to be
#'   conservative. When the individual group sizes of an unequal-n previous study
#'   are known, we highly encourage entering these explicitly, especially if the
#'   sample sizes are quite discrepant, as this allows for the greatest precision
#'   in estimating an appropriate planned study n.
#'
#' @param t.observed Observed \eqn{t}-value from a previous study used to plan
#'   sample size for a planned study
#' @param n Per group sample size (if equal) or the two group sample sizes of
#'   the previous study (enter either a single number or a vector of length 2)
#' @param N Total sample size of the previous study, assumed equal across groups
#'   if specified
#' @param alpha.prior Alpha-level \eqn{\alpha} for the previous study or the
#'   assumed statistical significance necessary for publishing in the field; to
#'   assume no publication bias, a value of 1 can be entered 
#' @param alpha.planned Alpha-level (\eqn{\alpha}) assumed for the planned study
#' @param assurance Desired level of assurance, or the long run proportion of
#'   times that the planned study power will reach or surpass desired level
#'   (assurance > .5 corrects for uncertainty; assurance < .5 not recommended)
#' @param power Desired level of statistical power for the planned study
#' @param step Value used in the iterative scheme to determine the noncentrality
#'   parameter necessary for sample size planning (0 < step < .01) (users
#'   should not generally need to change this value; smaller values lead to more
#'   accurate sample size planning results, but unnecessarily small values will
#'   add unnecessary computational time)
#'
#' @return Suggested per-group sample size for planned study
#' Publication bias and uncertainty- adjusted prior study noncentrality parameter
#'
#' @export
#' @import stats
#'
#' @examples
#' ss.power.it(t.observed=3, n=20, alpha.prior=.05, alpha.planned=.05,
#' assurance=.80, power=.80, step=.001)
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

ss.power.it <- function(t.observed, n, N, alpha.prior=.05, alpha.planned=.05,
                        assurance=.80, power=.80, step=.001)
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

  if(!missing(N))
  {
    if(!is.null(N))
    {
      if(N <= 2) stop("Your total sample size is too small")

      if(!missing(n)) stop("Because you specified 'N' you should not specify 'n'.")
      if(missing(n))
      {
        n.ru <- ceiling(N/2)
        N.ru <- 2*n.ru
        DF.ru <- N.ru-2

        n.rd <- floor(N/2)
        N.rd <- 2*n.rd
        DF.rd <- N.rd-2
      }
    }
  }

  if(missing(N))
  {
    if(!(length(n) %in% c(1, 2))) stop("The value of 'n' should be a vector of length two or a single value (for equal group sample sizes)")

    if(length(n)==2)
    {
      n.1 <- n[1]
      n.2 <- n[2]
      n.ru <- ceiling(2/((1/n.1)+(1/n.2)))
      N.ru <- 2*n.ru
      DF.ru <- N.ru-2
      n.rd <-  floor(2/((1/n.1)+(1/n.2)))
      N.rd <- 2*n.rd
      DF.rd <- N.rd-2
    }

    if(length(n)==1)
    {
      n.ru <- n
      N.ru <- 2*n.ru
      DF.ru <- N.ru-2
      n.rd <- n
      N.rd <- 2*n.rd
      DF.rd <- N.rd-2
    }
  }

  NCP <- seq(from=0, to=100, by=step)

  ## Rounding up

  d.density.ru <- dt(t.observed, df=DF.ru, ncp=NCP)

  value.critical.ru <- qt(1-alpha.prior/2, df=DF.ru)

  if(t.observed <= value.critical.ru) stop("Your observed t statistic is nonsignificant based on your specfied alpha of the prior study. Please increase 'alpha.prior' so 't.observed' exceeds the critical value")

  area.above.critical.value.ru <- 1 - pt(value.critical.ru, df=DF.ru, ncp=NCP)
  area.other.tail.ru <- pt(-1*value.critical.ru, df=DF.ru, ncp=NCP)
  power.values.ru <- area.above.critical.value.ru + area.other.tail.ru
  area.above.t.ru <- 1 - pt(t.observed, df=DF.ru, ncp = NCP)
  area.above.t.opp.ru <- pt(-1*t.observed, df = DF.ru, ncp = NCP)
  area.area.between.ru <- (area.above.critical.value.ru - area.above.t.ru) + (area.other.tail.ru - area.above.t.opp.ru)

  TM.ru <- area.area.between.ru/power.values.ru
  TM.Percentile.ru <- min(NCP[which(abs(TM.ru-assurance)==min(abs(TM.ru-assurance)))])

  if(TM.Percentile.ru==0) stop("The corrected noncentrality parameter is zero. Please either choose a lower value of assurance and/or a higher value of alpha for the prior study (e.g. accounting for less publication bias)")

  if(TM.Percentile.ru > 0)
  {

    nrep <- 2
    denom.df <- (2*nrep)-2
    diff.ru <- -1
    while (diff.ru < 0)
    {
      criticalT <- qt(1-alpha.planned/2, df = denom.df)
      powers1.ru <- 1 - pt(criticalT, df = denom.df, ncp = sqrt(nrep/n.ru)*TM.Percentile.ru)
      powers2.ru <- pt(-1*criticalT, df=denom.df, ncp = sqrt(nrep/n.ru)*TM.Percentile.ru)
      powers.ru <- powers1.ru + powers2.ru
      diff.ru <- powers.ru - power
      nrep <- nrep + 1
      denom.df = (2*nrep)-2
    }
  }
  repn.ru <- nrep-1
  ##
  ## Rounding up

  d.density.rd <- dt(t.observed, df=DF.rd, ncp=NCP)

  value.critical.rd <- qt(1-alpha.prior/2, df=DF.rd)

  if(t.observed <= value.critical.rd) stop("Your observed t statistic is nonsignificant based on your specfied alpha of the prior study. Please increase 'alpha.prior' so 't.observed' exceeds the critical value")

  area.above.critical.value.rd <- 1 - pt(value.critical.rd, df=DF.rd, ncp=NCP)
  area.other.tail.rd <- pt(-1*value.critical.rd, df=DF.rd, ncp=NCP)
  power.values.rd <- area.above.critical.value.rd + area.other.tail.rd
  area.above.t.rd <- 1 - pt(t.observed, df=DF.rd, ncp = NCP)
  area.above.t.opp.rd <- pt(-1*t.observed, df = DF.rd, ncp = NCP)
  area.area.between.rd <- (area.above.critical.value.rd - area.above.t.rd) + (area.other.tail.rd - area.above.t.opp.rd)

  TM.rd <- area.area.between.rd/power.values.rd
  TM.Percentile.rd <- min(NCP[which(abs(TM.rd-assurance)==min(abs(TM.rd-assurance)))])

  if(TM.Percentile.rd==0) stop("The corrected noncentrality parameter is zero. Please either choose a lower value of assurance and/or a higher value of alpha for the prior study (e.g. accounting for less publication bias)")

  if(TM.Percentile.rd > 0)
  {

    nrep <- 2
    denom.df <- (2*nrep)-2
    diff.rd <- -1
    while (diff.rd < 0)
    {
      criticalT <- qt(1-alpha.planned/2, df = denom.df)
      powers1.rd <- 1 - pt(criticalT, df = denom.df, ncp = sqrt(nrep/n.rd)*TM.Percentile.rd)
      powers2.rd <- pt(-1*criticalT, df=denom.df, ncp = sqrt(nrep/n.rd)*TM.Percentile.rd)
      powers.rd <- powers1.rd + powers2.rd
      diff.rd <- powers.rd - power
      nrep <- nrep + 1
      denom.df = (2*nrep)-2
    }
  }
  repn.rd <- nrep-1

  output.n <- max(repn.ru, repn.rd)
  return(list(output.n, min(TM.Percentile.ru, TM.Percentile.rd)))
}
