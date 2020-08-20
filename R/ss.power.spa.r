# Split plot ANOVA.

#' Necessary sample size to reach desired power for two-factor split-plot
#' (mixed) ANOVA using an uncertainty and publication bias correction procedure
#'
#' @description \code{ss.power.spa} returns the necessary per-group sample size
#'   to achieve a desired level of statistical power for a planned study testing
#'   an omnibus effect using a two-factor split-plot (mixed) ANOVA, based on
#'   information obtained from a previous study. The effect from the previous
#'   study can be corrected for publication bias and/or uncertainty to provide a
#'   sample size that will achieve more accurate statistical power for a planned
#'   study, when compared to approaches that use a sample effect size at face
#'   value or rely on sample size only. The bias and uncertainty adjusted previous
#'   study noncentrality parameter is also returned, which can be transformed to
#'   various effect size metrics.
#'
#' @details Researchers often use the sample effect size from a prior study as
#'   an estimate of the likely size of an expected future effect in sample size
#'   planning. However, sample effect size estimates should not usually be used
#'   at face value to plan sample size, due to both publication bias and
#'   uncertainty.
#'
#'   The approach implemented in \code{ss.power.spa} uses the observed
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
#'   2.1. and Anderson & Maxwell, 2017, for more details.)
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
#'   \code{ss.power.spa} assumes that the planned study will have equal n.
#'   Unequal n in the previous study is handled in the following way for split
#'   plot designs. If the user enters an N not equally divisible by the number
#'   of between-subjects cells, the function calculates n by dividing N by the
#'   number of cells and both rounding up and rounding down the result,
#'   effectively assuming equal n. The suggested sample size for the planned
#'   study is calculated using both of these values of n, and the function
#'   returns the larger of these two suggestions, to be conservative. The
#'   adjusted noncentrality parameter that is output is the lower of the two
#'   possibilities, again, to be conservative. Although equal-n previous
#'   studies are preferable, this approach will work well as long as the cell
#'   sizes are only slightly discrepant.
#'
#'   \code{ss.power.spa} assumes sphericity for the within-subjects effects.
#'
#' @param F.observed Observed F-value from a previous study used to plan sample
#'   size for a planned study
#' @param N Total sample size of the previous study
#' @param levels.between Number of levels for the between-subjects factor
#' @param levels.within Number of levels for the within-subjects factor
#' @param effect Effect most of interest to the planned study: between main
#'   effect (between), within main effect (within), interaction
#' @param alpha.prior Alpha-level \eqn{\alpha} for the previous study or the
#'   assumed statistical significance necessary for publishing in the field; to
#'   assume no publication bias, a value of 1 can be entered 
#' @param alpha.planned Alpha level assumed for the planned study
#' @param assurance Desired level of assurance, or the long run proportion of
#'   times that the planned study power will reach or surpass desired level
#'   (assurance > .5 corrects for uncertainty; assurance < .5 not recommended)
#' @param power Desired level of statistical power for the planned study
#' @param step Value used in the iterative scheme to determine the noncentral
#'   parameters necessary for sample size planning (0 < step < .01) (users
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
#' ss.power.spa(F.observed=5, N=60, levels.between=2, levels.within=3,
#' effect="within", alpha.prior=.05, alpha.planned=.05, assurance=.80,
#' power=.80, step=.001)
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


ss.power.spa <- function(F.observed, N, levels.between, levels.within, effect=c("between", "within", "interaction"), alpha.prior=.05, alpha.planned=.05, assurance=.80, power=.80, step=.001)
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

if(missing(levels.within)) stop("You must specify the number of levels of the within factor. If there is no within factor, use the between subjects appraoch.")
if(missing(levels.between)) stop("You must specify the number of levels of the between factor. If there is no between factor, use the within subjects appraoch.")


NCP <- seq(from=0, to=100, by=step) # sequence of possible values for the noncentral parameter.

### ROUND DOWN

if(effect=="between")
{
n.rd <- floor(N/levels.between) # To ensure that the between sample size is appropriate given specifications.
N.rd <- n.rd*levels.between

df.numerator <- levels.between - 1
df.denominator.rd <- N.rd-levels.between
}

if(effect=="within")
{
n.rd <- floor(N/levels.between) # To ensure that the per-cell sample size is equal. Rounds down (floor).
N.rd <- n.rd*levels.between

df.numerator <- levels.within - 1
df.denominator.rd <- (N.rd-levels.between)*(levels.within - 1)
}

if(effect=="interaction")
{
n.rd <- floor(N/levels.between) # To ensure that the per-cell sample size is equal. Rounds down (floor).
N.rd <- n.rd*levels.between

df.numerator <- (levels.between - 1)*(levels.within - 1)
df.denominator.rd <- (N.rd-levels.between)*(levels.within - 1)
}

f.density.rd <- df(F.observed, df1=df.numerator, df2=df.denominator.rd, ncp=NCP) # density of F using F observed
critF.rd <- qf(1-alpha.prior, df1=df.numerator, df2=df.denominator.rd)

if(F.observed <= critF.rd) stop("Your observed F statistic is nonsignificant based on your specfied alpha of the prior study. Please increase 'alpha.prior' so 'F.observed' exceeds the critical value")

power.values.rd <- 1 - pf(critF.rd, df1=df.numerator, df2=df.denominator.rd, ncp = NCP) # area above critical F
area.above.F.rd <- 1 - pf(F.observed, df1=df.numerator, df2=df.denominator.rd, ncp = NCP) # area above observed F
area.area.between.rd <- power.values.rd - area.above.F.rd

TM.rd <- area.area.between.rd/power.values.rd
TM.Percentile.rd <- min(NCP[which(abs(TM.rd-assurance)==min(abs(TM.rd-assurance)))])

if(TM.Percentile.rd==0) stop("The corrected noncentrality parameter is zero. Please either choose a lower value of assurance and/or a higher value of alpha for the prior study (e.g. accounting for less publication bias)")

if (TM.Percentile.rd > 0)
{
nrep <- 2

if(effect=="between")
{
denom.df <- nrep*levels.between-levels.between
}

if(effect=="within")
{
denom.df <- levels.between*(nrep-1)*(levels.within - 1)
}

if(effect=="interaction")
{
denom.df <- levels.between*(nrep-1)*(levels.within - 1)
}


diff.rd <- -1
while (diff.rd < 0 )
{
criticalF <- qf(1-alpha.planned, df1 = df.numerator, df2 = denom.df)
powers.rd <- 1 - pf(criticalF, df1 = df.numerator, df2 = denom.df, ncp = (nrep/n.rd)*TM.Percentile.rd)
diff.rd <- powers.rd - power
nrep <- nrep + 1

if(effect=="between")
{
denom.df <- (nrep*levels.between)-levels.between
}

if(effect=="within")
{
denom.df <- levels.between*(nrep-1)*(levels.within - 1)
}

if(effect=="interaction")
{
denom.df <- levels.between*(nrep-1)*(levels.within - 1)
}

}
}
repn.rd <- nrep-1

### ROUND UP

if(effect=="between")
{
  n.ru <- ceiling(N/levels.between) # To ensure that the between sample size is appropriate given specifications.
  N.ru <- n.ru*levels.between

  df.numerator <- levels.between - 1
  df.denominator.ru <- N.ru-levels.between
}

if(effect=="within")
{
  n.ru <- ceiling(N/levels.between) # To ensure that the per-cell sample size is equal. Rounds down (floor).
  N.ru <- n.ru*levels.between

  df.numerator <- levels.within - 1
  df.denominator.ru <- (N.ru-levels.between)*(levels.within - 1)
}

if(effect=="interaction")
{
  n.ru <- ceiling(N/levels.between) # To ensure that the per-cell sample size is equal. Rounds down (floor).
  N.ru <- n.ru*levels.between

  df.numerator <- (levels.between - 1)*(levels.within - 1)
  df.denominator.ru <- (N.ru-levels.between)*(levels.within - 1)
}

f.density.ru <- df(F.observed, df1=df.numerator, df2=df.denominator.ru, ncp=NCP) # density of F using F observed
critF.ru <- qf(1-alpha.prior, df1=df.numerator, df2=df.denominator.ru)

if(F.observed <= critF.ru) stop("Your observed F statistic is nonsignificant based on your specfied alpha of the prior study. Please increase 'alpha.prior' so 'F.observed' exceeds the critical value")

power.values.ru <- 1 - pf(critF.ru, df1=df.numerator, df2=df.denominator.ru, ncp = NCP) # area above critical F
area.above.F.ru <- 1 - pf(F.observed, df1=df.numerator, df2=df.denominator.ru, ncp = NCP) # area above observed F
area.area.between.ru <- power.values.ru - area.above.F.ru

TM.ru <- area.area.between.ru/power.values.ru
TM.Percentile.ru <- min(NCP[which(abs(TM.ru-assurance)==min(abs(TM.ru-assurance)))])

if(TM.Percentile.ru==0) stop("The corrected noncentrality parameter is zero. Please either choose a lower value of assurance and/or a higher value of alpha for the prior study (e.g. accounting for less publication bias)")

if (TM.Percentile.ru > 0)
{
  nrep <- 2

  if(effect=="between")
  {
    denom.df <- nrep*levels.between-levels.between
  }

  if(effect=="within")
  {
    denom.df <- levels.between*(nrep-1)*(levels.within - 1)
  }

  if(effect=="interaction")
  {
    denom.df <- levels.between*(nrep-1)*(levels.within - 1)
  }


  diff.ru <- -1
  while (diff.ru < 0 )
  {
    criticalF <- qf(1-alpha.planned, df1 = df.numerator, df2 = denom.df)
    powers.ru <- 1 - pf(criticalF, df1 = df.numerator, df2 = denom.df, ncp = (nrep/n.ru)*TM.Percentile.ru)
    diff.ru <- powers.ru - power
    nrep <- nrep + 1

    if(effect=="between")
    {
      denom.df <- (nrep*levels.between)-levels.between
    }

    if(effect=="within")
    {
      denom.df <- levels.between*(nrep-1)*(levels.within - 1)
    }

    if(effect=="interaction")
    {
      denom.df <- levels.between*(nrep-1)*(levels.within - 1)
    }

  }
}
repn.ru <- nrep-1

output.n <- max(repn.rd, repn.ru)
return(list(output.n, min(TM.Percentile.rd, TM.Percentile.ru)))
}
