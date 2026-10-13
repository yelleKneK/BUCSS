# Necessary Sample Size to Reach Desired Power for a Joint Test of Multiple Predictors in a Multiple Regression Using a Publication Bias and Uncertainty Correction Procedure

`ss_buc_reg_joint` returns the necessary total sample size to achieve a
desired level of statistical power for a joint test of multiple
predictors in a planned study using multiple regression, based on
information obtained from a previous study. The effect from the previous
study can be corrected for publication bias and/or uncertainty to
provide a sample size that will achieve more accurate statistical power
for a planned study, when compared to approaches that use a sample
effect size at face value or rely on sample size only. The bias and
uncertainty adjusted previous study noncentrality parameter is also
returned, together with the planning effect size, the effect size the
plan is built on, as Cohen's \\f^2\\.

## Usage

``` r
ss_buc_reg_joint(
  F_observed,
  N,
  p,
  p_joint,
  alpha_prior = 0.05,
  alpha_planned = 0.05,
  assurance = 0.8,
  desired_power = 0.8,
  random_predictors = FALSE
)
```

## Arguments

- F_observed:

  Observed *F* value from a previous study used to plan sample size for
  a planned study.

- N:

  Total sample size of the previous study.

- p:

  Number of predictors; be sure to include any product terms or
  polynomials that are in the model.

- p_joint:

  Number of predictors tested jointly for significance.

- alpha_prior:

  The significance level \\\alpha\\ a prior result had to reach to be
  published (more generally, to reach you): the correction assumes you
  see only results with *p* below `alpha_prior`. Enter a number above 0
  and at most 1. A value of 1, which may also be written
  `"No Publication Filter"`, means no publication filter: every *p*
  value is at most 1, so any result would have been published whatever
  its *p* value, there is no publication bias to correct, and only the
  correction for uncertainty (set by `assurance`) remains, as for a
  pilot study. Anderson, Kelley, and Maxwell (2017, p. 1552) note that
  this level, their \\\alpha_P\\, can be set equal to 1 when publication
  bias is thought not to apply. The string is read without regard to
  case, and the result echoes it as an `alpha_prior` of 1.

- alpha_planned:

  Alpha level \\\alpha\\ assumed for the planned study.

- assurance:

  Desired level of assurance: the long run proportion of planned studies
  whose statistical power at the true effect reaches or surpasses
  `desired_power`, when each is planned from a prior study published
  under the rule `alpha_prior` describes. The adjusted noncentrality
  parameter is the lower limit of a one-sided confidence interval at
  this level, computed from the publication-truncated distribution of
  the prior statistic, so `assurance = .5` plans from a median unbiased
  estimate (correcting for publication bias only), a larger value adds
  the correction for uncertainty, and a value below .5 is not
  recommended. The term follows Kelley and Maxwell (2003), who used
  assurance in sample size planning for the probability that a
  confidence interval will be no wider than desired. With
  `alpha_prior = 1` and `assurance = .80` the plan corresponds to
  safeguard power, which plans for the lower limit of a two-sided 60
  percent confidence interval around the prior effect size rather than
  for the effect size itself (see
  [`BUCSS-package`](https://yelleKneK.github.io/BUCSS/reference/BUCSS-package.md)
  for how, and for where the two differ). Enter it as a proportion in
  (0, 1) or as a percentage greater than 1 (for example, 80 is read as
  .80); a value of exactly 1 is read as 1 percent. A percentage is
  echoed in the result's planning-input rows as the coerced proportion
  (80 is echoed as .80); the same applies to `desired_power`.

- desired_power:

  Desired level of statistical power for the planned study: the
  probability that its test rejects the null hypothesis when the effect
  is present (Maxwell, Delaney, & Kelley, 2027). It must exceed
  `alpha_planned`: a test's power is never below its Type I error rate,
  so a value at or below `alpha_planned` is not allowable, and the
  function stops rather than return a plan. Enter it as a proportion in
  (0, 1) or as a percentage greater than 1; a value of exactly 1 is read
  as 1 percent (.01), which is then refused unless `alpha_planned` is
  below .01.

- random_predictors:

  `TRUE` if the predictors are sampled along with the outcome, as in
  most observational research, and `FALSE` if their values are fixed by
  design, as in an experiment that assigns them. The default is `FALSE`,
  the frame of BUCSS 1.x; see Details. With random predictors (`TRUE`)
  the planner uses the exact distribution for sampled predictors that
  the correlation planner also uses: with one predictor, tested
  (`p = 1`, `p_joint = 1`), it returns exactly what the correlation
  planner,
  [`ss_buc_correlation`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_correlation.md),
  returns. The choice is echoed in the result as a `random_predictors`
  row valued 1 (`TRUE`) or 0 (`FALSE`) and named in the printed footer.

## Value

An object of class `bucss_plan`: a `data.frame` with a character `term`
column and a numeric `value` column, in the same shape as the sibling
package `DMAR`'s planners. The design results come first: the necessary
total sample size for the planned study, named for its unit
(`necessary_n_per_group`, `necessary_n_per_cell`, or `necessary_N`); for
the per-group and per-cell designs, `total_N`, the implied total sample
size; `actual_power`, the power of the plan returned, that is, of the
returned sample size at the returned noncentrality parameter (see
Details); `ncp_adjusted`, the publication bias and uncertainty adjusted
prior study noncentrality parameter, on the scale of the *F* statistic
at the prior study's sample size; and the planning effect size,
`f2_planning`, Cohen's \\f^2\\ for the set of predictors tested,
\\\lambda/N\\, where \\\lambda\\ is `ncp_adjusted` and *N* the prior
sample size; the planned study's power is computed from it at the
planned sample size, in the frame `random_predictors` names (with
sampled predictors it is the population \\f^2\\ of that frame, and the
power statement that follows is in that frame).

A planning effect size (the row whose name ends in `_planning`) is an
adjusted effect size: the value, in the design's usual effect size
metric, used to plan the sample size given the situation at hand (the
prior study's result, the publication bias assumed through
`alpha_prior`, and the uncertainty allowed for through `assurance`). It
is `ncp_adjusted` restated in that metric and is what the plan is built
on, not an estimate of the population effect to report: above an
assurance of .5 it is deliberately conservative (a lower confidence
limit), and at .5 it is the median unbiased estimate. Because it
restates `ncp_adjusted`, `actual_power` is the ordinary power of a study
of the returned sample size when the population effect is the planning
effect size.

Rows echoing the planning inputs follow the design results, so the
assumptions travel with the result. The design, the unit the sample size
is counted in, the largest assurance this prior can support (the
assurance ceiling), and the planned test's degrees of freedom
(`df_effect`, and `df_error` where the test has a denominator) travel on
attributes. `print.bucss_plan` shows the design, for the regression
planners the predictor frame, the unit, and the ceiling beneath the
table; [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) give the
compact and wide views, the latter with every attribute as a column.

## Details

Researchers often use the sample effect size from a prior study as an
estimate of the likely size of an expected future effect in sample size
planning. However, sample effect size estimates should not usually be
used at face value to plan sample size, due to both publication bias and
uncertainty.

The approach implemented in `ss_buc_reg_joint` uses the observed *F*
value and sample size from a previous study to correct the noncentrality
parameter associated with the effect of interest for publication bias
and/or uncertainty. This new estimated noncentrality parameter is then
used to calculate the necessary total sample size to achieve the desired
level of power in the planned study.

**Fixed or random predictors.** Whether the predictors' values are fixed
by design, as in an experiment in which the researcher sets them, or
sampled along with the outcome, as in most observational research,
changes neither the estimates nor the distribution of the test statistic
when the tested coefficients are zero, but it does change the
statistic's distribution, and so the power, when they are not all zero
(Sampson, 1974; Gatsonis & Sampson, 1989; Kelley & Maxwell, 2003). With
`random_predictors = FALSE`, the default, the prior study's statistic
and the planned study's power come from the noncentral *F* distribution,
which is exact when the predictor values are fixed; this is the frame of
BUCSS 1.x and of Anderson (2021). With `random_predictors = TRUE` the
noncentrality parameter is itself random: it is proportional to the
sampled spread of the tested predictors once the other predictors are
accounted for, a chi square variable on *N* - 1 - (*p* - `p_joint`)
degrees of freedom, so the statistic follows a chi square mixture of
noncentral *F* distributions, which is more dispersed. (When every
predictor is tested, the mixture is the general sampling distribution of
the multiple correlation coefficient that Fisher, 1928, derived; for a
subset it is that distribution for the tested predictors once the other
predictors are accounted for.) The mixture has an exact closed form (the
Poisson weights of the noncentral *F* become negative binomial weights),
used for both the correction and the planned study's power. It is exact
when the tested predictors are multivariate normal around their linear
regression on the other predictors, which may themselves be fixed by
design (dummy codes for assigned groups, say). In observational research
in the behavioral sciences the predictors are typically random (Kelley &
Maxwell, 2003), and `random_predictors = TRUE` is then the frame to plan
in.

Treating sampled predictors as fixed overstates the planned study's
power, by about 1/*n* to 1.4/*n* for a plan of *n* participants when
every predictor is tested, and more the more predictors are left
untested: 2/*n* for two of four and 3.8/*n* for two of eight (at *n* =
100), so the assurance the plans deliver can fall short of the requested
`assurance`. The delivered assurance is what `assurance` promises: the
proportion of plans, each made from its own published prior study, whose
power at the true effect reaches `desired_power`, counting a refusal to
plan as meeting the target. In simulations with multivariate normal
predictors, an `assurance` of .80, and prior studies of 40 to 100
participants with power .3 to .8 at the true effect, testing two of four
or five predictors or three of six, the fixed frame delivered about .76
to .81, the lower values from smaller, better powered prior studies;
`random_predictors = TRUE` delivered .79 to .81 in the same
configurations.

Which frame plans the larger study depends on the untested predictors
and on how far the plan outgrows the prior study, because the fixed
frame scales the prior study's noncentrality parameter by *N* while the
tested predictors' sampled spread has *N* - 1 - (*p* - `p_joint`)
degrees of freedom. With every predictor tested the random-predictor
plan is usually a few participants larger (and at most a few percent
smaller). With predictors left untested the two plans are within a few
participants, the random-predictor plan usually the larger, when the
plan is no more than about twice the prior *N*; beyond that the
random-predictor plan is usually smaller, and the more the plan outgrows
the prior study, the closer the ratio of the two plans comes to (*N* -
1 - (*p* - `p_joint`))/*N*, so the random-predictor plan is smaller by
at most about (1 + *p* - `p_joint`)/*N*: 518 rather than 600
participants for *F* = 6 from a prior study of 30 testing two of six
predictors, where the limit is about 17 percent, and, with two of eight
predictors tested, a limit of 35 percent from a prior study of 20 and 7
percent from 100. The largest assurance a prior result can support is
the same in both frames, because at a zero effect both are the central
*F* distribution. With `p_joint = p` this function returns what
[`ss_buc_R2`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_R2.md)
returns, and with `p_joint = 1` what
[`ss_buc_reg_coef`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_coef.md)
returns for the square root of `F_observed`, in either frame. So with
one predictor, tested (`p = 1`, `p_joint = 1`), and
`random_predictors = TRUE`, it returns exactly what
[`ss_buc_correlation`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_correlation.md)
returns for the same prior *t* (the square root of `F_observed`) and *N*
(the same sample size, power, noncentrality parameter, and ceiling, with
the planning effect size stated as \\f^2\\ rather than as the
correlation), for any prior *N* of at least 4 (the correlation planner's
minimum).

The approach uses a likelihood function of a truncated noncentral *F*
distribution, where the truncation occurs due to small effect sizes
being unobserved due to publication bias. The numerator of the
likelihood function is the density of a noncentral *F* distribution
(with `random_predictors = TRUE`, of the mixture described above). The
denominator is the power of the test, which serves to truncate the
distribution. Thus, the ratio of the numerator and the denominator is a
truncated noncentral *F* distribution. (See Taylor & Muller, 1996,
Equation 2.1; Anderson, Kelley, & Maxwell, 2017, Equation 3; and
Anderson & Maxwell, 2017, for more details.)

Assurance is the proportion of times that power will be at or above the
desired level, if the experiment were to be reproduced many times. For
example, `assurance = .5` means that power will be above the desired
level half of the time, but below the desired level the other half of
the time. Selecting `assurance = .5` (selecting the noncentrality
parameter at the 50th percentile of the likelihood distribution) results
in a median unbiased estimate of the population noncentrality parameter
and does not correct for uncertainty. In order to correct for
uncertainty, `assurance > .5` can be selected, which corresponds to
selecting the noncentrality parameter associated with the (1 -
assurance) quantile of the likelihood distribution.

`alpha_prior` is the significance level a prior result had to reach to
be published: the correction assumes that only results with *p* values
below `alpha_prior` reach the literature, and so reach you. Different
\\\alpha\\ levels represent differing amounts of publication bias. For
example, setting `alpha_prior = .20` assumes that studies with *p*
values less than .20 are published, whereas those with larger *p* values
are not, which reflects less severe publication bias than the default of
.05. If the previous study of interest has not been subjected to
publication bias (e.g., a pilot study), set `alpha_prior = 1`, which may
also be written `alpha_prior = "No Publication Filter"`; Anderson,
Kelley, and Maxwell (2017, p. 1552) note that it can be set equal to 1
when publication bias is thought not to apply. Because no *p* value
exceeds 1, this means no publication filter at all: any result would
have been published whatever its *p* value, so there is no publication
bias to correct, and only the correction for uncertainty, set by
`assurance`, remains.

In some cases, the corrected noncentrality parameter for a given level
of assurance will be estimated to be zero. This is an indication that,
at the desired level of assurance, the previous study's effect cannot be
accurately estimated due to high levels of uncertainty and bias. When
this happens, subsequent sample size planning is not possible with the
chosen specifications, and the function stops with an informative error
rather than returning a result. The error reports the largest assurance
the prior result can support, which is \\1 - p/\alpha\\, where *p* is
the prior study's *p* value and \\\alpha\\ is `alpha_prior`; at any
higher assurance the corrected noncentrality parameter is driven to
zero. Two remedies are available. First, users can select a lower value
of assurance (e.g., .8 instead of .95), at or below the reported
ceiling. Second, users can reduce the influence of publication bias by
setting `alpha_prior` at a value greater than .05, which raises the
ceiling. Note that `alpha_prior` is not a property of the prior study,
which is fixed, but the analyst's assumption about the publication
threshold its literature faced. When the ceiling is at or below the
recommended floor of .5, lowering assurance cannot help and
`alpha_prior` is the only remaining lever. It is possible to correct for
uncertainty only by setting `alpha_prior = 1` and choosing the desired
level of assurance. We encourage users to make the adjustments as
minimal as possible.

The returned `actual_power` is the statistical power of the plan this
function reports: the power of a study of the returned sample size when
the effect is the returned bias and uncertainty adjusted noncentrality
parameter. Because a sample size must be a whole number, `actual_power`
is never below the `desired_power` that was requested and is usually a
little above it.

The search for the planned sample size starts at *p* + 3, the smallest
total that leaves the planned study two error degrees of freedom (as in
BUCSS 1.x), so no plan is smaller than that, even when the prior effect
is so large that a smaller study would already reach the desired power.
A plan at this floor signals a very large prior effect, not a
recommended design.

## References

Anderson, S. F., & Kelley, K. (2024). Sample size planning for
replication studies: The devil is in the design. *Psychological Methods,
29,* 844–867.
[doi:10.1037/met0000520](https://doi.org/10.1037/met0000520)

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28,* 1547–1562.
[doi:10.1177/0956797617723724](https://doi.org/10.1177/0956797617723724)

Anderson, S. F., & Maxwell, S. E. (2017). Addressing the "replication
crisis": Using original studies to design replication studies with
appropriate statistical power. *Multivariate Behavioral Research, 52,*
305–324.
[doi:10.1080/00273171.2017.1289361](https://doi.org/10.1080/00273171.2017.1289361)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8,* 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Taylor, D. J., & Muller, K. E. (1996). Bias in linear model power and
sample size calculation due to estimating noncentrality. *Communications
in Statistics: Theory and Methods, 25,* 1595–1610.
[doi:10.1080/03610929608831787](https://doi.org/10.1080/03610929608831787)

Anderson, S. F. (2021). Using prior information to plan appropriately
powered regression studies: A tutorial using BUCSS. *Psychological
Methods, 26,* 513–526.
[doi:10.1037/met0000366](https://doi.org/10.1037/met0000366)

Fisher, R. A. (1928). The general sampling distribution of the multiple
correlation coefficient. *Proceedings of the Royal Society of London.
Series A, Containing Papers of a Mathematical and Physical Character,
121,* 654–673.
[doi:10.1098/rspa.1928.0224](https://doi.org/10.1098/rspa.1928.0224)

Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact power
and sample size calculations. *Psychological Bulletin, 106,* 516–524.
[doi:10.1037/0033-2909.106.3.516](https://doi.org/10.1037/0033-2909.106.3.516)

Sampson, A. R. (1974). A tale of two regressions. *Journal of the
American Statistical Association, 69,* 682–689.
[doi:10.1080/01621459.1974.10480189](https://doi.org/10.1080/01621459.1974.10480189)

## Author

Ken Kelley (<kkelley@nd.edu>) and Samantha F. Anderson
(<samantha.f.anderson@asu.edu>)

## Examples

``` r
result <- ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
  alpha_prior = .05, alpha_planned = .05, assurance = .80, desired_power = .80)
result
#>  term              value  
#>  necessary_N       3960   
#>  actual_power      0.8    
#>  ncp_adjusted      0.365  
#>  f2_planning       0.00244
#>  F_observed        5      
#>  N                 150    
#>  p                 4      
#>  p_joint           2      
#>  random_predictors 0      
#>  alpha_prior       0.05   
#>  alpha_planned     0.05   
#>  assurance         0.8    
#>  desired_power     0.8    
#> 
#> Design: Multiple regression: joint test of predictors
#> Predictors: fixed by design
#> Sample size unit: total
#> Largest supportable assurance: .84

# Predictors sampled along with the outcome, as in most observational
# research:
ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
  random_predictors = TRUE)
#>  term              value  
#>  necessary_N       3898   
#>  actual_power      0.8    
#>  ncp_adjusted      0.372  
#>  f2_planning       0.00248
#>  F_observed        5      
#>  N                 150    
#>  p                 4      
#>  p_joint           2      
#>  random_predictors 1      
#>  alpha_prior       0.05   
#>  alpha_planned     0.05   
#>  assurance         0.8    
#>  desired_power     0.8    
#> 
#> Design: Multiple regression: joint test of predictors
#> Predictors: random (sampled along with the outcome)
#> Sample size unit: total
#> Largest supportable assurance: .84

# Requesting more assurance than the prior result can support stops with an
# informative error naming the largest workable assurance (here near .84):
try(ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
  assurance = .95))
#> Error : Sample size planning is not possible with these settings. After correcting the prior study's effect for publication bias and uncertainty at the requested level of assurance, the corrected noncentrality parameter is zero, which means there is not enough information in the prior result to plan the sample size as desired. With this prior result and 'alpha_prior', the largest 'assurance' that still yields a usable plan is about .84; re-run with 'assurance' at or below that value, or raise 'alpha_prior' to lift the ceiling. Re-running at a lower 'assurance' is not free: the plan it returns carries that lower assurance, not the one you first asked for, so the long run proportion of replications reaching your target power falls accordingly. 'alpha_prior' does not change the prior study, which is what it is; it is your assumption about the publication threshold the study's literature faced, so raising it toward 1 assumes less publication bias to undo (for example, 'alpha_prior = .10' assumes the literature's journals would have published a result significant at the .10 level, a more lenient policy than .05, and 'alpha_prior = 1' assumes no publication bias at all). See Anderson, Kelley, and Maxwell (2017) and Anderson and Kelley (2024).
```
