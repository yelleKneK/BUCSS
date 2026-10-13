# Necessary Sample Size to Reach Desired Power for a Two-Group Multivariate Comparison Using a Publication Bias and Uncertainty Correction Procedure

`ss_buc_manova` returns the necessary per-group sample size to achieve a
desired level of statistical power for a planned study comparing two
groups on several outcome variables at once (Hotelling's \\T^2\\,
equivalently a one-way multivariate analysis of variance with two
groups), based on information obtained from a previous study. The effect
from the previous study can be corrected for publication bias and/or
uncertainty to provide a sample size that will achieve more accurate
statistical power for a planned study, when compared to approaches that
use a sample effect size at face value or rely on sample size only. The
bias and uncertainty adjusted previous study noncentrality parameter is
also returned, together with the planning effect size, the effect size
the plan is built on, as a squared Mahalanobis distance.

## Usage

``` r
ss_buc_manova(
  F_observed,
  N,
  p_variables,
  T2_observed,
  alpha_prior = 0.05,
  alpha_planned = 0.05,
  assurance = 0.8,
  desired_power = 0.8
)
```

## Arguments

- F_observed:

  Observed *F* value for the multivariate comparison from a previous
  study. Supply either `F_observed` or `T2_observed`, not both.

- N:

  Total sample size of the previous study (both groups combined).

- p_variables:

  Number of outcome variables compared.

- T2_observed:

  Observed Hotelling's \\T^2\\ from the previous study, an alternative
  to `F_observed`.

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

## Value

An object of class `bucss_plan`: a `data.frame` with a character `term`
column and a numeric `value` column, in the same shape as the sibling
package `DMAR`'s planners. The design results come first: the necessary
per-group sample size for the planned study, named for its unit
(`necessary_n_per_group`, `necessary_n_per_cell`, or `necessary_N`); for
the per-group and per-cell designs, `total_N`, the implied total sample
size; `actual_power`, the power of the plan returned, that is, of the
returned sample size at the returned noncentrality parameter (see
Details); `ncp_adjusted`, the publication bias and uncertainty adjusted
prior study noncentrality parameter, on the scale of the *F* statistic,
at the prior study's sample size (for an odd prior *N*, at the balanced
sample size the planner substitutes for it; with one outcome variable it
is the square of the independent *t* planner's); and the planning effect
size, `D2_planning`, the squared Mahalanobis distance between the two
population mean vectors, \\2\lambda/n\\, where \\\lambda\\ is
`ncp_adjusted` and *n* the per-group size of the balanced prior design
it belongs to (see Details); with one outcome variable it is the square
of the `smd_planning` of
[`ss_buc_independent_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md).

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

**Scope.** The correction this package implements works on a single test
statistic with an exact noncentral *F* distribution. A multivariate
hypothesis has one only when its rank is one, which is the case for a
comparison of two groups: Hotelling's \\T^2\\ (the generalization of
Student's ratio; Hotelling, 1931) on *p* outcome variables with total
sample size *N* maps exactly to \\F = T^2 (N - p - 1) / (p (N - 2))\\ on
*p* and *N* - *p* - 1 degrees of freedom. Supply either the *F* or the
\\T^2\\. For a hypothesis of rank greater than one (three or more
groups, say) the multivariate noncentrality is a set of eigenvalues
rather than a single number, so the power of the usual tests depends on
more than one parameter and this correction, which solves for one, does
not apply; `ss_buc_manova` will not accept such a design.

The approach uses a likelihood function of a truncated noncentral *F*
distribution, where the truncation occurs due to small effect sizes
being unobserved due to publication bias. The numerator of the
likelihood function is the density of a noncentral *F* distribution. The
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

`ss_buc_manova` assumes that the planned study will have equal *n* per
group. If the prior *N* is odd, the per-group size is rounded both up
and down, the planned sample size is computed under each, and the larger
of the two is returned, with the smaller of the two adjusted
noncentrality parameters, to be conservative. The planning effect size
belongs to the same balanced design as the reported noncentrality
parameter, so an ordinary power analysis at it can return a sample size
a little smaller than the one reported here (never a larger one),
because the reported size is the larger of the two rounded designs'
requirements. For a prior study with unequal group sizes, or a total
sample size that does not divide evenly, the largest supportable
assurance is computed at the balanced design the planner plans from, so
it can be below \\1 - p/\alpha\\ at the study's own *p* value.

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

Hotelling, H. (1931). The generalization of Student's ratio. *The Annals
of Mathematical Statistics, 2,* 360–378.
[doi:10.1214/aoms/1177732979](https://doi.org/10.1214/aoms/1177732979)

## Author

Ken Kelley (<kkelley@nd.edu>) and Samantha F. Anderson
(<samantha.f.anderson@asu.edu>)

## Examples

``` r
result <- ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3,
  alpha_prior = .05, alpha_planned = .05, assurance = .80,
  desired_power = .80)
result
#>  term                  value 
#>  necessary_n_per_group 310   
#>  total_N               620   
#>  actual_power          0.801 
#>  ncp_adjusted          1.42  
#>  D2_planning           0.0709
#>  F_observed            4.5   
#>  N                     80    
#>  p_variables           3     
#>  alpha_prior           0.05  
#>  alpha_planned         0.05  
#>  assurance             0.8   
#>  desired_power         0.8   
#> 
#> Design: Two-group multivariate comparison (Hotelling's T squared)
#> Sample size unit: per group
#> Largest supportable assurance: .88

# The equivalent call from Hotelling's T squared:
ss_buc_manova(T2_observed = 4.5 * 3 * 78 / 76, N = 80, p_variables = 3)
#>  term                  value 
#>  necessary_n_per_group 310   
#>  total_N               620   
#>  actual_power          0.801 
#>  ncp_adjusted          1.42  
#>  D2_planning           0.0709
#>  T2_observed           13.9  
#>  N                     80    
#>  p_variables           3     
#>  alpha_prior           0.05  
#>  alpha_planned         0.05  
#>  assurance             0.8   
#>  desired_power         0.8   
#> 
#> Design: Two-group multivariate comparison (Hotelling's T squared)
#> Sample size unit: per group
#> Largest supportable assurance: .88

# Requesting more assurance than the prior result can support stops with an
# informative error naming the largest workable assurance (here near .88):
try(ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3,
  assurance = .95))
#> Error : Sample size planning is not possible with these settings. After correcting the prior study's effect for publication bias and uncertainty at the requested level of assurance, the corrected noncentrality parameter is zero, which means there is not enough information in the prior result to plan the sample size as desired. With this prior result and 'alpha_prior', the largest 'assurance' that still yields a usable plan is about .88; re-run with 'assurance' at or below that value, or raise 'alpha_prior' to lift the ceiling. Re-running at a lower 'assurance' is not free: the plan it returns carries that lower assurance, not the one you first asked for, so the long run proportion of replications reaching your target power falls accordingly. 'alpha_prior' does not change the prior study, which is what it is; it is your assumption about the publication threshold the study's literature faced, so raising it toward 1 assumes less publication bias to undo (for example, 'alpha_prior = .10' assumes the literature's journals would have published a result significant at the .10 level, a more lenient policy than .05, and 'alpha_prior = 1' assumes no publication bias at all). See Anderson, Kelley, and Maxwell (2017) and Anderson and Kelley (2024).
```
