# Necessary Sample Size to Reach Desired Power for a Two-Factor Split-Plot (Mixed) ANOVA Using a Publication Bias and Uncertainty Correction Procedure

`ss_buc_mixed_anova` returns the necessary per between-subjects cell
sample size to achieve a desired level of statistical power for a
planned study testing an omnibus effect using a two-factor split-plot
(mixed) ANOVA, based on information obtained from a previous study. The
effect from the previous study can be corrected for publication bias
and/or uncertainty to provide a sample size that will achieve more
accurate statistical power for a planned study, when compared to
approaches that use a sample effect size at face value or rely on sample
size only. The bias and uncertainty adjusted previous study
noncentrality parameter is also returned, together with the planning
effect size, the effect size the plan is built on, as Cohen's *f*.

## Usage

``` r
ss_buc_mixed_anova(
  F_observed,
  N,
  levels_between,
  levels_within,
  effect = c("between", "within", "interaction"),
  alpha_prior = 0.05,
  alpha_planned = 0.05,
  assurance = 0.8,
  desired_power = 0.8,
  epsilon = 1
)
```

## Arguments

- F_observed:

  Observed *F* value from a previous study used to plan sample size for
  a planned study.

- N:

  Total sample size of the previous study.

- levels_between:

  Number of levels for the between-subjects factor.

- levels_within:

  Number of levels for the within-subjects factor.

- effect:

  Effect of most interest to the planned study: between-subjects main
  effect (`between`), within-subjects main effect (`within`), or
  interaction (`interaction`). Case-insensitive shorthand is also
  accepted but not recommended in scripts: `"bs"` for `between`, `"ws"`
  for `within`, and `"BxW"` for `interaction`.

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

- epsilon:

  The sphericity parameter \\\epsilon\\ of the within-subjects
  contrasts, for the within-subjects and interaction effects, from
  1/(`levels_within` - 1) (its lower bound, for both effects) to 1 (the
  default, which assumes sphericity); with two within-subjects levels it
  takes only 1. For a prior study reported with Greenhouse-Geisser or
  Huynh-Feldt corrected degrees of freedom, it is the reported numerator
  degrees of freedom divided by the uncorrected ones. The same value
  describes the planned study, and it is echoed as the `epsilon` row of
  the result. It does not apply to the between-subjects effect, whose
  test is on subject means: there a value other than 1 is ignored with a
  warning, and the result has no `epsilon` row. See Details.

## Value

An object of class `bucss_plan`: a `data.frame` with a character `term`
column and a numeric `value` column, in the same shape as the sibling
package `DMAR`'s planners. The design results come first: the necessary
per-group sample size (per between-subjects cell; the row is
`necessary_n_per_group`) for the planned study, named for its unit
(`necessary_n_per_group`, `necessary_n_per_cell`, or `necessary_N`); for
the per-group and per-cell designs, `total_N`, the implied total sample
size; `actual_power`, the power of the plan returned, that is, of the
returned sample size at the returned noncentrality parameter (see
Details); `ncp_adjusted`, the publication bias and uncertainty adjusted
prior study noncentrality parameter, on the scale of the *F* statistic,
at the prior study's sample size (for a prior *N* that does not divide
evenly, at the balanced sample size the planner substitutes for it), and
on the sphericity scale for every `epsilon` (see Details); and the
planning effect size, `f_planning`, Cohen's *f* for the effect tested,
relative to its own error term, \\\sqrt{\lambda/(kN)}\\, where
\\\lambda\\ is `ncp_adjusted`, *N* the total of the balanced prior
design it belongs to, and *k* the error degrees of freedom per subject,
1 for the between-subjects effect and `levels_within` - 1 for the
within-subjects effect and the interaction, so \\f^2/(1 + f^2)\\ is the
population partial eta squared for the effect; like `ncp_adjusted` it is
on the sphericity scale, and for `epsilon` below 1 the power statement
that follows holds under the same Muller and Barton (1989)
approximation, with the corrected degrees of freedom and noncentrality
parameter; it is not the *f* that power functions for repeated measures
take (see Details).

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
assumptions travel with the result. The design, the effect tested, the
unit the sample size is counted in, the largest assurance this prior can
support (the assurance ceiling), and the planned test's degrees of
freedom (`df_effect` and `df_error`) travel on attributes.
`print.bucss_plan` shows the design, the effect, the unit, and the
ceiling beneath the table;
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) give the
compact and wide views, the latter with every attribute as a column.

## Details

Researchers often use the sample effect size from a prior study as an
estimate of the likely size of an expected future effect in sample size
planning. However, sample effect size estimates should not usually be
used at face value to plan sample size, due to both publication bias and
uncertainty.

The approach implemented in `ss_buc_mixed_anova` uses the observed *F*
value and sample size from a previous study to correct the noncentrality
parameter associated with the effect of interest for publication bias
and/or uncertainty. This new estimated noncentrality parameter is then
used to calculate the necessary per between-subjects cell sample size to
achieve the desired level of power in the planned study.

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

`ss_buc_mixed_anova` assumes that the planned study will have equal *n*.
Unequal *n* in the previous study is handled in the following way for
split-plot designs. If the user enters an *N* not equally divisible by
the number of between-subjects cells, the function calculates *n* by
dividing *N* by the number of cells and both rounding up and rounding
down the result, effectively assuming equal *n*. The suggested sample
size for the planned study is calculated using both of these values of
*n*, and the function returns the larger of these two suggestions, to be
conservative. The adjusted noncentrality parameter that is output is the
lower of the two possibilities, again, to be conservative. The planning
effect size belongs to the same balanced design as the reported
noncentrality parameter, so an ordinary power analysis at it can return
a sample size a little smaller than the one reported here (never a
larger one), because the reported size is the larger of the two rounded
designs' requirements. For a prior study with unequal group sizes, or a
total sample size that does not divide evenly, the largest supportable
assurance is computed at the balanced design the planner plans from, so
it can be below \\1 - p/\alpha\\ at the study's own *p* value. Although
equal-*n* previous studies are preferable, this approach will work well
as long as the cell sizes are only slightly discrepant.

**Sphericity.** The default `epsilon = 1` assumes sphericity, as the
uncorrected within-subjects test does. When a prior study reports
Greenhouse-Geisser (1959) or Huynh-Feldt (1976) corrected degrees of
freedom, supply `epsilon`, the reported numerator degrees of freedom
divided by the uncorrected ones. For example, *F*(1.6, 92.8) for the
within-subjects effect of a 2 x 3 design with 60 subjects gives
`epsilon = 1.6/2`, and the denominator agrees: 92.8/116 is also .8. The
ratio of the two reported degrees of freedom, 92.8/1.6 = 58, is the
number of subjects the analysis used less the number of groups (for the
interaction, that difference divided by one less than the number of
groups); when it points to a different number of subjects than the study
reports, that number is the *N* to supply. An implied `epsilon` above 1
or below its lower bound (see Arguments) means the reported degrees of
freedom are not a Greenhouse-Geisser or Huynh-Feldt correction of this
effect; a Huynh-Feldt estimate above 1 is reported as the uncorrected
degrees of freedom, which is `epsilon = 1`. `ss_buc_mixed_anova` then
applies the approximation of Muller and Barton (1989) to both studies:
the prior *F* is treated as noncentral *F* with both degrees of freedom
and the noncentrality parameter multiplied by `epsilon`, the publication
threshold is the critical value at the corrected degrees of freedom, and
the planned study's power is computed the same way. It applies to the
within-subjects effect and the interaction, which share the
within-subjects contrasts and so share one \\\epsilon\\; the
between-subjects test is on subject means and does not involve
sphericity. The returned `ncp_adjusted` stays on the sphericity scale,
so it is comparable across values of `epsilon`, and the planned test's
corrected degrees of freedom are the `df_effect` and `df_error`
attributes. One `epsilon` describes both studies, because it is a
property of the covariance matrix of the repeated measures, which the
planned study shares with the prior one; the reported value is the prior
study's estimate of it, and re-running across the values you find
plausible shows how much the plan depends on it. The value used is
echoed as the `epsilon` row of the result, 1 when sphericity is assumed.

**The correction under `epsilon` is an approximation.** Below sphericity
the distribution of a within-subjects *F* depends on the whole
covariance matrix of the repeated measures, and on how the effect lines
up with it, not on \\\epsilon\\ alone. In simulations of raw data
(one-way designs of 4 and 6 levels with \\\epsilon\\ from .45 to .86,
prior studies of 12 to 40 subjects published when their
Greenhouse-Geisser test reached *p* \< .05, and a 2 x 4 split-plot
interaction), planning from each published study with its own reported
Greenhouse-Geisser estimate at an `assurance` of .80 gave a plan
reaching .80 power for the Greenhouse-Geisser test in .70 to .80 of
replications when the effect followed the dominant pattern of variation
(a linear trend over autocorrelated occasions), and in .92 to .99 when
it lay in the least variable contrasts. Ignoring the departure
(`epsilon = 1`) reached it in .26 to .64 of replications for the
trend-like effects and in .70 to .87 for the effects in the least
variable contrasts. The shortfall does not shrink as the prior study
grows, and is not larger near the lower bound of \\\epsilon\\. The
planned study's power is approximated the same way. Separate simulations
of planned one-way studies of 10 to 30 subjects, with 4, 6, or 8 levels,
measured the real power of the Greenhouse-Geisser test, whose estimate
of \\\epsilon\\ is biased downward. For a trend-like effect (a linear
trend over occasions with first-order autoregressive correlations of .3
to .9, \\\epsilon\\ from .50 to .94) it fell below `actual_power` by
about .04 to .08 with four levels and .06 to .13 with six, the larger
shortfalls for the smaller plans; a planned study analyzed with the
Huynh-Feldt correction tracks `actual_power` more closely. Even when
sphericity holds and `epsilon = 1`, a planned study analyzed with the
Greenhouse-Geisser correction anyway, as some software does by default,
has less power than `actual_power`, the power of the uncorrected test:
about .02 to .06 less with four levels, .03 to .11 with six, and .04 to
.15 with eight, for plans of 30 down to 10 subjects (the Huynh-Feldt
test is within about .02). When a paper reports both corrections, supply
the Greenhouse-Geisser one: the Huynh-Feldt estimate is larger, and
planning from it attained about .02 to .09 less, the larger differences
from the smaller prior studies. When sphericity in fact holds but the
prior study reported a Greenhouse-Geisser correction anyway, supplying
it is conservative, with a median plan 7 to 34 percent larger.

**The planning effect size and other power software.** Here `f_planning`
is Cohen's *f* relative to the effect's own error term, whose variance
already reflects the correlation among the repeated measures, so
\\f^2/(1 + f^2)\\ is the population value of the partial eta squared a
paper reports for the effect. It is not the *f* that power functions for
repeated measures designs usually take, such as
`ss_power_split_plot_anova()` in the sibling package `DMAR`: those
standardize the effect by the standard deviation of the raw scores and
take the correlation among the measures as a separate argument, so
giving them `f_planning`, with any correlation, plans a different study.
The value that reproduces this plan in `ss_power_split_plot_anova()`,
for any of the three effects, with the correlation set to 0 and the same
`epsilon`, is `f_planning` times \\\sqrt{k/b}\\, where *b* is
`levels_within` and *k* is the error degrees of freedom per subject (1
for the between-subjects effect, *b* - 1 for the other two). The
between-subjects effect is included: its `f_planning` is relative to the
variance of subject means, not of raw scores.

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

Greenhouse, S. W., & Geisser, S. (1959). On methods in the analysis of
profile data. *Psychometrika, 24,* 95–112.
[doi:10.1007/BF02289823](https://doi.org/10.1007/BF02289823)

Huynh, H., & Feldt, L. S. (1976). Estimation of the Box correction for
degrees of freedom from sample data in randomized block and split-plot
designs. *Journal of Educational Statistics, 1,* 69–82.
[doi:10.3102/10769986001001069](https://doi.org/10.3102/10769986001001069)

Muller, K. E., & Barton, C. N. (1989). Approximate power for
repeated-measures ANOVA lacking sphericity. *Journal of the American
Statistical Association, 84,* 549–555.
[doi:10.1080/01621459.1989.10478802](https://doi.org/10.1080/01621459.1989.10478802)

## Author

Ken Kelley (<kkelley@nd.edu>) and Samantha F. Anderson
(<samantha.f.anderson@asu.edu>)

## Examples

``` r
result <- ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
  levels_within = 3, effect = "within", alpha_prior = .05,
  alpha_planned = .05, assurance = .80, desired_power = .80)
result
#>  term                  value 
#>  necessary_n_per_group 969   
#>  total_N               1938  
#>  actual_power          0.8   
#>  ncp_adjusted          0.299 
#>  f_planning            0.0499
#>  F_observed            5     
#>  N                     60    
#>  levels_between        2     
#>  levels_within         3     
#>  epsilon               1     
#>  alpha_prior           0.05  
#>  alpha_planned         0.05  
#>  assurance             0.8   
#>  desired_power         0.8   
#> 
#> Design: Two-factor split-plot (mixed) ANOVA (within)
#> Sample size unit: per between-subjects cell
#> Largest supportable assurance: .83

# A prior study reported with Greenhouse-Geisser corrected degrees of
# freedom of 1.6 and 92.8, a correction of 1.6/2: its plan is larger than
# the one that assumes sphericity (60 per group).
ss_buc_mixed_anova(F_observed = 7, N = 60, levels_between = 2,
  levels_within = 3, effect = "within", epsilon = .8)
#>  term                  value
#>  necessary_n_per_group 108  
#>  total_N               216  
#>  actual_power          0.803
#>  ncp_adjusted          3.17 
#>  f_planning            0.163
#>  F_observed            7    
#>  N                     60   
#>  levels_between        2    
#>  levels_within         3    
#>  epsilon               0.8  
#>  alpha_prior           0.05 
#>  alpha_planned         0.05 
#>  assurance             0.8  
#>  desired_power         0.8  
#> 
#> Design: Two-factor split-plot (mixed) ANOVA (within)
#> Sample size unit: per between-subjects cell
#> Largest supportable assurance: .93

# Requesting more assurance than the prior result can support stops with an
# informative error naming the largest workable assurance (here near .83):
try(ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2,
  levels_within = 3, effect = "within", assurance = .95))
#> Error : Sample size planning is not possible with these settings. After correcting the prior study's effect for publication bias and uncertainty at the requested level of assurance, the corrected noncentrality parameter is zero, which means there is not enough information in the prior result to plan the sample size as desired. With this prior result and 'alpha_prior', the largest 'assurance' that still yields a usable plan is about .83; re-run with 'assurance' at or below that value, or raise 'alpha_prior' to lift the ceiling. Re-running at a lower 'assurance' is not free: the plan it returns carries that lower assurance, not the one you first asked for, so the long run proportion of replications reaching your target power falls accordingly. 'alpha_prior' does not change the prior study, which is what it is; it is your assumption about the publication threshold the study's literature faced, so raising it toward 1 assumes less publication bias to undo (for example, 'alpha_prior = .10' assumes the literature's journals would have published a result significant at the .10 level, a more lenient policy than .05, and 'alpha_prior = 1' assumes no publication bias at all). See Anderson, Kelley, and Maxwell (2017) and Anderson and Kelley (2024).
```
