# Necessary sample size to reach desired power for a Welch (unequal variance) t test using a publication bias and uncertainty correction procedure

`ss_buc_welch_t` returns the necessary per-group sample size to achieve
a desired level of statistical power for a planned study using a Welch
*t* test, the two-group comparison that does not assume equal variances,
based on information obtained from a previous study. The effect from the
previous study can be corrected for publication bias and/or uncertainty
to provide a sample size that will achieve more accurate statistical
power for a planned study, when compared to approaches that use a sample
effect size at face value or rely on sample size only. The bias and
uncertainty adjusted previous study noncentrality parameter is also
returned, which can be transformed to various effect size metrics.

## Usage

``` r
ss_buc_welch_t(
  t_observed,
  n_1,
  n_2,
  sd_ratio = 1,
  sd_1,
  sd_2,
  var_1,
  var_2,
  alpha_prior = 0.05,
  alpha_planned = 0.05,
  assurance = 0.8,
  desired_power = 0.8
)
```

## Arguments

- t_observed:

  Observed Welch *t* value from a previous study used to plan sample
  size for a planned study. Either sign is accepted.

- n_1, n_2:

  Group sample sizes of the previous study.

- sd_ratio:

  Ratio of the second group's standard deviation to the first group's,
  in the previous study, assumed to hold in the planned study as well. A
  value of 1 means equal variances. Give this, or the two group standard
  deviations, or the two group variances, but not more than one of the
  three.

- sd_1, sd_2:

  Group standard deviations of the previous study, an alternative to
  `sd_ratio`. Only their ratio enters the computation, so the unit does
  not matter.

- var_1, var_2:

  Group variances of the previous study, a second alternative to
  `sd_ratio`.

- alpha_prior:

  Alpha level \\\alpha\\ for the previous study or the assumed
  statistical significance necessary for publishing in the field; to
  assume no publication bias, a value of 1 can be entered.

- alpha_planned:

  Alpha level \\\alpha\\ assumed for the planned study.

- assurance:

  Desired level of assurance, or the long run proportion of times that
  the planned study power will reach or surpass the desired level
  (assurance \> .5 corrects for uncertainty; assurance \< .5 is not
  recommended). Enter it as a proportion in (0, 1) or as a percentage
  greater than 1 (for example, 80 is read as .80); a value of exactly 1
  is read as 1 percent. A percentage is echoed in the result's
  planning-input rows as the coerced proportion (80 is echoed as .80);
  the same applies to `desired_power`.

- desired_power:

  Desired level of statistical power for the planned study. Enter it as
  a proportion in (0, 1) or as a percentage greater than 1; a value of
  exactly 1 is read as 1 percent.

## Value

An object of class `bucss_power`: a `data.frame` with a character `term`
column and a numeric `value` column, in the same shape as the sibling
package `DMAR`'s planners. The design results come first: the necessary
per-group sample size for the planned study, named for its unit
(`necessary_n_per_group`, `necessary_n_per_cell`, or `necessary_N`); for
the per-group and per-cell designs, `total_N`, the implied total sample
size; `actual_power`, the conservative achieved power at the returned
sample size (see Details); and `ncp_adjusted`, the publication bias and
uncertainty adjusted prior study noncentrality parameter. Rows echoing
the planning inputs follow, so the assumptions travel with the result.
The design, the unit the sample size is counted in, the largest
assurance this prior can support (the assurance ceiling), and the
planned test's degrees of freedom (`df_effect` and `df_error`) travel on
attributes; they are shown by `print.bucss_power`, and
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) give the
compact and wide views.

## Details

Researchers often use the sample effect size from a prior study as an
estimate of the likely size of an expected future effect in sample size
planning. However, sample effect size estimates should not usually be
used at face value to plan sample size, due to both publication bias and
uncertainty.

**This planner is an approximation, and needs one assumption the others
do not.** A Welch *t* statistic is not exactly distributed as a
noncentral *t*: its denominator degrees of freedom are estimated from
the data by the Welch-Satterthwaite formula, and they depend on the
ratio of the two population standard deviations. `ss_buc_welch_t` treats
the statistic as noncentral *t* on the Welch-Satterthwaite degrees of
freedom implied by the prior study's group sizes and the standard
deviation ratio you supply, which is the usual approximation in power
analysis for this test. In simulation that treatment is accurate: across
ratios from 1 to 3 the Type I error of the real Welch test stays between
.047 and .054 against a nominal .05, and the planned study's power is
within about .006 of what this function predicts, with the agreement no
worse at the extreme ratios.

**The assumption that matters is the ratio itself.** Planning assumes
the ratio you supply is the ratio the planned study will have, and
getting that wrong costs far more than the distributional approximation
does. For the example below, a plan built assuming equal variances
delivers about .43 power rather than .80 if the true ratio turns out to
be 2. Re-run across a range of plausible ratios and plan for the least
favorable one you find credible. The required sample size does not
always rise with the ratio: it falls when the group assumed to be more
variable is the smaller one.

When the ratio is 1 and the prior groups were the same size, this
planner returns exactly what
[`ss_buc_independent_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md)
returns. With unequal prior group sizes the two differ, and should: the
equal-variance planner reduces the prior study to an equivalent
equal-*n* design at the harmonic mean, while this one uses the actual
Welch-Satterthwaite degrees of freedom, which are smaller when the
groups are unbalanced. The same prior *t* is then less compelling and
the plan is larger.

The approach uses a likelihood function of a truncated noncentral *F*
distribution, where the truncation occurs due to small effect sizes
being unobserved due to publication bias. The numerator of the
likelihood function is the density of a noncentral *F* distribution. The
denominator is the power of the test, which serves to truncate the
distribution. In the two-group case, this formula reduces to the density
of a truncated noncentral *t* distribution. (See Taylor & Muller, 1996,
Equation 2.1, and Anderson & Maxwell, 2017, for more details.)

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

If the previous study of interest has not been subjected to publication
bias (e.g., a pilot study), `alpha_prior` can be set to 1 to indicate no
publication bias. Alternative \\\alpha\\ levels can also be accommodated
to represent differing amounts of publication bias. For example, setting
`alpha_prior = .20` would reflect less severe publication bias than the
default of .05. In essence, setting `alpha_prior` at .20 assumes that
studies with *p* values less than .20 are published, whereas those with
larger *p* values are not.

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

The observed *t* may be entered with either sign; the publication rule
the correction assumes is two-sided, so only the magnitude enters the
computation. The planned study is assumed to have equal group sizes,
which is the allocation that maximizes power when the variances are
equal and is a reasonable default otherwise. It is not the optimal
allocation when they are not: sampling in proportion to the standard
deviations would need about 11 percent fewer participants in total at a
ratio of 2, and about 25 percent fewer at a ratio of 3.

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

Anderson, S. F., & Maxwell, S. E. (2017). Addressing the 'replication
crisis': Using original studies to design replication studies with
appropriate statistical power. *Multivariate Behavioral Research, 52,*
305–324.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Taylor, D. J., & Muller, K. E. (1996). Bias in linear model power and
sample size calculation due to estimating noncentrality. *Communications
in Statistics: Theory and Methods, 25,* 1595–1610.

## See also

[`ss_buc_independent_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md)
for the equal-variance test.

## Author

Ken Kelley (<kkelley@nd.edu>)

## Examples

``` r
# A prior study with 40 and 55 participants whose second group was half
# again as variable as the first.
result <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55,
  sd_ratio = 1.5, alpha_prior = .05, alpha_planned = .05,
  assurance = .80, desired_power = .80)
result
#>  term                  value
#>  necessary_n_per_group 222  
#>  total_N               444  
#>  actual_power          0.801
#>  ncp_adjusted          1.32 
#>  t_observed            3    
#>  n_1                   40   
#>  n_2                   55   
#>  sd_ratio              1.5  
#>  alpha_prior           0.05 
#>  alpha_planned         0.05 
#>  assurance             0.8  
#>  desired_power         0.8  
#> 
#> Design: Welch (unequal variance) t test
#> Sample size unit: per group
#> Largest supportable assurance: .93

# The same study described by the group standard deviations, or by the
# group variances, rather than by their ratio.
ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_1 = 8.4, sd_2 = 12.6)
#>  term                  value
#>  necessary_n_per_group 222  
#>  total_N               444  
#>  actual_power          0.801
#>  ncp_adjusted          1.32 
#>  t_observed            3    
#>  n_1                   40   
#>  n_2                   55   
#>  sd_1                  8.4  
#>  sd_2                  12.6 
#>  alpha_prior           0.05 
#>  alpha_planned         0.05 
#>  assurance             0.8  
#>  desired_power         0.8  
#> 
#> Design: Welch (unequal variance) t test
#> Sample size unit: per group
#> Largest supportable assurance: .93
ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, var_1 = 70.56,
  var_2 = 158.76)
#>  term                  value
#>  necessary_n_per_group 222  
#>  total_N               444  
#>  actual_power          0.801
#>  ncp_adjusted          1.32 
#>  t_observed            3    
#>  n_1                   40   
#>  n_2                   55   
#>  var_1                 70.6 
#>  var_2                 159  
#>  alpha_prior           0.05 
#>  alpha_planned         0.05 
#>  assurance             0.8  
#>  desired_power         0.8  
#> 
#> Design: Welch (unequal variance) t test
#> Sample size unit: per group
#> Largest supportable assurance: .93

# How sensitive is the plan to the assumed spread?
ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1)
#>  term                  value
#>  necessary_n_per_group 213  
#>  total_N               426  
#>  actual_power          0.801
#>  ncp_adjusted          1.31 
#>  t_observed            3    
#>  n_1                   40   
#>  n_2                   55   
#>  sd_ratio              1    
#>  alpha_prior           0.05 
#>  alpha_planned         0.05 
#>  assurance             0.8  
#>  desired_power         0.8  
#> 
#> Design: Welch (unequal variance) t test
#> Sample size unit: per group
#> Largest supportable assurance: .92
ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 2)
#>  term                  value
#>  necessary_n_per_group 236  
#>  total_N               472  
#>  actual_power          0.801
#>  ncp_adjusted          1.31 
#>  t_observed            3    
#>  n_1                   40   
#>  n_2                   55   
#>  sd_ratio              2    
#>  alpha_prior           0.05 
#>  alpha_planned         0.05 
#>  assurance             0.8  
#>  desired_power         0.8  
#> 
#> Design: Welch (unequal variance) t test
#> Sample size unit: per group
#> Largest supportable assurance: .92

# Requesting more assurance than the prior result can support stops with an
# informative error naming the largest workable assurance (here near .93):
try(ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5,
  assurance = .95))
#> Error : Sample size planning is not possible with these settings. After correcting the prior study's effect for publication bias and uncertainty at the requested level of assurance, the corrected noncentrality parameter is zero, which means there is not enough information in the prior result to plan the sample size as desired. With this prior result and 'alpha_prior', the largest 'assurance' that still yields a usable plan is about .93; re-run with 'assurance' at or below that value, or raise 'alpha_prior' to lift the ceiling. Re-running at a lower 'assurance' is not free: the plan it returns carries that lower assurance, not the one you first asked for, so the long run proportion of replications reaching your target power falls accordingly. 'alpha_prior' does not change the prior study, which is what it is; it is your assumption about the publication threshold the study's literature faced, so raising it toward 1 assumes less publication bias to undo (for example, 'alpha_prior = .10' assumes the literature's journals would have published a result significant at the .10 level, a more lenient policy than .05, and 'alpha_prior = 1' assumes no publication bias at all). See Anderson, Kelley, and Maxwell (2017) and Anderson and Kelley (2024).
```
