# Necessary sample size to reach desired power for a nested model chi-square difference test using a publication bias and uncertainty correction procedure

`ss_buc_chisq_diff` returns the necessary total sample size to achieve a
desired level of statistical power for a planned study whose effect of
interest is tested by the chi-square difference between two nested
models (the likelihood ratio test used, for example, to test a
constrained path in a structural equation model), based on information
obtained from a previous study. The effect from the previous study can
be corrected for publication bias and/or uncertainty to provide a sample
size that will achieve more accurate statistical power for a planned
study, when compared to approaches that use a sample effect size at face
value or rely on sample size only. The bias and uncertainty adjusted
previous study noncentrality parameter is also returned.

## Usage

``` r
ss_buc_chisq_diff(
  chisq_observed,
  N,
  df_full,
  df_restricted,
  alpha_prior = 0.05,
  alpha_planned = 0.05,
  assurance = 0.8,
  desired_power = 0.8
)
```

## Arguments

- chisq_observed:

  Observed chi-square difference between the two nested models in the
  previous study.

- N:

  Total sample size of the previous study.

- df_full:

  Degrees of freedom of the full (less constrained) model. A saturated
  full model has zero, which identifies the comparison as an omnibus
  test of model fit rather than a difference test between two
  substantive models; the function refuses that case (see Details).

- df_restricted:

  Degrees of freedom of the restricted (more constrained) model. It
  exceeds `df_full` by the number of constraints imposed, which is the
  difference test's degrees of freedom.

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
total sample size for the planned study, named for its unit
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

**Scope, and one thing this function is not for.** The correction
applies to a test whose statistic is large when the effect is present,
so that publication selects large values. A chi-square difference test
between nested models is such a test: under the alternative it is
approximately noncentral chi-square with degrees of freedom equal to the
number of constraints released and a noncentrality parameter
proportional to sample size, and a paper reports the constrained path as
supported when the difference test is significant. The correction is
therefore built for the difference test.

It is *not* appropriate for an omnibus model fit chi-square, where a
publishable result is a *small* statistic (a model that fits). The
selection region is inverted there, so the truncated likelihood this
function builds does not describe it.

**The function refuses that case rather than only warning about it.**
This is why both models' degrees of freedom are required rather than
only their difference. An omnibus fit chi-square is the comparison of a
model against the *saturated* model, and a saturated model has zero
degrees of freedom, so `df_full = 0` identifies the misuse exactly.
Given only the difference, nothing distinguishes a fit test from a
genuine comparison of two substantive models, and the mistake would pass
silently. Supplying both is also a more faithful description of what was
tested, and the two values are what a paper reports for each model.

The correction uses a likelihood function of a truncated noncentral
chi-square distribution, the same construction the rest of the package
uses for the noncentral *F* (see Taylor & Muller, 1996, Equation 2.1,
and Anderson & Maxwell, 2017): the truncated area between the critical
value and the observed statistic, divided by the power of the test.
Because the noncentrality parameter of a likelihood ratio test is
proportional to sample size, the planned study's noncentrality parameter
is the corrected one scaled by the ratio of the planned to the prior
sample size.

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

The chi-square difference test is asymptotic, so both the correction and
the planned-study power inherit that approximation. The test is liberal
in small samples: in simulation its actual Type I error is nearer .06
than .05 at a total sample size of 50, and it reaches the nominal rate
only around 200. Treat a planned sample size under a hundred or so with
caution, and read a very small one as a signal that the prior effect was
large rather than as a serious recommendation.

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

## Author

Ken Kelley (<kkelley@nd.edu>)

## Examples

``` r
# A path constrained to zero: the full model has 42 degrees of freedom and
# the restricted model 43, so the difference test has 1.
result <- ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42,
  df_restricted = 43, alpha_prior = .05, alpha_planned = .05,
  assurance = .80, desired_power = .80)
result
#>  term           value
#>  necessary_N    726  
#>  actual_power   0.8  
#>  ncp_adjusted   2.71 
#>  chisq_observed 9.5  
#>  N              250  
#>  df_full        42   
#>  df_restricted  43   
#>  alpha_prior    0.05 
#>  alpha_planned  0.05 
#>  assurance      0.8  
#>  desired_power  0.8  
#> 
#> Design: Nested model chi-square difference test
#> Sample size unit: total
#> Largest supportable assurance: .95

# A saturated full model means the comparison is an omnibus test of model
# fit, which this method does not describe, and the function says so:
try(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 0,
  df_restricted = 43))
#> Error : 'df_full' is 0, so the full model is saturated and this comparison is the restricted model's omnibus test of model fit, not a difference test between two substantive models. This function does not apply to an omnibus fit chi-square: the correction assumes a literature publishes significant results, whereas a publishable fit chi-square is a small (nonsignificant) one, which inverts the selection region the method is built on. If you meant to compare two substantive nested models, give the full model's own degrees of freedom.

# Requesting more assurance than the prior result can support stops with an
# informative error naming the largest workable assurance (here near .95):
try(ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42,
  df_restricted = 43, assurance = .99))
#> Error : Sample size planning is not possible with these settings. After correcting the prior study's effect for publication bias and uncertainty at the requested level of assurance, the corrected noncentrality parameter is zero, which means there is not enough information in the prior result to plan the sample size as desired. With this prior result and 'alpha_prior', the largest 'assurance' that still yields a usable plan is about .95; re-run with 'assurance' at or below that value, or raise 'alpha_prior' to lift the ceiling. Re-running at a lower 'assurance' is not free: the plan it returns carries that lower assurance, not the one you first asked for, so the long run proportion of replications reaching your target power falls accordingly. 'alpha_prior' does not change the prior study, which is what it is; it is your assumption about the publication threshold the study's literature faced, so raising it toward 1 assumes less publication bias to undo (for example, 'alpha_prior = .10' assumes the literature's journals would have published a result significant at the .10 level, a more lenient policy than .05, and 'alpha_prior = 1' assumes no publication bias at all). See Anderson, Kelley, and Maxwell (2017) and Anderson and Kelley (2024).
```
