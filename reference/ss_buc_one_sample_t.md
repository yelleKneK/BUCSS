# Necessary sample size to reach desired power for a one-sample t test using a publication bias and uncertainty correction procedure

`ss_buc_one_sample_t` returns the necessary total sample size to achieve
a desired level of statistical power for a planned study using a
one-sample *t* test, based on information obtained from a previous
study. The effect from the previous study can be corrected for
publication bias and/or uncertainty to provide a sample size that will
achieve more accurate statistical power for a planned study, when
compared to approaches that use a sample effect size at face value or
rely on sample size only. The bias and uncertainty adjusted previous
study noncentrality parameter is also returned, which can be transformed
to various effect size metrics.

## Usage

``` r
ss_buc_one_sample_t(
  t_observed,
  N,
  alpha_prior = 0.05,
  alpha_planned = 0.05,
  assurance = 0.8,
  desired_power = 0.8
)
```

## Arguments

- t_observed:

  Observed *t* value from a previous study used to plan sample size for
  a planned study. Either sign is accepted: the publication rule is
  two-sided, so only the magnitude enters the correction.

- N:

  Total sample size of the previous study.

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

The one-sample *t* test and the dependent (paired) *t* test are the same
test applied to a single column of scores: the statistic has *N* - 1
degrees of freedom and a noncentrality parameter of \\d\sqrt{N}\\.
`ss_buc_one_sample_t` therefore performs exactly the computation
[`ss_buc_paired_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_paired_t.md)
performs, and reports it with a one-sample design label and a total
sample size rather than a number of pairs. See
[`ss_buc_paired_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_paired_t.md)
for the full description of the method, the meaning of `assurance`, the
role of `alpha_prior`, and what happens when the corrected noncentrality
parameter is zero.

The returned `actual_power` is the statistical power of the plan this
function reports: the power of a study of the returned sample size when
the effect is the returned bias and uncertainty adjusted noncentrality
parameter. Because a sample size must be a whole number, `actual_power`
is never below the `desired_power` that was requested and is usually a
little above it.

The observed *t* may be entered with either sign; the publication rule
the correction assumes is two-sided, so only the magnitude enters the
computation.

If you are working from a standardized mean (Cohen's \\d\\) rather than
a *t* statistic, convert it before calling: \\t = d\sqrt{N}\\.

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

[`ss_buc_paired_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_paired_t.md),
which performs the identical computation for a paired design.

## Author

Ken Kelley (<kkelley@nd.edu>) and Samantha F. Anderson
(<samantha.f.anderson@asu.edu>)

## Examples

``` r
result <- ss_buc_one_sample_t(t_observed = 3, N = 40, alpha_prior = .05,
  alpha_planned = .05, assurance = .80, desired_power = .80)
result
#>  term          value
#>  necessary_N   255  
#>  actual_power  0.801
#>  ncp_adjusted  1.12 
#>  t_observed    3    
#>  N             40   
#>  alpha_prior   0.05 
#>  alpha_planned 0.05 
#>  assurance     0.8  
#>  desired_power 0.8  
#> 
#> Design: One-sample t test
#> Sample size unit: total
#> Largest supportable assurance: .90

# Requesting more assurance than the prior result can support stops with an
# informative error naming the largest workable assurance (here near .90):
try(ss_buc_one_sample_t(t_observed = 3, N = 40, assurance = .95))
#> Error : Sample size planning is not possible with these settings. After correcting the prior study's effect for publication bias and uncertainty at the requested level of assurance, the corrected noncentrality parameter is zero, which means there is not enough information in the prior result to plan the sample size as desired. With this prior result and 'alpha_prior', the largest 'assurance' that still yields a usable plan is about .90; re-run with 'assurance' at or below that value, or raise 'alpha_prior' to lift the ceiling. Re-running at a lower 'assurance' is not free: the plan it returns carries that lower assurance, not the one you first asked for, so the long run proportion of replications reaching your target power falls accordingly. 'alpha_prior' does not change the prior study, which is what it is; it is your assumption about the publication threshold the study's literature faced, so raising it toward 1 assumes less publication bias to undo (for example, 'alpha_prior = .10' assumes the literature's journals would have published a result significant at the .10 level, a more lenient policy than .05, and 'alpha_prior = 1' assumes no publication bias at all). See Anderson, Kelley, and Maxwell (2017) and Anderson and Kelley (2024).
```
