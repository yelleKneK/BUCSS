# Necessary sample size to reach desired power for a between-subjects ANOVA with any number of factors using a publication bias and uncertainty correction procedure

`ss_buc_factorial_anova_general` returns the necessary per-cell sample
size to achieve a desired level of statistical power for a planned study
testing any type of effect (omnibus or contrast) using a fully
between-subjects ANOVA with any number of factors, based on information
obtained from a previous study. The effect from the previous study can
be corrected for publication bias and/or uncertainty to provide a sample
size that will achieve more accurate statistical power for a planned
study, when compared to approaches that use a sample effect size at face
value or rely on sample size only. The bias and uncertainty adjusted
previous study noncentrality parameter is also returned, which can be
transformed to various effect size metrics.

For the common one-way and two-way omnibus designs,
[`ss_buc_one_way_anova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_one_way_anova.md)
and
[`ss_buc_factorial_anova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_factorial_anova.md)
take the factor levels directly; this general form instead takes the
degrees of freedom, so it covers more factors and
single-degree-of-freedom contrasts.

## Usage

``` r
ss_buc_factorial_anova_general(
  F_observed,
  N,
  cells,
  df_numerator,
  df_denominator,
  alpha_prior = 0.05,
  alpha_planned = 0.05,
  assurance = 0.8,
  desired_power = 0.8
)
```

## Arguments

- F_observed:

  Observed *F* value from a previous study used to plan sample size for
  a planned study.

- N:

  Total sample size of the previous study.

- cells:

  Number of cells for the design (the product of the number of levels of
  each factor).

- df_numerator:

  Numerator degrees of freedom for the effect of interest.

- df_denominator:

  Denominator (error) degrees of freedom for the effect of interest in
  the prior study. This equals `N - cells` when the cell means are the
  only estimated parameters; supply a smaller value to account for
  covariates or other nuisance parameters. It cannot exceed `N - cells`.

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
per-cell sample size for the planned study, named for its unit
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

The approach implemented in `ss_buc_factorial_anova_general` uses the
observed *F* value and sample size from a previous study to correct the
noncentrality parameter associated with the effect of interest for
publication bias and/or uncertainty. This new estimated noncentrality
parameter is then used to calculate the necessary per-cell sample size
to achieve the desired level of power in the planned study.

The approach uses a likelihood function of a truncated noncentral *F*
distribution, where the truncation occurs due to small effect sizes
being unobserved due to publication bias. The numerator of the
likelihood function is the density of a noncentral *F* distribution. The
denominator is the power of the test, which serves to truncate the
distribution. Thus, the ratio of the numerator and the denominator is a
truncated noncentral *F* distribution. (See Taylor & Muller, 1996,
Equation 2.1, and Anderson & Maxwell, 2017, for more details.)

The denominator degrees of freedom you supply describe the prior study's
*F* test for the effect of interest. In a pure cell-means model the
residual degrees of freedom equal *N* minus the number of cells;
supplying a smaller value accommodates additional estimated parameters,
such as covariates or blocking factors. The difference,
`(N - cells) - df_denominator`, is treated as the number of such
nuisance parameters and is carried through to the planned study, whose
denominator degrees of freedom are reduced by the same amount. Supplying
`df_denominator = N - cells` reproduces the case with no such extra
parameters, and a smaller `df_denominator` calls for a larger sample
size. `df_denominator` cannot exceed `N - cells`.

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

`ss_buc_factorial_anova_general` assumes that the planned study will
have equal *n*. Unequal *n* in the previous study is handled in the
following way for between-subjects ANOVA designs. If the user enters an
*N* not equally divisible by the number of cells, the function
calculates *n* by dividing *N* by the number of cells and both rounding
up and rounding down the result, effectively assuming equal *n*. The
suggested sample size for the planned study is calculated using both of
these values of *n*, and the function returns the larger of these two
suggestions, to be conservative. The adjusted noncentrality parameter
that is output is the lower of the two possibilities, again, to be
conservative. Although equal-n previous studies are preferable, this
approach will work well as long as the cell sizes are only slightly
discrepant.

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

Ken Kelley (<kkelley@nd.edu>) and Samantha F. Anderson
(<samantha.f.anderson@asu.edu>)

## Examples

``` r
# With df_denominator = N - cells the general form reproduces the two-way
# omnibus result from ss_buc_factorial_anova (here 120 - 6 = 114).
result <- ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
  df_numerator = 2, df_denominator = 114, alpha_prior = .05,
  alpha_planned = .05, assurance = .80, desired_power = .80)
result
#>  term                 value
#>  necessary_n_per_cell 659  
#>  total_N              3954 
#>  actual_power         0.801
#>  ncp_adjusted         0.293
#>  F_observed           5    
#>  N                    120  
#>  cells                6    
#>  df_numerator         2    
#>  df_denominator       114  
#>  alpha_prior          0.05 
#>  alpha_planned        0.05 
#>  assurance            0.8  
#>  desired_power        0.8  
#> 
#> Design: Between-subjects ANOVA (general)
#> Sample size unit: per cell
#> Largest supportable assurance: .83

# Requesting more assurance than the prior result can support stops with an
# informative error naming the largest workable assurance (here near .83):
try(ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
  df_numerator = 2, df_denominator = 114, assurance = .95))
#> Error : Sample size planning is not possible with these settings. After correcting the prior study's effect for publication bias and uncertainty at the requested level of assurance, the corrected noncentrality parameter is zero, which means there is not enough information in the prior result to plan the sample size as desired. With this prior result and 'alpha_prior', the largest 'assurance' that still yields a usable plan is about .83; re-run with 'assurance' at or below that value, or raise 'alpha_prior' to lift the ceiling. Re-running at a lower 'assurance' is not free: the plan it returns carries that lower assurance, not the one you first asked for, so the long run proportion of replications reaching your target power falls accordingly. 'alpha_prior' does not change the prior study, which is what it is; it is your assumption about the publication threshold the study's literature faced, so raising it toward 1 assumes less publication bias to undo (for example, 'alpha_prior = .10' assumes the literature's journals would have published a result significant at the .10 level, a more lenient policy than .05, and 'alpha_prior = 1' assumes no publication bias at all). See Anderson, Kelley, and Maxwell (2017) and Anderson and Kelley (2024).
```
