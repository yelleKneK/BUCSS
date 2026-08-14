# Tidy and summarize a bias and uncertainty corrected sample size result

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
compact estimate view, mirroring the sibling package `DMAR`'s sample
size planners: one row with `term` (`"sample_size"`), `estimate` (the
necessary sample size in the design's unit), and `power` (the
conservative achieved power at that sample size; see the planner's help
page). [`glance()`](https://generics.r-lib.org/reference/glance.html)
returns the one-row wide view: the estimate plus one column per row of
the result (`total_N` where applicable, `ncp_adjusted`, and every echoed
planning input) together with the design, the effect tested, the unit,
the planned test's degrees of freedom, and the assurance ceiling.

## Usage

``` r
# S3 method for class 'bucss_power'
tidy(x, ...)

# S3 method for class 'bucss_power'
glance(x, ...)
```

## Arguments

- x:

  An object of class `bucss_power` returned by an `ss_buc_*` function.

- ...:

  Ignored.

## Value

A one-row `data.frame`. For
[`tidy()`](https://generics.r-lib.org/reference/tidy.html), the columns
`term`, `estimate`, and `power`. For
[`glance()`](https://generics.r-lib.org/reference/glance.html), the wide
summary described above.

## Details

These methods come from the generics package and need no extra setup;
the verbs are re-exported by BUCSS, so `tidy(result)` and
`glance(result)` work as soon as the package is loaded. The result is
also an ordinary `data.frame`, so the long view and the stored full
precision values remain available directly, for example
`result$value[result$term == "ncp_adjusted"]`.

## Author

Ken Kelley (<kkelley@nd.edu>)

## Examples

``` r
result <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
  levels_B = 3, effect = "factor_B")

# Compact estimate view, as in DMAR: the sample size and its power.
tidy(result)
#>          term estimate     power
#> 1 sample_size      659 0.8005714
tidy(result)$estimate
#> [1] 659

# One-row wide view: every quantity and echoed input as a column.
glance(result)
#>   necessary_n_per_cell total_N actual_power ncp_adjusted F_observed   N
#> 1                  659    3954    0.8005714    0.2930235          5 120
#>   levels_A levels_B alpha_prior alpha_planned assurance desired_power
#> 1        2        3        0.05          0.05       0.8           0.8
#>                           design   effect sample_size_unit df_effect df_error
#> 1 Two-way between-subjects ANOVA factor_B         per cell         2     3948
#>   assurance_ceiling
#> 1         0.8342054
glance(result)$total_N
#> [1] 3954
```
