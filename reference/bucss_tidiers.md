# Tidy and Summarize a Bias and Uncertainty Corrected Sample Size Result

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
compact estimate view, as the sibling package `DMAR`'s sample size
planners do: one row with `term` (`"sample_size"`), `estimate` (the
necessary sample size in the design's unit, which the result's
`sample_size_unit` attribute names), and `power` (the power of the plan
the result reports; see `actual_power` on the planner's help page).
[`glance()`](https://generics.r-lib.org/reference/glance.html) returns
the one-row wide view: the same three columns, then one column per row
of the result, namely the sample size row, `total_N` where applicable,
`actual_power`, `ncp_adjusted`, the planning effect size (the row whose
name ends in `_planning`, such as `smd_planning`), and every echoed
planning input, together with the design, the effect tested, the unit,
the planned test's degrees of freedom where the design has them, and the
assurance ceiling.

## Usage

``` r
# S3 method for class 'bucss_plan'
tidy(x, ...)

# S3 method for class 'bucss_plan'
glance(x, ...)
```

## Arguments

- x:

  An object of class `bucss_plan` returned by an `ss_buc_*` function.

- ...:

  Ignored.

## Value

A one-row `data.frame`. For
[`tidy()`](https://generics.r-lib.org/reference/tidy.html), the columns
`term`, `estimate`, and `power`. For
[`glance()`](https://generics.r-lib.org/reference/glance.html), the wide
summary described above. Both stop with an explanation if the object no
longer has the shape of a single plan (for example, several plans
stacked with `dplyr::bind_rows()`, or a plan with a row removed by
`dplyr::filter()`; both keep the class).

## Details

These methods come from the generics package and need no extra setup;
the verbs are re-exported by BUCSS, so `tidy(result)` and
`glance(result)` work as soon as the package is loaded. The result is
also an ordinary `data.frame`, so the long view and the stored full
precision values remain available directly, for example
`result$value[result$term == "ncp_adjusted"]`.

## Author

Ken Kelley (<kkelley@nd.edu>) and Samantha F. Anderson
(<samantha.f.anderson@asu.edu>)

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
#>          term estimate     power necessary_n_per_cell total_N actual_power
#> 1 sample_size      659 0.8005714                  659    3954    0.8005714
#>   ncp_adjusted f_planning F_observed   N levels_A levels_B alpha_prior
#> 1    0.2930235  0.0494152          5 120        2        3        0.05
#>   alpha_planned assurance desired_power                         design   effect
#> 1          0.05       0.8           0.8 Two-way between-subjects ANOVA factor_B
#>   sample_size_unit df_effect df_error assurance_ceiling
#> 1         per cell         2     3948         0.8342054
glance(result)$total_N
#> [1] 3954
```
