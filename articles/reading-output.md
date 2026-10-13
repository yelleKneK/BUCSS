# Reading and Using BUCSS Output

Once a planner has returned, how do you read the plan, pull the one
number you need out of it, or collect it into a table or a manuscript?
Every `ss_buc_*` planner returns the same kind of object: a tidy result
of class `bucss_plan`, holding a bias and uncertainty corrected plan
(Anderson, Kelley, & Maxwell, 2017). This vignette shows the three ways
to read it, which follow the same conventions as the **DMAR** package,
so the habits carry over.

We will use a two-way between-subjects ANOVA throughout.

``` r

result <- ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2,
                                 levels_B = 3, effect = "factor_B")
```

## 1. Read It as a Printed Summary

Printing the object shows the aligned `term`/`value` table, in the
**DMAR** style: the design results first (the necessary sample size,
named for its unit; the implied total sample size for per-group and
per-cell designs; the achieved power at the returned size; the bias and
uncertainty adjusted noncentrality parameter; and the planning effect
size, the effect size the plan is built on, in the design’s usual
metric), then the rows echoing the planning inputs you supplied. Footer
lines close with the design, for the regression planners the predictor
frame (`random_predictors`), the unit the sample size is counted in, and
the largest assurance this prior result can support.

``` r

result
```

| term                 | value  |
|:---------------------|:-------|
| necessary_n_per_cell | 659    |
| total_N              | 3954   |
| actual_power         | 0.801  |
| ncp_adjusted         | 0.293  |
| f_planning           | 0.0494 |
| F_observed           | 5      |
| N                    | 120    |
| levels_A             | 2      |
| levels_B             | 3      |
| alpha_prior          | 0.05   |
| alpha_planned        | 0.05   |
| assurance            | 0.8    |
| desired_power        | 0.8    |

Design: Two-way between-subjects ANOVA (factor_B); Sample size unit: per
cell; Largest supportable assurance: .83

Displayed values are rounded for reading (whole numbers print clean;
other values print at three significant figures, and a value below 1 is
given more figures rather than shown as 1); the stored values keep full
precision. Change the display with `options(bucss.digits = )`, and put
back whatever was set before when you are done:

``` r

old <- options(bucss.digits = 5)
result
```

| term                 | value    |
|:---------------------|:---------|
| necessary_n_per_cell | 659      |
| total_N              | 3954     |
| actual_power         | 0.80057  |
| ncp_adjusted         | 0.29302  |
| f_planning           | 0.049415 |
| F_observed           | 5        |
| N                    | 120      |
| levels_A             | 2        |
| levels_B             | 3        |
| alpha_prior          | 0.05     |
| alpha_planned        | 0.05     |
| assurance            | 0.8      |
| desired_power        | 0.8      |

Design: Two-way between-subjects ANOVA (factor_B); Sample size unit: per
cell; Largest supportable assurance: .83

``` r

options(old)  # restore the previous setting
```

## 2. Pull a Single Number Out

The object is an ordinary `data.frame` with a character `term` column
and a numeric `value` column, so you can index it like any data frame.
The design results come first (here `necessary_n_per_cell`, `total_N`,
`actual_power`, `ncp_adjusted`, and `f_planning`; the size row is
`necessary_n_per_group` or `necessary_N` in other designs, and the
planning effect size row is named for its metric), followed by one row
per planning input, so the assumptions travel with the result through
subsetting and
[`write.csv()`](https://rdrr.io/r/utils/write.table.html).

``` r

result$value[result$term == "necessary_n_per_cell"]
#> [1] 659
result$value[result$term == "total_N"]
#> [1] 3954
```

That works for any row. For the sample size and its power,
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) below gives
them directly, and
[`glance()`](https://generics.r-lib.org/reference/glance.html) puts
every row in a column of its own, so any single value can be pulled out
by name.

## 3. The Broom Verbs: `tidy()` and `glance()`

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
compact estimate view, the same shape as **DMAR**’s sample size
planners: one row with `term` (`"sample_size"`), `estimate` (the
necessary sample size in the design’s unit), and `power` (the achieved
power at that size).

``` r

tidy(result)
#>          term estimate     power
#> 1 sample_size      659 0.8005714
tidy(result)$estimate
#> [1] 659
```

[`glance()`](https://generics.r-lib.org/reference/glance.html) returns
the one-row wide view. It begins with the same `term`, `estimate`, and
`power` columns as
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) (as
**DMAR**’s planners do), then gives every quantity and echoed planning
input a column of its own, and ends with the design, the effect tested,
the unit, the planned test’s degrees of freedom, and the assurance
ceiling. It is the convenient shape for pulling one value out by name or
collecting many results into one table for reporting.

``` r

glance(result)
#>          term estimate     power necessary_n_per_cell total_N actual_power
#> 1 sample_size      659 0.8005714                  659    3954    0.8005714
#>   ncp_adjusted f_planning F_observed   N levels_A levels_B alpha_prior
#> 1    0.2930235  0.0494152          5 120        2        3        0.05
#>   alpha_planned assurance desired_power                         design   effect
#> 1          0.05       0.8           0.8 Two-way between-subjects ANOVA factor_B
#>   sample_size_unit df_effect df_error assurance_ceiling
#> 1         per cell         2     3948         0.8342054
glance(result)$ncp_adjusted
#> [1] 0.2930235
```

Both verbs come from the **generics** package and are re-exported by
BUCSS, so they work as soon as the package is loaded, with no extra
setup. The full precision values are always available in the object
itself, since only the printed display rounds. And for a manuscript’s
method section,
[`planning_sentence()`](https://yelleKneK.github.io/BUCSS/reference/planning_sentence.md)
writes the sentence:

``` r

planning_sentence(result)
#> [1] "A planned study with 659 participants per cell (total N = 3,954) attains 80% power with 80% assurance, after correcting the prior study's result for publication bias (assuming its literature published results significant at the .05 level) and for uncertainty. The plan is built on a planning effect size of 0.0494, expressed as Cohen's f."
```

## The Planning Effect Size

The noncentrality parameter is the natural scale for the correction, but
it is not the scale a paper reports an effect on, and it grows with the
sample size. So every plan also carries its **planning effect size**, an
adjusted effect size: the value, in the design’s usual effect size
metric, used to plan the sample size given the situation at hand (the
prior study’s result, the publication bias assumed through
`alpha_prior`, and the uncertainty allowed for through `assurance`). It
is `ncp_adjusted` restated in that metric, and it is what the plan is
built on, not an estimate of the population effect to report. It sits
directly after `ncp_adjusted`, and its name, which always ends in
`_planning`, says which metric it is:

| Planner | Row | Metric |
|----|----|----|
| `ss_buc_independent_t` (`ss_buc_smd`) | `smd_planning` | standardized mean difference |
| `ss_buc_welch_t` | `smd_planning` | standardized mean difference, standardized by the root mean square of the two standard deviations |
| `ss_buc_paired_t` (`ss_buc_smd_paired`) | `dz_planning` | Cohen’s (1988) *d_z*: the mean of the paired differences divided by their standard deviation (the *z* names the difference score, which Cohen writes Z = X - Y) |
| `ss_buc_one_sample_t` | `std_mean_planning` | standardized mean |
| the eight ANOVA and ANCOVA planners | `f_planning` | Cohen’s *f* for the effect, relative to its own error term |
| `ss_buc_reg_coef`, `ss_buc_reg_joint` | `f2_planning` | Cohen’s *f* squared |
| `ss_buc_R2` | `R2_planning` | squared multiple correlation |
| `ss_buc_correlation` | `rho_planning` | correlation |
| `ss_buc_manova` | `D2_planning` | squared Mahalanobis distance |
| `ss_buc_chisq_diff` | `rmsea_d_planning` | RMSEA_D, with D for difference: the root mean square error of approximation of the chi square difference test between the two nested models (Savalei, Brace, & Fouladi, 2024), not the RMSEA of either model’s fit |

For the example above, the published result was *F*(2, 114) = 5, whose
sample Cohen’s *f* is `sqrt(5 * 2/114)`, about 0.296. The plan is built
on a much smaller effect:

``` r

f_planning <- result$value[result$term == "f_planning"]
f_planning
#> [1] 0.0494152
f_planning^2 / (1 + f_planning^2)  # the same effect as a partial eta squared
#> [1] 0.002435914
```

“Built on” is meant literally. `actual_power` is the ordinary power of a
study of the returned size when the population effect is the planning
effect size, so any power calculator reproduces it from the effect
alone. Here the planned study has 3954 participants in 6 cells, and its
test of the effect has 2 and 3948 degrees of freedom:

``` r

total_N <- result$value[result$term == "total_N"]
1 - pf(qf(.95, 2, total_N - 6), 2, total_N - 6, ncp = total_N * f_planning^2)
#> [1] 0.8005714
result$value[result$term == "actual_power"]
#> [1] 0.8005714
```

A few things to know when you carry the planning effect size elsewhere:

- **It is the value the plan is built on, not a better estimate of the
  effect.** At `assurance = .5` it is the median unbiased estimate (with
  `alpha_prior` below 1, after the correction for publication bias);
  above .5 it is deliberately conservative, a lower confidence limit, an
  effect the plan must still be powered for, so report it as the
  planning effect size rather than as the size of the effect.
- **Within-subjects and split-plot designs.** There `f_planning` is
  relative to the effect’s own error term, which already reflects the
  correlation among the repeated measures. Power functions that take an
  *f* standardized by the raw scores together with a correlation (such
  as **DMAR**’s `ss_power_rm_anova()` and `ss_power_split_plot_anova()`)
  expect a different number; each planner’s help page gives the
  conversion. With `epsilon` below 1 the planning effect size stays on
  the sphericity scale, like `ncp_adjusted`, and `actual_power` is the
  power under the same sphericity correction.
- **Regression with sampled predictors.** With
  `random_predictors = TRUE`, `f2_planning` and `R2_planning` are the
  population values when the predictors are sampled along with the
  outcome, and `actual_power` is the power in that frame.
- **An unbalanced prior study.** When the planner reduces the prior
  study to a balanced design in two ways and reports the smaller
  noncentrality parameter, the planning effect size belongs to that same
  balanced design. An ordinary power analysis at it can then return a
  slightly smaller sample size than the plan, never a larger one,
  because the plan is the larger of the two designs’ requirements.

[`planning_sentence()`](https://yelleKneK.github.io/BUCSS/reference/planning_sentence.md)
names the planning effect size in its second sentence, as shown in the
previous section.

## What Lives Where

- Every planning result and every echoed input is a row of the `value`
  column: the design results (the size row named for its unit,
  `total_N`, `actual_power`, `ncp_adjusted`, and the planning effect
  size) and the planning inputs. Rows survive any further arithmetic,
  subsetting, or export you do on the object. A switch such as the
  regression planners’ `random_predictors` is echoed as 1 (`TRUE`) or 0
  (`FALSE`), since the value column is numeric, and the print method’s
  footer names it in words. The within-subjects and split-plot planners
  echo `epsilon`, the sphericity parameter the plan assumed (1 under
  sphericity), for every effect it applies to.
- Everything else travels on attributes: the design label, the effect
  tested, and the unit, which are text, and the assurance ceiling and
  the planned test’s degrees of freedom, which are numbers describing
  the prior result and the planned test rather than rows of the plan.
  All of it is surfaced by
  [`glance()`](https://generics.r-lib.org/reference/glance.html), and
  the design (with the effect tested), the unit, and the ceiling by the
  print method’s footer.

## References

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28*, 1547–1562. <https://doi.org/10.1177/0956797617723724>

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Lawrence Erlbaum Associates.

Savalei, V., Brace, J. C., & Fouladi, R. T. (2024). We need to change
how we compute RMSEA for nested model comparisons in structural equation
modeling. *Psychological Methods, 29*, 480–493.
<https://doi.org/10.1037/met0000537>
