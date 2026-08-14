# Worked Examples Across Designs

This vignette walks through one example of each `ss_buc_*` function. The
pattern is always the same: summarize the prior study with a single test
statistic and its design, and the function returns the necessary
planned-study sample size and the bias and uncertainty adjusted
noncentrality parameter (Anderson, Kelley, & Maxwell, 2017). Every
example uses 80% power and 80% assurance unless noted. For the meaning
of the `assurance` and `alpha_prior` settings, see
[`vignette("understanding-assurance", package = "BUCSS")`](https://yelleKneK.github.io/BUCSS/articles/understanding-assurance.md).

## Independent *t* Test

A prior two-group study reported *t* = 3.00 with 20 participants per
group. The returned sample size is per group.

``` r

ss_buc_independent_t(t_observed = 3, n = 20, alpha_prior = .05, alpha_planned = .05,
               assurance = .80, desired_power = .80)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 130   |
| total_N               | 260   |
| actual_power          | 0.801 |
| ncp_adjusted          | 1.11  |
| t_observed            | 3     |
| n                     | 20    |
| alpha_prior           | 0.05  |
| alpha_planned         | 0.05  |
| assurance             | 0.8   |
| desired_power         | 0.8   |

Design: Independent t test; Sample size unit: per group; Largest
supportable assurance: .90

`ss_buc_smd` is the same planner under an effect size name, for when you
are working from a standardized mean difference rather than a raw *t*;
prefer `ss_buc_independent_t` when the prior *t* was computed from
unstandardized data.

## Dependent (Paired) *t* Test

A prior within-pairs study of 40 pairs reported *t* = 3.00. The returned
sample size is the number of pairs.

``` r

ss_buc_paired_t(t_observed = 3, N = 40, alpha_prior = .05, alpha_planned = .05,
                assurance = .80, desired_power = .80)
```

| term          | value |
|:--------------|:------|
| necessary_N   | 255   |
| actual_power  | 0.801 |
| ncp_adjusted  | 1.12  |
| t_observed    | 3     |
| N             | 40    |
| alpha_prior   | 0.05  |
| alpha_planned | 0.05  |
| assurance     | 0.8   |
| desired_power | 0.8   |

Design: Dependent (paired) t test; Sample size unit: number of pairs;
Largest supportable assurance: .90

The alias `ss_buc_smd_paired` serves the standardized mean difference
framing of the same comparison.

## Between-Subjects ANOVA

A one-way fully between-subjects design takes the number of levels of
the single factor. Here a four-group design with total *N* = 120
reported *F* = 5.00 for the omnibus effect. The returned sample size is
per group.

``` r

ss_buc_one_way_anova(F_observed = 5, N = 120, levels_A = 4, alpha_prior = .05,
                     alpha_planned = .05, assurance = .80, desired_power = .80)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 89    |
| total_N               | 356   |
| actual_power          | 0.802 |
| ncp_adjusted          | 3.74  |
| F_observed            | 5     |
| N                     | 120   |
| levels_A              | 4     |
| alpha_prior           | 0.05  |
| alpha_planned         | 0.05  |
| assurance             | 0.8   |
| desired_power         | 0.8   |

Design: One-way between-subjects ANOVA; Sample size unit: per group;
Largest supportable assurance: .94

A two-way fully between-subjects design takes the levels of both factors
and the effect of interest (`factor_A`, `factor_B`, or `interaction`).
Here a 2 x 3 design with total *N* = 120 reported *F* = 5.00 for the
main effect of Factor B. The returned sample size is per cell.

``` r

ss_buc_factorial_anova(F_observed = 5, N = 120, levels_A = 2, levels_B = 3,
                       effect = "factor_B", alpha_prior = .05,
                       alpha_planned = .05, assurance = .80, desired_power = .80)
```

| term                 | value |
|:---------------------|:------|
| necessary_n_per_cell | 659   |
| total_N              | 3954  |
| actual_power         | 0.801 |
| ncp_adjusted         | 0.293 |
| F_observed           | 5     |
| N                    | 120   |
| levels_A             | 2     |
| levels_B             | 3     |
| alpha_prior          | 0.05  |
| alpha_planned        | 0.05  |
| assurance            | 0.8   |
| desired_power        | 0.8   |

Design: Two-way between-subjects ANOVA (factor_B); Sample size unit: per
cell; Largest supportable assurance: .83

When the effect of interest is not a standard main effect or
interaction, or the design has more than two factors, use the general
form and supply the degrees of freedom directly. The denominator degrees
of freedom equal *N* minus the number of cells when the cell means are
the only estimated parameters, so supplying `df_denominator = N - cells`
(here 120 - 6 = 114) reproduces the two-way result above. A smaller
`df_denominator` accounts for covariates or other nuisance parameters
and calls for a larger sample.

``` r

ss_buc_factorial_anova_general(F_observed = 5, N = 120, cells = 6,
                               df_numerator = 2, df_denominator = 114,
                               alpha_prior = .05, alpha_planned = .05,
                               assurance = .80, desired_power = .80)
```

| term                 | value |
|:---------------------|:------|
| necessary_n_per_cell | 659   |
| total_N              | 3954  |
| actual_power         | 0.801 |
| ncp_adjusted         | 0.293 |
| F_observed           | 5     |
| N                    | 120   |
| cells                | 6     |
| df_numerator         | 2     |
| df_denominator       | 114   |
| alpha_prior          | 0.05  |
| alpha_planned        | 0.05  |
| assurance            | 0.8   |
| desired_power        | 0.8   |

Design: Between-subjects ANOVA (general); Sample size unit: per cell;
Largest supportable assurance: .83

## Within-Subjects ANOVA

For a fully within-subjects design, supply the factor levels and the
effect. Here a 2 x 3 within-subjects design with *N* = 60 reported *F* =
5.00 for the main effect of Factor B. The returned sample size is the
total number of subjects.

``` r

ss_buc_rm_anova(F_observed = 5, N = 60, levels_A = 2, levels_B = 3,
                effect = "factor_B", alpha_prior = .05, alpha_planned = .05,
                assurance = .80, desired_power = .80)
```

| term          | value |
|:--------------|:------|
| necessary_N   | 1902  |
| actual_power  | 0.8   |
| ncp_adjusted  | 0.304 |
| F_observed    | 5     |
| N             | 60    |
| levels_A      | 2     |
| levels_B      | 3     |
| alpha_prior   | 0.05  |
| alpha_planned | 0.05  |
| assurance     | 0.8   |
| desired_power | 0.8   |

Design: One or two-way within-subjects ANOVA (factor_B); Sample size
unit: total; Largest supportable assurance: .83

The general form takes the numerator degrees of freedom directly and
works for any number of within-subjects factors. This example uses
`assurance = .50`, which corrects for publication bias only.

``` r

ss_buc_rm_anova_general(F_observed = 6.5, N = 80, df_numerator = 1,
                        alpha_prior = .05, alpha_planned = .05, assurance = .50,
                        desired_power = .80)
```

| term          | value |
|:--------------|:------|
| necessary_N   | 256   |
| actual_power  | 0.8   |
| ncp_adjusted  | 2.47  |
| F_observed    | 6.5   |
| N             | 80    |
| df_numerator  | 1     |
| alpha_prior   | 0.05  |
| alpha_planned | 0.05  |
| assurance     | 0.5   |
| desired_power | 0.8   |

Design: Within-subjects ANOVA (any number of factors); Sample size unit:
total; Largest supportable assurance: .74

## Split-Plot (Mixed) ANOVA

A split-plot design crosses a between-subjects factor with a
within-subjects factor. Supply the levels of each and the effect
(`between`, `within`, or `interaction`). Here a design with 2
between-subjects levels and 3 within-subjects levels and *N* = 60
reported *F* = 5.00 for the within-subjects effect. The returned sample
size is per between-subjects cell.

``` r

ss_buc_mixed_anova(F_observed = 5, N = 60, levels_between = 2, levels_within = 3,
                   effect = "within", alpha_prior = .05, alpha_planned = .05,
                   assurance = .80, desired_power = .80)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 969   |
| total_N               | 1938  |
| actual_power          | 0.8   |
| ncp_adjusted          | 0.299 |
| F_observed            | 5     |
| N                     | 60    |
| levels_between        | 2     |
| levels_within         | 3     |
| alpha_prior           | 0.05  |
| alpha_planned         | 0.05  |
| assurance             | 0.8   |
| desired_power         | 0.8   |

Design: Two-factor split-plot (mixed) ANOVA (within); Sample size unit:
per between-subjects cell; Largest supportable assurance: .83

The general form takes the numerator degrees of freedom and the number
of between-subjects groups directly. The effect is one of
`between_only`, `within_only`, or `between_within`; supply
`df_num_within` when the effect involves a within-subjects component.
The returned sample size is per group.

``` r

ss_buc_mixed_anova_general(F_observed = 5, N = 90, df_numerator = 2,
                           num_groups = 3, effect = "between_only",
                           alpha_prior = .05, alpha_planned = .05,
                           assurance = .80, desired_power = .80)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 1489  |
| total_N               | 4467  |
| actual_power          | 0.8   |
| ncp_adjusted          | 0.194 |
| F_observed            | 5     |
| N                     | 90    |
| df_numerator          | 2     |
| num_groups            | 3     |
| alpha_prior           | 0.05  |
| alpha_planned         | 0.05  |
| assurance             | 0.8   |
| desired_power         | 0.8   |

Design: Split-plot (mixed) ANOVA (any number of factors) (between_only);
Sample size unit: per group; Largest supportable assurance: .82

## Multiple Regression

Anderson (2021) is a step-by-step tutorial for planning regression
studies with BUCSS, and is the best starting point for the three
regression planners below.

For a single regression coefficient, supply the *t* statistic for that
coefficient, the prior sample size, and the number of predictors `p`. A
prior study with *N* = 150 and 3 predictors reported *t* = 3.00 for one
coefficient. The returned sample size is the total.

``` r

ss_buc_reg_coef(t_observed = 3, N = 150, p = 3, alpha_prior = .05,
                alpha_planned = .05, assurance = .80, desired_power = .80)
```

| term          | value |
|:--------------|:------|
| necessary_N   | 624   |
| actual_power  | 0.8   |
| ncp_adjusted  | 1.89  |
| t_observed    | 3     |
| N             | 150   |
| p             | 3     |
| alpha_prior   | 0.05  |
| alpha_planned | 0.05  |
| assurance     | 0.8   |
| desired_power | 0.8   |

Design: Multiple regression: single coefficient; Sample size unit:
total; Largest supportable assurance: .93

For the omnibus test of the model R-squared, supply the model *F*
statistic. Here a prior study with *N* = 150 and 4 predictors reported
*F* = 5.00.

``` r

ss_buc_R2(F_observed = 5, N = 150, p = 4, alpha_prior = .05,
          alpha_planned = .05, assurance = .80, desired_power = .80)
```

| term          | value |
|:--------------|:------|
| necessary_N   | 234   |
| actual_power  | 0.8   |
| ncp_adjusted  | 7.82  |
| F_observed    | 5     |
| N             | 150   |
| p             | 4     |
| alpha_prior   | 0.05  |
| alpha_planned | 0.05  |
| assurance     | 0.8   |
| desired_power | 0.8   |

Design: Multiple regression: omnibus test of model R2; Sample size unit:
total; Largest supportable assurance: .98

For a joint test of a subset of predictors, also supply `p_joint`, the
number of predictors tested together. Here 2 of the 4 predictors are
tested jointly.

``` r

ss_buc_reg_joint(F_observed = 5, N = 150, p = 4, p_joint = 2,
                 alpha_prior = .05, alpha_planned = .05, assurance = .80,
                 desired_power = .80)
```

| term          | value |
|:--------------|:------|
| necessary_N   | 3960  |
| actual_power  | 0.8   |
| ncp_adjusted  | 0.365 |
| F_observed    | 5     |
| N             | 150   |
| p             | 4     |
| p_joint       | 2     |
| alpha_prior   | 0.05  |
| alpha_planned | 0.05  |
| assurance     | 0.8   |
| desired_power | 0.8   |

Design: Multiple regression: joint test of predictors; Sample size unit:
total; Largest supportable assurance: .84

## Correlation

A prior study reporting a Pearson correlation is planned from the
correlation itself. A study with *N* = 100 reported *r* = .35. The
returned sample size is the total.

``` r

ss_buc_correlation(r_observed = .35, N = 100, alpha_prior = .05,
                   alpha_planned = .05, assurance = .80, desired_power = .80)
```

| term          | value |
|:--------------|:------|
| necessary_N   | 122   |
| actual_power  | 0.802 |
| ncp_adjusted  | 6.71  |
| r_observed    | 0.35  |
| N             | 100   |
| alpha_prior   | 0.05  |
| alpha_planned | 0.05  |
| assurance     | 0.8   |
| desired_power | 0.8   |

Design: Pearson correlation; Sample size unit: total; Largest
supportable assurance: .99

The test of a correlation is arithmetically the test of a simple
regression slope, so it is tempting to plan it as one. That is exact
only when the predictor is fixed by design. A correlation study samples
both variables, and then the statistic’s noncentrality parameter is
itself random, proportional to the sampled spread of the predictor. The
resulting distribution is more dispersed, and planning as though it were
not overstates the power. The planner uses the exact distribution for
sampled variables, so it asks for a few more participants than the
regression planner does:

``` r

t_equivalent <- .35 * sqrt(98 / (1 - .35^2))
c(correlation = ss_buc_correlation(t_observed = t_equivalent, N = 100)$value[1],
  regression = ss_buc_reg_coef(t_observed = t_equivalent, N = 100, p = 1)$value[1])
#> correlation  regression 
#>         122         119
```

If your predictor really is fixed by design, the analysis is a
regression slope, and
[`ss_buc_reg_coef()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_coef.md)
with `p = 1` is the planner for it.

## Analysis of Covariance

Covariates cost error degrees of freedom in both the prior study and the
planned one, so the planner takes their number directly. A prior study
with *N* = 120, three groups, and two covariates reported *F* = 6.00 for
the group effect. The returned sample size is per cell.

``` r

ss_buc_ancova(F_observed = 6, N = 120, cells = 3, n_covariates = 2,
              alpha_prior = .05, alpha_planned = .05, assurance = .80,
              desired_power = .80)
```

| term                 | value |
|:---------------------|:------|
| necessary_n_per_cell | 154   |
| total_N              | 462   |
| actual_power         | 0.802 |
| ncp_adjusted         | 2.53  |
| F_observed           | 6     |
| N                    | 120   |
| cells                | 3     |
| n_covariates         | 2     |
| df_numerator         | 2     |
| alpha_prior          | 0.05  |
| alpha_planned        | 0.05  |
| assurance            | 0.8   |
| desired_power        | 0.8   |

Design: Analysis of covariance; Sample size unit: per cell; Largest
supportable assurance: .93

With `n_covariates = 0` this reproduces
[`ss_buc_one_way_anova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_one_way_anova.md)
exactly. For a single-degree-of-freedom contrast among the adjusted
means, supply that contrast’s *F* and set `df_numerator = 1`.

## Multivariate Comparison of Two Groups

When two groups are compared on several outcome variables at once, the
prior study’s Hotelling *T*-squared (or the *F* it maps to) plans the
replication. A study with *N* = 80 comparing two groups on 3 outcomes
reported *F* = 4.50. The returned sample size is per group.

``` r

ss_buc_manova(F_observed = 4.5, N = 80, p_variables = 3, alpha_prior = .05,
              alpha_planned = .05, assurance = .80, desired_power = .80)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 310   |
| total_N               | 620   |
| actual_power          | 0.801 |
| ncp_adjusted          | 1.42  |
| F_observed            | 4.5   |
| N                     | 80    |
| p_variables           | 3     |
| alpha_prior           | 0.05  |
| alpha_planned         | 0.05  |
| assurance             | 0.8   |
| desired_power         | 0.8   |

Design: Two-group multivariate comparison (Hotelling’s T squared);
Sample size unit: per group; Largest supportable assurance: .88

This covers the multivariate hypothesis whose rank is one, which is what
a two-group comparison gives. A hypothesis of higher rank has a set of
noncentrality parameters rather than a single one, and the correction
does not apply to it.

## Nested Model Comparison

An effect tested by the chi-square difference between two nested models,
such as a constrained path in a structural equation model, is planned
from that difference. Supply both models’ degrees of freedom rather than
only their difference. A prior study with *N* = 250 fit a model with 42
degrees of freedom and a constrained version with 43, a difference of
9.50 on 1 degree of freedom.

``` r

ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 42, df_restricted = 43,
                  alpha_prior = .05, alpha_planned = .05, assurance = .80,
                  desired_power = .80)
```

| term           | value |
|:---------------|:------|
| necessary_N    | 726   |
| actual_power   | 0.8   |
| ncp_adjusted   | 2.71  |
| chisq_observed | 9.5   |
| N              | 250   |
| df_full        | 42    |
| df_restricted  | 43    |
| alpha_prior    | 0.05  |
| alpha_planned  | 0.05  |
| assurance      | 0.8   |
| desired_power  | 0.8   |

Design: Nested model chi-square difference test; Sample size unit:
total; Largest supportable assurance: .95

This is for the *difference* test, not for an omnibus model fit
chi-square: a publishable fit statistic is a small one, which inverts
the selection the correction assumes. Requiring both models’ degrees of
freedom is what lets the function enforce that. An omnibus fit
chi-square compares a model against the *saturated* model, which has
zero degrees of freedom, so `df_full = 0` identifies the misuse exactly
and the function refuses it:

``` r

ss_buc_chisq_diff(chisq_observed = 9.5, N = 250, df_full = 0, df_restricted = 43)
#> Error:
#> ! 'df_full' is 0, so the full model is saturated and this comparison is the restricted model's omnibus test of model fit, not a difference test between two substantive models. This function does not apply to an omnibus fit chi-square: the correction assumes a literature publishes significant results, whereas a publishable fit chi-square is a small (nonsignificant) one, which inverts the selection region the method is built on. If you meant to compare two substantive nested models, give the full model's own degrees of freedom.
```

## Unequal Variances

When the two groups’ variances differ, the Welch test is the appropriate
comparison, and planning it needs the groups’ relative spread in
addition to their sizes. A prior study with 40 and 55 participants,
whose second group was half again as variable as the first, reported *t*
= 3.00.

``` r

ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5,
               alpha_prior = .05, alpha_planned = .05, assurance = .80,
               desired_power = .80)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 222   |
| total_N               | 444   |
| actual_power          | 0.801 |
| ncp_adjusted          | 1.32  |
| t_observed            | 3     |
| n_1                   | 40    |
| n_2                   | 55    |
| sd_ratio              | 1.5   |
| alpha_prior           | 0.05  |
| alpha_planned         | 0.05  |
| assurance             | 0.8   |
| desired_power         | 0.8   |

Design: Welch (unequal variance) t test; Sample size unit: per group;
Largest supportable assurance: .93

Only the ratio enters the computation, but a paper usually reports the
group standard deviations or variances rather than their ratio, so those
may be given directly instead. All three calls describe the same study:

``` r

sizes <- c(
  ratio = ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = 1.5)$value[1],
  sds = ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_1 = 8.4, sd_2 = 12.6)$value[1],
  vars = ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, var_1 = 70.56,
                        var_2 = 158.76)$value[1])
sizes
#> ratio   sds  vars 
#>   222   222   222
```

The assumed ratio matters more than any other input here, so it is worth
seeing how the plan moves across the values you find plausible:

``` r

ratios <- c(1, 1.5, 2, 2.5, 3)
sapply(ratios, function(r) {
  res <- ss_buc_welch_t(t_observed = 3, n_1 = 40, n_2 = 55, sd_ratio = r)
  res$value[res$term == "necessary_n_per_group"]
})
#> [1] 213 222 236 248 257
```

Planning at equal variances and being wrong is expensive: if the true
ratio were 2, that plan would deliver about .43 power rather than .80.
Plan for the least favorable ratio you find credible.

## References

Anderson, S. F. (2021). Using prior information to plan appropriately
powered regression studies: A tutorial using BUCSS. *Psychological
Methods, 26*, 513–526. <https://doi.org/10.1037/met0000366>

Anderson, S. F., & Kelley, K. (2024). Sample size planning for
replication studies: The devil is in the design. *Psychological Methods,
29*, 844–867. <https://doi.org/10.1037/met0000520>

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28*, 1547–1562. <https://doi.org/10.1177/0956797617723724>
