# BUCSS <img src="man/figures/logo.png" align="right" height="139" alt="BUCSS hex logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/yelleKneK/BUCSS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yelleKneK/BUCSS/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/BUCSS)](https://CRAN.R-project.org/package=BUCSS)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/BUCSS)](https://CRAN.R-project.org/package=BUCSS)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

BUCSS is an R package for implementing Bias and Uncertainty Corrected Sample
Size planning. BUCSS implements a method of correcting for publication bias and
uncertainty when planning the sample size of a future study from an original
study. See [Anderson, Kelley, &amp; Maxwell (2017; *Psychological Science*, *28*, 1547--1562)](https://kenkelley.org/wp-content/uploads/articles/Anderson_Kelley_Maxwell_Psychological_Science_2017.pdf).

## Installing `BUCSS`

Install the released version from CRAN:

``` r
install.packages("BUCSS")
```

Then load the package:

``` r
library(BUCSS)
```

The `ss_buc_*` function family described below is the current BUCSS interface,
introduced in 2.0.0. The development version lives on GitHub:

``` r
# install.packages("remotes")
remotes::install_github("yelleKneK/BUCSS")
```

## Using `BUCSS`

Consider an original study in which two independent groups are used to test the
null hypothesis of no difference in the population means of the two groups
(e.g., a treatment and a control group). Suppose the independent samples
$t$ test in the original study reported $t=3.00$ based on sample sizes per group
of $n_1=50$ and $n_2=55$ with a Type I error rate of $\alpha=.05$. The planned
study seeks a statistical power of .80 with 90% assurance that the power will be
at least 80%. Assurance is the probability that a study planned with this method
will reach or surpass the desired level of statistical power. Note that an
assurance of .5 corrects for publication bias only, whereas assurance $> .5$
also corrects for uncertainty.

``` r
ss_buc_independent_t(t_observed = 3, n = c(50, 55), alpha_prior = .05,
               alpha_planned = .05, desired_power = .80, assurance = .90)
```

This yields a necessary sample size of $n_1=n_2=1485$ per group (i.e.,
$N=n_1+n_2=2970$). The 2017 paper reported 1482 for this example under the
original grid-search method; as of BUCSS 2.0.0 the bias and uncertainty adjusted
noncentrality parameter is found by direct root finding, which gives 1485 (see
`NEWS.md`).

The functions return a tidy object of class `bucss_power`, shaped like the
result tables in the [`DMAR`](https://github.com/yelleKneK/DMAR) package: an
ordinary `data.frame` of `term`/`value` rows holding the design results (the
size row is named for its unit, here `necessary_n_per_group`) followed by rows
echoing the planning inputs, so the assumptions travel with the result:

``` r
result <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), desired_power = .80,
                         assurance = .90)

result$value[result$term == "necessary_n_per_group"]  # 1485
result$value[result$term == "actual_power"]           # power at that size
result$value[result$term == "ncp_adjusted"]           # adjusted noncentrality
```

The broom verbs give the two convenient views, again exactly as in `DMAR`:
`tidy()` is the compact estimate view and `glance()` the one-row wide view
with every quantity and echoed input as a column:

``` r
tidy(result)            # term = "sample_size", estimate = 1485, power
tidy(result)$estimate   # 1485
glance(result)          # every row as a column, plus design metadata
planning_sentence(result)  # the sentence for a manuscript
```

## A Note on Function Names

BUCSS 2.0.0 renamed the function family to the `ss_buc_*` prefix (for example
`ss_buc_independent_t()` in place of the former `ss.power.it()`), so the names make
clear that these planners apply the bias and uncertainty correction, rather than
ordinary power analysis. The independent and paired *t* test planners also carry
effect size aliases (`ss_buc_smd()` and `ss_buc_smd_paired()`), and the
between-subjects ANOVA is now split into a one-way planner
(`ss_buc_one_way_anova()`) and a two-way planner (`ss_buc_factorial_anova()`).
Each function returns a tidy `bucss_power` object rather than a two-element list.

The family covers eighteen designs:

| Design | Function |
| --- | --- |
| One-sample, paired, independent, and Welch *t* tests | `ss_buc_one_sample_t()`, `ss_buc_paired_t()`, `ss_buc_independent_t()`, `ss_buc_welch_t()` |
| Between-subjects ANOVA and ANCOVA | `ss_buc_one_way_anova()`, `ss_buc_factorial_anova()`, `ss_buc_factorial_anova_general()`, `ss_buc_ancova()` |
| Within-subjects and split-plot ANOVA | `ss_buc_rm_anova()`, `ss_buc_rm_anova_general()`, `ss_buc_mixed_anova()`, `ss_buc_mixed_anova_general()` |
| Correlation and multiple regression | `ss_buc_correlation()`, `ss_buc_reg_coef()`, `ss_buc_R2()`, `ss_buc_reg_joint()` |
| Multivariate and latent variable | `ss_buc_manova()`, `ss_buc_chisq_diff()` |

`vignette("worked-examples", package = "BUCSS")` walks through one example of
each.

Calling `plot()` on a plan draws its assurance ceiling: the largest assurance
the prior result can support, the adjusted noncentrality parameter falling to
zero there, and the sample size growing asymptotically as it does. Assurance does not
degrade gracefully as it rises, it hits a wall, and that is easier to see than
to read about.

A plan can also be audited by simulation. `ss_buc_sensitivity()` takes the
result of any planner, simulates prior studies of that design from a true
effect you name, discards the ones a literature would not have published, runs
the same planner on each survivor, and reports how often a plan built this way
actually reaches the target power.

The old dotted names are retained for backward compatibility. As of BUCSS 2.0.0
each one (for example `ss.power.it()`) is **deprecated**: it still runs, but it
forwards to its `ss_buc_*` replacement and issues a one-time-per-session warning
naming the new function. The object it returns keeps the 1.x positional
extraction working, so legacy code that reads `result[[1]]` (the necessary
sample size) and `result[[2]]` (the adjusted noncentrality parameter) continues
to work unchanged. New code should call the `ss_buc_*` functions directly. See
the `NEWS.md` file for the full list of changes.

## Citing `BUCSS`

If you use BUCSS, please cite the package, which credits the software, its
authors, and the version you used:

> Kelley, K., Anderson, S. F., &amp; Maxwell, S. E. (2026). *BUCSS: Bias and
> Uncertainty Corrected Sample Size* (R package version 2.0.0).
> <https://CRAN.R-project.org/package=BUCSS>

For the full set of references, including the methodological articles the
correction relies on, run `citation("BUCSS")` in R.
