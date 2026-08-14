# Package index

## t Tests

One-sample, independent, dependent (paired), and unequal-variance
(Welch) t test designs, that is, the standardized mean difference.

- [`ss_buc_independent_t()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md)
  [`ss_buc_smd()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md)
  : Necessary sample size to reach desired power for an independent t
  test (or standardized mean difference) using a publication bias and
  uncertainty correction procedure
- [`ss_buc_paired_t()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_paired_t.md)
  [`ss_buc_smd_paired()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_paired_t.md)
  : Necessary sample size to reach desired power for a dependent
  (paired) t test (or standardized mean difference) using a publication
  bias and uncertainty correction procedure
- [`ss_buc_one_sample_t()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_one_sample_t.md)
  : Necessary sample size to reach desired power for a one-sample t test
  using a publication bias and uncertainty correction procedure
- [`ss_buc_welch_t()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_welch_t.md)
  : Necessary sample size to reach desired power for a Welch (unequal
  variance) t test using a publication bias and uncertainty correction
  procedure

## Between-Subjects ANOVA and ANCOVA

One-way, two-way, and general between-subjects designs, with or without
covariates, planned per group or per cell.

- [`ss_buc_one_way_anova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_one_way_anova.md)
  : Necessary sample size to reach desired power for a one-way
  between-subjects ANOVA using a publication bias and uncertainty
  correction procedure
- [`ss_buc_factorial_anova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_factorial_anova.md)
  : Necessary sample size to reach desired power for a two-way
  between-subjects ANOVA using a publication bias and uncertainty
  correction procedure
- [`ss_buc_factorial_anova_general()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_factorial_anova_general.md)
  : Necessary sample size to reach desired power for a between-subjects
  ANOVA with any number of factors using a publication bias and
  uncertainty correction procedure
- [`ss_buc_ancova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_ancova.md)
  : Necessary sample size to reach desired power for an analysis of
  covariance using a publication bias and uncertainty correction
  procedure

## Within-Subjects ANOVA

Repeated measures designs, planned on the total number of subjects.

- [`ss_buc_rm_anova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_rm_anova.md)
  : Necessary sample size to reach desired power for a one or two-way
  within-subjects ANOVA using a publication bias and uncertainty
  correction procedure
- [`ss_buc_rm_anova_general()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_rm_anova_general.md)
  : Necessary sample size to reach desired power for a within-subjects
  ANOVA with any number of factors using a publication bias and
  uncertainty correction procedure

## Split-Plot (Mixed) ANOVA

Designs crossing a between-subjects factor with a within-subjects
factor.

- [`ss_buc_mixed_anova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_mixed_anova.md)
  : Necessary sample size to reach desired power for a two-factor
  split-plot (mixed) ANOVA using a publication bias and uncertainty
  correction procedure
- [`ss_buc_mixed_anova_general()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_mixed_anova_general.md)
  : Necessary sample size to reach desired power for a split-plot
  (mixed) ANOVA with any number of factors using a publication bias and
  uncertainty correction procedure

## Correlation and Multiple Regression

A Pearson correlation, a single coefficient, the omnibus model R2, or a
joint subset of predictors.

- [`ss_buc_correlation()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_correlation.md)
  : Necessary sample size to reach desired power for a Pearson
  correlation using a publication bias and uncertainty correction
  procedure
- [`ss_buc_reg_coef()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_coef.md)
  : Necessary sample size to reach desired power for a single
  coefficient in a multiple regression using a publication bias and
  uncertainty correction procedure
- [`ss_buc_R2()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_R2.md)
  : Necessary sample size to reach desired power for a test of model
  \\R^2\\ in a multiple regression using a publication bias and
  uncertainty correction procedure
- [`ss_buc_reg_joint()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_joint.md)
  : Necessary sample size to reach desired power for a joint test of
  multiple predictors in a multiple regression using a publication bias
  and uncertainty correction procedure

## Multivariate and Latent Variable Designs

A two-group multivariate comparison, and the chi-square difference test
between nested models.

- [`ss_buc_manova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_manova.md)
  : Necessary sample size to reach desired power for a two-group
  multivariate comparison using a publication bias and uncertainty
  correction procedure
- [`ss_buc_chisq_diff()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_chisq_diff.md)
  : Necessary sample size to reach desired power for a nested model
  chi-square difference test using a publication bias and uncertainty
  correction procedure

## Seeing and Auditing a Plan

Draw the assurance ceiling a prior result can support, and simulate the
literature a plan came from to see how often a plan built this way
actually reaches the target power.

- [`plot(`*`<bucss_power>`*`)`](https://yelleKneK.github.io/BUCSS/reference/plot.bucss_power.md)
  : Plot the assurance ceiling of a plan
- [`ss_buc_sensitivity()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_sensitivity.md)
  : Simulate the literature a plan came from and audit how often it
  works
- [`print(`*`<bucss_sensitivity>`*`)`](https://yelleKneK.github.io/BUCSS/reference/print.bucss_sensitivity.md)
  : Print a simulated sensitivity audit of a plan

## The Result Object

Working with the tidy bucss_power object the planners return.

- [`print(`*`<bucss_power>`*`)`](https://yelleKneK.github.io/BUCSS/reference/print.bucss_power.md)
  : Print a bias and uncertainty corrected sample size result
- [`tidy(`*`<bucss_power>`*`)`](https://yelleKneK.github.io/BUCSS/reference/bucss_tidiers.md)
  [`glance(`*`<bucss_power>`*`)`](https://yelleKneK.github.io/BUCSS/reference/bucss_tidiers.md)
  : Tidy and summarize a bias and uncertainty corrected sample size
  result
- [`planning_sentence()`](https://yelleKneK.github.io/BUCSS/reference/planning_sentence.md)
  : The planning sentence for a manuscript

## Deprecated 1.x Functions

The dot-named 1.x API, kept working as thin wrappers around the
ss_buc\_\* planners.

- [`ss.power.it()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  [`ss.power.dt()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  [`ss.power.ba()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  [`ss.power.ba.general()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  [`ss.power.wa()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  [`ss.power.wa.general()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  [`ss.power.spa()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  [`ss.power.spa.general()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  [`ss.power.reg1()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  [`ss.power.reg.all()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  [`ss.power.reg.joint()`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md)
  : Deprecated 1.x functions in BUCSS

## Package Overview

The package-level help page, with the method and its references.

- [`BUCSS`](https://yelleKneK.github.io/BUCSS/reference/BUCSS-package.md)
  [`BUCSS-package`](https://yelleKneK.github.io/BUCSS/reference/BUCSS-package.md)
  : BUCSS: Bias and Uncertainty Corrected Sample Size
