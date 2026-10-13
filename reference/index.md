# Package index

## t Tests

One-sample, independent, dependent (paired), and unequal-variance
(Welch) t test designs, that is, the standardized mean difference.

- [`ss_buc_independent_t()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md)
  [`ss_buc_smd()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md)
  : Necessary Sample Size to Reach Desired Power for an Independent t
  Test (or Standardized Mean Difference) Using a Publication Bias and
  Uncertainty Correction Procedure
- [`ss_buc_paired_t()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_paired_t.md)
  [`ss_buc_smd_paired()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_paired_t.md)
  : Necessary Sample Size to Reach Desired Power for a Dependent
  (Paired) t Test (or Standardized Mean Difference) Using a Publication
  Bias and Uncertainty Correction Procedure
- [`ss_buc_one_sample_t()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_one_sample_t.md)
  : Necessary Sample Size to Reach Desired Power for a One-Sample t Test
  Using a Publication Bias and Uncertainty Correction Procedure
- [`ss_buc_welch_t()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_welch_t.md)
  : Necessary Sample Size to Reach Desired Power for a Welch (Unequal
  Variance) t Test Using a Publication Bias and Uncertainty Correction
  Procedure

## Between-Subjects ANOVA and ANCOVA

One-way, two-way, and general between-subjects designs, with or without
covariates, planned per group or per cell.

- [`ss_buc_one_way_anova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_one_way_anova.md)
  : Necessary Sample Size to Reach Desired Power for a One-Way
  Between-Subjects ANOVA Using a Publication Bias and Uncertainty
  Correction Procedure
- [`ss_buc_factorial_anova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_factorial_anova.md)
  : Necessary Sample Size to Reach Desired Power for a Two-Way
  Between-Subjects ANOVA Using a Publication Bias and Uncertainty
  Correction Procedure
- [`ss_buc_factorial_anova_general()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_factorial_anova_general.md)
  : Necessary Sample Size to Reach Desired Power for a Between-Subjects
  ANOVA With Any Number of Factors Using a Publication Bias and
  Uncertainty Correction Procedure
- [`ss_buc_ancova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_ancova.md)
  : Necessary Sample Size to Reach Desired Power for an Analysis of
  Covariance Using a Publication Bias and Uncertainty Correction
  Procedure

## Within-Subjects ANOVA

Repeated measures designs, planned on the total number of subjects.

- [`ss_buc_rm_anova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_rm_anova.md)
  : Necessary Sample Size to Reach Desired Power for a One- or Two-Way
  Within-Subjects ANOVA Using a Publication Bias and Uncertainty
  Correction Procedure
- [`ss_buc_rm_anova_general()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_rm_anova_general.md)
  : Necessary Sample Size to Reach Desired Power for a Within-Subjects
  ANOVA With Any Number of Factors Using a Publication Bias and
  Uncertainty Correction Procedure

## Split-Plot (Mixed) ANOVA

Designs crossing a between-subjects factor with a within-subjects
factor.

- [`ss_buc_mixed_anova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_mixed_anova.md)
  : Necessary Sample Size to Reach Desired Power for a Two-Factor
  Split-Plot (Mixed) ANOVA Using a Publication Bias and Uncertainty
  Correction Procedure
- [`ss_buc_mixed_anova_general()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_mixed_anova_general.md)
  : Necessary Sample Size to Reach Desired Power for a Split-Plot
  (Mixed) ANOVA With Any Number of Factors Using a Publication Bias and
  Uncertainty Correction Procedure

## Correlation and Multiple Regression

A Pearson correlation, a single coefficient, the omnibus model R2, or a
joint subset of predictors.

- [`ss_buc_correlation()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_correlation.md)
  : Necessary Sample Size to Reach Desired Power for a Pearson
  Correlation Using a Publication Bias and Uncertainty Correction
  Procedure
- [`ss_buc_reg_coef()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_coef.md)
  : Necessary Sample Size to Reach Desired Power for a Single
  Coefficient in a Multiple Regression Using a Publication Bias and
  Uncertainty Correction Procedure
- [`ss_buc_R2()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_R2.md)
  : Necessary Sample Size to Reach Desired Power for a Test of the
  Squared Multiple Correlation in a Multiple Regression Using a
  Publication Bias and Uncertainty Correction Procedure
- [`ss_buc_reg_joint()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_joint.md)
  : Necessary Sample Size to Reach Desired Power for a Joint Test of
  Multiple Predictors in a Multiple Regression Using a Publication Bias
  and Uncertainty Correction Procedure

## Multivariate and Latent Variable Designs

A two-group multivariate comparison, and the chi square difference test
between nested models.

- [`ss_buc_manova()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_manova.md)
  : Necessary Sample Size to Reach Desired Power for a Two-Group
  Multivariate Comparison Using a Publication Bias and Uncertainty
  Correction Procedure
- [`ss_buc_chisq_diff()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_chisq_diff.md)
  : Necessary Sample Size to Reach Desired Power for a Nested Model Chi
  Square Difference Test Using a Publication Bias and Uncertainty
  Correction Procedure

## Seeing and Auditing a Plan

Draw the assurance ceiling a prior result can support, and simulate the
literature a plan came from to see how often a plan built this way
actually reaches the target power.

- [`plot(`*`<bucss_plan>`*`)`](https://yelleKneK.github.io/BUCSS/reference/plot.bucss_plan.md)
  : Plot the Assurance Ceiling of a Plan
- [`ss_buc_sensitivity()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_sensitivity.md)
  : Simulate the Literature a Plan Came From and Audit How Often It
  Works
- [`print(`*`<bucss_sensitivity>`*`)`](https://yelleKneK.github.io/BUCSS/reference/print.bucss_sensitivity.md)
  : Print a Simulated Sensitivity Audit of a Plan

## The Result Object

Working with the tidy bucss_plan object the planners return.

- [`print(`*`<bucss_plan>`*`)`](https://yelleKneK.github.io/BUCSS/reference/print.bucss_plan.md)
  : Print a Bias and Uncertainty Corrected Sample Size Result
- [`tidy(`*`<bucss_plan>`*`)`](https://yelleKneK.github.io/BUCSS/reference/bucss_tidiers.md)
  [`glance(`*`<bucss_plan>`*`)`](https://yelleKneK.github.io/BUCSS/reference/bucss_tidiers.md)
  : Tidy and Summarize a Bias and Uncertainty Corrected Sample Size
  Result
- [`planning_sentence()`](https://yelleKneK.github.io/BUCSS/reference/planning_sentence.md)
  : The Planning Sentence for a Manuscript

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
  : Deprecated 1.x Functions in BUCSS

## Package Overview

The package-level help page, with the method and its references.

- [`BUCSS`](https://yelleKneK.github.io/BUCSS/reference/BUCSS-package.md)
  [`BUCSS-package`](https://yelleKneK.github.io/BUCSS/reference/BUCSS-package.md)
  : BUCSS: Bias and Uncertainty Corrected Sample Size
