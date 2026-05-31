# BUCSS 2.0.0

This is a major release that modernizes the user-facing API. The corrected
noncentrality and sample size computations are unchanged from 1.2.1 and are
verified against it bit-for-bit, with one deliberate exception: the general
between-subjects function's `df_denominator` argument, which was silently
ignored in 1.x, now functions (see Bug fixes). At `df_denominator = N - cells`
that function still reproduces the 1.2.1 result exactly.

BUCSS 2.0.0 is **not backward compatible** with the 1.x API: function names,
argument names, and the returned object all changed. To run scripts written for
BUCSS 1.x without modification, install version 1.2.1 from the
[CRAN archive](https://cran.r-project.org/src/contrib/Archive/BUCSS/).

## Breaking changes

* The exported functions are renamed to the `ss_buc_*` prefix, so the names make
  clear that these planners apply the bias and uncertainty correction rather
  than ordinary power analysis. The between-subjects ANOVA planner is split into
  a one-way and a two-way function, and the two *t* test planners gain
  effect-size aliases:

  | 1.x | 2.0.0 |
  | --- | --- |
  | `ss.power.it` | `ss_buc_indep_t` (alias `ss_buc_smd`) |
  | `ss.power.dt` | `ss_buc_paired_t` (alias `ss_buc_smd_paired`) |
  | `ss.power.ba` | `ss_buc_one_way_anova`, `ss_buc_factorial_anova` |
  | `ss.power.ba.general` | `ss_buc_factorial_anova_general` |
  | `ss.power.wa` | `ss_buc_rm_anova` |
  | `ss.power.wa.general` | `ss_buc_rm_anova_general` |
  | `ss.power.spa` | `ss_buc_mixed_anova` |
  | `ss.power.spa.general` | `ss_buc_mixed_anova_general` |
  | `ss.power.reg1` | `ss_buc_reg_coef` |
  | `ss.power.reg.all` | `ss_buc_R2` |
  | `ss.power.reg.joint` | `ss_buc_reg_joint` |

* The old dotted names are defunct. Calling one (for example `ss.power.it()`)
  raises an error naming its replacement and pointing to the CRAN archive for
  the old behavior; it no longer forwards to the new function.

* Arguments are renamed to snake_case (for example `t.observed` becomes
  `t_observed`, `alpha.prior` becomes `alpha_prior`, `df.numerator` becomes
  `df_numerator`).

* `effect` values are now snake_case: `"factor.A"` becomes `"factor_A"`,
  `"between.only"` becomes `"between_only"`, and so on.

* The functions now return a tidy object of class `bucss_power` rather than a
  printed two-element list. It is a `data.frame` with a character `term` column
  and a numeric `value` column whose two rows are `necessary_sample_size` and
  `ncp_adjusted`. The design, effect, sample size unit, and planning inputs
  travel on attributes and are shown by the print method. Recover the two
  quantities with, for example,
  `result$value[result$term == "necessary_sample_size"]`.

## New features

* The independent and paired *t* test planners are available under both
  test-specific names (`ss_buc_indep_t`, `ss_buc_paired_t`) and effect-size
  names (`ss_buc_smd`, `ss_buc_smd_paired`) for the standardized mean
  difference. The two names in each pair are the same function.

* Added `print.bucss_power()`, which formats the result for humans (design,
  effect of interest, necessary sample size with its unit, adjusted
  noncentrality parameter, and the planning inputs).

* Added a `testthat` (edition 3) test suite: regression oracles that pin the
  1.2.1 outputs exactly, object-structure smoke tests, input-validation tests,
  and tests confirming each defunct dotted name errors with a useful signpost.

## Bug fixes

* `df_denominator` now functions in `ss_buc_factorial_anova_general()`. In 1.x
  it was accepted but silently ignored, and the documented example supplied an
  impossible value (`df_denominator = 117` with `N = 120` and `cells = 6`, for
  which the maximum residual degrees of freedom is `120 - 6 = 114`). The
  difference `(N - cells) - df_denominator` is now treated as the number of
  nuisance parameters (covariates or blocking factors) and carried through to
  the planned study: supplying `df_denominator = N - cells` reproduces the
  no-nuisance case (and the 1.2.1 result), while a smaller value correctly calls
  for a larger sample. Values greater than `N - cells` are now rejected.

* The error shown when `N` is omitted in `ss_buc_reg_coef()`, `ss_buc_R2()`, and
  `ss_buc_reg_joint()` incorrectly described `N` as the number of pairs. It now
  correctly describes `N` as the total sample size of the original study.

* The nonsignificance error in `ss_buc_R2()` and `ss_buc_reg_joint()` referred
  to `t.observed`. These functions test an `F` statistic, so the message now
  refers to `F_observed`.

* Omitting `effect` in the between-subjects, within-subjects, and split-plot
  functions no longer fails with "the condition has length > 1" under
  R >= 4.2. `effect` now defaults to its first choice through `match.arg()`.

## Internal

* Removed dead density assignments and redundant always-true guards in the
  iterative noncentrality search. Output is unchanged.
