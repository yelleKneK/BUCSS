# BUCSS 2.0.0

This is a major release that modernizes the user-facing API and refines the
internal computation. It is the first CRAN release since 1.2.1.

The corrected noncentrality and sample size computations follow the same method
as 1.2.1; the rounded (three-decimal) adjusted noncentrality parameters match
the 1.2.1 values for every documented example. The one substantive numerical
change is the engine that finds the noncentrality parameter (see "Engine"
below), which shifts a small number of planned sample sizes by a few units and
makes them more accurate. The general between-subjects function's
`df_denominator` argument, silently ignored in 1.x, now functions (see "Bug
fixes").

## Engine: root finding instead of a fixed grid

* The bias and uncertainty adjusted noncentrality parameter is found by solving
  the truncated-likelihood equation directly with `stats::uniroot()`, replacing
  the search over the fixed grid `seq(0, 100, by = step)` used in 1.x. Two
  consequences:
  * The adjusted parameter is **no longer capped at 100**. Prior studies with
    very large effects, for which the implied parameter exceeds 100, previously
    returned a silently truncated (incorrect) result; they are now solved
    correctly.
  * The parameter is exact rather than snapped to the grid resolution. The
    rounded (three-decimal) parameters are unchanged from 1.2.1, but the extra
    precision shifts a small number of planned sample sizes by a few units. For
    the documented unequal-*n* replication example
    (`ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .90)`), the
    necessary sample size moves from 1482 (the value quoted by Anderson et al.,
    2017, under the grid method) to 1485, which is the smallest sample size that
    actually reaches the target power. Across the developer characterization
    grid, about 95% of sample sizes are unchanged and the rest move by a few
    units.
* There is no grid-resolution argument or option. The 1.x `step` argument and
  any grid resolution setting are gone, since root finding has no resolution
  parameter.

## API modernization

* The exported functions are renamed to the `ss_buc_*` prefix, so the names make
  clear that these planners apply the bias and uncertainty correction rather
  than ordinary power analysis. The between-subjects ANOVA planner is split into
  a one-way and a two-way function, and the two *t* test planners gain
  effect-size aliases:

  | 1.x | 2.0.0 |
  | --- | --- |
  | `ss.power.it` | `ss_buc_independent_t` (alias `ss_buc_smd`) |
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

* Arguments are renamed to snake_case (for example `t.observed` becomes
  `t_observed`, `alpha.prior` becomes `alpha_prior`, `df.numerator` becomes
  `df_numerator`), and `effect` values are snake_case (`"factor.A"` becomes
  `"factor_A"`, `"between.only"` becomes `"between_only"`, and so on).

* The functions now return a tidy object of class `bucss_power` rather than a
  printed two-element list. It is a `data.frame` with a character `term` column
  and a numeric `value` column whose rows are `necessary_sample_size` and
  `ncp_adjusted` (plus `total_N` for per-group and per-cell designs). The design,
  effect, sample size unit, assurance ceiling, and planning inputs travel on
  attributes and are shown by the print method.

## Backward compatibility: 1.x names are deprecated, not removed

* The dot-named 1.x functions (`ss.power.it`, `ss.power.dt`, `ss.power.ba`,
  `ss.power.ba.general`, `ss.power.wa`, `ss.power.wa.general`, `ss.power.spa`,
  `ss.power.spa.general`, `ss.power.reg1`, `ss.power.reg.all`,
  `ss.power.reg.joint`) still work. Each accepts the 1.x argument names,
  translates them, and calls the corresponding `ss_buc_*` function, so scripts
  written for BUCSS 1.x continue to run. Calling one issues a deprecation
  warning that names its replacement, once per function per session; the
  `ss_buc_*` functions themselves are silent. `ss.power.ba` dispatches to
  `ss_buc_one_way_anova` when `levels.B` is `NULL` and to
  `ss_buc_factorial_anova` otherwise; the dotted `effect` values are translated
  to snake_case; and the 1.x `step` argument is accepted and ignored.

* The deprecated functions return a `bucss_power` object, but a `[[` method
  preserves the 1.x positional extraction: `result[[1]]` is the necessary sample
  size and `result[[2]]` is the adjusted noncentrality parameter, the same two
  values the 1.x functions returned as an unnamed list. Column and name access
  (`result$value`, `result[["value"]]`) are unaffected.

## New features

* The independent and paired *t* test planners are available under both
  test-specific names (`ss_buc_independent_t`, `ss_buc_paired_t`) and effect-size
  names (`ss_buc_smd`, `ss_buc_smd_paired`) for the standardized mean
  difference. The two names in each pair are the same function.

* Added `print.bucss_power()`, which formats the result for humans (design,
  effect of interest, necessary sample size with its unit, the implied total
  sample size where applicable, the adjusted noncentrality parameter, the
  assurance ceiling the prior result supports, and the planning inputs). The
  printed noncentrality parameter respects `options(bucss.digits)`.

* The `effect` argument also accepts case-insensitive shorthand in the two-way
  between-subjects, within-subjects, and split-plot planners, for example `"A"`
  for `"factor_A"`, `"AxB"` for `"interaction"`, and `"bs"`/`"ws"` for
  `"between"`/`"within"`. The snake_case values remain the documented forms and
  the result always records the canonical one.

* Added `tidy()` and `glance()` methods for `bucss_power` results (from the
  `generics` package, as in the `DMAR` package). `tidy()` returns a one-row data
  frame with one column per quantity, so a single value is pulled out by name
  (`tidy(result)$necessary_sample_size`); `glance()` returns a one-row summary.
  Both report the planned study's test degrees of freedom (`df_effect`,
  `df_error`).

* When the correction drives the corrected noncentrality parameter to zero, the
  error reports the largest assurance the prior result can support (the
  closed-form ceiling `1 - p/alpha_prior`, where `p` is the prior study's *p*
  value) and states whether lowering assurance to at or below that ceiling can
  reach a usable plan, or whether raising `alpha_prior` is the only remaining
  lever. The message clarifies that `alpha_prior` is the analyst's assumption
  about the publication threshold the prior study's literature faced, not a
  change to the prior study itself. Each planner's documentation states this
  ceiling and illustrates it with a runnable `try()` example.

* On a successful call, the result also reports the largest assurance the prior
  study can support and, for per-group and per-cell designs, the implied total
  sample size (`total_N`). Both are shown by `print.bucss_power()` and travel on
  the object.

* Added an `inst/CITATION` file. `citation("BUCSS")` returns the package
  reference, so that the version and CRAN source are credited, together with the
  Anderson, Kelley, and Maxwell (2017) and Anderson and Kelley (2024) articles.

* Added a vignette on planning replication studies, which walks through why a
  literature built from underpowered studies inflates effect sizes and how the
  bias and uncertainty correction repairs the resulting sample size plan.

## Bug fixes

* An `assurance` or `power` of exactly 0 or 1 (after the percentage coercion,
  so also `power = 100`) is now rejected with a clear error. In 1.x these
  endpoint values were accepted and silently returned meaningless results: an
  `assurance` of 0 has no solution, so the reported noncentrality parameter was
  an artifact of the internal search bound, and a target power of exactly 1 is
  unattainable, so the returned sample size was where the computed power
  crossed 1 by floating-point error. The documented rule that an `assurance` or
  `power` of exactly 1 is read as 1 percent is unchanged.

* `df_denominator` now functions in `ss_buc_factorial_anova_general()`. In 1.x
  it was accepted but silently ignored, and the documented example supplied an
  impossible value (`df_denominator = 117` with `N = 120` and `cells = 6`, for
  which the maximum residual degrees of freedom is `120 - 6 = 114`). The
  difference `(N - cells) - df_denominator` is now treated as the number of
  nuisance parameters (covariates or blocking factors) and carried through to
  the planned study: supplying `df_denominator = N - cells` reproduces the
  no-nuisance case, while a smaller value correctly calls for a larger sample.
  Values greater than `N - cells` are rejected.

* The error shown when `N` is omitted in `ss_buc_reg_coef()`, `ss_buc_R2()`, and
  `ss_buc_reg_joint()` incorrectly described `N` as the number of pairs. It now
  correctly describes `N` as the total sample size of the original study.

* The nonsignificance error in `ss_buc_R2()` and `ss_buc_reg_joint()` referred
  to `t.observed`. These functions test an `F` statistic, so the message now
  refers to `F_observed`.

* Omitting `effect` in the between-subjects, within-subjects, and split-plot
  functions no longer fails with "the condition has length > 1" under
  R >= 4.2. `effect` now defaults to its first choice through `match.arg()`.

* `ss_buc_rm_anova()` no longer silently returns the single-factor result when
  `effect` is `"factor_B"` or `"interaction"` but `levels_B` is not supplied; it
  now stops and explains that `levels_B` is required for those effects.

* The planners give clearer errors for more input mistakes: a missing structural
  argument (`cells`, `df_numerator`, `num_groups`, `df_num_within`, `p`, or
  `p_joint`) and a prior `N` too small for the design (fewer than two
  observations per cell or group) are now reported with an informative message.

* Every planner now validates that its observed statistic is a single finite
  number and that its design counts (`n`, `N`, levels, `cells`, `num_groups`,
  `p`, `p_joint`, prior degrees of freedom) are single whole numbers, naming
  the offending argument. In 1.x a fractional count could flow through the
  search and return a non-integer "sample size" (`p = 2.5` returned 624.5),
  and `NA`, infinite, or vector inputs died inside base R with "missing value
  where TRUE/FALSE needed" or "the condition has length > 1". A within-subjects
  `N` of 1, or an independent-samples total `N` below 4, is likewise rejected
  up front instead of reaching `qt()`/`qf()` with zero degrees of freedom.

## Documentation

* Added the reference Maxwell, Delaney, and Kelley (2027, *Designing Experiments
  and Analyzing Data: A Model Comparison Perspective*, 4th ed.) to the package
  references, and ORCIDs for all three authors.

## Internal

* The noncentrality search, the truncated-likelihood construction, and the
  planned-sample-size search each live in a single internal helper
  (`.solve_ncp_assurance`, `.tm_f`/`.tm_t`, and `.smallest_n_for_power`), shared
  by all twelve planners, so the numerical core is maintained in one place.

* Factored the shared planning-input validation into a single internal helper,
  and moved the duplicated `@references`, planning-argument, and `@return`
  documentation into `man-roxygen` templates, so each is maintained in one place
  across the twelve planners. The `stats` functions are imported individually
  rather than wholesale, an `inst/WORDLIST` curates the domain terms for the
  spell check enabled by `Language: en-US`, and a GitHub Actions workflow runs
  R CMD check across platforms.

* Added a `testthat` (edition 3) test suite: regression oracles that pin the
  computed outputs, object-structure smoke tests, input-validation tests,
  coercion tests, and tests confirming each deprecated dotted name forwards to
  its `ss_buc_*` replacement.
