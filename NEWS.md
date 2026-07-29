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

## Engine: Root Finding Instead of a Fixed Grid

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
    units. Among the documented examples, five planned sample sizes move:
    1482 to 1485 and 1904 to 1902 (`ss_buc_independent_t`), 968 to 969
    (`ss_buc_mixed_anova`), 1491 to 1489 (`ss_buc_rm_anova`), and 3963 to 3960
    (`ss_buc_reg_joint`). Each new value is the smallest sample size that
    actually attains the target power at the exact noncentrality parameter;
    they are pinned in the package's regression tests.
* There is no grid-resolution argument or option. The 1.x `step` argument and
  any grid resolution setting are gone, since root finding has no resolution
  parameter.

## API Modernization

* The exported functions are renamed to the `ss_buc_*` prefix, so the names make
  clear that these planners apply the bias and uncertainty correction rather
  than ordinary power analysis. The between-subjects ANOVA planner is split into
  a one-way and a two-way function, and the two *t* test planners gain
  effect size aliases:

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
  `df_numerator`), the 1.x `power` argument is now `desired_power` (matching
  the `DMAR` package's planners and the `desired_power` row it echoes), and
  `effect` values are snake_case (`"factor.A"` becomes `"factor_A"`,
  `"between.only"` becomes `"between_only"`, and so on).

* The functions now return a tidy object of class `bucss_power` rather than a
  printed two-element list, shaped like the result tables of the sibling
  package `DMAR`. It is a `data.frame` with a character `term` column and a
  numeric `value` column. The design results come first: the necessary sample
  size, named for its unit (`necessary_n_per_group`, `necessary_n_per_cell`,
  or `necessary_N`); `total_N`, the implied total, for per-group and per-cell
  designs; `actual_power`, the conservative achieved power at the returned
  size; and `ncp_adjusted`. Rows echoing the planning inputs follow (a
  length-2 prior `n` echoes as `n_1`/`n_2`), so the assumptions travel with
  the result through subsetting and CSV export. Only non-numeric metadata
  (the design label, the effect tested, the unit) plus the assurance ceiling
  and the planned test's degrees of freedom travel on attributes.

## Backward Compatibility: 1.x Names Are Deprecated, Not Removed

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

## New Designs

Six designs join the twelve the package has always covered. Three are exact
within the published method and were already reachable by hand; three extend
the correction to statistics the package did not previously accept.

* `ss_buc_correlation()` plans from a prior Pearson correlation (or its *t*).
  A correlation study samples both variables, so its test statistic is not
  noncentral *F*: conditional on the sampled predictor it is, but the
  noncentrality parameter is itself random, proportional to the sampled spread
  of the predictor. The marginal distribution is a chi-square mixture of
  noncentral *F* distributions, and it is more dispersed. The planner uses
  that mixture, which has an exact closed form (the Poisson weights of the
  noncentral *F* become negative binomial weights), for both the correction
  and the planned study's power, so the returned sample size is exact for
  bivariate normal data. Treating the statistic as noncentral *F*, which is
  correct only when the predictor is fixed by design, overstates the planned
  study's power by up to about .04 and returns a sample size typically 2 to 5
  participants light. When the predictor really is fixed, the analysis is a
  regression slope and `ss_buc_reg_coef()` with `p = 1` is the planner for it.
  The largest supportable assurance is the same either way, because at a zero
  effect the mixture is the central *F* distribution.

* `ss_buc_one_sample_t()` plans a one-sample *t* test. This is exactly the
  computation `ss_buc_paired_t()` performs, reported with a one-sample design
  label and a total sample size instead of a number of pairs.

* `ss_buc_ancova()` plans an analysis of covariance from the number of
  covariates, rather than requiring the user to work out the error degrees of
  freedom for `ss_buc_factorial_anova_general()`. With no covariates it
  reproduces the corresponding ANOVA planner exactly. The omnibus group effect
  is the default; a single-degree-of-freedom contrast among the adjusted means
  is available through `df_numerator = 1`.

* `ss_buc_manova()` plans a two-group multivariate comparison (Hotelling's
  *T* squared), taking either the *F* or the *T* squared. This is the
  multivariate case whose hypothesis has rank one and therefore an exact
  noncentral *F*; hypotheses of higher rank have an eigenvalue spectrum rather
  than a single noncentrality parameter and are outside the method.

* `ss_buc_chisq_diff()` plans from a nested model chi-square difference
  (likelihood ratio) test, the test used for a constrained path in a
  structural equation model. Its help page is explicit that the function is
  *not* for an omnibus model fit chi-square, where a publishable result is a
  small statistic and the selection region is inverted.

* `ss_buc_welch_t()` plans a Welch (unequal variance) *t* test from the two
  prior group sizes and an assumed ratio of standard deviations. A Welch
  statistic is not exactly noncentral *t*, so this planner is an approximation
  and its help page says which assumption it needs and suggests re-running
  across a range of ratios.

## The Assurance Ceiling, Drawn

* Added `plot.bucss_power()`, so `plot()` on any plan draws its assurance
  ceiling. The method re-runs that plan's own planner across a range of
  assurance values and draws two panels sharing an assurance axis: the
  adjusted noncentrality parameter falling to zero at the ceiling, and the
  necessary sample size it implies running away on a logarithmic scale. The
  ceiling is marked, the plan itself is marked, and the region beyond the
  ceiling is shaded and labeled as one where no plan exists.

* The point of the figure is the shape. Assurance does not degrade gracefully
  as it rises, it hits a wall, and that is the hardest thing about the method
  to convey in prose. For a prior independent *t* test with *t* = 3.00 and 20
  per group, the plan runs from 25 per group at `assurance = .50` to over
  65,000 just short of the ceiling at .905, and is undefined above it.

* The grid of assurance values is spaced evenly in the logarithm of the
  remaining headroom rather than in assurance, so roughly half the plotted
  points fall in the last few hundredths below the ceiling. Spacing evenly in
  assurance would draw the wall as a single jump between two grid points and
  lose the shape.

* The plotted values are returned invisibly as an ordinary `data.frame`, one
  row per assurance, so the figure can be redrawn with any other graphics
  system. The method is written in base graphics; `graphics` joins `generics`
  and `stats` in `Imports`, and no plotting package is required.

## Auditing a Plan by Simulation

* Added `ss_buc_sensitivity()`, a Monte Carlo companion that takes a plan from
  any `ss_buc_*` function and exercises the claim the method actually makes.
  It reads the design and the planning inputs off the plan, simulates prior
  studies of that design from a true effect the user names, discards the ones
  a literature would not have published, runs the same planner on each
  survivor, and reports how often a plan built this way reaches the target
  power. There is deliberately no default for `true_ncp`: defaulting to the
  corrected estimate would ask the method to grade its own answer.

* The output distinguishes two attainment rates, and the distinction settles a
  question the method's documentation had left implicit. `attainment` counts
  only the prior studies for which a plan was issued, which is what a user
  experiences; `attainment_with_refusals` counts a refusal as attaining,
  because declining to plan is the limiting case of an arbitrarily
  conservative plan, and it is the rate the assurance guarantee is about. When
  the refusal rate is high the first sits well below `assurance`, which is the
  visible price of the method refusing rather than guessing rather than a
  defect.

## New Features

* The independent and paired *t* test planners are available under both
  test-specific names (`ss_buc_independent_t`, `ss_buc_paired_t`) and effect size
  names (`ss_buc_smd`, `ss_buc_smd_paired`) for the standardized mean
  difference. The two names in each pair are the same function.

* Added `print.bucss_power()`, which prints the aligned `term`/`value` table
  in the `DMAR` display style (whole numbers clean, other values at
  `options(bucss.digits)` significant figures, stored values never rounded)
  and closes with factual footer lines naming the design, the unit the sample
  size is counted in, and the largest assurance the prior result supports. A
  `knitr::knit_print` method renders the same table as a kable in R Markdown.

* The `effect` argument also accepts case-insensitive shorthand in the two-way
  between-subjects, within-subjects, and split-plot planners, for example `"A"`
  for `"factor_A"`, `"AxB"` for `"interaction"`, and `"bs"`/`"ws"` for
  `"between"`/`"within"`. The snake_case values remain the documented forms and
  the result always records the canonical one.

* Added `tidy()` and `glance()` methods for `bucss_power` results (from the
  `generics` package), with the same division of labor as `DMAR`'s planners:
  `tidy()` is the compact estimate view (`term = "sample_size"`, `estimate`,
  and `power`), and `glance()` is the one-row wide view with every quantity
  and echoed planning input as a column, plus the design metadata and the
  planned study's test degrees of freedom (`df_effect`, `df_error`).

* Every planner reports `actual_power`, the statistical power of the plan it
  reports: the power of a study of the returned sample size when the effect is
  the returned adjusted noncentrality parameter. Because a sample size must be
  a whole number, it is never below the requested `desired_power` and is
  usually a little above it. Each help page states this definition.

* Added `planning_sentence()`, which turns a result into the sentence an
  author writes in a manuscript's planning section (the analog of `DMAR`'s
  `results_sentence()`).

* The zero-noncentrality error now says plainly that re-running at a lower
  assurance is not free: the plan it returns carries that lower assurance, not
  the one originally requested, so the long run proportion of replications
  reaching the target power falls accordingly. Raising `alpha_prior` is a
  different kind of remedy, changing an assumption about the prior study's
  literature rather than lowering the guarantee. The assurance vignette
  develops the distinction.

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
  four articles behind the methods: Anderson, Kelley, and Maxwell (2017),
  Anderson and Kelley (2024), Anderson and Maxwell (2017), and the Anderson
  (2021) regression tutorial.

* Added a vignette on planning replication studies, which walks through why a
  literature built from underpowered studies inflates effect sizes and how the
  bias and uncertainty correction repairs the resulting sample size plan.

## Bug Fixes

* A negative observed *t* is now accepted by all three *t*-based planners
  (`ss_buc_independent_t()`, `ss_buc_paired_t()`, `ss_buc_reg_coef()`),
  which plan from its magnitude. The publication rule the correction
  assumes is two-sided and the truncated likelihood is symmetric in the
  sign of *t*, so `t_observed = -3` is the same evidence as
  `t_observed = 3`; the sign records which group was subtracted from
  which (or the direction of the paired difference, or the coefficient's
  coding), never the strength of the prior result. In 1.x and the
  earlier 2.0.0 drafts the two *t*-test planners rejected any negative
  *t* as "nonsignificant", however large its magnitude, with advice
  (raise `alpha_prior`) that could not help, while
  `ss.power.reg1`/`ss_buc_reg_coef` silently accepted the same input;
  the three planners now agree. The signed value is still echoed in the
  stored planning inputs, and a genuinely nonsignificant magnitude is
  still rejected. Scripts (and the designingexperiments.com apps) that
  previously received an error for a significant negative *t* now
  receive the plan its magnitude implies; no result that was previously
  returned has changed.

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

* `ss_buc_independent_t()` now echoes the prior study's sample size (`n` or
  `N`, whichever was supplied) in the result's planning-input rows, and hence
  in the printed table, as every other planner already did, so a saved or
  printed result can be traced back to the prior study.

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
  planned sample size search each live in a single internal helper
  (`.solve_ncp_assurance`, `.tm_f`/`.tm_t`, and `.smallest_n_for_power`), shared
  by all twelve planners, so the numerical core is maintained in one place.

* The planned sample size search brackets by doubling and then bisects instead
  of stepping by one, so its cost is logarithmic in the answer. Requests whose
  planned sample size runs to the millions (an `assurance` just below the
  prior's ceiling) now return instantly instead of iterating for minutes, and
  a huge answer can no longer overflow the integer step counter. Power is
  monotone in the sample size, so the returned value is the identical integer
  the incremental search produced (verified case for case on the developer
  characterization grid).

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
