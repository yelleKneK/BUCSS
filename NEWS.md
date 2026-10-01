# BUCSS 2.0.0

This is a major release that modernizes the user-facing API, refines the
internal computation, and adds planning options. It is the first CRAN release
since 1.2.1.

The bias and uncertainty correction and the sample size computations follow
the same method as 1.2.1. The new options default to what 1.2.1 assumed: the
regression planners' `random_predictors` defaults to `FALSE` (predictors fixed
by design) and the within-subjects and split-plot planners' `epsilon` defaults
to 1 (sphericity), so a plan that uses neither is computed in the frame 1.2.1
used. For every documented example that 1.2.1 can also compute, the rounded
(three-decimal) adjusted noncentrality parameter matches the 1.2.1 value. The
main numerical change is the engine that finds the noncentrality parameter
(see "Engine" below), which shifts some planned sample sizes, most by a few
units, and makes them more accurate. Two smaller deliberate changes are
described under "Bug Fixes": the general between-subjects function's
`df_denominator` argument, silently ignored in 1.x, now functions, and the
independent *t* planner rounds the harmonic mean of two prior group sizes in
integer arithmetic, which changes the plan (to a smaller, correct one) for
pairs whose harmonic mean is a whole number. Every plan also gains one row,
its planning effect size (see "Sampled Predictors, Sphericity, and the
Planning Effect Size" below).

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
    rounded (three-decimal) parameters are unchanged from 1.2.1 (apart from
    the harmonic-mean cases above), but the extra precision shifts some
    planned sample sizes. For the documented unequal-*n* replication example
    (`ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .90)`),
    the necessary sample size moves from 1482 (the value BUCSS 1.x reported
    under its grid search) to 1485, which is the smallest sample size that
    actually reaches the target power. Across the developer characterization grid,
    leaving out the cases the old cap affected, the `df_denominator` change,
    and the harmonic-mean correction, about 95% of sample sizes are
    unchanged; most of the rest move by a few units and none by more than 2
    percent, except close to the assurance ceiling, where the corrected
    parameter is small and the grid's resolution mattered most (up to about 5
    percent). Among the documented examples, seven planned sample sizes move:
    1482 to 1485 (`ss_buc_independent_t`); 1904 to 1902 (`ss_buc_rm_anova`);
    968 to 969 (`ss_buc_mixed_anova`); 1491 to 1489
    (`ss_buc_mixed_anova_general`); 3963 to 3960 (`ss_buc_reg_joint`); and,
    in the tables of the vignette "Assurance and the Publication Bias
    Correction", 5069 to 5042
    (`ss_buc_independent_t(t_observed = 3, n = 20, assurance = .90)`) and 915
    to 913 (`ss_buc_paired_t()` with `alpha_prior = 1` and `assurance = .80`,
    for a prior *t* with a two-sided *p* value of .15 from 30 pairs). Each new
    value is the smallest sample size that actually attains the target power
    at the exact noncentrality parameter; they are pinned in the package's
    regression tests.
* Every worked example in Anderson, Kelley, and Maxwell (2017) and its
  supplemental material is reproduced, and the values are pinned in the
  regression tests. The article's three researchers plan from effect sizes it
  prints rounded (they are simulation medians; its notes 18 to 20), and at the
  printed values four of its sample sizes differ: *d* = 0.68 gives 330 per
  group where the article prints 332, an *f* squared of .09 gives 1,400 where
  it prints 1,418, and a paired *d* of .50 gives 50 and 80 where it prints 49
  and 79. Each printed size is reproduced by a value that rounds to the
  printed one (*d* = 0.6799, an *f* squared of .08996, and a paired *d* of
  .502), in BUCSS 1.2.1 as in 2.0.0.
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

* The functions now return a tidy object of class `bucss_plan` rather than a
  printed two-element list, shaped like the result tables of the sibling
  package `DMAR`. It is a `data.frame` with a character `term` column and a
  numeric `value` column. The design results come first: the necessary sample
  size, named for its unit (`necessary_n_per_group`, `necessary_n_per_cell`,
  or `necessary_N`); `total_N`, the implied total, for per-group and per-cell
  designs; `actual_power`, the power of the plan returned; and `ncp_adjusted`,
  whose scale (that of the prior study's *t*, *F*, or chi square statistic)
  each help page states. Rows echoing the planning inputs follow (a
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

* The deprecated functions return a `bucss_plan` object, but a `[[` method
  preserves the 1.x positional extraction: `result[[1]]` is the necessary sample
  size and `result[[2]]` is the adjusted noncentrality parameter, the same two
  values the 1.x functions returned as an unnamed list. Column and name access
  (`result$value`, `result[["value"]]`) are unaffected.

## New Designs

Six designs join the twelve that 1.2.1 covered. Two re-present computations
the package already made (the one-sample *t* and the analysis of covariance);
four extend the correction to statistics or distributions it did not
previously handle.

* `ss_buc_correlation()` plans from a prior Pearson correlation (or its *t*).
  A correlation study samples both variables, so its test statistic is not
  noncentral *F*: conditional on the sampled predictor it is, but the
  noncentrality parameter is itself random, proportional to the sampled spread
  of the predictor. The marginal distribution is a chi square mixture of
  noncentral *F* distributions, and it is more dispersed. The planner uses
  that mixture, which has an exact closed form (the Poisson weights of the
  noncentral *F* become negative binomial weights), for both the correction
  and the planned study's power, so the returned sample size is exact for
  bivariate normal data. Treating the statistic as noncentral *F*, which is
  correct only when the predictor is fixed by design, overstates the planned
  study's power by roughly 1/*n* for a plan of *n* participants (about .01 at
  *n* = 100 and .05 at *n* = 20) and returns a sample size typically 2 to 5
  participants light. When the predictor really is fixed, the analysis is a
  regression slope and `ss_buc_reg_coef()` with `p = 1` and
  `random_predictors = FALSE` is the planner for it; with
  `random_predictors = TRUE` that call returns exactly this planner's result.
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

* `ss_buc_chisq_diff()` plans from a nested model chi square difference
  (likelihood ratio) test, the test used for a constrained path in a
  structural equation model. It is *not* for an omnibus model fit chi square,
  where a publishable result is a small statistic and the selection region is
  inverted, and the function **refuses** that case rather than only warning
  about it. That is why it takes both models' degrees of freedom, `df_full`
  and `df_restricted`, rather than only their difference: an omnibus fit
  chi square is the comparison of a model against the saturated model, which
  has zero degrees of freedom, so `df_full = 0` is its signature. Given only
  the difference, nothing distinguishes a fit test from a genuine comparison
  of two substantive models and the mistake would pass silently. (A
  just-identified full model also has zero degrees of freedom; the help page
  and the error say how to plan that case.) The help page states the
  finite-sample liberality of the likelihood ratio test for a named model.

* `ss_buc_welch_t()` plans a Welch (unequal variance) *t* test from the two
  prior group sizes and the groups' relative spread. Only the ratio of the
  standard deviations enters the computation, but a paper usually reports the
  group standard deviations or variances rather than their ratio, so the
  spread may be given any of three ways: `sd_ratio`, or `sd_1` and `sd_2`, or
  `var_1` and `var_2`. Exactly one form is accepted per call, and the form
  given is the one echoed back in the result. A Welch statistic is not exactly
  noncentral *t*, so this planner is an approximation, and its help page says
  which assumption it needs and suggests re-running across a range of ratios.
  A minimum planned sample size of 10 per group is imposed because of
  estimation accuracy issues at smaller sample sizes: below 10 per group the
  approximation overstates the real Welch test's power, badly so for the
  smallest, while from 10 per group it is within .006 of the real test's
  power for every plan returned.

## The Assurance Ceiling, Drawn

* Added `plot.bucss_plan()`, so `plot()` on any plan draws its assurance
  ceiling. The method re-runs that plan's own planner across a range of
  assurance values and draws two panels sharing an assurance axis: the
  adjusted noncentrality parameter falling to zero at the ceiling, and the
  necessary sample size it implies growing asymptotically on a logarithmic
  scale. The
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

* The package-wide verification ships with the package. The simulation harness
  that validated the planners, twenty-four checks covering assurance
  calibration across every sampling family, raw-data audits of the
  correlation planner, the Welch planner, and the regression coefficient
  planner with sampled predictors, and the finite-sample behavior behind the
  chi square caveat, is installed under `inst/tools/` and reruns in a few
  seconds (under a minute in its full mode) with
  `source(system.file("tools", "verify_by_simulation.R", package = "BUCSS"))`.
  It reports the version and location of the BUCSS it verified.

* `ss_buc_sensitivity()` reports the Monte Carlo standard error of both
  attainment rates. Only the planner's own refusals count as refusals; any
  other error stops the audit with its message, and the planner's warnings are
  summarized once rather than repeated for every simulated study. A plan given
  as Hotelling's *T* squared can be audited like any other.

* Neither leaves a trace in the session. `ss_buc_sensitivity()`'s `seed`
  defaults to `NULL`, which draws from the caller's own random number stream;
  a supplied seed is set through `withr::local_seed()`, so the caller's
  generator state and kind are restored on exit (or removed again, if the
  caller had none), and the package never reads or writes the random number
  state itself. `withr` (>= 3.0.3; earlier versions of `local_seed()` did not
  restore anything when called from inside a function) joins `Imports`. The
  harness runs entirely inside `local()` with the same seed handling, so
  sourcing it creates or overwrites nothing in the global environment, leaves
  the caller's random number state as it was, and reports a failed check as
  an error rather than quitting the session.

## Sampled Predictors, Sphericity, and the Planning Effect Size

* The three regression planners, `ss_buc_reg_coef()`, `ss_buc_R2()`, and
  `ss_buc_reg_joint()`, take a `random_predictors` argument (default `FALSE`,
  the fixed-predictor frame of 1.2.1) for predictors sampled along with the
  outcome, as is typical of observational research in the behavioral
  sciences (Kelley & Maxwell, 2003). Whether the predictors are fixed by
  design or sampled changes neither the estimates nor the distribution of
  the test statistic when the effect is zero, but it changes the statistic's
  distribution, and so the power, when the effect is not zero (Sampson, 1974;
  Gatsonis & Sampson, 1989). With `random_predictors = TRUE` the
  noncentrality parameter is itself random, proportional to the sampled
  spread of the tested predictors once the other predictors are accounted
  for (chi square on *N* - 1 - (*p* - *q*) degrees of freedom when *q* of *p*
  predictors are tested), so the statistic follows a chi square mixture of
  noncentral *F* distributions. The mixture has an exact closed form (the
  Poisson weights of the noncentral *F* become negative binomial weights),
  used for both the correction and the planned study's power; it is the
  distribution `ss_buc_correlation()` uses: with one predictor each of the
  three returns exactly what `ss_buc_correlation()` returns for the same prior
  study (`ss_buc_reg_coef(p = 1)` for the same *t* and *N*, and
  `ss_buc_R2(p = 1)` and `ss_buc_reg_joint(p = 1, p_joint = 1)` for an
  `F_observed` equal to *t* squared), with the planning effect size stated in
  each planner's own metric. Treating sampled predictors as fixed overstates
  the planned study's power, by about 1/*n* for a plan of *n* participants
  with one predictor and about 4/*n* for a coefficient among eight (1/*n* to
  1.4/*n* for the squared multiple correlation). Which frame plans the larger
  study depends on the other predictors and on how far the plan outgrows the
  prior study; the help pages give the pattern. The largest supportable
  assurance is the same in both frames. The frame is echoed as a
  `random_predictors` row (1 or 0) on every regression plan, named in the
  printed footer, and followed by `ss_buc_sensitivity()` and `plot()`. The
  deprecated `ss.power.reg1()`, `ss.power.reg.all()`, and
  `ss.power.reg.joint()` keep the fixed-predictor computation of 1.x.

* The within-subjects and split-plot planners, `ss_buc_rm_anova()`,
  `ss_buc_rm_anova_general()`, `ss_buc_mixed_anova()`, and
  `ss_buc_mixed_anova_general()`, take an `epsilon` argument (the last
  argument, default 1) for a prior study reported with Greenhouse-Geisser or
  Huynh-Feldt corrected degrees of freedom: `epsilon` is the reported
  numerator degrees of freedom divided by the uncorrected ones. The prior
  study's *F* and the planned study's power are both treated under the
  approximation of Muller and Barton (1989), with the degrees of freedom and
  the noncentrality parameter multiplied by `epsilon` (the convention of
  `DMAR`'s `ss_power_rm_anova()` and `ss_power_split_plot_anova()`), and
  `ncp_adjusted` stays on the sphericity scale. It applies to within-subjects
  effects and to the interaction of a between-subjects and a within-subjects
  factor, each with its own lower bound; the between-subjects effect of a
  split-plot design involves no sphericity, so a value other than 1 is
  ignored there, with a warning. The correction is an approximation, and the
  help pages report its accuracy in simulations of raw data: at an
  `assurance` of .80, plans reached the target power in .70 to .80 of
  replications for trend-like effects and in .92 to .99 for effects in the
  least variable contrasts, against .26 to .64 and .70 to .87, respectively,
  when the departure was ignored. The default, `epsilon = 1`, reproduces the
  plans that assume sphericity exactly. Every plan for an effect with a
  within-subjects component echoes `epsilon` as a row (1 under sphericity),
  and `ss_buc_sensitivity()`, `plot()`, and `planning_sentence()` follow it.

* Every plan reports its planning effect size, an adjusted effect size: the
  value, in the design's usual effect size metric, used to plan the sample
  size given the situation at hand (the prior study's result, the publication
  bias assumed through `alpha_prior`, and the uncertainty allowed for through
  `assurance`). It is `ncp_adjusted` restated in that metric, as a row
  directly after `ncp_adjusted`, and it is what the plan is built on, not an
  estimate of the population effect to report: above an assurance of .5 it is
  deliberately conservative (a lower confidence limit), and at .5 it is the
  median unbiased estimate. The rows are `smd_planning` (the independent *t*
  planner and its alias `ss_buc_smd()`, and the Welch *t* planner, which
  standardizes by the root mean square of the two standard deviations,
  following Bonett, 2008), `dz_planning` (the paired *t* planner and its
  alias; Cohen's *d_z*, the mean of the paired differences divided by their
  standard deviation, where the *z* names the difference score, which Cohen,
  1988, writes Z = X - Y), `std_mean_planning` (the one-sample *t*),
  `f_planning` (Cohen's *f* relative to the effect's own error term, for all
  eight analysis of variance and covariance planners, so that *f* squared over
  1 plus *f* squared is the population partial eta squared for the effect),
  `f2_planning` (Cohen's *f* squared, for a regression coefficient and for a
  subset of predictors), `R2_planning` (the squared multiple correlation),
  `rho_planning` (the correlation), `D2_planning` (the squared Mahalanobis
  distance, for the two-group multivariate comparison), and `rmsea_d_planning`
  (RMSEA_D, with D for difference: the root mean square error of approximation
  of the chi square difference test between the two nested models, as distinct
  from the RMSEA of either model's fit, following Savalei, Brace, & Fouladi,
  2024). `actual_power` is exactly the ordinary power of a study of the
  returned sample size when the population effect is the planning effect size,
  in the plan's own frame (with sampled predictors, the random-predictor
  power; with `epsilon` below 1, the Muller and Barton approximation, the
  effect staying on the sphericity scale). For an unbalanced prior study the
  row belongs to the same balanced design as `ncp_adjusted`.
  `planning_sentence()` names it in a second sentence and `glance()` gives it
  a column; `tidy()` and the 1.x `[[1]]`/`[[2]]` positions are unchanged. For
  the within-subjects and split-plot planners, `f_planning` is not the
  raw-score *f* that repeated measures power functions such as `DMAR`'s
  `ss_power_rm_anova()` take; the help pages give the conversion.

## New Features

* The independent and paired *t* test planners are available under both
  test-specific names (`ss_buc_independent_t`, `ss_buc_paired_t`) and effect size
  names (`ss_buc_smd`, `ss_buc_smd_paired`) for the standardized mean
  difference. The two names in each pair are the same function.

* Added `print.bucss_plan()`, which prints the aligned `term`/`value` table
  in the `DMAR` display style (whole numbers clean, other values at
  `options(bucss.digits)` significant figures, stored values never rounded)
  and closes with factual footer lines naming the design, for the regression
  planners the predictor frame, the unit the sample size is counted in, and
  the largest assurance the prior result supports. A `knitr::knit_print`
  method renders the same table as a kable in R Markdown.

* The `effect` argument also accepts case-insensitive shorthand in the two-way
  between-subjects, within-subjects, and split-plot planners, for example `"A"`
  for `"factor_A"`, `"AxB"` for `"interaction"`, and `"bs"`/`"ws"` for
  `"between"`/`"within"`. The snake_case values remain the documented forms and
  the result always records the canonical one.

* Added `tidy()` and `glance()` methods for `bucss_plan` results (from the
  `generics` package), with the same division of labor as `DMAR`'s planners:
  `tidy()` is the compact estimate view (`term = "sample_size"`, `estimate`,
  and `power`), and `glance()` is the one-row wide view that begins with the
  same three columns and then has every quantity and echoed planning input as
  a column, plus the design metadata and the planned study's test degrees of
  freedom (`df_effect`, and `df_error` where the test has a denominator). Both
  stop with an explanation, rather than reporting one plan's numbers for
  another, when given several plans stacked into one object or a plan with
  one of its rows removed.

* Every planner reports `actual_power`, the statistical power of the plan it
  reports: the power of a study of the returned sample size when the effect is
  the returned adjusted noncentrality parameter (equivalently, the planning
  effect size). Because a sample size must be a whole number, it is never
  below the requested `desired_power` and is usually a little above it. Each
  help page states this definition.

* Added `planning_sentence()`, which turns a result into the sentence an
  author writes in a manuscript's planning section (the analog of `DMAR`'s
  `results_sentence()`). The sentence names only the corrections the plan
  actually made: the publication bias correction, with its assumed threshold,
  when `alpha_prior` is below 1, and the uncertainty correction when
  `assurance` is above .5, and it notes a departure from sphericity when
  `epsilon` is below 1. A second sentence names the planning effect size and
  its metric.

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
  sample size (`total_N`). Both are shown by `print.bucss_plan()` and travel on
  the object.

* Added an `inst/CITATION` file. `citation("BUCSS")` returns the package
  reference, so that the version and CRAN source are credited, together with the
  four articles behind the methods: Anderson, Kelley, and Maxwell (2017),
  Anderson and Kelley (2024), Anderson and Maxwell (2017), and the Anderson
  (2021) regression tutorial.

* Added a vignette on planning replication studies, which walks through why a
  literature built from underpowered studies inflates effect sizes and how the
  bias and uncertainty correction repairs the resulting sample size plan.

* Added a vignette on planning a correlation study, which explains why a
  correlation is not a fixed-predictor regression slope, when
  `ss_buc_correlation()` and `ss_buc_reg_coef()` each apply, and how to audit a
  correlation plan with `ss_buc_sensitivity()` and read its two attainment
  rates.

## Bug Fixes

* A negative observed *t* is now accepted by every planner that takes a *t*
  (`ss_buc_independent_t()`, `ss_buc_paired_t()`, `ss_buc_one_sample_t()`,
  `ss_buc_welch_t()`, `ss_buc_reg_coef()`, and `ss_buc_correlation()`), each
  of which plans from its magnitude. The publication rule the correction
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
  so also `power = 100`) is now rejected with a clear error. In 1.x an
  `assurance` of 0 and a `power` of 0 or 100 percent were accepted and
  silently returned meaningless results: an `assurance` of 0 has no solution,
  so the reported noncentrality parameter was an artifact of the internal
  search bound, and a target power of 100 percent is unattainable, so the
  returned sample size was where the computed power crossed 1 by
  floating-point error. An `assurance` of 100 percent (`assurance = 100`)
  already stopped in 1.x, but with the misleading message that the corrected
  noncentrality parameter was zero. The documented rule that an `assurance` or
  `power` of exactly 1 is read as 1 percent is unchanged.

* When the two prior group sizes given to `ss_buc_independent_t()` (or
  `ss.power.it`) have a harmonic mean that is a whole number, the planner no
  longer plans from a spurious second rounding one unit off. The mean was
  computed in floating point, where a whole number can land a hair away from
  itself: with `t_observed = 3` and the other arguments at their defaults,
  `n = c(49, 49)` planned 222 per group where `n = 49` plans 219, and
  `n = c(10, 15)`, whose harmonic mean is exactly 12, planned 177 where 12
  plans 154. The mean is now rounded in integer arithmetic. Only such inputs
  change, and the change is in the direction of a smaller, correct plan.

* For a prior study whose group sizes are unequal, the between-subjects and
  split-plot planners (and the independent *t* and MANOVA planners) reduce it
  to a balanced design rounded both down and up. When a statistic is
  significant at the prior study's own degrees of freedom but not at the
  rounded-down design, the error now says exactly that, rather than calling
  the prior result nonsignificant, and names an `alpha_prior` at which a plan
  can be made. And when the correction reaches zero, the error names the
  smaller of the two rounded designs' ceilings, so re-running at the
  assurance it names now produces a plan; before, it could name the larger
  one, which failed again.

* `ss.power.ba` and `ss.power.wa` refuse `effect = "factor.B"` (as well as the
  interaction) when `levels.B` is not given. 1.2.1 refused only the
  interaction, and silently returned the Factor A plan for `"factor.B"`.

* `ss_buc_factorial_anova_general()` and `ss_buc_ancova()` reject a
  `df_numerator` larger than `cells - 1`, and `ss_buc_mixed_anova_general()`
  rejects numerator degrees of freedom no design of its kind can have and
  warns that `df_num_within` is ignored except for the between-within
  interaction. Each used to return an ordinary-looking plan. The 1.x wrappers
  `ss.power.ba.general` and `ss.power.spa.general` inherit the same checks, so
  such 1.x calls now stop with an explanation (the wrapper passes
  `df.num.within` on only for the interaction, as 1.x used it, so 1.x calls
  that supply it with another effect do not warn).

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
  correctly describes `N` as the total sample size of the previous study.

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
  search and return a non-integer "sample size"
  (`ss.power.reg1(t.observed = 3, N = 150, p = 2.5)` returned 624.5), and
  `NA`, infinite, or vector inputs died inside base R with "missing value
  where TRUE/FALSE needed" or "the condition has length > 1". A within-subjects
  `N` of 1, or an independent-samples total `N` below 4, is likewise rejected
  up front instead of reaching `qt()`/`qf()` with zero degrees of freedom.

## Documentation

* Added the reference Maxwell, Delaney, and Kelley (2027, *Designing Experiments
  and Analyzing Data: A Model Comparison Perspective*, 4th ed.) to the package
  references, DOIs for Anderson and Maxwell (2017) and Taylor and Muller
  (1996), and ORCIDs for all three authors.

* The package help page, every planner's `assurance` argument, and
  `vignette("understanding-assurance", package = "BUCSS")` spell out how
  assurance is used. It is the long run proportion of planned studies whose
  power at the true effect reaches `desired_power`, when each is planned from
  a prior study published under the rule `alpha_prior` describes, and the
  planners deliver it by planning from the lower limit of a one-sided
  confidence interval at level `assurance`, computed from the
  publication-truncated distribution of the prior statistic. The term follows
  Kelley and Maxwell (2003), who used assurance in sample size planning for
  the probability that a confidence interval will be no wider than desired;
  that reference joins every planner's help page.

* Safeguard power (Perugini, Gallucci, & Costantini, 2014) guards against a
  prior effect size that is too large by planning for the lower limit of a
  two-sided 60% confidence interval around it rather than for the effect size
  itself; it allows for uncertainty but not for publication bias. It is shown
  to be a special case of the method: with `alpha_prior = 1` and
  `assurance = .80`, the adjusted noncentrality parameter is the lower limit
  of the two-sided 60% confidence interval that safeguard power plans from,
  and in general a two-sided confidence level *C* corresponds to
  `assurance = (1 + C)/2`. The package help page and the assurance vignette
  explain how, and where the two differ: every difference makes the plan here
  the same as the safeguard plan or larger, never smaller (for the regression
  planners, within the frame `random_predictors` names).
  `tests/testthat/test-safeguard-power.R` pins the equivalence.

* Every planner's `alpha_prior` argument and Details say plainly what the
  value means: the significance level a prior result had to reach to be
  published, so the correction assumes you see only results with *p* below
  it. `alpha_prior = 1` means no publication filter, since no *p* value
  exceeds 1, so there is no publication bias to correct and only the
  correction for uncertainty, set by `assurance`, remains. The planners
  compute with .999 in place of exactly 1, as BUCSS 1.x did, so plans match
  earlier versions; the argument's help and the package help page say so,
  and that the difference is negligible.

* The regression planners' help pages describe both predictor frames, what
  treating sampled predictors as fixed costs, and which frame plans the
  larger study (see "Sampled Predictors, Sphericity, and the Planning Effect
  Size" above). The regression and chi square difference pages also state
  where the planned sample size search starts (a total of *p* + 3, and the
  number of constraints plus 2), so no plan is smaller than that.

* The package help page and the Getting Started vignette list all eighteen
  planners.

## Internal

* The noncentrality search, the truncated-likelihood construction, and the
  planned sample size search each live in a single internal helper
  (`.solve_ncp_assurance`, `.tm_f`/`.tm_t`, and `.smallest_n_for_power`), shared
  by every planner, so the numerical core is maintained in one place.

* The zero-noncentrality and nonsignificance stops carry the condition classes
  `bucss_zero_ncp` and `bucss_nonsignificant`, so code that loops over plans
  (including `plot()` and `ss_buc_sensitivity()`) can catch a refusal without
  catching anything else. `plot()` now draws the plans of a prior result whose
  ceiling is below .5 instead of treating the warning for such an assurance as
  a refusal.

* R's warning that a noncentral *t* probability "may not have been computed to
  full precision", which the root finder triggered with large observed *t*
  statistics, is muffled where it arises; the solved parameter agrees with
  independent quadrature in every such case examined.

* The planned sample size search stops with a message rather than looping
  forever once the answer exceeds about 4.5e15 (2^52). Doubles hold every
  whole number exactly only up to 2^53, beyond which the bisection could not
  finish; answers below the bound are found exactly. Only an assurance a hair
  below the ceiling gets there. The random-predictor series (the correlation
  planner's, and the regression planners' with `random_predictors = TRUE`) is
  one internal function, summed over a window around its mode, so it stays
  exact for a strong correlation in a very large prior study (an *r* of .9
  with an *N* in the hundreds of thousands), where a fixed number of terms
  would have dropped most of its mass.

* Seeds are set through `withr` (>= 3.0.3), and R 3.6.0 is now the stated
  minimum, which the package's delayed S3 registration and its dependencies
  already required.

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
  across every planner. The `stats` functions are imported individually
  rather than wholesale, an `inst/WORDLIST` curates the domain terms for the
  spell check enabled by `Language: en-US`, and a GitHub Actions workflow runs
  R CMD check across platforms.

* Added a `testthat` (edition 3) test suite: regression oracles that pin the
  computed outputs, object-structure smoke tests, input-validation tests,
  coercion tests, and tests confirming each deprecated dotted name forwards to
  its `ss_buc_*` replacement. The new options have their own files
  (`test-random-predictors.R`, `test-epsilon.R`, `test-effect-rows.R`, and
  `test-safeguard-power.R`), which check each against an independent
  computation: the random-predictor mixture against quadrature and against
  the correlation planner, `epsilon` against the corrected-degrees-of-freedom
  power, the planning effect size against the ordinary power of the returned
  sample size, and safeguard power against confidence limits computed
  directly.
