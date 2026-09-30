## Release summary

BUCSS 2.0.0 is a major update of a package already on CRAN (1.2.1). It
modernizes the user-facing API: the planners are renamed to the `ss_buc_*`
prefix with snake_case arguments, and each returns a tidy `bucss_plan` object.
The dotted 1.x function names are retained but deprecated: they still run,
forwarding to the new functions and issuing a one-time-per-session warning, and
the returned object preserves the 1.x positional extraction (`result[[1]]`,
`result[[2]]`). Some 1.x calls now stop with an explanation instead of
returning a number, each because the number 1.x returned was wrong or
meaningless: degrees of freedom no design can have (`ss.power.ba.general`
with a `df.denominator` above `N - cells` or a `df.numerator` above
`cells - 1`, and the corresponding split-plot inputs to
`ss.power.spa.general`), `ss.power.ba` or `ss.power.wa` asked for Factor B
without `levels.B` (1.x silently returned the Factor A plan), an assurance or
power of 0 or of 100 percent, and fractional or non-scalar counts. `NEWS.md`
lists each.

## Why some results differ from 1.2.1

The bias and uncertainty adjusted noncentrality parameter is now found by
direct root finding (`stats::uniroot`) rather than a fixed grid, so it is no
longer capped at 100 and is exact rather than grid-snapped. For the documented
examples that 1.2.1 can also compute, the rounded parameters are unchanged
from 1.2.1; some planned sample sizes shift, most by a few units, as
documented in `NEWS.md`. Two other numerical changes are corrections: the
`df_denominator` argument of `ss_buc_factorial_anova_general`, which 1.x
accepted but silently ignored, now functions, and the independent *t* planner
rounds the harmonic mean of two prior group sizes in integer arithmetic (1.x
could plan from a spurious rounding one unit off). The regression tests in
`tests/testthat/` pin each documented example whose planned sample size
differs from 1.2.1.

The update also repairs the released version on current R: in 1.2.1 every call
to `ss.power.ba`, `ss.power.wa`, `ss.power.spa`, or `ss.power.spa.general`
that omits `effect` fails under R >= 4.2 with the length-greater-than-one
condition error, and a prior study whose implied noncentrality parameter
exceeds the old grid's cap of 100 silently returns an incorrect result. Both
are fixed in 2.0.0.

## New in this version

* New planners: `ss_buc_correlation`, `ss_buc_one_sample_t`, `ss_buc_ancova`,
  `ss_buc_manova`, `ss_buc_chisq_diff`, and `ss_buc_welch_t`, joining the
  twelve designs 1.2.1 covered.
* Other new exported functions: `ss_buc_sensitivity` (a Monte Carlo audit of
  a plan), a `plot` method, `planning_sentence`, and the re-exported `tidy()`
  and `glance()`.
* New arguments: `random_predictors` on `ss_buc_reg_coef`, `ss_buc_R2`, and
  `ss_buc_reg_joint` (predictors sampled rather than fixed by design), and
  `epsilon` on `ss_buc_rm_anova`, `ss_buc_rm_anova_general`,
  `ss_buc_mixed_anova`, and `ss_buc_mixed_anova_general` (a sphericity
  correction). Each is the last argument, and each defaults to the
  assumption 1.2.1 made (`random_predictors = FALSE`, `epsilon = 1`); the
  deprecated 1.x functions do not take them.
* New output rows: every plan reports its planning effect size, the adjusted
  noncentrality parameter restated in the design's usual effect size metric
  (for example `smd_planning`, `f_planning`, or `R2_planning`).

Each function has tests, runnable examples, and a help page.

## Dependencies

`Imports` changes from `stats` alone to `generics` (for `tidy()` and
`glance()`), `graphics` (a base package, for the `plot` method), `stats`, and
`withr` (>= 3.0.3), through which `ss_buc_sensitivity()` sets a supplied seed
and restores the caller's random number state on exit; the package never reads
or writes `.Random.seed` itself. Version 3.0.3 is required because in earlier
versions `withr::local_seed()` called inside a function does not restore the
state. `Suggests` gains `testthat`, `knitr`, and `rmarkdown`, and `Depends` is
now R (>= 3.6.0), which the package's delayed S3 registration already required.

The package also installs `inst/tools/verify_by_simulation.R` (about 25 KB), a
simulation check of the planners that users can rerun with
`source(system.file("tools", "verify_by_simulation.R", package = "BUCSS"))`.
It runs inside `local()`, leaves the caller's workspace and random number state
as they were, and is not run by the examples, tests, or vignettes.

## Test environments

* local: macOS 27 (arm64), R 4.6.1 (2026-06-24), 2026-09-30.

## R CMD check results

`R CMD check --as-cran` on the local machine:

    Status: 1 NOTE

* "Skipping checking HTML validation: 'tidy' doesn't look like recent enough
  HTML Tidy." The local copy of HTML Tidy is too old; this note is specific to
  the local machine.

## Reverse dependencies

There are no reverse dependencies on CRAN.
