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

As of 2.0.0 the bias and uncertainty adjusted noncentrality parameter is found
by direct root finding (`stats::uniroot`) rather than a fixed grid, so it is no
longer capped and is exact rather than grid-snapped. The rounded parameters are
unchanged from 1.2.1; some planned sample sizes shift, most by a few units, as
documented in `NEWS.md` with the five documented examples that move. The other
deliberate numerical changes from 1.x are the `df_denominator` argument of
`ss_buc_factorial_anova_general`, which 1.x accepted but silently ignored and
which now functions, and the harmonic mean of two prior group sizes in the
independent *t* planner, which 1.x computed in floating point and which now
rounds in integer arithmetic (it could plan from a spurious rounding one unit
off). The regression tests in `tests/testthat/` pin every documented example.

The update also repairs the released version on current R: in 1.2.1 every call
to `ss.power.ba`, `ss.power.wa`, `ss.power.spa`, or `ss.power.spa.general`
that omits `effect` fails under R >= 4.2 with the length-greater-than-one
condition error, and a prior study whose implied noncentrality parameter
exceeds the old grid's cap of 100 silently returns an incorrect (capped or
grid-noise) result. Both are fixed by the 2.0.0 engine, so users on current R
are better served by this submission than by the version CRAN currently hosts.

## New in this version

Six designs join the twelve 1.2.1 covered, taking the family to eighteen:
`ss_buc_correlation`, `ss_buc_one_sample_t`, `ss_buc_ancova`, `ss_buc_manova`
(the rank-one multivariate case, Hotelling's *T* squared), `ss_buc_chisq_diff`
(a nested model chi square difference test), and `ss_buc_welch_t`. New
exported non-planners are `ss_buc_sensitivity`, a Monte Carlo companion that
simulates the published literature a plan came from and reports how often such
a plan attains its target power; a `plot` method that draws the largest
assurance a prior result can support; `planning_sentence`, which writes the
planning sentence for a manuscript; and the re-exported `tidy()` and
`glance()`. Each new function carries its own tests, examples, and help page,
and each approximation is stated in the help page rather than left implicit.
The Welch planner's approximation overstates power for very small planned
studies, so it never returns fewer than 10 per group, where it is accurate to
within .006.

One point deserves the reviewer's attention because it is a statistical
correction rather than a new feature. `ss_buc_correlation` does not treat the
correlation's test statistic as noncentral *F*. That is exact only when the
predictor is fixed by design, whereas a correlation study samples both
variables, and then the noncentrality parameter is itself random (proportional
to the sampled spread of the predictor, which is chi square distributed). The
marginal distribution is a chi square mixture of noncentral *F* distributions
and is more dispersed. The planner uses that mixture, which has an exact
closed form in `stats::dnbinom` and `stats::pbeta`, for both the correction
and the planned study's power. Treating the statistic as noncentral *F*
overstates the planned study's power by roughly 1/*n* for a plan of *n*
participants (about .01 at *n* = 100 and .05 at *n* = 20). This is verified by
simulation, and the package documents that a genuinely fixed predictor is a
regression slope and belongs in `ss_buc_reg_coef`.

`Imports` changes from `stats` alone to `generics` (for the `tidy()` and
`glance()` verbs), `graphics` (a base package, for the `plot` method), `stats`,
and `withr` (>= 3.0.3), through which `ss_buc_sensitivity()` sets a supplied
seed so the caller's random number state is restored on exit; the package
never reads or writes `.Random.seed` itself. Version 3.0.3 is required because
in earlier versions `withr::local_seed()` called inside a function does not
restore the state. `Suggests` gains `testthat`, `knitr`, and `rmarkdown`, and
`Depends` is now R (>= 3.6.0), which the package's delayed S3 registration
already required.

The simulation verification of the planners ships with the package
(`inst/tools/verify_by_simulation.R`, about 18 KB): twenty seeded checks that
any user can rerun in a few seconds (about half a minute in its full mode)
with `source(system.file("tools", "verify_by_simulation.R", package =
"BUCSS"))`, covering assurance calibration across every sampling family the
package uses and audits of the correlation and Welch planners against raw
simulated data. It runs inside `local()`, so sourcing it creates nothing in the
global environment and leaves the caller's random number state as it was; a
failed check is reported as an error. Nothing in the package's examples,
tests, or vignettes executes it, so it adds no check time.

## Test environments

* local: macOS 27 (arm64), R 4.6.1 (2026-06-24), 2026-09-29.

## R CMD check results

`R CMD check --as-cran` on the local machine:

    Status: 2 NOTEs

* "Found the following (possibly) invalid URLs: URL:
  https://yelleKneK.github.io/BUCSS/ ... Status: 404". This is the package's
  pkgdown site, listed in DESCRIPTION (URL), which is published from the
  package's GitHub repository before this submission.
* "Skipping checking HTML validation: 'tidy' doesn't look like recent enough
  HTML Tidy." The local copy of HTML Tidy is too old; this note is specific to
  the local machine.

## Reverse dependencies

There are no reverse dependencies on CRAN.
