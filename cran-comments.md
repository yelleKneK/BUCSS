## Release summary

BUCSS 2.0.0 is a major update of a package already on CRAN (1.2.1). It
modernizes the user-facing API: the planners are renamed to the `ss_buc_*`
prefix with snake_case arguments, and each returns a tidy `bucss_power` object.
The dotted 1.x function names are retained but deprecated: they still run,
forwarding to the new functions and issuing a one-time-per-session warning, so
scripts written for 1.x keep working. The returned object preserves the 1.x
positional extraction (`result[[1]]`, `result[[2]]`).

As of 2.0.0 the bias and uncertainty adjusted noncentrality parameter is found
by direct root finding (`stats::uniroot`) rather than a fixed grid, so it is no
longer capped and is exact rather than grid-snapped. The rounded parameters are
unchanged from 1.2.1; a small number of planned sample sizes shift by a few
units, documented in `NEWS.md`. The other deliberate change from 1.x is the
`df_denominator` argument of `ss_buc_factorial_anova_general`, which 1.x accepted
but silently ignored and which now functions. The regression tests in
`tests/testthat/` pin every documented example.

The update also repairs the released version on current R: in 1.2.1 the
one-way paths of `ss.power.ba` and `ss.power.wa` fail under R >= 4.2 with the
length-greater-than-one condition error, and a prior study whose implied
noncentrality parameter exceeds the old grid's cap of 100 silently returns an
incorrect (capped or grid-noise) result. Both are fixed by the 2.0.0 engine,
so users on current R are better served by this submission than by the
version CRAN currently hosts.

## Test environments

* local: macOS (arm64), R 4.5.2

Before submission this will also be checked on win-builder (devel and release)
and R-hub.

## R CMD check results

0 errors | 0 warnings | NOTEs as described below.

* "Found the following (possibly) invalid URLs:
  https://yelleKneK.github.io/BUCSS ... Status: 404". This is the package's
  pkgdown documentation website, listed in DESCRIPTION (URL) and referenced by
  the pkgdown configuration. It is hosted on GitHub Pages and is deployed from
  the package repository. The site is published before this submission, so the
  URL resolves on the CRAN check machines; a local check run before the site is
  deployed will report the 404, which is expected and not a package defect.
* "checking HTML version of manual": the local copy of HTML Tidy is too old and
  'V8' is unavailable, so the optional HTML validation is skipped. Environmental;
  it does not appear on a CRAN check machine.
* "checking for future file timestamps ... unable to verify current time" may
  also appear when the local machine cannot reach a time server. Environmental;
  it does not appear on a CRAN check machine.

There are no other NOTEs. The package was previously on CRAN as version 1.2.1.

## Reverse dependencies

There are no reverse dependencies on CRAN.
