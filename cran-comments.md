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

## Test environments

* local: macOS, R 4.5.2

(Before submission, also check on win-builder (devel and release) and R-hub.)

## R CMD check results

0 errors | 0 warnings | 2 notes

Both NOTEs are environmental rather than package issues and do not appear on a
clean check machine:

* "checking for future file timestamps ... unable to verify current time": the
  local machine could not reach a time server.
* "checking HTML version of manual": the system copy of HTML Tidy is too old and
  'V8' is unavailable, so the optional HTML validation and math-rendering checks
  are skipped.

## Reverse dependencies

There are no reverse dependencies on CRAN.
