## Release summary

BUCSS 2.0.0 is a major update of a package already on CRAN (1.2.1). It
modernizes the user-facing API: the planners are renamed to the `ss_buc_*`
prefix with snake_case arguments, and each returns a tidy `bucss_power` object
instead of an unnamed list. The release is deliberately not backward compatible
with the 1.x API; the dotted 1.x function names are defunct, and their error
messages direct users to the CRAN archive for version 1.2.1.

With a single documented exception (the `df_denominator` argument of
`ss_buc_factorial_anova_general`, which 1.x accepted but silently ignored and
which now functions), the corrected noncentrality and sample size computations
are unchanged from 1.2.1 and are verified against it by the regression tests in
`tests/testthat/`.

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
