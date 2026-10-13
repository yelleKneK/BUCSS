# The Planning Sentence for a Manuscript

Turns a `bucss_plan` result into the sentence an author writes in a
manuscript's planning or method section: the necessary sample size in
its design's unit, the desired power, the assurance, and a clause
stating which corrections the plan made. The clause is built from the
plan's own settings: the publication bias correction is named, with its
assumed threshold, only when `alpha_prior` is below 1, and the
uncertainty correction only when `assurance` is above .5 (at .5 the plan
uses the median estimate). A within-subjects or split-plot plan made
with `epsilon` below 1 adds the departure from sphericity it allowed
for. A second sentence names the planning effect size (the result's row
whose name ends in `_planning`), the adjusted effect size the plan is
built on, in the design's usual metric: a standardized mean difference,
Cohen's *f*, a correlation, and so on, as the planner's help page
describes. It is the value used to plan given the prior result,
`alpha_prior`, and `assurance`, not an estimate of the effect to report.
It is shown at `getOption("bucss.digits", 3)` significant figures (more
when fewer would show a value below 1 as 1), and without a leading zero
when the metric cannot exceed 1. This mirrors the `results_sentence()`
helper in the sibling package `DMAR`.

## Usage

``` r
planning_sentence(x)
```

## Arguments

- x:

  An object of class `bucss_plan`, as returned by any `ss_buc_*`
  function.

## Value

A single character string holding the two sentences.

## Author

Ken Kelley (<kkelley@nd.edu>) and Samantha F. Anderson
(<samantha.f.anderson@asu.edu>)

## Examples

``` r
result <- ss_buc_independent_t(t_observed = 3, n = c(50, 55),
  assurance = .90)
planning_sentence(result)
#> [1] "A planned study with 1,485 participants per group (total N = 2,970) attains 80% power with 90% assurance, after correcting the prior study's result for publication bias (assuming its literature published results significant at the .05 level) and for uncertainty. The plan is built on a planning effect size of 0.103, expressed as a standardized mean difference."

# With no publication bias assumed (a pilot study, say), the sentence
# names the uncertainty correction alone.
planning_sentence(ss_buc_independent_t(t_observed = 3, n = 20,
  alpha_prior = 1))
#> [1] "A planned study with 37 participants per group (total N = 74) attains 80% power with 80% assurance, after correcting the prior study's result for uncertainty, with no publication bias assumed. The plan is built on a planning effect size of 0.661, expressed as a standardized mean difference."
```
