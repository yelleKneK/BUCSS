# The planning sentence for a manuscript

Turns a `bucss_power` result into the sentence an author writes in a
manuscript's planning or method section: the necessary sample size in
its design's unit, the desired power, the assurance, and a clause noting
the publication bias and uncertainty correction. This mirrors the
`results_sentence()` helper in the sibling package `DMAR`.

## Usage

``` r
planning_sentence(x)
```

## Arguments

- x:

  An object of class `bucss_power`, as returned by any `ss_buc_*`
  function.

## Value

A single character string.

## Author

Ken Kelley (<kkelley@nd.edu>)

## Examples

``` r
result <- ss_buc_independent_t(t_observed = 3, n = c(50, 55),
  assurance = .90)
planning_sentence(result)
#> [1] "A planned study with 1485 participants per group (total N = 2970) attains 80% power with 90% assurance, after correcting the prior study's result for publication bias and uncertainty."
```
