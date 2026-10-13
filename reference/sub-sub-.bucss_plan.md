# Positional Extraction From a Bias and Uncertainty Corrected Sample Size Result

Backward-compatible positional extraction for
[`bucss_plan`](https://yelleKneK.github.io/BUCSS/reference/print.bucss_plan.md)
objects. The deprecated 1.x functions (see
[`bucss-deprecated`](https://yelleKneK.github.io/BUCSS/reference/bucss-deprecated.md))
returned an unnamed two-element list, `list(sample_size, ncp)`, so code
written for BUCSS 1.x reads `result[[1]]` and `result[[2]]`. The
`ss_buc_*` functions return a richer object, but `result[[1]]` still
yields the necessary sample size and `result[[2]]` the bias and
uncertainty adjusted noncentrality parameter, so that 1.x extraction
code keeps working.

Any other subscript (a column name, or a length other than one) falls
through to the ordinary `data.frame` method, so `result[["value"]]` and
`result$value` are unaffected.

The legacy reading applies only to an object that still has the shape of
a single plan, as its planner returned it. An object that keeps the
class without that shape (several plans stacked together, or a plan with
rows removed, which some data manipulation tools leave in the class) is
read as the plain `data.frame` it is, so `[[1]]` and `[[2]]` give its
`term` and `value` columns rather than one plan's sample size and
noncentrality parameter.

## Usage

``` r
# S3 method for class 'bucss_plan'
x[[i, ...]]
```

## Arguments

- x:

  A `bucss_plan` object.

- i:

  A subscript. A single double `1` or `2`, as 1.x code writes it,
  selects the legacy scalar; anything else, including the integers `1L`
  and `2L`, is passed to the `data.frame` method as a column index.

- ...:

  Passed to the `data.frame` method.

## Value

For `i` equal to the double `1` or `2` and an object with the shape of a
single plan, a single numeric value; otherwise whatever the `data.frame`
`[[` method returns.
