# Positional extraction from a bias and uncertainty corrected sample size result

Backward-compatible positional extraction for
[`bucss_power`](https://yelleKneK.github.io/BUCSS/reference/print.bucss_power.md)
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

## Usage

``` r
# S3 method for class 'bucss_power'
x[[i, ...]]
```

## Arguments

- x:

  A `bucss_power` object.

- i:

  A subscript. A single `1` or `2` selects the legacy scalar; anything
  else is passed to the `data.frame` method.

- ...:

  Passed to the `data.frame` method.

## Value

For `i` equal to `1` or `2`, a single numeric value; otherwise whatever
the `data.frame` `[[` method returns.
