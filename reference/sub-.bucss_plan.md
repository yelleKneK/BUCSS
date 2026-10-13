# Subset a Bias and Uncertainty Corrected Sample Size Result

Row and column subsetting returns a plain `data.frame`. The class is
dropped first so that `[.data.frame`'s internal column extraction (which
uses `[[` on columns 1 and 2) cannot collide with the legacy positional
`[[` contract (see
[`[[.bucss_plan`](https://yelleKneK.github.io/BUCSS/reference/sub-sub-.bucss_plan.md));
without this, `head(result)` and `result[i, ]` would silently return the
legacy scalars in place of the `term` and `value` columns. A subset is
no longer the full planning result (its attributes and print contract no
longer apply), so the return is an ordinary `data.frame`.

## Usage

``` r
# S3 method for class 'bucss_plan'
x[...]

# S3 method for class 'bucss_plan'
x[i, j] <- value

# S3 method for class 'bucss_plan'
x[[i, j]] <- value

# S3 method for class 'bucss_plan'
x$name <- value

# S3 method for class 'bucss_plan'
names(x) <- value

# S3 method for class 'bucss_plan'
rbind(..., deparse.level = 1)
```

## Arguments

- x:

  A `bucss_plan` object.

- ...:

  Passed to the `data.frame` method.

- i, j, value:

  Passed to the `data.frame` method.

- deparse.level:

  Passed to the `data.frame` method.

## Value

A `data.frame` (or vector, under the usual `drop` rules).

## Details

Modifying or stacking a result does the same. Assignment through `[<-`,
`[[<-`, `$<-`, and `names<-` (and so `colnames<-` and
[`setNames()`](https://rdrr.io/r/stats/setNames.html)), and
[`rbind()`](https://rdrr.io/r/base/cbind.html), all return an ordinary
`data.frame` without the print,
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html), or
positional `[[` behavior. Assigning to `result[[1]]` or `result[[2]]` is
refused with an error, because reading those subscripts gives the legacy
scalars while assigning to them would overwrite the `term` or `value`
column.
