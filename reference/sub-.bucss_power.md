# Subset a bias and uncertainty corrected sample size result

Row and column subsetting returns a plain `data.frame`. The class is
dropped first so that `[.data.frame`'s internal column extraction (which
uses `[[` on columns 1 and 2) cannot collide with the legacy positional
`[[` contract above; without this, `head(result)` and `result[i, ]`
would silently return the legacy scalars in place of the `term` and
`value` columns. A subset is no longer the full planning result (its
attributes and print contract no longer apply), so the honest return is
an ordinary `data.frame`.

## Usage

``` r
# S3 method for class 'bucss_power'
x[...]

# S3 method for class 'bucss_power'
x[i, j] <- value

# S3 method for class 'bucss_power'
x[[i, j]] <- value

# S3 method for class 'bucss_power'
x$name <- value

# S3 method for class 'bucss_power'
rbind(..., deparse.level = 1)
```

## Arguments

- x:

  A `bucss_power` object.

- ...:

  Passed to the `data.frame` method.

- i, j, value:

  Passed to the `data.frame` method.

- deparse.level:

  Passed to the `data.frame` method.

## Value

A `data.frame` (or vector, under the usual `drop` rules).
