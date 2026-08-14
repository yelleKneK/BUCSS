# Print a bias and uncertainty corrected sample size result

Prints the aligned `term`/`value` table (the planned-study results
followed by the rows echoing the planning inputs), then factual footer
lines naming the design, the unit the sample size is counted in, and the
largest assurance the prior result can support. Whole numbers print
clean; other values print at `getOption("bucss.digits", 3)` significant
figures. Only the display is rounded; the stored values keep full
precision.

## Usage

``` r
# S3 method for class 'bucss_power'
print(x, digits = getOption("bucss.digits", 3L), ...)
```

## Arguments

- x:

  An object of class `bucss_power`.

- digits:

  Significant figures for non-integer values; defaults to
  `getOption("bucss.digits", 3)`.

- ...:

  Further arguments, ignored.

## Value

`x`, invisibly.

## Examples

``` r
result <- ss_buc_independent_t(t_observed = 3, n = c(50, 55),
  assurance = .90)
result
#>  term                  value
#>  necessary_n_per_group 1485 
#>  total_N               2970 
#>  actual_power          0.803
#>  ncp_adjusted          0.526
#>  t_observed            3    
#>  n_1                   50   
#>  n_2                   55   
#>  alpha_prior           0.05 
#>  alpha_planned         0.05 
#>  assurance             0.9  
#>  desired_power         0.8  
#> 
#> Design: Independent t test
#> Sample size unit: per group
#> Largest supportable assurance: .93

# Display only: the stored values keep full precision.
print(result, digits = 6)
#>  term                  value   
#>  necessary_n_per_group 1485    
#>  total_N               2970    
#>  actual_power          0.802709
#>  ncp_adjusted          0.526243
#>  t_observed            3       
#>  n_1                   50      
#>  n_2                   55      
#>  alpha_prior           0.05    
#>  alpha_planned         0.05    
#>  assurance             0.9     
#>  desired_power         0.8     
#> 
#> Design: Independent t test
#> Sample size unit: per group
#> Largest supportable assurance: .93
result$value[result$term == "ncp_adjusted"]
#> [1] 0.5262426
```
