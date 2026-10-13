# Print a Simulated Sensitivity Audit of a Plan

Prints the aligned `term`/`value` table, then factual footer lines
naming the design, the sample size unit, and the requested `assurance`
and `desired_power` the simulated rates should be read against.

## Usage

``` r
# S3 method for class 'bucss_sensitivity'
print(x, digits = getOption("bucss.digits", 3L), ...)
```

## Arguments

- x:

  An object of class `bucss_sensitivity`.

- digits:

  Significant figures for non-integer values; defaults to
  `getOption("bucss.digits", 3)`.

- ...:

  Further arguments, ignored.

## Value

`x`, invisibly.

## Examples

``` r
plan <- ss_buc_paired_t(t_observed = 3, N = 30, assurance = .80)
ss_buc_sensitivity(plan, true_ncp = 2, replications = 200, seed = 113)
#>  term                          value 
#>  true_ncp                      2     
#>  ncp_adjusted                  0.975 
#>  replications                  200   
#>  publication_rate              0.49  
#>  refusal_rate                  0.51  
#>  attainment                    0.541 
#>  attainment_mcse               0.0503
#>  attainment_with_refusals      0.775 
#>  attainment_with_refusals_mcse 0.0295
#>  size_at_true_effect           61    
#>  size_q25                      45    
#>  size_median                   70.5  
#>  size_q75                      187   
#> 
#> Design: Dependent (paired) t test
#> Sample size unit: number of pairs
#> Requested assurance: .80
#> Target power: .80
```
