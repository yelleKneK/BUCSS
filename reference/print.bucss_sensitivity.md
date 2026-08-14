# Print a simulated sensitivity audit of a plan

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
ss_buc_sensitivity(plan, true_ncp = 2, replications = 200, seed = 2017)
#>  term                     value 
#>  true_ncp                 2     
#>  ncp_adjusted             0.975 
#>  replications             200   
#>  publication_rate         0.49  
#>  refusal_rate             0.505 
#>  attainment               0.616 
#>  attainment_mcse          0.0489
#>  attainment_with_refusals 0.81  
#>  size_at_true_effect      61    
#>  size_q25                 48.5  
#>  size_median              85    
#>  size_q75                 456   
#> 
#> Design: Dependent (paired) t test
#> Sample size unit: number of pairs
#> Requested assurance: .80
#> Target power: .80
```
