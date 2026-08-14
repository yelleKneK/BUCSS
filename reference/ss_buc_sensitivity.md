# Simulate the literature a plan came from and audit how often it works

`ss_buc_sensitivity` takes a plan returned by any `ss_buc_*` function
and asks the question the method's assurance claim is really about: if
the true effect were a particular value, and a literature published only
the prior studies that reached significance, how often would a plan
built this way actually reach the target power?

## Usage

``` r
ss_buc_sensitivity(object, true_ncp, replications = 1000L, seed = NULL)
```

## Arguments

- object:

  A plan returned by any `ss_buc_*` function.

- true_ncp:

  The assumed true noncentrality parameter of the prior study's design,
  on the same scale as the `ncp_adjusted` row of `object`. Required;
  there is deliberately no default.

- replications:

  Number of published prior studies to simulate. The default of 1000
  gives a Monte Carlo standard error near .013 on the attainment rate;
  10000 gives about .004.

- seed:

  Optional random number seed. `NULL` (the default) leaves the random
  number stream alone. When supplied, the caller's random number state
  is restored on exit.

## Value

An object of class `bucss_sensitivity`, a `data.frame` with columns
`term` and `value` carrying the assumed true parameter and the plan's
own `ncp_adjusted` for comparison, the number of replications, the
publication and refusal rates, the two attainment rates and the Monte
Carlo standard error of the first, the sample size the true effect
requires, and the quartiles of the sample sizes the simulated plans
called for. The design, the sample size unit, the requested `assurance`,
and the requested `desired_power` travel as attributes.

## Details

The bias and uncertainty correction rests on a claim about the long run.
At `assurance = .80`, eight replications in ten should reach the target
power. That claim is a statement about a whole literature, and a user
planning one study never sees it exercised.

This function exercises it on the user's own numbers. It reads the
design and the planning inputs off the plan you pass it, simulates prior
studies of that design from a true effect *you* name, discards the ones
a literature would not have published (those failing to reach
significance at `alpha_prior`), runs the same planner on each survivor,
and reports what happened.

**You must name the true effect.** There is no default, deliberately.
Defaulting to the corrected estimate would be convenient and circular:
it would ask the method to grade its own answer. `true_ncp` is on the
same scale as the `ncp_adjusted` the plan reports, so the plan's own
value is the natural anchor. Run the function at that value, at a
fraction of it, and at a multiple of it; the interesting reading is how
fast the attainment rate falls when the truth is smaller than the
corrected estimate.

**Reading the output.** A plan attains the target when its sample size
is at least the sample size the true effect actually requires, which is
reported as `size_at_true_effect`. Two attainment rates are given, and
the difference between them is not a rounding detail:

- `attainment` counts only the prior studies for which the method issued
  a plan. This is what a user experiences.

- `attainment_with_refusals` counts a refusal (the zero noncentrality
  parameter stop) as attaining, because declining to plan is the
  limiting case of an arbitrarily conservative plan. This is the
  quantity the assurance guarantee is about, and it should land on
  `assurance`.

When the refusal rate is high the first number sits well below
`assurance`, and that is not a defect: it is the visible price of the
method refusing rather than guessing. Its expected value is
`(assurance - refusal_rate)/(1 - refusal_rate)`.

**What this can and cannot detect.** The prior studies are drawn from
the sampling distribution the method assumes for that design, so this is
a check that the correction and the planned sample size are internally
right, and a demonstration of what publication bias does. It cannot
detect a violated distributional assumption, since the simulation shares
that assumption. For
[`ss_buc_welch_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_welch_t.md),
in particular, it holds the variance ratio fixed at the value you
supplied.

## References

Anderson, S. F., & Kelley, K. (2024). Sample size planning for
replication studies: The devil is in the design. *Psychological Methods,
29,* 844–867.
[doi:10.1037/met0000520](https://doi.org/10.1037/met0000520)

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28,* 1547–1562.
[doi:10.1177/0956797617723724](https://doi.org/10.1177/0956797617723724)

Anderson, S. F., & Maxwell, S. E. (2017). Addressing the 'replication
crisis': Using original studies to design replication studies with
appropriate statistical power. *Multivariate Behavioral Research, 52,*
305–324.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Taylor, D. J., & Muller, K. E. (1996). Bias in linear model power and
sample size calculation due to estimating noncentrality. *Communications
in Statistics: Theory and Methods, 25,* 1595–1610.

## See also

[`ss_buc_independent_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md)
and the other `ss_buc_*` planners, whose results this function takes as
input.

## Author

Ken Kelley (<kkelley@nd.edu>)

## Examples

``` r
plan <- ss_buc_independent_t(t_observed = 3, n = c(50, 55), assurance = .80)
plan$value[plan$term == "ncp_adjusted"]
#> [1] 1.338019

# Take the corrected estimate at face value as the truth.
ss_buc_sensitivity(plan, true_ncp = 1.9, replications = 200, seed = 2017)
#>  term                     value 
#>  true_ncp                 1.9   
#>  ncp_adjusted             1.34  
#>  replications             200   
#>  publication_rate         0.469 
#>  refusal_rate             0.535 
#>  attainment               0.602 
#>  attainment_mcse          0.0508
#>  attainment_with_refusals 0.815 
#>  size_at_true_effect      116   
#>  size_q25                 89    
#>  size_median              154   
#>  size_q75                 459   
#> 
#> Design: Independent t test
#> Sample size unit: per group
#> Requested assurance: .80
#> Target power: .80

# And ask what happens if the truth is smaller than the correction assumed.
ss_buc_sensitivity(plan, true_ncp = 1.2, replications = 200, seed = 2017)
#>  term                     value 
#>  true_ncp                 1.2   
#>  ncp_adjusted             1.34  
#>  replications             200   
#>  publication_rate         0.221 
#>  refusal_rate             0.605 
#>  attainment               0.519 
#>  attainment_mcse          0.0562
#>  attainment_with_refusals 0.81  
#>  size_at_true_effect      289   
#>  size_q25                 124   
#>  size_median              316   
#>  size_q75                 1370  
#> 
#> Design: Independent t test
#> Sample size unit: per group
#> Requested assurance: .80
#> Target power: .80
```
