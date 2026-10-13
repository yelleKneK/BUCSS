# Planning a Correlation Study

A Pearson correlation is one of the most commonly reported effects in
the social sciences, and one of the most commonly replicated. Suppose a
published study of 100 participants reported a correlation of *r* = .35,
and you want to plan a study powered to find that relationship again.
How large should it be? As with every planner in this package, the
reported correlation should not be taken at face value: it was selected
for publication in part because it was large, so planning from it
without correction tends to yield a study that is too small.
[`ss_buc_correlation()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_correlation.md)
applies the bias and uncertainty correction of Anderson, Kelley, and
Maxwell (2017), which reads the prior result as a draw from a
distribution truncated at the publication threshold (Taylor & Muller,
1996), and returns the total sample size the replication needs.

``` r

first_plan <- ss_buc_correlation(r_observed = .35, N = 100, alpha_prior = .05,
                                 alpha_planned = .05, assurance = .80,
                                 desired_power = .80)
first_plan
```

| term          | value |
|:--------------|:------|
| necessary_N   | 122   |
| actual_power  | 0.802 |
| ncp_adjusted  | 6.71  |
| rho_planning  | 0.251 |
| r_observed    | 0.35  |
| N             | 100   |
| alpha_prior   | 0.05  |
| alpha_planned | 0.05  |
| assurance     | 0.8   |
| desired_power | 0.8   |

Design: Pearson correlation; Sample size unit: total; Largest
supportable assurance: .99

The returned `necessary_N` is the total sample size and `ncp_adjusted`
the bias and uncertainty adjusted noncentrality parameter.
`rho_planning` is the planning effect size, the adjusted correlation the
plan is built on given the prior result, `alpha_prior`, and `assurance`
(not an estimate of the population correlation to report),
`sqrt(ncp_adjusted / (ncp_adjusted + N))` with the prior study’s *N*:
here .251, against the published .35. A study of 122 participants has
exactly the returned `actual_power` when the population correlation is
.251. For the meaning of `assurance` and `alpha_prior`, and how each
moves the recommended sample size, see
[`vignette("understanding-assurance", package = "BUCSS")`](https://yelleKneK.github.io/BUCSS/articles/understanding-assurance.md).
Briefly, raising assurance buys protection against an unlucky prior
estimate at the cost of a larger study:

``` r

assurances <- c(.5, .7, .8, .9)
sizes <- sapply(assurances, function(a) {
  plan <- ss_buc_correlation(r_observed = .35, N = 100, assurance = a)
  plan$value[plan$term == "necessary_N"]
})
data.frame(assurance = assurances, necessary_N = sizes)
#>   assurance necessary_N
#> 1       0.5          64
#> 2       0.7          92
#> 3       0.8         122
#> 4       0.9         205
```

## A Correlation Is Not a Regression Slope

The test of whether a correlation differs from zero is arithmetically
the test of the slope in a simple regression, so it is tempting to plan
a correlation study with the single-coefficient regression planner,
treating the predictor as fixed. Doing so gives a slightly different,
and slightly too small, answer:

``` r

t_equivalent <- .35 * sqrt(98 / (1 - .35^2))
by_correlation <- ss_buc_correlation(t_observed = t_equivalent, N = 100)
by_regression  <- ss_buc_reg_coef(t_observed = t_equivalent, N = 100, p = 1,
                                  random_predictors = FALSE)
c(correlation = by_correlation$value[by_correlation$term == "necessary_N"],
  regression  = by_regression$value[by_regression$term == "necessary_N"])
#> correlation  regression 
#>         122         119
```

The gap is small here, three participants, but it is systematic (only
near the assurance ceiling does it reverse, and then by a few percent,
which can be hundreds of participants because plans there are very
large), and the reason is worth understanding because it decides which
planner is correct for your design.

The regression planner, with `random_predictors = FALSE`, is exact when
the predictor is **fixed by design**: the researcher chooses the values
of the predictor, as in an experiment with assigned doses. A correlation
study is different. It **samples both variables**, and the sampled
spread of the predictor is itself random from study to study.
Conditional on the spread that happens to be drawn, the test statistic
is a noncentral *F*, but its noncentrality parameter rises and falls
with that spread, which follows a chi square distribution on *N* - 1
degrees of freedom. The distribution that actually governs the statistic
is therefore a mixture of noncentral *F* distributions, not a single
one, and a mixture is more spread out than its center. The two frames
give the same estimates, and the same distribution when there is no
effect, but different power when there is one (Sampson, 1974; Gatsonis &
Sampson, 1989; Kelley & Maxwell, 2003). Planning as though the statistic
were a plain noncentral *F* overstates the planned study’s power by
roughly 1/*n* for a plan of *n* participants: about .01 at *n* = 100,
.02 at 50, .05 at 20, and .10 at 10 (for the regression plan of 119
above, about .009). It therefore returns a sample size that is a few
participants short; the correct plan is typically 2 to 5 participants
larger, and 1 to 3 larger for the smallest plans.

[`ss_buc_correlation()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_correlation.md)
uses the mixture. It has an exact closed form (the Poisson weights that
build the noncentral *F* become negative binomial weights when mixed
over the chi square), so the returned sample size is exact for bivariate
normal data rather than an approximation. The rule of thumb is simple:

- If both variables are sampled, as in an observational correlation
  study, use
  [`ss_buc_correlation()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_correlation.md).
- If the predictor is fixed by the design, the effect is a regression
  slope, and
  [`ss_buc_reg_coef()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_coef.md)
  with `p = 1` and `random_predictors = FALSE` is the planner for it.
  With `random_predictors = TRUE` the same call returns exactly what
  [`ss_buc_correlation()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_correlation.md)
  returns, and the multiple regression planners use the same mixture for
  several sampled predictors.

You may supply the prior result either as the correlation `r_observed`
or as its *t* statistic `t_observed`; the two are equivalent.

## Assuming Less Publication Bias

`alpha_prior` is the significance level a prior result had to reach to
be published, so it encodes how strongly publication bias is assumed to
have shaped the prior literature. The default of .05 assumes the field
publishes results significant at the .05 level; a larger value assumes a
more permissive literature and therefore a smaller correction. A prior
study that escaped the file drawer entirely, such as a registered report
or a pilot, is described by `alpha_prior = 1`: no *p* value exceeds 1,
so that means no publication filter at all, and only the correction for
uncertainty, set by `assurance`, remains.

``` r

alphas <- c(.05, .20, 1)
sizes <- sapply(alphas, function(ap) {
  plan <- ss_buc_correlation(r_observed = .35, N = 100, alpha_prior = ap)
  plan$value[plan$term == "necessary_N"]
})
data.frame(alpha_prior = alphas, necessary_N = sizes)
#>   alpha_prior necessary_N
#> 1        0.05         122
#> 2        0.20         108
#> 3        1.00         103
```

## Auditing the Plan by Simulation

The plan above rests on a promise about the long run: at
`assurance = .80`, about 80% of studies planned this way should reach
the target power. That promise is easy to state and hard to picture.
[`ss_buc_sensitivity()`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_sensitivity.md)
makes it concrete by simulating the very process the correction models.
It draws many prior studies from a true effect you name, keeps only the
ones a literature would have published, re-plans from each survivor, and
reports how often those plans actually attain the target.

You must name the true effect; there is no default, because letting the
method assume the corrected estimate is true would be asking it to grade
its own work. A natural first choice is the plan’s own corrected
parameter, which asks “if the correction happened to be exactly right,
how would plans built this way behave?”

``` r

plan <- ss_buc_correlation(r_observed = .35, N = 100, assurance = .80)
ncp <- plan$value[plan$term == "ncp_adjusted"]
audit <- ss_buc_sensitivity(plan, true_ncp = ncp, replications = 2000,
                            seed = 2017)
audit
#>  term                          value  
#>  true_ncp                      6.71   
#>  ncp_adjusted                  6.71   
#>  replications                  2000   
#>  publication_rate              0.718  
#>  refusal_rate                  0.324  
#>  attainment                    0.708  
#>  attainment_mcse               0.0124 
#>  attainment_with_refusals      0.803  
#>  attainment_with_refusals_mcse 0.00889
#>  size_at_true_effect           122    
#>  size_q25                      110    
#>  size_median                   212    
#>  size_q75                      564    
#> 
#> Design: Pearson correlation
#> Sample size unit: total
#> Requested assurance: .80
#> Target power: .80
```

Two attainment rates are reported, and the difference between them is
the point of the exercise.

- **`attainment_with_refusals`** is the fraction of published prior
  studies for which a plan built this way reaches the target, counting a
  refusal (a prior result too weak to plan from, which the function
  declines rather than guesses at) as a success, because declining is
  the most cautious plan of all. This is the quantity the assurance
  guarantee is about, and it lands on the requested `assurance` of .80
  within simulation error: here it is .803, with a Monte Carlo standard
  error of .0089.
- **`attainment`** is the same fraction among only the studies that
  received a plan. It is lower, .708 here (Monte Carlo standard error
  .012), because it excludes the refusals, and it is what a researcher
  who did get a plan actually experiences.

The two are linked by the refusal rate reported alongside them:

`attainment_with_refusals = refusal_rate + (1 - refusal_rate) * attainment`

Here that is .324 + (1 - .324) times .708, or .803; the refusal rate’s
own Monte Carlo standard error is .01. The two rates coincide when
nothing is refused and separate as refusals rise. That separation is not
a flaw in the method; it is the visible price of a method that refuses
to plan from a prior result too weak to support the request, rather than
returning a confident but hollow number. `size_at_true_effect` is the
sample size the true effect really requires, the benchmark each
simulated plan is measured against, and the quartiles show how widely
the plans scatter around it.

The most instructive use of the function is to vary the assumed truth.
Ask what happens if the real correlation is smaller than the corrected
estimate, which is the case a cautious planner most wants to guard
against:

``` r

smaller <- ss_buc_sensitivity(plan, true_ncp = ncp * 0.6, replications = 2000,
                              seed = 2017)
rows <- c("attainment_with_refusals", "attainment", "refusal_rate",
          "size_at_true_effect")
data.frame(term = rows, value = smaller$value[match(rows, smaller$term)])
#>                       term      value
#> 1 attainment_with_refusals   0.803500
#> 2               attainment   0.637788
#> 3             refusal_rate   0.457500
#> 4      size_at_true_effect 200.000000
```

A smaller truth needs a larger study, 200 participants rather than 122,
yet `attainment_with_refusals` holds at .804 (Monte Carlo standard error
.0089): the guarantee does not depend on the corrected estimate being
right. What changes is how it is kept. More of the published prior
results are now too weak to plan from, so the refusal rate rises to .458
(Monte Carlo standard error .011), and among the plans actually issued
only .638 reach the target (Monte Carlo standard error .015).

## When the Prior Result Cannot Support the Request

Every prior result has a largest assurance it can support, the
closed-form value `1 - p/alpha_prior`, where *p* is the prior study’s
*p* value. Ask for more and the correction drives the noncentrality
parameter to zero: the prior effect cannot be distinguished from nothing
at that level of confidence, so no sample size attains the target and
the function stops rather than returning a meaningless answer. The
strong prior used throughout this vignette supports assurance up to
about .99, so this does not arise for it, but a weaker prior reaches its
ceiling quickly.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on any plan
draws the ceiling and the way the required sample size climbs toward it;
see
[`vignette("understanding-assurance", package = "BUCSS")`](https://yelleKneK.github.io/BUCSS/articles/understanding-assurance.md)
for that figure and a fuller discussion. A refusal is informative in its
own right: it says the prior study cannot anchor a replication aimed at
detecting the effect, and Anderson and Kelley (2024) describe the other
goals a replication can have, each with its own sample size logic.

## References

Anderson, S. F., & Kelley, K. (2024). Sample size planning for
replication studies: The devil is in the design. *Psychological Methods,
29*, 844–867. <https://doi.org/10.1037/met0000520>

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28*, 1547–1562. <https://doi.org/10.1177/0956797617723724>

Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact power
and sample size calculations. *Psychological Bulletin, 106*, 516–524.
<https://doi.org/10.1037/0033-2909.106.3.516>

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*, 305–321.
<https://doi.org/10.1037/1082-989X.8.3.305>

Sampson, A. R. (1974). A tale of two regressions. *Journal of the
American Statistical Association, 69*, 682–689.
<https://doi.org/10.1080/01621459.1974.10480189>

Taylor, D. J., & Muller, K. E. (1996). Bias in linear model power and
sample size calculation due to estimating noncentrality. *Communications
in Statistics: Theory and Methods, 25*, 1595–1610.
<https://doi.org/10.1080/03610929608831787>
