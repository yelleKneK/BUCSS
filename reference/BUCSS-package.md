# BUCSS: Bias and Uncertainty Corrected Sample Size

Implements sample size planning for a future study from the summary
statistics of a prior study, correcting the prior effect for publication
bias and uncertainty so the planned study attains its target statistical
power more reliably than methods that take a sample effect size at face
value. The correction builds on the truncated noncentral distribution of
Taylor and Muller (1996) and is described in Anderson, Kelley, and
Maxwell (2017).

## Details

Each `ss_buc_*` function plans one analysis. The *t* tests:
[`ss_buc_one_sample_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_one_sample_t.md),
[`ss_buc_paired_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_paired_t.md),
[`ss_buc_independent_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md),
and
[`ss_buc_welch_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_welch_t.md).
Between-subjects analysis of variance and covariance:
[`ss_buc_one_way_anova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_one_way_anova.md),
[`ss_buc_factorial_anova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_factorial_anova.md),
[`ss_buc_factorial_anova_general`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_factorial_anova_general.md),
and
[`ss_buc_ancova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_ancova.md).
Within-subjects and split-plot analysis of variance:
[`ss_buc_rm_anova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_rm_anova.md),
[`ss_buc_rm_anova_general`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_rm_anova_general.md),
[`ss_buc_mixed_anova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_mixed_anova.md),
and
[`ss_buc_mixed_anova_general`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_mixed_anova_general.md).
Regression and correlation:
[`ss_buc_reg_coef`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_coef.md),
[`ss_buc_R2`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_R2.md),
[`ss_buc_reg_joint`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_joint.md),
and
[`ss_buc_correlation`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_correlation.md).
And
[`ss_buc_manova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_manova.md)
for a two-group multivariate comparison and
[`ss_buc_chisq_diff`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_chisq_diff.md)
for the chi square difference between nested models. The experimental
designs among these, their tests, and the statistical power of those
tests are described in Maxwell, Delaney, and Kelley (2027).
[`ss_buc_sensitivity`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_sensitivity.md)
simulates the literature a plan came from and reports how often such a
plan attains its target power, and
[`plot`](https://yelleKneK.github.io/BUCSS/reference/plot.bucss_plan.md)
draws the largest assurance a prior result can support.

The planners carried over from BUCSS 1.x are also available, through a
web interface, as Shiny apps at <https://designingexperiments.com>; the
designs added in 2.0.0 are available in R only.

BUCSS follows the design and naming conventions of the DMAR package
(Design, Measurement, and Analysis in R): its `ss_buc_*` planners
parallel the `ss_aipe_*` (accuracy in parameter estimation) and
`ss_power_*` (power analysis) sample size planning families in DMAR,
returning a tidy result object in the same house style. DMAR is a
natural companion for effect size estimation, confidence intervals, and
sample size planning across additional frameworks.

The bias and uncertainty adjusted noncentrality parameter is found by
solving the truncated-likelihood equation directly with
[`uniroot`](https://rdrr.io/r/stats/uniroot.html), so it is no longer
capped at 100 (the 1.x grid bound) and needs no resolution argument.

**How assurance is used.** Every planner takes an `assurance`: the long
run proportion of planned studies whose power at the true effect reaches
or surpasses the desired power, when each is planned from a prior study
published under the rule `alpha_prior` describes. A planned study
reaches its target power whenever the noncentrality parameter it was
planned from is at or below the true one, because power rises with the
effect (and, apart from rounding the sample size up to a whole number
and from the smallest sample size a planner returns, the design's own
minimum or the Welch planner's 10 per group, only then). The planners
therefore plan from the lower limit of a one-sided confidence interval
for the prior study's noncentrality parameter, at confidence level
`assurance`, computed from the distribution of the prior statistic
truncated at the publication threshold: the value at which a published
prior study would exceed the observed statistic with probability 1 -
`assurance` (what Anderson, Kelley, and Maxwell, 2017, call the (1 -
`assurance`) quantile of the likelihood distribution). At
`assurance = .5` this is the median unbiased estimate, which corrects
for publication bias only; above .5 it adds the correction for
uncertainty. The term follows Kelley and Maxwell (2003), who used
assurance in sample size planning for the probability that a confidence
interval will be no wider than desired; here it is the probability that
the planned study's power will be at least the desired power.
[`ss_buc_sensitivity`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_sensitivity.md)
checks the proportion by simulation for any plan, counting a refusal to
plan as meeting the target, as the guarantee does, and also reporting
the rate among the plans issued.

**The planning effect size.** Every plan also reports a planning effect
size, in a row whose name ends in `_planning` (`smd_planning`,
`f_planning`, `rho_planning`, and so on). It is an adjusted effect size:
the value, in the design's usual effect size metric, used to plan the
sample size given the situation at hand (the prior study's result, the
publication bias assumed through `alpha_prior`, and the uncertainty
allowed for through `assurance`). It is `ncp_adjusted` restated in that
metric and is what the plan is built on, not an estimate of the
population effect to report: above an assurance of .5 it is deliberately
conservative (a lower confidence limit), and at .5 it is the median
unbiased estimate. Each planner's Value section names its metric.

**Safeguard power is a special case.** Safeguard power (Perugini,
Gallucci, & Costantini, 2014) is a way of planning a new study from a
prior study's effect size that guards against the prior estimate being
too large. Instead of planning for the effect size the prior study
observed, the researcher plans for the lower limit of a two-sided 60
percent confidence interval around it (a limit that, in repeated
sampling, falls below the true effect 80 percent of the time), so the
new study keeps its intended power even if the true effect is somewhat
smaller than the one observed. It allows for uncertainty in the prior
estimate but not for publication bias.

It is the special case `alpha_prior = 1` and `assurance = .80` of the
method here, that is, BUCSS without the publication bias correction, at
an assurance of .80. `alpha_prior` is the significance level a prior
result had to reach to be published, and no *p* value exceeds 1, so
`alpha_prior = 1`, which may also be written
`alpha_prior = "No Publication Filter"`, means no publication filter:
any result would have been published, the truncation is removed, and the
prior statistic has its ordinary distribution. (Anderson, Kelley, and
Maxwell, 2017, p. 1552, note that this level can be set equal to 1 when
publication bias is thought not to apply.) At `assurance = .80` the
adjusted noncentrality parameter is then the lower limit of a one-sided
80 percent confidence interval, which is the lower limit of the
two-sided 60 percent interval (either way, 20 percent lies below it),
and the plan is the smallest sample size whose ordinary power at that
value reaches `desired_power`. In general, the lower limit of a
two-sided interval at confidence level *C* is `assurance = (1 + C)/2`,
so a two-sided 80 percent interval is `assurance = .90`.

For an *F* or chi square statistic the correspondence is exact: the
adjusted noncentrality parameter is the lower limit of that interval,
computed from the statistic's ordinary noncentral distribution. Two
differences from safeguard power remain, and each can only make the plan
here the same as the safeguard plan or larger. First, for the *t* tests
and every other effect on one degree of freedom, safeguard power
computes its interval from the signed *t*, while the interval here
bounds the size of the effect whatever its sign, because the publication
rule is two-sided. That difference is negligible for a convincing prior
result and grows as the prior result weakens. Second, an unbalanced
prior design is read as a balanced one, and the Welch planner returns at
least 10 per group.
[`vignette("understanding-assurance", package = "BUCSS")`](https://yelleKneK.github.io/BUCSS/articles/understanding-assurance.md)
measures both. For the regression planners the equivalence holds in the
frame `random_predictors` names. With `alpha_prior` below 1 the same
construction is applied to the publication-truncated distribution, which
is what the bias correction adds to safeguard power.

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

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8,* 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Perugini, M., Gallucci, M., & Costantini, G. (2014). Safeguard power as
a protection against imprecise power estimates. *Perspectives on
Psychological Science, 9,* 319–332.
[doi:10.1177/1745691614528519](https://doi.org/10.1177/1745691614528519)

Taylor, D. J., & Muller, K. E. (1996). Bias in linear model power and
sample size calculation due to estimating noncentrality. *Communications
in Statistics: Theory and Methods, 25,* 1595–1610.
[doi:10.1080/03610929608831787](https://doi.org/10.1080/03610929608831787)

## See also

DMAR

## Author

Ken Kelley (<kkelley@nd.edu>), Samantha F. Anderson
(<samantha.f.anderson@asu.edu>), and Scott E. Maxwell
