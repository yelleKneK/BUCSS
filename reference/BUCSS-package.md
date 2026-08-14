# BUCSS: Bias and Uncertainty Corrected Sample Size

Implements sample size planning for a future study from the summary
statistics of a prior study, correcting the prior effect for publication
bias and uncertainty so the planned study attains its target statistical
power more reliably than methods that take a sample effect size at face
value. The correction builds on the truncated noncentral distribution of
Taylor and Muller (1996) and is described in Anderson, Kelley, and
Maxwell (2017).

## Details

The functions are also available, through a web interface, as Shiny apps
at <https://designingexperiments.com>.

BUCSS follows the design and naming conventions of the DMAR package
(Design, Measurement, and Analysis of Human-Centered Research): its
`ss_buc_*` planners parallel the `ss_aipe_*` (accuracy in parameter
estimation) and `ss_power_*` (power analysis) sample size planning
families in DMAR, returning a tidy result object in the same house
style. DMAR is a natural companion for effect size estimation,
confidence intervals, and sample size planning across additional
frameworks.

The bias and uncertainty adjusted noncentrality parameter is found by
solving the truncated-likelihood equation directly with
[`uniroot`](https://rdrr.io/r/stats/uniroot.html), so it is no longer
capped at 100 (the 1.x grid bound) and needs no resolution argument.

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

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

DMAR

## Author

Ken Kelley (<kkelley@nd.edu>), Samantha F. Anderson
(<samantha.f.anderson@asu.edu>), and Scott E. Maxwell
