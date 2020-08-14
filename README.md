# BUCSS
BUCSS is an R package for implementing Bias- and Uncertainty-Corrected Sample Size planning. BUCSS implements a method of correcting for publication bias and uncertainty when planning sample sizes in a future study from an original study. See [Anderson, Kelley, &amp; Maxwell (2017; *Psychological Science*, *28*, 1547-1562)](https://www3.nd.edu/~kkelley/publications/articles/Anderson_Kelley_Maxwell_Psychological_Science_2017.pdf).

# using `BUCSS`
After initially installing BUCSS on your R system `install.packages("BUCSS")`, load the package.	
``` r	
library(BUCSS)
```

Consider an original study in which two independent groups are used to test the null hypothesis of no difference in population means of the two groups (e.g., treatment and control group). The independent samples $t$-test in which the original study reported $t=3.00$ based on sample size per group of $n_1=50$ and $n_2=55$ with a Type I error rate of $\alpha=.05.$. The desired study seeks to have a statistical power of .80 and to have 90\% assurance that the power will be at least 80\%. The desired assurance is the probability that a planned study done using this method will reach or surpass desired level of statistical power. Note that an assurance of .5 corrects for publication bias only, whereas assurance $> .5$ corrects for uncertainty. 

``` r	
ss.power.it(t.observed=3, n=c(50, 55), alpha.prior=.05, alpha.planned=.05, power=.80, assurance=.90)
```
Which yeilds a necessary sample size of $n_1=n_2=1482$ per group (i.e., $N=n_1+n_2=2964$). 
