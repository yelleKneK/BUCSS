# Deprecated 1.x Functions in BUCSS

The dot-named functions from BUCSS 1.x are **deprecated** as of BUCSS
2.0.0 and are no longer actively developed. They are retained so that
scripts written for BUCSS 1.x continue to run. Each one translates its
arguments to the 2.x convention and calls the corresponding `ss_buc_*`
function, returning that function's
[`bucss_plan`](https://yelleKneK.github.io/BUCSS/reference/print.bucss_plan.md)
object.

For backward compatibility the returned object still supports the 1.x
positional extraction: `result[[1]]` is the necessary sample size and
`result[[2]]` is the bias and uncertainty adjusted noncentrality
parameter, the same two values the 1.x functions returned as an unnamed
list (this is provided by a `[[` method for `bucss_plan`). Only those
two subscripts are kept: code that treats the 1.x result as a list in
other ways (`unlist(result)`, `result[[1L]]`, or printing it through
[`rbind()`](https://rdrr.io/r/base/cbind.html)) now sees the data frame.
The richer named view (`result$value`, `tidy(result)`, `glance(result)`)
is also available.

Calling one of these functions issues a deprecation warning that names
its replacement, once per function per session; the `ss_buc_*` functions
themselves never issue this deprecation warning. New code should call
the `ss_buc_*` functions directly.

## Usage

``` r
ss.power.it(
  t.observed,
  n,
  N,
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)

ss.power.dt(
  t.observed,
  N,
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)

ss.power.ba(
  F.observed,
  N,
  levels.A,
  levels.B = NULL,
  effect = c("factor.A", "factor.B", "interaction"),
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)

ss.power.ba.general(
  F.observed,
  N,
  cells,
  df.numerator,
  df.denominator,
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)

ss.power.wa(
  F.observed,
  N,
  levels.A,
  levels.B = NULL,
  effect = c("factor.A", "factor.B", "interaction"),
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)

ss.power.wa.general(
  F.observed,
  N,
  df.numerator,
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)

ss.power.spa(
  F.observed,
  N,
  levels.between,
  levels.within,
  effect = c("between", "within", "interaction"),
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)

ss.power.spa.general(
  F.observed,
  N,
  df.numerator,
  num.groups,
  effect = c("between.only", "within.only", "between.within"),
  df.num.within,
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)

ss.power.reg1(
  t.observed,
  N,
  p,
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)

ss.power.reg.all(
  F.observed,
  N,
  p,
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)

ss.power.reg.joint(
  F.observed,
  N,
  p,
  p.joint,
  alpha.prior = 0.05,
  alpha.planned = 0.05,
  assurance = 0.8,
  power = 0.8,
  step = 0.001
)
```

## Arguments

- t.observed, F.observed:

  The observed *t* or *F* statistic from the prior study (the 1.x
  spelling of `t_observed` / `F_observed`).

- n:

  Per-group sample size(s) of the prior study, for `ss.power.it`.

- N:

  Total sample size of the prior study (number of pairs for
  `ss.power.dt`).

- alpha.prior, alpha.planned:

  The prior-study and planned-study Type I error rates (the 1.x spelling
  of `alpha_prior` / `alpha_planned`). `alpha.prior` also accepts
  `"No Publication Filter"`, which means 1.

- assurance:

  The desired assurance.

- power:

  The desired statistical power of the planned study, which must exceed
  `alpha.planned`.

- step:

  Ignored. Accepted only so that 1.x calls that set it do not error.

- levels.A, levels.B:

  Number of levels of the first and (optional) second factor, for
  `ss.power.ba` and `ss.power.wa`.

- effect:

  The effect to plan for, in the 1.x dotted spelling (translated
  internally to the 2.x form).

- cells:

  Number of cells (groups) in the between-subjects design, for
  `ss.power.ba.general`.

- df.numerator, df.denominator:

  Numerator and denominator degrees of freedom for the effect of
  interest in the prior study, for the `.general` planners.

- levels.between, levels.within:

  Number of levels of the between- and within-subjects factors, for
  `ss.power.spa`.

- num.groups:

  Number of between-subjects groups, for `ss.power.spa.general`.

- df.num.within:

  Numerator degrees of freedom of the within-subjects component, for
  `ss.power.spa.general` with `effect = "between.within"`.

- p, p.joint:

  Number of predictors, and the number tested jointly, for the
  regression planners.

## Value

A
[`bucss_plan`](https://yelleKneK.github.io/BUCSS/reference/print.bucss_plan.md)
object, as returned by the corresponding `ss_buc_*` function.

## Details

The replacements are:

|  |  |
|----|----|
| **Deprecated** | **Replacement** |
| `ss.power.it` | [`ss_buc_independent_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md) (alias [`ss_buc_smd`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_independent_t.md)) |
| `ss.power.dt` | [`ss_buc_paired_t`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_paired_t.md) (alias [`ss_buc_smd_paired`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_paired_t.md)) |
| `ss.power.ba` | [`ss_buc_one_way_anova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_one_way_anova.md) (when `levels.B` is `NULL`) or [`ss_buc_factorial_anova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_factorial_anova.md) |
| `ss.power.ba.general` | [`ss_buc_factorial_anova_general`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_factorial_anova_general.md) |
| `ss.power.wa` | [`ss_buc_rm_anova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_rm_anova.md) |
| `ss.power.wa.general` | [`ss_buc_rm_anova_general`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_rm_anova_general.md) |
| `ss.power.spa` | [`ss_buc_mixed_anova`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_mixed_anova.md) |
| `ss.power.spa.general` | [`ss_buc_mixed_anova_general`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_mixed_anova_general.md) |
| `ss.power.reg1` | [`ss_buc_reg_coef`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_coef.md) |
| `ss.power.reg.all` | [`ss_buc_R2`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_R2.md) |
| `ss.power.reg.joint` | [`ss_buc_reg_joint`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_joint.md) |

Because BUCSS 2.0.0 finds the noncentrality parameter by direct root
finding rather than the fixed 0 to 100 grid used in 1.x, the values
these wrappers return can differ slightly from 1.x and are no longer
capped when the implied parameter exceeds 100. The `step` argument is
accepted for call compatibility and ignored.

Some 1.x behaviors deliberately changed, so a 1.x call can now error
where it used to return a number. `ss.power.ba.general` rejects a
`df.denominator` greater than `N - cells` (1.x accepted such a value and
ignored the argument entirely, including in its own documented example).
As in 1.x, an `assurance` or `power` of 1 or more is read as a
percentage, so exactly 1 means 1 percent; after that reading, an
`assurance` or `power` of 0, or of 100 (which means 1), is rejected
rather than returning an artifact of the internal search. A `power` at
or below `alpha.planned` after that reading is not allowable and stops
too, so `power = 1` (which means .01) stops unless `alpha.planned` is
below .01; 1.x returned the smallest sample size it considers for such a
target, with no notice. An `alpha.prior` of 1 is computed as exactly 1,
where 1.x computed with .999 in its place: most such plans are
unchanged, and the rest are slightly smaller. `ss.power.ba` and
`ss.power.wa` refuse `effect = "factor.B"` or `"interaction"` when
`levels.B` is not given (1.x refused only the interaction, and silently
returned the Factor A plan for `"factor.B"`). `ss.power.ba.general` and
`ss.power.spa.general` reject numerator degrees of freedom no design of
their kind can have (`df.numerator` above `cells - 1`; for the split
plot, above `num.groups - 1` for the between-subjects effect, or not a
multiple of `df.num.within` within that bound for the interaction).
Non-scalar, non-finite, and fractional inputs are rejected with a
message naming the argument. Most errors raised through these wrappers
come from the `ss_buc_*` planners and name their arguments (`t_observed`
for `t.observed`, and so on); an unrecognized `effect` is rejected by
the wrapper itself. The regression wrappers keep the 1.x fixed-predictor
computation (`random_predictors = FALSE`); call
[`ss_buc_reg_coef`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_coef.md),
[`ss_buc_R2`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_R2.md),
or
[`ss_buc_reg_joint`](https://yelleKneK.github.io/BUCSS/reference/ss_buc_reg_joint.md)
directly for predictors sampled along with the outcome. See `NEWS` for
the full accounting.

## See also

The `ss_buc_*` replacements listed above.
