# Deprecated 1.x functions in BUCSS

The dot-named functions from BUCSS 1.x are **deprecated** as of BUCSS
2.0.0 and are no longer actively developed. They are retained so that
scripts written for BUCSS 1.x continue to run unchanged. Each one
translates its arguments to the 2.x convention and calls the
corresponding `ss_buc_*` function, returning that function's
[`bucss_power`](https://yelleKneK.github.io/BUCSS/reference/print.bucss_power.md)
object.

For backward compatibility the returned object still supports the 1.x
positional extraction: `result[[1]]` is the necessary sample size and
`result[[2]]` is the bias and uncertainty adjusted noncentrality
parameter, the same two values the 1.x functions returned as an unnamed
list (this is provided by a `[[` method for `bucss_power`). The richer
named view (`result$value`, `tidy(result)`, `glance(result)`) is also
available.

Calling one of these functions issues a deprecation warning that names
its replacement, once per function per session; the `ss_buc_*` functions
themselves never warn. New code should call the `ss_buc_*` functions
directly.

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
  of `alpha_prior` / `alpha_planned`).

- assurance:

  The desired assurance.

- power:

  The desired statistical power of the planned study.

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
[`bucss_power`](https://yelleKneK.github.io/BUCSS/reference/print.bucss_power.md)
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

Three further 1.x behaviors deliberately changed, so a 1.x call can now
error where it used to return a number. `ss.power.ba.general` rejects a
`df.denominator` greater than `N - cells` (1.x accepted such a value and
ignored the argument entirely, including in its own documented example).
An `assurance` or `power` of exactly 0 or 1, and a `power` of 100, are
rejected rather than returning an artifact of the internal search.
Non-scalar, non-finite, and fractional inputs are rejected with a
message naming the argument. See `NEWS` for the full accounting.

## See also

The `ss_buc_*` replacements listed above.
