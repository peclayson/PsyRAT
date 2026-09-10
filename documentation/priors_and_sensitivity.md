# Priors and prior sensitivity (test-retest occasion effects)

Methodology note for the variance-component priors used in CmdStan estimation,
the SCI-004 change to the test-retest occasion priors, and a sensitivity
analysis quantifying its effect. This note may be folded into future user
documentation.

## Measurement-scale assumption

The default priors are weakly informative **on the ERP microvolt scale**. They
are not scale-free. Scores on a very different scale (for example standardized
values, or other modalities) can make these priors informative. Until per-design
priors are user-configurable (SCI-003) and documented for other modalities,
treat the defaults as appropriate for ERP amplitude data.

## SCI-004: widened test-retest occasion priors

The crossed test-retest model (Persons x Trials x Occasions) estimates seven
variance components. Two of the occasion-related standard deviations previously
carried near-zero half-Cauchy priors that forced them toward zero:

| Component | Meaning | Old prior | New prior |
| --- | --- | --- | --- |
| `sig_occ` | occasion main effect (sigma_o) | `cauchy(0, 0.05)` | `cauchy(0, 5)` |
| `sig_trlxocc` | trial x occasion (sigma_oi) | `cauchy(0, 0.05)` | `cauchy(0, 1)` |

The old priors effectively encoded the assumption that occasion-to-occasion
variance (beyond person x occasion) was negligible. For **absolute-decision**
(dependability) test-retest coefficients, the occasion main effect enters the
absolute error term, so pinning it near zero removed a real source of
measurement error and inflated those coefficients. The widened priors let the
data inform these components.

## Sensitivity analysis

Estimated on the `trt_time` integration scenario (derived from a development
fixture that is not distributed with the release), 4 chains, 1000 warmup,
1000 sampling, seed 12345.
Coefficients evaluated at 30 trials and 2 occasions. This is an illustrative
single dataset; the magnitude of the shift depends on how much real occasion
variance is present.

Posterior means of the affected components moved from near zero to substantial,
while the other components changed only modestly:

| Component | Old mean | New mean |
| --- | --- | --- |
| `sig_occ` | 0.16 | 3.24 |
| `sig_trlxocc` | 0.20 | 1.37 |
| `sig_id` | 2.16 | 2.26 |
| `sig_trl` | 1.65 | 1.53 |
| `sig_occxid` | 1.50 | 1.33 |
| `sig_err` | 4.69 | 4.46 |
| `sig_trlxid` | 1.02 | 1.08 |

Reliability coefficients (point estimate; 95% credible interval), old vs new:

| Coefficient | Dependability old | Dependability new | Generalizability old | Generalizability new |
| --- | --- | --- | --- | --- |
| CE (equivalence) | 0.827 [0.278, 0.992] | 0.824 [0.270, 0.993] | 0.854 [0.323, 0.994] | 0.859 [0.344, 0.995] |
| CS (stability) | 0.582 [0.019, 0.984] | **0.415** [0.004, 0.971] | 0.589 [0.020, 0.985] | 0.627 [0.026, 0.989] |
| CES (trial equiv. + stability) | 0.550 [0.005, 0.979] | **0.399** [0.001, 0.966] | 0.571 [0.005, 0.984] | 0.606 [0.008, 0.989] |

Interpretation:

- The **dependability** coefficients that include the occasion main effect in
  absolute error (CS, CES) dropped substantially (about -0.15 to -0.17). CE
  dependability, which does not include the occasion main effect, was essentially
  unchanged.
- The **generalizability** coefficients, which use relative error and exclude the
  occasion main effect, were essentially unchanged (small increases attributable
  to modest posterior shifts in the person-interaction components).

The net effect is that absolute-decision test-retest reliability is no longer
inflated by an artificially suppressed occasion variance. Generalizability
(relative-decision) reliability is unaffected by design.

## Configuring priors (SCI-003)

The population-level priors are user-configurable for every builder except the subject-level
(sserr) family, whose intercept, between-person and residual log-SD constants stay fixed and
which exposes only
`priors.sserr.sig_trl`.
The four-event difference-of-differences model (analysis 9,
`psyrat_build_stan_dodiff`, and its native HMC twin) used to hard-code brms-style
constants that no `'priors'` override could reach; since 2026-09-05 (owner
ruling) they are exposed as `priors.dod` at the same values, so a default run
emits byte-identical Stan: `b_cell ~ normal(0, b_cell)` with `b_cell = 10`,
`b_sigma_cell ~ student_t(10, 0, b_sigma_cell)` with `b_sigma_cell = 2`,
`sd_id ~ student_t(10, 0, sd_id)` with `sd_id = 2` (the eight person SDs of the
location and scale block), `sd_trial ~ student_t(10, 0, sd_trial)` with
`sd_trial = 2` (the eight trial SDs), and `lkj_corr_cholesky(lkj)` with
`lkj = 2` on both correlation matrices. Those scales are on the microvolt scale
of the four cells; the student-t degrees of freedom stay fixed at 10. A result
saved before the family existed still prints `Priors: PsyRAT defaults` in its
export header, because a family absent from a stored prior set is at its
default by construction. For every design,
`psyrat_default_priors`
returns the documented ERP-scale defaults; the values are exposed on
`psyrat_prefs.proc.priors` and accepted as a `'priors'` name-value by
`psyrat_run` and by `psyrat_computevarcomp` (merged over the defaults, so
partial overrides are safe). A run started through `psyrat_run` records the
merged prior set in its run-config sidecar, so a sidecar replay reproduces the
same model. For example, to widen the between-person prior for a larger-scale
measure:

```matlab
pr = psyrat_default_priors;
pr.single.sig_u = 200;      % wider half-Cauchy scale
psyrat_run('file', 'mydata.csv', 'priors', pr, ...
    'savepath', pwd, 'savename', 'wide_priors.psyrat');
```

The same struct is accepted by a direct `psyrat_computevarcomp` call
(`'priors', pr`) for callers working below the scripted runner.

With the defaults, the generated Stan model is byte-identical to the previous
hard-coded models, so default-prior results are unchanged. The structural
non-centered `normal(0,1)` priors are not exposed. The Gaussian-family LKJ
correlation priors are fixed, except the difference-of-differences shape
`priors.dod.lkj`; the gamma family exposes its two LKJ shapes as
`priors.gamma.lkj` and `priors.gamma.modular_lkj`.

## Reproduction

The component draws for the old (pre-SCI-004) priors are frozen in the W0
baseline `tests/baselines/trt_time.mat`. The comparison re-runs the scenario
with the current priors and computes the coefficients with `psyrat_rel_trt`.
Because priors are now configurable, the same comparison can run both settings
programmatically by passing different `'priors'` structs.

## Splits reliability priors (nonparallel observed-design models)

The nonparallel data-splits designs (Rocha Tables 4/5; estimation cases 23
`ic_splits` and 24 `trt_splits`) reuse their single-occasion and test-retest
donor priors unchanged: `priors.single` for case 23 and `priors.trt` for case
24. Each split-mean model adds one explicit interaction component on top of the
donor, exposed under `priors.splits`:

| Component | Meaning | Prior | Case |
| --- | --- | --- | --- |
| `sig_splitxid` | person x split, sigma_ps | `cauchy(0, 5)` (half) | 23 (p x (i:s)) |
| `sig_posxid` | person x occasion x split, sigma_pos | `cauchy(0, 5)` (half) | 24 (p x o x (i:s)) |

Both default to a half-Cauchy scale of 5, matching the microvolt-scale crossed
interaction defaults `priors.trt.sig_trlxid` / `priors.trt.sig_occxid`. They are
merged over the defaults like every other family, so partial overrides via the
`'priors'` name-value are safe (see "Configuring priors" above). The parallel
splits path (`splits=2`) carries no splits-specific priors: it reuses the
one-facet / test-retest models unchanged, so the standard `priors.single` /
`priors.trt` apply.

**Modeling assumption (SD-scaling residual).** The input is one row per split,
where `meas` is the mean of `weight` (= n_i) items. Each split mean of n_i items
carries sampling variance sigma_err^2 / n_i, baked into the model as
`sig_obs = sig_err ./ sqrt(weight)`. The per-item residual `sig_err` and the
per-row person-by-split (or person-by-occasion-by-split) interaction are
identified **only** by variation in n_i across splits; near-equal n_i collapses
them. The preflight therefore errors on all-equal n_i and warns on low n_i
coefficient of variation under `splits=3`.

**Generalizability is a downward-biased lower bound.** Because only split means
are observed, the split-level item main effect collapses into `sig_err` and
cannot be separated from the person-involving residual. Absolute-decision
dependability (D) is exact, but relative-decision generalizability (G) is an
exact lower bound and the relative SEM an upper bound, by exactly the
sigma_(i:s) share. Every relative coefficient carries a `gbias` flag in the
output. The Rocha Table 4/5 formulas are summarized in the in-app Formula
Guide; the full derivation is maintained in the development repository.

## Gamma location-scale priors (`gammascale = 2`)

> **MIGRATION NOTE (2026-08-07): `gammascale = 2` is now the DEFAULT** on every
> design that supports it — analyses 6, 7, 8, 9, 10, 11, 12, 13, 14, 25, 26 and
> 27 — by owner rulings 33–36. It was previously `1` everywhere.
>
> **What this means for you.** A gamma run on one of those designs that does not
> name `gammascale` now fits the location-scale model and therefore a *different
> estimand* than the same script did before this date. The two are not
> comparable. **To reproduce earlier results, pass `gammascale = 1` explicitly**
> (CLI), or select "Log-nu (relative dispersion / CV)" in the GUI's new
> **Gamma scale parameterization** control.
>
> **What is unaffected.** Analyses 1–5 have no location-scale variant and still
> default to log-nu, as do the concurrent (correlated-residual) difference
> designs, which reject `gammascale = 2` outright. Gaussian runs are untouched.
>
> **Replaying an old run-config sidecar — read this before you do.** Runs made
> from **2026-08-03** onward stamp the *resolved* `gammascale` into both the saved
> result and the sidecar, so replaying one of those reproduces its own run. **A
> sidecar written before that date carries no `gammascale` at all**, and replaying
> it today re-derives the *current* default — so a movable gamma design that
> originally ran log-nu will silently replay as location-scale. **Pass
> `gammascale = 1` explicitly when replaying a pre-2026-08-03 config.** The same
> applies to a sidecar produced by calling `psyrat_write_runconfig` directly
> rather than through the normal run path, which records `[]` by design.
>
> **Why now.** The change was gated on every supported design being validated
> first (ruling 31); case 13 completed that in 2026-08. The rationale is that a
> log-nu contrast blends amplitude with residual SD and so cannot support the
> absolute-precision reading most users take from output labeled "error variance"
> or "dependability".
> sections G, H and K.

The location-scale parameterization replaces the log-nu scale submodel with a
direct log residual-SD submodel. **Every `gammascale = 2` design draws on this
block** — scope was widened beyond the original two dynamic designs by owner
ruling H (2026-08-04), and it now covers analyses 6, 7, 8, 9, 10, 11, 12, 13, 14,
25, 26 and 27. The static
designs (subject-level 6/25, difference 7/8/9/10) use `sigma_center`, `sigma_sd`
and `sig_sd` but carry no dimension slope, so `b_sigma` is simply unused there
rather than taking a separate default. Priors are **not** invariant under that
reparameterization, so four new defaults were derived rather than copied, using
the exact relation

```
log sigma = log mu + 0.5*log(2) - 0.5*log(nu)
```

| Prior | Default | How it was obtained |
|---|---|---|
| `gamma.sigma_center` | `log(5.5) + 0.5*log(2/18)` = `log(1.8333)` | reproduces exactly `nu = 18` at the reference point, matching `nu_center` |
| `gamma.sigma_sd` | 0.5 | `nu_sd = 0.75` implies ~0.375 after the 0.5 coefficient; widened to 0.5 to stay weakly informative and match `mu_sd` |
| `gamma.sig_sd` | 0.5 | `sig_nu = 1.0` halved by the same coefficient |
| `gamma.b_sigma` | 0.25 | **not** a translation, see below |

`b_sigma` cannot be derived from `b_nu`. The slope relation `b_nu = 2*b - 2*b_sigma`
is not invertible into a prior, because `b` and `b_nu` are specified
independently above (and the induced value would be ~0.56). The default is chosen
for interpretability on its own scale: `N(0, 0.25)` allows a one-SD change in the
predictor to shift log-sigma by up to ~+/-0.5 in the 95% prior interval, i.e. a
residual-SD change of up to ~1.65x per z-SD. That is generous for a precision
effect without being unbounded.

**Status: all four are owner-unconfirmed and no sensitivity analysis has been run
for them.** `b_sigma` is the one most worth checking, because it is a judgment
call rather than a derivation and it is the slope the parameterization exists to
make interpretable. The owner's Stan reference bundle uses tighter values for the
comparable quantities (`student_t(3, 0, 0.35)` for the person log-sigma SD and
`normal(0, 0.35)` for the scale slopes), so these defaults are the more
conservative choice in both cases. Override per run via the `priors` input.

## Stan console messages during estimation (verbose mode)

With "Verbose Stan output" enabled, CmdStan echoes its raw console log. Two
message families commonly appear near the start of sampling. Both are expected,
and neither reflects a problem with the model or the reported estimates. They do
not appear unless verbose output is selected.

**LKJ correlation warmup rejections.** Runs of the location-scale, difference-
score, dynamic-reliability, and residual-correlation designs may print repeated
lines of the form:

```
Informational Message: The current Metropolis proposal is about to be rejected ...
Exception: lkj_corr_cholesky_lpdf: Random variable[k] is 0, but must be positive!
```

These designs estimate a correlation matrix through its Cholesky factor. CmdStan
initializes parameters by default from a uniform draw on the unconstrained scale
(`init=2`), and an occasional early-warmup draw pushes a diagonal element of that
factor to underflow to zero, which the LKJ prior density rejects. The sampler
discards that single proposal and continues; adaptation quickly moves
initialization away from the boundary. As CmdStan's own text states, sporadic
rejections of this kind (here confined to early warmup) mean the sampler is fine;
rejections that persisted throughout sampling would instead point to an
ill-conditioned or misspecified model. The messages have no effect on the
reported variance components, generalizability or dependability coefficients,
SEMs, or convergence diagnostics (R-hat, effective sample size), all of which are
computed from the retained post-warmup draws. Plain one-facet and test-retest
models do not carry a correlation Cholesky and do not produce these lines. A
future opt-in initialization setting will allow starting from the identity
correlation matrix, which removes the messages; it is not the default because a
fixed common start reduces the diagnostic value of dispersed multi-chain
initialization.

**RDump data-format deprecation.** CmdStan may print:

```
Warning: file '.../temp.data.R' is being read as an 'RDump' file.
This format is deprecated and will not receive new features.
Consider saving your data in JSON format instead.
```

The bundled MatlabStan writes each model's data to CmdStan in the legacy R
"dump" format. CmdStan still reads this format across PsyRAT's supported range
(CmdStan 2.26 through the pinned 2.38) and returns identical results; the message
is a forward-looking deprecation notice, not an error, and has no effect on the
estimates. A migration of the data writer to JSON is planned (tracked in the
development repository) and will remove the notice.
