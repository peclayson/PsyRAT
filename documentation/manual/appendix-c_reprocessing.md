# Appendix C. Reprocessing saved gamma results

`psyrat_reprocess_gamma_result` brings a saved gamma (scaled chi-square)
reliability result up to the `gamma_target_version` 2 residual-target contract
adopted in September 2026. It is explicit and versioned: nothing is
reprocessed automatically when a file is opened. If you do not run this tool, an older result keeps its
original meaning and is labeled as such. The pooling convention behind the
contract is described in [Chapter 15](15_advanced-designs.md), Section 15.4.

## C.1 What changed and why it needs a versioned migration

Two families of gamma results are affected.

**Static subject-level test-retest (analysis 25, `trt_sserrvar`).** Earlier
results encoded each participant's residual as the observation variance alone
(`C_p`, with `C_typ` for the typical person). The current contract adds the
population highest-order (person x occasion x trial) mean-surface remainder
`R`, so the residual is `E_p = R + C_p` and the reference is `E_typ = R + C_typ`.
Because the stored `pop_sdlog` / `ind_sdlog` are log-SD encodings of the
residual, the saved draws must be re-encoded once. Doing it twice would add
`R` twice, which is why the result carries a version stamp and the tool
refuses to pool again. Because `R` is nonnegative, a migrated participant's
coefficient is lower than before wherever `R` is positive once you recompute
the summary. At an unchanged cutoff, a participant who was included can
therefore be excluded, and the Save IDs lists can change. Recompute the
summary with `psyrat_relsummary` before you apply a cutoff to the migrated
file.

**Dynamic non-difference designs (analyses 11, 14, 26, 27).** The saved draws
already contain everything needed. The surface target (population-marginal
observation variance, or the typical person with zero person effects) is now
an explicit, labeled choice, and the two-facet subject-level design (27) uses
`R + C_p` only when the result is stamped. Migration validates the draws and
stamps the target; the raw draws are not changed, and no refit is needed. What
the next summary derives from them can change. A stamped analysis-27 result is
summarized with `R + C_p`, so its per-participant coefficients are lower than
before the migration wherever `R` is positive, and their comparison with any
cutoff you apply can change (the dynamic designs store no inclusion flag). An
explicit `surface_target`, whether given at the first migration or at a later
one, can change the dynamic surface of any of these designs; the migration record
says so.

## C.2 Usage

```matlab
% saved file -> fresh file (the source is never modified)
[updated, record] = psyrat_reprocess_gamma_result('study_gamma25.psyrat', ...
    'output', 'study_gamma25_v2.psyrat');

% in-memory result (nothing is written unless 'output' is given)
[updated, record] = psyrat_reprocess_gamma_result(psyrat_data);

% dynamic result with an explicit surface target
[updated, record] = psyrat_reprocess_gamma_result(psyrat_data, ...
    'surface_target', 'typical_person');
```

The first argument, `source`, may be a full `psyrat_data` struct, a bare `REL`
struct, or a `.psyrat`/`.mat` file holding `psyrat_data`. Two options apply.

- `'output'`: the destination file. File input requires it; the destination
  must end in `.psyrat` or `.mat`, must not exist, and must differ from the
  source. Every other variable stored in the source file is carried into the
  destination unchanged. Struct input writes nothing unless `'output'` is
  given, and a bare `REL` struct cannot be written (there is no wrapper to
  preserve).
- `'surface_target'`: `population_marginal` (default) or `typical_person`. It
  applies to the dynamic designs only; for analysis 25 it is ignored, and the
  record shows `not_applicable`. If the result already carries an explicit
  target and none is given, it is kept. The historical summary settings
  `marginal = 0` and `marginal = 1` (the latter was ignored) both meant the old
  population surface; that flag is never read as a target choice.

Spreadsheet or text exports (`.xlsx`, `.xls`, `.csv`, `.txt`, `.tsv`) are
refused. They hold summaries, not posterior draws, so no reliability can be
rebuilt from them. Reprocess the `.psyrat` result they came from, or refit.

Versions. The only stamp version accepted is exactly `2`. An ABSENT
`gamma_target` field is the historical version-1 encoding, the only legacy
form ever produced. Any other explicit version (`1`, `99`, non-numeric) is
refused, and so is a `gamma_target` that is present but empty or malformed:
partial provenance is never treated as "unstamped" (that would add `R` a
second time to an already-migrated result). A `gamma_resid` block without a
stamp is likewise refused as inconsistent provenance, because that schema
exists only under version 2.

Storage. The per-stratum draw blocks (`gro_sds`, `b`, `b_nu`/`b_sigma`,
`chol_corrmat`, `ind_bs`, `ind_nu`, `ind_sdlog`) must be stored as the
estimator writes them, as cell arrays of scalar draws. A block stored as a
plain numeric array is refused (`unsupported_storage`) rather than converted,
because the summary, table, and provenance code could not open the migrated
result. Refit such a result.

Shapes. The `ssinfo` columns `z1`, `z2`, `trls`, and `id2` must hold exactly
one finite value per row; `id` must be a one-per-row numeric, cellstr, string,
or categorical column (char matrices and multi-column IDs are refused). IDs
must be unique within a stratum; the same ID may recur across strata.

Parameterization. `REL.gammascale` must be exactly `1` (log-nu) or `2`
(location-scale) when present. Every other identifier of the scale
parameterization (`gamma_disp` raw fields, raw `b_nu`/`b_sigma` fields, the
stamp, `gamma_resid`) must agree; a contradiction is refused rather than
restamped. With `gammascale` absent, only unambiguous provenance is accepted
(for example location-scale `gamma_disp.log_sigma`/`s_sd`); nothing defaults to
log-nu.

## C.3 What the tool does for analysis 25

For each stratum and each draw it reconstructs `R` from the retained
observed-scale quantities, without a refit and without any person effect,
in this normalized form:

```
a = (sd_P/mu)^2,  b = (sd_O/mu)^2,  c = (sd_I/mu)^2,  h_PO = (sd_PO/mu)^2 (h_PI, h_OI likewise)
q_P = a, q_O = b, q_I = c
q_PO = (h_PO - a b) / ((1+a)(1+b))      (q_PI, q_OI likewise)
v_A  = log1p(q_A),  alpha = log(mu) - sum(v)/2
R    = population highest-order remainder from (alpha, sqrt(v))
```

where `sd_P = gro_sds(:,1)`, `sd_O = sig_occ`, `sd_I = sig_trl`,
`sd_PO = sig_occxid`, `sd_PI = sig_trlxid`, `sd_OI = sig_trlxocc`.

Before any pair variance is formed, each pair deficit `h - a b` is classified
against a declared per-draw error bound (unit roundoff `2^-53`, a 32-operation
budget for the legacy forward arithmetic, a condition factor
`1 + 2|log mu| + sum|v_raw|` (with `v_raw` the pair inverse before any
canonicalization), and the pair/main/inverse terms listed in the function
header). A deficit above the bound is a genuine positive pair effect
and is retained, however small. A deficit below minus the bound means the
stored components cannot have come from the gamma log-link model (for example
a missing or corrupted interaction component) and the tool refuses, naming the
component, stratum, draw, deficit and bound. A deficit inside the bound is
recorded as roundoff-unresolved, and the record marks that stratum as not
exactly reconstructed; a negative deficit inside the bound is also
canonicalized to zero, while a positive one is kept. The bound is a
conservative model of the arithmetic in the PsyRAT version that wrote the
version-1 result, not a certified MATLAB bound, and it is not a fixed absolute
floor. Overflow or underflow anywhere in the reconstruction is refused as an
unsupported numeric regime. The tool then re-encodes the residual:

```
E_typ = R + exp(2*pop_sdlog_v1)
E_p   = R + exp(2*(pop_sdlog_v1 + ind_sdlog_v1))
pop_sdlog_v2 = .5*log(E_typ)
ind_sdlog_v2 = .5*log(E_p) - pop_sdlog_v2
```

Column order of the per-participant blocks and the `id_matches` rows are
preserved exactly. The raw pieces are stored per stratum in
`REL.out.gamma_resid{s}` (`mean_remainder`, `obs_var_typ`, `obs_var_ss`,
`scale0`, `ind_scale`, `parameterization`). Under location-scale the fitted
log-sigma draws are recovered exactly from the version-1 encoding; under
log-nu the raw `Intercept_nu`/`ind_nu` draws were not retained by version 1
and are recorded as absent, never manufactured. `REL.out.gamma_disp` (which
holds only the last stratum) keeps its raw fields; its `sigma_pot_e2_pop`
becomes `E_typ` of the last stratum and the version-1 value is kept beside it.
A cached `psyrat_data.relsummary` computed under version 1 is removed from the
new result (recorded); recompute it with `psyrat_relsummary`.

When `gamma_disp` retains the last stratum's raw log-scale pieces, the
normalized components they imply are compared with the stored ones under
twice the same per-draw bound; a discrepancy is refused.

A result that already carries version 2 is validated in full before it is
returned unchanged: alignment, IDs, finite positive `C_typ`/`C_p`, finite
nonnegative `R`, finite positive `E`, finite encodings that decode to finite
positive values, and the encoding identity checked elementwise in the log
domain,

```
|log(R + C) - 2 (pop + ind)| <= gamma_8 (1 + |log(R + C)| + 2|pop| + 2|ind|),   gamma_8 = 8u/(1 - 8u), u = 2^-53
```

For the typical reference, the right side uses `2|pop|` alone. The bound
counts the producer's
sum, log, scaling and subtraction and the validator's sum, log and
re-summation, each contributing at most `u` times its operand magnitudes;
`pop` and `ind` enter because their sum cancels when both are large. It is
not an absolute floor: a 100% mismatch such as stored `E = 1e-12` against
`R + C = 2e-12` (`|log 2| = .69`) is refused. Nothing is pooled twice.

### Numeric regime

The error model assumes finite NORMAL binary64 arithmetic. Every stored input
(`mu` and the six SDs), its square, the normalized ratios `(sd/mu)^2`, the
prefactor `B = mu^2`, the reconstructed `q` and `v`, the decoded `C`, the
residual `E`, and the encodings must have magnitude in `[realmin, realmax]`
(`2.2e-308` to `1.8e308`). A component may be exactly zero only where it is
structurally absent, that is, where the stored source value itself is zero
(or a negative pair deficit inside the bound was canonicalized). A strictly
positive source that maps to zero, a subnormal, or `Inf` under any transform
(its square, its normalization by `mu`, the pair quotient, `log1p`, or the
remainder itself) is positive-to-zero underflow or overflow and is refused
with a message naming the transform; it is never carried forward as a
structural zero.

In addition, every positive term of the remainder (a product of `q` over a
factor subset with support union `{P,O,I}`) is checked in log space, including
its left-to-right partial products and its product with `B`. A tiny-but-normal
`q` whose product would underflow inside the remainder evaluation is refused as
`unsupported_numeric_regime` instead of producing a truncated `R`. A complete
zero facet still yields the legitimate `R = 0`, and moderate small magnitudes
(for example `mu = 1e-100` with SDs `1e-101`) remain supported. Subnormal
inputs, overflowing squares, and underflowing decoded variances are refused the
same way; nothing is fabricated.

Consider two examples. With `mu = 1e100`, `sd_P = 1e-100`, `sd_OI = 1e150`,
and the other SDs 0, the squares `1e-200` and `1e300` are normal, but
`(sd_P/mu)^2 = 1e-400` underflows to 0, which would wrongly make every
P-supported remainder term look absent (`R = 0` instead of the true
`1e-100`); the tool refuses. With `mu = 1e50` and all six SDs `1e-40`, the true
`R ~ 6e-260` is normal, but every product of two or three `q ~ 1e-180`
underflows, and the tool refuses for the same reason.

## C.4 Records

`record` (also stored at `REL.out.gamma_target.migration_history(end)` when
the status is `migrated` or `target_restamped`) lists: version from/to, status
(`migrated`, `already_current`, `target_restamped`), source and output, time,
analysis, parameterization, method, surface target, fields read, written,
derived, and absent, per-stratum transforms (draw and participant counts,
roundoff count, cross-check flag, `R` summary), count source, ID mapping,
invalidation, limitations, and notes (for example, an ignored historical
`marginal` flag). Running the tool on an already-current result returns it
unchanged with status `already_current`.

`REL.out.gamma_target` carries: `gamma_target_version` (2), `surface_target`,
`reference_target`, `subject_target`, `remainder_scheme`, `parameterization`,
`count_source`, `surface_count_source`, `analysis`, `migrated_from_version`,
`migration_history`. For the dynamic designs `reference_target` is
`not_applicable`, `surface_count_source` is `default_stratum_median`, and
`count_source` is `not_applicable` (11, 14) or
`median_trials_per_occupied_occasion` (26, 27). For analysis 25
`reference_target` is `typical_person_zero_person_effects`, `count_source` is
`mean_trials_per_occupied_occasion`, and `surface_target` and
`surface_count_source` are `not_applicable`. A stamped result whose
identifiers disagree with this contract is refused; only `surface_target`
may be changed, and only by an explicit request.

## C.5 Dynamic designs: what is validated

The tool checks that the population draws and slopes are finite and that the
facet and heterogeneity SDs are finite and nonnegative. Each stratum must have
exactly one label and one positive count (`ntrials`, plus `nocc` for two
facets). `REL.ndim` must be exactly 1 or 2, with `dim_names` of matching
length, `dimz` finite with `ndim` columns, and a slope basis of width 1
(`ndim = 1`) or 3 (`ndim = 2`). A complete draw-aligned `[ndraws x 2 x 2]`
Cholesky correlation block is required per stratum under both
parameterizations and both surface choices, because the summary, variance
table, and provenance code read it unconditionally, and a reusable result must
be able to switch to the population-marginal primary surface. For subject-level
designs, the tool also requires finite own `z1`/`z2`, positive `trls`,
non-missing participant IDs unique within a stratum (the same ID may recur
across strata), `id2` a complete permutation, and own-effect blocks of size
draws x participants.

## C.6 Precise refusals

The tool stops with an identifier of the form
`psyrat_reprocess_gamma_result:<reason>` and a message naming the field,
stratum, draw or dimension involved: `missing_field`, `misaligned_draws`,
`misaligned_subjects`, `invalid_id_matches`, `invalid_ids`, `invalid_ssinfo`,
`invalid_strata`, `invalid_dimensions`, `invalid_correlation_draws`,
`invalid_source_values`, `invalid_inverse`, `unsupported_numeric_regime`,
`unsupported_storage`, `inconsistent_provenance`,
`inverse_crosscheck_failed`, `invalid_v2_claim`, `contradictory_stamp`,
`unsupported_target_version`, `unknown_parameterization`,
`contradictory_parameterization`, `ambiguous_parameterization`,
`unsupported_family`, `unsupported_analysis`, `summary_only_input`,
`output_required`, `output_exists`, `output_is_source`, `output_not_allowed`,
`bad_output`, `no_psyrat_data`, `bad_option` (an unknown or unpaired option
name), `bad_surface_target` (a value other than `population_marginal` or
`typical_person`), `bad_source`, `bad_source_file` (a source file that is not
`.psyrat` or `.mat`), and `source_not_found`. An output path matching the source
is refused as `output_is_source`; another existing destination is refused as
`output_exists`. Both refusals preserve the existing files. A result whose draws are
incomplete for its
design needs a refit; the tool does not guess missing effects, correlations
or a default parameterization.

## C.7 Out of scope

Only analyses 11, 14, 25, 26, and 27 are reprocessed. Gaussian results are
refused with `unsupported_family`; the difference-score, copula,
difference-of-differences (DoD), and splits designs, and the static gamma
designs other than analysis 25 (analyses 1 through 6), are refused with
`unsupported_analysis`. The tool does not compute reliability itself. After
migration, recompute an analysis-25 summary with `psyrat_relsummary` (analysis
key `trt_sserr`; see [Chapter 14](14_scripting.md)); for the dynamic designs,
open the migrated file in Dynamic Reliability: View Results or pass it to
`psyrat_report`, which summarizes the draws directly.
