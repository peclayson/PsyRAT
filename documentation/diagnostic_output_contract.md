# Diagnostic output contract — batches 1A/1B

The existing `out.conv.converged` flag remains an **early-workflow** check:
Rhat < 1.1 and effective sample size >= 10 × requested chains. It is not a
final-inference certificate. The owner approved preserving this policy on
2026-09-12; final-inference readiness is explicitly unassessed. Divergences
remain separately reported and do not change this flag. Native HMC acceptance
below .1 still fails; invalid, nonreal or greater-than-one acceptance also fails.
REML uses optimizer validity and the approved bootstrap count floor, not MCMC
diagnostics. Its uncertainty-calibration limits are separate.

## Saved schema and compatibility

Convergence data retain names, effective sample size and Rhat in columns 1–3.
Data columns 4–6 contain bulk ESS, tail ESS and mean MCSE. **Header cell 4 is
still the divergence count**, not a column label. Native acceptance ratios keep
their existing header positions. A header struct with schema
`psyrat_diagnostics/1` records source, precision, degradation, monitored names,
model-backed constant dispositions and scope. Read it using
`psyrat_diagnostic_metadata(table)`, rather than assuming a fixed table width.
`psyrat_checkconv` retains per-fit metadata in `out.conv.diagnostics`, plus
excluded/invalid names, thresholds, assessment scope and final-readiness state.
Saved results, runconfig sidecars, viewer notices and exported provenance carry
that distinction. Scope-only viewer notices appear after setup without blocking.
Failed checks, positive divergences, saved exclusions and recorded degraded
diagnostic formats require acknowledgment. Legacy exclusions are described as
recorded exclusions, without retroactively claiming model-backed proof.
DoD derived views retain the selected source fit's diagnostics;
they do not thereby assess the Monte Carlo error of the derived coefficient.

CmdStan's normal named CSV route remains at eight significant figures. Modern
CSV and recognized modern text use min(bulk ESS, tail ESS). Legacy N_Eff remains
in column 2, with unmeasured bulk/tail fields NaN. Printed fallback is explicitly
degraded at two significant figures; unknown positional text is labeled as
ambiguous. A malformed precise source cannot be cleared by a healthy rounded
fallback. Required rows come independently from every chain CSV header; missing
rows or unavailable/mismatched inventories fail closed. This inventory covers
the existing model-specific monitored prefixes, not every Stan latent variable.
Both formats reject nonreal diagnostic tokens. Text fallback retains malformed
or duplicate monitored rows as failing rows: an invalid string cannot become
the double-NaN signature of a permitted constant. Literal `nan`, `-nan` and
`+nan` remain valid source tokens, subject to the model-backed constant rule.

Only deterministic entries supported by the actual emitted source and data may
exclude double-NaN diagnostics. Recognized rules cover declared correlation
diagonals, structural Cholesky entries, bounded principal-block aliases and the
reviewed disabled case-8 residual-coupling graph. Unexpected stochastic NaNs,
single NaNs, infinities, malformed rows, all-constant tables and missing rows
fail. Unknown source syntax cannot authorize an exemption. The classifier is a
bounded recognizer, not a general Stan interpreter. Model/data hashes and
relevant flags accompany its evidence.
Nested and chained index assignments to a candidate alias invalidate its
exemption; nested reads alone do not. This conservative write check does not
attempt to prove that an overwrite preserves the original constant.
Aliases must be complete, unconditional top-level principal-block assignments.
A correlation parameter is available immediately; a transformed/generated
parent must have exactly one write: a complete unconditional assignment that
precedes the copy. A copy made from a temporary stochastic parent
does not inherit the parent's final constraint. The final typed parent's own
constant entries remain permitted.
The case-8 rule matches its data/output declarations and complete reviewed
top-level loops, including the order of calculations. Matching text in another
function, a prefixed identifier, a conditional/partial loop or a print string
does not prove the graph. Stan string/comment masking is internal to recognition;
saved source hashes still identify the original unmodified source.

## Native diagnostic calculations and scope

Native HMC now uses rank-normalized and folded split Rhat, bulk ESS, 5/95% tail
ESS and raw-scale mean MCSE. Odd chains omit their middle draw only for splitting;
saved samples are unchanged. Original frozen chains, nonfinite draws and fewer
than three observations per half fail. Six-draw and one-chain numerical results
are computational controls, not evidence of adequate inference. The original
sampler, log densities, priors, RNG and draw ordering are unchanged.

The independent reference is R posterior 1.6.0 on identical binary doubles.
Original decimal-reference artifacts remain preserved: R/MATLAB decimal decoding
differed by one ULP for some inputs, altering folded ties in two fixtures. No
reference algorithm or comparison tolerance was changed. The exact-zero native
guard deliberately avoids posterior's absolute `range < eps` guard, which can
declare moving tiny-scale quantities constant merely because of units. All
numerically assessed binary reference fixtures agree at the declared tolerances.
The same/different-constant, mixed-constant-chain and nonfinite fixtures instead
assert PsyRAT's separate hard-failure guards; those outputs are not required to
equal the reference's numeric estimates. CmdStan 2.38's tail
ESS differs on a scale-separated fixture with constant tail indicators; that
independent discrepancy is documented, not hidden by an equality assertion.

All required returned fields and dimensions are checked against constructor
schemas for the 20 native keys. Only declared zero-width dimension coefficients
may be absent numerically. All returned quantities are diagnosed with their
original contiguous chain coordinates. Subject log-SDs are additionally
diagnosed as `exp(2*er_var_ss)`; location-scale subject residual variances use
`exp(2*(Intercept_sigma+ind_sd))` at dimension zero. The test–retest constructor
now also retains `ind_sd` for this purpose. This is an additive native-fit field,
not a new reliability estimand.

Known DoD correlation diagonals are mathematical constants despite floating
roundoff. Raw values are retained. A gamma_32 forward-error bound, derived from
the C=4 constructor's operations, checks the unit-diagonal invariant before
exclusion; it is not a variance tolerance or an empirical constant detector.
Covariance diagonals and nonfirst Cholesky diagonals remain stochastic.

Dynamic surfaces away from the reference dimension, DoD subject/trial conditional
residual combinations, G/D/ICC ratios, interval endpoints and trial-count cutoffs
remain unassessed by these engine diagnostics. Numeric mean MCSE for a Cauchy
fixture does not establish a finite population mean. Five/ninety-five-percent
tail ESS does not establish the precision of 2.5/97.5% interval endpoints.

## Reuse and pending decisions

Preserved native HMC draws can be re-diagnosed if chain ordering and the required
quantities survive. A saved renamed/subset result may omit latent quantities
needed for the expanded coverage; those cases require a refit. CmdStan chain
CSVs can regenerate precise diagnostics and source-backed classifications;
rounded summaries alone cannot. Existing numeric coefficients need not change
merely because a diagnostic check changes, but a newly failed check requires
investigation before scientific use. Native ML and warning-contaminated or
incompletely retained bootstrap results have the refitting requirements in
`native_bootstrap_contract.md`.
Loading a historical table or calling `psyrat_checkconv` on it does not rerun
the parser or reconstruct its constant contract: it uses the stored ESS/Rhat
and permitted-constant metadata. Historical passes can persist until diagnostics
are rebuilt from the original evidence. This patch does not automatically
migrate or retroactively certify those saved tables.

Final-inference MCSE tolerances remain an owner decision for each reported
quantity. For context, mean MCSE is about .158 SD at ESS 40 versus .05 SD at
ESS 400 under finite variance; this does not specify endpoint/cutoff precision.
An alternative stricter screen (Rhat <1.01, bulk/tail ESS >=100 per chain,
zero divergences and precise modern diagnostics) would still leave those
quantity-specific requirements open. The approved option in this patch keeps
the existing early threshold, precise diagnostics and explicit unassessed final
state. No stricter gate or new scientific target was silently selected.

The unit suite (`runAllUnitTests`) exercises this contract, including the
independent R posterior reference arrays in `tests/fixtures/rank_diagnostics/`;
the user-visible changes are summarized in `CHANGELOG.md`. Long repeated-dataset
REML coverage and failure-selection validation has not been run and is not
implied by regression checks; see `native_bootstrap_contract.md`.
