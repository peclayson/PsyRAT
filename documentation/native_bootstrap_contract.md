# Native REML bootstrap results

The native `fitlme` engine fits by REML and simulates a parametric bootstrap on
the observed incidence design. Its intervals are approximate frequentist
intervals. A reported mean of bootstrap coefficient draws is not necessarily
the coefficient evaluated at the original REML point estimate.

Original fits and refits use the same validity rule: recognized optimizer
warnings or nonfinite/complex estimates invalidate the fit. Zero variance
components remain admissible. MATLAB's **last** warning ID and message are
recorded for each fit/refit; this is not a complete optimizer warning history.
Thrown refit errors are recorded too. Failed refits are omitted, never replaced
with the fitted point. All components within a fit retain the same mask.

The approved adequacy floor is unchanged: at least 90% of requested refits must
succeed, and at least 200 successes are required if 200 or more were requested.
The floor applies to each independent fit. It does not establish interval
coverage or remove failure-selection bias. An invalid original fit remains
invalid even if its subsequent bootstrap passes the count floor.

## Draw access and compatibility

For analyses 1–5 and 23/24, independent output columns remain numeric when their
successful-draw counts are equal. When they differ, each cell contains the
complete numeric vector for one fit. Read either representation with:

```matlab
person_sd = psyrat_draw_column(REL.out.sig_u, column);
trial_sd = psyrat_draw_column(REL.out.sig_trl, column);
residual_sd = psyrat_draw_column(REL.out.sig_e, column);
```

Joint components within one fit are aligned by successful replicate index.
Equal row numbers in different independent fits do not establish joint
cross-fit uncertainty. No values are padded, resampled, or discarded to obtain
a rectangular array. External scripts directly indexing `REL.out.x(:,k)` need
the accessor for unequal-size results. Existing equal-size numeric results and
their scientific summaries keep their representation.

## Evidence in saved results and exports

`REL.engine_meta.fits(k)` records the output label, design, seed, formula,
original point/validity, requested/successful/failed counts, success mask,
retained indices and per-replicate warning/error records. Top-level counts sum
over fits. `draw_layout` and `draw_dependence` identify the storage contract.

The JSON sidecar retains `estimation.bootstrap_reps` as the **per-fit replay
setting**. `bootstrap_requested_total`, `bootstrap_ok` and `bootstrap_failed`
describe all fits, as stated by `bootstrap_counts_scope`; `bootstrap_fits`
retains individual records. Viewer notices and exported headers report counts
by fit. MCMC trace plots are unavailable for bootstrap results; bootstrap
replicates must not be reshaped as chains.

Exported provenance calls these **bootstrap intervals** and identifies the
approximate frequentist parametric-bootstrap method. It records Bayesian priors
and MCMC chain/warmup/sampling settings as not used. Sidecars keep those
configured values for later Bayesian replay, with `priors_applied=false` and
`mcmc_settings_applied=false`. Saved diagnostic thresholds describe original-fit
validity and the per-fit 90%/conditional-200 success floor; the internal Rhat/ESS
compatibility row is not a measured MCMC diagnostic. Plot window titles use the
requested interval percentage. Bayesian controls retain credible-interval
wording. These presentation changes do not alter interval endpoints or draws.
The fitlme `REL.out.conv.thresholds` record replaces `rhat_less_than` and
`ess_at_least` with the bootstrap criteria; external consumers of those fields
must use the engine-appropriate threshold record.

For historical native results, complete original/replicate validity records and
all aligned draws, if available, can support reassessment and re-summary under
the same validity rule and success floor. Refit when necessary records or
discarded draws cannot be recovered, or the original fit or success floor fails.
Saved draws alone do not identify which earlier refits emitted warnings, and
absent warning records cannot be reconstructed. Historical native ML results
still require REML refitting.

Re-display/re-summary/table export updates presentation labels. Explicitly
regenerate the run-config sidecar with `psyrat_write_runconfig` to add its
application flags; table exports do not write that sidecar. To refresh historical
threshold/scope metadata, assign the result of `psyrat_checkconv`, resave it,
then regenerate its sidecar.
That recheck does not recover missing original/refit validity history or rerun
diagnostic parsing/classification; it is not new scientific validation.

Repeated-dataset coverage and failure-selection validation remain separate from
this code correction and have not been run. It is therefore unknown whether the
native REML bootstrap intervals attain their nominal coverage or whether omitting
failed refits biases them, and passing regression tests does not settle either
question. As the user manual's chapter 9 (Viewing results) states, the native
engines are for getting a pipeline working; the reporting workflow remains
CmdStan, which carries its own diagnostic and uncertainty validation
requirements.
