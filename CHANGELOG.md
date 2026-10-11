# Changelog

All notable changes to PsyRAT are recorded here. This file is the **only** release record.

Everything before `0.1.0-beta` belongs to the ERA Toolbox lineage, not to PsyRAT. That history was
kept in `documentation/VersionHistory.pdf` (ERA 0.4.4 through 0.5.2, last released 30 October 2020),
which was **retired on 2026-08-22 (D1)**: it was binary rather than diffable, had no in-repo source,
and described a toolbox this one descends from rather than this one. It is not lost — recover it
with `git show 0a1f694:documentation/VersionHistory.pdf` from a clone of the development repository
(the release archive carries no git history) — but it is not maintained and nothing here depends on it.

Version strings come from `psyrat_defineversion` — the single source of truth. The same string is
printed at startup, stored in `psyrat_prefs.ver`, and written into every exported table header and
`*_runconfig.json` sidecar. Between releases the string is shared by every commit, so the
sidecar also records the git source revision (`source_revision`); report the version plus that
revision to identify the code that produced your estimates.

**Format for new version entries.** The retired ERA history used a shape worth keeping, and the
maintainer asked for it on 2026-08-22: per version, a release date, then `### Major changes` and
`### Minor changes`, each a flat bullet list of what actually changed. Apply it to entries added
from here on; existing entries were not restructured.

---

## 0.1.1-beta — 2026-10-06

Correction release from the post-beta scientific review (the Q, PBA, G, GUI and F items named
below, from the development trackers). No prior, likelihood or sampling target changed, and
CmdStan's posterior draws are unchanged in every design. What changes is listed below: the
fitted variance components of the native fitlme engine, which now fits by restricted maximum
likelihood; the nonparallel-splits summary under unequal split counts per participant
(analysis 23) and, under unequal split counts across the observed person-by-occasion cells or a
missing occasion, the multi-occasion splits summary (analysis 24); the per-participant
coefficients of the gamma two-facet subject-level designs (analyses 25 and 27), in a new fit
directly and in a saved result only after an explicit versioned reprocessing; the
dynamic-reliability surface of gamma non-difference runs, whose Estimand selector was inert;
the retained and excluded participant rosters where a participant is missing from an event or
an ID is reused across groups, and the cutoff projections beside them; and the
difference-of-differences grid display, whose admissible-cell arithmetic is unchanged.
Diagnostics, labels, displays and exports change as listed.

### Major changes

- **The three gamma two-facet designs (analyses 14, 25 and 27) are reclassified from Supported
  to Restricted and disclosed at Analyze, in the viewers, on saved figures and in every export
  header (owner rulings 2-4, 2026-09-21; queue rows 3a and 29).**
  Location-scale and log-nu alike. With three or four occasion levels in each reviewed capture of
  these designs, the population variance components and the pooled residual reported as
  posterior means are not recovered at level. As variance-scale ratios of posterior mean to
  planted value, the five re-evaluated scope captures put the occasion component at 2.4 to 12.1
  times, the highest-order remainder at 2.2 to 7.2 times and the pooled residual at 1.1 to 1.9
  times, with several 95% intervals excluding the truth; the retained analysis-27
  location-scale capture, reviewed as a control, read 1.42, 1.66 and 1.02, which is not a new
  pass. Two prospective analysis-27 fits put the same components at 60.1 and 2.45, 31.6 and
  3.35, and 1.83 and 2.95 times, the second of them non-converged by the toolbox's own
  criterion. The pattern is consistent with raw-scale components that are exponentials of
  squared log-scale SDs whose posteriors are right-skewed with few occasions.
  G and D stay within about .03 and .06 of the truth for analyses 14 and 25; for analysis 27
  both G and D are biased downward in the converged fits (medians across participants of the
  posterior-mean coefficient minus its planted value: G about -.05 location-scale at four
  occasions and -.10 log-nu at three, D -.13 in the converged prospective fit, -.10 log-nu and
  -.02 in the retained capture; the non-converged second prospective fit gave D +.08) while
  participant ordering is kept (r >= .86). Wherever production was
  compared with an independent enumeration it agreed to 2e-12 or better, which shows
  implementation consistency, not recovery, calibration, identifiability or sampling
  validation; the reclassification itself changes no estimate. New: a modal notice at Analyze before the save dialog (Cancel
  there backs out), a note in the results viewers (`psyrat_startview_dynrel` for 14 and 27,
  `psyrat_startview_trt` for 25), a shorter note above the panels of the surface figure (14,
  27) and of the two subject plots (25), kept in a saved PNG, PDF or FIG with a
  `design_maturity` record in the FIG's UserData, and a `Design maturity: RESTRICTED` line in
  every exported table header, all from one text source (`psyrat_gamma_twofacet_restriction`). The
  predeclared LS27 recovery assertion (`tests/integration/TestGammaTwoFacetDynrelSserrLsRecovery.m:215`)
  keeps its bound and its recorded FAIL, which the release checklist now names as the
  adjudicated expected outcome. The design maturity table below and manual Section 5.4 carry
  the new row; the one-facet gamma dynamic designs (11, 26) do not show the pattern and stay
  Supported. The Gaussian arms are unaffected.

- **Native bootstrap validity and independent-fit output are now auditable (PBA-05, batch 1A).**
  Original/refit warning and finite-real checks agree; valid zero boundaries and the approved
  success floor are preserved. Unequal successful counts retain per-fit vectors, with masks,
  indices and warning/error evidence in results and sidecars, and counts in notices/exports.
  Failed refits are never padded, and independent successful draws are not discarded to force
  rectangular arrays. The fitlme trace action no longer treats bootstrap replicates as chains.
  See `documentation/native_bootstrap_contract.md` for saved-result compatibility and refitting.

- **The native fitlme engine now fits by restricted maximum likelihood (REML), as its label and
  the manual state.** Every `fitlme` call in the engine omitted `FitMethod`, so v0.1.0-beta fitted
  by maximum likelihood (ML) under a REML label, and the direct-fit reference tests repeated the
  omission. Variance components from the native fitlme engine (analyses 1-4; 5 when it is the
  chosen engine; 23/24 and parallel splits) change; CmdStan and native HMC results are unchanged.
  Rerun any analysis that used the native fitlme engine. (Post-beta review PBA-04, Q1.)
- **Failed parametric-bootstrap refits of the native fitlme engine are dropped, disclosed and
  floored.** A refit that failed used to be replaced by the original point estimate, an
  artificial point mass that narrowed the interval, with a warning only above a 10% failure
  rate. Now every failed refit is dropped and warned about; the run is refused
  (`PsyRAT:fitlmeBootstrapFloor`) when fewer than 90% of the requested replicates succeed, or
  fewer than 200 when 200 or more were requested; and the achieved count is recorded: exported
  table headers gain `Bootstrap replicates: N of M succeeded`, the runconfig sidecar gains
  `estimation.bootstrap_ok` and `bootstrap_failed`, the results viewer's native-engine notice
  names them, and the saved result carries the engine record as `rel.engine_meta`. (PBA-05, Q3.)
- **CmdStan convergence diagnostics are read from the stansummary CSV by column name at full
  precision.** They were parsed positionally from the printed table at two significant figures,
  so the stored effective sample size was `N_Eff` on CmdStan before 2.33 and the tail ESS from
  2.33 on, and a 1.01 R-hat criterion could never have been applied to the stored value. The
  gate's `n_eff` is now the smaller of bulk and tail ESS; bulk ESS, tail ESS and MCSE are stored
  beside it for every monitored parameter (data columns 4 to 6 of the convergence table; header
  cell 4 remains the divergence count, so their labels sit in header cells 5 to 7); the printed
  table is unchanged. A fit whose smaller ESS is below ten times the chain count now reads not
  converged where the tail value alone passed. When the CSV cannot be produced, the previous text
  parse is used and a warning says so. (PBA-01, Q4; the prerequisite for the G9 gate tightening.)
- **The native HMC engine's convergence diagnostics fail on frozen chains, and a chain's
  acceptance ratio now gates the fit.** The engine's R-hat returned 1 when every chain had zero
  variance and its effective sample size returned the full draw count when the autocovariance
  was zero, so four chains stuck at four different values reported R-hat 1 and passed the
  check; the effective sample size also centered each chain on its own mean, so chains sitting
  at different values reported a full sample size. R-hat is now infinite for frozen chains, the
  effective sample size is computed on the split chains against the between-chain-aware
  variance (Vehtari et al., 2021; batch 1B adds that reference's rank normalization to this
  statistic, and its folded R-hat companion, next bullet) and collapses when chains
  disagree, the fraction of accepted proposals is recorded per chain in the convergence table,
  and a chain that accepted fewer than 10% of its proposals marks the fit not converged. Native
  HMC runs of v0.1.0-beta that passed the check with a stuck chain would have reported
  unreliable components; the draws themselves are unchanged. (PBA-03, Q5.)
- **Native HMC convergence diagnostics follow the rank-normalized and folded conventions, a
  parameter is excluded from the gate as a constant only on model evidence, and the source of
  every CmdStan diagnostic is recorded (batch 1B).** The native HMC engine's split R-hat is now
  the larger of the rank-normalized and the folded statistic, its effective sample size is
  reported as bulk and tail ESS with the gate reading the smaller, and the mean MCSE is stored
  beside them in the saved convergence table. The kernels were checked against R `posterior`
  1.6.0 on identical binary inputs, with one deliberate departure: the guard that declares a
  chain constant is stricter than that reference's absolute range test, which can call a moving
  quantity on a tiny scale constant. Every required returned field and dimension is checked
  against its constructor's schema, and every returned quantity is diagnosed in its original
  chain coordinates, so a missing, nonfinite or frozen quantity fails instead of passing unseen.
  A parameter may be excluded from the gate as a constant only when the emitted model and data
  show it is one (declared correlation diagonals, structural Cholesky entries, bounded
  principal-block aliases, the reviewed disabled residual coupling of analysis 8); an unexpected
  NaN on a parameter that should vary fails closed, and a saved result whose exclusions predate
  this rule is described as carrying recorded exclusions, not certified. Which parser produced
  the CmdStan diagnostics (the named stansummary CSV at eight significant figures, degraded to
  the legacy `N_Eff` with no bulk or tail ESS on CmdStan before 2.33, or the degraded
  two-significant-figure printed fallback, whose modern, legacy and ambiguous positional layouts
  are distinguished), and whether it degraded, persists in the saved result, the runconfig
  sidecar, the viewer notice and the export header; a degraded format or a saved exclusion asks
  for acknowledgment in the viewer. The existing R-hat < 1.1 and ESS >= 10 x
  chains threshold is labeled an early-workflow check; final-inference readiness is explicitly
  unassessed. Sampling targets, priors and the draws are unchanged. See
  `documentation/diagnostic_output_contract.md`. (PBA-01/Q4, PBA-02/Q6c, PBA-03/Q5; batch 1B.)
- **Nonparallel-splits dependability averages each participant's own error term.** The
  single-occasion splits summary (analysis 23) used the median number of splits per participant
  and a harmonic mean of the items per split pooled over the whole stratum. That is exact only
  when every participant contributes the same number of splits; with unequal counts it
  misweighted participants (a two-participant probe reported a dependability of 0.683 for a
  true 0.578). The coefficients now use each participant's own split count and item counts,
  averaged over participants, through the effective design the unchanged Table 4 closed form
  is evaluated at; balanced designs, including the chapter 13 tutorial, are unchanged. The
  reported number of splits is then the harmonic mean of the counts, preflight warns when the
  counts differ, and the formula guide, the viewer note and the design note say that
  dependability is exact for the person-averaged error and that the relative coefficient's
  exclusion of the split main effect is approximate under unequal counts. The multi-occasion
  design (analysis 24) is corrected separately, next. (PBA-06, Q7.)
- **The multi-occasion nonparallel-splits design (analysis 24) scores each participant's available
  occasions and averages the universe and error variances over participants.** The summary applied
  the same stratum-level reduction as analysis 23 (a median number of splits and a pooled harmonic
  mean of items per split), which is exact only when every participant contributes the same
  schedule; a participant with a missing occasion, or unequal split counts across the observed
  person-by-occasion cells, was misweighted. Each split now carries the weight 1 / (occasions
  available to that participant x splits in that occasion), a missing occasion carries no weight,
  the universe-score and error variances are averaged equally over participants within each
  group-by-event stratum and posterior draw, and the coefficients and SEMs follow from those
  averages (neither an average of per-participant coefficients nor the reliability of a group
  mean). Balanced designs are unchanged (the balanced fixture's CES dependability stays
  .528085483838); an irregular fixture moves from .598354525056 to .547429646737 and a
  missing-occasion fixture from .528357962520 to .428789173073. A numeric split or occasion count
  typed into the viewer or passed to `psyrat_report` remains a balanced reference projection and
  is labeled as one; the observed score is labeled as conditional on the fitted design. Preflight
  warns (`preflight:splitsUnequalCounts`) when cells hold unequal split counts or participants
  differ in available occasions; the viewer's banner names both the available-occasion design and
  the count projection, and the coefficients window states which of the two is shown; the export
  headers and `psyrat_report` carry the design note and each stratum's incidence and
  projection metadata, and the multi-occasion splits windows were made large enough for the design
  and limitation notes to read whole; analysis 23's own numerics are untouched by this correction.
  An existing analysis-24 fit needs re-summary and re-export from the saved fit and its scored
  rows, not a refit; reopening an old summary does not reconstruct the new score. (PBA-06, Q7;
  batch 1C.)
- **Retained and excluded participant rosters are partitioned per group from the observed
  IDs, and a valid cutoff projection beyond the observed trial counts stays numeric and
  labeled.** The roster of a group is now the union of the IDs observed in that group in any
  event, and the excluded IDs are that union minus the participants retained in every
  required event; before, a participant missing from the first event could disappear from the
  excluded list, and an ID reused across groups could pool its support across them, so a
  cutoff above either group's own maximum was not flagged. Per-cell support is recorded from
  the original observed counts before any retention filter. When a valid reliability cutoff
  lies beyond every participant's observed trials, the projected coefficient is kept as a
  number and labeled a projection against original counts, where before it could be unmarked
  or replaced by -1; when nobody is retained, the ordinary coefficient and trial-count cells
  are dashed while the counts stay numeric, the raw values and the full-sample reference are
  kept in the saved summary, and the criterion export headers and `psyrat_report` say which cells
  are projections or full-sample references; the GUI trial-cutoff and overall tables carry the
  distinction in their labeled or dashed cells. The criterion-score reliability figures, tables and
  dialogs carry the same distinctions and qualify their fitted-model Mu, and the manual's
  zero-retention examples were corrected to match. The Save IDs export keeps its two leading
  columns and appends the group of each ID. Single-session common difference cutoffs do not
  re-filter the event-qualified roster, and the test-retest difference summary keeps its
  full-sample reference calculation, now labeled as such. (G84, G88; batch 2A.)
- **Per-participant coefficients of the gamma two-facet subject-level designs (analyses 25 and 27)
  pool the population highest-order mean remainder into each participant's residual, and a saved
  result acquires that change only through an explicit, versioned reprocessing.** Under both gamma
  parameterizations each participant's residual was the conditional observation-level variance
  alone; the pooled population person x occasion x trial mean-surface remainder was omitted. The
  residual is now that remainder, evaluated at the participant's own dimension value in the
  dynamic analysis 27 and at the population intercept in the static analysis 25, plus the
  participant's own observation variance, with all six population mean components and the existing
  count conventions kept (analysis 25 the mean, analysis 27 the median, of trials per occupied
  occasion; default occasion count 1) and the remainder evaluated as a positive product sum rather
  than a cancellation-prone subtraction; this is a mixed population/person target, not a fully
  person-conditioned decomposition. Per-participant coefficients fall (in the recorded fixture, G
  .808441254 to .799748348 under location-scale and .777698535 to .769650901 under log-nu), so at
  an unchanged cutoff an analysis-25 participant can move from included to excluded, and an
  analysis-27 coefficient can fall below a cutoff it met before (analysis 27 stores no inclusion
  flag). Nothing is reprocessed
  when a saved file is opened: an unstamped 25/27 result keeps its historical observation-only
  scheme and is labeled as such. `psyrat_reprocess_gamma_result` writes the updated result to a
  NEW destination from the saved draws, preserves the original, validates draw, stratum and ID
  alignment, records the conversion so a repeat cannot pool the remainder twice, and refuses
  summary-only files, missing or inconsistent draws and unsupported storage regimes with the
  reason (refit those). Current results carry surface, subject, remainder, count and
  parameterization identifiers in the saved file and the export header; a malformed or incomplete
  marker is refused rather than reinterpreted. See Appendix C of the user manual
  (`documentation/manual/appendix-c_reprocessing.md`). Priors, likelihoods and the sampler's
  posterior draws are unchanged; for analysis 25, whose stored draws are log-SD encodings of the
  residual itself, the migration re-encodes them exactly once into the new file. (PBA-08, Q8;
  batch 2C.)
- **The Estimand selector on gamma dynamic-reliability runs selects between two real, labeled
  surfaces: the population-marginal surface (primary) and the typical-person surface
  (companion).** The selector was inert on those runs. The primary surface adds the population
  observation variance, marginalized over the person effects before each draw's coefficient is
  formed, to the pooled mean-surface term (person x trial for one facet, the highest-order
  remainder for two); the companion sets the person effects to zero (location-scale: the
  person log-SD effect w = 0; log-nu: both u and v = 0) with the same population mean
  components. Neither is a mean of per-participant coefficients, and for log-nu the two
  observation variances stand in the ratio exp(Var(2u - v) / 2), so they coincide only when
  that variance is zero. The selected target and the surface trial-count override change the
  surface only; the per-participant coefficient table is fixed. The target is named in the
  viewer note, the saved figure caption (a graphics annotation, so it survives a PNG or PDF
  save), the saved result, `psyrat_report` and every export header, and a result saved under
  the historical `marginal = 0` keeps its old population-surface meaning. The selector stays
  disabled on the difference variants, where only the typical-person reference surface is
  computed. Changing the target on a saved result is a re-summary, not a new fit. (GUI-1, Q9;
  batch 2C.)
- **Analysis 8's per-participant difference-score rows are now exported (G107).** On the
  subject-level difference design the difference-score flag decides retention and the Save IDs
  lists, yet neither subject export carried those rows. Save Subject Coefficients now appends one
  labeled section per group after the unchanged event sections, and `psyrat_report` returns
  `tables.sscoeffs_diff`, unrounded. Each row is one participant's first-event-minus-second-event
  score at their own two trial counts: the selected coefficient with its interval, the ICC at one
  observation per event (derived, not a Rocha et al. Table 6 entry), the SEM of the selected
  error, the relative error variance of the mean difference (labeled apart from the event rows'
  single-trial error variance) and the difference flag; the header says what each column means.
  Excluded and NaN rows are kept, and the values are the stored summary's, unchanged. A missing
  or malformed stored table stops the export before anything is written. Export a saved result
  again to get the rows; one whose stored summary is refused is summarized again from its saved
  draws, with no refit. Separately, a zero-row table that `psyrat_report` writes to .xlsx, in any
  design, now keeps its header above the column-name row.
- **Analysis 8's subject-level difference-score plots draw the group's difference-score
  coefficient as their reference line, and the red markers follow it (G108).** The line came
  from a single-event formula, which leaves out the residual and trial covariances between the
  events, each event's own trial count and, for dependability, the trial terms; in a check where
  every participant's coefficient was 10/11, it sat at 20/21. It is now the difference-score
  coefficient of the group-level model, already stored in the summary, evaluated at the selected
  mean or median of the retained participants' trial counts in each event (of the full roster
  when nobody is retained); on the ICC plot it is that model's difference-score ICC at one
  observation per event. A caption under the panels states the basis. Participants are drawn in
  red when their interval lies wholly above or below this line, so the set of red markers can change.
  Each participant's coefficient and interval are unchanged, and the line is not the overall
  table's mean of the retained participants' coefficients. Regenerate the figures from the saved
  result to pick this up; no refit. A figure drawn directly from a summary stored before this
  change still uses the stored coefficient, and its caption says the count basis was not
  recorded; re-summarizing adds it.
- **On analysis 7, Plot: Trials vs Reliability follows the Difference-Score Coefficient the
  viewer shows (RC-31 residual).** That control already set the summary and its tables, but the
  per-event curves were drawn at the stored G-Theory Coefficient preference, which this viewer
  does not show, so a visible Dependability choice could draw a Generalizability curve and the
  reverse. The curves now use the coefficient resolved from the visible control. Plain one-facet
  runs keep their own G-Theory Coefficient, and the two selectors of the test-retest
  difference-score viewer are unchanged. Regenerate the figures from the saved result to pick
  this up; no refit.
- **A viewer's .xlsx save replaces the chosen workbook instead of writing into it, so no cells of
  the old file survive (G110).** On an existing workbook the viewers' .xlsx saves changed only the
  cells they wrote, so rows below the new table, columns beyond it and extra worksheets of the
  older file stayed in the saved file. Every .xlsx save from a viewer (table windows, Save IDs,
  the criterion table, and the dynamic-reliability and data-splits screens) now writes a new
  workbook in the same folder, reads it back, requires one sheet holding exactly the intended
  cells, and only then moves it onto the chosen name. An existing file is first copied and the
  copy checked byte for byte; the moved file is read back, and if the move failed and changed the
  original, the original is restored from the copy and checked, or the copy is kept and named
  when that fails. Dialogs describe the cleanup as observed and name any leftover file; a
  leftover after a successful save is warned about without calling the save failed. A folder
  name is refused and Cancel writes nothing. The replacement is not atomic and does not keep the
  old file's permissions or other file-system metadata. Some text is refused rather than written:
  text longer than the 32,767 characters an Excel cell holds; a carriage return, which the
  current workbook writer stores in a form standard readers return as a line feed; and some
  rarer text the check cannot confirm (manual Section 9.1 lists it). An infinite value is refused
  rather than written as a large number; NaN is still an empty cell. Long design lines, such as
  those of large multi-occasion data-splits exports, and text with line breaks can be saved as
  .csv instead; for other refused text, .csv can be tried, but check the exported text before
  relying on it.
  `psyrat_report` is not part of this change: it still deletes an existing .xlsx before writing
  and keeps no copy. A workbook that an earlier version saved over an existing file may hold stale
  cells; export it again from the saved result, with no refit.

### Minor changes

- **Trial-cutoff table headers fit in both single-session and test-retest windows.**
  The cutoff column now has a fixed width, with wider ordinary and projection windows
  that also leave room for a vertical scrollbar. Trial/split labels, coefficient values,
  projection annotations and exported table contents are unchanged. (G115.)
- **Saved-result reprocessing reports the specific refusal for in-place migration.**
  Choosing the source file as the output now raises `output_is_source` before the
  general `output_exists` check; another occupied destination still raises
  `output_exists`. Both refusals leave existing files untouched. (G116.)
- **Returning to the home screen no longer repeats the dependency and release checks.** The
  version checks (two or three `stanc --version` spawns, one behind a whole-second wait) and the
  GitHub release check (an 8 s timeout when offline) now run once per MATLAB session; a relaunch
  from a Back button or a cancelled dialog reuses the verdicts and prints one line saying so,
  and the relaunch after the installer or the updater always rechecks. The presence gate stays
  live, a launch that fails it reuses nothing, and the reuse is keyed to the toolbox, CmdStan,
  MatlabStan and MatlabProcessManager actually found, so a change in which of those the toolbox
  finds, a hand-installed dependency included, reruns the checks automatically and says why.
  `psyrat_startgate_cache('clear')` then `psyrat_start` repeats them on demand; restarting
  MATLAB always does. The update prompts and every fit-time check are unchanged. (G91.)
- **The difference-of-differences trial-count heatmaps no longer clip.** A cell whose projected
  coefficient is inadmissible (denominator <= 0) or outside [0,1] is masked in grey and named in a
  caption, and the color axis is fixed at [0,1]. Before, MATLAB's two-input `min`/`max` dropped
  the NaN of an inadmissible cell and painted it 1.0, and out-of-range values were painted 0 or 1.
  The coefficient table beside the heatmaps is draw-level and was never clipped. (PBA-07, Q2.)
- **The difference-of-differences heatmaps also mask a cell whose variance components are
  inadmissible, and name which summary each display uses.** A cell whose plug-in components
  (universe score, scaled between-trial, scaled residual, single-observation between-trial
  contrast) are negative or not finite is now grey with its raw values and the reason kept, even
  when the resulting ratio happened to fall inside [0,1] (an unequal-count example with a
  between-trial component of -.98 and a dependability of .330033 was painted valid). The window
  and the saved figure metadata name the heatmap as a plug-in ratio of posterior mean components
  with no interval and the table as the posterior mean of draw-level coefficients with its
  interval (two summaries of one posterior that can differ: in the recorded two-draw example the
  table reports .7 and the heatmap .833333), and the table's negative-component warning survives a
  FIG or image save instead of living only in a popup. The arithmetic of every admissible cell is
  unchanged (41,718 grid values bit-identical before and after), and the manual's heatmap image
  was regenerated from its saved result. Regenerate DoD displays and saved figures from saved
  results to pick this up; no refit. (G11, the PBA-07 residual; batch 2B.)
- **The convergence check names the parameters it excluded as constants.** Rows whose R-hat and
  effective sample size are both `nan` (correlation-matrix diagonals fixed at 1, disabled
  placeholders fixed at 0) are excluded from the gate; the result now records their names, the
  results viewer's diagnostics notice lists them, and export headers gain `Recorded diagnostic
  exclusions: N (...)` when any were. The source comment and test docstring that said chains stuck
  at different constants print an Inf R-hat were corrected: the tested CmdStan prints `nan` for
  both diagnostics in that case, so such a parameter would be excluded like a constant; the
  disclosure makes that visible. (PBA-02, Q6a and Q6b.)
- **A native fitlme run's intervals are labeled bootstrap intervals, and its saved settings
  record which of them were applied.** The viewer text, the plot-axis titles, the export header
  and `psyrat_report`'s interval note all said "credible interval" whatever engine produced the
  result, although the REML path's interval is a parametric bootstrap. They now follow the saved
  engine; a result saved before the engine was recorded keeps the Bayesian label. The reliability
  trials plots also print the requested interval width instead of a hardcoded 95%, the export
  header's MCMC line reads "not applicable (REML parametric bootstrap)" on a REML run, the
  runconfig sidecar records `priors_applied` and `mcmc_settings_applied` as false there (the
  requested configuration is still kept for replay), and the saved convergence record's
  threshold and scope fields publish the bootstrap rule actually applied (original fit valid, per
  fit, at least 90% and at least 200 successful replicates when 200 or more were requested)
  instead of the R-hat and ESS thresholds; the internal R-hat and ESS compatibility row remains
  and is not a measured diagnostic. Estimates are unchanged. (PBA-05, Q3; batch 1A.)
- **The native-engine notices no longer call a CmdStan run definitive.** They advised re-running
  with CmdStan "for a definitive analysis"; a CmdStan run has its own diagnostic and uncertainty
  requirements, so they now say that final-inference readiness is unassessed and that a CmdStan
  run also requires diagnostic and uncertainty validation. The native-versus-CmdStan differences
  they describe are unchanged. (Batch 1A.)
- **The user manual gains an Acknowledgments section** in its front matter: the beta testers,
  NIMH support under award number R01MH128208 with the NIH disclaimer, and the AI tools used to
  write code and documentation.
- **The README carries the edits made on the public repository after 0.1.0-beta, and its
  validation section moves to `documentation/validation_and_testing.md`.** The description names
  psychophysiological measurements (ERP component scores, oscillatory EEG, skin conductance,
  cardiac data) rather than ERP scores alone, and `CITATION.cff`'s abstract states the same
  scope; the citation guidance says "for now"; the bug-report paragraph asks what was expected;
  the README gains the NIMH Acknowledgments paragraph. The four suites and the evidence-tier
  table are now on their own page, which the manual's tutorials and the design-maturity table
  below cite; its JUnit note names the three runners that take the option (the workflow sweep
  writes its own CSV and MAT report instead). The copyright line keeps 2016-2026. (G105.)
- **The toolbox names itself PsyRAT, not "PsyRAT Toolbox".** The README heading (edited on the
  public repository after 0.1.0-beta), the startup, installer and updater messages and dialogs,
  the command-line input-check errors, the View Results file filter, `psyrat_printversions`, the
  function help text and the manual now say "PsyRAT", and the updater installs a new release into
  a folder named `PsyRAT_v<version>` (was `PsyRAT_Toolbox_v<version>`). Where an export header
  names the toolbox and its version, in the GUI exports and in `psyrat_report`'s files, that line
  now reads `PsyRAT v` and the version, which is chosen as before (`PsyRAT v0.1.1-beta` for a
  result processed with this release). The header keeps its lines and their order, and nothing in
  PsyRAT depends on the wording of that line; a script of your own that matches the old text needs
  the new one. (G127.)
- **The test-retest Output Prefs screen's controls fit the window.** The right-column fields
  and the Save and Back buttons took their widths and positions from the window height rather
  than its width, so on that taller screen they clipped and the Back button ran past the right
  edge.
- **First-party copyright notices read 2016-2026** in every file that carries one: the root
  scripts and runners, every `.m` file under `subroutines/`, the test files that carry a
  notice, and the README. Bundled third-party notices, the GPL text and reference-paper dates
  are untouched. (Q10, G18.)
- **Setting a facet popup back to `(none)` after using its Select button no longer crashes
  Analyze.** Opening Select Groups, Select Events or Select Occasions under a real column (or a
  failed Analyze under one) remembered that column; setting the popup back to `(none)` and
  clicking Analyze then read the sentinel as a table variable and threw "Unrecognized table
  variable name '(none)'". The facet is now simply not used. (G92, beta-tester report.)
- **Multi-measurement runs are named `<stem>_<measurement>.psyrat`.** On macOS the save dialog
  supplies the `.psyrat` extension, and the namer kept it inside the stem, so two measurements
  saved as `<stem>.psyrat_<measurement>.psyrat`. The sidecar and export stems follow. (G93,
  beta-tester report.)
- **The per-measurement error report is named `<stem>_error.txt`** (it was `<stem>.psyrat.txt`).
  The `_error` suffix keeps the report out of the input folder's `.txt` namespace, so a failed run
  can never overwrite an input file of the same stem. (G94.)
- **The Output Prefs screen refuses an empty list selection.** A Cmd or Ctrl click on the
  selected row of one of its three lists cleared the selection on a beta tester's installation
  (MATLAB release not recorded); saving that stored an empty preference, Generate then failed
  with "depcentmeas should be 1 (mean) or 2 (median)", and the reopened screen showed nothing
  selected. Save now raises a dialog naming the list and puts the previous choice back, on both
  the single-session and the test-retest screens. (G95, beta-tester report.)
- **The missing-cells preflight message names the cells, and the processing screen reports every
  preflight error.** "Input data include missing cells" now continues with each column that has
  missing cells, how many, the first row of the prepared table and that row's id, and says that
  a value which cannot be read as a number (NA, text, an empty cell) is imported as missing. When a file trips several preflight
  errors, the Analyze route reports all of them instead of the first alone, so an Inf score is no
  longer hidden behind a missing cell. (G96, beta-tester report.)
- **The overall summary's difference-score rows print n Excluded.** The difference-score row
  dashed that cell while the per-event rows beside it printed the complement; the diff row now
  prints the roster minus the retained count too. Manual section 12.4 names the window, table and columns as they
  appear on screen and marks its True D column as computed by hand. (G97, beta-tester report.)
- **The overall summary's Label column has a fixed width** so labels such as `GroupA - inc_err`
  render whole; the window is wider to hold it. Windows stay fixed-size. (G98, beta-tester
  report.)
- **An empty selection can no longer be saved on the Select Groups, Select Events and Select
  Occasions screens.** Saving with every level deselected stored an empty selection that the
  loader read as "no filter", so every level was processed while the button's tooltip claimed all
  were selected and the reopened list showed none. Each screen now refuses with its own modal
  dialog (No Groups Selected, No Events Selected, No Occasions Selected) and stays open; Back
  leaves the current selection as it was. (F1, a follow-up to the beta-tester fixes.)
- **Confirming the scoring-method or scoring-spec picker with nothing selected no longer
  throws.** MATLAB's list dialog returns an empty index on OK when its selection was cleared; the
  import and scoring screen treated only Cancel, so the next line failed on the console and the
  button looked dead. An empty OK now returns like Cancel. (F2.)
- **Cancelling a table or ID save no longer runs on silently, and a name outside .xlsx/.csv is
  refused.** The Save Table and Save IDs buttons of the overall, trial-cutoff, sources-of-variance
  and all-standard-deviations windows continued after a cancelled dialog and wrote nothing without
  saying so, and a typed extension the dialog's filter did not supply fell through the same way.
  They now return on Cancel and raise a modal Unsupported File Type dialog naming the file for any
  other extension; the two difference-of-differences table windows gain the same dialog, and the
  extension test is case-insensitive at every site. The one exception is the View All Standard
  Deviations window of a one-facet difference-score run, which offers only CSV and adds .csv to
  any other name (G109, below). (F3.)
- **The installer's copies of the bundled dependents keep their full version in the folder
  name.** MatlabProcessManager-0.5.1 was copied as MatlabProcessManager-0.5 (the last version
  component was dropped as if it were a file extension). Every lookup is by content, so nothing
  else changed, and existing installs are left as they are. (F4.)
- **Fixed-width table windows now hold their columns.** The one-facet sources-of-variance table,
  the all-standard-deviations and difference-components windows, the observed
  difference-of-differences table, the formula guide, the overall reliability tables (which now
  carry a width for every column, so Max # of Trials no longer clips behind a horizontal
  scrollbar) and the test-retest table's Estimate New Reliability Coefficient window are sized
  for their columns plus grid slack and, where the row count can call for one, a vertical
  scrollbar. (F5.)
- **Headers on the nine fixed-width table windows now read whole.** Long headers there
  (SEM (dependability), Max # of Trials, Within Person (Error), the difference-components and
  observed difference-of-differences headers, among others) were cut with an ellipsis inside
  their columns. They now break onto two lines on screen; the CSV and XLSX headers keep the
  one-line names, and the observed difference-of-differences table's four mean-trial columns
  read Mean Trials / ERP 1 to ERP 4 instead of four identical cuts. A table window wider than
  the primary display now opens 80 px narrower than that display, is moved on screen, and
  scrolls horizontally instead of running off-screen. (G104.)
- **The two-facet dynamic-reliability viewer's Reliability coefficient popup fits in its
  window.** The popup ended past the right edge of the 560 px window, and its longest option was
  cut short. Two-facet dynamic-reliability viewers now open 660 px wide with a 300 px popup,
  which shows the longest option at the default font; one-facet viewers are unchanged. (G106.)
- **The per-participant window of the gamma subject-level dynamic designs shows its whole
  target note.** On analyses 26 and 27 the note under the per-participant table names the
  conditioning set, the residual scheme and the count sources. It wraps to six to nine lines,
  depending on the design, the scheme and the font size, but its fixed box showed about three and
  a half at the default font, so the remainder-scheme identifier and the count sources were cut
  off on screen. The
  export headers were always complete. The box is now sized from the note's own line count and
  the window grows to hold it, so the table keeps its size; windows without the note are
  unchanged. (G111.)
- **A result without a recorded measurement name opens in the viewers.** The single-session,
  test-retest, difference-of-differences and criterion setup viewers print the name of the
  measurement column. A result that lacks it stopped with an "Unrecognized field name" error while
  the viewer was being built, leaving a partial window, and a result with an empty name showed a
  blank Measurement line. Shipped fitting always records the name, so only a hand-built or
  externally produced result could lack it. Such a result now opens with the measurement shown as
  "(not recorded)", the label the loader already gives a missing version stamp; results that
  record a measurement name are unchanged. (G112.)
- **Three CSV-only saves ask before replacing a file the save dialog never named, and write
  through a temporary file (G109).** Save Subject Coefficients, the Save Table of the one-facet
  difference-score View All Standard Deviations window and the criterion table's CSV branch add
  .csv to a name without it; an existing file under that name is replaced only after a question
  that names it, and a folder of that name is refused. Each writes a temporary file in the
  target's folder and moves it onto the target only after it was written and closed without a
  reported error. A failure before the move replaces nothing and says whether the temporary file
  was removed; a failed move is reported without claiming the old file is unchanged. In the
  Save Subject Coefficients and View All files, IDs, group and event names and labels that
  contain a comma, double quote or line break are quoted, so each stays one field. Header lines
  of the viewers' CSV exports that contain a carriage return or line feed are now quoted the same
  way, as lines with a comma or double quote already were, and so is such a label in the observed
  difference-of-differences table's CSV. Written values are unchanged.
- **An output name the save folder cannot hold is refused before the fit instead of failing after
  it.** Before estimation, the GUI and `psyrat_run` now check, on the actual disk, that a file can
  be created in the save folder under the output's name, its `_runconfig.json` sidecar's and its
  `_error.txt` report's. How long a name may be depends on the disk, not on a fixed number of
  characters, and names are never altered, so accented and other non-English characters stay as
  typed. A name the check finds the folder cannot hold, or a folder already carrying the output's
  or sidecar's name, stops that request before its fit with a `psyrat_data:savename` error naming
  the file and the remedy; in a batch only that measurement is refused and the others run. Before,
  a long measurement header in a batch, or a long `psyrat_run` savename, ran the whole fit and
  then failed at the save. To decide, the check briefly creates and removes a test folder with a
  name of the same length. If that folder is refused, it tries a test file of the same length,
  created only where no file of that name exists, and refuses the name only if that file fails
  too while a short test file made the same way succeeds. Existing files are not touched, fit
  seeds are unaffected, and a test entry that cannot be removed is named on the console and never
  deleted with content in it. When the check cannot decide (for example, no Java runtime after a
  test folder is refused), it refuses nothing and the run keeps its usual handling of a failed
  save. Passing the check does not guarantee the later save. The Processing Complete box names the
  error reports actually written (the first five, with a count of the rest; the console lists
  them all). No estimate, prior, seed or saved-result format changes. (Queue row 31, LA-2.)

---

## 0.1.0-beta — 2026-09-09

First version to carry an explicit maturity label. `-beta` states that the API and the shape of the
outputs are still moving. It is not a statement that the estimates are provisional: see the design
maturity table below, which is per-design.

### Design maturity

These labels describe **verification status in this repository**, not a judgment about the
underlying statistics. Read them together with the evidence-tier table in
`documentation/validation_and_testing.md`, which explains what a green run of each suite does and does not prove.

| Maturity | Designs | What the label means |
|---|---|---|
| **Supported** | One-facet and test-retest reliability (Gaussian, and gamma where offered); subject-level error variances (except analysis 25 under gamma, Restricted below); two-event difference scores, non-concurrent; the four-event difference-of-differences (analysis 9, non-concurrent by construction); reliability from splits (Gaussian only); the dynamic/conditional family except analysis 20 under gamma and except the gamma two-facet arms of analyses 14 and 27 (Restricted below); the person-specific dynamic non-concurrent DoD designs (analyses 28/29, single-occasion and retest, in both likelihood families - the Gaussian arms shipped 2026-08-20; gamma is location-scale only) | Covered by the unit and accuracy suites, and by live CmdStan recovery against known variance components where a recovery test exists (the gaps are listed below the table). |
| **Restricted** (2026-09-21) | Three gamma two-facet designs: dynamic test-retest (analysis 14), subject-level test-retest (analysis 25) and subject-level dynamic test-retest (analysis 27), under gamma in both parameterizations (location-scale and log-nu). The Gaussian arms of the same designs are Supported. | **The coefficients are usable with a caveat; the component estimates are not recovered at level with few occasions.** Each reviewed capture of these designs had three or four occasion levels. As variance-scale ratios of posterior mean to planted value, the five re-evaluated scope captures put the occasion component at 2.4 to 12.1 times, the highest-order remainder at 2.2 to 7.2 times and the pooled residual at 1.1 to 1.9 times, with several 95% intervals excluding the truth; the retained analysis-27 location-scale capture, reviewed as a control, read 1.42, 1.66 and 1.02, which is not a new pass. Two prospective analysis-27 fits: occasion component 60.1 and 2.45 times, highest-order remainder 31.6 and 3.35 times, pooled residual 1.83 and 2.95 times, the second of them non-converged by the toolbox's own criterion. The pattern is consistent with raw-scale components that are exponentials of squared log-scale SDs whose posteriors are right-skewed with few occasions, which makes their means tail-dominated. G and D stay within about .03 and .06 of the truth for analyses 14 and 25; for analysis 27 both G and D are biased downward in the converged fits (medians across participants of the posterior-mean coefficient minus its planted value: G about -.05 location-scale at four occasions and -.10 log-nu at three, D -.13 in the converged prospective fit, -.10 log-nu and -.02 in the retained capture; the non-converged second prospective fit gave D +.08) while participant ordering is kept (r >= .86). Wherever production was compared with an independent enumeration it agreed to 2e-12 or better, which shows implementation consistency, not recovery, calibration, identifiability or sampling validation. Disclosed at Analyze (a notice to acknowledge before the save dialog), in the results viewer, on the saved surface and subject figures and in every export header. Report the coefficients with the caveat; reading the component estimates against their intervals does not remedy the weak evidence. The predeclared LS27 recovery assertion keeps its recorded FAIL as the adjudicated outcome (release checklist section 6). |
| **Experimental** | Concurrent (correlated-residual) difference scores under the Gamma family with per-person ν (`dispersion=1`; analyses 7, 8, 10; the only gamma arm analysis 8 has) | **Convergence is not validated in this repository.** A live subject-level run reached R̂ ≈ 2.5 with ESS ≈ 2 on one fit. These designs need long warmup and, in practice, threaded compute. The runtime convergence gate (R̂ > 1.1) warns you about a specific fit, and the three live-recovery tests *skip* rather than pass when a fit has not converged. Treat any estimate from these designs as unverified until you have checked its diagnostics yourself. |
| **Tractability checked, robustness not** | The group-level concurrent difference designs under global fixed dispersion (`dispersion=2`, the default for analyses 7 and 10; analysis 8 has no such arm, since global dispersion is a group-level estimand) | One simulated known-truth dataset converged cleanly (max R̂ = 1.000; recovered G = 0.941 against a true 0.939). That establishes tractability at one set of sampler settings on one dataset. It is **not** a robustness result: there is no sweep over sample size, trial count, or real data. |
| **Development-only** | Analysis 20 (two-facet subject-level concurrent dynamic difference) under Gamma, location-scale | Ships and runs, but no empirical anchor exists for its two-facet cut-copula configuration and no user-facing label says so on screen or in an export; treat its estimates as unverified. |
| **Not implemented** | Analyses 15/16, 17/18 and 21/22 under Gamma; every splits design under Gamma; the exact-shared covariance count arm (not user-reachable); concurrent static DoD ("9 + rescor") | Barred pending an estimand decision. Requesting them errors rather than silently substituting a different model. |

Validation coverage behind the Supported label is not uniform, and the gaps are recorded here
rather than closed for this release (owner ruling 2026-09-05): the Gaussian four-event
difference-of-differences (analysis 9) has a live planted-truth recovery only through the native
HMC engine, not through CmdStan; the Gaussian subject-level dynamic designs (analyses 26 and 27)
share the fits of analyses 11 and 14 and have no recovery test of their own; the gamma arm of
analysis 19 is checked only by a replay of an external reference fit that runs no sampler; and
the splits designs (analyses 23 and 24) have no independent-oracle row. The CmdStan recoveries
for those arms are filed as post-beta work.

The dispersion parameterization actually used is recorded in the provenance block of every copula
export, so a saved table states which of the two rows above it belongs to.

### Changed

- **The Gaussian difference-of-differences priors are configurable.** Analysis 9 was the one
  design with every prior fixed inside its Stan builder and its native HMC log-posterior, so a
  `'priors'` override never reached it (the subject-level family still keeps its intercept,
  between-person and residual log-SD constants fixed and exposes only `priors.sserr.sig_trl`). They are now `priors.dod` (`b_cell`, `b_sigma_cell`,
  `sd_id`, `sd_trial`, `lkj`; owner ruling 2026-09-05) at the former constants, so the generated
  Stan and the native fit are unchanged unless a value is set. A result saved before the family
  existed still prints `Priors: PsyRAT defaults`, because a family absent from a stored prior set
  is at its default by construction.
- **Exported table headers now write the reliability cutoff with four decimal places instead of
  two** (now `Dependability Cutoff: 0.8000`; previously `Dependability Cutoff: 0.80`). **This changes the
  header of every table that carries the cutoff line** (the difference-of-differences and variance
  tables carry none). The cutoff is a recorded analytic *input*: it is
  consumed at full precision to pick the minimum trial count that defines the retained sample, and
  the cutoff field accepts any value in (0, 1). At two decimals a cutoff carrying a third or fourth
  decimal — as a sensitivity grid stepping .025 produces — was written to file as a different
  threshold than the one applied; a run at 0.9999 recorded `1.00`. Estimates in the table *body*
  (coefficients, credible intervals, SEM) are unchanged at two decimals; this is the input/estimate
  distinction, not a change to the reporting convention. Four decimals matches the format the
  toolbox already used for its other user-typed threshold (`Criterion Cutoff`). Filed as B36.
- **Version string is `0.1.0-beta`** (was `0.1.0`, which understated the toolbox against a Version
  History reaching 0.5.2). Filed as B4.
- **Log-nu subject-level gamma residuals now pool the population person-by-trial term (analyses 6,
  8, and 26 under `gammascale = 1`, the "Log-nu (relative dispersion / CV)" popup or
  `'gammascale', 1`).** Each participant's residual was the conditional observation-level variance
  alone; it omitted the induced person-by-trial term of the log-linked mean surface that the
  group-level residual carries and that the location-scale parameterization (the default since
  2026-08-07) already pooled. Every per-participant coefficient on the log-nu path was therefore
  optimistically biased, by an amount that depends on the person and trial log-SDs (about 1% to
  about 50% of the residual across the design points examined). Per-participant coefficients from
  those runs change, and so do the affected runs' own group-level residual, within-person SD, and
  ICC rows, which read the same pooled population reference; runs of other designs and every
  location-scale run are unchanged. The export header of an affected run now carries a "POOLED"
  residual line. Filed as S20.

### Added

- **Gaussian arms for the person-specific dynamic non-concurrent difference-of-differences
  (analyses 28/29).** Until 2026-08-20 these two designs were the only routed analyses with no
  Gaussian implementation ("no Gaussian model of the estimand has been derived"); an owner ruling
  reopened that closed decision and the arms shipped with their own derivation, recorded in the
  project's internal formula audit. Identity link on the mean (negative ERP amplitudes
  are in scope, unlike the gamma arm), log link on the person-specific residual SD, the same
  nonconcurrent counterfactual estimand and disclosures, and the same read-outs minus the
  chi-square-specific per-person nu columns. CmdStan-only for now (`native_key = ''`), like
  analyses 25-27. Validation is in-repo by design: no external reference bundle exists for a
  Gaussian DoD, so the live planted-truth recovery tests for both analyses are the load-bearing
  gates, alongside closed-form longhand checks, the identity-link theorem pins, and a full-space
  tie to the static Gaussian DoD kernel.
- **`CITATION.cff`** at the repository root. Previously the only DOI in `README.md` belonged to the
  2017 ERA Toolbox paper, so a user following the README to cite PsyRAT would have cited a
  different tool. The citation file lists Rocha et al. (2026), the source of the estimators the
  toolbox implements, first under `references` and asks for the software version alongside it; it
  deliberately carries no `preferred-citation` key, because in CFF 1.2.0 that key makes citation
  tooling emit the paper instead of the software. It also lists Rast & Clayson (in press) and the
  2021 ERA Toolbox test-retest paper under `references`, and carries an inline note to revisit that
  decision when the dedicated PsyRAT software publication is released. A matching "How to cite" section was added to `README.md`.
- **A user manual**, `documentation/manual/`: 18 chapters plus a glossary and a reference list,
  written for this toolbox rather than inherited from the ERA Toolbox. Four worked tutorials run on
  the simulated datasets in `test_data/`, and every reference value they print came from a real run
  at the settings its chapter states, transcribed with provenance in
  `documentation/manual/tutorial_expected_values.md`. The same content ships as one generated PDF,
  `documentation/manual/psyrat_manual.pdf`, rebuilt from the chapters at each release.
- This changelog.

### Fixed

- **Subject-level plots on an events-only or groups-only run were untitled.** Every panel of the
  caterpillar plot fell through to an empty title, so a reader could not tell which panel belonged
  to which event; the panels now carry the event or group name. Groups-by-events runs were already
  titled "group: event" and are unchanged.
- **Closing the "Chains did not converge" dialog with the window's close box errored after the fit
  had finished.** It now counts as Do Not Rerun: the non-converged result is saved as it is when
  that button is clicked (council G13).
- **Save Table wrote a half-integer median trial count in scientific notation** (`2.450000e+01`) in
  both the single-session and test-retest overall tables. It now prints as a decimal; integer
  medians are byte-identical to before.
- **The input checks' non-fatal warnings were computed and never shown.** Fewer than 20
  participants, a median trial count below 10 in a cell, and an incomplete event-by-occasion layout
  now print in the Command Window before the fit (the only channel for `psyrat_run`) and, in the
  GUI, appear in one modal dialog titled Preflight warnings that waits for OK (council G8).

- **Every table export and the headless report failed to load** on the pre-release tip: the S20
  disclosure commit left a stray statement (`, added 2026-08-04 (S16).`) in
  `psyrat_provenance_lines.m`, which every table writer and `psyrat_report` call, so the file did
  not parse. Found by three independent review lanes of the 2026-09-05 pre-beta audit; removed.
- **The convergence gate now fails closed on real pathologies without rejecting constants.** A NaN
  R-hat or effective sample size on a parameter that varies, and a convergence table with no
  monitored parameter, both counted as converged, because `NaN >= 1.1` is false and an empty
  comparison flags nothing; both now read as not converged, as does an Inf R-hat (chains stuck at
  different constants). Parameters that are constant by construction (`stansummary` prints NaN for
  both statistics; correlation-matrix diagonals fixed at 1, disabled covariance and
  residual-correlation terms fixed at 0) are excluded from the gate, which the old comparisons did
  by accident: an earlier form of this change flagged them, and no design that monitors one could
  converge. Finite values are judged exactly as before.
- **Preflight refuses three inputs that silently changed the estimand or the exported numbers**, on
  both the GUI and the scripted route: a non-finite (Inf) score (`meas:notfinite`); a label
  containing a comma, double quote, or line break, which the exported CSV tables would have written
  unquoted so every later column of that row shifted (`labels:reservedcharacters`); and an occasion
  column with a single level, whose presence routed the run to the two-facet test-retest design with
  an unidentified occasion facet (`preflight:singleOccasionLevel`). Three new warnings, not errors:
  a participant id listed under more than one group, which the estimation core counts as one person
  per group (`preflight:groupWithinParticipant`; a warning rather than a refusal because ids
  numbered within each group are a legitimate convention, and the message says when a within-person
  condition belongs in the event column instead); labels that differ only in letter case or
  surrounding spaces (`preflight:labelCaseVariants`); and exact duplicate rows
  (`preflight:duplicateRows`). The gamma concurrent designs are now refused before CmdStan compiles
  when the Statistics and Machine Learning Toolbox is missing
  (`preflight:statisticsToolboxRequired`) instead of failing after sampling.
- **The sidecar's input-file hash is taken from the file's full path.** Both entry routes record
  `raw.sourcepath`; hashing from the bare file name only worked when the file sat in the current
  folder or on the MATLAB path, and could hash a same-named file found there. The sidecar gains
  `data.source_path` and `estimation.threads_per_chain` (the resolved within-chain thread count,
  1 when threading was off; empty only for a result saved before the field existed), and a
  threaded run's export header now says so, because threaded models use `reduce_sum` and are not
  bit-reproducible across runs.
- **`psyrat_relsummary` rejects an unknown `analysis` key** with `varargin:analysis` and prints the
  mapping from stored `rel.analysis` values to dispatch keys; an unknown key used to fall through
  both dispatch switches and fail with MATLAB's generic output-argument error.
- **`psyrat_run` reads the input file through the toolbox reader**, gaining the supported-extension
  check and delimiter cycling the GUI already had. A file whose extension the reader does not
  support is now refused with `ext:filetype` before anything runs, where a bare `readtable` used
  to attempt it. The reader's documented `delimiter` option was unreachable because its parser
  looked for the key `idcol`.
- **The installer warns when it falls back to downloading an unpatched upstream MatlabStan or
  MatlabProcessManager** instead of copying the bundled, patched copies; the fallback used to be
  silent until the next launch's patch-level check.
- **The release check no longer calls an official beta "unreleased"**: when only prerelease
  versions have been published it says so.
- `help psyrat_startview` showed a `.mat` example path against a picker that accepts `.psyrat`;
  the dynamic-reliability plot printed a user-supplied dimension column name through the TeX
  interpreter (now literal, like event and group names).
- **Headless exports explain their `-1` cutoff cells.** When no trial count reached the reliability
  threshold, or the cutoff found lay beyond every participant's observed trials, the D-study stored
  `-1` in the cutoff table and the GUI raised a dialog; `psyrat_report` exported the same `-1` cells
  with no explanation. It now returns `report.cutoff_note` and writes that sentence into the header
  of every exported table (council G1).
- **`psyrat_relsummary` accepts an explicit `'interactive'` option.** The flag was read only from the
  saved `.psyrat`, which GUI-produced files never carry, so a batch loop over them could stall on
  the modal threshold-failure dialog with no way to say otherwise. The option wins over the stored
  flag; `psyrat_run` output is unchanged (council G3).
- **The in-app toolbox updater could not have worked.** The Yes button on the "old version"
  prompt fetched the GitHub releases web page and parsed a heading whose markup GitHub no longer
  emits and whose "Version x.y.z" title the publishing runbook never wrote, then downloaded a
  `PsyRAT.zip` no release carried, after every PsyRAT window had already been closed. It now
  resolves the latest stable release through the same GitHub releases API the startup check
  uses, downloads the `PsyRAT.zip` asset every release now attaches, and, when no stable release
  exists or the release has no such asset, reopens the home screen with a dialog naming the
  GitHub page to download from instead of failing (audit G87, owner ruling 2026-09-09). A release
  check that throws is routed into the same dialog. Prereleases are never installed by it.

### Known gaps in this release

Recorded so they are not mistaken for oversights:

- **There is no continuous integration.** All verification is local and manual: the four MATLAB
  suites (unit, accuracy, feature/integration, and the workflow sweep) are run by the maintainer
  before a release. The accuracy suite's regression class (`TestCalculationAccuracyOracle`)
  duplicates the production formulas, so its green run demonstrates "unchanged", not "correct";
  the suite's other two classes, the independent oracle that restates Rocha et al. (2026) Tables
  2, 3 and 6 and the external benchmark, are the corroboration
  (`documentation/validation_and_testing.md`). Live CmdStan
  recovery is likewise run manually by the maintainer, as is the check that re-derives every
  reference value the manual's tutorials print and diffs it against the transcription file, which
  is run before any release that touches estimation or calculation code.
- **The legacy 2020 ERA-era manual is not distributed.** It predates most of the current
  functionality and has no editable source.
- **The GUI requires CmdStan.** `psyrat_start` stops at its dependency check when CmdStan cannot be
  located and offers the guided installer. Without CmdStan the native MATLAB engines (REML via
  `fitlme`, and the lightweight HMC sampler) are reachable only through `psyrat_run`; once CmdStan
  is installed the GUI's Estimation engine popup offers them as well; see the manual's chapter 18.
- **Windows is not in the test matrix, and Linux installation is experimental.** The toolbox is
  developed and verified on macOS; the retired CI lane ran on Ubuntu. The Windows installer checks
  for GNU make and g++ but no Windows run has been verified by the maintainer, and the Linux
  auto-installer announces its experimental status with a warning when it starts. Details in
  `documentation/dependencies_support_matrix.md`.
- **No DOI or archived release yet.** `CITATION.cff` will gain a DOI when a release is tagged and
  archived, and the preferred citation will be revisited if a PsyRAT software publication appears.
- **The updater's download-and-install arm is unexercised live.** It can only run once a stable
  release carrying the `PsyRAT.zip` asset exists, and it is not driven on the maintainer's machine
  because it replaces and saves the MATLAB path. The release resolution and the asset lookup are
  unit-tested against mocked API responses, and the refusal dialog was verified on screen.
