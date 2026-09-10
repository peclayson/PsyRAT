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

## 0.1.0-beta — 2026-09-09

First version to carry an explicit maturity label. `-beta` states that the API and the shape of the
outputs are still moving. It is not a statement that the estimates are provisional: see the design
maturity table below, which is per-design.

### Design maturity

These labels describe **verification status in this repository**, not a judgment about the
underlying statistics. Read them together with the evidence-tier table in `README.md`
(Validation and Testing), which explains what a green run of each suite does and does not prove.

| Maturity | Designs | What the label means |
|---|---|---|
| **Supported** | One-facet and test-retest reliability (Gaussian, and gamma where offered); subject-level error variances; two-event difference scores, non-concurrent; the four-event difference-of-differences (analysis 9, non-concurrent by construction); reliability from splits (Gaussian only); the dynamic/conditional family except analysis 20 under gamma; the person-specific dynamic non-concurrent DoD designs (analyses 28/29, single-occasion and retest, in both likelihood families - the Gaussian arms shipped 2026-08-20; gamma is location-scale only) | Covered by the unit and accuracy suites, and by live CmdStan recovery against known variance components where a recovery test exists (the gaps are listed below the table). |
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
  2, 3 and 6 and the external benchmark, are the corroboration (`README.md`, Validation and
  Testing). Live CmdStan
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
