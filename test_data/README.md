# PsyRAT test and tutorial data

This directory holds the small data files used by the PsyRAT test suite and the
user manual's tutorials. The `manual_ern_*` files are simulated analysis
datasets meant to be opened by readers; everything else is a test fixture.

## Tutorial datasets (simulated; ship in releases)

All five are fully simulated: no real participants are involved, every
generating parameter (means, variance components, seeds) is documented in the
generator headers, and `tests/TestManualDatasets.m` verifies that each
committed CSV is exactly what its generator produces. Free to use and
redistribute under the repository's license (GPL-3).

| File | Generator (`tests/helpers/`) | Design |
|---|---|---|
| `manual_ern_singlesession.csv` | `psyrat_make_manual_singlesession.m` | 80 participants, 2 between-person groups, 4 events (congruency x accuracy), unequal trial counts; columns `id`, `group`, `event`, `meas` |
| `manual_ern_testretest.csv` | `psyrat_make_manual_testretest.m` | 30 participants, 2 occasions, 2 events, balanced trials; columns `id`, `time`, `event`, `meas` |
| `manual_ern_splits.csv` | `psyrat_make_manual_splits.m` | 36 participants x 16 split means with items-per-split weights (4 to 40); columns `id`, `meas`, `weight` |
| `manual_ern_dynrel.csv` | `psyrat_make_manual_dynrel.m` | 80 participants, 2 between-person groups of 40, 2 events, plus one questionnaire-style `worry` score per participant as the dimension covariate; columns `id`, `group`, `event`, `meas`, `worry` |
| `manual_power_concurrent.csv` | `psyrat_make_manual_power_concurrent.m` | 40 participants, single-trial theta-band and alpha-band power from the same error trials (strictly positive, equal counts within participant), plus the same `worry` covariate; columns `id`, `event`, `meas`, `worry` |

The user manual (`documentation/manual/`) walks through these files and prints
expected results produced from them at the toolbox's default seed.

## `legacy_headers_fixture.csv`

Synthetic fixture (generator: `tests/helpers/psyrat_make_legacy_headers_fixture.m`)
whose point is its non-canonical header row, `subjid, group, event, ern`. It
exercises the column auto-detect fallback in `psyrat_startproc` (see
`tests/TestStartprocColumnAutodetect.m`). Not a tutorial dataset.

## `testdata.csv` (development checkouts only)

Legacy long-format fixture (columns `subjid`, `group`, `event`, `ern`) used by
older unit and workflow tests via `PsyRATCmdStanIntegrationFactory`. It is
development-only and is not distributed in public releases; tests that derive
scenario tables from it skip when it is absent. It is not an analysis dataset:
its `group` labels are within-participant task labels, and several
participant x cell combinations are intentionally absent.

## `ep_toolkit_real_subset.mat` (development checkouts only)

A minimal slice derived from a real EP Toolkit single-trial ERP export, used to
test the EP Toolkit `.mat` import/scoring path against genuine recorded data.
Development-only; not distributed in public releases (the consuming test skips
when the file is absent).

- **Extent:** 8 channels, the first 5 single trials, and a -100 to 200 ms window
  (a tiny fraction of the source recording; the full (~260 MB) `.ept` source is
  intentionally not committed).
- **De-identification:** the file contains only numeric EEG values, channel
  names, and time points. It carries no participant name, date of birth, or other
  identifying information.
- **Regeneration:** produced by
  `tests/helpers/psyrat_make_ep_real_subset_fixture.m`, which requires a local
  copy of the source EP export (not distributed). The generator has no default
  source path and errors if one is not supplied.

Redistribution of this derived slice is governed by the data-sharing terms of the
originating study. Confirm the applicable consent / IRB approval and data-use
terms before reusing it outside this test suite.
