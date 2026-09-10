# Tutorial expected values

This file is the transcription record for every expected value the manual
states. Each slug below appears in exactly one chapter as an `EV` placeholder;
the value recorded here is the value the chapter prints. Values come only from
real producer runs on the committed `test_data/manual_ern_*.csv` files at the
settings each section records, on the capture machine (the expected-values
protocol in the production ledger governs). Slugs are written bare here, without
the placeholder wrapper, so that a count of live placeholders in the chapters
stays exact.

A blank Value means the producer run has not happened yet. Do not fill a value
from memory, from a test fixture, or from a run at different settings.

Reading a value: reliability coefficients are reported to two decimals with no
leading zero, standard deviations and SEMs in microvolts to two decimals, and
credible intervals as bracketed pairs, matching the chapters' checkpoint
policy (exact reference values, a "within about .02" acceptance band, and the
troubleshooting chapter's section 18.3 for what to do when a checkpoint
misses).

## Provenance (fill once per capture session)

| Field | Value |
|---|---|
| Capture date | 2026-09-03 (regenerated datasets; the 2026-09-02 values were superseded) |
| Toolbox version | the version printed in the startup banner at capture (the literal is single-sourced in index.md and not repeated here; the production ledger records the development revision) |
| CmdStan version | 2.38.0 |
| MATLAB version | R2025b Update 5 (25.2.0.3177638) |
| Estimation engine | cmdstan |
| Operating system | macOS 26.3.1 (Apple silicon) |

## Run QS: quick start (chapter 2)

Dataset `test_data/manual_ern_singlesession.csv`; event `inc_err` only, both
groups; all processing preferences at defaults (4 chains, 5,000 warmup, 5,000
sampling, seed 12345); the viewer left at its default cutoff (.80), at which
neither group retains anyone, so the coefficient is evaluated over all
participants.

| Slug | Value |
|---|---|
| quickstart_dep_pt | .72 |
| quickstart_dep_ci | [.56, .83] |
| quickstart_nincl80_grpA | 0 |
| quickstart_cut80_grpA_n | 38 |
| quickstart_cut80_grpB_n | 55 |

## Observed timings (chapter 3)

Wall-clock observations on the capture machine, recorded once.

| Slug | Value |
|---|---|
| first_compile_time | 8 s per generated model (measured 7.5 s wall for stanc + C++ compile + link on the capture machine; a four-event run compiles four models) |

## Import walkthrough (chapter 8)

Read from the Trial Browser on the example `.set` file the chapter walks
through (the file is local to the author and is not committed).

| Slug | Value |
|---|---|
| import-ern-nevents | 2 |
| import-ern-ntrials | 398 |
| import-ern-nincluded | 398 |
| import-ern-nexcluded | 0 |
| import-ern-trial1-value | 3.41604 |

## Run T1: single-session tutorial (chapter 10)

Dataset `test_data/manual_ern_singlesession.csv`; the two incongruent events
(`inc_err`, `inc_cor`, as section 10.4 selects), both groups; defaults (4
chains, 5,000 warmup, 5,000 sampling, seed 12345). Cutoff blocks read from the viewer
at dependability cutoff .70 and again at .80; the plot block reads the
dependability-versus-trials figure at 50 trials.

| Slug | Value |
|---|---|
| t1_rhat_max | 1.00 |
| t1_ndivergent | 0 |
| t1_overall70_grpA_nincl | 23 |
| t1_overall70_grpA_nexcl | 17 |
| t1_overall70_grpA_incerr_dep | .75 |
| t1_overall70_grpA_incerr_meantrl | 28.6 |
| t1_overall70_grpA_inccor_dep | .80 |
| t1_overall70_grpA_inccor_meantrl | 31.3 |
| t1_overall70_grpB_nincl | 5 |
| t1_overall70_grpB_nexcl | 35 |
| t1_overall70_grpB_incerr_dep | .72 |
| t1_overall70_grpB_incerr_meantrl | 34.4 |
| t1_overall70_grpB_inccor_dep | .78 |
| t1_overall70_grpB_inccor_meantrl | 30.0 |
| t1_overall80_grpA_nincl | 0 |
| t1_overall80_grpB_nincl | 0 |
| t1_cut70_grpA_incerr_dep | .70 [.55, .83] |
| t1_cut70_grpA_incerr_n | 22 |
| t1_cut70_grpA_inccor_dep | .70 [.58, .81] |
| t1_cut70_grpA_inccor_n | 18 |
| t1_cut70_grpB_incerr_dep | .70 [.52, .84] |
| t1_cut70_grpB_incerr_n | 32 |
| t1_cut70_grpB_inccor_dep | .70 [.57, .82] |
| t1_cut70_grpB_inccor_n | 20 |
| t1_cut80_grpA_incerr_n | 38 (extrapolated: -1 sentinel in the coefficient cell) |
| t1_cut80_grpA_inccor_n | 31 |
| t1_cut80_grpB_incerr_n | 55 (extrapolated: -1 sentinel in the coefficient cell) |
| t1_cut80_grpB_inccor_n | 34 |
| t1_plot_grpA_incerr_dep50 | .84 |
| t1_plot_grpA_inccor_dep50 | .87 |
| t1_plot_grpB_incerr_dep50 | .79 |
| t1_plot_grpB_inccor_dep50 | .85 |
| t1_var_grpA_incerr_betsd | 1.85 |
| t1_var_grpA_incerr_witsd | 5.48 |
| t1_var_grpA_incerr_trlsd | 0.43 |
| t1_var_grpA_incerr_icc | .10 |
| t1_var_grpA_incerr_sem | 1.03 |
| t1_var_grpA_inccor_betsd | 1.88 |
| t1_var_grpA_inccor_witsd | 4.98 |
| t1_var_grpA_inccor_trlsd | 1.01 |
| t1_var_grpA_inccor_icc | .12 |
| t1_var_grpA_inccor_sem | 0.91 |
| t1_var_grpB_incerr_betsd | 2.17 |
| t1_var_grpB_incerr_witsd | 7.59 |
| t1_var_grpB_incerr_trlsd | 1.14 |
| t1_var_grpB_incerr_icc | .08 |
| t1_var_grpB_incerr_sem | 1.31 |
| t1_var_grpB_inccor_betsd | 1.84 |
| t1_var_grpB_inccor_witsd | 5.09 |
| t1_var_grpB_inccor_trlsd | 1.16 |
| t1_var_grpB_inccor_icc | .11 |
| t1_var_grpB_inccor_sem | 0.95 |
| t1_rec_grpA_incerr_betsd_pt | 1.85 |
| t1_rec_grpA_incerr_betsd_ci | [1.30, 2.54] |
| t1_rec_grpA_incerr_witsd_pt | 5.48 |
| t1_rec_grpA_incerr_witsd_ci | [5.23, 5.75] |
| t1_rec_grpA_incerr_trlsd_pt | 0.43 |
| t1_rec_grpA_incerr_trlsd_ci | [0.02, 1.04] |
| t1_rec_grpA_inccor_betsd_pt | 1.88 |
| t1_rec_grpA_inccor_betsd_ci | [1.41, 2.47] |
| t1_rec_grpA_inccor_witsd_pt | 4.98 |
| t1_rec_grpA_inccor_witsd_ci | [4.79, 5.19] |
| t1_rec_grpB_incerr_betsd_pt | 2.17 |
| t1_rec_grpB_incerr_betsd_ci | [1.43, 3.06] |
| t1_rec_grpB_incerr_witsd_pt | 7.59 |
| t1_rec_grpB_incerr_witsd_ci | [7.24, 7.95] |
| t1_rec_grpB_inccor_betsd_pt | 1.84 |
| t1_rec_grpB_inccor_betsd_ci | [1.37, 2.45] |
| t1_rec_grpB_inccor_witsd_pt | 5.09 |
| t1_rec_grpB_inccor_witsd_ci | [4.89, 5.30] |

## Run T2: test-retest tutorial (chapter 11)

Dataset `test_data/manual_ern_testretest.csv`; both events, both occasions;
4 chains, 10,000 warmup, 10,000 sampling, seed 12345 (the chapter raises the
iteration counts and says so). Cutoff blocks read from the viewer at the
criterion settings the chapter names (CE70, CES50, CS50).

| Slug | Value |
|---|---|
| t2_rhat_max | 1.00 |
| t2_ndivergent | 909 |
| t2_overallCE70_nincl | 30 |
| t2_overallCE70_nexcl | 0 |
| t2_overallCE70_err_rel | .91 |
| t2_overallCE70_err_meantrl | 24 |
| t2_overallCE70_cor_rel | .93 |
| t2_overallCE70_cor_meantrl | 40 |
| t2_overallCES50_err_rel | .55 |
| t2_overallCES50_cor_rel | .51 |
| t2_overallCS50_err_rel | .56 |
| t2_overallCS50_cor_rel | .51 |
| t2_genCE70_err_rel | .91 |
| t2_genCE70_cor_rel | .93 |
| t2_cutCE70_err_rel | .72 [.62, .82] |
| t2_cutCE70_err_n | 6 |
| t2_cutCE70_cor_rel | .73 [.63, .82] |
| t2_cutCE70_cor_n | 8 |
| t2_cutCES50_err_n | 11 |
| t2_cutCES50_cor_n | 32 |
| t2_cutCS50_err_n | 9 |
| t2_cutCS50_cor_n | 29 |
| t2_var70_err_betsd | 3.47 |
| t2_var70_err_witsd | 5.42 |
| t2_var70_err_icc | .30 |
| t2_var70_err_sem | 1.20 |
| t2_var70_cor_betsd | 2.61 |
| t2_var70_cor_witsd | 5.05 |
| t2_var70_cor_icc | .25 |
| t2_var70_cor_sem | 0.83 |
| t2_comp_err_bp | 3.47 |
| t2_comp_err_bt | 1.25 |
| t2_comp_err_bo | 2.58 |
| t2_comp_err_txp | 1.73 |
| t2_comp_err_oxp | 1.70 |
| t2_comp_err_txo | 0.28 |
| t2_comp_err_err | 5.42 |
| t2_comp_cor_bp | 2.61 |
| t2_comp_cor_bt | 1.12 |
| t2_comp_cor_bo | 2.16 |
| t2_comp_cor_txp | 0.75 |
| t2_comp_cor_oxp | 1.56 |
| t2_comp_cor_txo | 0.22 |
| t2_comp_cor_err | 5.05 |

## Run T3: difference-score tutorial (chapter 12)

Dataset `test_data/manual_ern_singlesession.csv`; incongruent events, both
groups; two-event difference scores on (subject-level error variances and
residual covariance both left at their default, off, as the chapter states), all
other preferences at defaults (4 chains, 5,000 warmup, 5,000 sampling, seed
12345).

| Slug | Value |
|---|---|
| t3_rhat_max | 1.00 |
| t3_ndivergent | 0 |
| t3_overall70_grpA_incerr_dep | .75 |
| t3_overall70_grpA_incerr_meantrl | 28.6 |
| t3_overall70_grpA_inccor_dep | .80 |
| t3_overall70_grpA_inccor_meantrl | 31.3 |
| t3_overall70_grpA_diff_rel | .40 |
| t3_overall70_grpB_incerr_dep | .72 |
| t3_overall70_grpB_incerr_meantrl | 34.4 |
| t3_overall70_grpB_inccor_dep | .78 |
| t3_overall70_grpB_inccor_meantrl | 30.0 |
| t3_overall70_grpB_diff_rel | .46 |
| t3_gen_grpA_diff_rel | .41 |
| t3_gen_grpB_diff_rel | .47 |
| t3_cut70_grpA_incerr_dep | .70 [.55, .82] |
| t3_cut70_grpA_incerr_n | 22 |
| t3_cut70_grpA_inccor_dep | .71 [.59, .81] |
| t3_cut70_grpA_inccor_n | 19 |
| t3_cut70_grpA_diff_rel | .70 [.39, .87] |
| t3_cut70_grpA_diff_n | 114 |
| t3_cut70_grpB_incerr_dep | .71 [.53, .83] |
| t3_cut70_grpB_incerr_n | 32 |
| t3_cut70_grpB_inccor_dep | .71 [.58, .82] |
| t3_cut70_grpB_inccor_n | 21 |
| t3_cut70_grpB_diff_rel | .70 [.40, .87] |
| t3_cut70_grpB_diff_n | 97 |
| t3_alldiff_grpA_bpcov | 1.66 |
| t3_alldiff_grpA_errcov | 0.00 |
| t3_alldiff_grpA_icc1 | .02 |
| t3_alldiff_grpB_bpcov | 1.68 |
| t3_alldiff_grpB_icc1 | .03 |

## Run T4: splits tutorial (chapter 13)

Dataset `test_data/manual_ern_splits.csv`; weight column `weight`, split type
Nonparallel (unequal n_i); defaults (4 chains, 5,000 warmup, 5,000 sampling,
seed 12345). The two Splits n' values are the Dependability column of the
coefficients table after typing 3 and then 12 into Splits n' and clicking
Show coefficients (section 13.6), at the viewer's default 0.95 credible
interval.

| Slug | Value |
|---|---|
| t4_rhat_max | 1.00 |
| t4_ndivergent | 0 |
| t4_nsplits | 16 |
| t4_nocc | 1 |
| t4_dep_pt | .97 |
| t4_dep_low | .96 |
| t4_dep_high | .98 |
| t4_gen_pt | .97 |
| t4_gen_low | .96 |
| t4_gen_high | .99 |
| t4_sem_abs | 0.39 |
| t4_sem_rel | 0.37 |
| t4_vc_sigma_p_est | 5.52 |
| t4_vc_sigma_p_ci | [3.37, 9.03] |
| t4_vc_sigma_s_est | 0.24 |
| t4_vc_sigma_s_ci | [0.08, 0.56] |
| t4_vc_sigma_ps_est | 0.71 |
| t4_vc_sigma_ps_ci | [0.28, 1.13] |
| t4_vc_sigma_item_est | 21.49 |
| t4_vc_sigma_item_ci | [13.78, 30.85] |
| t4_dep_ns3 | .87 |
| t4_dep_ns12 | .96 |

## Scripting parity (chapter 14)

One sentence recording the outcome of the parity check: a GUI run and a
`psyrat_run` call at identical seed and settings produce identical draws.
Verified once before any transcription above.

| Slug | Value |
|---|---|
| script_parity_note | Verified 2026-09-02 on the capture machine, on the pre-regeneration realization of manual_ern_singlesession.csv with all four events selected: a GUI run (Process New Data, Analyze) and a psyrat_run call at identical settings (both groups, 4 chains, 5,000 warmup, 5,000 sampling, seed 12345) produced bit-identical posterior draws in all 32 group-by-event-by-component cells (maximum absolute difference 0). Both surfaces call the same estimation function, so the result is a property of the code path rather than of the data; the regenerated file and the two-event selection were not re-checked. |
