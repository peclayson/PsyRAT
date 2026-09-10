# Estimation baselines (W0)

Frozen seeded posterior-draw references used to prove that later changes
(especially the W3 prior edits) alter estimation output only where intended.

The `.mat` files are **not committed**. They are bit-reproducible only on the
machine and CmdStan version that generated them, so they serve as local
regression references rather than portable golden values. Regenerate them with
`PsyRATBaselineCapture`.

## Capture configuration (fixed)

4 chains, 1000 warmup, 1000 sampling, seed 12345. Reproducibility is exact at
any length; the W3 prior-change check re-runs at this same config and diffs, so
production length is not needed for the intentional-change comparison.

## Regenerate

```matlab
root = pwd;  % PsyRAT project root
addpath(genpath(fullfile(root,'subroutines')));
addpath(genpath(fullfile(root,'bundled_dependents')));
addpath(genpath(fullfile(root,'tests','helpers')));

scenarios = {'ic_base','trt_time','ic_diff_event','ic_sserr_base'};
for i = 1:numel(scenarios)
    PsyRATBaselineCapture.capture(root, scenarios{i}, ...
        fullfile(root,'tests','baselines',[scenarios{i} '.mat']));
end
```

## Verify reproducibility (sanity check)

```matlab
rpt = PsyRATBaselineCapture.compare(root, 'trt_time', ...
    fullfile(root,'tests','baselines','trt_time.mat'));
% rpt.bitIdentical should be true and rpt.maxAbsDiff 0 on the capture machine.
```

As of the W0 capture, all four scenarios verified `bitIdentical = true`,
`maxAbsDiff = 0`.

## Scenarios captured

| Baseline | Design family | Notes |
| --- | --- | --- |
| `ic_base` | single-occasion (basic) | person + residual model |
| `trt_time` | test-retest | 7-component crossed model; the SCI-004 occasion priors live here |
| `ic_diff_event` | difference score | bivariate location-scale model |
| `ic_sserr_base` | subject-level | per-subject residual (location-scale) |

Not captured separately: difference-of-differences (no 4-event scenario in the
integration data factory yet) and cut-score (derived from the single-occasion /
test-retest variance components, not a separate estimation run).
