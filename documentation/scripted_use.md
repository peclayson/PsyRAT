# Scripted (command-line) use of PsyRAT

PsyRAT is GUI-first, but every analysis the Process New Data workflow runs can
also be run from a script with `psyrat_run`. The scripted path is not a second
implementation: it builds the same `psyrat_prefs`/`psyrat_data` scaffold the GUI
builds and calls the same shared estimation code. Given the same inputs, the same
settings, and the same seed, a scripted run produces the same variance
components, the same convergence diagnostics, and the same saved `.psyrat` file
as the GUI. Two settings differ by default: `psyrat_run` runs headless
(`'interactive', false`), so on non-convergence it doubles the iterations and
refits automatically where the GUI asks first, and it never renders trace plots
(Chapter 14 of the manual covers both).

**Chapter 14 of the user manual is the scripting reference.** This page is the
quickstart; [Chapter 14](manual/14_scripting.md) is the full treatment, and
where the two ever disagree, the manual is correct.

## Requirements

Identical to the GUI: a MATLAB release the toolbox runs on (tested on R2025b;
no minimum release is enforced in code, and R2022a is the effective floor), a
working CmdStan for the full Bayesian
path, and the patched MatlabStan bundled in `bundled_dependents/`. See
[Chapter 3](manual/03_installation.md) for installation.

Get the toolbox on the path by running `psyrat_start` once in the session. That
adds the toolbox root, `subroutines/`, and `bundled_dependents/`. Do not
`addpath(genpath(...))` the whole repository: stale copies in tooling or
artifact folders can shadow the real files.

## A first run

This runs the single-session ERN tutorial dataset, the same one
[Chapter 10](manual/10_tutorial-single-session.md) walks through the GUI.

```matlab
% Locate the toolbox root from a file inside it rather than hardcoding a path.
psyratroot = fileparts(which('psyrat_run'));
datafile   = fullfile(psyratroot,'test_data','manual_ern_singlesession.csv');

% The output directory must already exist, and CmdStan rejects whitespace
% anywhere in the output path or filename.
outdir = fullfile(tempdir,'psyrat_demo');
if ~exist(outdir,'dir'); mkdir(outdir); end

psyrat_data = psyrat_run( ...
    'file',     datafile, ...
    'idcol',    'id', ...        % the default; shown so the mechanism is visible
    'meascol',  'meas', ...      % also the default; most real files need both
    'eventcol', 'event', ...
    'groupcol', 'group', ...           % chapter 10 keeps the two groups separate
    'whichevents', {'inc_err','inc_cor'}, ... % and processes these two events
    'savepath', outdir, ...
    'savename', 'ern_demo.psyrat', ...
    'seed',     12345);          % fix the seed, and report it
```

Estimation stops at the variance components. Two further headless calls turn
those into the reliability tables the GUI viewer shows: `psyrat_relsummary`
runs the D-study, then `psyrat_report` assembles and optionally exports the
tables. Both steps, the `analysis` dispatch key each design requires, and the
families that skip `psyrat_relsummary` entirely are covered in
[Chapter 14](manual/14_scripting.md).

## Where to read next

- [Chapter 14, Scripting with `psyrat_run`](manual/14_scripting.md) covers the
  two-step reporting chain, the `analysis` dispatch table, headless behavior and
  batch loops, replaying a run from its `*_runconfig.json` sidecar, engine
  choice, and a full option reference for both `psyrat_run` and `psyrat_report`.
- [Chapter 15, Advanced designs](manual/15_advanced-designs.md) covers the gamma
  family, dynamic reliability, and setting priors.
- [Priors and sensitivity](priors_and_sensitivity.md) documents the default
  priors, which are not all on one scale. The Gaussian-family defaults are
  weakly informative on the ERP microvolt scale and are not scale-free, so data
  on another scale needs values matched to it. The gamma-family defaults instead
  live on log scales and are calibrated to a time-frequency-power reference.
- `help psyrat_run` is the authoritative option list.
