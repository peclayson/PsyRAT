# Contributing to PsyRAT

Thank you for your interest in improving PsyRAT.

## Reporting bugs

Please open an issue on the
[issue tracker](https://github.com/peclayson/PsyRAT/issues) and include:

- the **PsyRAT version string** (printed in the startup banner, returned by
  `psyrat_defineversion`, and written into every exported table header);
- the **`*_runconfig.json` sidecar** saved with the analysis, when one exists —
  it records the exact configuration and makes the problem reproducible;
- what you expected, what happened instead, and the full text of any error.

## How this repository is maintained

Development happens in a separate working repository; this public repository
receives complete, tested release snapshots. Two practical consequences:

- **Pull requests are welcome as proposals**, but they are not merged directly
  into this repository's history. Accepted changes are applied in the
  development repository, validated there, and arrive here with the next
  release. Your contribution is credited in the changelog.
- **There is no continuous integration.** All verification is local and manual:
  the four MATLAB suites (`runAllUnitTests`, `runAccuracyValidationSuite`,
  `runFeatureValidationSuite`, `runExhaustiveWorkflowBugSweep`) are the gate,
  and changes touching Bayesian estimation additionally require live CmdStan
  recovery checks on the maintainer's machine. Please run at least the unit and
  accuracy suites against your proposed change and say so in the PR.

## Ground rules for changes

- **Scientific accuracy first.** Changes to reliability formulas or estimation
  code must cite the source (paper, table, and equation) they implement.
- **Behavior preservation.** Estimation is seeded and bit-reproducible
  (default seed 12345); a change that alters outputs needs an explicit
  justification, not just green tests.
- **Minimal diffs.** Especially in the two large core files
  (`subroutines/estimation/psyrat_computevarcomp.m` and
  `subroutines/calculation/psyrat_relsummary.m`), keep changes as small as
  reviewable correctness allows and do not refactor for style.

## Questions

For usage questions, please also use the issue tracker; questions that recur
become documentation.
