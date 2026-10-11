# Validation and testing

This page describes the four MATLAB test suites that gate a PsyRAT release and,
for the accuracy suite, what a green run of each of its classes does and does
not prove. The manual's tutorials and the design-maturity table in `CHANGELOG.md`
refer here when they say that PsyRAT's correctness claims rest on its validation
suites.

## Running the suites

There is no continuous integration. The four MATLAB suites below are run by the maintainer before every release, and you can run them yourself from the toolbox root:

```matlab
[allPassed, results] = runAllUnitTests('Verbose', true);            % core unit tests (non-integration)
[allPassed, results] = runAccuracyValidationSuite('Verbose', true); % deterministic calculation accuracy checks
[allPassed, results] = runFeatureValidationSuite('Verbose', true);  % integration/feature validation (needs CmdStan; slow)
[allPassed, report]  = runExhaustiveWorkflowBugSweep('Verbose', true); % exhaustive workflow bug sweep
```

Optional JUnit output is available from the unit, accuracy, and feature runners using `'ReportJUnit', true` and `'JUnitFile', '<path>.xml'`. The workflow sweep has no JUnit option; it writes its own timestamped `workflow_bug_sweep_*.csv` and `.mat` report files to the toolbox root (`'SaveReport'` and `'ReportPrefix'` control them; `help runExhaustiveWorkflowBugSweep` lists the options).

## What a green accuracy run proves

The accuracy suite spans three test classes with **different evidentiary weight**, and the distinction matters when reporting results:

| class | what a green run proves |
|---|---|
| `TestCalculationAccuracyOracle` | **Regression only.** `PsyRATAccuracyOracle` duplicates the production formulas, so this shows they are *unchanged*, not that they are correct. |
| `TestIndependentOracleAgreement` | **External corroboration.** `PsyRATIndependentOracle` states the published formulas of Rocha et al. (2026) in the paper's own symbols. Its expressions have been checked term by term against Tables 2, 3 and 6 of the published article (doi:10.1016/j.ijpsycho.2026.113321), and they agree. That verified agreement with the printed tables is what makes agreement with production evidence the formulas are *right*, not merely unchanged. The class header records that its Table 6 block has weaker provenance than Tables 2 and 3. |
| `TestExternalBenchmark` | **External numeric anchor.** Validates against trial-count tables from a separate analysis; covers one-facet, one-facet difference, and difference-of-differences only. |

Do not cite a green run of this suite as evidence of accuracy for a design family none of
the three covers. `TestIndependentOracleAgreement` restates Tables 2 and 3 for the one-facet,
test-retest, and subject-level test-retest coefficients and Table 6 for the difference-score
coefficients; the splits designs and the dynamic/conditional family have no independent-oracle
row, and the difference-of-differences designs rest on the external benchmark and the in-repo
recovery tests alone.

Live CmdStan recovery against simulated data with known variance components is run by the
maintainer where a recovery test exists (`tests/integration/`), and the check that re-derives every
reference value the manual's tutorials print is run before a release; the design-maturity table
in `CHANGELOG.md` states, design by design, which of these checks each design has and has not had.
