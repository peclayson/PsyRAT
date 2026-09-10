# PsyRAT Toolbox

The Psychophysiologist's Reliability Analysis Toolbox (PsyRAT) uses generalizability theory to evaluate the psychometric reliability of psychophysiological measurements, such as event-related potential (ERP) component scores, oscillatory EEG, skin conductance, and cardiac data. Variance components are estimated in a Bayesian framework via [CmdStan](https://mc-stan.org/users/interfaces/cmdstan), and the toolbox reports dependability (absolute-error) and generalizability (relative-error) coefficients for designs with any number of events, groups, or occasions. It also estimates how dependability changes with the number of trials and recommends a trial-count cutoff for including a participant's data, based on the stability of measurement as trials accumulate in a participant's average for a given group and event. Beyond single-session internal consistency, PsyRAT estimates test-retest reliability, the reliability of difference scores, subject-level (per-participant) reliability, reliability that varies along a participant-level dimension (dynamic or conditional reliability), and reliability from data splits. PsyRAT succeeds the [ERP Reliability Analysis (ERA) Toolbox](https://github.com/peclayson/ERA_Toolbox); the estimators it implements are described in [Rocha et al. (2026)](https://doi.org/10.1016/j.ijpsycho.2026.113321).

See the [user manual](documentation/manual/index.md) for full documentation of the toolbox. It can be read as chapters on GitHub, starting from the contents page, or as [one PDF](documentation/manual/psyrat_manual.pdf).

The current release is `0.1.0-beta`. The beta label means that the interface and the shape of the outputs are still moving; it is not a statement that the estimates are provisional. [`CHANGELOG.md`](CHANGELOG.md) records, design by design, what has been verified in this repository and what has not.

## Getting started

1. **Get the toolbox.** Download `PsyRAT.zip` from the [Releases page](https://github.com/peclayson/PsyRAT/releases) and unzip it, or clone this repository, into a folder whose full path contains no spaces (CmdStan does not accept whitespace in paths). In MATLAB, change into the PsyRAT folder and run `psyrat_start`. The launcher adds the toolbox to your MATLAB path and saves the path, so later sessions find it from anywhere. Do not `addpath(genpath(...))` the whole repository; the launcher adds exactly the folders it needs.
2. **Install CmdStan.** If the Stan dependencies are missing, `psyrat_start` stops before the home screen and offers the guided installer, which you can also run directly:

   ```matlab
   psyrat_installdependents
   ```

   It locates or downloads CmdStan and puts the bundled MatlabStan and MatlabProcessManager on the path.
3. **Run your first analysis.** [Chapter 2 of the manual](documentation/manual/02_quick-start.md) runs one analysis end to end on a simulated dataset that ships in `test_data/`.

### Requirements

- **MATLAB.** Tested on R2025b. No minimum release is enforced in code, but the plotting and export code calls functions introduced in R2022a (`clim`) and R2020a (`exportgraphics`, `writetable` with `WriteMode`) with no fallback, so R2022a is the effective floor.
- **CmdStan** 2.26 or newer (2.38.0 is the tested release), installed separately, plus the C++ toolchain CmdStan needs to compile each model (the Xcode Command Line Tools on macOS; GNU make and g++ on Windows and Linux). The GUI does not start without CmdStan. Without it, the native MATLAB engines (REML via `fitlme`, and a lightweight HMC sampler) for the simpler designs are reachable through `psyrat_run` only; with CmdStan installed they also appear in the GUI's Estimation engine popup.
- **Bundled dependencies**, shipped in `bundled_dependents/` with nothing to install: MatlabStan 2.15.1.0 and MatlabProcessManager 0.5.1. The bundled MatlabStan is a patched fork of a package its developer no longer maintains, and the patches are load-bearing; use it rather than the upstream release.
- **MATLAB toolboxes.** The Statistics and Machine Learning Toolbox is required by the native engines (`fitlme`, `hmcSampler`) and, on the default CmdStan path, by the gamma-family concurrent difference-score designs (analyses 7, 8, and 10), where the failure without it occurs only after sampling completes. The Deep Learning Toolbox is optional: native HMC uses it for automatic-differentiation gradients, and a numeric fallback runs without it.
- **Platforms.** The toolbox is developed and verified on macOS. Windows is not in the test matrix (no Windows run has been verified), and Linux installation is experimental.

[Chapter 3 of the manual](documentation/manual/03_installation.md) covers installation, and [documentation/dependencies_support_matrix.md](documentation/dependencies_support_matrix.md) gives the full, version-pinned dependency matrix and platform notes.

## Documentation

- The [user manual](documentation/manual/index.md): 18 chapters covering installation, generalizability theory for ERP scores, every screen of the interface, four worked tutorials on the simulated datasets in `test_data/`, an analysis reference, reporting guidance, and troubleshooting. Also available as [a single PDF](documentation/manual/psyrat_manual.pdf).
- [documentation/manual/tutorial_expected_values.md](documentation/manual/tutorial_expected_values.md): where every reference value printed in the tutorials came from. The tutorials print values produced by real runs at the settings each chapter states.
- [documentation/scripted_use.md](documentation/scripted_use.md): the quickstart for running analyses from scripts.
- [documentation/model_availability_matrix.md](documentation/model_availability_matrix.md): which designs can be fitted under which likelihood family and estimation engine.
- [documentation/priors_and_sensitivity.md](documentation/priors_and_sensitivity.md): the variance-component priors, with a sensitivity analysis of the test-retest occasion priors.
- [documentation/native_engine_validation.md](documentation/native_engine_validation.md): the validation record for the native MATLAB engines against CmdStan.

## Scripted use

Every analysis available in the GUI can also be run from a script with `psyrat_run`, which uses the same estimation pipeline and reproducibility contract as the GUI: the same inputs and seed produce the same variance components when within-chain parallelization is off, its default (with it on, the threaded models' draws can differ in the trailing digits between runs; see [Chapter 18](documentation/manual/18_troubleshooting.md)). [Chapter 14 of the manual](documentation/manual/14_scripting.md) is the scripting reference, including an option table for `psyrat_run` and `psyrat_report`; `help psyrat_run` is the authoritative option list.


## Bug reports, questions, and contributing

Please report bugs and ask questions on the [issue tracker](https://github.com/peclayson/PsyRAT/issues). Including the PsyRAT version string (printed in the startup banner and written into every exported table header) and, when possible, the `*_runconfig.json` sidecar saved with the analysis makes problems much faster to reproduce. [`CONTRIBUTING.md`](CONTRIBUTING.md) explains how this repository is maintained and what a proposed change needs. Be clear about what problem occurred and what you expected to happen.

## How to cite

For now, cite Rocha et al. (2026), the source of the estimators the toolbox implements (its Tables 2, 3 and 6 give the one-facet, test-retest, and difference-score expressions), together with the GitHub release of the toolbox version you ran. Report the version returned by `psyrat_defineversion`, which is printed in the startup banner and written into every exported table header and `*_runconfig.json` sidecar; with the source revision recorded in the sidecar, it identifies the code that produced your estimates. The current version is `0.1.0-beta`. Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).

&nbsp;

Rocha, H. A., Holbrook, A., Hajcak, G., Keil, A., Rast, P., Thayer, J. F., Verona, E., Vispoel, W. P., & Clayson, P. E. (2026). [Beyond classical metrics: Generalizability theory across psychophysiological modalities](https://doi.org/10.1016/j.ijpsycho.2026.113321). _International Journal of Psychophysiology_, _222_, Article 113321. doi: 10.1016/j.ijpsycho.2026.113321

&nbsp;

The Bayesian location-scale estimation the toolbox is built on is described in

Rast, P., & Clayson, P. E. (in press). [Enhancing generalizability theory with mixed-effects models for heteroscedasticity in psychological measurement: A theoretical introduction with an application from EEG data](https://doi.org/10.1111/bmsp.70026). _British Journal of Mathematical and Statistical Psychology_. doi: 10.1111/bmsp.70026

&nbsp;

The papers below describe the ERA Toolbox, PsyRAT's predecessor, and the framework PsyRAT carries forward: generalizability theory for ERP scores, test-retest reliability, subject-level reliability, and the reliability of difference scores. Clayson and Miller (2017) is the citation for ERA, not for PsyRAT.

&nbsp;

Clayson, P. E., & Miller, G. A. (2017). [ERP Reliability Analysis (ERA) Toolbox: An open-source toolbox for analyzing the reliability of event-related brain potentials](https://doi.org/10.1016/j.ijpsycho.2016.10.012). _International Journal of Psychophysiology_, _111_, 68-79. doi: 10.1016/j.ijpsycho.2016.10.012

&nbsp;

Clayson, P. E., Carbine, K. A., Baldwin, S. A., Olsen, J. A., & Larson, M. J. (2021). [Using generalizability theory and the ERP Reliability Analysis (ERA) Toolbox for assessing test-retest reliability of ERP scores part 1: Algorithms, framework, and implementation](https://doi.org/10.1016/j.ijpsycho.2021.01.006). _International Journal of Psychophysiology_, _166_, 174-187. doi: 10.1016/j.ijpsycho.2021.01.006 ([preprint](https://psyarxiv.com/kcven/))

&nbsp;

Clayson, P. E., Brush, C. J., & Hajcak, G. (2021). [Data quality and reliability metrics for event-related potentials (ERPs): The utility of subject-level reliability](https://doi.org/10.1016/j.ijpsycho.2021.04.004). _International Journal of Psychophysiology_, _165_, 121-136. doi: 10.1016/j.ijpsycho.2021.04.004 ([preprint](https://psyarxiv.com/ja6bw/))

&nbsp;

Clayson, P. E., Baldwin, S. A., & Larson, M. J. (2021). [Evaluating the internal consistency of subtraction-based and residualized difference scores: Considerations for psychometric reliability analyses of event-related potentials](https://doi.org/10.1111/psyp.13762). _Psychophysiology_, _58_(4), Article e13762. doi: 10.1111/psyp.13762 ([preprint](https://psyarxiv.com/nqwz6))

&nbsp;

***

Copyright (C) 2026 Peter E. Clayson
 
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program ([gpl.txt](gpl.txt)). If not, see 
  <http://www.gnu.org/licenses/>.
