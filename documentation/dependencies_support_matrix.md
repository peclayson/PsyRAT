# PsyRAT dependencies and support matrix

Developer and installation reference for PsyRAT's external dependencies: the
versions it bundles or targets, the version constraints that actually matter,
where it looks for CmdStan, and the PsyRAT-specific patches carried in the
bundled MatlabStan fork. Companion to the user manual (planned).

Every claim below cites the file it was verified against. Items that are not
yet pinned are marked "Open."

## Bundled dependency versions

Single source of truth: `subroutines/installation/psyrat_dependentsversions.m`.

| Dependency | Bundled / targeted version | Role |
|---|---|---|
| CmdStan | 2.38.0 | Compiles and samples the generated Stan models |
| MatlabStan (PsyRAT fork) | 2.15.1.0 | MATLAB interface to CmdStan, bundled at `bundled_dependents/MatlabStan-2.15.1.0` |
| MatlabProcessManager | 0.5.1 | Spawns and monitors the CmdStan subprocess, bundled at `bundled_dependents/MatlabProcessManager-0.5.1` |

MatlabStan and MatlabProcessManager are vendored in the repository. CmdStan is
installed separately (see "CmdStan discovery" below).

## CmdStan version support

- Minimum: 2.26. The generated Stan models use the `array[N] T name;`
  declaration syntax introduced in Stan 2.26. It appears in every generated
  model (for example `array[NOBS] int<lower=0, upper=NSUB> id;`, captured in
  `tests/golden/stan/ic_base.stan`). CmdStan older than 2.26 cannot compile
  these models. This minimum is a property of the generated code, not an
  enforced runtime check.
- Targeted and tested: 2.38.0, the version declared in
  `psyrat_dependentsversions.m` and used for local CmdStan integration and
  accuracy runs.
- Version detection: `subroutines/installation/psyrat_checkversionsofdeps.m`
  reads the installed CmdStan version two ways depending on the release. Older
  CmdStan ships `cmdstan-guide.tex`; the document layout changed around CmdStan
  2.22, so newer releases are queried via `bin/stanc --version`. This logic
  supports update reporting, not minimum enforcement.
- Preflight (enforced): `psyrat_require_min_cmdstan`
  (`subroutines/installation/`) detects the installed CmdStan version and
  reports a clear message when it is below 2.26 instead of letting CmdStan emit
  a cryptic compile error. It runs at startup in `psyrat_start` (warn mode, so
  the GUI still opens) and again before the first compile in `psyrat_run_stan`
  (error mode). The check is fail-open: if the version cannot be determined it
  does not block.

## Estimation engines (CmdStan and native fallbacks)

PsyRAT can estimate variance components with three engines. CmdStan is the
default and the only one required for full functionality; the native engines are
optional fallbacks for when CmdStan cannot be installed.

| Engine | Code | Requires | Notes |
|---|---|---|---|
| CmdStan | 1 | CmdStan >= 2.26 (+ bundled MatlabStan). Statistics and Machine Learning Toolbox additionally required for the gamma concurrent-difference designs (analyses 7, 8, 10 with `diffrescor=2`) | Default. Full Bayesian model; the only engine that covers all designs. See the toolbox exception in the notes below this table. |
| Native HMC | 3 | Statistics and Machine Learning Toolbox (`hmcSampler`); Deep Learning Toolbox (`dlarray`) for automatic-differentiation gradients | Faithful Bayesian fallback for the location-scale / multivariate designs. Wired for 18 design cases (5-22: test-retest, subject-level/location-scale, dynamic/conditional reliability, difference scores, and the difference-of-differences family). |
| Native fitlme | 2 | Statistics and Machine Learning Toolbox (`fitlme`) | Frequentist (REML) fallback for the simple plain-Gaussian designs. Intervals are APPROXIMATE (parametric bootstrap), NOT Bayesian credible intervals. Wired for 7 design cases (1-4: one-facet and group/event internal consistency; 5: test-retest; 23-24: nonparallel splits). |

- Engine preference: `psyrat_prefs.proc.engine` (`0` auto, `1` cmdstan,
  `2` fitlme, `3` hmc; default `0`). Exposed on the CLI as `psyrat_run(...,
  'engine', ...)` accepting the codes or the names `auto`/`cmdstan`/`fitlme`/`hmc`.
- Auto cascade (engine `0`): tries CmdStan, then native HMC, then native fitlme,
  taking the first engine that is both installed and implements the requested
  design. When CmdStan is present the run is unchanged; whenever the cascade
  falls back off CmdStan it emits a clear, non-suppressible warning (fitlme
  fallbacks additionally warn that the intervals are approximate).
- Capability detection (single source of truth):
  `subroutines/installation/psyrat_estimation_engines_available.m` reports which
  engines are runnable (CmdStan present; Statistics and Machine Learning Toolbox;
  Deep Learning Toolbox). The GUI dropdown, the cascade, and forced-engine
  validation all consult it, so an engine is never offered unless its
  prerequisites are installed.
- Optional dependency, with one documented exception: the native-engine toolboxes
  are otherwise needed only to use the corresponding native engine.
- **Exception (Statistics and Machine Learning Toolbox, gamma concurrent
  differences).** Three designs require that toolbox even on the default CmdStan
  path: the gamma (scaled chi-square) **concurrent** difference scores, i.e.
  `family='gamma'` together with `diffrescor=2`, for analyses **7** (two-event
  difference), **8** (its subject-level variant) and **10** (two-facet test-retest
  difference). Their Gaussian-copula extractors call `normcdf` and `gaminv`, which
  are toolbox functions with no base-MATLAB equivalent. Both dispersion
  parameterizations are affected (per-person log-nu and global fixed nu).
  Where: the five `psyrat_gamma_extract_diff_*copula*` local functions in
  `subroutines/estimation/psyrat_computevarcomp.m`. Non-concurrent gamma
  differences (`diffrescor=1`) do not reach them, and no other family does.
- The failure mode is late and expensive: extraction runs **after** sampling has
  completed, so a user without the toolbox loses the finished fit to
  `MATLAB:UndefinedFunction` rather than being stopped up front. **There is no
  CI to catch it: the GitHub Actions validation lane was retired 2026-08-20
  (`2fbd5d7`), and Actions is disabled at the repo level.** While it ran it did
  not catch this either, because the copula test classes skip via `assumeTrue`
  when the toolbox is absent. A developer machine typically has the toolbox
  installed, so a green local run is not evidence either.
- Provenance: every run records the engine preference and the engine actually
  used in the `*_runconfig.json` sidecar (`psyrat_write_runconfig`), and the
  estimated `REL.engine` field.

## MATLAB support

- Tested: MATLAB R2025b (Update 5) for the unit, accuracy, and CmdStan
  integration runs, all on the maintainer's local machine. **There is no CI:
  the GitHub Actions lane was retired 2026-08-20 (`2fbd5d7`), so every release
  claim in this document rests on a local run.**
- No minimum MATLAB release is enforced in code. The toolbox relies on `table`,
  `string`, `categorical`, and the `matlab.unittest` framework, all available
  in current releases.
- Open: the precise minimum supported MATLAB release has not been pinned.

## Operating systems

- Developed and tested on macOS (Darwin). The retired CI lane ran on `ubuntu-latest`.
- CmdStan additionally requires a C++ toolchain (make plus a C++ compiler) to
  build models on any platform. The Windows installer checks for both `make`
  (GNU Make) and `g++` before attempting a CmdStan build.
- Linux auto-install is EXPERIMENTAL and has not been verified by the
  maintainer. `psyrat_installdependents` includes a Linux branch
  (`psyrat_linuxdepsinstall`) that mirrors the macOS path: it checks for
  `make`/`g++`, downloads the CmdStan release tarball, and builds it with
  `make`. If it fails, install CmdStan manually and point the toolbox at it via
  the `STAN_HOME` or `CMDSTAN` environment variable. Estimation itself runs on
  Linux once CmdStan is discoverable (this is how the retired `ubuntu-latest` CI
  lane ran).
- Open: Windows is not part of the current test matrix; it follows the general
  MatlabStan and CmdStan support but is unverified here.

## CmdStan discovery

PsyRAT locates CmdStan through the bundled fork's `+mstan/stan_home.m`
(PsyRAT-patched). Resolution order (the FIRST candidate that is a valid CmdStan
root wins; a valid root is a directory containing a `makefile` and a `bin/`
subdirectory):

1. Environment variables, in order: `PSYRAT_TEST_CMDSTAN_PATH`, `STAN_HOME`,
   `CMDSTAN_HOME`, `CMDSTAN`. Every non-empty variable is considered, in that
   order.
2. `cmdstan-*` globs under home-directory locations, in order:
   `~/Documents/MATLAB/PsyRATDependents`, `~/Documents/PsyRATDependents`,
   `~/MATLAB/PsyRATDependents`, and `~` itself. Within one glob, candidates
   are ordered newest parseable `cmdstan-X.Y.Z` version first (G43, patch
   level 2; previously directory-listing order, under which `cmdstan-2.38.0`
   beat `cmdstan-2.40.0` digit-by-digit), with unparseable names after the
   parseable ones in listing order. A machine that must run an OLDER install
   than its newest one should pin it with an environment variable.
3. MATLAB-path fallback, consulted only when nothing above resolved:
   `which('runCmdStanTests.py')`, then `which('stanc')`, then
   `which('makefile')`, then a scan of every entry on the MATLAB path for
   `runCmdStanTests.py`.

Since G43 (2026-08-29) the toolbox side consults ONE discovery point,
`subroutines/installation/psyrat_locate_cmdstan.m`, which delegates to
`mstan.stan_home` -- so the startup dependency gate (`psyrat_start`), the
version check (`psyrat_checkversionsofdeps`), the version/availability
probe (`psyrat_detect_cmdstan_version`), and the stan_home argument injected
at the `stan()` call (`psyrat_find_cmdstan_home`, which consults the
resolver ahead of its legacy tiers) all see the install a fit would use.
Two behavior changes came with that: an environment-variable-only CmdStan
(not on the MATLAB path) now passes the startup gate instead of hitting the
install dialog, and its version is read correctly instead of tripping the
"Unable to determine CmdStan version" warning. Of the former divergent
lookups, exactly two legacy fallbacks remain, both consulted only when the
resolver finds nothing: the injection helper's own tiers (which skip
`STAN_HOME` entirely), and the version check's `cmdstan-guide.tex` read for
pre-~2.22 installs. The probe's path-then-env lookup is gone outright, and
the repair template inside `psyrat_computevarcomp` now emits the canonical
resolver order rather than its former which()-before-globs order.

The test harness honors the same `PSYRAT_TEST_CMDSTAN_PATH` override, plus
`PSYRAT_TEST_MATLABSTAN_PATH` and `PSYRAT_TEST_MPM_PATH`, resolved in
`tests/helpers/PsyRATCmdStanDependencySetup.m`.

## Bundled MatlabStan fork patches

The bundled MatlabStan differs from upstream `brian-lau/MatlabStan` 2.15.1.0 in
six PsyRAT-specific ways. The first five are covered by `patch_level = 1` in
`bundled_dependents/MatlabStan-2.15.1.0/PSYRAT_PATCHLEVEL.txt` (the stamp
postdates them); the sixth is `patch_level = 2`.

1. Within-chain threading (`threads_per_chain`). Adds a `threads_per_chain`
   option (scalar `>= 1`, or `'auto'`, default `1` meaning off) that drives
   CmdStan's `reduce_sum` within-chain parallelism via the `STAN_NUM_THREADS`
   environment variable. `'auto'` resolves to `floor(available_cores/chains)`
   with a minimum of 1 and warns on oversubscription. Touched files: `stan.m`,
   `StanModel.m` (see `resolve_threads_per_chain`), `README.md`, and
   `Tests/TestStanModel.m`. Upstream 2.15.1.0 has no threading support.
2. CmdStan auto-discovery. `+mstan/stan_home.m` adds `PSYRAT_TEST_CMDSTAN_PATH`
   to the environment-variable precedence list and adds the home-directory
   globs above (the three `PsyRATDependents/cmdstan-*` locations plus
   `~/cmdstan-*`), so a bundled-style CmdStan install is found without manual
   configuration.
3. Permuted-extract indexing repair (S19, 2026-08-10). `mcmc.m:180` replaces an
   upstream out-of-bounds draw index (`1:max(sz)`) with
   `self.permute_index(1:sz(1))`, so the permuted-extract path returns instead
   of erroring. Bit-identical on every previously working extraction; the
   change turns an erroring path into a returning one.
4. `nansum` removal (B39, 2026-08-21). `StanFit.m:272` and `:375` replace
   `nansum` with `sum(...,'omitnan')`, so the default CmdStan path no longer
   requires the Statistics and Machine Learning Toolbox.
5. `prctile` removal in PSIS (B40, 2026-08-21). `psis/psislw.m:38` replaces
   `prctile` with a toolbox-free equivalent, inside the separately GPL-licensed
   `psis/` subpackage.
6. CmdStan home-glob version sort (G43, 2026-08-29; patch level 2).
   `+mstan/stan_home.m` orders each managed `cmdstan-*` location newest
   parseable version first (see "CmdStan discovery" above), so a machine
   holding several managed installs resolves the newest one instead of the
   lexicographic first match.

On the PsyRAT side (not in the fork), `psyrat_run_stan` in
`subroutines/estimation/psyrat_computevarcomp.m` adds robustness around
MatlabStan: it repairs a malformed `+mstan/stan_home.m` and retries once on a
`StanModel:stan_home:InputFormat` error before failing.
