# PsyRAT model availability matrix

What PsyRAT can fit, across the three things a user chooses independently: the **design**
(facets, score type, error level, moderation, residual coupling), the **observation
likelihood family** (Gaussian, or scaled chi-square in one of two parameterizations), and
the **estimation engine** (CmdStan, native HMC, native fitlme).

This file records **capability**, not status. The status and sequence authority, the
findings ledger, and the derivation roadmap for the unbuilt gamma designs are all
maintained in the development repository and are not distributed with the release.
Nothing here restates what those own.

> **Line numbers in this file drift.** The estimator grows with every design increment.
> **Grep for the quoted construct, not the line.** The four constructs that decide
> everything below are greppable in
> `subroutines/estimation/psyrat_computevarcomp.m`:
>
> - `if strcmp(family,'gamma') && ~any(analysis ==` sets the family whitelist
> - `if gammascale == 2 && ~any(analysis ==` sets the location-scale accept list
> - `switch analysis` dispatches the design
> - `if strcmp(family,'gamma')` immediately followed by `engine == 2 || engine == 3` is
>   the CmdStan pin
>
> **Read the guards, not this document, for current scope.** Where this file and a guard
> disagree, the guard is right and this file is stale.

---

## 1. The design axes

A run's design is fully determined by six inputs, resolved into an integer `analysis`
code by the selector chain in `psyrat_computevarcomp.m` (grep `analysis = 2` to find it;
it runs from the `splits == 3` test to the `sserrvar == 1 && (ntime > 0)` catch-all).

| Axis | Input | Values |
|---|---|---|
| Facets | `ntime` (derived from a time/occasion column) | one-facet (trials) / two-facet (trials x occasions) |
| Score type | `diffest` | 1 single score / 2 two-event difference / 3 four-event difference-of-differences |
| Error level | `sserrvar` | 1 group-level / 2 subject-level (per-participant residual) |
| Moderation | `dynrel` | 1 static / 2 dynamic, i.e. dimension-conditional (Rast & Clayson) |
| Residual coupling | `diffrescor` (+ `ssrescor`) | non-concurrent / concurrent with a population rho / concurrent with a per-subject rho |
| Splits | `splits` | 1 single-trial / 2 parallel / 3 nonparallel |

`ngroup` and `nevent` sub-divide the one-facet static group-level cell into analyses 1-4
and change nothing about family or engine availability. Parallel splits (`splits == 2`)
get **no analysis number of their own**; they reuse analyses 1-5 and relabel "trial" as
"split".

**29 analyses are routed.** Two share a `case` arm with their group-level twin
(`case {11, 26}` and `case {14, 27}`). The switch has no `otherwise` arm; the selector
chain is exhaustive, so an unrouted value is unreachable.

---

## 2. Legend

| Mark | Meaning |
|---|---|
| **@** | Shipped and supported |
| **~** | Shipped, but carries a status caveat (EXPERIMENTAL or DEVELOPMENT-ONLY); see section 5 |
| **o** | Not built. Open work, or simply unwired; no decision bars it |
| **X** | Barred by a closed decision (scientific, estimand, or structural); see section 4 |
| **-** | Not applicable: the axis is meaningless for this cell |

Column headings: **Gau** = Gaussian; **g-nu** = gamma, log-nu parameterization
(`gammascale=1`); **g-LS** = gamma, location-scale (`gammascale=2`); **Stan** = CmdStan;
**HMC** = native MATLAB HMC; **flme** = native `fitlme`.

**An engine column says which engine IMPLEMENTS a design, not what the host machine needs to
run it.** On the default CmdStan path the Statistics and Machine Learning Toolbox is required
for exactly three designs: the gamma (scaled chi-square) **concurrent** difference scores
(`family='gamma'` with `diffrescor=2`) for analyses **7**, **8**, and **10**, whose
Gaussian-copula extractors call `normcdf` and `gaminv`. Every other design runs on base MATLAB
via CmdStan. (**B39**, which found the toolbox was needed for every design, was RESOLVED
2026-08-21: the two `nansum` call sites in the bundled MatlabStan `StanFit.m` were replaced with
`sum(...,'omitnan')`, and `prctile` was removed from `psis/psislw.m` under **B40**.) Read a `@`
in the Stan column as "this design has a CmdStan implementation"; see
`dependencies_support_matrix.md` for the full per-machine requirement table.

---

## 3. The matrix

### 3A. Single score: 8 cells, all shipped

The residual-coupling axis does not apply: there is only one event.

| Facets | Level | Moderation | # | Gau | g-nu | g-LS | Stan | HMC | flme |
|---|---|---|---|---|---|---|---|---|---|
| one | group | static | 1-4 | @ | @ | X [a] | @ | o | @ |
| one | subject | static | 6 | @ | @ | @ | @ | @ | X [b] |
| one | group | dynamic | 11 | @ | @ | @ | @ | @ | X [b] |
| one | subject | dynamic | 26 | @ | @ | @ | @ | o [c] | X [b] |
| two | group | static | 5 | @ | @ | X [a] | @ | @ | @ |
| two | subject | static | 25 | @ | @ | @ | @ | o [c] | X [b] |
| two | group | dynamic | 14 | @ | @ | @ | @ | @ | X [b] |
| two | subject | dynamic | 27 | @ | @ | @ | @ | o [c] | X [b] |

### 3B. Two-event difference score: 24 cells, 16 shipped

| Facets | Level | Moderation | Coupling | # | Gau | g-nu | g-LS | Stan | HMC | flme |
|---|---|---|---|---|---|---|---|---|---|---|
| one | group | static | none | 7 | @ | @ | @ | @ | @ | X [b] |
| one | group | static | population rho | 7* | @ | @ / ~ [d] | X [e] | @ | @ | X [b] |
| one | group | static | per-subject rho | - | - | - | - | - | - | - |
| one | subject | static | none | 8 | @ | @ | @ | @ | @ | X [b] |
| one | subject | static | population rho | - | o [f] | - | - | - | - | - |
| one | subject | static | per-subject rho | 8* | @ | ~ [g] | X [e] | @ | @ | X [b] |
| one | group | dynamic | none | 12 | @ | @ | @ | @ | @ | X [b] |
| one | group | dynamic | population rho | 17 | @ | o | o | @ | @ | X [b] |
| one | group | dynamic | per-subject rho | - | - | - | - | - | - | - |
| one | subject | dynamic | none | 13 | @ | X [h] | @ | @ | @ | X [b] |
| one | subject | dynamic | population rho | 19 | @ | X [h] | @ [i] | @ | @ [j] | X [b] |
| one | subject | dynamic | per-subject rho | 21 | @ | o | o | @ | @ | X [b] |
| two | group | static | none | 10 | @ | @ | @ | @ | @ | X [b] |
| two | group | static | population rho | 10* | @ | @ / ~ [d] | X [e] | @ | @ | X [b] |
| two | group | static | per-subject rho | - | - | - | - | - | - | - |
| two | subject | static | none | - | X [k] | X [k] | X [k] | - | - | - |
| two | subject | static | population rho | - | X [k] | X [k] | X [k] | - | - | - |
| two | subject | static | per-subject rho | - | X [k] | X [k] | X [k] | - | - | - |
| two | group | dynamic | none | 15 | @ | o | o | @ | @ | X [b] |
| two | group | dynamic | population rho | 18 | @ | o | o | @ | @ | X [b] |
| two | group | dynamic | per-subject rho | - | - | - | - | - | - | - |
| two | subject | dynamic | none | 16 | @ | o | o | @ | @ | X [b] |
| two | subject | dynamic | population rho | 20 | @ | X [h] | ~ [l] | @ | @ [j] | X [b] |
| two | subject | dynamic | per-subject rho | 22 | @ | o | o | @ | @ | X [b] |

`7*`, `8*`, `10*` are the same analysis number with `diffrescor=2`; concurrency selects
a different builder inside the case rather than a new analysis number. In the dynamic
branch concurrency *does* change the number (12 -> 17, 13 -> 19, 15 -> 18, 16 -> 20, and
`ssrescor=2` gives 19 -> 21, 20 -> 22).

### 3C. Four-event difference-of-differences: 16 cells, 3 shipped

| Facets | Level | Moderation | Coupling | # | Gau | g-nu | g-LS | Stan | HMC | flme |
|---|---|---|---|---|---|---|---|---|---|---|
| one | group | static | none | 9 | @ | @ | @ | @ | @ | X [b] |
| one | group | static | concurrent | - | X [m] | X [m] | X [m] | - | - | - |
| one | subject | static | either | - | X [n] | X [n] | X [n] | - | - | - |
| two | group | static | either | - | X [k] | X [k] | X [k] | - | - | - |
| two | subject | static | either | - | X [k,n] | X [k,n] | X [k,n] | - | - | - |
| one | group | dynamic | either | - | X [o] | X [o] | X [o] | - | - | - |
| one | subject | dynamic | none | 28 | @ | X [h] | @ | @ | o [c] | X [b] |
| one | subject | dynamic | concurrent | - | X [r] | X [r] | X [r] | - | - | - |
| two | group | dynamic | either | - | X [o] | X [o] | X [o] | - | - | - |
| two | subject | dynamic | none | 29 | @ | X [h] | @ | @ | o [c] | X [b] |
| two | subject | dynamic | concurrent | - | X [r] | X [r] | X [r] | - | - | - |

### 3D. Split scores: 4 cells, all shipped, Gaussian only

Splits are a Stage-1 design: non-difference, group-level, non-dynamic only (the guard
rejects any other combination before an analysis number is assigned).

| Splits | Facets | # | Gau | g-nu | g-LS | Stan | HMC | flme |
|---|---|---|---|---|---|---|---|---|
| parallel | one | 1-4 relabeled | @ | X [s] | X [s] | @ | o | @ |
| parallel | two | 5 relabeled | @ | X [s] | X [s] | @ | @ | @ |
| nonparallel | one | 23 | @ | X [s] | X [s] | @ | o | @ |
| nonparallel | two | 24 | @ | X [s] | X [s] | @ | o | @ |

Analysis 23 additionally rejects group/event facets; analysis 24 is group/event aware.
Neither offers a D-study projection, because they are observed designs.

### Footnotes

- **[a]** Analyses 1-5 have no location-scale variant **by decision, not as unbuilt
  work**: they are single-mean designs whose coefficients are unchanged by the
  parameterization. Stated in the accept-list guard's own error text.
- **[b]** Not `fitlme`-doable. `fitlme` fits single-residual plain-Gaussian mixed models;
  every cell marked here needs a per-subject heteroscedastic residual, correlated random
  effects, or a multivariate residual. Recorded per design key in
  `subroutines/estimation/psyrat_native_supported.m`.
- **[c]** CmdStan-only "for now": the case passes `native_key = ''`, which refuses a
  native fallback rather than silently fitting a different model. Analyses 25, 26, 27,
  and the Gaussian arms of 28/29 (owner ruling 2026-08-20: the native HMC port of the
  four-margin model is a separately scoped follow-up).
- **[d]** Group-level concurrent gamma has two copula parameterizations, selected by
  `dispersion`: `2` = global fixed nu (the default for this design, tractable) and `1` =
  per-person log-nu (**EXPERIMENTAL**, S11).
- **[e]** Location-scale is barred for the concurrent static designs: the observed-scale
  residual covariance of the Gaussian-copula scaled-chi-square model has not been defined
  for the direct residual-SD parameterization. Analyses 19/20 escape this via the modular
  cut-copula carve-out.
- **[f]** The static subject-level one-facet difference offers **no population-rho**
  concurrent variant, because its concurrent arm is per-subject rho by construction. The dynamic
  sibling offers both (19 population, 21 per-subject). Asymmetry, not a defect; see
  section 6.
- **[g]** Subject-level reliability requires per-person nu, because global fixed
  dispersion is a group-level estimand. This design therefore has **no non-experimental
  gamma parameterization**.
- **[h]** Location-scale only. The log-nu dispersion model is an **absent derivation**,
  not an unbuilt branch: the per-participant conversion these designs need is defined on
  a direct residual SD.
- **[i]** Modular cut-copula: stage-1 Stan fits Gamma margins with no copula, the residual
  coupling is estimated post hoc and converted by quadrature. The result is a **CUT
  posterior**, disclosed through the provenance carriers.
- **[j]** Native HMC covers the **Gaussian** arm only. The gamma arm sets
  `native_key = ''`, and gamma is pinned to CmdStan upstream regardless.
- **[k]** Rejected by preflight in every family: a two-facet subject-level static
  difference, and any two-facet static difference-of-differences.
- **[l]** **DEVELOPMENT-ONLY**: the two-facet cut-copula configuration has no empirical
  anchor. See section 5.
- **[m]** The four-event difference-of-differences supports only non-concurrent events
  (residual covariance defined as exactly 0), **in every family**.
- **[n]** Rejected by preflight: difference-of-differences with subject-level error
  variance is supported only together with `dynrel=2` (analyses 28/29).
- **[o]** The dynamic difference-of-differences is **person-specific only**, so the routing
  gate rejects a group-level request.
- **[p]** Retired 2026-08-20. This letter recorded "no Gaussian model of the estimand has
  been derived"; the owner ruling of that date reopened the decision and the Gaussian
  arms shipped (recorded in the development repository's formula audit, section 26). Kept
  as a tombstone so older references to [p] resolve; no cell carries it.
- **[q]** Retired 2026-08-20 with [p]: its "a gamma-only design can never reach a native
  engine" consequence dissolved when the Gaussian arms shipped. The native columns for
  28/29 now read from [c] (Gaussian: CmdStan-only for now) and [b] (not fitlme-doable).
- **[r]** The dynamic difference-of-differences is **nonconcurrent by estimand**: all six
  pairwise residual covariances are fixed to exactly zero, so concurrency cannot be
  enabled.
- **[s]** The gamma family models **single-trial** dispersion (Var = 2*mu^2/nu). A split
  score is the MEAN of n_i items, so the scaled-chi-square mean-variance relationship does
  not hold on split means. **Closed decision, not a deferral. Do not reopen.**

---

## 4. Gap register

### 4.1 Gamma gaps that are buildable

**Analyses 15, 16, 17, 18, 21, 22.** These are the only cells where "not available under
chi-square" means open work rather than a closed decision. All six already exist and are
fully tested under the Gaussian family, so none is a new design. Each needs the standing
seven-step increment recipe, whose scientific core is the observed-scale conversion.

Scope, blockers, difficulty ratings, and the derivation roadmap are maintained
in the development repository. Summary of what blocks each:

| # | Design | Blocker |
|---|---|---|
| 15 | two-facet group-level dynamic difference | no empirical anchor; owner decision pending |
| 16 | two-facet subject-level dynamic difference | 15, plus the conditional two-facet difference estimand |
| 17 | one-facet group-level concurrent dynamic difference | copula viability (S11); no group-level conversion route |
| 18 | two-facet group-level concurrent dynamic difference | same as 17 |
| 21 | one-facet subject-level per-subject-rho difference | S11; no estimand defined anywhere; won't-do candidate |
| 22 | two-facet subject-level per-subject-rho difference | same as 21 |

Also open, and listed there rather than here: the concurrent four-event
difference-of-differences ("9 + rescor"), which has **no template in either family** and
needs an estimand decision first.

### 4.2 Closed decisions

These are **not** backlog items. Each was decided on the merits and should not be reopened
without a new owner ruling.

| Cell | Decision |
|---|---|
| Splits under gamma (23, 24, and parallel splits) | Permanently excluded. A split mean is not a single trial; the scaled-chi-square mean-variance relationship does not hold on it |
| Location-scale for analyses 1-5 | No LS variant. Single-mean designs; coefficients are parameterization-invariant |
| Log-nu for analyses 13, 19, 20, 28, 29 | Absent derivation. These conversions are defined on a direct residual SD |
| Location-scale for concurrent 7, 8, 10 | The copula model's observed-scale residual covariance is undefined for the direct residual-SD parameterization |
| Gamma under any native engine | The native engines are family-blind; the pin exists so a gamma request cannot silently be fitted as Gaussian |

One former row is gone from this table: "Gaussian for analyses 28, 29 - no Gaussian model
of the estimand has been derived" was reopened by owner ruling 2026-08-20 and both Gaussian
arms shipped (recorded in the development repository's formula audit, section 26;
CmdStan-only per footnote [c]). Its old footnotes [p]/[q] are kept as tombstones.

### 4.3 Barred in every family

Design cells with no implementation in Gaussian **or** gamma, rejected upstream:

- Two-facet subject-level **static** two-event difference (preflight)
- Two-facet **static** difference-of-differences (preflight)
- Subject-level **static** difference-of-differences (preflight)
- **Concurrent** difference-of-differences (the case-9 guard, in every family)
- **Group-level dynamic** difference-of-differences (the routing gate; person-specific only)
- Per-subject rho outside the subject-level concurrent **dynamic** difference design

---

## 5. Status flags

Four shipped cells carry a caveat that a bare availability mark would hide: the gamma
concurrent arms of analyses 7, 8 and 10, and analysis 20. They fall into three kinds.

**EXPERIMENTAL: per-person log-nu copula (`dispersion=1`), analyses 7, 8, 10 concurrent.**
Convergence has **not** been established in this repository: one in-repo subject-level fit
reached R-hat of about 2.5 with ESS of about 2, and did not improve when simulated
participant and trial counts were increased. Analysis 8 has no non-experimental
alternative, because subject-level reliability requires per-person nu. Disclosed in
`psyrat_provenance_lines` and in the export headers; **there is no on-screen label**.

The `dispersion=2` (global fixed nu) default for analyses 7 and 10 is the tractable
alternative but is itself hedged in its own disclosure text as establishing
"tractability, not robustness". It rests on one dataset at one choice of sampler settings, with the
two-facet check still producing divergent transitions.

**DEVELOPMENT-ONLY: analysis 20 (two-facet subject-level concurrent dynamic difference,
gamma location-scale).** The two-facet configuration has no empirical anchor: the
reference bundle's own two-facet runs were prepared against a simulated second occasion
and were never sampled, so unlike the one-facet design there is no completed external fit
to replay against. **This status currently has no user-facing carrier.**

**CUT POSTERIOR: analyses 19 and 20 under gamma.** The modular cut-copula workflow fits
the margins without the copula and estimates the residual coupling post hoc, so the result
is a cut posterior, not a full Bayesian one. Analyses 28/29 are explicitly **not** cut
posteriors, because they have no dependence module at all: their residual covariances are
zero by estimand.

The toolbox version carries `-beta` **because of** the copula designs: the API and outputs
are stated to be still moving while analyses 7/8/10 ship experimental with convergence not
validated in-repo.

---

## 6. Known asymmetries

**Analysis 8's concurrent arm fits a different coupling structure in each family.** Under
Gaussian it estimates a **per-subject** correlation, `rho[s] = tanh(rescor_mu +
sd_rescor * z_rescor[s])`. Under gamma it delegates to analysis 7's **population**
Gaussian-copula model and folds the per-person copula residual cross-covariance into the
subject-level residual downstream. Same user-facing design, two coupling models. This is
deliberate and disclosed in the case comment; it is recorded here so nobody "fixes" one
side toward the other (Golden Rule 1).

**Only one design has two native engines.** Analysis 5 (`trt`) is the sole design both
native engines implement, so the HMC-before-fitlme preference branch of the auto cascade is
exercised by exactly one design in production.

**`native_key = ''` is overloaded.** It marks both "CmdStan-only by design" and "typo or
unregistered key", and both resolve through the same empty-engine-set path, so the overload
itself still stands.

**The testing gap it created is CLOSED (2026-08-20).**
`TestStanModelSnapshots.testNativeKeyMatchesCaseSite` now captures the key each analysis
actually passes to `psyrat_run_stan` and pins it three ways: against a golden map covering
every scenario, against an allowlist of the analyses permitted to pass an empty key, and
against the `psyrat_native_supported` registry. A mistyped key now fails the unit lane
instead of silently degrading that analysis to CmdStan-only. This closes the case-site leg
of the three-way contract `TestNativeEngine.testRegistryDispatchConsistency` documents.

---

## 7. Validation evidence

Most shipped cells have live CmdStan known-truth recovery coverage in `tests/integration/`.
The exceptions, recorded rather than closed for the beta (owner ruling 2026-09-05): the
Gaussian four-event difference-of-differences (analysis 9) is recovered live only through
the native HMC engine (`TestDodNativeHmcRecovery`; the integration factory carries no
four-event scenario); the Gaussian arms of analyses 26 and 27 have no recovery test of their
own (the qualification below); and the gamma arm of analysis 19 is checked only by
`TestModularCutCopulaExternalReplay`, a replay of an external reference fit that runs no
sampler and skips without the bundle. Coverage is otherwise complete for the gamma family:
analyses 1-4 (marginal estimand), 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 20, 25, 26, 27, 28 and
29 each have at least one gamma recovery test, with both parameterizations covered
separately wherever both are legal (the `*LsRecovery` twins).

The GAUSSIAN arms of 28/29 (shipped 2026-08-20) have their own live recovery tests,
`TestDodDynrelRecovery` and `TestDodDynrelTrtRecovery`, both load-bearing and both observed
PASS on 2026-08-20 (28: max R-hat 1.0000, G(0) recovered 0.857 vs true 0.829, per-person
rank correlation 0.881 at 2x750/750; 29: max R-hat 1.0000 over 161 parameters, CES G(0)
0.861 vs true 0.835, rank correlation 0.916 at the gamma-29 4x1500/1000 budget; occasion
SDs are prior-dominated at nocc=2 and stay printed diagnostics, the known weak
identification of occasion-level components). The gamma-28 recovery was re-run in the same
session and passed, pinning live non-regression of the gamma arm across the increment.

**Do not audit this coverage by filename.** Several two-facet analyses are tested inside
the file named for their one-facet sibling, because the two designs share a test class:
analysis 18 lives in `TestDiffDynrelRescorRecovery`, 20 in `TestSubjDiffDynrelRescorRecovery`,
and 22 in `TestSubjDiffDynrelSsRescorRecovery`. A filename sweep reports those three as
uncovered when they are not.

**The qualification: analyses 26 and 27 have no Gaussian recovery test of their own.** They
reuse the case-11 and case-14 Stan models verbatim, so the FIT is covered by those analyses'
recovery tests and the draws are identical; what 26/27 add is a downstream per-participant
read-out, which is covered offline rather than by a live run. That is a reasonable division,
but it means a change to the read-out alone would not be caught by any live CmdStan test.

Three cells carry **external** grounding, meaning a replay against the owner's completed reference
fits rather than against simulated truth. That is stronger evidence than self-consistent
recovery, and it covers: analysis 13, analysis 19, and analysis 28 — in every case the GAMMA
arm. **No external anchor exists for the Gaussian DoD arms, and none can be minted**: the
reference bundle is gamma end to end, so where the gamma family was anchored at 28 and
unanchored at 29, the Gaussian family is unanchored at both. Their in-repo evidence ladder
(owner ruling 2026-08-20) is live planted-truth recovery for BOTH analyses as the
load-bearing gate, plus constant-draw longhand closed forms, the identity-link theorem
pins, and a full-space tie to the vetted static Gaussian DoD kernel.

Two caveats on reading a green suite, both standing:

- The accuracy suite's **regression class is regression-only**: `TestCalculationAccuracyOracle`'s
  oracle duplicates the production formulas, so its green means "unchanged", not
  "scientifically correct". The suite's other two classes (the independent oracle restating
  Rocha et al. (2026) Tables 2, 3 and 6, and the external benchmark) are the corroboration; the
  splits and dynamic designs have no independent-oracle row.
- Native-engine agreement with CmdStan is validated only for the designs recorded in
  `native_engine_validation.md`, namely analyses 5, 6, 7 and 10. The remaining native HMC designs
  have recovery tests but no CmdStan-vs-native posterior comparison.

---

## 8. Where to look next

| Question | File |
|---|---|
| Which engine will actually run, and what does it need? | `dependencies_support_matrix.md` |
| How close is native HMC to CmdStan? | `native_engine_validation.md` |
| What do the split designs estimate? | the in-app Formula Guide (Rocha Table 4/5 entries) |
| How do I run one of these designs? | [the user manual](manual/index.md) |
| What changed, and what is known to be incomplete? | `../CHANGELOG.md` |

The development trackers that record current status, open findings, and the derivation
behind each formula are internal working documents and are not distributed with the
release.
