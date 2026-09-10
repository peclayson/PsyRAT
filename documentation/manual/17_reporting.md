# 17. Interpreting and reporting results

A dependability estimate earns its place in a manuscript at the point where it
changes something: how many trials a task presents, which participants enter an
analysis, or whether a score is treated as adequate for the question being
asked. Psychometric reporting is part of the analysis rather than an appendix to
it, and the case for treating it that way has been made at length for
psychophysiology (Clayson, 2024). The material below covers the decision rules
that a reliability estimate can legitimately support, paragraph templates for a
methods section, and the record a reviewer needs in order to check the work.

## 17.1 Decision rules

### What counts as adequate dependability

For internal consistency, a coefficient of .80 is the working standard once a
measure has proven promising, with group-difference research held to it, and
.70 is the floor for early-stage or exploratory work (Clayson & Miller, 2017b,
p. 62). The reason is not convention. Reliability bounds the correlation an
observed score can show with an external variable, and it lowers statistical
power as it falls. Against a criterion whose own reliability is .90, a score
with a reliability of .70 cannot correlate above about .79, whereas a score
with a reliability of .80 can reach about .85 (Clayson & Miller, 2017b, p. 61,
Table 1; the attenuation bound is defined on relative-error reliability, so
substituting a dependability estimate gives a conservative bound). A study
whose scores sit near the floor is therefore spending sample size on
measurement error, and a null result
from such a study is difficult to interpret.

The exception runs upward, not downward. When a decision is made about an
individual score rather than about a group mean, the standard summarized by
Clayson and Miller (2017b, p. 62) sets .90 as a bare minimum and prefers .95.
Screening and treatment selection both fall under that heading, and the gap
between a coefficient adequate for a group contrast and one adequate for an
individual decision is the reason psychometric readiness is a precondition for
biomarker claims rather than a consequence of them (Clayson, 2025).

The illegitimate reason to move a threshold is that too few participants survive
it. Clayson and Miller (2017b, p. 62) note that the .70 standard is frequently
cited to justify a single cutoff for every stage of research, which the source
of that standard does not support. The threshold belongs to the question, and it
should be fixed and recorded before any effect of interest is inspected. A
threshold chosen after inspecting how much data it costs is a researcher degree
of freedom, and reporting it as though it had been set in advance misdescribes
the analysis.

### Planning trial counts from the D-study projection

The dependability-versus-trials plot and the trial-cutoff table project the
coefficient forward under the variance components estimated from the data in
hand, and they answer a design question: how many trials per participant would
be needed to reach the stated threshold. The projection holds the estimated
components fixed, so it describes what dependability would be if scores at the
larger trial count behaved like the scores observed. It does not certify that a
future sample will reach the projected value, because reliability is a property
of scores in a context rather than of the measure (Clayson & Miller, 2017b), and
a dependability estimate obtained in one sample can fail to transfer to another
(Baldwin et al., 2015).

A trial cutoff reported as -1 means that no trial count within the projected
range reached the specified threshold. That result is informative rather than
broken. It says the between-person variance in these scores is small relative
to the within-person variance, and that the remedy lies in the task or the
scoring pipeline rather than in collecting more trials of the same kind. The
projection is also what lets you satisfy the third of the three reporting
guidelines Clayson and Miller (2017b, p. 62) set out, reporting how the trial
cutoff for data exclusion was justified, from your own data rather than from a
number inherited from a prior study; applied test-retest work has followed the
same practice (Carbine et al., 2021).

### Participant inclusion and the reporting of exclusions

Subject-level estimates give a dependability coefficient for each participant
rather than one coefficient for the group, so they support an inclusion rule
that a group-level estimate cannot: a participant whose own scores are too
imprecise for the intended use can be identified as such, even when the group
coefficient is adequate (Clayson, Brush, & Hajcak, 2021). The toolbox reports the
retained and excluded identifiers together with the coefficient that produced
the split, and both belong in the manuscript.

Three things should be reported, and reported in this order relative to the
analysis. The threshold and the rule come first, fixed before the effects of
interest are examined. The counts come second: how many participants met the
rule, how many did not, and what the trial distribution looked like on each side
of the line. The estimates come third. Reporting the exclusions without the
threshold, or the threshold without the count, leaves a reader unable to tell
whether the inclusion rule was a measurement decision or an outcome-dependent
one. Where the analysis was preregistered, the threshold is the element most
worth registering, because it is the one whose post hoc adjustment is hardest to
detect from the published record.

## 17.2 Methods paragraphs

Each block below is a complete paragraph with bracketed slots for the reader's
own numbers. Replace every bracketed span, delete the alternatives that do not
apply, and check the resulting coefficient name against the label printed in the
exported table header, which states the estimand the numbers belong to. Every
version number, chain count, and seed value should be copied from the run rather
than from this page: the toolbox version appears in the startup banner and in
the header of every exported table, and the CmdStan version and the sampler
settings appear in the run-configuration sidecar described in Section 17.3.

**Single-session (internal consistency, one facet).**

```
Internal consistency of [COMPONENT] scores was evaluated within a
generalizability-theory framework treating trials as the sole facet of
measurement (Rocha et al., 2026). Variance components were estimated in a
Bayesian framework with the PsyRAT Toolbox (version [X.Y.Z];
https://github.com/peclayson/PsyRAT), which fit the models in CmdStan
[CMDSTAN VERSION] using [4] chains with [5,000] warmup and [5,000] sampling
iterations per chain and a random seed of [12345]. Dependability (the
absolute-error coefficient) of [COMPONENT] scores averaged over [X] trials was
[.XX], 95% CrI [[.XX], [.XX]]. A decision study indicated that [X] trials were
required to reach the dependability threshold of [.XX], which was specified
before the data were examined. Applying that trial minimum retained [N] of [N]
participants; [N] participants were excluded for contributing too few usable
trials.
```

**Test-retest (two facets: trials and occasions).**

```
Test-retest reliability of [COMPONENT] scores across [N] sessions separated by
[INTERVAL] was evaluated within a generalizability-theory framework treating
trials and occasions as facets (Clayson, Carbine, et al., 2021; Rocha et al.,
2026). Variance components were estimated in a Bayesian framework with the PsyRAT
Toolbox (version [X.Y.Z]; https://github.com/peclayson/PsyRAT), which fit the
models in CmdStan [CMDSTAN VERSION] using [4] chains with [5,000] warmup and
[5,000] sampling iterations per chain and a random seed of [12345]. The
coefficient reported is the coefficient of [equivalence and stability (CES) /
equivalence (CE) / stability (CS)] under the [absolute-error (dependability) /
relative-error (generalizability)] definition, evaluated for a single occasion
(n'_o = 1). For scores averaged over [X] trials per occasion, the coefficient
was [.XX], 95% CrI [[.XX], [.XX]].
```

**Difference scores.**

```
Internal consistency of the [CONDITION A] minus [CONDITION B] difference score
was estimated with the generalizability-theory difference-score expressions of
Rocha et al. (2026), which model the two constituent scores jointly and carry
their covariance rather than treating them as independent (Clayson, Baldwin,
& Larson, 2021). Variance components were estimated in a Bayesian framework with
the PsyRAT Toolbox (version [X.Y.Z]; https://github.com/peclayson/PsyRAT), which
fit the models in CmdStan [CMDSTAN VERSION] using [4] chains with [5,000]
warmup and [5,000] sampling iterations per chain and a random seed of [12345].
Dependability of the difference score at [X] trials per condition was [.XX],
95% CrI [[.XX], [.XX]]; dependability of the constituent [CONDITION A] and
[CONDITION B] scores was [.XX], 95% CrI [[.XX], [.XX]], and [.XX], 95% CrI
[[.XX], [.XX]], respectively. The intraclass correlation reported alongside the
difference score is a single-observation (n' = 1) quantity derived from the same
variance components rather than a published difference-score expression.
```

**Reliability from data splits.**

```
Reliability was estimated from [N] splits of [X] items each rather than from
single trials, with the number of items per split supplied as a weight, using
the observed-design generalizability-theory expressions of Rocha et al. (2026);
see Vispoel, Xu, and Schneider (2022) for splits rather than items as the unit
of analysis in generalizability theory, and Vispoel et al. (2018) for the
broader framework. The splits were [parallel (equal items per split) / nonparallel
(unequal items per split)]. Variance components were estimated in a Bayesian
framework with the PsyRAT Toolbox (version [X.Y.Z];
https://github.com/peclayson/PsyRAT), which fit the models in CmdStan
[CMDSTAN VERSION] using [4] chains with [5,000] warmup and [5,000] sampling
iterations per chain and a random seed of [12345]. Dependability (absolute
error) was [.XX], 95% CrI [[.XX], [.XX]]. Because only split means were
observed, the item main effect within splits could not be separated from the
person-by-split residual; the generalizability (relative-error) coefficient is
therefore a lower bound and the relative standard error of measurement an upper
bound, whereas dependability is exact.
```

**Dynamic (conditional) reliability.**

```
Reliability was modeled as a function of [DIMENSION] using the mixed-effects
location-scale extension of generalizability theory described by Rast and
Clayson (in press), in which the residual standard deviation varies
systematically across measurement conditions and persons rather than being held
constant. [DIMENSION] was standardized across the full sample, so a value of
zero corresponds to the sample mean. Variance components were estimated in a
Bayesian framework with the PsyRAT Toolbox (version [X.Y.Z];
https://github.com/peclayson/PsyRAT), which fit the models in CmdStan
[CMDSTAN VERSION] using [4] chains with [5,000] warmup and [5,000] sampling
iterations per chain and a random seed of [12345]. At [X] trials, dependability
was [.XX], 95% CrI [[.XX], [.XX]], at one standard deviation below the mean of
[DIMENSION]; [.XX], 95% CrI [[.XX], [.XX]], at the mean; and [.XX], 95% CrI
[[.XX], [.XX]], at one standard deviation above. The posterior mean of the
residual-scale slope on [DIMENSION] was [X.XX], 95% CrI [[X.XX], [X.XX]].
```

Two habits make these paragraphs harder to misread. Name the coefficient rather
than calling it reliability, because dependability and generalizability differ
by whether the trial main effect enters absolute error, and a reader cannot
recover which was reported from the number alone. Report the interval with the
point estimate in every case, since a dependability estimate from a small sample
can sit above a threshold with an interval that spans it.

## 17.3 Reproducibility for reviewers

Four items constitute the archive for a reliability analysis. The first is the
set of exported tables with their headers intact. Each exported table carries
the toolbox version, the generation date, the estimation seed, the analysis
label, the estimation engine, whether the priors were the toolbox defaults or
customized, the convergence verdict, the divergent-transition count, and any
design-specific disclosure that applies to the estimand. Tables built from a
summary also record the credible-interval width the summary was built at; the
splits, dynamic-reliability, and difference-of-differences families build none
and carry no width line. Stripping the header
rows before archiving removes the only record that travels with the numbers. The
second is the `*_runconfig.json` sidecar written next to the saved result, which
records the seed, chains, warmup and sampling counts, any sampler overrides, the
column mapping, the design flags, the coefficient selectors in effect, the full
prior specification, the detected CmdStan version, the MATLAB version, and
content hashes of the input data. The third is the saved `.psyrat` result file,
which allows the posterior summaries to be regenerated without refitting. The
fourth is the analysis dataset itself, in the long format the toolbox reads,
where participant-level data can be shared. Where it cannot, the sidecar's data
hashes and row count still document what entered estimation.

The seed guarantees less than it appears to. Holding the data, the settings,
and the CmdStan build fixed, the same seed reproduces the same posterior draws
exactly, which is what makes a rerun on the same machine a genuine check.
Outside that combination, on a different CmdStan build or a different machine,
the draws are not guaranteed to be bit-identical, and the toolbox warns rather
than stops when a replay detects a CmdStan version different from the one
recorded. Two runs under different CmdStan builds should agree closely, but a
version mismatch is a reason to re-estimate rather than a license to compare
numbers across the two runs. A run produced by a native MATLAB engine instead
of CmdStan is Bayesian but is not bit-identical to a CmdStan run, and the
engine that produced a result is recorded in both the table header and the
sidecar. Reporting the seed without the CmdStan version, the engine, the
sampler settings, and the within-chain parallelization setting therefore claims
a reproducibility the seed alone does not deliver.

Three toolbox conventions account for most of the questions a careful reader
raises, and each is a choice the software makes rather than a result it found.
Test-retest coefficients default to a single occasion, n'_o = 1, so the reported
coefficient describes a score from one occasion rather than a composite
averaged over occasions; the multi-occasion composite is available and is labeled as such
in the output, and [Chapter 9](09_viewing-results.md) covers the selector.
Intraclass correlations reported alongside a difference score are derived from
the fitted variance components at n' = 1 rather than taken from a published
difference-score expression; the table column is labeled ICC (n'=1), and the
derivation is documented here and in the tutorial rather than in the header;
[Chapter 12](12_tutorial-difference-scores.md) and
[Chapter 15](15_advanced-designs.md) give the context. Splits analyses report an
exact dependability coefficient and a generalizability coefficient that is a
downward-biased lower bound, because split means do not identify the item main
effect within splits separately from the person-involving residual; the relative
standard error of measurement is correspondingly an upper bound, and
[Chapter 13](13_tutorial-splits.md) works through the design. Disclosing these
three in a supplement costs a paragraph and prevents a reviewer from reading a
convention as an error.

## 17.4 Reporting the toolbox itself

Cite two things: Rocha et al. (2026), which is the source of the estimators the
toolbox implements, and the tagged GitHub release of the version that produced
the estimates. The paper describes the method, and the version identifies the
code. The version string is returned by `psyrat_defineversion`, printed in the
startup banner, written into the header of every exported table, and recorded in
the run-configuration sidecar, so the correct value is available from any saved
output rather than from memory. Machine-readable metadata for the software entry
is in `CITATION.cff` in the repository, and the release page supplies the year a
reference style may require. Analyses that use the dynamic or conditional
designs should additionally cite Rast and Clayson (in press), which is the basis
for the location-scale estimation those designs rest on.

Clayson and Miller (2017a) describes the ERA Toolbox, PsyRAT's predecessor. It
is the right citation for ERA and the wrong one for PsyRAT, and the distinction
matters because the two toolboxes do not implement the same set of designs.

Bug reports and questions go to the issue tracker at
`https://github.com/peclayson/PsyRAT/issues`. A report that includes the version
string and, where the analysis can be shared, the `*_runconfig.json` sidecar is
reproducible on another machine; a report without them usually is not, because
the sidecar carries the seed, the priors, the design flags, and the data hashes
that determine the result. [Chapter 18](18_troubleshooting.md) covers the errors
that have a known cause and a documented fix, and is worth checking first.
