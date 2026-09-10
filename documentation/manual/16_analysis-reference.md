# 16. Analysis reference

One card per routed analysis, plus a card for parallel splits, which reuses
other analyses rather than receiving its own number. Use
[Chapter 5](05_choosing-an-analysis.md) to find the number; use this chapter to
check what it needs, what fits it, and where the results open.

Every card follows the same shape. Status is the maturity label from
[Chapter 5](05_choosing-an-analysis.md). Data lists the columns the design
requires and the constraints it enforces. GUI names the exact on-screen
controls that reach it, on the Specify Inputs screen and the Specify Processing
Preferences screen behind the Preferences button. Scripted lists the
`psyrat_run` name-value pairs that differ from the defaults; every call also
needs a data source, `savepath`, `savename`, and the column-mapping arguments
(`idcol` and `meascol` when your headers differ from their defaults of `id`
and `meas`, and whichever of `eventcol`, `groupcol`, `timecol`, `dim1col`,
`dim2col`, `weightcol` the design uses). Engines names which fitters
implement the design. Opens in names the results window. Notes and See carry
the caveats and the source of the formulas.

Three conventions run through every card. `ssubjrel` is the current name for
subject-level reliability, and `sserrvar` is accepted as a legacy alias.
`diffrescor` and `diffwpcov` are one setting under two names, so pass either.
The gamma family runs only on CmdStan: the Automatic engine setting is pinned
to it, and a forced Native HMC or Native fitlme selection is rejected with an
error rather than redirected.

Subject-level analyses are not collected in a section of their own. Each one
sits with the design it modifies, and its Use when line says that it reports a
coefficient per participant.

## 16.1 One-facet reliability, single occasion

These are the internal-consistency designs: persons crossed with trials, all
data from one session. They are the oldest and best-covered part of the
toolbox, and the only difference among the first four is which facet columns
you selected.

### Analysis 1: One-facet reliability

**Status:** Supported.

**Use when:** All the data come from one session, there are no group or
condition splits, and you want one dependability coefficient for the sample.

**Data:** `id` and `meas`, one row per single-trial score. No group, event, or
occasion column.

**GUI:** Specify Inputs, Participant ID and Measurement set; Group ID, Event
Type, Occasion ID, Dimension 1 (dynamic), and Items per split (n_i) all left on
none. Preferences: Estimate internal consistency of difference scores set to
No, Estimate subject-specific error variances set to No.

**Scripted:** No design flags at all. The defaults (`diffest` 1, `ssubjrel` 1,
`dynrel` 1, `splits` 1, `family` `'gaussian'`) select this analysis.

**Engines:** CmdStan, or native fitlme. There is no native HMC implementation.

**Opens in:** View Results Setup, the single-session viewer.

**Notes:** This is the design the D-study machinery was built for. You get the
dependability-against-trials curve, the trial-count cutoff, and the
participant-inclusion table. Gamma is available; the location-scale
parameterization is not, because a single-mean design gives the same
coefficients under either one.

**See:** Rocha et al. (2026), Table 2.

### Analysis 2: One-facet reliability by group

**Status:** Supported.

**Use when:** You have the analysis 1 design plus a grouping variable, and you
want a separate coefficient for each group.

**Data:** `id`, `meas`, and `group`. No event or occasion column.

**GUI:** As analysis 1, with Group ID set to the grouping column. Select Groups
restricts which groups are processed.

**Scripted:** Pass `groupcol`. Add `whichgroups` to restrict the levels.

**Engines:** CmdStan, or native fitlme.

**Opens in:** View Results Setup, the single-session viewer.

**Notes:** Each group is fitted as a separate stratum rather than as a fixed
effect, so the groups do not borrow strength from one another and a small group
gets wide intervals of its own. Everything else matches analysis 1.

**See:** Rocha et al. (2026), Table 2.

### Analysis 3: One-facet reliability by event

**Status:** Supported.

**Use when:** You have two or more conditions and want each one's reliability
separately, not the reliability of their difference.

**Data:** `id`, `meas`, and `event`. No group or occasion column.

**GUI:** As analysis 1, with Event Type set to the condition column. Select
Events restricts which events are processed. Leave Estimate internal
consistency of difference scores on No.

**Scripted:** Pass `eventcol`. Add `whichevents` to restrict the levels.

**Engines:** CmdStan, or native fitlme.

**Opens in:** View Results Setup, the single-session viewer.

**Notes:** Any number of events is allowed here, unlike the difference designs,
which require exactly two or exactly four. If you want the difference rather
than the constituents, you want analysis 7.

**See:** Rocha et al. (2026), Table 2.

### Analysis 4: One-facet reliability by group and event

**Status:** Supported.

**Use when:** Both a grouping variable and a condition variable are present and
you want a coefficient for each combination.

**Data:** `id`, `meas`, `group`, and `event`. No occasion column.

**GUI:** As analysis 1, with both Group ID and Event Type set.

**Scripted:** Pass both `groupcol` and `eventcol`.

**Engines:** CmdStan, or native fitlme.

**Opens in:** View Results Setup, the single-session viewer.

**Notes:** The number of strata is the product of the two level counts, so a
design with three groups and four events fits twelve separate models and takes
roughly twelve times as long. Trim with Select Groups and Select Events before
committing to a long run.

**See:** Rocha et al. (2026), Table 2.

### Analysis 6: Subject-level one-facet reliability

**Status:** Supported.

**Use when:** You want a dependability coefficient for each participant rather
than one for the sample, from a single session.

**Data:** `id` and `meas`; optional `group` and `event`. No occasion column.

**GUI:** As analysis 1, with Preferences, Estimate subject-specific error
variances set to Yes.

**Scripted:** `'ssubjrel', 2`.

**Engines:** CmdStan, or native HMC. There is no fitlme implementation,
because the model gives each participant their own residual standard
deviation.

**Opens in:** View Results Setup, the single-session viewer.

**Notes:** The fit is slower than analysis 1 by more than the parameter count
suggests, since per-participant scale is harder to identify than a pooled one.
Participants with few trials get wide per-person intervals, and that is the
honest answer rather than a defect.

**See:** Rocha et al. (2026), Table 2, applied per participant (the toolbox's
convention for the per-person read-out).

## 16.2 Test-retest reliability

Adding an occasion column crosses persons with trials and occasions. Both cards
here are static and non-difference; the difference and dynamic two-facet
designs live in later sections.

### Analysis 5: Test-retest reliability

**Status:** Supported.

**Use when:** Participants were measured on two or more occasions and you want
coefficients that generalize over occasions as well as trials.

**Data:** `id`, `meas`, and `time`; optional `group` and `event`.

**GUI:** Specify Inputs, Occasion ID set to the occasion column. Select
Occasions restricts which occasions are processed. Preferences: difference
scores No, subject-specific error variances No.

**Scripted:** Pass `timecol`. Add `whichtimes` to restrict the levels.

**Engines:** CmdStan, native HMC, or native fitlme. This is the only design
both native engines implement, and the automatic cascade prefers HMC.

**Opens in:** View Results Setup, the test-retest viewer.

**Notes:** The viewer reports coefficients of equivalence, of stability, and of
equivalence and stability, and you pick which one the tables use. Occasion-level
variance components are estimated from as many levels as you have occasions, so
at two occasions they stay wide and are prior-influenced. Group and event
columns are fitted as separate strata, as in analyses 2 to 4.

**See:** Rocha et al. (2026), Table 3.

### Analysis 25: Subject-level test-retest reliability

**Status:** Supported.

**Use when:** You want a per-participant coefficient from a multi-occasion
design, without difference scores and without a dimension.

**Data:** `id`, `meas`, and `time`; optional `group` and `event`.

**GUI:** As analysis 5, with Preferences, Estimate subject-specific error
variances set to Yes.

**Scripted:** Pass `timecol` and `'ssubjrel', 2`.

**Engines:** CmdStan only. No native engine implements it yet, and the toolbox
refuses a native fallback rather than fitting a different model quietly.

**Opens in:** View Results Setup, the test-retest viewer.

**Notes:** This is the direct two-facet analog of analysis 6. It inherits the
subject-level participant-inclusion contract rather than the test-retest one,
which matters if you drive the summary layer from a script. Criterion Score
Outputs is not offered for this design (Section 9.7).

**See:** Rocha et al. (2026), Table 3, applied per participant (the toolbox's
convention for the per-person read-out).

## 16.3 Difference scores

Two conditions are modeled together, and the reported coefficient belongs to
their difference. The covariance between the constituents enters the
universe-score variance, which is why a difference is often much less reliable
than either score behind it.

### Analysis 7: One-facet difference-score reliability

**Status:** Supported under the Gaussian family. Under the gamma family with
correlated residuals, the per-person dispersion arm (`'dispersion', 1`) is
Experimental; the default single-dispersion arm is tractability-checked but
not robustness-checked.

> EXPERIMENTAL under gamma with per-participant dispersion (`dispersion` 1):
> convergence has not been validated in this release, and exported tables from
> that arm carry an EXPERIMENTAL line in their headers. The default
> global-dispersion arm (`dispersion` 2) carries a header line saying its one
> check established tractability, not robustness. Gaussian runs of this design
> are Supported and carry no such line.

**Use when:** One session, exactly two conditions, and the score of interest is
the difference between them.

**Data:** `id`, `meas`, and `event` with exactly two unique event types;
optional `group`. No occasion column.

**GUI:** Specify Inputs, Event Type set to the condition column, with exactly
two events selected under Select Events. Preferences: Estimate internal
consistency of difference scores set to Yes: 2-event difference. Set Estimate
residual covariance for co-occurring events to Yes when both conditions come
from the same trials.

**Scripted:** `'diffest', 2` (and `'diffrescor', 2` for the concurrent
variant). Under gamma with residual covariance on, `'dispersion', 1` selects
the experimental per-person parameterization; the default for that variant is
`2`. The `dispersion` setting does not apply to the non-concurrent variant.

**Engines:** CmdStan, or native HMC.

**Opens in:** View Results Setup, the single-session viewer.

**Notes:** The concurrent variant keeps the same analysis number and changes
the model behind it, so a run with residual covariance on is not comparable to
one with it off. Under gamma, the concurrent variant additionally requires the
Statistics and Machine Learning Toolbox even on the CmdStan path, and the
missing-function failure lands after sampling has finished.

**See:** Rocha et al. (2026), Table 6.

### Analysis 8: Subject-level one-facet difference-score reliability

**Status:** Supported under the Gaussian family. Experimental under the gamma
family with correlated residuals; the non-concurrent gamma variant is
Supported.

> EXPERIMENTAL under gamma with correlated residuals: that variant's only gamma
> arm uses per-participant dispersion, convergence has not been validated in
> this release, and its exported tables carry an EXPERIMENTAL line in their
> headers. The non-concurrent gamma variant and Gaussian runs of this design
> are Supported and carry no such line.

**Use when:** You want each participant's own difference-score reliability from
a single session.

**Data:** `id`, `meas`, and `event` with exactly two unique event types;
optional `group`. No occasion column.

**GUI:** As analysis 7, with Preferences, Estimate subject-specific error
variances set to Yes.

**Scripted:** `'diffest', 2, 'ssubjrel', 2` (and `'diffrescor', 2` for the
concurrent variant).

**Engines:** CmdStan, or native HMC.

**Opens in:** View Results Setup, the single-session viewer.

**Notes:** Under the Gaussian family the concurrent variant estimates a
residual correlation for each participant, unlike analysis 7's single
population correlation; under gamma it delegates to analysis 7's population
copula and folds the per-person cross-covariance in downstream. The gamma
concurrent variant also requires the Statistics and Machine Learning Toolbox
even on the CmdStan path, and a missing toolbox surfaces only after sampling
finishes. This design has no non-experimental gamma parameterization: a
per-participant coefficient needs a per-participant dispersion parameter, and a
single sample-wide dispersion is a group-level estimand. Run it under the
Gaussian family unless you have a specific reason not to.

**See:** Rocha et al. (2026), Table 6.

### Analysis 10: Test-retest difference-score reliability

**Status:** Supported under the Gaussian family. Under the gamma family with
correlated residuals, the per-person dispersion arm (`'dispersion', 1`) is
Experimental; the default single-dispersion arm is tractability-checked but
not robustness-checked.

> EXPERIMENTAL under gamma with per-participant dispersion (`dispersion` 1):
> convergence has not been validated in this release, and exported tables from
> that arm carry an EXPERIMENTAL line in their headers. The default
> global-dispersion arm (`dispersion` 2) carries a header line saying its one
> check established tractability, not robustness. Gaussian runs of this design
> are Supported and carry no such line.

**Use when:** Two conditions measured on two or more occasions, and the score
of interest is the difference.

**Data:** `id`, `meas`, `time`, and `event` with exactly two unique event
types; optional `group`.

**GUI:** As analysis 7, with an Occasion ID column assigned on the input
screen (subject-level reliability off).

**Scripted:** Pass `timecol` with `'diffest', 2` (and `'diffrescor', 2` for the
concurrent variant; under gamma the same `'dispersion'` arms and the same
Statistics and Machine Learning Toolbox requirement as analysis 7 apply).

**Engines:** CmdStan, or native HMC.

**Opens in:** View Results Setup, the test-retest viewer.

**Notes:** Six crossed variance components are estimated for the difference,
and the facet main effects (occasions especially) are the weakly identified
ones at realistic occasion counts. Read the person and person-by-facet
components first. Criterion Score Outputs is not offered for this design
(Section 9.7).

**See:** Rocha et al. (2026), Table 6, two-facet form.

## 16.4 Difference of differences

One four-event contrast, (ERP_1 - ERP_2) - (ERP_3 - ERP_4), reported as a
single score. The static form is here; the person-specific dynamic forms have
their own section.

### Analysis 9: Difference-of-differences reliability

**Status:** Supported.

**Use when:** The quantity of interest is a contrast between two difference
scores, such as a condition effect compared across two tasks or two electrode
sites.

**Data:** `id`, `meas`, and `event` with exactly four unique event types;
optional `group`. No occasion column. The four events must be assigned to
ERP_1 through ERP_4 before the run.

**GUI:** Preferences, Estimate internal consistency of difference scores set to
Yes: 4-event DoD, and Estimate residual covariance for co-occurring events set
to No. On Specify Inputs, Event Type set to the condition column with exactly
four events selected, then click DoD Contrast... and assign four unique events
on the Configure DoD Contrast screen.

**Scripted:** `'diffest', 3` plus `'dodmap', {...}`, a 1-by-4 cell of event
labels in ERP_1 to ERP_4 order.

**Engines:** CmdStan, or native HMC.

**Opens in:** DoD Reliability Setup.

**Notes:** Residual covariance is unavailable in every family, and a concurrent
request is rejected rather than ignored: the contrast is defined with all four
constituents treated as measured on distinct trials. The default contrast
weights are 1, -1, -1, 1. The intraclass correlations this design reports are
derived quantities rather than entries in the published difference-score table,
and the exported headers label them that way.

**See:** Toolbox convention, built on the Rocha et al. (2026) Table 6
difference-score partitioning.

## 16.5 Dynamic reliability

Selecting a Dimension 1 column makes the mean and the residual scale functions
of a standardized between-person covariate, so the result is a reliability
surface over that dimension rather than a single number. The dimension is
standardized once across the whole sample before any group or event stratum is
split off, so zero is the overall mean everywhere. All fourteen cards in this
section open in the same viewer.

### Analysis 11: Dynamic one-facet reliability

**Status:** Supported.

**Use when:** One session, one score, and you want reliability reported as a
function of a participant characteristic.

**Data:** `id`, `meas`, and a continuous between-person `dim1` (one value per
participant); optional `dim2`, `group`, and `event`. No occasion column.

**GUI:** Specify Inputs, Dimension 1 (dynamic) set to the covariate column, and
optionally Dimension 2 (dynamic). Preferences: difference scores No,
subject-specific error variances No.

**Scripted:** Pass `dim1col` (and `dim2col`) with `'dynrel', 2`.

**Engines:** CmdStan, or native HMC.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Dimension 2 adds a second predictor and its interaction with the
first, which turns the reported surface two-dimensional and roughly doubles
what the viewer has to render. Group and event columns are fitted as separate
strata.

**See:** Rast and Clayson (in press).

### Analysis 26: Subject-level dynamic one-facet reliability

**Status:** Supported.

**Use when:** You want each participant's own conditional reliability from the
one-facet dynamic model rather than one population surface.

**Data:** As analysis 11.

**GUI:** As analysis 11, with Preferences, Estimate subject-specific error
variances set to Yes.

**Scripted:** `'dynrel', 2, 'ssubjrel', 2` plus `dim1col`.

**Engines:** CmdStan only.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** The model and the posterior draws are identical to analysis 11. What
changes is the read-out: each participant's own residual heterogeneity becomes
their own conditional coefficient instead of being pooled into one surface. If
you want both views, one fit gives you both.

**See:** Rast and Clayson (in press).

### Analysis 14: Dynamic test-retest reliability

**Status:** Supported.

**Use when:** Multi-occasion data, one score, and reliability reported as a
function of a participant characteristic.

**Data:** `id`, `meas`, `time`, and `dim1`; optional `dim2`, `group`, and
`event`.

**GUI:** As analysis 11, with Occasion ID also set.

**Scripted:** Pass `timecol` and `dim1col` with `'dynrel', 2`.

**Engines:** CmdStan, or native HMC.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Five crossed mean components sit on top of the one-facet dynamic
model. Occasion-level components stay weakly identified at two occasions, and
the viewer prints them as diagnostics rather than as findings.

**See:** Rast and Clayson (in press).

### Analysis 27: Subject-level dynamic test-retest reliability

**Status:** Supported.

**Use when:** You want per-participant conditional reliability from the
two-facet dynamic model.

**Data:** As analysis 14.

**GUI:** As analysis 14, with Preferences, Estimate subject-specific error
variances set to Yes.

**Scripted:** Pass `timecol` and `dim1col` with `'dynrel', 2, 'ssubjrel', 2`.

**Engines:** CmdStan only.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Same relationship to analysis 14 that analysis 26 has to analysis
11: one model, one set of draws, a per-participant read-out instead of a pooled
surface.

**See:** Rast and Clayson (in press).

### Analysis 12: Dynamic difference-score reliability

**Status:** Supported.

**Use when:** One session, a two-event difference score, and reliability
reported as a function of a participant characteristic.

**Data:** `id`, `meas`, `dim1`, and `event` with exactly two unique event
types; optional `dim2` and `group`. No occasion column.

**GUI:** Specify Inputs, Dimension 1 (dynamic) set; Preferences, Estimate
internal consistency of difference scores set to Yes: 2-event difference, with
Estimate residual covariance for co-occurring events set to No.

**Scripted:** `'dynrel', 2, 'diffest', 2` plus `dim1col`.

**Engines:** CmdStan, or native HMC.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** The dimension carries event-specific slopes on both the mean and the
per-event residual scale, so a dimension can move one condition's error more
than the other's. Both gamma parameterizations are available here, which is not
true of the subject-level sibling.

**See:** Rast and Clayson (in press), with the difference-score partitioning of
Rocha et al. (2026), Table 6.

### Analysis 13: Subject-level dynamic difference-score reliability

**Status:** Supported.

**Use when:** You want each participant's conditional reliability for a
two-event difference score, from one session.

**Data:** As analysis 12.

**GUI:** As analysis 12, with Preferences, Estimate subject-specific error
variances set to Yes.

**Scripted:** `'dynrel', 2, 'diffest', 2, 'ssubjrel', 2` plus `dim1col`.

**Engines:** CmdStan, or native HMC.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Under gamma this design accepts only the location-scale
parameterization, and an explicit request for the other one is rejected rather
than silently substituted. The Gaussian family has no such restriction.

**See:** Rast and Clayson (in press), with the difference-score partitioning of
Rocha et al. (2026), Table 6.

### Analysis 17: Dynamic difference-score reliability, correlated residuals

**Status:** Supported.

**Use when:** As analysis 12, and the two conditions are recorded on the same
trials so their residuals are coupled.

**Data:** As analysis 12.

**GUI:** As analysis 12, with Estimate residual covariance for co-occurring
events set to Yes.

**Scripted:** `'dynrel', 2, 'diffest', 2, 'diffrescor', 2` plus `dim1col`.

**Engines:** CmdStan, or native HMC.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Gaussian only. The gamma family has no implementation of this
design, and a gamma request is rejected before the fit starts. Unlike the
static difference designs, switching residual covariance on here changes the
analysis number, so the two runs are labeled differently in the output.

**See:** Rast and Clayson (in press), with the difference-score partitioning of
Rocha et al. (2026), Table 6.

### Analysis 19: Subject-level dynamic difference-score reliability, correlated residuals

**Status:** Supported.

**Use when:** As analysis 13, with co-occurring conditions, and one residual
correlation shared across the sample.

**Data:** As analysis 12.

**GUI:** As analysis 13, with Estimate residual covariance for co-occurring
events set to Yes, and Specify Inputs, Residual correlation left on Population
(shared).

**Scripted:** `'dynrel', 2, 'diffest', 2, 'ssubjrel', 2, 'diffrescor', 2` plus
`dim1col`.

**Engines:** CmdStan for either family; native HMC for the Gaussian arm only.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Under gamma this design is fitted by a modular cut workflow: the
margins are fitted first without the coupling, and the coupling is estimated
afterward, so the result is a cut posterior rather than a full Bayesian one.
The exported provenance says so on every table. The Gaussian arm is an ordinary
joint fit.

**See:** Rast and Clayson (in press), with the difference-score partitioning of
Rocha et al. (2026), Table 6.

### Analysis 21: Subject-level dynamic difference-score reliability, per-participant residual correlation

**Status:** Supported.

**Use when:** As analysis 19, but the residual coupling itself is expected to
differ across participants.

**Data:** As analysis 12.

**GUI:** As analysis 19, with Specify Inputs, Residual correlation set to
Per-subject. That popup stays greyed out until difference scores,
subject-specific error variances, and residual covariance are all on.

**Scripted:** `'dynrel', 2, 'diffest', 2, 'ssubjrel', 2, 'diffrescor', 2,
'ssrescor', 2` plus `dim1col`.

**Engines:** CmdStan, or native HMC.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Gaussian only; the gamma family has no implementation. A separate
correlation per participant is the most demanding thing this family asks of the
sampler, so budget warmup accordingly and read the diagnostics before the
estimates.

**See:** Rast and Clayson (in press), with the difference-score partitioning of
Rocha et al. (2026), Table 6.

### Analysis 15: Dynamic test-retest difference-score reliability

**Status:** Supported.

**Use when:** As analysis 12, with an occasion facet.

**Data:** `id`, `meas`, `time`, `dim1`, and `event` with exactly two unique
event types; optional `dim2` and `group`.

**GUI:** As analysis 12, with Occasion ID also set.

**Scripted:** Pass `timecol` and `dim1col` with `'dynrel', 2, 'diffest', 2`.

**Engines:** CmdStan, or native HMC.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Gaussian only; the gamma family has no implementation of the
two-facet dynamic difference designs. Six crossed cross-condition blocks are
estimated on top of the dimension slopes, which makes this one of the slower
fits in the toolbox.

**See:** Rast and Clayson (in press), with the difference-score partitioning of
Rocha et al. (2026), Table 6.

### Analysis 16: Subject-level dynamic test-retest difference-score reliability

**Status:** Supported.

**Use when:** As analysis 13, with an occasion facet.

**Data:** As analysis 15.

**GUI:** As analysis 15, with Preferences, Estimate subject-specific error
variances set to Yes.

**Scripted:** Pass `timecol` and `dim1col` with `'dynrel', 2, 'diffest', 2,
'ssubjrel', 2`.

**Engines:** CmdStan, or native HMC.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Gaussian only. The dynamic branch is the only route by which a
subject-level difference score reaches an occasion facet (this analysis and
its concurrent siblings, analyses 20 and 22); the static equivalent is
rejected before the fit starts.

**See:** Rast and Clayson (in press), with the difference-score partitioning of
Rocha et al. (2026), Table 6.

### Analysis 18: Dynamic test-retest difference-score reliability, correlated residuals

**Status:** Supported.

**Use when:** As analysis 15, with co-occurring conditions.

**Data:** As analysis 15.

**GUI:** As analysis 15, with Estimate residual covariance for co-occurring
events set to Yes.

**Scripted:** Pass `timecol` and `dim1col` with `'dynrel', 2, 'diffest', 2,
'diffrescor', 2`.

**Engines:** CmdStan, or native HMC.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Gaussian only. The bivariate residual is evaluated row by row, since
its scale moves with the dimension, which is the main reason this design costs
more per draw than analysis 15.

**See:** Rast and Clayson (in press), with the difference-score partitioning of
Rocha et al. (2026), Table 6.

### Analysis 20: Subject-level dynamic test-retest difference-score reliability, correlated residuals

**Status:** Supported under the Gaussian family. Development-only under the
gamma family.

**Use when:** As analysis 19, with an occasion facet.

**Data:** As analysis 15.

**GUI:** As analysis 16, with Estimate residual covariance for co-occurring
events set to Yes and Residual correlation left on Population (shared).

**Scripted:** Pass `timecol` and `dim1col` with `'dynrel', 2, 'diffest', 2,
'ssubjrel', 2, 'diffrescor', 2`.

**Engines:** CmdStan for either family; native HMC for the Gaussian arm only.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** The gamma arm carries no user-facing status label anywhere in the
software, so this card is the disclosure: its two-facet configuration has no
empirical anchor, nothing outside the toolbox has been reproduced with it, and
it should not be reported from. Like analysis 19 under gamma, it is a cut
posterior. The Gaussian arm is ordinary and supported.

**See:** Rast and Clayson (in press), with the difference-score partitioning of
Rocha et al. (2026), Table 6.

### Analysis 22: Subject-level dynamic test-retest difference-score reliability, per-participant residual correlation

**Status:** Supported.

**Use when:** As analysis 21, with an occasion facet.

**Data:** As analysis 15.

**GUI:** As analysis 20, with Residual correlation set to Per-subject.

**Scripted:** Pass `timecol` and `dim1col` with `'dynrel', 2, 'diffest', 2,
'ssubjrel', 2, 'diffrescor', 2, 'ssrescor', 2`.

**Engines:** CmdStan, or native HMC.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** Gaussian only. This is the most demanding model in the two-event
dynamic family: a per-participant location and scale block, five crossed
cross-condition blocks, dimension slopes, and a per-participant residual
correlation. Treat a first run as a pilot and check the diagnostics before
planning around the numbers.

**See:** Rast and Clayson (in press), with the difference-score partitioning of
Rocha et al. (2026), Table 6.

## 16.6 Dynamic difference of differences

These two designs report each participant's own conditional reliability for a
four-event contrast. They are person-specific by definition, so subject-level
error variances must be on, and they are non-concurrent by definition, so
residual covariance cannot be turned on.

### Analysis 28: Dynamic difference-of-differences reliability

**Status:** Supported.

**Use when:** One session, a four-event contrast, and you want each
participant's own conditional reliability for it.

**Data:** `id`, `meas`, `dim1`, and `event` with exactly four unique event
types, assigned to ERP_1 through ERP_4; optional `dim2` and `group`. No
occasion column. Every participant needs at least one trial in all four cells.

**GUI:** Specify Inputs, Dimension 1 (dynamic) set, Event Type set with exactly
four events selected, and the DoD Contrast... assignment saved. Preferences:
difference scores Yes: 4-event DoD, subject-specific error variances Yes,
residual covariance No.

**Scripted:** `'dynrel', 2, 'diffest', 3, 'ssubjrel', 2` plus `dim1col` and
`dodmap`.

**Engines:** CmdStan only, in both likelihood families.

**Opens in:** Dynamic Reliability: View Results, not the static DoD viewer.

**Notes:** All six pairwise residual covariances are fixed at zero by the
estimand, so this is a full Bayesian posterior rather than a cut one. The
per-participant coverage check runs before the fit: a participant missing a
cell entirely would make their coefficient undefined, and finding that out
after a multi-hour run is worse than being stopped at the start. Expect a large
saved result, since the posterior carries four margins and a per-participant
block.

**See:** Toolbox convention, built on Rast and Clayson (in press) and the Rocha
et al. (2026) Table 6 difference-score partitioning.

### Analysis 29: Dynamic test-retest difference-of-differences reliability

**Status:** Supported.

**Use when:** As analysis 28, with an occasion facet.

**Data:** As analysis 28, plus `time`.

**GUI:** As analysis 28, with Occasion ID also set.

**Scripted:** Pass `timecol` and `dim1col` with `'dynrel', 2, 'diffest', 3,
'ssubjrel', 2` and `dodmap`.

**Engines:** CmdStan only, in both likelihood families.

**Opens in:** Dynamic Reliability: View Results.

**Notes:** The largest and slowest routed design. Occasion standard deviations
are prior-dominated at two occasions and stay printed as diagnostics. One
observed run of this design took about half an hour and wrote a result file of
roughly 300 MB, and loading a file that size back into the viewer takes
noticeably longer than the fit's last iteration might lead you to expect.

**See:** Toolbox convention, built on Rast and Clayson (in press) and the Rocha
et al. (2026) Table 6 difference-score partitioning.

## 16.7 Reliability from splits

Each row is the mean of several items rather than one trial, and a weight
column says how many items went into it. Declaring the splits parallel reuses
the ordinary machinery; declaring them nonparallel fits a weighted
observed-design model with no projection.

### Parallel splits (no analysis number)

**Status:** Supported.

**Use when:** Your rows are split means and every split holds the same number
of items, as with an even split of the trials into halves or quarters.

**Data:** `id`, `meas` (the split mean), and a required items-per-split
`weight` column; optional `group`, `event`, and `time`.

**GUI:** Specify Inputs, Items per split (n_i) set to the weight column, and
Split type set to Parallel (equal n_i). Everything else as for the analysis it
reuses. Preferences must keep difference scores No and subject-specific error
variances No.

**Scripted:** Pass `weightcol` with `'splits', 2`.

**Engines:** Whatever the reused analysis offers, since the model is unchanged.

**Opens in:** View Results Setup, the single-session or test-retest viewer,
with trial relabeled as split throughout the controls and tables.

**Notes:** No analysis number of its own. A single-occasion run reuses analyses
1 to 4 and a multi-occasion run reuses analysis 5, with the split mean as the
score and the split count in place of the trial count. Unequal counts are
allowed but warned about, on the grounds that the nonparallel model probably
fits your data better. Gamma is permanently unavailable for any split design:
the family models single-trial dispersion, and a split mean is not a single
trial.

**See:** Rocha et al. (2026), Tables 4 and 5 for the split framing; the
coefficients computed are the Table 2 and Table 3 forms with the split count in
place of the trial count.

### Analysis 23: Nonparallel splits, one occasion

**Status:** Supported.

**Use when:** Your rows are split means from one session and the splits hold
different numbers of items.

**Data:** `id`, `meas`, and `weight`, with genuine variation in `weight` across
splits. No group, event, or occasion column.

**GUI:** Specify Inputs, Items per split (n_i) set and Split type set to
Nonparallel (unequal n_i), with Group ID, Event Type, Occasion ID, and
Dimension 1 (dynamic) all left on none.

**Scripted:** Pass `weightcol` with `'splits', 3`.

**Engines:** CmdStan, or native fitlme.

**Opens in:** Data-Splits Reliability: View Results.

**Notes:** Equal item counts are rejected outright, because the per-item
residual and the person-by-split interaction are identified only through the
variation in the counts. Low variation is allowed with a warning, and the
warning is worth heeding. Dependability is exact for this design;
generalizability is a downward-biased lower bound (and its standard error a
corresponding upper bound), because a split mean has already lost the
item-level partition the relative coefficient needs. Every relative coefficient
in the output carries a flag saying so. There is no D-study projection and no
trial-count cutoff. A group or event column is rejected rather than pooled.

**See:** Rocha et al. (2026), Table 4.

### Analysis 24: Nonparallel splits, test-retest

**Status:** Supported.

**Use when:** As analysis 23, with two or more occasions.

**Data:** `id`, `meas`, `weight`, and `time`; optional `group` and `event`.

**GUI:** As analysis 23, with Occasion ID also set. Group ID and Event Type are
allowed here.

**Scripted:** Pass `weightcol` and `timecol` with `'splits', 3`.

**Engines:** CmdStan, or native fitlme.

**Opens in:** Data-Splits Reliability: View Results.

**Notes:** Unlike analysis 23, this design is group and event aware and fits
each stratum separately. It reports the three test-retest coefficients at the
observed design only, with no projection and no cutoff. The same exactness
asymmetry applies: dependability is exact, generalizability is a lower bound,
and the flag travels with every relative coefficient.

**See:** Rocha et al. (2026), Table 5.
