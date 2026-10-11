# Appendix A. Glossary

Terms as the toolbox uses them. Where a term has a broader meaning elsewhere in
the psychometric literature, the definition here is the one that governs the
output, the GUI labels, and [Chapter 4](04_gtheory.md).

**chains.** Independent runs of the sampler over the same model and data,
started from different initial values. Agreement between chains is what makes
R-hat informative, so the toolbox requires at least 3 and uses 4 by default.

**coefficient of equivalence.** The test-retest coefficient that generalizes
over trials while holding the occasion fixed. Trial-related terms enter the
error; the occasion main effect does not.

**coefficient of equivalence and stability (CES).** The test-retest coefficient
that generalizes over both trials and occasions at once; only person variance
is in the numerator, and it cannot exceed either of the other two. Selected as
Equivalence + Stability in the viewer.

**coefficient of stability.** The test-retest coefficient that generalizes over
occasions while holding the trial sample fixed. Occasion-related terms enter the
error.

**credible interval.** The interval between two quantiles of the posterior
distribution of a coefficient, at the width the viewer uses: 0.95, or the value
entered in the Credible interval control on the viewers that offer one. Reported
coefficients are posterior means with such an interval around them, except in
the test-retest window Point Estimates of All Standard Deviation Components,
which prints point estimates only. A credible interval is a statement about
the posterior, not a frequentist confidence statement. A run on the native
fitlme engine reports a parametric-bootstrap interval instead; the viewers,
plot axes and export headers label it Bootstrap interval, and it is a
frequentist statement, not a posterior one.

**cut score (criterion-referenced reliability).** A threshold on the measurement
scale, entered as Criterion Cutoff, used to ask how reliably scores are
classified relative to that value. Cut score coefficients are absolute-error
coefficients only; a relative-error version is not defined and is rejected.

**dependability.** The absolute-error reliability coefficient, chosen as
Dependability in a viewer's coefficient popup or as `gcoeff = 1` from a script. Its error term includes the trial main effect divided by
the projected number of trials, so it speaks to how well an observed score
estimates a participant's own universe score on the measurement scale.
Dependability and generalizability are not interchangeable, and neither is a
synonym for reliability.

**difference score.** A score formed by subtracting one condition's score from
another's, selected with Estimate internal consistency of difference scores.
The two constituents are modeled jointly, so their between-person covariance
enters the difference's universe-score variance. The four-event contrast is a
difference of two difference scores.

**difference of differences (DoD).** The four-event contrast
(ERP_1 - ERP_2) - (ERP_3 - ERP_4), selected with Yes: 4-event DoD and mapped on
the Configure DoD Contrast screen. Its residual covariances are fixed at zero
in every family.

**dimension.** A continuous between-person covariate. Selecting a column in
the Dimension 1 (dynamic) popup turns on dynamic (conditional) reliability, and
the Dimension 2 (dynamic) popup adds an optional second covariate; Dimension 2
cannot be used without Dimension 1. Each
covariate must carry one value per participant, repeated across that
participant's rows; the toolbox standardizes it across participants before
estimation, so a dimension value of 0 is the sample mean.

**divergent transition.** A sampler iteration whose trajectory could not be
followed accurately, reported by CmdStan and counted by the toolbox. Divergences
never trigger a rerun and never change the convergence verdict, but a nonzero
count flags possible bias in the posterior. More warmup, a higher Adapt delta,
or stronger priors are the remedies.

**D-study.** The step that projects estimated variance components to a chosen
number of trials (n'_i) and occasions (n'_o), producing the coefficient at
those counts rather than at the numbers the sample happened to contain. Trial-count curves
and trial cutoffs are D-study output.

**dynamic (conditional) reliability.** Reliability reported as a function of one
or two standardized between-person dimensions rather than as one number, from
a model in which the residual scale varies with the dimension (Rast & Clayson,
in press). You select it by choosing a column in the Dimension 1 (dynamic)
popup.

**effective sample size.** The number of independent draws a correlated chain is
worth for a given parameter. The convergence check fails when any monitored
parameter's effective sample size falls below ten times the number of chains.
For CmdStan 2.33 or newer and for the native HMC engine, the effective sample
size checked is the smaller of the bulk and tail values; an older CmdStan
reports only one, and the toolbox records its diagnostics as degraded.

**ERN (error-related negativity).** The negative-going, response-locked ERP
component elicited by errors, and the manual's running example. The tutorial
datasets are simulated ERN-like scores.

**ERP (event-related potential).** An EEG waveform time-locked to an event and
averaged over trials. PsyRAT's reliability analyses take the single-trial
scores that go into such an average; the Import + Score ERP workflow reads the
waveforms to produce those scores.

**error variance (absolute and relative).** The denominator term that separates
the two coefficients. Absolute error is the variance of the difference between
an observed score and the universe score, and it includes the trial main effect
divided by the projected number of trials (n'_i). Relative error concerns only
the ordering of participants, so the trial main effect drops out. For a
single-score design, absolute error is therefore at least as large as relative
error, and dependability is at most generalizability (Section 4.7 notes the
difference-score exception).

**estimation engine.** The fitter that estimates the variance components:
CmdStan (the reference path), native HMC (Hamiltonian Monte Carlo), or native
fitlme (restricted maximum likelihood with a bootstrap). In the Estimation
engine popup, Automatic (recommended) uses CmdStan whenever it is installed and
falls back to a native engine only when it is not (a CmdStan older than the
minimum version stops the run instead of falling back); CmdStan only, Native
HMC, and Native fitlme each force one engine. The engine that ran is recorded
in every export header and in the run-config sidecar.

**event.** A condition, such as correct or error trials, taken from the column
assigned to Event Type. Difference scores are formed across two events; the
four-event contrast is formed across four.

**facet.** A source of measurement error over which a score is generalized. The
trial facet is present in every design; the occasion facet is added by assigning
an Occasion ID column. Participants are the object of measurement rather than a
facet.

**gamma family.** The scaled chi-square (Gamma with a log link) observation
likelihood, selected with Observation likelihood family, for strictly positive,
right-skewed scores such as single-trial time-frequency power. It runs only on
CmdStan and covers a subset of the designs.

**generalizability.** The relative-error reliability coefficient, chosen as
Generalizability in a viewer's coefficient popup or as `gcoeff = 2` from a
script. Its error term excludes the trial main effect, so it speaks to how
well the ordering of participants is reproduced. See dependability.

**G-study.** The step that estimates the variance components of the design from
the observed data. In this toolbox the G-study is the Bayesian fit, and its
output is a set of posterior draws for each component.

**ICC (intraclass correlation coefficient).** The proportion of total variance
attributable to between-person differences at a single observation. It is the
reliability coefficient evaluated at one trial, so it does not depend on a
projected trial count and is not comparable with a coefficient reported at n'
trials.

**location-scale.** A model in which the residual standard deviation (the
scale) is modeled alongside the mean (the location), per participant or as a
function of a dimension. The subject-level and dynamic designs are
location-scale models. Under the gamma family, Location-scale is also one of
the options in the Gamma scale parameterization popup.

**maturity label.** The per-design status printed on every Chapter 16 card:
Supported, Restricted, Experimental, Tractability checked (robustness not),
Development-only, or Not implemented. Section 5.4 defines each one.

**MCSE (Monte Carlo standard error).** The standard error of a posterior
summary that comes from using a finite number of draws. The toolbox stores the
MCSE of each monitored parameter's mean beside its effective sample size. It
does not set an MCSE tolerance for the reported coefficients, which is why
every notice and header reports final-inference readiness as unassessed.

**modular cut posterior.** The result of fitting a model in two stages: the
margins first, without a coupling term, and the coupling afterward,
conditional on the retained margin draws, with nothing fed back. Under the
gamma family, the concurrent subject-level dynamic difference designs
(analyses 19 and 20 in [Chapter 16](16_analysis-reference.md)) are fitted this
way. Their intervals are not full joint-posterior credible intervals, and the
export header says so.

**nonparallel splits.** Splits holding unequal numbers of items, declared as
Nonparallel (unequal n_i). They are fitted with a weighted observed-design
model that is identified only by the variation in n_i and report the observed
design by default, with no projection over item counts. On a single occasion
(analysis 23) generalizability is
reported as a lower bound; on several occasions (analysis 24) it is reported
under the conventional relative partition with no bound asserted.

**occasion.** A level of the occasion facet, taken from the column assigned to
Occasion ID. Its presence makes a design test-retest. Test-retest coefficients
default to a single-occasion score (n'_o = 1); the Multi-occasion composite
(n'_o = k) reports a score averaged over k occasions and is a different
estimand.

**parallel splits.** Splits holding equal numbers of items, declared as Parallel
(equal n_i). They reuse the one-facet and test-retest machinery with the split
mean as the score and the split count in place of the trial count.

**person variance.** The between-person variance component, the variance of
participants' universe scores. It is the numerator of the one-facet
coefficients.

**posterior.** The distribution of a parameter given the data, the priors, and
the model, represented in PsyRAT by the draws retained from every chain.
Reported point estimates are the means of those draws, and credible intervals
are pairs of their quantiles.

**prior.** The distribution assumed for a parameter before the data are seen.
The Gaussian-family defaults are weakly informative on the ERP microvolt scale
and are not scale-free. Priors are set from a script, never from a screen, and
the merged set a run used is recorded in its sidecar.

**reliability.** The umbrella term for the psychometric quality of scores.
Dependability and generalizability are its two precise generalizability-theory
forms in this toolbox, and output that names one of them means that one.

**R-hat (potential scale reduction factor).** A convergence diagnostic
comparing variation between chains with variation within them; values near 1
indicate the chains describe the same distribution. The convergence check fails
when any monitored parameter has R-hat of 1.1 or higher.

**run-config sidecar.** The `<name>_runconfig.json` file written beside every
`.psyrat` result, recording the seed, sampler settings, column mapping, design
flags, priors, and versions; `psyrat_run('config', ...)` replays it against a
data source.

**seed.** The integer that initializes the sampler's random number generator,
12345 by default. The same data, settings, and seed reproduce the same draws on
the same platform and CmdStan version.

**SEM (standard error of measurement).** The square root of the error variance
at the projected number of observations, reported on the measurement scale. The
absolute SEM pairs with dependability and includes the trial main effect; the
relative SEM pairs with generalizability and excludes it.

**sentinel.** A value that marks the absence of something rather than measuring
it. `(none)` in a Specify Inputs popup means that no column is assigned to
that role (no groups, no events, no occasion facet, or no dimension). In the
trial-cutoff table, -1 in the Trial Cutoff column means no cutoff could be
found at the threshold in the evaluated range, and the coefficient cells of
that row read -1 as well. A cutoff beyond the observed trial counts is not a
sentinel: it stays numeric, and its coefficient carries the suffix "(projection
beyond observed data)" (a summary saved by an earlier version may still show
-1 there until its tables are generated again). `(none)` is not a column name,
and -1 is not a count.

**split.** One row of a splits input, holding the mean of several items rather
than a single trial. A splits analysis treats the split as the unit, and output
that would otherwise say "trial" says "split".

**stratum.** One separately modeled cell of the design: a group, an event, or a
group-by-event combination. Variance components and coefficients are estimated
within each stratum, and dynamic-reliability figures draw one panel per
stratum.

**subject-level (subject-specific) error variance.** A residual variance
estimated separately for each participant, rather than one residual shared
across the sample. Turning it on yields a reliability estimate for each
participant, so inclusion decisions can rest on each participant's own
reliability rather than on a trial count compared with the group's cutoff.

**trial cutoff.** The smallest number of trials at which the projected
coefficient reaches the value entered in Reliability Cutoff. A cutoff above the
largest observed trial count is an extrapolation and is flagged as one. A cutoff
of -1 is not a count: it records that no cutoff could be found at that
threshold.

**universe score.** The score a participant would obtain averaged over the
universe of admissible observations. A reliability coefficient asks how well
an observed score estimates it. Its variance, the universe-score variance, is
the numerator of every coefficient in this toolbox, and it is a property of
the score rather than of the construct the score indexes.

**warmup.** The sampler iterations used for adaptation and discarded before the
posterior summaries are computed. Warmup and sampling iterations are set
separately; warmup draws never enter a reported estimate.

**weight (n_i).** The items-per-split column required by a splits analysis,
giving the number of items averaged into each split mean. It must be a positive
whole number. Variation in n_i across splits is what identifies the per-item
residual separately from the person-by-split interaction, so the nonparallel
model rejects all-equal weights.
