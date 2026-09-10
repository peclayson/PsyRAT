# 4. Generalizability theory for ERP scores

An ERP score is an average. A participant's ERN amplitude is the mean of some
number of single-trial values, each of which carries the participant's standing
on the construct plus whatever else was happening in that epoch. Classical test
theory splits such a score into true score and error, which forces every source
of inconsistency into one undifferentiated term and makes the resulting
coefficient depend on which source the chosen method happens to capture.
Generalizability theory replaces that single error term with a set of estimated
variance components, one for each source the design can identify, and then
builds reliability coefficients out of them (Rocha et al., 2026). PsyRAT
implements that framework. This chapter defines the quantities the toolbox
estimates, states the formulas it evaluates, and maps each symbol onto the
column the toolbox prints.

Readers who want the shortest path to a working analysis should read
[Chapter 2](02_quick-start.md) first and return here when the output needs
interpreting. Readers deciding which analysis their design supports should read
this chapter and then [Chapter 5](05_choosing-an-analysis.md).

## 4.1 Persons crossed with trials

Consider a researcher who records ERN from 80 participants in a flanker task and
scores mean amplitude on every error trial. Each participant contributes many
trials, and every participant is measured on the same kind of trial, so persons
and trials are crossed. A facet is a characteristic of the measurement situation
that varies and that could be sampled differently: here, trials. The object of
measurement is the participant, and it is not a facet. The universe score is the
participant's expected score across all conditions in the universe of
generalization, which is the generalizability-theory analog of the true score
in classical test theory (Rocha et al., 2026, Table 1, p. 3).

A single-trial score decomposes into a grand mean and three deviations:

$$
X_{pi} \;=\; \mu \;+\; (\mu_p - \mu) \;+\; (\mu_i - \mu) \;+\; (X_{pi} - \mu_p - \mu_i + \mu)
\qquad (1)
$$

In words, the score of person $p$ on trial $i$ equals the mean across persons and
trials, plus a person effect, plus a trial effect, plus everything left over.
The leftover term is the person-by-trial interaction confounded with residual
error, because each person is observed on each trial only once, and the two
cannot be separated.

The three deviation terms are uncorrelated in the random-effects model that
generalizability theory assumes, so the variance of a single-trial score
partitions into three components:

$$
\sigma^2_X \;=\; \sigma^2_p \;+\; \sigma^2_i \;+\; \sigma^2_{pi,e}
\qquad (2)
$$

That is, the total variance of single-trial scores is the sum of variance
between persons, variance between trials, and person-by-trial residual variance.
The first component is the signal for individual-differences research. The
second reflects mean differences among trials, which shift every participant's
score in the same direction and therefore do not reorder participants. The third
is what averaging over trials reduces.

A facet observed at only one level is hidden: its variance component is not
identifiable, and the reliability that results is conditional on the single level
observed (Rocha et al., 2026, Table 1, p. 3). A single-session study hides the
occasion facet, so its coefficients confound stable person differences with
whatever was specific to that session. This is not a defect in the analysis, but
it does bound the claim the analysis supports.

*Where this is derived: Rocha et al. (2026), Tables 1 and 2 (pp. 3, 7); Baldwin,
Larson, and Clayson (2015), Equations 1 and 2 (p. 792).*

## 4.2 Adding occasions

If the same participants return for a second session, occasion becomes an
estimable facet and the design becomes persons crossed with trials crossed with
occasions. The partition then carries seven components:

$$
\sigma^2_X \;=\; \sigma^2_p + \sigma^2_i + \sigma^2_o + \sigma^2_{pi} + \sigma^2_{po} + \sigma^2_{io} + \sigma^2_{pio,e}
\qquad (3)
$$

The three main effects are persons, trials, and occasions. The three two-way
interactions are person-by-trial, person-by-occasion, and trial-by-occasion. The
final term is the three-way interaction confounded with residual error.

Two of these components carry most of the interpretive weight in retest studies.
The person-by-trial component, $\sigma^2_{pi}$, is specific-factor error: a
systematic source that is not the construct and that reflects idiosyncratic
responding across trials. The person-by-occasion component, $\sigma^2_{po}$, is
transient error: variation unique to the individual that changes across sessions,
often reflecting state fluctuations unrelated to the construct (Rocha et al.,
2026, pp. 5 to 6). Clayson, Carbine, et al. (2021) report an N2 dataset in
which $\sigma^2_{po}$ was the largest error component apart from the residual,
roughly 40% the size of $\sigma^2_p$, which placed a ceiling on retest
reliability that no number of additional trials could raise.

*Where this is derived: Rocha et al. (2026), Table 3 (p. 8); Clayson, Carbine,
et al. (2021), pp. 178 to 183.*

## 4.3 Reading the variance components in PsyRAT output

Selecting Table: Sources of Variance on the View Results Setup screen produces a
table of standard deviations, standard errors of measurement, and intraclass
correlations, each with a point estimate and a credible interval at the width
selected in the viewer (95% by default). The toolbox reports these components on
the standard-deviation scale rather than the variance scale, because standard
deviations are in microvolts and can be compared directly against the ERP scores
themselves. Every formula in this chapter uses variances, so each printed value
is squared before it enters a coefficient.

The mapping between the symbols above and the printed columns is the single most
useful thing to memorize about the output. For a one-facet analysis:

| Symbol | Meaning | Printed column |
|---|---|---|
| $\sigma_p$ | between-person standard deviation | Between Std Dev |
| $\sigma_i$ | trial main effect | Trial Std Dev |
| $\sigma_{pi,e}$ | person-by-trial residual | Within Std Dev |
| $\sqrt{\hat\sigma^2(\Delta)}$ or $\sqrt{\hat\sigma^2(\delta)}$ | standard error of measurement | SEM (dependability) or SEM (generalizability) |
| $\hat\rho_{\text{ICC}}$ | intraclass correlation | ICC (dependability) or ICC (generalizability) |

Two features of that table are easy to misread. The Within Std Dev column holds
the residual $\sigma_{pi,e}$ alone, not a combined within-person term, so the
trial main effect is reported separately in its own column rather than folded
into it. The SEM and ICC headers carry a parenthetical decision type because
both quantities change value with the G-Theory Coefficient setting, and an
archived export would otherwise be impossible to reconstruct. If the analysis
uses data splits rather than single trials, the Trial Std Dev header reads Split
Std Dev, and the analysis unit is the split throughout.

Test-retest analyses print a five-column version of the same table, in which
Between Std Dev is $\sigma_p$ and Within Std Dev is the three-way residual
$\sigma_{pio,e}$. All seven components are listed in the window reached through
the View All Standard Deviations button, two of which repeat the five-column
table:

| Symbol | Printed column |
|---|---|
| $\sigma_p$ | Between Person |
| $\sigma_o$ | Between Session |
| $\sigma_i$ | Between Trial |
| $\sigma_{po}$ | Person x Session |
| $\sigma_{pi}$ | Person x Trial |
| $\sigma_{io}$ | Session x Trial |
| $\sigma_{pio,e}$ | Within Person (Error) |

![Sources of Variance table for a one-facet analysis, showing point estimates and credible intervals for each component.](images/tbl_variance.png)

Comparing the components against one another is the first diagnostic step in any
generalizability study, because their relative sizes indicate which sources of
variance dominate the reliability estimate (Clayson, Carbine, et al., 2021,
p. 180). Which redesign would help follows from where the dominant component
sits in Equation 5, Equation 6, or Equation 9. A large $\sigma_{pi,e}$ relative
to $\sigma_p$ means more trials will help. A large $\sigma_{po}$ relative to
$\sigma_p$ means more trials will not.

*Where this is derived: Rocha et al. (2026), Tables 2 and 3 (pp. 7 to 8);
Clayson, Carbine, et al. (2021), p. 180.*

## 4.4 Two coefficients, two decisions

Reliability in generalizability theory is a ratio of universe-score variance to
universe-score variance plus error variance, and the two coefficients the
toolbox reports differ only in what counts as error. For a one-facet design with
scores averaged over $n'_i$ trials:

$$
\hat\sigma^2(\delta) \;=\; \frac{\sigma^2_{pi,e}}{n'_i}
\qquad\qquad
\hat\sigma^2(\Delta) \;=\; \frac{\sigma^2_{pi,e}}{n'_i} + \frac{\sigma^2_i}{n'_i}
\qquad (4)
$$

Relative error variance, $\hat\sigma^2(\delta)$, contains only the component that
varies idiosyncratically across persons. Absolute error variance,
$\hat\sigma^2(\Delta)$, adds the trial main effect, because a shift common to all
participants still moves each participant's score away from its universe value
even though it leaves the ordering of participants intact.

The generalizability coefficient uses relative error:

$$
E\rho^2 \;=\; \frac{\sigma^2_p}{\sigma^2_p + \dfrac{\sigma^2_{pi,e}}{n'_i}}
\qquad (5)
$$

The dependability coefficient uses absolute error:

$$
\phi \;=\; \frac{\sigma^2_p}{\sigma^2_p + \dfrac{\sigma^2_{pi,e}}{n'_i} + \dfrac{\sigma^2_i}{n'_i}}
\qquad (6)
$$

In words, generalizability is the proportion of variance in an averaged score
that reflects stable person differences once sources that reorder participants
are treated as error, and dependability is the same proportion once every
measured source is treated as error. For a single-score design, dependability
is never larger than generalizability, and the two are equal when
$\sigma^2_i = 0$, which happens when trials do not differ systematically in
their means. (For difference scores, where estimated trial-effect covariances
enter the absolute error, the ordering can reverse in unusual configurations;
Section 4.7 returns to this.)

A worked contrast makes the size of the gap concrete. The manual's tutorial
dataset is simulated, so its generating values are known: for GroupA
incongruent-error trials (the inc_err cell), $\sigma_p = 2.0$,
$\sigma_{pi,e} = 5.5$, and $\sigma_i = 1.0$ microvolts, giving
$\sigma^2_p = 4.00$, $\sigma^2_{pi,e} = 30.25$, and $\sigma^2_i = 1.00$.
Averaged over 17 error trials,

$$
E\rho^2 = \frac{4.00}{4.00 + 30.25/17} = \frac{4.00}{5.779} = 0.692
$$

The trial main effect enters absolute error only, so dependability sits below
that: the absolute error is $(30.25 + 1.00)/17 = 1.838$, and
$\phi = 4.00/5.838 = 0.685$. With this dataset's modest trial main effect the
gap is small. Now suppose a researcher analyzed real data in which the trial
main effect was estimated at $\sigma_i = 2.0$ microvolts, so
$\sigma^2_i = 4.00$. Relative error is unchanged at 1.779, but absolute error
becomes $(30.25 + 4.00)/17 = 2.015$, and

$$
\phi = \frac{4.00}{4.00 + 2.015} = 0.665
$$

Generalizability of 0.692 and dependability of 0.665 describe the same scores.
They differ because they quantify different error definitions, and the choice
between them is a research decision, not a statistical one. If the inference concerns the
relative standing of participants, as it does for correlations, group
comparisons, and most individual-differences work, the generalizability
coefficient is the relevant one. If the inference concerns a participant's
absolute level, as it does for cut scores, clinical thresholds, and any claim
about the magnitude of a score, the dependability coefficient is the relevant
one. The toolbox exposes this as the G-Theory Coefficient control on the View
Results Setup screen, with the options Dependability and Generalizability. Cut
scores are always computed with absolute error regardless of that setting,
because a threshold decision is by definition an absolute one.

*Where this is derived: Rocha et al. (2026), Table 2 (p. 7) and their
Section 3.3.1 (p. 5).*

## 4.5 Standard error of measurement

The two error variances have matching standard errors of measurement, which are
simply their square roots:

$$
\text{SEM}_{\delta} = \sqrt{\frac{\sigma^2_{pi,e}}{n'_i}}
\qquad\qquad
\text{SEM}_{\Delta} = \sqrt{\frac{\sigma^2_{pi,e}}{n'_i} + \frac{\sigma^2_i}{n'_i}}
\qquad (7)
$$

The standard error of measurement is expressed in the units of the score, which
for ERP amplitude means microvolts, and it quantifies the typical distance
between an observed averaged score and the universe score behind it.
Continuing the example above, the relative standard error at 17 trials is
$\sqrt{1.779} = 1.33$ microvolts and the absolute standard error is
$\sqrt{1.838} = 1.36$ microvolts. The between-person standard deviation for
that cell is 2.0 microvolts, so an individual averaged score carries roughly
1.4 microvolts of measurement error against 2.0 microvolts of true
between-person spread, which is the comparison the standard error of
measurement exists to make available, and which is why the dependability there
is only .685.

The toolbox evaluates Equation 7 at the trial count actually available rather
than at a single trial. That count is the mean number of retained trials across
participants who met the trial cutoff, or the median if the median was selected
in the viewer preferences. The reported standard error is therefore the standard
error of the scores in hand, and it moves with the G-Theory Coefficient setting
in step with the coefficient printed beside it.

*Where this is derived: Rocha et al. (2026), Table 2 (p. 7).*

## 4.6 The decision study: projecting reliability onto a design

Generalizability theory separates two phases. A generalizability study estimates
the variance components. A decision study uses those components to project
reliability for a measurement procedure that was not necessarily the one run,
most often by varying the number of trials retained for averaging (Baldwin et
al., 2015, p. 793). The variance components do not change between the two
phases. Only $n'_i$ changes.

Because $n'_i$ appears only in the denominator of Equations 5 and 6, the
projection from those two equations is monotonic: more trials yield a higher
coefficient, with diminishing returns. (The coefficient of stability of
Equation 9 is the exception: when the person-by-trial variance is large
relative to the person-by-occasion terms, adding trials can lower it.) For the
tutorial values above, dependability rises from 0.390 at 5 trials to 0.561 at 10, 0.685 at 17, 0.793 at 30, and 0.865 at 50.
Inverting Equation 6 gives the count a researcher needs for a target coefficient
$\rho$:

$$
n'_i \;=\; \frac{\rho}{1-\rho} \cdot \frac{\sigma^2_{pi,e} + \sigma^2_i}{\sigma^2_p}
\qquad (8)
$$

In words, the required number of trials is the target coefficient expressed as
odds, multiplied by the ratio of error variance to between-person variance. With
$\sigma^2_i = 0$, reaching a dependability of 0.80 requires $4 \times (30.25/4.00) = 30.25$ trials, so 31
trials. Adding the hypothetical trial main effect of $\sigma^2_i = 4.00$ raises
the requirement to 34.25, so 35 trials. That is the practical cost of a trial
main effect: four additional error trials per participant, in a paradigm where
error trials are the scarce resource.

![Dependability as a function of the number of trials retained for averaging, one panel per event and one curve per group, with the reliability cutoff drawn as a dotted horizontal line.](images/fig_rel_vs_trials.png)

The plot produced by the toolbox draws that curve directly, one panel per
event, with the number of observations on the horizontal axis and the
coefficient on the vertical axis, bounded above at 1. The dotted horizontal
line is the reliability cutoff, and the point where a curve crosses it is the
trial count that design requires.

Choosing that cutoff is a judgment the software cannot make. Clayson and Miller
(2017b, p. 62) recommend a threshold of at least 0.70 for exploratory ERP
research and 0.80 in most contexts. The reasoning behind those two values is that
0.70 buys enough precision to decide whether a measure is worth developing
further, while 0.80 is what group-difference research needs before unreliability
begins attenuating the effects under test. The toolbox default is 0.80. The
exception runs upward, not downward: when a decision is made about an individual participant's score, as in clinical classification,
0.90 is a minimum and 0.95 is preferred. The illegitimate reason to lower the
threshold is that the data did not reach it. A threshold chosen after the
coefficient is known is not a threshold, and the honest response to a coefficient
below the stated cutoff is to report it and to state what the analysis therefore
cannot support.

For a two-facet test-retest design the projection has two arguments, $n'_i$ and
$n'_o$. The dependability coefficient of stability, which is the retest analog,
takes the form

$$
\phi_{CS} \;=\; \frac{\sigma^2_p + \dfrac{\sigma^2_{pi}}{n'_i}}{\sigma^2_p + \dfrac{\sigma^2_{pi}}{n'_i} + \dfrac{\sigma^2_{po}}{n'_o} + \dfrac{\sigma^2_{pio,e}}{n'_i n'_o} + \dfrac{\sigma^2_o}{n'_o} + \dfrac{\sigma^2_{io}}{n'_i n'_o}}
\qquad (9)
$$

The person-by-trial component sits in the numerator here, because a coefficient
of stability treats consistency across trials as part of the universe score and
asks only whether scores hold up across occasions. The coefficient of
equivalence exchanges those roles, placing $\sigma^2_{po}$ in the numerator and
$\sigma^2_{pi}$ in the error term. The coefficient of trial equivalence and
stability keeps $\sigma^2_p$ alone in the numerator and treats both interactions
as error, which makes it the most demanding of the three.

By the toolbox's convention, every test-retest coefficient defaults to
$n'_o = 1$. The estimand is then the reliability of a score collected at one
occasion and generalized over the occasion facet, which is the score a researcher
actually has in hand after a single visit. Setting $n'_o = k$ instead reports the
reliability of a composite averaged over $k$ occasions, a different and generally
higher quantity, and the toolbox labels the two distinctly wherever they appear
so that a coefficient is never paired with the other estimand's standard error.
For the coefficient of stability the $n'_o = 1$ default matches the practice
described by Clayson, Carbine, et al. (2021, p. 183). Selecting the multi-occasion
composite is an explicit choice, described in [Chapter 9](09_viewing-results.md).

*Where this is derived: Rocha et al. (2026), Table 3 (p. 8); Baldwin et al.
(2015), Equation 4 (p. 793); Clayson and Miller (2017b), p. 62; Clayson,
Carbine, et al. (2021), p. 183.*

## 4.7 Difference scores

Many ERP measures are contrasts. ERN is frequently scored as error minus correct,
and lateralized indices are scored as one hemisphere minus the other. The
reliability of such a score is not the reliability of either constituent, and it
is usually lower than both.

Write $X$ and $Y$ for the two constituent measures. The universe-score variance
of the difference and its two error terms are

$$
\sigma^2_{p(X-Y)} = \sigma^2_{X(p)} + \sigma^2_{Y(p)} - 2\sigma_{XY(p)}
\qquad (10)
$$

$$
\hat\sigma^2(\delta) = \frac{\sigma^2_{X(pi,e)}}{n'_{iX}} + \frac{\sigma^2_{Y(pi,e)}}{n'_{iY}} - \frac{2\sigma_{XY(pi,e)}}{\ddot{n}'_i}
\qquad\qquad
\hat\sigma^2(\Delta) = \hat\sigma^2(\delta) + \frac{\sigma^2_{X(i)}}{n'_{iX}} + \frac{\sigma^2_{Y(i)}}{n'_{iY}} - \frac{2\sigma_{XY(i)}}{\ddot{n}'_i}
\qquad (11)
$$

Each coefficient places Equation 10 over Equation 10 plus one of the error terms
in Equation 11: the relative error term for generalizability, the absolute one
for dependability. In words, the numerator is the between-person variance the two
conditions do not share, and the denominator adds the measurement error of both
conditions minus whatever error they share. Covariance terms are divided by the
harmonic mean of the two trial counts, written $\ddot{n}'_i$, the convention
Rocha et al. (2026, Table 6) adopt; since $2/\ddot{n}'_i = 1/n'_{iX} +
1/n'_{iY}$, this is exactly the scaling that makes the difference-score error
variance vanish when the two conditions' components are identical.

The subtraction in the numerator is why difference scores are typically less
dependable than their constituents. Whatever between-person variance the two
conditions share is removed from the signal, while the error variances add rather
than cancel. The tutorial dataset makes the arithmetic explicit. Error-trial and
correct-trial person effects were generated with a correlation of $r = 0.75$,
with $\sigma_p = 2.0$ for error trials and $\sigma_p = 1.7$ for correct trials,
so $\sigma_{XY(p)} = 0.75 \times 2.0 \times 1.7 = 2.55$. The universe-score
variance of the difference is $4.00 + 2.89 - 2(2.55) = 1.79$, which is smaller
than either constituent. Error trials are noisier and scarcer than correct
trials, so at 17 error trials and 32 correct trials the relative error variance is
$30.25/17 + 25.00/32 = 2.561$, and

$$
E\rho^2_{X-Y} = \frac{1.79}{1.79 + 2.561} = 0.411
$$

The generating trial main effects ($\sigma_i = 1.0$ per cell) are drawn
independently across conditions, so their cross-condition covariance is zero
and dependability differs from generalizability only by the two
$\sigma^2_i/n'_i$ terms: absolute error is $31.25/17 + 26.00/32 = 2.651$ and
$\phi_{X-Y} = 1.79/4.441 = 0.403$. The constituent generalizability
coefficients at those same trial counts are 0.692 for error trials and 0.787
for correct trials. The difference score is markedly less reliable than
either, and the correlation between conditions is the reason. Had
the person effects been uncorrelated, the numerator would have been 6.89 instead
of 1.79, and the same error variance would have produced a coefficient of 0.729.

Two further points govern how the toolbox handles these designs. First, the
residual covariance $\sigma_{XY(pi,e)}$ is estimable only when the two measures
are recorded concurrently, as in an ipsilateral-versus-contralateral contrast.
When they come from separate trials, as error and correct trials do, the residual
covariance is set to zero (Rocha et al., 2026, p. 15). Second, a low
difference-score coefficient does not by itself indicate a flawed measure. When
true between-person differences in the contrast are genuinely small, or when
participants differ but all to the same degree, the numerator is small or zero
and the coefficient is correspondingly low even for a perfect instrument; the
toolbox reports a zero numerator as an exact zero with a diagnostic rather
than treating it as an error (Rocha et al., 2026, Section 5.1, p. 15).
Correlated measurement error runs the other way and can inflate the
coefficient. One further property of these estimates deserves stating: the
difference-score error terms are differences of estimated components, so a
legitimate posterior draw can push a coefficient outside the zero-to-one range,
and the toolbox reports the value as evaluated together with a diagnostic
rather than silently repairing it.

The difference-score row of the All Standard Deviation Components window,
reached from the View All Standard Deviations button, carries an ICC computed
at a single observation, labeled ICC (n'=1). It is a derived quantity for the
difference score rather than an entry in the published difference-score
formulas; the label records the n'=1 convention, and the derivation is
documented here and in [Chapter 12](12_tutorial-difference-scores.md).

*Where this is derived: Rocha et al. (2026), Table 6 and their Section 5.1
(pp. 13, 15).*

## 4.8 The intraclass correlation

The intraclass correlation is the coefficient evaluated at a single
observation. It addresses a different question from a dependability
coefficient: not how reliable an averaged score is, but how much of the
variance in one trial is between-person variance.

$$
\rho_{\text{ICC},\delta} = \frac{\sigma^2_p}{\sigma^2_p + \sigma^2_{pi,e}}
\qquad\qquad
\rho_{\text{ICC},\Delta} = \frac{\sigma^2_p}{\sigma^2_p + \sigma^2_{pi,e} + \sigma^2_i}
\qquad (12)
$$

The absolute version has the plainest reading. An absolute ICC of 0.40
indicates that 40% of the variance in single-trial scores is stable
between-person variance and the remaining 60% counts as error for an absolute
decision. The relative version removes the trial main effect from the
denominator, so its percentage is the share of variance among only those sources
that affect the relative standing of participants. The two are equal when
$\sigma^2_i = 0$.

Because the ICC ignores the number of trials, it is a trial-independent index of
how much person-level signal each observation carries (Rocha et al., 2026, p. 16).
That property is what makes it a useful companion to the dependability
coefficient rather than a substitute for it. The tutorial values illustrate the
gap. For GroupA error trials the relative ICC is $4.00/34.25 = 0.117$ and the
absolute ICC, which keeps the trial main effect in the denominator, is
$4.00/35.25 = 0.113$: fewer than 12% of the variance in a single error trial is
between-person variance either way. The dependability of the same scores
averaged over 17 trials is 0.685. Averaging is
what converts a weak per-trial signal into a usable score, and a design that
yields more usable trials can produce more dependable scores from an
equally noisy signal.

![ICC Estimates plot, one point estimate and credible interval per event and group.](images/fig_icc_forest.png)

The toolbox plots these estimates as a forest plot, with the events listed down
the left, the intraclass correlation running along the bottom, and one point
estimate with its credible interval for each group and event. Reading it means
comparing the horizontal positions of the points, which indicate per-trial
signal, and the widths of the intervals, which indicate how well that signal
was estimated. Overlap between two intervals is not by itself evidence that two events carry the
same per-trial signal, since that claim concerns the posterior distribution of
their difference and not the two marginal intervals. Events with similar ICCs can
still differ substantially in dependability whenever their trial counts differ,
which is the point of reporting both.

*Where this is derived: Rocha et al. (2026), Table 2 and their Section 8
(pp. 7, 16).*

## 4.9 Why the estimates are Bayesian

Every formula in this chapter takes variance components as input. How those
components are estimated determines whether the formulas can be evaluated at all.
PsyRAT estimates them with Bayesian multilevel models fitted in CmdStan, for
reasons that bear directly on the structure of ERP data.

| Problem with the usual approach | What the Bayesian estimates provide |
|---|---|
| Random-effects ANOVA on unbalanced data can return negative variance components, which makes the coefficients uninterpretable | Variance components are constrained to be positive by the model's parameter space, so a negative component estimate cannot occur (derived difference-score error terms can still leave the usual range; Section 4.7) |
| Restricted maximum likelihood sets negative estimates to zero during iteration, which blocks standard errors and intervals | Every component carries a full posterior distribution, so credible intervals accompany the coefficients, standard errors, and intraclass correlations in the main output tables (a few summary windows print point estimates only) |
| Conventional inference assumes homoscedastic residuals, and departures distort variance-component estimates in the small, unbalanced samples common to psychophysiology | In the subject-level and dynamic designs, residual variance is modeled explicitly rather than assumed constant |
| Reliability is reported for the group as a whole | The same posterior draws support subject-level coefficients, which the toolbox reports for the designs that estimate them |

Trial counts in ERP research are unbalanced by construction, because artifact
rejection removes different numbers of epochs from different participants, and
error trials are scarce. That is the condition under which the first two problems
actually bite (Rocha et al., 2026, Section 7, pp. 15 to 16).

### Priors

The toolbox ships weakly informative default priors, and those defaults are
calibrated on the ERP microvolt scale. They are not scale-free. A measure
recorded on a very different scale, such as standardized scores or a modality
whose values are orders of magnitude larger or smaller, can make the same priors
informative, which would let the prior rather than the data drive the estimates.
Researchers analyzing scores on another scale should set the priors deliberately
rather than accepting the defaults. [Chapter 15](15_advanced-designs.md) covers
the priors the toolbox exposes and how to override them.

### Sampling controls

Four processing preferences govern the fit. Each is set on the Specify Processing
Preferences screen, described in [Chapter 7](07_processing.md).

| Control | What it does | Default |
|---|---|---|
| Number of Chains | independent sampling runs started from different points; agreement among them is the evidence of convergence | 4 |
| Warmup iterations | iterations used to tune the sampler and discarded before estimation | 5000 |
| Sampling iterations | iterations retained per chain to form the posterior | 5000 |
| Random seed | fixes the random number stream so a rerun reproduces the same draws | 12345 |

The seed is the reproducibility control. Holding the seed, the settings, and the
input file fixed reproduces a run's draws exactly on the same machine and CmdStan
version, with within-chain parallelization off ([Chapter 7](07_processing.md)),
which is what makes an archived analysis checkable. Warmup and sampling
iterations are the controls to raise when a model does not converge.

### Convergence

The toolbox checks convergence automatically after every fit and applies two
rules to every monitored parameter. A run is marked as not converged if any
parameter has a potential scale reduction factor (R-hat) of 1.1 or higher, or
if any parameter has an effective sample size below ten times the number of chains,
which is 40 at the default of four chains. Passing both rules is a necessary
condition for trusting the estimates, not a sufficient one.

When a fit fails either rule, the toolbox offers to rerun with the warmup and
sampling iterations doubled. Accepting is usually the right response, because
non-convergence most often reflects too short a run rather than a
misspecified model. The estimates from a non-converged fit should not be
reported.

Divergent transitions are counted separately and reported, but they do not
trigger a rerun. A nonzero count is a warning that the variance components may be
biased even when both convergence rules pass, and the response is more warmup or
stronger priors. [Chapter 18](18_troubleshooting.md) covers both cases.

### Credible intervals are not confidence intervals

The intervals PsyRAT prints are credible intervals, and the distinction from
confidence intervals is not cosmetic. A 95% credible interval is a statement
about the parameter given the data, the prior, and the model: the posterior
probability that the dependability coefficient lies between the reported bounds
is 0.95, up to the Monte Carlo error of a finite sample of draws. A 95%
confidence interval is a statement about the procedure: intervals constructed
this way would contain the fixed true value in 95% of hypothetical replications,
and the particular interval in hand either contains that value or does not. The
everyday reading that researchers apply to a confidence interval is the one that
is actually licensed for a credible interval. Reporting a credible interval in
frequentist language, or the reverse, misstates what was estimated.

Two details of construction belong in a methods section. The toolbox forms each
interval from the empirical quantiles of the posterior draws, taking equal tails,
so a 95% interval runs from the 2.5th to the 97.5th percentile. The point
estimate beside it is the posterior mean of the same draws, not the median or the
mode.

*Where this is derived: Rocha et al. (2026), Section 7 (pp. 15 to 16). The
convergence rules and the priors are the toolbox's own settings.*

## 4.10 Where the framework goes next

The designs in this chapter assume that residual variance is a single quantity
shared by all participants and all conditions. Rast and Clayson (in press)
integrate generalizability theory with multilevel location-scale models, which
relax that assumption and let residual variance vary systematically across
persons and conditions. The result is a reliability function rather than a single
coefficient, and it supports allocating trials to the participants or conditions
where precision is lowest. PsyRAT implements this as its dynamic and conditional
reliability analyses, described in [Chapter 15](15_advanced-designs.md).

For readers who want a broader treatment of generalizability theory outside
psychophysiology, including its relations to classical test theory and structural
equation modeling, Vispoel et al. (2018) is the standard reference. Rocha et al.
(2026) is the source of record for the formulas the toolbox evaluates and is the
preferred citation for a PsyRAT analysis.

[Chapter 5](05_choosing-an-analysis.md) maps these designs onto the analyses the
toolbox offers. [Chapter 17](17_reporting.md) covers what to report.
Definitions of the terms used here are collected in
[Appendix A](appendix-a_glossary.md), and full references in
[Appendix B](appendix-b_references.md).
