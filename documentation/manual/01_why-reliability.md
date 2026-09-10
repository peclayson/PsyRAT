# 1. Why analyze the reliability of ERP scores?

Reliability is a property of scores (the data in hand), not a property of
measures. Everything this toolbox does follows from that one sentence, and most
of this chapter is spent defending it against a habit the event-related
potential (ERP) literature has held for decades. If you are already convinced, [Chapter 2](02_quick-start.md) puts a
number on screen in one sitting.

## 1.1 Reliability is a property of scores, not measures

Most published ERP work cites reliability rather than estimating it. A methods
section reports that internal consistency for error-related negativity (ERN)
scores has been established elsewhere, supplies a citation, and moves on. Clayson and Miller (2017b) named
this the common misunderstanding that reliability is a fundamental property of a
measure. It is not. A reliability estimate describes one particular set of
scores, obtained from one sample, with one task and one processing pipeline.

The evidence against borrowing is not subtle. A preregistered reliability
generalization study pooled 189 internal-consistency estimates of ERN scores
drawn from 68 samples nested in 43 studies (n = 4,499 participants). The pooled
estimate was .68, 95% CI [.63, .72], and about 90% of the variation among those
estimates was not attributable to sampling error (Clayson, 2020). An average
that heterogeneous is not a number anyone can inherit.

Population is one reason. Baldwin, Larson, and Clayson (2015) estimated the
dependability of ERN and error positivity (Pe) scores in 34 participants meeting criteria for major
depression, 29 meeting criteria for an anxiety disorder, and 319 controls.
Within-person variance exceeded between-person variance in every group, and the
controls needed fewer trials than the clinical groups to reach the same
dependability. Under a fixed eight-trial inclusion rule, internal consistency
ran from .52 to .60 in psychopathology samples against .67 to .69 in controls
(Clayson, 2025). Same component, same rule, different scores, different
reliability.

Task is another. ERN scores recorded during flanker, Stroop, and go/no-go tasks
show both shared and unique variance, and their internal consistency differs by
task (Clayson & Miller, 2017b). Then there is the processing pipeline, which is
where the variability is worst. A multiverse analysis computed 3,456 ERN scores
per person across defensible choices about filtering, ocular correction,
reference scheme, baseline window, sensor, and scoring approach. Those choices
changed both the mean amplitudes and the variability of the scores (Clayson,
2024). Reliability moves with all of it. Your pipeline is not the cited study's
pipeline, and there is no reason to expect its estimate to survive the move.

Trial-count rules of thumb are the most common form the borrowing takes. Olvet
and Hajcak (2009) reported that six to eight error trials could support internal
consistency near .50 in undergraduates, and six to eight became a field-wide
inclusion rule (Clayson, 2025). Across the ERN literature as a whole, the
reliability generalization study put the requirement at 9 trials for .70 and 16
trials for .80. After influential estimates were dropped it was 10 and 17
(Clayson, 2020). Neither pair is your number. Both were averaged over samples,
tasks, and pipelines that are not yours, and the second pair was estimated from
the very heterogeneity that makes the first pair unsafe. Computing reliability on the
data in hand replaces a borrowed cutoff with one your own data justify, which is
the entire reason PsyRAT reports a trial cutoff at all.

One kind of reliability also does not stand in for another. Carbine et al.
(2021) recorded four food-related ERP scores from 132 participants
twice, two weeks apart. Coefficients of equivalence (internal consistency)
exceeded .96 for every score. Coefficients of stability (test-retest) in the same
data ranged from .48 to above .65. Scores that separate people cleanly inside one
session need not preserve that ordering across two weeks, and a paper reporting
only the first coefficient has not addressed the second.

Now the costs. Poor reliability caps the correlations you can observe. The
maximum correlation between two measures is the square root of the product of
their reliabilities, so nothing else about the hypothesis matters until that
ceiling is high enough. If ERN scores are reliable at .50 and the second measure
at .90, no observed correlation can exceed .67 in absolute value. Raise the ERN
score reliability to .90 and the ceiling rises to .90 (Clayson & Miller, 2017b,
Table 1). A null individual-differences result under the first condition is weak
evidence about the relationship and strong evidence about the measurement.

Low reliability also costs statistical power. Take a simulated
independent-samples t test at alpha = .05 with a true standardized difference of
0.5. Power was .65 when scores were reliable at .70, .74 at .85, and .80 under
perfect reliability (reported in Clayson & Miller, 2017b). Holding true-score
variance constant, better reliability shrinks error variance, and shrinking error
variance is what buys the power.

Group comparisons are where it gets genuinely misleading. Suppose two groups
differ by the same true amount on two tasks and one task is measured with less
error. The effect size will be larger for the better-measured task. Simple-effects
tests following a significant interaction then track the measurement difference
rather than the construct (Clayson & Miller, 2017b). High measurement error can
also produce sign errors, so a clinical group can appear to show larger scores
than controls when the true difference runs the other way (Clayson, 2024).
Comparing groups whose scores differ in reliability means comparing measurement
quality alongside everything else. The only way to find out whether that is
happening is to estimate reliability separately in each group.
PsyRAT does that by construction: it estimates separate variance components
for every group and event, and every coefficient it reports is labeled with
the cell it came from.

None of this requires heroics. Report the reliability threshold you adopted and
the estimate you obtained, report how you estimated it, and justify the trial
cutoff you used for exclusions (Clayson & Miller, 2017b, Table 2). On the
threshold itself, .70 is a floor for exploratory work and .80 is more appropriate
in most contexts (Clayson & Miller, 2017b). Decisions about individual people, of
the kind a clinical application would require, are usually held to .90 or higher,
and ERN scores in clinical samples do not generally reach that yet (Clayson,
2025). [Chapter 17](17_reporting.md) works through what to write down.

## 1.2 What PsyRAT does and does not do

PsyRAT takes a long-format table of single-trial ERP scores and estimates the
reliability of the averaged scores those trials produce. Variance components come
from a Bayesian multilevel model, fit in CmdStan whenever it is installed, so
every coefficient arrives as a point estimate with an interval rather than as
a bare number (a native fallback engine covers a subset of designs with
approximate intervals; [Chapter 5](05_choosing-an-analysis.md) explains). The
models tolerate unequal trial counts across participants and cells, which is the
ordinary situation after artifact rejection and which coefficient alpha cannot
handle without truncating the data (Clayson, 2024).

From a single session, PsyRAT reports internal consistency for each group and
event, the number of trials needed to reach a reliability threshold you set, and
which participants clear it. Asking for subject-specific error variances adds
per-participant (subject-level) estimates, one coefficient per person rather than
one for the sample. It handles two-event difference scores and the four-event
difference-of-differences contrast. Supply an occasion column and it estimates
test-retest reliability. Split means with an items-per-split count are analyzed
as splits rather than as trials. Supply a continuous between-person dimension
and it estimates reliability as a function of that dimension. Reliability can then be reported conditionally rather than as one
number for a whole sample (Rast & Clayson, in press).
[Chapter 5](05_choosing-an-analysis.md) helps you choose, and
[Chapter 16](16_analysis-reference.md) lists every design.

PsyRAT also imports epoched ERP files (EEGLAB `.set` and EP Toolkit exports;
ERPLAB `.erp` files are averaged and import as reference only) and scores
single-trial component amplitudes into the table the analysis needs. That path
is [Chapter 8](08_import-score.md).

The toolbox is GUI-first, and this manual is organized around its screens. Every
analysis the Process New Data workflow runs is also available from
`psyrat_run`, which uses the same estimation pipeline and the same seed
contract, so the same inputs at the same seed on the same CmdStan build give
the same variance components. Neither route is the lesser one, and I
use both. [Chapter 14](14_scripting.md) is the scripting reference.

There are four things PsyRAT does not do. It does not record or preprocess EEG:
filtering, re-referencing, ocular correction, and epoching happen in your
existing pipeline before PsyRAT sees anything. It does not decide what counts as
adequate reliability, because you type the threshold and the toolbox applies it.
It does not establish validity, since reliability limits validity without
supplying it, and a dependable score can still index something other than the
construct it was named after. And it does not deliver a biomarker: these
coefficients say how precisely you measured, which is a precondition for a
clinical claim rather than a substitute for one (Clayson, 2025).

### Migrating from the ERA Toolbox

PsyRAT descends from the ERP Reliability Analysis (ERA) Toolbox, which is
documented in Clayson and Miller (2017a). That paper describes the predecessor,
not this toolbox, so cite it for ERA. For PsyRAT, cite Rocha et al. (2026) for
the estimators together with the release you ran; `README.md` carries the current
citation guidance.

If you used ERA, expect renaming rather than a port. Toolbox function names now
begin with `psyrat_` where they began with `era_`, so the launcher is
`psyrat_start`. Processed results are written with the `.psyrat` extension, and
the View Results file dialog lists only those. Read
[Chapter 6](06_preparing-data.md) through [Chapter 9](09_viewing-results.md)
rather than assuming ERA habits carry over: the input contract and every output
are documented there from scratch. I make no compatibility promise in either
direction, and ERA result files and ERA scripts are not supported here.

## 1.3 The three coefficient terms used everywhere

Three words do a great deal of work in this manual, and they are not
interchangeable.

Reliability is the umbrella term. Use it for the general idea that scores
consistently distinguish among people, or within a person over time. The manual
uses it that way and never as a stand-in for either coefficient below.

Dependability (the absolute-error coefficient, written phi) asks how well a
participant's observed score reflects their absolute level on the construct
across the conditions being generalized over. It charges systematic differences
among trials or occasions to error (Rocha et al., 2026). If the value
of the score matters (a cut score, a change score, a comparison against a norm),
dependability is the coefficient you want. It is the toolbox's default.

Generalizability (the relative-error coefficient) asks how well the scores
preserve relative differences among people, and it leaves facet main effects out
of the error term (Rocha et al., 2026). If a systematic shift between occasions
is expected and irrelevant (habituation, say), and all that matters is whether
the same participants stay high, generalizability is the coefficient that ignores
the shift.

The two coincide when facet means do not differ. If scores do not systematically
differ across trials or occasions, absolute and relative error variance are the
same quantity, and the two coefficients are equal (Rocha et al., 2026).
Otherwise dependability is the smaller of the two, because absolute error
variance contains every term relative error variance contains and adds the facet
main effects on top.

Both coefficients are ratios of universe-score variance to universe-score
variance plus error, which is why both live on the same 0-to-1 scale as a
classical reliability coefficient and are read the same way. What differs between
them is only what got counted as error. [Chapter 4](04_gtheory.md) gives the
formulas and restates each of them in words. Everywhere else, the label tells you
which one you are looking at. The G-Theory Coefficient control on the viewer
selects between Dependability and Generalizability, and the exported overall
and cutoff tables name the coefficient that produced them.
