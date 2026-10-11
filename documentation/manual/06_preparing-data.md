# 6. Preparing your data

PsyRAT reads one kind of file, and everything downstream depends on getting it
right. Keep this chapter open while you export from EEGLAB, ERPLAB, or your own
scripts. It is the reference the tutorials in Part IV assume you have read.

## 6.1 Long format: one row per single-trial score

The input is a long-format table with one row per single-trial score. Two
columns are required: a participant identifier and a numeric measurement. Three
more are optional: group, which splits the sample into separately modeled
strata, and event (condition) and occasion, which enter the model as facets.

These ten rows are lifted verbatim from `test_data/manual_ern_singlesession.csv`,
the file Tutorial 1 uses:

| id | group | event | meas |
|---|---|---|---|
| 101 | GroupA | inc_err | -12.041412 |
| 101 | GroupA | inc_err | -6.130699 |
| 101 | GroupA | inc_err | -1.701634 |
| 101 | GroupA | inc_err | -1.552103 |
| 101 | GroupA | inc_err | -3.576956 |
| 101 | GroupA | inc_cor | -3.551382 |
| 101 | GroupA | inc_cor | -2.451203 |
| 101 | GroupA | inc_cor | -11.254837 |
| 101 | GroupA | inc_cor | -4.363838 |
| 101 | GroupA | inc_cor | -4.467972 |

Read the columns left to right. `id` is the participant, `group` is the
between-person grouping, `event` is the condition label, and `meas` is the
single-trial amplitude in microvolts (the ERN on error trials). That last
column is the only one the model treats as a score. The window shown here straddles a condition boundary: the first five
rows are error trials on incongruent stimuli and the next five are correct
trials on the same stimuli. Nothing marks the boundary except the `event` value,
and nothing needs to.

Rows may be interleaved across participants and conditions. Within a
participant and condition, though, row order matters. PsyRAT numbers each
participant's trials in the order the rows appear, and that number is the
trial level shared across participants in the trial main effect (a
participant's third row within a condition is that condition's trial 3 for
everyone). The concurrent difference designs pair trials the same way, the
k-th trial of one event with the k-th trial of the other. Keep every
participant's rows in the order the trials occurred, and do not shuffle rows
within a participant.

Trial counts do not have to match. Participants routinely differ in how many
artifact-free error trials they contribute, and the trial-count projection in
the D-study exists precisely because they do. What must match across
participants is the number of distinct facet levels, which Section 6.5 covers.

Do not average before loading. A file of participant means carries no
within-person replicates, so the residual variance the model needs is not
identified and the reliability estimate has nothing to be built from. The one
exception is the splits format below, where each row is a mean of a known number
of items and you declare how many.

## 6.2 Your column names, not mine

Header names are free. The GUI maps columns to roles through the popup menus on
the Specify Inputs screen ([Chapter 7](07_processing.md)), so a file whose
headers read `subjid`, `cond`, and `ern_uv` works exactly as well as one written
to the canonical headers `id`, `meas`, `group`, `event`, and `time`. Scripted
callers pass the same mapping as arguments (`idcol`, `meascol`, `groupcol`,
`eventcol`, `timecol`, and for the optional roles `dim1col`, `dim2col`, and
`weightcol`); see [Chapter 14](14_scripting.md).

The canonical names buy one thing: they are what auto-detection finds first.
When the Specify Inputs screen opens, PsyRAT preselects each role by matching
headers case-insensitively against a short candidate list. Participant ID matches
Subject, ID, Participant, SubjID, or Subj. Measurement matches Measurement or
Meas. Group ID matches Group. Event Type matches Event or Type. Occasion ID
matches Time, Occasion, or Session. A header has to match one of those
candidates exactly (apart from case), so `subj_id` and `ERN` match nothing.

If no header matches Measurement, PsyRAT falls back to the LAST numeric column
that no other role has already claimed. In a long-format export the score is the
melted value column and usually does come last, so the fallback is right more
often than it is wrong. It is still a guess. Confirm the Measurement popup before
you click Analyze, because a run on the wrong column completes normally and
reports the reliability of whatever you pointed at.

Dimension columns and the Items per split (n_i) column are never auto-detected. They
stay on `(none)` until you pick them, which is deliberate: selecting either one
changes which analysis runs.

One header caveat, and it will surprise you the first time. MATLAB coerces table
headers into valid variable names as it reads the file, so `ERN (uV)` arrives as
`ERN_uV_` (the space is dropped and each parenthesis becomes an underscore), and
that coerced name is what the popups show. Your data are untouched. If you
want the popups to read cleanly, export headers that start with a letter and
contain only letters, digits, and underscores.

## 6.3 The splits format: one row per split, with a weight

Reliability from data splits takes a different input, and the difference is not
cosmetic. Each row is one SPLIT rather than one trial, the measurement column
holds the split MEAN, and an Items per split (n_i) column (the weight column,
for short) reports how many items were averaged into that split. That column
is required, and selecting it on the Items per split (n_i) popup is what
declares a splits analysis. Without it the split-level model
cannot be identified, because the sampling variance of a split mean of `n_i`
items is the per-item residual variance divided by `n_i`.

These ten rows come verbatim from `test_data/manual_ern_splits.csv`:

| id | meas | weight |
|---|---|---|
| 301 | -4.171339 | 24 |
| 301 | -6.509574 | 40 |
| 301 | -4.995126 | 28 |
| 301 | -4.531105 | 28 |
| 301 | -4.46002 | 6 |
| 301 | -2.156673 | 8 |
| 301 | -4.495167 | 26 |
| 301 | -5.342279 | 20 |
| 301 | -8.377559 | 30 |
| 301 | -3.916375 | 35 |

Participant 301 contributes sixteen splits; the ten rows shown are its first
ten, with weights of 24, 40, 28, 28, 6, 8, 26, 20, 30, and 35 items. The mean of
the 24 items in the first split is -4.17 microvolts. Because those counts are
unequal, this file is declared Nonparallel on the Split type popup; a file whose splits all held the same number of items
would be declared Parallel instead. The whole shipped file runs 4 to 40 items per
split across 36 participants.

Every weight must be a positive whole number. Fractions, zeros, negatives, and
non-finite entries are refused when you click Analyze, along with the other
data checks in Section 6.5. Group, event, and occasion columns
load here exactly as they do for single-trial input, and Parallel splits use
them as single-trial input does. One restriction applies to Nonparallel splits:
their single-occasion model is fitted as one pooled model over all
participants, so a run with a Group ID or Event Type column selected stops with
a message saying so rather than fitting one model per stratum. The fix the
message names is to remove the column or to add an Occasion ID column, which
routes to the group- and event-aware two-facet split design
([Chapter 13](13_tutorial-splits.md)).

The two split types are not interchangeable declarations of taste.
Nonparallel splits are identified ONLY by variation in the weights, so a
nonparallel run on all-equal weights is refused outright, and a run whose weights
barely vary is allowed but flagged. Parallel splits assume equal weights and
reuse the standard one-facet machinery with the split mean as the score. Declare
what your splitting procedure actually produced. [Chapter 13](13_tutorial-splits.md)
walks the workflow end to end.

## 6.4 Dimension columns for dynamic reliability

Dynamic (conditional) reliability reports reliability as a function of one or two
continuous between-person covariates, and you enable it by selecting a Dimension
1 column. Dimension 2 is optional and adds a second predictor with its
interaction. Selecting Dimension 2 without Dimension 1 is refused.

A dimension column must be numeric, and it must be constant within participant.
One value per person, repeated down that person's rows, is what the format
expects (a questionnaire total or a symptom score, for example). Fill the value
on every row: a blank Dimension cell anywhere stops the run with the
missing-cells message in Section 6.5. If a participant's rows disagree, the
run stops and names the participant and the Dimension (as dim1 or dim2). A
column with no between-participant variance cannot be standardized and is also
refused.

Nothing else about the file changes. The dimension column sits alongside the
usual `id`, `meas`, and facet columns, one more column in the same long-format
table. [Chapter 15](15_advanced-designs.md) covers what the analysis reports.

## 6.5 What Analyze checks, and what to do about it

The data checks do not run when the file picker closes. Opening the file only
reads it into a table. The checks below run when you click Analyze, once per
selected measurement, and (this trips people up) AFTER you have chosen a name
and folder for the output. The one exception is the setup-screen wording quoted
in the non-numeric row below, which fires once for all selected measurements
before the save dialog opens. [Chapter 7](07_processing.md), Section 7.5,
follows that sequence in order.

The table gives the message fragment you will see and the fix. Messages are
quoted from what the toolbox prints, trimmed to the identifying phrase.

| What you will see | How to fix |
|---|---|
| "Input data include missing cells, which CmdStan cannot process." (followed by each column with missing cells, how many, the first such row of the prepared table, and that row's id) | Remove or impute the affected rows so every row has an id, a measurement, and a value for every facet you selected. A value that cannot be read as a number, such as `NA`, text, or an empty cell, is imported as missing, so the message also catches non-numeric entries in the measurement column. |
| "Participants differ in the number of distinct events" (with a printed table of ids and `Number_of_Events`) | Read the printed table, find the ids whose count differs from the rest, and either restore the missing condition or drop those participants. |
| "Participants differ in the number of distinct occasion cells" (with a printed table of ids and `Number_of_Occasions`) | Same remedy, applied to the event-by-occasion grid (or, when no event column is selected, to occasions alone; the check runs either way). |
| "The data in the measurement column is not numeric" (or, from the setup screen, "The following selected measurement columns are not numeric:") | Point Measurement at a numeric column. A text or date column selected by mistake is the usual cause. |
| "Column 'group' is a logical (true/false) array, which cannot be used as a label column." | Convert the column to numeric codes, a categorical, or text before loading, for example with `double(t.group)` or `categorical(t.group)`. |
| "has numeric values that convert to the SAME text label, so two or more distinct levels would be merged into one" | Supply that label column as text yourself, rounding or naming the levels the way you want them reported. |
| "The 'id' column contains missing values" (or 'group', 'event', 'time') | Give every row a label in that column. A missing label cannot be filled in safely, because every missing entry would collapse into one level. |
| "The dataset is a table with column headers but no data rows." | Re-export the file with its single-trial rows. |
| "Measurement column '...' contains N non-finite (Inf or -Inf) value(s)" | Remove or correct the infinite scores. Every score must be a finite number. |
| "Column '...' has N label(s) containing a comma, a double quote, or a line break" | Rename those labels. The exported CSV tables would otherwise misplace every later column of such a row. |
| "The occasion column 'time' has only one level" | An occasion column with one level routes the run to the two-facet test-retest design, whose occasion components cannot be identified from a single occasion. Set Occasion ID back to (none), or supply a second occasion. |

Three of those rows deserve more than a line.

The two participant-count checks compare COUNTS, not membership, and the error
text says so itself. A participant with four distinct events passes even if they
are not the same four events everybody else has. That is a real limit, not a
subtlety I am glossing: uniform counts are what the check can enforce cheaply on
a long-format table.

A genuinely incomplete design still runs. When both an
event and an occasion column are selected, PsyRAT warns rather than blocking,
reporting how many participant-by-event-by-occasion cells hold no observations
and naming the first few. In a single-session design there is no such warning,
so a mismatched-but-equal-count event set runs with no notice. Estimation is valid
with empty cells, because variance-component estimation does not require a
complete design. The estimates for an empty cell are model-based rather than
observed, and you should check that each empty cell reflects an intended
exclusion, such as artifact rejection, and not a data-preparation error.

Two more warnings belong to the same non-blocking class. "The dataset has only N
participant(s), fewer than the recommended minimum of 20" fires below 20
participants (variance components are unstable and prior-dominated with so
few); when a group column is selected the check runs per group and reads
"Group 'GroupA' has only N participant(s)..." instead. "Median trials per cell
is X, fewer than the recommended minimum of 10" fires when the median trial
count per cell (each participant crossed with whichever group, event, and
occasion columns are selected) falls below 10; splits runs skip it. None of
the three stops the run. In the GUI they appear together in one modal dialog
titled Preflight warnings, raised after the save dialog and before each fit,
and the run waits until you click OK; the same lines print in the Command
Window, which is the only channel on the `psyrat_run` route
([Chapter 14](14_scripting.md)). The same dialog carries the other advisory
checks: labels that differ only in letter case or surrounding spaces, exact
duplicate rows, a participant id listed under more than one group, and, for
splits, weights that barely vary, unequal weights declared Parallel, or
participants with different numbers of splits.

The numeric-label collision is the check people meet unprepared. MATLAB renders
numbers with about five significant digits, so two event codes that differ only
beyond that point convert to the same text label and would silently merge into
one level. PsyRAT refuses instead. Consider this example: two numeric event
codes of `0.123456789` and `0.123456111` both convert to the label `0.12346`,
so two conditions that look distinct in the spreadsheet would arrive as one.
Writing label columns as text (`t1` and `t2`, `cor` and `err`) in the first
place avoids the whole question.

## 6.6 File types PsyRAT reads

The file picker offers a combined All Readable Files filter covering `.xlsx`,
`.xls`, `.csv`, `.dat`, `.txt`, and `.ods`, plus one filter per type. Underneath,
the reader also accepts `.xlsb`, `.xlsm`, `.xltm`, and `.xltx`, so those open if
you type the name or widen the filter to All Files. Anything else is refused with
a message about unsupported file types. The extension must be lowercase: a file
saved as `DATA.CSV` is refused, so rename it `DATA.csv` first.

Delimited text does not need a declaration. PsyRAT first asks MATLAB to read the
file on the strength of its extension, and if that succeeds but yields only one
column it retries with a comma, then a space, then a tab, then a semicolon, then
a pipe, stopping at the first delimiter that produces more than one column. A file that
survives none of those raises a loading error. Comma-separated files are the
safest choice, and they are what the tutorial datasets ship as.

## 6.7 Before you load

Run down this list once against a new export. It takes a minute and saves the
round trip through a save dialog and a failed run.

1. One row per single-trial score, unaveraged.
2. A participant column and a numeric measurement column, both present in every row.
3. No blank cells in any column you intend to use, and no infinite values in the measurement column.
4. Facet labels as text where you can manage it (`t1`, `err`), not as bare numbers, and with no commas, double quotes, or line breaks in them.
5. The same number of distinct events, and of event-by-occasion cells, for every participant.
6. An occasion column with at least two levels, or none at all.
7. Plain headers: a letter first, then letters, digits, underscores.
8. For splits, a positive whole-number weight on every row, and a decision about Parallel or Nonparallel.
9. For dynamic reliability, one numeric Dimension value per participant, constant down that participant's rows.
10. A save location whose full path contains no spaces (a CmdStan restriction, covered in Section 7.5).

The last item is not about the data at all, and it is one of the mistakes that
most often stop a first run. With the ten satisfied you have a file PsyRAT
will load, and [Chapter 7](07_processing.md) picks up at the file picker.
