# PsyRAT User Manual

This manual documents PsyRAT 0.1.1-beta, the Psychophysiologist's Reliability
Analysis Toolbox: a MATLAB toolbox that uses generalizability theory to
evaluate the psychometric reliability of event-related potential (ERP) scores,
with Bayesian variance-component estimation via CmdStan.

The version above is the only PsyRAT version number printed in this manual.
Your installed version appears in the startup banner when you run
`psyrat_start`. The header of every table the toolbox exports names the
version that produced the results, which is your installed version for a fresh
run but an earlier one for a results file made with an earlier release (a file
with no version stamp shows "(not recorded)").

## How to read this manual

- New to the toolbox? Start with [Chapter 2, Quick start](02_quick-start.md),
  then work through the tutorials (Chapters 10 to 13).
- New to generalizability theory? Read
  [Chapter 1](01_why-reliability.md) and [Chapter 4](04_gtheory.md) first.
- Hit an error message? Go straight to
  [Chapter 18, Troubleshooting](18_troubleshooting.md).
- Migrating from the ERP Reliability Analysis (ERA) Toolbox?
  [Chapter 1](01_why-reliability.md) notes what changed, and Chapters 6 to 9
  (Part III) are the new reference.

## Contents

### Part I. Orientation
1. [Why analyze the reliability of ERP scores?](01_why-reliability.md)
2. [Quick start: your first analysis](02_quick-start.md)
3. [Installation and dependencies](03_installation.md)

### Part II. Concepts
4. [Generalizability theory for ERP scores](04_gtheory.md)
5. [Choosing an analysis](05_choosing-an-analysis.md)

### Part III. The toolbox, screen by screen
6. [Preparing your data](06_preparing-data.md)
7. [Processing data, screen by screen](07_processing.md)
8. [Importing and scoring ERP data](08_import-score.md)
9. [Viewing results](09_viewing-results.md)

### Part IV. Tutorials and scripting
10. [Tutorial 1: Single-session ERN reliability](10_tutorial-single-session.md)
11. [Tutorial 2: Test-retest reliability](11_tutorial-test-retest.md)
12. [Tutorial 3: Difference-score reliability](12_tutorial-difference-scores.md)
13. [Tutorial 4: Reliability from data splits](13_tutorial-splits.md)
14. [Scripting with psyrat_run](14_scripting.md)

### Part V. Advanced use
15. [Advanced designs](15_advanced-designs.md)
16. [Analysis reference](16_analysis-reference.md)

### Part VI. Reporting and help
17. [Interpreting and reporting results](17_reporting.md)
18. [Troubleshooting and FAQ](18_troubleshooting.md)

### Appendices
- [Appendix A. Glossary](appendix-a_glossary.md)
- [Appendix B. References](appendix-b_references.md)
- [Appendix C. Reprocessing saved gamma results](appendix-c_reprocessing.md)

## Conventions used in this manual

- GUI control names appear with their exact on-screen capitalization: click
  Analyze, select the Preferences button.
- Dialog and window titles are given exactly as the toolbox sets them, in
  plain text.
- In running text, file names and MATLAB commands appear in `code font`.
- The tutorial datasets ship with the toolbox in `test_data/`, and the
  expected results printed in the tutorials were produced from those exact
  files at the settings each tutorial states.
- Reliability is the umbrella term. Dependability refers to the
  absolute-error coefficient and generalizability to the relative-error
  coefficient; the two are not interchangeable, and the manual always names
  the one it means.

macOS file panels do not display a title. On a Mac, a file dialog the manual
describes as titled Data opens as an untitled panel.

## Acknowledgments

I thank Gavin Heindorf and Bohyun Park for their help with beta testing.

Development of PsyRAT was supported by the National Institute of Mental
Health of the National Institutes of Health under award number R01MH128208.
I am grateful to NIMH for its generous support of this work. The content is
solely the responsibility of the author and does not necessarily represent
the official views of the National Institutes of Health.

Claude (Anthropic) and ChatGPT (OpenAI) were used to write code and
documentation for PsyRAT, including this manual. Responsibility for the
resulting code and text rests with the author.
