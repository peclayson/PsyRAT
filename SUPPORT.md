# Getting help with PsyRAT

## Reporting bugs

Please open an issue on the
[issue tracker](https://github.com/peclayson/PsyRAT/issues) and include:

- the **PsyRAT version string** (printed in the startup banner, returned by
  `psyrat_defineversion`, and written into every exported table header);
- the **`*_runconfig.json` sidecar** saved with the analysis, when one exists
  (it records the seed and every estimation setting, and is usually enough to
  reproduce the problem without your data);
- what you expected, what happened instead, and the full text of any error.

The bug-report form on the issue tracker asks for these items, along with the
MATLAB release, the operating system, and the dependency versions.

## Before you post

The issue tracker is public. Do not post or attach participant-level data,
participant IDs, `.psyrat` result files, exported tables that list
participants, or files from a leftover `psyrat_stan_*` folder, which can hold
your data and per-participant estimates.

Several files and messages also carry local paths. The `*_runconfig.json`
sidecar records the full paths of the input file and the output folder. An
`_error.txt` report and the Command Window output can name the save folder in
an error message, and the report's stack trace lists full paths to the toolbox
files. These paths can show your user name and folder names, so delete the
paths in the sidecar and shorten each path in an error report so that it starts
at the PsyRAT folder (Section 7.6 of the
[user manual](documentation/manual/07_processing.md) describes both files).

Participant IDs can appear in output as well. Some error messages name a
participant, some checks print a table of participant IDs before the error,
and the Import + Score status log names each ERP file, whose name often
contains a participant ID. Replace the IDs, or delete the table and say how
many participants it listed. When possible, reproduce the problem with one of
the example data files in [`test_data/`](test_data/).

## Questions

For usage questions, please also use the issue tracker; questions that recur
become documentation.
[Chapter 18 of the user manual](documentation/manual/18_troubleshooting.md)
covers the most common error messages and questions.

## Pull requests and feature requests

PsyRAT does not accept pull requests or other code contributions. It is
developed in a separate repository, and each release is published at
https://github.com/peclayson/PsyRAT. To suggest a new analysis or option, use
the feature-request form on the issue tracker.
