# 3. Installation and dependencies

PsyRAT itself is plain MATLAB code, and MATLAB code needs no installing. What
needs installing is CmdStan, which does the Bayesian estimation, and the C++
toolchain CmdStan needs to compile each model. That is the whole problem. This
chapter is deliberately short on click-by-click steps, because the ones I wrote
for the predecessor toolbox ran to pages and then rotted with the operating
systems they described. Requirements keep. Dialogs do not.

If the toolbox is already working on your machine, skip to
[Section 3.4](#34-verify-the-installation) and confirm it.

## 3.1 What you need

| Component | Version | Where it comes from |
|---|---|---|
| MATLAB | Tested on R2025b (Update 5) | You install it |
| CmdStan | 2.26 or newer; 2.38.0 is the targeted and tested release | The guided installer, or you |
| MatlabStan (PsyRAT fork) | 2.15.1.0 | Ships with the toolbox |
| MatlabProcessManager | 0.5.1 | Ships with the toolbox |
| C++ toolchain | make and a C++ compiler | Your platform (see below) |
| Statistics and Machine Learning Toolbox | Conditional, see below | MathWorks |
| Deep Learning Toolbox | Optional | MathWorks |

**MATLAB.** The tested release is R2025b (Update 5), and every release claim for
the toolbox rests on runs I made on my own machine. No minimum release is
enforced anywhere in the code, and the precise minimum has not been pinned.
The effective floor is R2022a: the dynamic-reliability surface plot and the
difference-of-differences heatmaps call `clim`, introduced in that release, and the export code calls `exportgraphics`
and `writetable` with `WriteMode`, introduced in R2020a, all without fallback.
Beyond those the toolbox relies on `table`, `string`, and `categorical` (the
test suite also uses the `matlab.unittest` framework), all of which are old
features by now. Between
R2022a and R2025b, the honest answer is that nobody has checked.

**CmdStan.** CmdStan is not bundled and must be installed separately. The
generated Stan models use the `array[N] T name;` declaration syntax introduced in
Stan 2.26, so anything older cannot compile them. Version 2.38.0 is the version
the guided installer fetches and the version used for the accuracy and
integration runs. If you already have a CmdStan on your machine, it is checked at
startup: an installed version below 2.26 produces a warning when the toolbox
opens and a hard error before the first model is compiled. The check is
deliberately permissive, so a CmdStan whose version cannot be read is allowed
through rather than blocked. CmdStan before 2.33 reports a single effective
sample size rather than separate bulk and tail values; PsyRAT still runs on
it, but it records those diagnostics as degraded, and the viewer asks you to
acknowledge that each time you open such a result, so use 2.33 or newer,
ideally the tested 2.38.0. To replace an old CmdStan that is on the MATLAB
path, answer Yes to the update offer described in Section 3.4, which removes
the old copy and installs the tested one. Running `psyrat_installdependents`
directly keeps any CmdStan the toolbox can already find, so an old install
that an environment variable names, or that sits in one of the searched
folders, has to be removed (or its variable cleared) first.

**The two bundled dependencies.** MatlabStan and MatlabProcessManager ship inside
the toolbox in `bundled_dependents/`, and there is nothing to install for either.
The bundled MatlabStan is a patched fork, not the upstream release, and the
patches are load-bearing. They add CmdStan auto-discovery and within-chain
threading, and they repair an indexing bug in the draw extraction. They also
remove three Statistics and Machine Learning Toolbox calls: two in the fit
reader, which would otherwise have made that toolbox mandatory on the default
path, and one in a separately licensed subpackage. Use the
bundled copy. If a different MatlabStan is on your MATLAB path ahead of it, the
startup check says so.

**C++ toolchain.** CmdStan compiles every model with make and a C++ compiler, so
one has to be present before CmdStan can be built or used. On macOS that is the
Xcode Command Line Tools, and the installer starts Apple's installer for you if
they are missing. On Windows it is GNU make plus g++, both on the PATH; Rtools provides
them, and the installer cannot install them for you. On Linux it is make plus
g++, from `build-essential` on Debian and Ubuntu or the Development Tools group
on Fedora and RHEL.

**Statistics and Machine Learning Toolbox.** For most analyses you do not need
it. Three designs need it even on the default CmdStan path, all of them the
gamma (scaled chi-square) family combined with concurrent difference scores:
analysis 7 (two-event difference), analysis 8 (its subject-level variant), and
analysis 10 (two-facet test-retest difference). Those three read their results
through a Gaussian copula that calls `normcdf` and `gaminv`, and neither
function has a base-MATLAB equivalent. Both gamma scale parameterizations
(location-scale and log-nu; see [Chapter 5](05_choosing-an-analysis.md)) are
affected. Without
the toolbox the run is refused before anything is compiled, with a message
saying that the gamma family with residual covariance needs the Statistics and
Machine Learning Toolbox; the fixes it offers are to install the toolbox, set
residual covariance to No, or use the Gaussian family. If you plan to use those
designs (see [Chapter 15](15_advanced-designs.md)), confirm the toolbox is
installed first. Whether other designs have a similar dependency has not been settled.
The toolbox is also required by both native estimation engines, which are the
CmdStan-free fallbacks described in [Chapter 5](05_choosing-an-analysis.md).

**Deep Learning Toolbox.** Optional. It supplies automatic-differentiation
gradients for the native HMC engine, which runs without it on numeric gradients
instead.

## 3.2 Install

FIRST, get the toolbox. Download `PsyRAT.zip` from
[the Releases page](https://github.com/peclayson/PsyRAT/releases) and unzip
it, or clone the repository.

Put it somewhere whose full path contains no spaces. CmdStan does not accept
whitespace in paths or filenames, and the toolbox enforces that on the files it
writes: choose an output folder or an output filename containing a space and you
get a Whitespace detected dialog rather than a failed run. Nothing enforces it on
the folder you unpack the toolbox into, or on the folder CmdStan is built in, so
that one is on you. `~/Documents/MATLAB/PsyRAT` is a safe choice.
`~/My Lab Stuff/PsyRAT` is not.

Then run this in MATLAB, once, from inside the toolbox folder:

```matlab
psyrat_start
```

That is also the installation step. `psyrat_start` adds its own directories to
the MATLAB path (the toolbox root, `subroutines/`, and `bundled_dependents/`) and
saves the path, so later sessions find it from anywhere. On a shared or
centrally managed MATLAB where the path cannot be saved, you get a warning
instead, and the toolbox is on the path for that session only; add the toolbox
folder (and, once the installer has run, `addpath(genpath(...))` of your
`PsyRATDependents` folder) to your `startup.m` so they load in every session.
Do not `addpath(genpath(...))` the whole repository. Doing so also adds tooling and
artifact folders, and a stale copy of a function inside one of those will shadow
the real file with no warning at all.

If the Stan dependencies are missing, `psyrat_start` prints a warning that
begins "Warning: Dependencies for PsyRAT were not found.", stops
before the home screen, and raises a small
window offering two buttons: Install Dependents and Exit. Selecting Install
Dependents runs the guided installer, which you can also run directly with
`psyrat_installdependents`.

The installer asks one question first, in a save dialog headed "Where would you
like to save the directory for the dependents?". Answer it with a folder whose
path has no spaces. The installer creates a folder called `PsyRATDependents`
there and builds CmdStan inside it. Choose `~/Documents/MATLAB` if you have no
preference, because that is one of the locations the toolbox searches on its
own later ([Section 3.3](#33-installing-or-pointing-at-cmdstan-yourself)), and
a `PsyRATDependents` folder somewhere else is then found only through the saved
MATLAB path.

From there the installer works down a checklist and shows its progress in a
window named Installation Summary (its heading reads Installation Status).

![The Installation Summary window after a completed install. Each row is one dependency; during the install its status reads Not Installed (Not Built for the build row) until that step finishes, then Installed (Built).](images/install_summary.png)

The window lists five rows: Command Line Tools, CmdStan Files, CmdStan Build,
Matlab Process Manager, and Matlab Stan. Each carries a status of Installed or
Not Installed (Built or Not Built for the build row), and each is redrawn as the
step behind it completes. The same status text also goes to the Command Window,
which is where to look if something fails, because the failure message is printed
there and not in the window. The line at the bottom reading "Installation may
take a while..." is not a figure of speech. Building CmdStan compiles a large
amount of C++, takes several minutes on a current laptop, and uses all but one
of your cores while it does (all of them on a two-core machine). A failure
stops the installer where it is: the Installation Summary window stays open,
the home screen does not come back, and the MATLAB current folder may have been
left inside `PsyRATDependents` or the CmdStan folder. Close the window, fix
what the message names, and run `psyrat_start` again from the toolbox folder.

What the installer actually does: it checks for a C++ toolchain, downloads the
CmdStan release tarball for the targeted version from the Stan project's GitHub
releases, unpacks it into `PsyRATDependents`, and builds it with `make`. It
then copies the bundled MatlabStan and MatlabProcessManager in beside it when
neither is already on the path (in the flow above they are, because
`psyrat_start` added `bundled_dependents/` first, so those two rows read
Installed at once), adds the whole `PsyRATDependents` tree to the MATLAB path,
which is how the toolbox finds the new CmdStan when the folder is not one of
the locations listed in Section 3.3, saves the path, and relaunches
`psyrat_start`.

On macOS it inserts one extra step. If the Xcode Command Line Tools are not
present, the installer starts Apple's own installer and then waits behind a
window reading "Please indicate when the XCode installation is complete", with a
single button labeled Installation Complete. Click it when Apple's installer has
finished, not when it has started. The toolbox rechecks the install location
rather than taking the button press as proof, so clicking early produces an
error telling you to install the command line tools manually.

On Windows the installer cannot set up the toolchain. If GNU make and g++ are not
already on the PATH it stops with a message pointing you at Rtools, and you have
to install that yourself before rerunning. Windows is not part of the test
matrix: MatlabStan and CmdStan both support Windows, and the toolbox relies on
that support, but I have not verified a Windows run.

On Linux the installer runs, and it is EXPERIMENTAL. That branch mirrors the
macOS path (check for make and g++, download the CmdStan tarball, build it), and
I have not verified it on a Linux machine. It announces this with
a warning when it starts. Estimation itself is fine on Linux once CmdStan is
discoverable. If the automatic install fails, install CmdStan yourself and point
the toolbox at it as described next.

## 3.3 Installing or pointing at CmdStan yourself

You do not have to use the guided installer. A CmdStan you built yourself, or one
that came with another Stan interface, works as long as the toolbox can find it
and it meets the minimum version.

Discovery happens in a fixed order, and the first candidate that looks like a
real CmdStan root (a `makefile` and a `bin` directory) is the one used:

1. The environment variables `STAN_HOME`, `CMDSTAN_HOME`, and `CMDSTAN`, in that
   order. (One further name, checked ahead of these, is reserved for the test
   harness.)
2. A `cmdstan-*` folder inside `PsyRATDependents` under your home directory:
   `Documents/MATLAB/PsyRATDependents`, `Documents/PsyRATDependents`, or
   `MATLAB/PsyRATDependents`. A `cmdstan-*` folder directly in your home
   directory is checked too. Where one location holds several versioned
   folders, the newest version wins; pin an older one with an environment
   variable if you need it.
3. Whatever is reachable on the MATLAB path.

An environment variable therefore always beats a toolbox-managed install, which
is what you want when you are deliberately testing against a different CmdStan.
The home-directory install beats a MATLAB path entry, which is also what you
want, since it keeps an unrelated CmdStan that happens to be on the path from
displacing the one built for the toolbox.

The startup check in `psyrat_start` uses this same discovery order, so any
CmdStan the list above can find, including one named only by an environment
variable, satisfies it. There is no need to add the CmdStan root to your
MATLAB path just to get past the launch check. (Earlier versions checked the
MATLAB path alone, and an environment-variable-only setup produced a spurious
missing-dependencies report at launch; that behavior is gone.)

For building CmdStan itself, follow the Stan project's own instructions at
<https://mc-stan.org/users/interfaces/cmdstan> rather than anything reproduced
here. Those steps are maintained by people who test them on every release, which
is more than can be said for a printed manual.

## 3.4 Verify the installation

Run `psyrat_start` and read the Command Window before you touch the buttons. A
working installation prints six lines:

```
Ensuring dependents are found in the Matlab path
PsyRAT files found
CmdStan, MatlabProcessManager, and MatlabStan files found
CmdStan version is current
MatlabProcessManager version is current
MatlabStan version is current
```

Three blank lines later comes the banner, reading "PsyRAT Version" followed by
the version string of your installation. Report that string in any bug report.
The first time you run `psyrat_start` from a new toolbox folder, the second
line instead reads "Directories for PsyRAT files added to Matlab
path". That is the toolbox putting itself on the path, and later launches
print the line shown above, unless the path could not be saved.

Those checks run once per MATLAB session, and so does the release check
described below. Returning to the home screen with a Back button, canceling a
file dialog, or running `psyrat_start` again in the same session prints the
first three lines as before. Next comes one line saying that the dependency
and release checks were reused from this session, with the time of the most
recent fresh check and the release-check outcome, and then the three version
lines again from that check. Nothing is spawned or downloaded on a reused
launch. If that check found a CmdStan below the minimum PsyRAT requires, the
reused launch does not repeat the full warning; a line beginning "Reminder:
the CmdStan found earlier this session is below the minimum PsyRAT requires"
appears in its place, and estimation still refuses to compile until a newer
CmdStan is installed. To repeat the checks without restarting MATLAB, run
`psyrat_startgate_cache('clear')` and then `psyrat_start`. The toolbox also
repeats them on its own, and says why, when the toolbox, CmdStan, MatlabStan,
or MatlabProcessManager it found has changed since that check. The relaunch
after the guided installer or the updater always starts from a fresh check,
and restarting MATLAB always repeats them.

![The MATLAB Command Window immediately after psyrat_start, showing the dependency check lines and the version banner.](images/console_banner.png)

The three version lines are advisory. For CmdStan and MatlabProcessManager, a
dependency behind the version the toolbox targets produces a line saying a new
version is available with an offer to update, and one ahead of it produces a
warning that it is newer than the version the toolbox was tested with. The
offer is a window naming the dependents that are behind, reading "It is
recommended that you update the dependents. Would you like to do so now?"
with a warning that updating deletes the old dependents directories, and it
offers Yes and No. If the release check has also found a newer PsyRAT release,
the toolbox update prompt takes precedence on every launch until the toolbox is
current, and the dependents offer appears only after that. The MatlabStan line
reports something different: whether the copy in use carries the patches
PsyRAT ships, with a third outcome when the patch level cannot be determined. A
copy without the patches (the line reads "The MatlabStan in use is not the
patched version PsyRAT ships") is named in the same update offer. PsyRAT's
bundled MatlabStan carries a patch-level stamp; when the copy in use is that
bundle inside the PsyRAT folder and its stamp is older than this release
expects, or missing, the line instead reads "MatlabStan: PsyRAT's bundled
pinned copy is in use (updated with the toolbox, not the dependents)": the
bundle cannot be replaced by the dependents updater, so no update is offered
for it, and updating the toolbox itself is what refreshes it. Estimation is
not blocked in any of these cases.

You may also see a message from the release check, and on a beta build you should
expect one. The toolbox asks GitHub for the latest published release and compares
it against the version you are running. When there is no published stable release
to compare against, or when your version is ahead of the newest one, it says so
and recommends caution. That message is informational, and it reports nothing
about whether your installation works.

The last check is a real fit. Run the analysis in
[Chapter 2](02_quick-start.md) and let it finish. If you need to stop it,
interrupting MATLAB with Ctrl+C stops MATLAB waiting for the sampler, not the
sampler itself, because the models run in C++ and not in MATLAB; to stop a run
for real, end each CmdStan process yourself (see "Ctrl+C does not stop a run"
in [Chapter 18](18_troubleshooting.md)). Every estimation starts by compiling
the generated Stan model to a binary before sampling, which adds about 8 s per
generated model on top of the sampling itself (the toolbox names each
generated model uniquely, so the compile happens on every run, not just the
first). The binaries land in a folder named `PsyRAT_CmdStan` on your Desktop,
or in the current MATLAB directory when there is no Desktop folder; if that
path contains a space, they go to a temporary folder instead, and the Command
Window prints a warning naming it. The folder is safe to empty between
sessions. If that first run finishes and produces a results window, everything
downstream of this chapter is working.

Success! CmdStan is installed and discoverable, the bundled MatlabStan is talking
to it, and the toolbox is on your MATLAB path and, unless the path could not be
saved, will still be there in your next session. If any of the checks above failed instead, the symptom-by-symptom
fixes are in [Chapter 18](18_troubleshooting.md).
