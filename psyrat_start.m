function psyrat_start
%Open the launcher gui for the Psychophysiologist's Reliability Analysis Toolbox (PsyRAT).
%
%For scripted (non-gui) analyses that follow the same estimation pipeline
% and reproducibility contract, see psyrat_run,
% documentation/scripted_use.md (quickstart), and chapter 14 of the user
% manual (the full scripting reference).
%

%The PsyRAT toolbox uses generalizability theory as a method for evaluating
%reliability of ERP data. Dependability estimates (generalizability-theory
%analog to reliability) can be computed for any number of groups or events.
%The influence of number of trials on dependability of measurements will
%also be determined, and a recommended cutoff for inclusion of ERP data
%will be provided based on the stability of measurement as the number of
%trials included in a single-subject average for a given group and event
%increases.
%
%A description of how to apply generalizability theory to ERP data can be
% found in
%
% Clayson, P. E., & Miller, G. A. (2017). ERP Reliability Analysis
% (ERA) Toolbox: An open-source toolbox for analyzing the reliability of
% event-related potentials. International Journal of Psychophysiology, 111,
% 68-79. doi: 10.1016/j.ijpsycho.2016.10.012
%
% Baldwin, S. A., Larson, M. J., & Clayson, P. E. (2015). The dependability
% of electrophysiological measurements of performance monitoring in a
% clinical sample: A generalizability and decision analysis of the ERN and
% Pe. Psychophysiology, 52, 790-800. doi: 10.1111/psyp.12401
%
% Clayson, P. E., Carbine, K. A., Baldwin, S. A., Olsen, J. A., &
% Larson, M. J. (2021). Using generalizability theory and the ERP
% Reliability Analysis (ERA) Toolbox for assessing test-retest reliability
% of ERP scores part 1: Algorithms, framework, and implementation.
% International Journal of Psychophysiology, 166, 174-187.
%
% Clayson, P. E., Brush, C. J., & Hajcak, G. (2021). Data quality
% and reliability metrics for event-related potentials (ERPs): The utility
% of subject-level reliability. International Journal of Psychophysiology,
% 165, 121-136.
%
%The notion of reporting estimates of reliability in all ERP studies and
% this toolbox are specifically discussed in
%
% Clayson, P. E., & Miller, G. A. (2017). Psychometric
% considerations in the measurement of event-related brain potentials:
% Guidelines for measurement and reporting. International Journal of
% Psychophysiology, 111, 57-67. doi: 10.1016/j.ijpsycho.2016.09.005
%
%
%Input
% There are no required inputs to execute this script.
% However, paths to the data to be processed or viewed will need to be
%  provided.
%
%Output
% This script will not output any variables to the workspace. However, if
%  the user chooses, various files (data, tables, plots) could be saved
%  from the guis initiated by this script.
%
%High-level matlab files contained in PsyRAT toolbox
%
% Entry points
%  psyrat_start - this launcher gui (buttons for processing, importing and
%   scoring, and viewing data)
%  psyrat_run - scripted (command-line) analyses using the same estimation
%   pipeline as the gui; documented in chapter 14 of the user manual, with a
%   quickstart in documentation/scripted_use.md
%
% Gui-related files
%  psyrat_startproc - gui to specify inputs for the analysis of data
%   (runs psyrat_loadfile and psyrat_computevarcomp)
%  psyrat_startimportscore - gui to import ERP files and score single-trial
%   data into a PsyRAT-ready table
%  psyrat_startview - gui to specify inputs for displaying processed data
%   (runs psyrat_relfigures)
%
% Data analysis files
%  psyrat_loadfile - function to load data files and prepare them for
%   processing
%  psyrat_computevarcomp - function to estimate variance components (via CmdStan
%   or the native fitlme/HMC engines) for dependability and generalizability
%   analyses
%  psyrat_relfigures - function to display information about dependability
%   (figures and tables)
%  psyrat_defaults - function to define default settings for processing and
%   viewing data
%
%Dependencies
% MatlabStan and MatlabProcessManager are bundled with the toolbox in
%  bundled_dependents/ and are added to the Matlab path automatically.
% CmdStan is NOT bundled and must be installed separately; run
%  psyrat_installdependents for guided installation, and see
%  documentation/dependencies_support_matrix.md for supported versions.
% CmdStan is required for the default estimation engine.
%

% Copyright (C) 2016-2025 Peter E. Clayson
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program (gpl.txt). If not, see
%     <http://www.gnu.org/licenses/>.
%

%If you have any questions or comments or if you run into any bugs, please
%report them on the GitHub issue tracker at
%https://github.com/peclayson/PsyRAT/issues. Including the PsyRAT version
%string (shown at startup) and, when available, the runconfig sidecar saved
%with the analysis makes problems much faster to reproduce and fix.

fprintf('\nEnsuring dependents are found in the Matlab path\n');
if exist('psyrat_startproc.m','file') ~= 2 || ...
        exist('psyrat_startview.m','file') ~= 2 || ...
        exist('psyrat_computevarcomp.m','file') ~= 2 || ...
        exist('psyrat_loadfile.m','file') ~= 2 || ...
        exist('psyrat_relfigures.m','file') ~= 2 || ...
        exist('psyrat_checkconv.m','file') ~= 2 || ...
        exist('psyrat_readtable.m','file') ~= 2 || ...
        exist('psyrat_reruncheck.m','file') ~= 2 || ...
        exist('psyrat_updatecheck.m','file') ~= 2 || ...
        exist('psyrat_checkgithubreleases.m','file') ~= 2 || ...
        exist('psyrat_defaults.m','file') ~= 2 || ...
        exist('psyrat_depvtrialsplot.m','file') ~= 2 || ...
        exist('psyrat_ptintervalplot.m','file') ~= 2 || ...
        exist('psyrat_installdependents.m','file') ~= 2 || ...
        exist('psyrat_checkversionsofdeps.m','file') ~= 2 || ...
        exist('psyrat_locate_cmdstan.m','file') ~= 2
    %psyrat_locate_cmdstan.m is in the sentinel list because the dependency
    %gate below calls it unguarded: a stale saved path from a pre-G43
    %checkout could satisfy the other sentinels while lacking the resolver,
    %and the gate would then die on an undefined function (adversarial
    %wave, 2026-08-29).

    %find where the PsyRAT files are located and add the directory and
    %sub-directories
    psyrat_path = which('psyrat_start');
    udir = fileparts(psyrat_path);
    %Add only the toolbox's KNOWN runtime locations to the path: the toolbox
    %root (entry points + runners), the 'subroutines' tree (all toolbox code),
    %and the bundled 'bundled_dependents' tree (MatlabStan). Do NOT genpath the
    %whole repository: genpath over the root also pulls in version-control
    %internals (.git), tooling folders (.claude), git worktrees, and throwaway
    %artifact folders nested under the toolbox. A stale or mock copy of a
    %function inside any of those (for example a test 'uiputfile' stub) would
    %then SHADOW the real file on the path, and adding '.git' subfolders also
    %prints spurious "nonexistent or not a directory" warnings at startup.
    %Test-only folders (e.g. 'tests') are added by the test runners themselves.
    addpath(udir);
    addpath(genpath(fullfile(udir,'subroutines')));
    addpath(genpath(fullfile(udir,'bundled_dependents')));
    psyrat_savepath_safe;
    fprintf('Directories for PsyRAT Toolbox files added to Matlab path\n');

else
    fprintf('PsyRAT Toolbox files found\n');
end

%look for the MatlabStan/MatlabProcessManager files on the Matlab path, and
%locate CmdStan through the single toolbox resolver (G43, owner ruling
%2026-08-29). CmdStan does NOT need to be on the Matlab path: an
%environment-variable (CMDSTAN/STAN_HOME/CMDSTAN_HOME) or toolbox-managed
%home install is found the same way estimation's own resolver
%(mstan.stan_home, via psyrat_locate_cmdstan) finds it, so those setups no
%longer hit this install dialog while estimation would have run fine. The
%companion version check (psyrat_checkversionsofdeps) consults the same
%resolver, so a passing gate is read by the same discovery and the old
%class of spurious unable-to-determine-version warnings (install present
%but invisible to a path-only lookup) is gone. A genuinely broken stanc
%under a valid-looking root can still warn there, correctly -- validity
%here is makefile + bin/, not a successful stanc run.
if exist('mcmc.m','file') ~= 2 || ...
        exist('stan.m','file') ~= 2 || ...
        exist('StanFit.m','file') ~= 2 || ...
        exist('StanModel.m','file') ~= 2 || ...
        exist('processManager.m','file') ~= 2 || ...
        exist('processState.m','file') ~= 2 || ...
        isempty(psyrat_locate_cmdstan)


    dlg = {'Warning: Dependencies for the PsyRAT toolbox were not found.';...
        'Please include the folders containing the scripts for';...
        'MatlabProcessManager and MatlabStan in your Matlab path, and';...
        'install CmdStan where the toolbox can find it (on the Matlab';...
        'path, in a CMDSTAN/STAN_HOME environment variable, or via the';...
        'guided installer).';...
        'Run psyrat_installdependents for guided setup, or see README.md';...
        'for manual installation instructions'};

    for i = 1:length(dlg)
        fprintf('%s\n',dlg{i});
    end
    fprintf('\n\n');
    psyrat_ask2install;
    return;

else

    fprintf('CmdStan, MatlabProcessManager, and MatlabStan files found\n');

    %check if the most up-to-date releases are being used for each of the
    %dependents

    depvercheck = psyrat_checkversionsofdeps;

    switch depvercheck.cmdstan
        case 0
            fprintf('There is a new version of CmdStan available\n');
        case 1
            fprintf('CmdStan version is current\n');
        case 2
            warning(strcat('Installed CmdStan version is newer than the',...
                ' version tested for the Toolbox.',...
                ' Toolbox may not work properly as a result'));
    end

    %warn (rather than block the GUI) if the installed CmdStan is below the
    %minimum PsyRAT can compile against. Estimation re-checks this and errors
    %before compiling, so the user still gets a clear message at run time.
    psyrat_require_min_cmdstan('2.26.0','warn');

    switch depvercheck.matlabprocessmanager
        case 0
            fprintf('There is a new version of MatlabProcessManager available\n');
        case 1
            fprintf('MatlabProcessManager version is current\n');
        case 2
            warning(strcat('Installed MatlabProcessManager version is newer than the',...
                ' version tested for the Toolbox.',...
                ' Toolbox may not work properly as a result'));
    end

    switch depvercheck.matlabstan
        case 0
            fprintf(['The MatlabStan in use is not the patched version ' ...
                'PsyRAT ships\n']);
        case 1
            fprintf('MatlabStan version is current\n');
        case 2
            warning(strcat('Installed MatlabStan version is newer than the',...
                ' version tested for the Toolbox.',...
                ' Toolbox may not work properly as a result'));
        case 3
            %G56d: the in-tree bundled copy with a would-be-old stamp. A
            %dependents update cannot replace the bundle (its removal is
            %refused by design), so the update dialog is not raised for it;
            %this line is the honest console signal that distinguishes the
            %pinned state from a genuinely current stamp.
            fprintf(['MatlabStan: PsyRAT''s bundled pinned copy is in use ' ...
                '(updated with the toolbox, not the dependents)\n']);
        otherwise
            %-1: no MatlabStan on the path, or a stamp that would not parse.
            fprintf('Could not determine the MatlabStan version in use\n');
    end

end

%pull version for the PsyRAT Toolbox
psyratver = psyrat_defineversion;

%Output info about PsyRAT Toolbox
fprintf('\n\n\nPsyRAT Version %s\n\n',psyratver);

%check whether running the newest release of the toolbox
if exist('psyrat_checkgithubreleases','file') == 2
    releaseinfo = psyrat_checkgithubreleases('currentversion',psyratver,...
        'verbose',true);
    psyratvercheck = releaseinfo.statusCode;
else
    psyratvercheck = psyrat_updatecheck(psyratver);
end

fprintf('\n');

%load default preferences for processing and viewing data
psyrat_prefs = psyrat_defaults;

%attach the current version number to psyrat_prefs
psyrat_prefs.ver = psyratver;

%create the gui for indicating whether the user wants to process data or view
%previously processed data

%define parameters for figure position
figwidth = 560;
figheight = 220;
psyrat_prefs.guis.fsize = psyrat_guifsize;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%initialize gui
psyrat_gui= psyrat_newfigure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'resize','off');

movegui(psyrat_gui,'center');

%Write text
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','What would you like to do?',...
    'Position',[0 row figwidth 25]);

%Create a button that will take the user to the gui for setting the inputs
%to process data
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',{'Process';'New Data'},...
    'Position', [20 25 160 95],...
    'TooltipString','Load a table of ERP measurements and run PsyRAT processing in CmdStan.',...
    'Callback',{@psyrat_startproc,psyrat_gui,'psyrat_prefs',psyrat_prefs});

%Create button that will start the ERP import + scoring workflow.
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',{'Import +';'Score ERP'},...
    'Position', [200 25 160 95],...
    'TooltipString','Import ERP waveform files, define component scoring, and export a validated handoff table.',...
    'Callback',{@psyrat_startimportscore,psyrat_gui,'psyrat_prefs',psyrat_prefs});

%Create button that will take the user to the gui for setting the inputs
%for viewing the data
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','View Results',...
    'Position', [380 25 160 95],...
    'TooltipString','Open an existing processed PsyRAT (.psyrat) file and view figures/tables.',...
    'Callback',{@psyrat_startview,psyrat_gui,'psyrat_prefs',psyrat_prefs});

%tag gui
psyrat_gui.Tag = 'psyrat_gui';

%display gui
set(psyrat_gui,'Visible','on');

%if applicable, ask the user whether PsyRAT Toolbox should be updated
if psyratvercheck == 0

    psyrat_ask2updatetoolbox;

%if applicable, ask the user whether dependents should be updated
elseif depvercheck.cmdstan == 0 || depvercheck.matlabstan == 0 || ...
        depvercheck.matlabprocessmanager == 0

    psyrat_ask2updatedeps(depvercheck);

end

end

function psyrat_ask2install
%gui to ask the user whether the PsyRAT Toolbox dependents should be installed

%define parameters for figure position
figwidth = 400;
figheight = 200;
fsize = psyrat_guifsize;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%initialize gui
psyrat_gui= psyrat_newfigure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'resize','off');

movegui(psyrat_gui,'center');

str = ['PsyRAT Toolbox dependents were not located. Would you like'...
    ' to install the dependents?'];

%Write text
uicontrol(psyrat_gui,'Style','text','fontsize',fsize+2,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 40]);

%Create a button that will install the dependents
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String',{'Install';'Dependents'},...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@installdeps});

%Create button to quit
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','Exit',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@giveup});

%tag gui
psyrat_gui.Tag = 'psyrat_gui';

%display gui
set(psyrat_gui,'Visible','on');

end

function psyrat_ask2updatetoolbox
%gui to ask the user whether the PsyRAT Toolbox dependents should be installed

%define parameters for figure position
figwidth = 400;
figheight = 200;
fsize = psyrat_guifsize;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%initialize gui
%
%B44, second half (ruled 2026-08-22): MODAL, same mechanism as the dependents
%prompt below. The justification here is NOT blast radius -- psyrat_updatepsyrat
%websaves, unzips and rmpaths, and never deletes the old toolbox directory --
%it is that a user can CONFUSE the two prompts. What they actually share,
%measured: the closing question "Would you like to do so now?", still
%verbatim in both; a 200 px figure height; button positions that are the
%identical figwidth-relative formula; the Tag psyrat_gui_update (invisible to
%a user, but it is why the findobj lookups need care); and -- until B45 --
%the same second sentence. What they do NOT share: the opening sentence, the
%WARNING clause (dependents only), and the text control size (400x40 vs
%500x75). This comment claimed they were "near-identical in wording, layout
%and Tag" until B50; that was the overstatement B45 had already refuted, and
%it survived here because nobody grepped for what cited it. Confusability is
%all the argument needs: leaving one prompt dismissable-by-burial and the
%other not is worse than either state on its own.
%
%The B20 precondition holds the same way: the main psyrat_gui already exists,
%and both callbacks close this figure before doing anything else (updatepsyrat
%runs psyrat_closepsyratfigs, which closes every 'psyrat_'-tagged window
%including this one).
%
%KNOWN COST, accepted with the ruling. This prompt fires whenever the release
%check returns 0 -- every launch AND every relaunch, and psyrat_start is
%invoked from EIGHTEEN sites across nine files including the viewer Back
%buttons (this said "ten" until 2026-08-22; the figure was never derived). A modal traps keyboard and mouse
%across MATLAB while it is visible, so a user several releases behind must
%now answer it every time instead of ignoring it. Being ignorable was the
%defect, so this is the intended direction, but it is a real increase in how
%often the user is interrupted and it lands hardest on the half whose reason
%is consistency rather than blast radius.
%
%psyrat_ask2install is deliberately NOT converted: it is raised at :172 and
%psyrat_start returns at :173 without ever building the home screen, so there
%is no window for a modal to bind and the B20 precondition does not hold.
psyrat_gui_update = psyrat_newfigure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'WindowStyle','modal',...
  'resize','off');

movegui(psyrat_gui_update,'center');

str = ['You are using an old version of the PsyRAT Toolbox.'...
    ' It is recommended that you update the toolbox.'...
    ' Would you like to do so now?'];

%Write text
uicontrol(psyrat_gui_update,'Style','text','fontsize',fsize+2,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 40]);

%Create a button that will install the dependents
uicontrol(psyrat_gui_update,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','Yes',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@updatepsyrat});

%Create button to quit
uicontrol(psyrat_gui_update,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','No',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@donotupdatepsyrat});

%tag gui
psyrat_gui_update.Tag = 'psyrat_gui_update';

%display gui
set(psyrat_gui_update,'Visible','on');

end

function psyrat_ask2updatedeps(depvercheck)
%gui to ask the user whether the PsyRAT Toolbox dependents should be installed

%define parameters for figure position
figwidth = 500;
figheight = 200;
fsize = psyrat_guifsize;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*4;

%initialize gui
%
%B44 (ruled 2026-08-22, after the click-through OBSERVED the defect): this
%prompt is MODAL. It confirms a destructive action -- Yes deletes the old
%dependency directories -- and non-modal it could be buried while still open
%and live. That was measured, not inferred: a figure census taken while the
%screen showed no dialog at all returned it as vis=on style=normal, sitting
%behind the MATLAB desktop window, with the main menu still clickable
%underneath a pending delete confirmation.
%
%MODALITY ALONE, NO uiwait -- and NOT for the reason it is tempting to write
%down. An earlier draft of this comment said nothing runs after this call
%because it is the last statement in psyrat_start. That is true only INSIDE
%this function, and it was refuted before commit: psyrat_start returns to its
%callers, and TWELVE call sites continue immediately with a modal dialog of
%their own -- each a bare "psyrat_start;" followed by a modal errordlg: the
%'File Error' cancel rejections in psyrat_startview.m:73,
%psyrat_startproc.m:83 and psyrat_updatepsyrat (three uiputfile-cancel arms
%plus, since the 2026-09-09 G87 ruling, its 'Update unavailable' refusal;
%line citations drifted twice and are not restated) and
%psyrat_installdependents (three uiputfile-cancel arms; the old 112/150/183
%line citation had drifted and is not restated), plus the three G56 'Update
%location not found' guards in psyrat_installdependents.
%(Corrected 2026-08-22 under B46: this said "two of them", naming only the
%two GUI sites. The undercount mattered because the closing clause below
%made a claim about those dialogs' modality while listing a fraction of
%them. That claim has since been measured and corrected -- see below -- and
%the count was re-derived from eight to eleven under G56.)
%
%uiwait is still wrong, for a different reason: it would hold psyrat_start
%open until this prompt is answered, so the caller's own "No data file was
%selected" message would be withheld behind an unrelated update prompt. B28
%needed uiwait because psyrat_startproc raced on into its measurement loop
%with nothing to stop it; here the follow-on dialogs are themselves modal and
%MATLAB gives focus to the most recently created modal, so the ordering
%resolves itself.
%
%MEASURED 2026-08-22, correcting what this comment used to assert. It said
%the resolution DEPENDS on those follow-on dialogs staying modal -- that a
%site which lost its 'modal' would be created on top of this prompt and be
%UNREACHABLE until this one is answered. The reachability half is WRONG.
%MATLAB modality binds only the windows that EXIST when the modal is created,
%so anything created afterwards escapes it. Probed on screen with a modal
%figure up (WindowStyle read from MATLAB, not from the accessibility layer):
%a plain non-modal figure and a bare errordlg (third argument omitted, which
%is exactly the failure mode described above) each took focus -- gcf moved to
%the new window -- and fired on the FIRST click, with a button on the modal
%itself as the positive control to prove clicks were being delivered at all.
%What was measured is FOCUS and INTERACTIVITY, not occlusion: the two probe
%figures did not overlap, so "created on top" is neither confirmed nor
%contradicted here, and it is not the half that mattered. The ordering
%resolves itself whether or not the follow-on dialogs keep their 'modal'.
%
%This is the same rule psyrat_startview.m:69-72 states. The two comments used
%to contradict each other; this one was the wrong half.
%
%Creation order is already correct for B20: the main psyrat_gui is created
%and made visible before this figure, so the modal binds it. Both callbacks
%close this figure before doing anything else (updatedeps runs
%psyrat_closepsyratfigs first), so no modal survives into the installer.
%THAT GUARANTEE IS WHY THE TEARDOWN IS PREFIX-SCOPED AND NOT AN EXPLICIT
%LIST: this figure's Tag, psyrat_gui_update, must stay inside whatever
%psyrat_closepsyratfigs matches, or MATLAB is left trapped behind a modal
%for the whole of the installer run (B48).
psyrat_gui_update = psyrat_newfigure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'WindowStyle','modal',...
  'resize','off');

movegui(psyrat_gui_update,'center');

deps2update = '';

%These are three separate tests, not an if/elseif chain. rmolddeps removes
%EVERY dependent flagged 0, so naming only the first one understated what was
%about to be deleted: a dialog that said "CmdStan" could remove MatlabStan and
%MatlabProcessManager too.
outofdate = {};
if depvercheck.cmdstan == 0
    outofdate{end+1} = 'CmdStan';
end
if depvercheck.matlabstan == 0
    outofdate{end+1} = 'MatlabStan';
end
if depvercheck.matlabprocessmanager == 0
    outofdate{end+1} = 'MatlabProcessManager';
end

if ~isempty(outofdate)
    deps2update = [' ' strjoin(outofdate,', ')];
end

%B45 (ruled 2026-08-22; found on screen during the B42 click-through): this
%sentence names THE DEPENDENTS, and it is deliberately NOT shared with the
%sibling psyrat_ask2updatetoolbox above, where "the toolbox" is correct. The
%two prompts shared this SENTENCE until B45 and still share their Tag, so the
%sentence itself did not tell the user which prompt it belonged to -- and the
%one it was wrong on is the one whose Yes DELETES the dependency directories.
%They are NOT otherwise identical, and claiming so would overstate the defect:
%the opening sentence differs, only this prompt carries the WARNING clause,
%and the text control is 500x75 HERE vs 400x40 in the sibling.
%
%"the dependents" rather than "them" because it reads the same in the
%one-dependent and three-dependent renderings, and it reuses the noun the
%WARNING sentence below already uses. Pinned PER FUNCTION in
%TestModalRejectionDialogs, since a file-level count cannot tell the two
%prompts apart.
str = ['You are using old version(s) of' deps2update...
    '. It is recommended that you update the dependents.'...
    ' Would you like to do so now? WARNING: Doing so will delete the '...
    'old directories for the dependents to avoid confusion.'];

%Write text
uicontrol(psyrat_gui_update,'Style','text','fontsize',fsize+2,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 75]);

%Create a button that will install the dependents
uicontrol(psyrat_gui_update,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','Yes',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@updatedeps,depvercheck});

%Create button to quit
uicontrol(psyrat_gui_update,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','No',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@donotupdatedeps});

%tag gui
psyrat_gui_update.Tag = 'psyrat_gui_update';

%display gui
set(psyrat_gui_update,'Visible','on');

end

function donotupdatepsyrat(varargin)
%if the user does not what the PsyRAT toolbox updated, just continue on
%check if psyrat_gui_update is open.
psyrat_gui_update = findobj('Tag','psyrat_gui_update');

if ~isempty(psyrat_gui_update)
    close(psyrat_gui_update);
end

fprintf('\nPsyRAT Toolbox will not be updated\n\n');

end

function donotupdatedeps(varargin)
%if the user does not what the PsyRAT toolbox updated, just continue on
%check if psyrat_gui_update is open.
psyrat_gui_update = findobj('Tag','psyrat_gui_update');

if ~isempty(psyrat_gui_update)
    close(psyrat_gui_update);
end

fprintf('\nDependents will not be updated\n\n');

end

function updatepsyrat(varargin)
%if the user wants the toolbox updated, begin updating
psyrat_closepsyratfigs;
psyrat_updatepsyrat;
end

function updatedeps(varargin)
%if the user wants the dependents updated, begin updating
psyrat_closepsyratfigs;
depvercheck = varargin{3};
psyrat_installdependents('depvercheck',depvercheck);
end

function giveup(varargin)
%if the user does not want the dependents installed, close PsyRAT's own
%windows and stop. This button is the "no thanks" route, so it must not
%destroy anything the user owns (B48).
psyrat_closepsyratfigs;
end

function installdeps(varargin)
%if the user wants the dependents installed, begin the installation
psyrat_closepsyratfigs;
psyrat_installdependents;
end

function psyrat_closepsyratfigs
%Tear down PsyRAT's own windows and leave every other figure alone (B48).
%
%The four startup callbacks above used to call a bare `close all`, which
%closes EVERY figure in the session whose HandleVisibility is on, whoever
%owns it. All four fire during startup prompts -- exactly when a user is
%likely to have EEGLAB/ERPLAB figures open from earlier work -- so accepting
%an update, or even declining the install prompt via giveup, silently
%destroyed their unsaved figures.
%
%Scoped by Tag PREFIX rather than by an explicit list, because the set of
%PsyRAT windows open at these moments is not fixed: psyrat_start is invoked
%from EIGHTEEN sites across nine files (counted, not inherited -- an earlier
%draft of this comment said "ten", copied from the note at
%psyrat_ask2updatetoolbox, and that figure is wrong), including the viewer
%Back buttons, so psyrat_output tables and plots, psyrat_import_gui,
%psyrat_gui_trialbrowser and the selector guis can all still be up.
%TestModalRejectionDialogs pins the Tags of the windows that were untagged
%when this was written, so they cannot be reintroduced untagged and silently
%stranded here.
%
%A KNOWN AND ACCEPTED NARROWING -- READ THIS BEFORE "FIXING" IT. `close all`
%also closed PsyRAT's OWN msgbox/warndlg/errordlg windows; this does not.
%MATLAB tags those 'Msgbox_<title>' with HandleVisibility 'callback', and
%these four sites ARE uicontrol callbacks, so in a genuine callback context
%they were inside what `close all` swept. MEASURED in R2025b, not argued: in
%a real callback context `get(0,'Children')` returned the errordlg and
%`close all` destroyed it, while this selector does not match it. (A probe
%that merely fevals a Callback property does NOT reproduce that -- gcbo is
%empty there and the dialog stays hidden, which is how the first measurement
%got it backwards.)
%
%IT IS ACCEPTED BECAUSE THE ALTERNATIVE REINTRODUCES THE BUG. 'Msgbox_' is
%MATLAB's own prefix, not PsyRAT's -- EEGLAB and ERPLAB dialogs carry it too
%-- so matching it here would once again close windows PsyRAT does not own,
%which is exactly what B48 exists to stop. The residue is a lingering
%NON-MODAL PsyRAT advisory outliving the teardown. A modal one cannot be in
%that state: MATLAB modality is creation-ordered, so a dialog raised after
%these prompts blocks them and must be dismissed before either button can be
%clicked at all.
%
%LOAD-BEARING, both directions. It must keep closing psyrat_gui_update: that
%figure is WindowStyle 'modal' (B44), so leaving it up would trap keyboard
%and mouse across MATLAB while the installer runs. It must also keep closing
%psyrat_instgui, because psyrat_guistatus (in psyrat_installdependents)
%REUSES an existing one by Tag and draws fresh uicontrols into it without
%clearing, so a stale one surviving would stack new status text on old.
%Named rather than cited by line: this file's own edits move those numbers.
%
%findall (not findobj) so hidden-handle PsyRAT figures are reached too, and
%'Type','figure' so only windows match -- this is the same idiom the viewers
%already use for their "Close PsyRAT Outputs" buttons.
psyratfigs = findall(0,'Type','figure','-regexp','Tag','^psyrat_');

if ~isempty(psyratfigs)
    close(psyratfigs);
end

end
