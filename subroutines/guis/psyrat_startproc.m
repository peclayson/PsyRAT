function psyrat_startproc(varargin)
%Initiate Matlab gui to begin processing data in Stan
%
%
%
%
%Input
% There are no required inputs to execute this script. The user will be
%  asked to identify the location of the data to processed and the columns
%  of interest for analysis
%
%Output
% This script will not output any variables to the workspace. A .mat with
%  processed data from Stan will be saved in a user specified locaiton.
%  This .mat will automatically be loaded into psyrat_startview for viewing
%  the data
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

psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    %Grab preferences if they exist
    if ~isempty(varargin)

        %check for psyrat_prefs
        ind = find(strcmp('psyrat_prefs',varargin),1);
        if ~isempty(ind)
            psyrat_prefs = varargin{ind+1};
        else
            %load default preferences for processing and viewing data
            psyrat_prefs = psyrat_defaults;

            %attach the current version number to psyrat_prefs
            psyrat_prefs.ver = psyrat_defineversion;

            %define parameters for figure position
            psyrat_prefs.guis.fsize = psyrat_guifsize;
        end
    end

    close(psyrat_gui);

elseif isempty(psyrat_gui)

    %load default preferences for processing and viewing data
    psyrat_prefs = psyrat_defaults;

    %attach the current version number to psyrat_prefs
    psyrat_prefs.ver = psyrat_defineversion;

    %define parameters for figure position
    psyrat_prefs.guis.fsize = psyrat_guifsize;

end

%ask the user to identify the data file to be loaded
[filepart, pathpart] = uigetfile({'*.xlsx;*.xls;*.csv;*.dat;*.txt;*.ods',...
    'All Readable Files'; '*.xlsx','Excel File (.xlsx)';'*.xls',...
    'Excel File 97-2003 (.xls)';'*.csv',...
    'Comma-Separated Value File (.csv)';'*.dat',...
    'Tab-Delimited Text (.dat)';'*.txt',...
    'Raw Text File (.txt)';'*.ods',...
    'OpenDocument Spreadsheet (.ods)'},'Data');

%if the user does not select a file, then take the user back to psyrat_start
if filepart == 0
    psyrat_start;
    errordlg({'No data file was selected.';...
        'How to fix: click "Process New Data" again and choose an input file.'},...
        'File Error','modal');
    return;
end

fprintf('\n\nLoading Data...');
fprintf('\nThis may take awhile depending on the amount of data...\n\n');

%load file
dataraw = psyrat_readtable('file',fullfile(pathpart,filepart));

%pull the headernames from the file
collist = dataraw.Properties.VariableNames;

%A "not used" entry is added so it will be presented as an option to the user
%(e.g., if there are no groups in the dataset then groups will be ignored by
%selecting it from a drop-down menu).
%
%The label is '(none)', not 'none', and the parentheses are load-bearing.
%psyrat_readtable calls readtable WITHOUT 'VariableNamingRule','preserve', so
%MATLAB coerces every header to a valid identifier -- and '(none)' is not one.
%A data column can therefore NEVER carry this label, which makes the sentinel
%unambiguous in the popup by construction rather than by convention. Under the
%old 'none' label a file with a column literally named none showed the word
%twice, indistinguishably (B25; observed live 2026-08-17).
%
%This is DISPLAY ONLY. Nothing reads the shown string back: the GUI reads the
%popup's 'value' (an index) and 'UserData', never 'String', and every
%facet-absent guard keys on the sentinel INDEX
%(selection == length(proc.collist)). Selecting it writes '' downstream, which
%is what psyrat_loadfile and psyrat_write_runconfig consume. Changing the label
%is therefore safe; changing the INDEX convention is not.
%
%Cost, accepted by owner ruling 2026-08-18: documentation/UserManual.pdf still
%says to select 'none' in three places and has no source in this repo to
%regenerate from. Tracked as B32.
vcollist = collist;
vcollist{end+1} = '(none)';

%psyrat_prefs may have been defined in the function call not using the gui or
%after the chains did not converge
if ~exist('psyrat_prefs','var')
    [psyrat_prefs,~] = psyrat_findprefsdata(varargin);
elseif isempty(psyrat_prefs)
    %load default preferences for processing and viewing data
    psyrat_prefs = psyrat_defaults;

    %attach the current version number to psyrat_prefs
    psyrat_prefs.ver = psyrat_defineversion;

    %define parameters for figure position
    psyrat_prefs.guis.fsize = psyrat_guifsize;
end

%Insert the information into a data structure for holding the PsyRAT data
psyrat_data = struct;
psyrat_data.ver = psyrat_prefs.ver;
psyrat_data.raw.filename = filepart;
psyrat_data.raw.filepath = pathpart;
%full input path for the sidecar's file hash (raw.filename is the bare name the
%screens display; see psyrat_write_runconfig, audit 2026-09-05)
psyrat_data.raw.sourcepath = fullfile(pathpart,filepart);
psyrat_data.raw.data = dataraw;
psyrat_data.raw.colnames = collist;
psyrat_data.proc.collist = vcollist;

%use the loaded dateset in the gui for setting up the data to be run
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

end

function psyrat_startproc_gui(varargin)
%
%Inputs
% psyrat_prefs - preferences for PsyRAT toolbox
% psyrat_data - PsyRAT toolbox data structure
%
%Optional Inputs:
% procprefs - prefences for processing in stan
%  .nchains - number of chains
%  .nwarmup - number of warmup iterations
%  .nsampling - number of post-warmup sampling iterations
%  .niter - total number of iterations (legacy alias)
% inpchoices - choices for the assignment of table columns to data inputs
%  (e.g., id, measurement). Used when the preferences executes psyrat_gui.
%
%Output
% There are no direct outputs to the Matlab workspace. The user will have
%  the option (via the gui) to go back to the psyrat_start window or pass the
%  data to psyrat_startproc to process the data in Stan
%

%somersault through varargin inputs to check for psyrat_prefs and psyrat_data
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%define parameters for figure position
figwidth = 600;
%taller than the legacy 470 to fit the two optional dimension-column pickers
%plus the residual-correlation parameterization picker (dynamic reliability)
%and the two reliability-with-splits pickers (items-per-split column + split
%type), all added after the occasion row.
figheight = 645;

%if the assignment of which columns belong to which category has already
%been made, then use those choices. Otherwise, load the defaults and try to
%figure out the data.
if ~isfield(psyrat_prefs.proc,'inp')

    %help the user set things up by checking if any generic names are present
    %in the headers that identify the columns
    psyrat_prefs.proc.inp.id = 1;
    psyrat_prefs.proc.inp.meas = 2;
    psyrat_prefs.proc.inp.group = length(psyrat_data.proc.collist);
    psyrat_prefs.proc.inp.event = length(psyrat_data.proc.collist);
    psyrat_prefs.proc.inp.time = length(psyrat_data.proc.collist);
    %dimension columns for the dynamic-reliability analysis (Rast & Clayson);
    %default 'none'. Selecting a Dimension 1 column enables the analysis.
    psyrat_prefs.proc.inp.dim1 = length(psyrat_data.proc.collist);
    psyrat_prefs.proc.inp.dim2 = length(psyrat_data.proc.collist);
    %items-per-split (n_i) column for reliability with data splits; default
    %'none'. Selecting a column declares that each row is a split mean.
    psyrat_prefs.proc.inp.weight = length(psyrat_data.proc.collist);
    psyrat_prefs.proc.inp.whichgroups = '';
    psyrat_prefs.proc.inp.whichgroupscol = [];
    psyrat_prefs.proc.inp.whichevents = '';
    psyrat_prefs.proc.inp.whicheventscol = [];
    psyrat_prefs.proc.inp.whichtimes = '';
    psyrat_prefs.proc.inp.whichtimescol = [];

    %check for a participant header
    %
    %ANCHOR -- this loop and the measurement loop below walk proc.collist,
    %which is raw.colnames plus the trailing sentinel, so they iterate ncol+1
    %entries. inp.id and inp.meas, though, are consumed against raw.colnames,
    %which has only ncol: the Participant ID popup further down is built with
    %'String',psyrat_data.raw.colnames. Entries 1:ncol are identical in both
    %lists (the index-space note below the measurement fallback works this
    %through), so the extra iteration is harmless ONLY because no candidate
    %list in this block contains the sentinel's label. If one ever did, the
    %match would set the role index to ncol+1 and the popup's Value would
    %exceed its String. TestStartprocNoneSentinel pins that precondition.
    %
    %The group/event/time loops below carry no such coupling -- their indices
    %address proc.collist, which really does have ncol+1 entries, and the
    %sentinel index is exactly how they encode "facet not used" (B25).
    poss = {'Subject' 'ID' 'Participant' 'SubjID' 'Subj'};
    for i = 1:length(psyrat_data.proc.collist)
        ind = strcmpi(psyrat_data.proc.collist(i),poss);
        if any(ind)
            psyrat_prefs.proc.inp.id = i;
        end
    end

    %check for a mesurement header
    %'Meas' is the toolbox's own canonical column name (the I/O contract, and
    %psyrat_run's meascol default), so a file already written to the contract
    %auto-detected nothing before it was listed here.
    %measfound only RECORDS whether a header matched. The candidate list and the
    %comparison are unchanged; the flag exists so the fallback further down knows
    %whether the seeded positional default is still in place.
    %Same ncol+1 index-space coupling as the participant loop above: this
    %walks the sentinel-extended list, but inp.meas addresses raw.colnames.
    measfound = false;
    poss = {'Measurement' 'Meas'};
    for i = 1:length(psyrat_data.proc.collist)
        ind = strcmpi(psyrat_data.proc.collist(i),poss);
        if any(ind)
            psyrat_prefs.proc.inp.meas = i;
            measfound = true;
        end
    end

    %check for a group header
    poss = {'Group'};
    for i = 1:length(psyrat_data.proc.collist)
        ind = strcmpi(psyrat_data.proc.collist(i),poss);
        if any(ind)
            psyrat_prefs.proc.inp.group = i;
        end
    end

    %check for an event type header
    poss = {'Event' 'Type'};
    for i = 1:length(psyrat_data.proc.collist)
        ind = strcmpi(psyrat_data.proc.collist(i),poss);
        if any(ind)
            psyrat_prefs.proc.inp.event = i;
        end
    end

    %check for a occasion header
    %the comparison below is strcmpi (exact, case-insensitive), so the leading
    %space this list used to carry on ' Session' meant a column named Session
    %could never match.
    poss = {'Time' 'Occasion' 'Session'};
    for i = 1:length(psyrat_data.proc.collist)
        ind = strcmpi(psyrat_data.proc.collist(i),poss);
        if any(ind)
            psyrat_prefs.proc.inp.time = i;
        end
    end

    %fallback for the measurement column when no header matched the candidate
    %names above. It sits here, after the group/event/time loops, because it
    %needs the FINAL role assignments.
    %
    %The seeded positional default (column 2) is a guess about long-format
    %layout, and a dangerous one. In the toolbox's own example file
    %(test_data/testdata.csv, headers subjid,group,event,ern) the score column
    %'ern' matches neither candidate, so Measurement stayed on column 2 -- the
    %GROUP column -- with two bad outcomes. (a) Group was itself auto-detected,
    %so one column sat in two roles and the run died at the Run button on the
    %measurement/facet overlap dialog below; the shipped example was unrunnable
    %as preselected. (b) Had the facet column been numerically coded and matched
    %no literal (say 'grp' holding 1/2), Group would have stayed 'none', nothing
    %would have overlapped, the numeric-measurement guard below would have
    %passed, and the run would have silently reported the reliability of a group
    %code.
    %
    %So prefer the last NUMERIC column that no other role already claimed. In
    %long-format ERP exports the score is the dependent variable and comes last,
    %and a column already serving as ID/Group/Event/Occasion can never also be
    %the measurement, which removes both failure modes by construction.
    %
    %This is a heuristic, not a universal improvement. It costs accuracy on a
    %file whose score is NOT last (say subjid,ern,trialnum): the preselection
    %moves off 'ern' onto 'trialnum'. But that is still a numeric column no role
    %claims, which the user retargets in one click, whereas the two failures
    %above were a dead end and a silent wrong answer. "Last" is the tie-break
    %because in long-format data the score is the melted value column.
    %
    %Mind the two index spaces: inp.meas indexes raw.colnames, while
    %inp.group/event/time index proc.collist, whose extra final entry is the
    %'none' sentinel. Entries 1:ncol are identical in both, so a role index above
    %ncol is 'none' and claims no column. Every facet-absent guard in this file
    %keys on that INDEX, never on the selected column's NAME, so a real data
    %column literally named 'none' stays selectable and legitimate (B25).
    %
    %If nothing qualifies, the positional default is deliberately left alone:
    %there is then no defensible measurement column in the file, and the two
    %Run-button guards give the user a specific dialog. Guessing would hide that.
    %
    %Leaving inp.meas unset is not an option: psyrat_coercemeasureindices maps
    %empty to 1, so an unset value silently preselects the ID column instead.
    %
    %This changes a preselection only. The user still confirms the Measurement
    %popup before running, and a stored column mapping skips this whole block.
    if ~measfound && isfield(psyrat_data,'raw') && ...
            isfield(psyrat_data.raw,'colnames') && ...
            isfield(psyrat_data.raw,'data') && istable(psyrat_data.raw.data)
        ncol = length(psyrat_data.raw.colnames);
        roleinds = [psyrat_prefs.proc.inp.id psyrat_prefs.proc.inp.group ...
            psyrat_prefs.proc.inp.event psyrat_prefs.proc.inp.time];
        claimed = roleinds(roleinds >= 1 & roleinds <= ncol);
        for i = ncol:-1:1
            if any(i == claimed)
                continue;
            end
            if ~isnumeric(psyrat_data.raw.data.(psyrat_data.raw.colnames{i}))
                continue;
            end
            psyrat_prefs.proc.inp.meas = i;
            break;
        end
    end
end

%initialize DoD contrast mapping fields if they do not exist
if ~isfield(psyrat_prefs.proc,'dodmap') || isempty(psyrat_prefs.proc.dodmap)
    psyrat_prefs.proc.dodmap = {};
end
if ~isfield(psyrat_prefs.proc,'dodmapcol')
    psyrat_prefs.proc.dodmapcol = [];
end
if ~isfield(psyrat_prefs.proc,'dodset') || isempty(psyrat_prefs.proc.dodset)
    psyrat_prefs.proc.dodset = 0;
end
if ~isfield(psyrat_prefs.proc,'diffest') || isempty(psyrat_prefs.proc.diffest)
    psyrat_prefs.proc.diffest = 1;
end

%normalize stored input selections so older preference files still open
%cleanly with the new multi-measure workflow.
psyrat_prefs.proc.inp.id = psyrat_coercescalarindex(psyrat_prefs.proc.inp.id,...
    length(psyrat_data.raw.colnames),1);
psyrat_prefs.proc.inp.meas = psyrat_coercemeasureindices(psyrat_prefs.proc.inp.meas,...
    length(psyrat_data.raw.colnames));
psyrat_prefs.proc.inp.group = psyrat_coercescalarindex(psyrat_prefs.proc.inp.group,...
    length(psyrat_data.proc.collist),length(psyrat_data.proc.collist));
psyrat_prefs.proc.inp.event = psyrat_coercescalarindex(psyrat_prefs.proc.inp.event,...
    length(psyrat_data.proc.collist),length(psyrat_data.proc.collist));
psyrat_prefs.proc.inp.time = psyrat_coercescalarindex(psyrat_prefs.proc.inp.time,...
    length(psyrat_data.proc.collist),length(psyrat_data.proc.collist));
%dimension selections (back-compat default 'none' for older preference files)
if ~isfield(psyrat_prefs.proc.inp,'dim1'); psyrat_prefs.proc.inp.dim1 = length(psyrat_data.proc.collist); end
if ~isfield(psyrat_prefs.proc.inp,'dim2'); psyrat_prefs.proc.inp.dim2 = length(psyrat_data.proc.collist); end
psyrat_prefs.proc.inp.dim1 = psyrat_coercescalarindex(psyrat_prefs.proc.inp.dim1,...
    length(psyrat_data.proc.collist),length(psyrat_data.proc.collist));
psyrat_prefs.proc.inp.dim2 = psyrat_coercescalarindex(psyrat_prefs.proc.inp.dim2,...
    length(psyrat_data.proc.collist),length(psyrat_data.proc.collist));
%items-per-split (n_i) column selection for reliability with data splits
%(back-compat default 'none' for older preference files).
if ~isfield(psyrat_prefs.proc.inp,'weight'); psyrat_prefs.proc.inp.weight = length(psyrat_data.proc.collist); end
psyrat_prefs.proc.inp.weight = psyrat_coercescalarindex(psyrat_prefs.proc.inp.weight,...
    length(psyrat_data.proc.collist),length(psyrat_data.proc.collist));
%split type for reliability with data splits: 1 = parallel (equal n_i, default),
%2 = nonparallel (unequal n_i). Inert unless an n_i column is selected. Back-compat
%default for older preference files, clamped to a valid 1/2 mode.
if ~isfield(psyrat_prefs.proc,'splitsmode') || isempty(psyrat_prefs.proc.splitsmode) || ...
        ~any(psyrat_prefs.proc.splitsmode == [1 2])
    psyrat_prefs.proc.splitsmode = 1;
end
%residual-correlation parameterization for the subject-level concurrent dynamic
%difference variants (cases 21/22): 1 = Population (shared across subjects,
%default), 2 = Per-subject. Inert for every other run; the engine combination
%guard is the backstop. Back-compat default for older preference files, and
%clamp to a valid 1/2 mode so a stray value cannot break the popup Value.
if ~isfield(psyrat_prefs.proc,'ssrescor') || isempty(psyrat_prefs.proc.ssrescor) || ...
        ~any(psyrat_prefs.proc.ssrescor == [1 2])
    psyrat_prefs.proc.ssrescor = 1;
end
%estimation engine preference: 0 = Automatic (CmdStan, else native fallback,
%default), 1 = CmdStan only, 2 = native fitlme, 3 = native HMC. Drives
%psyrat_resolve_engine via the warp (which reads proc.engine). Back-compat
%default for older preference files, clamped to a valid 0-3 code so a stray value
%cannot break the popup Value mapping.
if ~isfield(psyrat_prefs.proc,'engine') || isempty(psyrat_prefs.proc.engine) || ...
        ~any(psyrat_prefs.proc.engine == [0 1 2 3])
    psyrat_prefs.proc.engine = 0;
end

%define space between rows and first row location
rowspace = 35;
row = figheight - rowspace*2;

%define locations of column 1 and 2 for the gui
lcol = 30;
rcol = (figwidth/2);

%allow room for subset action buttons and summary text
subset_pop_w = 170;
subset_btn_w = 120;
subset_btn_x = rcol + subset_pop_w + 5;

%create the basic psyrat_gui
psyrat_gui= psyrat_newfigure('unit','pix',...
    'position',[400 400 figwidth figheight],...
    'menub','no',...
    'name','Specify Inputs',...
    'numbertitle','off',...
    'resize','off');

%Print the name of the loaded dataset
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' psyrat_data.raw.filename],...
    'Position',[0 row figwidth 25]);

%next row
row = row - rowspace*1.5;

%Print the gui headers
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Variable',...
    'Position', [figwidth/8 row figwidth/4 25]);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Input',...
    'Position',[5*(figwidth/8) row figwidth/4 25]);

%next row
row = row - rowspace*1.5;

%print the name of the variables and place a listbox with possible options
%start with participant ID row
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Participant ID:',...
    'Position', [lcol row figwidth/4 25]);

inplists(1) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_data.raw.colnames,'Value',psyrat_prefs.proc.inp.id,...
    'Position', [rcol row subset_pop_w 25]);

%next row
row = row - rowspace;

%measurement row
measrow = row;
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Measurement:',...
    'Position', [lcol row figwidth/4 25]);

[meas_pop_str,meas_pop_val,meas_pop_userdata] = ...
    psyrat_measurepopupdisplay(psyrat_data,psyrat_prefs);
inplists(2) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',meas_pop_str,'Value',meas_pop_val,'UserData',meas_pop_userdata,...
    'Position', [rcol row subset_pop_w 25]);

%next row
row = row - rowspace;
grouprow = row;

%group row
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Group ID:',...
    'Position', [lcol row figwidth/4 25]);

inplists(3) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_data.proc.collist,'Value',psyrat_prefs.proc.inp.group,...
    'Position', [rcol row subset_pop_w 25]);

%next row
row = row - rowspace;
eventrow = row;

%event row
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Event Type:',...
    'Position', [lcol row figwidth/4 25]);

inplists(4) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_data.proc.collist,'Value',psyrat_prefs.proc.inp.event,...
    'Position', [rcol row subset_pop_w 25]);

%next row
row = row - rowspace;
timerow = row;

%time row
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Occasion ID:',...
    'Position', [lcol row figwidth/4 25]);

inplists(5) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_data.proc.collist,'Value',psyrat_prefs.proc.inp.time,...
    'Position', [rcol row subset_pop_w 25]);

%next row: dimension columns for dynamic/conditional reliability (Rast &
%Clayson). These are continuous, between-person covariates. Selecting a
%Dimension 1 column enables the dynamic-reliability analysis; Dimension 2 is
%optional (adds a second predictor and its interaction). Leave as 'none' for a
%standard analysis.
row = row - rowspace;
dim_help = ['Optional. Continuous, between-person covariate (one value per ',...
    'participant, e.g., a questionnaire total). Selecting Dimension 1 runs ',...
    'the Rast & Clayson dynamic-reliability model, reporting reliability as a ',...
    'function of the standardized dimension(s). Works with single-session or ',...
    'test-retest data, and with non-difference or two-event difference scores ',...
    '(including subject-level error variances).'];
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Dimension 1 (dynamic):','Tooltip',dim_help,...
    'Position', [lcol row figwidth/4 25]);

inplists(6) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_data.proc.collist,'Value',psyrat_prefs.proc.inp.dim1,...
    'Tooltip',dim_help,'Position', [rcol row subset_pop_w 25]);

%next row
row = row - rowspace;
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Dimension 2 (dynamic):','Tooltip',dim_help,...
    'Position', [lcol row figwidth/4 25]);

inplists(7) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_data.proc.collist,'Value',psyrat_prefs.proc.inp.dim2,...
    'Tooltip',dim_help,'Position', [rcol row subset_pop_w 25]);

%next row: residual-correlation parameterization for the subject-level concurrent
%dynamic difference variants (cases 21/22). Population shares one residual
%correlation across participants (default); Per-subject estimates a separate
%residual correlation per participant. Only meaningful for a subject-level
%concurrent difference (difference scores + subject-specific error variances +
%within-person residual covariance, all set in Processing Preferences); enabled
%only for that combination and disabled otherwise. The Dimension-1 requirement
%(dynrel = 2) is the user's final selection on this window and is enforced at run
%time by the engine combination guard ('varargin:ssrescor'), which ignores this
%setting for any non-dynamic-subject-level-concurrent run. Gating on the three
%Processing Preferences (which redraw this window on save) rather than on the
%Dimension-1 popup (which does not redraw on change) keeps Per-subject reachable.
row = row - rowspace;
ssrescor_help = ['Residual-correlation parameterization for the subject-level ',...
    'concurrent dynamic difference variants. Population shares one residual ',...
    'correlation across participants; Per-subject estimates a separate residual ',...
    'correlation per participant. Takes effect only for a dynamic (Dimension 1 ',...
    'selected) subject-level concurrent difference run (difference scores + ',...
    'subject-specific error variances + within-person residual covariance).'];
%this gate keys on diffwpcov==2 (the user-facing "within-person covariance"
%preference), which the pipeline keeps equivalent to diffrescor==2: psyrat_procprefs
%and psyrat_prefs_save both run psyrat_canonicalize_proc, which forces
%diffrescor = (diffwpcov==2 ? 2 : 1). The warp's canonicalizer call and the
%engine ssrescor guard key on diffrescor==2, so all three checks describe the
%same "subject-level concurrent difference" combination. Keep that coupling in mind
%if diffwpcov/diffrescor are ever decoupled.
ssrescor_applies = ...
    isfield(psyrat_prefs.proc,'diffest') && psyrat_prefs.proc.diffest == 2 && ...
    isfield(psyrat_prefs.proc,'ssubjrel') && psyrat_prefs.proc.ssubjrel == 2 && ...
    isfield(psyrat_prefs.proc,'diffwpcov') && psyrat_prefs.proc.diffwpcov == 2;
if ssrescor_applies; ssrescor_enable = 'on'; else; ssrescor_enable = 'off'; end
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Residual correlation:','Tooltip',ssrescor_help,...
    'Position', [lcol row figwidth/4 25]);

inplists(8) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',{'Population (shared)','Per-subject'},'Value',psyrat_prefs.proc.ssrescor,...
    'Tooltip',ssrescor_help,'Enable',ssrescor_enable,...
    'Position', [rcol row subset_pop_w 25]);

%next row: reliability with data splits (Rocha Tables 4-5). Selecting an
%items-per-split (n_i) column declares that each row is a split mean (the
%number of items/trials averaged into that split) rather than a single trial,
%enabling the splits analysis. Leave as 'none' for a standard single-trial
%analysis. The Split type popup below declares whether the splits are parallel
%(equal n_i) or nonparallel (unequal n_i).
row = row - rowspace;
weight_help = ['Optional. Number of items/trials averaged into each split ',...
    'mean (n_i). Selecting this column declares that each row is a split mean ',...
    'rather than a single trial, running the reliability-with-splits analysis. ',...
    'Leave as none for a standard single-trial analysis.'];
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Items per split (n_i):','Tooltip',weight_help,...
    'Position', [lcol row figwidth/4 25]);

inplists(9) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_data.proc.collist,'Value',psyrat_prefs.proc.inp.weight,...
    'Tooltip',weight_help,'Position', [rcol row subset_pop_w 25]);

%next row: split type (parallel vs nonparallel). Only consulted when an n_i
%column is selected above; ignored for a standard single-trial analysis.
%Parallel (equal n_i) reuses the standard one-facet/test-retest machinery (the
%split mean is the score, n' = number of splits). Nonparallel (unequal n_i)
%fits the weighted observed-design model, which is identified only by n_i
%variation across splits.
row = row - rowspace;
splittype_help = ['Declares whether the data splits are parallel (equal n_i ',...
    'across splits) or nonparallel (unequal n_i). Only applies when an Items ',...
    'per split (n_i) column is selected above. Parallel reuses the standard ',...
    'one-facet/test-retest model; nonparallel fits the weighted observed-design ',...
    'model and requires variation in n_i across splits.'];
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Split type:','Tooltip',splittype_help,...
    'Position', [lcol row figwidth/4 25]);

inplists(10) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',{'Parallel (equal n_i)','Nonparallel (unequal n_i)'},...
    'Value',psyrat_prefs.proc.splitsmode,...
    'Tooltip',splittype_help,'Position', [rcol row subset_pop_w 25]);

%next row: estimation engine. Default Automatic uses CmdStan when installed and
%otherwise falls back to the first native engine that supports the design (native
%HMC preferred over native fitlme). The other choices FORCE an engine. The native
%options are annotated with a "(toolbox not installed)" hint when their required
%toolbox is missing; a forced engine that is unavailable or cannot handle the
%design errors clearly at run time (psyrat_resolve_engine). The popup is listed in
%cascade-preference order (Auto, CmdStan, native HMC, native fitlme), so the popup
%index maps to the engine code via engine_codes, NOT directly.
row = row - rowspace;
av_engines = psyrat_estimation_engines_available();
%single source of truth for the popup option order and the index->engine-code
%map (shared with the selection read site via local_engine_popup_map so the two
%cannot drift; a drift would silently store the wrong engine code)
[engine_strings,engine_codes] = local_engine_popup_map();
%annotate the native options when their required toolbox is missing
if ~av_engines.hmc
    engine_strings{3} = [engine_strings{3} ' (toolbox not installed)'];
end
if ~av_engines.fitlme
    engine_strings{4} = [engine_strings{4} ' (toolbox not installed)'];
end
engine_value = find(engine_codes == psyrat_prefs.proc.engine,1);
if isempty(engine_value)
    engine_value = 1;   % default to Automatic for any stray stored code
end
engine_help = ['Estimation engine. Automatic uses CmdStan when it is installed ',...
    'and otherwise falls back to a native MATLAB engine that supports the design ',...
    '(native HMC, a faithful Bayesian fit, is preferred over native fitlme, whose ',...
    'intervals are approximate frequentist bootstrap intervals). The other ',...
    'choices force a specific engine; a forced engine errors clearly if its ',...
    'toolbox is not installed or it cannot handle the requested design. Native ',...
    'engines are not bit-identical to CmdStan.'];
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Estimation engine:','Tooltip',engine_help,...
    'Position', [lcol row figwidth/4 25]);

inplists(11) = uicontrol(psyrat_gui,'Style','pop','fontsize',psyrat_prefs.guis.fsize,...
    'String',engine_strings,'Value',engine_value,...
    'Tooltip',engine_help,'Position', [rcol row subset_pop_w 25]);

%define summary strings for subset selections
meas_status = psyrat_measurestatus(psyrat_data,psyrat_prefs);
group_status = psyrat_subsetstatus(psyrat_data,psyrat_prefs,'group');
event_status = psyrat_subsetstatus(psyrat_data,psyrat_prefs,'event');
time_status = psyrat_subsetstatus(psyrat_data,psyrat_prefs,'time');
dod_status = psyrat_dodcontraststatus(psyrat_data,psyrat_prefs);
if psyrat_prefs.proc.diffest == 3
    dod_btn_enable = 'on';
else
    dod_btn_enable = 'off';
end

%Since these buttons use inplists as an input, it needed to be specified
%after all of inputs had been placed in inplists
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Select Measurements',...
    'TooltipString',['Select one or more measurement columns to process (' meas_status ')'],...
    'Position', [subset_btn_x measrow+2 subset_btn_w 27],...
    'Callback',{@selectmeasurements_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'inplists',inplists});

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Select Groups',...
    'TooltipString',['Select a subset of groups to process (' group_status ')'],...
    'Position', [subset_btn_x (grouprow+2) subset_btn_w 27],...
    'Callback',{@selectgroups_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'inplists',inplists});

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Select Events',...
    'TooltipString',['Select a subset of events to process (' event_status ')'],...
    'Position', [subset_btn_x eventrow+2 subset_btn_w 27],...
    'Callback',{@selectevents_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'inplists',inplists});

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','DoD Contrast...',...
    'Enable',dod_btn_enable,...
    'TooltipString',['Configure ERP_1..ERP_4 mapping for ',...
    '(ERP_1 - ERP_2) - (ERP_3 - ERP_4). Current: ' dod_status],...
    'Position', [subset_btn_x timerow-rowspace+2 subset_btn_w 27],...
    'Callback',{@dodcontrast_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'inplists',inplists});

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Select Occasions',...
    'TooltipString',['Select a subset of occasions to process (' time_status ')'],...
    'Position', [subset_btn_x timerow+2 subset_btn_w 27],...
    'Callback',{@selecttime_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'inplists',inplists});

%next row with extra space. Anchor the action buttons BELOW the last input row
%(Estimation engine, the current value of `row`), not to timerow. The
%dynamic-reliability rows (Dimension 1/2 and Residual correlation), the splits
%rows, and the Estimation engine row were all added after the Occasion (time)
%row, so anchoring the buttons to timerow placed them on top of those controls.
%The engine row consumed one rowspace, so the gap here is 1.5 (was 2.5) to keep
%the action buttons at their original position.
row = row - rowspace*1.5;

%Create a back button that will take the user back to psyrat_start
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back to Home',...
    'Position', [figwidth/8 row figwidth/5 50],...
    'Callback',{@bb_call,psyrat_gui});

%Create button that will check the inputs and begin processing the data
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Analyze',...
    'Position', [3*figwidth/8 row figwidth/5 50],...
    'Callback',{@psyrat_exec,'psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data,...
    'inplists',inplists});

%Create button that will display preferences
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Preferences',...
    'Position', [5*figwidth/8 row figwidth/5 50],...
    'Callback',{@psyrat_procprefs,'psyrat_prefs',psyrat_prefs,'psyrat_data',...
    psyrat_data,'inplists',inplists});

%tag gui
psyrat_gui.Tag = 'psyrat_gui';

end


function bb_call(varargin)
%if back button is pushed, go back to psyrat_start

%close the psyrat_gui
close(varargin{3});

%go back to psyrat_start
psyrat_start;

end

function psyrat_procprefs(varargin)
%
%Input
% psyrat_prefs - toolbox preferences
% psyrat_data - toolbox data
%
%Output
% There are no direct outputs to the Matlab workspace. The user will have
%  the option (via the gui) to go back to the psyrat_gui or save the specified
%  preferences to be used for stan
%

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%reconcile both alias pairs (diffrescor/diffwpcov, ssubjrel/sserrvar) in the
%shared canonicalizer so the preference widgets render canonical values.
%'quiet' preserves this callback's historical silence -- the inline blocks it
%replaced never warned, and this runs on every Preferences traversal.
psyrat_prefs.proc = psyrat_canonicalize_proc(psyrat_prefs.proc,'quiet',true);
%legacy preference files may not include newer processing options
if ~isfield(psyrat_prefs.proc,'withinchain') || isempty(psyrat_prefs.proc.withinchain)
    psyrat_prefs.proc.withinchain = 1;
end
%legacy preference files may not include the CmdStan RNG seed
if ~isfield(psyrat_prefs.proc,'seed') || isempty(psyrat_prefs.proc.seed)
    psyrat_prefs.proc.seed = 12345;
end
%legacy preference files may not include the optional sampler overrides. Empty
%([]) is the intended default (use per-analysis delta / CmdStan max_depth=10), so
%only add the field when it is missing; a user-set value is preserved.
if ~isfield(psyrat_prefs.proc,'adapt_delta')
    psyrat_prefs.proc.adapt_delta = [];
end
if ~isfield(psyrat_prefs.proc,'max_treedepth')
    psyrat_prefs.proc.max_treedepth = [];
end
%legacy preference files may not include the observation-likelihood family
if ~isfield(psyrat_prefs.proc,'family') || isempty(psyrat_prefs.proc.family)
    psyrat_prefs.proc.family = 'gaussian';
end
%legacy preference files may not include the gamma scale parameterization. It
%backfills to [] (= Automatic), NOT to 1: [] is what makes the warp omit the
%argument so the ENGINE applies its context-aware default. Backfilling a concrete
%1 here would silently pin every legacy GUI user to log-nu and defeat rulings
%33-36, and backfilling a 2 would force it on designs that cannot take it.
if ~isfield(psyrat_prefs.proc,'gammascale')
    psyrat_prefs.proc.gammascale = [];
end
psyrat_prefs.proc = psyrat_sync_samplingprefs(psyrat_prefs.proc);

%find inplists
ind = find(strcmp('inplists',varargin),1);
if ~isempty(ind)
    inplists = varargin{ind+1};
    psyrat_prefs = psyrat_store_inputchoices(psyrat_prefs,psyrat_data,inplists);
end

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
minw = 980;
%height fits every processing-option row, including the four added independently
%on separate branches: the two sampler-tuning rows (adapt_delta, max_treedepth),
%the observation-likelihood family selector, and the gamma scale parameterization
%selector. Historically this was tracked as "pre-feature 840 + N rows * rowspace
%+ 8 px of slack", and 840 + 4*56 + 8 = 1072 still lands on the right number.
%
%THE ACTUAL CONSTRAINT, since the historical formula is folklore that has drifted:
%the panel lays out from row = height - 80, then -1.75*rowspace (header gap), then
%13 * -rowspace (one per input row after the first), then -1.5*rowspace (button
%gap). So the Save/Back buttons sit at y = height - 990 and are 50 tall. They are
%fully visible for any height >= 990; at 1072 the bottom margin is 82 px. Each
%added row costs one more rowspace (56).
%
%THIS MUST GROW WITH EVERY ADDED ROW. A headless geometry check compares computed
%positions and passes at equality, so it does NOT catch a control drawn below the
%window - only opening the window does. Do not add a row without both bumping this
%and looking at the result.
minh = 1072;
if ~isempty(psyrat_gui)
    psyrat_prefs.guis.pos = psyrat_gui.Position;
    psyrat_prefs.guis.pos(3) = max(psyrat_prefs.guis.pos(3),minw);
    psyrat_prefs.guis.pos(4) = max(psyrat_prefs.guis.pos(4),minh);
    close(psyrat_gui);
else
    psyrat_prefs.guis.pos=[400 200 minw minh];
end

%define layout constants for robust spacing across displays and font scaling
rowspace = 56;
labelh = 44;
row = psyrat_prefs.guis.pos(4) - 80;
lcol = 45;
labelw = 430;
colgap = 35;
rcol = lcol + labelw + colgap;
inputw = 320;
btnw = 220;
headery = row + 5;

%create the basic psyrat_prefs
psyrat_gui= psyrat_newfigure('unit','pix',...
    'position',psyrat_prefs.guis.pos,...
    'menub','no',...
    'name','Specify Processing Preferences',...
    'numbertitle','off',...
    'resize','off');

%print the gui headers
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Preferences',...
    'Position', [lcol headery labelw 25]);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Input',...
    'Position',[rcol headery inputw 25]);

%next row
row = row - rowspace*1.75;

%ask for an input for the number of stan chains
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Number of Chains:',...
    'Tooltip','Number of chains to be run in Stan (must specify at least 3)',...
    'Position', [lcol row labelw labelh]);

newprefs.nchains = uicontrol(psyrat_gui,'Style','edit','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_prefs.proc.nchains,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying the number of warmup iterations to run
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Warmup iterations:',...
    'Tooltip','Number of warmup (burn-in) iterations run in Stan',...
    'Position', [lcol row labelw labelh]);

newprefs.nwarmup = uicontrol(psyrat_gui,'Style','edit','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_prefs.proc.nwarmup,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying the number of sampling iterations to run
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Sampling iterations:',...
    'Tooltip','Number of post-warmup draws to save per chain in Stan',...
    'Position', [lcol row labelw labelh]);

newprefs.nsampling = uicontrol(psyrat_gui,'Style','edit','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_prefs.proc.nsampling,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying the CmdStan RNG seed (controls reproducibility)
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Random seed:',...
    'Tooltip','CmdStan RNG seed. Reusing the same seed reproduces the same estimates (default 12345).',...
    'Position', [lcol row labelw labelh]);

newprefs.seed = uicontrol(psyrat_gui,'Style','edit','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_prefs.proc.seed,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying the NUTS target acceptance rate (adapt_delta). Leave blank
%to keep each model's tuned default; a value in (0,1) overrides it.
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Adapt delta (blank = model default):',...
    'Tooltip',['NUTS target acceptance rate passed to CmdStan. Leave blank to ',...
    'use each model''s tuned default (0.90-0.98). Enter a value between 0 and 1 ',...
    '(e.g., 0.95 or 0.99) to reduce divergent transitions at the cost of speed.'],...
    'Position', [lcol row labelw labelh]);

newprefs.adapt_delta = uicontrol(psyrat_gui,'Style','edit','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_prefs.proc.adapt_delta,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying the NUTS maximum tree depth (max_treedepth). Leave blank to
%keep CmdStan's default of 10; a positive integer overrides it.
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Max tree depth (blank = 10):',...
    'Tooltip',['NUTS maximum tree depth passed to CmdStan. Leave blank to use ',...
    'CmdStan''s default of 10. Enter a positive integer (e.g., 12-15) to avoid ',...
    'premature trajectory termination on hard posteriors.'],...
    'Position', [lcol row labelw labelh]);

newprefs.max_treedepth = uicontrol(psyrat_gui,'Style','edit','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_prefs.proc.max_treedepth,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying whether to use verbose stan output
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Verbose Stan output (print each iteration)',...
    'Tooltip','Displays Stan output in the Matlab command window',...
    'Position', [lcol row labelw labelh]);

newprefs.verbose = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,'String',{'No';'Yes'},...
    'Value',psyrat_prefs.proc.verbose,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying whether to view trace plots
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','View trace plots prior to saving Stan output',...
    'Tooltip','Displays trace plots to give the user the option to rerun',...
    'Position', [lcol row labelw labelh]);

newprefs.traceplots = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,'String',{'No';'Yes'},...
    'Value',psyrat_prefs.proc.traceplots,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying whether to estimate subject-specific error variance
expl_str = ['Subject-specific error variances can be used to calculate'...
    ' subject-level reliability estimates. With an occasion/time column they'...
    ' are supported for the non-difference test-retest design and the'...
    ' dynamic-reliability designs; they are not supported for non-dynamic'...
    ' DIFFERENCE scores in a test-retest workflow.'];

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Estimate subject-specific error variances',...
    'Tooltip',expl_str,...
    'Position', [lcol row labelw labelh]);

newprefs.sserrvar = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,'String',{'No';'Yes'},...
    'Value',psyrat_prefs.proc.ssubjrel,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying the observation-likelihood family (gaussian/gamma)
expl_str = ['Observation likelihood for the single-trial scores. Gaussian '...
    '(default) suits ERP amplitudes. Gamma / scaled chi-square suits strictly '...
    'positive, right-skewed data such as single-trial time-frequency power, and '...
    'is available for the one-facet internal-consistency analyses, the crossed '...
    'test-retest design at the group or subject level, the single-occasion '...
    'subject-level design, the '...
    'two-event and four-event difference scores, the one-facet or '...
    'trial+occasion two-facet dynamic/conditional (dimensional) reliability at '...
    'the group or subject level, '...
    'the one-facet dynamic/conditional two-event difference score at the group '...
    'or subject level, its concurrent (correlated-residual) subject-level form '...
    'with a population residual correlation at one or two facets, '...
    'and the person-specific dynamic difference-of-differences designs '...
    '(single-occasion or test-retest). '...
    'Gamma is NOT available for the trial+occasion two-facet dynamic DIFFERENCE '...
    'designs, the group-level concurrent dynamic difference designs, the '...
    'per-subject-correlation variants, or any split design, and it always runs '...
    'on CmdStan (the native engines fit Gaussian designs only). '...
    'Non-positive data are rejected under Gamma.'];

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Observation likelihood family',...
    'Tooltip',expl_str,...
    'Position', [lcol row labelw labelh]);

%single source of truth for the popup option order and the index->family-key map
%(shared with the read site via local_family_popup_map so the two cannot drift)
[family_labels,family_keys] = local_family_popup_map();
family_value = find(strcmpi(family_keys,char(psyrat_prefs.proc.family)),1);
if isempty(family_value)
    family_value = 1;   % default to Gaussian for any stray stored value
end
newprefs.family = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,'String',family_labels,...
    'Value',family_value,'Tooltip',expl_str,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying the gamma scale parameterization (gammascale). Applies only
%under the Gamma family, whose dynamic model is the one this selects for. On the
%Gaussian path Automatic ([]) and Log-nu (1) are inert, but Location-scale (2) is
%NOT: psyrat_computevarcomp rejects gammascale==2 outside the gamma family and
%aborts the run. Nothing gates this popup on the selected family, so the tooltip
%has to say so.
expl_str = ['Which quantity the Gamma dispersion submodel is built on. '...
    'Automatic (recommended) lets the analysis choose: location-scale for every '...
    'design that supports it, log-nu for the rest. '...
    'Log-nu models RELATIVE dispersion (constant coefficient of variation), so '...
    'an event, person or dimension contrast on it blends amplitude with error. '...
    'Location-scale models the ABSOLUTE residual SD directly, which is what you '...
    'want when the output is read as precision, individual error variance, or '...
    'is compared against a Gaussian analysis. '...
    'The two are different estimands, not two views of one fit, so results are '...
    'not comparable across the setting. Under the Gaussian family, Automatic and '...
    'Log-nu are ignored, but Location-scale is REJECTED and stops the run - the '...
    'Gaussian dynamic model is already location-scale, so leave this on Automatic.'];

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Gamma scale parameterization',...
    'Tooltip',expl_str,...
    'Position', [lcol row labelw labelh]);

%shared map keeps the build/read order from drifting (see local_family_popup_map)
[gscale_labels,gscale_values] = local_gammascale_popup_map();
gscale_value = find(cellfun(@(v) isequal(v, psyrat_prefs.proc.gammascale), ...
    gscale_values),1);
if isempty(gscale_value)
    gscale_value = 1;   % default to Automatic for any stray stored value
end
newprefs.gammascale = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,'String',gscale_labels,...
    'Value',gscale_value,'Tooltip',expl_str,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying whether to estimate difference-score reliability
expl_str = ['Estimate either a two-event difference score (X - Y) or a ',...
    'four-event difference-of-differences contrast ',...
    '(ERP_1 - ERP_2) - (ERP_3 - ERP_4). ',...
    'Two-event difference scores support an occasion facet: with ',...
    'subject-level reliability off (the static two-facet test-retest ',...
    'design), or with it on when a dimension column is selected (the ',...
    'dynamic two-facet designs). ',...
    'Difference-of-differences with an occasion facet requires dynamic ',...
    '(dimensional) reliability with subject-level reliability on (the ',...
    'person-specific retest design); otherwise it is single-session only.'];

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Estimate internal consistency of difference scores',...
    'Tooltip',expl_str,...
    'Position', [lcol row labelw labelh]);

newprefs.diffest = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',{'No';'Yes: 2-event difference';'Yes: 4-event DoD'},...
    'Value',psyrat_prefs.proc.diffest,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

expl_str = ['Enable CmdStan within-chain parallelization for all Stan ',...
    'models in PsyRAT. Threads per chain are selected automatically ',...
    'so at least 2 CPU cores remain free after chains*threads_per_chain.'];

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Use within-chain parallelization',...
    'Tooltip',expl_str,...
    'Position', [lcol row labelw labelh]);

newprefs.withinchain = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,'String',{'No';'Yes'},...
    'Value',psyrat_prefs.proc.withinchain,...
    'Position', [rcol row+5 inputw 34]);

%next row
row = row - rowspace;

%input for specifying whether to estimate residual covariance for
%co-occurring difference-score events
expl_str = ['For co-occurring events (e.g., ipsi/contra), estimate the ',...
    'within-person residual covariance. Set to No to force wp cov = 0.'];

%widths use the shared labelw/inputw constants like every other row on this page.
%They previously derived from the figure HEIGHT (psyrat_prefs.guis.pos(4)/2),
%which silently widened this row every time minh grew: at the merged minh=1016 the
%label ran 39 px into the input column and the popup was clipped 34 px past the
%980 px figure width. Decoupling it from the height fixes that permanently.
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Estimate residual covariance for co-occurring events',...
    'Tooltip',expl_str,...
    'Position', [lcol row labelw 30]);

newprefs.diffwpcov = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',{'No';'Yes'},...
    'Value',psyrat_prefs.proc.diffwpcov,...
    'Position', [rcol row+5 inputw 30]);

%next row with extra space
row = row - rowspace*1.5;

%Create a back button that will save inputs for preferences
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [lcol+60 row btnw 50],...
    'Callback',{@psyrat_prefs_save,'psyrat_prefs',psyrat_prefs,'psyrat_data',...
    psyrat_data,'newprefs',newprefs});

%Create button that will go back to psyrat_gui without saving
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [rcol+50 row btnw 50],...
    'Callback',{@psyrat_prefs_back,'psyrat_prefs',psyrat_prefs,'psyrat_data',...
    psyrat_data});

%tag gui
psyrat_gui.Tag = 'psyrat_gui';

%recommend that the user use at least 10k iterations if running retest
%reliability analyses. Raised LAST, after every control above exists: the
%advisory previously fired before the figure was built, so the screen buried
%it by creation order (a live advisory instance of the B20 hazard); moving it
%here removes that guaranteed burial, and the on-screen look is the live
%click-through's to verify. Deliberately NON-modal per the 2026-08-17
%notification ruling: the user should see the advice but may proceed.
if psyrat_prefs.proc.inp.time ~= length(psyrat_data.proc.collist) && ...
        psyrat_prefs.proc.niter < 10000
    str = {'Test-retest workflows usually require at least 10,000 iterations.'};
    str{end+1} = 'You currently selected fewer iterations, which may prevent convergence.';
    str{end+1} = 'How to fix: open Preferences and increase warmup and/or sampling iterations.';
    errordlg(str,'Recommended Iteration Increase');
end

end

function psyrat_prefs_save(varargin)
%if the preferences are to be saved

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find newprefs
ind = find(strcmp('newprefs',varargin),1);
newprefs = varargin{ind+1};

%validate chain and iteration inputs before leaving this gui.
%
%Every check below leads with ~isreal because str2double parses complex
%literals ('3+2i'), and none of the usual guards catch one: isnan is false,
%isfinite is true, and ordering comparisons use only the real part. Three
%different things happened, depending on the field:
%
%  nwarmup/nsampling/seed/max_treedepth -- the mod() term here errored
%    outright, since mod is undefined for complex arguments ("Argument must
%    be real"), so the user got an unhandled callback error instead of the
%    errordlg.
%  nchains -- no mod() term here, so a complex value left this function and
%    died later, in psyrat_computevarcomp's own nchains guard, which IS
%    mod-based; the error therefore surfaced only after estimation setup had
%    begun. It never reached CmdStan.
%  adapt_delta -- the only field with no realness check anywhere on its path.
%    A complex value cleared this check, cleared psyrat_computevarcomp's
%    adapt_delta guard (no mod term, and its comparisons use the real part),
%    cleared MatlabStan's validateattributes, and reached CmdStan as the
%    malformed token "delta=0.5+0.5i" via num2str. Silent end to end.
%
%Leading with ~isreal short-circuits before mod is ever reached and keeps
%every field on the documented errordlg path.
nchains = str2double(newprefs.nchains.String);
nwarmup = str2double(newprefs.nwarmup.String);
nsampling = str2double(newprefs.nsampling.String);

if ~isreal(nchains) || isnan(nchains) || nchains < 3
    errordlg({'Number of chains is invalid.';...
        'How to fix: set "Number of Chains" to an integer value of 3 or higher.'},...
        'Invalid Number of Chains','modal');
    if ishghandle(newprefs.nchains)
        uicontrol(newprefs.nchains);
    end
    return;
end

if ~isreal(nwarmup) || isnan(nwarmup) || nwarmup < 0 || mod(nwarmup,1) ~= 0
    errordlg({'Warmup iterations are invalid.';...
        'How to fix: set "Warmup iterations" to an integer value of 0 or higher.'},...
        'Invalid Warmup Iterations','modal');
    if ishghandle(newprefs.nwarmup)
        uicontrol(newprefs.nwarmup);
    end
    return;
end

if ~isreal(nsampling) || isnan(nsampling) || nsampling <= 0 || mod(nsampling,1) ~= 0
    errordlg({'Sampling iterations are invalid.';...
        'How to fix: set "Sampling iterations" to an integer value greater than 0.'},...
        'Invalid Sampling Iterations','modal');
    if ishghandle(newprefs.nsampling)
        uicontrol(newprefs.nsampling);
    end
    return;
end

seed = str2double(newprefs.seed.String);
if ~isreal(seed) || isnan(seed) || seed < 0 || mod(seed,1) ~= 0
    errordlg({'Random seed is invalid.';...
        'How to fix: set "Random seed" to an integer value of 0 or higher (default 12345).'},...
        'Invalid Random Seed','modal');
    if ishghandle(newprefs.seed)
        uicontrol(newprefs.seed);
    end
    return;
end

%validate the optional sampler overrides. A blank field means "use the default"
%and is stored as []; a supplied value must be well formed (the same checks
%psyrat_computevarcomp enforces).
adapt_delta_str = strtrim(char(newprefs.adapt_delta.String));
if isempty(adapt_delta_str)
    adapt_delta = [];
else
    adapt_delta = str2double(adapt_delta_str);
    if ~isreal(adapt_delta) || isnan(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1
        errordlg({'Adapt delta is invalid.';...
            ['How to fix: leave "Adapt delta" blank to use the model default, ',...
            'or set a value between 0 and 1 (e.g., 0.95).']},...
            'Invalid Adapt Delta','modal');
        if ishghandle(newprefs.adapt_delta)
            uicontrol(newprefs.adapt_delta);
        end
        return;
    end
end

max_treedepth_str = strtrim(char(newprefs.max_treedepth.String));
if isempty(max_treedepth_str)
    max_treedepth = [];
else
    max_treedepth = str2double(max_treedepth_str);
    if ~isreal(max_treedepth) || isnan(max_treedepth) || max_treedepth <= 0 || mod(max_treedepth,1) ~= 0
        errordlg({'Max tree depth is invalid.';...
            ['How to fix: leave "Max tree depth" blank to use CmdStan''s ',...
            'default of 10, or set a positive integer (e.g., 12).']},...
            'Invalid Max Tree Depth','modal');
        if ishghandle(newprefs.max_treedepth)
            uicontrol(newprefs.max_treedepth);
        end
        return;
    end
end

%pull the preferences and store them in psyrat_prefs
psyrat_prefs.proc.nchains = nchains;
psyrat_prefs.proc.nwarmup = nwarmup;
psyrat_prefs.proc.nsampling = nsampling;
psyrat_prefs.proc.niter = nwarmup + nsampling;
psyrat_prefs.proc.seed = seed;
psyrat_prefs.proc.adapt_delta = adapt_delta;
psyrat_prefs.proc.max_treedepth = max_treedepth;
psyrat_prefs.proc.verbose = newprefs.verbose.Value;
psyrat_prefs.proc.traceplots = newprefs.traceplots.Value;
psyrat_prefs.proc.ssubjrel = newprefs.sserrvar.Value;
psyrat_prefs.proc.sserrvar = psyrat_prefs.proc.ssubjrel;
%map the family popup index back to the family key string (default Gaussian for
%any stray value); shared map keeps the build/read order from drifting
[~,family_keys] = local_family_popup_map();
fidx = newprefs.family.Value;
if fidx >= 1 && fidx <= numel(family_keys)
    psyrat_prefs.proc.family = family_keys{fidx};
else
    psyrat_prefs.proc.family = 'gaussian';
end
%map the gamma scale popup index back to its value ([] = Automatic, i.e. let the
%engine apply its context-aware default). Any stray value falls back to Automatic
%rather than to a concrete parameterization, so a corrupted preference cannot
%silently pin an estimand the user did not choose.
[~,gscale_values] = local_gammascale_popup_map();
gidx = newprefs.gammascale.Value;
if gidx >= 1 && gidx <= numel(gscale_values)
    psyrat_prefs.proc.gammascale = gscale_values{gidx};
else
    psyrat_prefs.proc.gammascale = [];
end
psyrat_prefs.proc.diffest = newprefs.diffest.Value;
psyrat_prefs.proc.withinchain = newprefs.withinchain.Value;
psyrat_prefs.proc.diffwpcov = newprefs.diffwpcov.Value;
%reconcile both alias pairs after the widget writes: diffrescor follows the
%widget-fresh diffwpcov (the ssubjrel pair was already written consistently
%above). 'quiet' is mandatory here, not stylistic -- diffrescor still holds
%the PREVIOUS setting's mirror when the user changes the popup, so a warning
%policy at this site would announce a "conflict" on any Save that changes it.
psyrat_prefs.proc = psyrat_canonicalize_proc(psyrat_prefs.proc,'quiet',true);

%find psyrat_gui
psyrat_gui = findobj('Tag','psyrat_gui');

%close psyrat_gui
close(psyrat_gui);

%execute psyrat_startproc_gui with the new preferences
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

end

function psyrat_prefs_back(varargin)
%if the back button was pressed then the inputs will not be saved

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);
%reconcile ONLY the ssubjrel/sserrvar pair, as this callback always has.
%The ssubjrel-only restriction is load-bearing, not cosmetic: canonicalizing
%the diff pair here too would change what downstream readers see for
%hand-built structs (see the canonicalizer's help). Do not quote the option
%literal in comments -- TestCanonicalizeProc counts it over this source.
psyrat_prefs.proc = psyrat_canonicalize_proc(psyrat_prefs.proc, ...
    'quiet',true,'pairs','ssubjrel');

%find psyrat_gui
psyrat_gui = findobj('Tag','psyrat_gui');

%close psyrat_gui
close(psyrat_gui);

%execute psyrat_startproc_gui with the old preferences
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

end

function proc = psyrat_sync_samplingprefs(proc)
%Ensure processing preferences expose warmup/sampling plus legacy total.

if ~isfield(proc,'niter') || isempty(proc.niter)
    proc.niter = 10000;
end

if ~isfield(proc,'nwarmup') || isempty(proc.nwarmup) || ...
        ~isfield(proc,'nsampling') || isempty(proc.nsampling)
    total_iters = max(round(proc.niter),2);
    proc.nwarmup = max(floor(total_iters/2),1);
    proc.nsampling = total_iters - proc.nwarmup;
end

proc.nwarmup = max(0,round(proc.nwarmup));
proc.nsampling = max(1,round(proc.nsampling));
proc.niter = proc.nwarmup + proc.nsampling;
end

function selectmeasurements_call(varargin)
%select measurement columns to be processed in the dataset

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find inplists
ind = find(strcmp('inplists',varargin),1);
if ~isempty(ind)
    inplists = varargin{ind+1};
    psyrat_prefs = psyrat_store_inputchoices(psyrat_prefs,psyrat_data,inplists);
end

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    psyrat_prefs.guis.pos = psyrat_gui.Position;
    close(psyrat_gui);
else
    psyrat_prefs.guis.pos=[400 400 600 400];
end

%pull a list of measurement columns from the loaded data
mnames = psyrat_data.raw.colnames;
mind = psyrat_coercemeasureindices(psyrat_prefs.proc.inp.meas,length(mnames));

%define parameters for figure position
figwidth = 500;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

psyrat_gui_whichmeasurements = psyrat_newfigure('unit','pix','Visible','off',...
    'position',[400 400 figwidth figheight],...
    'menub','no',...
    'name','Specify Which Measurements to Process',...
    'numbertitle','off',...
    'resize','off');

movegui(psyrat_gui_whichmeasurements,'center');

str = {'Select one or more measurement columns to process:'};

%Write text
uicontrol(psyrat_gui_whichmeasurements,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 25]);

%bump down to next row
row = row - rowspace*8.5;

mlist = uicontrol(psyrat_gui_whichmeasurements,'style','list',...
    'min',0,'max',length(mnames),...
    'Value',mind,...
    'Position',[.5*figwidth/4 row 3*figwidth/4 figheight/2],...
    'string',mnames);

%Create a button that will go back to psyrat_proc without saving
uicontrol(psyrat_gui_whichmeasurements,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@psyrat_whichmeasurements_back_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data});

%Create button that will save measurement inputs
uicontrol(psyrat_gui_whichmeasurements,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [4.5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@psyrat_whichmeasurements_save_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'mlist',mlist});

%display gui
set(psyrat_gui_whichmeasurements,'Visible','on');

%tag the gui
psyrat_gui_whichmeasurements.Tag = 'psyrat_gui_whichmeasurements';
end

function psyrat_whichmeasurements_back_call(varargin)
%go back to psyrat_startproc and do not save measurement selections

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find psyrat_gui_whichmeasurements
psyrat_gui_whichmeasurements = findobj('Tag','psyrat_gui_whichmeasurements');

%close psyrat_gui
close(psyrat_gui_whichmeasurements);

%execute psyrat_startproc_gui with the old preferences
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

end

function psyrat_whichmeasurements_save_call(varargin)
%go back to psyrat_startproc and save measurement selections

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find mlist
ind = find(strcmp('mlist',varargin),1);
mlistraw = varargin{ind+1}.Value;
if isempty(mlistraw)
    errordlg({'No measurement columns were selected.';...
        'How to fix: choose one or more columns and click Save.'},...
        'Measurement Not Defined','modal');
    return;
end
mlist = psyrat_coercemeasureindices(mlistraw,length(psyrat_data.raw.colnames));

if isempty(mlist)
    errordlg({'No measurement columns were selected.';...
        'How to fix: choose one or more columns and click Save.'},...
        'Measurement Not Defined','modal');
    return;
end

psyrat_prefs.proc.inp.meas = mlist;

%find psyrat_gui_whichmeasurements
psyrat_gui_whichmeasurements = findobj('Tag','psyrat_gui_whichmeasurements');

%close psyrat_gui
close(psyrat_gui_whichmeasurements);

%execute psyrat_startproc_gui with the new preferences
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
end

function selectgroups_call(varargin)
%select groups to be processed in the dataset
%
%Input
% psyrat_prefs - toolbox preferences
% psyrat_data - toolbox data
%
%Output
% There are no direct outputs to the Matlab workspace. The inputs will be
%  stored for later processing
%

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find inplists
ind = find(strcmp('inplists',varargin),1);
if ~isempty(ind)
    inplists = varargin{ind+1};
    psyrat_prefs = psyrat_store_inputchoices(psyrat_prefs,psyrat_data,inplists);
end

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    psyrat_prefs.guis.pos = psyrat_gui.Position;
    close(psyrat_gui);
else
    psyrat_prefs.guis.pos=[400 400 600 400];
end

%ensure the 'none' sentinel (the last collist entry) is not selected. If it
%is, spit out an error and go back to psyrat_startproc_gui. The check keys on
%the sentinel INDEX, not the name, so a real column named 'none' stays
%selectable (B25).
if psyrat_prefs.proc.inp.group == length(psyrat_data.proc.collist)
    psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
    errordlg({'No group column is selected.';...
        'How to fix: choose a data column for Group before opening group selection.'},...
        'Group Column Not Defined','modal');
    return;
end

%pull a list of groups from the file based on the input from psyrat_gui
gnames = unique(psyrat_data.raw.data.(psyrat_data.proc.collist{...
    psyrat_prefs.proc.inp.group}));

%check whether this function has been called before and whether the column
%assigned to group has changed
if isempty(psyrat_prefs.proc.inp.whichgroupscol) || ...
        (~isempty(psyrat_prefs.proc.inp.whichgroupscol) && ...
        psyrat_prefs.proc.inp.whichgroupscol ~= psyrat_prefs.proc.inp.group)
    psyrat_prefs.proc.inp.whichgroupscol = psyrat_prefs.proc.inp.group;
    psyrat_prefs.proc.inp.whichgroups = gnames;
    gind = 1:length(gnames);

    %if the column assigned to group is the same, then pull the previous inputs
elseif (~isempty(psyrat_prefs.proc.inp.whichgroupscol) && ...
        psyrat_prefs.proc.inp.whichgroupscol == psyrat_prefs.proc.inp.group)
    gind = zeros(1,length(psyrat_prefs.proc.inp.whichgroups));
    for ii = 1:length(psyrat_prefs.proc.inp.whichgroups)
        if ~isnumeric(gnames)
            gind(ii) = find(strcmp(psyrat_prefs.proc.inp.whichgroups{ii},...
                gnames));
        elseif isnumeric(gnames)
            gind(ii) = find(gnames(:) == ...
                psyrat_prefs.proc.inp.whichgroups{ii});
        end
    end
end

%define parameters for figure position
figwidth = 500;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

psyrat_gui_whichgroups = psyrat_newfigure('unit','pix','Visible','off',...
    'position',[400 400 figwidth figheight],...
    'menub','no',...
    'name','Specify Which Groups to Process',...
    'numbertitle','off',...
    'resize','off');

movegui(psyrat_gui_whichgroups,'center');

str = {'Select which groups you would like to be processed:'};

%Write text
uicontrol(psyrat_gui_whichgroups,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 25]);

%bump down to next row
row = row - rowspace*8.5;

glist = uicontrol(psyrat_gui_whichgroups,'style','list',...
    'min',0,'max',length(gnames),...
    'Value',gind,...
    'Position',[.5*figwidth/4 row 3*figwidth/4 figheight/2],...
    'string',gnames);

%Create a button that will go back to psyrat_proc without saving
uicontrol(psyrat_gui_whichgroups,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@psyrat_whichgroups_back_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data});

%Create button that will save group inputs
uicontrol(psyrat_gui_whichgroups,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [4.5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@psyrat_whichgroups_save_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'glist',glist,'gnames',gnames});

%display gui
set(psyrat_gui_whichgroups,'Visible','on');

%tag the gui
psyrat_gui_whichgroups.Tag = 'psyrat_gui_whichgroups';
end

function psyrat_whichgroups_back_call(varargin)
%go back to psyrat_startproc and do not save which groups should be processed

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find psyrat_gui_whichgroups
psyrat_gui_whichgroups = findobj('Tag','psyrat_gui_whichgroups');

%close psyrat_gui
close(psyrat_gui_whichgroups);

%execute psyrat_startproc_gui with the old preferences
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

end

function psyrat_whichgroups_save_call(varargin)
%go back to psyrat_startproc and save which groups should be processed

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find glist and gnames
ind = find(strcmp('glist',varargin),1);
glist = varargin{ind+1}.Value;
ind = find(strcmp('gnames',varargin),1);
gnames = varargin{ind+1};

groups = cell(1,length(glist));

if ~isnumeric(gnames)
    for ii = 1:length(glist)
        groups{ii} = gnames{glist(ii)};
    end
elseif isnumeric(gnames)
    for ii = 1:length(glist)
        groups{ii} = gnames(glist(ii));
    end
end

psyrat_prefs.proc.inp.whichgroups = groups;

%find psyrat_gui_whichgroups
psyrat_gui_whichgroups = findobj('Tag','psyrat_gui_whichgroups');

%close psyrat_gui
close(psyrat_gui_whichgroups);

%execute psyrat_startproc_gui with the old preferences
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
end

function selectevents_call(varargin)
%select events to be processed in the dataset
%
%Input
% psyrat_prefs - toolbox preferences
% psyrat_data - toolbox data
%
%Output
% There are no direct outputs to the Matlab workspace. The inputs will be
%  stored for later processing
%

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find inplists
ind = find(strcmp('inplists',varargin),1);
if ~isempty(ind)
    inplists = varargin{ind+1};
    psyrat_prefs = psyrat_store_inputchoices(psyrat_prefs,psyrat_data,inplists);
end

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    psyrat_prefs.guis.pos = psyrat_gui.Position;
    close(psyrat_gui);
else
    psyrat_prefs.guis.pos=[400 400 600 400];
end

%ensure the 'none' sentinel (the last collist entry) is not selected. If it
%is, spit out an error and go back to psyrat_startproc_gui. The check keys on
%the sentinel INDEX, not the name, so a real column named 'none' stays
%selectable (B25).
if psyrat_prefs.proc.inp.event == length(psyrat_data.proc.collist)
    psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
    errordlg({'No event column is selected.';...
        'How to fix: choose a data column for Event before opening event selection.'},...
        'Event Column Not Defined','modal');
    return;
end

%pull a list of events from the file based on the input from psyrat_gui
enames = unique(psyrat_data.raw.data.(psyrat_data.proc.collist{...
    psyrat_prefs.proc.inp.event}));

%check whether this function has been called before and whether the column
%assigned to event has changed
if isempty(psyrat_prefs.proc.inp.whicheventscol) || ...
        (~isempty(psyrat_prefs.proc.inp.whicheventscol) && ...
        psyrat_prefs.proc.inp.whicheventscol ~= psyrat_prefs.proc.inp.event)
    psyrat_prefs.proc.inp.whicheventscol = psyrat_prefs.proc.inp.event;
    psyrat_prefs.proc.inp.whichevents = enames;
    eind = 1:length(enames);

    %if the column assigned to event is the same, then pull the previous inputs
elseif (~isempty(psyrat_prefs.proc.inp.whicheventscol) && ...
        psyrat_prefs.proc.inp.whicheventscol == psyrat_prefs.proc.inp.event)
    eind = zeros(1,length(psyrat_prefs.proc.inp.whichevents));
    for ii = 1:length(psyrat_prefs.proc.inp.whichevents)
        if ~isnumeric(enames)
            eind(ii) = find(strcmp(psyrat_prefs.proc.inp.whichevents{ii},...
                enames));
        elseif isnumeric(enames)
            eind(ii) = find(enames(:) == ...
                psyrat_prefs.proc.inp.whichevents{ii});
        end
    end
end

%define parameters for figure position
figwidth = 500;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

psyrat_gui_whichevents = psyrat_newfigure('unit','pix','Visible','off',...
    'position',[400 400 figwidth figheight],...
    'menub','no',...
    'name','Specify Which Events to Process',...
    'numbertitle','off',...
    'resize','off');

movegui(psyrat_gui_whichevents,'center');

str = {'Select which events you would like to be processed:'};

%Write text
uicontrol(psyrat_gui_whichevents,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 25]);

%bump down to next row
row = row - rowspace*8.5;

elist = uicontrol(psyrat_gui_whichevents,'style','list',...
    'min',0,'max',length(enames),...
    'Value',eind,...
    'Position',[.5*figwidth/4 row 3*figwidth/4 figheight/2],...
    'string',enames);

%Create a button that will go back to psyrat_proc without saving
uicontrol(psyrat_gui_whichevents,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@psyrat_whichevents_back_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data});

%Create button that will save events inputs
uicontrol(psyrat_gui_whichevents,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [4.5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@psyrat_whichevents_save_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'elist',elist,'enames',enames});

%display gui
set(psyrat_gui_whichevents,'Visible','on');

%tag the gui
psyrat_gui_whichevents.Tag = 'psyrat_gui_whichevents';
end

function psyrat_whichevents_back_call(varargin)
%go back to psyrat_startproc and do not save which events should be processed

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find psyrat_gui_whichevents
psyrat_gui_whichevents = findobj('Tag','psyrat_gui_whichevents');

%close psyrat_gui
close(psyrat_gui_whichevents);

%execute psyrat_startproc_gui with the old preferences
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

end

function psyrat_whichevents_save_call(varargin)
%go back to psyrat_startproc and save which events should be processed

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find elist and enames
ind = find(strcmp('elist',varargin),1);
elist = varargin{ind+1}.Value;
ind = find(strcmp('enames',varargin),1);
enames = varargin{ind+1};

events = cell(1,length(elist));

if ~isnumeric(enames)
    for ii = 1:length(elist)
        events{ii} = enames{elist(ii)};
    end
elseif isnumeric(enames)
    for ii = 1:length(elist)
        events{ii} = enames(elist(ii));
    end
end

%B20 (owner-ruled 2026-08-16): preserve a configured DoD contrast when the
%saved selection leaves the event set unchanged. Every Save used to wipe
%dodset/dodmap unconditionally, and the only on-screen signal was the DoD
%Contrast button's tooltip flipping to "Not configured" -- invisible without
%a deliberate hover -- while re-opening the configurator seeded a DIFFERENT,
%immediately saveable mapping (data first-appearance order, made alphabetical
%by this very Save). Set comparison, not positional: the mapping references
%event NAMES, so the same set in any order leaves every entry valid. When
%the set genuinely changes the mapping may reference removed events, so the
%wipe stays -- but it becomes VISIBLE (modal, per the B21 modality ruling)
%whenever it discards a configured contrast.
prevevents = {};
if isfield(psyrat_prefs.proc.inp,'whichevents')
    prevevents = psyrat_prefs.proc.inp.whichevents;
end
selectionchanged = ~isequal(psyrat_eventset_key(prevevents), ...
    psyrat_eventset_key(events));

%set-equality alone is NOT sufficient to preserve: switching the Event
%COLUMN rewrites inp.whichevents outside this callback without wiping the
%map (psyrat_store_inputchoices writes no dod fields), so the comparison
%above can be made against a selection the mapping was never configured
%under (caught by the 2026-08-16 refuter pass). A mapping is preservable
%only if it was configured for the CURRENT event column -- the same triple
%condition psyrat_get_dodmap validates.
mapcurrent = isfield(psyrat_prefs.proc,'dodmapcol') && ...
    ~isempty(psyrat_prefs.proc.dodmapcol) && ...
    isequal(psyrat_prefs.proc.dodmapcol, psyrat_prefs.proc.inp.event);

psyrat_prefs.proc.inp.whichevents = events;
warnwipe = false;
if selectionchanged || ~mapcurrent
    hadmap = isfield(psyrat_prefs.proc,'dodset') && ...
        isequal(psyrat_prefs.proc.dodset,1);
    psyrat_prefs.proc.dodset = 0;
    psyrat_prefs.proc.dodmap = {};
    %the dialog itself is fired BELOW, after the setup screen relaunch:
    %MATLAB modality binds only the windows that exist when the dialog is
    %created, so a warning raised here would be buried by -- and would not
    %block -- the relaunched screen (observed live 2026-08-16)
    warnwipe = hadmap;
end

%find psyrat_gui_whichevents
psyrat_gui_whichevents = findobj('Tag','psyrat_gui_whichevents');

%close psyrat_gui
close(psyrat_gui_whichevents);

%execute psyrat_startproc_gui with the new preferences
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

%B20: a discarded configured contrast must be VISIBLE. Created last, the
%modal dialog sits above the relaunched setup screen and blocks it until
%acknowledged.
if warnwipe
    warndlg(['The event selection or event column changed, so the '...
        'configured difference-of-differences contrast was cleared. '...
        'Reconfigure it with the DoD Contrast button before running.'],...
        'DoD contrast cleared','modal');
end
end

function dodcontrast_call(varargin)
%configure ERP_1..ERP_4 mapping for DoD contrast on processing screen

[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%sync input assignments from process screen controls
ind = find(strcmp('inplists',varargin),1);
if ~isempty(ind)
    inplists = varargin{ind+1};
    psyrat_prefs = psyrat_store_inputchoices(psyrat_prefs,psyrat_data,inplists);
end

if ~isfield(psyrat_prefs.proc,'diffest') || psyrat_prefs.proc.diffest ~= 3
    errordlg({'DoD contrast setup is only used when diffest = 3.';...
        'How to fix: open Preferences and set "Estimate internal consistency of difference scores" to "Yes: 4-event DoD".'},...
        'DoD mode not enabled','modal');
    return;
end

if psyrat_prefs.proc.inp.event == length(psyrat_data.proc.collist)
    errordlg({'DoD contrast setup requires an event column.';...
        'How to fix: choose an Event Type column first, then configure DoD contrast.'},...
        'Event column not specified','modal');
    return;
end

events = psyrat_get_dod_events(psyrat_data,psyrat_prefs);
if numel(events) ~= 4
    errordlg({sprintf('DoD contrast setup requires exactly 4 selected events; current selection has %d.',numel(events));...
        'How to fix: click "Select Events" and choose exactly four events before configuring DoD contrast.'},...
        'Need exactly 4 selected events','modal');
    return;
end

[dodmap, map_valid] = psyrat_get_dodmap(psyrat_prefs,events);
if ~map_valid
    dodmap = events(:)';
end

eind = zeros(1,4);
for i = 1:4
    eind(i) = find(strcmp(dodmap{i},events),1);
end

%layout
figwidth = 560;
figheight = 420;
rowspace = 46;
row = figheight - 70;
lcol = 35;
rcol = 230;

psyrat_gui_dod = psyrat_newfigure('unit','pix','Visible','off',...
    'position',[400 400 figwidth figheight],...
    'menub','no',...
    'name','Configure DoD Contrast',...
    'numbertitle','off',...
    'resize','off');
movegui(psyrat_gui_dod,'center');

uicontrol(psyrat_gui_dod,'Style','text','fontsize',15,...
    'HorizontalAlignment','center',...
    'String','Configure: (ERP_1 - ERP_2) - (ERP_3 - ERP_4)',...
    'Position',[0 row figwidth 26]);
row = row - rowspace;

uicontrol(psyrat_gui_dod,'Style','text','fontsize',11,...
    'HorizontalAlignment','left',...
    'String','ERP_1 (first minuend):',...
    'Position',[lcol row 180 26]);
plist(1) = uicontrol(psyrat_gui_dod,'Style','pop','fontsize',11,...
    'String',events(:),'Value',eind(1),...
    'Position',[rcol row 270 28]);
row = row - rowspace;

uicontrol(psyrat_gui_dod,'Style','text','fontsize',11,...
    'HorizontalAlignment','left',...
    'String','ERP_2 (first subtrahend):',...
    'Position',[lcol row 180 26]);
plist(2) = uicontrol(psyrat_gui_dod,'Style','pop','fontsize',11,...
    'String',events(:),'Value',eind(2),...
    'Position',[rcol row 270 28]);
row = row - rowspace;

uicontrol(psyrat_gui_dod,'Style','text','fontsize',11,...
    'HorizontalAlignment','left',...
    'String','ERP_3 (second minuend):',...
    'Position',[lcol row 180 26]);
plist(3) = uicontrol(psyrat_gui_dod,'Style','pop','fontsize',11,...
    'String',events(:),'Value',eind(3),...
    'Position',[rcol row 270 28]);
row = row - rowspace;

uicontrol(psyrat_gui_dod,'Style','text','fontsize',11,...
    'HorizontalAlignment','left',...
    'String','ERP_4 (second subtrahend):',...
    'Position',[lcol row 180 26]);
plist(4) = uicontrol(psyrat_gui_dod,'Style','pop','fontsize',11,...
    'String',events(:),'Value',eind(4),...
    'Position',[rcol row 270 28]);
row = row - rowspace*1.2;

uicontrol(psyrat_gui_dod,'Style','text','fontsize',10,...
    'HorizontalAlignment','left',...
    'String','Each ERP must be unique. This mapping is required before Analyze when diffest=3.',...
    'Position',[lcol row figwidth-(2*lcol) 22]);

uicontrol(psyrat_gui_dod,'Style','push','fontsize',12,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position',[figwidth/8 22 figwidth/3 72],...
    'Callback',{@psyrat_dodcontrast_back_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data});

uicontrol(psyrat_gui_dod,'Style','push','fontsize',12,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position',[4.5*figwidth/8 22 figwidth/3 72],...
    'Callback',{@psyrat_dodcontrast_save_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'plist',plist,'events',events});

set(psyrat_gui_dod,'Visible','on');
psyrat_gui_dod.Tag = 'psyrat_gui_dodcontrast';
end

function psyrat_dodcontrast_back_call(varargin)
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);
psyrat_gui_dod = findobj('Tag','psyrat_gui_dodcontrast');
if ~isempty(psyrat_gui_dod)
    close(psyrat_gui_dod);
end
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
end

function psyrat_dodcontrast_save_call(varargin)
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

ind = find(strcmp('plist',varargin),1);
plist = varargin{ind+1};
ind = find(strcmp('events',varargin),1);
events = varargin{ind+1};

vals = zeros(1,4);
for i = 1:4
    vals(i) = plist(i).Value;
end

if numel(unique(vals)) ~= 4
    errordlg({'ERP_1 to ERP_4 must map to four unique events.';...
        'How to fix: set each ERP to a different event label.'},...
        'DoD mapping requires unique events','modal');
    return;
end

dodmap = cell(1,4);
for i = 1:4
    dodmap{i} = events{vals(i)};
end

psyrat_prefs.proc.dodmap = dodmap;
psyrat_prefs.proc.dodmapcol = psyrat_prefs.proc.inp.event;
psyrat_prefs.proc.dodset = 1;

psyrat_gui_dod = findobj('Tag','psyrat_gui_dodcontrast');
if ~isempty(psyrat_gui_dod)
    close(psyrat_gui_dod);
end

psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
end

function selecttime_call(varargin)
%select occasions to be processed in the dataset
%
%Input
% psyrat_prefs - toolbox preferences
% psyrat_data - toolbox data
%
%Output
% There are no direct outputs to the Matlab workspace. The inputs will be
%  stored for later processing
%

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find inplists
ind = find(strcmp('inplists',varargin),1);
if ~isempty(ind)
    inplists = varargin{ind+1};
    psyrat_prefs = psyrat_store_inputchoices(psyrat_prefs,psyrat_data,inplists);
end

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    psyrat_prefs.guis.pos = psyrat_gui.Position;
    close(psyrat_gui);
else
    psyrat_prefs.guis.pos=[400 400 600 400];
end

%ensure the 'none' sentinel (the last collist entry) is not selected. If it
%is, spit out an error and go back to psyrat_startproc_gui. The check keys on
%the sentinel INDEX, not the name, so a real column named 'none' stays
%selectable (B25).
if psyrat_prefs.proc.inp.time == length(psyrat_data.proc.collist)
    psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
    errordlg({'No occasion/time column is selected.';...
        'How to fix: choose a data column for Occasion before opening occasion selection.'},...
        'Occasion Column Not Defined','modal');
    return;
end

%pull a list of occasions from the file based on the input from psyrat_gui
tnames = unique(psyrat_data.raw.data.(psyrat_data.proc.collist{...
    psyrat_prefs.proc.inp.time}));

%check whether this function has been called before and whether the column
%assigned to event has changed
if isempty(psyrat_prefs.proc.inp.whichtimescol) || ...
        (~isempty(psyrat_prefs.proc.inp.whichtimescol) && ...
        psyrat_prefs.proc.inp.whichtimescol ~= psyrat_prefs.proc.inp.time)
    psyrat_prefs.proc.inp.whichtimescol = psyrat_prefs.proc.inp.time;
    psyrat_prefs.proc.inp.whichtimes = tnames;
    tind = 1:length(tnames);

    %if the column assigned to occasion is the same, then pull the previous inputs
elseif (~isempty(psyrat_prefs.proc.inp.whichtimescol) && ...
        psyrat_prefs.proc.inp.whichtimescol == psyrat_prefs.proc.inp.time)
    tind = zeros(1,length(psyrat_prefs.proc.inp.whichtimes));
    for ii = 1:length(psyrat_prefs.proc.inp.whichtimes)
        if ~isnumeric(tnames)
            tind(ii) = find(strcmp(psyrat_prefs.proc.inp.whichtimes{ii},...
                tnames));
        elseif isnumeric(tnames)
            tind(ii) = find(tnames(:) == ...
                psyrat_prefs.proc.inp.whichtimes{ii});
        end
    end
end

%define parameters for figure position
figwidth = 500;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

psyrat_gui_whichtimes = psyrat_newfigure('unit','pix','Visible','off',...
    'position',[400 400 figwidth figheight],...
    'menub','no',...
    'name','Specify Which Occasions to Process',...
    'numbertitle','off',...
    'resize','off');

movegui(psyrat_gui_whichtimes,'center');

str = {'Select which occasions you would like to be processed:'};

%Write text
uicontrol(psyrat_gui_whichtimes,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 25]);

%bump down to next row
row = row - rowspace*8.5;

tlist = uicontrol(psyrat_gui_whichtimes,'style','list',...
    'min',0,'max',length(tnames),...
    'Value',tind,...
    'Position',[.5*figwidth/4 row 3*figwidth/4 figheight/2],...
    'string',tnames);

%Create a button that will go back to psyrat_proc without saving
uicontrol(psyrat_gui_whichtimes,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@psyrat_whichtimes_back_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data});

%Create button that will save events inputs
uicontrol(psyrat_gui_whichtimes,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [4.5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@psyrat_whichtimes_save_call,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'tlist',tlist,'tnames',tnames});

%display gui
set(psyrat_gui_whichtimes,'Visible','on');

%tag the gui
psyrat_gui_whichtimes.Tag = 'psyrat_gui_whichtimes';
end

function psyrat_whichtimes_back_call(varargin)
%go back to psyrat_startproc and do not save which occasions should be processed

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find psyrat_gui_whichevents
psyrat_gui_whichtimes = findobj('Tag','psyrat_gui_whichtimes');

%close psyrat_gui
close(psyrat_gui_whichtimes);

%execute psyrat_startproc_gui with the old preferences
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

end

function psyrat_whichtimes_save_call(varargin)
%go back to psyrat_startproc and save which events should be processed

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find tlist and tnames
ind = find(strcmp('tlist',varargin),1);
tlist = varargin{ind+1}.Value;
ind = find(strcmp('tnames',varargin),1);
tnames = varargin{ind+1};

times = cell(1,length(tlist));

if ~isnumeric(tnames)
    for ii = 1:length(tlist)
        times{ii} = tnames{tlist(ii)};
    end
elseif isnumeric(tnames)
    for ii = 1:length(tlist)
        times{ii} = tnames(tlist(ii));
    end
end

psyrat_prefs.proc.inp.whichtimes = times;

%find psyrat_gui_whichevents
psyrat_gui_whichtimes = findobj('Tag','psyrat_gui_whichtimes');

%close psyrat_gui
close(psyrat_gui_whichtimes);

%execute psyrat_startproc_gui with the new preferences
psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
end

function psyrat_exec(varargin)
%if execute button is pushed, parse input to run in psyrat_computevarcompwarp
%
%Input
% psyrat_prefs - toolbox preferences
% psyrat_data - toolbox data
% inplists - list of input choices from psyrat_startproc_gui
%
%Output
% No variables outputted to the Matlab workspace
% A directory will be created where the temporary files will be saved for
%  running the Stan model (stan creates various temp files). After the
%  model is done running (executed with psyrat_computevarcomp), the temporary
%  directory and files will be removed.
% The user will also be asked to specify the location of where the .mat
%  file that contains the Stan results should be saved. This .mat file will
%  be used by psyrat_startview, which will be automatically pulled up after
%  stan has finished.

%pull psyrat_prefs and psyrat_data from varargin
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);
%reconcile ONLY the ssubjrel/sserrvar pair, as this callback always has
%(see psyrat_prefs_back's twin call for why the restriction is load-bearing)
psyrat_prefs.proc = psyrat_canonicalize_proc(psyrat_prefs.proc, ...
    'quiet',true,'pairs','ssubjrel');

%find inplists
ind = find(strcmp('inplists',varargin),1);
inplists = varargin{ind+1};
psyrat_prefs = psyrat_store_inputchoices(psyrat_prefs,psyrat_data,inplists);

measinds = psyrat_prefs.proc.inp.meas;
if isempty(measinds)
    errordlg({'No measurement columns were selected.';...
        'How to fix: select one or more numeric data columns for Measurement.'},...
        'Measurement Not Defined','modal');
    if ishghandle(inplists(2))
        uicontrol(inplists(2));
    end
    return;
end

%parse inputs
psyrat_data.proc.idheader = psyrat_data.proc.collist{psyrat_prefs.proc.inp.id};
measheaders = psyrat_data.proc.collist(measinds);
psyrat_data.proc.measheaders = measheaders;
psyrat_data.proc.measheader = measheaders{1};

if psyrat_prefs.proc.inp.group ~= length(psyrat_data.proc.collist)
    psyrat_data.proc.groupheader = psyrat_data.proc.collist{psyrat_prefs.proc.inp.group};
elseif psyrat_prefs.proc.inp.group == length(psyrat_data.proc.collist)
    psyrat_data.proc.groupheader = '';
end

if psyrat_prefs.proc.inp.event ~= length(psyrat_data.proc.collist)
    psyrat_data.proc.eventheader = psyrat_data.proc.collist{psyrat_prefs.proc.inp.event};
elseif psyrat_prefs.proc.inp.event == length(psyrat_data.proc.collist)
    psyrat_data.proc.eventheader = '';
end

if psyrat_prefs.proc.inp.time ~= length(psyrat_data.proc.collist)
    psyrat_data.proc.timeheader = psyrat_data.proc.collist{psyrat_prefs.proc.inp.time};
elseif psyrat_prefs.proc.inp.time == length(psyrat_data.proc.collist)
    psyrat_data.proc.timeheader = '';
end

%dimension columns (dynamic/conditional reliability, Rast & Clayson). The last
%collist entry is the 'none' sentinel. Selecting a Dimension 1 column enables
%the analysis (psyrat_prefs.proc.dynrel = 2); Dimension 2 is optional.
nonecol = length(psyrat_data.proc.collist);
if psyrat_prefs.proc.inp.dim1 ~= nonecol
    psyrat_data.proc.dim1header = psyrat_data.proc.collist{psyrat_prefs.proc.inp.dim1};
else
    psyrat_data.proc.dim1header = '';
end
if psyrat_prefs.proc.inp.dim2 ~= nonecol
    psyrat_data.proc.dim2header = psyrat_data.proc.collist{psyrat_prefs.proc.inp.dim2};
else
    psyrat_data.proc.dim2header = '';
end

if ~isempty(psyrat_data.proc.dim1header)
    psyrat_prefs.proc.dynrel = 2;
else
    %no Dimension 1: a Dimension 2 selection on its own is invalid because the
    %dynamic-reliability model requires Dimension 1. Error here rather than
    %silently dropping the user's Dimension 2 choice. (This is the check the
    %later dynrel==2 guard could never reach, since dynrel==2 already implies a
    %non-empty Dimension 1 header.)
    if ~isempty(psyrat_data.proc.dim2header)
        errordlg({'Dimension 2 was selected without Dimension 1.';...
            'How to fix: select a Dimension 1 column (or set both to none).'},...
            'Dynamic reliability','modal');
        return;
    end
    %no dimensions selected -> standard (non-dynamic) analysis.
    psyrat_prefs.proc.dynrel = 1;
    psyrat_data.proc.dim1header = '';
    psyrat_data.proc.dim2header = '';
end

%guard the GUI against unsupported dynamic-reliability combinations before the
%run starts (the estimation core also enforces these, but a friendly dialog is
%clearer). The dynamic-reliability family (Rast & Clayson, cases 11-22) now
%supports difference scores, the occasion/time facet, and subject-level error
%variances; the engine routes diffest/sserrvar/time/diffwpcov/ssrescor to a
%specific case (11-22). Only a few combinations remain unimplemented and are
%blocked here. The facet-specific requirements that the standard difference path
%already enforces (exactly two events; the diff+time and sserrvar+time
%retest-workflow blocks) are handled by the shared validation below, which now
%exempts dynamic-reliability runs (dynrel == 2). The sserrvar+time block is also
%exempt for the NON-difference, non-dynamic case (trt_sserrvar, analysis 25).
if psyrat_prefs.proc.dynrel == 2
    %(the Dimension 2 without Dimension 1 case is handled above, where the dim
    %headers are resolved, so it is not re-checked here.)
    %difference-of-differences (DoD, diffest == 3) under dynamic reliability is
    %implemented ONLY as the person-specific NONCONCURRENT design (analyses
    %28/29): subject-level reliability on. BOTH likelihood families support it
    %(Gaussian, and Gamma location-scale; owner ruling 2026-08-20), so the
    %family is no longer part of this gate. Any other DoD+dynrel combination
    %gets the friendly early error; the engine enforces the same requirements
    %with its own guards.
    if isfield(psyrat_prefs.proc,'diffest') && psyrat_prefs.proc.diffest == 3
        dodDynrelOk = isfield(psyrat_prefs.proc,'ssubjrel') && ...
            psyrat_prefs.proc.ssubjrel == 2;
        if ~dodDynrelOk
            errordlg({['Dynamic reliability supports difference-of-differences '...
                '(DoD) scores only as the person-specific nonconcurrent design: '...
                'subject-level reliability must be Yes.'];...
                ['How to fix: set subject-level reliability to Yes, or set '...
                'difference scores to No or to Yes with exactly two events, '...
                'or clear the Dimension columns.']},...
                'Dynamic reliability','modal');
            return;
        end
    end
    %a non-difference run with subject-level error variances (diffest == 1,
    %ssubjrel == 2) is now supported: it routes to the subject-level dynamic
    %reliability estimators (analysis 26 one-facet / 27 two-facet), which report
    %each participant's conditional reliability phi_s(z) from the per-subject
    %residual the location-scale model already estimates. No guard needed here;
    %the shared preflight below allows the two-facet (occasion) case.
end

%reliability with data splits (Rocha Tables 4-5). Selecting an items-per-split
%(n_i) column declares that each row is a split mean. The Split type popup then
%declares parallel (equal n_i; splits = 2, which reuses the standard one-facet/
%test-retest machinery) vs nonparallel (unequal n_i; splits = 3, the weighted
%observed-design model). 'none' -> standard single-trial analysis (splits = 1).
%The n_i column reaches the loader via psyrat_data.proc.weightheader, mirroring
%the dimension headers resolved above. nonecol is the 'none' sentinel set above.
if isfield(psyrat_prefs.proc.inp,'weight') && psyrat_prefs.proc.inp.weight ~= nonecol
    psyrat_data.proc.weightheader = psyrat_data.proc.collist{psyrat_prefs.proc.inp.weight};
    if isfield(psyrat_prefs.proc,'splitsmode') && psyrat_prefs.proc.splitsmode == 2
        psyrat_prefs.proc.splits = 3; %nonparallel
    else
        psyrat_prefs.proc.splits = 2; %parallel
    end
else
    psyrat_data.proc.weightheader = '';
    psyrat_prefs.proc.splits = 1;
end

%reliability with data splits is a Stage-1 design: non-difference, group-level,
%and not combined with dynamic reliability. Guard the unsupported combinations
%with a friendly dialog before the run (the estimation core also hard-errors on
%these, but a clear dialog beats a batch-loop failure report).
if psyrat_prefs.proc.splits > 1
    if psyrat_prefs.proc.dynrel == 2
        errordlg({'Reliability with data splits cannot be combined with dynamic reliability.';...
            'How to fix: clear the Dimension columns, or set Items per split (n_i) to none.'},...
            'Reliability with splits','modal');
        return;
    end
    if isfield(psyrat_prefs.proc,'diffest') && psyrat_prefs.proc.diffest ~= 1
        errordlg({'Reliability with data splits does not support difference scores.';...
            'How to fix: set difference scores to No, or set Items per split (n_i) to none.'},...
            'Reliability with splits','modal');
        return;
    end
    if isfield(psyrat_prefs.proc,'ssubjrel') && psyrat_prefs.proc.ssubjrel == 2
        errordlg({'Reliability with data splits does not support subject-level error variances.';...
            'How to fix: set subject-specific error variances to No, or set Items per split (n_i) to none.'},...
            'Reliability with splits','modal');
        return;
    end
end

%check whether particular events were specified to process. If not, process
%all event types. Also, if event was changed, overwrite the old and replace
%with all the new event types.
if (isempty(psyrat_prefs.proc.inp.whicheventscol) && ...
        psyrat_prefs.proc.inp.event ~= length(psyrat_data.proc.collist)) || ...
        (~isempty(psyrat_prefs.proc.inp.whicheventscol) && ...
        psyrat_prefs.proc.inp.whicheventscol ~= psyrat_prefs.proc.inp.event)

    %pull a list of events from the file based on the input from psyrat_gui
    enames = unique(psyrat_data.raw.data.(psyrat_data.proc.collist{...
        psyrat_prefs.proc.inp.event}));

    psyrat_prefs.proc.inp.whicheventscol = psyrat_prefs.proc.inp.event;
    psyrat_prefs.proc.inp.whichevents = enames;

end

%check whether particular group were specified to process. If not, process
%all group types. Also, if group was changed, overwrite the old and replace
%with all the new group types.
if (isempty(psyrat_prefs.proc.inp.whichgroupscol) &&...
        psyrat_prefs.proc.inp.group ~= length(psyrat_data.proc.collist)) || ...
        (~isempty(psyrat_prefs.proc.inp.whichgroupscol) && ...
        psyrat_prefs.proc.inp.whichgroupscol ~= psyrat_prefs.proc.inp.group)

    %pull a list of groups from the file based on the input from psyrat_gui
    gnames = unique(psyrat_data.raw.data.(psyrat_data.proc.collist{...
        psyrat_prefs.proc.inp.group}));

    psyrat_prefs.proc.inp.whichgroupscol = psyrat_prefs.proc.inp.group;
    psyrat_prefs.proc.inp.whichgroups = gnames;

end

%check whether particular occasions were specified to process. If not, process
%all occasions. Also, if occasion was changed, overwrite the old and replace
%with all the new occasion types.
if (isempty(psyrat_prefs.proc.inp.whichtimescol) &&...
        psyrat_prefs.proc.inp.time ~= length(psyrat_data.proc.collist)) || ...
        (~isempty(psyrat_prefs.proc.inp.whichtimescol) && ...
        psyrat_prefs.proc.inp.whichtimescol ~= psyrat_prefs.proc.inp.time)

    %pull a list of occasions from the file based on the input from psyrat_gui
    tnames = unique(psyrat_data.raw.data.(psyrat_data.proc.collist{...
        psyrat_prefs.proc.inp.time}));

    psyrat_prefs.proc.inp.whichtimescol = psyrat_prefs.proc.inp.time;
    psyrat_prefs.proc.inp.whichtimes = tnames;

end

%put the information about which events, groups, and occasions to process
%in psyrat_data
if psyrat_prefs.proc.inp.event ~= length(psyrat_data.proc.collist)
    psyrat_data.proc.whichevents = psyrat_prefs.proc.inp.whichevents;
else
    psyrat_data.proc.whichevents = '';
end
if psyrat_prefs.proc.inp.group ~= length(psyrat_data.proc.collist)
    psyrat_data.proc.whichgroups = psyrat_prefs.proc.inp.whichgroups;
else
    psyrat_data.proc.whichgroups = '';
end
if psyrat_prefs.proc.inp.time ~= length(psyrat_data.proc.collist)
    psyrat_data.proc.whichtimes = psyrat_prefs.proc.inp.whichtimes;
else
    psyrat_data.proc.whichtimes = '';
end

%check if there are duplicates (other than 'none') for non-measurement
%inputs (e.g., group and event do not refer to the same column in the data)
noneind = length(psyrat_data.proc.collist);
%dimension columns (dynamic reliability) must also be unique from every other
%role and from the measurement columns, so include them in the duplicate and
%overlap checks below. Read defensively in case an older saved config predates
%the dimension pickers.
%the items-per-split (n_i) column (reliability with data splits) must also be
%unique from every other role and from the measurement columns, so include it in
%the duplicate and overlap checks. Read defensively for older saved configs.
dim1sel = noneind; dim2sel = noneind; weightsel = noneind;
if isfield(psyrat_prefs.proc.inp,'dim1'); dim1sel = psyrat_prefs.proc.inp.dim1; end
if isfield(psyrat_prefs.proc.inp,'dim2'); dim2sel = psyrat_prefs.proc.inp.dim2; end
if isfield(psyrat_prefs.proc.inp,'weight'); weightsel = psyrat_prefs.proc.inp.weight; end
facetinds = [psyrat_prefs.proc.inp.id psyrat_prefs.proc.inp.group ...
    psyrat_prefs.proc.inp.event psyrat_prefs.proc.inp.time dim1sel dim2sel weightsel];
facetkeep = [true psyrat_prefs.proc.inp.group ~= noneind ...
    psyrat_prefs.proc.inp.event ~= noneind ...
    psyrat_prefs.proc.inp.time ~= noneind ...
    dim1sel ~= noneind dim2sel ~= noneind weightsel ~= noneind];
facetused = facetinds(facetkeep);
[facetuniq,~,facetbin] = unique(facetused);
facetcounts = accumarray(facetbin(:),1);
dupfacets = facetuniq(facetcounts > 1);
if ~isempty(dupfacets)
    probcol = psyrat_data.proc.collist(dupfacets);
    dlg = {'The same source column was assigned to multiple PsyRAT inputs:'};
    for i = 1:length(probcol)
        dlg{end+1} = probcol{i}; %#ok<AGROW>
    end
    dlg{end+1} = 'How to fix: assign unique columns for ID, Group, Event, Occasion, Dimension, and Items-per-split inputs.';
    errordlg(dlg, 'Unique variable names not provided','modal');
    roleinds = [psyrat_prefs.proc.inp.id psyrat_prefs.proc.inp.group ...
        psyrat_prefs.proc.inp.event psyrat_prefs.proc.inp.time];
    rolecontrols = [1 3 4 5];
    focusrole = find(ismember(roleinds,dupfacets),1);
    if ~isempty(focusrole) && ishghandle(inplists(rolecontrols(focusrole)))
        uicontrol(inplists(rolecontrols(focusrole)));
    end
    return;
end

%check for overlap between measurement columns and other required inputs.
overlapinds = intersect(measinds,facetused);
if ~isempty(overlapinds)
    probcol = psyrat_data.proc.collist(overlapinds);
    dlg = {'Measurement columns must be unique from ID, Group, Event, Occasion, Dimension, and Items-per-split columns.';...
        'The following selected measurement columns conflict with other inputs:'};
    for i = 1:length(probcol)
        dlg{end+1} = ['  ' probcol{i}]; %#ok<AGROW>
    end
    dlg{end+1} = 'How to fix: remove overlapping columns from Measurement selection.';
    errordlg(dlg, 'Measurement column conflict','modal');
    if ishghandle(inplists(2))
        uicontrol(inplists(2));
    end
    return;
end

%make sure every selected measurement column is numeric.
nonnumeric = {};
for i = 1:length(measheaders)
    if ~isnumeric(psyrat_data.raw.data.(measheaders{i}))
        nonnumeric{end+1} = measheaders{i}; %#ok<AGROW>
    end
end
if ~isempty(nonnumeric)
    dlg = {'The following selected measurement columns are not numeric:'};
    for i = 1:length(nonnumeric)
        dlg{end+1} = ['  ' nonnumeric{i}]; %#ok<AGROW>
    end
    dlg{end+1} = 'PsyRAT requires numeric measurement values (for example ERP amplitudes).';
    dlg{end+1} = 'How to fix: choose only numeric columns for Measurement and rerun.';
    errordlg(dlg, 'Measurement data not numeric','modal');
    if ishghandle(inplists(2))
        uicontrol(inplists(2));
    end
    return;
end

%Subject-level error variance combined with an occasion facet. This is a
%friendly EARLY copy of the check psyrat_preflight_validate performs at
%:244-250; that shared validator is the authority and runs below (before
%CmdStan), but it runs after the data file is loaded, so the dialog here is
%clearer. THE TWO CONDITIONS MUST AGREE -- TestStartprocPreflightAgreement pins
%that, because their disagreement is exactly what caused RC-24.
%
%Supported combinations that must NOT be blocked (each has a named exemption in
%preflight):
%  - dynrel == 2, diffest == 2  : subject-level two-facet difference (16/20/22)
%  - dynrel == 2, non-difference: subject-level dynamic two-facet (27)
%  - non-dynamic, non-difference: trt_sserrvar, analysis 25
%
%Analysis 25 is the reason this condition changed. The block used to reject it
%outright: it was inherited from the ERA Toolbox port (2026-02), predating the
%analysis-25 engine by ~4.5 months, and the commit that added that engine relaxed
%psyrat_preflight_validate but never this file. What remains blocked here is the
%non-dynamic SUBJECT-LEVEL difference designs (ssubjrel == 2 with diffest 2 or 3)
%with a time facet; the plain two-facet difference (diffest == 2, ssubjrel == 1,
%analysis 10) passes this gate and the diffest gate below.
if (psyrat_prefs.proc.inp.time ~= length(psyrat_data.proc.collist) && ...
        psyrat_prefs.proc.ssubjrel == 2 && psyrat_prefs.proc.dynrel ~= 2 && ...
        any(psyrat_prefs.proc.diffest == [2 3]))
    dlg = {'Single-subject error variance estimation is not supported for';...
        'DIFFERENCE scores when an occasion/time facet is present';...
        '(test-retest workflow).'; ...
        'How to fix: set subject-level reliability off (ssubjrel = 1),';...
        'remove the occasion/time facet, or turn off difference scores';...
        'for this run.'};
    errordlg(dlg, 'Single-subject error variance not supported for retest','modal');
    return;
end

%make sure that if the user asked to estimate difference-score reliability,
%required events are provided and the difference/occasion combination is one the
%engine ships
if any(psyrat_prefs.proc.diffest == [2 3])
    %dynamic-reliability runs (dynrel == 2) support the occasion facet: the
    %two-facet difference variants (cases 15/16/18/20/22) are diffest == 2 with a
    %time facet, and the person-specific nonconcurrent DoD retest design
    %(analysis 29) is diffest == 3 with a time facet. The dynrel gate earlier in
    %this function permits DoD + dynrel when subject-level reliability is on
    %(analyses 28/29, in either likelihood family since the 2026-08-20
    %Gaussian arm), so BOTH classes reach this exemption
    %and are routed by the engine. (B17: this comment previously claimed DoD
    %with dynrel was blocked above; that has been false since the gate learned
    %the 28/29 design, and the dynrel ~= 2 exemption below is precisely what
    %lets analysis 29 through this retest block.)
    %
    %The static two-facet difference (diffest == 2 with ssubjrel == 1, analysis
    %10) is ALSO exempt: it is shipped, routed by the engine selector, and
    %explicitly permitted by psyrat_preflight_validate's blockDiffWithTime --
    %"the supported two-facet (test-retest) difference pipeline". This gate
    %blocked it anyway from 2026-02 until 2026-08-29 (finding G34), the same
    %inherited-port over-block class as the ssubjrel gate above (RC-24): the
    %commit that taught preflight the design never relaxed this file. THE TWO
    %CONDITIONS MUST AGREE with preflight; the exemption below mirrors
    %blockDiffWithTime exactly. What stays blocked here is DoD (diffest == 3)
    %with a time facet outside the dynamic 29 design; the subject-level
    %diffest == 2 combination (ssubjrel == 2) is already rejected by the gate
    %above, and the term below keeps this gate correct even if that ordering
    %ever changes.
    if psyrat_prefs.proc.inp.time ~= length(psyrat_data.proc.collist) && ...
            psyrat_prefs.proc.dynrel ~= 2 && ...
            ~(psyrat_prefs.proc.diffest == 2 && psyrat_prefs.proc.ssubjrel == 1)
        dlg = {'This difference-score configuration is not supported when';...
            'an occasion/time facet is present (test-retest workflow).'; ...
            'Two-event difference scores (diffest = 2) are supported;';...
            'difference-of-differences (diffest = 3) is not.';...
            'How to fix: set diffest = 1 or 2 for test-retest runs,';...
            'or remove the occasion/time facet and run single-session DoD.'};
        errordlg(dlg, 'Difference score reliability is not supported for retest','modal');
        return;
    end

    if psyrat_prefs.proc.inp.event == length(psyrat_data.proc.collist)
        dlg = {'Difference-score reliability requires an event/condition column.';...
            'How to fix: choose the column that stores event labels'; ...
            '(for example, correct/error) and rerun.'};
        errordlg(dlg, 'Event column not specified','modal');
        return;

    end

    nevents = size(psyrat_prefs.proc.inp.whichevents,1) * ...
        size(psyrat_prefs.proc.inp.whichevents,2);

    if psyrat_prefs.proc.diffest == 2
        if nevents < 2
            dlg = {'Difference-score reliability requires exactly two event types.';...
                'Only one event type is currently selected/found.'; ...
                'How to fix: select two unique events to process.'};
            errordlg(dlg, 'Only 1 event','modal');
            return;
        elseif nevents > 2
            dlg = {'Difference-score reliability requires exactly two event types.';...
                'More than two event types are currently selected/found.'; ...
                'How to fix: click "Select Events" and choose exactly two events.'};
            errordlg(dlg, 'Only 2 events for difference score','modal');
            return;
        end
    elseif psyrat_prefs.proc.diffest == 3
        if nevents < 4
            dlg = {'Difference-of-differences reliability requires exactly four event types.';...
                'Fewer than four event types are currently selected/found.'; ...
                'How to fix: click "Select Events" and choose exactly four events.'};
            errordlg(dlg, 'Need 4 events for difference-of-differences','modal');
            return;
        elseif nevents > 4
            dlg = {'Difference-of-differences reliability requires exactly four event types.';...
                'More than four event types are currently selected/found.'; ...
                'How to fix: click "Select Events" and choose exactly four events.'};
            errordlg(dlg, 'Need exactly 4 events for difference-of-differences','modal');
            return;
        end

        %Residual covariance is not available for the four-event
        %difference-of-differences design in ANY family:
        %psyrat_build_stan_dodiff has no residual-covariance term, the case-9
        %body never reads the flag, and psyrat_dod_build_virtual_rel hard-codes
        %rescor = 0. TWO engine guards enforce that, and the DoD branch alone
        %covers both: varargin:familyconcurrent (grep the id in
        %psyrat_computevarcomp.m; line numbers drift with every increment)
        %for the static group-level DoD (analysis 9), and the varargin:dynrel
        %nonconcurrent guard for the person-specific dynamic DoD (28/29), whose
        %estimand fixes all six pairwise residual covariances to zero.
        %
        %THOSE GUARDS REMAIN THE AUTHORITY -- this is only a friendly EARLY
        %copy, and it exists because of WHERE their error surfaces: psyrat_exec
        %reaches the engine through the per-measurement try/catch below, which
        %is after the uiputfile save dialog, so without this block the user
        %picks an output location and only then learns the run is rejected.
        %Same doctrine as the sserrvar+time early copy earlier in this
        %function. Because no DoD route ACCEPTS residual covariance, this
        %cannot over-block, which is the direction that matters: over-blocking
        %would make a supported design unreachable from the GUI (RC-24).
        %(A psyrat_preflight_validate mirror was tried and reverted --
        %psyrat_computevarcomp calls preflight at :900, 747 lines BEFORE the
        %guard at :1647, so the copy shadowed the guard and changed the raised
        %id out from under TestDodConcurrentGuard/TestGammaScaleOption.)
        %
        %KEYED ON diffwpcov ALONE, matching the ssrescor gate earlier in this
        %file: the two names are ONE setting. The rule that decides is the
        %ENGINE's own normalization at psyrat_computevarcomp.m:461-467 (NOT the
        %warp, which resolves the pair in whichever direction the present field
        %dictates; the engine's version also holds for direct calls):
        %  diffwpcov==2                    -> diffrescor:=2, engine REJECTS
        %  diffwpcov==1                    -> diffrescor:=1, engine ACCEPTS
        %  diffwpcov absent, diffrescor==2 -> diffwpcov:=2,  engine REJECTS
        %So diffwpcov is decisive in every state this block can see. Only the
        %third state would need diffrescor read here, and it cannot arise at
        %this point: psyrat_exec dereferences psyrat_prefs.proc.diffwpcov
        %unguarded at :2981, so an absent field throws there regardless.
        %Reading diffrescor as well would instead OVER-BLOCK the second state,
        %refusing a run the engine accepts. The GUI cannot decouple the pair
        %anyway -- diffest==3 is settable only in psyrat_prefs_save at :1382,
        %two lines above the diffwpcov widget write, whose value the
        %canonicalizer call just below it mirrors into diffrescor.
        %
        %The block's EXECUTION is not headlessly testable -- psyrat_startproc
        %is a long GUI callback -- so the rejection itself is verified by
        %click-through. Its MODALITY is pinned headlessly, however:
        %TestModalRejectionDialogs counts the modal-dialog sites in this file
        %and had to be updated from 38 to 39 for this dialog. (Do NOT write
        %that needle literally in a comment -- the count is a plain substring
        %search over the source, so a comment quoting it is itself counted.)
        dod_diffcov = isfield(psyrat_prefs.proc,'diffwpcov') && ...
            ~isempty(psyrat_prefs.proc.diffwpcov) && ...
            psyrat_prefs.proc.diffwpcov == 2;
        if dod_diffcov
            %The gamma-scale line INSTRUCTS Automatic rather than assuming it:
            %gammascale persists across runs, so a user who once chose Log-nu
            %still has it. Naming a concrete value would misroute one of the
            %two populations that reach this dialog -- a static concurrent
            %two-event gamma run is analyses 7/8/10 under log-nu, but a
            %subject-level dynamic one is 19/20, which :1465 rejects unless
            %gammascale==2. Automatic resolves each correctly (:1370-1372).
            dlg = {'Residual covariance is not available for difference-of-differences.';...
                'The four-event contrast is estimated with the within-person';...
                'residual covariance defined to be zero, in every family.'; ...
                'How to fix: in Preferences, set "Estimate residual covariance';...
                'for co-occurring events" to No. Co-occurring events ARE';...
                'supported for the two-event difference score, so you can';...
                'instead switch "Estimate internal consistency of difference';...
                'scores" to "Yes: 2-event difference" and set the gamma scale';...
                'to "Automatic (recommended)".'};
            errordlg(dlg,'Residual covariance not supported for DoD','modal');
            return;
        end

        events = psyrat_get_dod_events(psyrat_data,psyrat_prefs);
        [dodmap,map_valid] = psyrat_get_dodmap(psyrat_prefs,events);
        if ~map_valid
            dlg = {'Difference-of-differences reliability requires an explicit contrast mapping.';...
                'How to fix: click "DoD Contrast..." and assign unique events to ERP_1..ERP_4 before Analyze.'};
            errordlg(dlg, 'Configure DoD contrast before Analyze','modal');
            return;
        end

        %store a normalized mapping to ensure downstream code receives
        %stable char labels in ERP_1..ERP_4 order.
        psyrat_prefs.proc.dodmap = dodmap;
        psyrat_prefs.proc.dodset = 1;
        psyrat_prefs.proc.dodmapcol = psyrat_prefs.proc.inp.event;
    end


end

%cmdstan cannot handle paths with white space. User will be required to
%provide a path to working directory that does not include a space.
%space will be checked for each selection
while true

    %prompt the user to indicate where the output from stan should be saved
    [savename, savepath] = uiputfile(fullfile(psyrat_data.raw.filepath,'*.psyrat'),...
        'Save output file as');

    if isequal(savename,0) || isequal(savepath,0)
        errordlg({'Output file location was not selected.';...
            'How to fix: choose a filename and folder to continue.'},...
            'File Error','modal');
        return;
    end

    if any(isspace(savename)) || any(isspace(savepath))
        str = {};
        str{end+1} = 'The selected filename/path contains whitespace.';
        str{end+1} = 'CmdStan does not support whitespace in filenames or paths.';
        str{end+1} = 'How to fix: choose a filename and path with no spaces.';
        errordlg(str,'Whitespace detected','modal');
        continue;
    end

    psyrat_data.proc.savepath = savepath;
    break;

end

run_measheaders = measheaders;
if psyrat_prefs.proc.diffest == 3 && numel(run_measheaders) > 1
    run_measheaders = run_measheaders(1);
    psyrat_data.proc.measheaders = run_measheaders;
    psyrat_data.proc.measheader = run_measheaders{1};
    %B28 ruling (2026-08-18): the user must ACKNOWLEDGE this before the run
    %starts, so it is uiwait + modal. It is not a notification about output
    %already on screen -- it discloses that the requested analysis SCOPE was
    %silently narrowed to the first measurement, which is the
    %precondition/operation class.
    %
    %Non-modal it was provably unseeable. Fired here, BEFORE the
    %close(psyrat_gui) below, it survived the close, survived the end-of-batch
    %psyrat_startproc_gui relaunch buried under the new screen, and escaped
    %that dialog's modality because it was created first (the B20
    %creation-order mechanism). Confirmed live 2026-08-18 by a handle census
    %that found it still open, vis=on, after the modal was dismissed -- while
    %the visible dialog read "All 1 selected measurements failed" to a user
    %who had selected two.
    %
    %uiwait IS THE LOAD-BEARING HALF; 'modal' alone would NOT fix this, and an
    %earlier draft of this comment wrongly claimed it did. warndlg RETURNS
    %IMMEDIATELY whatever its CreateMode (only uiwait blocks), and this file
    %has no uiwait/drawnow/pause anywhere else, so execution would run straight
    %through close(psyrat_gui) into the measurement loop without ever yielding
    %to the event queue -- the user would get no opportunity to dismiss the box
    %before the run began, and it could still be open at the relaunch. What
    %'modal' buys on its own is only that the box is raised above the setup GUI
    %at creation. uiwait is what makes the acknowledgment real and therefore
    %what prevents the orphaning. Same pattern, same rationale, as
    %psyrat_startview.m's engine and convergence notices.
    uiwait(warndlg({'Difference-of-differences mode currently processes one measurement at a time.';...
        ['Using only the first selected measurement: ' run_measheaders{1}];...
        'How to fix: select a single measurement column when diffest is set to "Yes: 4-event DoD".'},...
        'DoD mode: single measurement only','modal'));
end

outputsavename = psyrat_measurement_savenames(savename,run_measheaders);

%all input checks have passed; close the gui before starting analysis
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

runsuccess = 0;
runfail = 0;
lastsuccess = [];
nmeas = length(run_measheaders);

for i = 1:nmeas
    run_data = psyrat_data;
    run_data.proc.measheader = run_measheaders{i};
    run_data.proc.savename = outputsavename{i};
    run_data.proc.batchmode = nmeas > 1;
    fprintf('\nProcessing measurement %d/%d: %s\n',i,nmeas,run_data.proc.measheader);

    try
        %set up the run-specific data table with the selected measurement.
        run_data.proc.data = psyrat_loadfile('psyrat_prefs',psyrat_prefs,...
            'psyrat_data',run_data);

        %run consolidated validation before launching CmdStan. dynrel must be
        %forwarded here (as in the warp's own preflight) so the dynamic-reliability
        %relaxation applies: the subject-level two-facet difference variants (cases
        %16/20/22; sserrvar=2 + time facet) need dynrel and would otherwise be
        %rejected by the test-retest guards before CmdStan launches. sserrvar=2 +
        %time is ALSO valid without dynrel when the design is non-difference
        %(trt_sserrvar, analysis 25); diffest is what separates the two.
        %Backfill the observation-likelihood family for legacy prefs so the gamma
        %positivity guard fires here, the earliest GUI gate before CmdStan.
        if ~isfield(psyrat_prefs.proc,'family') || isempty(psyrat_prefs.proc.family)
            psyrat_prefs.proc.family = 'gaussian';
        end
        preflight = psyrat_preflight_validate(...
            'datatable', run_data.proc.data, ...
            'sserrvar', psyrat_prefs.proc.ssubjrel, ...
            'diffest', psyrat_prefs.proc.diffest, ...
            'diffwpcov', psyrat_prefs.proc.diffwpcov, ...
            'diffrescor', psyrat_prefs.proc.diffrescor, ...
            'dynrel', psyrat_prefs.proc.dynrel, ...
            'splits', psyrat_prefs.proc.splits, ...
            'family', psyrat_prefs.proc.family, ...
            'savename', run_data.proc.savename, ...
            'savepath', run_data.proc.savepath, ...
            'requireSaveName', true, ...
            'requireSavePath', true);
        if ~preflight.ok
            diag = preflight.errors(1);
            error(diag.id,'%s',diag.message);
        end

        %Non-fatal preflight warnings were dropped here before the pre-beta
        %batch (council G8: only errors(1) was read). Raise them in one modal
        %dialog and wait for acknowledgment, so a multi-hour run does not
        %start on a design the user has never seen flagged. Same construction
        %as the DoD notice below (uiwait on the figure warndlg returns).
        if ~isempty(preflight.warnings)
            wmsgs = arrayfun(@(w) w.message, preflight.warnings, ...
                'UniformOutput', false);
            uiwait(warndlg(wmsgs(:), 'Preflight warnings', 'modal'));
        end

        %pass the data to be setup for processing.
        run_data = psyrat_computevarcompwarp('psyrat_prefs',psyrat_prefs,...
            'psyrat_data',run_data);

        if ~isfield(run_data,'rel')
            error('psyrat:missingReliabilityOutput',...
                ['No reliability output was generated for measurement ''%s''. ',...
                'How to fix: review convergence and CmdStan logs for this run.'],...
                run_data.proc.measheader);
        end

        runsuccess = runsuccess + 1;
        lastsuccess = run_data;

    catch runME
        runfail = runfail + 1;
        reportfile = psyrat_write_measurement_error(run_data.proc.savepath,...
            run_data.proc.savename,run_data.proc.measheader,runME);
        fprintf('\nMeasurement failed: %s\n',run_data.proc.measheader);
        fprintf('Reason: %s\n',runME.message);
        fprintf('Wrote error report: %s\n',reportfile);
    end
end

if runsuccess == 0
    %relaunch FIRST, create the dialog LAST (B20): modality binds only the
    %windows that exist at the dialog's creation, so a dialog fired before
    %the relaunch would be covered by the new screen AND escape its modality
    psyrat_startproc_gui('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
    errordlg({...
        sprintf('All %d selected measurements failed to process.',nmeas);...
        ['How to fix: inspect the generated .txt error reports in ' ...
        psyrat_data.proc.savepath]},...
        'All Measurements Failed','modal');
    return;
end

if nmeas == 1 && runsuccess == 1
    %single-measurement workflow keeps the existing auto-open behavior.
    psyrat_startview('psyrat_prefs',psyrat_prefs,'psyrat_data',lastsuccess);
else
    str = {sprintf('Completed %d measurement runs.',nmeas); ...
        sprintf('Successful: %d',runsuccess); ...
        sprintf('Failed: %d',runfail); ...
        ['Outputs saved in: ' psyrat_data.proc.savepath]};
    if runfail > 0
        str{end+1} = 'Failed runs produced .txt error reports with matching filenames.'; %#ok<AGROW>
    end
    msgbox(str,'Processing Complete');
end

end

function [haspairs,badpairs] = psyrat_checkdiffpairs(datatable)
%verify matched event counts so rows can be paired for co-occurring events

haspairs = true;
badpairs = table;

if ~all(ismember({'id','event'},datatable.Properties.VariableNames))
    haspairs = false;
    return;
end

eventnames = unique(datatable.event(:));
if length(eventnames) ~= 2
    haspairs = false;
    return;
end

hasgroup = any(strcmpi(datatable.Properties.VariableNames,'group'));
idvals = string(datatable.id);
eventvals = string(datatable.event);
eventnames = string(eventnames);

if hasgroup
    groupvals = string(datatable.group);
end

if hasgroup
    keytab = unique(datatable(:,{'group','id'}),'rows');
    bad_group = {};
    bad_id = {};
    bad_n1 = [];
    bad_n2 = [];
else
    keytab = unique(datatable(:,{'id'}),'rows');
    bad_id = {};
    bad_n1 = [];
    bad_n2 = [];
end

for ii = 1:height(keytab)
    if hasgroup
        ind = (idvals == string(keytab.id(ii))) & ...
            (groupvals == string(keytab.group(ii)));
    else
        ind = idvals == string(keytab.id(ii));
    end

    n1 = sum(eventvals(ind) == eventnames(1));
    n2 = sum(eventvals(ind) == eventnames(2));

    if n1 ~= n2
        haspairs = false;
        if hasgroup
            bad_group{end+1,1} = char(string(keytab.group(ii))); %#ok<AGROW>
        end
        bad_id{end+1,1} = char(string(keytab.id(ii))); %#ok<AGROW>
        bad_n1(end+1,1) = n1; %#ok<AGROW>
        bad_n2(end+1,1) = n2; %#ok<AGROW>
    end
end

if ~haspairs
    if hasgroup
        badpairs = table(bad_group,bad_id,bad_n1,bad_n2,...
            'VariableNames',{'group','id','event1_count','event2_count'});
    else
        badpairs = table(bad_id,bad_n1,bad_n2,...
            'VariableNames',{'id','event1_count','event2_count'});
    end
end

end

function subsetstr = psyrat_subsetstatus(psyrat_data, psyrat_prefs, facet)
%build a short label describing selected subsets for groups/events/occasions

switch lower(facet)
    case 'group'
        colind = psyrat_prefs.proc.inp.group;
        prevcol = psyrat_prefs.proc.inp.whichgroupscol;
        selected = psyrat_prefs.proc.inp.whichgroups;
    case 'event'
        colind = psyrat_prefs.proc.inp.event;
        prevcol = psyrat_prefs.proc.inp.whicheventscol;
        selected = psyrat_prefs.proc.inp.whichevents;
    case 'time'
        colind = psyrat_prefs.proc.inp.time;
        prevcol = psyrat_prefs.proc.inp.whichtimescol;
        selected = psyrat_prefs.proc.inp.whichtimes;
    otherwise
        subsetstr = 'Selected: all';
        return;
end

%facet unused when the index names the sentinel (last collist entry) or is a
%stale out-of-range preference; keyed by index, not name (B25)
if isempty(colind) || colind >= length(psyrat_data.proc.collist)
    subsetstr = 'Not used';
    return;
end

allvals = unique(psyrat_data.raw.data.(psyrat_data.proc.collist{colind}));
ntotal = numel(allvals);

if isempty(prevcol) || prevcol ~= colind || isempty(selected)
    nsel = ntotal;
else
    nsel = numel(selected);
end

nsel = min(nsel,ntotal);
subsetstr = sprintf('Selected: %d/%d',nsel,ntotal);

end

function subsetstr = psyrat_measurestatus(psyrat_data,psyrat_prefs)
%build a short label describing selected measurement columns.

measinds = psyrat_coercemeasureindices(psyrat_prefs.proc.inp.meas,...
    length(psyrat_data.raw.colnames));
subsetstr = sprintf('Selected: %d/%d',numel(measinds),...
    length(psyrat_data.raw.colnames));

end

function [popupstr,popupval,popupuserdata] = ...
        psyrat_measurepopupdisplay(psyrat_data,psyrat_prefs)
%define display text for the measurement popup in the main input gui.

measinds = psyrat_coercemeasureindices(psyrat_prefs.proc.inp.meas,...
    length(psyrat_data.raw.colnames));
if numel(measinds) > 1
    popupstr = {'multiple measurements'};
    popupval = 1;
    popupuserdata = measinds;
else
    popupstr = psyrat_data.raw.colnames;
    popupval = measinds(1);
    popupuserdata = [];
end

end

function psyrat_prefs = psyrat_store_inputchoices(psyrat_prefs,psyrat_data,inplists)
%pull all GUI input choices and store them in processing preferences.

choices = psyrat_extract_inputchoices(inplists,psyrat_data);
psyrat_prefs.proc.inp.id = choices.id;
psyrat_prefs.proc.inp.meas = choices.meas;
psyrat_prefs.proc.inp.group = choices.group;
psyrat_prefs.proc.inp.event = choices.event;
psyrat_prefs.proc.inp.time = choices.time;
psyrat_prefs.proc.inp.dim1 = choices.dim1;
psyrat_prefs.proc.inp.dim2 = choices.dim2;
%items-per-split (n_i) column for reliability with data splits.
psyrat_prefs.proc.inp.weight = choices.weight;
%residual-correlation parameterization (subject-level concurrent dynamic
%difference, cases 21/22). Stored at proc level (not proc.inp) to match the other
%dynamic-reliability mode flags read by the engine warp.
psyrat_prefs.proc.ssrescor = choices.ssrescor;
%split type (parallel/nonparallel) for reliability with data splits. Stored at
%proc level to match the other mode flags; the n_i column index lives in proc.inp.
psyrat_prefs.proc.splitsmode = choices.splitsmode;
%estimation-engine preference (0 auto / 1 cmdstan / 2 fitlme / 3 hmc). Stored at
%proc level; the warp (psyrat_computevarcompwarp) reads proc.engine and forwards
%it to psyrat_resolve_engine.
psyrat_prefs.proc.engine = choices.engine;

end

function choices = psyrat_extract_inputchoices(inplists,psyrat_data)
%read GUI inputs while handling scalar popup controls plus measurement
%multi-select list values.

vals = get(inplists(:),'value');
if ~iscell(vals)
    vals = {vals};
end

measud = get(inplists(2),'UserData');
if isempty(measud)
    measchoice = vals{2};
else
    measchoice = measud;
end

choices.id = psyrat_coercescalarindex(vals{1},...
    length(psyrat_data.raw.colnames),1);
choices.meas = psyrat_coercemeasureindices(measchoice,...
    length(psyrat_data.raw.colnames));
choices.group = psyrat_coercescalarindex(vals{3},...
    length(psyrat_data.proc.collist),length(psyrat_data.proc.collist));
choices.event = psyrat_coercescalarindex(vals{4},...
    length(psyrat_data.proc.collist),length(psyrat_data.proc.collist));
choices.time = psyrat_coercescalarindex(vals{5},...
    length(psyrat_data.proc.collist),length(psyrat_data.proc.collist));
%dimension pickers (inplists 6 and 7); default to 'none' when absent so older
%layouts without these controls still extract cleanly.
nonecol = length(psyrat_data.proc.collist);
if numel(vals) >= 6
    choices.dim1 = psyrat_coercescalarindex(vals{6},nonecol,nonecol);
else
    choices.dim1 = nonecol;
end
if numel(vals) >= 7
    choices.dim2 = psyrat_coercescalarindex(vals{7},nonecol,nonecol);
else
    choices.dim2 = nonecol;
end
%residual-correlation parameterization popup (inplists 8). This is a 2-state mode
%(1 = Population, 2 = Per-subject), NOT a column index, so it is clamped to [1 2]
%directly rather than coerced against collist. Default to Population (1) when the
%control is absent (older layouts) or holds a stray value.
if numel(vals) >= 8 && isscalar(vals{8}) && any(vals{8} == [1 2])
    choices.ssrescor = vals{8};
else
    choices.ssrescor = 1;
end
%items-per-split (n_i) column picker (inplists 9); default 'none' when absent
%(older layouts) so the standard single-trial analysis still extracts cleanly.
if numel(vals) >= 9
    choices.weight = psyrat_coercescalarindex(vals{9},nonecol,nonecol);
else
    choices.weight = nonecol;
end
%split type popup (inplists 10): 1 = Parallel, 2 = Nonparallel. A 2-state mode,
%not a column index, so clamp to [1 2] directly. Default Parallel (1) when the
%control is absent (older layouts) or holds a stray value.
if numel(vals) >= 10 && isscalar(vals{10}) && any(vals{10} == [1 2])
    choices.splitsmode = vals{10};
else
    choices.splitsmode = 1;
end
%estimation-engine popup (inplists 11). The popup is listed in cascade-preference
%order (Auto, CmdStan, native HMC, native fitlme), so the popup Value (1-4) maps
%to the engine code via engine_codes, matching the layout in psyrat_startproc.
%Default Automatic (0) when the control is absent (older layouts) or holds a stray
%value. The resolver re-validates a forced engine at run time.
[~,engine_codes] = local_engine_popup_map();
if numel(vals) >= 11 && isscalar(vals{11}) && any(vals{11} == 1:numel(engine_codes))
    choices.engine = engine_codes(vals{11});
else
    choices.engine = 0;
end

end

function [labels,keys] = local_family_popup_map()
%Single source of truth for the observation-likelihood-family popup: the option
%labels and the popup-index -> family-key mapping. Both the popup BUILD site and
%the selection READ site call this so the option order and the keys cannot drift
%apart (a drift would silently store the wrong family, e.g. selecting Gamma but
%persisting 'gaussian'). The keys are the strings psyrat_computevarcomp accepts.
labels = {'Gaussian';'Gamma'};
keys   = {'gaussian','gamma'};   % popup index -> family key
end

function [labels,values] = local_gammascale_popup_map()
%Single source of truth for the gamma scale-parameterization popup: the option
%labels and the popup-index -> gammascale-value mapping. Both the popup BUILD
%site and the selection READ site call this so the order and the values cannot
%drift apart, exactly as local_family_popup_map does.
%
%WHY THE FIRST OPTION IS [] AND NOT 1. [] means "unset", which makes
%psyrat_computevarcompwarp OMIT the argument, which lets psyrat_computevarcomp
%apply its context-aware default: location-scale on the designs that support it,
%log-nu everywhere else (rulings 33-36). Storing a concrete value here instead
%would force the GUI to replicate the engine's movable-design list, and the two
%would drift - the engine is deliberately the single authority on which designs
%accept which parameterization.
%
%Options 2 and 3 are what ruling 32 requires: log-nu must remain reachable as an
%EXPLICIT choice, and before this control existed it was not reachable from the
%GUI in either direction.
labels = {'Automatic (recommended)'; ...
    'Log-nu (relative dispersion / CV)'; ...
    'Location-scale (absolute residual SD)'};
values = {[], 1, 2};   % popup index -> gammascale value
end

function [labels,codes] = local_engine_popup_map()
%Single source of truth for the estimation-engine popup: the base option labels
%and the popup-index -> engine-code mapping, in cascade-preference order
%(Auto, CmdStan, native HMC, native fitlme). Both the popup BUILD site and the
%selection READ site call this so the option order and the codes cannot drift
%apart - a drift would silently store the wrong engine (e.g. selecting native HMC
%but persisting the fitlme code), mislabeling an approximate frequentist run as a
%Bayesian one. The labels are the base text; the build site annotates the native
%options with a "(toolbox not installed)" hint when their toolbox is missing.
labels = {'Automatic (recommended)','CmdStan only','Native HMC','Native fitlme'};
codes  = [0 1 3 2];   % popup index -> engine code (0 auto/1 cmdstan/3 hmc/2 fitlme)
end

function out = psyrat_coercescalarindex(value,maxind,defaultind)
%coerce GUI selection to a valid scalar index.

if nargin < 3 || isempty(defaultind)
    defaultind = 1;
end

if isempty(value)
    out = defaultind;
else
    out = value(1);
end

if ~isnumeric(out) || ~isfinite(out)
    out = defaultind;
end

out = round(out);
if out < 1 || out > maxind
    out = defaultind;
end

end

function out = psyrat_coercemeasureindices(value,maxind)
%coerce measurement listbox selections to valid unique column indices.

if maxind < 1
    out = [];
    return;
end

if isempty(value)
    value = 1;
end

out = value(:)';
out = out(isfinite(out));
out = unique(round(out),'stable');
out = out(out >= 1 & out <= maxind);

if isempty(out)
    out = 1;
end

end

function savenames = psyrat_measurement_savenames(basesavename,measheaders)
%build one output filename per measurement, adding _measurement suffixes
%when multiple measurement columns are processed.

if ~iscell(measheaders)
    measheaders = cellstr(measheaders);
end

nmeas = length(measheaders);
if nmeas <= 1
    savenames = {basesavename};
    return;
end

[basename,~,ext] = fileparts(basesavename);
if isempty(basename)
    basename = basesavename;
end
if isempty(ext)
    ext = '.psyrat';
end

savenames = cell(1,nmeas);
used = struct;
for i = 1:nmeas
    token = psyrat_measurementtoken(measheaders{i});
    basecandidate = [basename '_' token];
    key = matlab.lang.makeValidName(lower(basecandidate));
    if isfield(used,key)
        used.(key) = used.(key) + 1;
        basecandidate = sprintf('%s_%d',basecandidate,used.(key));
    else
        used.(key) = 1;
    end
    savenames{i} = [basecandidate ext];
end

end

function token = psyrat_measurementtoken(measheader)
%sanitize measurement header text so it is safe in output filenames.

token = char(string(measheader));
token = strtrim(token);
token = regexprep(token,'\s+','_');
token = regexprep(token,'[^A-Za-z0-9_-]','_');
token = regexprep(token,'_+','_');
token = regexprep(token,'^_+|_+$','');

if isempty(token)
    token = 'measurement';
end

end

function reportfile = psyrat_write_measurement_error(savepath,savename,measheader,runME)
%write a per-measurement text report when processing fails.

[basename,~,~] = fileparts(savename);
if isempty(basename)
    basename = savename;
end
if isempty(basename)
    basename = ['psyrat_' psyrat_measurementtoken(measheader)];
end

reportfile = fullfile(savepath,[basename '.txt']);
fid = fopen(reportfile,'w');
if fid == -1
    reportfile = fullfile(savepath,...
        ['psyrat_error_' datestr(now,'yyyymmdd_HHMMSS') '.txt']);
    fid = fopen(reportfile,'w');
    if fid == -1
        warning('psyrat:errorReportWriteFailed',...
            'Could not write error report for measurement ''%s'' in %s.',...
            measheader,savepath);
        return;
    end
end

cleanupobj = onCleanup(@() fclose(fid));

fprintf(fid,'PsyRAT measurement processing error\n');
fprintf(fid,'Measurement: %s\n',measheader);
fprintf(fid,'Requested output file: %s\n',savename);
fprintf(fid,'Timestamp: %s\n',datestr(now,31));
fprintf(fid,'Error identifier: %s\n',runME.identifier);
fprintf(fid,'\nError message:\n%s\n',runME.message);

if ~isempty(runME.stack)
    fprintf(fid,'\nStack trace:\n');
    for i = 1:length(runME.stack)
        fprintf(fid,'  %s (line %d)\n',runME.stack(i).file,runME.stack(i).line);
    end
end

clear cleanupobj

end

function events = psyrat_get_dod_events(psyrat_data,psyrat_prefs)
%return currently selected event labels as a cellstr vector (stable order)

events = {};

if psyrat_prefs.proc.inp.event == length(psyrat_data.proc.collist)
    return;
end

colind = psyrat_prefs.proc.inp.event;
rawvals = unique(psyrat_data.raw.data.(psyrat_data.proc.collist{colind}),'stable');
allvals = cellstr(string(rawvals(:)));

if isfield(psyrat_prefs.proc,'inp') && ...
        ~isempty(psyrat_prefs.proc.inp.whicheventscol) && ...
        psyrat_prefs.proc.inp.whicheventscol == colind && ...
        ~isempty(psyrat_prefs.proc.inp.whichevents)
    sel = psyrat_prefs.proc.inp.whichevents;
    events = cell(1,numel(sel));
    for i = 1:numel(sel)
        events{i} = char(string(sel{i}));
    end
else
    events = allvals(:)';
end

events = unique(events,'stable');
end

function keys = psyrat_eventset_key(v)
%normalize an event selection (cell of char/string/numeric entries, a bare
%array, or empty) to a sorted string vector, so two selections can be
%compared as SETS regardless of order or storage type (B20: the DoD mapping
%references event names, so an order-only difference changes nothing it
%depends on)

if isempty(v)
    keys = string.empty(1,0);
    return;
end
if ~iscell(v)
    v = num2cell(v);
end
keys = strings(1,numel(v));
for k = 1:numel(v)
    keys(k) = string(v{k});
end
keys = sort(keys);
end

function [dodmap,map_valid] = psyrat_get_dodmap(psyrat_prefs,events)
%resolve currently stored DoD mapping against available event labels

map_valid = false;
dodmap = {};

if ~isfield(psyrat_prefs.proc,'dodset') || psyrat_prefs.proc.dodset ~= 1
    return;
end
if ~isfield(psyrat_prefs.proc,'dodmap') || numel(psyrat_prefs.proc.dodmap) ~= 4
    return;
end
if ~isfield(psyrat_prefs.proc,'dodmapcol') || ...
        isempty(psyrat_prefs.proc.dodmapcol) || ...
        psyrat_prefs.proc.dodmapcol ~= psyrat_prefs.proc.inp.event
    return;
end

dodmap = cell(1,4);
for i = 1:4
    dodmap{i} = char(string(psyrat_prefs.proc.dodmap{i}));
end

if numel(unique(dodmap)) ~= 4
    dodmap = {};
    return;
end

if ~all(ismember(dodmap,events))
    dodmap = {};
    return;
end

map_valid = true;
end

function dodstr = psyrat_dodcontraststatus(psyrat_data,psyrat_prefs)
%build status message for DoD contrast setup on process screen

if ~isfield(psyrat_prefs.proc,'diffest') || psyrat_prefs.proc.diffest ~= 3
    dodstr = 'Off (enable diffest = 3 in Preferences)';
    return;
end

if psyrat_prefs.proc.inp.event == length(psyrat_data.proc.collist)
    dodstr = 'Needs event column';
    return;
end

events = psyrat_get_dod_events(psyrat_data,psyrat_prefs);
if numel(events) ~= 4
    dodstr = sprintf('Needs exactly 4 selected events (current: %d)',numel(events));
    return;
end

[dodmap,map_valid] = psyrat_get_dodmap(psyrat_prefs,events);
if ~map_valid
    dodstr = 'Not configured';
    return;
end

dodstr = sprintf('Configured: (%s-%s)-(%s-%s)',...
    dodmap{1},dodmap{2},dodmap{3},dodmap{4});
end
