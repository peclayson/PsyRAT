function psyrat_startimportscore(varargin)
%Start PsyRAT ERP import + component scoring workflow.
%
% psyrat_startimportscore
% psyrat_startimportscore('psyrat_prefs',psyrat_prefs)
%
%Optional inputs
% psyrat_prefs - PsyRAT preferences structure (default: psyrat_defaults).
%
%Opens the gui for importing ERP data files (EEGLAB .set, ERPLAB .erp, or
% EP Toolkit exported .mat/.ept), scoring single-trial component amplitudes,
% and exporting a long-format table ready for reliability analysis in
% psyrat_startproc or psyrat_run.
%
%Output
% No variables are returned to the workspace; results are saved through the
% gui.

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

%Ensure import-workflow functions are discoverable even if MATLAB path
%predates this new module.
localEnsureImportWorkflowPath();

psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

psyrat_prefs = [];
if ~isempty(varargin)
    ind = find(strcmp('psyrat_prefs',varargin),1);
    if ~isempty(ind)
        psyrat_prefs = varargin{ind+1};
    end
end
if isempty(psyrat_prefs)
    psyrat_prefs = psyrat_defaults;
    psyrat_prefs.ver = psyrat_defineversion;
    psyrat_prefs.guis.fsize = psyrat_guifsize;
end
if ~isfield(psyrat_prefs,'guis') || ~isstruct(psyrat_prefs.guis)
    psyrat_prefs.guis = struct();
end
if ~isfield(psyrat_prefs.guis,'fsize') || isempty(psyrat_prefs.guis.fsize)
    psyrat_prefs.guis.fsize = psyrat_guifsize;
end

figwidth = 920;
figheight = 560;
fsize = psyrat_prefs.guis.fsize;

hFig = psyrat_newfigure('unit','pix','Visible','off', ...
    'position',[250 150 figwidth figheight], ...
    'menub','no','numbertitle','off','resize','off', ...
    'Name','PsyRAT: Import + Score ERP');
movegui(hFig,'center');
hFig.Tag = 'psyrat_import_gui';

state = struct();
state.project = [];
state.project_file = '';
state.last_import_report = [];
state.last_score_report = [];
state.last_handoff = [];
guidata(hFig,state);

uicontrol(hFig,'Style','text','fontsize',fsize+3,'HorizontalAlignment','left', ...
    'String','PsyRAT: Import + Score ERP (single-trial authoritative workflow)', ...
    'Position',[20 figheight-45 figwidth-40 30]);

uicontrol(hFig,'Style','text','fontsize',fsize,'HorizontalAlignment','left', ...
    'String',['Safety policy: scoring/handoff use single-trial data only. ',...
    'Imported averaged ERP data are stored as reference and flagged in warnings.'], ...
    'Position',[20 figheight-75 figwidth-40 22]);

statusBox = uicontrol(hFig,'Style','listbox','fontsize',fsize,'HorizontalAlignment','left', ...
    'Max',2,'Min',0,'Position',[20 90 figwidth-320 figheight-190], ...
    'String',{'Ready. Create or load a .psyratproj project to begin.'});

uicontrol(hFig,'Style','push','fontsize',fsize,'String','New Project', ...
    'Position',[figwidth-270 figheight-110 240 34], ...
    'TooltipString','Create a new resumable .psyratproj project before importing data.', ...
    'Callback',@onNewProject);

uicontrol(hFig,'Style','push','fontsize',fsize,'String','Open Project', ...
    'Position',[figwidth-270 figheight-150 240 34], ...
    'TooltipString','Open an existing .psyratproj file to continue a prior import/scoring session.', ...
    'Callback',@onOpenProject);

uicontrol(hFig,'Style','push','fontsize',fsize,'String','Save Project', ...
    'Position',[figwidth-270 figheight-190 240 34], ...
    'TooltipString','Save the current import/scoring state so work can be resumed later.', ...
    'Callback',@onSaveProject);

uicontrol(hFig,'Style','push','fontsize',fsize,'String','Import Source', ...
    'Position',[figwidth-270 figheight-230 240 34], ...
    'TooltipString','Import one ERP source file (.set, .erp, .mat, .ept) into the current project.', ...
    'Callback',@onImportSource);

uicontrol(hFig,'Style','push','fontsize',fsize,'String','Add Scoring Spec', ...
    'Position',[figwidth-270 figheight-270 240 34], ...
    'TooltipString','Define a component measurement method and parameters for single-trial scoring.', ...
    'Callback',@onAddSpec);

uicontrol(hFig,'Style','push','fontsize',fsize,'String','Run Scoring', ...
    'Position',[figwidth-270 figheight-310 240 34], ...
    'TooltipString','Compute scores for the active spec; failed trials are excluded from handoff.', ...
    'Callback',@onRunScoring);

uicontrol(hFig,'Style','push','fontsize',fsize,'String','Export Handoff Table', ...
    'Position',[figwidth-270 figheight-350 240 34], ...
    'TooltipString','Build a validated process-data table from included scored trials and send it to workspace.', ...
    'Callback',@onExportHandoff);

uicontrol(hFig,'Style','push','fontsize',fsize,'String','Browse Trials', ...
    'Position',[figwidth-270 figheight-390 240 34], ...
    'TooltipString',['Open a trial browser (next/prev/slider) to cycle through many trials ',...
    'and inspect event/cell metadata and scores.'], ...
    'Callback',@onBrowseTrials);

uicontrol(hFig,'Style','push','fontsize',fsize,'String','Cells / Events Help', ...
    'Position',[figwidth-270 figheight-430 240 34], ...
    'TooltipString','Explain how event/cell labels are defined from imported source data.', ...
    'Callback',@onCellHelp);

uicontrol(hFig,'Style','push','fontsize',fsize,'String','Back to Home', ...
    'Position',[figwidth-270 95 240 40], ...
    'TooltipString','Return to the PsyRAT start screen.', ...
    'Callback',@(~,~) localBackToStart());

set(hFig,'Visible','on');

    function onNewProject(~,~)
        answer = inputdlg({'Project name:'},'New PsyRAT Project',1,{'Untitled PsyRAT ERP Project'});
        if isempty(answer)
            return;
        end
        st = guidata(hFig);
        st.project = psyrat_project_new('name',answer{1});
        st.project_file = '';
        st.last_import_report = [];
        st.last_score_report = [];
        st.last_handoff = [];
        guidata(hFig,st);
        appendStatus(sprintf('Created project: %s',st.project.name));
    end

    function onOpenProject(~,~)
        [f,p] = uigetfile({'*.psyratproj','PsyRAT Project (*.psyratproj)'},'Open PsyRAT Project');
        if isequal(f,0)
            return;
        end
        infile = fullfile(p,f);
        try
            st = guidata(hFig);
            st.project = psyrat_project_load('file',infile);
            st.project_file = infile;
            guidata(hFig,st);
            appendStatus(sprintf('Opened project: %s',infile));
        catch ME
            errordlg({ME.message},'Open Project Error','modal');
        end
    end

    function onSaveProject(~,~)
        st = guidata(hFig);
        if isempty(st.project)
            errordlg({'No project is loaded.'; ...
                'How to fix: click New Project or Open Project first.'},'Save Error','modal');
            return;
        end

        outfile = st.project_file;
        if isempty(outfile)
            [f,p] = uiputfile({'*.psyratproj','PsyRAT Project (*.psyratproj)'}, ...
                'Save PsyRAT Project', [st.project.project_id '.psyratproj']);
            if isequal(f,0)
                return;
            end
            outfile = fullfile(p,f);
        end

        try
            psyrat_project_save('project',st.project,'file',outfile);
            st.project_file = outfile;
            guidata(hFig,st);
            appendStatus(sprintf('Saved project: %s',outfile));
        catch ME
            errordlg({ME.message},'Save Project Error','modal');
        end
    end

    function onImportSource(~,~)
        st = guidata(hFig);
        if isempty(st.project)
            errordlg({'No project is loaded.'; ...
                'How to fix: click New Project or Open Project first.'},'Import Error','modal');
            return;
        end

        [f,p] = uigetfile({'*.set;*.erp;*.mat;*.ept','Supported ERP Sources (*.set, *.erp, *.mat, *.ept)'; ...
            '*.set','EEGLAB (*.set)'; '*.erp','ERPLAB (*.erp)'; '*.mat','MATLAB (*.mat)'; '*.ept','EP Toolkit (.ept)'}, ...
            'Import ERP Source(s)','MultiSelect','on');
        if isequal(f,0)
            return;
        end
        % uigetfile returns a char for a single file and a cell for several;
        % normalize to a cell so one loop handles both cases.
        if ischar(f)
            f = {f};
        end

        % Import each selected file in turn. Sequential imports are merged into
        % the canonical single-trial store (psyrat_project_import_file +
        % localMergeCanonical) with strict time-axis/fs/channel compatibility
        % checks, so a per-file failure is reported and skipped rather than
        % aborting the whole batch.
        nFiles = numel(f);
        nImported = 0;
        nFailed = 0;
        for k = 1:nFiles
            infile = fullfile(p,f{k});
            try
                [st.project, rpt] = psyrat_project_import_file('project',st.project,'file',infile);
                st.last_import_report = rpt;
                nImported = nImported + 1;
                appendStatus(sprintf('Imported source (%d/%d): %s',k,nFiles,infile));
                appendStatus(sprintf('Detected format: %s (confidence %.2f)', ...
                    rpt.detect.display_name, rpt.detect.confidence));
                if ~isempty(rpt.warnings)
                    for iWarn = 1:numel(rpt.warnings)
                        appendStatus(['Warning: ' rpt.warnings{iWarn}]);
                    end
                end
            catch ME
                nFailed = nFailed + 1;
                appendStatus(sprintf('FAILED import (%d/%d): %s -- %s',k,nFiles,infile,ME.message));
            end
        end

        guidata(hFig,st);
        appendCellSummary(st.project);
        appendStatus(sprintf('Import batch complete: %d imported, %d failed.',nImported,nFailed));
        if nFailed > 0
            appendStatus(['Note: failed files were skipped; common cause is a ',...
                'time-axis/sampling-rate/channel mismatch with already-imported data.']);
        end
    end

    function onAddSpec(~,~)
        st = guidata(hFig);
        if isempty(st.project)
            errordlg({'No project is loaded.'; ...
                'How to fix: click New Project or Open Project first.'},'Spec Error','modal');
            return;
        end

        methods = psyrat_scoring_methods();   % shared with spec-create validation
        [idx,ok] = listdlg('PromptString','Select scoring method:', ...
            'SelectionMode','single','ListString',methods);
        if ~ok
            return;
        end
        method = methods{idx};

        % Build the form dynamically: common fields, then only the parameter(s)
        % the chosen method actually uses. keys label the answers so they can be
        % mapped back regardless of which fields appear. Single trials are
        % assumed to be already baseline-adjusted upstream, so no baseline field.
        keys    = {'name','window_start_ms','window_end_ms','channel','polarity'};
        prompt  = {'Spec name:','Window start (ms):','Window end (ms):', ...
            'Channel selector (index/label; comma-separated ROI allowed):', ...
            'Polarity (positive/negative/both):'};
        defaults = {'Spec Name','0','500','1','positive'};

        switch method
            case 'local_peak_amplitude'
                keys{end+1} = 'neighborhood';
                prompt{end+1} = 'Neighborhood samples (local peak):';
                defaults{end+1} = '3';
            case 'centroid_latency'
                keys{end+1} = 'weight_mode';
                prompt{end+1} = 'Centroid weight mode (positive_only/negative_only/absolute):';
                defaults{end+1} = 'absolute';
            case 'fractional_area_latency'
                keys{end+1} = 'area_fraction';
                prompt{end+1} = 'Area fraction (0-1; fractional area latency):';
                defaults{end+1} = '0.5';
        end

        answ = inputdlg(prompt,'Scoring Spec',1,defaults);
        if isempty(answ)
            return;
        end

        % Map answers to a struct of string fields for the pure spec builder.
        inputs = struct();
        for iKey = 1:numel(keys)
            inputs.(keys{iKey}) = answ{iKey};
        end

        try
            spec = psyrat_scoring_spec_from_inputs(method,inputs);
            st.project = psyrat_project_add_scoring_spec('project',st.project,'spec',spec);
            guidata(hFig,st);
            appendStatus(sprintf('Added scoring spec: %s (%s)',spec.spec_id,spec.method));
        catch ME
            errordlg({ME.message},'Spec Error','modal');
        end
    end

    function onRunScoring(~,~)
        st = guidata(hFig);
        if isempty(st.project)
            errordlg({'No project is loaded.'; ...
                'How to fix: click New Project or Open Project first.'},'Scoring Error','modal');
            return;
        end
        specId = localChooseSpec(st.project,'Select scoring spec to run:');
        if isempty(specId)
            return;
        end
        try
            st.project.ui_state.active_spec_id = specId;
            [st.project, rpt] = psyrat_project_score('project',st.project,'spec_id',specId);
            st.last_score_report = rpt;
            guidata(hFig,st);
            appendStatus(rpt.message);
        catch ME
            errordlg({ME.message},'Scoring Error','modal');
        end
    end

    function onExportHandoff(~,~)
        st = guidata(hFig);
        if isempty(st.project)
            errordlg({'No project is loaded.'; ...
                'How to fix: click New Project or Open Project first.'},'Handoff Error','modal');
            return;
        end
        specId = localChooseSpec(st.project,'Select scoring spec to export:');
        if isempty(specId)
            return;
        end
        try
            st.project.ui_state.active_spec_id = specId;
            [tbl, meta, st.project] = psyrat_project_export_proc_table('project',st.project, ...
                'spec_id',specId);
            st.last_handoff.proc_table = tbl;
            st.last_handoff.meta = meta;
            guidata(hFig,st);
            assignin('base','psyrat_import_proc_table',tbl);
            assignin('base','psyrat_import_handoff_meta',meta);
            appendStatus(meta.notification);
            appendStatus('Exported process-data table to workspace variable: psyrat_import_proc_table');
        catch ME
            errordlg({ME.message},'Handoff Error','modal');
        end
    end

    function onBrowseTrials(~,~)
        st = guidata(hFig);
        if isempty(st.project)
            errordlg({'No project is loaded.'; ...
                'How to fix: click New Project or Open Project first.'},'Trial Browser Error','modal');
            return;
        end
        if ~isfield(st.project,'canonical') || ~st.project.canonical.has_single_trial
            errordlg({'This project has no single-trial data to browse.'; ...
                'How to fix: import a single-trial ERP source first.'},'Trial Browser Error','modal');
            return;
        end

        try
            psyrat_import_trialbrowser('project',st.project, ...
                'spec_id',st.project.ui_state.active_spec_id);
        catch ME
            errordlg({ME.message},'Trial Browser Error','modal');
        end
    end

    function onCellHelp(~,~)
        dlg = { ...
            'How PsyRAT defines cells/events in this workflow:'; ...
            ''; ...
            '- "event(cell)" means the trial-level condition label stored in canonical trial_table.event_id.'; ...
            '- EEGLAB imports: event(cell) is the time-locking event (epoch eventtype at latency ~0; first listed event used if latency is unavailable).'; ...
            '- EP Toolkit imports: event(cell) comes from EP cell names (e.g., EPdata.cellNames).'; ...
            '- ERPLAB .erp imports are averaged/reference-only (no single-trial event cells for scoring).'; ...
            ''; ...
            'Handoff mapping:'; ...
            '- subject_id -> id'; ...
            '- event_id -> event (if present)'; ...
            '- occasion_id -> time (if present)'; ...
            '- group_id -> group (if present)'; ...
            ''; ...
            'Scoring always uses single-trial data. Averaged data are reference-only.'};
        msgbox(dlg,'Cells / Events Help','help');
    end

    function appendStatus(msg)
        old = get(statusBox,'String');
        if ischar(old)
            old = {old};
        end
        old{end+1} = sprintf('[%s] %s',datestr(now,'HH:MM:SS'),msg);
        set(statusBox,'String',old,'Value',numel(old));
        drawnow;
    end

    function appendCellSummary(project)
        if ~isfield(project,'canonical') || ~project.canonical.has_single_trial
            return;
        end
        tt = project.canonical.trial_table;
        nTrials = height(tt);
        eventVals = string(tt.event_id);
        eventVals = eventVals(strlength(eventVals)>0);
        uev = unique(eventVals,'stable');
        nEvents = numel(uev);
        if nEvents > 0
            preview = strjoin(cellstr(uev(1:min(5,nEvents))),', ');
            if nEvents > 5
                preview = [preview ', ...'];
            end
            appendStatus(sprintf(['Cell/Event definition: %d unique event(cell) labels across %d trials ',...
                '(event_id). Example labels: %s'],nEvents,nTrials,preview));
        else
            appendStatus(sprintf(['Cell/Event definition: event(cell) labels are currently empty ',...
                'for %d trials.'],nTrials));
        end
    end

end

function localBackToStart()
psyrat_import_gui = findobj('Tag','psyrat_import_gui');
if ~isempty(psyrat_import_gui)
    close(psyrat_import_gui);
end
psyrat_start;
end

function localEnsureImportWorkflowPath()
%Add importworkflow directories to path when this module is launched from
%an environment with stale saved paths.

if exist('psyrat_project_new','file') == 2
    return;
end

thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));
importRoot = fullfile(repoRoot,'subroutines','importworkflow');

if exist(importRoot,'dir') == 7
    addpath(genpath(importRoot));
end
end

function specId = localChooseSpec(project,prompt)
%Pick a scoring spec for Run Scoring / Export. Returns '' when cancelled or
%when no specs exist. With a single spec, returns it directly (no dialog) so
%the common case is friction-free; with several, shows a picker preselecting
%the current active spec.

specId = '';
if ~isfield(project,'scoring_specs') || isempty(project.scoring_specs)
    errordlg({'No scoring specs are defined.'; ...
        'How to fix: click Add Scoring Spec first.'},'Spec Selection','modal');
    return;
end

specs = project.scoring_specs;
if numel(specs) == 1
    specId = specs(1).spec_id;
    return;
end

labels = cell(numel(specs),1);
for i = 1:numel(specs)
    labels{i} = sprintf('%s  [%s]',specs(i).name,specs(i).method);
end

preIdx = 1;
if isfield(project,'ui_state') && isfield(project.ui_state,'active_spec_id') ...
        && ~isempty(project.ui_state.active_spec_id)
    hit = find(strcmp({specs.spec_id},project.ui_state.active_spec_id),1,'first');
    if ~isempty(hit)
        preIdx = hit;
    end
end

[idx,ok] = listdlg('PromptString',prompt,'SelectionMode','single', ...
    'ListString',labels,'InitialValue',preIdx);
if ~ok
    return;
end
specId = specs(idx).spec_id;
end
