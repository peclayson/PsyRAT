function [proc_table, handoff_meta, project] = psyrat_project_export_proc_table(varargin)
%Export validated process-data table from single-trial scoring results.
%
% [proc_table, handoff_meta, project] = psyrat_project_export_proc_table(...)

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

opts = struct('project',[],'spec_id','','run_id','', ...
    'validate_with_preflight',true,'sserrvar',1,'diffest',1);

if mod(length(varargin),2)
    error('varargin:incomplete','Inputs are incomplete. Provide name/value pairs.');
end

for i = 1:2:length(varargin)
    key = lower(char(string(varargin{i})));
    val = varargin{i+1};
    switch key
        case 'project'
            opts.project = val;
        case 'spec_id'
            opts.spec_id = char(string(val));
        case 'run_id'
            opts.run_id = char(string(val));
        case 'validate_with_preflight'
            opts.validate_with_preflight = logical(val);
        case 'sserrvar'
            opts.sserrvar = val;
        case 'diffest'
            opts.diffest = val;
        otherwise
            error('varargin:unknown','Unknown option ''%s''.',key);
    end
end

if isempty(opts.project)
    error('psyrat_handoff:missingProject','Input ''project'' is required.');
end
project = opts.project;

if isempty(opts.spec_id)
    opts.spec_id = project.ui_state.active_spec_id;
end
if isempty(opts.spec_id)
    error('psyrat_handoff:missingSpec','No active scoring spec selected for handoff.');
end

if isempty(project.scoring_results) || ~istable(project.scoring_results)
    error('psyrat_handoff:noScores', ...
        'Project contains no scoring results to export.');
end

rows = project.scoring_results(strcmp(project.scoring_results.spec_id,opts.spec_id),:);
if isempty(rows)
    error('psyrat_handoff:noScoresForSpec', ...
        'No scoring results found for spec_id %s', opts.spec_id);
end

if any(strcmp(rows.Properties.VariableNames,'is_stale'))
    rows = rows(~rows.is_stale,:);
end

if isempty(rows)
    error('psyrat_handoff:staleOnly', ...
        'All scoring results for spec_id %s are stale. Recompute before handoff.', opts.spec_id);
end

%Default run: prefer the project's active run, but only when it belongs to the
%SELECTED spec's surviving rows. active_run_id is spec-agnostic (scoring ANY
%spec overwrites it), so after scoring spec B a default export of spec A used
%to filter A's rows by B's run id and die runNotFound below (finding G35). The
%membership test routes that case to the per-spec-latest fallback instead --
%rows is already spec- and stale-filtered here, so its last row is the latest
%surviving run of the chosen spec (the same per-spec-latest rule the trial
%browser applies). An explicitly passed 'run_id' is untouched and still errors
%runNotFound when absent, which is the caller asking for that exact run.
if isempty(opts.run_id)
    if isfield(project.ui_state,'active_run_id') && ~isempty(project.ui_state.active_run_id) && ...
            any(strcmp(rows.run_id,project.ui_state.active_run_id))
        opts.run_id = project.ui_state.active_run_id;
    else
        opts.run_id = rows.run_id{end};
    end
end

rows = rows(strcmp(rows.run_id,opts.run_id),:);
if isempty(rows)
    error('psyrat_handoff:runNotFound', ...
        'No scoring results found for run_id %s', opts.run_id);
end

% Protect against hash drift.
if any(~strcmp(rows.data_hash,project.canonical.single_trial.data_hash))
    error('psyrat_handoff:dataHashMismatch', ...
        'Scoring results were generated from different canonical data. Recompute scores before handoff.');
end

trials = project.canonical.trial_table(:,{'trial_uid','subject_id','event_id','occasion_id','group_id'});
joined = outerjoin(rows,trials,'Keys','trial_uid','MergeKeys',true,'Type','left');

excluded_mask = joined.excluded_from_handoff | ~strcmp(joined.status,'ok');
excluded_rows = joined(excluded_mask,:);
included = joined(~excluded_mask,:);

if isempty(included)
    error('psyrat_handoff:noIncludedTrials', ...
        ['All trials were excluded by scoring failures. ',...
        'Adjust scoring parameters and recompute before handoff.']);
end

proc_table = table();
proc_table.id = cellstr(string(included.subject_id));
proc_table.meas = included.value;

if any(strlength(string(included.group_id)) > 0)
    proc_table.group = cellstr(string(included.group_id));
end
if any(strlength(string(included.event_id)) > 0)
    proc_table.event = cellstr(string(included.event_id));
end
if any(strlength(string(included.occasion_id)) > 0)
    proc_table.time = cellstr(string(included.occasion_id));
end

if any(~isfinite(proc_table.meas)) || any(ismissing(proc_table.meas))
    error('psyrat_handoff:missingMeasurements', ...
        ['Final handoff dataset contains missing/non-finite measurements. ',...
        'PsyRAT requires complete meas values for process-data handoff.']);
end

preflight = struct('ok',true,'errors',[],'warnings',[]);
if opts.validate_with_preflight
    preflight = psyrat_preflight_validate('datatable',proc_table, ...
        'sserrvar',opts.sserrvar,'diffest',opts.diffest);
    if ~preflight.ok
        error(preflight.errors(1).id,'%s',preflight.errors(1).message);
    end
end

handoff_meta = struct();
handoff_meta.project_id = project.project_id;
handoff_meta.spec_id = opts.spec_id;
handoff_meta.run_id = opts.run_id;
handoff_meta.generated_utc = psyrat_iso8601_utc;
handoff_meta.n_input_trials = height(joined);
handoff_meta.n_included_trials = height(included);
handoff_meta.n_excluded_trials = height(excluded_rows);
handoff_meta.exclude_reasons = localCountReasons(excluded_rows.status_reason);
handoff_meta.notification = sprintf([ ...
    'Handoff dataset prepared with %d included trials and %d excluded trials. ',...
    'Excluded trials were removed from final dataset and logged.'], ...
    handoff_meta.n_included_trials,handoff_meta.n_excluded_trials);
handoff_meta.preflight = preflight;

project.ui_state.active_handoff_spec_id = opts.spec_id;
project.updated_utc = handoff_meta.generated_utc;
project.audit_log = localAppendAudit(project.audit_log,'handoff_exported',struct( ...
    'spec_id',opts.spec_id, ...
    'run_id',opts.run_id, ...
    'n_included',handoff_meta.n_included_trials, ...
    'n_excluded',handoff_meta.n_excluded_trials));

end

function counts = localCountReasons(reasons)
counts = struct('reason',{{}},'count',[]);
if isempty(reasons)
    return;
end

vals = string(reasons(:));
vals = vals(strlength(vals)>0);
if isempty(vals)
    return;
end

u = unique(vals,'stable');
counts.reason = cellstr(u);
counts.count = zeros(numel(u),1);
for i = 1:numel(u)
    counts.count(i) = sum(vals == u(i));
end
end

function out = localAppendAudit(in,action,details)
entry = struct('timestamp_utc',psyrat_iso8601_utc, ...
    'action',action, ...
    'details',details);
if isempty(in)
    out = entry;
else
    out = in;
    out(end+1) = entry;
end
end
