function report = psyrat_project_validate(project)
%Validate PsyRAT ERP import/scoring project schema and invariants.

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

report = struct();
report.ok = true;
report.errors = struct('id',{},'what',{},'howto',{},'message',{});
report.warnings = struct('id',{},'what',{},'howto',{},'message',{});

requiredTop = {'schema_name','schema_version','project_id','name', ...
    'created_utc','updated_utc','psyrat_version','source_manifest', ...
    'raw_source_refs','canonical','scoring_specs','scoring_results', ...
    'validation_reports','audit_log','ui_state'};

for i = 1:numel(requiredTop)
    if ~isfield(project, requiredTop{i})
        report = localAddError(report,'psyrat_project:missingField', ...
            sprintf('Project is missing required field ''%s''.',requiredTop{i}), ...
            'Load a valid .psyratproj file or recreate the project.');
    end
end

if ~isfield(project,'schema_name') || ~strcmp(project.schema_name,'PsyRATImportProject')
    report = localAddError(report,'psyrat_project:schemaName', ...
        'Project schema_name is invalid.', ...
        'Use projects created by psyrat_project_new or migrate legacy files.');
end

if isfield(project,'schema_version') && ~localIsVersion(project.schema_version)
    report = localAddError(report,'psyrat_project:schemaVersion', ...
        sprintf('Project schema_version is invalid: ''%s''.',char(string(project.schema_version))), ...
        'Use semantic versions like ''1.0.0''.');
end

if isfield(project,'canonical')
    c = project.canonical;
    cReq = {'has_single_trial','has_average','single_trial','trial_table'};
    for i = 1:numel(cReq)
        if ~isfield(c,cReq{i})
            report = localAddError(report,'psyrat_project:canonicalField', ...
                sprintf('Canonical data are missing field ''%s''.',cReq{i}), ...
                'Re-import source data into a new project.');
        end
    end

    if isfield(c,'has_single_trial') && c.has_single_trial
        stReq = {'signals','time_ms','fs_hz','channel_labels','data_hash'};
        for i = 1:numel(stReq)
            if ~isfield(c.single_trial,stReq{i})
                report = localAddError(report,'psyrat_project:singleTrialField', ...
                    sprintf('single_trial is missing field ''%s''.',stReq{i}), ...
                    'Re-import source data so canonical single_trial data are complete.');
            end
        end

        if isfield(c.single_trial,'signals') && ~isempty(c.single_trial.signals) && ndims(c.single_trial.signals) ~= 3
            report = localAddError(report,'psyrat_project:singleTrialShape', ...
                'single_trial.signals must be trial x channel x sample.', ...
                'Normalize imported data to 3-D trial/channel/sample format.');
        end

        if isfield(c,'trial_table') && (~istable(c.trial_table) || isempty(c.trial_table))
            report = localAddError(report,'psyrat_project:trialTable', ...
                'trial_table is missing or empty while has_single_trial is true.', ...
                'Ensure trial metadata are generated during import.');
        end
    end
end

if isfield(project,'scoring_results') && ~isempty(project.scoring_results)
    if ~istable(project.scoring_results)
        report = localAddError(report,'psyrat_project:scoringResultsType', ...
            'scoring_results must be a table.', ...
            'Recompute scoring results or recreate project.');
    else
        needed = {'run_id','spec_id','trial_uid','value','status','excluded_from_handoff'};
        vars = project.scoring_results.Properties.VariableNames;
        for i = 1:numel(needed)
            if ~any(strcmp(vars,needed{i}))
                report = localAddError(report,'psyrat_project:scoringResultsField', ...
                    sprintf('scoring_results is missing column ''%s''.',needed{i}), ...
                    'Recompute scores so results table has required fields.');
            end
        end
    end
end

report.ok = isempty(report.errors);

end

function tf = localIsVersion(v)
if ~(ischar(v) || (isstring(v) && isscalar(v)))
    tf = false;
    return;
end
v = char(string(v));
tf = ~isempty(regexp(v,'^\d+\.\d+\.\d+$','once')); %#ok<RGXP1>
end

function report = localAddError(report,id,what,howto)
entry = struct('id',id,'what',what,'howto',howto, ...
    'message',sprintf('What failed: %s\nHow to fix: %s',what,howto));
if isempty(report.errors)
    report.errors = entry;
else
    report.errors(end+1) = entry; %#ok<AGROW>
end
end
