function project = psyrat_project_new(varargin)
%Create a new resumable PsyRAT ERP import/scoring project.
%
% project = psyrat_project_new('name','My Project')

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

opts = struct();
opts.name = '';
opts.created_by = '';

if mod(length(varargin),2)
    error('varargin:incomplete', ...
        'Inputs are incomplete. Provide name/value pairs.');
end

for i = 1:2:length(varargin)
    key = lower(char(string(varargin{i})));
    val = varargin{i+1};
    switch key
        case 'name'
            opts.name = char(string(val));
        case 'created_by'
            opts.created_by = char(string(val));
        otherwise
            error('varargin:unknown', ...
                'Unknown option ''%s''.', key);
    end
end

if isempty(opts.name)
    opts.name = 'Untitled PsyRAT ERP Project';
end

if exist('psyrat_defineversion','file') == 2
    psyrat_version = psyrat_defineversion;
else
    psyrat_version = 'unknown';
end

created_utc = psyrat_iso8601_utc;

project = struct();
project.schema_name = 'PsyRATImportProject';
project.schema_version = '1.0.0';
project.project_id = psyrat_uuid('proj');
project.name = opts.name;
project.created_by = opts.created_by;
project.created_utc = created_utc;
project.updated_utc = created_utc;
project.psyrat_version = psyrat_version;
project.matlab_version = version;
project.source_manifest = struct([]);
project.raw_source_refs = struct([]);
project.canonical = localEmptyCanonical();
project.scoring_specs = struct([]);
project.scoring_results = table();
project.validation_reports = struct([]);
project.audit_log = localAuditEntry('project_created', ...
    struct('name',opts.name,'created_by',opts.created_by));
project.ui_state = struct('active_spec_id','', ...
    'active_run_id','', ...
    'active_handoff_spec_id','', ...
    'notes','');

end

function canonical = localEmptyCanonical()
canonical = struct();
canonical.has_single_trial = false;
canonical.has_average = false;
canonical.averaged_reference_origin = '';
canonical.warning_flags = struct('averaged_not_used_for_scoring_handoff',false);
canonical.single_trial = struct( ...
    'signals', [], ...
    'time_ms', [], ...
    'fs_hz', [], ...
    'amplitude_unit', '', ...
    'channel_labels', {{}}, ...
    'data_hash', '');
canonical.averaged_reference = struct( ...
    'signals', [], ...
    'time_ms', [], ...
    'fs_hz', [], ...
    'channel_labels', {{}}, ...
    'source', '', ...
    'data_hash', '');
canonical.trial_table = table();
end

function entry = localAuditEntry(action,details)
entry = struct('timestamp_utc',psyrat_iso8601_utc, ...
    'action',action, ...
    'details',details);
end
