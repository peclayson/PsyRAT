function project = psyrat_project_add_scoring_spec(varargin)
%Add a validated scoring spec to project state.
%
% project = psyrat_project_add_scoring_spec('project',project,'spec',spec)

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

if mod(length(varargin),2)
    error('varargin:incomplete','Inputs are incomplete. Provide name/value pairs.');
end

project = [];
spec = [];

for i = 1:2:length(varargin)
    key = lower(char(string(varargin{i})));
    val = varargin{i+1};
    switch key
        case 'project'
            project = val;
        case 'spec'
            spec = val;
        otherwise
            error('varargin:unknown','Unknown option ''%s''.',key);
    end
end

if isempty(project)
    error('psyrat_scoring:missingProject','Input ''project'' is required.');
end
if isempty(spec)
    error('psyrat_scoring:missingSpec','Input ''spec'' is required.');
end

if ~isstruct(spec) || ~isfield(spec,'spec_id') || ~isfield(spec,'spec_hash')
    error('psyrat_scoring:invalidSpec', ...
        'Spec must be created with psyrat_scoring_spec_create.');
end

existing = {};
if ~isempty(project.scoring_specs)
    existing = {project.scoring_specs.spec_id};
end
if any(strcmp(existing,spec.spec_id))
    error('psyrat_scoring:duplicateSpec', ...
        'Spec ID already exists in project: %s', spec.spec_id);
end

if isempty(project.scoring_specs)
    project.scoring_specs = spec;
else
    project.scoring_specs(end+1) = spec;
end

project.ui_state.active_spec_id = spec.spec_id;
project.updated_utc = psyrat_iso8601_utc;
project.audit_log = localAppendAudit(project.audit_log,'scoring_spec_added', ...
    struct('spec_id',spec.spec_id,'name',spec.name,'method',spec.method));

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
