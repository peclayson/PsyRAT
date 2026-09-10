function project = psyrat_project_load(varargin)
%Load PsyRAT ERP import/scoring project from .psyratproj file.
%
% project = psyrat_project_load('file','/tmp/my.psyratproj')

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

infile = '';

for i = 1:2:length(varargin)
    key = lower(char(string(varargin{i})));
    val = varargin{i+1};
    switch key
        case 'file'
            infile = char(string(val));
        otherwise
            error('varargin:unknown','Unknown option ''%s''.',key);
    end
end

if isempty(infile)
    error('psyrat_project:missingFile', ...
        'Input ''file'' is required.');
end

if exist(infile,'file') ~= 2
    error('psyrat_project:fileNotFound', ...
        'Project file was not found: %s', infile);
end

s = load(infile,'-mat');
if ~isfield(s,'project')
    error('psyrat_project:missingVariable', ...
        'Project file does not contain variable ''project'': %s', infile);
end

project = s.project;
project = localMigrateLegacySpecs(project);

report = psyrat_project_validate(project);
if ~report.ok
    error(report.errors(1).id,'%s',report.errors(1).message);
end

localRequireCompatibleSchema(project.schema_version);

end

function localRequireCompatibleSchema(schemaVersion)
parts = sscanf(char(string(schemaVersion)),'%d.%d.%d')';
if numel(parts) ~= 3
    error('psyrat_project:schemaVersion', ...
        'Project schema_version is invalid: %s', char(string(schemaVersion)));
end

current = [1 0 0];
if parts(1) ~= current(1)
    error('psyrat_project:schemaIncompatible', ...
        ['Project schema major version (%d) is incompatible with this PsyRAT ',...
        'import workflow schema (%d).'], parts(1), current(1));
end
end

function project = localMigrateLegacySpecs(project)
% Migrate scoring specs saved by older builds. The opt-in baseline-correction
% feature was removed, so drop the obsolete baseline_policy field. Without this,
% adding a new (field-free) spec to a loaded legacy project would fail
% struct-array concatenation in psyrat_project_add_scoring_spec ("Subscripted
% assignment between dissimilar structures").
if isfield(project,'scoring_specs') && ~isempty(project.scoring_specs) && ...
        isfield(project.scoring_specs,'baseline_policy')
    project.scoring_specs = rmfield(project.scoring_specs,'baseline_policy');
end
end
