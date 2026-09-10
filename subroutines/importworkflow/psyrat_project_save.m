function psyrat_project_save(varargin)
%Save PsyRAT ERP import/scoring project to .psyratproj file.
%
% psyrat_project_save('project',project,'file','/tmp/my.psyratproj')

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
outfile = '';

for i = 1:2:length(varargin)
    key = lower(char(string(varargin{i})));
    val = varargin{i+1};
    switch key
        case 'project'
            project = val;
        case 'file'
            outfile = char(string(val));
        otherwise
            error('varargin:unknown','Unknown option ''%s''.',key);
    end
end

if isempty(project)
    error('psyrat_project:missingProject', ...
        'Input ''project'' is required.');
end

if isempty(outfile)
    error('psyrat_project:missingFile', ...
        'Output ''file'' is required.');
end

[folder,~,ext] = fileparts(outfile);
if isempty(ext)
    outfile = [outfile '.psyratproj'];
    [folder,~,~] = fileparts(outfile);
end

if isempty(folder)
    folder = pwd;
    outfile = fullfile(folder,outfile);
end

if exist(folder,'dir') ~= 7
    error('psyrat_project:savePathMissing', ...
        'Output folder does not exist: %s',folder);
end

project.updated_utc = psyrat_iso8601_utc;
report = psyrat_project_validate(project);
if ~report.ok
    error(report.errors(1).id,'%s',report.errors(1).message);
end

save(outfile,'project','-v7.3');

end
