function [removed,reason] = psyrat_safe_rmdepdir(target,label)
%PSYRAT_SAFE_RMDEPDIR Recursively remove an installed dependency directory.
%
%   [removed,reason] = psyrat_safe_rmdepdir(target,label)
%
%   Wraps rmdir(target,'s') with the checks that make it safe to point at a
%   directory resolved from the Matlab path. Refuses rather than errors, so a
%   caller updating several dependents removes the ones it may and reports the
%   ones it may not.
%
%Input
% target - directory to remove, normally fileparts(which(<a file in it>))
% label  - dependency name used in messages (e.g. 'MatlabStan')
%
%Output
% removed - true when the directory was removed
% reason  - '' on success, otherwise why the removal was refused
%
%WHY THIS EXISTS
% rmolddeps in psyrat_installdependents deletes a dependency by resolving one of
% its files on the Matlab path and removing the containing folder:
%
%     rmdir(fileparts(which('mcmc.m')),'s')
%
% That target is whatever sits first on the path, which is NOT necessarily an
% installed copy. PsyRAT adds its own vendored bundle_dependents tree to the path
% (psyrat_start) and persists it with savepath, and the installed copy carries an
% IDENTICAL folder name, so on a normal developer or repository checkout the
% in-repo bundle wins the lookup. Removing it would delete PsyRAT's own vendored
% source - every local patch with it - out of the working tree, and 'stale
% dependency' is never a reason to do that.
%
% The guard is therefore a boundary check, not a confirmation prompt: a directory
% that lives inside the PsyRAT toolbox tree is source, not an install, and is
% never removable here no matter what a version check reported.
%
%SEE ALSO psyrat_installdependents, psyrat_checkversionsofdeps

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

removed = false;

if nargin < 2 || isempty(label)
    label = 'dependency';
end

%An unresolved which() yields '', and fileparts('') is '' - rmdir('','s') must
%never be reached with that.
if isempty(target) || ~ischar(target) && ~isstring(target)
    reason = sprintf(['No directory could be resolved for %s, so nothing ' ...
        'was removed.'],label);
    return;
end
target = char(target);

if exist(target,'dir') ~= 7
    reason = sprintf(['Directory for %s does not exist (%s), so nothing was ' ...
        'removed.'],label,target);
    return;
end

%Refuse a filesystem root or a bare drive letter outright. Cheap, and the one
%mistake that would be unrecoverable.
canon = local_canonical(target);
if isempty(canon) || strcmp(canon,filesep) || ...
        ~isempty(regexp(canon,'^[A-Za-z]:\\?$','once'))
    reason = sprintf(['Refusing to remove %s: the resolved path (%s) is a ' ...
        'filesystem root.'],label,target);
    return;
end

%The check that matters: never delete PsyRAT's own tree.
toolboxroot = local_canonical(local_toolboxroot());
if ~isempty(toolboxroot) && local_is_within(canon,toolboxroot)
    reason = sprintf(['Refusing to remove %s: the resolved path (%s) is ' ...
        'inside the PsyRAT toolbox directory (%s), so it is PsyRAT''s own ' ...
        'vendored copy rather than an installed one. Nothing was removed. ' ...
        'The bundled copy is used directly and does not need reinstalling.'], ...
        label,target,toolboxroot);
    return;
end

rmdir(target,'s');
removed = true;
reason = '';

end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function root = local_toolboxroot()
%LOCAL_TOOLBOXROOT The PsyRAT installation directory, or '' if not resolvable.
%   psyrat_start lives at the toolbox root; this file lives two levels below it.
%   The which() form is preferred because it reports where the RUNNING toolbox
%   is, but it is empty when the root is not on the path, so fall back to this
%   file's own location.

entry = which('psyrat_start');
if ~isempty(entry)
    root = fileparts(entry);
else
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

end

function p = local_canonical(p)
%LOCAL_CANONICAL Absolute path with a single trailing filesep, for prefix tests.
%   Uses what the filesystem reports where possible so that symlinks, '..'
%   segments and relative inputs cannot be used to sidestep the boundary check.

if isempty(p)
    p = '';
    return;
end
p = char(p);

d = dir(p);
if ~isempty(d) && isfield(d,'folder') && ~isempty(d(1).folder)
    %For a directory, dir() reports the parent in .folder and '.' in .name, so
    %rebuilding from those two yields the resolved absolute path.
    if strcmp(d(1).name,'.')
        p = d(1).folder;
    end
end

if ~endsWith(p,filesep)
    p = [p filesep];
end

end

function tf = local_is_within(child,parent)
%LOCAL_IS_WITHIN True when child is parent or sits underneath it.
%   Both arguments must already be canonical with a trailing filesep, so this is
%   a plain prefix test and 'PsyRAT_test_other' cannot match 'PsyRAT_test'.
%   Case-insensitive on Windows and macOS, where the filesystem is too.

if isempty(child) || isempty(parent)
    tf = false;
    return;
end

if ispc || ismac
    tf = strncmpi(child,parent,numel(parent));
else
    tf = strncmp(child,parent,numel(parent));
end

end
