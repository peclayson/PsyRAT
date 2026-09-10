function info = psyrat_matlabstan_bundleinfo(msdir)
%PSYRAT_MATLABSTAN_BUNDLEINFO Identify the MatlabStan that is actually loaded.
%
%   info = psyrat_matlabstan_bundleinfo()
%   info = psyrat_matlabstan_bundleinfo(msdir)
%
%   Reads the PSYRAT_PATCHLEVEL.txt stamp that PsyRAT ships inside its vendored
%   MatlabStan bundle, and reports whether the copy on the MATLAB path is the
%   one this release of PsyRAT expects.
%
%Input
% msdir - OPTIONAL directory to inspect. When omitted, the directory is resolved
%         from the MatlabStan that MATLAB would actually load.
%
%Output
% info - struct describing the resolved bundle
%   dir        - directory inspected ('' when no MatlabStan could be found)
%   stampfile  - full path to the stamp that was looked for ('' when no dir)
%   found      - true when the stamp exists and both fields parsed
%   version    - upstream MatlabStan release recorded in the stamp ('' if none)
%   patchlevel - PsyRAT patch level recorded in the stamp ('' if none)
%   status     - verdict, extending psyrat_checkversionsofdeps' scale:
%                 0 - old (stamp missing, or older than expected) -- reported
%                     only for copies the containment test does not place
%                     inside the resolved toolbox root
%                 1 - up to date
%                 2 - newer than the version PsyRAT was tested against
%                 3 - the bundled pinned copy with a would-be-old stamp
%                     (G56d): the copy lies inside the PsyRAT toolbox tree,
%                     which the dependents-update flow can never replace
%                     (psyrat_safe_rmdepdir refuses it by design), so the
%                     every-launch update offer is suppressed and
%                     psyrat_start prints the pinned-copy console line
%                     instead; only a toolbox update can change the bundle
%                -1 - could not determine (no MatlabStan on the path)
%   detail     - one-line human-readable explanation of status
%
%WHY THIS ANCHORS ON which('StanModel.m') RATHER THAN ON THE STAMP NAME
% Two MatlabStan copies can sit on the path at once: the in-repo bundle added by
% psyrat_start, and an installed copy added by psyrat_installdependents. Looking
% the stamp up BY NAME would return whichever copy comes first on the path, so a
% freshly patched in-repo stamp could answer on behalf of a stale installed copy
% and report "up to date" while the unpatched files are the ones being loaded.
% Anchoring on StanModel.m finds the copy MATLAB will actually execute, and the
% stamp is then read from that same directory by absolute path. The stamp is a
% plain text file specifically so that a name lookup is not even available.
%
%SEE ALSO psyrat_checkversionsofdeps, psyrat_dependentsversions

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

info = struct('dir','','stampfile','','found',false,'version','', ...
    'patchlevel','','status',-1,'detail','');

%Resolve the MatlabStan that MATLAB would load, unless the caller named one.
if nargin < 1 || isempty(msdir)
    stanmodel = which('StanModel.m');
    if isempty(stanmodel)
        info.detail = 'No MatlabStan found on the Matlab path.';
        return;
    end
    msdir = fileparts(stanmodel);
end
info.dir = msdir;
info.stampfile = fullfile(msdir,'PSYRAT_PATCHLEVEL.txt');

%G56d. A copy living inside the PsyRAT toolbox tree is the vendored bundle:
%psyrat_safe_rmdepdir refuses to remove it by design, so the dependents-update
%flow can never act on an "old" verdict for it. Before this check, a stale (or
%missing) stamp on the bundle re-raised the launch update offer every session
%with an offer that could not succeed. Every would-be-old verdict below is
%therefore reported as status 3 (bundled, pinned) when the inspected directory
%is in-tree -- a DISTINCT value rather than a bare 1, because the launch
%console is the only surface a user sees and "version is current" would be
%false for a stale stamp; psyrat_start prints a dedicated line for 3. An
%out-of-tree stale copy still reports 0 -- that offer IS actionable, because
%its removal is permitted. The containment test mirrors psyrat_safe_rmdepdir's
%boundary check so the verdict is keyed on the remover's own condition (the
%mirror is a copy, welded byte-identical by TestMatlabStanPatchLevel; it lacks
%the remover's exist() early-out, which cannot matter for the which()-resolved
%directories this function inspects).
intree = local_intree(msdir);

depvers = psyrat_dependentsversions;

%A MatlabStan directory with no stamp predates stamping, so it cannot contain
%PsyRAT's local patches. That is "old", not "unknown" - reporting it as unknown
%would leave the user with no signal at all, which is the gap this closes.
if exist(info.stampfile,'file') ~= 2
    if intree
        info.status = 3;
        info.detail = sprintf(['MatlabStan at %s is PsyRAT''s bundled pinned ' ...
            'copy; its patch-level stamp is missing (reinstalling the ' ...
            'toolbox restores it), but a dependents update cannot act on ' ...
            'the bundle.'],msdir);
    else
        info.status = 0;
        info.detail = sprintf(['MatlabStan at %s carries no PsyRAT patch-level ' ...
            'stamp, so it predates PsyRAT''s local patches.'],msdir);
    end
    return;
end

stamp = fileread(info.stampfile);
info.version    = local_read_field(stamp,'upstream_version');
info.patchlevel = local_read_field(stamp,'patch_level');

%A stamp we cannot parse is worse than no stamp, because it suggests the file was
%edited by hand. Say so rather than guessing a verdict from a partial read.
if isempty(info.version) || isempty(info.patchlevel)
    info.status = -1;
    info.detail = sprintf(['MatlabStan stamp at %s could not be parsed ' ...
        '(expected an upstream_version and a patch_level line).'],info.stampfile);
    return;
end
info.found = true;

%The upstream release dominates: a different MatlabStan release is a bigger
%difference than a patch level within one release, and the patch levels of two
%different releases are not comparable.
vercmp = local_compare_dotted(depvers.matlabstan,info.version);
if vercmp ~= 0
    info.status = local_verdict(vercmp);
    info.detail = sprintf('MatlabStan release %s found; PsyRAT expects %s.', ...
        info.version,depvers.matlabstan);
    if info.status == 0 && intree
        %G56d: would-be-old, but it is the pinned bundle (see above).
        info.status = 3;
        info.detail = sprintf(['MatlabStan release %s is PsyRAT''s bundled ' ...
            'pinned copy (PsyRAT expects %s; update the PsyRAT toolbox, ' ...
            'not the dependents).'],info.version,depvers.matlabstan);
    end
    return;
end

patchcmp = local_compare_dotted(depvers.matlabstanpatch,info.patchlevel);
info.status = local_verdict(patchcmp);
switch info.status
    case 0
        if intree
            %G56d: would-be-old, but it is the pinned bundle (see above).
            info.status = 3;
            info.detail = sprintf(['MatlabStan %s patch level %s is ' ...
                'PsyRAT''s bundled pinned copy (PsyRAT expects patch level ' ...
                '%s; update the PsyRAT toolbox, not the dependents).'], ...
                info.version,info.patchlevel,depvers.matlabstanpatch);
        else
            info.detail = sprintf(['MatlabStan %s patch level %s found; PsyRAT ' ...
                'expects patch level %s.'],info.version,info.patchlevel, ...
                depvers.matlabstanpatch);
        end
    case 2
        info.detail = sprintf(['MatlabStan %s patch level %s found; PsyRAT ' ...
            'was tested against patch level %s.'],info.version, ...
            info.patchlevel,depvers.matlabstanpatch);
    otherwise
        info.detail = sprintf('MatlabStan %s patch level %s is current.', ...
            info.version,info.patchlevel);
end

end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function value = local_read_field(text,name)
%LOCAL_READ_FIELD Pull "name = value" out of the stamp.
%   Anchored to the start of a line so that prose elsewhere in the stamp - which
%   discusses these field names by name - cannot be mistaken for the payload.
%   'once' makes the first matching line authoritative if the file is malformed.

value = '';
pattern = ['(?m)^\s*' name '\s*=\s*(\S+)\s*$'];
tok = regexp(text,pattern,'tokens','once');
if ~isempty(tok)
    value = strtrim(tok{1});
end

end

function cmp = local_compare_dotted(expected,found)
%LOCAL_COMPARE_DOTTED Compare two dot-separated numeric version strings.
%   Returns +1 when expected is HIGHER than found (i.e. found is old), -1 when
%   expected is lower than found, and 0 when they are equal.
%
%   Written to take an arbitrary number of components rather than the fixed
%   three the CmdStan and MatlabProcessManager checks assume: MatlabStan's
%   release string is '2.15.1.0', which has four. Missing trailing components
%   are treated as zero, so '2.15.1' and '2.15.1.0' compare equal.

e = local_parts(expected);
f = local_parts(found);
n = max(numel(e),numel(f));
e(end+1:n) = 0;
f(end+1:n) = 0;

cmp = 0;
for ii = 1:n
    if e(ii) > f(ii)
        cmp = 1;
        return;
    elseif e(ii) < f(ii)
        cmp = -1;
        return;
    end
end

end

function parts = local_parts(str)
%LOCAL_PARTS Dot-separated string to a numeric row vector.
%   Non-numeric components become 0 rather than erroring, so a hand-edited stamp
%   degrades to a comparison rather than an exception in a startup path.

parts = zeros(1,0);
if isempty(str)
    return;
end
chunks = strsplit(strtrim(char(str)),'.');
parts = zeros(1,numel(chunks));
for ii = 1:numel(chunks)
    v = str2double(chunks{ii});
    if ~isnan(v)
        parts(ii) = v;
    end
end

end

function status = local_verdict(cmp)
%LOCAL_VERDICT Map a comparison to the 0/1/2 scale the toolbox already uses.

if cmp > 0
    status = 0;      %expected is higher: what is installed is old
elseif cmp < 0
    status = 2;      %expected is lower: newer than PsyRAT was tested against
else
    status = 1;
end

end

function tf = local_intree(msdir)
%LOCAL_INTREE True when msdir sits inside the PsyRAT toolbox tree.
%   Same semantics as psyrat_safe_rmdepdir's boundary check, which is the
%   remover this verdict must agree with: resolve the running toolbox root
%   (which('psyrat_start'), falling back to this file's own location two
%   levels below the root), canonicalize both paths through the filesystem,
%   then prefix-test with a trailing filesep -- case-insensitively where the
%   filesystem is. A false negative keeps the old nag for one setup; a false
%   positive hides an actionable update, so the test stays a strict mirror of
%   the remover's rather than a looser heuristic.

tf = false;

entry = which('psyrat_start');
if ~isempty(entry)
    root = fileparts(entry);
else
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

root = local_canonical(root);
target = local_canonical(msdir);
if isempty(root) || isempty(target)
    return;
end

if ispc || ismac
    tf = strncmpi(target,root,numel(root));
else
    tf = strncmp(target,root,numel(root));
end

end

function p = local_canonical(p)
%LOCAL_CANONICAL Absolute path with a single trailing filesep, for prefix tests.
%   Mirrors psyrat_safe_rmdepdir's helper: uses what the filesystem reports
%   where possible so that symlinks, '..' segments and relative inputs cannot
%   sidestep the containment test.

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
