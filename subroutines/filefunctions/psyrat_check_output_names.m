function [problems, leftovers, notes] = psyrat_check_output_names(savepath, savename)
%Check, before estimation, that every file a run derives from its output name
%can be created in the save folder (row 31, LA-2).
%
%[problems, leftovers, notes] = psyrat_check_output_names(savepath, savename)
%
%A run writes three files named after SAVENAME: the output itself (saved by
%psyrat_computevarcompwarp; like MATLAB's save, a name without an extension
%gets .mat), the run-configuration sidecar <stem>_runconfig.json
%(psyrat_write_runconfig), and, only if the run fails, its error report
%<stem>_error.txt (psyrat_startproc's report writer). The first write happens
%after the fit, so a name the folder cannot hold used to cost the whole fit.
%This function answers the question beforehand, from the actual volume, with
%no fixed byte or character cap: the limit differs between filesystems. On the
%one APFS volume where it was measured (a no-fit probe, 2026-09-23), every
%accepted name had at most 255 UTF-16 units and every refused one at least
%256, for one- to four-byte UTF-8 characters alike: three-byte names of 735 to
%751 bytes were accepted, and four-byte names (two units each) were refused one
%character past 131 to 135 code points. That volume's limit was therefore not
%255 UTF-8 bytes or 255 code points. Other volumes may differ, and none of
%these figures is used anywhere in this function.
%
%Input
% savepath - the existing save folder
% savename - the output name as the run will use it (no path separators)
%
%Output
% problems  - struct array with fields id ('psyrat_data:savename'), what,
%   howto and component; empty when nothing was found. At most one problem is
%   reported (the first component that fails). It states that estimation for
%   the requested output was refused, names the failing file, and gives the
%   remedy.
% leftovers - full paths of test entries (folders or files) this function
%   created but did not remove. Each one is also announced with warning
%   'psyrat:outputNameProbeLeft'. Normally empty.
% notes     - one entry per component examined, with fields component, method
%   ('none' for an existing entry, 'directory', or 'file' when a refused test
%   folder was followed by a test file), outcome ('exists', 'folder',
%   'created', 'refused' or 'indeterminate') and reason. It says what evidence
%   each decision rests on; callers may ignore it.
%
%What counts as a problem
% - A FOLDER already occupies the output or the sidecar path. No file can be
%   written there. The folder and its contents are left as they are.
% - A file with that name cannot be created in the save folder although the
%   folder itself accepts new entries (typically: the name is too long).
%What does NOT count, on purpose
% - An existing FILE at a component path. The name is representable, and the
%   callers' own rules apply unchanged: the save dialog's consent to overwrite
%   the output, the sidecar overwrite, the batch allocator's reservations and
%   the report writer's free-name suffixing.
% - A folder at the report path. The report writer already moves on to a free
%   or timestamped name (and then to the temp folder).
% - A save folder that refuses new entries altogether (read-only, no
%   permission). That is not a name problem. It is left to the run's existing
%   handling of that folder, which reports the failure through the usual error
%   path, so nothing new is refused for it.
% - A name containing a path separator, or '.' or '..'. The existing calling
%   contract is kept, nothing is checked, and no other folder is touched.
% - A case the check cannot decide (no disjoint test name could be built,
%   every test name tried was already taken, or a test folder was refused and
%   no test file can confirm it: no Java runtime, the save folder's absolute
%   location cannot be read, or MATLAB does not see the test file under its
%   name). The name is
%   not refused, the
%   reason is in NOTES, and the run keeps its existing late handling.
%
%How the check works, and what it never does
% For each component that does not exist yet, a test entry of EXACTLY the
% component's length is created in the save folder and removed again. Its name
% keeps the user's characters verbatim and replaces only an ASCII tail (the
% fixed suffix, extended over trailing ASCII characters when it is shorter than
% 7) with '.' plus hexadecimal characters of the same length. The length is
% therefore identical in UTF-8 bytes, UTF-16 units and code points, and the
% test entry is never the target itself.
% - The test entry is a FOLDER first. It is ours only if mkdir reports success
%   with no message identifier (an existing folder answers
%   'MATLAB:MKDIR:DirectoryExists' and is never treated as ours).
% - A refused test folder is not taken as proof. MATLAB's mkdir refuses some
%   long multibyte paths ('MATLAB:FILE:PATHMAX'; on the one macOS machine where
%   it was measured, from about 470 to 510 bytes of full path depending on the
%   folder) that save and fopen still write, so the refusal is
%   then checked with a test FILE of the same name, created by
%   java.io.File.createNewFile. That call is documented as an atomic
%   check-and-create that answers false when the name already exists, so a
%   file it reports as created is ours when it is created; it does not stop
%   another process from touching the file afterwards. The name is refused only
%   if that file cannot be created either.
% - A name that happens to exist is retried, at most 3 times, with a new tail
%   of the same length (no suffix growth).
% - A test folder is removed with a NON-recursive rmdir. A test file is removed
%   only while it is still empty (content in it can only come from someone
%   else). Either way the entry is checked to be gone afterwards; if it is
%   still there, its path goes to LEFTOVERS and to a visible warning, and
%   removal is never claimed without that check.
% - When a component cannot be created, a short control test FILE is made
%   once, the same way. If that fails too, the folder refuses new files
%   altogether or the file test cannot reach it; neither is shown to be a
%   name problem (see above), so nothing is refused.
% - The hexadecimal characters come from tempname's UUID. The global random
%   number generator is not touched, so fit seeds are unaffected.
% - Nothing existing is opened, overwritten, truncated or deleted.
% A passing check is evidence at this moment only, not a promise: a disk that
% fills or a folder whose permissions change before the save can still fail.

% Copyright (C) 2016-2026 Peter E. Clayson
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

problems = struct('id', {}, 'what', {}, 'howto', {}, 'component', {});
leftovers = {};
notes = struct('component', {}, 'method', {}, 'outcome', {}, 'reason', {});

savepath = char(string(savepath));
savename = char(string(savename));

%nothing to check (the caller's own save-path and save-name checks report
%these), or a name that is not a single component: keep the existing contract
if isempty(savename) || isempty(savepath) || ~isfolder(savepath) || ...
        any(savename == '/') || any(savename == '\') || ...
        any(strcmp(savename, {'.', '..'}))
    return;
end

%the three names, derived exactly as their writers derive them
[~, stem, ext] = fileparts(savename);
if isempty(ext)
    outname = [savename '.mat'];    %save appends .mat to a bare name
    outsuffix = '.mat';
else
    outname = savename;
    outsuffix = ext;
end
reportstem = stem;                  %psyrat_write_measurement_error falls back
if isempty(reportstem)              %to the whole name when the stem is empty
    reportstem = savename;
end
%component name, the ASCII suffix the writer appends, the role for messages,
%and whether an existing folder at that path blocks the write
comps = { ...
    outname,                         outsuffix,         'output',                    true; ...
    [stem '_runconfig.json'],        '_runconfig.json', 'run-configuration sidecar', true; ...
    [reportstem '_error.txt'],       '_error.txt',      'error-report',              false};

controlResult = [];   %the short control entry is tried at most once

for c = 1:size(comps, 1)
    component = comps{c, 1};
    target = fullfile(savepath, component);

    if isfolder(target)
        notes(end + 1) = localNote(component, 'none', 'folder', ''); %#ok<AGROW>
        if comps{c, 4}
            problems(end + 1) = localProblem(savename, savepath, component, ...
                sprintf('a folder named ''%s'' already occupies the name of its %s file', ...
                component, comps{c, 3}), ...
                ['Choose a different output name, or move or rename that folder, then ' ...
                'rerun. The folder was not changed, and nothing was estimated for this output.']); %#ok<AGROW>
            return;
        end
        continue;   %report path: the report writer moves on by itself
    end
    if isfile(target)
        notes(end + 1) = localNote(component, 'none', 'exists', ''); %#ok<AGROW>
        continue;   %an existing file: the name is representable here
    end

    [outcome, cause, left, method] = localProbe(savepath, component, comps{c, 2});
    leftovers = [leftovers, left]; %#ok<AGROW>
    if ~strcmp(outcome, 'refused')
        notes(end + 1) = localNote(component, method, outcome, cause); %#ok<AGROW>
        continue;   %creatable, or no evidence either way: no restriction
    end

    %the name could not be created: a name problem only if a short test
    %file, made the same way, can be
    if isempty(controlResult)
        [controlResult, left] = localControl(savepath);
        leftovers = [leftovers, left]; %#ok<AGROW>
    end
    if ~controlResult
        notes(end + 1) = localNote(component, method, 'indeterminate', ...
            ['a short test file made the same way could not be created either (' ...
            cause '), so this is not shown to be a name problem']); %#ok<AGROW>
        return;     %no name problem shown: nothing is refused
    end
    notes(end + 1) = localNote(component, method, 'refused', cause); %#ok<AGROW>
    problems(end + 1) = localProblem(savename, savepath, component, ...
        sprintf('its %s file name ''%s'' (%d characters) cannot be created in that folder (%s)', ...
        comps{c, 3}, component, numel(component), cause), ...
        ['Choose a shorter output name, or shorten the measurement column header ' ...
        'the name is built from, then rerun. Nothing was estimated for this output.']); %#ok<AGROW>
    return;
end

end

function p = localProblem(savename, savepath, component, why, howto)
%The lead states the refusal of the REQUESTED output, not a save failure: when
%only the sidecar or report name is the problem, the output file itself could
%have been written, so "the output cannot be saved" would be false.
p = struct('id', 'psyrat_data:savename', ...
    'what', sprintf('Estimation for the requested output ''%s'' in ''%s'' was refused before the fit: %s.', ...
    savename, savepath, why), ...
    'howto', howto, 'component', component);
end

function n = localNote(component, method, outcome, reason)
n = struct('component', component, 'method', method, 'outcome', outcome, 'reason', reason);
end

function [outcome, cause, left, method] = localProbe(savepath, component, suffix)
%Test COMPONENT with an equal-length, disjoint test entry. outcome is
%'created', 'refused' or 'indeterminate' (no disjoint name could be built, or
%every name tried was taken, or the refusal could not be confirmed).
outcome = 'indeterminate';
cause = '';
left = {};
method = 'none';

%the replaced region: the writer's ASCII suffix, extended over trailing
%ASCII characters of the user's part until it holds 7 characters
keep = component(1:end - numel(suffix));
repl = suffix;
while numel(repl) < 7 && ~isempty(keep) && keep(end) < 128
    repl = [keep(end) repl]; %#ok<AGROW>
    keep = keep(1:end - 1);
end
if numel(repl) < 7 || any(repl >= 128)
    cause = 'no ASCII tail of 7 characters to replace, so no disjoint test name';
    return;
end

for attempt = 1:3
    probe = [keep '.' localHex(numel(repl) - 1)];
    if strcmpi(probe, component)
        continue;   %never the target itself (cannot happen; guarded anyway)
    end
    [outcome, cause, left, method] = localCreate(savepath, probe);
    if ~strcmp(outcome, 'exists')
        return;
    end
end
%every attempt met an existing name: no evidence either way, no restriction
outcome = 'indeterminate';
cause = 'every test name tried was already taken';
end

function [ok, left] = localControl(savepath)
%Whether the save folder accepts a short test FILE made exactly the way the
%refused one was (a refusal always rests on the file test). Success shows
%that the folder accepts new files AND that the Java route reaches the save
%folder under the name MATLAB uses, so the refusal is about the name. A
%short test FOLDER would show neither for the Java route.
ok = false;
left = {};
for attempt = 1:3
    name = ['psyrat_np_' localHex(8)];
    [outcome, ~, left] = localCreateFile(fullfile(savepath, name), ...
        localJavaPath(savepath, name), 'not tried for the control');
    if strcmp(outcome, 'exists')
        continue;
    end
    ok = strcmp(outcome, 'created');
    return;
end
end

function [outcome, cause, left, method] = localCreate(savepath, name)
%Create NAME as a test entry, record the outcome, and remove it again.
%'exists' means something else already carries NAME: it is never ours.
left = {};
cause = '';
method = 'directory';
p = fullfile(savepath, name);
[ok, msg, msgid] = mkdir(savepath, name);
if ok && isempty(msgid)
    %created by this call: it is ours, so remove it and prove it is gone
    outcome = 'created';
    left = localRemove(p);
    return;
end
if ok || isfolder(p) || isfile(p)
    outcome = 'exists';     %something already carries this test name
    return;
end
%The test folder was refused. That is not proof the writers' FILE would be:
%MATLAB's mkdir refuses some long multibyte paths that save and fopen write.
dircause = strtrim(msg);
if isempty(dircause)
    dircause = 'the filesystem refused the name';
end
method = 'file';
[outcome, cause, left] = localCreateFile(p, localJavaPath(savepath, name), dircause);
end

function [outcome, cause, left] = localCreateFile(p, javapath, dircause)
%Confirm a refused test folder with an equal-length test FILE, created with
%java.io.File.createNewFile: an atomic check-and-create that answers false
%when the name exists, so the file is ours only if this call created it. A
%file open ('w' or append) cannot give that guarantee: it may reuse an empty
%file another process created a moment earlier.
left = {};
cause = '';
if ~usejava('jvm')
    outcome = 'indeterminate';
    cause = sprintf(['the test folder was refused (%s), and without a Java runtime no ' ...
        'test file can confirm whether a file of that name can be created'], dircause);
    return;
end
if isempty(javapath)
    outcome = 'indeterminate';
    cause = sprintf(['the test folder was refused (%s), and the save folder''s ' ...
        'absolute location could not be read to address a test file'], dircause);
    return;
end
f = java.io.File(javapath);
try
    made = f.createNewFile();
catch javaME
    outcome = 'refused';
    cause = localJavaCause(javaME, dircause);
    return;
end
if ~made
    outcome = 'exists';     %something already carries this name: never ours
    return;
end
%this call created the file, so it is ours; MATLAB must see it under the same
%name, or the test would not be about the name the writers will use
if ~isfile(p)
    left = localRemoveFile(f);
    outcome = 'indeterminate';
    cause = 'the test file could not be matched to the name MATLAB uses';
    return;
end
outcome = 'created';
left = localRemoveFile(f);
end

function cause = localJavaCause(javaME, dircause)
%the operating system's reason, without the Java stack trace
cause = '';
try
    if isprop(javaME, 'ExceptionObject') && ~isempty(javaME.ExceptionObject)
        cause = strtrim(char(javaME.ExceptionObject.getMessage()));
    end
catch
    cause = '';
end
if isempty(cause)
    tok = regexp(javaME.message, 'Exception:\s*([^\r\n]*)', 'tokens', 'once');
    if ~isempty(tok)
        cause = strtrim(tok{1});
    end
end
if isempty(cause)
    cause = dircause;
end
end

function javapath = localJavaPath(savepath, name)
%NAME in the save folder as an absolute path, for the Java file test only.
%java.io.File documents that a relative path is resolved against the user.dir
%system property (getAbsolutePath uses it), which need not be MATLAB's current
%folder. Where a relative path is actually created is up to the runtime: on
%one macOS runtime (2026-09-23), getAbsolutePath reported user.dir while
%createNewFile created the file in MATLAB's current folder. An absolute path
%removes the dependence on either behavior.
%dir reports the absolute folder MATLAB itself resolves SAVEPATH to. Only
%the fallback calls this, so the ordinary folder test never lists the folder.
%'' when the folder cannot be listed, or when its path holds '*', which dir
%would read as a pattern (and so report another folder).
javapath = '';
if any(savepath == '*')
    return;
end
d = dir(savepath);
if ~isempty(d)
    javapath = fullfile(d(1).folder, name);
end
end

function left = localRemove(probepath)
%Remove a test folder THIS function created: non-recursive, so anything that
%appeared inside it is never deleted, and verified afterwards.
left = {};
[~, msg] = rmdir(probepath);
if isfolder(probepath)
    left = {probepath};
    %the cause is unknown (something may have been put inside it, or a
    %permission changed), so the advice never tells anyone to destroy content
    %the check did not create
    warning('psyrat:outputNameProbeLeft', ...
        ['The output-name check created the test folder %s and could not ' ...
        'remove it (%s). Nothing else was touched; check its contents before ' ...
        'deleting it.'], probepath, strtrim(msg));
end
end

function left = localRemoveFile(f)
%Remove a test file THIS function created, only while it is still empty (any
%content can only have come from someone else), and verify afterwards. Nothing
%here can stop another process from replacing or reopening the file between
%the check and the removal.
left = {};
fpath = char(f.getAbsolutePath());
why = '';
if f.length() == 0
    if ~f.delete()
        why = 'removal failed';
    end
else
    why = 'it is no longer empty, so it was not removed';
end
if f.exists()
    if isempty(why)
        why = 'it is still there';
    end
    left = {fpath};
    warning('psyrat:outputNameProbeLeft', ...
        ['The output-name check created the test file %s and it remains (%s). ' ...
        'Nothing else was touched; check its contents before deleting it.'], fpath, why);
end
end

function h = localHex(n)
%n hexadecimal characters from tempname's UUID. tempname builds the name from
%matlab.lang.internal.uuid (R2025b tempname.m) and never calls the random
%number generator; a folder argument keeps it from consulting tempdir.
h = '';
while numel(h) < n
    [~, t] = fileparts(tempname('.'));
    t = strrep(t(3:end), '_', '');   %drop the 'tp' prefix and separators
    h = [h t]; %#ok<AGROW>
end
h = h(1:n);
end
