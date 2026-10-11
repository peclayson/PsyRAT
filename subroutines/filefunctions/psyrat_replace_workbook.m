function psyrat_replace_workbook(target, cells, datatable, namesrow)
%Save a GUI XLSX export as a fresh, checked workbook, and only then replace the
%chosen file (G110).
%
% psyrat_replace_workbook(target, cells)
% psyrat_replace_workbook(target, cells, datatable, namesrow)
%
%The GUI XLSX Save callbacks used to write their header and table with
%writecell/writetable straight into the file the save dialog returned. On an
%existing workbook those writers change only the cells they cover, so the old
%workbook's trailing rows, wider columns and extra sheets survived a Replace.
%This writes the same cells, with the same calls, into a NEW workbook beside
%TARGET, reads that workbook back and compares every cell with what the caller
%asked for, and only then moves it onto TARGET. TARGET is never deleted first.
%
%Input
% target    - full path of the .xlsx file the user chose in the save dialog.
% cells     - the header block, written from A1 with writecell exactly as the
%             callbacks always wrote it (a column of header lines, or the whole
%             Save IDs sheet as one cell array).
% datatable - optional table, written as writetable(datatable, ..., 'Range',
%             ['A' namesrow]) always wrote it: variable names on NAMESROW and
%             the rows below. A table with NO rows gets its variable names on
%             NAMESROW by writecell instead, because writetable writes an empty
%             table's names at A1, over the header.
% namesrow  - the row of the variable names; each callback passes the row it
%             always used.
%
%Checked BEFORE the chosen file is touched. Any failure stops the save, shows a
%"File Not Saved" dialog and raises an error, so the caller reports no success.
% * The payload must be one this save can write and check exactly: finite
%   numbers, logicals, empty cells, and text of at most 32,767 characters
%   (longer text would be cut). NaN is stored as an EMPTY CELL - that is how the
%   XLSX writer serializes it - and the check accepts an empty cell for NaN and
%   nothing else. Inf and -Inf are refused: the writer would store 65535. Any
%   other kind of value is refused as unsupported rather than written.
% * Some text is refused before anything is written because this save cannot
%   confirm it, although the format may hold it: a carriage return (the writer
%   stores it in a form XML readers return as a line feed), any other control
%   character except tab and line feed (the writer escapes it, or drops
%   character 0), and "_xHHHH_" with a lower-case x and four hexadecimal
%   digits (the writer escapes it, and xlsread, the reader of the check,
%   returns the escaped form). These are limits of the writer and reader used
%   here, not of XLSX; CSV keeps such text. "_XHHHH_" with an upper-case X is
%   not an escape and is saved as it is.
% * The new workbook must have exactly one sheet, and reading it back gives
%   exactly the payload in its places. It is read from A1 over the payload's
%   extent (an explicit range keeps leading empty rows and columns, so a payload
%   written anywhere but A1 fails) and over the sheet's used range (so nothing
%   may lie outside the payload). The comparison is exact: the same text
%   (including leading and trailing spaces, which is why it is read with xlsread
%   in 'basic' mode; readcell trims them), the same finite values exactly (0 and
%   -0 compare equal), the same logicals, and empty cells where the payload is
%   empty. Text the check's reader does not return exactly, such as characters
%   beyond the Basic Multilingual Plane, fails here, after writing, and the save
%   stops.
%
%Temporary files: psyrat_stage_<k>.xlsx (the new workbook) and
%psyrat_backup_<k>.xlsx (the copy of the existing file) in TARGET's folder,
%where k is the first number whose name no file or folder there uses, compared
%without regard to case, and which is not TARGET's own name. The names do not
%grow with TARGET's name. Only files this call created are removed. Whether a
%removal worked is checked afterwards, and a temporary file still there is
%reported by its full path, never as removed.
%
%Replacing. When TARGET exists, its bytes are first written to a backup file
%and read back to check the copy. The new workbook is then moved onto TARGET
%and TARGET is read back. If the move fails, or TARGET does not hold the
%checked workbook afterwards, TARGET is compared with the original bytes and,
%if it differs or is gone, copied back from the backup and checked again. The
%backup is removed only after TARGET has been shown to hold the new workbook or
%the original. If the original cannot be put back, the backup is kept and the
%message says where it is. If the save succeeds but the backup cannot be
%removed, the save stands, and a dialog and a warning name the backup. Every
%file is read by its size and checked (count, ferror, close), so a partial read
%is never taken for the file. This is not an atomic replacement, and nothing
%here can guard against a failure of the computer or the file system itself.

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

if nargin ~= 2 && nargin ~= 4
    error('psyrat:replaceWorkbook:args', ...
        'Call psyrat_replace_workbook(target, cells) or (target, cells, datatable, namesrow).');
end
hastable = nargin == 4;
target = char(target);
[folder, name, ext] = fileparts(target);
shown = [name ext];
nothing = 'Nothing was written.';

%1. Refusals that need no writing.
if ~strcmpi(ext, '.xlsx')
    local_fail(shown, 'psyrat:replaceWorkbook:refused', 'The name does not end in .xlsx.', nothing);
end
if isfolder(target)
    local_fail(shown, 'psyrat:replaceWorkbook:refused', 'A folder has that name.', nothing);
end
if isempty(folder) || ~isfolder(folder)
    local_fail(shown, 'psyrat:replaceWorkbook:refused', 'Its folder does not exist.', nothing);
end
if hastable
    [expected, problem] = local_expected(cells, datatable, namesrow);
else
    [expected, problem] = local_expected(cells, [], []);
end
if ~isempty(problem)
    local_fail(shown, 'psyrat:replaceWorkbook:refused', problem, nothing);
end

%2. Write the new workbook beside TARGET and check it cell by cell.
stage = local_free_name(folder, 'psyrat_stage_', shown);
if isempty(stage)
    local_fail(shown, 'psyrat:replaceWorkbook:refused', ...
        'No free name for a temporary workbook was found in its folder.', nothing);
end
try
    writecell(cells, stage);
    if hastable
        if height(datatable) > 0
            writetable(datatable, stage, 'Range', sprintf('A%d', namesrow));
        else
            writecell(datatable.Properties.VariableNames, stage, 'Range', sprintf('A%d', namesrow));
        end
    end
    local_check(stage, expected);
    newbytes = local_bytes(stage);
catch ME
    [sleft, sgone] = local_remove(stage);
    local_fail(shown, 'psyrat:replaceWorkbook:check', ...
        ['The new workbook did not pass the check made before replacing the file: ' ME.message], ...
        local_sentences(local_stage_text(sleft, sgone), sprintf('%s was not replaced.', shown)));
end

%3. Keep a checked copy of the existing file.
hadtarget = isfile(target);
original = [];
backup = '';
if hadtarget
    try
        original = local_bytes(target);
        backup = local_free_name(folder, 'psyrat_backup_', shown);
        if isempty(backup)
            error('psyrat:replaceWorkbook:noname', 'no free name for the copy was found');
        end
        local_write_bytes(backup, original);
        if ~local_same(backup, original)
            error('psyrat:replaceWorkbook:copy', 'the copy did not read back byte for byte');
        end
    catch ME
        bleft = local_remove(backup);
        [sleft, sgone] = local_remove(stage);
        local_fail(shown, 'psyrat:replaceWorkbook:backup', ...
            ['A checked copy of the existing file could not be made: ' ME.message], ...
            local_sentences(local_stage_text(sleft, sgone), local_backup_text(bleft), ...
            sprintf('%s was not replaced.', shown)));
    end
end

%4. Move the new workbook onto TARGET and read TARGET back.
try
    [moved, why] = movefile(stage, target, 'f');
catch ME
    moved = false;
    why = ME.message;
end
if moved && local_same(target, newbytes)
    bleft = local_remove(backup);
    if ~isempty(bleft)
        local_report_leftover(shown, bleft);
    end
    return;
end
if moved
    why = 'the file in place afterwards is not the checked workbook';
end
[sleft, sgone] = local_remove(stage);
stagetext = local_stage_text(sleft, sgone);

%5. The replacement failed: say exactly what is known about TARGET.
failed = sprintf('The new workbook could not be moved into place (%s).', why);
if ~hadtarget
    placed = sprintf('No file named %s is in that folder.', shown);
    if isfile(target)
        placed = sprintf(['A file named %s now exists in that folder; it did not pass ', ...
            'the check, so do not rely on it.'], shown);
    end
    local_fail(shown, 'psyrat:replaceWorkbook:replace', failed, local_sentences(stagetext, placed));
end
if local_same(target, original)
    bleft = local_remove(backup);
    local_fail(shown, 'psyrat:replaceWorkbook:replace', failed, local_sentences( ...
        sprintf('%s was read back afterwards and is byte for byte the same as before.', shown), ...
        stagetext, local_backup_text(bleft)));
end
try
    [restored, why2] = copyfile(backup, target, 'f');
catch ME
    restored = false;
    why2 = ME.message;
end
if restored && local_same(target, original)
    bleft = local_remove(backup);
    local_fail(shown, 'psyrat:replaceWorkbook:replace', ...
        sprintf('Moving the new workbook into place failed (%s) and changed the existing file.', why), ...
        local_sentences(sprintf('The original %s was put back from its copy and checked byte for byte.', ...
        shown), stagetext, local_backup_text(bleft)));
end
if restored || isempty(why2)
    why2 = 'the file in place afterwards is not the original';
end
%the backup is deliberately not removed here: it is the only checked copy
local_fail(shown, 'psyrat:replaceWorkbook:restore', ...
    sprintf(['Moving the new workbook into place failed (%s) and changed or removed the ', ...
    'existing file, and putting the original back failed (%s).'], why, why2), ...
    local_sentences(sprintf(['A byte-for-byte copy of the original %s is kept at %s. Keep that ', ...
    'file: it is the only checked copy of the original.'], shown, backup), stagetext));

end

% -- payload ---------------------------------------------------------------------

function [E, problem] = local_expected(cells, T, namesrow)
%Every cell the workbook must hold, from the payload alone: [] for an empty
%cell, char for text, double for a number, logical for a logical.
E = {};
problem = '';
if ~iscell(cells) || ~ismatrix(cells) || isempty(cells)
    problem = 'The header is not a non-empty two-dimensional cell array.';
    return;
end
[R, C] = size(cells);
head = cell(R, C);
for r = 1:R
    for c = 1:C
        [head{r, c}, what] = local_value(cells{r, c});
        if ~isempty(what)
            problem = sprintf('Cell %s holds %s.', local_address(r, c), what);
            return;
        end
    end
end
if isempty(T) && ~istable(T)
    E = head;
else
    if ~istable(T) || ~isnumeric(namesrow) || ~isscalar(namesrow) || namesrow < 1 || ...
            namesrow ~= round(namesrow)
        problem = 'The table or the row of its variable names is not valid.';
        return;
    end
    names = T.Properties.VariableNames;
    [h, w] = size(T);
    for j = 1:w
        v = T.(names{j});
        if ~ischar(v) && size(v, 2) ~= 1
            problem = sprintf('The table variable %s has more than one column.', names{j});
            return;
        end
    end
    for r = namesrow:R
        for c = 1:C
            if ~isempty(head{r, c})
                problem = sprintf('Header cell %s lies where the table is written.', local_address(r, c));
                return;
            end
        end
    end
    E = repmat({[]}, max(R, namesrow + h), max(C, w));
    E(1:R, 1:C) = head;
    for j = 1:w
        [E{namesrow, j}, what] = local_value(names{j});
        if ~isempty(what) || ~ischar(E{namesrow, j})
            if isempty(what)
                what = 'no text';
            end
            problem = sprintf('The table variable name %s holds %s.', names{j}, what);
            return;
        end
    end
    for i = 1:h
        for j = 1:w
            [E{namesrow + i, j}, what] = local_value(T{i, j});
            if ~isempty(what)
                problem = sprintf('Cell %s holds %s.', local_address(namesrow + i, j), what);
                return;
            end
        end
    end
end
if isempty(E{1, 1})
    problem = 'Cell A1 is empty, so the new workbook cannot be checked against it.';
end
end

function [x, what] = local_value(v)
%One payload value in the form the check compares, or WHAT describes why this
%save cannot write and check it exactly.
x = [];
what = '';
if iscell(v) && isscalar(v)
    v = v{1};
end
if iscell(v) || isstruct(v) || (isobject(v) && ~isstring(v) && ~iscategorical(v))
    what = sprintf('a %s value, which this export does not support', class(v));
elseif isstring(v) || iscategorical(v)
    if ~isscalar(v) || ismissing(v)
        what = 'a missing or multi-element text value';
    else
        [x, what] = local_value(char(string(v)));
    end
elseif ischar(v)
    if isempty(v)
        x = [];
    elseif size(v, 1) ~= 1
        what = 'text with more than one row';
    elseif numel(v) > 32767
        what = sprintf('text of %d characters; XLSX keeps at most 32,767 and would cut it', numel(v));
    elseif any(v < 32 & v ~= 9 & v ~= 10 & v ~= 13)
        what = ['a control character, which this save cannot write and read back exactly ', ...
            '(the writer escapes it, or drops character 0)'];
    elseif any(v == 13)
        what = ['a carriage return, which the XLSX writer used here stores in a form that ', ...
            'XML readers return as a line feed; save as CSV to keep it'];
    elseif ~isempty(regexp(v, '_x[0-9A-Fa-f]{4}_', 'once'))
        %lower-case x only: that is the escape form; "_X0041_" is plain text
        what = sprintf(['the text %s, which the XLSX writer used here stores escaped and ', ...
            'the check reads back in escaped form, so this save cannot confirm it; save as ', ...
            'CSV to keep it'], regexp(v, '_x[0-9A-Fa-f]{4}_', 'match', 'once'));
    else
        x = v;
    end
elseif islogical(v)
    if isscalar(v)
        x = v;
    else
        what = 'several logical values in one cell';
    end
elseif isnumeric(v)
    if isempty(v)
        x = [];
    elseif ~isscalar(v) || ~isreal(v)
        what = 'several or complex numbers in one cell';
    elseif isinteger(v) && abs(double(v)) > flintmax
        what = 'an integer too large to store exactly';
    elseif isnan(v)
        x = [];   %the XLSX writer stores NaN as an empty cell
    elseif isinf(v)
        what = 'an infinite value, which the XLSX writer used here would save as 65535';
    else
        x = double(v);
    end
else
    what = sprintf('a %s value, which this export does not support', class(v));
end
end

% -- check -----------------------------------------------------------------------

function local_check(stage, E)
%Raise an error unless STAGE has one sheet holding exactly E from A1. It is
%read twice with xlsread 'basic' (which keeps leading and trailing spaces):
%over E's own extent from A1, an explicit range that keeps leading empty rows
%and columns (a used-range read drops them, so a payload written one row or
%column away would read back alike), and over the sheet's used range, which
%starts at A1 once A1 holds E's first cell (never empty), so a cell outside E
%is seen too.
sheets = sheetnames(stage);
if numel(sheets) ~= 1
    error('psyrat:replaceWorkbook:sheets', 'it has %d sheets instead of 1', numel(sheets));
end
[R, C] = size(E);
[~, ~, anchored] = xlsread(stage, '', ['A1:' local_address(R, C)], 'basic'); %#ok<XLSRD>
local_compare(anchored, E);
[~, ~, used] = xlsread(stage, '', '', 'basic'); %#ok<XLSRD>
local_compare(used, E);
end

function local_compare(raw, E)
%Raise an error at the first cell where RAW, read from A1, differs from E;
%every cell beyond either one must be empty.
[R, C] = size(E);
[R2, C2] = size(raw);
for r = 1:max(R, R2)
    for c = 1:max(C, C2)
        want = [];
        if r <= R && c <= C
            want = E{r, c};
        end
        got = NaN;
        if r <= R2 && c <= C2
            got = raw{r, c};
        end
        if ~local_match(want, got)
            error('psyrat:replaceWorkbook:mismatch', 'cell %s reads %s where the export holds %s', ...
                local_address(r, c), local_show(got), local_show(want));
        end
    end
end
end

function tf = local_match(want, got)
%Exact agreement: text character for character, finite values exactly (0 and
%-0 agree, as they do under ==), logicals by value, empty with empty.
blank = (isnumeric(got) && (isempty(got) || (isscalar(got) && isnan(got)))) || ...
    (ischar(got) && isempty(got));
if isnumeric(want) && isempty(want)
    tf = blank;
elseif ischar(want)
    tf = ischar(got) && strcmp(got, want);
elseif islogical(want)
    tf = islogical(got) && isscalar(got) && got == want;
else
    tf = isnumeric(got) && ~islogical(got) && isscalar(got) && isreal(got) && ...
        double(got) == want;
end
end

function s = local_show(v)
if (isnumeric(v) && (isempty(v) || (isscalar(v) && isnan(v)))) || (ischar(v) && isempty(v))
    s = 'an empty cell';
elseif ischar(v)
    s = ['"' v(1:min(end, 60)) '"'];
elseif islogical(v) && isscalar(v)
    s = mat2str(v);
elseif isnumeric(v) && isscalar(v)
    s = sprintf('%.17g', v);
else
    s = sprintf('a %s value', class(v));
end
end

% -- files -----------------------------------------------------------------------

function p = local_free_name(folder, prefix, shown)
%The first <prefix><k>.xlsx that no file or folder in FOLDER uses, compared
%without regard to case, and that is not the chosen file's own name.
entries = dir(folder);
taken = lower({entries.name});
p = '';
for k = 1:999
    candidate = sprintf('%s%d.xlsx', prefix, k);
    if strcmpi(candidate, shown) || any(strcmp(taken, lower(candidate)))
        continue;
    end
    full = fullfile(folder, candidate);
    if exist(full, 'file') || isfolder(full)
        continue;
    end
    p = full;
    return;
end
end

function b = local_bytes(file)
%Every byte of FILE, read by its size and checked, or an error: the size is
%taken from the open file, exactly that many bytes are read (reading to Inf
%ends at end of file, which ferror then reports), and the count, ferror, a
%one-byte look past the end and fclose must all show a whole, clean read. A
%partial read is never returned. The file is closed on every path.
fid = fopen(file, 'r');
if fid < 0
    error('psyrat:replaceWorkbook:read', 'could not open %s', file);
end
closer = onCleanup(@() local_close(fid)); %#ok<NASGU>
n = -1;
if fseek(fid, 0, 'eof') == 0
    n = ftell(fid);
end
if n < 0 || fseek(fid, 0, 'bof') ~= 0
    error('psyrat:replaceWorkbook:read', 'could not find the size of %s', file);
end
ferror(fid, 'clear');   %so the ferror below reports on the read alone
[b, count] = fread(fid, n, '*uint8');
msg = ferror(fid);
more = fread(fid, 1, '*uint8');
status = fclose(fid);
if count ~= n || numel(b) ~= n || ~isempty(msg) || ~isempty(more) || status ~= 0
    why = sprintf('read %d of its %d bytes', count, n);
    if ~isempty(msg)
        why = sprintf('%s (%s)', why, msg);
    elseif ~isempty(more)
        why = 'it grew while it was read';
    elseif status ~= 0
        why = 'it could not be closed after reading';
    end
    error('psyrat:replaceWorkbook:read', 'could not read all of %s: %s', file, why);
end
end

function local_write_bytes(file, b)
%Write B to FILE and check the count and the close; the file is closed on
%every path, also when fwrite raises.
fid = fopen(file, 'w');
if fid < 0
    error('psyrat:replaceWorkbook:write', 'could not create %s', file);
end
closer = onCleanup(@() local_close(fid)); %#ok<NASGU>
count = fwrite(fid, b, 'uint8');
status = fclose(fid);
if count ~= numel(b) || status ~= 0
    error('psyrat:replaceWorkbook:write', 'could not write all of %s', file);
end
end

function local_close(fid)
%The guard behind each fopen here: close FID if it is still open, which it is
%only when a read or write raised before its own checked fclose. A failure here
%is not raised, so it cannot replace the error already on its way out.
try
    if ~isempty(fopen(fid))
        fclose(fid);
    end
catch
end
end

function tf = local_same(file, b)
tf = false;
try
    tf = isfile(file) && isequal(local_bytes(file), b);
catch
end
end

function [left, gone] = local_remove(file)
%Remove a file this call created and report what is observed afterwards: LEFT
%is FILE when it is still there ('' otherwise), GONE is true when it existed
%and is now gone. delete can refuse with no more than a warning, so the result
%is checked, never assumed.
left = '';
gone = false;
if isempty(file) || ~isfile(file)
    return;
end
try
    delete(file);
catch
end
if isfile(file)
    left = file;
else
    gone = true;
end
end

function s = local_stage_text(left, gone)
%What became of the new workbook, as observed by local_remove.
s = '';
if ~isempty(left)
    s = sprintf('The new workbook could not be removed and is still at %s; it can be deleted.', left);
elseif gone
    s = 'The new workbook was removed.';
end
end

function s = local_backup_text(left)
%A copy of the existing file that could not be removed after a failed save
%('' when it is gone or was never made). TARGET is intact in every branch that
%reports this, so the copy is not needed.
s = '';
if ~isempty(left)
    s = sprintf(['The copy made of the existing file could not be removed and is still at %s; ', ...
        'it can be deleted.'], left);
end
end

function s = local_sentences(varargin)
%The non-empty sentences, in order, joined by spaces.
parts = varargin(~cellfun(@isempty, varargin));
s = strjoin(parts, ' ');
end

function local_report_leftover(shown, backup)
%The save succeeded but its copy of the previous file could not be removed: the
%save stands, and the copy is named, by its full path, in a dialog and in a
%warning. Neither can turn the successful save into an error: a caller who has
%made this warning an error gets the same text on the console instead.
lines = {sprintf('%s was saved.', shown); ...
    sprintf(['The copy of the previous %s that this save made could not be removed and is ', ...
    'still at %s. It holds the previous contents; delete it when it is no longer needed.'], ...
    shown, backup)};
try
    warndlg(lines, 'File Saved', 'modal');
catch
end
msg = strjoin(lines', ' ');
try
    warning('psyrat:replaceWorkbook:cleanup', '%s', msg);
catch
    fprintf(2, 'Warning: %s\n', msg);
end
end

function a = local_address(r, c)
letters = '';
while c > 0
    m = mod(c - 1, 26);
    letters = [char(double('A') + m) letters]; %#ok<AGROW>
    c = floor((c - 1) / 26);
end
a = sprintf('%s%d', letters, r);
end

function local_fail(shown, id, reason, outcome)
%One dialog and one error with the same words, so the caller shows no success.
lines = {sprintf('%s was not saved.', shown); reason; outcome};
errordlg(lines, 'File Not Saved', 'modal');
error(id, '%s', strjoin(lines', ' '));
end
