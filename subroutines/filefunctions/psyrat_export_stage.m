function varargout = psyrat_export_stage(action, varargin)
%Write a table export to a temporary file in the target's folder, and move it
%onto the target only after it was written and closed without a reported
%error (G109).
%
% [fid, stage] = psyrat_export_stage('open', target)
% psyrat_export_stage('close', fid, stage)
% psyrat_export_stage('commit', fid, stage, target)
% psyrat_export_stage('abort', fid, stage, target)
%
%'open' creates a file in TARGET's folder that no file or folder occupies
%(psyrat_partial_<k><ext>, a short name so a target at the folder's name-length
%limit can still be staged, in the same folder so the final move never crosses a
%volume) and opens it for writing. The open is in append mode and a file that
%turns out to hold bytes is closed and skipped, the collision guard of
%psyrat_fopen_unused, so a file someone else created is never truncated. If no
%such file can be opened it tells the user that nothing was saved and returns
%fid = -1 and stage = ''; TARGET is not touched.
%
%'close' closes the stage file and raises an error if ferror reports a write
%error for the stream at that point or the close itself fails. The file is
%closed either way. A write failure that neither reports is not detected here.
%A writer that closes the stage before another writer appends to it (the
%criterion table's writetable) uses this instead of a bare fclose; what that
%second writer does afterwards is checked only through the errors it raises.
%
%'commit' runs the same checked close if the stage is still open, refuses a
%TARGET that has become a folder, and then moves the stage file onto TARGET.
%That move is the only step that can change an existing TARGET. When the move
%fails, by returning status 0 or by raising an error, and when TARGET has
%become a folder, commit tries to remove the stage file, tells the user what
%happened and whether the temporary file was removed or is still there, and
%raises an error. It does not say TARGET is unchanged: a replacement the
%operating system abandoned part-way is not observable from here.
%
%'abort' is the cleanup a caller registers with onCleanup right after 'open'.
%After commit has moved the file or reported its own failure (a failed move, or
%a TARGET that became a folder) it does nothing. After any other failure,
%including a write or close failure raised by commit's checked close, it closes
%the stage file if it is still open, tries to remove it, and tells the user
%that nothing was saved, whether the temporary file was removed or is still
%there, and that no existing file was replaced; the error itself still reaches
%the Command Window.
%
%Used by the CSV Save callbacks that write with fprintf (Save Subject
%Coefficients, the View All difference table, and the criterion table's CSV
%branch). The ext of TARGET is kept on the stage file, so writetable can append
%to it in the same format.

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

%stage files a commit has settled (moved, or failed and reported), so the
%caller's onCleanup abort neither repeats commit's dialog nor touches a file
%that now carries the stage's name but is not this task's
persistent settled
if isempty(settled)
    settled = {};
end

switch action
    case 'open'
        target = varargin{1};
        [folder,~,ext] = fileparts(target);
        fid = -1;
        stage = '';
        %the first free name; bounded so a folder full of leftovers cannot spin
        %forever. A refused open is not retried under another name: the folder
        %is refusing new files, not this one.
        for k = 1:1000
            candidate = fullfile(folder, sprintf('psyrat_partial_%d%s', k, ext));
            %a target that itself carries a stage name is never its own stage
            %(compared without case: the default macOS volume ignores it)
            if isfile(candidate) || isfolder(candidate) || strcmpi(candidate, target)
                continue;
            end
            %append mode creates a missing file but never truncates one that
            %appeared since the check; a file that already holds bytes belongs
            %to someone else, so it is closed untouched and the next name tried
            trial = fopen(candidate, 'a');
            if trial < 0
                break;
            end
            fseek(trial, 0, 'eof');
            if ftell(trial) > 0
                fclose(trial);
                continue;
            end
            fid = trial;
            stage = candidate;
            break;
        end
        if fid < 0
            local_notsaved(target, {['A temporary file could not be created ', ...
                'in that folder, so nothing was written and no existing file ', ...
                'was replaced.']; 'How to fix: save to a folder you can write to.'});
        else
            settled(strcmp(settled, stage)) = [];
        end
        varargout = {fid, stage};

    case 'close'
        [fid, stage] = varargin{1:2};
        local_checkedclose(fid, stage);

    case 'commit'
        [fid, stage, target] = varargin{1:3};
        %a write or close failure raises here and is reported by abort
        local_checkedclose(fid, stage);
        [~,tname,tx] = fileparts(target);
        tfile = [tname tx];
        if isfolder(target)
            %the destination was checked before the stage opened; this catches
            %a folder created with the target's name since then
            settled{end+1} = stage;
            local_notsaved(target, [{sprintf(['"%s" is now a folder, so the ', ...
                'new file was not moved there.'], tfile)}; ...
                local_removestage(stage); ...
                {'How to fix: save again under a different name.'}]);
            error('psyrat:exportStage:folder', ...
                'The finished file was not moved onto %s, which is a folder.', target);
        end
        %a move that raises instead of returning status 0 is reported the same
        %way: whether it changed TARGET is not observable here either
        try
            [done, msg] = movefile(stage, target, 'f');
        catch ME
            done = false;
            msg = ME.message;
        end
        settled{end+1} = stage;
        if ~done
            local_notsaved(target, [{sprintf(['The new file was written, but it ', ...
                'could not replace %s (%s).'], tfile, msg)}; ...
                {sprintf(['If %s already existed, check it before relying on it: ', ...
                'whether the failed replacement changed it cannot be confirmed.'], ...
                tfile)}; ...
                local_removestage(stage)]);
            error('psyrat:exportStage:replace', ...
                'The finished file could not replace %s: %s', target, msg);
        end

    case 'abort'
        [fid, stage, target] = varargin{1:3};
        if isempty(stage)
            return;
        end
        hit = strcmp(settled, stage);
        if any(hit)
            settled(hit) = [];
            return;
        end
        local_forceclose(fid, stage);
        if isfile(stage)
            local_notsaved(target, [{['The new file could not be written ', ...
                'completely, so it did not replace any existing file.']}; ...
                local_removestage(stage); ...
                {'The cause is shown in the Command Window.'}]);
        end

    otherwise
        error('psyrat:exportStage:action', 'Unknown action ''%s''.', action);
end

end

function local_checkedclose(fid, stage)
%Close the stage file if FID is still that file (a caller that already closed it
%through 'close' leaves nothing to do here), then raise an error if the stream
%recorded a write error or the close failed. The close happens either way, so a
%failure never leaves the file open.
if fid < 0 || ~strcmp(fopen(fid), stage)
    return;
end
werr = ferror(fid);
status = fclose(fid);
if ~isempty(werr)
    error('psyrat:exportStage:write', 'Writing %s failed: %s', stage, werr);
end
if status ~= 0
    error('psyrat:exportStage:write', 'Closing %s failed.', stage);
end
end

function local_forceclose(fid, stage)
%Close the stage file if FID is still that file, without checking the result:
%the export has already failed and the file is about to be removed.
if fid >= 0 && strcmp(fopen(fid), stage)
    try
        fclose(fid);
    catch
        %the removal below reports what is left on disk
    end
end
end

function line = local_removestage(stage)
%Remove the stage file and describe the result as observed on disk afterwards:
%nothing when there was no stage file, otherwise whether it is gone.
line = {};
if ~isfile(stage)
    return;
end
try
    delete(stage);
catch
    %reported below from what is on disk
end
[~,sname,sext] = fileparts(stage);
if isfile(stage)
    line = {sprintf(['The temporary file %s could not be removed; delete it ', ...
        'yourself.'], [sname sext])};
else
    line = {sprintf('The temporary file %s was removed.', [sname sext])};
end
end

function local_notsaved(target, lines)
[folder,name,ext] = fileparts(target);
errordlg([{sprintf('%s was not saved.', [name ext])}; lines(:); ...
    {sprintf('Folder: %s', folder)}], 'File Not Saved', 'modal');
end
