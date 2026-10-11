function [target, ok] = psyrat_export_destination(savepath, savename, ext)
%Resolve the file a CSV-only export will actually write, and ask before an
%existing file the save dialog never named is replaced (G109).
%
% [target, ok] = psyrat_export_destination(savepath, savename, ext)
%
%The CSV-only table exports append EXT to a name that does not already end in
%it, so a typed or forced .xlsx name cannot receive CSV bytes. That keeps the
%file type honest, but the appended name is one the save dialog never showed,
%and the dialog's own replace confirmation covered only the name it returned.
%When the appended name is already a file, this asks about that exact file
%and returns ok = false unless the user chooses Replace. A returned name that
%already ends in EXT is used as is, with no second question: the dialog
%confirmed it. Either name is refused, with a dialog, when it is an existing
%folder. Nothing is opened or written here; callers write through
%psyrat_export_stage, whose commit checks for a folder again.
%
%Input
% savepath, savename - the folder and name uiputfile returned, after the
%   caller's Cancel guard.
% ext - the extension the export requires, e.g. '.csv'.
%
%Output
% target - full path of the file the export should write.
% ok     - false when the user declined, or when TARGET is a folder.

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

ok = true;
[~,~,given] = fileparts(savename);
exact = strcmpi(given, ext);
if exact
    target = fullfile(savepath, savename);
else
    target = fullfile(savepath, [savename ext]);
end

%a folder is refused under either name, before anything is opened: the move
%that finishes the export would otherwise put the file inside it
[~,fname,fext] = fileparts(target);
if isfolder(target)
    ok = false;
    errordlg({sprintf('"%s" is a folder, so nothing was saved.', [fname fext]); ...
        'How to fix: save again under a different name.'}, ...
        'File Not Saved', 'modal');
elseif exact
    %the save dialog confirmed this exact name
    return;
elseif isfile(target)
    answer = questdlg(sprintf(['This table is saved as %s, so it would be written to ', ...
        '"%s", a file the save dialog did not ask about. That file already exists ', ...
        'in %s. Replace it?'], ext, [fname fext], savepath), ...
        'Replace Existing File?', 'Replace', 'Cancel', 'Cancel');
    ok = strcmp(answer, 'Replace');
end

end
