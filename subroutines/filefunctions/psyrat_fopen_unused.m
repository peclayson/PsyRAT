function [fid, usedpath] = psyrat_fopen_unused(folder,stem,ext)
%Open a NEW file for writing under the first unused name in a folder.
%
%   [fid, usedpath] = psyrat_fopen_unused(folder,stem,ext)
%
%   Tries <stem><ext>, then <stem>_2<ext>, <stem>_3<ext>, ... in folder and
%   opens the first name that no file occupies. An existing file is never
%   truncated: the candidate is skipped when it exists, and the open itself
%   is in append mode, so even a file created between the check and the
%   open keeps its bytes (a non-empty one is then treated as occupied and
%   the next suffix is tried; an empty one receives the content). Returns
%   fid = -1 and usedpath = '' when a free name cannot be opened, which is
%   how a refusing folder shows itself (a FOLDER that carries a candidate
%   name has the same effect: isfile is false and the open fails), so the
%   caller can move to another destination.
%
%   The free-name convention is unchanged: a name that is free is used as
%   is; only collisions get a suffix. Introduced for the per-measurement
%   error reports (G94 residual, 2026-09-20): the report's stem is derived
%   from the output name the user typed in the input file's folder, so the
%   loaded input, or an earlier report, can already carry the report's name.

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

fid = -1;
usedpath = '';

%bound the search so a folder full of reports cannot spin forever; the
%bound is far beyond any real number of collisions
maxtries = 10000;

for n = 1:maxtries
    if n == 1
        candidate = fullfile(folder,[stem ext]);
    else
        candidate = fullfile(folder,sprintf('%s_%d%s',stem,n,ext));
    end

    %an existing file (a loaded input, an earlier report) is occupied
    if isfile(candidate)
        continue;
    end

    %append mode ('a') creates a missing file but never truncates an
    %existing one, unlike 'w'; that is what makes the check-then-open safe
    trial = fopen(candidate,'a');
    if trial == -1
        %a free name that cannot be opened means the folder refuses the
        %write; the caller decides on another destination
        return;
    end

    %race guard: if the name was taken between isfile and fopen, the file
    %already holds bytes; leave them alone and try the next suffix
    fseek(trial,0,'eof');
    if ftell(trial) > 0
        fclose(trial);
        continue;
    end

    fid = trial;
    usedpath = candidate;
    return;
end

end
