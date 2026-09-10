function note = psyrat_admissibility_collect(relsummary)
%Gather the admissibility records scattered across a summary into one note
%
%note = psyrat_admissibility_collect(relsummary)
%
%Inputs
% relsummary - the relsummary struct built by psyrat_relsummary (or by
%  psyrat_dynrel_summary / psyrat_splits_summary)
%
%Outputs
% note - char row explaining what was inadmissible anywhere in the summary, or
%  '' when everything is fine
%
%Why this exists
% The coefficient kernels attach an admissibility record to each difference
% score they return, but a summary holds one per group and, on the subject-level
% paths, one per participant. A user wants to be told once, with the totals,
% rather than once per cell. This walks the summary, collects every record it
% can find, and hands the lot to psyrat_admissibility_note.
%
% Only the notable case produces text. A clean summary returns '', so callers
% can append the result unconditionally and a run with nothing wrong looks
% exactly as it did before this diagnostic existed -- the same convention
% psyrat_provenance_lines uses for its caveat lines.
%
% The note is NOT persisted anywhere, and does not need to be: a coefficient is
% not stored in the .psyrat file either. The summary is rebuilt every time
% results are viewed or reported, so the diagnostic is rebuilt with it and
% always describes the numbers actually on screen.

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
%

note = '';

%fail open: a diagnostic must never be the reason a viewer or export breaks
if isempty(relsummary) || ~isstruct(relsummary) || ~isfield(relsummary,'group')
    return;
end

adm_all = [];

for gloc = 1:numel(relsummary.group)
    grp = relsummary.group(gloc);

    if ~isfield(grp,'diffscore') || isempty(grp.diffscore)
        continue;
    end

    adm_all = local_gather(adm_all, grp.diffscore);

    %the cutoff row is a second, independently evaluated difference score and
    %carries its own record
    if isfield(grp.diffscore,'trlcutoff') && isstruct(grp.diffscore.trlcutoff)
        adm_all = local_gather(adm_all, grp.diffscore.trlcutoff);
    end
end

note = psyrat_admissibility_note(adm_all);

end

function adm_all = local_gather(adm_all, s)
%Append one struct's admissibility record, if it has one.
if isstruct(s) && isfield(s,'admissibility') && ~isempty(s.admissibility)
    adm_all = [adm_all, s.admissibility];
end
end
