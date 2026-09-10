function note = psyrat_dod_admissibility_note(dodiffs)
%Gather the admissibility records from difference-of-differences results
%
%note = psyrat_dod_admissibility_note(dodiffs)
%
%Inputs
% dodiffs - one struct returned by psyrat_dodiffrel, or an array or cell array
%  of them (one per evaluation, e.g. one per group)
%
%Outputs
% note - char row explaining what was inadmissible, or '' when everything is
%  fine
%
%Why this exists
% psyrat_dodiffrel attaches an admissibility record to every coefficient it
% returns, but the difference-of-differences route builds no relsummary, so
% psyrat_admissibility_collect -- which walks relsummary.group(g).diffscore --
% cannot reach it. The record was written and read by nothing: on this route a
% negative variance produced psyrat_admissibility_warn's console line and
% nothing else, no report block and no dialog, which is exactly the surfacing a
% GUI user never sees. This is the DoD counterpart of that collector (RC-43).
%
% Only the notable case produces text. A clean run returns '', so callers can
% append or test the result unconditionally and a run with nothing wrong looks
% exactly as it did before this diagnostic existed -- the same convention
% psyrat_admissibility_collect and psyrat_provenance_lines already use.
%
%PASS ONE RECORD PER EVALUATION, NOT BOTH HALVES OF A dep/gen PAIR
% Every caller computes dependability and generalizability together from the
% same draws. Every recorded quantity -- u, bt_scaled, er_scaled, bt_unscaled --
% is built before psyrat_dodiffrel branches on est, so the two records are
% bit-identical. psyrat_admissibility_note rolls up by unique({label}) and SUMS
% nneg and ndraw, so passing both would not check anything twice; it would
% double every count inside one line and report that twice as many draws were
% inspected as actually were. Hand in either one and drop the other.
%
% Records from DIFFERENT evaluations are the opposite case and SHOULD all be
% passed: one group's counts genuinely add to another's, which is how
% psyrat_admissibility_collect aggregates across groups.

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
if nargin < 1 || isempty(dodiffs)
    return;
end

%accept a single struct, a struct array, or a cell array of structs, so callers
%can hand over one evaluation or a whole per-group sweep without reshaping
if ~iscell(dodiffs)
    dodiffs = num2cell(dodiffs);
end

adm_all = [];

for k = 1:numel(dodiffs)
    d = dodiffs{k};

    %psyrat_dod_observedt substitutes a NaN summary when a trial count is zero,
    %and that stand-in carries no record -- skip rather than fail
    if ~isstruct(d) || ~isfield(d,'admissibility') || isempty(d.admissibility)
        continue;
    end

    adm_all = [adm_all, d.admissibility]; %#ok<AGROW>
end

note = psyrat_admissibility_note(adm_all);

end
