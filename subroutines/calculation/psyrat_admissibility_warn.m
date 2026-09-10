function psyrat_admissibility_warn(adm)
%Warn once about inadmissible quantities recorded by psyrat_admissibility
%
%psyrat_admissibility_warn(adm)
%
%Inputs
% adm - diagnostic struct (or array of them) from psyrat_admissibility
%
%Emits the shared warning identifier psyrat:negerrvar with the text built by
%psyrat_admissibility_note, so the console says exactly what the viewer, the
%exports and the report say. Does nothing when there is nothing to report.
%
%Only a NEGATIVE quantity warrants a warning. An exactly-zero universe-score
%variance is a legitimate result (a difference score with no between-person
%variance) and is carried in the note without raising a warning, so ordinary
%runs stay quiet.
%
%Callers that invoke a coefficient kernel in a loop -- once per participant, or
%once per point on a conditional-reliability surface -- should suppress this
%identifier around the loop, accumulate the per-call records, and call this
%once with the aggregate. Otherwise a single table emits one block per row.

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

if ~psyrat_admissibility_hasneg(adm)
    return;
end

msg = psyrat_admissibility_note(adm);

if ~isempty(msg)
    warning('psyrat:negerrvar','%s',msg);
end

end
