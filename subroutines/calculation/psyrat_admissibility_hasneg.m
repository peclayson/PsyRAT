function tf = psyrat_admissibility_hasneg(adm)
%True when an admissibility record contains a negative quantity
%
%tf = psyrat_admissibility_hasneg(adm)
%
%Inputs
% adm - diagnostic struct (or array of them) from psyrat_admissibility
%
%Outputs
% tf - logical scalar; false for empty, malformed or clean records
%
%This is the single test for "is there something here the user must be told
%about". It is separate from psyrat_admissibility_note because several callers
%need to branch on the condition without building the text, and because a
%record that contains only exactly-zero entries is NOT a fault: it is a
%difference score with no between-person variance, which is a real result.
%
%Fails open, like the other readers: anything unexpected returns false rather
%than throwing, so a diagnostic can never be the reason a viewer or an export
%breaks.

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

tf = false;

if isempty(adm) || ~isstruct(adm)
    return;
end

for k = 1:numel(adm)
    if isfield(adm(k),'entries') && ~isempty(adm(k).entries)
        if any([adm(k).entries.nneg] > 0)
            tf = true;
            return;
        end
    end
end

end
