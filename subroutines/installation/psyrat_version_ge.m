function tf = psyrat_version_ge(a,b)
%Compare two 'major.minor.patch' version strings: true if a >= b.
%
%tf = psyrat_version_ge(a,b)
%
%Inputs:
% a, b - version strings of the form 'major.minor.patch' (e.g. '2.26.0').
%
%Output:
% tf - logical; true if a is greater than or equal to b. FAIL-OPEN: if either
%  string cannot be parsed into three integer parts, returns true (do not block
%  on an unparseable version).
%
%This is the SINGLE source of truth for CmdStan version gating. It is shared by
%psyrat_estimation_engines_available (which sets the .cmdstan availability flag)
%and psyrat_require_min_cmdstan (the run-time preflight) so the two cannot drift:
%the availability dropdown and the preflight always agree on whether a detected
%CmdStan meets the minimum. The fail-open behavior matches both callers' prior
%stance (an unreadable version never blocks an otherwise valid run).
%
%See also: psyrat_estimation_engines_available, psyrat_require_min_cmdstan,
% psyrat_detect_cmdstan_version

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
%

pa = sscanf(a,'%d.%d.%d')';
pb = sscanf(b,'%d.%d.%d')';
if numel(pa) < 3 || numel(pb) < 3
    tf = true;   % unparseable -> fail-open, do not block
    return;
end

%lexicographic major.minor.patch comparison
for k = 1:3
    if pa(k) < pb(k)
        tf = false;
        return;
    elseif pa(k) > pb(k)
        tf = true;
        return;
    end
end
tf = true;   % all three parts equal -> a >= b
end
