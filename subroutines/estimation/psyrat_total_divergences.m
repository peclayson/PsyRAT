function ndiv = psyrat_total_divergences(convdata)
%Total the HMC divergent-transition counts that psyrat_storeconv stamped into
%the header cell (column 4) of each fit's convergence table.
%
%ndiv = psyrat_total_divergences(convdata)
%
%Input:
% convdata - the REL.out.conv.data field. It takes the same two shapes
%  psyrat_checkconv handles: a 1xN cell whose elements are per-fit convergence
%  tables (one fit per measurement/occasion/group), or a single convergence
%  table (an Mx3-or-4 cell) when there is only one fit.
%
%Output:
% ndiv - the summed divergent-transition count across all fits, or NaN when no
%  fit recorded a count. NaN means "not assessed": native engines (MATLAB's
%  hmcSampler exposes no divergence diagnostic; fitlme is not MCMC), legacy
%  results files predating this field, and older Stan versions all yield NaN,
%  and the count is reported only when it is a real, non-negative number.
%
%Divergences are report-only. This total is surfaced to the user (console plus
%the results viewer via psyrat_convergence_notice) but never feeds the
%converged flag set by psyrat_checkconv.

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

ndiv = NaN;

if ~iscell(convdata) || isempty(convdata)
    return;
end

%mirror psyrat_checkconv's shape test: a row count of 1 means convdata is a
%1xN container of per-fit tables; otherwise convdata is a single table.
[a,nloops] = size(convdata);
if a == 1
    counts = nan(1,nloops);
    for i = 1:nloops
        counts(i) = local_cell_count(convdata{i});
    end
else
    counts = local_cell_count(convdata);
end

%report the sum only if at least one fit recorded a real count; otherwise leave
%NaN so downstream surfacing stays silent.
real = counts(~isnan(counts));
if ~isempty(real)
    ndiv = sum(real);
end

end

function c = local_cell_count(tbl)
%Read the divergence count from a single convergence table's header cell,
%returning NaN unless it is present and a real numeric scalar.
c = NaN;
if iscell(tbl) && size(tbl,2) >= 4
    v = tbl{1,4};
    if isnumeric(v) && isscalar(v)
        c = v;
    end
end
end
