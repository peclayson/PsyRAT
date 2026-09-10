function [maxRhat, nParams] = psyrat_conv_maxrhat(convdata)
%Summarize the potential scale reduction factors (r_hat) stored in a
%REL.out.conv.data field: return the worst (largest) r_hat and the number of
%real sampled parameters it was computed over.
%
%[maxRhat, nParams] = psyrat_conv_maxrhat(convdata)
%
%Input:
% convdata - the REL.out.conv.data field. It takes the same two shapes
%  psyrat_checkconv and psyrat_total_divergences handle: a 1xN cell whose
%  elements are per-fit convergence tables (one fit per measurement / occasion /
%  group), or a single convergence table (an Mx3-or-4 cell {name,n_eff,r_hat})
%  when there is only one fit. Row 1 of each table is the header; data rows are
%  2:end.
%
%Outputs:
% maxRhat - the largest r_hat across every data row of every fit table. A failed
%  extraction stores the {'convergence_unavailable',0,Inf} sentinel (see
%  psyrat_storeconv), so an unextractable summary yields Inf. An empty /
%  header-only table also yields Inf. maxRhat is therefore never a small,
%  reassuring number unless real r_hats were actually read.
% nParams - the count of data rows carrying a FINITE r_hat (i.e. real sampled
%  parameters). The Inf sentinel and NaN rows are excluded, so nParams == 0 flags
%  "no usable r_hat was found" and can gate a test independently of maxRhat.
%
%This is a read-only diagnostic helper for the live-recovery integration tests:
%those tests call psyrat_computevarcomp directly (bypassing the
%psyrat_computevarcompwarp -> psyrat_checkconv path), so they need to surface the
%worst r_hat themselves before trusting a recovery. The pass/fail decision still
%uses psyrat_checkconv (same r_hat>=1.1 + n_eff criteria as production); this
%helper only supplies the number to report and a non-empty guard.

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

%default to the safe, non-reassuring answer: if nothing usable is found, callers
%see Inf (worst possible r_hat) and 0 parameters, which reads as "not converged".
maxRhat = Inf;
nParams = 0;

if ~iscell(convdata) || isempty(convdata)
    return;
end

%mirror psyrat_checkconv's shape test: a row count of 1 means convdata is a
%1xN container of per-fit tables; otherwise convdata is itself a single table.
[a,nloops] = size(convdata);
if a == 1
    tables = cell(1,nloops);
    for i = 1:nloops
        tables{i} = convdata{i};
    end
else
    tables = {convdata};
end

%collect every data-row r_hat across all fit tables
rhats = [];
for k = 1:numel(tables)
    tbl = tables{k};
    if ~iscell(tbl) || size(tbl,1) < 2 || size(tbl,2) < 3
        continue;   %empty / header-only / malformed table contributes nothing
    end
    col = tbl(2:end,3);   %the r_hat column, data rows only
    for r = 1:numel(col)
        v = col{r};
        if isnumeric(v) && isscalar(v)
            rhats(end+1) = v; %#ok<AGROW>
        end
    end
end

if isempty(rhats)
    return;   %keep the [Inf, 0] default
end

%worst r_hat keeps Inf (unextractable sentinel) but ignores NaN entries; the
%usable-parameter count is the number of FINITE r_hats.
maxRhat = max(rhats,[],'omitnan');
if isempty(maxRhat) || isnan(maxRhat)
    maxRhat = Inf;   %all-NaN column: no usable r_hat
end
nParams = sum(isfinite(rhats));

end
