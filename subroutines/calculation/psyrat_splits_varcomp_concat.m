function T = psyrat_splits_varcomp_concat(summ)
%Stack the per-stratum data-splits variance-component tables into one table.
%
% T = psyrat_splits_varcomp_concat(summ)
%
%Single source of truth for the data-splits variance-component export: each
%stratum's variance-component table, prepended with the stratum label and its
%observed design (n_s splits, n_o occasions, harmonic-mean items-per-split
%n_i = nbar_i), so the export records the design each component was reported at.
%Used by BOTH the splits viewer (psyrat_startview_splits) and the headless
%report (psyrat_report's local_splits_tables) so the two stay identical.
%
%Input
% summ - the struct returned by psyrat_splits_summary (per-stratum .varcomp
%        table plus .nsplit/.nocc/.nbar_i design fields).
%
%Output
% T - the stacked table (stratum, n_splits, n_occasions, nbar_items, then the
%     original variance-component columns).

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

parts = cell(numel(summ.strata),1);
for s = 1:numel(summ.strata)
    st = summ.strata(s);
    vc = st.varcomp;
    h = height(vc);
    vc.stratum = repmat(string(st.label),h,1);
    vc.n_splits = repmat(st.nsplit,h,1);
    vc.n_occasions = repmat(st.nocc,h,1);
    vc.nbar_items = repmat(st.nbar_i,h,1);
    parts{s} = movevars(vc,{'stratum','n_splits','n_occasions','nbar_items'},...
        'Before','component');
end
T = vertcat(parts{:});

end
