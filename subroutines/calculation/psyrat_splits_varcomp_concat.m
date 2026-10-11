function T = psyrat_splits_varcomp_concat(summ)
%Stack the per-stratum data-splits variance-component tables into one table.
%
% T = psyrat_splits_varcomp_concat(summ)
%
%Single source of truth for the data-splits variance-component export: each
%stratum's variance-component table, prepended with the stratum label and its
%design. Case23 retains n_s/n_o/nbar_i columns; case24 explicitly names scalar
%reference counts and adds the observed incidence and projection metadata.
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
    if strcmp(summ.analysis,'trt_splits')
        % Reference counts are not effective incidence counts. Preserve the
        % complete observed weights/counts and explicit projection provenance
        % in a portable JSON column shared by GUI and headless exports.
        vc.Properties.VariableNames{'n_splits'} = 'reference_n_splits';
        vc.Properties.VariableNames{'n_occasions'} = 'reference_n_occasions';
        vc.Properties.VariableNames{'nbar_items'} = 'pooled_harmonic_items';
        vc.design_mode = repmat(string(st.design.mode),h,1);
        vc.score_rule = repmat(string(st.design.score_rule),h,1);
        vc.design_metadata = repmat(string(jsonencode(st.design)),h,1);
        parts{s} = movevars(vc,{'stratum','reference_n_splits',...
            'reference_n_occasions','pooled_harmonic_items','design_mode',...
            'score_rule','design_metadata'},'Before','component');
    else
        parts{s} = movevars(vc,{'stratum','n_splits','n_occasions','nbar_items'},...
            'Before','component');
    end
end
T = vertcat(parts{:});

end
