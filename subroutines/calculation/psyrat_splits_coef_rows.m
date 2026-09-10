function [data, cols, widths] = psyrat_splits_coef_rows(summ, L)
%Flatten a nonparallel data-splits summary into reliability-coefficient rows.
%
% [data, cols, widths] = psyrat_splits_coef_rows(summ, L)
%
%Single source of truth for the data-splits coefficient table: one row per
%stratum x reported coefficient. Used by BOTH the splits results viewer
%(psyrat_startview_splits) and the headless report (psyrat_report's
%local_splits_tables), so the on-screen table, the GUI export, and the headless
%export are identical by construction.
%
%Dependability D and generalizability G are each [ll pt ul]; G and the relative
%SEM are the downward-biased lower/upper bounds (the column headers say so).
%
%Input
% summ - the struct returned by psyrat_splits_summary (per-stratum .coef array).
% L    - the analysis-unit label struct from psyrat_unitlabels (its .units field
%        names the per-person count column 'n splits'/'n trials').
%
%Output
% data   - ncoef-by-12 cell array of row values (Stratum first).
% cols   - 1-by-12 cell array of column header strings.
% widths - 1-by-12 cell array of uitable column widths (GUI display only; the
%          headless caller ignores it).
%
%Example
% L = psyrat_unitlabels(REL.splits);
% [data,cols] = psyrat_splits_coef_rows(summ, L);
% T = cell2table(data,'VariableNames',matlab.lang.makeValidName(cols));

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

%preallocate (each stratum contributes numel(.coef) rows) rather than growing
ncoef = 0;
for s = 1:numel(summ.strata)
    ncoef = ncoef + numel(summ.strata(s).coef);
end
data = cell(ncoef, 12);
r = 0;
for s = 1:numel(summ.strata)
    st = summ.strata(s);
    for c = 1:numel(st.coef)
        k = st.coef(c);
        r = r + 1;
        data(r,:) = {
            char(string(st.label)), k.label, st.nsplit, st.nocc, ...
            round(k.D(2),3), round(k.D(1),3), round(k.D(3),3), ...
            round(k.G(2),3), round(k.G(1),3), round(k.G(3),3), ...
            round(k.SEM_abs(2),3), round(k.SEM_rel(2),3)};
    end
end
cols = {'Stratum','Coefficient',['n ' L.units],'n occasions',...
    'Dependability','D low','D high',...
    'Generalizability (lower bound)','G low','G high',...
    'SEM (abs)','SEM (rel, upper bound)'};
widths = {110 200 70 80 100 60 60 150 60 60 80 130};

end
