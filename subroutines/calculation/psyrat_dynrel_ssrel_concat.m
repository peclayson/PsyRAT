function T = psyrat_dynrel_ssrel_concat(summ)
%Stack the per-stratum dynamic-reliability per-participant tables into one.
%
% T = psyrat_dynrel_ssrel_concat(summ)
%
%Single source of truth for the dynamic-reliability per-participant
%(subject-level) export: each stratum's ssrel_table, prepended with the stratum
%label. Only the subject-level variants populate ssrel_table; strata without one
%are skipped. Used by BOTH the dynrel viewer (psyrat_startview_dynrel) and the
%headless report (psyrat_report's local_dynrel_tables) so the two match.
%
%Input
% summ - the struct returned by psyrat_dynrel_summary (per-stratum .ssrel_table,
%        empty for non-subject-level variants).
%
%Output
% T - the stacked table, or [] when no stratum has a per-participant table (so
%     callers can test isempty and omit the table). isempty([]) and isempty of
%     an empty table both hold, so existing isempty guards are unaffected.

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
keep = false(numel(summ.strata),1);
for s = 1:numel(summ.strata)
    st = summ.strata(s);
    sr = st.ssrel_table;
    if isempty(sr); continue; end
    sr.stratum = repmat(string(st.label),height(sr),1);
    parts{s} = movevars(sr,'stratum','Before',sr.Properties.VariableNames{1});
    keep(s) = true;
end
if any(keep)
    T = vertcat(parts{keep});
else
    T = [];
end

end
