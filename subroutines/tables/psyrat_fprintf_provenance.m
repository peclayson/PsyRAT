function psyrat_fprintf_provenance(fid, psyrat_data)
%Write the standard run-provenance lines to an open file id (CSV export header).
%
% psyrat_fprintf_provenance(fid, psyrat_data)
%
%Convenience wrapper that prints psyrat_provenance_lines(psyrat_data) -- the
%seed, estimation engine, prior status, convergence verdict, and divergent-
%transition count -- one line per row to an already-open text file. The XLSX
%export paths append the same psyrat_provenance_lines output to their header
%cell directly; this keeps the CSV (fprintf) paths a single call so the two
%formats carry identical provenance. The underlying helper is fail-open (an
%absent field yields a "not recorded"/"not assessed" line, never an error), so
%adding this call can never break a table export. Comma-bearing lines are
%conditionally quoted (B19, psyrat_quote_csv_header) so each provenance line
%stays ONE spreadsheet field; comma-free lines stay byte-identical.
%
%Input
% fid        - file id from fopen on the open CSV being written.
% psyrat_data- toolbox data structure (passed through to psyrat_provenance_lines).

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

prov = psyrat_provenance_lines(psyrat_data);
for k = 1:numel(prov)
    fprintf(fid, '%s\n', psyrat_quote_csv_header(prov{k}));
end

end
