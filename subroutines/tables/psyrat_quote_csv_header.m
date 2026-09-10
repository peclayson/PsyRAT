function line = psyrat_quote_csv_header(line)
%Quote one CSV header/provenance line so it survives a spreadsheet round-trip
%as a single field (B19).
%
% line = psyrat_quote_csv_header(line)
%
%A comma inside an unquoted CSV line splits it across spreadsheet columns, so
%a comma- or quote-bearing header line (the 'Chains: 3, Iterations: 600' line,
%a dataset filename with a comma, the comma-bearing provenance disclosures) is
%wrapped in double quotes with internal quotes doubled (RFC-4180 escaping).
%Comma- and quote-free lines are returned untouched, so every existing
%comma-free header line stays byte-identical. This mirrors the inline CSV
%convention the dynrel viewer shipped first (psyrat_startview_dynrel.m, grep
%B19); this shared copy serves the table writers, the criterion export, the
%splits viewer, and psyrat_fprintf_provenance so the quoting rule has one
%definition. The XLSX branches (writecell) already quote at write time and
%never call this. Data rows and column-name rows must NOT be routed through
%this function: their commas are field separators.
%
%Input
% line - char row vector, one header line WITHOUT its trailing newline.
%
%Output
% line - the same line, quoted per RFC 4180 only when it contains a comma or
%        a double quote; otherwise returned byte-identical.

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

if contains(line, ',') || contains(line, '"')
    line = ['"' strrep(line, '"', '""') '"'];
end

end
