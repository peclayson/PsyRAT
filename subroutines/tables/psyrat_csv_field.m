function value = psyrat_csv_field(value)
%Quote one CSV text field so it survives a spreadsheet round-trip as a single
%field (G109).
%
% value = psyrat_csv_field(value)
%
%A participant ID, group or event name, or other text written as one CSV field
%(a value in a data row, or a line that is a single field, such as a section
%label or an empty-status line that names a group) is wrapped in double quotes,
%with internal quotes doubled (RFC 4180), when it contains a comma, a double
%quote, a carriage return or a line feed. Any other value is returned
%unchanged, so plain IDs and labels stay byte-identical. This is the rule
%psyrat_depoverallt's Save IDs writer already applies to its identity fields.
%psyrat_quote_csv_header applies the same convention (comma, double quote, CR
%or LF) to whole header lines.
%
%Input
% value - char, string scalar, or a 1x1 cell holding one.
%
%Output
% value - char row, quoted only when it needs to be.

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

value = char(string(value));
if contains(value,',') || contains(value,'"') || contains(value,newline) || ...
        contains(value,char(13))
    value = ['"' strrep(value,'"','""') '"'];
end

end
