function psyrat_write_table(filepath, filehead, datatable)
%Write a header block followed by a data table to .xlsx or .csv, no dialogs.
%
%psyrat_write_table(filepath, filehead, datatable)
%
%This is the headless counterpart of the per-table GUI save callbacks: it takes
%an already-built header block and table and writes them to a caller-supplied
%path, so a scripted/batch report can export tables without any uiputfile
%prompt. It is used by psyrat_report. The interactive GUI exporters keep their
%own (table-specific) write paths; this writer provides a single consistent
%format for the programmatic path.
%
%Input
% filepath  - full output path; the extension selects the format ('.xlsx' or
%             '.csv', case-insensitive).
% filehead  - cell array of char rows written above the table (provenance /
%             header lines). May be empty ({}).
% datatable - a MATLAB table written below the header block.
%
%The table is written via writetable, so the column names form its header row;
%the header block is written above it via writecell. A zero-row table written
%to .xlsx gets its column-name row from writecell instead (see the .xlsx
%branch), so it still sits directly below the header block. An unsupported
%extension errors.

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
%

if nargin < 2 || isempty(filehead)
    filehead = {};
end
% Normalize the header block to a cell column of char rows.
filehead = filehead(:);

[~,~,ext] = fileparts(filepath);
ext = lower(ext);

switch ext
    case '.xlsx'
        % writecell/writetable overwrite cells top-down but do NOT clear cells
        % below the new content, so a re-run producing FEWER rows than a prior
        % run into the same path would leave the prior run's stale trailing rows
        % in the sheet. Delete any existing file first so the overwrite fully
        % truncates. (The .csv branch is safe: writecell recreates the file.)
        if isfile(filepath)
            delete(filepath);
        end
        if height(datatable) == 0 && width(datatable) > 0
            % A zero-row table is only its column-name row. writetable ignores
            % 'Range' for such a table (observed on R2025b) and puts the names
            % in A1, over the first header line, so the header lines are
            % written first and the names are written explicitly one row below
            % them (in A1 when there is no header). No data row is written.
            if ~isempty(filehead)
                writecell(filehead, filepath);
            end
            writecell(local_column_names(datatable), filepath, ...
                'Range', sprintf('A%d', numel(filehead)+1));
        elseif ~isempty(filehead)
            writecell(filehead, filepath);
            % Start the table one row below the last header row so nothing is
            % overwritten.
            writetable(datatable, filepath, 'Range', sprintf('A%d', numel(filehead)+1));
        else
            writetable(datatable, filepath);
        end
    case '.csv'
        if ~isempty(filehead)
            writecell(filehead, filepath);
            % Append the table beneath the header block. WriteMode 'append'
            % defaults WriteVariableNames to false, so force it on to keep the
            % column-name header row.
            writetable(datatable, filepath, 'WriteMode', 'append', ...
                'WriteVariableNames', true);
        else
            writetable(datatable, filepath);
        end
    otherwise
        error('psyrat_write_table:ext', ...
            ['Unsupported table output extension ''%s''. ', ...
            'How to fix: use a path ending in .xlsx or .csv.'], ext);
end
end

function names = local_column_names(datatable)
%The column-name row writetable writes for DATATABLE, as a 1xN cell: one name
%per column, with the columns of a multicolumn variable named <name>_1,
%<name>_2, ... as writetable names them.
vn = datatable.Properties.VariableNames;
names = {};
for k = 1:numel(vn)
    x = datatable.(vn{k});
    ncol = size(x, 2);
    if ischar(x) || ncol == 1
        names{end+1} = vn{k}; %#ok<AGROW>
    else
        for c = 1:ncol
            names{end+1} = sprintf('%s_%d', vn{k}, c); %#ok<AGROW>
        end
    end
end
end
