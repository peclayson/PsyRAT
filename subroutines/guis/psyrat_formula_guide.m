function psyrat_formula_guide(varargin)
%Display formulas from Tables 2, 3, and 6 with runnable usage templates.
%
%psyrat_formula_guide
%psyrat_formula_guide('analysis','sing')

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

analysis = 'all';

if ~isempty(varargin)
    for i = 1:length(varargin)-1
        if (ischar(varargin{i}) || (isstring(varargin{i}) && isscalar(varargin{i}))) && ...
                strcmpi(char(varargin{i}),'analysis')
            analysis = varargin{i+1};
        end
    end
end

catalog = psyrat_formula_catalog('analysis',analysis);

fsize = psyrat_guifsize;
fig = psyrat_newfigure('unit','pix',...
    'position',[120 120 1260 720],...
    'menub','none',...
    'name','Formula Guide (Tables 2, 3, and 6)',...
    'numbertitle','off',...
    'resize','off');
fig.Tag = 'psyrat_output';

uicontrol(fig,'Style','text',...
    'fontsize',fsize+2,...
    'HorizontalAlignment','left',...
    'String','Formula Guide: formulas, implementation status, and runnable templates',...
    'Position',[20 680 1080 24]);

uicontrol(fig,'Style','text',...
    'fontsize',fsize-1,...
    'HorizontalAlignment','left',...
    'String','Select a row, then use Show Usage / Run Formula. Cut-score formulas use criterion as ''cut''.',...
    'Position',[20 658 1220 18]);

tableData = table2cell(catalog(:,{'paper_table','design','coefficient',...
    'decision_type','formula_ascii','psyrat_function','gui_status'}));

hTable = uitable(fig,...
    'Data',tableData,...
    'ColumnName',{'Table','Design','Coefficient','Decision Type',...
    'Formula','Function','GUI Status'},...
    'ColumnWidth',{50 190 170 150 390 150 230},...
    'RowName',[],...
    'Position',[20 280 1220 365],...
    'CellSelectionCallback',@psyrat_selectrow);

hTable.ColumnFormat = {'char','char','char','char','char','char','char'};

uicontrol(fig,'Style','push',...
    'fontsize',max(fsize-1,8),...
    'String','Show Usage',...
    'Position',[20 235 120 34],...
    'Callback',@psyrat_showusage);

uicontrol(fig,'Style','push',...
    'fontsize',max(fsize-1,8),...
    'String','Run Formula',...
    'Tooltip','Opens a small input form and runs the selected formula using entered values',...
    'Position',[150 235 120 34],...
    'Callback',@psyrat_runformula);

uicontrol(fig,'Style','push',...
    'fontsize',max(fsize-1,8),...
    'String','Copy Command',...
    'Tooltip','Copy the current usage command to clipboard',...
    'Position',[280 235 130 34],...
    'Callback',@psyrat_copycmd);

missing_impl = strcmp(catalog.implementation_status,'Needs implementation');
if any(missing_impl)
    implSummary = sprintf('Missing implementation rows: %d',sum(missing_impl));
else
    implSummary = 'Implementation coverage: all formulas in this view are implemented.';
end

uicontrol(fig,'Style','text',...
    'fontsize',fsize-1,...
    'HorizontalAlignment','left',...
    'String',implSummary,...
    'Position',[430 240 810 24]);

usageBox = uicontrol(fig,'Style','edit',...
    'Min',0,...
    'Max',2,...
    'Enable','inactive',...
    'HorizontalAlignment','left',...
    'fontsize',fsize-1,...
    'String','Select a formula row and click Show Usage.',...
    'Position',[20 20 1220 205]);

state.selected_row = 1;
state.last_command = '';
setappdata(fig,'formula_state',state);
setappdata(fig,'formula_catalog',catalog);
setappdata(fig,'formula_usage_box',usageBox);

%initialize usage text for first row
psyrat_showusage([],[]);

    function psyrat_selectrow(~,evt)
        if isempty(evt.Indices)
            return;
        end
        state = getappdata(fig,'formula_state');
        state.selected_row = evt.Indices(1);
        setappdata(fig,'formula_state',state);
    end

    function psyrat_showusage(~,~)
        catalogLocal = getappdata(fig,'formula_catalog');
        usageHandle = getappdata(fig,'formula_usage_box');
        state = getappdata(fig,'formula_state');
        
        row = min(max(state.selected_row,1),height(catalogLocal));
        fid = catalogLocal.formula_id{row};
        [tmpl,titleText] = psyrat_formula_template('formula_id',fid);
        
        msg = sprintf('%s\n\n%s',titleText,tmpl);
        usageHandle.String = msg;
        state.last_command = tmpl;
        setappdata(fig,'formula_state',state);
    end

    function psyrat_copycmd(~,~)
        state = getappdata(fig,'formula_state');
        if isempty(state.last_command)
            psyrat_showusage([],[]);
            state = getappdata(fig,'formula_state');
        end
        
        try
            clipboard('copy',state.last_command);
        catch
            errordlg('Could not access clipboard on this system.','Clipboard Error','modal');
        end
    end

    function psyrat_runformula(~,~)
        catalogLocal = getappdata(fig,'formula_catalog');
        usageHandle = getappdata(fig,'formula_usage_box');
        state = getappdata(fig,'formula_state');
        
        row = min(max(state.selected_row,1),height(catalogLocal));
        fid = catalogLocal.formula_id{row};
        
        try
            out = psyrat_formula_prompt('formula_id',fid);
        catch ME
            errordlg(ME.message,'Formula Error','modal');
            return;
        end
        if ~out.ok
            return;
        end
        
        usageHandle.String = out.summary;
        state.last_command = out.command;
        setappdata(fig,'formula_state',state);
    end

end
