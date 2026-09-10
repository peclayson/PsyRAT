function psyrat_autofitbuttontext(hObject,~)
%Normalize uicontrol labels and shrink font size so text fits control bounds.

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

if ~ishghandle(hObject,'uicontrol')
    return;
end

style = lower(get(hObject,'Style'));
fitStyles = {'pushbutton','push','text','checkbox','radiobutton'};
if ~any(strcmp(style,fitStyles))
    return;
end

oldUnits = 'pixels';
try
    oldUnits = get(hObject,'Units');
    set(hObject,'Units','pixels');

    pos = get(hObject,'Position');
    if pos(3) <= 0 || pos(4) <= 0
        set(hObject,'Units',oldUnits);
        return;
    end

    str = psyrat_cleanupcontrolstring(get(hObject,'String'));
    set(hObject,'String',str);
    if isempty(str) || (iscell(str) && all(cellfun(@isempty,str)))
        set(hObject,'Units',oldUnits);
        return;
    end

    origSize = get(hObject,'FontSize');
    if isempty(origSize) || ~isnumeric(origSize) || ~isfinite(origSize) || origSize <= 0
        origSize = get(0,'DefaultUicontrolFontSize');
        if isempty(origSize) || ~isnumeric(origSize) || ~isfinite(origSize) || origSize <= 0
            origSize = 11;
        end
        set(hObject,'FontSize',origSize);
    end

    minSize = max(8,floor(origSize*0.65));
    targetW = max(pos(3) - 8, 1);
    targetH = max(pos(4) - 4, 1);

    for s = origSize:-0.5:minSize
        set(hObject,'FontSize',s);
        ext = get(hObject,'Extent');
        if ext(3) <= targetW && ext(4) <= targetH
            break;
        end
    end

    set(hObject,'Units',oldUnits);
catch
    if ishghandle(hObject,'uicontrol') && exist('oldUnits','var')
        try
            set(hObject,'Units',oldUnits);
        catch
            % If units cannot be reset, preserve callback safety.
        end
    end
end

end

function out = psyrat_cleanupcontrolstring(in)
%Convert legacy HTML uicontrol labels to plain text and preserve line breaks.

if iscell(in)
    parts = cellfun(@psyrat_coercetext,in,'UniformOutput',false);
    txt = strjoin(parts,sprintf('\n'));
else
    txt = psyrat_coercetext(in);
end

if contains(lower(txt),'<html>')
    txt = regexprep(txt,'(?i)<br\\s*/?>',sprintf('\n'));
    txt = regexprep(txt,'(?i)</p>|</tr>|</td>',sprintf('\n'));
    txt = regexprep(txt,'(?i)<[^>]+>','');
    txt = strrep(txt,'&nbsp;',' ');
end

txt = strtrim(txt);

if contains(txt,sprintf('\n'))
    lines = regexp(txt,'\n','split');
    lines = lines(~cellfun(@isempty,lines));
    if isempty(lines)
        out = '';
    elseif numel(lines) == 1
        out = lines{1};
    else
        out = lines(:);
    end
else
    out = txt;
end

end

function out = psyrat_coercetext(in)
%Convert numeric/string inputs to a character vector for cleanup.

if isempty(in)
    out = '';
elseif isnumeric(in)
    out = num2str(in);
elseif isstring(in)
    if isscalar(in)
        out = char(in);
    else
        out = strjoin(cellstr(in(:)),sprintf('\n'));
    end
elseif ischar(in)
    if size(in,1) > 1
        out = strjoin(cellstr(in),sprintf('\n'));
    else
        out = in;
    end
else
    try
        out = char(in);
    catch
        out = '';
    end
end

end
