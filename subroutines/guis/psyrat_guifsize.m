function fsize = psyrat_guifsize
%Default GUI font size for PsyRAT windows.
%
%Returns:
% fsize - text/control font size in points

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

baseSize = get(0,'DefaultUicontrolFontSize');

if isempty(baseSize) || ~isnumeric(baseSize) || ~isfinite(baseSize) || baseSize <= 0
    baseSize = get(0,'DefaultTextFontSize');
end

%Cap sizes to avoid clipping in fixed-position legacy GUIs while still
%respecting larger defaults where possible.
fsize = min(max(round(baseSize),11),13);

end
