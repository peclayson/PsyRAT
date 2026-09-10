function h = psyrat_addhline(y,varargin)
%Add a horizontal reference line at y using base MATLAB graphics.
%
% h = psyrat_addhline(y)
% h = psyrat_addhline(y,'Name',Value,...)
%
% Input:
% y - y-value for horizontal line
%
% Output:
% h - line object handle

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

ax = gca;
xl = xlim(ax);
was_hold = ishold(ax);
hold(ax,'on');
h = line(ax,xl,[y y]);
if ~isempty(varargin)
    set(h,varargin{:});
end
if ~was_hold
    hold(ax,'off');
end

end
