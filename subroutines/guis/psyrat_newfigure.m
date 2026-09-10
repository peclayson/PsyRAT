function h = psyrat_newfigure(varargin)
%Create a PsyRAT figure with automatic pushbutton text fitting.
%
%This keeps button labels readable across operating systems and display
%scales by shrinking the button font only when needed.

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

h = figure(varargin{:},'DefaultUicontrolCreateFcn',@psyrat_autofitbuttontext);

end
