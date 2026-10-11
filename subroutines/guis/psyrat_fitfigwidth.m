function figwidth = psyrat_fitfigwidth(figwidth)
%Cap a window width at the display so a wide table window opens on screen.
%
%   figwidth = psyrat_fitfigwidth(figwidth) returns figwidth unchanged, or
%   the primary display's width minus a 40 px margin on each side when
%   figwidth is wider than that.
%
%Why: the table windows are resize-off and compute every control from
%figwidth, so a capped window keeps its layout while its uitable becomes
%narrower than its columns, and the classic uitable then scrolls
%horizontally (G104, owner ruling 2026-09-17: cap at the display and scroll
%rather than run off it). Pair it with movegui(fig,'onscreen') right after
%the figure is created, which slides a window inside the display; movegui
%would also clamp a window larger than the screen, which this cap prevents
%for the width. Without a real display (a ScreenSize of a few pixels in some
%headless sessions) nothing is capped.

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

margin = 40;
scr = get(groot, 'ScreenSize');
if numel(scr) == 4 && scr(3) > 2 * margin + 200
    figwidth = min(figwidth, scr(3) - 2 * margin);
end

end
