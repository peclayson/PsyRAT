function h = psyrat_restricted_figure_note(fig, analysis, family)
%PSYRAT_RESTRICTED_FIGURE_NOTE Put the RESTRICTED design-maturity note on an
% output figure of the three gamma two-facet designs (analyses 14, 25 and 27
% under the gamma family, location-scale and log-nu alike), and return its
% handle. For every other figure it returns [] and changes nothing.
%
%   h = psyrat_restricted_figure_note(fig, analysis, family)
%
% fig      - the output figure: the dynamic surface of psyrat_dynrelplot (14,
%            27) or a subject plot of psyrat_ssrelplot (25)
% analysis - the result's rel.analysis label
% family   - the result's rel.family; only 'gamma' qualifies
%
% WHAT IT ADDS.
%  (1) A band ABOVE the panels holding a graphics TEXTBOX annotation, Tag
%      'psyrat_restricted_disclosure', with the saved-figure note R.figure of
%      psyrat_gamma_twofacet_restriction (the one text source). It is an
%      annotation, not a uicontrol, because exportgraphics drops UI components
%      from PNG and PDF exports; print/saveas and savefig keep annotations.
%  (2) UserData.design_maturity, a struct with label, analysis, family, text
%      and source, added beside whatever UserData the figure already carries
%      (the dynamic figure's dynrel/summary fields stay as they are), so a
%      saved FIG keeps the note in machine-readable form. A figure whose
%      UserData is not a scalar struct keeps it untouched.
%
% GEOMETRY. The figure grows by the band's height with its top edge fixed on
% screen. Every axes and every existing annotation keeps its pixel size and
% its distance from the figure's bottom edge, and uicontrols (pixel units, such
% as the surface figure's Save button) do not move. So panels, titles, axis
% labels, captions and buttons keep their layout, and the band sits above
% them. Colorbars follow their axes. The band is sized for the default figure
% sizes (950 x 520 and 900 x 450 px) at 10 pt; long text, another screen DPI or
% a narrower window can crowd it, so the rendered figure is read in the native
% look.
%
% GATE. Exact label AND the gamma family. The analysis labels nest as prefixes
% of one another (see psyrat_gamma_twofacet_restriction), so a prefix or
% contains test would put the note on a design it does not describe; the
% Gaussian arms of the same designs, the one-facet gamma dynamic designs and
% every other design are left unchanged.
%
% See also PSYRAT_GAMMA_TWOFACET_RESTRICTION, PSYRAT_DYNRELPLOT, PSYRAT_SSRELPLOT.

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

h = [];

%---- gate: gamma family AND an exact Restricted label -------------------------
if ~(ischar(analysis) || (isstring(analysis) && isscalar(analysis))) || ...
        ~(ischar(family) || (isstring(family) && isscalar(family)))
    return;
end
analysis = char(analysis);
R = psyrat_gamma_twofacet_restriction();
if ~strcmpi(strtrim(char(family)), 'gamma') || ~any(strcmp(analysis, R.designs))
    return;
end

%---- wrap the note to the figure width ----------------------------------------
oldunits = get(fig, 'Units');
set(fig, 'Units', 'pixels');
restoreUnits = onCleanup(@() set(fig, 'Units', oldunits));
pos = get(fig, 'Position');
maxchars = max(40, floor(0.96 * pos(3) / 8));   %about 8 px per character at 10 pt
lines = local_wrap(R.figure, maxchars);
linepx = 17;                                   %line pitch for 10-pt text
band = numel(lines) * linepx + 12;             %pixels added above the panels

%---- grow the figure upward and keep every normalized object's pixel geometry --
oldh = pos(4);
newh = oldh + band;
axs = findall(fig, 'Type', 'axes');
for k = 1:numel(axs)
    set(axs(k), 'Units', 'normalized');
    p = get(axs(k), 'Position');
    set(axs(k), 'Position', [p(1) p(2)*oldh/newh p(3) p(4)*oldh/newh]);
end
everything = findall(fig);
notes = everything(arrayfun(@(o) isa(o, 'matlab.graphics.shape.TextBox'), everything));
for k = 1:numel(notes)
    set(notes(k), 'Units', 'normalized');
    p = get(notes(k), 'Position');
    set(notes(k), 'Position', [p(1) p(2)*oldh/newh p(3) p(4)*oldh/newh]);
end
set(fig, 'Position', [pos(1) pos(2)-band pos(3) newh]);

%---- the note itself ------------------------------------------------------------
h = annotation(fig, 'textbox', [0.02 (oldh + 4)/newh 0.96 (band - 8)/newh], ...
    'String', lines(:), 'Tag', 'psyrat_restricted_disclosure', 'FontSize', 10, ...
    'Interpreter', 'none', 'EdgeColor', 'none', 'Margin', 2, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FitBoxToText', 'off');

%---- the machine-readable record for a saved FIG ------------------------------
ud = get(fig, 'UserData');
if isempty(ud)
    ud = struct();
end
if isstruct(ud) && isscalar(ud)
    ud.design_maturity = struct('label', R.label, 'analysis', analysis, ...
        'family', 'gamma', 'text', R.figure, ...
        'source', 'psyrat_gamma_twofacet_restriction');
    set(fig, 'UserData', ud);
end
end

function out = local_wrap(str, maxchars)
%Greedy word wrap by character count. A character count is not a pixel width,
%so the rendered figure is read in the native look.
words = strsplit(str, ' ');
out = {};
cur = '';
for k = 1:numel(words)
    if isempty(cur)
        cand = words{k};
    else
        cand = [cur ' ' words{k}];
    end
    if numel(cand) > maxchars && ~isempty(cur)
        out{end + 1} = cur; %#ok<AGROW>
        cur = words{k};
    else
        cur = cand;
    end
end
if ~isempty(cur)
    out{end + 1} = cur;
end
end
