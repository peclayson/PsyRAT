classdef LS27FigureFixture
    %LS27FIGUREFIXTURE Fixtures, figure routes and comparison utilities for the
    %Restricted disclosure on SAVED FIGURES (optional queue row 29: analyses 14,
    %25 and 27 under the gamma family, location-scale and log-nu alike). Used by
    %TestRestrictedFigureDisclosure and the local artifact writer.
    %
    %THE CONTRACT THE TESTS ENCODE (the names psyrat_restricted_figure_note
    %writes, kept in one place):
    %  TAG   Tag of the disclosure object on a figure of a Restricted arm. It is
    %        a graphics TEXTBOX annotation, never a uicontrol, because
    %        exportgraphics drops UI components from PNG/PDF exports while
    %        exportgraphics, print (saveas) and savefig all keep annotations.
    %  META  Top-level UserData field on the same figure: a struct with label
    %        'RESTRICTED', analysis (the exact rel.analysis label) and text
    %        (the disclosure text), so a saved FIG carries the disclosure in
    %        machine-readable form. Every other figure keeps its UserData.
    %  The disclosure text comes from the one text source,
    %  psyrat_gamma_twofacet_restriction: its figure field when one exists,
    %  otherwise its viewer-note field.
    %
    %THE SIX ARMS AND THEIR REAL ROUTES.
    %  14 and 27: psyrat_startview_dynrel -> "Plot reliability surface" ->
    %    psyrat_dynrelplot (one figure; for 27 the participants are overlaid);
    %    saved through the figure's own "Save figure (PNG/PDF)" button
    %    (uiputfile + exportgraphics) or as a FIG.
    %  25: psyrat_startview_trt -> Generate -> psyrat_relfigures('analysis',
    %    'trt') -> psyrat_relsummary 'trt_sserr' -> psyrat_ssrelplot 'dep' and
    %    'icc' (two figures, the only figures this design draws); saved with
    %    MATLAB's figure menu (saveas/print) or as a FIG.
    %
    %FIXTURES are the repository's deterministic draws, not new data:
    %LS27WordingFixture's gamma analyses 14, 27 and 26 (log-nu and
    %location-scale), PsyRATTestDataFactory's Gaussian 14 and 27, and its
    %analysis-25 result labeled gamma with gammascale 1 or 2. The 'trt_sserr'
    %summary branch and psyrat_ssrelplot read only the shared slots, never the
    %family, so the family label is the only input the gate can use there.
    %twoDimData extends the same four Restricted dynamic arms to two
    %dimensions and two strata (see its help), so the filled-contour surface
    %with its colorbars is exercised without any file outside the repository.
    %
    %LAYOUT PRESERVATION (the helper's GEOMETRY contract). The note sits in a
    %band added above the panels, so a Restricted figure is its silent twin
    %(the same numbers under a non-Restricted label) plus that band: every
    %axes, colorbar, existing caption and button keeps its pixel box measured
    %from the figure's lower-left corner (box, GEOM_TOL_PX), and the figure
    %grows by the band alone.

    properties (Constant)
        TAG = 'psyrat_restricted_disclosure';
        META = 'design_maturity';
        TARGET_TAG = 'psyrat_dynrel_target_caption';
        %minimum clearance (pixels) between the disclosure and any axes'
        %title, tick or axis-label region, and between it and other objects
        CLEAR_PX = 4;
        %tolerance (pixels) when a Restricted figure's layout is compared with
        %its silent twin's. Both are measured from the figure's lower-left
        %corner, and the helper rescales normalized positions by oldh/newh,
        %which is exact up to floating-point rounding (about 1e-13 px). One
        %pixel allows MATLAB to snap a position it computes during layout (a
        %colorbar's) to whole pixels differently in two figures of different
        %heights. It is far below any real shift: a figure that failed to grow
        %by its 63-px band would squeeze a 380-px panel by about 40 px, and the
        %smallest designed gap (CLEAR_PX) is four times larger.
        GEOM_TOL_PX = 1;
    end

    methods (Static)

        % ------------------------------------------------------------------
        % Text source and labels
        % ------------------------------------------------------------------

        function txt = expectedText()
            %The disclosure text from the one text source.
            R = psyrat_gamma_twofacet_restriction();
            if isfield(R, 'figure')
                txt = LS27FigureFixture.joinText(R.figure);
            else
                txt = LS27FigureFixture.joinText(R.screen);
            end
        end

        function lab = label(arm)
            %The exact rel.analysis label of an arm name ('g27ls', 'g25lognu', ...).
            switch arm(1:3)
                case 'g14'; lab = 'ic_dynrel_trt';
                case 'g27'; lab = 'ic_dynrel_sserrvar_trt';
                case 'g25'; lab = 'trt_sserrvar';
                otherwise; error('LS27FigureFixture:arm', 'Not a Restricted arm: %s', arm);
            end
        end

        % ------------------------------------------------------------------
        % Dynamic designs (14, 27 and their neighbors)
        % ------------------------------------------------------------------

        function pd = dynrelData(arm)
            %Raw dynamic result for an arm: g14ls, g14lognu, g27ls, g27lognu
            %(Restricted); g26ls, g26lognu, g11ls (gamma one-facet neighbors);
            %n14, n27 (Gaussian arms of the Restricted designs).
            switch arm
                case 'g14ls';    pd = LS27WordingFixture.dyn14('ls');
                case 'g14lognu'; pd = LS27WordingFixture.dyn14('lognu');
                case 'g27ls';    pd = LS27WordingFixture.dyn27('ls', 'constant');
                case 'g27lognu'; pd = LS27WordingFixture.dyn27('lognu', 'constant');
                case 'g26ls';    pd = LS27WordingFixture.dyn26('ls');
                case 'g26lognu'; pd = LS27WordingFixture.dyn26('lognu');
                case 'g11ls'
                    pd = LS27WordingFixture.dyn26('ls');
                    pd.rel.analysis = 'ic_dynrel';
                    drop = intersect(fieldnames(pd.rel.out), {'ssinfo'; 'ind_bs'; 'ind_nu'; 'ind_sdlog'});
                    if ~isempty(drop); pd.rel.out = rmfield(pd.rel.out, drop); end
                case 'n14'; pd = PsyRATTestDataFactory.makeDynrelTrtSummaryData();
                case 'n27'; pd = PsyRATTestDataFactory.makeDynrelSubjNondiffTrtSummaryData();
                otherwise; error('LS27FigureFixture:arm', 'Unknown dynamic arm: %s', arm);
            end
        end

        function summ = dynrelSummary(arm)
            pd = LS27FigureFixture.dynrelData(arm);
            summ = psyrat_dynrel_summary(pd.rel, 'CI', .95, 'ngrid', 5);
        end

        function pd = twoDimData(arm)
            %A two-dimension, two-stratum variant of a Restricted dynamic arm
            %('g14ls', 'g14lognu', 'g27ls', 'g27lognu'), built from the same
            %repository draws as dynrelData. It reaches the two layouts the
            %one-dimension arms never do: the filled-contour surface with one
            %colorbar per panel (psyrat_dynrelplot.m:109-131), and a 1 x 2
            %panel grid.
            %  Dimension 2: a second standardized value per participant, with
            %    slopes on z2 and on z1*z2 (the design row is [z1 z2 z1*z2])
            %    for the log-mean and for the scale submodel (b_sigma under
            %    location-scale, b_nu under log-nu), so the surface varies
            %    over both axes.
            %  Stratum B: every per-stratum slot duplicated, with the log-mean
            %    intercept raised by log 1.3 and the reference trial count set
            %    to 11, so the two panels differ; for analysis 27 its three
            %    participants are new IDs in another row order with other
            %    trial counts.
            %The draws stay unstamped, as in dynrelData.
            if ~any(strcmp(arm, {'g14ls', 'g14lognu', 'g27ls', 'g27lognu'}))
                error('LS27FigureFixture:arm', 'Not a Restricted dynamic arm: %s', arm);
            end
            pd = LS27FigureFixture.dynrelData(arm);
            o = pd.rel.out;
            nd = size(o.mu, 1);
            %---- dimension 2, stratum A --------------------------------------
            z2 = [-0.7; 0.8; 0.2];
            o.dimz{1} = [o.dimz{1}(:), z2];
            if isfield(o, 'ssinfo'); o.ssinfo{1}.z2 = z2; end
            o.b{1} = num2cell([cell2mat(o.b{1}), repmat([0.07 -0.05], nd, 1)]);
            if isfield(o, 'b_sigma'); scale = 'b_sigma'; else; scale = 'b_nu'; end
            o.(scale){1} = num2cell([cell2mat(o.(scale){1}), repmat([0.10 0.04], nd, 1)]);
            %---- stratum B: every per-stratum slot duplicated ---------------
            %(cell slots hold one cell per stratum; numeric slots one column)
            names = fieldnames(o);
            for k = 1:numel(names)
                x = o.(names{k});
                if iscell(x) || isnumeric(x)
                    o.(names{k}) = [x x];
                end
            end
            o.labels = {'Stratum A', 'Stratum B'};
            o.mu(:, 2) = o.mu(:, 2) + log(1.3);
            o.ntrials(2) = 11;
            if isfield(o, 'ssinfo')
                o.ssinfo{2} = o.ssinfo{1}([3 1 2], :);
                o.ssinfo{2}.id = [703; 701; 702];
                o.ssinfo{2}.trls = [7; 19; 12];
            end
            pd.rel.out = o;
            pd.rel.ndim = 2;
            pd.rel.dim_names = {'Dimension 1', 'Dimension 2'};
        end

        function summ = twoDimSummary(arm)
            pd = LS27FigureFixture.twoDimData(arm);
            summ = psyrat_dynrel_summary(pd.rel, 'CI', .95, 'ngrid', 5);
        end

        function s = relabel(s, lab)
            %The same summary under another analysis label: the only design
            %identifiers a dynamic summary carries are gamma_target.analysis and,
            %if an implementation adds one, a top-level analysis field. Every
            %number is untouched, so a figure drawn from it can differ from the
            %original only through a label-keyed element.
            s.gamma_target.analysis = lab;
            if isfield(s, 'analysis'); s.analysis = lab; end
        end

        % ------------------------------------------------------------------
        % Static analysis 25 and its neighbors
        % ------------------------------------------------------------------

        function pd = trt25(fam)
            %Analysis 25 (trt_sserrvar): 'ls' (gamma, gammascale 2), 'lognu'
            %(gamma, gammascale 1), 'gaussian', or 'nofamily' (a legacy result).
            pd = PsyRATTestDataFactory.makeTRTSSErrRelData();
            switch fam
                case 'ls';       pd.rel.family = 'gamma'; pd.rel.gammascale = 2;
                case 'lognu';    pd.rel.family = 'gamma'; pd.rel.gammascale = 1;
                case 'gaussian'; pd.rel.family = 'gaussian';
                case 'nofamily'
                otherwise; error('LS27FigureFixture:fam', 'Unknown family: %s', fam);
            end
        end

        function summ = trt25Summary(fam)
            %The 'trt_sserr' summary the analysis-25 subject plots draw from.
            [summ, ~] = psyrat_relsummary('psyrat_data', LS27FigureFixture.trt25(fam), ...
                'analysis', 'trt_sserr', 'gcoeff', 1, 'reltype', 3, ...
                'depcutoff', 0.7, 'meascutoff', 2, 'depcentmeas', 1, 'noccmode', 1);
        end

        function prefs = trtViewPrefs()
            %Test-retest viewer preferences with only the two subject plots on,
            %so psyrat_relfigures draws exactly the analysis-25 figures.
            v = struct();
            v.gcoeff = 1; v.reltype = 3; v.diffgcoeff = 1;
            v.noccmode = 1; v.nocc = [];
            v.relvalue = 0.5; v.relcentmeas = 1; v.meascutoff = 2;
            v.plotrel = 0; v.plotrelline = 2;
            v.ploticc = 0; v.inctrltable = 0; v.overalltable = 0;
            v.ntrials = 20; v.showstddevt = 0; v.showstddevf = 0;
            v.plotssrel = 1; v.plotssicc = 1; v.tablessrel = 0;
            prefs = struct('view', v, 'ver', '0-test');
        end

        function fig = subjectPlot(summ, stat)
            %One analysis-25 subject plot ('dep' or 'icc'); psyrat_ssrelplot
            %returns no handle, so the new output figure is found by difference.
            before = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
            psyrat_ssrelplot('psyrat_data', summ, 'stat', stat);
            fig = setdiff(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output'), before);
        end

        % ------------------------------------------------------------------
        % Reading a figure
        % ------------------------------------------------------------------

        function h = disclosure(fig)
            %Annotations live in the annotation pane: findall, not findobj.
            h = findall(fig, 'Tag', LS27FigureFixture.TAG);
        end

        function t = textOf(h)
            t = LS27FigureFixture.joinText(get(h, 'String'));
        end

        function t = joinText(s)
            %Lines joined with single spaces (a wrapped caption compares equal
            %to its one-line source text).
            if iscell(s); s = strjoin(cellfun(@char, s(:)', 'UniformOutput', false), ' '); end
            if isstring(s); s = strjoin(s(:)', ' '); end
            t = strtrim(regexprep(char(s), '\s+', ' '));
        end

        function tf = hasMeta(fig)
            ud = get(fig, 'UserData');
            tf = isstruct(ud) && isfield(ud, LS27FigureFixture.META);
        end

        function sig = signature(fig)
            %What a reader sees in the panels, independent of layout: the figure
            %name and, per axes in creation order, the title, axis labels, axis
            %limits and the data and colors of every plotted object. Annotations
            %are excluded on purpose, so two figures that differ only by a
            %caption or a disclosure compare equal.
            sig = struct('name', get(fig, 'Name'), 'axes', {{}});
            props = {'XData', 'YData', 'ZData', 'YNegativeDelta', 'YPositiveDelta', ...
                'CData', 'Vertices', 'Faces', 'SizeData', 'Color', ...
                'MarkerFaceColor', 'MarkerEdgeColor', 'LineStyle', 'LineWidth'};
            axs = findobj(fig, 'Type', 'axes');
            for k = 1:numel(axs)
                ax = axs(k);
                a = struct('title', get(get(ax, 'Title'), 'String'), ...
                    'xlabel', get(get(ax, 'XLabel'), 'String'), ...
                    'ylabel', get(get(ax, 'YLabel'), 'String'), ...
                    'xlim', get(ax, 'XLim'), 'ylim', get(ax, 'YLim'), 'objs', {{}});
                ch = get(ax, 'Children');
                for j = 1:numel(ch)
                    o = struct('type', get(ch(j), 'Type'));
                    for p = props
                        if isprop(ch(j), p{1}); o.(p{1}) = get(ch(j), p{1}); end
                    end
                    a.objs{end + 1} = o;
                end
                sig.axes{end + 1} = a;
            end
        end

        % ------------------------------------------------------------------
        % Geometry (pixels, figure-relative)
        % ------------------------------------------------------------------

        function b = box(h, fig)
            %[left bottom width height] in figure pixels, measured from the
            %figure's lower-left corner, for an annotation or a colorbar
            %(read in normalized figure units and scaled by the figure's pixel
            %size), an axes or a uicontrol. A position MATLAB computes during
            %layout (a colorbar's) is current only after drawnow.
            fp = getpixelposition(fig);
            if isa(h, 'matlab.graphics.shape.TextBox') || strcmp(get(h, 'Type'), 'colorbar')
                u = get(h, 'Units');
                if ~strcmp(u, 'normalized'); set(h, 'Units', 'normalized'); end
                p = get(h, 'Position');
                if ~strcmp(u, 'normalized'); set(h, 'Units', u); end
                b = [p(1) * fp(3), p(2) * fp(4), p(3) * fp(3), p(4) * fp(4)];
            else
                b = getpixelposition(h, true);
            end
        end

        function b = colorbarReadableBox(cb, fig)
            %A colorbar plus the space its tick labels and label take: the
            %bar's box, widened to the figure's right edge (the tick labels
            %and the rotated label of an eastoutside bar sit to its right) and
            %extended by one tick-label line above and below (the end tick
            %labels are centred on the bar's ends, so they reach half a line
            %past them). The rotated label is centred on the bar and, in these
            %figures, shorter than it, so the bar's height bounds it.
            fp = getpixelposition(fig);
            b = LS27FigureFixture.box(cb, fig);
            linepx = ceil(get(cb, 'FontSize') * get(groot, 'ScreenPixelsPerInch') / 72);
            b = [b(1), b(2) - linepx, fp(3) - b(1), b(4) + 2 * linepx];
        end

        function b = readableBox(ax)
            %An axes plus its title, x label and y label: the region the
            %disclosure must not cover.
            b = getpixelposition(ax, true);
            for t = [get(ax, 'Title'), get(ax, 'XLabel'), get(ax, 'YLabel')]
                s = get(t, 'String');
                if isempty(s) || (iscell(s) && all(cellfun(@isempty, s))); continue; end
                u = get(t, 'Units');
                set(t, 'Units', 'pixels');
                e = get(t, 'Extent');
                set(t, 'Units', u);
                eb = [b(1) + e(1), b(2) + e(2), e(3), e(4)];
                lo = min([b(1:2); eb(1:2)], [], 1);
                hi = max([b(1:2) + b(3:4); eb(1:2) + eb(3:4)], [], 1);
                b = [lo, hi - lo];
            end
        end

        function tf = overlaps(a, b, pad)
            %True when boxes a and b, each grown by pad pixels, intersect.
            if nargin < 3; pad = 0; end
            tf = a(1) - pad < b(1) + b(3) && b(1) - pad < a(1) + a(3) && ...
                a(2) - pad < b(2) + b(4) && b(2) - pad < a(2) + a(4);
        end

        % ------------------------------------------------------------------
        % Saved images
        % ------------------------------------------------------------------

        function n = ink(file)
            %Number of near-black pixels in a saved PNG: black text (captions,
            %labels, the disclosure) counts; the colored curves, bands and
            %markers do not.
            img = imread(file);
            if size(img, 3) == 1; img = repmat(img, [1 1 3]); end
            img = double(img);
            if max(img(:)) > 255; img = img / 257; end
            n = nnz(all(img < 96, 3));
        end

        function tf = isPdf(file)
            %A nonempty file that starts with the PDF signature.
            tf = false;
            fid = fopen(file, 'r');
            if fid < 0; return; end
            c = onCleanup(@() fclose(fid));
            head = fread(fid, 5, '*char')';
            tf = strcmp(head, '%PDF-');
        end

    end
end
