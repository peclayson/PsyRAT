classdef TestDodHeatmapMask < matlab.unittest.TestCase
    %TESTDODHEATMAPMASK The DoD trial-count heatmaps must not clip: a cell whose
    %projected coefficient is inadmissible (denominator <= 0, NaN) or outside
    %[0,1] is masked, not painted as 0 or 1.
    %
    % Finding (post-beta review 2026-09-10, PBA-07; work queue Q2): the grid
    % was clipped with max(min(z,1),0). MATLAB's two-input min and max ignore
    % NaN, so an inadmissible cell rendered as 1.0 (maximum reliability), and a
    % value outside [0,1] rendered as 0 or 1 with no mark. The adjacent
    % coefficient table is draw-level through psyrat_dodiffrel and never
    % clipped, so the two disagreed silently.
    %
    % Headless, on the real screens, the way TestPlotLabelInterpreter drives
    % the Generate DoD Outputs button. The fixture inflates group 1's event-2
    % trial covariance (variance 36, covariance 5.4 with event 1; still
    % positive semidefinite) so the harmonic-scaled pairwise term is strongly
    % negative where event-1 trials are few and event-2 trials are many, which
    % drives the projected dependability denominator below zero in that corner
    % of the first heatmap's grid.
    %
    % 2B / G11 (2026-09-13): the mask is now the cell's component-level
    % admissibility (the four quantities psyrat_dodiffrel records, strict
    % computed sign) with the range rule kept as a second layer, so it is a
    % SUPERSET of the range mask. The first test therefore asserts that every
    % range-bad cell is still masked and that the mask equals an INDEPENDENT
    % hand-algebra component mask (localIndependentMask), rather than the old
    % exact equality with the range rule. On this fixture the two differ, e.g.
    % at (n1,n2) = (2,15) the scaled trial contrast is about -.06 with D about
    % .86: in range, yet not a variance. TestDodGridAdmissibility carries the
    % frozen anchors and reasons.

    methods (Test)

        function testInadmissibleCellsAreMaskedNotClipped(testCase)
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = localFixture();
            pd.rel.out.trl_varcov{1}(:, 1, 1) = 1;
            pd.rel.out.trl_varcov{1}(:, 2, 2) = 36;
            pd.rel.out.trl_varcov{1}(:, 1, 2) = 5.4;
            pd.rel.out.trl_varcov{1}(:, 2, 1) = 5.4;
            [imgs, out] = localHeatmapImages(testCase, pd);
            anyBad = false;
            anyComponentOnly = false;
            for k = 1:numel(imgs)
                z = imgs(k).CData;
                bad = ~isfinite(z) | z < 0 | z > 1;
                anyBad = anyBad || any(bad(:));
                testCase.verifyEqual(size(imgs(k).AlphaData), size(z), ...
                    'the heatmap must carry a per-cell mask (AlphaData) the size of its data');
                if isequal(size(imgs(k).AlphaData), size(z))
                    alpha = double(imgs(k).AlphaData);
                    testCase.verifyTrue(all(alpha(bad) == 0), ...
                        'every inadmissible (NaN) or out-of-range cell must be masked');
                    expected = localIndependentMask(pd, imgs(k), out);
                    testCase.verifyEqual(alpha, double(expected), ...
                        ['the mask must equal the independent component-level rule ' ...
                        '(every recorded variance >= 0 and finite, denominator > 0, value in [0,1]; ' ...
                        'zero-count cells visible)']);
                    anyComponentOnly = anyComponentOnly || any(~expected(:) & ~bad(:));
                end
                ax = ancestor(imgs(k), 'axes');
                testCase.verifyEqual(ax.CLim, [0 1], ...
                    'the color axis must span [0,1] without the data being clipped to it');
            end
            testCase.verifyTrue(anyBad, ['the fixture must produce at least one ' ...
                'inadmissible or out-of-range cell, or this test is vacuous']);
            testCase.verifyTrue(anyComponentOnly, ['the fixture must produce at least one ' ...
                'in-range cell with a negative component, or the component rule is untested here']);
        end

        function testHealthyFixtureIsFullyVisible(testCase)
            % positive control: the unmodified factory fixture has every cell
            % finite in [0,1] and nothing masked
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = localFixture();
            imgs = localHeatmapImages(testCase, pd);
            for k = 1:numel(imgs)
                z = imgs(k).CData;
                testCase.verifyTrue(all(isfinite(z(:)) & z(:) >= 0 & z(:) <= 1), ...
                    'the healthy fixture must have no inadmissible or out-of-range cell');
                testCase.verifyEqual(size(imgs(k).AlphaData), size(z));
                testCase.verifyTrue(all(imgs(k).AlphaData(:) == 1), ...
                    'nothing may be masked on the healthy fixture');
            end
        end

        function testMaskBackgroundIsGreyNotWhite(testCase)
            % Live GUI look, Item 9 (2026-09-12): the masked cells rendered
            % WHITE. imagesc runs newplot, which resets an axes whose NextPlot
            % is 'replace' to its default properties and so discards a Color
            % given when the axes was created. The mask must show through as a
            % neutral grey that is neither white nor the figure background, or
            % a masked cell cannot be told from an empty part of the window.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            [imgs, out] = localHeatmapImages(testCase, localInflatedFixture());
            for k = 1:numel(imgs)
                ax = ancestor(imgs(k), 'axes');
                c = ax.Color;
                testCase.assertTrue(isnumeric(c) && numel(c) == 3, ...
                    'the axes color must be an RGB triplet');
                testCase.verifyTrue(all(abs(c - c(1)) < 1e-9), ...
                    'the mask background must be a neutral grey');
                testCase.verifyLessThan(c(1), 0.95, ...
                    'the mask background must not be white');
                testCase.verifyLessThan(c(1), out.Color(1) - 0.05, ...
                    'the mask background must be darker than the figure background');
            end
        end

        function testMaskCaptionClearsAxisLabels(testCase)
            % Live GUI look, Item 9 (2026-09-12): the third caption sat at
            % y = 28 and covered the "Trials <event>" x-axis labels, whose boxes
            % span roughly y = 40 to 56 in figure pixels. The caption must not
            % intersect either label's box or the two existing captions, and it
            % must lie inside the figure. uicontrols paint above axes, so an
            % intersection is an occlusion, not a harmless overlap.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            [imgs, out] = localHeatmapImages(testCase, localInflatedFixture());
            caps = findall(out, 'Type', 'uicontrol', 'Style', 'text');
            isGrey = arrayfun(@(t) startsWith(localString(t.String), 'Grey cells'), caps);
            testCase.assertEqual(nnz(isGrey), 1, 'expected exactly one "Grey cells" caption');
            greyBox = localBox(caps(isGrey).Position);
            testCase.verifyGreaterThanOrEqual(greyBox(2), 0, ...
                'the caption must lie inside the figure');
            others = caps(~isGrey);
            for k = 1:numel(others)
                testCase.verifyFalse(localOverlaps(greyBox, localBox(others(k).Position)), ...
                    sprintf('the caption must not overlap the text "%s"', ...
                    localString(others(k).String)));
            end
            for k = 1:numel(imgs)
                ax = ancestor(imgs(k), 'axes');
                testCase.verifyFalse(localOverlaps(greyBox, localXLabelBox(ax)), ...
                    sprintf('the caption must not cover the x-axis label "%s"', ...
                    localString(ax.XLabel.String)));
            end
        end

    end
end

function pd = localFixture()
pd = PsyRATTestDataFactory.makeDoDRelData();
if ~isfield(pd, 'ver'); pd.ver = '0-test'; end
if ~isfield(pd, 'proc') || ~isfield(pd.proc, 'measheader')
    pd.proc.measheader = 'meas';
end
end

function pd = localInflatedFixture()
% group 1's event-2 trial covariance inflated (variance 36, covariance 5.4
% with event 1; still positive semidefinite): the first heatmap's grid gets
% inadmissible and out-of-range cells where event-1 trials are few and
% event-2 trials are many
pd = localFixture();
pd.rel.out.trl_varcov{1}(:, 1, 1) = 1;
pd.rel.out.trl_varcov{1}(:, 2, 2) = 36;
pd.rel.out.trl_varcov{1}(:, 1, 2) = 5.4;
pd.rel.out.trl_varcov{1}(:, 2, 1) = 5.4;
end

function s = localString(s)
% a uicontrol or label String as one char row (cellstr rows joined by a space)
if iscell(s); s = strjoin(s, ' '); end
s = char(string(s));
end

function b = localBox(pos)
% [x0 y0 x1 y1] corners from a [left bottom width height] position
b = [pos(1), pos(2), pos(1) + pos(3), pos(2) + pos(4)];
end

function tf = localOverlaps(a, b)
% true when two [x0 y0 x1 y1] boxes share any interior area
tf = a(1) < b(3) && b(1) < a(3) && a(2) < b(4) && b(2) < a(4);
end

function b = localXLabelBox(ax)
% The x-axis label's box in figure pixels: a text object's Extent is measured
% from the axes' lower-left corner, so add the axes position (both in pixels).
au = ax.Units; ax.Units = 'pixels';
restoreAxes = onCleanup(@() set(ax, 'Units', au)); %#ok<NASGU>
lab = ax.XLabel; lu = lab.Units; lab.Units = 'pixels';
restoreLabel = onCleanup(@() set(lab, 'Units', lu)); %#ok<NASGU>
e = lab.Extent;
b = [ax.Position(1) + e(1), ax.Position(2) + e(2), ...
    ax.Position(1) + e(1) + e(3), ax.Position(2) + e(2) + e(4)];
end

function expected = localIndependentMask(pd, img, out)
% Independent component-level mask for one heatmap of the outputs window,
% from hand-written contrast algebra on the posterior MEAN blocks of the
% fixture's group 1 (the screen opens on group 1 with the recorded mapping
% A,B,C,D). The entered counts and the metric are read from the on-screen
% caption, the grid size from the image itself, and which slots the image
% varies from its x-axis label -- none of it from the figure UserData.
events = cellstr(string(pd.rel.events(:)));
bp = squeeze(mean(pd.rel.out.id_varcov{1}, 1));
bt = squeeze(mean(pd.rel.out.trl_varcov{1}, 1));
er = mean(exp(pd.rel.out.b_sigma{1}).^2, 1);
caps = findall(out, 'Type', 'uicontrol', 'Style', 'text');
obs = [];
metric = 1;
for k = 1:numel(caps)
    s = localString(caps(k).String);
    tok = regexp(s, 'Entered trials: \[(\d+) (\d+) (\d+) (\d+)\]', 'tokens', 'once');
    if ~isempty(tok)
        obs = cellfun(@str2double, tok);
        if contains(s, 'Heatmap metric: Generalizability'); metric = 2; end
    end
end
assert(numel(obs) == 4, 'the entered-trials caption was not found');
xl = localString(ancestor(img, 'axes').XLabel.String);
first = endsWith(xl, events{1});
ng = size(img.CData, 1);
c = [1; -1; -1; 1];
expected = false(ng, ng);
for xi = 1:ng
    for yi = 1:ng
        if first
            n = [xi - 1, yi - 1, obs(3), obs(4)];
        else
            n = [obs(1), obs(2), xi - 1, yi - 1];
        end
        if any(n <= 0)
            expected(yi, xi) = true;      % zero-count convention: visible 0
            continue;
        end
        u = c' * bp * c;
        s = zeros(4);
        for i = 1:4
            for j = 1:4
                if i == j
                    s(i, j) = bt(i, j) / n(i);
                else
                    s(i, j) = bt(i, j) / (2 / (1 / n(i) + 1 / n(j)));
                end
            end
        end
        bts = c' * s * c;
        ers = sum((c .^ 2) .* (er(:) ./ n(:)));
        btu = c' * bt * c;
        if metric == 1
            den = u + bts + ers;
        else
            den = u + ers;
        end
        v = u / den;
        expected(yi, xi) = all(isfinite([u bts ers btu])) && u >= 0 && bts >= 0 && ...
            ers >= 0 && btu >= 0 && isfinite(den) && den > 0 && isfinite(v) && v >= 0 && v <= 1;
    end
end
end

function [imgs, out] = localHeatmapImages(testCase, pd)
% Open the setup screen, fire its own Generate DoD Outputs callback, and
% return the two heatmap image objects and the outputs window.
psyrat_startview_dodiff('psyrat_prefs', psyrat_defaults, 'psyrat_data', pd);
setup = findall(groot, 'Type', 'figure', 'Name', 'DoD Reliability Setup');
testCase.assertNotEmpty(setup, 'the DoD setup screen did not open');
btns = findall(setup(1), 'Type', 'uicontrol', 'Style', 'pushbutton');
gen = [];
for k = 1:numel(btns)
    s = get(btns(k), 'String');
    if iscell(s); s = strjoin(s, ' '); end
    if contains(s, 'Generate DoD'); gen = btns(k); end
end
testCase.assertNotEmpty(gen, 'Generate DoD Outputs button not found');
cb = get(gen, 'Callback');
feval(cb{1}, gen, [], cb{2:end});
out = findall(groot, 'Type', 'figure', 'Name', 'DoD Reliability Outputs');
testCase.assertNotEmpty(out, 'the DoD outputs window did not open');
out = out(1);
imgs = findall(out, 'Type', 'image');
testCase.assertEqual(numel(imgs), 2, 'expected the two heatmap images');
end
