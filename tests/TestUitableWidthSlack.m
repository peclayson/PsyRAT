classdef TestUitableWidthSlack < PsyRATTestBase
    %Sweep follow-up F5 (ruled 2026-09-16): five all-numeric ColumnWidth
    %lists met or exceeded their uitable's width (psyrat_variancet's three
    %windows, psyrat_dod_observedt, psyrat_formula_guide), and the two
    %overall-table lists (psyrat_depoverallt, psyrat_trt_reloverallt) had
    %eight entries for nine columns, so the last column fell back to a default
    %width and clipped at 940 px (seen live, GUI look Item 11 E). Every one
    %of these windows is resize-off, so the user cannot widen around it.
    %
    %THE RULE (findings ledger, the owner-confirmed 2026-07-31 note
    %"leaves 20 px of slack" under the trial-cutoff table finding, and
    %psyrat_criterionfigures' in-source precedent): a rendered table needs its
    %column widths plus 20 px of grid slack inside the table width, and a
    %window whose row count is data-driven also needs the 18 px a vertical
    %scrollbar takes; a check that passes at equality clips on screen. This
    %class builds each window with gui = 1 and reads the BUILT uitable handle.
    %It cannot see a scrollbar; the look is live GUI Item 12.
    %
    %G104 (2026-09-17, live GUI Item 12's residual): header TEXT must fit
    %its column too. The classic uitable draws the header bold at about 9 pt
    %(the handle reports 8 pt normal) and truncates it with an ellipsis once
    %the column is narrower than the text plus a 16-17 px cell inset; a '|'
    %in a ColumnName is a line break and each line is measured alone. The
    %proxy is a hidden bold 9 pt text uicontrol's Extent, which is at most
    %3.3 px narrower and at most 7 px wider than the rendered text over 70
    %measured header strings, 65 of them reachable from the nine windows
    %(R2025b Update 5, macOS, exportapp at 1x; the inset from a 1 px grid on
    %41 headers). HeaderPad = 17 px inset + 4 px proxy error + 3 px headroom.
    %A failure on another platform is a re-calibration of these two
    %constants, not a regression. A window capped at the display width
    %(psyrat_fitfigwidth) is meant to scroll, so the column-sum rule is an
    %assumption there, not a failure.

    properties (Access = private)
        Slack = 20;
        ScrollbarAllowance = 18;   % a vertical scrollbar takes width from the columns on a data-driven table
    end

    properties (Constant)
        HeaderPad = 24;        % px a header needs beyond its bold 9 pt extent (G104 calibration, class comment)
        HeaderFontSize = 9;    % the classic uitable header renders bold at about this size, not the reported 8 pt
    end

    methods (TestMethodSetup)
        function closeOutputs(~)
            delete(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output'));
        end
    end

    methods (TestMethodTeardown)
        function closeOutputsAgain(~)
            delete(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output'));
            delete(findall(groot, 'Type', 'figure', 'Name', 'Formula Guide (Tables 2, 3, and 6)'));
        end
    end

    methods (Test)

        function testSingleSessionCutoffLayoutsFit(testCase)
            localAssertCutoffLayouts(testCase, @psyrat_depcutofft, localIcData());
        end

        function testTestRetestCutoffLayoutsFit(testCase)
            localAssertCutoffLayouts(testCase, @psyrat_trt_relcutofft, localTrtData());
        end

        function testVarianceTableIcLayoutFits(testCase)
            t = localTable(testCase, @() psyrat_variancet('psyrat_data', localIcData(), 'gui', 1));
            localAssertFits(testCase, t, 'psyrat_variancet ic (six columns)', testCase.ScrollbarAllowance);
        end

        function testVarianceTableFiveColumnLayoutFits(testCase)
            %the trt layout has five columns and already fit; a regression pin
            t = localTable(testCase, @() psyrat_variancet('psyrat_data', localTrtData(), 'gui', 1));
            localAssertFits(testCase, t, 'psyrat_variancet trt (five columns)', testCase.ScrollbarAllowance);
        end

        function testAllStandardDeviationsWindowFits(testCase)
            fig = localFig(testCase, @() psyrat_variancet('psyrat_data', localTrtData(), 'gui', 1));
            win = localOpenChild(testCase, fig, 'psyrat_viewvarcomp');
            localAssertFits(testCase, findall(win, 'Type', 'uitable'), 'psyrat_viewvarcomp (trt all SDs)', testCase.ScrollbarAllowance);
        end

        function testDifferenceComponentsWindowFits(testCase)
            %the equality case: the sum used to equal the table width exactly
            fig = localFig(testCase, @() psyrat_variancet('psyrat_data', localIcDiffData(), 'gui', 1));
            win = localOpenChild(testCase, fig, 'psyrat_viewdiffcomp');
            localAssertFits(testCase, findall(win, 'Type', 'uitable'), 'psyrat_viewdiffcomp (ic_diff all SDs)', testCase.ScrollbarAllowance);
        end

        function testObservedDodTableFits(testCase)
            pd = localAddHeaderFields(PsyRATTestDataFactory.makeDoDRelData());
            map = localDodMap(pd.rel);
            t = localTable(testCase, @() psyrat_dod_observedt('psyrat_data', pd, 'groupidx', 1, 'map', map, 'gui', 1));
            localAssertFits(testCase, t, 'psyrat_dod_observedt');
        end

        function testFormulaGuideTableFits(testCase)
            %28 rows in a 365 px table: the vertical scrollbar is always there
            t = localTable(testCase, @() psyrat_formula_guide());
            localAssertFits(testCase, t, 'psyrat_formula_guide', testCase.ScrollbarAllowance);
        end

        function testOverallTableRetainedBranchFits(testCase)
            t = localTable(testCase, @() psyrat_depoverallt('psyrat_data', localIcData(), 'gui', 1));
            localAssertFits(testCase, t, 'psyrat_depoverallt (retained)', testCase.ScrollbarAllowance);
        end

        function testOverallTableFullSampleBranchFits(testCase)
            %the 350 px coefficient column of the full-sample reference rows
            pd = localIcDiffData();
            pd.relsummary.group(1).diffscore.summary_population = 'full_sample_reference';
            t = localTable(testCase, @() psyrat_depoverallt('psyrat_data', pd, 'gui', 1));
            d = get(t, 'Data');
            testCase.assertTrue(any(contains(string(d(:, 4)), 'full-sample reference')), ...
                'the fixture must take the full-sample branch (positive control)');
            localAssertFits(testCase, t, 'psyrat_depoverallt (full-sample reference)', testCase.ScrollbarAllowance);
        end

        function testTrtOverallTableRetainedBranchFits(testCase)
            t = localTable(testCase, @() psyrat_trt_reloverallt('psyrat_data', localTrtData(), 'gui', 1));
            localAssertFits(testCase, t, 'psyrat_trt_reloverallt (retained)', testCase.ScrollbarAllowance);
        end

        function testNewReliabilityWindowWidthIsPinnedAtSource(testCase)
            %the trt overall table's recalculation window (psyrat_shownewrel):
            %three columns of 151 + 120 + 100 = 371 px (G104 widened the
            %coefficient column from 110 for 'Generalizability') in a window
            %that was 415 px wide before F5 (365 px of table, 4 px of slack).
            %The window recomputes from a fitted rel.out that the summary
            %fixtures do not carry, so it cannot be built here; its width is
            %pinned at source against the same arithmetic (371 + 2 x 25
            %margins + 20 slack = 441 <= 450).
            seg = localNewRelSegment(testCase);
            testCase.assertEqual(numel(strfind(seg, "set(t,'ColumnWidth',{151 120 100});")), 1, ...
                'the recalculation table must still carry its three 371 px columns (positive control)');
            testCase.verifyEqual(numel(strfind(seg, 'figwidth = 450;')), 1, ...
                'the recalculation window must be 450 px wide (371 px of columns + margins + 20 px slack)');
            testCase.verifyEqual(numel(strfind(seg, 'figwidth = 415;')), 0, ...
                'the 415 px window left 4 px of slack, below the 20 px rule');
        end

        function testTrtOverallTableFullSampleBranchFits(testCase)
            pd = localTrtData();
            for g = 1:numel(pd.relsummary.group)
                for e = 1:numel(pd.relsummary.group(g).event)
                    pd.relsummary.group(g).event(e).summary_population = 'full_sample_reference';
                end
            end
            t = localTable(testCase, @() psyrat_trt_reloverallt('psyrat_data', pd, 'gui', 1));
            d = get(t, 'Data');
            testCase.assertTrue(any(contains(string(d(:, 4)), 'full-sample reference')), ...
                'the fixture must take the full-sample branch (positive control)');
            localAssertFits(testCase, t, 'psyrat_trt_reloverallt (full-sample reference)', testCase.ScrollbarAllowance);
        end

        function testHeaderCheckSeesATruncatedHeader(testCase)
            %positive control for the header term: 'SEM (dependability)' in a
            %130 px column is exactly what Item 12 saw truncated on screen,
            %and the same words on two lines fit that column
            [t, cleanup] = localScratchTable(); %#ok<ASGLU>
            [need1, ext1] = localHeaderNeed(t, 'SEM (dependability)');
            [need2, ext2] = localHeaderNeed(t, 'SEM|(dependability)');
            testCase.verifyGreaterThan(need1, 130, sprintf( ...
                'the header term must flag "SEM (dependability)" in 130 px (extent %.1f, need %d)', ext1, need1));
            testCase.verifyLessThanOrEqual(need2, 130, sprintf( ...
                'the header term must accept "SEM|(dependability)" in 130 px (longest line %.1f, need %d)', ext2, need2));
        end

        function testVarianceTableIcLayoutFitsGeneralizability(testCase)
            %the relative-error headers are the longer variant of columns 5-6
            pd = localIcData();
            pd.relsummary.gcoeff = 2; pd.relsummary.gcoeff_name = 'gen';
            t = localTable(testCase, @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
            cn = cellstr(get(t, 'ColumnName'));
            testCase.assertTrue(any(contains(cn, 'generalizability')), ...
                'the fixture must render the generalizability headers (positive control)');
            localAssertFits(testCase, t, 'psyrat_variancet ic (generalizability)', testCase.ScrollbarAllowance);
        end

        function testVarianceTableFiveColumnLayoutFitsGeneralizability(testCase)
            pd = localTrtData();
            pd.relsummary.gcoeff = 2; pd.relsummary.gcoeff_name = 'gen';
            t = localTable(testCase, @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
            cn = cellstr(get(t, 'ColumnName'));
            testCase.assertTrue(any(contains(cn, 'generalizability')), ...
                'the fixture must render the generalizability headers (positive control)');
            localAssertFits(testCase, t, 'psyrat_variancet trt (generalizability)', testCase.ScrollbarAllowance);
        end

        function testSserrVarianceTableFits(testCase)
            %the subject-level layout carries the longest sources-of-variance
            %headers ('Group SEM (generalizability)') in the same five-column
            %widths as the trt and difference layouts
            pd = localSserrData(testCase);
            t = localTable(testCase, @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
            cn = cellstr(get(t, 'ColumnName'));
            testCase.assertTrue(any(contains(cn, 'Group SEM')), ...
                'the subject-level headers must be on the table (positive control)');
            localAssertFits(testCase, t, 'psyrat_variancet ic_sserrvar (generalizability)', testCase.ScrollbarAllowance);
        end

        function testSserrOverallTableFits(testCase)
            %'Subject Level Generalizability' is the longest coefficient header
            pd = localSserrData(testCase);
            t = localTable(testCase, @() psyrat_depoverallt('psyrat_data', pd, 'gui', 1));
            cn = cellstr(get(t, 'ColumnName'));
            testCase.assertTrue(any(contains(cn, 'Subject Level')), ...
                'the subject-level coefficient header must be on the table (positive control)');
            localAssertFits(testCase, t, 'psyrat_depoverallt ic_sserrvar (generalizability)', testCase.ScrollbarAllowance);
        end

        function testVarianceTableIcLayoutFitsSplits(testCase)
            %the parallel-splits run relabels the unit ('Split Std Dev')
            pd = localIcData(); pd.rel.splits = 2;
            t = localTable(testCase, @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
            cn = cellstr(get(t, 'ColumnName'));
            testCase.assertTrue(any(contains(cn, 'Split')), 'the fixture must render the split labels (positive control)');
            localAssertFits(testCase, t, 'psyrat_variancet ic (splits)', testCase.ScrollbarAllowance);
        end

        function testAllStandardDeviationsWindowFitsSplits(testCase)
            %'Between Split', 'Person x Split', 'Session x Split'
            pd = localTrtData(); pd.rel.splits = 2;
            fig = localFig(testCase, @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
            win = localOpenChild(testCase, fig, 'psyrat_viewvarcomp');
            t = findall(win, 'Type', 'uitable');
            cn = cellstr(get(t, 'ColumnName'));
            testCase.assertTrue(any(contains(cn, 'Split')), 'the fixture must render the split labels (positive control)');
            localAssertFits(testCase, t, 'psyrat_viewvarcomp (trt all SDs, splits)', testCase.ScrollbarAllowance);
        end

        function testOverallTableRetainedBranchFitsSplits(testCase)
            %'Mean # of Splits', 'Std Dev of Splits' and the other split forms
            pd = localIcData(); pd.rel.splits = 2;
            t = localTable(testCase, @() psyrat_depoverallt('psyrat_data', pd, 'gui', 1));
            cn = cellstr(get(t, 'ColumnName'));
            testCase.assertTrue(any(contains(cn, 'Splits')), 'the fixture must render the split labels (positive control)');
            localAssertFits(testCase, t, 'psyrat_depoverallt (retained, splits)', testCase.ScrollbarAllowance);
        end

        function testTrtOverallTableRetainedBranchFitsGeneralizability(testCase)
            %'Generalizability' is the longer coefficient header of the trt table
            pd = localTrtData();
            pd.relsummary.gcoeff = 2; pd.relsummary.gcoeff_name = 'gen';
            t = localTable(testCase, @() psyrat_trt_reloverallt('psyrat_data', pd, 'gui', 1));
            cn = cellstr(get(t, 'ColumnName'));
            testCase.assertTrue(any(strcmp(cn, 'Generalizability')), 'the fixture must render the generalizability header (positive control)');
            localAssertFits(testCase, t, 'psyrat_trt_reloverallt (retained, generalizability)', testCase.ScrollbarAllowance);
        end

        function testNewReliabilityWindowHeadersFitAtSource(testCase)
            %the recalculation window cannot be built (see the width pin
            %above); its coefficient column shows 'Dependability' or
            %'Generalizability', so the longer one is measured against the
            %width literal in psyrat_shownewrel, and '# of Trials' against its
            seg = localNewRelSegment(testCase);
            tok = regexp(seg, "set\(t,'ColumnWidth',\{151 (\d+) 100\}\);", 'tokens', 'once');
            testCase.assertNotEmpty(tok, 'the recalculation table must still carry {151 <coefficient> 100} widths');
            coeffwidth = str2double(tok{1});
            [t, cleanup] = localScratchTable(); %#ok<ASGLU>
            [need, ext] = localHeaderNeed(t, 'Generalizability');
            testCase.verifyLessThanOrEqual(need, coeffwidth, sprintf( ...
                'the recalculation window''s coefficient column must hold "Generalizability" (%d px: bold 9 pt extent %.1f + %d) but is %d px', ...
                need, ext, TestUitableWidthSlack.HeaderPad, coeffwidth));
            [need, ext] = localHeaderNeed(t, '# of Trials');
            testCase.verifyLessThanOrEqual(need, 100, sprintf( ...
                'the recalculation window''s trial column must hold "# of Trials" (%d px: extent %.1f + %d) in 100 px', ...
                need, ext, TestUitableWidthSlack.HeaderPad));
        end

        function testHeaderMapInsertsOneBreak(testCase)
            %G104: the on-screen line break for the headers Item 12 saw
            %truncated; one break each, both lines non-empty, keys unique
            in = {'SEM (dependability)', 'ICC (generalizability)', 'Group SEM (dependability)', ...
                'Within Person (Error)', 'Error Covariance Std Dev', 'Max # of Trials', ...
                'Std Dev of Splits', 'Subject Level Generalizability', 'n Excluded'};
            expected = {'SEM|(dependability)', 'ICC|(generalizability)', 'Group SEM|(dependability)', ...
                'Within Person|(Error)', 'Error Covariance|Std Dev', 'Max|# of Trials', ...
                'Std Dev of|Splits', 'Subject Level|Generalizability', 'n|Excluded'};
            testCase.verifyEqual(psyrat_uitable_headers(in), expected);
            map = psyrat_uitable_headers();
            testCase.assertTrue(iscell(map) && size(map, 2) == 2, 'the zero-argument call returns the two-column map');
            testCase.verifyEqual(numel(unique(map(:, 1))), size(map, 1), 'map keys must be unique');
            testCase.verifyFalse(any(contains(map(:, 1), '|')), 'map keys are the one-line display names');
            for k = 1:size(map, 1)
                parts = strsplit(map{k, 2}, '|');
                testCase.verifyEqual(numel(parts), 2, sprintf('"%s" must break into exactly two lines', map{k, 1}));
                testCase.verifyTrue(all(~cellfun(@isempty, strtrim(parts))), ...
                    sprintf('"%s" must not break into an empty line', map{k, 1}));
            end
        end

        function testHeaderMapLeavesUnmappedNamesAlone(testCase)
            in = {'Label', 'Between Std Dev', 'Trial Std Dev', 'ICC (n''=1)', 'Dependability', 'Table'};
            testCase.verifyEqual(psyrat_uitable_headers(in), in);
            s = psyrat_uitable_headers(["Label", "Max # of Trials"]);
            testCase.verifyClass(s, 'string');
            testCase.verifyEqual(s, ["Label", "Max|# of Trials"]);
        end

        function testHeaderMapNamesTheDodColumnsReadably(testCase)
            %the observed DoD table shows its CSV variable names; on screen
            %the counts break into two readable lines and the four mean-trial
            %columns become distinguishable ('Mean Trials' / 'ERP n')
            in = {'Label', 'n_Included', 'n_Excluded', 'Dependability', 'Generalizability', ...
                'ICC_Dep', 'ICC_Gen', 'MeanTrials_ERP_1', 'MeanTrials_ERP_2', 'MeanTrials_ERP_3', 'MeanTrials_ERP_4'};
            expected = {'Label', 'n|Included', 'n|Excluded', 'Dependability', 'Generalizability', ...
                'ICC_Dep', 'ICC_Gen', 'Mean Trials|ERP 1', 'Mean Trials|ERP 2', 'Mean Trials|ERP 3', 'Mean Trials|ERP 4'};
            testCase.verifyEqual(psyrat_uitable_headers(in), expected);
        end

        function testFitFigWidthCapsAtTheDisplay(testCase)
            %G104: a window wider than the display opens 40 px inside each
            %edge and its table scrolls; a narrower one keeps its width
            scr = get(groot, 'ScreenSize');
            testCase.assumeTrue(scr(3) > 280, 'needs a display wider than the helper''s 280 px floor');
            testCase.verifyEqual(psyrat_fitfigwidth(100), 100);
            testCase.verifyEqual(psyrat_fitfigwidth(scr(3) - 80), scr(3) - 80);
            testCase.verifyEqual(psyrat_fitfigwidth(100000), scr(3) - 80, ...
                'a window wider than the display must open 40 px inside each edge');
        end

    end
end


function localAssertCutoffLayouts(testCase, builder, pd)
%Exercise both header vocabularies and both coefficient widths through the
%real builders. Enough groups to require a vertical scrollbar, with table
%data checked against the same builder's headless return.
pd.rel.groups = arrayfun(@(n) sprintf('Group %02d', n), 1:20, 'UniformOutput', false);
pd.relsummary.group = repmat(pd.relsummary.group(1), 1, 20);
for splits = [1 2]
    pd.rel.splits = splits;
    L = psyrat_unitlabels(splits);
    for coeff = {'dep', 'gen'}
        pd.relsummary.gcoeff_name = coeff{1};
        for projection = [false true]
            for g = 1:numel(pd.relsummary.group)
                for e = 1:numel(pd.relsummary.group(g).event)
                    pd.relsummary.group(g).event(e).cutoff_extrapolated = projection;
                end
            end
            expected = builder('psyrat_data', pd, 'gui', 0);
            t = localTable(testCase, @() builder('psyrat_data', pd, 'gui', 1));
            testCase.verifyEqual(get(t, 'Data'), table2cell(expected));
            headers = cellstr(get(t, 'ColumnName'));
            testCase.verifyEqual(headers{2}, L.Cutoff);
            testCase.verifyEqual(all(contains(string(expected{:, 3}), ...
                'projection beyond observed data')), projection);
            localAssertFits(testCase, t, sprintf('%s %s %s projection=%d', ...
                func2str(builder), L.Cutoff, coeff{1}, projection), 18);
        end
    end
end
end


function localAssertFits(testCase, t, what, extra)
%The guard: every column has a numeric width, and the widths plus 20 px of
%slack (plus a scrollbar allowance where the caller says one is always
%present) fit inside the table's pixel width.
if nargin < 4; extra = 0; end
testCase.assertNumElements(t, 1, sprintf('%s: expected one uitable', what));
cw = get(t, 'ColumnWidth'); cn = get(t, 'ColumnName');
testCase.verifyEqual(numel(cw), numel(cn), sprintf( ...
    '%s: every column needs a width entry (%d widths for %d columns)', what, numel(cw), numel(cn)));
testCase.verifyTrue(all(cellfun(@isnumeric, cw)), sprintf( ...
    '%s: widths must be numeric so the sum is checkable (an ''auto'' column has no pinned width)', what));
if ~all(cellfun(@isnumeric, cw)); return; end
localAssertHeadersFit(testCase, t, what);
total = sum(cellfun(@(x) x, cw));
pos = get(t, 'Position');
width = pos(3);
if strcmp(get(t, 'Units'), 'normalized')
    parent = get(t, 'Parent');
    pp = get(parent, 'Position');
    width = width * pp(3);
end
%G104: a window capped at the display width (psyrat_fitfigwidth) is meant to
%be narrower than its columns and to scroll, so on such a display the sum
%rule is an assumption, not a failure
figpos = get(ancestor(t, 'figure'), 'Position');
if total + 20 + extra > width && figpos(3) >= psyrat_fitfigwidth(1e6)
    testCase.assumeFail(sprintf( ...
        '%s: the window is capped at the display width (%g px), so the column sum cannot be judged here', what, figpos(3)));
end
testCase.verifyLessThanOrEqual(total + 20 + extra, width, sprintf( ...
    '%s: %d px of columns + 20 px slack%s must fit the %g px table', what, total, ...
    localAllowanceNote(extra), width));
end


function s = localAllowanceNote(extra)
s = '';
if extra > 0; s = sprintf(' + %d px scrollbar', extra); end
end


function localAssertHeadersFit(testCase, t, what)
%G104: every header text fits its column (class comment for the calibration);
%a '|' is a line break, at most two lines, each measured on its own
cw = get(t, 'ColumnWidth'); cn = cellstr(get(t, 'ColumnName'));
for j = 1:min(numel(cn), numel(cw))   % a width-count mismatch is reported by the caller
    lines = strsplit(cn{j}, '|');
    testCase.verifyLessThanOrEqual(numel(lines), 2, sprintf( ...
        '%s: column %d "%s" may break into two lines at most', what, j, cn{j}));
    testCase.verifyTrue(all(~cellfun(@isempty, strtrim(lines))), sprintf( ...
        '%s: column %d "%s" must not carry an empty line', what, j, cn{j}));
    [need, ext] = localHeaderNeed(t, cn{j});
    testCase.verifyGreaterThanOrEqual(cw{j}, need, sprintf( ...
        '%s: column %d "%s" needs %d px (bold %g pt extent %.1f + %d px header inset) but has %d', ...
        what, j, cn{j}, need, TestUitableWidthSlack.HeaderFontSize, ext, TestUitableWidthSlack.HeaderPad, cw{j}));
end
end


function [need, ext] = localHeaderNeed(t, name)
%the pixels a header needs: the widest of its lines measured as a hidden
%bold 9 pt text uicontrol on the table's own figure, plus HeaderPad. The
%explicit empty CreateFcn and the 4000 px box keep psyrat_newfigure's
%autofit (psyrat_autofitbuttontext shrinks overflowing text controls) out
%of the measurement.
fig = ancestor(t, 'figure');
lines = strsplit(char(name), '|');
ext = 0;
for k = 1:numel(lines)
    c = uicontrol(fig, 'Style', 'text', 'String', lines{k}, ...
        'FontName', get(t, 'FontName'), 'FontUnits', 'points', ...
        'FontSize', TestUitableWidthSlack.HeaderFontSize, 'FontWeight', 'bold', ...
        'Units', 'pixels', 'Position', [1 1 4000 40], 'Visible', 'off', 'CreateFcn', '');
    e = get(c, 'Extent');
    delete(c);
    ext = max(ext, e(3));
end
need = ceil(ext) + TestUitableWidthSlack.HeaderPad;
end


function [t, cleanup] = localScratchTable()
%a uitable on a PsyRAT figure (psyrat_newfigure, the real windows' uicontrol
%defaults) for measuring header text on its own
fig = psyrat_newfigure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 400 200], 'Tag', 'psyrat_output');
t = uitable('Parent', fig, 'Position', [20 20 360 150]);
cleanup = onCleanup(@() delete(fig));
end


function seg = localNewRelSegment(testCase)
%the psyrat_shownewrel local function's source text (the recalculation
%window recomputes from a fitted rel.out the fixtures do not carry)
src = fileread(fullfile(testCase.projectRoot(), 'subroutines', 'tables', 'psyrat_trt_reloverallt.m'));
idef = regexp(src, '\nfunction psyrat_shownewrel\(', 'once');
testCase.assertNotEmpty(idef, 'psyrat_shownewrel must still be defined');
seg = src(idef + 1:end);
nextf = strfind(seg, [newline 'function ']);
if ~isempty(nextf); seg = seg(1:nextf(1)); end
end


function fig = localFig(testCase, builderThunk)
delete(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output'));
builderThunk();
fig = findobj(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
testCase.assertNotEmpty(fig, 'the builder created no output window');
fig = fig(1);
end


function t = localTable(testCase, builderThunk)
fig = localFig(testCase, builderThunk);
t = findall(fig, 'Type', 'uitable');
testCase.assertNumElements(t, 1, 'expected one uitable on the window');
end


function win = localOpenChild(testCase, fig, callbackName)
%fire the button whose Callback function is callbackName and return the
%window it opened
btns = findall(fig, 'Style', 'pushbutton');
before = findobj(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
for b = btns(:)'
    cb = b.Callback;
    if iscell(cb) && isa(cb{1}, 'function_handle') && strcmp(func2str(cb{1}), callbackName)
        feval(cb{1}, b, [], cb{2:end});
        after = findobj(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
        win = after(~ismember(after, before));
        testCase.assertNumElements(win, 1, sprintf('%s must open one window', callbackName));
        return;
    end
end
testCase.assertFail(sprintf('no button on the window calls %s', callbackName));
end


function map = localDodMap(rel)
events = cellstr(string(rel.events(:)));
mapnames = cellstr(string(rel.dod_map(:)));
[~, map] = ismember(mapnames, events);
map = map(:)';
end


function pd = localTrtData()
pd = localAddHeaderFields(PsyRATTestDataFactory.makeTRTSummaryData());
end


function pd = localSserrData(testCase)
%the subject-level (ic_sserrvar) fixture summarized the way the sserr tests
%do (TestSummaryAndPresentationFunctions; pure MATLAB, no CmdStan) at the
%relative-error coefficient, so psyrat_variancet renders the 'Group ...'
%headers and psyrat_depoverallt the 'Subject Level Generalizability' header
pd = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
[pd, relerr] = psyrat_relsummary('psyrat_data', pd, 'analysis', 'sing_sserr', ...
    'depcutoff', 0.01, 'meascutoff', 2, 'depcentmeas', 1, 'gcoeff', 2, 'CI', 0.95);
testCase.assertEqual(relerr.nogooddata, 0, 'the sserr summary must keep its data (positive control)');
pd.relsummary.gcoeff = 2;
pd.relsummary.gcoeff_name = 'gen';
pd = localAddHeaderFields(pd);
end


function pd = localIcData()
pd = localAddHeaderFields(PsyRATTestDataFactory.makeICSummaryDataWithTrial());
end


function pd = localIcDiffData()
%the factory's ic_diff summary (one ungrouped pair of events, a group-level
%diffscore) is what psyrat_depoverallt reads; the variance windows also read
%per-event variance structs (icc, betsd, trlsd, witsd, sem, each m/ll/ul) and
%the difference components (bt1/bt2 variances, the covariances, the ICC), so
%add them with deterministic values
pd = localAddHeaderFields(PsyRATTestDataFactory.makeICDiffSummaryData());
for g = 1:numel(pd.relsummary.group)
    for e = 1:numel(pd.relsummary.group(g).event)
        base = 1 + 0.1 * e;
        names = {'icc', 'betsd', 'trlsd', 'witsd', 'sem'};
        for k = 1:numel(names)
            v = base * k / 5;
            pd.relsummary.group(g).event(e).(names{k}) = struct('m', v, 'll', v * 0.8, 'ul', v * 1.2);
        end
    end
    ds = pd.relsummary.group(g).diffscore;
    names = {'bt1_var', 'bt2_var', 'bp_cov', 'bt_cov', 'icc'};
    vals = [2.1 1.7 0.9 0.4 0.55];
    for k = 1:numel(names)
        ds.([names{k} '_pt']) = vals(k);
        ds.([names{k} '_ll']) = vals(k) * 0.8;
        ds.([names{k} '_ul']) = vals(k) * 1.2;
    end
    ds.wp_cov = 0.3;
    pd.relsummary.group(g).diffscore = ds;
end
end


function pd = localAddHeaderFields(pd)
pd.ver = 'test';
pd.rel.filename = 'width_slack_fixture.psyrat';
pd.rel.nchains = 2;
pd.rel.niter = 600;
pd.rel.seed = 12345;
pd.rel.engine = 'cmdstan';
pd.rel.priors = psyrat_default_priors;
pd.rel.out.conv.converged = 1;
pd.rel.out.conv.ndivergent = 0;
end
