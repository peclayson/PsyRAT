classdef TestDynrelPopupBounds < PsyRATTestBase
    %G106 (P3 cosmetic): on the two-facet dynamic viewer the "Reliability
    %coefficient" popup overran its window. psyrat_startview_dynrel places
    %the right-column controls at rcol = 5/9 of the 560 px window (311.1 px)
    %and gave this popup 260 px, so its right edge sat at 571.1 and the
    %dropdown arrow was drawn on or past the window edge (seen on screen on
    %a Gaussian analysis-14 viewer and a gamma analysis-27 viewer). Every
    %other popup ends at 541.1 or less.
    %
    %Pinned for each two-facet arm (non-subject 14; subject-level 27; and
    %the gamma 27 Restricted layout, which adds the on-screen note):
    % 1. the reliability-coefficient popup ends at least MinRightMargin px
    %    inside the window (RED on the unfixed viewer: 571.1 in 560);
    % 2. every control lies horizontally inside the window, and the
    %    interactive ones keep the same margin (RED on the unfixed viewer,
    %    through that popup alone);
    % 3. the fix keeps the controls: four popups and three edit fields; the
    %    coefficient, reliability-coefficient and occasion option lists in
    %    order with their defaults, which buildSummary maps straight to
    %    gcoeff, psyrat_rel_trt's reltype and the nocc option; each
    %    two-facet selector still on its label's row, right of the label
    %    (control on the unfixed viewer).
    %Pinned on the two Gaussian arms:
    % 4. the resized popup still drives the surface: reltype 1 and 3 draw
    %    different curves, and returning to 3 redraws the first curve
    %    exactly (control).
    %Counterfactual:
    % 5. a one-facet viewer has no reliability popup and keeps its 560 px
    %    width, the native width of the manual's viewer_dynrel.png (an
    %    analysis-11 capture), so the fix stays on the two-facet arms
    %    (control). A deliberate change there means recapturing that image.
    %
    %LIMITS. This is a structural read-back, not a look. Arithmetic passes at
    %equality while the rendered arrow can still touch the frame, which is
    %why MinRightMargin exists, and nothing here measures whether the
    %longest option is shown whole. An exportapp read of the fixed windows
    %is still owed. No restriction wording is pinned; the note's Tag only
    %identifies the arm. The viewers build from PsyRATTestDataFactory
    %fixtures: no CmdStan, no fit.
    %
    %RED and control are the expectations on the unfixed viewer; nothing
    %here was executed by the author.

    properties (Constant)
        ViewerName = 'Dynamic Reliability: View Results';
        %slack beyond bare containment for popups, edit fields and buttons.
        %Every right-column control other than the G106 popup keeps 18.9 px
        %or more on the unfixed viewer.
        MinRightMargin = 10;
        CoeffOptions = {'Generalizability G(z)', 'Dependability D(z)'};
        RelOptions = {'Equivalence (occasion fixed)', 'Stability (trial fixed)', ...
            'Equivalence + stability (both random)'};
        OccOptions = {'1 occasion', 'Observed occasions'};
        %native width of the manual's one-facet viewer capture
        OneFacetWidth = 560;
    end

    properties (TestParameter)
        %{factory method, whether the arm renders the Restricted note}
        TwoFacetArm = struct( ...
            'nonsubject14', {{'makeDynrelTrtSummaryData', false}}, ...
            'subject27', {{'makeDynrelSubjNondiffTrtSummaryData', false}}, ...
            'gammaRestricted27', {{'makeGammaDynrelSubjNondiffTrtSummaryData', true}});
        GaussianTwoFacetArm = struct( ...
            'nonsubject14', 'makeDynrelTrtSummaryData', ...
            'subject27', 'makeDynrelSubjNondiffTrtSummaryData');
    end

    methods (Test)

        function testReliabilityPopupEndsInsideWindow(testCase, TwoFacetArm)
            fig = testCase.openViewer(TwoFacetArm{1});
            pop = testCase.reltypePopup(fig);
            testCase.assertEqual(pop.Units, 'pixels', 'premise: popup geometry in pixels');
            pos = pop.Position;
            figw = fig.Position(3);
            testCase.verifyGreaterThanOrEqual(pos(1), 0, ...
                'the Reliability coefficient popup starts left of the window');
            testCase.verifyLessThanOrEqual(pos(1) + pos(3), ...
                figw - testCase.MinRightMargin, sprintf( ...
                ['G106: the Reliability coefficient popup must end at least ' ...
                '%d px inside the %.1f px window; its right edge is at %.1f'], ...
                testCase.MinRightMargin, figw, pos(1) + pos(3)));
        end

        function testEveryControlStaysInsideWindow(testCase, TwoFacetArm)
            fig = testCase.openViewer(TwoFacetArm{1});
            figw = fig.Position(3);
            u = findall(fig, 'Type', 'uicontrol');
            testCase.assertNotEmpty(u, 'the viewer holds no uicontrols');
            for h = u(:)'
                testCase.assertEqual(h.Units, 'pixels', 'premise: control geometry in pixels');
                pos = h.Position;
                name = localLabel(h);
                testCase.verifyGreaterThanOrEqual(pos(1), 0, sprintf( ...
                    '%s "%s" starts left of the window (x = %.1f)', h.Style, name, pos(1)));
                %text rows may span the full width (the centered Dataset and
                %Dimensions lines do, at equality); interactive controls keep
                %the margin
                limit = figw;
                if ismember(h.Style, {'popupmenu', 'edit', 'pushbutton'})
                    limit = figw - testCase.MinRightMargin;
                end
                testCase.verifyLessThanOrEqual(pos(1) + pos(3), limit, sprintf( ...
                    '%s "%s" ends at %.1f, past %.1f in a %.1f px window', ...
                    h.Style, name, pos(1) + pos(3), limit, figw));
            end
        end

        function testSelectorsKeepTheirOptionsAndDefaults(testCase, TwoFacetArm)
            fig = testCase.openViewer(TwoFacetArm{1});
            testCase.verifyNumElements(findall(fig, 'Style', 'popupmenu'), 4, ...
                'coefficient, estimand, reliability coefficient and occasions popups');
            testCase.verifyNumElements(findall(fig, 'Style', 'edit'), 3, ...
                'trials, credible interval and reference cutoff fields');

            coeff = testCase.onePopupWithOption(fig, testCase.CoeffOptions{1});
            testCase.verifyEqual(localOptions(coeff), testCase.CoeffOptions);
            testCase.verifyEqual(coeff.Value, 1, 'the coefficient default is G(z)');

            rel = testCase.reltypePopup(fig);
            testCase.verifyEqual(localOptions(rel), testCase.RelOptions, ...
                'buildSummary passes the popup Value as reltype, so the order is load-bearing');
            testCase.verifyEqual(rel.Value, 3, ...
                'the reliability-coefficient default is reltype 3 (both facets random)');
            testCase.verifyLabelLeftOf(fig, rel, 'Reliability coefficient:');

            occ = testCase.onePopupWithOption(fig, testCase.OccOptions{1});
            testCase.verifyEqual(localOptions(occ), testCase.OccOptions);
            testCase.verifyEqual(occ.Value, 1, 'the occasion default is 1 occasion');
            testCase.verifyLabelLeftOf(fig, occ, 'Occasions n'':');

            %arm identity only (the wording is not pinned here)
            testCase.verifyNumElements( ...
                findall(fig, 'Tag', 'psyrat_dynrel_restriction_note'), ...
                double(TwoFacetArm{2}), ...
                'the Restricted note must appear exactly on the gamma two-facet arm');
        end

        function testReliabilityPopupStillDrivesTheSurface(testCase, GaussianTwoFacetArm)
            fig = testCase.openViewer(GaussianTwoFacetArm);
            rel = testCase.reltypePopup(fig);
            testCase.assertEqual(rel.Value, 3, 'premise: the viewer opens on reltype 3');
            y3 = testCase.plottedCurve(fig);
            rel.Value = 1;
            y1 = testCase.plottedCurve(fig);
            rel.Value = 3;
            y3again = testCase.plottedCurve(fig);
            %positive control first: the same selection redraws the same
            %curve, so a difference below comes from the popup
            testCase.verifyEqual(y3again, y3, ...
                'reltype 3 -> 1 -> 3 must redraw the first curve exactly');
            testCase.verifyNotEqual(y1, y3, ...
                'the popup no longer reaches the surface: reltype 1 and 3 drew the same curve');
        end

        function testOneFacetViewerIsUnchanged(testCase)
            fig = testCase.openViewer('makeDynrelSubjNondiffSummaryData');
            testCase.verifyEmpty(testCase.popupsWithOption(fig, testCase.RelOptions{3}), ...
                'premise: a one-facet viewer carries no reliability-coefficient popup');
            testCase.verifyEqual(fig.Position(3), testCase.OneFacetWidth, ...
                ['the G106 fix is scoped to the two-facet arms: the one-facet ' ...
                'viewer keeps the 560 px width of the manual''s viewer_dynrel.png']);
        end

    end

    methods (Access = private)

        function fig = openViewer(testCase, factoryMethod)
            %build the viewer from a PsyRATTestDataFactory fixture
            make = str2func(['PsyRATTestDataFactory.' factoryMethod]);
            pd = make();
            close(findall(0, 'Type', 'figure', 'Name', testCase.ViewerName));
            psyrat_startview_dynrel('psyrat_prefs', psyrat_defaults(), 'psyrat_data', pd);
            fig = findall(0, 'Type', 'figure', 'Name', testCase.ViewerName);
            testCase.addTeardown(@() localClose(fig));
            testCase.assertNumElements(fig, 1, 'expected one dynamic-reliability viewer');
            testCase.assertEqual(fig.Units, 'pixels', 'premise: window geometry in pixels');
        end

        function h = popupsWithOption(~, fig, option)
            %every popup menu whose option list contains OPTION
            pops = findall(fig, 'Style', 'popupmenu');
            keep = false(size(pops));
            for k = 1:numel(pops)
                keep(k) = any(strcmp(localOptions(pops(k)), option));
            end
            h = pops(keep);
        end

        function h = onePopupWithOption(testCase, fig, option)
            h = testCase.popupsWithOption(fig, option);
            testCase.assertNumElements(h, 1, sprintf( ...
                'expected one popup offering "%s"', option));
        end

        function h = reltypePopup(testCase, fig)
            %found by its longest option, the one G106 is about
            h = testCase.onePopupWithOption(fig, testCase.RelOptions{3});
        end

        function verifyLabelLeftOf(testCase, fig, ctrl, labelText)
            %the selector stays on its label's row, clear of the label box
            lab = findall(fig, 'Style', 'text', 'String', labelText);
            testCase.assertNumElements(lab, 1, sprintf('expected one "%s" label', labelText));
            testCase.verifyEqual(lab.Position(2), ctrl.Position(2), sprintf( ...
                'the "%s" selector left its label''s row', labelText));
            testCase.verifyLessThanOrEqual(lab.Position(1) + lab.Position(3), ...
                ctrl.Position(1), sprintf( ...
                'the "%s" selector now overlaps its label box', labelText));
        end

        function y = plottedCurve(testCase, fig)
            %press Plot reliability surface through its own callback and
            %return the YData of the one curve it draws
            btn = findall(fig, 'Style', 'pushbutton', 'String', 'Plot reliability surface');
            testCase.assertNumElements(btn, 1, 'expected one Plot reliability surface button');
            before = findall(0, 'Type', 'figure', 'Tag', 'psyrat_output');
            cb = btn.Callback;
            cb(btn, []);
            out = setdiff(findall(0, 'Type', 'figure', 'Tag', 'psyrat_output'), before);
            testCase.addTeardown(@() localClose(out));
            testCase.assertNumElements(out, 1, 'the plot button must open one figure');
            curve = findobj(out, 'Type', 'line');
            testCase.assertNumElements(curve, 1, 'expected one reliability curve');
            y = curve.YData;
        end

    end
end

function c = localOptions(h)
%a popup's option list as a 1-by-n cellstr, whatever form String takes
c = reshape(cellstr(h.String), 1, []);
end

function s = localLabel(h)
%a control's String as one line, for failure messages
s = h.String;
if iscell(s)
    s = strjoin(s, ' ');
elseif size(s, 1) > 1
    s = strjoin(cellstr(s), ' ');
end
s = char(s);
end

function localClose(h)
%close whichever of H is still open
h = h(isgraphics(h));
if ~isempty(h)
    close(h, 'force');
end
end
