classdef TestRestrictedFigureDisclosure < PsyRATTestBase
    %TESTRESTRICTEDFIGUREDISCLOSURE The Restricted disclosure on SAVED FIGURES
    %of the three gamma two-facet designs (analyses 14, 25 and 27, both gamma
    %parameterizations): the optional LS27 disclosure extension (queue row 29).
    %The export headers, the viewer notes and the Analyze notice carried it
    %before; a saved PNG, PDF or FIG of the surface or subject plots did not
    %(four of the first six methods below failed on that code; the two
    %control methods passed, as intended).
    %
    %WHAT IS PINNED.
    %  - Each of the six arms, drawn through its real route, carries exactly
    %    one disclosure: an exportable textbox annotation with the shared
    %    saved-figure text (R.figure), and the same facts in the figure's
    %    UserData.
    %  - The existing gamma target caption stays a separate, unchanged object,
    %    and the disclosure covers no panel, label, caption or button.
    %  - The saved PNG shows the disclosure (its ink disappears when the
    %    disclosure is deleted), the PDF is written, and a reopened FIG keeps the
    %    disclosure and its metadata beside the existing carriers.
    %  - CONTROLS: no disclosure on the Gaussian arms, the gamma one-facet
    %    designs, the gamma occasion designs outside the ruling (5, 10, 29) or
    %    label cousins (15, 20, 6, 8); and a Restricted figure differs from an
    %    otherwise identical non-Restricted figure only by the disclosure.
    %  - LAYOUT PRESERVATION (review N3): against its silent twin, a
    %    Restricted figure keeps every axes, colorbar, existing caption and
    %    button at the same pixel box from its lower-left corner, keeps its
    %    width, and grows by exactly the band that holds the disclosure.
    %  - TWO DIMENSIONS, TWO STRATA (review N4): the filled-contour surface
    %    with one colorbar per panel carries one disclosure, clear of every
    %    colorbar and its tick labels, for both designs, both
    %    parameterizations and both coefficients.
    %
    %WHAT IS NOT PINNED. The wording itself (it is the row-27 text source, whose
    %content rules live in TestRestrictedDisclosureWording), rendered layout on
    %a real display (a native look; the geometry here is MATLAB's own layout
    %of invisible figures, not pixels on a screen), PDF text (no reader in base
    %MATLAB; the local artifact writer hands the files to a rendered read), and
    %the scripted route (psyrat_run draws no figure; the disclosure reaches a
    %script through psyrat_report's provenance lines).
    %
    %Fixtures and routes: LS27FigureFixture. No fit, no file outside tempdir,
    %no skip.

    properties (Access = private)
        ShadowDir
        TmpDir
    end

    methods (TestMethodSetup)
        function setupShadowAndTmp(testCase)
            %The surface figure's "Save figure (PNG/PDF)" button asks uiputfile
            %for a path; a shadow answers from a global (the TestExportCancelGuards
            %pattern), so the button's own callback does the writing.
            testCase.ShadowDir = tempname;
            mkdir(testCase.ShadowDir);
            fid = fopen(fullfile(testCase.ShadowDir, 'uiputfile.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function [filename,pathname,filterindex] = uiputfile(varargin)', ...
                '% TEST SHADOW: answers the save dialog from a global; 0 is Cancel.', ...
                'global PSYRAT_TEST_UIPUTFILE', ...
                'if isempty(PSYRAT_TEST_UIPUTFILE) || isequal(PSYRAT_TEST_UIPUTFILE, 0)', ...
                '    filename = 0; pathname = 0; filterindex = 0; return;', ...
                'end', ...
                '[p,n,e] = fileparts(PSYRAT_TEST_UIPUTFILE);', ...
                'filename = [n e]; pathname = [p filesep]; filterindex = 1;', ...
                'end');
            fclose(fid);
            addpath(testCase.ShadowDir, '-begin');
            rehash;
            testCase.assertTrue(startsWith(which('uiputfile'), testCase.ShadowDir), ...
                'the uiputfile shadow must be the resolved uiputfile');
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);
            testCase.addTeardown(@() localCleanup(testCase.ShadowDir, testCase.TmpDir));
        end
    end

    methods (Test)

        function testSixArmsCarryOneExportableDisclosureWithMetadata(testCase)
            % The six arms through their real routes: the dynamic viewer's
            % "Plot reliability surface" button for 14 and 27, and the
            % test-retest viewer's Generate route (psyrat_relfigures) for 25,
            % whose two subject plots are the only figures that design draws.
            expected = LS27FigureFixture.expectedText();
            for arm = {'g14ls', 'g14lognu', 'g27ls', 'g27lognu'}
                fig = testCase.surfaceFromViewer(arm{1});
                testCase.verifyDisclosure(fig, arm{1}, expected, LS27FigureFixture.label(arm{1}));
            end
            for fam = {'ls', 'lognu'}
                figs = testCase.trt25FromGenerate(fam{1});
                testCase.verifyNumElements(figs, 2, sprintf( ...
                    'analysis 25 %s: expected the two subject plots', fam{1}));
                for k = 1:numel(figs)
                    testCase.verifyDisclosure(figs(k), sprintf('g25%s %s', fam{1}, ...
                        get(figs(k), 'Name')), expected, 'trt_sserrvar');
                end
            end
        end

        function testTargetCaptionsStayBesideTheDisclosure(testCase)
            % The gamma target caption (surface target and subject scheme) is kept
            % as its own object with its own text; the disclosure is a second
            % object that covers no panel, title, axis label, caption or button.
            for arm = {'g27ls', 'g14lognu'}
                fh = psyrat_dynrelplot(LS27FigureFixture.dynrelSummary(arm{1}), 'gcoeff', 2);
                cap = findall(fh, 'Tag', LS27FigureFixture.TARGET_TAG);
                testCase.verifyNumElements(cap, 1, [arm{1} ': the target caption must stay one object']);
                if numel(cap) == 1
                    testCase.verifyFalse(contains(lower(LS27FigureFixture.textOf(cap)), 'restricted'), ...
                        [arm{1} ': the disclosure must not be merged into the target caption']);
                end
                testCase.verifyPlacement(fh, arm{1});
            end
            summ = LS27FigureFixture.trt25Summary('lognu');
            for stat = {'dep', 'icc'}
                fig = LS27FigureFixture.subjectPlot(summ, stat{1});
                testCase.verifyPlacement(fig, ['g25lognu ' stat{1}]);
            end
        end

        function testNoDisclosureOutsideTheSixArms(testCase)
            % CONTROL (guard for the implementation). Exact labels under the gamma
            % family only: the Gaussian arms of 14, 25 and 27, a legacy analysis-25
            % result without a family, the gamma one-facet designs 11 and 26, the
            % gamma occasion designs outside the ruling (5, 10, 29) and the label
            % cousins 15, 20, 6 and 8 carry no disclosure and no metadata. The
            % relabeled cases reuse Restricted data under another label, so a
            % silent result there is the gate at work, not missing data.
            for arm = {'n14', 'n27', 'g26ls', 'g26lognu', 'g11ls'}
                fh = psyrat_dynrelplot(LS27FigureFixture.dynrelSummary(arm{1}), 'gcoeff', 2);
                testCase.verifySilent(fh, arm{1});
            end
            s27 = LS27FigureFixture.dynrelSummary('g27ls');
            for lab = {'ic_dodiff_dynrel_sserrvar_trt', 'ic_diff_dynrel_trt', ...
                    'ic_diff_dynrel_sserrvar_trt_rescor', 'ic_dynrel_sserrvar', 'ic_dynrel'}
                fh = psyrat_dynrelplot(LS27FigureFixture.relabel(s27, lab{1}), 'gcoeff', 2);
                testCase.verifySilent(fh, ['analysis-27 data relabeled ' lab{1}]);
            end
            for fam = {'gaussian', 'nofamily'}
                figs = testCase.trt25FromGenerate(fam{1});
                testCase.verifyNumElements(figs, 2, sprintf('analysis 25 %s: two subject plots', fam{1}));
                for k = 1:numel(figs)
                    testCase.verifySilent(figs(k), ['analysis 25 ' fam{1}]);
                end
            end
            s25 = LS27FigureFixture.trt25Summary('ls');
            for lab = {'trt', 'trt_diff', 'ic_sserrvar', 'ic_diff_sserrvar'}
                other = s25;
                other.rel.analysis = lab{1};
                for stat = {'dep', 'icc'}
                    fig = LS27FigureFixture.subjectPlot(other, stat{1});
                    testCase.verifySilent(fig, sprintf('analysis-25 data relabeled %s (%s)', lab{1}, stat{1}));
                end
            end
        end

        function testDisclosureLeavesPlottedContentUnchanged(testCase)
            % CONTROL (guard for the implementation). A Restricted figure and a
            % figure drawn from the same numbers under a non-Restricted label have
            % the same panels, data, colors, titles, axis labels and limits, the
            % same target caption and the same existing UserData; the surface
            % figure keeps its save button. Today they are identical; after the
            % change they may differ only by the disclosure.
            for arm = {'g27ls', 'g14lognu'}
                s = LS27FigureFixture.dynrelSummary(arm{1});
                twin = LS27FigureFixture.relabel(s, 'ic_dynrel_sserrvar');
                fa = psyrat_dynrelplot(s, 'gcoeff', 1);
                fb = psyrat_dynrelplot(twin, 'gcoeff', 1);
                testCase.verifyEqual(LS27FigureFixture.signature(fa), LS27FigureFixture.signature(fb), ...
                    [arm{1} ': plotted content must not depend on the disclosure']);
                ca = findall(fa, 'Tag', LS27FigureFixture.TARGET_TAG);
                cb = findall(fb, 'Tag', LS27FigureFixture.TARGET_TAG);
                testCase.assertNumElements(ca, 1, [arm{1} ': target caption missing']);
                testCase.assertNumElements(cb, 1, [arm{1} ': twin target caption missing']);
                testCase.verifyEqual(LS27FigureFixture.textOf(ca), LS27FigureFixture.textOf(cb), ...
                    [arm{1} ': the target caption text must be unchanged']);
                ua = get(fa, 'UserData'); ub = get(fb, 'UserData');
                testCase.assertTrue(isstruct(ua) && isfield(ua, 'dynrel') && isfield(ua, 'summary'));
                da = ua.dynrel; db = ub.dynrel;
                da.gamma_target = rmfield(da.gamma_target, 'analysis');
                db.gamma_target = rmfield(db.gamma_target, 'analysis');
                testCase.verifyEqual(da, db, [arm{1} ': UserData.dynrel must be unchanged']);
                testCase.verifyEqual(ua.summary.strata, ub.summary.strata, ...
                    [arm{1} ': the stored summary must be unchanged']);
                btn = findall(fa, 'Style', 'pushbutton');
                testCase.verifyTrue(any(contains(string(get(btn, 'String')), 'Save figure')), ...
                    [arm{1} ': the Save figure button must stay']);
            end
            g = LS27FigureFixture.trt25Summary('ls');
            n = g;
            n.rel.family = 'gaussian';
            for stat = {'dep', 'icc'}
                fa = LS27FigureFixture.subjectPlot(g, stat{1});
                fb = LS27FigureFixture.subjectPlot(n, stat{1});
                testCase.verifyEqual(LS27FigureFixture.signature(fa), LS27FigureFixture.signature(fb), ...
                    ['analysis 25 ' stat{1} ': plotted content must not depend on the disclosure']);
            end
        end

        function testSavedPngAndPdfCarryTheDisclosure(testCase)
            % The real save routes. Surface figure: its own "Save figure
            % (PNG/PDF)" button. Subject plot: saveas, the route MATLAB's figure
            % menu takes. The PNG must show the disclosure: deleting the
            % disclosure and saving again removes its ink. A positive control
            % shows the ink measure sees a caption (the target caption) or a
            % label (the subject plot's x label). The PDF must be written; its
            % text is read on the rendered artifacts, not here.
            global PSYRAT_TEST_UIPUTFILE %#ok<GVMIS>
            fh = psyrat_dynrelplot(LS27FigureFixture.dynrelSummary('g27ls'), 'gcoeff', 2);
            btn = findall(fh, 'Style', 'pushbutton');
            testCase.assertNumElements(btn, 1, 'the surface figure must have its one save button');
            pressSave = @(f) localPressSave(btn, f);
            pngA = fullfile(testCase.TmpDir, 'surface_A.png');
            pdfA = fullfile(testCase.TmpDir, 'surface_A.pdf');
            pressSave(pngA); pressSave(pdfA);
            PSYRAT_TEST_UIPUTFILE = [];
            testCase.verifyTrue(isfile(pngA), 'button route: PNG not written');
            testCase.verifyTrue(LS27FigureFixture.isPdf(pdfA), 'button route: PDF not written');
            testCase.verifyTrue(strcmp(char(get(btn, 'Visible')), 'on'), ...
                'the save button must be visible again after saving');
            d = LS27FigureFixture.disclosure(fh);
            testCase.verifyNumElements(d, 1, 'surface figure: no disclosure to find in the saved PNG');
            if numel(d) == 1
                delete(d);
                pngB = fullfile(testCase.TmpDir, 'surface_B.png');
                pressSave(pngB);
                testCase.verifyGreaterThan(LS27FigureFixture.ink(pngA) - LS27FigureFixture.ink(pngB), 300, ...
                    'surface PNG: the disclosure is not visible in the saved image');
            end
            % positive control for the ink measure: the target caption
            cap = findall(fh, 'Tag', LS27FigureFixture.TARGET_TAG);
            testCase.assertNumElements(cap, 1);
            pngC = fullfile(testCase.TmpDir, 'surface_C.png');
            pngD = fullfile(testCase.TmpDir, 'surface_D.png');
            pressSave(pngC); delete(cap); pressSave(pngD);
            PSYRAT_TEST_UIPUTFILE = [];
            testCase.verifyGreaterThan(LS27FigureFixture.ink(pngC) - LS27FigureFixture.ink(pngD), 300, ...
                'control: deleting the target caption must remove ink from the PNG');

            fig = LS27FigureFixture.subjectPlot(LS27FigureFixture.trt25Summary('lognu'), 'dep');
            pngA = fullfile(testCase.TmpDir, 'subject_A.png');
            pdfA = fullfile(testCase.TmpDir, 'subject_A.pdf');
            saveas(fig, pngA); saveas(fig, pdfA);
            testCase.verifyTrue(LS27FigureFixture.isPdf(pdfA), 'menu route: PDF not written');
            d = LS27FigureFixture.disclosure(fig);
            testCase.verifyNumElements(d, 1, 'subject plot: no disclosure to find in the saved PNG');
            if numel(d) == 1
                delete(d);
                pngB = fullfile(testCase.TmpDir, 'subject_B.png');
                saveas(fig, pngB);
                testCase.verifyGreaterThan(LS27FigureFixture.ink(pngA) - LS27FigureFixture.ink(pngB), 300, ...
                    'subject PNG: the disclosure is not visible in the saved image');
            end
            ax = findobj(fig, 'Type', 'axes');
            pngC = fullfile(testCase.TmpDir, 'subject_C.png');
            pngD = fullfile(testCase.TmpDir, 'subject_D.png');
            saveas(fig, pngC);
            set(get(ax(1), 'XLabel'), 'String', '');
            saveas(fig, pngD);
            testCase.verifyGreaterThan(LS27FigureFixture.ink(pngC) - LS27FigureFixture.ink(pngD), 100, ...
                'control: deleting the x label must remove ink from the PNG');
        end

        function testSavedFigKeepsDisclosureMetadataAndExistingCarriers(testCase)
            % A FIG round trip keeps the disclosure (tag, class, text) and its
            % UserData record, beside the carriers the figure already had: the
            % target caption and the stored summary on the surface figure.
            expected = LS27FigureFixture.expectedText();
            summ = LS27FigureFixture.dynrelSummary('g27lognu');
            fh = psyrat_dynrelplot(summ, 'gcoeff', 2);
            capText = LS27FigureFixture.textOf(findall(fh, 'Tag', LS27FigureFixture.TARGET_TAG));
            f1 = fullfile(testCase.TmpDir, 'surface.fig');
            savefig(fh, f1);
            close(fh, 'force');
            fh2 = openfig(f1, 'invisible');
            testCase.verifyDisclosure(fh2, 'reopened g27lognu FIG', expected, 'ic_dynrel_sserrvar_trt');
            cap2 = findall(fh2, 'Tag', LS27FigureFixture.TARGET_TAG);
            testCase.verifyNumElements(cap2, 1, 'reopened FIG: the target caption was lost');
            if numel(cap2) == 1
                testCase.verifyEqual(LS27FigureFixture.textOf(cap2), capText);
            end
            ud = get(fh2, 'UserData');
            testCase.verifyTrue(isstruct(ud) && isfield(ud, 'summary') && isequal(ud.summary, summ), ...
                'reopened FIG: the stored summary must be unchanged');

            fig = LS27FigureFixture.subjectPlot(LS27FigureFixture.trt25Summary('ls'), 'icc');
            f2 = fullfile(testCase.TmpDir, 'subject.fig');
            savefig(fig, f2);
            close(fig, 'force');
            fig2 = openfig(f2, 'invisible');
            testCase.verifyDisclosure(fig2, 'reopened g25ls ICC FIG', expected, 'trt_sserrvar');
        end

        function testBandKeepsEveryExistingObjectAtItsPixelPosition(testCase)
            % Review N3 (layout preservation; no science). The helper promises
            % that the note sits in a band ADDED above the panels, so a
            % Restricted figure is its silent twin (the same numbers under a
            % non-Restricted label) plus that band. Measured in pixels from
            % each figure's lower-left corner, every axes, colorbar, existing
            % caption and button keeps its box within GEOM_TOL_PX, the figure
            % keeps its width, and it grows by exactly the band. A helper that
            % rescaled the panels but did not grow the figure would squeeze
            % every panel by about 11%, which no other method here measures:
            % they check content and non-overlap, not size. Covers
            % the one-dimension surfaces (14 and 27, both parameterizations
            % across the pair), the two-dimension two-stratum surfaces
            % (colorbars, a 1 x 2 grid) and both analysis-25 subject plots.
            for arm = {'g27ls', 'g14lognu'}
                s = LS27FigureFixture.dynrelSummary(arm{1});
                fr = psyrat_dynrelplot(s, 'gcoeff', 2);
                ft = psyrat_dynrelplot(LS27FigureFixture.relabel(s, 'ic_dynrel_sserrvar'), 'gcoeff', 2);
                testCase.verifyKeepsTwinGeometry(fr, ft, [arm{1} ' surface']);
            end
            for arm = {'g14ls', 'g14lognu', 'g27ls', 'g27lognu'}
                s = LS27FigureFixture.twoDimSummary(arm{1});
                fr = psyrat_dynrelplot(s, 'gcoeff', 2);
                ft = psyrat_dynrelplot(LS27FigureFixture.relabel(s, 'ic_dynrel_sserrvar'), 'gcoeff', 2);
                testCase.assertNumElements(findall(ft, 'Type', 'colorbar'), 2, ...
                    [arm{1} ' two-dimension twin: expected one colorbar per panel']);
                testCase.verifyKeepsTwinGeometry(fr, ft, [arm{1} ' two-dimension surface']);
            end
            for fam = {'ls', 'lognu'}
                g = LS27FigureFixture.trt25Summary(fam{1});
                n = g;
                n.rel.family = 'gaussian';
                for stat = {'dep', 'icc'}
                    fr = LS27FigureFixture.subjectPlot(g, stat{1});
                    ft = LS27FigureFixture.subjectPlot(n, stat{1});
                    testCase.verifyKeepsTwinGeometry(fr, ft, sprintf('g25%s %s', fam{1}, stat{1}));
                end
            end
        end

        function testTwoDimensionTwoStratumSurfacesKeepNoteClearOfColorbars(testCase)
            % Review N4. The two-dimension surface (a filled contour with one
            % colorbar per panel) in a two-stratum grid, for analyses 14 and
            % 27, both parameterizations and both coefficients: exactly one
            % disclosure with its metadata, the target caption kept as its own
            % object, and the disclosure clear of every panel, label, caption,
            % button and colorbar with its tick labels. The silent twin keeps
            % both colorbars, so they belong to the plot, not to the note.
            expected = LS27FigureFixture.expectedText();
            for arm = {'g14ls', 'g14lognu', 'g27ls', 'g27lognu'}
                s = LS27FigureFixture.twoDimSummary(arm{1});
                testCase.assertEqual(s.ndim, 2, [arm{1} ': the fixture must be two-dimensional']);
                testCase.assertNumElements(s.strata, 2, [arm{1} ': the fixture must have two strata']);
                if startsWith(arm{1}, 'g27')
                    testCase.assertNotEmpty(s.strata(2).ssrel_table, ...
                        [arm{1} ': analysis 27 must overlay its participants']);
                end
                for gc = [1 2]
                    what = sprintf('%s two-dimension, gcoeff %d', arm{1}, gc);
                    fh = psyrat_dynrelplot(s, 'gcoeff', gc);
                    testCase.verifyNumElements(findobj(fh, 'Type', 'axes'), 2, [what ': expected two panels']);
                    testCase.verifyNumElements(findall(fh, 'Type', 'colorbar'), 2, ...
                        [what ': expected one colorbar per panel']);
                    testCase.verifyNumElements(findall(fh, 'Tag', LS27FigureFixture.TARGET_TAG), 1, ...
                        [what ': the target caption must stay one object']);
                    testCase.verifyDisclosure(fh, what, expected, LS27FigureFixture.label(arm{1}));
                    testCase.verifyPlacement(fh, what);
                end
                tw = psyrat_dynrelplot(LS27FigureFixture.relabel(s, 'ic_dynrel_sserrvar'), 'gcoeff', 2);
                testCase.verifySilent(tw, [arm{1} ' two-dimension data relabeled ic_dynrel_sserrvar']);
                testCase.verifyNumElements(findall(tw, 'Type', 'colorbar'), 2, ...
                    [arm{1} ' two-dimension twin: expected one colorbar per panel']);
            end
        end

    end

    methods (Access = private)

        function fig = surfaceFromViewer(testCase, arm)
            %The surface figure through the dynamic viewer's own button.
            name = 'Dynamic Reliability: View Results';
            close(findall(0, 'Type', 'figure', 'Name', name), 'force');
            psyrat_startview_dynrel('psyrat_prefs', psyrat_defaults(), ...
                'psyrat_data', LS27FigureFixture.dynrelData(arm));
            viewer = findall(0, 'Type', 'figure', 'Name', name);
            testCase.assertNotEmpty(viewer, [arm ': the dynamic viewer did not open']);
            before = findall(0, 'Type', 'figure', 'Tag', 'psyrat_output');
            btns = findall(viewer(1), 'Style', 'pushbutton');
            isPlot = false(size(btns));
            for k = 1:numel(btns)
                s = get(btns(k), 'String');
                if iscell(s); s = strjoin(s, ' '); end
                isPlot(k) = contains(string(s), 'Plot reliability surface');
            end
            hit = btns(isPlot);
            testCase.assertNumElements(hit, 1, [arm ': Plot reliability surface button not found']);
            cb = get(hit, 'Callback');
            cb(hit, []);
            fig = setdiff(findall(0, 'Type', 'figure', 'Tag', 'psyrat_output'), before);
            testCase.assertNumElements(fig, 1, [arm ': expected one new surface figure']);
        end

        function figs = trt25FromGenerate(testCase, fam)
            %The analysis-25 figures through psyrat_relfigures, the call the
            %test-retest viewer's Generate button makes (label 'trt').
            before = findall(0, 'Type', 'figure', 'Tag', 'psyrat_output');
            psyrat_relfigures('psyrat_data', LS27FigureFixture.trt25(fam), ...
                'psyrat_prefs', LS27FigureFixture.trtViewPrefs(), 'analysis', 'trt');
            figs = setdiff(findall(0, 'Type', 'figure', 'Tag', 'psyrat_output'), before);
            testCase.assertNotEmpty(figs, ['analysis 25 ' fam ': Generate drew no figure']);
        end

        function verifyDisclosure(testCase, fig, what, expected, analysis)
            %One disclosure: exportable, visible, the shared text, stating
            %analysis 27's G and D caveat; and the same facts in UserData.
            h = LS27FigureFixture.disclosure(fig);
            testCase.verifyNumElements(h, 1, sprintf( ...
                '%s: expected one Restricted disclosure (Tag %s)', what, LS27FigureFixture.TAG));
            if numel(h) == 1
                testCase.verifyTrue(isa(h, 'matlab.graphics.shape.TextBox'), sprintf( ...
                    '%s: the disclosure must be a textbox annotation (a uicontrol is dropped from PNG/PDF)', what));
                testCase.verifyTrue(strcmp(char(get(h, 'Visible')), 'on'), [what ': disclosure hidden']);
                txt = LS27FigureFixture.textOf(h);
                testCase.verifyEqual(txt, expected, [what ': the disclosure must be the shared text']);
                testCase.verifyTrue(contains(lower(txt), 'restricted'), [what ': the label word is missing']);
                c = LS27WordingRules.clauses(LS27WordingRules.normalize(txt));
                testCase.verifyTrue(any(cellfun(@LS27WordingRules.states27Bias, c)), ...
                    [what ': the disclosure must state analysis 27''s G and D bias']);
            end
            testCase.verifyTrue(LS27FigureFixture.hasMeta(fig), sprintf( ...
                '%s: UserData.%s missing (the saved FIG would lack the disclosure record)', ...
                what, LS27FigureFixture.META));
            if LS27FigureFixture.hasMeta(fig)
                m = fig.UserData.(LS27FigureFixture.META);
                testCase.verifyEqual(m.label, 'RESTRICTED', [what ': metadata label']);
                testCase.verifyEqual(m.analysis, analysis, [what ': metadata analysis label']);
                if numel(h) == 1
                    testCase.verifyEqual(LS27FigureFixture.joinText(m.text), ...
                        LS27FigureFixture.textOf(h), [what ': metadata text']);
                end
            end
        end

        function verifySilent(testCase, fig, what)
            testCase.verifyEmpty(LS27FigureFixture.disclosure(fig), ...
                [what ': must not carry the Restricted disclosure']);
            testCase.verifyFalse(LS27FigureFixture.hasMeta(fig), ...
                [what ': must not carry the Restricted metadata']);
        end

        function verifyPlacement(testCase, fig, what)
            %The disclosure lies inside the figure and covers no axes region
            %(axes plus title and axis labels), no other annotation, no
            %visible uicontrol and no colorbar with its tick labels. drawnow
            %first, so the colorbars MATLAB places during layout are where a
            %saved figure shows them.
            drawnow;
            h = LS27FigureFixture.disclosure(fig);
            testCase.verifyNumElements(h, 1, [what ': no disclosure to place']);
            if numel(h) ~= 1; return; end
            fp = getpixelposition(fig);
            b = LS27FigureFixture.box(h, fig);
            pad = LS27FigureFixture.CLEAR_PX;
            testCase.verifyTrue(b(1) >= 0 && b(2) >= 0 && b(1) + b(3) <= fp(3) + 1 && ...
                b(2) + b(4) <= fp(4) + 1, [what ': the disclosure extends past the figure']);
            axs = findobj(fig, 'Type', 'axes');
            for k = 1:numel(axs)
                testCase.verifyFalse(LS27FigureFixture.overlaps(b, ...
                    LS27FigureFixture.readableBox(axs(k)), pad), ...
                    sprintf('%s: the disclosure covers axes %d (panel, title or labels)', what, k));
            end
            everything = findall(fig);
            others = everything(arrayfun(@(o) isa(o, 'matlab.graphics.shape.TextBox'), everything));
            others = others(others ~= h);
            for k = 1:numel(others)
                testCase.verifyFalse(LS27FigureFixture.overlaps(b, ...
                    LS27FigureFixture.box(others(k), fig), pad), ...
                    sprintf('%s: the disclosure covers another annotation (%s)', what, get(others(k), 'Tag')));
            end
            ui = findall(fig, 'Type', 'uicontrol', 'Visible', 'on');
            for k = 1:numel(ui)
                testCase.verifyFalse(LS27FigureFixture.overlaps(b, ...
                    LS27FigureFixture.box(ui(k), fig), pad), ...
                    sprintf('%s: the disclosure covers a button or control', what));
            end
            cbs = findall(fig, 'Type', 'colorbar');
            for k = 1:numel(cbs)
                testCase.verifyFalse(LS27FigureFixture.overlaps(b, ...
                    LS27FigureFixture.colorbarReadableBox(cbs(k), fig), pad), ...
                    sprintf('%s: the disclosure covers colorbar %d or its tick labels', what, k));
            end
        end

        function verifyKeepsTwinGeometry(testCase, fr, ft, what)
            %Review N3. FR is a Restricted figure and FT its silent twin. After
            %both are laid out (drawnow, so the colorbars have been placed),
            %measured in pixels from each figure's lower-left corner: FR is as
            %wide as FT and taller by exactly the band holding the disclosure
            %(the note starts at least CLEAR_PX above FT's height, ends at
            %least CLEAR_PX below FR's top, and the growth is no more than the
            %note plus those two clearances), and every axes, colorbar,
            %existing caption and button of FR sits at the box of its twin,
            %matched in creation order, within GEOM_TOL_PX.
            drawnow;
            tol = LS27FigureFixture.GEOM_TOL_PX;
            pad = LS27FigureFixture.CLEAR_PX;
            testCase.verifyEmpty(LS27FigureFixture.disclosure(ft), [what ': the twin must be silent']);
            h = LS27FigureFixture.disclosure(fr);
            testCase.assertNumElements(h, 1, [what ': no disclosure to measure']);
            pr = getpixelposition(fr);
            pt = getpixelposition(ft);
            b = LS27FigureFixture.box(h, fr);
            growth = pr(4) - pt(4);
            testCase.verifyEqual(pr(3), pt(3), 'AbsTol', tol, [what ': the figure width must not change']);
            testCase.verifyGreaterThan(growth, 0, sprintf( ...
                '%s: the figure must grow to hold the band (twin %g px tall, Restricted %g px)', ...
                what, pt(4), pr(4)));
            testCase.verifyGreaterThanOrEqual(b(2), pt(4) + pad - tol, sprintf( ...
                ['%s: the disclosure starts %.1f px from the twin''s top edge; the band must ', ...
                'lie wholly above the preserved layout'], what, b(2) - pt(4)));
            testCase.verifyLessThanOrEqual(b(2) + b(4), pr(4) - pad + tol, ...
                [what ': the disclosure must end inside the grown figure']);
            testCase.verifyLessThanOrEqual(growth, b(4) + 2 * pad + tol, sprintf( ...
                '%s: the figure grew %.1f px for a %.1f px disclosure; the growth must be the band', ...
                what, growth, b(4)));
            kinds = { ...
                'axes',     findobj(fr, 'Type', 'axes'),     findobj(ft, 'Type', 'axes'); ...
                'colorbar', findall(fr, 'Type', 'colorbar'), findall(ft, 'Type', 'colorbar'); ...
                'caption',  findall(fr, 'Tag', LS27FigureFixture.TARGET_TAG), ...
                            findall(ft, 'Tag', LS27FigureFixture.TARGET_TAG); ...
                'button',   findall(fr, 'Type', 'uicontrol'), findall(ft, 'Type', 'uicontrol')};
            for g = 1:size(kinds, 1)
                a = kinds{g, 2};
                c = kinds{g, 3};
                testCase.verifyNumElements(a, numel(c), sprintf( ...
                    '%s: the number of %s objects differs from the twin', what, kinds{g, 1}));
                for k = 1:min(numel(a), numel(c))
                    testCase.verifyEqual(LS27FigureFixture.box(a(k), fr), ...
                        LS27FigureFixture.box(c(k), ft), 'AbsTol', tol, sprintf( ...
                        '%s: %s %d moved or changed size against the twin (pixels from the lower-left corner)', ...
                        what, kinds{g, 1}, k));
                end
            end
        end

    end
end

function localPressSave(btn, file)
%Answer the save dialog with FILE and press the figure's own save button.
global PSYRAT_TEST_UIPUTFILE %#ok<GVMIS>
PSYRAT_TEST_UIPUTFILE = file;
cb = get(btn, 'Callback');
cb(btn, []);
end

function localCleanup(shadowDir, tmpDir)
clear global PSYRAT_TEST_UIPUTFILE
if any(strcmp(strsplit(path, pathsep), shadowDir)); rmpath(shadowDir); end
if isfolder(shadowDir); rmdir(shadowDir, 's'); end
if isfolder(tmpDir); rmdir(tmpDir, 's'); end
end
