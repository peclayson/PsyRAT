classdef TestDifferencePresentationContracts < PsyRATTestBase
    % Two presentation contracts on the single-session difference designs
    % (G108 and the RC-31 residual).
    %
    % G108 (analysis 8, ic_diff_sserrvar, summarized as sing_diff_sserr). The
    % difference caterpillar's black reference line, and the red/blue marker rule
    % that reads it, is the STORED group target: group(g).diffscore.pt on the
    % coefficient panel and group(g).diffscore.icc_pt on the ICC panel. The plot
    % used to take its one-facet fallback here, which on a difference table
    % ignores the contrast covariance, the per-event divisors and the trial block,
    % and divides by the full roster's mean MINIMUM paired count. The count basis
    % of the stored target (the selected mean or median of the retained paired
    % counts, or of the full roster when nobody is retained) is recorded beside it
    % in diffscore.obsbasis and stated by an exportable caption. A summary without
    % that record gets a neutral caption. The line is never the overall table's
    % mean of retained participants' coefficients, a different quantity.
    %
    % RC-31 residual (analysis 7, ic_diff, summarized as sing_diff). The viewer
    % shows only the Difference-Score Coefficient popup, and its Generate callback
    % leaves the stored prefs.view.gcoeff alone. The per-event trials curve follows
    % the resolved summary coefficient (relsummary.gcoeff, which the sing_diff
    % branch sets from diffgcoeff), never that hidden preference. Plain 'sing' and
    % the two independent trt_diff selectors are controls and must not move.
    %
    % REFERENCE VALUES. Every expected coefficient, interval, count basis, marker
    % mask and curve comes from G108DiffFixture, which computes it from the raw
    % draw arrays with test-owned arithmetic; no production kernel and no
    % PsyRATAccuracyOracle call supplies an expectation. Comparisons against
    % psyrat_report are PARITY checks against a second production surface and are
    % labelled as such.
    %
    % TAGS
    %   [CONTROL]  an unchanged route that must not move
    %   [CALLBACK] fires the real uicontrol callback with feval on an invisible
    %              figure. This proves the wiring and the drawn values; it is not
    %              a click and not a native look at the rendered figure.
    %
    % The diffscore field name and the caption Tag are constants in
    % G108DiffFixture; the caption phrase rules are in localCaptionRules.

    methods (Test)

        % ================= G108: the stored target and its record =================

        function testStoredGroupTargetMatchesIndependentArithmetic(testCase)
            % [CONTROL] The target the line adopts is the quantity the contract
            % names. For every scenario: the group target's point, interval and
            % ICC; every participant's coefficient, ICC and interval; both paired
            % counts; the inclusion flags; the retained IDs. trlinfo still
            % describes the FULL roster, which is why it cannot supply the count
            % basis.
            [pd, truth] = G108DiffFixture.build(G108DiffFixture.groupA());
            S = localScenarios();
            for k = 1:numel(S)
                sc = S(k);
                out = localSummary(pd, sc);
                ex = G108DiffFixture.expected(truth, sc.gc, sc.cut, sc.meas, sc.cent, 0.95);
                testCase.assertEqual(ex.nretained, sc.nret, sprintf(['%s: the fixture no ' ...
                    'longer yields the retention it was designed for; recheck ' ...
                    'G108DiffFixture.groupA'], sc.label));
                ds = out.relsummary.group(1).diffscore;
                localVerifyTarget(testCase, ds, ex.target, 1e-12, sc.label);
                localVerifyParticipants(testCase, ds.ssrel_table, ex, 1e-12, sc.label);
                testCase.verifyEqual(ds.trlinfo.mean, [mean(ex.n1) mean(ex.n2)], ...
                    'AbsTol', 1e-12, sprintf('%s: trlinfo.mean no longer describes the full roster', sc.label));
                localVerifySameIds(testCase, out.relsummary.group(1).goodids, ...
                    ex.ids(ex.included), sprintf('%s: retained IDs', sc.label));
            end
        end

        function testFixtureExtractorAgreesWithBuilderTruth(testCase)
            % [CONTROL] G108DiffFixture.truthFromRel reads any ic_diff_sserrvar
            % REL (it is how saved results are checked). Prove that extractor
            % against the builder's own arrays, ungrouped and grouped.
            cases = {G108DiffFixture.groupA(), ...
                {G108DiffFixture.groupA('GA'), G108DiffFixture.groupB('GB')}};
            flds = {'n1','n2','u','W1','W2','C','b11','b22','b12','w1','w2','c','popwp'};
            for c = 1:numel(cases)
                [pd, truth] = G108DiffFixture.build(cases{c});
                got = G108DiffFixture.truthFromRel(pd.rel);
                testCase.assertNumElements(got, numel(truth));
                for g = 1:numel(truth)
                    testCase.verifyEqual(got(g).ids, truth(g).ids);
                    testCase.verifyEqual(got(g).name, truth(g).name);
                    for f = 1:numel(flds)
                        testCase.verifyEqual(got(g).(flds{f}), truth(g).(flds{f}), ...
                            'AbsTol', 1e-12, sprintf('case %d group %d field %s', c, g, flds{f}));
                    end
                end
            end
        end

        function testSummaryRetainsTargetBasisMetadata(testCase)
            % The selected counts, the central tendency and the retained/full-
            % roster fallback are recorded in diffscore.(META), and the stored
            % coefficient is exactly the target at the recorded counts. Grouped:
            % one group retained, one on the full roster.
            [pd, truth] = G108DiffFixture.build(G108DiffFixture.groupA());
            S = localScenarios();
            for k = 1:numel(S)
                sc = S(k);
                out = localSummary(pd, sc);
                ex = G108DiffFixture.expected(truth, sc.gc, sc.cut, sc.meas, sc.cent, 0.95);
                localVerifyMeta(testCase, out.relsummary.group(1).diffscore, ex, truth, ...
                    sc.gc, sc.label);
            end

            [pd2, truth2] = G108DiffFixture.build( ...
                {G108DiffFixture.groupA('GA'), G108DiffFixture.groupB('GB')});
            sc = S(1);
            out2 = localSummary(pd2, sc);
            ex2 = G108DiffFixture.expected(truth2, sc.gc, sc.cut, sc.meas, sc.cent, 0.95);
            testCase.assertEqual({ex2.population}, {'retained', 'full_roster'}, ...
                'fixture design: GA retained, GB on the full roster');
            for g = 1:2
                localVerifyMeta(testCase, out2.relsummary.group(g).diffscore, ex2(g), ...
                    truth2(g), sc.gc, sprintf('grouped %s', ex2(g).name));
            end
        end

        % ================= G108: the drawn line and the markers =================

        function testReferenceLineIsStoredGroupTarget(testCase)
            % On both panels the drawn line is the stored group target, bit for
            % bit, and equals the independent arithmetic. The scenarios vary the
            % coefficient, the inclusion measure, the central tendency and the
            % retention (partial, changed, zero). In the primary scenario the
            % fixture is first shown to separate the target from every
            % wrong-but-plausible reference, so a line matching any of them fails.
            [pd, truth] = G108DiffFixture.build(G108DiffFixture.groupA());
            S = localScenarios();
            for k = 1:numel(S)
                sc = S(k);
                out = localSummary(pd, sc);
                ds = out.relsummary.group(1).diffscore;
                ex = G108DiffFixture.expected(truth, sc.gc, sc.cut, sc.meas, sc.cent, 0.95);
                if k == 1
                    localAssertSeparates(testCase, ex, sc.gc);
                end
                spec = {'dep', ds.pt,     ex.target.pt,     ex.alt.prefix; ...
                        'icc', ds.icc_pt, ex.target.icc_pt, ex.alt.prefix_icc};
                for s = 1:2
                    lbl = sprintf('%s (%s)', sc.label, spec{s,1});
                    fig = localPlot(out, spec{s,1});
                    y = localRefLine(testCase, localPanel(testCase, fig, 'diff score'));
                    testCase.verifyEqual(y, [spec{s,2} spec{s,2}], ...
                        sprintf('%s: the reference line is not the stored group target', lbl));
                    testCase.verifyEqual(y(1), spec{s,3}, 'AbsTol', 1e-12, ...
                        sprintf('%s: the reference line is not the independently computed target', lbl));
                    if abs(spec{s,3} - spec{s,4}) > 1e-3
                        testCase.verifyGreaterThan(abs(y(1) - spec{s,4}), 1e-6, ...
                            sprintf('%s: the line is still the one-facet fallback expression', lbl));
                    end
                    close(fig);
                end
            end
        end

        function testMarkerMasksFollowTheDrawnTarget(testCase)
            % Red markers are exactly the participants whose interval lies wholly
            % above or wholly below the target; the rest turn blue; with nobody
            % flagged, one default-colored series. The expected masks use the
            % independently computed participant intervals and target. In the
            % primary scenario the one-facet fallback line is first shown to give
            % a different mask on both panels.
            [pd, truth] = G108DiffFixture.build(G108DiffFixture.groupA());
            S = localScenarios();
            pick = [1 3 8];   % D partial mean (primary), D zero mean, G partial median
            for k = pick
                sc = S(k);
                out = localSummary(pd, sc);
                ex = G108DiffFixture.expected(truth, sc.gc, sc.cut, sc.meas, sc.cent, 0.95);
                P = ex.participants;
                testCase.assertEqual(numel(unique(round(P.pt, 12))), numel(P.pt), ...
                    'fixture: tied point estimates would make the plot order ambiguous');
                testCase.assertEqual(numel(unique(round(P.icc_pt, 12))), numel(P.icc_pt), ...
                    'fixture: tied ICC point estimates would make the plot order ambiguous');
                spec = {'dep', P.pt,     P.ll,     P.ul,     ex.target.pt,     ex.alt.prefix; ...
                        'icc', P.icc_pt, P.icc_ll, P.icc_ul, ex.target.icc_pt, ex.alt.prefix_icc};
                for s = 1:2
                    lbl = sprintf('%s (%s)', sc.label, spec{s,1});
                    [pt, ll, ul, T, L] = spec{s, 2:6};
                    testCase.assertGreaterThan(G108DiffFixture.margin(ll, ul, T), 1e-9, ...
                        sprintf('%s: an interval endpoint sits on the target; the mask is ambiguous', lbl));
                    [red, ord] = G108DiffFixture.redRanks(pt, ll, ul, T);
                    if k == pick(1)
                        testCase.assertNotEqual(G108DiffFixture.redRanks(pt, ll, ul, L), red, ...
                            sprintf('%s: the fixture cannot tell the target mask from the fallback mask', lbl));
                    end

                    fig = localPlot(out, spec{s,1});
                    mk = localMarkers(localPanel(testCase, fig, 'diff score'));
                    n = numel(pt);
                    testCase.verifyEqual(sort(mk.all), 1:n, ...
                        sprintf('%s: every participant must be drawn exactly once', lbl));
                    testCase.verifyEqual(sort(mk.red), red, sprintf('%s: red markers', lbl));
                    if isempty(red)
                        testCase.verifyEqual(mk.n, 1, ...
                            sprintf('%s: a panel with nobody flagged draws one series', lbl));
                    else
                        testCase.verifyEqual(sort(mk.blue), setdiff(1:n, red), ...
                            sprintf('%s: blue markers', lbl));
                    end
                    %the red markers sit at their own participants' estimates
                    [xs, i] = sort(mk.red);
                    ys = mk.redY(i);
                    sorted = pt(ord);
                    testCase.verifyEqual(ys, reshape(sorted(xs), 1, []), 'AbsTol', 1e-12, ...
                        sprintf('%s: red marker heights', lbl));
                    close(fig);
                end
            end
        end

        function testConstantDrawProbeResolves(testCase)
            % The constant-draw control: equal event variances of 1, residual
            % variances of 1, zero covariances and trial effects, ten trials per
            % event. Every participant's generalizability is 10/11 and ICC 1/2,
            % and so is the group target (u = 2, relative error 2/10,
            % single-observation error 2). The one-facet fallback drew 20/21 and
            % 2/3.
            p = PsyRATTestDataFactory.makeICDiffSSErrRelData();
            n = 5;
            V = zeros(n, 2, 2);
            V(:,1,1) = 1;
            V(:,2,2) = 1;
            p.rel.out.id_varcov = {num2cell(V)};
            p.rel.out.trl_varcov = {num2cell(zeros(n, 2, 2))};
            p.rel.out.b_sigma = {num2cell(zeros(n, 2))};
            p.rel.out.err_varcov = {num2cell(V)};
            p.rel.out.er_var_ss1 = {zeros(n, 3)};
            p.rel.out.er_var_ss2 = {zeros(n, 3)};
            p.rel.out.wp_cov_ss = {zeros(n, 3)};
            ids = repelem({'s1'; 's2'; 's3'}, 20);
            ev = repmat([repmat({'loss'}, 10, 1); repmat({'gain'}, 10, 1)], 3, 1);
            p.rel.data = table(ids, ones(60, 1), ev, 'VariableNames', {'id','meas','event'});
            out = psyrat_relsummary('psyrat_data', p, 'analysis', 'sing_diff_sserr', ...
                'depcutoff', 0.5, 'meascutoff', 2, 'depcentmeas', 1, 'diffgcoeff', 2, ...
                'CI', 0.95, 'interactive', false);
            d = out.relsummary.group(1).diffscore;
            testCase.assertEqual(d.ssrel_table.dep_pt, repmat(10/11, 3, 1), 'AbsTol', 1e-12);
            testCase.assertEqual(d.ssrel_table.icc_pt, repmat(0.5, 3, 1), 'AbsTol', 1e-12);
            testCase.verifyEqual(d.pt, 10/11, 'AbsTol', 1e-12);
            testCase.verifyEqual(d.icc_pt, 0.5, 'AbsTol', 1e-12);

            spec = {'dep', 10/11, 20/21; 'icc', 0.5, 2/3};
            for s = 1:2
                fig = localPlot(out, spec{s,1});
                y = localRefLine(testCase, localPanel(testCase, fig, 'diff score'));
                testCase.verifyEqual(y(1), spec{s,2}, 'AbsTol', 1e-12, sprintf( ...
                    '%s: constant-draw control (the one-facet fallback drew %.12f)', spec{s,1}, spec{s,3}));
                close(fig);
            end
        end

        function testTwoGroupPanelsEachDrawOwnGroupTarget(testCase)
            % Each group's panel draws its OWN stored target: GA at its retained
            % counts, GB (nobody retained) at its full roster. The two targets
            % differ, so a copied or swapped line fails.
            [pd, truth] = G108DiffFixture.build( ...
                {G108DiffFixture.groupA('GA'), G108DiffFixture.groupB('GB')});
            S = localScenarios();
            sc = S(1);
            out = localSummary(pd, sc);
            ex = G108DiffFixture.expected(truth, sc.gc, sc.cut, sc.meas, sc.cent, 0.95);
            testCase.assertEqual({ex.population}, {'retained', 'full_roster'}, ...
                'fixture design: GA retained, GB on the full roster');
            testCase.assertGreaterThan(abs(ex(1).target.pt - ex(2).target.pt), 1e-3, ...
                'fixture: the two group targets must differ');
            testCase.assertGreaterThan(abs(ex(1).target.icc_pt - ex(2).target.icc_pt), 1e-3, ...
                'fixture: the two group ICC targets must differ');
            for s = {'dep', 'icc'}
                fig = localPlot(out, s{1});
                for g = 1:2
                    ds = out.relsummary.group(g).diffscore;
                    if strcmp(s{1}, 'dep')
                        stored = ds.pt;
                        indep = ex(g).target.pt;
                    else
                        stored = ds.icc_pt;
                        indep = ex(g).target.icc_pt;
                    end
                    lbl = sprintf('%s panel %s', s{1}, ex(g).name);
                    y = localRefLine(testCase, localPanel(testCase, fig, [ex(g).name ': diff score']));
                    testCase.verifyEqual(y, [stored stored], ...
                        sprintf('%s: not its own stored group target', lbl));
                    testCase.verifyEqual(y(1), indep, 'AbsTol', 1e-12, ...
                        sprintf('%s: not the independently computed target', lbl));
                end
                close(fig);
            end
        end

        % ========================= G108: the caption =========================

        function testCaptionStatesTheTargetBasis(testCase)
            % An exportable caption (a graphics text object or textbox annotation,
            % never a uicontrol, which exportgraphics drops) tagged CAPTION_TAG
            % says what the line is: the group-model difference-score coefficient
            % at the selected mean or median of the retained paired counts, or of
            % the full roster when nobody is retained; on the ICC panel, one
            % observation per event. The caption sits below every panel. Grouped:
            % both bases on one figure.
            [pd, truth] = G108DiffFixture.build(G108DiffFixture.groupA());
            S = localScenarios();
            for k = [1 2 3 4]   % partial mean, partial median, zero mean, zero median
                sc = S(k);
                out = localSummary(pd, sc);
                ex = G108DiffFixture.expected(truth, sc.gc, sc.cut, sc.meas, sc.cent, 0.95);
                for s = {'dep', 'icc'}
                    lbl = sprintf('%s (%s)', sc.label, s{1});
                    fig = localPlot(out, s{1});
                    txt = localCaptionText(testCase, fig, lbl);
                    localCaptionRules(testCase, txt, s{1}, {ex.population}, ex.centmeas, lbl);
                    localVerifyCaptionClearsPanels(testCase, fig, lbl);
                    close(fig);
                end
            end

            [pd2, truth2] = G108DiffFixture.build( ...
                {G108DiffFixture.groupA('GA'), G108DiffFixture.groupB('GB')});
            sc = S(1);
            out2 = localSummary(pd2, sc);
            ex2 = G108DiffFixture.expected(truth2, sc.gc, sc.cut, sc.meas, sc.cent, 0.95);
            for s = {'dep', 'icc'}
                lbl = sprintf('grouped (%s)', s{1});
                fig = localPlot(out2, s{1});
                txt = localCaptionText(testCase, fig, lbl);
                localCaptionRules(testCase, txt, s{1}, {ex2.population}, 'mean', lbl);
                localVerifyCaptionClearsPanels(testCase, fig, lbl);
                close(fig);
            end
        end

        function testCaptionReadsStoredMetadataNotTrlinfo(testCase)
            % The line and caption come from the stored target and its stored
            % count basis, never from trlinfo (the full roster). Poison trlinfo:
            % nothing drawn may change or show the poison. Then flip only the
            % stored central-tendency record: the caption must follow the record.
            % (Preferences cannot reach psyrat_ssrelplot, which takes none; the
            % Generate test covers the preference-to-summary route.)
            [pd, ~] = G108DiffFixture.build(G108DiffFixture.groupA());
            S = localScenarios();
            out = localSummary(pd, S(1));   % D partial mean
            ds = out.relsummary.group(1).diffscore;
            out.relsummary.group(1).diffscore.trlinfo.mean = [999 999];
            out.relsummary.group(1).diffscore.trlinfo.med = [999 999];

            fig = localPlot(out, 'dep');
            y = localRefLine(testCase, localPanel(testCase, fig, 'diff score'));
            testCase.verifyEqual(y, [ds.pt ds.pt], ...
                'poisoned trlinfo: the line must remain the stored group target');
            txt = localCaptionText(testCase, fig, 'poisoned trlinfo');
            testCase.verifyEmpty(regexp(txt, '999', 'once'), 'the caption read trlinfo');
            close(fig);

            meta = G108DiffFixture.META;
            testCase.assertTrue(isfield(ds, meta), sprintf('diffscore.%s is absent', meta));
            out.relsummary.group(1).diffscore.(meta).centmeas = 'median';
            fig = localPlot(out, 'dep');
            txt = localCaptionText(testCase, fig, 'flipped record');
            testCase.verifyNotEmpty(regexp(txt, '\<median\>', 'once'), ...
                'the caption did not follow the stored central-tendency record');
            close(fig);
        end

        function testLegacySummaryWithoutMetadataGetsNeutralCaption(testCase)
            % A summary built before the count basis was recorded carries the
            % stored target but not its counts. The line is still the stored
            % target; the coefficient caption makes no count-basis claim it
            % cannot support (no retained, median or full roster) and says the
            % basis was not recorded or to re-summarize. The ICC needs no record
            % (one observation per event by definition) and keeps its caption.
            meta = G108DiffFixture.META;
            [pd, ~] = G108DiffFixture.build(G108DiffFixture.groupA());
            S = localScenarios();
            out = localSummary(pd, S(2));   % D partial median: the richest claim to suppress
            d = out.relsummary.group(1).diffscore;
            if isfield(d, meta)
                d = rmfield(d, meta);
            end
            out.relsummary.group(1).diffscore = d;

            fig = localPlot(out, 'dep');
            y = localRefLine(testCase, localPanel(testCase, fig, 'diff score'));
            testCase.verifyEqual(y, [d.pt d.pt], ...
                'legacy summary: the line must still be the stored group target');
            txt = localCaptionText(testCase, fig, 'legacy (dep)');
            testCase.verifyEmpty(regexp(txt, '\<retained\>|\<median\>|full roster', 'once'), ...
                'legacy caption asserts a count basis the summary does not record');
            testCase.verifyNotEmpty(regexp(txt, 'not recorded|re-?summari[sz]', 'once'), ...
                'legacy caption must say the basis was not recorded, or to re-summarize');
            close(fig);

            fig = localPlot(out, 'icc');
            y = localRefLine(testCase, localPanel(testCase, fig, 'diff score'));
            testCase.verifyEqual(y, [d.icc_pt d.icc_pt], ...
                'legacy summary: the ICC line must still be the stored ICC target');
            txt = localCaptionText(testCase, fig, 'legacy (icc)');
            localCaptionRules(testCase, txt, 'icc', {}, '', 'legacy (icc)');
            close(fig);
        end

        function testSavedFigureKeepsTargetLineAndCaption(testCase)
            % A FIG round trip keeps the target line and the caption (tag and
            % text). PNG/PDF pixels need a look at the rendered image.
            [pd, ~] = G108DiffFixture.build(G108DiffFixture.groupA());
            S = localScenarios();
            out = localSummary(pd, S(2));   % D partial median
            ds = out.relsummary.group(1).diffscore;
            tmp = [tempname '.fig'];
            cleanFile = onCleanup(@() localDeleteIfExists(tmp)); %#ok<NASGU>

            fig = localPlot(out, 'dep');
            savefig(fig, tmp);
            close(fig);
            fig2 = openfig(tmp, 'invisible');
            y = localRefLine(testCase, localPanel(testCase, fig2, 'diff score'));
            testCase.verifyEqual(y, [ds.pt ds.pt], 'saved FIG: reference line');
            txt = localCaptionText(testCase, fig2, 'saved FIG');
            testCase.verifyNotEmpty(regexp(txt, '\<median\>', 'once'), ...
                'saved FIG: the caption lost its count basis');
            close(fig2);
        end

        function testSavedSyntheticResultResummarizesWithTargetBasis(testCase)
            % A saved result carries draws, not a summary. Save a synthetic
            % two-group analysis-8 result, reopen it, summarize it and plot it:
            % the record, the lines and the caption all come from the fresh
            % summary and match the independent arithmetic.
            [pd, truth] = G108DiffFixture.build( ...
                {G108DiffFixture.groupA('GA'), G108DiffFixture.groupB('GB')});
            re = localSaveAndReload(testCase, pd);
            S = localScenarios();
            sc = S(2);   % D partial median
            out = localSummary(re, sc);
            ex = G108DiffFixture.expected(truth, sc.gc, sc.cut, sc.meas, sc.cent, 0.95);
            testCase.assertEqual({ex.population}, {'retained', 'full_roster'}, ...
                'fixture design: GA retained, GB on the full roster');
            for g = 1:2
                localVerifyMeta(testCase, out.relsummary.group(g).diffscore, ex(g), ...
                    truth(g), sc.gc, sprintf('reopened %s', ex(g).name));
            end
            for s = {'dep', 'icc'}
                fig = localPlot(out, s{1});
                for g = 1:2
                    if strcmp(s{1}, 'dep')
                        indep = ex(g).target.pt;
                    else
                        indep = ex(g).target.icc_pt;
                    end
                    y = localRefLine(testCase, localPanel(testCase, fig, [ex(g).name ': diff score']));
                    testCase.verifyEqual(y(1), indep, 'AbsTol', 1e-12, ...
                        sprintf('reopened %s panel %s', s{1}, ex(g).name));
                end
                txt = localCaptionText(testCase, fig, sprintf('reopened (%s)', s{1}));
                localCaptionRules(testCase, txt, s{1}, {ex.population}, 'median', ...
                    sprintf('reopened (%s)', s{1}));
                close(fig);
            end
        end

        % ====================== G108: the callback route ======================

        function testGenerateCallbackDrawsSelectedMedianTarget(testCase)
            % [CALLBACK] The real analysis-8 viewer Generate callback. Median is
            % selected in the preferences, and the hidden G-Theory preference is
            % set opposite the visible Difference-Score choice. Each group panel
            % draws that group's target at the median of ITS retained counts (GA)
            % or of its full roster (GB); the caption says median.
            [pd, truth] = G108DiffFixture.build( ...
                {G108DiffFixture.groupA('GA'), G108DiffFixture.groupB('GB')});
            ex = G108DiffFixture.expected(truth, 1, 0.50, 2, 2, 0.95);
            testCase.assertEqual({ex.population}, {'retained', 'full_roster'}, ...
                'fixture design: GA retained, GB on the full roster');

            prefs = localViewPrefs(2, 1, 0.50, 20);
            prefs.view.depcentmeas = 2;
            prefs.view.ploticc = 1;
            viewer = localOpenViewer(testCase, pd, prefs);
            localFireGenerate(testCase, viewer);

            figD = localOneFigure(testCase, 'Subject-Level Dependability Estimates', '');
            figI = localOneFigure(testCase, 'Subject-Level ICCs', '');
            for g = 1:2
                ttl = [ex(g).name ': diff score'];
                y = localRefLine(testCase, localPanel(testCase, figD, ttl));
                testCase.verifyEqual(y(1), ex(g).target.pt, 'AbsTol', 1e-12, ...
                    sprintf('Generate: %s coefficient line', ex(g).name));
                y = localRefLine(testCase, localPanel(testCase, figI, ttl));
                testCase.verifyEqual(y(1), ex(g).target.icc_pt, 'AbsTol', 1e-12, ...
                    sprintf('Generate: %s ICC line', ex(g).name));
            end
            txt = localCaptionText(testCase, figD, 'Generate (dep)');
            localCaptionRules(testCase, txt, 'dep', {ex.population}, 'median', 'Generate (dep)');
            txt = localCaptionText(testCase, figI, 'Generate (icc)');
            localCaptionRules(testCase, txt, 'icc', {ex.population}, 'median', 'Generate (icc)');
        end

        % ===================== G108: unchanged-route controls =====================

        function testAnalysis6OneFacetReferenceUnchanged(testCase)
            % [CONTROL] Analysis 6 (ic_sserrvar) keeps its one-facet reference:
            % bp_var/(bp_var + err/mean(trls)) on the coefficient panel and
            % bp_var/(bp_var + err) on the ICC panel, err = pop_errvar plus
            % trl_var under dependability. A regression pin on the existing route,
            % read from the table the plot reads, not a scientific reference. No
            % difference caption appears on this design.
            base = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
            gn = cellstr(string(base.rel.groups(:)));
            drawn = struct();
            for gc = 1:2
                out = psyrat_relsummary('psyrat_data', base, 'analysis', 'sing_sserr', ...
                    'depcutoff', 0.01, 'meascutoff', 2, 'depcentmeas', 1, 'gcoeff', gc, ...
                    'CI', 0.95, 'interactive', false);
                for s = {'dep', 'icc'}
                    fig = localPlot(out, s{1});
                    testCase.verifyEmpty(findall(fig, 'Tag', G108DiffFixture.CAPTION_TAG), ...
                        'analysis 6 must not carry the difference-score caption');
                    for g = 1:numel(gn)
                        t = out.relsummary.group(g).event(1).ssrel_table;
                        err = t.pop_errvar(1);
                        if gc == 1
                            err = err + t.trl_var(1);
                        end
                        if strcmp(s{1}, 'dep')
                            want = t.bp_var(1) / (t.bp_var(1) + err / mean(t.trls));
                        else
                            want = t.bp_var(1) / (t.bp_var(1) + err);
                        end
                        y = localRefLine(testCase, localPanel(testCase, fig, gn{g}));
                        testCase.verifyEqual(y(1), want, 'AbsTol', 1e-12, ...
                            sprintf('analysis 6, %s, gcoeff %d, %s', s{1}, gc, gn{g}));
                        drawn.(sprintf('%s%d%d', s{1}, gc, g)) = y(1);
                    end
                    close(fig);
                end
            end
            testCase.verifyGreaterThan(abs(drawn.dep11 - drawn.dep21), 1e-6, ...
                'anti-vacuity: the trial main effect must separate the two coefficients');
        end

        function testAnalysis25TwoFacetReferenceUnchanged(testCase)
            % [CONTROL] Analysis 25 (trt_sserrvar) keeps its gro_rel / gro_icc
            % route, here through the real trt_sserr summary.
            base = PsyRATTestDataFactory.makeTRTSSErrRelData();
            for gc = 1:2
                out = psyrat_relsummary('psyrat_data', base, 'analysis', 'trt_sserr', ...
                    'gcoeff', gc, 'reltype', 3, 'depcutoff', 0.01, 'meascutoff', 2, ...
                    'depcentmeas', 1, 'noccmode', 1, 'CI', 0.95, 'interactive', false);
                t = out.relsummary.group(1).event(1).ssrel_table;
                testCase.assertTrue(all(ismember({'gro_rel', 'gro_icc'}, t.Properties.VariableNames)), ...
                    'the trt_sserr table must carry gro_rel and gro_icc');
                spec = {'dep', t.gro_rel(1); 'icc', t.gro_icc(1)};
                for s = 1:2
                    fig = localPlot(out, spec{s,1});
                    axs = findobj(fig, 'Type', 'axes');
                    testCase.assertNumElements(axs, 1);
                    y = localRefLine(testCase, axs(1));
                    testCase.verifyEqual(y, [spec{s,2} spec{s,2}], ...
                        sprintf('analysis 25, %s, gcoeff %d', spec{s,1}, gc));
                    testCase.verifyEmpty(findall(fig, 'Tag', G108DiffFixture.CAPTION_TAG), ...
                        'analysis 25 must not carry the difference-score caption');
                    close(fig);
                end
            end
        end

        % ===================== RC-31 residual (analysis 7) =====================

        function testSingDiffCurveFollowsVisibleControlThroughGenerate(testCase)
            % [CALLBACK] Both opposed preference pairs through the real Generate
            % callback, in two forms: the stored preferences already disagree when
            % the viewer opens, and the user changes the visible popup after it
            % opens. The per-event curve, its y-axis label and its window name
            % follow the VISIBLE Difference-Score Coefficient, and the curve
            % equals the independently computed per-event curve (nonzero trial
            % variance).
            [pd, arr] = G108DiffFixture.rc31Build();
            N = 20;
            want = localRc31Expected(testCase, arr, N);
            forms = {'stored preferences disagree', 'popup changed after opening'};
            for v = 1:2
                h = 3 - v;
                for f = 1:numel(forms)
                    if f == 1
                        viewer = localOpenViewer(testCase, pd, localViewPrefs(h, v, 0.30, N));
                        pop = localOnePopup(testCase, viewer);
                    else
                        viewer = localOpenViewer(testCase, pd, localViewPrefs(h, h, 0.30, N));
                        pop = localOnePopup(testCase, viewer);
                        pop.Value = v;
                    end
                    testCase.assertEqual(pop.Value, v, ...
                        'the visible popup must show the choice under test');
                    localFireGenerate(testCase, viewer);
                    lbl = sprintf('visible %s, hidden %s (%s)', G108DiffFixture.coefName(v), ...
                        G108DiffFixture.coefName(h), forms{f});
                    fig = localOneFigure(testCase, ' v Number of Trials', 'Difference Score');
                    localVerifyRc31Figure(testCase, fig, arr.events, want, v, N, lbl);
                    close('all', 'force');
                end
            end
        end

        function testSingDiffGenerateMatchesHeadlessReport(testCase)
            % [CALLBACK] Figure/export parity: the curve Generate draws equals the
            % headless report's curve for the same summary. The report is itself
            % checked against the independent curve first, then the figure against
            % the report, to the report's four-decimal rounding.
            [pd, arr] = G108DiffFixture.rc31Build();
            N = 20;
            want = localRc31Expected(testCase, arr, N);
            for v = 1:2
                h = 3 - v;
                coef = G108DiffFixture.coefName(v);
                viewer = localOpenViewer(testCase, pd, localViewPrefs(h, v, 0.30, N));
                localFireGenerate(testCase, viewer);
                fig = localOneFigure(testCase, ' v Number of Trials', 'Difference Score');

                q = psyrat_relsummary('psyrat_data', pd, 'analysis', 'sing_diff', ...
                    'depcutoff', 0.30, 'meascutoff', 2, 'depcentmeas', 1, ...
                    'gcoeff', h, 'diffgcoeff', v, 'CI', 0.95, 'interactive', false);
                rep = psyrat_report(q, 'ntrials', N);
                crv = rep.tables.curve;
                testCase.verifyEqual(crv.Properties.VariableNames{5}, [coef '_PointEstimate'], ...
                    'control: the headless report must follow the visible coefficient');
                for e = 1:2
                    rows = strcmp(cellstr(string(crv.Event)), arr.events{e});
                    ref = reshape(crv{rows, 5}, 1, []);
                    testCase.verifyEqual(ref, want{v, e}, 'AbsTol', 5e-5 + 1e-12, ...
                        'control: the headless report disagrees with the independent curve');
                    ys = localCurveY(localPanel(testCase, fig, arr.events{e}), N);
                    testCase.assertNumElements(ys, 1);
                    testCase.verifyEqual(ys{1}, ref, 'AbsTol', 5e-5 + 1e-12, sprintf( ...
                        'Generate curve vs headless report (%s, visible %s)', arr.events{e}, coef));
                end
                close('all', 'force');
            end
        end

        function testSingDiffSavedFigureCarriesVisibleCoefficient(testCase)
            % [CALLBACK] The exported figure: the Generate figure saved as FIG and
            % reopened carries the visible coefficient's curve, y-axis label and
            % window name. PNG/PDF pixels need a look at the rendered image.
            [pd, arr] = G108DiffFixture.rc31Build();
            N = 20;
            want = localRc31Expected(testCase, arr, N);
            tmp = [tempname '.fig'];
            cleanFile = onCleanup(@() localDeleteIfExists(tmp)); %#ok<NASGU>
            for v = 1:2
                h = 3 - v;
                viewer = localOpenViewer(testCase, pd, localViewPrefs(h, h, 0.30, N));
                pop = localOnePopup(testCase, viewer);
                pop.Value = v;
                localFireGenerate(testCase, viewer);
                fig = localOneFigure(testCase, ' v Number of Trials', 'Difference Score');
                savefig(fig, tmp);
                close('all', 'force');
                fig2 = openfig(tmp, 'invisible');
                localVerifyRc31Figure(testCase, fig2, arr.events, want, v, N, ...
                    sprintf('saved FIG, visible %s', G108DiffFixture.coefName(v)));
                close(fig2);
            end
        end

        function testSavedSyntheticIcDiffResultFollowsVisibleCoefficient(testCase)
            % [CALLBACK] The same contract on a saved result: a synthetic
            % analysis-7 result saved and reopened, then viewed with the stored
            % preferences disagreeing, draws the visible coefficient's curve.
            [pd, arr] = G108DiffFixture.rc31Build();
            re = localSaveAndReload(testCase, pd);
            N = 20;
            want = localRc31Expected(testCase, arr, N);
            for v = 1:2
                viewer = localOpenViewer(testCase, re, localViewPrefs(3 - v, v, 0.30, N));
                localFireGenerate(testCase, viewer);
                fig = localOneFigure(testCase, ' v Number of Trials', 'Difference Score');
                localVerifyRc31Figure(testCase, fig, arr.events, want, v, N, ...
                    sprintf('reopened result, visible %s', G108DiffFixture.coefName(v)));
                close('all', 'force');
            end
        end

        function testRelfiguresSingDiffUsesResolvedSummaryCoefficient(testCase)
            % The contract location. psyrat_relfigures itself, given preferences
            % whose gcoeff and diffgcoeff disagree, passes the resolved
            % relsummary.gcoeff (= diffgcoeff on sing_diff) to the per-event
            % curve. A fix made only in the viewer callback would fail here.
            [pd, arr] = G108DiffFixture.rc31Build();
            N = 15;
            want = localRc31Expected(testCase, arr, N);
            for v = 1:2
                close('all', 'force');
                psyrat_relfigures('psyrat_data', pd, ...
                    'psyrat_prefs', localViewPrefs(3 - v, v, 0.30, N), 'analysis', 'sing_diff');
                fig = localOneFigure(testCase, ' v Number of Trials', 'Difference Score');
                localVerifyRc31Figure(testCase, fig, arr.events, want, v, N, ...
                    sprintf('psyrat_relfigures, diffgcoeff %d, gcoeff %d', v, 3 - v));
            end
        end

        function testPlainSingCurveStillFollowsGcoeffPopup(testCase)
            % [CONTROL][CALLBACK] Plain one-facet 'sing' keeps following its own
            % visible G-Theory Coefficient popup (prefs.view.gcoeff), whatever the
            % difference-score preference holds. The REL comes from the DoD
            % single-cell builder (production, as the DoD viewer route uses it);
            % the expected curve comes from the source draw arrays.
            src = PsyRATTestDataFactory.makeDoDRelData();
            pd = psyrat_dod_build_virtual_rel('psyrat_data', src, 'mode', 'single', ...
                'groupidx', 1, 'eventidx', 1);
            ev = cellstr(string(pd.rel.events(:)));
            testCase.assertNumElements(ev, 1);
            idv = G108DiffFixture.unwrap(src.rel.out.id_varcov{1});
            trv = G108DiffFixture.unwrap(src.rel.out.trl_varcov{1});
            bs = G108DiffFixture.unwrap(src.rel.out.b_sigma{1});
            N = 20;
            [yD, yG] = G108DiffFixture.onefacetCurves(idv(:,1,1), exp(bs(:,1)).^2, trv(:,1,1), N);
            testCase.assertGreaterThan(max(abs(yD - yG)), 1e-3, ...
                'fixture: dependability and generalizability curves must differ');
            want = {yD; yG};
            for v = 1:2
                viewer = localOpenViewer(testCase, pd, localViewPrefs(3 - v, 3 - v, 0.30, N));
                pop = localOnePopup(testCase, viewer);
                pop.Value = v;
                localFireGenerate(testCase, viewer);
                fig = localOneFigure(testCase, ' v Number of Trials', 'Difference Score');
                coef = G108DiffFixture.coefName(v);
                testCase.verifyTrue(startsWith(localStr(get(fig, 'Name')), [coef ' v ']), ...
                    sprintf('plain sing: window must name %s', coef));
                ax = localPanel(testCase, fig, ev{1});
                testCase.verifyEqual(localStr(ax.YLabel.String), coef);
                ys = localCurveY(ax, N);
                testCase.assertNumElements(ys, 1);
                testCase.verifyEqual(ys{1}, want{v}, 'AbsTol', 1e-12, ...
                    sprintf('plain sing curve, visible gcoeff %d', v));
                close('all', 'force');
            end
        end

        function testTrtDiffSelectorsUnchanged(testCase)
            % [CONTROL] trt_diff keeps its two independent selectors: the
            % per-condition curve follows gcoeff and the difference-score curve
            % follows diffgcoeff, in both opposed pairs, and each figure agrees
            % with its headless table (parity). The primary pins remain
            % TestTrtDiffCurveWiring; this one guards the edit site in
            % psyrat_relfigures.
            pd = localTrtDiffData();
            N = 12;
            pairs = [2 1; 1 2];   % [gcoeff diffgcoeff]
            cond = cell(1, 2);
            for k = 1:size(pairs, 1)
                gc = pairs(k, 1);
                dg = pairs(k, 2);
                close('all', 'force');
                psyrat_relfigures('psyrat_data', pd, 'psyrat_prefs', localTrtPrefs(gc, dg, N), ...
                    'analysis', 'trt_diff');
                figC = localOneFigure(testCase, ' v Number of Trials', 'Difference Score');
                figD = localOneFigure(testCase, 'Difference Score', '');
                testCase.verifyTrue(startsWith(localStr(get(figC, 'Name')), ...
                    [G108DiffFixture.coefName(gc) ' v ']), 'per-condition curve must follow gcoeff');
                testCase.verifyTrue(contains(localStr(get(figD, 'Name')), ...
                    G108DiffFixture.coefName(dg)), 'difference-score curve must follow diffgcoeff');

                out = psyrat_relsummary('psyrat_data', pd, 'analysis', 'trt_diff', ...
                    'gcoeff', gc, 'diffgcoeff', dg, 'reltype', 3, 'relcutoff', 0.5, ...
                    'meascutoff', 2, 'relcentmeas', 1, 'noccmode', 1, 'CI', 0.95, ...
                    'interactive', false);
                rep = psyrat_report(out, 'ntrials', N);
                for e = 1:2
                    ename = pd.rel.events{e};
                    rows = strcmp(cellstr(string(rep.tables.curve.Event)), ename);
                    ys = localCurveY(localPanel(testCase, figC, ename), N);
                    testCase.assertNumElements(ys, 1);
                    testCase.verifyEqual(ys{1}, reshape(rep.tables.curve{rows, 5}, 1, []), ...
                        'AbsTol', 5e-5 + 1e-12, 'per-condition figure vs headless curve');
                    cond{gc}(e, :) = ys{1};
                end
                axD = findobj(figD, 'Type', 'axes');
                testCase.assertNumElements(axD, 1);
                ysD = localCurveY(axD(1), N);
                testCase.assertNumElements(ysD, 1);
                testCase.verifyEqual(ysD{1}, reshape(rep.tables.diffcurve{:, 4}, 1, []), ...
                    'AbsTol', 5e-5 + 1e-12, 'difference figure vs headless difference curve');
            end
            testCase.verifyGreaterThan(max(abs(cond{1}(:) - cond{2}(:))), 1e-3, ...
                'anti-vacuity: gcoeff 1 and 2 must draw different per-condition curves');
        end

    end
end

% =========================================================================
% Local helpers
% =========================================================================

function S = localScenarios()
%G108 scenarios on G108DiffFixture.groupA: coefficient (1 dependability,
%2 generalizability), cutoff, inclusion measure (1 lower limit, 2 point,
%3 upper limit), central tendency (1 mean, 2 median), and the retention the
%fixture was designed to produce (asserted against the independent arithmetic
%in the control test, so a drifted fixture is reported as such).
S = struct( ...
    'label', {'D partial mean', 'D partial median', 'D zero mean', 'D zero median', ...
              'D changed retention', 'D lower-limit rule', 'G partial mean', ...
              'G partial median', 'G zero mean'}, ...
    'gc',    {1, 1, 1, 1, 1, 1, 2, 2, 2}, ...
    'cut',   {0.50, 0.50, 0.95, 0.95, 0.40, 0.50, 0.62, 0.62, 0.995}, ...
    'meas',  {2, 2, 2, 2, 2, 1, 2, 2, 2}, ...
    'cent',  {1, 2, 1, 2, 1, 1, 1, 2, 1}, ...
    'nret',  {4, 4, 0, 0, 5, 2, 4, 4, 0});
end

function out = localSummary(pd, sc)
%The real analysis-8 summary route, headless.
out = psyrat_relsummary('psyrat_data', pd, 'analysis', 'sing_diff_sserr', ...
    'depcutoff', sc.cut, 'meascutoff', sc.meas, 'depcentmeas', sc.cent, ...
    'diffgcoeff', sc.gc, 'CI', 0.95, 'interactive', false);
end

function re = localSaveAndReload(testCase, pd)
%Write the result wrapper as a .psyrat (MAT) file in a temporary location and
%read it back, as the viewer's loader would. The file is removed afterwards.
tmp = [tempname '.psyrat'];
testCase.addTeardown(@() localDeleteIfExists(tmp));
psyrat_data = pd; %#ok<NASGU> saved by name below
save(tmp, 'psyrat_data', '-mat');
z = load(tmp, '-mat');
re = z.psyrat_data;
testCase.assertFalse(isfield(re, 'relsummary'), 'a saved result carries draws, not a summary');
end

function localVerifyTarget(testCase, ds, T, tol, label)
flds = {'pt', 'll', 'ul', 'icc_pt', 'icc_ll', 'icc_ul'};
for f = 1:numel(flds)
    testCase.verifyEqual(ds.(flds{f}), T.(flds{f}), 'AbsTol', tol, sprintf( ...
        '%s: stored group %s differs from the raw-draw arithmetic', label, flds{f}));
end
end

function localVerifyParticipants(testCase, tbl, ex, tol, label)
testCase.assertEqual(height(tbl), numel(ex.ids), sprintf('%s: participant rows', label));
tid = cellstr(string(tbl.id));
P = ex.participants;
for j = 1:numel(ex.ids)
    r = find(strcmp(tid, ex.ids{j}));
    testCase.assertNumElements(r, 1, sprintf('%s: participant %s', label, ex.ids{j}));
    who = sprintf('%s, %s', label, ex.ids{j});
    testCase.verifyEqual([tbl.trls1(r) tbl.trls2(r)], [ex.n1(j) ex.n2(j)], ...
        sprintf('%s: paired counts', who));
    testCase.verifyEqual([tbl.dep_pt(r) tbl.dep_ll(r) tbl.dep_ul(r)], ...
        [P.pt(j) P.ll(j) P.ul(j)], 'AbsTol', tol, sprintf('%s: coefficient and interval', who));
    testCase.verifyEqual([tbl.icc_pt(r) tbl.icc_ll(r) tbl.icc_ul(r)], ...
        [P.icc_pt(j) P.icc_ll(j) P.icc_ul(j)], 'AbsTol', tol, sprintf('%s: ICC and interval', who));
    testCase.verifyEqual(logical(tbl.ind2include(r)), ex.included(j), ...
        sprintf('%s: inclusion flag', who));
end
end

function localVerifySameIds(testCase, a, b, label)
a = cellstr(string(a));
b = cellstr(string(b));
testCase.verifyEqual(numel(a), numel(b), sprintf('%s: count', label));
testCase.verifyEmpty(setxor(a, b), sprintf('%s: membership', label));
end

function localVerifyMeta(testCase, ds, ex, truth, gc, label)
%The record of the stored target's count basis (G108DiffFixture.META).
meta = G108DiffFixture.META;
testCase.verifyTrue(isfield(ds, meta), sprintf(['%s: diffscore.%s is absent, so the ' ...
    'count basis of the stored group target is not recorded'], label, meta));
if ~isfield(ds, meta)
    return;
end
m = ds.(meta);
testCase.verifyEqual(reshape(double(m.obs), 1, []), ex.obs, 'AbsTol', 1e-12, ...
    sprintf('%s: recorded counts', label));
testCase.verifyEqual(char(m.centmeas), ex.centmeas, sprintf('%s: recorded central tendency', label));
testCase.verifyEqual(char(m.population), ex.population, sprintf('%s: retained vs full roster', label));
testCase.verifyEqual(double(m.nretained), ex.nretained, sprintf('%s: retained count', label));
%the recorded counts are the ones the stored coefficient was evaluated at
R = G108DiffFixture.ratios(truth.u, truth.W1, truth.W2, truth.C, truth.b11, truth.b22, ...
    truth.b12, reshape(double(m.obs), 1, []), gc);
testCase.verifyEqual(ds.pt, mean(R), 'AbsTol', 1e-12, ...
    sprintf('%s: the stored coefficient is not the target at the recorded counts', label));
end

function localAssertSeparates(testCase, ex, gc)
%Fixture validity for the primary scenario: the target must sit clearly apart
%from every wrong-but-plausible reference, or an equality check could pass for
%the wrong reason. Asserted, not verified: everything after depends on it.
T = ex.target.pt;
alts = {'prefix',     'the one-facet fallback line'; ...
        'fullroster', 'the target at full-roster counts'; ...
        'othercent',  'the target at the other central tendency'; ...
        'nocov',      'the target without the covariance terms'; ...
        'mincount',   'the target at the mean minimum paired count'; ...
        'subjmean',   'the mean of retained participants'' coefficients'};
if gc == 1
    alts(end+1, :) = {'notrial', 'dependability without the trial block'};
end
for a = 1:size(alts, 1)
    testCase.assertGreaterThan(abs(T - ex.alt.(alts{a,1})), 1e-3, ...
        sprintf('fixture cannot tell the target from %s', alts{a,2}));
end
testCase.assertGreaterThan(abs(ex.target.icc_pt - ex.alt.prefix_icc), 1e-3, ...
    'fixture cannot tell the ICC target from the one-facet fallback ICC line');
end

function fig = localPlot(out, stat)
%The real caterpillar, on a clean figure set.
close('all', 'force');
psyrat_ssrelplot('psyrat_data', out, 'stat', stat);
fig = gcf;
end

function ax = localPanel(testCase, fig, want)
%The one axes on FIG whose title reads WANT.
axs = findobj(fig, 'Type', 'axes');
titles = arrayfun(@(a) localStr(a.Title.String), axs, 'UniformOutput', false);
hit = axs(strcmp(titles, want));
testCase.assertNumElements(hit, 1, sprintf('expected one panel titled "%s"', want));
ax = hit(1);
end

function y = localRefLine(testCase, ax)
%The caterpillar's reference line is the only Line at LineWidth 1.5; the
%participant bars are errorbar objects, not lines.
h = findobj(ax, 'Type', 'line', 'LineWidth', 1.5);
testCase.assertNumElements(h, 1, 'expected exactly one reference line (LineWidth 1.5)');
y = reshape(get(h, 'YData'), 1, []);
end

function mk = localMarkers(ax)
%Plot positions of every drawn participant, of the red ones (with their
%heights) and of the blue ones. Red and blue are the plot's own RGB values.
eb = findobj(ax, 'Type', 'errorbar');
mk = struct('all', zeros(1, 0), 'red', zeros(1, 0), 'redY', zeros(1, 0), ...
    'blue', zeros(1, 0), 'n', numel(eb));
for k = 1:numel(eb)
    x = reshape(get(eb(k), 'XData'), 1, []);
    y = reshape(get(eb(k), 'YData'), 1, []);
    col = get(eb(k), 'Color');
    mk.all = [mk.all x];
    if isequal(col, [1 0 0])
        mk.red = [mk.red x];
        mk.redY = [mk.redY y];
    elseif isequal(col, [0 0 1])
        mk.blue = [mk.blue x];
    end
end
end

function txt = localCaptionText(testCase, fig, label)
%Lower-cased text of every caption object tagged CAPTION_TAG, after checking
%each is exportable (exportgraphics drops UI components).
cap = findall(fig, 'Tag', G108DiffFixture.CAPTION_TAG);
testCase.assertNotEmpty(cap, sprintf(['%s: no caption tagged %s, so the figure does ' ...
    'not say which group target the reference line is'], label, G108DiffFixture.CAPTION_TAG));
parts = strings(0, 1);
for k = 1:numel(cap)
    testCase.verifyFalse(isa(cap(k), 'matlab.ui.control.UIControl'), ...
        sprintf('%s: the caption is a uicontrol, which exportgraphics drops', label));
    testCase.verifyTrue(isa(cap(k), 'matlab.graphics.shape.TextBox') || ...
        isa(cap(k), 'matlab.graphics.primitive.Text'), ...
        sprintf('%s: the caption must be a graphics text object or a textbox annotation', label));
    s = string(get(cap(k), 'String'));
    parts = [parts; s(:)]; %#ok<AGROW>
end
txt = lower(char(strjoin(parts.', ' ')));
end

function localVerifyCaptionClearsPanels(testCase, fig, label)
%Every panel's box sits at least 20 pixels above every caption textbox. This is
%box geometry only: whether the rendered axis-label text clears the caption is
%a matter for a look at the rendered figure.
cap = findall(fig, 'Tag', G108DiffFixture.CAPTION_TAG);
axs = findobj(fig, 'Type', 'axes');
for c = 1:numel(cap)
    if ~isa(cap(c), 'matlab.graphics.shape.TextBox')
        continue;
    end
    oldunits = get(cap(c), 'Units');
    set(cap(c), 'Units', 'pixels');
    cp = get(cap(c), 'Position');
    set(cap(c), 'Units', oldunits);
    for k = 1:numel(axs)
        ap = getpixelposition(axs(k));
        testCase.verifyGreaterThanOrEqual(ap(2), cp(2) + cp(4) + 20, ...
            sprintf('%s: a panel overlaps the caption band', label));
    end
end
end

function localCaptionRules(testCase, txt, stat, populations, centmeas, label)
%The content the caption must carry. Loose on wording, strict on content.
%MATLAB regexp: \< and \> are word boundaries (\b is a backspace).
has = @(pat) ~isempty(regexp(txt, pat, 'once'));
testCase.verifyTrue(has('group'), ...
    sprintf('%s: the caption does not say the line is a group-level target', label));
testCase.verifyTrue(has('difference'), ...
    sprintf('%s: the caption does not name the difference score', label));
if strcmp(stat, 'icc')
    testCase.verifyTrue(has('(one|single|1) observation (per|in each|for each) event'), ...
        sprintf('%s: the ICC caption must say one observation per event', label));
    return;
end
testCase.verifyTrue(has(['\<' centmeas '\>']), ...
    sprintf('%s: the caption must name the %s', label, centmeas));
if strcmp(centmeas, 'mean')
    testCase.verifyFalse(has('\<median\>'), ...
        sprintf('%s: the caption claims a median on a mean run', label));
end
if any(strcmp(populations, 'retained'))
    testCase.verifyTrue(has('\<retained\>'), ...
        sprintf('%s: the caption must say the counts are the retained participants''', label));
end
if any(strcmp(populations, 'full_roster'))
    testCase.verifyTrue(has('full roster'), ...
        sprintf('%s: the caption must say the full roster was used (nobody retained)', label));
else
    testCase.verifyFalse(has('full roster'), ...
        sprintf('%s: the caption claims a full-roster basis on a retained run', label));
end
end

function prefs = localViewPrefs(gcoeff, diffgcoeff, depvalue, ntrials)
%Shipped defaults, then only the fields under test, with every other output off
%so the figure search finds one window. gcoeff and diffgcoeff are parameters:
%RC-31 is only observable when they disagree.
prefs = psyrat_defaults;
prefs.view.gcoeff = gcoeff;
prefs.view.diffgcoeff = diffgcoeff;
prefs.view.depvalue = depvalue;
prefs.view.ntrials = ntrials;
prefs.view.meascutoff = 2;
prefs.view.depcentmeas = 1;
prefs.view.plotdepline = 2;
prefs.view.plotdep = 1;
prefs.view.ploticc = 0;
prefs.view.inctrltable = 0;
prefs.view.tablessrel = 0;
prefs.view.overalltable = 0;
prefs.view.showstddevt = 0;
prefs.view.showstddevf = 0;
end

function viewer = localOpenViewer(testCase, pd, prefs)
close('all', 'force');
psyrat_startview_sing('psyrat_data', pd, 'psyrat_prefs', prefs);
viewer = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui');
testCase.assertNumElements(viewer, 1, 'the single-session viewer did not open one setup window');
end

function pop = localOnePopup(testCase, viewer)
pop = findall(viewer, 'Style', 'popupmenu');
testCase.assertNumElements(pop, 1, 'the viewer should show exactly one coefficient popup');
end

function localFireGenerate(testCase, viewer)
%Fire the viewer's own Generate Figures/Tables callback (feval of the
%uicontrol's Callback cell). Not a click.
btns = findall(viewer, 'Style', 'pushbutton');
hit = [];
for k = 1:numel(btns)
    if contains(localStr(get(btns(k), 'String')), 'Generate')
        hit = btns(k);
        break;
    end
end
testCase.assertNotEmpty(hit, 'no Generate Figures/Tables button on the viewer');
cb = get(hit, 'Callback');
feval(cb{1}, hit, [], cb{2:end});
end

function fig = localOneFigure(testCase, include, exclude)
%The single open figure whose Name contains INCLUDE (and not EXCLUDE). Dialogs
%and the setup window are excluded by name.
figs = findall(groot, 'Type', 'figure');
keep = false(size(figs));
for k = 1:numel(figs)
    nm = localStr(get(figs(k), 'Name'));
    keep(k) = contains(nm, include) && (isempty(exclude) || ~contains(nm, exclude));
end
testCase.assertEqual(nnz(keep), 1, sprintf('expected exactly one figure named like "%s"', include));
fig = figs(find(keep, 1));
end

function ys = localCurveY(ax, N)
%Every curve on AX drawn over n = 1..N (the flat cutoff line has two points).
ys = {};
ls = findobj(ax, 'Type', 'line');
for k = 1:numel(ls)
    x = reshape(get(ls(k), 'XData'), 1, []);
    if numel(x) == N && isequal(x, 1:N)
        ys{end+1} = reshape(get(ls(k), 'YData'), 1, []); %#ok<AGROW>
    end
end
end

function want = localRc31Expected(testCase, arr, N)
%Independent per-event curves, want{coefficient, event}, and the two
%non-vacuity conditions: nonzero trial variance, and dependability and
%generalizability curves that differ materially.
testCase.assertTrue(all(arr.si2(:) > 0), 'fixture: the trial variance must be nonzero');
want = cell(2, 2);
for e = 1:2
    [yD, yG] = G108DiffFixture.onefacetCurves(arr.su2(:,e), arr.se2(:,e), arr.si2(:,e), N);
    want{1, e} = yD;
    want{2, e} = yG;
    testCase.assertGreaterThan(max(abs(yD - yG)), 1e-3, ...
        'fixture: dependability and generalizability curves must differ');
end
end

function localVerifyRc31Figure(testCase, fig, events, want, v, N, label)
%Window name, y-axis label and curve values of an analysis-7 trials figure
%against the visible coefficient V.
coef = G108DiffFixture.coefName(v);
nm = localStr(get(fig, 'Name'));
testCase.verifyTrue(startsWith(nm, [coef ' v ']), sprintf('%s: window reads "%s"', label, nm));
for e = 1:numel(events)
    ax = localPanel(testCase, fig, events{e});
    testCase.verifyEqual(localStr(ax.YLabel.String), coef, ...
        sprintf('%s: %s y-axis label', label, events{e}));
    ys = localCurveY(ax, N);
    testCase.assertNumElements(ys, 1, sprintf('%s: %s should carry one curve', label, events{e}));
    testCase.verifyEqual(ys{1}, want{v, e}, 'AbsTol', 1e-12, ...
        sprintf('%s: %s curve is not the %s curve', label, events{e}, coef));
end
end

function prefs = localTrtPrefs(gcoeff, diffgcoeff, N)
%trt viewer preferences with only the trials curves on (after
%TestTrtDiffCurveWiring's localViewPrefsStruct).
v = struct();
v.gcoeff = gcoeff;
v.reltype = 3;
v.diffgcoeff = diffgcoeff;
v.noccmode = 1;
v.nocc = [];
v.relvalue = 0.5;
v.relcentmeas = 1;
v.meascutoff = 2;
v.plotrel = 1;
v.plotrelline = 2;
v.ploticc = 0;
v.inctrltable = 0;
v.overalltable = 0;
v.ntrials = N;
v.showstddevt = 0;
v.showstddevf = 0;
v.plotssrel = 0;
v.plotssicc = 0;
v.tablessrel = 0;
prefs = struct('view', v, 'ver', '0-test');
end

function pd = localTrtDiffData()
%Synthetic non-concurrent trt_diff REL with constant-across-draws blocks and
%unequal event variances (the pattern of TestTrtDiffCurveWiring.localBuildREL,
%smaller).
ndraw = 40;
nSub = 6;
nOcc = 3;
nTrl = 4;
cov2 = @(s1, s2, r) [s1^2, r * s1 * s2; r * s1 * s2, s2^2];
wrap = @(M) {num2cell(repmat(reshape(M, [1 2 2]), [ndraw 1 1]))};
sigma = [2.0 2.2];

REL = struct();
REL.groups = 'none';
REL.events = {'A', 'B'};
REL.time = arrayfun(@(o) sprintf('o%d', o), 1:nOcc, 'UniformOutput', false);
REL.analysis = 'trt_diff';
REL.diffrescor = 1;
REL.nchains = 1;
REL.niter = ndraw;
REL.filename = 'synthetic_trt_diff';

ids = {};
meas = [];
event = {};
time = {};
for s = 1:nSub
    for o = 1:nOcc
        for e = 1:2
            for t = 1:nTrl
                ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
                event{end+1, 1} = REL.events{e}; %#ok<AGROW>
                time{end+1, 1} = REL.time{o}; %#ok<AGROW>
                meas(end+1, 1) = 0; %#ok<AGROW>
            end
        end
    end
end
REL.data = table(ids, meas, event, time, 'VariableNames', {'id', 'meas', 'event', 'time'});

REL.out = struct();
REL.out.id_varcov  = wrap(cov2(3.0, 2.6, 0.6));
REL.out.tid_varcov = wrap(cov2(2.0, 1.8, 0.4));
REL.out.oid_varcov = wrap(cov2(1.2, 1.1, 0.3));
REL.out.trl_varcov = wrap(cov2(0.9, 0.7, 0.3));
REL.out.occ_varcov = wrap(cov2(0.6, 0.5, 0.3));
REL.out.to_varcov  = wrap(cov2(0.4, 0.3, 0.2));
REL.out.b_sigma    = {num2cell(repmat(log(sigma), ndraw, 1))};
REL.out.err_varcov = wrap(cov2(sigma(1), sigma(2), 0));
REL.out.rescor     = {num2cell(zeros(ndraw, 1))};
REL.out.elabels = {REL.events(:)};
REL.out.glabels = {{'none'}};

pd = struct('rel', REL, 'ver', '0-test');
end

function s = localStr(v)
%Any title/label/name value (char, cellstr, string) as one char row.
v = string(v);
if isempty(v)
    s = '';
else
    s = char(strjoin(v(:).', ' '));
end
end

function localDeleteIfExists(f)
if isfile(f)
    delete(f);
end
end
