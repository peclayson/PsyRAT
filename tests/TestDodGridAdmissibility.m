classdef TestDodGridAdmissibility < matlab.unittest.TestCase
    %TESTDODGRIDADMISSIBILITY The DoD trial-count heatmaps must carry the
    %established component-level admissibility of every cell, label the grid
    %as the plug-in ratio of posterior MEAN variance components (distinct from
    %the table's posterior mean of draw-level ratios), and preserve the raw
    %grid, its components and its reasons in the outputs figure (2B / G11).
    %
    % Finding (PBA-07 / G11, classical-alignment review of 2026-09-12): the grid
    % is masked on range only (finite and inside [0,1]). psyrat_dodiffrel
    % records four quantities that must be variances -- the universe-score
    % contrast variance, the harmonic-scaled between-trial contrast variance,
    % the scaled residual contrast variance and the unscaled (ICC) between-trial
    % contrast variance -- and a negative one does not have to push the final
    % ratio out of [0,1]. The frozen anchor is a PSD trial block [[1,2],[2,4]]
    % at counts [1,100,1,1]: the scaled trial contrast is -.98, the dependability
    % is 100/303 = .330033, and the current mask shows the cell as valid while
    % the adjacent table's kernel flags the same evaluation.
    %
    % Every expected number here is declared independently of the production
    % code: the anchors are frozen in tests/fixtures/dod_grid_reference.json and
    % the whole-grid expectations are closed forms derived by hand from the
    % contrast algebra (see localClosedForms). The fixture factory is used only
    % as a structural carrier; its draws are overwritten with the anchors.
    %
    % Headless, on the real screens: the setup screen's own controls are set
    % and its Generate DoD Outputs callback is fired, the way
    % TestDodHeatmapMask and TestDodRefreshReadsControls do. A callback check is
    % not a native click; the live look is root's.
    %
    % Failure messages are prefixed so the two kinds of red can be told apart:
    %   BEHAVIOR:    the callback produces the wrong value or mask
    %   MISSING API: the figure UserData / caption contract is absent
    % (Phase 1 RED at 025258a: every numeric value control passed; the
    % failures were masks, labels and the absent UserData.)
    %
    % Contract (psyrat_dod_plugin_cell + psyrat_startview_dodiff, 2B / G11):
    %   fig.UserData.dodgrid with fields
    %     grid_estimand, table_estimand, metric ('dep'|'gen'), grid, obs, cvec,
    %     events, group, components (bp_mean, bt_mean, er_mean),
    %     uni_var, bt_unscaled_var, mask_rule, heat(1:2) each with
    %       xevent, yevent, fixed_obs, value, bt_scaled_var, er_scaled_var,
    %       denominator, admissible (logical == AlphaData), zero_count
    %       (logical; components NaN and a 'not computed' marker there),
    %       reason (cell of char; '' for a computed admissible cell)
    %   a text caption in the outputs window containing 'plug-in' and
    %   'posterior mean' (heatmap) and 'draw-level' (table), and the existing
    %   'Grey cells' caption naming the component rule.
    %   fig.UserData.dodtable (Phase 3): the raw dep/gen psyrat_dodiffrel
    %   results as built for the table (or the zero-count stand-ins), the
    %   dep-only admissibility note, invalid = psyrat_admissibility_hasneg on
    %   the dep record, estimand/obs/events/group/record metadata; and a
    %   visible 'Table warning' caption when invalid. The grid mask is a
    %   deterministic check at the component MEANS and the table's record
    %   counts DRAWS, so the two are not interchangeable
    %   (testDrawMixtureWarnsTableWhileGridStaysVisible).
    %
    % Floating-boundary control (root, Phase 1 RED review). On the PSD anchor
    % the scaled trial contrast is exactly 0 in exact arithmetic where
    % n2 = 2 n1, and the existing loop's operation order returns a tiny
    % negative at some of those cells (observed -1.3877787807814457e-17 at
    % (29,58) and -6.938893903907228e-18 at (49,98) on a grid of 100). The
    % scientific target is unchanged: the STRICT computed sign with no
    % tolerance and no clip. So the whole-grid checks compare the mask to the
    % closed-form sign rule AWAY from that boundary, and ON it verify that the
    % retained raw component is ~0 and that the mask equals the strict computed
    % sign (localVerifyComponentMask). Grid 12 and the mapping test's grid 8
    % have no such negatives; the Generalizability test at grid 100 does.

    methods (Test)

        function testAnchorCellIsFlaggedDespiteInRangeRatio(testCase)
            % THE failing case. PSD block, counts [1,100,1,1]: D = 100/303 in
            % [0,1], scaled trial contrast -.98. The value must be the exact
            % formula (no clipping), the cell must be masked, and the reason
            % must name the between-trial contrast variance.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            ref = localReference();
            a = localAnchor(ref, 'psd_block_unequal_counts');
            pd = localAnchorFixture('psd');
            out = localGenerate(testCase, pd, struct('obs', a.obs, 'gridmax', 100, ...
                'metric', 1, 'group', [], 'map', []));

            % the kernel behind the table flags this same evaluation today:
            % positive control that "established admissibility" fires here
            dlg = findall(groot, 'Type', 'figure', 'Name', ...
                'Coefficient is not interpretable as a reliability');
            testCase.assertNotEmpty(dlg, ['precondition: psyrat_dodiffrel must flag ' ...
                'the evaluation at [1 100 1 1] (the table-side note dialog)']);

            [r, c] = localCell(1, 100);
            z = out.im12.CData;
            testCase.verifyEqual(z(r, c), a.dependability, 'AbsTol', 1e-12, ...
                'BEHAVIOR: the cell must be the exact plug-in formula 100/303, unclipped');
            testCase.verifyEqual(z(r, c), 100/303, 'AbsTol', 1e-12, ...
                'reference cross-check: the frozen anchor must equal 100/303');
            testCase.verifyEqual(double(out.im12.AlphaData(r, c)), 0, ...
                ['BEHAVIOR: the cell at trials [1,100] has a negative between-trial ' ...
                'contrast variance (-.98) and must be masked even though D = .330033']);

            % equal-count control on the same screen: [1,1,1,1] is a variance
            % everywhere, D = 1/6, visible
            b = localAnchor(ref, 'psd_block_equal_counts_control');
            [r1, c1] = localCell(1, 1);
            testCase.verifyEqual(z(r1, c1), b.dependability, 'AbsTol', 1e-12, ...
                'control: the equal-count cell must be 1/6');
            testCase.verifyEqual(double(out.im12.AlphaData(r1, c1)), 1, ...
                'control: the equal-count cell must stay visible');

            g = localGrid(testCase, out.fig);
            h = localHeat(testCase, g, 1);
            testCase.verifyEqual(h.bt_scaled_var(r, c), a.bt_scaled, 'AbsTol', 1e-12, ...
                'MISSING API: the scaled between-trial contrast variance of the cell must be retained (-.98)');
            testCase.verifyEqual(h.er_scaled_var(r, c), a.er_scaled, 'AbsTol', 1e-12, ...
                'MISSING API: the scaled residual contrast variance of the cell must be retained (3.01)');
            testCase.verifyEqual(h.denominator(r, c), a.dependability_denominator, 'AbsTol', 1e-12, ...
                'MISSING API: the denominator of the cell must be retained (3.03)');
            testCase.verifyFalse(logical(h.admissible(r, c)), ...
                'MISSING API: the retained admissibility flag must be false for the anchor cell');
            testCase.verifyTrue(contains(h.reason{r, c}, 'between-trial contrast variance'), ...
                sprintf(['MISSING API: the reason must name the between-trial contrast ' ...
                'variance; got "%s"'], h.reason{r, c}));
            testCase.verifyTrue(logical(h.admissible(r1, c1)) && isempty(h.reason{r1, c1}), ...
                'MISSING API: the equal-count control must be admissible with an empty reason');
            testCase.verifyEqual(g.uni_var, a.universe, 'AbsTol', 1e-12, ...
                'MISSING API: the universe-score contrast variance must be retained');
        end

        function testWholeGridMatchesClosedFormsDependability(testCase)
            % Whole-grid check against hand-derived closed forms (no production
            % formula involved). First heatmap: D(n1,n2) = n2/(3 n2 + 3),
            % admissible iff n2 <= 2 n1. Second heatmap (events 1-2 fixed at
            % [1,100]): D = 1/(1.03 + 1/n3 + 1/n4), every cell inadmissible.
            % Zero-count rows/columns are 0 by the existing convention and
            % stay visible.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            gridmax = 12;
            pd = localAnchorFixture('psd');
            out = localGenerate(testCase, pd, struct('obs', [1 100 1 1], 'gridmax', gridmax, ...
                'metric', 1, 'group', [], 'map', []));
            cf = localClosedForms('psd', gridmax);

            testCase.verifyEqual(out.im12.CData, cf.z12_dep, 'AbsTol', 1e-12, ...
                'BEHAVIOR: first heatmap values must equal the closed form n2/(3 n2 + 3) (0 at a zero count)');
            testCase.verifyEqual(out.im34.CData, cf.z34_dep, 'AbsTol', 1e-12, ...
                'BEHAVIOR: second heatmap values must equal the closed form 1/(1.03 + 1/n3 + 1/n4)');
            testCase.verifyEqual(double(out.im34.AlphaData), double(cf.mask34), ...
                'BEHAVIOR: second heatmap must mask every nonzero-count cell (scaled trial contrast -.98 throughout)');

            % this fixture must actually separate the range rule from the
            % component rule, or the test is vacuous
            rangeOnly12 = isfinite(cf.z12_dep) & cf.z12_dep >= 0 & cf.z12_dep <= 1;
            testCase.assertTrue(all(rangeOnly12(:)) && ~all(cf.mask12(:)), ...
                'fixture precondition: every first-heatmap value is in [0,1] yet some cells are inadmissible');

            g = localGrid(testCase, out.fig);
            h1 = localHeat(testCase, g, 1);
            h2 = localHeat(testCase, g, 2);
            localVerifyComponentMask(testCase, out.im12, h1, cf.mask12, cf.boundary12, ...
                'first heatmap (n2 <= 2 n1, zero-count cells visible)');
            testCase.verifyEqual(h1.value, out.im12.CData, ...
                'MISSING API: retained first-heatmap values must be the plotted CData');
            testCase.verifyEqual(logical(h1.admissible), logical(out.im12.AlphaData), ...
                'MISSING API: retained first-heatmap admissibility must be the plotted AlphaData');
            nz = ~cf.zero12;
            testCase.verifyEqual(h1.bt_scaled_var(nz), cf.bt12(nz), 'AbsTol', 1e-12, ...
                'MISSING API: retained scaled between-trial contrast variance must be -1/n1 + 2/n2');
            testCase.verifyEqual(logical(h1.zero_count), cf.zero12, ...
                'MISSING API: zero-count cells must be marked as the convention, not as computed coefficients');
            % a zero-count cell shows 0 by convention; nothing was computed, so
            % its components are NaN and its marker says so (never a
            % calculated error variance of 0)
            testCase.verifyTrue(all(isnan(h1.bt_scaled_var(cf.zero12))) && ...
                all(isnan(h1.er_scaled_var(cf.zero12))) && all(isnan(h1.denominator(cf.zero12))), ...
                'MISSING API: zero-count cells must carry NaN components, not placeholder zeros');
            testCase.verifyTrue(all(cellfun(@(r) contains(r, 'not computed'), h1.reason(cf.zero12))), ...
                'MISSING API: zero-count cells must carry an explicit "not computed" marker');
            testCase.verifyTrue(all(h1.value(cf.zero12) == 0) && all(out.im12.AlphaData(cf.zero12) == 1), ...
                'zero-count cells keep the visible-0 convention on the plotted image');
            testCase.verifyTrue(all(cellfun(@isempty, h1.reason(nz & cf.mask12 & ~cf.boundary12))), ...
                'MISSING API: computed admissible cells (off the floating boundary) carry an empty reason');
            testCase.verifyTrue(all(cellfun(@(r) contains(r, 'between-trial contrast variance'), ...
                h1.reason(nz & ~cf.mask12))), ...
                'MISSING API: every masked first-heatmap cell must name the between-trial contrast variance');
            testCase.verifyEqual(h2.value, out.im34.CData, ...
                'MISSING API: retained second-heatmap values must be the plotted CData');
            testCase.verifyEqual(logical(h2.admissible), logical(out.im34.AlphaData), ...
                'MISSING API: retained second-heatmap admissibility must be the plotted AlphaData');
            testCase.verifyEqual(g.grid, 0:gridmax, ...
                'MISSING API: the grid axis must be retained');
            testCase.verifyEqual(g.obs, [1 100 1 1], ...
                'MISSING API: the entered counts must be retained');
            testCase.verifyEqual(strtrim(lower(g.metric)), 'dep', ...
                'MISSING API: the metric must be recorded on the retained grid');
        end

        function testGeneralizabilityMetricKeepsEvaluationFlag(testCase)
            % Both metric modes. psyrat_dodiffrel builds and records every
            % quantity BEFORE it branches on est, so the record is a property
            % of the evaluation (bp, bt, er, obs), not of the metric; the
            % table-side dialog fires from the dep record whatever heatmap
            % metric is selected. The grid follows the same contract: the
            % G values follow the closed form 1/(3 + 1/n1 + 1/n2), the mask is
            % the same "n2 <= 2 n1", and the reason makes the metric explicit
            % (the flagged quantity does not enter the G denominator).
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            % gridmax 100 so the declared anchor cell (1,100) lies on the grid
            % (harness correction before RED; the closed forms scale with it)
            gridmax = 100;
            pd = localAnchorFixture('psd');
            out = localGenerate(testCase, pd, struct('obs', [1 100 1 1], 'gridmax', gridmax, ...
                'metric', 2, 'group', [], 'map', []));
            cf = localClosedForms('psd', gridmax);
            ref = localReference();
            a = localAnchor(ref, 'psd_block_unequal_counts');

            testCase.verifyEqual(out.im12.CData, cf.z12_gen, 'AbsTol', 1e-12, ...
                'BEHAVIOR: first heatmap G values must equal 1/(3 + 1/n1 + 1/n2)');
            testCase.verifyEqual(out.im34.CData, cf.z34_gen, 'AbsTol', 1e-12, ...
                'BEHAVIOR: second heatmap G values must equal 1/(2.01 + 1/n3 + 1/n4)');
            testCase.verifyEqual(double(out.im34.AlphaData), double(cf.mask34), ...
                'BEHAVIOR: under Generalizability the second heatmap is still masked at every nonzero-count cell');
            [r, c] = localCell(1, 100);
            testCase.verifyEqual(double(out.im12.AlphaData(r, c)), 0, ...
                'BEHAVIOR: under Generalizability the anchor cell (scaled trial contrast -.98) must be masked');
            testCase.verifyEqual(out.im12.CData(r, c), a.generalizability, 'AbsTol', 1e-12, ...
                'BEHAVIOR: the anchor cell G must be 1/4.01');

            caps = localCaptions(out.fig);
            testCase.verifyTrue(any(contains(caps, 'Heatmap metric: Generalizability')), ...
                'the window must state the metric the heatmaps were drawn for');

            g = localGrid(testCase, out.fig);
            testCase.verifyEqual(strtrim(lower(g.metric)), 'gen', ...
                'MISSING API: the retained grid must record the generalizability metric');
            h1 = localHeat(testCase, g, 1);
            % whole-grid mask: closed-form rule away from n2 = 2 n1, strict
            % computed sign on it (this grid has the observed tiny negatives)
            localVerifyComponentMask(testCase, out.im12, h1, cf.mask12, cf.boundary12, ...
                'first heatmap under Generalizability');
            testCase.verifyTrue(contains(h1.reason{r, c}, 'between-trial contrast variance'), ...
                'MISSING API: the reason must still name the flagged quantity under G');
            testCase.verifyTrue(contains(lower(h1.reason{r, c}), 'dependability') || ...
                contains(lower(h1.reason{r, c}), 'generalizability'), ...
                ['MISSING API: the reason must make the metric context explicit (the flagged ' ...
                'between-trial quantity does not enter the G denominator)']);
        end

        function testStrictlyPositiveDefiniteBlockIsFlaggedByProjection(testCase)
            % Robustness: a strictly PD raw block cannot be caught by any PSD
            % test of the raw components; only the projected (harmonic-scaled)
            % component check sees it. bt = [[1,1.9],[1.9,4]] plus .1 on the
            % other two diagonals; at [1,100,1,1] the scaled trial contrast is
            % -.679 and D = 1/3.331. Equal counts are valid (scaled 1.4).
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            ref = localReference();
            a = localAnchor(ref, 'strict_pd_block_unequal_counts');
            pd = localAnchorFixture('strictpd');
            bt = squeeze(pd.rel.out.trl_varcov{1}(1, :, :));
            testCase.assertTrue(all(eig((bt + bt') / 2) > 0), ...
                'fixture precondition: the raw between-trial block must be strictly positive definite');
            out = localGenerate(testCase, pd, struct('obs', a.obs, 'gridmax', 100, ...
                'metric', 1, 'group', [], 'map', []));
            [r, c] = localCell(1, 100);
            [r1, c1] = localCell(1, 1);
            z = out.im12.CData;
            testCase.verifyEqual(z(r, c), a.dependability, 'AbsTol', 1e-12, ...
                'BEHAVIOR: the strict-PD anchor cell must be 1/3.331, unclipped');
            testCase.verifyEqual(double(out.im12.AlphaData(r, c)), 0, ...
                'BEHAVIOR: the strict-PD anchor cell (scaled trial contrast -.679) must be masked');
            testCase.verifyEqual(double(out.im12.AlphaData(r1, c1)), 1, ...
                'control: the equal-count cell of the strict-PD block must stay visible');
            g = localGrid(testCase, out.fig);
            h1 = localHeat(testCase, g, 1);
            testCase.verifyEqual(h1.bt_scaled_var(r, c), a.bt_scaled, 'AbsTol', 1e-12, ...
                'MISSING API: retained scaled between-trial contrast variance must be -.679');
            testCase.verifyEqual(h1.bt_scaled_var(r1, c1), a.equal_count_bt_scaled, 'AbsTol', 1e-12, ...
                'MISSING API: retained equal-count scaled between-trial contrast variance must be 1.4');
            testCase.verifyTrue(g.bt_unscaled_var > 0, ...
                'MISSING API: the unscaled (ICC) contrast variance of a PD block is positive and must be retained as such');
        end

        function testNonPsdRawBlockFlagsThroughUnscaledRecord(testCase)
            % The fourth kernel record: the unscaled (single-observation, ICC)
            % contrast variance. Raw block [[1,3],[3,4]] (det -5) gives
            % c' Sb c = -1 at every cell, while the scaled form -2/n1 + 1/n2 is
            % POSITIVE when n1 > 2 n2. The cell (n1,n2) = (10,1) is therefore
            % flagged by the unscaled record alone; (1,1) by both.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = localAnchorFixture('nonpsd');
            out = localGenerate(testCase, pd, struct('obs', [1 1 1 1], 'gridmax', 10, ...
                'metric', 1, 'group', [], 'map', []));
            [r10, c10] = localCell(10, 1);
            [r1, c1] = localCell(1, 1);
            z = out.im12.CData;
            % (10,1): u = 1, scaled trial = -2/10 + 1 = .8, residual = 1/10 + 3 = 3.1
            testCase.verifyEqual(z(r10, c10), 1 / (1 + .8 + 3.1), 'AbsTol', 1e-12, ...
                'BEHAVIOR: (10,1) must be the exact formula 1/4.9');
            % (1,1): scaled trial = 5 - 6 = -1, residual = 4 -> 1/(1 - 1 + 4)
            testCase.verifyEqual(z(r1, c1), .25, 'AbsTol', 1e-12, ...
                'BEHAVIOR: (1,1) must be the exact formula .25');
            testCase.verifyEqual(double(out.im12.AlphaData(r1, c1)), 0, ...
                'BEHAVIOR: an equal-count cell of a non-PSD raw block must be masked (scaled contrast -1)');
            testCase.verifyEqual(double(out.im12.AlphaData(r10, c10)), 0, ...
                'BEHAVIOR: (10,1) must be masked by the unscaled (ICC) record although its scaled contrast is +.8');
            g = localGrid(testCase, out.fig);
            testCase.verifyEqual(g.bt_unscaled_var, -1, 'AbsTol', 1e-12, ...
                'MISSING API: the unscaled contrast variance must be retained (-1)');
            h1 = localHeat(testCase, g, 1);
            testCase.verifyTrue(contains(h1.reason{r10, c10}, 'single observation'), ...
                'MISSING API: the (10,1) reason must name the single-observation (ICC) quantity');
            testCase.verifyEqual(h1.bt_scaled_var(r10, c10), .8, 'AbsTol', 1e-12, ...
                'MISSING API: the (10,1) scaled between-trial contrast variance must be retained as +.8');
        end

        function testNonfiniteAndOutOfRangeCellsKeepRawValuesAndReasons(testCase)
            % The existing range rule stays as a second layer and its reasons
            % are honest: bp = diag([.5 0 0 0]), PSD block, residual .01 each.
            % At (1,3): scaled trial -1/3, residual .0333.., denominator .2,
            % D = 2.5 (finite, > 1). At (1,5): denominator -.068 -> NaN.
            % At (1,1): D = .5/1.54, valid. Nothing is clipped.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = localAnchorFixture('psd');
            pd = localConstantDraws(pd, 1, diag([.5 0 0 0]), localPsdBlock(), [.01 .01 .01 .01]);
            out = localGenerate(testCase, pd, struct('obs', [1 5 1 1], 'gridmax', 5, ...
                'metric', 1, 'group', [], 'map', []));
            z = out.im12.CData;
            [r3, c3] = localCell(1, 3);
            [r5, c5] = localCell(1, 5);
            [r1, c1] = localCell(1, 1);
            testCase.verifyEqual(z(r3, c3), 2.5, 'AbsTol', 1e-12, ...
                'BEHAVIOR: an out-of-range cell must keep its raw value 2.5 (no clipping)');
            testCase.verifyTrue(isnan(z(r5, c5)), ...
                'BEHAVIOR: a non-positive denominator must stay NaN (no clipping)');
            testCase.verifyEqual(z(r1, c1), .5 / 1.54, 'AbsTol', 1e-12, ...
                'control: the equal-count cell must be .5/1.54');
            testCase.verifyEqual(double(out.im12.AlphaData([r3 r5], c3)), [0; 0], ...
                'BEHAVIOR: out-of-range and nonfinite cells must be masked');
            testCase.verifyEqual(double(out.im12.AlphaData(r1, c1)), 1, ...
                'control: the valid cell must be visible');
            g = localGrid(testCase, out.fig);
            h1 = localHeat(testCase, g, 1);
            testCase.verifyEqual(h1.value(r3, c3), 2.5, 'AbsTol', 1e-12, ...
                'MISSING API: the retained raw value must be 2.5');
            testCase.verifyTrue(isnan(h1.value(r5, c5)), ...
                'MISSING API: the retained raw value must be NaN where the denominator is <= 0');
            testCase.verifyEqual(h1.denominator(r5, c5), .5 - .6 + .032, 'AbsTol', 1e-12, ...
                'MISSING API: the retained denominator must be the raw negative number, not NaN');
            testCase.verifyTrue(contains(h1.reason{r3, c3}, 'between-trial contrast variance') && ...
                contains(h1.reason{r3, c3}, '[0,1]'), ...
                'MISSING API: the (1,3) reason must carry both the negative component and the range failure');
            testCase.verifyTrue(contains(h1.reason{r5, c5}, 'between-trial contrast variance') && ...
                contains(h1.reason{r5, c5}, 'denominator'), ...
                'MISSING API: the (1,5) reason must carry both the negative component and the denominator failure');
        end

        function testZeroUniverseVarianceIsValidZero(testCase)
            % Valid zero-component boundary (numerator): bp = 0 gives u = 0 and
            % a coefficient of exactly 0 -- a legitimate result, not a fault
            % (psyrat_admissibility: zero is counted separately from negative).
            % A no-trial-variance, PSD fixture keeps every other quantity a
            % variance, so every nonzero-count cell is a visible 0.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = localAnchorFixture('psd');
            pd = localConstantDraws(pd, 1, zeros(4), diag([1 1 1 1]), [1 1 1 1]);
            out = localGenerate(testCase, pd, struct('obs', [3 3 3 3], 'gridmax', 6, ...
                'metric', 1, 'group', [], 'map', []));
            testCase.verifyTrue(all(out.im12.CData(:) == 0), ...
                'BEHAVIOR: a zero universe-score variance must give a coefficient of exactly 0 everywhere');
            testCase.verifyTrue(all(out.im12.AlphaData(:) == 1), ...
                'BEHAVIOR: a zero universe-score variance is valid, so nothing may be masked');
            g = localGrid(testCase, out.fig);
            testCase.verifyEqual(g.uni_var, 0, ...
                'MISSING API: the retained universe-score variance must be exactly 0');
            h1 = localHeat(testCase, g, 1);
            [r, c] = localCell(3, 3);
            testCase.verifyTrue(logical(h1.admissible(r, c)) && isempty(h1.reason{r, c}), ...
                'MISSING API: a valid zero must carry no reason');
        end

        function testZeroErrorComponentsGiveValidOne(testCase)
            % Valid zero-component boundary (denominator): no trial variance and
            % a residual SD of exactly zero (log-SD -Inf, the log(0) encoding)
            % drop every error term and the coefficient is exactly 1. Zero in
            % the denominator is a real (extreme) result, not inadmissible.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = localAnchorFixture('psd');
            pd = localConstantDraws(pd, 1, diag([1 0 0 0]), zeros(4), [0 0 0 0]);
            out = localGenerate(testCase, pd, struct('obs', [2 2 2 2], 'gridmax', 4, ...
                'metric', 1, 'group', [], 'map', []));
            z = out.im12.CData;
            nz = localNonzeroCountCells(4);
            testCase.verifyTrue(all(z(nz) == 1), ...
                'BEHAVIOR: with all error components exactly zero the coefficient must be exactly 1');
            testCase.verifyTrue(all(z(~nz) == 0), ...
                'BEHAVIOR: zero-count cells report 0 by the existing convention');
            testCase.verifyTrue(all(out.im12.AlphaData(:) == 1), ...
                'BEHAVIOR: a zero error variance is valid (extreme, not inadmissible), so nothing may be masked');
            g = localGrid(testCase, out.fig);
            h1 = localHeat(testCase, g, 1);
            [r, c] = localCell(2, 2);
            testCase.verifyEqual(h1.er_scaled_var(r, c), 0, ...
                'MISSING API: the retained scaled residual contrast variance must be exactly 0');
            testCase.verifyTrue(logical(h1.admissible(r, c)), ...
                'MISSING API: a zero error component must not be recorded as inadmissible');
        end

        function testPlugInGridIsLabeledDistinctFromDrawLevelTable(testCase)
            % Two equally weighted admissible draws (universe 1 and 9, residual
            % 1 at one trial, no trial variance): the table's posterior mean of
            % draw-level ratios is .7; the grid's plug-in ratio of posterior
            % mean components is 5/6. Both are correct under their own label.
            % Nothing is masked (every component is a variance) -- the cure is
            % the label, not the mask.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            ref = localReference();
            a = localAnchor(ref, 'two_draw_nonlinear_summary');
            pd = localTwoDrawFixture();
            out = localGenerate(testCase, pd, struct('obs', [1 1 1 1], 'gridmax', 3, ...
                'metric', 1, 'group', [], 'map', []));

            tab = findall(out.fig, 'Type', 'uitable');
            testCase.assertEqual(numel(tab), 1, 'expected the one coefficient table');
            data = tab(1).Data;
            testCase.assertTrue(iscell(data) && size(data, 1) >= 2 && size(data, 2) >= 2);
            testCase.assertTrue(strcmp(data{1, 1}, 'Dependability'), 'row 1 must be Dependability');
            testCase.verifyEqual(data{1, 2}, a.table_posterior_mean_of_ratios, 'AbsTol', 1e-9, ...
                'the table must be the posterior mean of draw-level ratios (.7)');

            [r, c] = localCell(1, 1);
            testCase.verifyEqual(out.im12.CData(r, c), a.grid_plug_in_ratio_of_means, 'AbsTol', 1e-12, ...
                'the grid cell at [1,1,1,1] must be the plug-in ratio of posterior mean components (5/6)');
            testCase.verifyEqual(out.im12.CData(r, c), 5/6, 'AbsTol', 1e-12, ...
                'reference cross-check: the frozen anchor must equal 5/6');
            testCase.verifyTrue(all(out.im12.AlphaData(:) == 1) && all(out.im34.AlphaData(:) == 1), ...
                'control: every component is a variance, so no cell may be masked');

            caps = localCaptions(out.fig);
            hasPlugIn = any(contains(lower(caps), 'plug-in') & contains(lower(caps), 'posterior mean'));
            hasDrawLevel = any(contains(lower(caps), 'draw-level'));
            testCase.verifyTrue(hasPlugIn, ...
                'MISSING API: the outputs window must label the heatmaps as the plug-in ratio of posterior mean components');
            testCase.verifyTrue(hasDrawLevel, ...
                'MISSING API: the outputs window must label the table as the posterior mean of draw-level ratios');

            g = localGrid(testCase, out.fig);
            testCase.verifyTrue(ischar(g.grid_estimand) && contains(lower(g.grid_estimand), 'plug-in'), ...
                'MISSING API: the retained grid must carry its own estimand label');
            testCase.verifyTrue(ischar(g.table_estimand) && contains(lower(g.table_estimand), 'draw-level'), ...
                'MISSING API: the retained grid must name the neighboring table''s estimand');
            % posterior mean bp = mean(.25 I, 2.25 I) = 1.25 I, so u = 4 * 1.25 = 5
            testCase.verifyEqual(g.components.bp_mean, 1.25 * eye(4), 'AbsTol', 1e-12, ...
                'MISSING API: the plug-in components (posterior mean bp) must be retained');
            testCase.verifyEqual(g.uni_var, 5, 'AbsTol', 1e-12, ...
                'MISSING API: the plug-in universe-score contrast variance must be 5');
        end

        function testGreyCaptionNamesComponentRule(testCase)
            % The on-screen carrier of the masking rule must say the mask is
            % also driven by a negative variance component, not only by range.
            % TestDodHeatmapMask requires exactly one caption starting with
            % "Grey cells"; that prefix is kept.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = localAnchorFixture('psd');
            out = localGenerate(testCase, pd, struct('obs', [1 100 1 1], 'gridmax', 5, ...
                'metric', 1, 'group', [], 'map', []));
            caps = localCaptions(out.fig);
            isGrey = startsWith(caps, 'Grey cells');
            testCase.assertEqual(nnz(isGrey), 1, 'expected exactly one "Grey cells" caption');
            grey = lower(caps{isGrey});
            testCase.verifyTrue(contains(grey, 'component') || contains(grey, 'negative'), ...
                'MISSING API: the Grey cells caption must name the negative-component rule');
        end

        function testMappingMovesTheFlaggedBlockToTheSecondHeatmap(testCase)
            % Map control: with ERP_1..4 = events C, D, A, B the PSD block sits
            % in contrast slots 3-4 (c3*c4 = -1, the same sign as c1*c2), so the
            % second heatmap follows "n4 <= 2 n3" and the first heatmap, whose
            % slots 3-4 are fixed at [1,100], is inadmissible throughout.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            gridmax = 8;
            pd = localAnchorFixture('psd');
            out = localGenerate(testCase, pd, struct('obs', [1 1 1 100], 'gridmax', gridmax, ...
                'metric', 1, 'group', [], 'map', [3 4 1 2]));
            cf = localClosedForms('psd', gridmax);
            testCase.verifyEqual(out.im34.CData, cf.z12_dep, 'AbsTol', 1e-12, ...
                'BEHAVIOR: after remapping, the second heatmap must carry the n2/(3 n2 + 3) surface');
            testCase.verifyEqual(out.im12.CData, cf.z34_dep, 'AbsTol', 1e-12, ...
                'BEHAVIOR: after remapping, the first heatmap must carry the 1/(1.03 + 1/n + 1/m) surface');
            testCase.verifyEqual(double(out.im12.AlphaData), double(cf.mask34), ...
                'BEHAVIOR: after remapping, the first heatmap must be masked at every nonzero-count cell');
            g = localGrid(testCase, out.fig);
            h2 = localHeat(testCase, g, 2);
            localVerifyComponentMask(testCase, out.im34, h2, cf.mask12, cf.boundary12, ...
                'second heatmap after remapping (n4 <= 2 n3)');
            testCase.verifyEqual(g.events(:)', {'ERP_C', 'ERP_D', 'ERP_A', 'ERP_B'}, ...
                'MISSING API: the retained grid must record the mapped event names');
            testCase.verifyEqual(h2.fixed_obs, [1 1], ...
                'MISSING API: the second heatmap must record the fixed counts of slots 1-2');
        end

        function testHealthySecondGroupIsFullyVisibleUnderComponentRule(testCase)
            % Group control: group 2 keeps the factory's small positive
            % covariances (each cross term is at most .10 against diagonals of
            % at least .88), so every projected quantity is a variance at every
            % count and nothing may be masked, while group 1 carries the
            % anchor. Cross-checked by an independent component evaluation of
            % group 2's posterior mean blocks.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            gridmax = 20;
            pd = localAnchorFixture('psd');
            out = localGenerate(testCase, pd, struct('obs', [3 17 5 2], 'gridmax', gridmax, ...
                'metric', 1, 'group', 2, 'map', []));
            testCase.verifyTrue(all(out.im12.AlphaData(:) == 1) && all(out.im34.AlphaData(:) == 1), ...
                'group 2 (healthy) must have nothing masked');
            testCase.verifyTrue(all(isfinite(out.im12.CData(:))) && all(isfinite(out.im34.CData(:))), ...
                'group 2 (healthy) must be finite everywhere');
            bp = squeeze(mean(pd.rel.out.id_varcov{2}, 1));
            bt = squeeze(mean(pd.rel.out.trl_varcov{2}, 1));
            er = mean(exp(pd.rel.out.b_sigma{2}).^2, 1);
            ok = true;
            for n1 = 1:gridmax
                for n2 = 1:gridmax
                    q = localIndependentComponents(bp, bt, er, [n1 n2 5 2]);
                    ok = ok && q.u >= 0 && q.bt_scaled >= 0 && q.er_scaled >= 0 && q.bt_unscaled >= 0;
                end
            end
            testCase.assertTrue(ok, 'fixture precondition: group 2 is admissible at every first-heatmap cell');
            caps = localCaptions(out.fig);
            testCase.verifyTrue(any(contains(caps, 'Group: Clinical')), ...
                'the window must name the selected group');
            % healthy table: no invalidity warning on the window, none in the
            % retained record, and the record still carries the draw-level structs
            testCase.verifyFalse(any(startsWith(caps, 'Table warning')), ...
                'a healthy draw-level table must carry no invalidity warning');
            t = localTable(testCase, out.fig);
            testCase.verifyFalse(logical(t.invalid), 'healthy control: dodtable.invalid must be false');
            testCase.verifyTrue(isempty(t.adm_note), 'healthy control: the stored note must be empty');
            testCase.verifyTrue(isfield(t.dep, 'admissibility') && isfield(t.dep, 'rel_draws'), ...
                'the retained dep struct must be the raw psyrat_dodiffrel result');
        end

        function testDrawMixtureWarnsTableWhileGridStaysVisible(testCase)
            % The grid mask and the table's draw-level warning are not
            % interchangeable. Equally weighted trial-contrast draws of -.5 and
            % 1.5 (bt = diag(-.5,0,0,0) on half the draws, diag(1.5,0,0,0) on
            % the other half) have a valid MEAN of .5, so every grid cell is
            % drawn, while the table's record counts the negative draws and the
            % table is not interpretable. Universe diag(1,0,0,0), residual 1
            % each, counts [1,1,1,1]: raw table dep = mean(1/4.5, 1/6.5) =
            % 22/117; plug-in cell = 1/(1 + .5 + 4) = 2/11. Nothing recomputed.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = localDrawMixtureFixture();
            out = localGenerate(testCase, pd, struct('obs', [1 1 1 1], 'gridmax', 3, ...
                'metric', 1, 'group', [], 'map', []));
            [r, c] = localCell(1, 1);
            testCase.verifyEqual(out.im12.CData(r, c), 2/11, 'AbsTol', 1e-12, ...
                'the plug-in cell must be 2/11 (mean trial contrast .5)');
            testCase.verifyTrue(all(out.im12.AlphaData(:) == 1) && all(out.im34.AlphaData(:) == 1), ...
                'every grid cell is drawn: the plug-in components are all variances');
            tab = findall(out.fig, 'Type', 'uitable');
            testCase.assertEqual(numel(tab), 1);
            testCase.verifyEqual(tab(1).Data{1, 2}, round(22/117 * 1e4) / 1e4, 'AbsTol', 1e-9, ...
                'the table must show the draw-level mean 22/117 (rounded to 4 dp)');
            caps = localCaptions(out.fig);
            testCase.verifyTrue(any(startsWith(caps, 'Table warning')), ...
                'the table''s draw-level warning must be visible on the window although no grid cell is masked');
            t = localTable(testCase, out.fig);
            testCase.verifyEqual(t.dep.pt, 22/117, 'AbsTol', 1e-12, ...
                'the retained raw dep result must be the unmodified psyrat_dodiffrel mean');
            testCase.verifyTrue(logical(t.invalid), 'dodtable.invalid must be true with negative draws');
            testCase.verifyTrue(contains(t.adm_note, 'between-trial contrast variance'), ...
                'the stored note must name the quantity that went negative');
            e = t.dep.admissibility.entries;
            k = find(strcmp({e.label}, 'between-trial contrast variance'), 1);
            testCase.assertNotEmpty(k, 'the dep record must carry the between-trial entry');
            testCase.verifyEqual(e(k).nneg, size(pd.rel.out.trl_varcov{1}, 1) / 2, ...
                'half the draws are negative, and the count is stored once (dep only), not doubled');
            testCase.verifyEqual(t.obs, [1 1 1 1], 'the entered counts must be recorded with the table');
            testCase.verifyTrue(contains(lower(t.estimand), 'draw-level'), ...
                'the table record must state its draw-level estimand');
        end

        function testZeroCountStandInTableCarriesNoWarning(testCase)
            % An entered count of 0 makes local_dod_summary return the 0-by-
            % convention stand-in with no draws and no admissibility record.
            % That existing convention must be carried as such: no invalidity
            % warning, an empty note, and a record line saying it is a stand-in.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = localAnchorFixture('psd');
            out = localGenerate(testCase, pd, struct('obs', [0 5 1 1], 'gridmax', 5, ...
                'metric', 1, 'group', [], 'map', []));
            caps = localCaptions(out.fig);
            testCase.verifyFalse(any(startsWith(caps, 'Table warning')), ...
                'a zero-count stand-in table must carry no invalidity warning');
            testCase.verifyTrue(any(contains(caps, 'reported as 0')), ...
                'the existing red zero-count note must still be on the window');
            t = localTable(testCase, out.fig);
            testCase.verifyFalse(logical(t.invalid));
            testCase.verifyTrue(isempty(t.adm_note));
            testCase.verifyTrue(contains(lower(t.record), 'stand-in'), ...
                'the record must say the table is the zero-count stand-in');
            testCase.verifyEqual(t.dep.pt, 0, 'the stand-in dep point is 0 by convention');
        end

        function testSavedFigRetainsRawGridReasonsAndCaptions(testCase)
            % The outputs window has no export route of its own, so the saved
            % figure IS the saved output. A real savefig / openfig round trip
            % through a test-owned temporary folder must bring back the raw
            % grid record (values, components, reasons, mask, metadata), the
            % plotted CData / AlphaData, and the captions. This is the FIG
            % carrier only: an image export keeps the rendered mask and
            % captions and nothing else, and the setup screen's table and
            % variance-component exports keep their draw-level meaning.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = localAnchorFixture('psd');
            out = localGenerate(testCase, pd, struct('obs', [1 100 1 1], 'gridmax', 12, ...
                'metric', 1, 'group', [], 'map', []));
            g = localGrid(testCase, out.fig);
            t = localTable(testCase, out.fig);
            % the table's own draw-level warning is on the window before saving
            testCase.assertTrue(any(startsWith(localCaptions(out.fig), 'Table warning')), ...
                'the PSD anchor table must show the draw-level invalidity warning on the window');
            testCase.verifyTrue(logical(t.invalid) && contains(t.adm_note, 'between-trial contrast variance'), ...
                'dodtable must carry the invalidity flag and the full note');
            cdata12 = out.im12.CData;   alpha12 = out.im12.AlphaData;
            cdata34 = out.im34.CData;   alpha34 = out.im34.AlphaData;
            caps = sort(localCaptions(out.fig));
            names = out.names;

            tmp = tempname;
            mkdir(tmp);
            cleanupDir = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
            file = fullfile(tmp, 'dod_outputs.fig');
            savefig(out.fig, file);
            close(out.fig);
            testCase.assertTrue(isfile(file), 'savefig must write the .fig');

            fig2 = openfig(file, 'new', 'invisible');
            ud = get(fig2, 'UserData');
            testCase.assertTrue(isstruct(ud) && isfield(ud, 'dodgrid'), ...
                'the reopened figure must still carry UserData.dodgrid');
            testCase.verifyTrue(isequaln(ud.dodgrid, g), ...
                'the raw grid record (values, components, NaNs, reasons, mask, metadata) must survive the FIG round trip');
            testCase.assertTrue(isfield(ud, 'dodtable'), 'the reopened figure must still carry UserData.dodtable');
            testCase.verifyTrue(isequaln(ud.dodtable, t), ...
                'the raw dep/gen results, note and invalidity flag must survive the FIG round trip');
            testCase.verifyTrue(any(startsWith(localCaptions(fig2), 'Table warning')), ...
                'the table warning must be on the reopened figure (it is also in any image export)');
            % the grid here is 0:12, so the masked cell read back is (1,12)
            [r, c] = localCell(1, 12);
            testCase.verifyTrue(contains(ud.dodgrid.heat(1).reason{r, c}, 'between-trial contrast variance'), ...
                'a masked cell''s reason must be readable from the reopened figure');
            testCase.verifyEqual(ud.dodgrid.heat(1).bt_scaled_var(r, c), -1 + 2/12, 'AbsTol', 1e-12, ...
                'the raw scaled trial contrast must be readable from the reopened figure');

            [im12b, im34b] = localImages(testCase, fig2, names);
            testCase.verifyTrue(isequaln(im12b.CData, cdata12) && isequaln(double(im12b.AlphaData), double(alpha12)), ...
                'first heatmap CData / AlphaData must survive the FIG round trip');
            testCase.verifyTrue(isequaln(im34b.CData, cdata34) && isequaln(double(im34b.AlphaData), double(alpha34)), ...
                'second heatmap CData / AlphaData must survive the FIG round trip');
            testCase.verifyEqual(sort(localCaptions(fig2)), caps, ...
                'every caption (estimand labels, Grey cells rule, metric line) must survive the FIG round trip');
        end

    end
end

%% ---------------------------------------------------------------- fixtures

function ref = localReference()
% The frozen anchors, declared before any production change.
root = fileparts(mfilename('fullpath'));
ref = jsondecode(fileread(fullfile(root, 'fixtures', 'dod_grid_reference.json')));
end

function a = localAnchor(ref, name)
% jsondecode returns a cell array (not a struct array) when the anchor
% objects do not share one field set, which these deliberately do not
anchors = ref.anchors;
if ~iscell(anchors); anchors = num2cell(anchors); end
names = cellfun(@(s) s.name, anchors, 'UniformOutput', false);
k = find(strcmp(names, name), 1);
assert(~isempty(k), 'anchor "%s" not found in dod_grid_reference.json', name);
a = anchors{k};
a.obs = a.obs(:)';
end

function pd = localFixture()
pd = PsyRATTestDataFactory.makeDoDRelData();
if ~isfield(pd, 'ver'); pd.ver = '0-test'; end
if ~isfield(pd, 'proc') || ~isfield(pd.proc, 'measheader')
    pd.proc.measheader = 'meas';
end
end

function bt = localPsdBlock()
% [[1,2],[2,4]] on events 1-2 (rank one, PSD), all other trial entries 0
bt = zeros(4);
bt(1:2, 1:2) = [1 2; 2 4];
end

function pd = localAnchorFixture(kind)
% Group 1 of the factory fixture with its draws overwritten by one of the
% frozen anchor component sets (every draw identical, so the posterior mean
% IS the anchor). Universe diag([1 0 0 0]) gives u = c' Sp c = 1; residual
% variances 1 each. Group 2 is left as the factory made it.
pd = localFixture();
switch kind
    case 'psd'
        bt = localPsdBlock();
    case 'strictpd'
        % det [[1,1.9],[1.9,4]] = .39 > 0; the .1 diagonals keep the full
        % 4 x 4 block strictly positive definite
        bt = diag([0 0 .1 .1]);
        bt(1:2, 1:2) = [1 1.9; 1.9 4];
    case 'nonpsd'
        % det [[1,3],[3,4]] = -5: not a covariance matrix
        bt = zeros(4);
        bt(1:2, 1:2) = [1 3; 3 4];
    otherwise
        error('unknown anchor kind %s', kind);
end
pd = localConstantDraws(pd, 1, diag([1 0 0 0]), bt, [1 1 1 1]);
end

function pd = localConstantDraws(pd, g, bp4, bt4, ervar4)
% Overwrite group g's draws with constant matrices. The kernel and the grid
% both square exp(b_sigma) to get residual variances, so a residual VARIANCE
% v is stored as the log-SD .5*log(v) (log(0) = -Inf encodes exactly zero).
nd = size(pd.rel.out.id_varcov{g}, 1);
pd.rel.out.id_varcov{g} = repmat(reshape(bp4, 1, 4, 4), nd, 1, 1);
pd.rel.out.trl_varcov{g} = repmat(reshape(bt4, 1, 4, 4), nd, 1, 1);
pd.rel.out.b_sigma{g} = repmat(.5 * log(ervar4(:)'), nd, 1);
end

function pd = localTwoDrawFixture()
% Two equally weighted admissible draws: universe contrast variance 1 (bp =
% .25 I, u = 4 * .25) and 9 (bp = 2.25 I), residual SD .5 each (variance .25,
% contrast residual 1 at one trial), no trial variance. Half the factory's
% draws take each value, so the posterior mean bp is 1.25 I (u = 5).
pd = localFixture();
nd = size(pd.rel.out.id_varcov{1}, 1);
assert(mod(nd, 2) == 0, 'fixture precondition: an even number of draws');
bp = zeros(nd, 4, 4);
for d = 1:nd
    if d <= nd / 2
        bp(d, :, :) = .25 * eye(4);
    else
        bp(d, :, :) = 2.25 * eye(4);
    end
end
pd.rel.out.id_varcov{1} = bp;
pd.rel.out.trl_varcov{1} = zeros(nd, 4, 4);
pd.rel.out.b_sigma{1} = repmat(log(.5), nd, 4);
end

function pd = localDrawMixtureFixture()
% Universe diag(1,0,0,0) on every draw, residual variance 1 each, and a
% between-trial block that is diag(-.5,0,0,0) on half the draws and
% diag(1.5,0,0,0) on the other half: posterior MEAN trial contrast .5 (valid
% plug-in), yet half the DRAWS carry a negative trial contrast (table warning).
pd = localFixture();
nd = size(pd.rel.out.id_varcov{1}, 1);
assert(mod(nd, 2) == 0, 'fixture precondition: an even number of draws');
pd = localConstantDraws(pd, 1, diag([1 0 0 0]), zeros(4), [1 1 1 1]);
bt = zeros(nd, 4, 4);
for d = 1:nd
    if d <= nd / 2
        bt(d, 1, 1) = -.5;
    else
        bt(d, 1, 1) = 1.5;
    end
end
pd.rel.out.trl_varcov{1} = bt;
end

%% ------------------------------------------------------- independent algebra

function q = localIndependentComponents(bp, bt, er, obs)
% Hand-written contrast algebra for c = [1 -1 -1 1], independent of the
% production point metric: quadratic forms in the raw blocks, the pairwise
% harmonic-mean scaler on the trial block, own-count scaling on the residual.
c = [1; -1; -1; 1];
q.u = c' * bp * c;
q.bt_unscaled = c' * bt * c;
s = zeros(4);
for i = 1:4
    for j = 1:4
        if i == j
            s(i, j) = bt(i, j) / obs(i);
        else
            s(i, j) = bt(i, j) / (2 / (1 / obs(i) + 1 / obs(j)));
        end
    end
end
q.bt_scaled = c' * s * c;
q.er_scaled = sum((c .^ 2) .* (er(:) ./ obs(:)));
end

function cf = localClosedForms(kind, gridmax)
% Whole-grid expectations for the PSD anchor, derived by hand (see the
% closed_forms block of dod_grid_reference.json). Rows index the y-axis
% count, columns the x-axis count, both from 0; a zero count reports 0 by the
% existing convention and is visible.
assert(strcmp(kind, 'psd'), 'closed forms are declared for the PSD anchor only');
n = 0:gridmax;
[NX, NY] = meshgrid(n, n);          % NX(row,col) = x count, NY(row,col) = y count
nz = NX > 0 & NY > 0;
cf.zero12 = ~nz;
% first heatmap: slots 3-4 at one trial each
%   bt_scaled = 1/n1 + 4/n2 - 4/H(n1,n2) = 1/n1 + 4/n2 - 2(1/n1 + 1/n2) = -1/n1 + 2/n2
%   er_scaled = 1/n1 + 1/n2 + 2
%   D = 1 / (1 + (-1/n1 + 2/n2) + (1/n1 + 1/n2 + 2)) = 1 / (3 + 3/n2) = n2 / (3 n2 + 3)
%   G = 1 / (1 + 1/n1 + 1/n2 + 2)
cf.bt12 = zeros(size(NX));
cf.z12_dep = zeros(size(NX));
cf.z12_gen = zeros(size(NX));
cf.bt12(nz) = -1 ./ NX(nz) + 2 ./ NY(nz);
cf.z12_dep(nz) = NY(nz) ./ (3 * NY(nz) + 3);
cf.z12_gen(nz) = 1 ./ (3 + 1 ./ NX(nz) + 1 ./ NY(nz));
cf.mask12 = ~nz | (NY <= 2 * NX);   % admissible iff n2 <= 2 n1; zero-count cells visible
% theoretical-zero boundary of the scaled trial contrast (exactly 0 in exact
% arithmetic): the floating loop may return a tiny negative there, so the
% mask is compared to the strict computed sign on these cells, not to the
% closed form (see the class comment)
cf.boundary12 = nz & (NY == 2 * NX);
% second heatmap: slots 1-2 fixed at [1,100]
%   bt_scaled = -.98 everywhere; er_scaled = 1 + .01 + 1/n3 + 1/n4
cf.z34_dep = zeros(size(NX));
cf.z34_gen = zeros(size(NX));
cf.z34_dep(nz) = 1 ./ (1.03 + 1 ./ NX(nz) + 1 ./ NY(nz));
cf.z34_gen(nz) = 1 ./ (2.01 + 1 ./ NX(nz) + 1 ./ NY(nz));
cf.mask34 = ~nz;                    % every nonzero-count cell inadmissible
end

function localVerifyComponentMask(testCase, im, h, cfmask, boundary, what)
% Mask against the hand-derived rule, boundary-aware. Away from the
% theoretical-zero boundary the plotted mask must equal the closed-form sign
% rule. ON the boundary the retained raw scaled trial contrast must be ~0 and
% the mask must equal the STRICT computed sign combined with the other
% established rules (residual >= 0, finite value in [0,1]) -- no tolerance and
% no clip anywhere. This is a floating-boundary control on the existing
% arithmetic, not a changed scientific target; the values there must still
% match the closed form (checked by the caller), which proves nothing was
% clipped.
alpha = double(im.AlphaData);
off = ~boundary;
testCase.verifyEqual(alpha(off), double(cfmask(off)), ...
    sprintf('BEHAVIOR: %s: off the boundary the mask must follow the closed-form rule', what));
if any(boundary(:))
    raw = h.bt_scaled_var(boundary);
    testCase.verifyTrue(all(abs(raw) < 1e-12), ...
        sprintf('%s: on the n2 = 2 n1 boundary the raw scaled trial contrast must be ~0 (retained, not rounded)', what));
    v = h.value(boundary);
    strict = double(raw >= 0 & h.er_scaled_var(boundary) >= 0 & isfinite(v) & v >= 0 & v <= 1);
    testCase.verifyEqual(alpha(boundary), strict, ...
        sprintf('%s: on the boundary the mask must equal the STRICT computed sign (no tolerance, no clip)', what));
    % a boundary cell that came out negative must say so with its raw value,
    % never as -0
    neg = find(raw < 0);
    rs = h.reason(boundary);
    for k = 1:numel(neg)
        testCase.verifyTrue(contains(rs{neg(k)}, 'between-trial contrast variance') && ...
            contains(rs{neg(k)}, sprintf('%.15g', raw(neg(k)))) && ~contains(rs{neg(k)}, '(-0)'), ...
            sprintf('%s: a tiny negative boundary cell must report its raw value (%.17g)', what, raw(neg(k))));
    end
end
end

function [r, c] = localCell(nx, ny)
% CData index of the cell at x count nx and y count ny (grid starts at 0)
r = ny + 1;
c = nx + 1;
end

function nz = localNonzeroCountCells(gridmax)
n = 0:gridmax;
[NX, NY] = meshgrid(n, n);
nz = NX > 0 & NY > 0;
end

%% ------------------------------------------------------------ screen driving

function out = localGenerate(testCase, pd, opts)
% Open the setup screen, set its controls, fire its own Generate DoD Outputs
% callback, and return the outputs window with its two heatmaps identified by
% the event names on their x-axis labels (first heatmap = slots 1-2).
psyrat_startview_dodiff('psyrat_prefs', psyrat_defaults, 'psyrat_data', pd);
setup = findall(groot, 'Type', 'figure', 'Name', 'DoD Reliability Setup');
testCase.assertEqual(numel(setup), 1, 'exactly one DoD Reliability Setup screen must be open');
setup = setup(1);
ctl = localControls(testCase, setup, pd);
if ~isempty(opts.group)
    testCase.assertNotEmpty(ctl.group, 'group popup not found');
    ctl.group.Value = opts.group;
end
if ~isempty(opts.map)
    for k = 1:4
        ctl.erp(k).Value = opts.map(k);
    end
end
for k = 1:4
    ctl.trials(k).String = num2str(opts.obs(k));
end
ctl.gridmax.String = num2str(opts.gridmax);
ctl.metric.Value = opts.metric;

btn = localButton(testCase, setup, 'Generate DoD');
cb = get(btn, 'Callback');
feval(cb{1}, btn, [], cb{2:end});

fig = findall(groot, 'Type', 'figure', 'Name', 'DoD Reliability Outputs');
testCase.assertNotEmpty(fig, 'the DoD outputs window did not open');
out.fig = fig(1);

events = cellstr(string(pd.rel.events(:)));
if isempty(opts.map)
    map = arrayfun(@(h) h.Value, ctl.erp);
else
    map = opts.map;
end
out.names = events(map);
[out.im12, out.im34] = localImages(testCase, out.fig, out.names);
end

function [im12, im34] = localImages(testCase, fig, names)
% The two heatmap images of an outputs window, identified by the event name
% on their x-axis labels (first heatmap = contrast slots 1-2).
imgs = findall(fig, 'Type', 'image');
testCase.assertEqual(numel(imgs), 2, 'expected the two heatmap images');
im12 = gobjects(0);
im34 = gobjects(0);
for k = 1:2
    ax = ancestor(imgs(k), 'axes');
    xl = localString(ax.XLabel.String);
    if endsWith(xl, names{1})
        im12 = imgs(k);
    elseif endsWith(xl, names{3})
        im34 = imgs(k);
    end
end
testCase.assertTrue(isscalar(im12) && isscalar(im34), ...
    'could not identify the two heatmaps by their x-axis event names');
testCase.assertEqual(size(im12.AlphaData), size(im12.CData), ...
    'the first heatmap must carry a per-cell mask the size of its data');
testCase.assertEqual(size(im34.AlphaData), size(im34.CData), ...
    'the second heatmap must carry a per-cell mask the size of its data');
end

function ctl = localControls(testCase, fig, pd)
% The setup screen's controls: group popup (may be absent), the four ERP
% popups top to bottom, the metric popup, the four trial edits (by the row
% they share with an ERP popup) and the one remaining edit (grid max).
groups = cellstr(string(pd.rel.groups(:)));
events = cellstr(string(pd.rel.events(:)));
pops = findall(fig, 'Style', 'popupmenu');
ctl.group = gobjects(0);
ctl.metric = gobjects(0);
erp = gobjects(1, 0);
for k = 1:numel(pops)
    s = cellstr(pops(k).String);
    if isequal(s(:), groups(:))
        ctl.group = pops(k);
    elseif isequal(s(:), events(:))
        erp(end+1) = pops(k); %#ok<AGROW>
    elseif isequal(s(:), {'Dependability'; 'Generalizability'})
        ctl.metric = pops(k);
    end
end
testCase.assertEqual(numel(erp), 4, 'four ERP popups expected');
testCase.assertNotEmpty(ctl.metric, 'metric popup not found');
[~, order] = sort(arrayfun(@(h) h.Position(2), erp), 'descend');
ctl.erp = erp(order);
edits = findall(fig, 'Style', 'edit');
ys = arrayfun(@(h) h.Position(2), edits);
used = false(size(edits));
ctl.trials = gobjects(1, 4);
for k = 1:4
    onrow = abs(ys - ctl.erp(k).Position(2)) < 1e-6;
    testCase.assertEqual(nnz(onrow), 1, sprintf('ERP row %d must carry exactly one edit field', k));
    ctl.trials(k) = edits(onrow);
    used = used | onrow;
end
rest = edits(~used);
testCase.assertEqual(numel(rest), 1, 'exactly one grid max edit field expected');
ctl.gridmax = rest;
end

function btn = localButton(testCase, fig, label)
btns = findall(fig, 'Style', 'pushbutton');
btn = gobjects(0);
for k = 1:numel(btns)
    s = localString(btns(k).String);
    if contains(s, label)
        btn = btns(k);
    end
end
testCase.assertNotEmpty(btn, sprintf('button "%s" not found', label));
end

function caps = localCaptions(fig)
% every text uicontrol's String in the outputs window, one char row each
h = findall(fig, 'Type', 'uicontrol', 'Style', 'text');
caps = arrayfun(@(t) localString(t.String), h, 'UniformOutput', false);
caps = caps(:);
end

function s = localString(s)
if iscell(s); s = strjoin(s, ' '); end
s = char(string(s));
end

%% ------------------------------------------------------- planned contract

function g = localGrid(testCase, fig)
% The retained grid record. Its absence is reported as MISSING API so root
% can separate it from a BEHAVIOR failure on the plotted objects.
ud = get(fig, 'UserData');
testCase.assertTrue(isstruct(ud) && isfield(ud, 'dodgrid') && isstruct(ud.dodgrid), ...
    'MISSING API: the outputs figure must carry UserData.dodgrid (saved by savefig)');
g = ud.dodgrid;
needed = {'grid_estimand', 'table_estimand', 'metric', 'grid', 'obs', 'cvec', 'events', ...
    'group', 'components', 'uni_var', 'bt_unscaled_var', 'mask_rule', 'heat'};
for k = 1:numel(needed)
    testCase.assertTrue(isfield(g, needed{k}), ...
        sprintf('MISSING API: UserData.dodgrid must carry the field "%s"', needed{k}));
end
testCase.assertEqual(numel(g.heat), 2, 'MISSING API: dodgrid.heat must describe the two heatmaps');
end

function t = localTable(testCase, fig)
% The retained draw-level table record: the raw dep/gen psyrat_dodiffrel
% results (or the zero-count stand-ins), the dep-only admissibility note,
% the invalidity flag and the entered-count metadata.
ud = get(fig, 'UserData');
testCase.assertTrue(isstruct(ud) && isfield(ud, 'dodtable') && isstruct(ud.dodtable), ...
    'MISSING API: the outputs figure must carry UserData.dodtable');
t = ud.dodtable;
needed = {'estimand', 'obs', 'events', 'group', 'dep', 'gen', 'adm_note', 'invalid', 'record'};
for k = 1:numel(needed)
    testCase.assertTrue(isfield(t, needed{k}), ...
        sprintf('MISSING API: UserData.dodtable must carry the field "%s"', needed{k}));
end
end

function h = localHeat(testCase, g, k)
h = g.heat(k);
needed = {'xevent', 'yevent', 'fixed_obs', 'value', 'bt_scaled_var', 'er_scaled_var', ...
    'denominator', 'admissible', 'zero_count', 'reason'};
for j = 1:numel(needed)
    testCase.assertTrue(isfield(h, needed{j}), ...
        sprintf('MISSING API: dodgrid.heat(%d) must carry the field "%s"', k, needed{j}));
end
testCase.assertTrue(iscell(h.reason) && isequal(size(h.reason), size(h.value)), ...
    'MISSING API: dodgrid.heat.reason must be a cell array the size of the grid');
end
