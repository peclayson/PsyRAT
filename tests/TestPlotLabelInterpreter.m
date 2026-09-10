classdef TestPlotLabelInterpreter < PsyRATTestBase
    %Every plotting function that prints a user-supplied event or group name
    %must render it with the text interpreter set to 'none'. MATLAB's default
    %TeX interpreter turns the underscore in a name such as inc_cor into a
    %subscript ("inc" with a small "c"), which is exactly what the manual's
    %first plot captures showed (2026-09-03) in panel titles, legends, and
    %the forest plot's event tick labels.
    %
    %Each test builds a fixture whose group and event names carry an
    %underscore, draws the figure through the production plotter, and reads
    %the Interpreter property back from every text object that prints one of
    %those names. Literal titles ('Measurement', 'diff score') are not under
    %test; they hold no user text.
    %
    %The positive control proves the assertions can fail: a bare title(),
    %legend, and axes on a fresh figure all read 'tex', so an assertion of
    %'none' is falsifiable and a plotter that drops the property would trip
    %it.

    properties (Constant)
        Groups = {'grp_a'; 'grp_b'};
        Events = {'inc_cor'; 'inc_err'};
    end

    methods (Test)

        function testPositiveControlDefaultIsTex(testCase)
            %Anti-vacuity: without the property, MATLAB defaults to 'tex'.
            fig = figure('Visible', 'off');
            ax = axes(fig);
            plot(ax, 1:3, 1:3);
            t = title(ax, 'a_b');
            lg = legend(ax, 'c_d');
            testCase.verifyEqual(t.Interpreter, 'tex');
            testCase.verifyEqual(lg.Interpreter, 'tex');
            testCase.verifyEqual(ax.TickLabelInterpreter, 'tex');
            testCase.verifyEqual(ax.Legend, lg, ...
                'the axes Legend property must reach the legend the plotters set');
        end

        function testSingleSessionCurvePlot(testCase)
            pd = localExpand(PsyRATTestDataFactory.makeICSummaryData(), ...
                testCase.Groups, testCase.Events);
            fig = psyrat_depvtrialsplot('psyrat_data', pd, ...
                'trials', [1 6], 'depcutoff', 0.8, 'depline', 2, 'CI', 0.95);
            testCase.verifyTitlesAndLegends(fig, testCase.Events, testCase.Groups);
        end

        function testTestRetestCurvePlot(testCase)
            pd = localExpand(PsyRATTestDataFactory.makeTRTSummaryData(), ...
                testCase.Groups, testCase.Events);
            fig = psyrat_trt_relvtrialsplot('psyrat_data', pd, ...
                'trials', [1 6], 'relcutoff', 0.8, 'relline', 2, 'CI', 0.95);
            testCase.verifyTitlesAndLegends(fig, testCase.Events, testCase.Groups);
        end

        function testCriterionCurvePlot(testCase)
            pd = localExpand(PsyRATTestDataFactory.makeICSummaryData(), ...
                testCase.Groups, testCase.Events);
            prefs = struct('view', struct( ...
                'criterioncutoff', 0.5, ...
                'plotcriterion', 1, ...
                'tablecriterion', 0), 'ver', '0-test');
            before = findall(groot, 'Type', 'figure');
            psyrat_criterionfigures('psyrat_prefs', prefs, ...
                'psyrat_data', pd, 'analysis', 'sing');
            figs = setdiff(findall(groot, 'Type', 'figure'), before);
            testCase.assertNotEmpty(figs, 'the criterion route drew no figure');
            testCase.verifyTitlesAndLegends(figs, testCase.Events, testCase.Groups);
        end

        function testPointIntervalPlotsTwoEvents(testCase)
            %Two groups x two events takes the errorbar-matrix branch, whose
            %legend is built afterwards from the group names.
            pd = localExpand(PsyRATTestDataFactory.makeICSummaryData(), ...
                testCase.Groups, testCase.Events);
            for stat = {'icc', 'bet', 'wit'}
                h = psyrat_ptintervalplot('psyrat_data', pd, 'stat', stat{1});
                testCase.verifyForestAxes(h(1).Parent, testCase.Events, testCase.Groups);
            end
        end

        function testPointIntervalPlotSingleEvent(testCase)
            %Two groups x ONE event takes the '-DynamicLegend' branch, where
            %the legend is created from the series DisplayNames inside the
            %loop; the fix reaches it through the axes Legend property.
            pd = localExpand(PsyRATTestDataFactory.makeICSummaryData(), ...
                testCase.Groups, testCase.Events(1));
            h = psyrat_ptintervalplot('psyrat_data', pd, 'stat', 'icc');
            testCase.verifyForestAxes(h(1).Parent, testCase.Events(1), testCase.Groups);
        end

        function testSubjectLevelPlots(testCase)
            %psyrat_ssrelplot is the production caterpillar (the viewer's
            %Subject-Level route). psyrat_depssplot carries the same title
            %sites but has no caller and cannot draw two groups (it indexes a
            %scalar struct per group, psyrat_depssplot.m:231), so it is not
            %exercised here; see the 2026-09-03 finding in the trackers.
            pd = localExpand(PsyRATTestDataFactory.makeICSummaryData(), ...
                testCase.Groups, testCase.Events);
            %built with sprintf rather than strcat, which strips the trailing
            %space from ': ' and would never match the plotter's titles
            [gi, ei] = ndgrid(1:numel(testCase.Groups), 1:numel(testCase.Events));
            expected = arrayfun(@(g, e) sprintf('%s: %s', ...
                testCase.Groups{g}, testCase.Events{e}), gi(:), ei(:), ...
                'UniformOutput', false);

            before = findall(groot, 'Type', 'figure');
            psyrat_ssrelplot('psyrat_data', pd, 'stat', 'dep');
            figs = setdiff(findall(groot, 'Type', 'figure'), before);
            testCase.assertNotEmpty(figs, 'psyrat_ssrelplot drew no figure');
            testCase.verifyTitlesAndLegends(figs, expected, {});
        end

        function testTracePlotYLabels(testCase)
            %psyrat_checktraceplots (estimation/, not plotting/) prints each
            %group x event cell's name as the y-label of its mu trace; the
            %manual's first trace-plot capture (2026-09-04) showed the same
            %subscripted rendering the plotters had. The 'ic' branch is the
            %one the manual shows; the 'trt' branch carries the identical
            %ylabel site and received the identical edit.
            REL = localBuildTraceREL(testCase.Groups, testCase.Events);
            expected = {};
            for e = 1:numel(testCase.Events)
                for g = 1:numel(testCase.Groups)
                    expected{end+1} = sprintf('%s - %s', ...
                        testCase.Events{e}, testCase.Groups{g}); %#ok<AGROW>
                end
            end
            before = findall(groot, 'Type', 'figure');
            psyrat_checktraceplots(REL);
            figs = setdiff(findall(groot, 'Type', 'figure'), before);
            testCase.addTeardown(@() close(figs));
            testCase.assertNotEmpty(figs, 'psyrat_checktraceplots drew no figure');
            axs = findall(figs, 'Type', 'axes');
            nchecked = 0;
            for k = 1:numel(axs)
                yl = cellstr(axs(k).YLabel.String);
                if isempty(yl) || isempty(yl{1}) || ~any(strcmp(yl{1}, expected))
                    continue;
                end
                nchecked = nchecked + 1;
                testCase.verifyEqual(axs(k).YLabel.Interpreter, 'none', ...
                    sprintf('trace y-label "%s" must use Interpreter none', yl{1}));
            end
            testCase.verifyEqual(nchecked, numel(expected), ...
                'expected one labelled mu trace per group x event cell');
        end

        function testDodHeatmapLabels(testCase)
            %The DoD Reliability Outputs window (psyrat_startview_dodiff in
            %guis/, drawn by its local function local_plot_dod_outputs) titles
            %its two heatmaps and labels their axes with the four event names.
            %The factory fixture's events carry underscores (ERP_A ... ERP_D).
            %Built the way a user reaches it: the setup screen, then its own
            %Generate DoD Outputs callback (owner ruling 2026-09-04).
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = PsyRATTestDataFactory.makeDoDRelData();
            if ~isfield(pd, 'ver'); pd.ver = '0-test'; end
            if ~isfield(pd, 'proc') || ~isfield(pd.proc, 'measheader')
                pd.proc.measheader = 'meas';
            end
            psyrat_startview_dodiff('psyrat_prefs', psyrat_defaults, ...
                'psyrat_data', pd);
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
            axs = findall(out(1), 'Type', 'axes');
            isheat = false(size(axs));
            for k = 1:numel(axs)
                isheat(k) = contains(char(string(axs(k).Title.String)), 'Heatmap');
            end
            heat = axs(isheat);
            testCase.assertEqual(numel(heat), 2, 'expected the two heatmap axes');
            for k = 1:numel(heat)
                t = char(string(heat(k).Title.String));
                testCase.verifyEqual(heat(k).Title.Interpreter, 'none', ...
                    sprintf('heatmap title "%s" must use Interpreter none', t));
                testCase.verifyEqual(heat(k).XLabel.Interpreter, 'none', ...
                    sprintf('x-label under "%s" must use Interpreter none', t));
                testCase.verifyEqual(heat(k).YLabel.Interpreter, 'none', ...
                    sprintf('y-label under "%s" must use Interpreter none', t));
            end
        end

        function testTestRetestDifferenceCurvePlot(testCase)
            %The contrast title carries both event names. This fixture has one
            %group, so no legend is drawn; the legend site in this plotter
            %shares its set(leg, ...) form with the two curve plots above.
            pd = localBuildTrtDiffREL(testCase.Events);
            out = psyrat_relsummary('psyrat_data', pd, 'analysis', 'trt_diff', ...
                'gcoeff', 1, 'diffgcoeff', 1, 'reltype', 3, ...
                'relcutoff', 0.5, 'meascutoff', 2, 'relcentmeas', 1, ...
                'CI', 0.95, 'noccmode', 1);
            fig = psyrat_trt_diffvtrialsplot('psyrat_data', out, ...
                'trials', [1 10], 'relcutoff', 0.5, 'relline', 2, 'CI', 0.95);
            ax = findall(fig, 'Type', 'axes');
            testCase.assertNumElements(ax, 1);
            ttl = cellstr(ax.Title.String);
            testCase.verifySubstring(ttl{1}, ...
                [testCase.Events{1} ' - ' testCase.Events{2}]);
            testCase.verifyEqual(ax.Title.Interpreter, 'none', ...
                'the difference-score contrast title must print event names literally');
        end

    end

    methods (Access = private)

        function verifyTitlesAndLegends(testCase, figs, titles, groups)
            %Every axes whose title is one of TITLES must carry Interpreter
            %'none' on that title and, when a legend exists, on the legend,
            %and the legend must list GROUPS. Requires at least one such axes
            %so a plotter that silently draws nothing cannot pass.
            axs = findall(figs, 'Type', 'axes');
            nchecked = 0;
            for k = 1:numel(axs)
                ax = axs(k);
                t = cellstr(ax.Title.String);
                if isempty(t) || ~any(strcmp(t{1}, titles))
                    continue;
                end
                nchecked = nchecked + 1;
                testCase.verifyEqual(ax.Title.Interpreter, 'none', ...
                    sprintf('title "%s" must use Interpreter none', t{1}));
                if ~isempty(ax.Legend)
                    testCase.verifyEqual(ax.Legend.Interpreter, 'none', ...
                        sprintf('legend under title "%s" must use Interpreter none', t{1}));
                    if ~isempty(groups)
                        testCase.verifyEqual(cellstr(ax.Legend.String(:)), groups(:), ...
                            'the legend must list the group names as given');
                    end
                elseif ~isempty(groups)
                    testCase.verifyFail(sprintf( ...
                        'axes titled "%s" has no legend, so the group-name site was not exercised', t{1}));
                end
            end
            testCase.verifyEqual(nchecked, numel(titles), ...
                'expected one titled axes per event/stratum name');
        end

        function verifyForestAxes(testCase, ax, events, groups)
            %The forest plots print event names as tick labels and group
            %names in the legend.
            testCase.verifyEqual(ax.TickLabelInterpreter, 'none', ...
                'forest-plot tick labels must use Interpreter none');
            testCase.verifyEqual(cellstr(ax.XTickLabel(:)), events(:), ...
                'the tick labels must be the event names as given');
            testCase.assertNotEmpty(ax.Legend, 'the forest plot must carry a legend');
            testCase.verifyEqual(ax.Legend.Interpreter, 'none', ...
                'forest-plot legend must use Interpreter none');
            testCase.verifyEqual(cellstr(ax.Legend.String(:)), groups(:), ...
                'the legend must list the group names as given');
        end

    end
end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function pd = localExpand(pd, groups, events)
%Replicate the single (group, event) cell of a factory summary fixture into
%numel(groups) x numel(events) cells that carry the given names. The plotters
%read the names from psyrat_data.rel.groups / rel.events and the per-cell
%summaries from relsummary.group(g).event(e) and relsummary.data.g(g).e(e),
%so the replicated values are what the fixture already provides; only the
%labels matter here.
pd.rel.groups = groups(:);
pd.rel.events = events(:);
g1 = pd.relsummary.group(1);
e1 = g1.event(1);
d1 = pd.relsummary.data.g(1).e(1);
for g = 1:numel(groups)
    pd.relsummary.group(g) = g1;
    pd.relsummary.group(g).name = groups{g};
    for e = 1:numel(events)
        pd.relsummary.group(g).event(e) = e1;
        pd.relsummary.group(g).event(e).name = events{e};
        pd.relsummary.data.g(g).e(e) = d1;
    end
end
end

function pd = localBuildTrtDiffREL(events)
%Synthetic psyrat_data.rel for the 'trt_diff' relsummary branch with constant
%known components, mirroring the builder in TestTrtDiffCurveWiring (kept
%local there for the same reason it is kept local here: the two files pin
%different properties of the branch). Non-concurrent (rho = 0) and one group,
%which is all the contrast title needs.
ndraw = 200;
nSub = 20;
nOcc = 4;
nTrl = 8;
rho = 0;

cov2 = @(s1, s2, r) [s1^2, r*s1*s2; r*s1*s2, s2^2];
Sp  = cov2(3.0, 3.0, 0.6);
Spt = cov2(2.0, 2.0, 0.4);
Spo = cov2(1.2, 1.2, 0.3);
St  = cov2(0.6, 0.6, 0.3);
So  = cov2(0.6, 0.6, 0.3);
Sto = cov2(0.4, 0.4, 0.2);
sigma = [2.0 2.0];

REL = struct();
REL.groups = 'none';
REL.events = events(:)';
REL.time = arrayfun(@(o) sprintf('o%d', o), 1:nOcc, 'UniformOutput', false);
REL.analysis = 'trt_diff';
REL.diffrescor = 1 + double(rho ~= 0);
REL.nchains = 2; REL.niter = 1; REL.filename = 'synthetic';

ids = {}; meas = []; event = {}; time = {};
for s = 1:nSub
    for o = 1:nOcc
        for e = 1:2
            for t = 1:nTrl %#ok<*AGROW>
                ids{end+1, 1} = sprintf('s%02d', s);
                event{end+1, 1} = REL.events{e};
                time{end+1, 1} = REL.time{o};
                meas(end+1, 1) = 0;
            end
        end
    end
end
REL.data = table(ids, meas, event, time, ...
    'VariableNames', {'id', 'meas', 'event', 'time'});

REL.out = struct();
REL.out.id_varcov  = {localConstCov(Sp,  ndraw)};
REL.out.tid_varcov = {localConstCov(Spt, ndraw)};
REL.out.oid_varcov = {localConstCov(Spo, ndraw)};
REL.out.trl_varcov = {localConstCov(St,  ndraw)};
REL.out.occ_varcov = {localConstCov(So,  ndraw)};
REL.out.to_varcov  = {localConstCov(Sto, ndraw)};
REL.out.b_sigma    = {num2cell(repmat(log(sigma), ndraw, 1))};
REL.out.err_varcov = {localConstCov(cov2(sigma(1), sigma(2), rho), ndraw)};
REL.out.rescor     = {num2cell(rho * ones(ndraw, 1))};
REL.out.elabels = {REL.events(:)};
REL.out.glabels = {{'none'}};

pd = struct();
pd.rel = REL;
pd.ver = '0-test';
end

function REL = localBuildTraceREL(groups, events)
%Minimal one-facet REL for psyrat_checktraceplots' 'ic' branch: per-chain
%draw blocks stacked down the rows (nsampling rows per chain) and one column
%per group x event cell, labelled event_;_group as psyrat_computevarcomp
%emits them (see the case-4 label-order note in the trackers).
nsampling = 5;
nchains = 2;
ncell = numel(groups) * numel(events);
REL = struct();
REL.analysis = 'ic';
REL.groups = groups(:);
REL.events = events(:);
REL.nchains = nchains;
REL.niter = 2 * nsampling;
REL.nsampling = nsampling;
REL.out = struct();
REL.out.mu = rand(nsampling * nchains, ncell);
REL.out.sig_u = rand(nsampling * nchains, ncell);
REL.out.sig_e = rand(nsampling * nchains, ncell);
labels = cell(1, ncell);
k = 0;
for e = 1:numel(events)
    for g = 1:numel(groups)
        k = k + 1;
        labels{k} = [events{e} '_;_' groups{g}];
    end
end
REL.out.labels = labels;
end

function c = localConstCov(S, ndraw)
A = repmat(reshape(S, [1 2 2]), [ndraw 1 1]);
c = num2cell(A);
end
