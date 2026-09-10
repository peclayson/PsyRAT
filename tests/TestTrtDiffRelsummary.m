classdef TestTrtDiffRelsummary < PsyRATTestBase
    % Relsummary-level wiring test for the two-facet (test-retest)
    % difference-score analysis ('trt_diff'). Builds a synthetic REL with KNOWN,
    % constant-across-draws cross-condition components (the structure
    % psyrat_computevarcomp case 10 produces) and runs the real psyrat_relsummary
    % path. Because draws are constant, the difference-score point estimate is the
    % closed form, so it must equal PsyRATAccuracyOracle.diffrelTrt exactly. This
    % guards the wiring (parse case + summary branch + helper), not just the
    % underlying psyrat_diffrel_trt function -- the F1/F7 lesson. No CmdStan.

    methods (Test)
        function testNonConcurrentSingleOccasionCES(testCase)
            [pd, truth] = localBuildREL(testCase, 0, 4, 8);  % rescor 0
            out = psyrat_relsummary('psyrat_data', pd, 'analysis', 'trt_diff', ...
                'gcoeff', 1, 'diffgcoeff', 1, 'reltype', 3, ...
                'noccmode', 1, 'relcutoff', 0.5, 'meascutoff', 2, ...
                'relcentmeas', 1, 'CI', 0.95);
            rs = out.relsummary;

            ds = rs.group(1).diffscore;
            oracle = localOracle(truth, 3, 'dep', [8 8], [1 1]);
            testCase.verifyLessThan(abs(ds.pt - oracle.pt), 1e-9, ...
                sprintf('CES dep wiring mismatch: relsummary %.6f vs oracle %.6f', ...
                ds.pt, oracle.pt));
            testCase.verifyGreaterThan(ds.pt, 0);
            testCase.verifyLessThan(ds.pt, 1);

            % single-occasion estimand recorded
            testCase.verifyEqual(rs.nocc, 1);
            testCase.verifyEqual(rs.nocc_name, 'single');

            % per-event scaffolding the difference tables consume exists
            for eloc = 1:2
                ev = rs.group(1).event(eloc);
                testCase.verifyTrue(isfield(ev, 'trlinfo') && isfinite(ev.trlinfo.mean));
                testCase.verifyTrue(isfield(ev, 'betsd') && isfinite(ev.betsd.m));
                testCase.verifyTrue(isfield(ev, 'witsd') && isfinite(ev.witsd.m));
                testCase.verifyTrue(isfield(ev, 'icc') && isfinite(ev.icc.m));
                testCase.verifyTrue(isfield(ev, 'goodn'));
            end
            testCase.verifyTrue(isfield(ds, 'trlcutoff'));
        end

        function testConcurrentSingleOccasionCEgen(testCase)
            [pd, truth] = localBuildREL(testCase, 0.4, 4, 8);  % rescor 0.4
            out = psyrat_relsummary('psyrat_data', pd, 'analysis', 'trt_diff', ...
                'gcoeff', 2, 'diffgcoeff', 2, 'reltype', 1, ...
                'noccmode', 1, 'relcutoff', 0.5, 'meascutoff', 2, ...
                'relcentmeas', 1, 'CI', 0.95);
            ds = out.relsummary.group(1).diffscore;
            oracle = localOracle(truth, 1, 'gen', [8 8], [1 1]);
            testCase.verifyLessThan(abs(ds.pt - oracle.pt), 1e-9, ...
                sprintf('CE gen (concurrent) wiring mismatch: %.6f vs %.6f', ...
                ds.pt, oracle.pt));
            % residual covariance should be carried through (non-zero)
            testCase.verifyGreaterThan(ds.wp_cov, 0);
        end

        function testMultiOccasionCompositeCSdep(testCase)
            [pd, truth] = localBuildREL(testCase, 0, 4, 8);
            out = psyrat_relsummary('psyrat_data', pd, 'analysis', 'trt_diff', ...
                'gcoeff', 1, 'diffgcoeff', 1, 'reltype', 2, ...
                'noccmode', 2, 'nocc', 4, 'relcutoff', 0.5, 'meascutoff', 2, ...
                'relcentmeas', 1, 'CI', 0.95);
            rs = out.relsummary;
            ds = rs.group(1).diffscore;
            oracle = localOracle(truth, 2, 'dep', [8 8], [4 4]);
            testCase.verifyLessThan(abs(ds.pt - oracle.pt), 1e-9, ...
                sprintf('CS dep (multi-occasion k=4) wiring mismatch: %.6f vs %.6f', ...
                ds.pt, oracle.pt));
            testCase.verifyEqual(rs.nocc, 4);
            testCase.verifyEqual(rs.nocc_name, 'multi');
        end

        function testOutputTablesAndPlotRenderHeadlessly(testCase)
            % Exercise the full stage-(b) output wiring headlessly: the
            % difference tables (gui=0 returns the table) and the trials-vs-
            % reliability plot must build for 'trt_diff' without error and carry
            % the difference-score row / two-facet columns.
            [pd, ~] = localBuildREL(testCase, 0, 4, 8);
            out = psyrat_relsummary('psyrat_data', pd, 'analysis', 'trt_diff', ...
                'gcoeff', 1, 'diffgcoeff', 1, 'reltype', 3, ...
                'noccmode', 1, 'relcutoff', 0.5, 'meascutoff', 2, ...
                'relcentmeas', 1, 'CI', 0.95);

            cutoffTbl = psyrat_depcutofft('psyrat_data', out, 'gui', 0);
            overallTbl = psyrat_depoverallt('psyrat_data', out, 'gui', 0);
            varianceTbl = psyrat_variancet('psyrat_data', out, 'gui', 0);

            testCase.verifyGreaterThanOrEqual(height(cutoffTbl), 1);
            testCase.verifyGreaterThanOrEqual(height(overallTbl), 1);
            testCase.verifyGreaterThanOrEqual(height(varianceTbl), 1);

            % difference-score row present in the cutoff + overall tables
            testCase.verifyTrue(any(contains(string(overallTbl.Label),'diff score')));
            testCase.verifyTrue(any(contains(string(cutoffTbl.Label),'diff score')));

            % variance table uses the difference (two-facet) column layout.
            % ICC and SEM carry the decision type (RC-11): their value changes
            % with the coefficient selector, so the header must say which one.
            % This run is diffgcoeff = 1 (dependability, absolute error).
            testCase.verifyEqual(varianceTbl.Properties.VariableNames, ...
                {'Label','Between_StdDev','Within_StdDev','SEM_Dep','ICC_Dep'});

            fig = psyrat_depvtrialsplot('psyrat_data', out, 'trials', [1 20], ...
                'depline', 2, 'gcoeff', 1, 'depcutoff', 0.5);
            testCase.verifyTrue(isgraphics(fig,'figure'));
            close(fig);
        end

        function testWorkflowRunsMapsTrtDiff(testCase)
            runs = psyrat_workflow_runs('trt_diff');
            testCase.verifyEqual(numel(runs), 1);
            testCase.verifyEqual(runs{1}.analysis, 'trt_diff');
            testCase.verifyEqual(runs{1}.diffgcoeff, [1 2]);
        end

        function testStartviewTrtShowsDiffgcoeffControlOnlyForTrtDiff(testCase)
            % The test-retest viewer must expose the difference-score coefficient
            % selector for 'trt_diff' and hide it for plain 'trt'. Built off-screen.
            prefs = localStartviewPrefs();
            rel = struct('analysis','trt_diff','filename','synthetic', ...
                'time',{{'o1','o2'}},'groups','none','events',{{'A','B'}}, ...
                'nchains',2,'niter',1);
            pd = struct('rel',rel,'proc',struct('measheader','meas'),'ver','0-test');
            cleanup = onCleanup(@() close('all','force')); %#ok<NASGU>

            pd.rel.analysis = 'trt_diff';
            close('all','force');
            psyrat_startview_trt('psyrat_prefs',prefs,'psyrat_data',pd);
            testCase.verifyTrue(localHasDiffControl(), ...
                'Difference-Score Coefficient control missing for trt_diff.');

            pd.rel.analysis = 'trt';
            close('all','force');
            psyrat_startview_trt('psyrat_prefs',prefs,'psyrat_data',pd);
            testCase.verifyFalse(localHasDiffControl(), ...
                'Difference-Score Coefficient control should not appear for plain trt.');
        end

        function testRelfiguresDispatchRendersTrtDiff(testCase)
            % Exercise the GUI entry point (psyrat_relfigures) end-to-end for
            % 'trt_diff': prefs extraction -> relsummary -> difference tables +
            % trials-vs-reliability plot, all off-screen. Validates the dispatch
            % wiring (not just the individual render functions).
            [pd, ~] = localBuildREL(testCase, 0.4, 4, 8);   % concurrent
            prefs = localViewPrefs();
            cleanup = onCleanup(@() close('all','force')); %#ok<NASGU>

            nBefore = numel(findall(groot,'Type','figure'));
            psyrat_relfigures('psyrat_data', pd, 'psyrat_prefs', prefs, ...
                'analysis', 'trt_diff');
            nAfter = numel(findall(groot,'Type','figure'));

            testCase.verifyGreaterThan(nAfter, nBefore, ...
                'psyrat_relfigures(trt_diff) did not produce any output figures/tables.');
        end

        function testTrtDiffIccIsCompositeInvariant(testCase)
            % Per-condition ICC is a single-observation quantity, so it must be
            % identical regardless of the multi-occasion composite estimand
            % (noccmode). Regression for the code-review blocker where the
            % trt_diff per-condition ICC used nocc=effnocc instead of nocc=1: a
            % noccmode=2 run divided the occasion variance by k, so the same
            % condition's ICC diverged from a plain 'trt' run and was no longer
            % single-observation.
            [pd, ~] = localBuildREL(testCase, 0, 4, 8);   % 4 occasions, balanced
            common = {'psyrat_data', pd, 'analysis', 'trt_diff', 'gcoeff', 1, ...
                'diffgcoeff', 1, 'reltype', 3, 'relcutoff', 0.5, ...
                'meascutoff', 2, 'relcentmeas', 1, 'CI', 0.95};
            single = psyrat_relsummary(common{:}, 'noccmode', 1);
            multi  = psyrat_relsummary(common{:}, 'noccmode', 2, 'nocc', 4);

            % confirm the multi-occasion estimand actually took effect (guards a
            % vacuous pass); only the difference coefficient scales with it.
            testCase.verifyEqual(single.relsummary.nocc, 1);
            testCase.verifyEqual(multi.relsummary.nocc, 4);

            for eloc = 1:2
                s = single.relsummary.group(1).event(eloc).icc;
                m = multi.relsummary.group(1).event(eloc).icc;
                testCase.verifyEqual(m.m,  s.m,  'AbsTol', 1e-12, ...
                    sprintf('event %d ICC point differs across noccmode', eloc));
                testCase.verifyEqual(m.ll, s.ll, 'AbsTol', 1e-12);
                testCase.verifyEqual(m.ul, s.ul, 'AbsTol', 1e-12);
                testCase.verifyGreaterThan(s.m, 0);
                testCase.verifyLessThan(s.m, 1);
            end
        end

        function testTrtDiffHandlesUnbalancedOccasionTrialCounts(testCase)
            % Unbalanced occasion trial counts make the per-participant mean
            % trial count fractional; the reliability-vs-trials range bound must
            % stay an integer. Regression for the code-review blocker where
            % ntrials = max(percell)+1000 was fractional, so psyrat_rel_trt built
            % zeros(1,1007.5) and errored (MATLAB:NonIntegerInput) on the common
            % unbalanced-occasion dataset.
            [pd, ~] = localBuildREL(testCase, 0, 2, 8);   % 2 occasions, 8 trials
            d = pd.rel.data;
            % drop the last 'o2' trial of every (id,event) -> o1=8, o2=7 so the
            % per-participant mean trial count is 7.5 (fractional).
            drop = false(height(d), 1);
            ids = unique(d.id, 'stable');
            for k = 1:numel(ids)
                for e = {'A','B'}
                    rows = find(strcmp(d.id, ids{k}) & ...
                        strcmp(d.time, 'o2') & strcmp(d.event, e{1}));
                    drop(rows(end)) = true;
                end
            end
            pd.rel.data = d(~drop, :);

            out = psyrat_relsummary('psyrat_data', pd, 'analysis', 'trt_diff', ...
                'gcoeff', 1, 'diffgcoeff', 1, 'reltype', 3, 'noccmode', 1, ...
                'relcutoff', 0.5, 'meascutoff', 2, 'relcentmeas', 1, 'CI', 0.95);
            rs = out.relsummary;
            for eloc = 1:2
                ev = rs.group(1).event(eloc);
                % fractional mean confirms the unbalanced setup took effect
                testCase.verifyEqual(ev.trlinfo.mean, 7.5, 'AbsTol', 1e-9);
                % the located cutoff is an integer index or the -1 sentinel
                testCase.verifyTrue(ev.trlcutoff == -1 || ...
                    (ev.trlcutoff > 0 && ev.trlcutoff == round(ev.trlcutoff)));
            end
        end
    end
end

function prefs = localStartviewPrefs()
% A complete psyrat_prefs (guis + view) for building the test-retest viewer.
v = localViewPrefsStruct();
prefs = struct('guis', struct('fsize', 10), 'view', v, 'ver', '0-test');
end

function tf = localHasDiffControl()
% True when a 'Difference-Score Coefficient' text label is present among the
% currently open figures' uicontrols.
txts = findall(groot, 'Style', 'text');
if isempty(txts)
    tf = false;
    return;
end
labels = cellstr(string(get(txts, 'String')));
tf = any(cellfun(@(s) contains(s, 'Difference-Score Coefficient'), labels));
end

function v = localViewPrefsStruct()
v = struct();
v.gcoeff = 1; v.reltype = 3; v.diffgcoeff = 1;
v.noccmode = 1; v.nocc = [];
v.relvalue = 0.5; v.relcentmeas = 1; v.meascutoff = 2;
v.plotrel = 1; v.plotrelline = 2;
v.ploticc = 1; v.inctrltable = 1; v.overalltable = 1;
v.ntrials = 20; v.showstddevt = 1; v.showstddevf = 1;
v.plotssrel = 0; v.plotssicc = 0; v.tablessrel = 0;
end

function prefs = localViewPrefs()
% A complete psyrat_prefs.view for the test-retest/difference view flow.
prefs = struct('view', localViewPrefsStruct(), 'ver', '0-test');
end

function [pd, truth] = localBuildREL(testCase, rho, nOcc, nTrl) %#ok<INUSD>
% Build a synthetic psyrat_data.rel for the 'trt_diff' path with KNOWN,
% constant-across-draws components. nTrl trials per (id, occasion, event) so the
% per-occasion trial count resolves deterministically to nTrl.
ndraw = 200;
nSub = 20;

cov2 = @(s1,s2,r) [s1^2, r*s1*s2; r*s1*s2, s2^2];
truth = struct();
truth.Sp  = cov2(3.0, 3.0, 0.6);
truth.Spt = cov2(2.0, 2.0, 0.4);
truth.Spo = cov2(1.2, 1.2, 0.3);
truth.St  = cov2(0.6, 0.6, 0.3);
truth.So  = cov2(0.6, 0.6, 0.3);
truth.Sto = cov2(0.4, 0.4, 0.2);
truth.sigma = [2.0 2.0];
truth.rho = rho;

REL = struct();
REL.groups = 'none';
REL.events = {'A','B'};
REL.time = arrayfun(@(o) sprintf('o%d',o), 1:nOcc, 'UniformOutput', false);
REL.analysis = 'trt_diff';
REL.diffrescor = 1 + double(rho ~= 0);   % 1 non-concurrent, 2 concurrent
REL.nchains = 2; REL.niter = 1; REL.filename = 'synthetic';

% balanced long data table for trial-count resolution
ids = {}; meas = []; event = {}; time = {};
for s = 1:nSub
    for o = 1:nOcc
        for e = 1:2
            for t = 1:nTrl %#ok<*AGROW>
                ids{end+1,1} = sprintf('s%02d', s);
                event{end+1,1} = REL.events{e};
                time{end+1,1} = REL.time{o};
                meas(end+1,1) = 0;
            end
        end
    end
end
REL.data = table(ids, meas, event, time, ...
    'VariableNames', {'id','meas','event','time'});

% constant-across-draws component cells (num2cell of [ndraw x 2 x 2])
REL.out = struct();
REL.out.id_varcov  = {localConstCov(truth.Sp,  ndraw)};
REL.out.tid_varcov = {localConstCov(truth.Spt, ndraw)};
REL.out.oid_varcov = {localConstCov(truth.Spo, ndraw)};
REL.out.trl_varcov = {localConstCov(truth.St,  ndraw)};
REL.out.occ_varcov = {localConstCov(truth.So,  ndraw)};
REL.out.to_varcov  = {localConstCov(truth.Sto, ndraw)};
REL.out.b_sigma    = {num2cell(repmat(log(truth.sigma), ndraw, 1))};

errc = cov2(truth.sigma(1), truth.sigma(2), rho);
REL.out.err_varcov = {localConstCov(errc, ndraw)};
REL.out.rescor     = {num2cell(repmat(rho, ndraw, 1))};
REL.out.elabels = {REL.events(:)};
REL.out.glabels = {{'none'}};

pd = struct();
pd.rel = REL;
pd.ver = '0-test';
testCase.assertEqual(REL.diffrescor, 1 + double(rho ~= 0));
end

function c = localConstCov(S, ndraw)
A = repmat(reshape(S, [1 2 2]), [ndraw 1 1]);
c = num2cell(A);
end

function out = localOracle(truth, reltype, est, obs, nocc)
errc = truth.rho * truth.sigma(1) * truth.sigma(2);
out = PsyRATAccuracyOracle.diffrelTrt( ...
    reshape(truth.Sp,[1 2 2]),  reshape(truth.Spt,[1 2 2]), ...
    reshape(truth.Spo,[1 2 2]), reshape(truth.St,[1 2 2]), ...
    reshape(truth.So,[1 2 2]),  reshape(truth.Sto,[1 2 2]), ...
    log(truth.sigma), errc, obs, nocc, reltype, est, 0.95);
end
