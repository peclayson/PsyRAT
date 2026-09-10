classdef TestTrtDiffCurveWiring < PsyRATTestBase
    % RC-31. The two-facet difference design ('trt_diff') reports TWO
    % independently selected error types: gcoeff for the per-condition output and
    % diffgcoeff for the difference score (the contract the toolbox states in its
    % own words at psyrat_trt_coeflabel.m:94-99). The reliability-vs-trials curve
    % was the one surface that broke it: the curve is per-condition -- it calls
    % psyrat_rel_trt, the plain two-facet kernel, on the per-condition components
    % -- but both the figure (psyrat_relfigures) and the headless export
    % (psyrat_report) fed it diffgcoeff, and the export then labeled it
    % "(difference score)".
    %
    % These tests pin the two halves of the fix:
    %   (a) each curve now follows the selector it names, and
    %   (b) the genuine difference-score curve that was missing is correct
    %       against the EXTERNAL oracle, not merely unchanged.
    %
    % The synthetic REL is built with KNOWN, constant-across-draws components
    % (the localBuildREL pattern from TestTrtDiffRelsummary), so every point on
    % the difference curve has a closed form and can be checked against
    % PsyRATAccuracyOracle. No CmdStan.

    methods (Test)

        function testDiffCurveMatchesTheOracleAtEveryN(testCase)
            % ACCURACY, not regression. The oracle is an independent
            % implementation, so this is an external anchor: it would catch the
            % curve being built from the wrong components, at the wrong occasion
            % estimand, or with the trial count applied to only one condition.
            [out, truth] = localSummary(testCase, 0, 3, 1, 1);
            rep = psyrat_report(out, 'ntrials', 12);

            crv = rep.tables.diffcurve;
            testCase.verifyEqual(height(crv), 12, ...
                'One row per trial count for the single group.');

            for n = [1 2 5 12]
                oracle = localOracle(truth, 3, 'dep', [n n], [1 1]);
                row = crv(crv.n_trials == n, :);
                testCase.verifyEqual(row.Dependability_PointEstimate, ...
                    oracle.pt, 'AbsTol', 1e-4, ...
                    sprintf('Difference curve at n = %d disagrees with the oracle.', n));
                testCase.verifyEqual(row.Dependability_LowerLimit, ...
                    oracle.ll, 'AbsTol', 1e-4);
                testCase.verifyEqual(row.Dependability_UpperLimit, ...
                    oracle.ul, 'AbsTol', 1e-4);
            end
        end

        function testDiffCurveHonorsTheOccasionEstimand(testCase)
            % The difference score is a two-facet quantity, so the curve must be
            % drawn at the SELECTED n'_o, not silently at 1 -- the RC-10 defect,
            % on the surface RC-10 did not reach. Checked against the oracle at
            % the composite k rather than merely "different", so a curve that
            % moved for the wrong reason still fails.
            [out, truth] = localSummary(testCase, 0, 3, 1, 1, ...
                'noccmode', 2, 'nocc', 4);
            testCase.assertEqual(out.relsummary.nocc, 4);
            rep = psyrat_report(out, 'ntrials', 8);
            crv = rep.tables.diffcurve;

            for n = [1 8]
                oracle = localOracle(truth, 3, 'dep', [n n], [4 4]);
                row = crv(crv.n_trials == n, :);
                testCase.verifyEqual(row.Dependability_PointEstimate, ...
                    oracle.pt, 'AbsTol', 1e-4);
            end

            % non-vacuity: n'_o = 4 must actually move the curve off n'_o = 1
            single = psyrat_report(localSummary(testCase, 0, 3, 1, 1), 'ntrials', 8);
            testCase.verifyGreaterThan( ...
                max(abs(crv.Dependability_PointEstimate - ...
                    single.tables.diffcurve.Dependability_PointEstimate)), 1e-3, ...
                'The composite occasion estimand must change the difference curve.');
        end

        function testEachCurveFollowsItsOwnCoefficient(testCase)
            % The core RC-31 assertion. Under the shipped default both selectors
            % are dependability, so the fix is a strict no-op there and a test
            % that does not vary them INDEPENDENTLY proves nothing. All four
            % combinations are built; the per-condition curve must track gcoeff
            % and ignore diffgcoeff, and the difference curve the reverse.
            dd = localPointCurves(testCase, 1, 1);
            dg = localPointCurves(testCase, 1, 2);   % only diffgcoeff moves
            gd = localPointCurves(testCase, 2, 1);   % only gcoeff moves
            gg = localPointCurves(testCase, 2, 2);

            % non-vacuity first: dependability and generalizability must be
            % materially different on these components, or every check below
            % passes for the wrong reason.
            testCase.assertGreaterThan(max(abs(dd.cond - gd.cond)), 1e-3, ...
                'Per-condition dep and gen are indistinguishable; test is vacuous.');
            testCase.assertGreaterThan(max(abs(dd.diff - dg.diff)), 1e-3, ...
                'Difference-score dep and gen are indistinguishable; test is vacuous.');

            % per-condition curve: follows gcoeff, invariant to diffgcoeff
            testCase.verifyEqual(dg.cond, dd.cond, 'AbsTol', 1e-12, ...
                'The per-condition curve moved with diffgcoeff (the RC-31 defect).');
            testCase.verifyEqual(gg.cond, gd.cond, 'AbsTol', 1e-12);

            % difference curve: follows diffgcoeff, invariant to gcoeff
            testCase.verifyEqual(gd.diff, dd.diff, 'AbsTol', 1e-12, ...
                'The difference-score curve moved with gcoeff.');
            testCase.verifyEqual(gg.diff, dg.diff, 'AbsTol', 1e-12);
        end

        function testCurveNoteNoLongerCallsThePerConditionSeriesADifference(testCase)
            % The exported header is the only place a standalone file records
            % which quantity the curve is. It claimed "(difference score)" for a
            % per-condition series.
            out = localSummary(testCase, 0, 3, 1, 1);
            rep = psyrat_report(out, 'ntrials', 5);

            testCase.verifySubstring(rep.curve_note, 'per condition');
            testCase.verifyFalse( ...
                contains(rep.curve_note, '(difference score)'), ...
                'The per-condition curve is still labeled a difference score.');
            testCase.verifySubstring(rep.diffcurve_note, 'Difference-score');
        end

        function testFigureAndHeadlessDiffCurveAreTheSameEstimand(testCase)
            % Figure/table parity, extending the invariant TestTrtMultiOccasion
            % already pins for the plain test-retest curve. The plot and the
            % export are two renderings of one quantity; if they drift, the user
            % reads a required trial count off a figure that the file beside it
            % contradicts.
            ntrials = 15;
            out = localSummary(testCase, 0.4, 3, 1, 2);   % concurrent, gen diff
            rep = psyrat_report(out, 'ntrials', ntrials);

            y = testCase.plotDiffCurveY(out, ntrials);
            testCase.verifyEqual(y, ...
                rep.tables.diffcurve.Generalizability_PointEstimate, ...
                'AbsTol', 1e-4, ...
                'The difference figure and its headless export disagree.');
        end

        function testCriterionRouteInputsAbsentFromTrtDiffSummary(testCase)
            % Root-cause pin for the trt viewer suppressing Criterion Score
            % Outputs on trt_diff data (G38): psyrat_criterionfigures passes
            % 'mu',d.mu.raw on every trt-family route it accepts, but the
            % trt_diff D-study stores no mu on relsummary.data.g().e(), so the
            % criterion route dies on a missing field with nothing it could
            % compute. If mu ever appears here, revisit that suppression.
            out = localSummary(testCase, 0, 3, 1, 1);
            d = out.relsummary.data.g(1).e(1);
            testCase.verifyTrue(isfield(d, 'sig_id'), ...
                'anti-vacuity: the component store must be present at this path');
            testCase.verifyFalse(isfield(d, 'mu'), ...
                ['trt_diff stores no universe-score mean, which is why the ' ...
                'trt viewer does not offer Criterion Score Outputs for it']);
        end

        function testRelfiguresDrawsBothCurvesForTrtDiff(testCase)
            % Dispatch wiring: the viewer's single "trials vs reliability"
            % checkbox must now produce BOTH the per-condition figure and the
            % difference-score figure.
            [pd, ~] = localBuildREL(testCase, 0, 4, 8);
            prefs = localViewPrefs(1, 1);
            cleanup = onCleanup(@() close('all','force'));

            close('all','force');
            psyrat_relfigures('psyrat_data', pd, 'psyrat_prefs', prefs, ...
                'analysis', 'trt_diff');

            names = localFigureNames();
            testCase.verifyTrue( ...
                any(contains(names,'Difference Score')), ...
                'psyrat_relfigures(trt_diff) did not draw the difference-score curve.');
        end

        function testRelfiguresGivesEachFigureItsOwnCoefficient(testCase)
            % THE GUI HALF OF RC-31, and the defect as a user met it. The two
            % figures psyrat_relfigures draws must be labeled by the selectors
            % they actually name: the per-condition curve by gcoeff, the
            % difference curve by diffgcoeff.
            %
            % Driven at gcoeff = 2 / diffgcoeff = 1 because that is the state a
            % user reaches WITHOUT touching either popup: prefs.view.gcoeff is
            % written on every test-retest view and persists, so opening a
            % trt_diff result after a plain trt run at Generalizability lands
            % exactly here. Before the fix the per-condition figure was drawn at
            % diffgcoeff and so came up titled "Dependability" while the
            % per-condition rows in the table beside it were generalizability.
            %
            % Asserted on the window titles rather than the plotted data because
            % the title is what the reader uses to identify the quantity, and
            % psyrat_depvtrialsplot derives it from the same gcoeff it computes
            % with -- so a title check cannot pass while the curve is wrong.
            [pd, ~] = localBuildREL(testCase, 0, 4, 8);
            prefs = localViewPrefs(2, 1);    % per-condition gen, difference dep
            cleanup = onCleanup(@() close('all','force'));

            close('all','force');
            psyrat_relfigures('psyrat_data', pd, 'psyrat_prefs', prefs, ...
                'analysis', 'trt_diff');
            names = localFigureNames();

            cond = names(~contains(names,'Difference Score') & ...
                contains(names,'v Number of Trials'));
            diff = names(contains(names,'Difference Score'));
            testCase.assertNotEmpty(cond, 'No per-condition trials figure drawn.');
            testCase.assertNotEmpty(diff, 'No difference-score trials figure drawn.');

            testCase.verifyTrue(all(contains(cond,'Generalizability')), ...
                sprintf(['The per-condition curve must follow gcoeff (= 2, ', ...
                'generalizability). Title was: %s'], cond{1}));
            testCase.verifyTrue(all(contains(diff,'Dependability')), ...
                sprintf(['The difference-score curve must follow diffgcoeff ', ...
                '(= 1, dependability). Title was: %s'], diff{1}));

            % and the mirror image, so neither assertion can pass by both
            % figures happening to carry the same label
            close('all','force');
            psyrat_relfigures('psyrat_data', pd, 'psyrat_prefs', ...
                localViewPrefs(1, 2), 'analysis', 'trt_diff');
            names = localFigureNames();
            cond = names(~contains(names,'Difference Score') & ...
                contains(names,'v Number of Trials'));
            diff = names(contains(names,'Difference Score'));
            testCase.verifyTrue(all(contains(cond,'Dependability')));
            testCase.verifyTrue(all(contains(diff,'Generalizability')));
        end

        function testDiffCurveWarnsOnceNotOncePerTrialCount(testCase)
            % The kernel warns per CALL and the curve calls it once per trial
            % count. At obs = [n n] every adjusted block reduces to raw/n, so an
            % inadmissible input is inadmissible at EVERY point -- without
            % aggregation this emits ntrials byte-identical blocks.
            %
            % Reaching it needs a block that is not positive semi-definite (a
            % real covariance matrix cannot make a raw variance of a difference
            % negative), which is the gamma-family failure mode of RC-23: the
            % between-person covariance below exceeds sqrt(v1*v2).
            ntrials = 20;
            out = localSummary(testCase, 0, 3, 1, 1);
            nd = size(out.relsummary.group(1).diffcomp.er_var, 1);
            one = ones(nd,1);
            out.relsummary.group(1).diffcomp.bp = ...
                cat(3, [1*one 1.5*one], [1.5*one 1*one]);   % r > 1, not PSD

            state = warning('on','psyrat:negerrvar');
            restore = onCleanup(@() warning(state));

            % The banner is the first line psyrat_admissibility_note emits, one
            % per warning block, so counting it counts blocks. Matched on the
            % literal text the note actually prints -- an earlier version of this
            % test matched a lowercase variant that never appears, so it counted
            % zero and passed no matter how many blocks were emitted.
            txt = evalc('psyrat_report(out, ''ntrials'', ntrials);');
            nblocks = numel(strfind(txt, 'NOT INTERPRETABLE AS A RELIABILITY'));

            % Non-vacuity: the mutated components must actually be inadmissible,
            % or "few blocks" is trivially true. Asserted, not verified, because
            % every claim below depends on it.
            testCase.assertGreaterThanOrEqual(nblocks, 1, ...
                ['The non-PSD components did not trip the admissibility guard, ' ...
                'so this test cannot detect a missing aggregation.']);

            testCase.verifyEqual(nblocks, 1, ...
                sprintf(['The difference curve emitted %d admissibility blocks ' ...
                'for a %d-point curve. At obs = [n n] an inadmissible input is ' ...
                'inadmissible at every n, so without aggregation this is one ' ...
                'block per trial count.'], nblocks, ntrials));
        end

    end

    methods (Access = private)

        function y = plotDiffCurveY(~, psyrat_data, ntrials)
            %Draw the difference-score figure and read the plotted curve back.
            %The figure also carries a flat cutoff reference line, so the curve
            %is identified by its length rather than by draw order -- the same
            %approach TestTrtMultiOccasion uses for the per-condition curve.
            f = psyrat_trt_diffvtrialsplot('psyrat_data', psyrat_data, ...
                'trials', [1 ntrials], 'relcutoff', 0.5, 'relline', 2, ...
                'CI', 0.95);
            closer = onCleanup(@() close(f));
            y = [];
            lines = findall(f, 'Type', 'line');
            for k = 1:numel(lines)
                yy = get(lines(k), 'YData');
                if numel(yy) == ntrials
                    y = yy(:);
                end
            end
            assert(~isempty(y), 'Could not find the plotted difference curve.');
        end

    end
end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function [out, truth] = localSummary(testCase, rho, reltype, gcoeff, diffgcoeff, varargin)
%Build the synthetic REL and run the real psyrat_relsummary trt_diff branch.
[pd, truth] = localBuildREL(testCase, rho, 4, 8);
args = {'psyrat_data', pd, 'analysis', 'trt_diff', ...
    'gcoeff', gcoeff, 'diffgcoeff', diffgcoeff, 'reltype', reltype, ...
    'relcutoff', 0.5, 'meascutoff', 2, 'relcentmeas', 1, 'CI', 0.95};
if isempty(varargin)
    args = [args, {'noccmode', 1}];
else
    args = [args, varargin];
end
out = psyrat_relsummary(args{:});
end

function c = localPointCurves(testCase, gcoeff, diffgcoeff)
%Point-estimate columns of both curves for one (gcoeff, diffgcoeff) pair.
%Extracted by column position rather than by name because the column NAMES
%encode the coefficient ('Dependability_*' vs 'Generalizability_*') and so
%change with the very selectors under test.
out = localSummary(testCase, 0, 3, gcoeff, diffgcoeff);
rep = psyrat_report(out, 'ntrials', 10);
c = struct();
c.cond = rep.tables.curve{:, 5};        % Group, Event, n_trials, LL, PT, UL
c.diff = rep.tables.diffcurve{:, 4};    % Group, n_trials, LL, PT, UL
end

function v = localViewPrefsStruct(gcoeff, diffgcoeff)
%The two coefficient selectors are parameters, not constants: RC-31 is only
%observable when they DIFFER, and under the shipped default (both dependability)
%the whole fix is a strict no-op.
v = struct();
v.gcoeff = gcoeff; v.reltype = 3; v.diffgcoeff = diffgcoeff;
v.noccmode = 1; v.nocc = [];
v.relvalue = 0.5; v.relcentmeas = 1; v.meascutoff = 2;
v.plotrel = 1; v.plotrelline = 2;
%only the trials curves are wanted; the other outputs open extra windows this
%test would then have to filter out of the figure-name list.
v.ploticc = 0; v.inctrltable = 0; v.overalltable = 0;
v.ntrials = 20; v.showstddevt = 0; v.showstddevf = 0;
v.plotssrel = 0; v.plotssicc = 0; v.tablessrel = 0;
end

function prefs = localViewPrefs(gcoeff, diffgcoeff)
prefs = struct('view', localViewPrefsStruct(gcoeff, diffgcoeff), 'ver', '0-test');
end

function names = localFigureNames()
%Names of every currently open figure, as a cellstr (get returns a bare char
%when exactly one figure is open, which breaks contains(...) over a cell).
names = get(findall(groot,'Type','figure'), 'Name');
if ~iscell(names)
    names = {names};
end
names = cellfun(@(s) string(s), names, 'UniformOutput', false);
names = cellfun(@char, names, 'UniformOutput', false);
end

function [pd, truth] = localBuildREL(testCase, rho, nOcc, nTrl)
% Synthetic psyrat_data.rel for the 'trt_diff' path with KNOWN,
% constant-across-draws components, so every curve point has a closed form.
% Mirrors the builder in TestTrtDiffRelsummary; kept local because the two files
% pin different properties of the same branch and a shared fixture would couple
% them.
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
