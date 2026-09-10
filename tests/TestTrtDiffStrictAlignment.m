classdef TestTrtDiffStrictAlignment < PsyRATTestBase
    %G45 (owner ruling 2026-08-29: FULL alignment with the G39 strict-fold
    %conventions). The trt_diff branch was strict BEFORE G39, but its
    %strictness had a different face: an empty both-events intersection set
    %the run-global relerr.nogooddata flag, which raised the modal threshold
    %dialog and made every caller abort with no output -- one empty group
    %blocked every group -- and its per-event goodn stayed
    %numel(eventgoodids) while psyrat_depoverallt pairs goodn with the
    %GROUP-level badids count, so printed rows did not complement under
    %partial retention. Both diverged from what G39 standardized everywhere
    %else, and the route is GUI-reachable since G34 un-gated analysis 10.
    %
    %Post-G45: the empty intersection warns on the console per group and the
    %summary is returned in full (nogooddata stays 0), and goodn is rekeyed
    %to the group carrier on both event rows, exactly like the three case-4
    %rekeys G39 applied. The per-event carriers keep the informative
    %per-event partition.
    %
    %The fixture stages retention through per-id trial counts against a
    %component-derived cutoff (percell >= trlcutoff): event A gives everyone
    %rich counts; event B gives rich counts only to s01-s05 (partial) or to
    %nobody (zero). Every staged state is welded by anti-vacuity assertions
    %on the per-event carriers, so a drifted cutoff fails loudly instead of
    %passing against the wrong state.

    methods (Test)

        function testPartialRetentionRowsComplement(testCase)
            %The goodn rekey, at the state where it is observable: event A
            %retains everyone, event B retains 5, so pre-fix row A printed
            %n Included 20 beside n Excluded 15 (badids holds the 15 outside
            %the intersection) -- 35 participants in a 20-participant group.
            pd = localBuildREL(2, 24, 24, 8, 5);
            [out, relerr] = localSummary(pd, 0.80);

            g = out.relsummary.group(1);
            %anti-vacuity welds: the staged per-event partition is real
            testCase.assertEqual(numel(g.event(1).eventgoodids), 20, ...
                'weld: event A must retain everyone at this cutoff');
            testCase.assertEqual(numel(g.event(2).eventgoodids), 5, ...
                'weld: event B must retain exactly the five rich ids');

            testCase.verifyEqual(relerr.nogooddata, 0);
            testCase.verifyEqual(sort(g.goodids), ...
                sort(arrayfun(@(s) sprintf('s%02d', s), 1:5, ...
                'UniformOutput', false)'), ...
                'the group carrier must be the both-events intersection');
            for eloc = 1:2
                testCase.verifyEqual(g.event(eloc).goodn, 5, sprintf( ...
                    ['event %d: goodn must be the GROUP-level intersection ' ...
                    'count, not the per-event count'], eloc));
            end
            testCase.verifyEqual(numel(g.badids), 15, ...
                'badids must complement the intersection over the roster');
        end

        function testOverallTableComplementsOnPartialRetention(testCase)
            %The user-facing half of the rekey: psyrat_depoverallt pairs
            %goodn with numel(badids), so every rendered row must complement
            %to the 20-participant roster.
            pd = localBuildREL(2, 24, 24, 8, 5);
            out = localSummary(pd, 0.80);

            %anti-vacuity welds on the per-event carriers, same as the
            %sibling tests: the staged partial state is real, so the 5s
            %below are the intersection count and not a coincidence
            g = out.relsummary.group(1);
            testCase.assertEqual(numel(g.event(1).eventgoodids), 20, ...
                'weld: event A must retain everyone at this cutoff');
            testCase.assertEqual(numel(g.event(2).eventgoodids), 5, ...
                'weld: event B must retain exactly the five rich ids');

            tbl = psyrat_depoverallt('psyrat_data', out, 'gui', 0);
            vn = tbl.Properties.VariableNames;
            %rows: event A, event B, diff score (case 3: one group, 2 events)
            testCase.assertEqual(height(tbl), 3);
            for r = 1:2
                testCase.verifyEqual(tbl.(vn{2}){r}, 5, sprintf( ...
                    'row %d: n Included must be the intersection count', r));
                testCase.verifyEqual(tbl.(vn{3}){r}, 15, sprintf( ...
                    'row %d: n Excluded must complement to the roster', r));
            end
            testCase.verifyEqual(tbl.(vn{2}){3}, 5, ...
                'the diff-score row carries the same group-level count');
        end

        function testEmptyIntersectionContinuesZeroed(testCase)
            %The abort replacement: event B retains nobody (its counts sit
            %below the cutoff), so the intersection is empty. Pre-G45 this
            %set nogooddata and every caller aborted with no output;
            %post-G45 the run continues with a zeroed, complementing group,
            %a per-group console notice, and informative per-event rows.
            pd = localBuildREL(2, 24, 8, 8, 0);
            txt = evalc('[out, relerr] = localSummary(pd, 0.80);');

            g = out.relsummary.group(1);
            testCase.assertEqual(numel(g.event(1).eventgoodids), 20, ...
                'weld: event A must retain everyone at this cutoff');
            testCase.assertTrue(iscell(g.event(2).eventgoodids) && ...
                isempty(g.event(2).eventgoodids), ...
                'weld: event B must retain nobody at this cutoff');

            testCase.verifyEqual(relerr.nogooddata, 0, ...
                'an empty intersection must no longer set the abort flag');
            testCase.verifyTrue(iscell(g.goodids) && isempty(g.goodids));
            testCase.verifyEqual(numel(g.badids), 20, ...
                'the zeroed group must complement over the full roster');
            testCase.verifyEqual(g.event(1).goodn, 0);
            testCase.verifyEqual(g.event(2).goodn, 0);
            testCase.verifyGreaterThan(g.event(1).relcutoff.m, 0, ...
                'the surviving event''s cutoff row must stay informative');
            testCase.verifyTrue(contains(txt, 'both events'), ...
                'the per-group console notice must name the condition');
        end

        function testRelfiguresDrawsOutputForAZeroedGroup(testCase)
            %The caller half of the abort replacement: psyrat_relfigures
            %returns on nogooddata ahead of every plot/table block, so
            %pre-G45 a zeroed group produced nothing. Post-G45 the flag
            %stays 0 and the outputs are drawn.
            pd = localBuildREL(2, 24, 8, 8, 0);
            out = localSummary(pd, 0.80);

            prefs = struct('view', localViewPrefsStruct(), 'ver', '0-test');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            close('all', 'force');
            nBefore = numel(findall(groot, 'Type', 'figure'));
            psyrat_relfigures('psyrat_data', out, 'psyrat_prefs', prefs, ...
                'analysis', 'trt_diff');
            nAfter = numel(findall(groot, 'Type', 'figure'));
            testCase.verifyGreaterThan(nAfter, nBefore, ...
                'a zeroed trt_diff group must still produce output figures');
        end

        function testRetentionIsMonotoneInTheCutoff(testCase)
            %The G39 monotonicity property, held on this branch too: a lax
            %cutoff retains everyone, the staged cutoff retains the
            %intersection, and a count-unreachable cutoff retains nobody.
            %
            %Honesty note (adversarial wave, 2026-08-29): this is a
            %REGRESSION pin of a property trt_diff already had, not a pin of
            %the G45 change -- the intersection/goodids assignments it reads
            %are byte-identical pre- and post-G45, so it passes against the
            %pre-change code. It exists so the property cannot be lost while
            %this branch is being aligned with the strict-fold conventions;
            %the change-detecting pins are the four tests above.
            pd = localBuildREL(2, 24, 24, 8, 5);
            nLax = numel(localSummary(pd, 0.70).relsummary.group(1).goodids);
            nMid = numel(localSummary(pd, 0.80).relsummary.group(1).goodids);
            pdZero = localBuildREL(2, 24, 8, 8, 0);
            nZero = numel(localSummary(pdZero, 0.80).relsummary.group(1).goodids);

            testCase.assertEqual(nLax, 20, ...
                'anti-vacuity: the lax cutoff must retain the full roster');
            testCase.assertEqual(nMid, 5, ...
                'anti-vacuity: the staged cutoff must retain the intersection');
            testCase.verifyEqual(nZero, 0);
            testCase.verifyGreaterThanOrEqual(nLax, nMid);
            testCase.verifyGreaterThanOrEqual(nMid, nZero);
        end

    end
end

%--------------------------------------------------------------------------
function [out, relerr] = localSummary(pd, relcutoff)
[out, relerr] = psyrat_relsummary( ...
    'psyrat_data', pd, ...
    'analysis', 'trt_diff', ...
    'gcoeff', 1, ...
    'diffgcoeff', 1, ...
    'reltype', 3, ...
    'relcutoff', relcutoff, ...
    'meascutoff', 2, ...
    'relcentmeas', 1, ...
    'noccmode', 1);
end

function v = localViewPrefsStruct()
%One plot plus the overall table: the two surfaces the pre-G45 abort blanked.
v = struct();
v.gcoeff = 1; v.reltype = 3; v.diffgcoeff = 1;
v.noccmode = 1; v.nocc = [];
v.relvalue = 0.80; v.relcentmeas = 1; v.meascutoff = 2;
v.plotrel = 1; v.plotrelline = 2;
v.ploticc = 0; v.inctrltable = 0; v.overalltable = 1;
v.ntrials = 20; v.showstddevt = 0; v.showstddevf = 0;
v.plotssrel = 0; v.plotssicc = 0; v.tablessrel = 0;
end

function pd = localBuildREL(nOcc, nTrlA, nTrlBRich, nTrlBPoor, nRich)
%Synthetic trt_diff REL (constant-across-draws components, the
%TestTrtDiffRelsummary pattern) with PER-ID trial-count knobs: event A gives
%every participant nTrlA trials per occasion; event B gives the first nRich
%participants nTrlBRich and the rest nTrlBPoor. Retention keys on
%percell >= trlcutoff, so the knobs stage full/partial/zero intersections
%against one component-derived cutoff.
ndraw = 200;
nSub = 20;

cov2 = @(s1,s2,r) [s1^2, r*s1*s2; r*s1*s2, s2^2];
Sp  = cov2(3.0, 3.0, 0.6);
Spt = cov2(2.0, 2.0, 0.4);
Spo = cov2(1.2, 1.2, 0.3);
St  = cov2(0.6, 0.6, 0.3);
So  = cov2(0.6, 0.6, 0.3);
Sto = cov2(0.4, 0.4, 0.2);
sigma = [2.0 2.0];

REL = struct();
REL.groups = 'none';
REL.events = {'A','B'};
REL.time = arrayfun(@(o) sprintf('o%d', o), 1:nOcc, 'UniformOutput', false);
REL.analysis = 'trt_diff';
REL.diffrescor = 1;   %non-concurrent
REL.nchains = 2; REL.niter = 1; REL.filename = 'synthetic';

ids = {}; meas = []; event = {}; time = {};
for s = 1:nSub
    for o = 1:nOcc
        for e = 1:2
            if e == 1
                n = nTrlA;
            elseif s <= nRich
                n = nTrlBRich;
            else
                n = nTrlBPoor;
            end
            for t = 1:n %#ok<*AGROW>
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
REL.out.id_varcov  = {localConstCov(Sp,  ndraw)};
REL.out.tid_varcov = {localConstCov(Spt, ndraw)};
REL.out.oid_varcov = {localConstCov(Spo, ndraw)};
REL.out.trl_varcov = {localConstCov(St,  ndraw)};
REL.out.occ_varcov = {localConstCov(So,  ndraw)};
REL.out.to_varcov  = {localConstCov(Sto, ndraw)};
REL.out.b_sigma    = {num2cell(repmat(log(sigma), ndraw, 1))};
errc = cov2(sigma(1), sigma(2), 0);
REL.out.err_varcov = {localConstCov(errc, ndraw)};
REL.out.rescor     = {num2cell(repmat(0, ndraw, 1))};
REL.out.elabels = {REL.events(:)};
REL.out.glabels = {{'none'}};

pd = struct();
pd.rel = REL;
pd.ver = '0-test';
pd.proc = struct('interactive', false);
end

function c = localConstCov(S, ndraw)
A = repmat(reshape(S, [1 2 2]), [ndraw 1 1]);
c = num2cell(A);
end
