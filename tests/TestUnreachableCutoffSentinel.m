classdef TestUnreachableCutoffSentinel < PsyRATTestBase
    %G54 regression pins. When no trial count in the projected range reaches
    %the reliability threshold, sing case 1 and trt case 1 store the -1
    %sentinel -- and before G54 both then fell through to the found-cutoff
    %arm, where GroupCount >= -1 marked EVERY participant good. That
    %contradicted cases 2-4 in both families, the trt_diff helper's
    %sentinel arm, and the manual ch. 9 promise that an unreachable cutoff
    %places every participant in Bad IDs (trt_sserr has no cutoff search
    %and no sentinel, so it is not part of that comparison). G54 respelled
    %the second guard on the sentinel VALUE (trlcutoff == -1), reviving the
    %everyone-bad arms.
    %
    %Every unreachable-threshold pin here fails decisively under the pre-G54
    %code (which retained everyone: goodn = 3, goodids = the full roster), and
    %each family carries a reachable-threshold POSITIVE CONTROL on the same
    %fixture so the everyone-bad assertions cannot pass for a fixture that is
    %simply doomed at any threshold (the differential-assertion rule).
    %
    %The cell-shaped variants pin the revived arms' REL.data normalization:
    %psyrat_computevarcomp stores a bare table on most cases and a per-event
    %cell on its event-array cases, so the revived arms normalize both
    %shapes like their live siblings (their pre-G54 dead bodies read
    %REL.data{1} bare in sing, and wrapped a plain assignment -- which
    %cannot throw -- in try/catch in trt). A cell paired with events='none'
    %is not a shipped shape; the variants pin the normalization itself.

    methods (Test)

        function testSingCaseOneUnreachableCutoffMarksEveryoneBad(testCase)
            pd = localSingRelData('table');
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'sing', 'depcutoff', 0.9999, 'meascutoff', 2, ...
                'depcentmeas', 1, 'CI', 0.95);
            localVerifySentinelAllBad(testCase, out, relerr);
        end

        function testSingCaseOneUnreachableCutoffCellData(testCase)
            pd = localSingRelData('cell');
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'sing', 'depcutoff', 0.9999, 'meascutoff', 2, ...
                'depcentmeas', 1, 'CI', 0.95);
            localVerifySentinelAllBad(testCase, out, relerr);
        end

        function testSingCaseOneReachableCutoffRetainsEveryone(testCase)
            %Positive control: the same fixture at a reachable threshold
            %retains the full roster, so the all-bad pins above are driven by
            %the sentinel state, not by the fixture.
            pd = localSingRelData('table');
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'sing', 'depcutoff', 0.30, 'meascutoff', 2, ...
                'depcentmeas', 1, 'CI', 0.95);
            localVerifyEveryoneRetained(testCase, out, relerr);
        end

        function testTrtCaseOneUnreachableCutoffMarksEveryoneBad(testCase)
            pd = localTrtRelData('table');
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'trt', 'gcoeff', 1, 'reltype', 3, ...
                'relcutoff', 0.9999, 'meascutoff', 2, 'relcentmeas', 1, ...
                'CI', 0.95);
            localVerifySentinelAllBad(testCase, out, relerr);
            %The icc/betsd/witsd block moved OUTSIDE the cutoff if/else at
            %G54 so the sentinel path also produces these fields (matching
            %case 2 and the sing sibling, which compute them on common
            %paths). Before the move the trt sentinel arm returned without
            %them and any consumer read would have crashed.
            ev = out.relsummary.group(1).event(1);
            %assert (not verify) before dereferencing: a reverted block-move
            %must fail with the named diagnostic, not error mid-test
            testCase.assertTrue(isfield(ev, 'icc'), ...
                'trt sentinel path must carry icc');
            testCase.verifyTrue(all(isfinite([ev.icc.ll ev.icc.m ev.icc.ul])), ...
                'icc summaries must be finite on the sentinel path');
            testCase.assertTrue(isfield(ev, 'betsd') && isfield(ev, 'witsd'), ...
                'trt sentinel path must carry betsd and witsd');
            testCase.verifyTrue(isfinite(ev.betsd.m) && isfinite(ev.witsd.m), ...
                'sd summaries must be finite on the sentinel path');
        end

        function testTrtCaseOneUnreachableCutoffCellData(testCase)
            pd = localTrtRelData('cell');
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'trt', 'gcoeff', 1, 'reltype', 3, ...
                'relcutoff', 0.9999, 'meascutoff', 2, 'relcentmeas', 1, ...
                'CI', 0.95);
            localVerifySentinelAllBad(testCase, out, relerr);
        end

        function testTrtCaseOneReachableCutoffRetainsEveryone(testCase)
            %Positive control for the trt family (see the sing twin).
            pd = localTrtRelData('table');
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'trt', 'gcoeff', 1, 'reltype', 3, ...
                'relcutoff', 0.20, 'meascutoff', 2, 'relcentmeas', 1, ...
                'CI', 0.95);
            localVerifyEveryoneRetained(testCase, out, relerr);
        end

    end
end

function localVerifySentinelAllBad(testCase, out, relerr)
%Shared G54 assertions for the unreachable-cutoff state: the sentinel is
%stored, the failure is reported through relerr.trlcutoff, and the ID
%partition is everyone-bad -- the exact inverse of the pre-G54 behavior.
ev = out.relsummary.group(1).event(1);
grp = out.relsummary.group(1);

testCase.verifyEqual(relerr.trlcutoff, 1, ...
    'An unreachable threshold must set relerr.trlcutoff');
testCase.verifyEqual(ev.trlcutoff, -1, ...
    'The stored trial cutoff must be the -1 sentinel');

testCase.verifyTrue(iscell(ev.eventgoodids) && isempty(ev.eventgoodids), ...
    'Sentinel state: eventgoodids must be the empty cell carrier');
testCase.verifyEqual(sort(ev.eventbadids), {'s1'; 's2'; 's3'}, ...
    'Sentinel state: every participant lands in eventbadids');

testCase.verifyTrue(iscell(grp.goodids) && isempty(grp.goodids), ...
    'Sentinel state: group goodids must be the empty cell carrier');
testCase.verifyEqual(sort(grp.badids), {'s1'; 's2'; 's3'}, ...
    'Sentinel state: every participant lands in the group badids');

testCase.verifyEqual(ev.goodn, 0, ...
    'Sentinel state: goodn must report zero retained participants');

%The revived arms' isempty(badids) backstop must not fire (badids holds
%the entire roster there), and trial info must come out real. On this
%complete-by-design fixture the sentinel arm's bad-roster join and the
%found arm's good-roster join cover the same subjects, so the trlinfo
%check is a LIVENESS pin (the join ran and produced finite stats), not a
%which-join-ran discriminator.
testCase.verifyEqual(relerr.nogooddata, 0, ...
    'The revived arm''s backstop must not fire when data exist');
testCase.verifyTrue(isfinite(ev.trlinfo.max) && ev.trlinfo.max > 0, ...
    'Trial info must be computed and finite on the sentinel path');
end

function localVerifyEveryoneRetained(testCase, out, relerr)
%Positive control: at a reachable threshold the same fixtures retain the
%full roster, so the sentinel pins cannot pass vacuously.
ev = out.relsummary.group(1).event(1);
grp = out.relsummary.group(1);

testCase.verifyEqual(relerr.trlcutoff, 0, ...
    'A reachable threshold must not flag the cutoff failure');
testCase.verifyGreaterThan(ev.trlcutoff, 0, ...
    'A reachable threshold must store a real cutoff');
testCase.verifyEqual(sort(ev.eventgoodids), {'s1'; 's2'; 's3'}, ...
    'Positive control: everyone qualifies at the reachable threshold');
testCase.verifyTrue(iscell(ev.eventbadids) && isempty(ev.eventbadids), ...
    'Positive control: nobody lands in eventbadids');
testCase.verifyEqual(sort(grp.goodids), {'s1'; 's2'; 's3'}, ...
    'Positive control: the group carrier retains everyone');
testCase.verifyEqual(ev.goodn, 3, ...
    'Positive control: goodn reports the full roster');
end

function psyrat_data = localSingRelData(shape)
%Minimal PRE-relsummary one-facet struct (rel.out draws plus a rel.data trial
%table) routed to sing analysis 1 (groups and events both 'none'). Three
%participants x five trials. With between-person sig_u ~2 against residual
%sig_e ~5, the projected dependability at the search ceiling (max trials +
%1000) sits near 0.994, so a 0.9999 threshold is unreachable while 0.30 is
%met by the third projected trial -- inside every participant's five.
subs = {'s1','s2','s3'};
ntrl = 5;
ids = {}; meas = [];
for s = 1:numel(subs)
    for t = 1:ntrl
        ids{end+1,1} = subs{s}; %#ok<AGROW>
        meas(end+1,1) = numel(meas) + 1; %#ok<AGROW>
    end
end

rel = struct();
rel.analysis = 'sing';
rel.groups = 'none';
rel.events = 'none';
rel.filename = 'g54_sing_sentinel_fixture.psyrat';
rel.nchains = 2;
rel.niter = 100;
rel.data = table(ids, meas, 'VariableNames', {'id','meas'});
if strcmp(shape, 'cell')
    rel.data = {rel.data};
end

rel.out = struct();
rel.out.labels = {'measure'};
rel.out.mu    = [0.10; 0.12; 0.11; 0.09; 0.10];
rel.out.sig_u = [2.00; 2.10; 1.90; 2.05; 1.95];
rel.out.sig_e = [4.90; 5.00; 5.10; 4.80; 5.20];
rel.out.sig_trl = [0.20; 0.25; 0.23; 0.21; 0.22];

psyrat_data = struct('rel', rel);
%headless: never block on a modal even if one regresses back in
psyrat_data.proc = struct('interactive', false, 'measheader', 'meas');
end

function psyrat_data = localTrtRelData(shape)
%Minimal PRE-relsummary test-retest struct routed to trt analysis 1. Three
%participants x two occasions x five trials. The component draws are a copy
%of the RC-33 criterion fixture's (large sig_err against sig_id ~2). The
%thresholds are derived from the formula, not from RC-33 (whose committed
%test exercises only 0.9999, and its first call at a different occasion
%estimand): at reltype 3 with nocc 1 the occasion terms put a hard
%asymptote near 0.94, so 0.9999 is unreachable at any trial count, and the
%mean coefficient crosses 0.20 by the second projected trial, inside every
%participant's five.
subs = {'s1','s2','s3'};
occ = {'t1','t2'};
ntrl = 5;
ids = {}; times = {}; meas = [];
for s = 1:numel(subs)
    for o = 1:numel(occ)
        for t = 1:ntrl
            ids{end+1,1} = subs{s}; %#ok<AGROW>
            times{end+1,1} = occ{o}; %#ok<AGROW>
            meas(end+1,1) = numel(meas) + 1; %#ok<AGROW>
        end
    end
end

rel = struct();
rel.analysis = 'trt';
rel.groups = 'none';
rel.events = 'none';
rel.filename = 'g54_trt_sentinel_fixture.psyrat';
rel.nchains = 2;
rel.niter = 100;
rel.time = occ(:);
rel.data = table(ids, meas, times, 'VariableNames', {'id','meas','time'});
if strcmp(shape, 'cell')
    rel.data = {rel.data};
end

rel.out = struct();
rel.out.labels = {'measure'};
rel.out.mu = [0.10; 0.12; 0.11; 0.09; 0.10];
rel.out.sig_id      = [2.00; 2.10; 1.90; 2.05; 1.95];
rel.out.sig_occ     = [0.30; 0.35; 0.28; 0.32; 0.31];
rel.out.sig_trl     = [0.20; 0.25; 0.23; 0.21; 0.22];
rel.out.sig_trlxid  = [0.20; 0.22; 0.18; 0.21; 0.19];
rel.out.sig_occxid  = [0.40; 0.42; 0.38; 0.41; 0.39];
rel.out.sig_trlxocc = [0.09; 0.10; 0.11; 0.09; 0.10];
rel.out.sig_err     = [4.90; 5.00; 5.10; 4.80; 5.20];

psyrat_data = struct('rel', rel);
%headless: never block on a modal even if one regresses back in
psyrat_data.proc = struct('interactive', false, 'measheader', 'meas');
end
