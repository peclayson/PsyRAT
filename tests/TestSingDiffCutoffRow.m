classdef TestSingDiffCutoffRow < PsyRATTestBase
    % Wiring tests for the one-facet difference-score cutoff row ('sing_diff').
    %
    % These run the real psyrat_relsummary path and the real cutoff table, not
    % the kernel in isolation, because the defect being guarded lived in the
    % wiring: psyrat_relsummary evaluated the difference score at the two
    % INDEPENDENTLY derived per-event trial cutoffs and psyrat_depcutofft then
    % printed '---' for the count. Two things follow from that, and only an
    % end-to-end test sees either:
    %
    %   1. the row reported a design nobody asked about, and
    %   2. when one event never reached its cutoff, its -1 sentinel was passed
    %      straight through to psyrat_harmmean, which rejects non-positive
    %      trial counts, so the whole summary threw.
    %
    % The fixture is the DoD difference-mode virtual REL, which emits a genuine
    % ic_diff structure rather than a hand-assembled struct.

    methods (Test)
        function testCutoffRowUsesACommonTrialCountFromTheDifferenceScore(testCase)
            out = localSummarize(localFixture(), 0.5);
            rs = out.relsummary;
            diffcut = rs.group(1).diffscore.trlcutoff;

            % a real, reported count rather than a suppressed one
            testCase.verifyTrue(isfield(diffcut, 'ntrials'));
            testCase.verifyGreaterThan(diffcut.ntrials, 0);
            testCase.verifyEqual(diffcut.ntrials, round(diffcut.ntrials));

            % and the score really is the difference score at that common
            % count, evaluated through the production kernel
            direct = localDiffrelAt(out, diffcut.ntrials);
            testCase.verifyEqual(diffcut.pt, direct.pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(diffcut.ll, direct.ll, 'AbsTol', 1e-12);
            testCase.verifyEqual(diffcut.ul, direct.ul, 'AbsTol', 1e-12);
        end

        function testCutoffRowCountIsTheSmallestThatClearsTheThreshold(testCase)
            cutoff = 0.5;
            out = localSummarize(localFixture(), cutoff);
            n = out.relsummary.group(1).diffscore.trlcutoff.ntrials;

            testCase.verifyGreaterThanOrEqual( ...
                out.relsummary.group(1).diffscore.trlcutoff.pt, cutoff);

            if n > 1
                below = localDiffrelAt(out, n-1);
                testCase.verifyLessThan(below.pt, cutoff);
            end
        end

        function testCutoffTablePrintsTheCountInsteadOfAPlaceholder(testCase)
            out = localSummarize(localFixture(), 0.5);
            t = psyrat_depcutofft('psyrat_data', out, 'gui', 0);

            row = strcmp(t.Label, 'diff score');
            testCase.verifyTrue(any(row), ...
                'The cutoff table did not contain a difference-score row.');

            printed = t.Trial_Cutoff{row};
            testCase.verifyTrue(isnumeric(printed), ...
                'The difference-score row printed a placeholder, not a count.');
            testCase.verifyEqual(printed, ...
                out.relsummary.group(1).diffscore.trlcutoff.ntrials);
        end

        function testOneEventWithoutACutoffNoLongerThrows(testCase)
            % Regression. Zeroing event 2's between-person variance in both
            % representations makes that event unreachable at any trial count,
            % so its per-event cutoff is the -1 sentinel while event 1 still
            % resolves. The old code handed [tc1 -1] to psyrat_diffrel, which
            % failed inside psyrat_harmmean ("Harmonic mean input must be
            % finite values greater than zero") and took the whole summary
            % with it. Verified against the pre-fix file: this exact input
            % threw at psyrat_relsummary's cutoff call.
            vd = localFixture();

            sdid = cell2mat(vd.rel.out.sd_id{:,1});
            sdid(:,2) = 0;
            vd.rel.out.sd_id = {num2cell(sdid)};

            iv = cell2mat(vd.rel.out.id_varcov{:,1});
            iv(:,2,2) = 0; iv(:,1,2) = 0; iv(:,2,1) = 0;
            vd.rel.out.id_varcov = {num2cell(iv)};

            out = localSummarize(vd, 0.5);
            rs = out.relsummary;

            % the state the regression is about really was reached
            testCase.verifyEqual(rs.group(1).event(2).trlcutoff, -1);
            testCase.verifyGreaterThan(rs.group(1).event(1).trlcutoff, 0);

            % and the difference score still resolves on its own terms
            testCase.verifyGreaterThan( ...
                rs.group(1).diffscore.trlcutoff.ntrials, 0);
        end

        function testUnreachableDifferenceScoreFallsBackToObservedTrials(testCase)
            % When no common trial count clears the threshold the row reports
            % the observed-trial difference score with the -1 sentinel the
            % per-event rows already use, mirroring the test-retest sibling.
            %
            % Simply raising the threshold does not reach this state: if no
            % EVENT reaches its cutoff the summary returns early on
            % nogooddata and never gets here. The state that does reach it is
            % a difference score that is unreliable while both events are
            % individually fine, which is what perfectly correlated events
            % produce -- the between-person variance of the difference is
            % zero, so no number of trials helps. Only the difference-score
            % block is touched; the per-event scaffolding (sd_id) is left
            % alone so both events still resolve their own cutoffs.
            vd = localFixture();
            iv = cell2mat(vd.rel.out.id_varcov{:,1});
            v = mean(iv(:,1,1));
            iv(:,1,1) = v; iv(:,2,2) = v; iv(:,1,2) = v; iv(:,2,1) = v;
            vd.rel.out.id_varcov = {num2cell(iv)};

            [out, relerr] = localSummarize(vd, 0.5);
            rs = out.relsummary;
            ds = rs.group(1).diffscore;

            % both events still resolved, so the early return did not fire
            testCase.verifyGreaterThan(rs.group(1).event(1).trlcutoff, 0);
            testCase.verifyGreaterThan(rs.group(1).event(2).trlcutoff, 0);

            testCase.verifyEqual(ds.trlcutoff.ntrials, -1);
            testCase.verifyEqual(ds.trlcutoff.pt, ds.pt, 'AbsTol', 0);

            % G60 wave: this arm ALIASES the parent evaluation, and it used
            % to carry the parent's admissibility record forward, which made
            % psyrat_admissibility_collect gather the identical record twice
            % (the found-cutoff arm's independent evaluation keeps its own).
            % The trt_diff fallback arm carries the identical two-line fix;
            % no fixture reaches it decisively, so this pin represents the
            % family.
            testCase.verifyFalse(isfield(ds.trlcutoff, 'admissibility'), ...
                ['the aliased fallback must not duplicate the parent''s ' ...
                'admissibility record']);

            % the -1 must be accompanied by the flag every other producer of
            % that sentinel raises, or the user sees it with no explanation
            testCase.verifyEqual(relerr.trlcutoff, 1);

            % the table prints the -1 sentinel, exactly as a per-event row
            % without a cutoff does
            t = psyrat_depcutofft('psyrat_data', out, 'gui', 0);
            testCase.verifyEqual(t.Trial_Cutoff{strcmp(t.Label,'diff score')}, -1);
        end

        function testExtrapolatedDifferenceCutoffRaisesTheTrialMaxFlag(testCase)
            % A derived count above every observed trial count is an
            % extrapolation beyond the data. The per-event rows already flag
            % that condition (relerr.trlmax, which raises the "cutoffs
            % represent an extrapolation" dialog); the difference row has to
            % as well, or it becomes the one row in the table that can print
            % an unreachable count while looking reachable.
            %
            % Raising the between-person covariance shrinks the variance of
            % the difference without touching either event's own reliability,
            % so both events still resolve small cutoffs while the difference
            % score needs more trials than anyone provided. That is the
            % ordinary situation for difference scores, not a contrived one.
            vd = localFixture();
            iv = cell2mat(vd.rel.out.id_varcov{:,1});
            iv(:,1,2) = 0.55; iv(:,2,1) = 0.55;   % PSD: 0.55 < sqrt(bp1*bp2)
            vd.rel.out.id_varcov = {num2cell(iv)};

            [out, relerr] = localSummarize(vd, 0.5);
            rs = out.relsummary;

            % the bound is the SMALLER of the two maxima: a common n' needs
            % that many trials in both conditions
            testCase.verifyGreaterThan( ...
                rs.group(1).diffscore.trlcutoff.ntrials, localObsMin(rs), ...
                'Fixture did not produce an extrapolated difference cutoff.');
            testCase.verifyEqual(relerr.trlmax, 1);
        end

        function testCutoffBetweenTheTwoEventMaximaIsStillAnExtrapolation(testCase)
            % Regression for the bound itself, and the case the other two
            % cannot see. When the two events have UNEQUAL maxima, a derived
            % n' that sits between them is unreachable -- no participant
            % supplied that many trials in the sparser condition -- but it is
            % not above the larger maximum. Comparing against the larger one
            % leaves this whole band unflagged, which for an unequal-trial
            % design (the case this row exists for) is most of the realistic
            % range.
            %
            % The fixture drops event 2's trial counts so the two maxima
            % differ, then shrinks the variance of the difference enough to
            % push the derived n' into the gap between them.
            vd = localFixture();

            enames = unique(vd.rel.data.event);
            drop = strcmp(vd.rel.data.event, enames{2}) & ...
                mod((1:height(vd.rel.data))', 3) ~= 0;
            vd.rel.data = vd.rel.data(~drop,:);

            % 0.40 puts the derived n' at 8, between the two maxima (4 and
            % 11). Raising it further pushes n' past the larger maximum, where
            % both the old and the new bound flag it and the test would prove
            % nothing.
            iv = cell2mat(vd.rel.out.id_varcov{:,1});
            iv(:,1,2) = 0.40; iv(:,2,1) = 0.40;
            vd.rel.out.id_varcov = {num2cell(iv)};

            [out, relerr] = localSummarize(vd, 0.5);
            rs = out.relsummary;
            lo = localObsMin(rs);
            hi = max([rs.group(1).event(1).trlinfo.max ...
                rs.group(1).event(2).trlinfo.max]);
            n = rs.group(1).diffscore.trlcutoff.ntrials;

            % verify, not assume: an assumption failure is reported as an
            % incomplete, which runAllUnitTests does not count as a failure, so
            % a fixture that stopped producing the state under test would
            % silently retire this test rather than failing it
            testCase.verifyGreaterThan(hi, lo, ...
                'Fixture did not produce unequal per-event maxima.');
            testCase.verifyGreaterThan(n, lo, ...
                'Derived n'' was not above the smaller maximum.');
            testCase.verifyLessThanOrEqual(n, hi, ...
                'Derived n'' was not inside the band between the two maxima.');

            % strictly inside (min, max], so this fails against a max() bound
            testCase.verifyEqual(relerr.trlmax, 1, ...
                sprintf(['A cutoff of %d is unreachable in the sparser ' ...
                'condition (max %d) but was not flagged as an ' ...
                'extrapolation.'], n, lo));
        end

        function testAReachableCutoffRaisesNeitherFlag(testCase)
            % The companion to the two above: when the difference score
            % reaches the threshold within the observed data, neither flag
            % fires, so the new flags cannot be passing by always being set.
            [out, relerr] = localSummarize(localFixture(), 0.5);
            rs = out.relsummary;

            testCase.verifyLessThanOrEqual( ...
                rs.group(1).diffscore.trlcutoff.ntrials, localObsMin(rs));
            testCase.verifyEqual(relerr.trlmax, 0);
            testCase.verifyEqual(relerr.trlcutoff, 0);
        end
    end
end

function vd = localFixture()
%A genuine one-facet difference-score REL, built from the DoD fixture's
%difference mode so the summary path is exercised rather than simulated.
pd = PsyRATTestDataFactory.makeDoDRelData();
vd = psyrat_dod_build_virtual_rel('psyrat_data', pd, ...
    'mode', 'diff', 'groupidx', 1, 'eventpair', [1 4]);
end

function lo = localObsMin(rs)
%The largest common trial count the data can support: a difference score at a
%common n' needs n' trials in BOTH conditions, so the binding limit is the
%SMALLER of the two events' observed maxima, not the larger.
lo = min([rs.group(1).event(1).trlinfo.max rs.group(1).event(2).trlinfo.max]);
end

function [out, relerr] = localSummarize(vd, depcutoff)
[out, relerr] = psyrat_relsummary('psyrat_data', vd, ...
    'analysis', 'sing_diff', ...
    'depcutoff', depcutoff, ...
    'meascutoff', 2, ...
    'depcentmeas', 1, ...
    'diffgcoeff', 1, ...
    'CI', 0.95);
end

function ds = localDiffrelAt(out, ntrl)
%Evaluate the production difference-score kernel at a common trial count,
%pulling the components the summary itself used.
REL = out.rel;
bsig = cell2mat(REL.out.b_sigma{:,1});
wp_cov = zeros(size(bsig,1),1);
if isfield(REL.out,'err_varcov') && ~isempty(REL.out.err_varcov) && ...
        ~isempty(REL.out.err_varcov{1,1})
    errvc = cell2mat(REL.out.err_varcov{:,1});
    wp_cov = errvc(:,2,1);
end

ds = psyrat_diffrel( ...
    'bp', cell2mat(REL.out.id_varcov{:,1}), ...
    'bt', cell2mat(REL.out.trl_varcov{:,1}), ...
    'er_var', bsig, ...
    'wp_cov', wp_cov, ...
    'obs', [ntrl ntrl], ...
    'est', 'dep', ...
    'CI', 0.95);
end
