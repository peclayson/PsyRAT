classdef TestSingDiffSserrStrictAlignment < PsyRATTestBase
    %G48 (owner ruling 2026-08-29: FULL alignment with the G39/G45
    %strict-fold conventions). The subject-level difference helper
    %(psyrat_relsummary_sing_diff_sserr) was the LAST route where an empty
    %inclusion set the run-global relerr.nogooddata flag inside its loop and
    %raised the one-per-run modal after it, so every caller aborted with no
    %output -- and because the flag is run-global, ONE empty group blocked
    %the output for every group. This is exactly the divergence the G45
    %ruling removed from trt_diff, on a GUI-reachable route
    %(ic_diff_sserrvar viewing).
    %
    %Post-G48: the empty state warns on the console per group and the
    %summary is returned in full (nogooddata stays 0 on every path); the
    %overall table renders the zeroed group through psyrat_depoverallt's
    %G44 dash-out arm with complementing counts. The G48-era claim that "no
    %goodn rekey is needed here" survived only in part: no shipped consumer
    %reads event.goodn on this route (that half stands), but depoverallt's
    %sserr branch used to print goodn as the event indicator summed over the
    %group-carrier rows, which under-counts when a subject is diff-included
    %but event-excluded -- the G53 filing. Since G53 the branch prints the
    %group-carrier count itself, so the rendered rows complement by
    %construction; the divergent-state pin below exercises exactly the
    %configuration the shared fixture cannot express.
    %
    %Retention is staged through the depcutoff against the fixture's
    %per-subject difference-score dependability draws, with cutoffs derived
    %in-test from a clean run so a drifted fixture fails the anti-vacuity
    %welds loudly instead of passing against the wrong state. The two-group
    %blanking pin stages the empty group via all-NaN per-subject residual
    %draws (a failed fit), which read as excluded under the G49 marker.

    methods (Test)

        function testEmptyInclusionNoticesAndContinues(testCase)
            %The abort replacement: nobody clears the staged cutoff, so the
            %difference-score inclusion is empty. Pre-G48 this set
            %nogooddata and every caller aborted with no output; post-G48
            %the run continues with a zeroed, complementing group, a
            %per-group console notice, and informative per-event rows.
            %
            %Honesty note (adversarial wave, 2026-08-29): the helper
            %RETURNED its summary in full pre-G48 too -- the abort vehicle
            %was the flag, consumed by the callers -- so the structure
            %asserts at the end are regression pins. The change-detecting
            %asserts here are nogooddata == 0, the console-notice text,
            %and (in the sibling class) the no-dialog check.
            pd = localFixture();
            [cutEmpty, ~] = localStagedCutoffs(testCase, pd);

            txt = evalc('[out, relerr] = localSummary(pd, cutEmpty);');

            g = out.relsummary.group(1);
            T = g.diffscore.ssrel_table;
            testCase.assertFalse(any(T.ind2include), ...
                'weld: nobody may clear the staged empty cutoff');

            testCase.verifyEqual(relerr.nogooddata, 0, ...
                'an empty inclusion must no longer set the abort flag');
            testCase.verifyTrue(iscell(g.goodids) && isempty(g.goodids));
            testCase.verifyEqual(sort(g.badids), sort({'s1'; 's2'; 's3'}), ...
                'the zeroed group must complement over the full roster');
            testCase.verifyTrue(contains(txt, 'difference score'), ...
                'the per-group console notice must name the condition');
            %the summary is returned in full: both events and the
            %difference-score block are present and populated
            testCase.verifyEqual(numel(g.event), 2);
            testCase.verifyTrue(istable(g.event(1).ssrel_table));
            testCase.verifyTrue(isfield(g.diffscore, 'trlinfo'));
        end

        function testOverallTableRendersZeroedGroupWithComplement(testCase)
            %The user-facing half: psyrat_depoverallt must render the zeroed
            %group -- per-event rows through the G44 dash-out arm, the diff
            %row from the marked diff table -- with every count pair
            %complementing to the 3-subject roster.
            %
            %Honesty note (adversarial wave, 2026-08-29): this is a
            %REGRESSION pin, not a G48 differential -- pre-G48 the helper
            %still published the summary (only callers aborted on the
            %flag), so this table rendered then too. It exists so the
            %G44/complement rendering cannot be lost while the route is
            %aligned; the change-detecting pins are the sibling tests'.
            pd = localFixture();
            [cutEmpty, ~] = localStagedCutoffs(testCase, pd);
            out = localSummary(pd, cutEmpty);
            testCase.assertTrue(isempty(out.relsummary.group(1).goodids), ...
                'weld: the staged group must retain nobody');

            tbl = psyrat_depoverallt('psyrat_data', out, 'gui', 0);
            testCase.assertEqual(height(tbl), 3, ...
                'expected two event rows plus the diff row');
            for r = 1:3
                testCase.verifyEqual(cell2mat(tbl{r,2}), 0, sprintf( ...
                    'row %d: n Included must be 0 for the zeroed group', r));
                testCase.verifyEqual(cell2mat(tbl{r,3}), 3, sprintf( ...
                    'row %d: n Excluded must complement over the roster', r));
            end
            %the per-event coefficient cells are dashed (G44: a mean over
            %zero retained rows is undefined)
            testCase.verifyEqual(tbl{1,4}{1}, '---');
            testCase.verifyEqual(tbl{2,4}{1}, '---');
            %G57: the DIFF row's coefficient and five trial cells dash too.
            %Pre-fix, this row fell back to the FULL table and printed a
            %coefficient plus trial statistics computed over every excluded
            %subject beside n Included = 0.
            for c = 4:9
                testCase.verifyEqual(tbl{3,c}{1}, '---', sprintf( ...
                    'diff row col %d: zero retention must dash, not render', c));
            end
        end

        function testPartialRetentionComplements(testCase)
            %Partial retention: the retained set and its complement must
            %partition the roster, with no modal and no abort flag.
            %
            %Honesty note (adversarial wave, 2026-08-29): a regression pin
            %-- partial retention never set the flag, so this passes
            %against pre-G48 code too. It also structurally CANNOT catch
            %the per-event/diff inclusion independence on this route (the
            %fixture's diff error variance exceeds either event's, so
            %event-included is a superset of diff-included by
            %construction); that independence was the separate ledger
            %filing G53, closed by the group-carrier rekey in
            %psyrat_depoverallt and pinned by
            %testDivergentInclusionRowsStillComplement, whose fixture
            %flips the between-person covariance to express the state this
            %fixture cannot.
            pd = localFixture();
            [~, cutPartial] = localStagedCutoffs(testCase, pd);
            [out, relerr] = localSummary(pd, cutPartial);

            g = out.relsummary.group(1);
            n = numel(g.goodids);
            testCase.assertTrue(n >= 1 && n <= 2, ...
                'weld: the staged cutoff must retain a strict subset');

            testCase.verifyEqual(relerr.nogooddata, 0);
            testCase.verifyEqual(numel(g.badids), 3 - n, ...
                'good/bad ids must complement over the roster');
            testCase.verifyEmpty(intersect(g.goodids, g.badids));

            tbl = psyrat_depoverallt('psyrat_data', out, 'gui', 0);
            for r = 1:height(tbl)
                testCase.verifyEqual( ...
                    cell2mat(tbl{r,2}) + cell2mat(tbl{r,3}), 3, sprintf( ...
                    'row %d: rendered counts must complement to the roster', r));
            end
        end

        function testDivergentInclusionRowsStillComplement(testCase)
            %G53 change pin. On ic_diff_sserrvar the group carrier derives
            %from the DIFFERENCE-score inclusion while each event table is
            %marked independently against the same cutoff -- two different
            %dep_pt quantities with no enforced ordering. With a negative
            %between-person covariance the difference's true-score variance
            %(bp1 + bp2 - 2*bp_cov) inflates past either event's, so a
            %subject can clear the cutoff on the difference score while
            %failing it in an event. Pre-G53, depoverallt's event rows then
            %summed short of the roster (goodn counted "diff-included AND
            %event-included" while badn counted roster minus diff-goodids),
            %so the complement loop below fails decisively against that
            %code.
            pd = localDivergentFixture();

            %derive the divergence-splitting cutoff from a clean run of the
            %SAME fixture (the class convention): for the subject with the
            %largest diff-over-event gap, a cutoff at the midpoint includes
            %the subject on the difference score and excludes it in the
            %weaker event.
            clean = localSummary(pd, 0.01);
            g = clean.relsummary.group(1);
            diffT = g.diffscore.ssrel_table;
            e1 = g.event(1).ssrel_table;
            e2 = g.event(2).ssrel_table;
            gap = -inf; cut = NaN;
            for i = 1:height(diffT)
                p1 = e1.dep_pt(strcmp(e1.id, diffT.id{i}));
                p2 = e2.dep_pt(strcmp(e2.id, diffT.id{i}));
                lo = min(p1, p2);
                if diffT.dep_pt(i) - lo > gap
                    gap = diffT.dep_pt(i) - lo;
                    cut = (diffT.dep_pt(i) + lo) / 2;
                end
            end
            testCase.assertGreaterThan(gap, 0, ['anti-vacuity: the ' ...
                'doctored covariance must push some subject''s difference ' ...
                'coefficient above one of its event coefficients']);

            [out, relerr] = localSummary(pd, cut);
            g2 = out.relsummary.group(1);
            testCase.verifyEqual(relerr.nogooddata, 0);
            testCase.assertGreaterThanOrEqual(numel(g2.goodids), 1, ...
                'weld: the staged cutoff must retain someone on the difference score');

            %positive control: the divergent state must actually have
            %materialized after marking -- without a diff-included but
            %event-excluded subject the complement loop cannot fail
            divergent = false;
            for e = 1:2
                T = g2.event(e).ssrel_table;
                rows = ismember(T.id, g2.goodids);
                if any(~T.ind2include(rows))
                    divergent = true;
                end
            end
            testCase.assertTrue(divergent, ['positive control: at least ' ...
                'one diff-included subject must be event-excluded, or ' ...
                'this pin cannot fail']);

            tbl = psyrat_depoverallt('psyrat_data', out, 'gui', 0);
            testCase.assertEqual(height(tbl), 3, ...
                'expected two event rows plus the diff row');
            for r = 1:3
                testCase.verifyEqual( ...
                    cell2mat(tbl{r,2}) + cell2mat(tbl{r,3}), 3, sprintf( ...
                    'row %d: rendered counts must complement to the roster', r));
            end
            %the G53 semantic: event rows print the GROUP-carrier count
            testCase.verifyEqual(cell2mat(tbl{1,2}), numel(g2.goodids), ...
                'event rows must print the group-carrier count');
            testCase.verifyEqual(cell2mat(tbl{2,2}), numel(g2.goodids), ...
                'event rows must print the group-carrier count');
        end

        function testDivergentPartialRetentionComplements(testCase)
            %G53 wave follow-up. The gap-splitting staging above lands in
            %the badn == 0 corner (its cutoff retains everyone on the
            %difference score), where the carrier assertions are
            %arithmetically identical to the complement ones. This staging
            %cuts between the two LARGEST diff coefficients instead, so
            %the diff retains a strict subset while a retained subject is
            %still event-excluded: the complement welds see a nonzero
            %excluded count, and the carrier assertions become independent
            %of them (pre-G53 the event rows here rendered the indicator
            %sum 0, not the carrier count).
            pd = localDivergentFixture();
            clean = localSummary(pd, 0.01);
            g = clean.relsummary.group(1);
            diffT = g.diffscore.ssrel_table;
            pts = sort(diffT.dep_pt, 'descend');
            testCase.assertGreaterThan(pts(1), pts(2), ...
                'weld: the top two diff coefficients must be distinct');
            cut = (pts(1) + pts(2)) / 2;

            %anti-vacuity: the subject this cutoff retains must fail it in
            %at least one event, or the staging expresses no divergence
            e1 = g.event(1).ssrel_table;
            e2 = g.event(2).ssrel_table;
            iTop = find(diffT.dep_pt == pts(1), 1);
            pTop = [e1.dep_pt(strcmp(e1.id, diffT.id{iTop})), ...
                e2.dep_pt(strcmp(e2.id, diffT.id{iTop}))];
            testCase.assertTrue(any(pTop < cut), ['anti-vacuity: the ' ...
                'retained subject must fail the cutoff in some event']);

            [out, relerr] = localSummary(pd, cut);
            g2 = out.relsummary.group(1);
            testCase.verifyEqual(relerr.nogooddata, 0);
            n = numel(g2.goodids);
            testCase.assertTrue(n >= 1 && n <= 2, ...
                'weld: the staged cutoff must retain a strict subset');

            %positive control, post-marking: a retained subject really is
            %event-excluded, so the pre-G53 indicator sum undercounts
            divergent = false;
            for e = 1:2
                T = g2.event(e).ssrel_table;
                rows = ismember(T.id, g2.goodids);
                if any(~T.ind2include(rows))
                    divergent = true;
                end
            end
            testCase.assertTrue(divergent, ['positive control: a ' ...
                'diff-included subject must be event-excluded here']);

            tbl = psyrat_depoverallt('psyrat_data', out, 'gui', 0);
            testCase.assertEqual(height(tbl), 3, ...
                'expected two event rows plus the diff row');
            for r = 1:3
                testCase.verifyEqual( ...
                    cell2mat(tbl{r,2}) + cell2mat(tbl{r,3}), 3, sprintf( ...
                    'row %d: rendered counts must complement to the table roster', r));
                testCase.verifyGreaterThan(cell2mat(tbl{r,3}), 0, sprintf( ...
                    'row %d: the excluded count must be nonzero here', r));
            end
            testCase.verifyEqual(cell2mat(tbl{1,2}), n, ...
                'event rows must print the group-carrier count');
            testCase.verifyEqual(cell2mat(tbl{2,2}), n, ...
                'event rows must print the group-carrier count');
        end

        function testOneEmptyGroupDoesNotBlankTheOther(testCase)
            %The run-global half of the old abort: ONE empty group used to
            %block the output for every group. Stage two groups, fail the
            %first via all-NaN per-subject residual draws (a failed fit,
            %excluded under the G49 marker), and require the second group's
            %retention and rendered rows to survive.
            %
            %Honesty note (adversarial wave, 2026-08-29): the blanking
            %lived in the CALLERS (psyrat_relfigures returns on the
            %run-global flag), so the change-detecting asserts here are
            %nogooddata == 0 and the notice; the table asserts pin that
            %both groups' rows render from the returned summary, which the
            %relfigures test below completes at the caller level.
            pd = localTwoGroupFixture();
            %fail GroupA only: NaN both events' per-subject residual draws
            pd.rel.out.er_var_ss1{1,1} = nan(size(pd.rel.out.er_var_ss1{1,1}));
            pd.rel.out.er_var_ss2{1,1} = nan(size(pd.rel.out.er_var_ss2{1,1}));

            txt = evalc('[out, relerr] = localSummary(pd, 0.01);');

            gA = out.relsummary.group(1);
            gB = out.relsummary.group(2);
            testCase.assertTrue(all(isnan(gA.diffscore.ssrel_table.dep_pt)), ...
                'weld: GroupA''s diff dep draws must genuinely be NaN');
            testCase.assertEqual(sort(gB.goodids), sort({'s1'; 's2'; 's3'}), ...
                'weld: GroupB must retain everyone at the lax cutoff');

            testCase.verifyEqual(relerr.nogooddata, 0, ...
                'one empty group must not raise the run-global abort flag');
            testCase.verifyTrue(iscell(gA.goodids) && isempty(gA.goodids));
            testCase.verifyEqual(sort(gA.badids), sort({'s1'; 's2'; 's3'}));
            testCase.verifyTrue(contains(txt, 'difference score'), ...
                'the notice must fire for the failed group');

            %the rendered table carries BOTH groups: 3 rows each, GroupA
            %zeroed and complementing, GroupB fully retained
            tbl = psyrat_depoverallt('psyrat_data', out, 'gui', 0);
            testCase.assertEqual(height(tbl), 6, ...
                'both groups must render: 2 event rows + 1 diff row each');
            testCase.verifyEqual(cell2mat(tbl{1,2}), 0);
            testCase.verifyEqual(cell2mat(tbl{1,3}), 3);
            testCase.verifyEqual(cell2mat(tbl{4,2}), 3, ...
                'GroupB''s first event row must show full retention');
            testCase.verifyEqual(cell2mat(tbl{4,3}), 0);

            %G57: GroupA's diff row (row 3) must dash its coefficient and
            %trial cells at zero retention while the counts stay numeric
            %and complement. Pre-fix, the all-NaN diff draws rendered
            %' NaN SD: NaN' here via the full-table fallback.
            testCase.verifyEqual(cell2mat(tbl{3,2}), 0, ...
                'GroupA diff row: n Included must be 0');
            testCase.verifyEqual(cell2mat(tbl{3,3}), 3, ...
                'GroupA diff row: n Excluded must complement');
            for c = 4:9
                testCase.verifyEqual(tbl{3,c}{1}, '---', sprintf( ...
                    'GroupA diff row col %d: zero retention must dash', c));
            end
            %positive control: GroupB's diff row (row 6) still renders a
            %real coefficient through the untouched else-arm
            testCase.verifyEqual(cell2mat(tbl{6,2}), 3);
            testCase.verifyEqual(cell2mat(tbl{6,3}), 0);
            testCase.verifyTrue(ischar(tbl{6,4}{1}) && ...
                contains(tbl{6,4}{1}, 'SD:'), ...
                'GroupB diff row must keep its numeric coefficient cell');
        end

        function testRelfiguresDrawsOutputForAZeroedGroup(testCase)
            %The caller half of the abort replacement: psyrat_relfigures
            %re-summarizes and returns on nogooddata ahead of every
            %plot/table block, so pre-G48 a zeroed run produced nothing.
            %Post-G48 the flag stays 0 and the outputs are drawn.
            pd = localFixture();
            [cutEmpty, ~] = localStagedCutoffs(testCase, pd);

            prefs = struct('view', localViewPrefsStruct(cutEmpty), ...
                'ver', '0-test');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            close('all', 'force');
            nBefore = numel(findall(groot, 'Type', 'figure'));
            psyrat_relfigures('psyrat_data', pd, 'psyrat_prefs', prefs, ...
                'analysis', 'sing_diff_sserr');
            nAfter = numel(findall(groot, 'Type', 'figure'));
            testCase.verifyGreaterThan(nAfter, nBefore, ...
                'a zeroed sing_diff_sserr run must still produce output');
        end

        function testZeroRetentionFallbackHonorsMedian(testCase)
            %G60 fix proper. The zero-retention obs fallback used to
            %hardcode MEAN, silently ignoring a median request
            %(depcentmeas = 2) while the ic_diff sibling honors it. On a
            %fixture whose event-1 trial counts are skewed (mean ~= median),
            %the stored group-level coefficient must therefore DIFFER
            %between the two central-tendency choices at zero retention.
            %Pre-fix both runs produced the mean-based value, so the
            %inequality below fails decisively on the old code.
            pd = localSkewedFixture(testCase);
            [cutEmpty, ~] = localStagedCutoffs(testCase, pd);

            outMean = localSummaryCent(pd, cutEmpty, 1);
            outMed = localSummaryCent(pd, cutEmpty, 2);

            %welds: both runs are genuinely zero-retention, and the skew is
            %real on the diff table the fallback reads (anti-vacuity: with
            %mean == median this pin could not fail on ANY implementation)
            testCase.assertTrue(isempty(outMean.relsummary.group(1).goodids));
            testCase.assertTrue(isempty(outMed.relsummary.group(1).goodids));
            dtab = outMean.relsummary.group(1).diffscore.ssrel_table;
            testCase.assertTrue( ...
                mean(dtab.trls1) ~= median(dtab.trls1) || ...
                mean(dtab.trls2) ~= median(dtab.trls2), ...
                'weld: the skewed fixture must separate mean from median');

            testCase.verifyNotEqual( ...
                outMean.relsummary.group(1).diffscore.pt, ...
                outMed.relsummary.group(1).diffscore.pt, ...
                ['the zero-retention fallback must evaluate the design ' ...
                'point with the requested central tendency']);
        end

        function testCutoffRowStructIsMinimal(testCase)
            %G60 wave. diffscore.trlcutoff on this route used to be a full
            %self-copy of diffscore, which dragged the admissibility record
            %along and made psyrat_admissibility_collect gather the
            %IDENTICAL record twice (its second gather exists for routes
            %whose cutoff row is independently evaluated). The producer now
            %stores exactly the three fields the cutoff table's diff row
            %reads; no admissibility field is what the collector's
            %local_gather guard keys on.
            pd = localFixture();
            clean = localSummary(pd, 0.01);
            ds = clean.relsummary.group(1).diffscore;

            testCase.assertTrue(isstruct(ds.trlcutoff));
            testCase.verifyEqual(sort(fieldnames(ds.trlcutoff)), ...
                sort({'ll'; 'pt'; 'ul'}), ...
                'the cutoff-row struct must carry exactly pt/ll/ul');
            testCase.verifyEqual(ds.trlcutoff.pt, ds.pt, ...
                'the cutoff row IS the observed-design coefficient');
            testCase.verifyEqual(ds.trlcutoff.ll, ds.ll);
            testCase.verifyEqual(ds.trlcutoff.ul, ds.ul);
            testCase.verifyFalse(isfield(ds.trlcutoff, 'admissibility'), ...
                ['no admissibility field: carrying one double-counts the ' ...
                'parent record in the admissibility note']);
        end

        function testCutoffTableRendersObservedDiffRowAtZeroRetention(testCase)
            %G60 surface pin (owner ruling: keep + document). The trial-
            %cutoff table is the one surface that renders the group-level
            %diff coefficient on this route, and no test exercised it. At
            %zero retention the diff row keeps a NUMERIC coefficient (it is
            %draw-based -- the G44 kept-coefficient rule, reaffirmed by the
            %G60 value ruling) beside the '---' count placeholder (this
            %route runs no cutoff search). The event rows' count cells are
            %numeric on either producer arm (min over included, or min over
            %the full table when that event retained nobody), so they carry
            %no assertion here -- no code path can dash them, and a control
            %that cannot fail is not a control.
            pd = localFixture();
            [cutEmpty, ~] = localStagedCutoffs(testCase, pd);
            out = localSummary(pd, cutEmpty);
            testCase.assertTrue(isempty(out.relsummary.group(1).goodids), ...
                'weld: the staged run must retain nobody');

            tbl = psyrat_depcutofft('psyrat_data', out, 'gui', 0);
            testCase.assertEqual(height(tbl), 3, ...
                'two event rows plus the diff row');
            testCase.verifyEqual(tbl.Label{3}, 'diff score');
            testCase.verifyEqual(tbl{3,2}{1}, '---', ...
                'no cutoff search on this route: the count cell dashes');
            testCase.verifyTrue(ischar(tbl{3,3}{1}) && ...
                contains(tbl{3,3}{1}, 'CI [') && ...
                ~contains(tbl{3,3}{1}, 'NaN'), ...
                ['the coefficient cell stays a real number (draw-based; ' ...
                'the renderer''s sprintf would happily print NaN, so the ' ...
                'absence of NaN is part of the pin)']);
        end

    end
end

%--------------------------------------------------------------------------
function [out, relerr] = localSummary(pd, depcutoff)
[out, relerr] = psyrat_relsummary( ...
    'psyrat_data', pd, ...
    'analysis', 'sing_diff_sserr', ...
    'depcutoff', depcutoff, ...
    'meascutoff', 2, ...
    'depcentmeas', 1, ...
    'diffgcoeff', 1, ...
    'CI', 0.95);
end

function pd = localFixture()
pd = PsyRATTestDataFactory.makeICDiffSSErrRelData();
pd.ver = '0-test';
pd.rel.filename = 'sing_diff_sserr_fixture.psyrat';
%headless: never block on a modal even if one regresses back in
pd.proc = struct('interactive', false, 'measheader', 'meas');
end

function [out, relerr] = localSummaryCent(pd, depcutoff, centmeas)
%localSummary with a caller-chosen central tendency (G60 median pin).
[out, relerr] = psyrat_relsummary( ...
    'psyrat_data', pd, ...
    'analysis', 'sing_diff_sserr', ...
    'depcutoff', depcutoff, ...
    'meascutoff', 2, ...
    'depcentmeas', centmeas, ...
    'diffgcoeff', 1, ...
    'CI', 0.95);
end

function pd = localSkewedFixture(testCase)
%The shared fixture's trial counts are symmetric (trls1 = [6 7 8],
%trls2 = [5 6 7]: mean == median in both events), which would make the G60
%median pin vacuous. Drop two of s3's first-event rows so trls1 becomes
%[6 7 6] (mean 6.33, median 6). The draws are per-subject, not per-trial,
%so removing data rows changes only the counts the idtable builders derive.
pd = localFixture();
evs = unique(pd.rel.data.event, 'stable');
idx = find(strcmp(pd.rel.data.id, 's3') & ...
    strcmp(pd.rel.data.event, evs{1}));
testCase.assertGreaterThanOrEqual(numel(idx), 4, ...
    'weld: the fixture must have rows to spare for the skew');
pd.rel.data(idx(1:2), :) = [];
end

function pd = localDivergentFixture()
%G53 divergent-state fixture: identical to the shared fixture except the
%between-person covariance is NEGATIVE. In the SHARED fixture the ordering
%comes from the error term (the diff error variance exceeds either
%event's), so its diff coefficient sits below either event's and
%event-included is a superset of diff-included by construction. Flipping
%bp_cov to -0.40 inflates the difference's true-score variance
%(bp1 + bp2 - 2*bp_cov) to ~1.94 against event variances of ~0.60/0.54,
%which pushes the diff coefficient above at least one event's -- the
%min(p1,p2) gap the tests actually weld -- so a midpoint cutoff separates
%them. The doctored matrix stays positive definite
%(det = 0.60*0.54 - 0.40^2 = 0.164 > 0 at drift 0, growing with drift).
pd = localFixture();
ndraw = 5;
idVarcov = zeros(ndraw, 2, 2);
for d = 1:ndraw
    drift = (d - 1) * 0.01;
    idVarcov(d,:,:) = [0.60 + drift, -0.40; -0.40, 0.54 + drift];
end
pd.rel.out.id_varcov = {num2cell(idVarcov)};
end

function [cutEmpty, cutPartial] = localStagedCutoffs(testCase, pd)
%Derive the staged cutoffs from a clean run of the SAME fixture, so the
%tests key on the fixture's actual difference-score dependability draws
%rather than hardcoded values that could silently drift.
clean = localSummary(pd, 0.01);
pts = sort(clean.relsummary.group(1).diffscore.ssrel_table.dep_pt);
testCase.assertEqual(numel(pts), 3, ...
    'weld: the fixture must carry three per-subject diff dep values');
testCase.assertEqual( ...
    numel(clean.relsummary.group(1).goodids), 3, ...
    'anti-vacuity: the lax cutoff must retain the full roster');
testCase.assertGreaterThan(pts(3), pts(1), ...
    'weld: the per-subject values must differ for a partial cutoff to exist');
cutEmpty = pts(3) + (1 - pts(3)) / 2;
testCase.assertTrue(cutEmpty > pts(3) && cutEmpty < 1, ...
    'weld: the empty cutoff must clear every subject yet stay in (0,1)');
cutPartial = (pts(1) + pts(2)) / 2;
testCase.assertTrue(cutPartial > pts(1) && cutPartial <= pts(3), ...
    'weld: the partial cutoff must split the roster');
end

function pd = localTwoGroupFixture()
%Two-group variant of the shared ic_diff_sserrvar fixture: duplicate the
%posterior columns and stamp a group column onto the data rows, so the
%helper's per-group loop runs twice on identical draws until one group is
%doctored.
pd = localFixture();
rel = pd.rel;
rel.groups = {'GroupA'; 'GroupB'};
flds = {'b_sigma', 'id_varcov', 'trl_varcov', 'err_varcov', ...
    'er_var_ss1', 'er_var_ss2', 'wp_cov_ss', 'id_matches', 'elabels'};
for f = 1:numel(flds)
    rel.out.(flds{f})(:, 2) = rel.out.(flds{f})(:, 1);
end
rel.out.glabels = {'GroupA'; 'GroupB'};
rel.out.conv.data{2} = rel.out.conv.data{1};
dataA = rel.data;
dataA.group = repmat({'GroupA'}, height(dataA), 1);
dataB = rel.data;
dataB.group = repmat({'GroupB'}, height(dataB), 1);
rel.data = [dataA; dataB];
pd.rel = rel;
end

function v = localViewPrefsStruct(depvalue)
%The overall table alone: the surface the pre-G48 abort blanked, rendered
%through the G44 dash-out arm for the zeroed group.
v = struct();
v.gcoeff = 1; v.diffgcoeff = 1;
v.depvalue = depvalue; v.depcentmeas = 1; v.meascutoff = 2;
v.plotdep = 0; v.plotdepline = 2; v.ploticc = 0;
v.inctrltable = 0; v.overalltable = 1;
v.ntrials = 20; v.showstddevt = 0; v.showstddevf = 0;
v.tablessrel = 0;
end
