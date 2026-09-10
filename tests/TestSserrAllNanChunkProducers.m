classdef TestSserrAllNanChunkProducers < PsyRATTestBase
    %G46 + G49: NaN handling at the sserr chunk producers.
    %
    %G46 (2026-08-29): the producers assigned the good/bad carriers inside
    %`if all(ind2exclude) ... elseif any(ind2include) ... end` with no else.
    %Both indicators came from independent >=/< comparisons, which NaN fails,
    %so an ALL-NaN chunk (a catastrophically failed fit) satisfied neither
    %arm and left the carriers unassigned. The fix added else-arms assigning
    %the canonical dead-event pair.
    %
    %G49 (this batch): the same >=/< independence dropped a PARTIAL-NaN
    %subject from BOTH carriers -- siblings finite, so the any-include arm
    %ran, and the NaN subject landed in neither eventgoodids nor
    %eventbadids, silently breaking the complement the G39/G45 conventions
    %guarantee. All five producer sites now derive the indicators through
    %the shared psyrat_relsummary_mark_inclusion (ind2exclude =
    %~ind2include), so NaN reads as excluded. Consequence for the all-NaN
    %state: it now satisfies all(ind2exclude) and flows the FIRST arm --
    %whose assignments are identical to the G46 else-arms, which remain as
    %documented backstops. The all-NaN tests below therefore weld to the
    %first-arm state (all excluded, none includable); the OUTCOME
    %assertions (dead pair, complement, no abort) are unchanged from G46.
    %
    %Coverage, stated honestly (re-derived after the adversarial wave
    %caught the earlier rationale asserting the all-NaN lanes cover the
    %partial-NaN arm -- they do not: all-NaN takes the all-excluded arm,
    %partial-NaN the any-include arm, DIFFERENT bodies): three of the five
    %sites are driven for all-NaN (case 4, case 2, the trt_sserr storage
    %helper) and two for partial-NaN (case 2 and the helper). The marker
    %itself is one shared function, so what the untested sites could get
    %wrong is only their any-include arm bodies; those at cases 1, 3, and
    %4 are byte-identical copies of case 2's tested body (verified by
    %diff in the wave, not by execution), and the clean-fixture control
    %runs execute them with finite draws.

    methods (Test)

        function testCase4AllNanEventFoldsAsDeadEvent(testCase)
            %The loud pre-fix path: case 4 (groups x events), one event's
            %draws all NaN, so the G39 strict fold used to die on the
            %unassigned eventgoodids (MATLAB:nonExistentField).
            pd = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroupEvent();
            pd.proc = struct('interactive', false);

            %counterfactual control: undoctored, GroupA retains someone at
            %this cutoff, so the zeroing below is attributable to the NaNs
            clean = localSummary(testCase, pd);
            testCase.assertNotEmpty(clean.relsummary.group(1).goodids, ...
                'anti-vacuity: the clean fixture must retain someone in GroupA');

            %doctor GroupA x err (column 2 by the factory's layout) to all-NaN
            testCase.assertEqual(pd.rel.out.labels{2}, 'GroupA_;_err', ...
                'weld to the factory column layout this doctoring relies on');
            sz = size(cell2mat(pd.rel.out.ind_sdlog{2}));
            pd.rel.out.ind_sdlog{2} = num2cell(nan(sz));

            [out, relerr] = localSummaryWithErr(testCase, pd);

            %weld to the post-G49 state: all dep draws NaN, so exclusion
            %(= ~inclusion) covers everyone and the FIRST arm owns the state
            gA = out.relsummary.group(1);
            errloc = find(strcmp({gA.event.name}, 'err'), 1);
            corloc = find(strcmp({gA.event.name}, 'cor'), 1);
            testCase.assertNotEmpty(errloc);
            testCase.assertNotEmpty(corloc);
            T = gA.event(errloc).ssrel_table;
            testCase.assertTrue(all(isnan(T.dep_pt)), ...
                'the doctored chunk must genuinely carry all-NaN dep draws');
            testCase.assertTrue(all(T.ind2exclude), ...
                'weld: post-G49, NaN reads as excluded, so the all-excluded arm covers this state');
            testCase.assertFalse(any(T.ind2include), ...
                'weld: nobody is includable on an all-NaN chunk');

            %the producer assigned the canonical dead-event pair
            testCase.verifyTrue(iscell(gA.event(errloc).eventgoodids) && ...
                isempty(gA.event(errloc).eventgoodids), ...
                'all-NaN chunk: eventgoodids must be the empty cell carrier');
            testCase.verifyEqual(sort(gA.event(errloc).eventbadids), ...
                sort({'s1'; 's2'}), ...
                'all-NaN chunk: eventbadids must carry the full roster');

            %and the G39 strict fold consumed it without aborting: the group
            %zeroes, counts complement, the sibling event stays informative
            testCase.verifyEqual(relerr.nogooddata, 0, ...
                'an all-NaN chunk must not abort the whole summary');
            testCase.verifyTrue(iscell(gA.goodids) && isempty(gA.goodids), ...
                'the dead all-NaN event must empty GroupA''s retained set');
            testCase.verifyEqual(sort(gA.badids), sort({'s1'; 's2'}), ...
                'GroupA good/bad ids must complement');
            testCase.verifyNotEmpty(gA.event(corloc).eventgoodids, ...
                'the sibling event''s own row must stay informative');

            %positive control inside the same run: GroupB is untouched
            testCase.verifyNotEmpty(out.relsummary.group(2).goodids, ...
                'GroupB must be unaffected by GroupA''s dead chunk');
        end

        function testCase2AllNanGroupCarriesDeadPairAndGroupCarriers(testCase)
            %Case 2 (groups only): the producer arms assign the GROUP carriers
            %directly (no fold runs at nevents == 1), so the else-arm must
            %assign those too.
            pd = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
            pd.proc = struct('interactive', false);

            clean = localSummary(testCase, pd);
            testCase.assertNotEmpty(clean.relsummary.group(1).goodids, ...
                'anti-vacuity: the clean fixture must retain someone in GroupA');

            sz = size(cell2mat(pd.rel.out.ind_sdlog{1}));
            pd.rel.out.ind_sdlog{1} = num2cell(nan(sz));

            [out, relerr] = localSummaryWithErr(testCase, pd);

            g1 = out.relsummary.group(1);
            T = g1.event(1).ssrel_table;
            testCase.assertTrue(all(isnan(T.dep_pt)), ...
                'the doctored chunk must genuinely carry all-NaN dep draws');
            testCase.assertTrue(all(T.ind2exclude) && ~any(T.ind2include), ...
                'weld: post-G49 the all-excluded arm covers the all-NaN state');

            testCase.verifyEqual(relerr.nogooddata, 0, ...
                'an all-NaN group must not abort the whole summary');
            testCase.verifyTrue(iscell(g1.event(1).eventgoodids) && ...
                isempty(g1.event(1).eventgoodids));
            testCase.verifyEqual(sort(g1.event(1).eventbadids), ...
                sort({'s1'; 's2'}));
            testCase.verifyTrue(iscell(g1.goodids) && isempty(g1.goodids), ...
                'case 2: the else-arm must assign the group carriers too');
            testCase.verifyEqual(sort(g1.badids), sort({'s1'; 's2'}));
            testCase.verifyEqual(g1.event(1).goodn, 0);

            testCase.verifyNotEmpty(out.relsummary.group(2).goodids, ...
                'GroupB must be unaffected by GroupA''s dead chunk');
        end

        function testTrtSserrHelperAssignsTheDeadPairOnAllNan(testCase)
            %The storage helper (local_store_ss_chunk, trt_sserr route) is
            %the one G46 site with a real BEHAVIOR change rather than a
            %missing-arm fix: its pre-G46 catch-all else assigned EMPTY
            %eventbadids on all-NaN draws (nobody fails the < test either),
            %so the group's good/bad ids did not complement. Post-G46 the
            %all-NaN chunk carries the canonical dead pair and the group
            %complements over the roster.
            pd = PsyRATTestDataFactory.makeTRTSSErrRelDataMultiGroup();
            pd.proc = struct('interactive', false);

            clean = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'trt_sserr', 'gcoeff', 1, 'reltype', 3, ...
                'depcutoff', 0.1, 'meascutoff', 2, 'depcentmeas', 1, ...
                'noccmode', 1);
            testCase.assertNotEmpty(clean.relsummary.group(1).goodids, ...
                'anti-vacuity: the clean fixture must retain someone in GroupA');

            %err_ss = exp(pop_sdlog + ind_sdlog): NaN ind_sdlog for GroupA
            %makes every per-subject dep draw NaN on that chunk
            sz = size(cell2mat(pd.rel.out.ind_sdlog{1}));
            pd.rel.out.ind_sdlog{1} = num2cell(nan(sz));

            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'trt_sserr', 'gcoeff', 1, 'reltype', 3, ...
                'depcutoff', 0.1, 'meascutoff', 2, 'depcentmeas', 1, ...
                'noccmode', 1);

            g1 = out.relsummary.group(1);
            T = g1.event(1).ssrel_table;
            testCase.assertTrue(all(isnan(T.dep_pt)), ...
                'the doctored chunk must genuinely carry all-NaN dep draws');
            testCase.assertTrue(all(T.ind2exclude) && ~any(T.ind2include), ...
                'weld: post-G49 the all-excluded arm covers the all-NaN state');

            testCase.verifyEqual(relerr.nogooddata, 0, ...
                'an all-NaN chunk must not abort the whole summary');
            testCase.verifyTrue(iscell(g1.event(1).eventgoodids) && ...
                isempty(g1.event(1).eventgoodids));
            testCase.verifyEqual(sort(g1.event(1).eventbadids), ...
                sort({'s1'; 's2'}), ...
                ['the dead pair''s full-roster eventbadids IS the behavior ' ...
                'change: the old catch-all assigned an empty list here']);
            testCase.verifyEqual(g1.event(1).goodn, 0);
            testCase.verifyTrue(iscell(g1.goodids) && isempty(g1.goodids));
            testCase.verifyEqual(sort(g1.badids), sort({'s1'; 's2'}), ...
                'the group must complement over the roster');

            testCase.verifyNotEmpty(out.relsummary.group(2).goodids, ...
                'GroupB must be unaffected by GroupA''s dead chunk');
        end

        function testCase2PartialNanSubjectLandsInBadIds(testCase)
            %G49: a subject whose dep draws are NaN while a sibling's are
            %finite used to fail BOTH >=/< comparisons -- the any-include arm
            %ran for the finite sibling and the NaN subject landed in NEITHER
            %carrier, so good/bad no longer complemented over the roster.
            %Post-G49 exclusion is ~inclusion, so the NaN subject must land
            %in the bad carrier.
            pd = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
            pd.proc = struct('interactive', false);

            clean = localSummary(testCase, pd);
            testCase.assertEqual(sort(clean.relsummary.group(1).goodids), ...
                sort({'s1'; 's2'}), ...
                'anti-vacuity: the clean fixture must retain BOTH subjects in GroupA');

            %NaN only s1's per-subject residual draws in GroupA (column 1 =
            %subject 1 by the factory's ndraw x nsubject layout)
            raw = cell2mat(pd.rel.out.ind_sdlog{1});
            raw(:, 1) = NaN;
            pd.rel.out.ind_sdlog{1} = num2cell(raw);

            [out, relerr] = localSummaryWithErr(testCase, pd);
            g1 = out.relsummary.group(1);
            T = g1.event(1).ssrel_table;
            testCase.assertTrue(all(isnan(T.dep_pt(strcmp(T.id, 's1')))) && ...
                all(~isnan(T.dep_pt(strcmp(T.id, 's2')))), ...
                'weld: exactly one subject''s dep draws are NaN');

            testCase.verifyEqual(relerr.nogooddata, 0, ...
                'a partial-NaN chunk must not abort the summary');
            testCase.verifyEqual(g1.event(1).eventgoodids, {'s2'}, ...
                'the finite sibling must stay retained');
            testCase.verifyEqual(g1.event(1).eventbadids, {'s1'}, ...
                ['the partial-NaN subject must land in eventbadids; an empty ' ...
                'list here is the pre-G49 dropped-from-both-carriers state']);
            testCase.verifyEqual(g1.goodids, {'s2'});
            testCase.verifyEqual(g1.badids, {'s1'}, ...
                'the group carriers must complement over the roster');
            testCase.verifyEqual(g1.event(1).goodn, 1);
        end

        function testTrtSserrHelperPartialNanComplements(testCase)
            %G49 on the trt_sserr storage helper route: same partial-NaN
            %drop, fixed through the same shared marker. The caller's
            %intersection must fold the NaN subject into badids, not lose it.
            pd = PsyRATTestDataFactory.makeTRTSSErrRelDataMultiGroup();
            pd.proc = struct('interactive', false);

            clean = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'trt_sserr', 'gcoeff', 1, 'reltype', 3, ...
                'depcutoff', 0.1, 'meascutoff', 2, 'depcentmeas', 1, ...
                'noccmode', 1);
            testCase.assertEqual(sort(clean.relsummary.group(1).goodids), ...
                sort({'s1'; 's2'}), ...
                'anti-vacuity: the clean fixture must retain BOTH subjects in GroupA');

            raw = cell2mat(pd.rel.out.ind_sdlog{1});
            raw(:, 1) = NaN;
            pd.rel.out.ind_sdlog{1} = num2cell(raw);

            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'trt_sserr', 'gcoeff', 1, 'reltype', 3, ...
                'depcutoff', 0.1, 'meascutoff', 2, 'depcentmeas', 1, ...
                'noccmode', 1);
            g1 = out.relsummary.group(1);
            T = g1.event(1).ssrel_table;
            testCase.assertTrue(all(isnan(T.dep_pt(strcmp(T.id, 's1')))) && ...
                all(~isnan(T.dep_pt(strcmp(T.id, 's2')))), ...
                'weld: exactly one subject''s dep draws are NaN');

            testCase.verifyEqual(relerr.nogooddata, 0, ...
                'a partial-NaN chunk must not abort the summary');
            testCase.verifyEqual(g1.event(1).eventgoodids, {'s2'});
            testCase.verifyEqual(g1.event(1).eventbadids, {'s1'}, ...
                'the partial-NaN subject must land in eventbadids');
            testCase.verifyEqual(g1.goodids, {'s2'});
            testCase.verifyEqual(g1.badids, {'s1'});
        end

    end
end

%--------------------------------------------------------------------------
function out = localSummary(testCase, pd) %#ok<INUSD>
out = psyrat_relsummary( ...
    'psyrat_data', pd, ...
    'analysis', 'sing_sserr', ...
    'depcutoff', 0.1, ...
    'meascutoff', 2, ...
    'depcentmeas', 1, ...
    'CI', 0.95);
end

function [out, relerr] = localSummaryWithErr(testCase, pd) %#ok<INUSD>
[out, relerr] = psyrat_relsummary( ...
    'psyrat_data', pd, ...
    'analysis', 'sing_sserr', ...
    'depcutoff', 0.1, ...
    'meascutoff', 2, ...
    'depcentmeas', 1, ...
    'CI', 0.95);
end
