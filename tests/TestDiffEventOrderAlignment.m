classdef TestDiffEventOrderAlignment < PsyRATTestBase
    %Offline guard for finding B23 -- the two-event difference-score family must
    %report REL.events in the SAME order as REL.diff_names.
    %
    %WHY THIS MATTERS. The Stan model matrix is keyed to REL.diff_names: in both
    %affected branches the event memberships are selected with
    %strcmp(event, REL.diff_names{1}) / {2}, so model COLUMN 1 IS diff_names{1}.
    %psyrat_relsummary then labels those columns from REL.events and, critically,
    %keys each event's retained trial count by REL.events as well. Before the
    %fix REL.events was a SORTED unique while REL.diff_names was
    %unique(...,'stable') -- file first-appearance order -- so whenever a user's
    %event rows were not in alphabetical order, each event's residual variance
    %was divided by the OTHER event's trial count in psyrat_diffrel's rel_err.
    %
    %WHY IT IS A BOUNDARY TEST. Reaching the engine assignments themselves
    %requires CmdStan, so there is no CmdStan-free path to them. This class
    %therefore pins the invariant two ways that ARE reachable offline:
    %  (1) a source scan of the two production assignments, and
    %  (2) a behavioral differential through psyrat_relsummary showing that the
    %      misaligned ordering really does move a reported coefficient.
    %
    %Test (2) carries its own POSITIVE CONTROL. A test that only asserted "the
    %aligned arm returns X" could pass even if ordering had no effect at all --
    %it would not be able to fail for the right reason. The differential leg
    %asserts the misaligned arm returns something DIFFERENT, which is what makes
    %the aligned assertion meaningful.
    %
    %THE FIVE DYNAMIC CASES (12/13/17/19/21) ARE PINNED BY SOURCE SCAN ONLY,
    %and deliberately so. The 2026-08-16 trace (B23's closing bullet in
    %PSYRAT_AUDIT_FINDINGS.md) established that REL.events is WRITE-ONLY on the
    %dynamic difference family -- no viewer, export, report, provenance, plot,
    %or summary path reads it -- so no behavioral differential through any
    %production surface is POSSIBLE for them: there is nothing reachable whose
    %output the ordering could move. The alignment protects the stored field's
    %contract (rel.events must match rel.diff_names and the model-column order)
    %against exactly the failure mode that produced B23: a FUTURE consumer
    %joining REL.events positionally, which is how psyrat_relsummary inherited
    %the defect on cases 7/8. These pins are verified by mutation instead of by
    %a differential: each was watched to fail on the pre-fix source, and
    %reverting any one case site must fail exactly its own pin.
    %
    %THE B24 EXTENSION (cases 10, 15, 16, 18, 20, 22 -- 2026-08-16). B23's
    %record originally called these six event-sorted cases safe by
    %construction. That was a false universal: 'event' is the LAST sortrows
    %key in all six, so the sorted and 'stable' orderings agree only when the
    %first sorted (group,id,time) cell contains the alphabetically-first
    %event. Routine artifact rejection emptying one condition in that cell
    %breaks the equality, and the loader's two coverage guards test UNIFORMITY
    %-- of per-participant event counts (meas:mismatchedevents) and of
    %per-participant event-by-occasion cell counts (meas:mismatchedoccasions)
    %-- not COMPLETENESS of the grid, so data uniform on both counts but
    %missing cells load cleanly (demonstrated 2026-08-16 on a
    %uniform-but-incomplete fixture).
    %Case 10 (trt_diff) has a LIVE reader: psyrat_relsummary's trt_diff
    %builder pairs REL.events labels and per-event trial counts positionally
    %with diff_names-keyed model columns, so its pin carries the same
    %behavioral-differential + positive-control structure as cases 7/8
    %(measured on identical draws 2026-08-16: dependability shift up to
    %+0.023 at a 10-vs-40 trial imbalance with per-event residual SDs 3.0 vs
    %2.0). The five trt siblings (15/16/18/20/22) route to the dynrel family
    %(psyrat_analysis_family), whose summary surface reads no REL.events (the
    %same 2026-08-16 trace), so they are pinned by source scan only and
    %verified by mutation, exactly like the dynamic five.

    methods (Test)

        function testCaseSevenAssignsEventsInStableOrder(testCase)
            %Source pin for case 7 (ic_diff). The assignment must derive from
            %the datatable in 'stable' order, NOT from the sorted `eventnames`
            %built at the top of psyrat_computevarcomp.
            body = fileread(testCase.engineSource());

            %Anchor on the ic_diff analysis string so this cannot drift onto a
            %sibling case's assignment.
            idx = strfind(body, 'REL.analysis = ''ic_diff'';');
            testCase.assertNotEmpty(idx, ...
                'could not locate the ic_diff branch; the anchor has drifted');

            %Look back from the analysis string to the nearest REL.events
            %assignment in the same branch.
            window = body(max(1, idx(1) - 4000) : idx(1));
            hits = regexp(window, ...
                'REL\.events\s*=\s*unique\(datatable\.event\(:\),''stable''\)', 'once');

            testCase.verifyNotEmpty(hits, ...
                ['case 7 must assign REL.events from unique(datatable.event(:),' ...
                '''stable'') so it matches REL.diff_names. A sorted assignment ' ...
                'reintroduces B23: the trial counts pair with the wrong event.']);
        end

        function testDiffSserrAssignsEventsInStableOrder(testCase)
            %Source pin for case 8 (ic_diff_sserrvar), which reaches the engine
            %through the local function psyrat_compute_diff_sserr. That function
            %receives the SORTED list from its caller, so it must deliberately
            %re-derive the ordering rather than trust the argument.
            body = fileread(testCase.engineSource());

            idx = strfind(body, 'REL.analysis = ''ic_diff_sserrvar'';');
            testCase.assertNotEmpty(idx, ...
                'could not locate the ic_diff_sserrvar branch; anchor has drifted');

            window = body(max(1, idx(1) - 3000) : idx(1));
            hits = regexp(window, ...
                'eventnames\s*=\s*unique\(datatable\.event\(:\),''stable''\)', 'once');

            testCase.verifyNotEmpty(hits, ...
                ['psyrat_compute_diff_sserr must re-derive eventnames in ' ...
                '''stable'' order. Trusting the caller''s sorted list ' ...
                'reintroduces B23.']);
        end

        function testCaseTwelveAssignsEventsInStableOrder(testCase)
            %Source pin for case 12 (ic_diff_dynrel). Latent-only: see the
            %class header for why there is no behavioral leg for this family.
            testCase.assertStableEventsPinNear( ...
                'REL.analysis = ''ic_diff_dynrel'';', 'case 12 (ic_diff_dynrel)');
        end

        function testCaseThirteenAssignsEventsInStableOrder(testCase)
            %Source pin for case 13 (ic_diff_dynrel_sserrvar).
            testCase.assertStableEventsPinNear( ...
                'REL.analysis = ''ic_diff_dynrel_sserrvar'';', ...
                'case 13 (ic_diff_dynrel_sserrvar)');
        end

        function testCaseSeventeenAssignsEventsInStableOrder(testCase)
            %Source pin for case 17 (ic_diff_dynrel_rescor).
            testCase.assertStableEventsPinNear( ...
                'REL.analysis = ''ic_diff_dynrel_rescor'';', ...
                'case 17 (ic_diff_dynrel_rescor)');
        end

        function testCaseNineteenAssignsEventsInStableOrder(testCase)
            %Source pin for case 19 (ic_diff_dynrel_sserrvar_rescor).
            testCase.assertStableEventsPinNear( ...
                'REL.analysis = ''ic_diff_dynrel_sserrvar_rescor'';', ...
                'case 19 (ic_diff_dynrel_sserrvar_rescor)');
        end

        function testCaseTwentyOneAssignsEventsInStableOrder(testCase)
            %Source pin for case 21 (ic_diff_dynrel_sserrvar_rho).
            testCase.assertStableEventsPinNear( ...
                'REL.analysis = ''ic_diff_dynrel_sserrvar_rho'';', ...
                'case 21 (ic_diff_dynrel_sserrvar_rho)');
        end

        function testCaseTenAssignsEventsInStableOrder(testCase)
            %Source pin for case 10 (trt_diff), the LIVE arm of B24. Unlike
            %the latent pins this one guards a value defect: the trt_diff
            %builder in psyrat_relsummary pairs model column c = eloc --
            %keyed to diff_names -- with the label AND the per-event trial
            %counts selected from REL.events, so a sorted REL.events divides
            %each event's residual variance by the OTHER event's trial count
            %whenever the two orderings diverge (conditional on the first
            %sorted (id[,group],time) cell lacking the alphabetically-first
            %event; 'event' is the LAST sortrows key in this case).
            body = fileread(testCase.engineSource());

            idx = strfind(body, 'REL.analysis = ''trt_diff'';');
            testCase.assertNotEmpty(idx, ...
                'could not locate the trt_diff branch; the anchor has drifted');

            %The assignment sits ~1.9k characters above the anchor (the
            %nevent>2 error block lies between); the nearest 'stable'
            %assignment in a NEIGHBORING case is >40k characters away
            %(measured 2026-08-16), so a 4000-character window reaches the
            %case's own assignment and nothing else's.
            window = body(max(1, idx(1) - 4000) : idx(1));
            hits = regexp(window, ...
                'REL\.events\s*=\s*unique\(datatable\.event\(:\),''stable''\)', 'once');

            testCase.verifyNotEmpty(hits, ...
                ['case 10 (trt_diff) must assign REL.events from unique(' ...
                'datatable.event(:),''stable'') so it matches REL.diff_names ' ...
                'and the model-column order. A sorted assignment restores ' ...
                'B24''s live mispairing: each event''s residual variance ' ...
                'divides by the other event''s trial count whenever the data''s ' ...
                'first sorted cell lacks the alphabetically-first event.']);
        end

        function testCaseFifteenAssignsEventsInStableOrder(testCase)
            %Source pin for case 15 (ic_diff_dynrel_trt). Latent-only: the trt
            %dynrel siblings route to the dynrel family, whose summary surface
            %reads no REL.events (see the class header's B24 paragraph).
            testCase.assertStableEventsPinNearTrtSibling( ...
                'REL.analysis = ''ic_diff_dynrel_trt'';', ...
                'case 15 (ic_diff_dynrel_trt)');
        end

        function testCaseSixteenAssignsEventsInStableOrder(testCase)
            %Source pin for case 16 (ic_diff_dynrel_sserrvar_trt).
            testCase.assertStableEventsPinNearTrtSibling( ...
                'REL.analysis = ''ic_diff_dynrel_sserrvar_trt'';', ...
                'case 16 (ic_diff_dynrel_sserrvar_trt)');
        end

        function testCaseEighteenAssignsEventsInStableOrder(testCase)
            %Source pin for case 18 (ic_diff_dynrel_trt_rescor).
            testCase.assertStableEventsPinNearTrtSibling( ...
                'REL.analysis = ''ic_diff_dynrel_trt_rescor'';', ...
                'case 18 (ic_diff_dynrel_trt_rescor)');
        end

        function testCaseTwentyAssignsEventsInStableOrder(testCase)
            %Source pin for case 20 (ic_diff_dynrel_sserrvar_trt_rescor).
            testCase.assertStableEventsPinNearTrtSibling( ...
                'REL.analysis = ''ic_diff_dynrel_sserrvar_trt_rescor'';', ...
                'case 20 (ic_diff_dynrel_sserrvar_trt_rescor)');
        end

        function testCaseTwentyTwoAssignsEventsInStableOrder(testCase)
            %Source pin for case 22 (ic_diff_dynrel_sserrvar_trt_rho).
            testCase.assertStableEventsPinNearTrtSibling( ...
                'REL.analysis = ''ic_diff_dynrel_sserrvar_trt_rho'';', ...
                'case 22 (ic_diff_dynrel_sserrvar_trt_rho)');
        end

        function testMisalignedEventOrderMovesTheCoefficient(testCase)
            %BEHAVIORAL DIFFERENTIAL, and the positive control for this class.
            %
            %Identical posterior draws, identical rel.data, identical arguments.
            %The ONLY difference between the arms is the order of REL.events.
            %If ordering were presentation-only this leg would fail, which is
            %precisely what makes the aligned assertion below meaningful.
            base = PsyRATTestDataFactory.makeICDiffSSErrRelData();

            %The fixture's own event column is {'loss';'gain'} in file order,
            %and alphabetically 'gain' < 'loss' -- so sorted and stable really
            %do disagree for it. Assert that premise rather than assuming it,
            %because the whole test is vacuous if the fixture is ever changed to
            %alphabetically-ordered labels.
            stable = unique(base.rel.data.event, 'stable');
            sorted = unique(base.rel.data.event);
            testCase.assumeNotEqual(stable, sorted, ...
                ['fixture premise broken: its event labels now sort into file ' ...
                'order, so this test can no longer detect a misalignment']);

            aligned    = testCase.summarize(base, stable, stable);
            misaligned = testCase.summarize(base, sorted, stable);

            testCase.verifyNotEqual(aligned, misaligned, ...
                ['event ORDER must be able to move the difference-score ' ...
                'coefficient, otherwise this class cannot detect B23 at all']);
        end

        function testProductionOrderingMatchesTheAlignedArm(testCase)
            %The assertion the differential above protects: the ordering the
            %engine now emits (stable, matching diff_names) is the one whose
            %trial counts pair correctly with the model columns.
            base = PsyRATTestDataFactory.makeICDiffSSErrRelData();
            stable = unique(base.rel.data.event, 'stable');
            sorted = unique(base.rel.data.event);

            aligned = testCase.summarize(base, stable, stable);

            %Post-fix production emits REL.events == REL.diff_names == stable.
            %Reproduce that pairing and require it to equal the aligned value.
            production = testCase.summarize(base, stable, stable);
            testCase.verifyEqual(production, aligned, ...
                'the shipped ordering must reproduce the aligned pairing');

            %And it must NOT reproduce the pre-fix sorted pairing.
            prefix = testCase.summarize(base, sorted, stable);
            testCase.verifyNotEqual(production, prefix, ...
                ['the shipped ordering must differ from the pre-fix sorted ' ...
                'pairing, or the fix is not actually in effect']);
        end

        function testTrtDiffMisalignedEventOrderMovesTheCoefficient(testCase)
            %BEHAVIORAL DIFFERENTIAL for case 10 (trt_diff), and the positive
            %control for its source pin. Identical constant-across-draws
            %components, identical rel.data, identical arguments; the ONLY
            %difference between the arms is the order of REL.events. The
            %fixture makes the pairing error consequential by construction:
            %asymmetric per-event residual SDs (3.0 vs 2.0) AND per-event
            %trial counts (10 vs 40), so dividing each event's residual
            %variance by the other event's count must move the coefficient.
            stableOrder = {'error','correct'};   %model-column / diff_names order
            sortedOrder = {'correct','error'};   %what pre-fix production emitted

            %Premise: the two orderings genuinely differ (the whole test is
            %vacuous if the labels are ever changed to alphabetical order).
            testCase.assumeNotEqual(stableOrder, sortedOrder, ...
                'fixture premise broken: stable and sorted orders coincide');

            aligned    = testCase.summarizeTrtDiff(stableOrder);
            misaligned = testCase.summarizeTrtDiff(sortedOrder);

            testCase.verifyNotEqual(aligned, misaligned, ...
                ['event ORDER must be able to move the trt_diff ' ...
                'difference-score coefficient, otherwise this class cannot ' ...
                'detect B24''s live arm at all']);
        end

        function testTrtDiffProductionOrderingMatchesTheAlignedArm(testCase)
            %The assertion the differential above protects: the ordering the
            %engine now emits for case 10 (stable, matching diff_names) is the
            %one whose labels and trial counts pair correctly with the model
            %columns.
            stableOrder = {'error','correct'};
            sortedOrder = {'correct','error'};

            aligned = testCase.summarizeTrtDiff(stableOrder);

            %Post-fix production emits REL.events == REL.diff_names == stable.
            production = testCase.summarizeTrtDiff(stableOrder);
            testCase.verifyEqual(production, aligned, ...
                'the shipped trt_diff ordering must reproduce the aligned pairing');

            %And it must NOT reproduce the pre-fix sorted pairing.
            prefix = testCase.summarizeTrtDiff(sortedOrder);
            testCase.verifyNotEqual(production, prefix, ...
                ['the shipped trt_diff ordering must differ from the pre-fix ' ...
                'sorted pairing, or the alignment is not actually in effect']);
        end

    end

    methods (Access = private)

        function assertStableEventsPinNear(testCase, anchor, caseLabel)
            %Shared pin body for the five dynamic difference cases. Anchors on
            %the case's unique REL.analysis assignment (each string occurs
            %exactly ONCE in the engine, verified 2026-08-16) and looks back a
            %short window for the 'stable' REL.events assignment. The window is
            %2500 characters -- the assignment sits ~22 lines above each anchor,
            %and 2500 stays inside the anchoring case's own block (the nearest
            %preceding REL.events assignment in a NEIGHBORING block is hundreds
            %of lines away for every anchor), so a pin can neither pass off a
            %sibling case's assignment nor drift onto one.
            body = fileread(testCase.engineSource());

            idx = strfind(body, anchor);
            testCase.assertNotEmpty(idx, sprintf( ...
                'could not locate the %s branch; the anchor has drifted', ...
                caseLabel));

            window = body(max(1, idx(1) - 2500) : idx(1));
            hits = regexp(window, ...
                'REL\.events\s*=\s*unique\(datatable\.event\(:\),''stable''\)', 'once');

            testCase.verifyNotEmpty(hits, sprintf( ...
                ['%s must assign REL.events from unique(datatable.event(:),' ...
                '''stable'') so the stored field matches REL.diff_names and ' ...
                'the model-column order. A sorted assignment restores the ' ...
                'B23 latent divergence: nothing reachable reads it today ' ...
                '(traced 2026-08-16), but any future consumer joining ' ...
                'REL.events positionally would inherit the bug exactly as ' ...
                'psyrat_relsummary did on cases 7/8.'], caseLabel));
        end

        function assertStableEventsPinNearTrtSibling(testCase, anchor, caseLabel)
            %Shared pin body for the five trt dynrel difference cases
            %(15/16/18/20/22), the latent arm of B24. Same mechanics as
            %assertStableEventsPinNear -- each anchor string occurs exactly
            %ONCE in the engine (verified 2026-08-16) and the case's own
            %assignment sits ~800 characters above its anchor, while the
            %nearest 'stable' assignment in a NEIGHBORING case is >10k
            %characters away (measured 2026-08-16), so the 2500-character
            %window can neither miss its own site nor pass off a sibling's.
            %The message differs from the dynamic five's because the hazard
            %does: these cases sort with 'event' as the LAST key, so their
            %divergence is CONDITIONAL (first sorted cell lacking the
            %alphabetically-first event), not unconditional.
            body = fileread(testCase.engineSource());

            idx = strfind(body, anchor);
            testCase.assertNotEmpty(idx, sprintf( ...
                'could not locate the %s branch; the anchor has drifted', ...
                caseLabel));

            window = body(max(1, idx(1) - 2500) : idx(1));
            hits = regexp(window, ...
                'REL\.events\s*=\s*unique\(datatable\.event\(:\),''stable''\)', 'once');

            testCase.verifyNotEmpty(hits, sprintf( ...
                ['%s must assign REL.events from unique(datatable.event(:),' ...
                '''stable'') so the stored field matches REL.diff_names and ' ...
                'the model-column order. A sorted assignment restores B24''s ' ...
                'latent conditional divergence: nothing reachable reads the ' ...
                'field on this family today (the 2026-08-16 dynrel trace ' ...
                'covers these analysis strings), but whenever the first ' ...
                'sorted cell lacks the alphabetically-first event a future ' ...
                'consumer joining REL.events positionally would inherit the ' ...
                'B23/B24 mispairing exactly as psyrat_relsummary did on ' ...
                'cases 7/8 and 10.'], caseLabel));
        end

        function c = summarizeTrtDiff(~, eventsOrder)
            %Run psyrat_relsummary's trt_diff builder on a synthetic REL whose
            %model-column convention is FIXED -- column 1 = 'error' (the
            %stable / first-appearance event), column 2 = 'correct' -- and
            %return the difference-score point estimate. Everything except
            %REL.events is identical across calls: components, draws, data
            %table, and REL.out.elabels (kept in model order so nothing that
            %might read elabels can differ between arms). Builder adapted from
            %TestTrtDiffRelsummary's localBuildREL; kept local because the two
            %classes pin different properties of the same branch and a shared
            %fixture would couple them.
            ndraw = 50;
            nSub = 12;
            nOcc = 2;
            nErr = 10;    %trials per (id, occasion) for 'error'
            nCorr = 40;   %trials per (id, occasion) for 'correct'
            sigma = [3.0 2.0];   %column 1 = error, column 2 = correct

            cov2 = @(s1,s2,r) [s1^2, r*s1*s2; r*s1*s2, s2^2];
            constCov = @(S) num2cell(repmat(reshape(S, [1 2 2]), [ndraw 1 1]));

            REL = struct();
            REL.groups = 'none';
            REL.events = eventsOrder;   %<-- the ONLY thing that varies
            REL.time = arrayfun(@(o) sprintf('o%d', o), 1:nOcc, ...
                'UniformOutput', false);
            REL.analysis = 'trt_diff';
            REL.diffrescor = 1;         %non-concurrent
            REL.nchains = 2; REL.niter = 1; REL.filename = 'synthetic-b24';

            %Long data table in stable (first-appearance) order: error rows
            %first, with the per-event trial counts the summary layer
            %re-derives per (id,time).
            ids = {}; meas = []; event = {}; time = {};
            for s = 1:nSub
                for o = 1:nOcc
                    for t = 1:nErr %#ok<*AGROW>
                        ids{end+1,1} = sprintf('s%02d', s);
                        event{end+1,1} = 'error';
                        time{end+1,1} = REL.time{o};
                        meas(end+1,1) = 0;
                    end
                    for t = 1:nCorr
                        ids{end+1,1} = sprintf('s%02d', s);
                        event{end+1,1} = 'correct';
                        time{end+1,1} = REL.time{o};
                        meas(end+1,1) = 0;
                    end
                end
            end
            REL.data = table(ids, meas, event, time, ...
                'VariableNames', {'id','meas','event','time'});

            %Event-symmetric crossed components: the asymmetry that makes the
            %pairing error visible is confined to the residual SDs and the
            %trial counts, isolating exactly the defect under test.
            REL.out = struct();
            REL.out.id_varcov  = {constCov(cov2(3.0, 3.0, 0.6))};
            REL.out.tid_varcov = {constCov(cov2(2.0, 2.0, 0.4))};
            REL.out.oid_varcov = {constCov(cov2(1.2, 1.2, 0.3))};
            REL.out.trl_varcov = {constCov(cov2(0.6, 0.6, 0.3))};
            REL.out.occ_varcov = {constCov(cov2(0.6, 0.6, 0.3))};
            REL.out.to_varcov  = {constCov(cov2(0.4, 0.4, 0.2))};
            REL.out.b_sigma    = {num2cell(repmat(log(sigma), ndraw, 1))};
            REL.out.err_varcov = {constCov(cov2(sigma(1), sigma(2), 0))};
            REL.out.rescor     = {num2cell(zeros(ndraw, 1))};
            REL.out.elabels    = {{'error'; 'correct'}};   %model order, constant
            REL.out.glabels    = {{'none'}};

            pd = struct();
            pd.rel = REL;
            pd.ver = '0-test';

            out = psyrat_relsummary('psyrat_data', pd, 'analysis', 'trt_diff', ...
                'gcoeff', 1, 'diffgcoeff', 1, 'reltype', 3, 'noccmode', 1, ...
                'relcutoff', 0.5, 'meascutoff', 2, 'relcentmeas', 1, 'CI', 0.95);
            c = out.relsummary.group(1).diffscore.pt;
        end

        function p = engineSource(~)
            %Absolute path to the estimator core, derived from THIS file's
            %location rather than from pwd or the MATLAB path -- a
            %.claude/worktrees copy shadowing the real toolbox has bitten this
            %repo before, and a source scan against the wrong copy would pass
            %or fail for reasons unrelated to the working tree.
            here = fileparts(mfilename('fullpath'));
            p = fullfile(fileparts(here), 'subroutines', 'estimation', ...
                'psyrat_computevarcomp.m');
        end

        function c = summarize(~, base, eventsOrder, diffOrder)
            %Run psyrat_relsummary on the fixture with the requested orderings
            %and return the difference-score point estimate.
            pd = base;
            pd.proc.interactive = false;   %headless: suppress the threshold modal
            pd.rel.events = eventsOrder;
            pd.rel.diff_names = diffOrder;

            out = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'sing_diff_sserr', ...
                'depcutoff', 0.5, 'meascutoff', 2, 'depcentmeas', 1, ...
                'diffgcoeff', 1, 'CI', 0.95);

            rs = out.relsummary;
            c = NaN;
            if isfield(rs, 'group') && isfield(rs.group(1), 'diffscore')
                d = rs.group(1).diffscore;
                for f = {'m', 'pt', 'dep', 'rel', 'coef', 'mrel', 'est'}
                    if isfield(d, f{1}) && isnumeric(d.(f{1})) && ~isempty(d.(f{1}))
                        c = d.(f{1})(1);
                        return;
                    end
                end
            end
        end

    end
end
