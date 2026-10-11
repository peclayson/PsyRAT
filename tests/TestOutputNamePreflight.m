classdef TestOutputNamePreflight < PsyRATTestBase
    %Row 31 (LA-2): a supported request whose derived file names the save
    %folder cannot hold passed every check, reached estimation, and failed
    %only at the save after the fit. A no-fit probe on 2026-09-23 established
    %this on one APFS volume. There, every accepted name had at most 255
    %UTF-16 units and every refused one at least 256, for one- to four-byte
    %UTF-8 characters alike: three-byte names of 735 to 751 bytes were
    %accepted, and four-byte names (two units each) were refused one character
    %past 131 to 135 code points. So that volume's limit was not 255 UTF-8
    %bytes or 255 code points. These are measurements on that volume, not
    %constants of the check. The four observed failures:
    %  - a 247-character measurement header under the typed name p.psyrat
    %    gives the 256-unit output p_<header>.psyrat;
    %  - psyrat_run with a 256-unit savename fails the same way;
    %  - a name whose output fits but whose _runconfig.json sidecar does not
    %    saved the output and lost the sidecar behind one console line;
    %  - the batch summary said the failed runs' reports were "named after
    %    their output files" when the report had fallen back to a timestamped
    %    name.
    %
    %CONTRACT UNDER TEST (root-approved 2026-09-23).
    %  - Before estimation starts for a request, the actual volume must be
    %    able to create its output, its <stem>_runconfig.json sidecar and its
    %    own-name <stem>_error.txt report. There is no universal byte cap, and
    %    no name is truncated or replaced.
    %  - A request that fails is refused with a psyrat_data:save* error
    %    naming the request, the failing file and the remedy. A folder
    %    occupying the output or sidecar path is refused the same way, and
    %    its contents are kept.
    %  - Refusal is per member, the existing batch policy: valid siblings
    %    still run.
    %  - Existing consent, reservation and fallback rules are unchanged: the
    %    consented overwrite target and its sidecar may be rewritten, G93
    %    reservations still apply, and a data error still comes first with
    %    the G94 report fallbacks. A folder at the report path uses those
    %    fallbacks and is no new restriction.
    %  - Nothing that existed before is removed or rewritten otherwise, and no
    %    probe entry is left behind.
    %  - The batch summary names the first five report files actually written
    %    and counts the rest, instead of claiming they are named after the
    %    outputs.
    %  - A passing check is pre-fit evidence only. It does not promise that
    %    every later disk failure is prevented.
    %The helper's own safeguards (probe ownership, cleanup honesty, random
    %number generator untouched, equal-length disjoint probe names) are
    %tested in TestCheckOutputNames.
    %
    %OBSERVED BASELINE (root, 2026-09-23, unfixed source): 5 PASS, 6 FAIL,
    %1 Incomplete. The six expected defects reproduced for their intended
    %reasons: estimation was reached; MATLAB's couldNotWriteFile and
    %filenameIsADirectory errors came back instead of the refusal; the
    %sidecar case wrote a surplus output; the summary named no report and
    %made a false naming claim. The Incomplete was an assertion on the
    %sidecar case's missing error, which is the defect itself. It is now a
    %verification with a guard, so the same source gives 6 FAIL and
    %0 Incomplete, with the same six failure reasons.
    %
    %GROUND TRUTH IS MEASURED, NEVER ASSUMED. The class measures this
    %volume's component limit with the pipeline's own writer forms in a
    %task-owned folder. Every request's validity is then asked of the
    %filesystem directly (OutputNameFixture.oracle), not taken from a
    %constant, so a byte-limited or unlimited volume changes the expected
    %outcome, not the test's truth. No test skips. If no limit exists up to
    %the search bound, the long-name tests check that the long names run and
    %save (every name is writable there). The folder-at-the-output-path test
    %is the deterministic, platform-independent refusal seam, and its
    %delivery control proves the writer really fails.
    %
    %NO FIT: the estimation core, the convergence check and the CmdStan
    %detector are shadowed, and which() proves it before every drive.
    %
    %On the unfixed code (RED as observed; every test is expected to pass once
    %psyrat_check_output_names is wired in):
    %  RED     testLongHeaderOutputIsRefusedBeforeEstimation       (GUI)
    %  RED     testScriptedOutputPastTheLimitIsRefusedBeforeEstimation
    %  RED     testSidecarPastTheLimitIsRefusedEvenWhenTheOutputFits
    %  RED     testFolderAtTheOutputPathIsRefusedBeforeEstimation
    %  RED     testSummaryNamesTheFallbackReport                    (GUI)
    %  RED     testSummaryNamesAnOwnNameReport                      (GUI)
    %  control testFiftyCharacterHeadersStillRun                    (GUI)
    %  control testNamesExactlyAtTheLimitStillRun                   (GUI)
    %  control testLongValidUnicodeNameIsSaved
    %  control testReservedAndPreexistingEntriesSurviveABatch       (GUI)
    %  control testConsentedTargetsAreUntouchedBeforeEstimation
    %Added after review (contract combinations, not RED requirements; not
    %run against the unfixed source):
    %          testDataErrorStaysFirstWhenTheSameMemberIsAlsoTooLong (GUI)
    %          testSummaryListsTheFirstFiveReportsAndCountsTheRest   (GUI)
    %          testMultibyteNamesFollowTheVolumeNotAByteCap          (GUI)
    %          (observed RED on the phase-3 source, 35 PASS, 1 FAIL: the
    %          check's test folder hit MATLAB:FILE:PATHMAX for a name the
    %          writers accept; fixed in phase 5 by confirming a refused test
    %          folder with a test file, truth unchanged)

    properties (Access = private)
        Lim
        OracleDir
    end

    methods (TestClassSetup)
        function measureTheVolume(testCase)
            testCase.assertFalse(any(isspace(tempdir)), ...
                ['these drives need a whitespace-free temp path: the save-dialog ' ...
                'shadow would answer the same whitespace path forever']);
            folder = tempname;
            mkdir(folder);
            testCase.addTeardown(@() rmdir(folder, 's'));
            testCase.OracleDir = folder;
            testCase.Lim = OutputNameFixture.measureLimits(folder);
            testCase.log(1, sprintf(['Measured in %s: output (save) component limit %d ' ...
                'characters (found %d), sidecar (fopen) %d (found %d).'], folder, ...
                testCase.Lim.save, testCase.Lim.saveFound, ...
                testCase.Lim.sidecar, testCase.Lim.sidecarFound));
        end
    end

    methods (Test)

        % ---- refusals (RED on the unfixed code) ------------------------------

        function testLongHeaderOutputIsRefusedBeforeEstimation(testCase)
            %The supported header route. m_short is an ordinary member; the
            %long header makes p_<header>.psyrat one character past the
            %measured output limit.
            long = OutputNameFixture.header(testCase.Lim.save - 8, 'b');
            members = [testCase.headerMember('m_short', 'm_short'), ...
                testCase.headerMember('long header', long)];
            testCase.premise(members, [true false], ...
                'the long member''s output is one character past the limit');
            [ctx, obs, before, after] = testCase.runGui({'m_short', long}, 0);
            OutputNameFixture.verifyOutcome(testCase, ctx, obs, members, before, after, 'gui');
            testCase.verifyTrue(contains(obs.dialogs.complete, testCase.counts(members)), ...
                'the batch must finish with every writable member saved and the refused one counted as failed');
        end

        function testScriptedOutputPastTheLimitIsRefusedBeforeEstimation(testCase)
            %The scripted twin: psyrat_run with a savename one character past
            %the measured output limit.
            savename = [repmat('r', 1, testCase.Lim.save - 6) '.psyrat'];
            m = testCase.scriptedMember('scripted output', savename);
            testCase.premise(m, false, 'the scripted output is one character past the limit');
            [ctx, obs, before, after] = testCase.runScripted(savename);
            OutputNameFixture.verifyOutcome(testCase, ctx, obs, m, before, after, 'script');
        end

        function testSidecarPastTheLimitIsRefusedEvenWhenTheOutputFits(testCase)
            %The output name fits; its <stem>_runconfig.json is one character
            %past the measured sidecar limit. Before the fix the run saved the
            %output and lost the sidecar behind a console line.
            stem = repmat('s', 1, testCase.Lim.sidecar - numel('_runconfig.json') + 1);
            savename = [stem '.psyrat'];
            [ok, detail] = OutputNameFixture.oracle(testCase.OracleDir, savename);
            m = OutputNameFixture.member('sidecar past the limit', savename, ok);
            if testCase.limitsFound()
                testCase.assertTrue(detail.output, 'premise: the output name alone must fit');
                testCase.assertFalse(detail.sidecar, 'premise: the sidecar must be one character past the limit');
            end
            [ctx, obs, before, after] = testCase.runScripted(savename);
            OutputNameFixture.verifyOutcome(testCase, ctx, obs, m, before, after, 'script');
            %wording (review L1): the output itself fits here, so the refusal
            %names the sidecar and must not claim the output cannot be saved
            if ~m.valid && ~isempty(obs.escaped)
                testCase.verifyTrue(contains(obs.escaped.message, [stem '_runconfig.json']), ...
                    'the refusal must name the sidecar, the file that cannot be created');
                testCase.verifyFalse(contains(obs.escaped.message, 'cannot be saved'), ...
                    'the output itself fits, so the refusal must not say it cannot be saved');
            end
        end

        function testFolderAtTheOutputPathIsRefusedBeforeEstimation(testCase)
            %The deterministic, platform-independent seam: a FOLDER occupies
            %the exact output path, so no writer can create the output there,
            %whatever the volume's name limit. The folder and its content must
            %survive, and estimation must not start.
            ctx = testCase.freshDrive();
            occupied = fullfile(ctx.save, 'occ.psyrat');
            mkdir(occupied);
            fid = fopen(fullfile(occupied, 'inside.txt'), 'w');
            fprintf(fid, 'a folder that happens to carry the output name, do not touch');
            fclose(fid);
            %delivery control: the writer really cannot write a file there
            probe = fopen(occupied, 'w');
            if probe >= 3; fclose(probe); end
            testCase.assertEqual(probe, -1, 'premise: a folder must refuse a file write under its own name');
            m = OutputNameFixture.member('folder at the output path', 'occ.psyrat', false);
            csv = OutputNameFixture.writeFixture(ctx, {'meas1'}, 0);
            before = OutputNameFixture.snapshot(ctx.root);
            obs = OutputNameFixture.driveScripted(testCase, ctx, csv, 'occ.psyrat');
            after = OutputNameFixture.snapshot(ctx.root);
            OutputNameFixture.verifyOutcome(testCase, ctx, obs, m, before, after, 'script');
        end

        function testSummaryNamesTheFallbackReport(testCase)
            %The refused member's own report name (<stem>_error.txt) is too
            %long as well, so its report falls back to psyrat_error_<time>.txt.
            %The summary must name that file and must not claim the reports
            %are named after the outputs.
            long = OutputNameFixture.header(testCase.Lim.save - 8, 'b');
            members = [testCase.headerMember('m_short', 'm_short'), ...
                testCase.headerMember('long header', long)];
            testCase.premise(members, [true false], 'one member is past the limit');
            [~, obs, before, after] = testCase.runGui({'m_short', long}, 0);
            reports = OutputNameFixture.newReports(before, after);
            if ~testCase.limitsFound()
                testCase.verifyEmpty(reports, 'no member failed, so no report may appear');
                return;
            end
            testCase.assertNotEmpty(reports, 'premise: the refused member must leave a report');
            txt = obs.dialogs.complete;
            testCase.assertNotEmpty(txt, 'premise: the batch must end in the Processing Complete summary');
            testCase.verifySummaryNamesReports(txt, reports);
        end

        function testSummaryNamesAnOwnNameReport(testCase)
            %A member that fails for its data (NA in its column only) with an
            %ordinary name: its report is <stem>_error.txt. The summary must
            %name it too; the listing is general, not a fallback special case.
            [~, obs, before, after] = testCase.runGui({'mA', 'mB'}, 2);
            reports = OutputNameFixture.newReports(before, after);
            testCase.assertNumElements(reports, 1, 'premise: only the member with NA fails');
            testCase.assertEqual(reports(1).name, 'p_mB_error.txt', ...
                'premise: the report takes its own name');
            testCase.assertTrue(contains(reports(1).identifier, 'missingcells'), ...
                'premise: the member failed for its data, in the preflight');
            txt = obs.dialogs.complete;
            testCase.assertNotEmpty(txt, 'premise: the batch must end in the Processing Complete summary');
            testCase.verifyTrue(contains(txt, 'Successful: 1 Failed: 1'), ...
                'the other member must still run');
            testCase.verifySummaryNamesReports(txt, reports);
        end

        % ---- approved-contract combinations (review T1, T2, T4) -------------
        %These add no product behavior: they pin clauses of the approved
        %contract that the tests above reach only one at a time.

        function testDataErrorStaysFirstWhenTheSameMemberIsAlsoTooLong(testCase)
            %T1. The long-header member also has NA in its own column. The
            %data error keeps the first identifier (the GUI raises every
            %preflight error under the first one's), the name refusal is
            %reported with it, estimation never starts for that member, and
            %the sibling still runs.
            long = OutputNameFixture.header(testCase.Lim.save - 8, 'b');
            short = testCase.headerMember('m_short', 'm_short');
            bad = testCase.headerMember('long header with NA', long);
            testCase.premise([short bad], [true false], ...
                'the long member''s output is one character past the limit');
            [ctx, obs, before, after] = testCase.runGui({'m_short', long}, 2);
            testCase.verifyEqual(obs.engine, {short.requested}, ...
                'estimation must start for the sibling only');
            testCase.verifyTrue(isfile(fullfile(ctx.save, short.requested)), ...
                'the sibling''s output must be saved');
            testCase.verifyTrue(isfile(fullfile(ctx.save, [short.stem '_runconfig.json'])), ...
                'and its sidecar');
            reports = OutputNameFixture.newReports(before, after);
            testCase.assertNumElements(reports, 1, 'exactly one report, for the member with NA');
            r = reports(1);
            testCase.verifyEqual(r.requested, bad.requested);
            testCase.verifyEqual(r.identifier, 'varargin:missingcells', ...
                'the data error must keep the first identifier');
            testCase.verifyTrue(contains(r.message, 'missing cells'), ...
                'the data error must be in the report');
            if testCase.limitsFound()
                testCase.verifyTrue(contains(r.message, bad.requested), ...
                    'the name refusal must be reported beside it');
            end
            testCase.verifyTrue(contains(obs.dialogs.complete, 'Successful: 1 Failed: 1'));
            OutputNameFixture.verifyNothingElse(testCase, before, after, ...
                {fullfile(ctx.save, short.requested), ...
                fullfile(ctx.save, [short.stem '_runconfig.json']), r.path});
        end

        function testSummaryListsTheFirstFiveReportsAndCountsTheRest(testCase)
            %T2. Seven members fail for their data and one runs. The summary
            %lists the first five reports in processing order and counts the
            %other two; the console names all seven.
            fails = arrayfun(@(k) sprintf('f%d', k), 1:7, 'UniformOutput', false);
            [ctx, obs, before, after] = testCase.runGui([{'mok'}, fails], 2:8);
            reports = OutputNameFixture.newReports(before, after);
            expected = cellfun(@(f) ['p_' f '_error.txt'], fails, 'UniformOutput', false);
            testCase.assertEqual(sort({reports.name}), sort(expected), ...
                'premise: each failing member writes its own-name report in the save folder');
            txt = obs.dialogs.complete;
            testCase.assertNotEmpty(txt, 'premise: the batch must end in the Processing Complete summary');
            testCase.verifyTrue(contains(txt, 'Successful: 1 Failed: 7'));
            listed = cellfun(@(n) contains(txt, n), expected);
            testCase.verifyEqual(listed, [true(1, 5) false(1, 2)], ...
                'the summary must list exactly the first five reports, in processing order');
            testCase.verifyTrue(contains(txt, 'and 2 more'), 'the remainder must be counted');
            for k = 1:numel(reports)
                testCase.verifyTrue(contains(obs.console, ['Wrote error report: ' reports(k).path]), ...
                    sprintf('the console must name every report (%s)', reports(k).name));
            end
            testCase.verifyEqual(obs.engine, {'p_mok.psyrat'}, 'only the member without NA is estimated');
            OutputNameFixture.verifyNothingElse(testCase, before, after, ...
                [{fullfile(ctx.save, 'p_mok.psyrat'), fullfile(ctx.save, 'p_mok_runconfig.json')}, ...
                {reports.path}]);
        end

        function testMultibyteNamesFollowTheVolumeNotAByteCap(testCase)
            %T4. A typed base of 200 two-byte characters. The member with a
            %short header is over 255 UTF-8 bytes but within this volume's
            %capacity and must be saved; the member whose header takes its
            %output one character past the measured limit must be refused
            %before its fit. Each name's validity is the oracle's, measured
            %on this volume, never a fixed byte count.
            %Observed RED (root, 2026-09-23, phase-3 source; 35 PASS, 1 FAIL):
            %the writable member was refused before its fit. Root's diagnostic
            %(adjudication H1) traced it to the check's test FOLDER: MATLAB's
            %mkdir refused the equal-length multibyte folder with
            %MATLAB:FILE:PATHMAX at this depth, while fopen and save wrote the
            %files. The truth of this test is unchanged; the check now
            %confirms a refused test folder with a test file.
            %Fixture (review F7). The names observed RED are kept wherever the
            %writers confirm them, as on the measured APFS volume, so that RED
            %is still detected there. Where they do not (a volume that does
            %not count e-acute the way the ASCII limit suggests), both members
            %come from the writers themselves: the longest base whose
            %short-header member they write, then the longest header they
            %write after it; one header character more is the invalid twin.
            %Either way the measured premise is logged and both members are
            %asserted before the drive.
            base = repmat(char(233), 1, 200);
            longh = OutputNameFixture.header(testCase.Lim.save - numel(base) - 7, 'b');
            valid = testCase.scriptedMember('multibyte, short header', [base '_mv.psyrat']);
            invalid = testCase.scriptedMember('multibyte, past the limit', [base '_' longh '.psyrat']);
            observed = ~testCase.limitsFound() || isequal([valid.valid invalid.valid], [true false]);
            branch = 'the names observed RED';
            if ~observed
                branch = 'names measured by the writers';
                n = OutputNameFixture.largestValid(testCase.OracleDir, ...
                    @(k) [repmat(char(233), 1, k) '_mv.psyrat'], 1, numel(base));
                testCase.assertGreaterThanOrEqual(n, 1, ...
                    'premise: the writers here must hold a short-header name with e-acute');
                base = repmat(char(233), 1, n);
                [h, found] = OutputNameFixture.largestValid(testCase.OracleDir, ...
                    @(k) [base '_' OutputNameFixture.header(k, 'b') '.psyrat'], 2, OutputNameFixture.Kmax);
                testCase.assertTrue(found, ...
                    'premise: the writers here must refuse some header length after this base');
                longh = OutputNameFixture.header(h + 1, 'b');
                valid = testCase.scriptedMember('multibyte, short header', [base '_mv.psyrat']);
                invalid = testCase.scriptedMember('multibyte, past the limit', [base '_' longh '.psyrat']);
            end
            testCase.log(1, sprintf(['T4 fixture (%s): base %d e-acute; valid output %d UTF-16 ' ...
                'units, %d UTF-8 bytes, writers %d; invalid twin %d units, %d bytes, writers %d.'], ...
                branch, numel(base), numel(valid.requested), ...
                numel(unicode2native(valid.requested, 'UTF-8')), valid.valid, ...
                numel(invalid.requested), numel(unicode2native(invalid.requested, 'UTF-8')), ...
                invalid.valid));
            if observed
                testCase.assertGreaterThan(numel(unicode2native(valid.requested, 'UTF-8')), 255, ...
                    'premise: the writable name is over 255 UTF-8 bytes');
            end
            testCase.premise([valid invalid], [true false], ...
                'the long member is one character past the limit');
            ctx = testCase.freshDrive();
            %Positive-control pin at the save folder's own depth. The class
            %oracle's folder is shallower (not canonicalized), so the same
            %question is asked again, independently of the product, in a
            %sibling of the save folder whose name has the same length; it
            %must agree with the class oracle, or the premise does not hold at
            %the depth the run writes to.
            twin = fullfile(ctx.root, 'twin');
            mkdir(twin);
            sameDepth = [OutputNameFixture.oracle(twin, valid.requested), ...
                OutputNameFixture.oracle(twin, invalid.requested)];
            rmdir(twin);
            testCase.assertEqual(sameDepth, [valid.valid invalid.valid], ...
                'premise: the oracle at the save folder''s depth must agree with the class oracle');
            headers = {'mv', longh};
            csv = OutputNameFixture.writeFixture(ctx, headers, 0);
            prefs = OutputNameFixture.guiPrefs(testCase, csv, headers, 0);
            before = OutputNameFixture.snapshot(ctx.root);
            obs = OutputNameFixture.driveGui(testCase, ctx, csv, prefs, [base '.psyrat']);
            after = OutputNameFixture.snapshot(ctx.root);
            OutputNameFixture.verifyOutcome(testCase, ctx, obs, [valid invalid], before, after, 'gui');
            testCase.verifyTrue(contains(obs.dialogs.complete, testCase.counts([valid invalid])));
        end

        % ---- controls (pass before and after the fix) -------------------------

        function testFiftyCharacterHeadersStillRun(testCase)
            %The accepted 50-character pin of the Analyze drive tests, driven
            %through estimation: both members run and save under their full
            %names. It is also the delivery control for the shadows: the
            %engine stub saw both names, and the detector stub answered.
            pin1 = ['m' repmat('a', 1, 48) 'b'];
            pin2 = ['m' repmat('a', 1, 48) 'c'];
            members = [testCase.headerMember('pin b', pin1), testCase.headerMember('pin c', pin2)];
            testCase.assertEqual([members.valid], [true true], ...
                'premise: 59-character names must be writable here');
            [ctx, obs, before, after] = testCase.runGui({pin1, pin2}, 0);
            OutputNameFixture.verifyOutcome(testCase, ctx, obs, members, before, after, 'gui');
            testCase.verifyTrue(contains(obs.dialogs.complete, 'Successful: 2 Failed: 0'));
            testCase.verifyGreaterThanOrEqual(obs.detect, 1, ...
                'delivery: the CmdStan detector stub must have answered (no stanc ran)');
        end

        function testNamesExactlyAtTheLimitStillRun(testCase)
            %Off-by-one guard: each sidecar is exactly at the measured limit,
            %so every derived name is legal and nothing may be refused.
            h = testCase.Lim.sidecar - numel('p_') - numel('_runconfig.json');
            hb = OutputNameFixture.header(h, 'b');
            hc = OutputNameFixture.header(h, 'c');
            members = [testCase.headerMember('edge b', hb), testCase.headerMember('edge c', hc)];
            if testCase.limitsFound()
                testCase.assertEqual([members.valid], [true true], 'premise: names at the limit are writable');
                [~, over] = OutputNameFixture.oracle(testCase.OracleDir, ...
                    ['p_' OutputNameFixture.header(h + 1, 'b') '.psyrat']);
                testCase.assertFalse(over.sidecar, ...
                    'premise: one more character must not fit, so this pair sits exactly at the limit');
            end
            [ctx, obs, before, after] = testCase.runGui({hb, hc}, 0);
            OutputNameFixture.verifyOutcome(testCase, ctx, obs, members, before, after, 'gui');
            testCase.verifyTrue(contains(obs.dialogs.complete, testCase.counts(members)));
        end

        function testLongValidUnicodeNameIsSaved(testCase)
            %No universal byte cap: 128 two-byte characters make an output of
            %more than 255 UTF-8 bytes but well under the unit limit. Where
            %the filesystem holds it, and APFS does, the run must save it.
            savename = [repmat(char(233), 1, 128) '_u.psyrat'];
            testCase.assertGreaterThan(numel(unicode2native(savename, 'UTF-8')), 255, ...
                'premise: the name is over 255 UTF-8 bytes');
            m = testCase.scriptedMember('multibyte name', savename);
            testCase.log(1, sprintf('multibyte request writable here: %d', m.valid));
            [ctx, obs, before, after] = testCase.runScripted(savename);
            OutputNameFixture.verifyOutcome(testCase, ctx, obs, m, before, after, 'script');
        end

        function testReservedAndPreexistingEntriesSurviveABatch(testCase)
            %G93 reservation around pre-existing entries: a zero-byte file, a
            %case variant of a derived name, a folder carrying an output name
            %(with content), an unrelated decoy, and the typed base itself.
            %Each must be untouched when estimation starts and after the run,
            %and the batch must take G93's documented suffixes.
            ctx = testCase.freshDrive();
            S = ctx.save;
            zero = {fullfile(S, 'p.psyrat'), fullfile(S, 'p_meas1.psyrat')};
            for i = 1:numel(zero)
                fid = fopen(zero{i}, 'w'); fclose(fid);
            end
            content = {fullfile(S, 'P_MEAS1_2_RUNCONFIG.JSON'), fullfile(S, 'p_meas2_2_error_2.txt')};
            for i = 1:numel(content)
                fid = fopen(content{i}, 'w'); fprintf(fid, 'planted %d, do not touch', i); fclose(fid);
            end
            folder = fullfile(S, 'p_meas2.psyrat');
            mkdir(folder);
            fid = fopen(fullfile(folder, 'inside.txt'), 'w'); fprintf(fid, 'inside a planted folder'); fclose(fid);
            planted = [zero, content, {folder, fullfile(folder, 'inside.txt')}];
            %G93: p_meas1 is taken (zero-byte output) and so is p_meas1_2 (its
            %sidecar, case-insensitively) -> p_meas1_3; p_meas2 is taken by
            %the folder -> p_meas2_2 (the decoy is not one of its names)
            members = [testCase.scriptedMember('meas1', 'p_meas1_3.psyrat'), ...
                testCase.scriptedMember('meas2', 'p_meas2_2.psyrat')];
            testCase.assertEqual([members.valid], [true true], 'premise: both allocated names are writable');
            csv = OutputNameFixture.writeFixture(ctx, {'meas1', 'meas2'}, 0);
            prefs = OutputNameFixture.guiPrefs(testCase, csv, {'meas1', 'meas2'}, 0);
            OutputNameFixture.watch(planted);
            before = OutputNameFixture.snapshot(ctx.root);
            obs = OutputNameFixture.driveGui(testCase, ctx, csv, prefs, 'p.psyrat');
            after = OutputNameFixture.snapshot(ctx.root);
            OutputNameFixture.verifyOutcome(testCase, ctx, obs, members, before, after, 'gui');
            testCase.verifyWatchedUntouchedAtEstimation();
        end

        function testConsentedTargetsAreUntouchedBeforeEstimation(testCase)
            %Single-request overwrite: the consented output exists as a
            %zero-byte file and a case variant of its sidecar exists. The run
            %may rewrite both afterwards (existing behavior), but nothing may
            %touch them, or the unrelated decoy, before estimation starts.
            ctx = testCase.freshDrive();
            S = ctx.save;
            target = fullfile(S, 'keep.psyrat');
            fid = fopen(target, 'w'); fclose(fid);
            sidecar = fullfile(S, 'KEEP_runconfig.json');
            fid = fopen(sidecar, 'w'); fprintf(fid, 'an earlier sidecar'); fclose(fid);
            decoy = fullfile(S, 'keep_error_2.txt');
            fid = fopen(decoy, 'w'); fprintf(fid, 'a decoy, do not touch'); fclose(fid);
            m = testCase.scriptedMember('consented target', 'keep.psyrat');
            testCase.assertTrue(m.valid, 'premise: the name is writable');
            csv = OutputNameFixture.writeFixture(ctx, {'meas1'}, 0);
            OutputNameFixture.watch({target, sidecar, decoy});
            before = OutputNameFixture.snapshot(ctx.root);
            obs = OutputNameFixture.driveScripted(testCase, ctx, csv, 'keep.psyrat');
            after = OutputNameFixture.snapshot(ctx.root);
            OutputNameFixture.verifyOutcome(testCase, ctx, obs, m, before, after, 'script', ...
                {target, sidecar});
            testCase.verifyWatchedUntouchedAtEstimation();
        end

    end

    methods (Access = private)

        function ctx = freshDrive(testCase)
            %a task-owned folder tree (canonical, so every path compares
            %exactly), the shadows, and teardown for both
            root = tempname;
            mkdir(root);
            testCase.addTeardown(@() rmdir(root, 's'));
            root = char(java.io.File(root).getCanonicalPath());
            ctx = OutputNameFixture.install(testCase, root);
        end

        function [ctx, obs, before, after] = runGui(testCase, headers, nancol)
            ctx = testCase.freshDrive();
            csv = OutputNameFixture.writeFixture(ctx, headers, nancol);
            prefs = OutputNameFixture.guiPrefs(testCase, csv, headers, nancol);
            before = OutputNameFixture.snapshot(ctx.root);
            obs = OutputNameFixture.driveGui(testCase, ctx, csv, prefs, 'p.psyrat');
            after = OutputNameFixture.snapshot(ctx.root);
        end

        function [ctx, obs, before, after] = runScripted(testCase, savename)
            ctx = testCase.freshDrive();
            csv = OutputNameFixture.writeFixture(ctx, {'meas1'}, 0);
            before = OutputNameFixture.snapshot(ctx.root);
            obs = OutputNameFixture.driveScripted(testCase, ctx, csv, savename);
            after = OutputNameFixture.snapshot(ctx.root);
        end

        function m = headerMember(testCase, label, header)
            %a batch member under the typed name p.psyrat in a fresh folder
            m = testCase.scriptedMember(label, ['p_' header '.psyrat']);
        end

        function m = scriptedMember(testCase, label, requested)
            m = OutputNameFixture.member(label, requested, ...
                OutputNameFixture.oracle(testCase.OracleDir, requested));
        end

        function tf = limitsFound(testCase)
            tf = testCase.Lim.saveFound && testCase.Lim.sidecarFound;
        end

        function premise(testCase, members, expectedValid, what)
            %The intended boundary, asserted where this volume has a limit.
            %Where it has none, every name is writable, and the drive checks
            %that outcome instead (the oracle already says so).
            if testCase.limitsFound()
                testCase.assertEqual([members.valid], expectedValid, ...
                    ['premise, measured on this volume: ' what]);
            else
                testCase.log(1, sprintf(['No component-name limit was found up to %d ' ...
                    'characters here, so "%s" cannot be built; the drive checks that ' ...
                    'these names run and save.'], OutputNameFixture.Kmax, what));
            end
        end

        function s = counts(~, members)
            s = sprintf('Successful: %d Failed: %d', nnz([members.valid]), nnz(~[members.valid]));
        end

        function verifySummaryNamesReports(testCase, txt, reports)
            for r = reports(:)'
                testCase.verifyTrue(contains(txt, r.name), sprintf( ...
                    'the summary must name the report the run wrote (%s)', r.name));
            end
            %a naming claim is allowed only when it is true of every report
            if contains(lower(txt), 'named after')
                for r = reports(:)'
                    [~, stem] = fileparts(r.requested);
                    testCase.verifyTrue(startsWith(r.name, [stem '_error']), sprintf( ...
                        'the summary says reports are named after their outputs, but %s is not', r.name));
                end
            end
        end

        function verifyWatchedUntouchedAtEstimation(testCase)
            W = OutputNameFixture.watched();
            for w = W(:)'
                testCase.assertNotEmpty(w.atengine, 'premise: estimation must have started');
                for k = 1:numel(w.atengine)
                    testCase.verifyEqual(w.atengine{k}, w.before, sprintf( ...
                        'changed before estimation started (call %d): %s', k, w.path));
                end
            end
        end

    end
end
