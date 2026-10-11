classdef TestWorkbookReplacement < PsyRATTestBase
    %G110 (queue row 26): a GUI XLSX save onto an existing workbook must
    %replace the whole workbook.
    %
    %FINDING. Every GUI XLSX Save callback wrote its header and table with
    %writecell/writetable straight into the file the save dialog returned. On
    %an existing workbook those writers update only the cells they cover, so
    %the old workbook's trailing rows, wider columns and extra sheets survived
    %a "Replace". Reproduced through the real Save IDs callback: a planted
    %OLD_PARTICIPANT at A100 and a sheet Old both survived, read back by
    %readcell/sheetnames and independently by openpyxl (2026-09-22,
    %reconfirmed 2026-09-23).
    %
    %APPROVED CORRECTION. Write a fresh sibling workbook in the destination
    %folder, complete header and body, check that the write completed, then
    %replace the consented target. A writer failure or a failed final
    %replacement leaves the original bytes untouched and removes only the
    %task's own scratch file. Not delete-first. No cross-platform atomicity is
    %promised, and unrelated old sheets are not kept.
    %
    %ROUTES. 16 button routes reach the 13 writer sites in 10 files; the
    %census lives in G110WorkbookHarness.route and is pinned against the
    %source by testCensusCoversEveryGuiXlsxSaveSite. Each route is built from
    %the shared synthetic fixtures and fired through its button's own
    %Callback, with the save dialog answered by a double (a path, or 0 for
    %Cancel). No fit, no CmdStan.
    %
    %EXPECTED ON THE ORIGINAL CODE (e8b6e51; unchanged at 8f629e7):
    % RED, every route: testSaveReplacesWholeWorkbook (the stale row, cells
    %   and sheet survive); testEveryWriterCallTargetsAFreshSibling (lists
    %   every in-place writer call by file:line, so both writers of a
    %   two-call callback are named, not just the first callback);
    %   testWriterFailureKeepsOriginalBytes (the in-place writer has already
    %   changed the workbook when the injected failure fires);
    %   testFinalReplacementFailureKeepsOriginalBytes (see that test).
    % RED, both ID writers: testIdWritersLongToShortZeroGoodGroupAware.
    % RED, every route (added 2026-09-23 after review of the helper):
    %   testSilentHeaderCorruptionIsNotSaved and
    %   testSilentBodyCorruptionIsNotSaved; and on the 12 routes whose
    %   table holds a number, testSilentNumericCorruptionIsNotSaved. The
    %   writer returns normally but leaves a valid workbook with one required
    %   header cell, body cell or number (x + eps(x)) replaced; the original
    %   code writes that damage into the chosen workbook and says nothing.
    %   Only an exact read-back of the written workbook against the
    %   callback's own payload can refuse it; waiting for the writers to
    %   return, opening the workbook, counting its rows or comparing numbers
    %   within a tolerance cannot.
    %The fresh-workbook helper itself (psyrat_replace_workbook: payload
    %refusals, whitespace, empty tables, scratch names, backup and restore)
    %is pinned directly in TestWorkbookReplacementHelper.
    % CONTROLS, passing before and after the correction: a writer failure
    %   before any byte is written, Cancel, the CSV branch, the stale
    %   detector, the I/O seams themselves (including the corruption seam),
    %   and the census pin.
    %
    %NOT COVERED. The native macOS save dialog and Excel's view of the result
    %(native look C: Replace and Cancel, owed and not replaceable by these
    %callback runs); OS-level atomicity; the headless psyrat_report writer.
    %psyrat_write_table deletes the target first and appears here only as a
    %no-stale success control, never as a model for failure safety.

    properties (TestParameter)
        Route = struct( ...
            'depoverallt_table', 'depoverallt_table', ...
            'depoverallt_ids', 'depoverallt_ids', ...
            'trt_reloverallt_table', 'trt_reloverallt_table', ...
            'trt_reloverallt_ids', 'trt_reloverallt_ids', ...
            'variancet_table', 'variancet_table', ...
            'variancet_allsds', 'variancet_allsds', ...
            'depcutofft_table', 'depcutofft_table', ...
            'trt_relcutofft_table', 'trt_relcutofft_table', ...
            'dod_observed_table', 'dod_observed_table', ...
            'dod_varcomp_export', 'dod_varcomp_export', ...
            'dod_varcomp_window', 'dod_varcomp_window', ...
            'criterion_table', 'criterion_table', ...
            'dynrel_varcomp', 'dynrel_varcomp', ...
            'dynrel_subject', 'dynrel_subject', ...
            'splits_coef', 'splits_coef', ...
            'splits_varcomp', 'splits_varcomp')
        IdRoute = struct( ...
            'depoverallt_ids', 'depoverallt_ids', ...
            'trt_reloverallt_ids', 'trt_reloverallt_ids')
        %the routes whose callback payload holds a finite number in the table
        %body (the Save IDs sheets and the two variance tables hold text only)
        NumericRoute = struct( ...
            'depoverallt_table', 'depoverallt_table', ...
            'trt_reloverallt_table', 'trt_reloverallt_table', ...
            'depcutofft_table', 'depcutofft_table', ...
            'trt_relcutofft_table', 'trt_relcutofft_table', ...
            'dod_observed_table', 'dod_observed_table', ...
            'dod_varcomp_export', 'dod_varcomp_export', ...
            'dod_varcomp_window', 'dod_varcomp_window', ...
            'criterion_table', 'criterion_table', ...
            'dynrel_varcomp', 'dynrel_varcomp', ...
            'dynrel_subject', 'dynrel_subject', ...
            'splits_coef', 'splits_coef', ...
            'splits_varcomp', 'splits_varcomp')
    end

    methods (Test)

        function testSaveReplacesWholeWorkbook(testCase, Route)
            %Long-to-short, wide-to-narrow, a stale extra sheet and a repeated
            %save, each through the route's real Save callback onto a planted
            %workbook. The expected result is the same export saved to a
            %brand-new file, whose own completeness is first checked against
            %expectations that do not come from the writer
            %(G110WorkbookHarness.completeness).
            [env, spec, fig] = localOpen(testCase, Route);
            name = G110WorkbookHarness.TargetName;

            ref = fullfile(env.ref, 'reference.xlsx');
            r0 = G110WorkbookHarness.fireSave(fig, spec.save, ref);
            testCase.assertEmpty(r0.errText, sprintf( ...
                '%s: precondition: a save to a new file must succeed: %s', Route, r0.errText));
            testCase.assertTrue(isfile(ref), sprintf( ...
                '%s: precondition: a save to a new file must write it', Route));
            refC = G110WorkbookHarness.content(ref);
            refRep = G110WorkbookHarness.staleReport(ref);
            testCase.assertLessThan(size(refC), ...
                [G110WorkbookHarness.LongRow, G110WorkbookHarness.WideCol], sprintf( ...
                '%s: precondition: the planted stale cells must lie outside the export', Route));
            testCase.assertTrue(size(refC, 2) < 2 || isempty(refC{1, 2}), sprintf( ...
                '%s: precondition: B1 must be empty in a fresh export', Route));
            testCase.assertEqual(refRep.sheets, {'Sheet1'}, sprintf( ...
                '%s: precondition: a fresh export holds one sheet', Route));
            p = G110WorkbookHarness.completeness(spec, fig, ref);
            testCase.verifyEmpty(p, sprintf( ...
                '%s: control: a save to a new file must hold the complete header and body:\n%s', ...
                Route, localJoin(p)));

            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            planted = G110WorkbookHarness.staleReport(target);
            testCase.assertTrue(planted.long && planted.wide && planted.near && planted.sheet, ...
                sprintf('%s: precondition: the planted workbook must carry every stale marker', Route));
            before = G110WorkbookHarness.snapshot(env.work);

            for k = 1:2
                what = sprintf('%s, save %d onto an existing workbook', Route, k);
                r = G110WorkbookHarness.fireSave(fig, spec.save, target);
                testCase.verifyEmpty(r.errText, sprintf('%s: the save raised an error: %s', what, r.errText));
                rep = G110WorkbookHarness.staleReport(target);
                testCase.verifyFalse(rep.long, sprintf( ...
                    '%s: long-to-short: the old cell at %s survived', what, G110WorkbookHarness.LongCell));
                testCase.verifyFalse(rep.wide, sprintf( ...
                    '%s: wide-to-narrow: the old cell at %s survived', what, G110WorkbookHarness.WideCell));
                testCase.verifyFalse(rep.near, sprintf( ...
                    '%s: wide-to-narrow: the old cell at %s, beside the first header line, survived', ...
                    what, G110WorkbookHarness.NearCell));
                testCase.verifyFalse(rep.sheet, sprintf( ...
                    '%s: the unrelated old sheet %s survived', what, G110WorkbookHarness.OldSheet));
                testCase.verifyEqual(rep.sheets, refRep.sheets, sprintf( ...
                    '%s: the workbook must hold only the export''s sheet', what));
                testCase.verifyEqual(rep.sheetnames, rep.sheets, sprintf( ...
                    '%s: harness: sheetnames and the workbook XML disagree', what));
                C = G110WorkbookHarness.content(target);
                testCase.verifyTrue(isequal(C, refC), sprintf( ...
                    '%s: the first sheet (%dx%d) differs from a save to a new file (%dx%d)', ...
                    what, size(C, 1), size(C, 2), size(refC, 1), size(refC, 2)));
                p = G110WorkbookHarness.completeness(spec, fig, target);
                testCase.verifyEmpty(p, sprintf('%s: incomplete or stale content:\n%s', ...
                    what, localJoin(p)));
                p = G110WorkbookHarness.compareSnapshots(before, ...
                    G110WorkbookHarness.snapshot(env.work), name, false);
                testCase.verifyEmpty(p, sprintf('%s: other files in the folder changed:\n%s', ...
                    what, localJoin(p)));
            end
        end

        function testEveryWriterCallTargetsAFreshSibling(testCase, Route)
            %Every writecell/writetable call of a successful save must write a
            %scratch workbook beside the chosen file (same folder, another
            %name), and the chosen file must then be replaced by a move or
            %copy onto it. On the original code every writer call lands in the
            %chosen file itself; the message lists each one with the
            %callback's file:line.
            [env, spec, fig] = localOpen(testCase, Route);
            name = G110WorkbookHarness.TargetName;
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);

            G110WorkbookHarness.armSeams(env, name);
            r = G110WorkbookHarness.fireSave(fig, spec.save, target);
            st = G110WorkbookHarness.disarmSeams(env);

            testCase.verifyEmpty(r.errText, sprintf('%s: the save raised an error: %s', Route, r.errText));
            w = G110WorkbookHarness.classifyWrites(st, env.work, name);
            testCase.assertGreaterThan(numel(w.inplace) + numel(w.sibling) + numel(w.elsewhere), 0, ...
                sprintf(['%s: harness: no writecell/writetable call reached the seam, so the ', ...
                'checks below would be vacuous'], Route));
            testCase.verifyEmpty(w.inplace, sprintf( ...
                '%s: %d writer call(s) wrote into the chosen workbook in place:\n%s', ...
                Route, numel(w.inplace), localJoin(w.inplace)));
            testCase.verifyEmpty(w.elsewhere, sprintf( ...
                '%s: the scratch workbook must be a sibling in the destination folder:\n%s', ...
                Route, localJoin(w.elsewhere)));
            testCase.verifyGreaterThanOrEqual(st.replaceCalls, 1, sprintf( ...
                '%s: no final replacement (movefile or copyfile onto %s) was made; calls seen:\n%s', ...
                Route, name, localJoin(G110WorkbookHarness.describeWrites(st))));
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, false);
            testCase.verifyEmpty(p, sprintf('%s: other files in the folder changed:\n%s', ...
                Route, localJoin(p)));
        end

        function testWriterFailureKeepsOriginalBytes(testCase, Route)
            %The first writer call into the destination folder writes and then
            %fails. The chosen workbook must keep its exact bytes, the partial
            %scratch workbook must be gone, other files must be untouched, and
            %the failure must be visible. RED on the original code because
            %that first writer call went into the chosen workbook itself, so
            %its bytes had already changed when the failure fired. The forced
            %failure is confirmed to have happened before anything is judged.
            [env, spec, fig] = localOpen(testCase, Route);
            localForcedFailure(testCase, env, spec, fig, {'Writer', 'after', 'FailAt', 1}, 'writer');
        end

        function testFinalReplacementFailureKeepsOriginalBytes(testCase, Route)
            %The first move or copy onto the chosen workbook fails, the way
            %movefile/copyfile report a failure. The workbook must keep its
            %exact bytes, the finished scratch workbook must be removed, other
            %files must be untouched, and the failure must be visible.
            %
            %WHY THIS IS RED ON THE ORIGINAL CODE AND NOT VACUOUS. The
            %original callbacks have no replacement step: the seam is armed
            %but never reached, so a check of "original preserved when the
            %replacement fails" alone could pass for a callback that wrote
            %nothing. Two things rule that out. First, the test requires that
            %a replacement onto the chosen file was attempted (the injected
            %failure really happened); on the original code that count is 0,
            %and the message lists the writer calls made instead. Second, the
            %same save still writes -- in place -- so the planted workbook's
            %bytes change (its A1 marker is overwritten by the header) and the
            %byte check fails independently. After the correction the count
            %check guarantees the failure path was exercised, so a pass means
            %preservation, not absence of a replacement.
            [env, spec, fig] = localOpen(testCase, Route);
            localForcedFailure(testCase, env, spec, fig, {'Replace', 'failFirst'}, 'replace');
        end

        function testWriterFailureBeforeFirstWriteLeavesNoTrace(testCase, Route)
            %Control. The first writer call into the destination folder fails
            %before writing anything (as a writer that refuses a path would).
            %Passes on the original code (nothing was written) and must keep
            %passing: a helper that creates its scratch file before writing
            %must still remove it.
            [env, spec, fig] = localOpen(testCase, Route);
            localForcedFailure(testCase, env, spec, fig, {'Writer', 'before', 'FailAt', 1}, 'writer');
        end

        function testSilentHeaderCorruptionIsNotSaved(testCase, Route)
            %A silent writer defect in the HEADER. Every writer call returns
            %normally, but after the call whose payload carries the export's
            %first header line (the callback's own literal, cell A1) that
            %cell of the file just written becomes G110_CORRUPT_HEADER. The
            %workbook still opens and has every row. The save must refuse it:
            %the chosen workbook keeps its exact bytes and never holds the
            %marker, only the task's scratch file is removed, the failure is
            %visible, and a viewer prints no "Saved ..." line. RED on the
            %original code, which wrote the damaged header into the chosen
            %workbook itself and reported nothing.
            [env, spec, fig] = localOpen(testCase, Route);
            localSilentCorruption(testCase, env, spec, fig, 'corruptHeader');
        end

        function testSilentBodyCorruptionIsNotSaved(testCase, Route)
            %A silent writer defect in the BODY: as above, but the replaced
            %cell is a required body value, the first text value of the first
            %body row as the callback's own payload (or, for Save IDs, the
            %fixture) holds it; the cell becomes G110_CORRUPT_BODY. Row and
            %column counts are unchanged, so only a cell-by-cell read-back
            %of the body can refuse it.
            [env, spec, fig] = localOpen(testCase, Route);
            localSilentCorruption(testCase, env, spec, fig, 'corruptBody');
        end

        function testSilentNumericCorruptionIsNotSaved(testCase, NumericRoute)
            %A silent writer defect in a NUMBER: after the call whose payload
            %holds the first finite number of the table (the callback's own
            %value), that cell becomes the next larger double, x + eps(x).
            %Text, row and column counts are unchanged and the two values
            %agree under any tolerance, so only an exact comparison of the
            %doubles can refuse it. Same outcome checks as above. RED on the
            %original code, which wrote the changed number into the chosen
            %workbook itself and reported nothing.
            [env, spec, fig] = localOpen(testCase, NumericRoute);
            localSilentCorruption(testCase, env, spec, fig, 'corruptNumeric');
        end

        function testCancelWritesNothing(testCase, Route)
            %Control. Cancel in the save dialog writes, moves and shows
            %nothing, and the planted workbook keeps its bytes. The dialog
            %double must have been asked once, so Cancel was really exercised.
            [env, spec, fig] = localOpen(testCase, Route);
            name = G110WorkbookHarness.TargetName;
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);

            G110WorkbookHarness.armSeams(env, name);
            r = G110WorkbookHarness.fireSave(fig, spec.save, 0);
            st = G110WorkbookHarness.disarmSeams(env);

            testCase.verifyEqual(r.saveCalls, 1, sprintf( ...
                '%s: the save dialog must be asked exactly once (Cancel not exercised)', Route));
            testCase.verifyEmpty(r.errText, sprintf('%s: Cancel raised an error: %s', Route, r.errText));
            testCase.verifyEmpty(r.dialogs, sprintf('%s: Cancel opened a window: %s', ...
                Route, localJoin(r.dialogs)));
            testCase.verifyEmpty(st.writerLog, sprintf('%s: Cancel wrote:\n%s', ...
                Route, localJoin(G110WorkbookHarness.describeWrites(st))));
            testCase.verifyEmpty(st.moveLog, sprintf('%s: Cancel moved or copied a file:\n%s', ...
                Route, localJoin(G110WorkbookHarness.describeWrites(st))));
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, sprintf('%s: Cancel changed the folder:\n%s', Route, localJoin(p)));
        end

        function testCsvExportStillReplacesTextFile(testCase, Route)
            %Control. The CSV branch of the same callback already recreates
            %its file; a longer, wider planted CSV must be fully replaced by
            %text identical to a save to a new file, no workbook bytes may
            %reach a .csv name, and other files stay untouched. This guards
            %the correction against rerouting CSV saves.
            [env, spec, fig] = localOpen(testCase, Route);
            name = G110WorkbookHarness.CsvTargetName;

            ref = fullfile(env.ref, 'reference.csv');
            r0 = G110WorkbookHarness.fireSave(fig, spec.save, ref);
            testCase.assertEmpty(r0.errText, sprintf( ...
                '%s: precondition: a CSV save to a new file must succeed: %s', Route, r0.errText));
            testCase.assertTrue(isfile(ref), sprintf('%s: precondition: the CSV was not written', Route));
            refText = G110WorkbookHarness.csvText(ref);
            testCase.assertNotEmpty(refText, sprintf('%s: precondition: the CSV is empty', Route));

            target = fullfile(env.work, name);
            G110WorkbookHarness.plantCsv(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            r = G110WorkbookHarness.fireSave(fig, spec.save, target);

            testCase.verifyEmpty(r.errText, sprintf('%s: the CSV save raised an error: %s', Route, r.errText));
            txt = G110WorkbookHarness.csvText(target);
            testCase.verifyFalse(startsWith(txt, 'PK'), sprintf( ...
                '%s: a .csv save wrote workbook (zip) bytes', Route));
            testCase.verifyFalse(contains(txt, 'G110_STALE_CSV'), sprintf( ...
                '%s: lines of the planted CSV survived', Route));
            testCase.verifyTrue(strcmp(txt, refText), sprintf( ...
                '%s: the CSV differs from a CSV save to a new file', Route));
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, false);
            testCase.verifyEmpty(p, sprintf('%s: other files in the folder changed:\n%s', ...
                Route, localJoin(p)));
        end

        function testIdWritersLongToShortZeroGoodGroupAware(testCase, IdRoute)
            %The ID writers through real data, not planted cells. First a
            %two-group roster with 10 Good IDs and 4 Bad IDs, where s01 and s02
            %belong to both groups, so only the group column tells the rows
            %apart; then, onto the same file, a roster with zero Good IDs and 4
            %Bad IDs. The expected ID blocks come from the fixtures' groups,
            %not from the writer. RED on the original code: rows 5 to 10 of
            %the first roster (Good IDs and their groups) survive under the
            %second.
            env = G110WorkbookHarness.install();
            testCase.addTeardown(@() G110WorkbookHarness.uninstall(env));
            spec = G110WorkbookHarness.route(IdRoute);
            target = fullfile(env.work, 'g110_ids.xlsx');
            header = {'Good IDs', 'Bad IDs', 'Good Group', 'Bad Group'};

            long = G110WorkbookHarness.idFixture(spec, 'long');
            figLong = localOpenWith(testCase, spec, long);
            r1 = G110WorkbookHarness.fireSave(figLong, spec.save, target);
            testCase.assertEmpty(r1.errText, sprintf('%s: the first save failed: %s', IdRoute, r1.errText));
            [blk1, hdr1] = G110WorkbookHarness.idBlock(target);
            expectedLong = G110WorkbookHarness.expectedIdBlock(long);
            testCase.assertEqual(size(expectedLong), [10 4], 'harness: the long roster has 10 rows');
            testCase.verifyEqual(hdr1, header, sprintf('%s: control: the ID header row', IdRoute));
            testCase.verifyTrue(isequal(blk1, expectedLong), sprintf( ...
                '%s: control: a group-aware 10-row roster must export exactly (%dx%d written)', ...
                IdRoute, size(blk1, 1), size(blk1, 2)));

            zero = G110WorkbookHarness.idFixture(spec, 'zeroGood');
            figZero = localOpenWith(testCase, spec, zero);
            refZero = fullfile(env.ref, 'reference_zero_good.xlsx');
            r0 = G110WorkbookHarness.fireSave(figZero, spec.save, refZero);
            testCase.assertEmpty(r0.errText, sprintf('%s: the zero-Good-IDs save to a new file failed: %s', ...
                IdRoute, r0.errText));
            expectedZero = G110WorkbookHarness.expectedIdBlock(zero);
            testCase.assertEqual(size(expectedZero), [4 4], 'harness: the zero-Good-IDs roster has 4 rows');
            testCase.verifyTrue(isequal(G110WorkbookHarness.idBlock(refZero), expectedZero), sprintf( ...
                '%s: control: a zero-Good-IDs roster saved to a new file must export exactly', IdRoute));

            r2 = G110WorkbookHarness.fireSave(figZero, spec.save, target);
            testCase.verifyEmpty(r2.errText, sprintf('%s: the second save raised an error: %s', IdRoute, r2.errText));
            [blk2, hdr2] = G110WorkbookHarness.idBlock(target);
            testCase.verifyEqual(hdr2, header, sprintf('%s: the ID header row after the second save', IdRoute));
            testCase.verifyTrue(isequal(blk2, expectedZero), sprintf( ...
                ['%s: long-to-short: after a 10-row roster was replaced by a zero-Good-IDs roster ', ...
                'the ID block is %dx%d; the new roster is 4x4'], IdRoute, size(blk2, 1), size(blk2, 2)));
            if size(blk2, 2) >= 1
                testCase.verifyTrue(all(cellfun(@isempty, blk2(:, 1))), sprintf( ...
                    '%s: Good IDs from the earlier roster remain under a zero-Good-IDs save', IdRoute));
            end
            testCase.verifyTrue(isequal(G110WorkbookHarness.content(target), ...
                G110WorkbookHarness.content(refZero)), sprintf( ...
                '%s: the replaced workbook differs from a zero-Good-IDs save to a new file', IdRoute));
        end

        function testStaleDetectorControls(testCase)
            %Control for the readers the tests above rely on. Positive: an
            %in-place writecell over a planted workbook (the defect's
            %mechanism) reaches A1 yet must leave every stale marker
            %detectable by the XML scan, readcell must see the same stale row,
            %and both sheet readers must agree. Negative: psyrat_write_table,
            %which deletes the file before writing, must be reported clean
            %with the expected cells. That writer is a no-stale SUCCESS
            %control only; it destroys the original when a write then fails,
            %so it is not a model for the correction. (The XML scan is not
            %used to prove that an overwritten cell is gone: a writer may keep
            %an orphaned shared string. The tests use it only on live stale
            %cells and on freshly written workbooks.)
            env = G110WorkbookHarness.install();
            testCase.addTeardown(@() G110WorkbookHarness.uninstall(env));

            inplace = fullfile(env.work, 'inplace.xlsx');
            G110WorkbookHarness.plantWorkbook(inplace);
            writecell({'new header'; 'new line'}, inplace);
            C = G110WorkbookHarness.content(inplace);
            testCase.verifyEqual(C{1, 1}, 'new header', 'positive control: the in-place write reached A1');
            rep = G110WorkbookHarness.staleReport(inplace);
            testCase.verifyTrue(rep.long && rep.wide && rep.near && rep.sheet, ...
                'positive control: an in-place write must leave every stale marker detectable');
            testCase.verifyEqual(readcell(inplace, 'Range', ...
                [G110WorkbookHarness.LongCell ':' G110WorkbookHarness.LongCell]), {'G110_STALE_LONG'}, ...
                'positive control: readcell sees the stale row the XML scan reports');
            testCase.verifyEqual(sort(rep.sheets), sort({'Sheet1', G110WorkbookHarness.OldSheet}), ...
                'positive control: the workbook XML lists both sheets');
            testCase.verifyEqual(rep.sheetnames, rep.sheets, 'sheetnames and the workbook XML agree');

            headless = fullfile(env.work, 'headless.xlsx');
            G110WorkbookHarness.plantWorkbook(headless);
            psyrat_write_table(headless, {'Header line'; ''}, ...
                table({'r1'; 'r2'}, [1; 2], 'VariableNames', {'Label', 'Value'}));
            rep = G110WorkbookHarness.staleReport(headless);
            testCase.verifyFalse(rep.long || rep.wide || rep.near || rep.sheet || rep.plant, ...
                'negative control: a delete-first rewrite must be reported clean');
            testCase.verifyEqual(rep.sheets, {'Sheet1'}, 'negative control: one sheet after a rewrite');
            testCase.verifyEqual(G110WorkbookHarness.content(headless), ...
                {'Header line', ''; '', ''; 'Label', 'Value'; 'r1', '1'; 'r2', '2'}, ...
                'negative control: the rewritten cells read back exactly');

            testCase.verifyEqual(G110WorkbookHarness.maskTime('Generated: 03-Sep-2026 09:05:07'), ...
                'Generated: <TIME>', 'the time mask covers the header''s datestr(clock) text');
        end

        function testSeamsInterceptWritersAndReplacement(testCase)
            %Control for the I/O seams, called straight from this file (code
            %outside MATLAB's toolboxes, like the callbacks). writecell and
            %writetable are logged in order with this file as the caller and
            %forwarded (the files appear); a 'before' failure writes nothing;
            %an 'after' failure writes and then raises; the first move onto
            %the watched name fails the way movefile reports a failure and the
            %second succeeds; copyfile called without outputs raises. If
            %writetable were dispatched as a table method the log would miss
            %it, and the writer results above would need rereading.
            env = G110WorkbookHarness.install();
            testCase.addTeardown(@() G110WorkbookHarness.uninstall(env));
            name = G110WorkbookHarness.TargetName;
            a = fullfile(env.work, 'a.xlsx');
            b = fullfile(env.work, 'b.xlsx');
            c = fullfile(env.work, 'c.xlsx');
            d = fullfile(env.work, 'd.xlsx');
            target = fullfile(env.work, name);

            G110WorkbookHarness.armSeams(env, name);
            writecell({'x'}, a);
            writetable(table(1, 'VariableNames', {'v'}), b);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.verifyEqual(reshape(st.writerLog(:, 1), 1, []), {'writecell', 'writetable'}, ...
                'both writers must be logged, in call order');
            testCase.verifyTrue(all(contains(st.writerLog(:, 3), 'TestWorkbookReplacement.m')), ...
                'the log must name the calling file');
            testCase.verifyTrue(isfile(a) && isfile(b), 'logged writes must be forwarded');

            %direct calls, not verifyError(@() ...): the seams find their
            %caller with dbstack, and a direct call keeps that caller this file
            G110WorkbookHarness.armSeams(env, name, 'Writer', 'before');
            raised = '';
            try
                writecell({'x'}, c);
            catch ME
                raised = ME.identifier;
            end
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.verifyEqual(raised, 'G110Test:injectedWriterFailure', 'the before-failure raises');
            testCase.verifyTrue(st.writerInjected, 'the before-failure was injected');
            testCase.verifyFalse(isfile(c), 'a before-failure writes nothing');

            G110WorkbookHarness.armSeams(env, name, 'Writer', 'after');
            raised = '';
            try
                writetable(table(2, 'VariableNames', {'v'}), d);
            catch ME
                raised = ME.identifier;
            end
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.verifyEqual(raised, 'G110Test:injectedWriterFailure', 'the after-failure raises');
            testCase.verifyTrue(st.writerInjected, 'the after-failure was injected');
            testCase.verifyTrue(isfile(d), 'an after-failure writes, then raises');

            G110WorkbookHarness.armSeams(env, name, 'Replace', 'failFirst');
            [ok1, msg1] = movefile(a, target, 'f');
            [ok2, ~] = movefile(a, target, 'f');
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.verifyFalse(logical(ok1), 'the first move onto the watched name reports failure');
            testCase.verifyNotEmpty(msg1, 'a failed move reports a message');
            testCase.verifyTrue(logical(ok2) && isfile(target) && ~isfile(a), ...
                'the second move is forwarded');
            testCase.verifyEqual(st.replaceCalls, 2, 'both moves onto the watched name are counted');

            G110WorkbookHarness.armSeams(env, name, 'Replace', 'failFirst');
            raised = '';
            try
                copyfile(b, target, 'f');
            catch ME
                raised = ME.identifier;
            end
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.verifyEqual(raised, 'G110Test:injectedReplaceFailure', ...
                'copyfile without outputs raises when its replacement fails');
            testCase.verifyTrue(st.replaceInjected, 'copyfile onto the watched name is watched too');
        end

        function testCorruptionSeamWritesThenAltersOneCell(testCase)
            %Positive forwarding control for the silent-corruption seam,
            %called straight from this file. Each write is forwarded (the
            %file holds what was written) and raises nothing; afterwards
            %exactly the needle's cell is replaced, at the address given by
            %the call's own Range origin plus, for writetable, the
            %variable-name row; a call whose payload lacks the needle, and a
            %write outside the destination folder, are left alone.
            env = G110WorkbookHarness.install();
            testCase.addTeardown(@() G110WorkbookHarness.uninstall(env));
            name = G110WorkbookHarness.TargetName;
            a = fullfile(env.work, 'a.xlsx');
            b = fullfile(env.work, 'b.xlsx');
            c = fullfile(env.work, 'c.xlsx');
            outside = fullfile(env.ref, 'outside.xlsx');

            G110WorkbookHarness.armSeams(env, name, 'Writer', 'corruptHeader', 'Needle', 'HEAD_1');
            raised = '';
            try
                writecell({'HEAD_1'; 'HEAD_2'; ''}, a);
                writecell({'HEAD_1'}, outside);
            catch ME
                raised = ME.identifier;
            end
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.verifyEmpty(raised, 'a silent corruption raises nothing');
            testCase.verifyEqual(st.corruptApplied, 1, ...
                'only the write into the destination folder is corrupted');
            testCase.verifyEqual(reshape(st.corruptLog(:, 3), 1, []), {'A1'}, ...
                'the header needle sits at A1 of an unranged writecell');
            testCase.verifyEqual(readcell(a, 'Range', 'A1:A2'), {'G110_CORRUPT_HEADER'; 'HEAD_2'}, ...
                'the write was forwarded, then only the needle cell was replaced');
            testCase.verifyEqual(readcell(outside, 'Range', 'A1:A1'), {'HEAD_1'}, ...
                'a write outside the destination folder is left alone');

            G110WorkbookHarness.armSeams(env, name, 'Writer', 'corruptBody', 'Needle', 'r2');
            raised = '';
            try
                writecell({'H'}, b);
                writetable(table({'r1'; 'r2'}, [1; 2], 'VariableNames', {'Label', 'Value'}), b, ...
                    'Range', 'A3');
                writecell({'q'; 'r2'}, c, 'Range', 'C4');
            catch ME
                raised = ME.identifier;
            end
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.verifyEmpty(raised, 'a silent corruption raises nothing');
            testCase.verifyEqual(reshape(st.corruptLog(:, 3), 1, []), {'A5', 'C5'}, ...
                'addresses follow each call''s Range origin and the variable-name row');
            testCase.verifyEqual(G110WorkbookHarness.content(b), ...
                {'H', ''; '', ''; 'Label', 'Value'; 'r1', '1'; 'G110_CORRUPT_BODY', '2'}, ...
                'the header and table were forwarded; only the needle cell of the body changed');
            testCase.verifyEqual(readcell(c, 'Range', 'C4:C5'), {'q'; 'G110_CORRUPT_BODY'}, ...
                'a ranged writecell is corrupted at its own offset');

            %numeric: the first finite number of the table payload becomes
            %the next larger double; the text beside it is untouched
            d = fullfile(env.work, 'd.xlsx');
            G110WorkbookHarness.armSeams(env, name, 'Writer', 'corruptNumeric');
            raised = '';
            try
                writecell({'H'}, d);
                writetable(table({'r1'}, 0.1, 'VariableNames', {'Label', 'Value'}), d, 'Range', 'A3');
            catch ME
                raised = ME.identifier;
            end
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.verifyEmpty(raised, 'a silent numeric corruption raises nothing');
            testCase.verifyEqual(reshape(st.corruptLog(:, 3), 1, []), {'B4'}, ...
                'the first finite number of the table payload sits at B4');
            got = readcell(d, 'Range', 'A4:B4');
            testCase.verifyEqual(got{1}, 'r1', 'the text beside the number is untouched');
            testCase.verifyTrue(isequal(got{2}, 0.1 + eps(0.1)) && ~isequal(got{2}, 0.1), ...
                'the number was replaced by the next larger double, stored exactly');
        end

        function testCensusCoversEveryGuiXlsxSaveSite(testCase)
            %Census pin (control). Every save dialog in the toolbox that
            %offers '*.xlsx' belongs to a GUI XLSX writer, and each must be
            %reached by at least one route above: per file, the number of such
            %dialogs must equal the number of distinct writer sites the routes
            %cover. A new XLSX save site without a route fails here. Source
            %census 2026-09-23: 13 sites in 10 files, all under subroutines/.
            root = testCase.projectRoot();
            keys = G110WorkbookHarness.routeKeys();

            mc = ?TestWorkbookReplacement;
            prop = mc.PropertyList(strcmp({mc.PropertyList.Name}, 'Route'));
            testCase.assertNumElements(prop, 1, 'the Route parameter property');
            testCase.verifyEqual(sort(fieldnames(prop.DefaultValue)), sort(keys(:)), ...
                'the Route parameters and G110WorkbookHarness.routeKeys must list the same routes');

            specs = cellfun(@G110WorkbookHarness.route, keys, 'UniformOutput', false);
            specs = [specs{:}];
            files = unique({specs.file});
            sites = unique(strcat({specs.file}, '#', {specs.site}));
            testCase.verifyNumElements(sites, 13, 'the routes must reach the 13 census writer sites');
            testCase.verifyNumElements(files, 10, 'the census writer sites lie in 10 files');

            needle = '{''*.xlsx''';
            listing = [dir(fullfile(root, 'subroutines', '**', '*.m')); dir(fullfile(root, '*.m'))];
            found = containers.Map('KeyType', 'char', 'ValueType', 'double');
            for k = 1:numel(listing)
                full = fullfile(listing(k).folder, listing(k).name);
                n = numel(strfind(fileread(full), needle));
                if n > 0
                    found(full(numel(root) + 2:end)) = n;
                end
            end
            testCase.verifyEqual(sort(found.keys()), sort(files), ...
                ['every file whose save dialog offers *.xlsx must be in the census, and every ', ...
                'census file must still offer it']);
            for k = 1:numel(files)
                if found.isKey(files{k})
                    nsites = numel(unique({specs(strcmp({specs.file}, files{k})).site}));
                    testCase.verifyEqual(found(files{k}), nsites, sprintf( ...
                        '%s: %d save dialogs offer *.xlsx but the routes reach %d writer sites', ...
                        files{k}, found(files{k}), nsites));
                end
            end
        end
    end
end

% -- local helpers -------------------------------------------------------------

function [env, spec, fig] = localOpen(testCase, key)
%Install the doubles for this test (removed again by the registered teardown,
%also when a verification fails) and open the route's window.
env = G110WorkbookHarness.install();
testCase.addTeardown(@() G110WorkbookHarness.uninstall(env));
spec = G110WorkbookHarness.route(key);
fig = localOpenWith(testCase, spec, []);
end

function fig = localOpenWith(testCase, spec, pd)
%A route that does not open is a harness failure for this test only.
fig = [];
try
    fig = G110WorkbookHarness.openRoute(spec, pd);
catch ME
    testCase.assertFail(sprintf('%s: harness: the route did not open: %s', spec.key, ME.message));
end
end

function s = localJoin(c)
%One line per entry; '' for an empty list (so no message depends on how strjoin
%treats an empty cell).
if isempty(c)
    s = '';
else
    s = strjoin(reshape(c, 1, []), newline);
end
end

function localForcedFailure(testCase, env, spec, fig, seamArgs, mode)
%Shared body of the three forced-failure tests. MODE says which injection must
%be confirmed ('writer' or 'replace') before the outcome is judged; the
%outcome checks run only when the failure really happened, so a message about
%preservation always describes an exercised failure.
name = G110WorkbookHarness.TargetName;
target = fullfile(env.work, name);
G110WorkbookHarness.plantWorkbook(target);
G110WorkbookHarness.plantDecoys(env.work, name);
before = G110WorkbookHarness.snapshot(env.work);
original = G110WorkbookHarness.fileBytes(target);

G110WorkbookHarness.armSeams(env, name, seamArgs{:});
r = G110WorkbookHarness.fireSave(fig, spec.save, target);
st = G110WorkbookHarness.disarmSeams(env);

key = spec.key;
switch mode
    case 'writer'
        injected = st.writerInjected;
        testCase.verifyTrue(injected, sprintf(['%s: harness: the forced writer failure was never ', ...
            'reached (no writecell/writetable call into the destination folder):\n%s'], ...
            key, localJoin(G110WorkbookHarness.describeWrites(st))));
    case 'replace'
        injected = st.replaceInjected;
        testCase.verifyGreaterThanOrEqual(st.replaceCalls, 1, sprintf(['%s: no final replacement ', ...
            'of %s was attempted, so the forced replacement failure never happened. Calls made ', ...
            'instead:\n%s'], key, name, localJoin(G110WorkbookHarness.describeWrites(st))));
end

testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(target), original), sprintf( ...
    '%s: the chosen workbook''s bytes changed although the save did not complete', key));
p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, false);
testCase.verifyEmpty(p, sprintf('%s: the failed save left or changed other files:\n%s', ...
    key, localJoin(p)));
if injected
    testCase.verifyTrue(~isempty(r.errText) || ~isempty(r.dialogs), sprintf( ...
        '%s: the failed save must not be silent (no error and no dialog)', key));
    if strcmp(spec.kind, 'viewer')
        testCase.verifyFalse(contains(r.out, 'Saved '), sprintf( ...
            '%s: a failed save must not print its "Saved ..." line', key));
    end
end
end

function localSilentCorruption(testCase, env, spec, fig, mode)
%Shared body of the two silent-corruption tests. The needle is one required
%cell of the export, taken from the callback's literal first header line
%(header) or from the callback's payload / fixture (body, see
%G110WorkbookHarness.bodyNeedle), never from the exporter's output. A save to
%a new file without seams first locates the column-name row, so the delivered
%corruption can be checked to have hit a header cell (A1) or a body cell
%(below that row). Outcome checks follow the forced-failure tests; the
%visibility checks run only when the corruption was delivered.
key = spec.key;
name = G110WorkbookHarness.TargetName;
seamArgs = {'Writer', mode};
switch mode
    case 'corruptHeader'
        part = 'header';
        needle = spec.first;
    case 'corruptBody'
        part = 'body';
        needle = G110WorkbookHarness.bodyNeedle(spec, fig);
    otherwise
        part = 'numeric body';
        needle = '(the first finite number of the payload)';
end
testCase.assertNotEmpty(needle, sprintf('%s: harness: no text %s cell to corrupt', key, part));
if ~strcmp(mode, 'corruptNumeric')
    seamArgs = [seamArgs, {'Needle', needle}];
end

ref = fullfile(env.ref, 'reference.xlsx');
r0 = G110WorkbookHarness.fireSave(fig, spec.save, ref);
testCase.assertEmpty(r0.errText, sprintf( ...
    '%s: precondition: a save to a new file must succeed: %s', key, r0.errText));
namesRow = G110WorkbookHarness.columnNameRowOf(ref);
testCase.assertNotEmpty(namesRow, sprintf( ...
    '%s: precondition: a save to a new file has a column-name row', key));

target = fullfile(env.work, name);
G110WorkbookHarness.plantWorkbook(target);
G110WorkbookHarness.plantDecoys(env.work, name);
before = G110WorkbookHarness.snapshot(env.work);
original = G110WorkbookHarness.fileBytes(target);

G110WorkbookHarness.armSeams(env, name, seamArgs{:});
r = G110WorkbookHarness.fireSave(fig, spec.save, target);
st = G110WorkbookHarness.disarmSeams(env);

%seam delivery: the damage really landed, and in a required cell of the part
delivered = st.corruptApplied >= 1;
testCase.verifyTrue(delivered, sprintf(['%s: harness: the silent %s corruption was never ', ...
    'delivered (no writer call into the destination folder carried "%s"):\n%s'], ...
    key, part, needle, localJoin(G110WorkbookHarness.describeWrites(st))));
if delivered
    addresses = st.corruptLog(:, 3);
    rows = cellfun(@(a) str2double(regexp(a, '\d+', 'match', 'once')), addresses);
    if strcmp(part, 'header')
        testCase.verifyTrue(all(strcmp(addresses, 'A1')), sprintf( ...
            '%s: harness: the header corruption must hit A1, not %s', key, strjoin(addresses', ', ')));
    else
        testCase.verifyTrue(all(rows > namesRow), sprintf(['%s: harness: the body corruption ', ...
            'must hit a row below the column-name row %d, not %s'], key, namesRow, ...
            strjoin(addresses', ', ')));
    end
end

%outcome: nothing damaged reaches the chosen workbook, and nothing else moves
testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(target), original), sprintf( ...
    '%s: the chosen workbook''s bytes changed although the workbook written for it had a corrupted %s cell', ...
    key, part));
if ~strcmp(mode, 'corruptNumeric')
    rep = G110WorkbookHarness.staleReport(target);
    testCase.verifyFalse(rep.corrupt, sprintf( ...
        '%s: the chosen workbook now holds the silently corrupted %s cell', key, part));
end
p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, false);
testCase.verifyEmpty(p, sprintf('%s: the refused save left or changed other files:\n%s', ...
    key, localJoin(p)));
if delivered
    testCase.verifyTrue(~isempty(r.errText) || ~isempty(r.dialogs), sprintf( ...
        '%s: a save whose written %s failed its read-back must not be silent (no error and no dialog)', ...
        key, part));
    if strcmp(spec.kind, 'viewer')
        testCase.verifyFalse(contains(r.out, 'Saved '), sprintf( ...
            '%s: a refused save must not print its "Saved ..." line', key));
    end
end
end
