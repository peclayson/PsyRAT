classdef TestWorkbookReplacementHelper < PsyRATTestBase
    %psyrat_replace_workbook (G110), the fresh-workbook helper the GUI XLSX
    %Save callbacks call, pinned directly on small payloads. The callback
    %routes themselves are covered in TestWorkbookReplacement.
    %
    %EXPECTATIONS come from the payloads and from what the XLSX writer is
    %known to store (round trip of 2026-09-23, reopened independently):
    %finite doubles exactly, denormals included; logicals written by
    %writecell; NaN as an EMPTY cell; Inf and -Inf as 65535, so they must be
    %refused; text over 32,767 characters cut, so it must be refused; leading
    %and trailing spaces stored exactly. readcell trims those spaces when it
    %reads, so the whitespace checks read the workbook's own XML; every
    %number is compared with readcell exactly, never within a tolerance.
    %
    %The helper's documented temporary names (psyrat_stage_<k>.xlsx and
    %psyrat_backup_<k>.xlsx, first free k without regard to case) are what
    %the decoys imitate: an empty file, a folder and a case variant of each.
    %
    %REPLACEMENT FAULTS use the G110WorkbookHarness movefile/copyfile seam:
    %the first move onto the chosen file fails leaving it intact, or damages
    %it or removes it and then fails, or (damageAll) every later copy onto it
    %fails too. The chosen file must end byte-identical to the original,
    %except in the last case, where the checked copy must be kept and the
    %message must say where, without claiming the original is unchanged.
    %
    %PHASE 3 (RED on the phase-2 helper, controls pass on it). The
    %G110WorkbookHarness fault doubles (fread, ferror, fwrite, delete) act only
    %on calls made directly from psyrat_replace_workbook.m on files in the test
    %folder; MATLAB's own I/O, writecell's included, passes straight through.
    %REDs: a temporary file whose removal is refused must be reported with its
    %full path and never called removed (after a writer failure, after a failed
    %checked copy, after a success); an original read only in part must stop
    %the save before anything is moved; a failed write into the checked copy
    %must not leave its file id open; a carriage return, which XML readers turn
    %into a line feed, must be refused with the CSV remedy. Controls: the
    %message capture reads the current message of a REUSED modal dialog; the
    %fault doubles act only where aimed and change nothing when unaimed; tab,
    %line feed and spaces keep saving.
    %
    %PHASE 3b. Open file ids are listed by G110WorkbookHarness.openFileIds
    %(openedFiles, which replaces fopen('all')). The fault-double control
    %compares every unaimed ferror with the built-in's own answer for the same
    %whole read, and covers a whole re-read of the aimed file on a reused id.
    %RED: a failed save through the DoD setup screen's Export button must show
    %ONE error dialog; the helper shows one and the button's catch shows a
    %second. errordlg CALLS are counted (G110WorkbookHarness errordlg
    %recorder), not windows, because a modal errordlg reuses an open window of
    %its title. Control: the same DoD writer saved from its own window (no
    %wrapping catch) shows one, and a successful DoD export shows none.
    %
    %PHASE 4 (RED on the phase-2 helper unless marked control). Reader facts
    %from the probe of 2026-09-23: xlsread and readcell return a sheet's used
    %range, dropping leading empty rows and columns, so a workbook written one
    %row lower or one column to the right reads back as the same block; an
    %explicit A1-anchored read tells them apart. REDs: every write moved down,
    %right, or both (the harness writer seam's shift mode) must fail the
    %check; text holding "_X0041_" with an upper-case X, which is not the XLSX
    %escape form and is stored and read back exactly, must save; text in the
    %escape form with a lower-case x stays refused before writing (the writer
    %stores it correctly, but xlsread 'basic' returns it escaped), and the
    %message must place that limit in this save's check, not in XLSX. Control:
    %the same seam with a zero shift saves and reads back from A1.
    %
    %PHASE 7 (pins of clauses the helper already enforces: they pass on it and
    %fail on a copy without the clause). The harness writer seam's extra mode
    %adds a second sheet, or one text cell right of or below the payload,
    %after each forwarded write; only the one-sheet check or the used-range
    %read can refuse that, each in its own words. Its control adds nothing,
    %and the payload's own empty cells (a blank beside the header, a NaN at
    %the end of the last row) must pass. The move seam's
    %wrongButReportsSuccess mode leaves other bytes under the chosen name and
    %reports success: with an original, the helper must put it back and say
    %so; with none, it must say the file now there is not to be relied on. The
    %shift pin now requires the refusal at A1, and the warning-as-error pin
    %also reads the copy's path from the console.
    %
    %No fit, no CmdStan, no desktop.

    properties (TestParameter)
        %how far every write of the new workbook is moved: [rows columns]
        Shift = struct('downOneRow', [1 0], 'rightOneColumn', [0 1], 'downAndRight', [1 1]);
        %what the extra-content pin adds after each write of the new workbook:
        %a second sheet, or one text cell just outside the payload (A1:B4)
        Extra = struct('secondSheet', 'sheet', 'cellRightOfPayload', 'C1', 'cellBelowPayload', 'A5');
    end

    methods (Test)

        function testSavesNumbersLogicalsAndEmptyCellsExactly(testCase)
            env = localEnv(testCase);
            f = fullfile(env.work, 'numbers.xlsx');
            head = {'Header line 1'; 'Header line 2'; ''};
            T = table({'a'; 'b'; 'c'}, [pi; 1/3; NaN], [0.1; -0; 5e-324], ...
                [2^53; realmax; 42], {'x'; ''; 3}, ...
                'VariableNames', {'Label', 'Irrational', 'Small', 'Large', 'Mixed'});
            r = localCall(f, head, T, 3);
            testCase.assertEmpty(r.errText, ['a representable table was refused: ' r.errText]);
            raw = readcell(f);
            testCase.verifyEqual(size(raw), [6 5], 'header rows, names row and three body rows');
            testCase.verifyEqual(raw(1:2, 1)', {'Header line 1', 'Header line 2'});
            testCase.verifyEqual(raw(3, :), T.Properties.VariableNames, 'names on the given row');
            testCase.verifyEqual(raw(4:6, 1)', {'a', 'b', 'c'});
            testCase.verifyTrue(isequal(raw{4, 2}, pi) && isequal(raw{5, 2}, 1/3), ...
                'irrational doubles are kept exactly');
            testCase.verifyTrue(ismissing(raw{6, 2}), 'NaN is stored as an empty cell');
            testCase.verifyTrue(isequal(raw{4, 3}, 0.1) && raw{5, 3} == 0 && isequal(raw{6, 3}, 5e-324), ...
                'small doubles, -0 and the smallest denormal are kept exactly');
            testCase.verifyTrue(isequal(raw{4, 4}, 2^53) && isequal(raw{5, 4}, realmax) && ...
                isequal(raw{6, 4}, 42), 'large doubles are kept exactly');
            testCase.verifyTrue(strcmp(raw{4, 5}, 'x') && ismissing(raw{5, 5}) && isequal(raw{6, 5}, 3), ...
                'a mixed cell column keeps its text, its empty cell and its number');
            testCase.verifyNumElements(sheetnames(f), 1, 'one sheet');

            g = fullfile(env.work, 'logicals.xlsx');
            r = localCall(g, {'Header'; true; false; []});
            testCase.assertEmpty(r.errText, ['a header with logicals was refused: ' r.errText]);
            raw = readcell(g);
            testCase.verifyEqual(raw(1:3)', {'Header', true, false}, 'logicals stay logicals');
            testCase.verifyEqual(localEntries(env.work), sort({'logicals.xlsx', 'numbers.xlsx'}), ...
                'no temporary file is left');
        end

        function testKeepsLeadingAndTrailingSpacesExactly(testCase)
            %The regression that must not happen: every overall table's
            %coefficient text begins with a space (' 0.85 CI [...]'). A check
            %that trimmed would refuse every normal export; one that ignored
            %spaces could not see a real change.
            env = localEnv(testCase);
            f = fullfile(env.work, 'spaces.xlsx');
            T = table({'  pad  '}, {' 0.85 CI [0.76 0.92]'}, ...
                'VariableNames', {'Label', 'Dependability'});
            r = localCall(f, {' Header with a leading space'; ''}, T, 2);
            testCase.assertEmpty(r.errText, ['padded text was refused: ' r.errText]);
            xml = G110WorkbookHarness.xmlText(f);
            testCase.verifyTrue(contains(xml, '>  pad  <'), 'two spaces on each side are stored');
            testCase.verifyTrue(contains(xml, '> 0.85 CI [0.76 0.92]<'), 'the leading space is stored');
            testCase.verifyTrue(contains(xml, '> Header with a leading space<'), ...
                'a header line keeps its leading space');
        end

        function testKeepsUnicodeAndLineBreaksExactly(testCase)
            %Both are stored exactly by the writer (round trip, 2026-09-23).
            %If this fails with a check message, the helper's reader does not
            %read them back exactly; that is a limitation to adjudicate, not a
            %reason to loosen the comparison.
            env = localEnv(testCase);
            f = fullfile(env.work, 'text.xlsx');
            symbol = char([955 32 8805 32 945]);
            note = sprintf('line one\nline two');
            T = table({symbol}, {note}, 'VariableNames', {'Symbol', 'Note'});
            r = localCall(f, {'Header'; ''}, T, 2);
            testCase.assertEmpty(r.errText, ['representable text was refused: ' r.errText]);
            raw = readcell(f);
            testCase.verifyEqual(raw{3, 1}, symbol, 'non-ASCII text is kept');
            testCase.verifyEqual(raw{3, 2}, note, 'a line break is kept');
        end

        function testKeepsLongTextUpToTheXlsxLimit(testCase)
            %Long provenance lines and the data-splits design metadata (JSON)
            %can pass 255 characters. Text up to 32,767 characters is stored
            %whole (the writer cut only the 32,768th character in the round
            %trip), so it must save; a reader that cut it would refuse normal
            %exports, which is a limitation to adjudicate, not to hide.
            env = localEnv(testCase);
            f = fullfile(env.work, 'long.xlsx');
            line300 = ['Header ' repmat('h', 1, 293)];
            json = ['{"design":"' repmat('j', 1, 1000) '"}'];
            longest = repmat('x', 1, 32767);
            T = table({json; longest}, 'VariableNames', {'Text'});
            r = localCall(f, {line300; ''}, T, 2);
            testCase.assertEmpty(r.errText, ['representable long text was refused: ' r.errText]);
            raw = readcell(f);
            testCase.verifyEqual(raw{1, 1}, line300, 'a 300-character header line is kept whole');
            testCase.verifyEqual(raw{3, 1}, json, 'a 1,000-character JSON cell is kept whole');
            testCase.verifyEqual(numel(raw{4, 1}), 32767, 'a 32,767-character cell is kept whole');
        end

        function testEmptyTableKeepsHeaderAndPlacesNames(testCase)
            %writetable(emptyTable, ..., 'Range', 'A3') writes the names at A1,
            %over the header (observed 2026-09-23). A typed table with no rows
            %must keep the header and put its names on the given row.
            env = localEnv(testCase);
            f = fullfile(env.work, 'empty.xlsx');
            T0 = table(cell(0, 1), zeros(0, 1), 'VariableNames', {'Label', 'Value'});
            r = localCall(f, {'Header 1'; 'Header 2'; ''}, T0, 3);
            testCase.assertEmpty(r.errText, ['an empty table was refused: ' r.errText]);
            raw = readcell(f);
            testCase.verifyEqual(size(raw), [3 2], 'two header rows and the names row, nothing else');
            testCase.verifyEqual(raw{1, 1}, 'Header 1', 'the first header line is not overwritten');
            testCase.verifyEqual(raw{2, 1}, 'Header 2');
            testCase.verifyEqual(raw(3, :), {'Label', 'Value'}, 'the names are on the given row');
        end

        function testRefusesValuesThisSaveCannotWriteAndCheckExactly(testCase)
            %Refused before anything is written: the chosen workbook keeps its
            %bytes, the folder is unchanged, no writer is called, and the
            %refusal is shown. Some of these the XLSX writer would change (Inf,
            %text over 32,767 characters) or this export does not support (a
            %datetime); the others the writer stores, but the check's reader
            %cannot return them exactly (a control character, "_xHHHH_" with a
            %lower-case x), a limit of this save's writer and reader, not of XLSX.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            cases = {
                'Inf in a table',        {'Header'; ''},                 table({'a'}, Inf, 'VariableNames', {'L', 'V'})
                '-Inf in a cell',        {'Header'; -Inf},               []
                '32,768 characters',     {'Header'; repmat('x', 1, 32768)}, []
                'a control character',   {'Header'; ['a' char(1) 'b']},  []
                'the _xHHHH_ sequence',  {'Header'; 'a_x0041_b'},        []
                'a datetime',            {'Header'; datetime(2026, 9, 23)}, []};
            for k = 1:size(cases, 1)
                what = cases{k, 1};
                G110WorkbookHarness.armSeams(env, name);
                if isempty(cases{k, 3})
                    r = localCall(target, cases{k, 2});
                else
                    r = localCall(target, cases{k, 2}, cases{k, 3}, 2);
                end
                st = G110WorkbookHarness.disarmSeams(env);
                testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:refused', ...
                    sprintf('%s must be refused before writing: %s', what, r.errText));
                testCase.verifyEmpty(st.writerLog, sprintf('%s: nothing may be written', what));
                testCase.verifyEmpty(st.moveLog, sprintf('%s: nothing may be moved', what));
                testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(target), original), ...
                    sprintf('%s: the chosen workbook keeps its bytes', what));
                p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), ...
                    name, true);
                testCase.verifyEmpty(p, sprintf('%s: the folder must be unchanged:\n%s', what, localJoin(p)));
                testCase.verifyNotEmpty(r.dialogText, sprintf('%s: the refusal must be shown', what));
            end
        end

        function testRefusesAFolderTarget(testCase)
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            mkdir(target);
            before = G110WorkbookHarness.snapshot(env.work);
            G110WorkbookHarness.armSeams(env, name);
            r = localCall(target, {'Header'});
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:refused', r.errText);
            testCase.verifyTrue(isfolder(target), 'the folder is left in place');
            testCase.verifyEmpty(st.writerLog, 'nothing may be written');
            testCase.verifyEmpty(st.moveLog, 'nothing may be moved');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, ['the folder must be unchanged: ' localJoin(p)]);
        end

        function testLongTargetNameUsesShortScratchNames(testCase)
            %A scratch name built from the chosen name would pass the
            %file-system limit for a chosen name near it. The longest name
            %the platform allows here: 240 characters, or on Windows what fits
            %the 260-character path limit.
            env = localEnv(testCase);
            len = 240;
            if ispc
                len = min(len, 259 - numel(env.work) - 1);
            end
            name = [repmat('n', 1, len - 5) '.xlsx'];
            target = fullfile(env.work, name);
            r = localCall(target, {'First'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            testCase.assertEmpty(r.errText, ['a long valid name was refused: ' r.errText]);

            G110WorkbookHarness.armSeams(env, name);
            r = localCall(target, {'Second'; ''}, table({'b'}, 2, 'VariableNames', {'L', 'V'}), 2);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.assertEmpty(r.errText, ['replacing a long valid name failed: ' r.errText]);
            w = G110WorkbookHarness.classifyWrites(st, env.work, name);
            testCase.verifyNotEmpty(w.sibling, 'the new workbook is written beside the chosen file');
            testCase.verifyEmpty([w.inplace, w.elsewhere], 'and nowhere else');
            for k = 1:size(st.writerLog, 1)
                [~, n, e] = fileparts(st.writerLog{k, 2});
                testCase.verifyLessThanOrEqual(numel([n e]), 30, ...
                    sprintf('the scratch name %s must not grow with the chosen name', [n e]));
            end
            testCase.verifyEqual(readcell(target, 'Range', 'A1:A1'), {'Second'}, 'the file was replaced');
            testCase.verifyEqual(localEntries(env.work), {name}, 'no temporary file is left');
        end

        function testScratchAndBackupNamesSkipDecoysWithoutRegardToCase(testCase)
            %Decoys hold the helper's own temporary names: an empty file, a
            %folder and an upper-case variant of each. None may be used,
            %written, moved or removed, and no temporary file may be left.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            decoys = lower(regexprep(setdiff({before.name}, {name}), '[\\/]$', ''));
            G110WorkbookHarness.armSeams(env, name);
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.assertEmpty(r.errText, ['the save failed: ' r.errText]);
            for k = 1:size(st.writerLog, 1)
                [~, n, e] = fileparts(st.writerLog{k, 2});
                testCase.verifyFalse(ismember(lower([n e]), decoys), ...
                    sprintf('%s is a decoy''s name', [n e]));
                testCase.verifyTrue(startsWith([n e], 'psyrat_stage_'), ...
                    sprintf('%s is not the documented scratch name', [n e]));
            end
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, false);
            testCase.verifyEmpty(p, ['decoys must be untouched and nothing left: ' localJoin(p)]);
        end

        function testSilentNumericCorruptionIsRefused(testCase)
            %The written number becomes the next larger double; the check must
            %compare exactly and refuse, removing only its own scratch file.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            G110WorkbookHarness.armSeams(env, name, 'Writer', 'corruptNumeric');
            r = localCall(target, {'Header'; ''}, table({'a'}, pi, 'VariableNames', {'L', 'V'}), 2);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.assertGreaterThanOrEqual(st.corruptApplied, 1, 'harness: the corruption was not delivered');
            testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:check', r.errText);
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, ['the chosen workbook and the folder must be unchanged: ' localJoin(p)]);
            testCase.verifyNotEmpty(r.dialogText, 'the refusal must be shown');
        end

        function testSuccessfulReplaceLeavesNoBackup(testCase)
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            testCase.assertEmpty(r.errText, ['the save failed: ' r.errText]);
            rep = G110WorkbookHarness.staleReport(target);
            testCase.verifyFalse(rep.long || rep.wide || rep.near || rep.sheet, 'the old workbook is gone');
            testCase.verifyEqual(readcell(target, 'Range', 'A1:A1'), {'Header'});
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, false);
            testCase.verifyEmpty(p, ['no backup or scratch file may be left: ' localJoin(p)]);
        end

        function testNewTargetMakesNoBackup(testCase)
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            r = localCall(fullfile(env.work, name), {'Header'});
            testCase.assertEmpty(r.errText, ['the save failed: ' r.errText]);
            after = G110WorkbookHarness.snapshot(env.work);
            testCase.verifyEqual(setdiff({after.name}, {before.name}), {name}, ...
                'only the chosen file is new: no backup, no scratch file');
            testCase.verifyEmpty(localChanged(before, after), 'every decoy is untouched');
        end

        function testFailedReplaceLeavesTheOriginalByteIdentical(testCase)
            %failFirst: the move fails and leaves the file alone. damageFirst /
            %removeFirst: the move changes or removes the file, then fails;
            %the original must be put back from the checked copy. Either way
            %the chosen file ends byte-identical, no temporary file is left,
            %and the failure is shown.
            for mode = {'failFirst', 'damageFirst', 'removeFirst'}
                env = localEnv(testCase);
                name = 'chosen.xlsx';
                target = fullfile(env.work, name);
                G110WorkbookHarness.plantWorkbook(target);
                G110WorkbookHarness.plantDecoys(env.work, name);
                before = G110WorkbookHarness.snapshot(env.work);
                G110WorkbookHarness.armSeams(env, name, 'Replace', mode{1});
                r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
                st = G110WorkbookHarness.disarmSeams(env);
                testCase.assertTrue(st.replaceInjected, sprintf('%s: harness: the fault was not delivered', mode{1}));
                testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:replace', ...
                    sprintf('%s: %s', mode{1}, r.errText));
                p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), ...
                    name, true);
                testCase.verifyEmpty(p, sprintf(['%s: the chosen file must be byte-identical to the ', ...
                    'original and nothing left:\n%s'], mode{1}, localJoin(p)));
                testCase.verifyNotEmpty(r.dialogText, sprintf('%s: the failure must be shown', mode{1}));
            end
        end

        function testRestoreFailureKeepsTheBackupAndSaysWhere(testCase)
            %damageAll: the move damages the chosen file and fails, and putting
            %the original back fails too. The checked copy must be kept,
            %byte-identical to the original; the message must give its full
            %path and must not claim the original is unchanged; the scratch
            %file must be gone and the decoys untouched.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            G110WorkbookHarness.armSeams(env, name, 'Replace', 'damageAll');
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.assertTrue(st.replaceInjected, 'harness: the fault was not delivered');
            testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:restore', r.errText);
            after = G110WorkbookHarness.snapshot(env.work);
            added = setdiff({after.name}, {before.name});
            testCase.assertNumElements(added, 1, ['exactly one new entry, the kept copy: ' localJoin(added)]);
            kept = fullfile(env.work, added{1});
            testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(kept), original), ...
                'the kept copy is byte-identical to the original');
            testCase.verifyTrue(contains(r.errText, kept), 'the message gives the kept copy''s full path');
            said = lower([r.errText ' ' r.dialogText]);
            testCase.verifyFalse(contains(said, 'unchanged') || contains(said, 'same as before'), ...
                'the message must not claim the original is unchanged');
            testCase.verifyEmpty(localChanged(before, after, name), 'every decoy is untouched');
        end

        % -- phase 3: capture control, leftover honesty, read/write faults, CR --

        function testDialogCaptureReadsTheCurrentMessage(testCase)
            %Control for the capture itself. A modal errordlg REUSES an open
            %dialog of the same title (msgbox.m), so no new window appears;
            %counting new windows lost every message after the first. Here a
            %dialog titled like the helper's is left open first, then two
            %refusals with different reasons follow. Each capture must hold its
            %own reason and neither the other reason nor the placeholder, and
            %the placeholder window must still be open, so the reuse path, the
            %one the old capture missed, is the one exercised.
            env = localEnv(testCase);
            target = fullfile(env.work, 'chosen.xlsx');
            placeholder = errordlg({'G110 placeholder: an earlier message'}, 'File Not Saved', 'modal');
            testCase.addTeardown(@() delete(placeholder(isvalid(placeholder))));
            rA = localCall(target, {'Header'; Inf});
            rB = localCall(target, {'Header'; repmat('x', 1, 32768)});
            testCase.verifyEqual({rA.errId, rB.errId}, ...
                {'psyrat:replaceWorkbook:refused', 'psyrat:replaceWorkbook:refused'}, 'both are refusals');
            testCase.verifyTrue(isvalid(placeholder), 'the open dialog was reused, not closed');
            testCase.verifyTrue(contains(rA.dialogText, 'infinite') && ...
                ~contains(rA.dialogText, '32,767') && ~contains(rA.dialogText, 'placeholder'), ...
                ['the first capture must hold the first reason only: ' rA.dialogText]);
            testCase.verifyTrue(contains(rB.dialogText, '32,767') && ...
                ~contains(rB.dialogText, 'infinite') && ~contains(rB.dialogText, 'placeholder'), ...
                ['the second capture must hold the second reason only: ' rB.dialogText]);
        end

        function testScratchLeftAfterWriterFailureIsReportedNotClaimedRemoved(testCase)
            %The writer writes the scratch workbook and fails; deleting it is
            %then refused the way MATLAB refuses it (a warning, the file stays).
            %Required: the chosen file keeps its bytes, the decoys are
            %untouched, the scratch file is the only new entry, and the error
            %gives its full path and does not say it was removed. RED on the
            %phase-2 helper: local_remove only warns, and the message still says
            %"The new workbook was removed".
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            G110WorkbookHarness.armSeams(env, name, 'Writer', 'after', 'FailAt', 1);
            G110WorkbookHarness.armFaults(env, 'RefuseDelete', {'psyrat_stage_'});
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            fst = G110WorkbookHarness.disarmFaults(env);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.assertTrue(st.writerInjected, 'harness: the writer failure was not delivered');
            testCase.assertNotEmpty(fst.deleteRefused, 'harness: no delete of the scratch file was refused');
            localVerifyLeftovers(testCase, env, before, name, original, r, fst.deleteRefused, ...
                'psyrat:replaceWorkbook:check');
        end

        function testLeftoversAfterBackupFailureAreReportedNotClaimedRemoved(testCase)
            %The checked copy of the chosen file fails part-way (the first
            %write into psyrat_backup_* writes half, then raises), and deleting
            %both temporary files is refused. Required as above, with both files
            %named. RED on the phase-2 helper: the message says "The new
            %workbook was removed" and names neither file.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            fids = G110WorkbookHarness.openFileIds();
            G110WorkbookHarness.armFaults(env, 'FailWrite', 'psyrat_backup_', ...
                'RefuseDelete', {'psyrat_stage_', 'psyrat_backup_'});
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            fst = G110WorkbookHarness.disarmFaults(env);
            %a leaked id is the next test's subject
            localCloseFiles(setdiff(G110WorkbookHarness.openFileIds(), fids));
            testCase.assertTrue(fst.writeFailed, 'harness: the backup write failure was not delivered');
            testCase.assertNotEmpty(fst.deleteRefused, 'harness: no delete was refused');
            localVerifyLeftovers(testCase, env, before, name, original, r, fst.deleteRefused, ...
                'psyrat:replaceWorkbook:backup');
        end

        function testBackupLeftAfterSuccessIsReported(testCase)
            %The save succeeds, but removing the checked copy afterwards is
            %refused. Required: no error (the new workbook is in place), the
            %copy is the only new entry and holds the original bytes, and its
            %full path is reported in a dialog or warning, not claimed removed.
            %RED on the phase-2 helper: local_remove passes on nothing but
            %MATLAB's own warning ("File not found or permission denied"), which
            %does not say which file, so the copy is left where nobody looks.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            G110WorkbookHarness.armFaults(env, 'RefuseDelete', {'psyrat_backup_'});
            lastwarn('');
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            [warned, ~] = lastwarn;
            fst = G110WorkbookHarness.disarmFaults(env);
            testCase.assertNotEmpty(fst.deleteRefused, 'harness: the delete of the copy was not refused');
            testCase.verifyEmpty(r.errText, ['the save itself succeeded: ' r.errText]);
            testCase.verifyEqual(readcell(target, 'Range', 'A1:A1'), {'Header'}, 'the new workbook is in place');
            after = G110WorkbookHarness.snapshot(env.work);
            added = setdiff({after.name}, {before.name});
            testCase.assertNumElements(added, 1, ['exactly one new entry, the copy: ' localJoin(added)]);
            kept = fullfile(env.work, added{1});
            testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(kept), original), ...
                'the copy left behind holds the original bytes');
            reported = [r.dialogText ' ' warned];
            testCase.verifyTrue(contains(reported, kept), ...
                ['the copy left behind must be reported with its full path; reported: ' reported]);
            testCase.verifyFalse(contains(reported, 'was removed'), 'a file still there is not "removed"');
            testCase.verifyEmpty(localChanged(before, after, name), 'every decoy is untouched');
        end

        function testUnreadableOriginalRefusesBeforeReplacing(testCase)
            %The helper's first read of the chosen file returns half its bytes
            %and ferror then reports an error (a read that failed part-way).
            %Required: the save stops before anything is moved or copied onto
            %the chosen file, the chosen file is byte-identical, no temporary
            %file is left, and the refusal is shown. RED on the phase-2 helper:
            %local_bytes asks neither ferror nor the file size, so the half file
            %becomes the "checked copy" and the file is replaced; had the move
            %failed, the half file is what would have been put back.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            G110WorkbookHarness.armSeams(env, name);   %counts moves and copies onto the chosen file
            G110WorkbookHarness.armFaults(env, 'CutRead', target);
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            fst = G110WorkbookHarness.disarmFaults(env);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.assertTrue(fst.readCut, 'harness: the partial read was not delivered');
            testCase.verifyNotEmpty(r.errId, 'an original that was not read in full must stop the save');
            testCase.verifyEqual(st.replaceCalls, 0, 'nothing may be moved or copied onto the chosen file');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, ['the chosen file and the folder must be unchanged: ' localJoin(p)]);
            testCase.verifyNotEmpty(r.dialogText, 'the refusal must be shown');
        end

        function testFailedBackupWriteLeavesNoOpenFile(testCase)
            %The first write into the checked copy writes half and raises.
            %Required: no file id is left open, the chosen file is
            %byte-identical, no temporary file is left, and the refusal is
            %shown. RED on the phase-2 helper: local_write_bytes has no cleanup,
            %so the copy's file id stays open.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            fids = G110WorkbookHarness.openFileIds();
            G110WorkbookHarness.armFaults(env, 'FailWrite', 'psyrat_backup_');
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            fst = G110WorkbookHarness.disarmFaults(env);
            leaked = setdiff(G110WorkbookHarness.openFileIds(), fids);
            stillOpen = arrayfun(@(f) fopen(f), leaked, 'UniformOutput', false);
            localCloseFiles(leaked);
            testCase.assertTrue(fst.writeFailed, 'harness: the backup write failure was not delivered');
            testCase.verifyEmpty(stillOpen, ['file ids left open: ' localJoin(stillOpen)]);
            testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:backup', r.errText);
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, ['the chosen file and the folder must be unchanged: ' localJoin(p)]);
            testCase.verifyNotEmpty(r.dialogText, 'the refusal must be shown');
        end

        function testFaultSeamsActOnlyWhereAimed(testCase)
            %Control for the fault doubles, aimed at THIS file. A read by
            %another file (the harness) is whole; the first read here of the
            %aimed file returns half and ferror reports the injected error; a
            %second, whole read of the aimed file (normally on the same, reused
            %file id) reports what the built-in reports, not the earlier
            %injected error; another file reads whole and ferror gives the
            %built-in's answer. The first write into an aimed name writes half
            %and raises. Deleting an aimed name leaves the file with MATLAB's
            %warning; another delete is done. "What the built-in reports" is
            %measured, not assumed: the same whole read of the same file is
            %made before the doubles are armed.
            env = localEnv(testCase);
            a = fullfile(env.work, 'g110_read_me.bin');
            other = fullfile(env.work, 'g110_other.bin');
            localWriteBytes(a, uint8(1:100));
            localWriteBytes(other, uint8(1:100));
            keep = fullfile(env.work, 'g110_keep_1.txt');
            gone = fullfile(env.work, 'g110_gone_1.txt');
            localWriteBytes(keep, uint8('k'));
            localWriteBytes(gone, uint8('g'));
            w = fullfile(env.work, 'g110_write_1.bin');
            %the built-in's own answer for a whole read, doubles not yet armed
            fid = fopen(other, 'r');
            bRef = fread(fid, inf, '*uint8');
            eRef = ferror(fid);
            fclose(fid);
            testCase.log(matlab.unittest.Verbosity.Terse, ...
                sprintf('built-in ferror after a whole read: "%s"', eRef));
            G110WorkbookHarness.armFaults(env, 'Caller', 'TestWorkbookReplacementHelper.m', ...
                'CutRead', a, 'FailWrite', 'g110_write_', 'RefuseDelete', {'g110_keep_'});
            viaHarness = G110WorkbookHarness.fileBytes(a);
            fid = fopen(a, 'r');
            fidCut = fid;
            b1 = fread(fid, inf, '*uint8');
            e1 = ferror(fid);
            fclose(fid);
            fid = fopen(a, 'r');
            fidWhole = fid;
            b2 = fread(fid, inf, '*uint8');
            e2 = ferror(fid);
            fclose(fid);
            fid = fopen(other, 'r');
            b3 = fread(fid, inf, '*uint8');
            e3 = ferror(fid);
            fclose(fid);
            fid = fopen(w, 'w');
            raised = '';
            try
                fwrite(fid, uint8(1:10), 'uint8');
            catch ME
                raised = ME.identifier;
            end
            fclose(fid);
            lastwarn('');
            delete(keep);
            [~, warnedId] = lastwarn;
            delete(gone);
            fst = G110WorkbookHarness.disarmFaults(env);
            testCase.verifyEqual(bRef, uint8(1:100)', 'harness: the built-in reads the whole reference file');
            testCase.verifyEqual(numel(viaHarness), 100, 'a read from another file is not cut');
            testCase.verifyEqual(numel(b1), 50, 'the first aimed read returns half the bytes');
            testCase.verifySubstring(e1, 'Injected read error', 'ferror reports the cut read');
            testCase.verifyEqual(b2, uint8(1:100)', 'only the first aimed read is cut');
            testCase.verifyEqual(e2, eRef, sprintf(['a whole re-read of the aimed file (file id %d; ', ...
                'the cut read used %d) must give the built-in''s ferror, not the earlier injected ', ...
                'error'], fidWhole, fidCut));
            testCase.verifyEqual(b3, uint8(1:100)', 'another file reads whole');
            testCase.verifyEqual(e3, eRef, 'ferror on another file is the built-in''s answer');
            testCase.verifyEqual(raised, 'G110Test:injectedWriteFailure', 'the aimed write raises');
            testCase.verifyEqual(G110WorkbookHarness.fileBytes(w), uint8(1:5)', 'after writing half');
            testCase.verifyTrue(isfile(keep), 'the refused delete leaves the file');
            testCase.verifyEqual(warnedId, 'MATLAB:DELETE:Permission', 'with MATLAB''s own warning');
            testCase.verifyFalse(isfile(gone), 'an unaimed delete is done');
            testCase.verifyTrue(fst.readCut && fst.writeFailed && numel(fst.deleteRefused) == 1, ...
                'the doubles record what they did');
        end

        function testFaultSeamsArmedButUnaimedChangeNothing(testCase)
            %Control: every fault double armed at names and a file the save
            %never touches. The helper's own reads, writes and deletes all
            %pass through the doubles to the built-ins, and the save is
            %normal: no error, the file replaced, nothing left, no id open.
            %Of that pass-through, only the ferror calls of the helper's reads
            %are checked: the gate counts them (two per read), and at least one
            %must arrive. Zero shows only that none reached the ferror double;
            %the fread, fwrite and delete calls are not counted.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            fids = G110WorkbookHarness.openFileIds();
            G110WorkbookHarness.armFaults(env, 'CutRead', fullfile(env.work, 'g110_never.xlsx'), ...
                'FailWrite', 'g110_never_', 'RefuseDelete', {'g110_never_'});
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            fst = G110WorkbookHarness.disarmFaults(env);
            leaked = setdiff(G110WorkbookHarness.openFileIds(), fids);
            localCloseFiles(leaked);
            testCase.assertGreaterThan(fst.ferrorCalls, 0, ['harness: none of the helper''s ferror ', ...
                'calls reached the ferror double']);
            testCase.verifyEmpty(r.errText, ['the save failed: ' r.errText]);
            testCase.verifyEqual(readcell(target, 'Range', 'A1:A1'), {'Header'}, 'the file was replaced');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, false);
            testCase.verifyEmpty(p, ['nothing may be left and the decoys untouched: ' localJoin(p)]);
            testCase.verifyFalse(fst.readCut || fst.writeFailed, 'no fault was delivered');
            testCase.verifyEmpty(fst.deleteRefused, 'no delete was refused');
            testCase.verifyEmpty(leaked, 'no file id is left open');
        end

        function testRefusesCarriageReturnAndNamesTheCsvRemedy(testCase)
            %The writer stores a carriage return as a literal character, and XML
            %readers (readcell, openpyxl, and so Excel) read it back as a line
            %feed; only xlsread 'basic' returns it, so the phase-2 check
            %accepted a workbook other readers see differently (round trip of
            %2026-09-23). Until the writer can store it, a carriage return must
            %be refused before writing, the message must name it and point to
            %CSV, and nothing may change. RED on the phase-2 helper, which
            %accepts it.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            cases = {
                'a carriage return in a header cell', {'Header'; ['a' char(13) 'b']}, []
                'a CR LF pair in a table cell',       {'Header'; ''}, ...
                    table({['a' char(13) char(10) 'b']}, 'VariableNames', {'Text'})};
            for k = 1:size(cases, 1)
                what = cases{k, 1};
                G110WorkbookHarness.armSeams(env, name);
                if isempty(cases{k, 3})
                    r = localCall(target, cases{k, 2});
                else
                    r = localCall(target, cases{k, 2}, cases{k, 3}, 2);
                end
                st = G110WorkbookHarness.disarmSeams(env);
                testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:refused', ...
                    sprintf('%s must be refused before writing: %s', what, r.errText));
                testCase.verifyEmpty(st.writerLog, sprintf('%s: nothing may be written', what));
                testCase.verifyEmpty(st.moveLog, sprintf('%s: nothing may be moved', what));
                p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), ...
                    name, true);
                testCase.verifyEmpty(p, sprintf('%s: the folder must be unchanged:\n%s', what, localJoin(p)));
                said = lower(r.errText);
                testCase.verifyTrue(contains(said, 'carriage return') && contains(said, 'csv'), ...
                    sprintf('%s: the message must name the character and the CSV remedy: %s', what, r.errText));
            end
        end

        function testTabLineFeedAndSpacesStayValid(testCase)
            %Control beside the carriage-return refusal: a tab, a line feed and
            %leading and trailing spaces are stored as they are and read back
            %alike by XML readers, so they keep saving.
            env = localEnv(testCase);
            f = fullfile(env.work, 'plain.xlsx');
            T = table({sprintf('a\tb')}, {sprintf('line\nnext')}, {'  two  '}, ...
                'VariableNames', {'Tab', 'Break', 'Spaces'});
            r = localCall(f, {'Header'; ''}, T, 2);
            testCase.assertEmpty(r.errText, ['tab, line feed or spaces were refused: ' r.errText]);
            raw = readcell(f);
            testCase.verifyEqual(raw{3, 1}, sprintf('a\tb'), 'a tab is kept');
            testCase.verifyEqual(raw{3, 2}, sprintf('line\nnext'), 'a line feed is kept');
            testCase.verifyTrue(contains(G110WorkbookHarness.xmlText(f), '>  two  <'), ...
                'leading and trailing spaces are stored');
        end

        % -- phase 3b: one error dialog per failed save, through the real DoD button --

        function testDodExportFailureShowsOneErrorDialog(testCase)
            %The DoD setup screen's "Export DoD Var Components" button
            %(psyrat_startview_dodiff) calls psyrat_dod_varcompt inside a try
            %whose catch shows its own errordlg with the error's message. The
            %helper has already shown "File Not Saved" before it raises, so one
            %failed save shows two error dialogs, the second repeating the first
            %under another title. The failure is delivered through the button's
            %own Callback: the first writer call into the new workbook writes
            %and then fails. Required: the chosen workbook keeps its bytes, the
            %folder is unchanged, and errordlg is called exactly once, naming
            %the file. RED on the phase-2 code: two calls.
            env = localEnv(testCase);
            [spec, fig] = localOpenRoute(testCase, 'dod_varcomp_export');
            name = G110WorkbookHarness.TargetName;
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            G110WorkbookHarness.armSeams(env, name, 'Writer', 'after', 'FailAt', 1);
            G110WorkbookHarness.armErrordlgRecorder(env);
            r = G110WorkbookHarness.fireSave(fig, spec.save, target);
            calls = G110WorkbookHarness.disarmErrordlgRecorder(env);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.assertEqual(r.saveCalls, 1, 'harness: the button did not ask the save dialog once');
            testCase.assertTrue(st.writerInjected, sprintf(['harness: the writer failure was not ', ...
                'delivered:\n%s'], localJoin(G110WorkbookHarness.describeWrites(st))));
            testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(target), original), ...
                'the chosen workbook keeps its bytes');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, ['the folder must be unchanged: ' localJoin(p)]);
            testCase.verifyNumElements(calls, 1, sprintf(['one failed save must show one error ', ...
                'dialog; errordlg was called %d times:\n%s'], numel(calls), localDescribeCalls(calls)));
            testCase.verifyTrue(any(contains({calls.message}, name)), ...
                sprintf('the error dialog must name %s:\n%s', name, localDescribeCalls(calls)));
        end

        function testErrordlgCountControls(testCase)
            %Controls for the count above. (1) The SAME writer
            %(psyrat_dod_varcompt's local_save_table) saved from its own
            %variance window, whose Save Table callback has no catch of its own,
            %with the same writer failure: the chosen workbook keeps its bytes,
            %the helper's error reaches the caller, and errordlg is called once,
            %by the helper, titled "File Not Saved" and naming the file. So the
            %recorder sees the helper's dialog, and the second call in the RED
            %comes from the Export button's catch. (2) A successful save through
            %the Export button calls errordlg zero times, and the planted
            %workbook is replaced (no planted or stale cell, one sheet, the
            %export's first header line in A1), so the zero is not a save that
            %never happened.
            env = localEnv(testCase);
            name = G110WorkbookHarness.TargetName;
            target = fullfile(env.work, name);

            %(1) the writer's own window, the same injected failure
            [spec, fig] = localOpenRoute(testCase, 'dod_varcomp_window');
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            G110WorkbookHarness.armSeams(env, name, 'Writer', 'after', 'FailAt', 1);
            G110WorkbookHarness.armErrordlgRecorder(env);
            r = G110WorkbookHarness.fireSave(fig, spec.save, target);
            calls = G110WorkbookHarness.disarmErrordlgRecorder(env);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.assertEqual(r.saveCalls, 1, 'harness: the window did not ask the save dialog once');
            testCase.assertTrue(st.writerInjected, sprintf(['harness: the writer failure was not ', ...
                'delivered:\n%s'], localJoin(G110WorkbookHarness.describeWrites(st))));
            testCase.verifyTrue(startsWith(r.errText, 'psyrat:replaceWorkbook:check'), ...
                ['the helper''s error reaches the Save Table callback: ' r.errText]);
            testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(target), original), ...
                'the chosen workbook keeps its bytes');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, ['the folder must be unchanged: ' localJoin(p)]);
            testCase.verifyNumElements(calls, 1, sprintf(['the failed save from the variance window ', ...
                'calls errordlg once; it was called %d times:\n%s'], numel(calls), localDescribeCalls(calls)));
            if numel(calls) == 1
                testCase.verifyEqual(calls(1).title, 'File Not Saved', 'the one dialog is the helper''s');
                testCase.verifySubstring(calls(1).caller, 'psyrat_replace_workbook.m', ...
                    'the one dialog is the helper''s');
                testCase.verifySubstring(calls(1).message, name, 'the one dialog names the file');
            end

            %(2) a successful save through the Export button onto the same file
            [spec2, fig2] = localOpenRoute(testCase, 'dod_varcomp_export');
            before2 = G110WorkbookHarness.snapshot(env.work);
            G110WorkbookHarness.armErrordlgRecorder(env);
            r2 = G110WorkbookHarness.fireSave(fig2, spec2.save, target);
            calls2 = G110WorkbookHarness.disarmErrordlgRecorder(env);
            testCase.assertEqual(r2.saveCalls, 1, 'harness: the button did not ask the save dialog once');
            testCase.verifyEmpty(r2.errText, ['the save failed: ' r2.errText]);
            testCase.verifyEmpty(calls2, sprintf('a successful save calls errordlg zero times:\n%s', ...
                localDescribeCalls(calls2)));
            rep = G110WorkbookHarness.staleReport(target);
            testCase.verifyFalse(rep.plant || rep.near || rep.wide || rep.long || rep.sheet, ...
                'the planted workbook was replaced: no planted or stale cell and no old sheet remain');
            testCase.verifyNumElements(rep.sheets, 1, 'the replaced workbook has one sheet');
            C = G110WorkbookHarness.content(target);
            testCase.assertNotEmpty(C, 'the replaced workbook reads back empty');
            testCase.verifyEqual(C{1, 1}, spec2.first, 'A1 is the export''s first header line');
            p2 = G110WorkbookHarness.compareSnapshots(before2, G110WorkbookHarness.snapshot(env.work), ...
                name, false);
            testCase.verifyEmpty(p2, ['nothing may be left and the decoys untouched: ' localJoin(p2)]);
        end

        % -- phase 4: anchoring at A1, and the _xHHHH_ text rule --------------

        function testShiftedWorkbookFailsTheCheck(testCase, Shift)
            %Every write of the new workbook lands the same distance from where
            %it was aimed: the harness writer seam moves each call's 'Range'
            %origin by Shift ([rows columns]). The values and their layout
            %relative to each other are unchanged, so a check that reads the
            %sheet's used range sees the same block (probe of 2026-09-23: all
            %four placements of one payload read back alike through xlsread
            %and readcell; only an explicit A1 read told them apart).
            %Required: the check refuses the new workbook, the chosen file keeps
            %its bytes, nothing is left, and the refusal is shown. RED on the
            %phase-2 helper, whose check reads the used range, so the shifted
            %workbook replaces the chosen file.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            [head, T, namesRow] = localAnchorPayload();
            G110WorkbookHarness.armSeams(env, name, 'Writer', 'shift', 'Shift', Shift);
            r = localCall(target, head, T, namesRow);
            st = G110WorkbookHarness.disarmSeams(env);
            writes = localJoin(G110WorkbookHarness.describeWrites(st));
            testCase.assertGreaterThanOrEqual(st.writerCalls, 1, ...
                'harness: no writer call reached the seam, so nothing was moved');
            testCase.assertEqual(st.shiftApplied, st.writerCalls, sprintf(['harness: every writer ', ...
                'call into the folder must be moved, the header and the table together:\n%s'], writes));
            testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:check', sprintf(['a workbook ', ...
                'that does not start at A1 must fail the check. Writes:\n%s\nError: %s'], writes, r.errText));
            %phase 7: the refusal must come from the anchored read, at A1, not
            %from any other failure of the step that writes and checks
            testCase.verifySubstring(r.errText, ...
                'cell A1 reads an empty cell where the export holds "G110 anchor header"', ...
                'the check must refuse at A1, where the payload''s first cell belongs');
            testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(target), original), ...
                'the chosen workbook keeps its bytes');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, ['the chosen file and the folder must be unchanged: ' localJoin(p)]);
            testCase.verifyNotEmpty(r.dialogText, 'the refusal must be shown');
        end

        function testUnshiftedWriteThroughTheSeamSavesAtA1(testCase)
            %Control for the shift RED: the same payload through the same seam
            %with a zero shift, so the seam rewrites each call's 'Range' origin
            %to itself. The save must succeed, and an explicit-range read
            %(which keeps leading empty rows and columns) must find the payload
            %from A1 with nothing beside or below it. A failure here means the
            %seam, not the helper, changed the write.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            [head, T, namesRow] = localAnchorPayload();
            G110WorkbookHarness.armSeams(env, name, 'Writer', 'shift', 'Shift', [0 0]);
            r = localCall(target, head, T, namesRow);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.assertGreaterThanOrEqual(st.writerCalls, 1, 'harness: no writer call reached the seam');
            testCase.assertEqual(st.shiftApplied, st.writerCalls, ...
                'harness: every writer call must go through the origin rewrite');
            testCase.verifyEmpty(r.errText, ['a write at A1 through the seam was refused: ' r.errText]);
            [~, ~, raw] = xlsread(target, '', 'A1:C4', 'basic'); %#ok<XLSRD> keeps leading empty rows and columns
            expected = {head{1}, NaN, NaN; 'Label', 'Value', NaN; 'p21', 2.5, NaN; NaN, NaN, NaN};
            testCase.verifyTrue(isequaln(raw, expected), ...
                ['A1:C4 must hold the payload from A1 and nothing else; it reads:' newline localShowCells(raw)]);
            rep = G110WorkbookHarness.staleReport(target);
            testCase.verifyFalse(rep.plant || rep.near || rep.wide || rep.long || rep.sheet, ...
                'the planted workbook was replaced');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, false);
            testCase.verifyEmpty(p, ['nothing may be left and the decoys untouched: ' localJoin(p)]);
        end

        function testUppercaseXTextIsSavedExactly(testCase)
            %"_X0041_" with an upper-case X is not the XLSX character escape,
            %whose form has a lower-case x. The writer stores it as written, and
            %xlsread 'basic', readcell and a strict XML parse all return it
            %exactly (probe of 2026-09-23: the shared string holds a_X0041_b).
            %Required: it saves, in a header line and in a table cell, and reads
            %back exactly. RED on the phase-2 helper, whose refusal before
            %writing matches the escape form without regard to case.
            env = localEnv(testCase);
            f = fullfile(env.work, 'upper.xlsx');
            text = 'a_X0041_b';
            T = table({text}, 'VariableNames', {'Text'});
            r = localCall(f, {'Header'; text; ''}, T, 3);
            testCase.verifyEmpty(r.errText, ['text with an upper-case _X0041_ was refused: ' r.errText]);
            if isempty(r.errText)
                [~, ~, raw] = xlsread(f, '', 'A1:A4', 'basic'); %#ok<XLSRD>
                testCase.verifyEqual(raw(:, 1)', {'Header', text, 'Text', text}, ...
                    'xlsread reads the header line and the table cell back exactly');
                testCase.verifyEqual(readcell(f, 'Range', 'A1:A4')', {'Header', text, 'Text', text}, ...
                    'readcell reads them back exactly');
                testCase.verifyTrue(contains(G110WorkbookHarness.xmlText(f), '>a_X0041_b<'), ...
                    'the workbook stores the text as written, not escaped');
            end
        end

        function testLowercaseEscapeTextIsRefusedAsAReaderLimit(testCase)
            %Text in the XLSX escape form, "_x" with four hex digits and "_"
            %(lower-case x), with lower- and upper-case hex digits. The writer
            %stores it correctly by escaping the underscore (a_x005f_x00e9_b),
            %and readcell and a strict XML parse return the original text; only
            %xlsread 'basic', the reader this check uses because readcell trims
            %spaces, returns the escaped form (probe of 2026-09-23). The policy
            %keeps a bounded refusal rather than a new decoder. Required: the
            %text is refused before anything is written, nothing changes, the
            %refusal is shown, and the message places the limit in this save's
            %check, not in the XLSX format or in XLSX readers generally. RED on
            %the phase-2 helper for the message only, which says XLSX readers
            %decode the text; the refusal itself holds.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            cases = {
                'lower-case hex in a header line', {'Header'; 'a_x00e9_b'}, []
                'upper-case hex in a header line', {'Header'; 'a_x000D_b'}, []
                'lower-case hex in a table cell',  {'Header'; ''}, ...
                    table({'a_x00e9_b'}, 'VariableNames', {'Text'})};
            for k = 1:size(cases, 1)
                what = cases{k, 1};
                G110WorkbookHarness.armSeams(env, name);
                if isempty(cases{k, 3})
                    r = localCall(target, cases{k, 2});
                else
                    r = localCall(target, cases{k, 2}, cases{k, 3}, 2);
                end
                st = G110WorkbookHarness.disarmSeams(env);
                testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:refused', ...
                    sprintf('%s must be refused before writing: %s', what, r.errText));
                testCase.verifyEmpty(st.writerLog, sprintf('%s: nothing may be written', what));
                testCase.verifyEmpty(st.moveLog, sprintf('%s: nothing may be moved', what));
                p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), ...
                    name, true);
                testCase.verifyEmpty(p, sprintf('%s: the folder must be unchanged:\n%s', what, localJoin(p)));
                testCase.verifyNotEmpty(r.dialogText, sprintf('%s: the refusal must be shown', what));
                said = lower(r.errText);
                testCase.verifyTrue(contains(said, 'check') || contains(said, 'read back') || ...
                    contains(said, 'reads back') || contains(said, 'verif'), sprintf(['%s: the message ', ...
                    'must say that this save''s check cannot confirm the text: %s'], what, r.errText));
                testCase.verifyFalse(contains(said, 'readers decode') || contains(said, 'xlsx cannot') || ...
                    contains(said, 'cannot be stored') || contains(said, 'cannot store'), sprintf(['%s: ', ...
                    'the message must not state a limit of XLSX or of XLSX readers generally: %s'], ...
                    what, r.errText));
            end
        end

        % -- phase 6: the leftover-copy warning set to 'error' by a caller -----

        function testBackupLeftAfterSuccessStaysSavedWhenItsWarningIsAnError(testCase)
            %The case of testBackupLeftAfterSuccessIsReported (the same refused
            %delete of the checked copy after a verified replacement), with ONLY
            %the helper's cleanup warning, psyrat:replaceWorkbook:cleanup, set to
            %'error', as a caller may set it. That warning is one way the copy
            %left behind is disclosed, so the setting must not turn the finished
            %save into a failure. Required: no error returned; the chosen file
            %holds the new checked payload from A1; the copy is the only new
            %entry, keeps the original bytes, and its full path is shown in the
            %File Saved dialog and on the console (the helper prints the notice
            %to fid 2 when its warning is an error; the call is captured with
            %evalc, after a probe shows that evalc captures fid 2 here); no File
            %Not Saved dialog; and the caller's warning settings are as they
            %were. The setting is restored right after the save and again at
            %teardown. RED on the phase-5 helper, whose warning call is not
            %guarded, so it raises after the move.
            id = 'psyrat:replaceWorkbook:cleanup';
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            %positive control for the console check: text written to fid 2
            %must reach an evalc capture at all
            probe = evalc('fprintf(2, ''G110 console probe\n'');');
            saved = warning('query', id);
            testCase.addTeardown(@() warning(saved));
            warning('error', id);
            armed = warning('query', id);
            testCase.assertEqual(armed.state, 'error', 'harness: the cleanup warning was not set to error');
            G110WorkbookHarness.armFaults(env, 'RefuseDelete', {'psyrat_backup_'});
            w = warning;
            allBefore = w(strcmp({w.identifier}, 'all')).state;
            r = [];
            doSave = @() localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            console = evalc('r = doSave();');
            stateNow = warning('query', id);
            w = warning;
            allAfter = w(strcmp({w.identifier}, 'all')).state;
            warning(saved);
            fst = G110WorkbookHarness.disarmFaults(env);
            testCase.assertNotEmpty(fst.deleteRefused, 'harness: the delete of the copy was not refused');
            testCase.verifyEmpty(r.errText, ['a verified save must not end in an error: ' r.errText]);
            testCase.verifyEqual(stateNow.state, 'error', ...
                'the caller''s setting for the cleanup warning must be left as it was');
            testCase.verifyEqual(allAfter, allBefore, 'the caller''s default warning state must be left as it was');
            testCase.verifySubstring(probe, 'G110 console probe', ...
                'harness: evalc must capture text written to fid 2, or the console check below means nothing');
            [~, ~, raw] = xlsread(target, '', 'A1:C4', 'basic'); %#ok<XLSRD>
            expected = {'Header', NaN, NaN; 'L', 'V', NaN; 'a', 1, NaN; NaN, NaN, NaN};
            testCase.verifyTrue(isequaln(raw, expected), ['the chosen file must hold the new ', ...
                'checked payload from A1; A1:C4 reads:' newline localShowCells(raw)]);
            rep = G110WorkbookHarness.staleReport(target);
            testCase.verifyFalse(rep.plant || rep.near || rep.wide || rep.long || rep.sheet, ...
                'the planted workbook was replaced');
            folderAfter = G110WorkbookHarness.snapshot(env.work);
            added = setdiff({folderAfter.name}, {before.name});
            testCase.verifyNumElements(added, 1, ['exactly one new entry, the copy: ' localJoin(added)]);
            if numel(added) == 1
                kept = fullfile(env.work, added{1});
                testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(kept), original), ...
                    'the copy left behind holds the original bytes');
                %the dialog wraps long lines, so compare without white space
                squash = @(s) regexprep(s, '\s', '');
                testCase.verifyTrue(any(strcmp(r.dialogs, 'File Saved')) && ...
                    contains(squash(r.dialogText), squash(kept)), sprintf(['the File Saved dialog ', ...
                    'must give the copy''s full path %s; windows: %s; text: %s'], kept, ...
                    localJoin(r.dialogs), r.dialogText));
                if contains(probe, 'G110 console probe')
                    testCase.verifyTrue(contains(squash(console), squash(kept)), sprintf(['the console ', ...
                        'must give the copy''s full path %s when the warning is an error; it shows: %s'], ...
                        kept, console));
                end
            end
            testCase.verifyFalse(any(strcmp(r.dialogs, 'File Not Saved')) || ...
                contains(r.dialogText, 'was not saved'), 'a verified save must not show File Not Saved');
            testCase.verifyEmpty(localChanged(before, folderAfter, name), 'every decoy is untouched');
        end

        % -- phase 7: pins of the one-sheet, used-range and moved-file checks --

        function testExtraSheetOrCellFailsTheCheck(testCase, Extra)
            %After each forwarded write of the new workbook, the harness writer
            %seam's extra mode adds a second sheet, or one text cell just
            %outside the payload (A1:B4): right of its widest column on row 1,
            %or on the row below its last row. The payload itself is written in
            %its places, so the explicit read from A1 over its extent passes;
            %only the one-sheet check (a second sheet) or the used-range read
            %(a cell outside the extent) can see the addition. Required: the
            %check refuses the new workbook in that clause's own words, the
            %chosen file keeps its bytes, nothing is left, and the refusal is
            %shown. A pin of clauses the helper already enforces: it passes on
            %the helper and fails on a copy without the sheet count or without
            %the used-range read.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            [head, T, namesRow] = localNeighborPayload();
            G110WorkbookHarness.armSeams(env, name, 'Writer', 'extra', 'Extra', Extra);
            r = localCall(target, head, T, namesRow);
            st = G110WorkbookHarness.disarmSeams(env);
            writes = localJoin(G110WorkbookHarness.describeWrites(st));
            testCase.assertGreaterThanOrEqual(st.writerCalls, 1, 'harness: no writer call reached the seam');
            testCase.assertEqual(st.extraApplied, st.writerCalls, sprintf(['harness: the addition must ', ...
                'follow every writer call into the folder:\n%s'], writes));
            if strcmp(Extra, 'sheet')
                testCase.assertEqual(st.extraSheets, 2, sprintf(['harness: the new workbook must ', ...
                    'have held two sheets:\n%s'], writes));
                wording = 'it has 2 sheets instead of 1';
            else
                wording = sprintf('cell %s reads "G110_EXTRA" where the export holds an empty cell', Extra);
            end
            testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:check', sprintf(['a new workbook ', ...
                'holding more than the payload must fail the check. Writes:\n%s\nError: %s'], ...
                writes, r.errText));
            testCase.verifySubstring(r.errText, wording, 'the check must name what it found');
            testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(target), original), ...
                'the chosen workbook keeps its bytes');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, ['the chosen file and the folder must be unchanged: ' localJoin(p)]);
            testCase.verifyNotEmpty(r.dialogText, 'the refusal must be shown');
        end

        function testExtraSeamAddingNothingKeepsValidEmptyNeighbors(testCase)
            %Control for the extra-content pin: the same payload through the
            %same seam and hook, with nothing added. The payload's own empty
            %cells, the blank beside its header line (B1) and the NaN that ends
            %its last row (B4, stored as an empty cell), lie next to where the
            %pin writes (C1, A5) and must not be taken for missing or extra
            %content. Required: the save succeeds with one sheet, and an
            %explicit read from A1 finds the payload and nothing else, the two
            %cells the pin writes included. A failure here means the seam or
            %the payload, not the helper, is at fault.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            [head, T, namesRow] = localNeighborPayload();
            G110WorkbookHarness.armSeams(env, name, 'Writer', 'extra', 'Extra', 'none');
            r = localCall(target, head, T, namesRow);
            st = G110WorkbookHarness.disarmSeams(env);
            testCase.assertGreaterThanOrEqual(st.writerCalls, 1, 'harness: no writer call reached the seam');
            testCase.assertEqual(st.extraCalls, st.writerCalls, ...
                'harness: every writer call must go through the addition hook');
            testCase.assertEqual(st.extraApplied, 0, 'harness: the control must add nothing');
            testCase.verifyEmpty(r.errText, ['a workbook holding exactly the payload was refused: ' r.errText]);
            [~, ~, raw] = xlsread(target, '', 'A1:C6', 'basic'); %#ok<XLSRD> keeps leading empty rows and columns
            expected = {head{1}, NaN, NaN; 'Label', 'Value', NaN; 'p21', 2.5, NaN; 'p31', NaN, NaN; ...
                NaN, NaN, NaN; NaN, NaN, NaN};
            testCase.verifyTrue(isequaln(raw, expected), ...
                ['A1:C6 must hold the payload from A1 and nothing else; it reads:' newline localShowCells(raw)]);
            rep = G110WorkbookHarness.staleReport(target);
            testCase.verifyFalse(rep.plant || rep.near || rep.wide || rep.long || rep.sheet, ...
                'the planted workbook was replaced');
            testCase.verifyNumElements(rep.sheets, 1, 'the new workbook has one sheet');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, false);
            testCase.verifyEmpty(p, ['nothing may be left and the decoys untouched: ' localJoin(p)]);
        end

        function testMoveReportingSuccessWithWrongBytesRestoresTheOriginal(testCase)
            %The first move of the new workbook onto the chosen file is carried
            %out and reported as a success, but what it leaves under that name
            %is other bytes (the harness move seam's wrongButReportsSuccess
            %mode). Only reading the chosen file back after the move can tell.
            %Required: no success; the message says the file in place is not
            %the checked workbook and that the original was put back from its
            %copy and checked; the chosen file ends byte-identical to the
            %original; nothing is left; the failure is shown as File Not Saved.
            %A pin of a clause the helper already enforces: it passes on the
            %helper and fails on a copy that trusts the move's reported success.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            G110WorkbookHarness.armSeams(env, name, 'Replace', 'wrongButReportsSuccess');
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            st = G110WorkbookHarness.disarmSeams(env);
            moves = localJoin(G110WorkbookHarness.describeWrites(st));
            testCase.assertTrue(st.replaceInjected, sprintf(['harness: no move onto the chosen file ', ...
                'reached the seam:\n%s'], moves));
            testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:replace', sprintf(['a move that left ', ...
                'other bytes must not count as a save:\n%s\nError: %s'], moves, r.errText));
            testCase.verifySubstring(r.errText, 'the file in place afterwards is not the checked workbook', ...
                'the message must say why the move was not accepted');
            testCase.verifySubstring(r.errText, sprintf(['The original %s was put back from its copy ', ...
                'and checked byte for byte.'], name), 'the message must say the original was put back and checked');
            testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(target), original), ...
                'the chosen file ends byte-identical to the original');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, ['the chosen file must be the original and nothing left: ' localJoin(p)]);
            testCase.verifyTrue(any(strcmp(r.dialogs, 'File Not Saved')), ...
                ['the failure must be shown as File Not Saved; windows: ' localJoin(r.dialogs)]);
        end

        function testMoveReportingSuccessWithWrongBytesToANewNameIsReported(testCase)
            %As above, but no file had the chosen name, so there is no original
            %to put back. Required: no success; the message says the file in
            %place is not the checked workbook and that a file of that name now
            %exists and must not be relied on; that file is the only new entry
            %and is still there; the decoys are untouched; the failure is shown
            %as File Not Saved.
            env = localEnv(testCase);
            name = 'chosen.xlsx';
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            G110WorkbookHarness.armSeams(env, name, 'Replace', 'wrongButReportsSuccess');
            r = localCall(target, {'Header'; ''}, table({'a'}, 1, 'VariableNames', {'L', 'V'}), 2);
            st = G110WorkbookHarness.disarmSeams(env);
            moves = localJoin(G110WorkbookHarness.describeWrites(st));
            testCase.assertTrue(st.replaceInjected, sprintf(['harness: no move onto the chosen name ', ...
                'reached the seam:\n%s'], moves));
            testCase.verifyEqual(r.errId, 'psyrat:replaceWorkbook:replace', sprintf(['a move that left ', ...
                'other bytes must not count as a save:\n%s\nError: %s'], moves, r.errText));
            testCase.verifySubstring(r.errText, 'the file in place afterwards is not the checked workbook', ...
                'the message must say why the move was not accepted');
            testCase.verifySubstring(r.errText, sprintf(['A file named %s now exists in that folder; it ', ...
                'did not pass the check, so do not rely on it.'], name), ...
                'the message must say the file now there is not to be relied on');
            after = G110WorkbookHarness.snapshot(env.work);
            testCase.verifyEqual(setdiff({after.name}, {before.name}), {name}, ...
                'the only new entry is the file under the chosen name');
            testCase.verifyTrue(isfile(target), 'the file the message warns about is there');
            testCase.verifyEmpty(localChanged(before, after), 'every decoy is untouched');
            testCase.verifyTrue(any(strcmp(r.dialogs, 'File Not Saved')), ...
                ['the failure must be shown as File Not Saved; windows: ' localJoin(r.dialogs)]);
        end
    end
end

% -- local helpers -------------------------------------------------------------

function env = localEnv(testCase)
env = G110WorkbookHarness.install();
testCase.addTeardown(@() G110WorkbookHarness.uninstall(env));
end

function r = localCall(varargin)
%Call the helper as a callback does, and record its error and the message it
%showed. A modal errordlg reuses an open dialog of the same title, so the
%message is read from the message objects created during this call (in a new
%or a reused window), never from the set of new windows alone; the windows
%this call opened are closed after they are read, so each call starts clean.
mark = G110WorkbookHarness.dialogMark();
r = struct('errId', '', 'errText', '', 'dialogText', '', 'dialogs', {{}});
try
    psyrat_replace_workbook(varargin{:});
catch ME
    r.errId = ME.identifier;
    r.errText = ME.message;
end
[r.dialogText, r.dialogs] = G110WorkbookHarness.dialogsSince(mark, true);
end

function names = localEntries(folder)
L = dir(folder);
names = sort(setdiff({L.name}, {'.', '..'}));
end

function p = localChanged(before, after, skip)
%Entries of BEFORE (other than SKIP) that are missing from AFTER or whose
%bytes changed.
if nargin < 3
    skip = '';
end
p = {};
for k = 1:numel(before)
    if strcmp(before(k).name, skip)
        continue;
    end
    j = find(strcmp({after.name}, before(k).name), 1);
    if isempty(j)
        p{end+1} = [before(k).name ' is gone']; %#ok<AGROW>
    elseif ~isequal(before(k).bytes, after(j).bytes)
        p{end+1} = [before(k).name ' changed']; %#ok<AGROW>
    end
end
end

function s = localJoin(c)
if isempty(c)
    s = '';
else
    s = strjoin(reshape(c, 1, []), newline);
end
end

function localVerifyLeftovers(testCase, env, before, name, original, r, refused, id)
%Shared outcome of the leftover tests: the refusal id; the chosen file and every
%decoy unchanged; the only new entries are the files whose removal was refused;
%the one error gives each one's full path, and nothing says it "was removed".
testCase.verifyEqual(r.errId, id, r.errText);
testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(fullfile(env.work, name)), original), ...
    'the chosen file keeps its bytes');
after = G110WorkbookHarness.snapshot(env.work);
added = sort(setdiff({after.name}, {before.name}));
left = sort(unique(cellfun(@G110WorkbookHarness.nameOf, refused, 'UniformOutput', false)));
testCase.verifyEqual(added, left, 'the only new entries are the files whose removal was refused');
testCase.verifyEmpty(localChanged(before, after, name), 'every decoy is untouched');
testCase.verifyNotEmpty(r.dialogText, 'the failure must be shown');
said = [r.errText ' ' r.dialogText];
testCase.verifyFalse(contains(said, 'was removed'), ...
    ['a file that is still there must not be reported as removed: ' r.errText]);
for k = 1:numel(left)
    testCase.verifyTrue(contains(r.errText, fullfile(env.work, left{k})), ...
        sprintf('the error must say where %s is: %s', left{k}, r.errText));
end
end

function localCloseFiles(fids)
%Close file ids a faulty call left open, so later tests start clean.
for k = 1:numel(fids)
    if ~isempty(fopen(fids(k)))
        fclose(fids(k));
    end
end
end

function localWriteBytes(file, bytes)
fid = fopen(file, 'w');
fwrite(fid, bytes, 'uint8');
fclose(fid);
end

function [spec, fig] = localOpenRoute(testCase, key)
%Open a route's window through the harness (the app's own screen and button),
%and close every window this test opened from here on (the route's windows and
%their dialogs) at teardown. A route that does not open fails this test only.
figsBefore = findall(groot, 'Type', 'figure');
testCase.addTeardown(@() localCloseNewFigures(figsBefore));
spec = G110WorkbookHarness.route(key);
fig = [];
try
    fig = G110WorkbookHarness.openRoute(spec);
catch ME
    testCase.assertFail(sprintf('%s: harness: the route did not open: %s', key, ME.message));
end
end

function localCloseNewFigures(figsBefore)
figs = findall(groot, 'Type', 'figure');
opened = figs(~ismember(figs, figsBefore));
delete(opened(isvalid(opened)));
end

function [head, T, namesRow] = localAnchorPayload()
%A small header and table passed the way the callbacks pass them: the header's
%last line is blank and the table's variable names go on that row, so the new
%workbook must hold the header line in A1, the names in A2:B2 and one row below.
head = {'G110 anchor header'; ''};
T = table({'p21'}, 2.5, 'VariableNames', {'Label', 'Value'});
namesRow = 2;
end

function [head, T, namesRow] = localNeighborPayload()
%A payload whose own empty cells sit at its edges, passed the way the callbacks
%pass one: B1 is blank beside the header line, and the last row ends in NaN,
%which the writer stores as an empty cell. The new workbook must hold the header
%line in A1, the names in A2:B2 and two rows below: A1:B4, with B1 and B4 empty.
head = {'G110 neighbor header'; ''};
T = table({'p21'; 'p31'}, [2.5; NaN], 'VariableNames', {'Label', 'Value'});
namesRow = 2;
end

function s = localShowCells(raw)
%A small cell array as text, one row per line, for a diagnostic.
lines = cell(1, size(raw, 1));
for r = 1:size(raw, 1)
    parts = cell(1, size(raw, 2));
    for c = 1:size(raw, 2)
        v = raw{r, c};
        if ischar(v)
            parts{c} = ['"' v '"'];
        elseif isnumeric(v) && isscalar(v)
            parts{c} = sprintf('%.17g', v);
        else
            parts{c} = ['<' class(v) '>'];
        end
    end
    lines{r} = strjoin(parts, ' | ');
end
s = localJoin(lines);
end

function s = localDescribeCalls(calls)
%One line per recorded errordlg call: where it was called, its title, its
%message on one line.
lines = cell(1, numel(calls));
for k = 1:numel(calls)
    lines{k} = sprintf('%d. %s, "%s": %s', k, calls(k).caller, calls(k).title, ...
        strrep(calls(k).message, newline, ' | '));
end
s = localJoin(lines);
end
