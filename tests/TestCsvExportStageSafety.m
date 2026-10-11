classdef TestCsvExportStageSafety < PsyRATTestBase
    %G109: psyrat_export_stage called directly, the way the three CSV Save
    %callbacks call it (localWriter: open, onCleanup abort, write, commit, as
    %in psyrat_variancet.m psyrat_savesscoeffs and psyrat_savealldifftable
    %and psyrat_criterionfigures.m local_save_table). The callback round
    %trips at all three sites are in TestSubjectDifferenceExport. These tests
    %pin the helper's own branches that no callback reaches on demand
    %(export-review.phase3 F1, F2, F3 and F5):
    %  - the final move failing by status 0 and by an exception;
    %  - a temporary file that cannot be removed after a failure;
    %  - a folder created with the target's name between open and commit;
    %  - a target named like a temporary file, in either letter case.
    %
    %SEAMS (G107ExportFixtures). installIoFaults puts recording fclose,
    %ferror, movefile and delete doubles first on the path. Each forwards to
    %its built-in and fails only a budgeted number of calls on files under
    %the test folder; movefile can raise instead of returning status 0, and
    %delete can decline to delete. installErrordlgRecorder records every
    %errordlg call, so "exactly one dialog" counts calls, not windows. The
    %doubles are removed before any file is read back. Nothing here opens a
    %PsyRAT window or fits a model.
    %
    %BEFORE THE CORRECTION. testFinalMoveFailureReportsOneReplaceFailure with
    %MoveFailure = thrown fails on the phase-5 helper (F1). Every other test
    %here passes on it; each comment names the guard it depends on.

    properties (TestParameter)
        %how the final move fails: status 0 (documented) or an exception
        MoveFailure = struct('statusFalse', false, 'thrown', true);
        %what failed before the temporary file could not be removed
        CleanupAfter = struct('refusedMove', 'refusedMove', 'writeFailure', 'writeFailure');
    end

    properties (SetAccess = private)
        TmpDir
        DialogDir
    end

    methods (TestMethodSetup)
        function setupFolderAndRecorder(testCase)
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);
            testCase.DialogDir = G107ExportFixtures.installErrordlgRecorder();
            testCase.assertTrue(startsWith(which('errordlg'), testCase.DialogDir), ...
                'the recording errordlg must be the resolved errordlg');
            delete(findall(groot, 'Type', 'figure'));
        end
    end

    methods (TestMethodTeardown)
        function teardownFolderAndRecorder(testCase)
            delete(findall(groot, 'Type', 'figure'));
            G107ExportFixtures.removeErrordlgRecorder(testCase.DialogDir);
            if isfolder(testCase.TmpDir)
                rmdir(testCase.TmpDir, 's');
            end
        end
    end

    methods (Test)

        function testFinalMoveFailureReportsOneReplaceFailure(testCase, MoveFailure)
            %The final move onto an existing kept.csv fails: movefile returns
            %status 0 (statusFalse, MATLAB's documented failure signal) or
            %raises (thrown, F1). The double moves nothing, so kept.csv must
            %keep its bytes. The temporary file must be removed (the delete
            %double at budget 0 forwards that delete, the control for that
            %seam), nothing may stay open, the helper must raise
            %psyrat:exportStage:replace, and exactly one dialog must say the
            %new file could not replace kept.csv, ask the user to check it,
            %report the observed removal, and say nothing the write-failure
            %dialog says. EXPECTED: statusFalse passes on the phase-5 helper.
            %thrown fails there: the exception leaves commit before its
            %report, so the caller's abort reports that the new file "could
            %not be written completely, so it did not replace any existing
            %file", and the error that reaches the caller is the double's.
            target = fullfile(testCase.TmpDir, 'kept.csv');
            G107ExportFixtures.plant(target, sprintf('PLANTED SENTINEL G109\n'));
            bytes0 = G107ExportFixtures.fileBytes(target);
            R = localRun(testCase, target, @localNoHook, struct('fclose', 0, ...
                'ferror', 0, 'movefile', 1, 'movefileThrows', MoveFailure, 'delete', 0));
            testCase.verifyNumElements(R.io.movefileFailed, 1, ...
                'seam: the final move was refused once');
            testCase.verifyNumElements(R.io.deleteSeen, 1, ...
                'seam: the helper removed its temporary file through delete');
            localVerifyError(testCase, R, 'psyrat:exportStage:replace');
            txt = localOneDialog(testCase, R);
            G107ExportFixtures.verifyReplaceFailureText(testCase, txt, 'kept.csv', 'final move');
            testCase.verifyTrue(contains(txt, 'was removed'), ...
                'the dialog reports the removal of the temporary file, which was observed');
            testCase.verifyEqual(G107ExportFixtures.fileBytes(target), bytes0, ...
                'kept.csv keeps its bytes (the double moved nothing)');
            tree = G107ExportFixtures.folderTree(testCase.TmpDir);
            testCase.verifyEqual(tree.paths, {'kept.csv'}, ...
                'the temporary file is gone and nothing else was written');
            testCase.verifyEmpty(R.leaked, 'nothing is left open');
        end

        function testStageCleanupFailureLeavesOwnStageAndSaysSoOnce(testCase, CleanupAfter)
            %After a failure the temporary file cannot be removed: the delete
            %double declines its one delete and deletes nothing.
            %refusedMove: the final move was refused (status 0).
            %writeFailure: the writer failed after writing, before commit.
            %The target is an existing kept.csv, beside the three stage-named
            %decoys (G107ExportFixtures.plantStageDecoys). Expected:
            %  - exactly one dialog, which names the leftover temporary file
            %    and says it could not be removed (never that it "was
            %    removed");
            %  - the leftover is this export's own file, still holding what
            %    was written, left for the user;
            %  - the decoys and kept.csv keep every byte, and nothing stays
            %    open;
            %  - the dialog is the replace-failure report after a refused
            %    move and the write-failure report after a writer failure.
            %refusedMove is the one branch where the helper's 'settled' list
            %decides the dialog count: the leftover is still on disk when the
            %caller's abort runs, so without the list abort would add a
            %second, contradictory dialog (F2). Passes on the phase-5 helper.
            target = fullfile(testCase.TmpDir, 'kept.csv');
            G107ExportFixtures.plant(target, sprintf('PLANTED SENTINEL G109\n'));
            decoys = G107ExportFixtures.plantStageDecoys(testCase.TmpDir);
            tree0 = G107ExportFixtures.folderTree(testCase.TmpDir);
            switch CleanupAfter
                case 'refusedMove'
                    hook = @localNoHook;
                    fails = struct('fclose', 0, 'ferror', 0, 'movefile', 1, 'delete', 1);
                    wantError = 'psyrat:exportStage:replace';
                case 'writeFailure'
                    hook = @localFailingWriter;
                    fails = struct('fclose', 0, 'ferror', 0, 'movefile', 0, 'delete', 1);
                    wantError = 'G109Test:writerFailed';
            end
            R = localRun(testCase, target, hook, fails);
            testCase.verifyNumElements(R.io.deleteFailed, 1, ...
                'seam: the one delete of the temporary file was declined');
            localVerifyError(testCase, R, wantError);
            tree1 = G107ExportFixtures.folderTree(testCase.TmpDir);
            testCase.assertEmpty(setdiff(tree0.paths, tree1.paths), ...
                'nothing that was there before was removed');
            added = setdiff(tree1.paths, tree0.paths);
            testCase.assertNumElements(added, 1, 'exactly one file is left over');
            leftover = added{1};
            testCase.verifyTrue(~isempty(regexp(leftover, '^psyrat_partial_\d+\.csv$', 'once')) && ...
                ~any(strcmp(leftover, decoys)), sprintf(['the leftover (%s) is a temporary ' ...
                'file of this export, not a decoy'], leftover));
            if ~isempty(R.io.deleteFailed)
                [~, n, e] = fileparts(R.io.deleteFailed{1});
                testCase.verifyEqual([n e], leftover, 'the declined delete was of the leftover file');
            end
            testCase.verifyTrue(contains(char(G107ExportFixtures.fileBytes( ...
                fullfile(testCase.TmpDir, leftover))'), 'G109 STAGE CONTENT'), ...
                'the leftover still holds what this export wrote, for the user');
            [~, keep] = ismember(tree0.paths, tree1.paths);
            testCase.verifyEqual(tree1.bytes(keep), tree0.bytes, ...
                'kept.csv and every decoy keep their bytes');
            txt = localOneDialog(testCase, R);
            testCase.verifyTrue(contains(txt, [leftover ' could not be removed']), ...
                'the dialog names the leftover and says it could not be removed');
            testCase.verifyFalse(contains(txt, 'was removed'), ...
                'the dialog claims no removal that did not happen');
            switch CleanupAfter
                case 'refusedMove'
                    G107ExportFixtures.verifyReplaceFailureText(testCase, txt, 'kept.csv', CleanupAfter);
                case 'writeFailure'
                    G107ExportFixtures.verifyWriteFailureText(testCase, txt, CleanupAfter);
            end
            testCase.verifyEmpty(R.leaked, 'nothing is left open');
        end

        function testFolderCreatedBeforeCommitIsRefusedOnce(testCase)
            %F5. A folder with the target's name appears between open and
            %commit. The destination check had passed because late.csv did
            %not exist; localWriter's hook makes the folder, with a sentinel
            %inside, after the data were written. The commit backstop must
            %refuse before any move:
            %  - the zero-budget movefile double sees no move;
            %  - the folder keeps only its sentinel;
            %  - the helper's own temporary file is removed (the delete
            %    double at budget 0 forwards it);
            %  - psyrat:exportStage:folder is raised and nothing stays open;
            %  - exactly one dialog says late.csv is now a folder and reports
            %    the removal, with none of the write-failure or
            %    replace-failure text.
            %No callback test reaches this branch, because the callbacks'
            %folder refusals stop in psyrat_export_destination. Without it
            %the move would put the finished file inside the folder. Passes
            %on the phase-5 helper.
            target = fullfile(testCase.TmpDir, 'late.csv');
            sentinel = sprintf('G109 FOLDER SENTINEL\n');
            R = localRun(testCase, target, @(~) localMakeFolder(target, sentinel), ...
                struct('fclose', 0, 'ferror', 0, 'movefile', 0, 'delete', 0));
            testCase.verifyEmpty(R.io.movefileSeen, 'no move was attempted onto the folder');
            testCase.verifyNumElements(R.io.deleteSeen, 1, ...
                'seam: the helper removed its temporary file through delete');
            localVerifyError(testCase, R, 'psyrat:exportStage:folder');
            txt = localOneDialog(testCase, R);
            testCase.verifyTrue(contains(txt, '"late.csv" is now a folder'), ...
                'the dialog says late.csv is now a folder');
            testCase.verifyTrue(contains(txt, 'was removed'), ...
                'the dialog reports the observed removal of the temporary file');
            testCase.verifyFalse(contains(txt, 'could not be written completely') || ...
                contains(txt, 'could not replace'), ...
                'the dialog is neither the write-failure nor the replace-failure report');
            tree = G107ExportFixtures.folderTree(testCase.TmpDir);
            testCase.verifyEqual(tree.paths, {['late.csv' filesep]; fullfile('late.csv', 'inside.txt')}, ...
                'nothing was written into or beside the folder, and the temporary file is gone');
            testCase.verifyEqual(G107ExportFixtures.fileBytes(fullfile(target, 'inside.txt')), ...
                uint8(sentinel(:)), 'the folder''s sentinel keeps its bytes');
            testCase.verifyEmpty(R.leaked, 'nothing is left open');
        end

        function testTargetNamedLikeATemporaryFileIsWritten(testCase)
            %F3. Two NEW targets, each in a folder of its own:
            %psyrat_partial_1.csv, then PSYRAT_PARTIAL_1.CSV. Each save must
            %write exactly that file, under the name given, with the content,
            %and leave nothing else, with no error and no dialog. The
            %zero-budget movefile double shows that the file moved onto the
            %target was not the target itself (compared without regard to
            %case).
            %This pins the target-equality skip in 'open'. Without it the
            %lower-case target becomes its own temporary file: the move onto
            %itself fails and the cleanup deletes the new file. The
            %upper-case target matters on a volume that ignores letter case
            %(the macOS default), where both spellings name one file. The
            %name comparison above holds on any volume, but what MATLAB's
            %movefile does with a case-only move of one file is not observed
            %here. Passes on the phase-5 helper.
            names = {'psyrat_partial_1.csv', 'PSYRAT_PARTIAL_1.CSV'};
            for k = 1:numel(names)
                what = names{k};
                sub = fullfile(testCase.TmpDir, sprintf('named_%d', k));
                mkdir(sub);
                target = fullfile(sub, names{k});
                R = localRun(testCase, target, @localNoHook, ...
                    struct('fclose', 0, 'ferror', 0, 'movefile', 0));
                testCase.verifyEmpty(R.out.errorText, [what ': the save must not throw']);
                testCase.verifyEmpty(R.dialogCalls, [what ': no dialog']);
                testCase.assertNumElements(R.io.movefileSeen, 1, [what ': one finished ' ...
                    'file was moved onto the target (seam fired)']);
                testCase.verifyFalse(strcmpi(R.io.movefileSeen{1}, target), ...
                    [what ': the target was not its own temporary file']);
                tree = G107ExportFixtures.folderTree(sub);
                testCase.verifyEqual(tree.paths, names(k), ...
                    [what ': exactly the target, under the name given, and nothing else']);
                testCase.verifyTrue(contains(char(G107ExportFixtures.fileBytes(target)'), ...
                    'G109 STAGE CONTENT'), [what ': the target holds the content']);
                testCase.verifyEmpty(R.leaked, [what ': nothing is left open']);
            end
        end

    end
end

% -- the caller, the hooks and the checks ---------------------------------------

function localWriter(target, hook)
%A caller written the way the three CSV Save callbacks use the helper: open,
%stop quietly on a refused open, register abort with onCleanup, write, HOOK
%(the moment between writing and commit), commit.
[fid, stage] = psyrat_export_stage('open', target);
if fid < 0
    return;
end
guard = onCleanup(@() psyrat_export_stage('abort', fid, stage, target)); %#ok<NASGU>
fprintf(fid, '%s\n', 'G109 STAGE CONTENT');
hook(stage);
psyrat_export_stage('commit', fid, stage, target);
end

function localNoHook(~)
%Nothing happens between writing and commit.
end

function localFailingWriter(~)
%The writer fails after its data were written, before commit.
error('G109Test:writerFailed', 'Injected writer failure after the data were written (G109 test).');
end

function localMakeFolder(folder, text)
%A folder with the target's name appears between open and commit, holding a
%sentinel file.
mkdir(folder);
G107ExportFixtures.plant(fullfile(folder, 'inside.txt'), text);
end

function R = localRun(testCase, target, hook, fails)
%Run localWriter(TARGET, HOOK) with the I/O doubles at the budgets in FAILS
%(G107ExportFixtures.installIoFaults) and the errordlg recorder reset. Returns
%what the run raised and showed (R.out, from G107ExportFixtures.capture), every
%errordlg call (R.dialogCalls), the doubles' log (R.io) and the file IDs left
%open (R.leaked; closed here once the doubles are gone).
doubled = intersect({'fclose', 'ferror', 'movefile', 'delete'}, fieldnames(fails)');
for f = doubled
    testCase.assertEqual(exist(f{1}, 'builtin'), 5, sprintf( ...
        'precondition: %s is a built-in the double can forward to', f{1}));
end
w = G107ExportFixtures.installIoFaults(testCase.TmpDir, fails);
testCase.addTeardown(@() G107ExportFixtures.removeIoFaults(w));
for f = doubled
    testCase.assertTrue(startsWith(which(f{1}), w), sprintf( ...
        'the %s double must be the resolved %s', f{1}, f{1}));
end
G107ExportFixtures.resetErrordlgCalls();
fids0 = G107ExportFixtures.openFileIDs();
R.out = G107ExportFixtures.capture(@() localWriter(target, hook));
R.dialogCalls = G107ExportFixtures.errordlgCalls();
R.io = G107ExportFixtures.ioFaultLog();
R.leaked = setdiff(G107ExportFixtures.openFileIDs(), fids0);
G107ExportFixtures.removeIoFaults(w);
for f = R.leaked(:)'
    fclose(f);   % after the doubles are gone; a leak must not outlive the test
end
end

function localVerifyError(testCase, R, id)
%The run raised, so the caller stops, and the error that reached the caller is
%ID.
testCase.verifyNotEmpty(R.out.error, 'the helper raises an error, so the caller stops');
if ~isempty(R.out.error)
    testCase.verifyEqual(R.out.error.identifier, id, 'the error that reaches the caller');
end
end

function txt = localOneDialog(testCase, R)
%Exactly one errordlg call. Returns the text of every recorded call
%(G107ExportFixtures.errordlgText), so a phrase check also sees a second
%dialog if there is one.
testCase.verifyNumElements(R.dialogCalls, 1, 'exactly one dialog');
txt = G107ExportFixtures.errordlgText(R.dialogCalls);
end
