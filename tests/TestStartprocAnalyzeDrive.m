classdef TestStartprocAnalyzeDrive < PsyRATTestBase
    %Beta-tester feedback of 2026-09-16 (G92 sentinel crash, G93 output
    %naming, G94 error-report naming, and the GUI half of G96): the first
    %tests to build the Specify Inputs screen headlessly and fire its Analyze
    %callback.
    %
    %ROUTE. psyrat_startproc_gui and psyrat_exec are local functions of
    %psyrat_startproc.m, so the screen is entered the way a user enters it:
    %psyrat_startproc with a uigetfile shadow that returns the fixture CSV.
    %psyrat_startproc honors a passed psyrat_prefs only while a figure tagged
    %psyrat_gui (the home screen) is open, so each drive opens a stand-in home
    %figure first; the entry closes it. A uiputfile shadow then answers the
    %save dialog: 0 (Cancel) stops psyrat_exec right after every input check
    %with a 'File Error' dialog, and a path lets the measurement loop run into
    %psyrat_loadfile and psyrat_preflight_validate, where a NaN score fails the
    %run inside the catch that writes the per-measurement error report. No
    %CmdStan is ever reached.
    %
    %OBSERVABLES. A crash is an error escaping the callback (the tester's
    %"Unrecognized table variable name '(none)'"); the cancel route leaves a
    %'File Error' dialog, which proves every check between the facet blocks
    %and the save dialog ran; the report route leaves .txt files in the save
    %folder whose names and 'Requested output file:' lines are the naming
    %contract under test.

    properties (Access = private)
        ShadowDir
        TmpDir
    end

    methods (TestMethodSetup)
        function setupShadowsAndTmp(testCase)
            testCase.ShadowDir = tempname;
            mkdir(testCase.ShadowDir);
            %uigetfile shadow: the input file comes from a global
            fid = fopen(fullfile(testCase.ShadowDir, 'uigetfile.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function [filename,pathname,filterindex] = uigetfile(varargin)', ...
                '% TEST SHADOW: answers the open dialog from a global.', ...
                'global PSYRAT_TEST_UIGETFILE', ...
                '[p,n,e] = fileparts(PSYRAT_TEST_UIGETFILE);', ...
                'filename = [n e]; pathname = [p filesep]; filterindex = 1;', ...
                'end');
            fclose(fid);
            %uiputfile shadow: 0 means Cancel, a path means that name/folder
            fid = fopen(fullfile(testCase.ShadowDir, 'uiputfile.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function [filename,pathname,filterindex] = uiputfile(varargin)', ...
                '% TEST SHADOW: answers the save dialog from a global; 0 is Cancel.', ...
                'global PSYRAT_TEST_UIPUTFILE', ...
                'if isempty(PSYRAT_TEST_UIPUTFILE) || isequal(PSYRAT_TEST_UIPUTFILE, 0)', ...
                '    filename = 0; pathname = 0; filterindex = 0; return;', ...
                'end', ...
                '[p,n,e] = fileparts(PSYRAT_TEST_UIPUTFILE);', ...
                'filename = [n e]; pathname = [p filesep]; filterindex = 1;', ...
                'end');
            fclose(fid);
            %tempdir shadow (G94 residual): the error-report writer's last
            %fallback destination is tempdir; a global redirects it to a
            %per-test folder so the fallback tests never touch the real
            %temp folder and an occupied fallback name is deterministic.
            %With the global empty it forwards to the real tempdir, so the
            %other tests see no difference.
            fid = fopen(fullfile(testCase.ShadowDir, 'tempdir.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function d = tempdir', ...
                '% TEST SHADOW: the writer''s fallback folder from a global; forwards when unset.', ...
                'global PSYRAT_TEST_TEMPDIR', ...
                'if ~isempty(PSYRAT_TEST_TEMPDIR)', ...
                '    d = PSYRAT_TEST_TEMPDIR; return;', ...
                'end', ...
                'here = fileparts(mfilename(''fullpath''));', ...
                'rmpath(here);', ...
                'restore = onCleanup(@() addpath(here, ''-begin'')); %#ok<NASGU>', ...
                'd = tempdir;', ...
                'end');
            fclose(fid);
            addpath(testCase.ShadowDir, '-begin');
            rehash;
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);
            testCase.addTeardown(@() teardownShadowsAndTmp(testCase));
            %psyrat_exec re-opens the save dialog on a path with whitespace;
            %the shadow answers the same way every time, so that would spin
            testCase.assumeFalse(any(isspace(testCase.TmpDir)), ...
                'the drive needs a whitespace-free temp path');
        end
    end

    methods (Test)

        function testAnalyzeSurvivesSentinelAfterGroupSubset(testCase)
            %G92: Select Groups was opened under a real Group column (which
            %writes whichgroupscol), then the popup was set back to (none).
            %Analyze must reach the save dialog, not index collist{end}.
            [prefs, data] = localFixture(testCase);
            prefs.proc.inp.whichgroupscol = localColumn(data, 'group');
            prefs.proc.inp.whichgroups = {'A'};
            prefs.proc.inp.group = numel(data.proc.collist);
            localDriveToCancel(testCase, prefs, data);
        end

        function testAnalyzeSurvivesSentinelAfterEventSubset(testCase)
            %G92, the event twin
            [prefs, data] = localFixture(testCase);
            prefs.proc.inp.whicheventscol = localColumn(data, 'event');
            prefs.proc.inp.whichevents = {'e1'};
            prefs.proc.inp.event = numel(data.proc.collist);
            localDriveToCancel(testCase, prefs, data);
        end

        function testAnalyzeSurvivesSentinelAfterOccasionSubset(testCase)
            %G92, the occasion twin
            [prefs, data] = localFixture(testCase);
            prefs.proc.inp.whichtimescol = localColumn(data, 'time');
            prefs.proc.inp.whichtimes = {'t1'};
            prefs.proc.inp.time = numel(data.proc.collist);
            localDriveToCancel(testCase, prefs, data);
        end

        function testAnalyzeWithSentinelAndNoSubsetReachesTheSaveDialog(testCase)
            %control: the sentinel alone never crashed (the first disjunct is
            %guarded), so this passes before and after the fix and proves the
            %drive itself reaches the save dialog
            [prefs, data] = localFixture(testCase);
            prefs.proc.inp.group = numel(data.proc.collist);
            localDriveToCancel(testCase, prefs, data);
        end

        function testAnalyzeSurvivesSentinelSeededByAPriorAnalyze(testCase)
            %G92, the route with no dialog visit: the facet block itself
            %writes whichgroupscol, and the all-measurements-failed relaunch
            %hands the mutated prefs back to the screen. A failed Analyze
            %under a real Group column, then (none), then Analyze.
            global PSYRAT_TEST_UIPUTFILE %#ok<GVMIS>
            [prefs, data] = localFixture(testCase, 'nanrow', 3);
            prefs.proc.inp.group = localColumn(data, 'group');
            PSYRAT_TEST_UIPUTFILE = fullfile(testCase.TmpDir, 'seed.psyrat');
            fig = localOpenSetup(testCase, prefs, data);
            localFireAnalyze(testCase, fig);
            %the relaunched screen carries the prefs psyrat_exec mutated
            fig = localSetupFigure(testCase);
            inplists = localInplists(testCase, fig);
            set(inplists(3), 'Value', numel(data.proc.collist));
            PSYRAT_TEST_UIPUTFILE = 0;
            localFireAnalyze(testCase, fig);
            testCase.verifyTrue(localDialogExists('File Error'), ...
                'Analyze must reach the save dialog after the popup returns to (none)');
        end

        function testMultiMeasurementOutputNamesAppendTheTokenOnce(testCase)
            %G93: uiputfile returns 'testing_splits.psyrat' (macOS supplies
            %the extension); two measurements must save as
            %testing_splits_meas1.psyrat and testing_splits_meas2.psyrat, not
            %testing_splits.psyrat_meas1.psyrat
            reports = localDriveToReports(testCase, 'meas2', true, 'nanrow', 3);
            testCase.assertNumElements(reports, 2, 'one report per failed measurement');
            txt = strjoin(cellfun(@fileread, reports, 'UniformOutput', false), newline);
            testCase.assertTrue(contains(txt, 'Error identifier: varargin:missingcells'), ...
                'the runs must have failed in the preflight, not earlier');
            testCase.verifyTrue(contains(txt, 'Requested output file: testing_splits_meas1.psyrat'), ...
                'the first measurement must be named <stem>_meas1.psyrat');
            testCase.verifyTrue(contains(txt, 'Requested output file: testing_splits_meas2.psyrat'), ...
                'the second measurement must be named <stem>_meas2.psyrat');
            testCase.verifyFalse(contains(txt, '.psyrat_meas'), ...
                'the typed extension must not survive inside the stem');
        end

        function testErrorReportNameDoesNotClobberARawInputOfTheSameStem(testCase)
            %G94 control: the save dialog opens in the input file's folder and
            %.txt is an accepted input format, so the report must never be
            %named <stem>.txt. Passes today (the doubled extension keeps the
            %report out of that namespace) and must keep passing.
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            sentinel = fullfile(savedir, 'testing_splits.txt');
            fid = fopen(sentinel, 'w'); fprintf(fid, 'raw input data, do not touch'); fclose(fid);
            localDriveToReports(testCase, 'nanrow', 3, 'savedir', savedir);
            %the sentinel first: under the naive <stem>.txt fix it IS the report
            testCase.verifyEqual(fileread(sentinel), 'raw input data, do not touch', ...
                'a failed run must not overwrite an input file named <stem>.txt');
            testCase.verifyTrue(isfile(fullfile(savedir, 'testing_splits_error.txt')), ...
                'the report must exist under its own name beside the untouched input');
        end

        function testErrorReportIsNamedFromTheStemWithAnErrorSuffix(testCase)
            %G94: the report is <stem>_error.txt, so it is neither
            %<stem>.psyrat.txt (the tester's l.psyrat.txt) nor <stem>.txt
            reports = localDriveToReports(testCase, 'nanrow', 3);
            [~, names, exts] = cellfun(@fileparts, reports, 'UniformOutput', false);
            testCase.verifyEqual(strcat(names, exts), {'testing_splits_error.txt'});
            reports = localDriveToReports(testCase, 'meas2', true, 'nanrow', 3);
            [~, names, exts] = cellfun(@fileparts, reports, 'UniformOutput', false);
            testCase.verifyEqual(sort(strcat(names, exts)), ...
                {'testing_splits_meas1_error.txt', 'testing_splits_meas2_error.txt'});
        end

        function testAllPreflightErrorsReachTheReport(testCase)
            %G96, GUI half: a file with a missing score AND an Inf score
            %produced two preflight errors but only the first was thrown, so
            %the Inf stayed hidden behind the missing-cells message
            reports = localDriveToReports(testCase, 'nanrow', 3, 'infrow', 5);
            testCase.assertNumElements(reports, 1);
            txt = fileread(reports{1});
            testCase.assertTrue(contains(txt, 'Error identifier: varargin:missingcells'), ...
                'the run must have failed in the preflight, not earlier');
            testCase.verifyTrue(contains(txt, 'missing cells'), 'the missing-cells error must be reported');
            testCase.verifyTrue(contains(txt, 'non-finite'), 'the non-finite error must be reported too');
        end

        %G94 residual (release review 2026-09-20): the _error suffix keeps
        %the report out of the <stem>.txt namespace, but it is not a
        %namespace of its own: the loaded input can itself be named
        %<stem>_error.txt, and an earlier report of that name can already
        %sit in the save folder. The writer must never truncate an existing
        %file at any destination; a taken name gets the first free numeric
        %suffix, a refusing folder sends the report to tempdir, and when every
        %destination fails the diagnostic goes to the console and the dialog
        %says so. Every case drives the REAL Analyze callback; no fit.

        function testReportNeverTruncatesAnInputNamedLikeTheReport(testCase)
            %(a) the input file IS <stem>_error.txt: it survives byte for
            %byte and the report lands at <stem>_error_2.txt
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            %the fixture CSV, copied under the report's own name and loaded
            %from there (the save dialog opens in the input's folder)
            [~, ~, csvfile] = localFixture(testCase, 'nanrow', 3);
            input = fullfile(savedir, 'testing_splits_error.txt');
            copyfile(csvfile, input);
            inputhash = psyrat_hash_file(input);
            localDriveToReports(testCase, 'nanrow', 3, 'savedir', savedir, 'input', input);
            testCase.verifyEqual(psyrat_hash_file(input), inputhash, ...
                'the loaded input named <stem>_error.txt must survive a failed run byte for byte');
            report = fullfile(savedir, 'testing_splits_error_2.txt');
            testCase.assertTrue(isfile(report), 'the report must land at the first free suffix');
            testCase.assertTrue(contains(fileread(report), 'Error identifier: varargin:missingcells'), ...
                'the run must have failed in the preflight, not earlier');
            testCase.verifyEqual(localTxtNames(savedir), ...
                {'testing_splits_error.txt', 'testing_splits_error_2.txt'});
        end

        function testReportSkipsAnOccupiedSuffix(testCase)
            %(b) <stem>_error.txt AND <stem>_error_2.txt are taken: both
            %survive and the report lands at _error_3
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            planted = localPlant(savedir, {'testing_splits_error.txt', 'testing_splits_error_2.txt'});
            localDriveToReports(testCase, 'nanrow', 3, 'savedir', savedir);
            localVerifyUntouched(testCase, planted);
            testCase.verifyEqual(localTxtNames(savedir), ...
                {'testing_splits_error.txt', 'testing_splits_error_2.txt', 'testing_splits_error_3.txt'});
        end

        function testReportSkipsAnOccupiedFallbackName(testCase)
            %(c) the save folder refuses the write AND the fallback name in
            %tempdir is taken: the planted fallback survives and the report
            %lands beside it at _error_2
            global PSYRAT_TEST_TEMPDIR %#ok<GVMIS>
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            PSYRAT_TEST_TEMPDIR = tempname(testCase.TmpDir);
            mkdir(PSYRAT_TEST_TEMPDIR);
            planted = localPlant(PSYRAT_TEST_TEMPDIR, {'testing_splits_error.txt'});
            localRefuseFolder(testCase, savedir);
            localDriveToReports(testCase, 'nanrow', 3, 'savedir', savedir);
            localVerifyUntouched(testCase, planted);
            report = fullfile(PSYRAT_TEST_TEMPDIR, 'testing_splits_error_2.txt');
            testCase.assertTrue(isfile(report), 'the report must take the first free suffix in the fallback folder');
            testCase.verifyTrue(contains(fileread(report), 'Error identifier: varargin:missingcells'), ...
                'the run must have failed in the preflight, not earlier');
            testCase.verifyEmpty(localTxtNames(savedir), 'nothing can be written into the refusing folder');
        end

        function testRepeatFailureKeepsBothReports(testCase)
            %(d) a second failed run into the same folder keeps the first
            %report byte for byte and writes its own at _error_2
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            localDriveToReports(testCase, 'nanrow', 3, 'savedir', savedir);
            first = fullfile(savedir, 'testing_splits_error.txt');
            testCase.assertTrue(isfile(first));
            %the report body carries a one-second timestamp, so a second
            %report written over the first within the same second could be
            %byte-identical; a marker makes the hash decisive
            fid = fopen(first, 'a'); fprintf(fid, '\nmarker: first report'); fclose(fid);
            firsthash = psyrat_hash_file(first);
            localDriveToReports(testCase, 'nanrow', 3, 'savedir', savedir);
            testCase.verifyEqual(psyrat_hash_file(first), firsthash, ...
                'the first report must survive the second failed run');
            testCase.verifyEqual(localTxtNames(savedir), ...
                {'testing_splits_error.txt', 'testing_splits_error_2.txt'});
        end

        function testRefusedFolderPersistsTheReportInTempdir(testCase)
            %(f1) the save folder refuses the write: the report is persisted
            %in tempdir and the console names that path (the free-name
            %control (e) is testErrorReportIsNamedFromTheStemWithAnErrorSuffix)
            global PSYRAT_TEST_TEMPDIR %#ok<GVMIS>
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            PSYRAT_TEST_TEMPDIR = tempname(testCase.TmpDir);
            mkdir(PSYRAT_TEST_TEMPDIR);
            localRefuseFolder(testCase, savedir);
            [~, ~, console] = localDriveToReports(testCase, 'nanrow', 3, 'savedir', savedir);
            report = fullfile(PSYRAT_TEST_TEMPDIR, 'testing_splits_error.txt');
            testCase.assertTrue(isfile(report), 'the report must be persisted in the fallback folder');
            testCase.verifyTrue(contains(fileread(report), 'Error identifier: varargin:missingcells'), ...
                'the run must have failed in the preflight, not earlier');
            testCase.verifyTrue(contains(console, ['Wrote error report: ' report]), ...
                'the console must name the fallback path');
            testCase.verifyEmpty(localTxtNames(savedir), 'nothing can be written into the refusing folder');
            %the closing dialog must not send the user to the refusing folder
            %for a report that is not there; 'refused' is the relocated
            %branch's own word (the unreported branch says 'could not be
            %written'), so a swap of the two branches fails here
            testCase.assertTrue(localDialogExists('All Measurements Failed'));
            testCase.verifyTrue(contains(localDialogText('All Measurements Failed'), 'refused'), ...
                'the dialog must say the save folder refused the report and point at the console');
        end

        function testRefusedFolderWithTheTempFolderNestedUnderIt(testCase)
            %(f1b, review 2026-09-21) the temp folder can lie BENEATH the
            %refusing save folder (a drive root or / chosen as the save
            %folder); the report is still relocated, so the dialog must not
            %infer "in the save folder" from a path prefix
            global PSYRAT_TEST_TEMPDIR %#ok<GVMIS>
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            PSYRAT_TEST_TEMPDIR = fullfile(savedir, 'nested_tmp');
            mkdir(PSYRAT_TEST_TEMPDIR);   %created before the parent refuses
            localRefuseFolder(testCase, savedir);
            [~, ~, console] = localDriveToReports(testCase, 'nanrow', 3, 'savedir', savedir);
            report = fullfile(PSYRAT_TEST_TEMPDIR, 'testing_splits_error.txt');
            testCase.assertTrue(isfile(report), 'the report must be persisted in the nested fallback folder');
            testCase.verifyTrue(contains(console, ['Wrote error report: ' report]), ...
                'the console must name the fallback path');
            testCase.verifyEmpty(localTxtNames(savedir), 'nothing can be written into the refusing folder');
            testCase.assertTrue(localDialogExists('All Measurements Failed'));
            testCase.verifyTrue(contains(localDialogText('All Measurements Failed'), 'refused'), ...
                'a report under a nested temp folder is still relocated; the dialog must say so');
        end

        function testEveryDestinationFailingKeepsTheDiagnosticOnTheConsole(testCase)
            %(f2) a controlled, injected failure of ALL destinations (the
            %allocation helper shadowed to refuse every name): the full
            %diagnostic is printed to the Command Window, the dialog says
            %that no report was written and points at the console, and no
            %pre-existing file changes
            global PSYRAT_TEST_TEMPDIR %#ok<GVMIS>
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            PSYRAT_TEST_TEMPDIR = tempname(testCase.TmpDir);
            mkdir(PSYRAT_TEST_TEMPDIR);
            planted = localPlant(savedir, {'testing_splits_error.txt'});
            localShadowAllocationFailure(testCase);
            [~, ~, console] = localDriveToReports(testCase, 'nanrow', 3, 'savedir', savedir);
            localVerifyUntouched(testCase, planted);
            testCase.verifyEqual(localTxtNames(savedir), {'testing_splits_error.txt'}, ...
                'no report can be written when every destination refuses');
            testCase.verifyEmpty(localTxtNames(PSYRAT_TEST_TEMPDIR), ...
                'no report can be written when every destination refuses');
            for needle = {'Measurement: meas1', 'Requested output file: testing_splits.psyrat', ...
                    'Error identifier: varargin:missingcells', 'Error message:', 'Stack trace:'}
                testCase.verifyTrue(contains(console, needle{1}), ...
                    sprintf('the console must carry the diagnostic line ''%s''', needle{1}));
            end
            testCase.verifyFalse(contains(console, 'Wrote error report:'), ...
                'the console must not claim a report was written');
            testCase.assertTrue(localDialogExists('All Measurements Failed'));
            dialogtxt = localDialogText('All Measurements Failed');
            testCase.verifyTrue(contains(dialogtxt, 'could not be written'), ...
                'the dialog must say no report was written (the unreported branch''s own words)');
            testCase.verifyTrue(contains(dialogtxt, 'Command Window'), ...
                'the dialog must say the diagnostic is on the console');
        end

        %G93 residual (release review 2026-09-20): the batch allocator keyed
        %collisions on makeValidName(lower(name)) and never registered the
        %suffixed name it emitted, so headers that readtable keeps distinct
        %(A, a, a_2) produced the same output name twice and a later
        %measurement's save replaced an earlier one's .psyrat and sidecar.
        %Independently, the save dialog confirmed only the BASE name, so a
        %batch output already in the folder (<name>.psyrat, its
        %_runconfig.json, or a <name>_error.txt from an earlier run) was
        %replaced with no consent at all. The allocator must reserve every
        %emitted name case-insensitively, treat those three derived names in
        %the save folder as occupied, and take the first free numeric suffix
        %(same convention as the error reports); a free name is used as is.
        %The report route (a NaN row fails every measurement in the preflight)
        %shows the allocated names in the reports' 'Requested output file:'
        %lines; the stub route (a warp shadow that saves and writes the
        %sidecar without a fit) shows the files themselves. No CmdStan.

        function testBatchNamesStayDistinctAcrossCaseAndSuffixCollisions(testCase)
            %(a) A, a, a_2: the second name collides with the first
            %case-insensitively and takes _2; the third then collides with
            %THAT emitted name and takes the next free suffix
            [reports, savedir] = localDriveToReports(testCase, 'measnames', {'A', 'a', 'a_2'}, 'nanrow', 3);
            testCase.assertNumElements(reports, 3, 'one report per failed measurement');
            testCase.verifyEqual(localRequestedNames(reports), ...
                sort({'testing_splits_A.psyrat', 'testing_splits_a_2.psyrat', 'testing_splits_a_2_2.psyrat'}), ...
                'every emitted name must be reserved, so a header equal to an emitted suffixed name gets its own suffix');
            testCase.verifyEqual(localTxtNames(savedir), ...
                sort({'testing_splits_A_error.txt', 'testing_splits_a_2_error.txt', 'testing_splits_a_2_2_error.txt'}), ...
                'the reports must follow the allocated names');
        end

        function testBatchNamesPreserveAnOrdinaryNameThatLooksLikeASuffix(testCase)
            %(b) A, a, a_2, a_3: a_3 collides with nothing and keeps its
            %ordinary name; only taken names are suffixed
            reports = localDriveToReports(testCase, 'measnames', {'A', 'a', 'a_2', 'a_3'}, 'nanrow', 3);
            testCase.assertNumElements(reports, 4, 'one report per failed measurement');
            testCase.verifyEqual(localRequestedNames(reports), ...
                sort({'testing_splits_A.psyrat', 'testing_splits_a_2.psyrat', ...
                'testing_splits_a_2_2.psyrat', 'testing_splits_a_3.psyrat'}));
        end

        function testBatchNamesStayDistinctAfterSanitizing(testCase)
            %(c) x_y, x__y, x_y_2: the sanitizer collapses runs of
            %underscores, so the first two share a token; the third equals
            %the suffixed name the second receives
            reports = localDriveToReports(testCase, 'measnames', {'x_y', 'x__y', 'x_y_2'}, 'nanrow', 3);
            testCase.assertNumElements(reports, 3, 'one report per failed measurement');
            testCase.verifyEqual(localRequestedNames(reports), ...
                sort({'testing_splits_x_y.psyrat', 'testing_splits_x_y_2.psyrat', 'testing_splits_x_y_2_2.psyrat'}));
        end

        function testBatchNamesPreserveLongDistinctHeaders(testCase)
            %(d) two 50-character headers that differ only in the last
            %character are distinct names and must both be emitted as they
            %are. A pin, not a defect on this release: the old key was
            %makeValidName, which truncates at namelengthmax; that is 2048 on
            %R2025b Update 7 (probed 2026-09-21) and was 63 on earlier
            %releases (the README's effective floor is R2022a), so the
            %spurious suffix appears only where the key is truncated; the fix
            %keys on the full name everywhere
            long1 = ['m' repmat('a', 1, 48) 'b'];
            long2 = ['m' repmat('a', 1, 48) 'c'];
            reports = localDriveToReports(testCase, 'measnames', {long1, long2}, 'nanrow', 3);
            testCase.assertNumElements(reports, 2, 'one report per failed measurement');
            testCase.verifyEqual(localRequestedNames(reports), ...
                sort({['testing_splits_' long1 '.psyrat'], ['testing_splits_' long2 '.psyrat']}), ...
                'distinct long names must not be suffixed');
        end

        function testSingleMeasurementKeepsTheConsentedBaseName(testCase)
            %(e) control, passes before and after: with one measurement the
            %output IS the name the save dialog confirmed, so an existing
            %file of that name is the user's consented target and the name
            %is never re-allocated (the report route touches nothing)
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            planted = localPlant(savedir, {'testing_splits.psyrat'});
            reports = localDriveToReports(testCase, 'nanrow', 3, 'savedir', savedir);
            testCase.verifyEqual(localRequestedNames(reports), {'testing_splits.psyrat'}, ...
                'a single measurement keeps the base name the dialog confirmed');
            localVerifyUntouched(testCase, planted);
        end

        function testBatchWritesOneOutputAndSidecarPerMeasurement(testCase)
            %(f) the stub route for (a): three distinct .psyrat files and
            %three sidecars must be on disk; on the unchanged code the third
            %save replaced the second
            localShadowEstimationStub(testCase);
            [~, savedir] = localDriveToReports(testCase, 'measnames', {'A', 'a', 'a_2'});
            localAssertBatchCompleted(testCase, 3);
            testCase.verifyEqual(localNames(savedir, '*.psyrat'), ...
                sort({'testing_splits_A.psyrat', 'testing_splits_a_2.psyrat', 'testing_splits_a_2_2.psyrat'}), ...
                'every measurement must have its own .psyrat output');
            testCase.verifyEqual(localNames(savedir, '*_runconfig.json'), ...
                sort({'testing_splits_A_runconfig.json', 'testing_splits_a_2_runconfig.json', ...
                'testing_splits_a_2_2_runconfig.json'}), ...
                'every measurement must have its own sidecar');
        end

        function testExistingBatchOutputIsNeverReplaced(testCase)
            %(g) <name>.psyrat from an earlier batch is in the folder: the
            %dialog confirmed only the base name, so it survives byte for
            %byte and the run takes the first free suffix, sidecar included
            localShadowEstimationStub(testCase);
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            planted = localPlant(savedir, {'testing_splits_meas1.psyrat'});
            localDriveToReports(testCase, 'meas2', true, 'savedir', savedir);
            localAssertBatchCompleted(testCase, 2);
            localVerifyUntouched(testCase, planted);
            testCase.verifyEqual(localNames(savedir, '*.psyrat'), ...
                sort({'testing_splits_meas1.psyrat', 'testing_splits_meas1_2.psyrat', 'testing_splits_meas2.psyrat'}), ...
                'the occupied name must move the run to the first free suffix');
            testCase.verifyEqual(localNames(savedir, '*_runconfig.json'), ...
                sort({'testing_splits_meas1_2_runconfig.json', 'testing_splits_meas2_runconfig.json'}), ...
                'the sidecar must follow the allocated name');
        end

        function testExistingSidecarOccupiesTheName(testCase)
            %(h) only <name>_runconfig.json is in the folder: it belongs to
            %an earlier run, survives, and the name is occupied
            localShadowEstimationStub(testCase);
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            planted = localPlant(savedir, {'testing_splits_meas1_runconfig.json'});
            localDriveToReports(testCase, 'meas2', true, 'savedir', savedir);
            localAssertBatchCompleted(testCase, 2);
            localVerifyUntouched(testCase, planted);
            testCase.verifyEqual(localNames(savedir, '*.psyrat'), ...
                sort({'testing_splits_meas1_2.psyrat', 'testing_splits_meas2.psyrat'}), ...
                'an existing sidecar must occupy the name');
            testCase.verifyEqual(localNames(savedir, '*_runconfig.json'), ...
                sort({'testing_splits_meas1_runconfig.json', 'testing_splits_meas1_2_runconfig.json', ...
                'testing_splits_meas2_runconfig.json'}));
        end

        function testExistingErrorReportOccupiesTheName(testCase)
            %(i) only <name>_error.txt is in the folder (an earlier failed
            %run): a stem is a run identity, so the name is occupied and the
            %successful re-run takes the next suffix (owner ruling 2026-09-21)
            localShadowEstimationStub(testCase);
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            planted = localPlant(savedir, {'testing_splits_meas1_error.txt'});
            localDriveToReports(testCase, 'meas2', true, 'savedir', savedir);
            localAssertBatchCompleted(testCase, 2);
            localVerifyUntouched(testCase, planted);
            testCase.verifyEqual(localNames(savedir, '*.psyrat'), ...
                sort({'testing_splits_meas1_2.psyrat', 'testing_splits_meas2.psyrat'}), ...
                'an existing error report must occupy the name');
        end

        function testExistingOutputOccupiesTheNameCaseInsensitively(testCase)
            %(j) <NAME>.psyrat differs from the batch name only by case: on a
            %case-insensitive filesystem the save would replace it, so the
            %probe is case-insensitive on every platform and the names come
            %out the same everywhere
            localShadowEstimationStub(testCase);
            savedir = tempname(testCase.TmpDir);
            mkdir(savedir);
            planted = localPlant(savedir, {'testing_splits_MEAS1.psyrat'});
            localDriveToReports(testCase, 'meas2', true, 'savedir', savedir);
            localAssertBatchCompleted(testCase, 2);
            localVerifyUntouched(testCase, planted);
            testCase.verifyEqual(localNames(savedir, '*.psyrat'), ...
                sort({'testing_splits_MEAS1.psyrat', 'testing_splits_meas1_2.psyrat', 'testing_splits_meas2.psyrat'}), ...
                'a case variant of the name must occupy it');
        end

        function testBatchNamesReserveDerivedNamesUnderAnotherExtension(testCase)
            %(k, refutation 2026-09-21) the save dialog lets any extension
            %through: saved as study.json, a header ERN_runconfig would name
            %its output exactly like ERN's sidecar, so each member reserves
            %its derived names too, not only its output
            reports = localDriveToReports(testCase, 'measnames', {'meas1', 'meas1_runconfig'}, ...
                'nanrow', 3, 'savename', 'testing_splits.json');
            testCase.assertNumElements(reports, 2, 'one report per failed measurement');
            testCase.verifyEqual(localRequestedNames(reports), ...
                sort({'testing_splits_meas1.json', 'testing_splits_meas1_runconfig_2.json'}), ...
                'a name equal to an earlier member''s sidecar must take a suffix');
        end


        function testRestrictedNoticeFiresForGammaTrtSserrvar(testCase)
            %Queue row 3a (owner rulings 2-4, 2026-09-21): the gamma two-facet
            %designs are RESTRICTED and the user must see it AT SELECTION,
            %before committing to a run. The GUI never sees the analysis
            %number (psyrat_computevarcomp resolves it), so psyrat_exec
            %mirrors the resolver from the settings it already reads and
            %raises the notice BEFORE the save dialog, where Cancel still
            %backs out. Analysis 25: gamma, an occasion column, subject-level
            %on, no dimension. The notice is acknowledged (uiwait), which
            %would park a headless drive, so uiwait is shadowed by a recorder
            %that logs the dialog's name and text and dismisses it.
            [prefs, data] = localFixture(testCase);
            prefs.proc.family = 'gamma';
            prefs.proc.inp.time = localColumn(data, 'time');
            prefs.proc.ssubjrel = 2; prefs.proc.sserrvar = 2;
            rec = localDriveToCancelRecordingNotices(testCase, prefs, data);
            R = psyrat_gamma_twofacet_restriction(25);
            k = find(strcmp(rec.names, R.dialogtitle));
            testCase.verifyNumElements(k, 1, ...
                'gamma 25: exactly one restricted-design notice before the save dialog');
            if ~isempty(k)
                testCase.verifyTrue(contains(rec.texts{k(1)}, 'analysis 25'), ...
                    'the notice must name the resolved design');
                testCase.verifyTrue(contains(rec.texts{k(1)}, 'RESTRICTED'));
                testCase.verifyTrue(contains(rec.texts{k(1)}, 'Cancel at the next'), ...
                    'the notice must tell the user how to back out');
            end
        end

        function testRestrictedNoticeFiresForGammaDynrelTrt(testCase)
            %analysis 14: gamma, an occasion column and a Dimension 1 column
            %(dynamic reliability), subject-level off
            [prefs, data] = localFixture(testCase, 'dimcol', true);
            prefs.proc.family = 'gamma';
            prefs.proc.inp.time = localColumn(data, 'time');
            prefs.proc.inp.dim1 = localColumn(data, 'dim');
            rec = localDriveToCancelRecordingNotices(testCase, prefs, data);
            R = psyrat_gamma_twofacet_restriction(14);
            k = find(strcmp(rec.names, R.dialogtitle));
            testCase.verifyNumElements(k, 1, ...
                'gamma 14: exactly one restricted-design notice before the save dialog');
            if ~isempty(k)
                testCase.verifyTrue(contains(rec.texts{k(1)}, 'analysis 14'), ...
                    'the notice must name the resolved design');
            end
        end

        function testRestrictedNoticeFiresForGammaDynrelSserrvarTrt(testCase)
            %analysis 27: gamma, occasion, Dimension 1 and subject-level on
            [prefs, data] = localFixture(testCase, 'dimcol', true);
            prefs.proc.family = 'gamma';
            prefs.proc.inp.time = localColumn(data, 'time');
            prefs.proc.inp.dim1 = localColumn(data, 'dim');
            prefs.proc.ssubjrel = 2; prefs.proc.sserrvar = 2;
            rec = localDriveToCancelRecordingNotices(testCase, prefs, data);
            R = psyrat_gamma_twofacet_restriction(27);
            k = find(strcmp(rec.names, R.dialogtitle));
            testCase.verifyNumElements(k, 1, ...
                'gamma 27: exactly one restricted-design notice before the save dialog');
            if ~isempty(k)
                testCase.verifyTrue(contains(rec.texts{k(1)}, 'analysis 27'), ...
                    'the notice must name the resolved design');
            end
        end

        function testRestrictedNoticeStaysSilentOutsideTheLabel(testCase)
            %the neighbors of the three designs must NOT raise the notice:
            %the Gaussian arms of 25 and 14 (Supported), gamma 5 (occasion,
            %population level), gamma 6 (subject-level, single occasion),
            %gamma 11 and 26 (one-facet dynamic; outside the label by ruling
            %4). Each drive still reaches the save dialog (asserted by the
            %helper), so a silent case is a real negative, not a drive that
            %stopped early.
            R = psyrat_gamma_twofacet_restriction();
            cases = { ...
                'gaussian 25', 'gaussian', true,  false, 2; ...
                'gaussian 14', 'gaussian', true,  true,  1; ...
                'gamma 5',     'gamma',    true,  false, 1; ...
                'gamma 6',     'gamma',    false, false, 2; ...
                'gamma 11',    'gamma',    false, true,  1; ...
                'gamma 26',    'gamma',    false, true,  2};
            for c = 1:size(cases, 1)
                [prefs, data] = localFixture(testCase, 'dimcol', true);
                prefs.proc.family = cases{c, 2};
                if cases{c, 3}; prefs.proc.inp.time = localColumn(data, 'time'); end
                if cases{c, 4}; prefs.proc.inp.dim1 = localColumn(data, 'dim'); end
                prefs.proc.ssubjrel = cases{c, 5}; prefs.proc.sserrvar = cases{c, 5};
                rec = localDriveToCancelRecordingNotices(testCase, prefs, data);
                testCase.verifyFalse(any(strcmp(rec.names, R.dialogtitle)), sprintf( ...
                    '%s is outside the RESTRICTED label and must raise no notice', cases{c, 1}));
            end
        end

    end

    methods (Access = private)

        function teardownShadowsAndTmp(testCase)
            clear global PSYRAT_TEST_UIGETFILE PSYRAT_TEST_UIPUTFILE PSYRAT_TEST_TEMPDIR
            %the base-class teardown may already have restored the original
            %path, so only rmpath when the shadow dir is still on it
            if any(strcmp(strsplit(path, pathsep), testCase.ShadowDir))
                rmpath(testCase.ShadowDir);
            end
            rehash;
            if isfolder(testCase.ShadowDir); rmdir(testCase.ShadowDir, 's'); end
            if isfolder(testCase.TmpDir); rmdir(testCase.TmpDir, 's'); end
        end

    end
end


function [prefs, data, csvfile] = localFixture(testCase, varargin)
%A balanced single-trial file: 6 ids x 2 groups x 2 events x 2 occasions x 3
%trials, deterministic scores. 'meas2' adds a second measurement column,
%'measnames' (a cellstr) names the measurement columns outright and selects
%them all (the G93 residual tests choose headers that collide), 'nanrow'
%writes NA in every measurement column of that data row, 'infrow' writes
%Inf. prefs carries a fully explicit proc.inp (as after the screen has been
%shown once), so the popups take their Values from it.
opts = struct('meas2', false, 'nanrow', 0, 'infrow', 0, 'measnames', {{}}, ...
    'dimcol', false);
for k = 1:2:numel(varargin); opts.(varargin{k}) = varargin{k+1}; end
if ~isempty(opts.measnames)
    measnames = cellstr(opts.measnames);
elseif opts.meas2
    measnames = {'meas1', 'meas2'};
else
    measnames = {'meas1'};
end
nmeas = numel(measnames);

csvfile = fullfile(testCase.TmpDir, 'drive.csv');
fid = fopen(csvfile, 'w');
%'dimcol' appends a numeric per-participant 'dim' column (a Dimension 1
%candidate for the dynamic-reliability routes; the restricted-design notice
%tests select it)
dimhead = '';
if opts.dimcol; dimhead = ',dim'; end
fprintf(fid, 'id,%s,group,event,time%s\n', strjoin(measnames, ','), dimhead);
row = 0;
for g = 1:2
    for s = 1:6
        for e = 1:2
            for t = 1:2
                for r = 1:3
                    row = row + 1;
                    %the first two columns keep the formulas the earlier
                    %tests were written against; further columns get their
                    %own deterministic series
                    vals = cell(1, nmeas);
                    for m = 1:nmeas
                        switch m
                            case 1; vals{m} = sprintf('%0.3f', mod(0.37 * row, 5) - 2.5);
                            case 2; vals{m} = sprintf('%0.3f', mod(0.53 * row, 4) - 2);
                            otherwise; vals{m} = sprintf('%0.3f', mod((0.29 * m + 0.11) * row, 3) - 1.5);
                        end
                        if row == opts.nanrow; vals{m} = 'NA'; end
                        if row == opts.infrow; vals{m} = 'Inf'; end
                    end
                    dimval = '';
                    if opts.dimcol; dimval = sprintf(',%0.2f', 0.3 * s - 1.05 + 0.5 * (g - 1)); end
                    fprintf(fid, 'g%ds%d,%s,%s,e%d,t%d%s\n', g, s, strjoin(vals, ','), char('A' + g - 1), e, t, dimval);
                end
            end
        end
    end
end
fclose(fid);

global PSYRAT_TEST_UIGETFILE %#ok<GVMIS>
PSYRAT_TEST_UIGETFILE = csvfile;

%collist as psyrat_startproc builds it: the raw headers plus the sentinel.
%Weld the route the report tests rely on: the tokens must import as a
%numeric column with NaN/Inf at the requested rows (readtable's typing is
%release-dependent; a text column would fail earlier, in the loader, and
%the naming tests would stay green against the wrong failure).
raw = psyrat_readtable('file', csvfile);
for m = 1:nmeas
    %the header must survive the import unchanged: the naming tests choose
    %headers readtable keeps as they are (valid names), and a modified
    %header would make the test drive a different collision than intended
    testCase.assertTrue(ismember(measnames{m}, raw.Properties.VariableNames), ...
        sprintf('fixture: readtable must keep the header %s unchanged', measnames{m}));
    col = raw.(measnames{m});
    testCase.assertTrue(isnumeric(col), sprintf('fixture: %s must import as numeric', measnames{m}));
    if opts.nanrow > 0
        testCase.assertTrue(isnan(col(opts.nanrow)), 'fixture: the NA token must import as NaN');
    end
    if opts.infrow > 0
        testCase.assertTrue(isinf(col(opts.infrow)), 'fixture: the Inf token must import as Inf');
    end
end
collist = raw.Properties.VariableNames;
data = struct();
data.raw.colnames = collist;
data.proc.collist = [collist {'(none)'}];
nonecol = numel(data.proc.collist);

prefs = psyrat_defaults;
prefs.ver = psyrat_defineversion;
prefs.guis.fsize = psyrat_guifsize;
inp = struct();
inp.id = localColumn(data, 'id');
inp.meas = cellfun(@(n) localColumn(data, n), measnames);
inp.group = nonecol;
inp.event = nonecol;
inp.time = nonecol;
inp.dim1 = nonecol;
inp.dim2 = nonecol;
inp.weight = nonecol;
inp.whichgroups = '';
inp.whichgroupscol = [];
inp.whichevents = '';
inp.whicheventscol = [];
inp.whichtimes = '';
inp.whichtimescol = [];
prefs.proc.inp = inp;
end


function k = localColumn(data, name)
k = find(strcmp(data.proc.collist, name), 1);
assert(~isempty(k), 'fixture column %s not found', name);
end


function fig = localOpenSetup(testCase, prefs, ~)
%Enter through psyrat_startproc: a stand-in home figure makes it honor the
%passed prefs; the uigetfile shadow supplies the fixture file.
close all force;
figure('Tag', 'psyrat_gui', 'Visible', 'off');
psyrat_startproc('psyrat_prefs', prefs);
fig = localSetupFigure(testCase);
end


function fig = localSetupFigure(testCase)
fig = findall(groot, 'Type', 'figure', 'Name', 'Specify Inputs');
testCase.assertNumElements(fig, 1, 'expected exactly one Specify Inputs screen');
end


function inplists = localInplists(testCase, fig)
btn = localAnalyzeButton(testCase, fig);
cb = btn.Callback;
k = find(strcmp('inplists', cb), 1);
testCase.assertNotEmpty(k, 'the Analyze callback must carry inplists');
inplists = cb{k + 1};
end


function btn = localAnalyzeButton(testCase, fig)
btn = findall(fig, 'Style', 'pushbutton', 'String', 'Analyze');
testCase.assertNumElements(btn, 1, 'Analyze button not found');
end


function localFireAnalyze(testCase, fig)
btn = localAnalyzeButton(testCase, fig);
cb = btn.Callback;
fcn = cb{1};
fcn(btn, [], cb{2:end});
end


function tf = localDialogExists(name)
tf = ~isempty(findall(groot, 'Type', 'figure', 'Name', name));
end


function localDriveToCancel(testCase, prefs, data)
%Fire Analyze with the save dialog answering Cancel: the callback must return
%normally (no escaping error) and leave the cancel dialog behind.
global PSYRAT_TEST_UIPUTFILE %#ok<GVMIS>
PSYRAT_TEST_UIPUTFILE = 0;
fig = localOpenSetup(testCase, prefs, data);
localFireAnalyze(testCase, fig);
testCase.verifyTrue(localDialogExists('File Error'), ...
    'Analyze must run every input check and reach the save dialog');
end


function rec = localDriveToCancelRecordingNotices(testCase, prefs, data)
%Fire Analyze with the save dialog answering Cancel, under a uiwait shadow
%that RECORDS every acknowledged dialog (name and text) and dismisses it.
%The restricted-design notice fires before the save dialog, so the drive
%can stop at Cancel; the 'File Error' dialog proves the drive reached the
%save dialog (every input check passed), which makes a silent case a real
%negative. Returns rec.names and rec.texts (row cells, {} when none).
global PSYRAT_TEST_UIPUTFILE PSYRAT_TEST_UIWAIT_NAMES PSYRAT_TEST_UIWAIT_TEXTS %#ok<GVMIS>
localShadowRecordingUiwait(testCase);
PSYRAT_TEST_UIPUTFILE = 0;
PSYRAT_TEST_UIWAIT_NAMES = {};
PSYRAT_TEST_UIWAIT_TEXTS = {};
fig = localOpenSetup(testCase, prefs, data);
localFireAnalyze(testCase, fig);
testCase.assertTrue(localDialogExists('File Error'), ...
    'Analyze must run every input check and reach the save dialog');
rec = struct('names', {reshape(PSYRAT_TEST_UIWAIT_NAMES, 1, [])}, ...
    'texts', {reshape(PSYRAT_TEST_UIWAIT_TEXTS, 1, [])});
end


function localShadowRecordingUiwait(testCase)
%Shadow uiwait with a recorder: the dialog's Name and its joined text go
%into globals, then the dialog is deleted (the user dismisses it). Removed
%at the end of the test, so the other tests see the real uiwait (or the
%dismiss-only shadow localShadowEstimationStub writes).
wait = fullfile(testCase.ShadowDir, 'uiwait.m');
fid = fopen(wait, 'w');
fprintf(fid, '%s\n', ...
    'function uiwait(h, varargin)', ...
    '% TEST SHADOW: records the acknowledged dialog, then dismisses it.', ...
    'global PSYRAT_TEST_UIWAIT_NAMES PSYRAT_TEST_UIWAIT_TEXTS', ...
    'if nargin >= 1 && all(ishghandle(h))', ...
    '    PSYRAT_TEST_UIWAIT_NAMES{end + 1} = get(h(1), ''Name'');', ...
    '    t = findall(h(1), ''Type'', ''text'');', ...
    '    parts = {};', ...
    '    for i = 1:numel(t); parts = [parts; cellstr(char(t(i).String))]; end', ...
    '    PSYRAT_TEST_UIWAIT_TEXTS{end + 1} = regexprep(strjoin(parts'', '' ''), ''\s+'', '' '');', ...
    '    delete(h);', ...
    'end', ...
    'end');
fclose(fid);
rehash;
testCase.addTeardown(@() localRemoveShadow(wait));
testCase.assertEqual(which('uiwait'), wait, 'the recording uiwait shadow is not first on the path');
end


function [reports, savedir, console] = localDriveToReports(testCase, varargin)
%Fire Analyze with the save dialog answering testing_splits.psyrat in the
%temp folder; the fixture's NaN (and optional Inf) fails each measurement
%inside the loop, so the catch writes one report per measurement. Returns
%the report paths (the .txt files the drive added to the folder; a planted
%.txt, a raw-input sentinel or an earlier report, is not one of them) and
%the captured console output of the drive.
global PSYRAT_TEST_UIPUTFILE PSYRAT_TEST_UIGETFILE %#ok<GVMIS>
%'savedir' names the save folder (a test that plants a raw-input sentinel
%chooses it); otherwise every drive gets a fresh folder, so a test that
%drives twice never reads the earlier drive's reports
k = find(strcmp('savedir', varargin), 1);
if isempty(k)
    savedir = tempname(testCase.TmpDir);
    mkdir(savedir);
else
    savedir = varargin{k + 1};
    varargin(k:k + 1) = [];
end
%'savename' is the name the save dialog answers with (the dialog lets any
%extension through, and one test relies on that)
k = find(strcmp('savename', varargin), 1);
savename = 'testing_splits.psyrat';
if ~isempty(k)
    savename = varargin{k + 1};
    varargin(k:k + 1) = [];
end
%'input' loads the run from that file instead of the fixture CSV (a copy of
%it the test has already made under another name)
k = find(strcmp('input', varargin), 1);
input = '';
if ~isempty(k)
    input = varargin{k + 1};
    varargin(k:k + 1) = [];
end
[prefs, data] = localFixture(testCase, varargin{:});
if ~isempty(input)
    PSYRAT_TEST_UIGETFILE = input;
end
prefs.proc.inp.group = localColumn(data, 'group');
prefs.proc.inp.event = localColumn(data, 'event');
PSYRAT_TEST_UIPUTFILE = fullfile(savedir, savename);
before = localTxtNames(savedir);
fig = localOpenSetup(testCase, prefs, data);
fire = @() localFireAnalyze(testCase, fig);
console = evalc('fire()');
names = reshape(setdiff(localTxtNames(savedir), before), 1, []);
reports = fullfile(savedir, names);
end


function names = localTxtNames(folder)
%sorted names of the .txt files in a folder (a row cell, {} when none)
d = dir(fullfile(folder, '*.txt'));
names = sort({d.name});
if isempty(names); names = {}; end
end


function planted = localPlant(folder, names)
%plant files with distinct content under the given names; returns a struct
%array of path and hash for localVerifyUntouched
planted = struct('path', {}, 'hash', {});
for i = 1:numel(names)
    p = fullfile(folder, names{i});
    fid = fopen(p, 'w');
    fprintf(fid, 'planted file %d: %s, do not touch', i, names{i});
    fclose(fid);
    planted(end + 1) = struct('path', p, 'hash', psyrat_hash_file(p)); %#ok<AGROW>
end
end


function localVerifyUntouched(testCase, planted)
for i = 1:numel(planted)
    testCase.verifyTrue(isfile(planted(i).path), ...
        sprintf('the pre-existing file must still exist: %s', planted(i).path));
    testCase.verifyEqual(psyrat_hash_file(planted(i).path), planted(i).hash, ...
        sprintf('the pre-existing file must be byte-identical: %s', planted(i).path));
end
end


function localRefuseFolder(testCase, folder)
%Make a folder refuse new files (owner write bit off) for the rest of the
%test, restoring it before the temp tree is removed. The refusal is proven
%before it is relied on: a differential assertion needs its positive
%control, and a root-run suite would not be refused at all.
testCase.addTeardown(@() fileattrib(folder, '+w', 'u'));
fileattrib(folder, '-w', 'u');
probe = fullfile(folder, 'refusal_probe.txt');
fid = fopen(probe, 'w');
if fid ~= -1
    fclose(fid);
    delete(probe);
end
testCase.assumeEqual(fid, -1, ...
    'this test needs a folder that refuses new files (not running as root)');
end


function localShadowAllocationFailure(testCase)
%Shadow the allocation helper so EVERY destination refuses: the writer must
%then keep the diagnostic on the console rather than pretend. The shadow is
%removed at the end of the test, so the other tests see the real helper.
shadow = fullfile(testCase.ShadowDir, 'psyrat_fopen_unused.m');
fid = fopen(shadow, 'w');
fprintf(fid, '%s\n', ...
    'function [fid, usedpath] = psyrat_fopen_unused(varargin)', ...
    '% TEST SHADOW: every destination refuses.', ...
    'fid = -1; usedpath = '''';', ...
    'end');
fclose(fid);
rehash;
testCase.addTeardown(@() localRemoveShadow(shadow));
testCase.assertEqual(which('psyrat_fopen_unused'), shadow, ...
    'the allocation shadow is not first on the path');
end


function localRemoveShadow(shadow)
if isfile(shadow); delete(shadow); end
rehash;
end


function txt = localDialogText(name)
%the text of a message dialog as one whitespace-normalized line:
%errordlg/msgbox draw the message as a text object and WRAP long lines
%(about 75 characters), so a phrase can straddle a wrap; joining the
%pieces with spaces lets contains() see the phrase as written
fig = findall(groot, 'Type', 'figure', 'Name', name);
t = findall(fig, 'Type', 'text');
parts = {};
for i = 1:numel(t)
    parts = [parts; cellstr(char(t(i).String))]; %#ok<AGROW>
end
txt = regexprep(strjoin(parts', ' '), '\s+', ' ');
end


function names = localRequestedNames(reports)
%the 'Requested output file:' line of each report, sorted: on the report
%route this is where the allocator's result is visible
names = cell(1, numel(reports));
for i = 1:numel(reports)
    tok = regexp(fileread(reports{i}), 'Requested output file: (\S+)', 'tokens', 'once');
    names{i} = tok{1};
end
names = sort(names);
if isempty(names); names = {}; end
end


function names = localNames(folder, pattern)
%sorted names of the files in a folder that match a pattern (a row cell,
%{} when none)
d = dir(fullfile(folder, pattern));
names = sort({d.name});
if isempty(names); names = {}; end
end


function localAssertBatchCompleted(testCase, nsuccess)
%positive control for the stub route: the batch ran through the stub to
%the Processing Complete box, whose own counts say every measurement
%succeeded (a planted .txt in the folder is not a report, so the folder
%listing is not the evidence here)
testCase.assertTrue(localDialogExists('Processing Complete'), ...
    'the batch must run through the stub to the Processing Complete box');
txt = localDialogText('Processing Complete');
testCase.assertTrue(contains(txt, sprintf('Successful: %d Failed: 0', nsuccess)), ...
    'every measurement must succeed on the stub route');
end


function localShadowEstimationStub(testCase)
%Shadow the estimation warp with a stub that skips CmdStan and does only
%what the warp's tail does after a fit (psyrat_computevarcompwarp.m, after
%the convergence check): attach a rel field, save the .psyrat output, write
%the sidecar. A batch then runs to Processing Complete without a fit and the
%files on disk are the observable. The preflight's warning gate (uiwait on a
%warndlg in psyrat_exec) fires for this fixture, whose 6 participants per
%group and 3 trials per cell are below its thresholds, and would park a
%headless drive, so uiwait is shadowed too: it dismisses the dialog it is
%given and returns, which is what the user does. Both shadows are removed at
%the end of the test, so the other tests see the real functions.
stub = fullfile(testCase.ShadowDir, 'psyrat_computevarcompwarp.m');
fid = fopen(stub, 'w');
fprintf(fid, '%s\n', ...
    'function psyrat_data = psyrat_computevarcompwarp(varargin)', ...
    '% TEST SHADOW: no estimation; saves the output and sidecar as the warp''s tail does.', ...
    'k = find(strcmp(''psyrat_data'', varargin), 1); psyrat_data = varargin{k + 1};', ...
    'k = find(strcmp(''psyrat_prefs'', varargin), 1); psyrat_prefs = varargin{k + 1};', ...
    'psyrat_data.rel = struct(''analysis'', ''stub'');', ...
    'save(fullfile(psyrat_data.proc.savepath, psyrat_data.proc.savename), ''psyrat_data'');', ...
    'try', ...
    '    psyrat_write_runconfig(psyrat_data, psyrat_prefs);', ...
    'catch sidecarME', ...
    '    fprintf(''\nRun-config sidecar could not be written: %s\n'', sidecarME.message);', ...
    'end', ...
    'end');
fclose(fid);
wait = fullfile(testCase.ShadowDir, 'uiwait.m');
fid = fopen(wait, 'w');
fprintf(fid, '%s\n', ...
    'function uiwait(h, varargin)', ...
    '% TEST SHADOW: the user dismisses the dialog; nothing blocks.', ...
    'if nargin >= 1 && all(ishghandle(h)); delete(h); end', ...
    'end');
fclose(fid);
rehash;
testCase.addTeardown(@() localRemoveShadow(stub));
testCase.addTeardown(@() localRemoveShadow(wait));
testCase.assertEqual(which('psyrat_computevarcompwarp'), stub, ...
    'the estimation stub is not first on the path');
testCase.assertEqual(which('uiwait'), wait, 'the uiwait shadow is not first on the path');
end
