classdef TestCheckOutputNames < PsyRATTestBase
    %Row 31 (LA-2): psyrat_check_output_names on its own. These are the
    %safeguards root required of its creation probe, each observed directly
    %rather than through a drive (the drives are TestOutputNamePreflight):
    %  - a refusal names the request, the failing file and the remedy, and
    %    only when this volume really cannot hold the name (measured here,
    %    never a fixed cap);
    %  - a folder at the output or sidecar path is refused and kept intact; a
    %    folder at the report path is no new restriction;
    %  - existing files are neither probed as targets nor changed;
    %  - a folder that refuses every new entry is not called a name problem;
    %  - a name with a path separator, '.' or '..' touches nothing anywhere;
    %  - no probe folder survives a check, and a probe folder that cannot be
    %    removed is reported with its exact path, never deleted recursively,
    %    and never claimed as removed;
    %  - the global random number generator is not touched;
    %  - probe names have exactly the component's length, keep the user's
    %    characters verbatim, and are never the target itself; a short ASCII
    %    extension borrows trailing ASCII characters for the 7-character tail,
    %    and where no such tail exists (a non-ASCII character before .mat)
    %    the output is not probed and the longer sidecar carries the check;
    %  - the refusal and leftover messages say only what is known (review L1,
    %    L2);
    %  - a refused test FOLDER is not proof (phase 5, root adjudication H1:
    %    MATLAB's mkdir refused a 210-unit multibyte name with
    %    MATLAB:FILE:PATHMAX at a depth where save and fopen wrote it). The
    %    refusal is confirmed with an equal-length test FILE created
    %    exclusively (java.io.File.createNewFile); only if that fails too is the
    %    name refused. A name someone else takes in the meantime is never
    %    claimed, a test file holding foreign data is never deleted, a test
    %    file that cannot be removed is reported with its exact path, and
    %    without a Java runtime the name is left unrefused with the reason.
    %    The control is a short test FILE made the same way, so a file test
    %    that cannot reach the save folder refuses nothing (counterfactual:
    %    a short test FOLDER would have turned that into a refusal).
    %  - review coverage (phase 7): a test file Java makes in another, owned
    %    folder is removed and decides nothing; a relative save folder is
    %    tested in that folder, not in Java's working folder; valid three- and
    %    four-byte names pass the file test. Multibyte fixtures that the ASCII
    %    limit cannot size on some volume are measured by the writers instead
    %    (OutputNameFixture.largestValid), never assumed.
    %
    %Tests that need to see inside the check shadow MATLAB functions for the
    %one call under test: mkdir (recordProbes: recording, optionally dropping
    %a file into each new test folder; runWithSeams: recording, or refusing
    %every or only the first test folder without creating anything), dir
    %(runWithSeams: pointing the Java route at a missing folder or at an
    %owned sibling of the save folder), isfile
    %(runWithSeams: planting a
    %name just seen free, writing foreign data into a test file, or locking
    %the folder, each simulating another process) and usejava (runWithSeams:
    %no Java runtime). Each shadow forwards to the real function by taking
    %itself off the path for the call, the same forwarding pattern as the
    %tempdir shadow of the Analyze drive tests. runWithSeams also installs a
    %warning recorder (phase 6), which forwards through builtin instead,
    %because re-joining the path raises a name-conflict warning that
    %overwrites the last warning; it records the check's leftover warning as
    %MATLAB reports it when it is issued. The deep-Unicode test uses no seam:
    %it establishes the real counterexample first.

    properties (Access = private)
        Lim
        OracleDir
        Root
        Folder
        ShadowDir
    end

    methods (TestClassSetup)
        function measureTheVolume(testCase)
            folder = tempname;
            mkdir(folder);
            testCase.addTeardown(@() rmdir(folder, 's'));
            testCase.OracleDir = folder;
            testCase.Lim = OutputNameFixture.measureLimits(folder);
        end
    end

    methods (TestMethodSetup)
        function freshFolder(testCase)
            %a task-owned root holding the save folder, so a '..' request
            %resolves into a folder the test owns and can inspect
            root = tempname;
            mkdir(root);
            testCase.addTeardown(@() rmdir(root, 's'));
            testCase.Root = char(java.io.File(root).getCanonicalPath());
            testCase.Folder = fullfile(testCase.Root, 'save');
            mkdir(testCase.Folder);
            %the shadow lives OUTSIDE the inspected tree, so writing it is
            %never mistaken for a change the check made
            shadow = tempname;
            mkdir(shadow);
            testCase.addTeardown(@() rmdir(shadow, 's'));
            testCase.ShadowDir = char(java.io.File(shadow).getCanonicalPath());
        end
    end

    methods (Test)

        function testWritableNameIsAcceptedAndLeavesNothing(testCase)
            before = OutputNameFixture.snapshot(testCase.Root);
            [problems, leftovers] = psyrat_check_output_names(testCase.Folder, 'plain.psyrat');
            testCase.verifyEmpty(problems, 'an ordinary name must pass');
            testCase.verifyEmpty(leftovers, 'no probe folder may be left');
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before, ...
                'the check must leave the folder exactly as it found it');
        end

        function testOutputPastTheLimitIsNamedAndLeavesNothing(testCase)
            savename = [repmat('o', 1, testCase.Lim.save - 6) '.psyrat'];
            [writable, detail] = OutputNameFixture.oracle(testCase.OracleDir, savename);
            before = OutputNameFixture.snapshot(testCase.Root);
            [problems, leftovers] = psyrat_check_output_names(testCase.Folder, savename);
            testCase.verifyEmpty(leftovers, 'no probe folder may be left');
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before, ...
                'a refusal must leave the folder exactly as it found it');
            if writable
                testCase.log(1, 'no name limit here: the long output is writable and must pass');
                testCase.verifyEmpty(problems);
                return;
            end
            testCase.assertFalse(detail.output, 'premise: the output is past the measured limit');
            testCase.assertNumElements(problems, 1, 'exactly one problem, the first failing file');
            testCase.verifyEqual(problems.id, 'psyrat_data:savename');
            testCase.verifyEqual(problems.component, savename, 'the failing file is the output itself');
            testCase.verifyTrue(contains(problems.what, savename), 'the request must be named');
            testCase.verifyTrue(contains(problems.what, testCase.Folder), 'the folder must be named');
            testCase.verifyNotEmpty(problems.howto, 'a remedy must be given');
        end

        function testSidecarPastTheLimitNamesTheSidecar(testCase)
            stem = repmat('s', 1, testCase.Lim.sidecar - numel('_runconfig.json') + 1);
            savename = [stem '.psyrat'];
            [writable, detail] = OutputNameFixture.oracle(testCase.OracleDir, savename);
            [problems, leftovers] = psyrat_check_output_names(testCase.Folder, savename);
            testCase.verifyEmpty(leftovers);
            if writable
                testCase.verifyEmpty(problems);
                return;
            end
            testCase.assertTrue(detail.output, 'premise: the output alone fits');
            testCase.assertNumElements(problems, 1);
            testCase.verifyEqual(problems.component, [stem '_runconfig.json'], ...
                'the failing file is the sidecar, not the output');
            testCase.verifyTrue(contains(problems.what, [stem '_runconfig.json']));
            testCase.verifyTrue(contains(problems.what, savename));
            %wording (review L1): the output itself fits, so the message states
            %the refusal and must not claim the output cannot be saved
            testCase.verifyTrue(contains(problems.what, 'refused'), 'the message must state the refusal');
            testCase.verifyFalse(contains(problems.what, 'cannot be saved'), ...
                'the output itself fits, so it must not be called unsavable');
        end

        function testShortAsciiExtensionProbesAnEqualLengthName(testCase)
            %Review T3 (1). A .json output (5 characters): the replaced region
            %extends over trailing ASCII characters of the user's part until
            %it holds 7, the probe keeps the rest of the name verbatim, and a
            %name one past the measured limit is refused as the output.
            savename = [repmat('a', 1, testCase.Lim.save - 4) '.json'];   %one past the limit
            keepPart = savename(1:end - 7);
            [writable, detail] = OutputNameFixture.oracle(testCase.OracleDir, savename);
            before = OutputNameFixture.snapshot(testCase.Root);
            [names, problems, leftovers] = testCase.recordProbes( ...
                @() psyrat_check_output_names(testCase.Folder, savename));
            testCase.verifyEmpty(leftovers, 'no probe folder may be left');
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before);
            testCase.assertNotEmpty(names, 'premise: the output was probed');
            first = names{1};
            testCase.verifyEqual(numel(first), numel(savename), 'the output probe has the output''s length');
            testCase.verifyTrue(startsWith(first, keepPart), ...
                'the probe keeps the user''s part before the 7-character tail');
            testCase.verifyFalse(strcmpi(first, savename), 'a probe is never the target');
            if writable
                testCase.verifyEmpty(problems, 'the volume holds the name, so nothing is refused');
            else
                testCase.assertFalse(detail.output, 'premise (oracle): the output name itself cannot be written');
                testCase.assertNumElements(problems, 1);
                testCase.verifyEqual(problems.component, savename, 'the output is the first failing file');
            end
        end

        function testNonAsciiBeforeAShortExtensionIsCoveredByTheSidecar(testCase)
            %Review T3 (2). A stem ending in e-acute with .mat: no ASCII tail of
            %7 exists before the extension, so the output is not probed. The
            %sidecar, always longer, is probed, so a name one past the limit
            %is still refused, with the sidecar as the failing file (the
            %dominance the helper relies on).
            stem = [repmat('a', 1, testCase.Lim.save - 4) char(233)];
            savename = [stem '.mat'];   %one past the limit
            [writable, detail] = OutputNameFixture.oracle(testCase.OracleDir, savename);
            before = OutputNameFixture.snapshot(testCase.Root);
            [names, problems, leftovers] = testCase.recordProbes( ...
                @() psyrat_check_output_names(testCase.Folder, savename));
            testCase.verifyEmpty(leftovers, 'no probe folder may be left');
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before);
            testCase.verifyFalse(any(cellfun(@numel, names) == numel(savename)), ...
                'the output cannot be probed without replacing a user character, so it is not');
            testCase.assertNotEmpty(names, 'premise: the sidecar was probed');
            testCase.verifyEqual(numel(names{1}), numel([stem '_runconfig.json']), ...
                'the first probe is the sidecar''s');
            testCase.verifyTrue(startsWith(names{1}, stem), 'the user''s characters are kept verbatim');
            if writable
                testCase.verifyEmpty(problems, 'the volume holds the names, so nothing is refused');
            else
                testCase.assertFalse(detail.output, 'premise (oracle): the output itself is past the limit');
                testCase.assertFalse(detail.sidecar, 'premise (oracle): so is the longer sidecar');
                testCase.assertNumElements(problems, 1, 'the request must still be refused');
                testCase.verifyEqual(problems.component, [stem '_runconfig.json'], ...
                    'the sidecar is the failing file that covers the unprobed output');
            end
        end

        function testFolderAtTheOutputOrSidecarPathIsRefusedAndKept(testCase)
            cases = {'occ.psyrat', 'occ.psyrat'; 'occ_runconfig.json', 'occ.psyrat'};
            for k = 1:size(cases, 1)
                blocker = fullfile(testCase.Folder, cases{k, 1});
                mkdir(blocker);
                fid = fopen(fullfile(blocker, 'inside.txt'), 'w');
                fprintf(fid, 'content of a folder the check must not touch');
                fclose(fid);
                before = OutputNameFixture.snapshot(testCase.Root);
                [problems, leftovers] = psyrat_check_output_names(testCase.Folder, cases{k, 2});
                testCase.verifyNumElements(problems, 1, sprintf('a folder at %s must be refused', cases{k, 1}));
                if isscalar(problems)
                    testCase.verifyEqual(problems.component, cases{k, 1});
                    testCase.verifyTrue(contains(problems.what, cases{k, 1}));
                end
                testCase.verifyEmpty(leftovers);
                testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before, ...
                    'the folder and its content must be left exactly as they were');
                rmdir(blocker, 's');   %the test's own folder, removed for the next case
            end
        end

        function testFolderAtTheReportPathIsNoNewRestriction(testCase)
            %the report writer already moves on to a free or timestamped name
            blocker = fullfile(testCase.Folder, 'rep_error.txt');
            mkdir(blocker);
            before = OutputNameFixture.snapshot(testCase.Root);
            [problems, leftovers] = psyrat_check_output_names(testCase.Folder, 'rep.psyrat');
            testCase.verifyEmpty(problems, 'a folder at the report path must not refuse the run');
            testCase.verifyEmpty(leftovers);
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before);
        end

        function testExistingFilesAreNeitherProbedNorChanged(testCase)
            %a zero-byte consented target and an earlier sidecar: their names
            %exist, so they are representable; nothing may open or change them
            fid = fopen(fullfile(testCase.Folder, 'keep.psyrat'), 'w'); fclose(fid);
            fid = fopen(fullfile(testCase.Folder, 'keep_runconfig.json'), 'w');
            fprintf(fid, 'an earlier sidecar'); fclose(fid);
            before = OutputNameFixture.snapshot(testCase.Root);
            names = testCase.recordProbes(@() psyrat_check_output_names(testCase.Folder, 'keep.psyrat'));
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before, ...
                'existing files must keep their bytes, and nothing may be left');
            testCase.verifyFalse(any(strcmpi(names, 'keep.psyrat')) || ...
                any(strcmpi(names, 'keep_runconfig.json')), ...
                'an existing target must never be used as a probe name');
            testCase.verifyNumElements(names, 1, ...
                'only the report name, which does not exist yet, needs a probe');
        end

        function testRefusingFolderIsNotANameProblem(testCase)
            %A save folder that refuses every new entry is outside this
            %check: calling it a name problem would be false, so nothing is
            %refused (the existing late path and its report handle it).
            %Delivery control first: the folder must really refuse. A run
            %with root privileges cannot establish that, and fails here
            %explicitly rather than skipping.
            testCase.addTeardown(@() fileattrib(testCase.Folder, '+w', 'u'));
            fileattrib(testCase.Folder, '-w', 'u');
            probe = fullfile(testCase.Folder, 'refusal_probe.txt');
            fid = fopen(probe, 'w');
            if fid >= 3; fclose(fid); delete(probe); end
            testCase.assertEqual(fid, -1, ...
                'premise: the save folder must refuse new files (this cannot be arranged when running as root)');
            long = [repmat('o', 1, testCase.Lim.save - 6) '.psyrat'];
            for savename = {'short.psyrat', long}
                [problems, leftovers] = psyrat_check_output_names(testCase.Folder, savename{1});
                testCase.verifyEmpty(problems, sprintf( ...
                    'a refusing folder is not a name problem (%d characters)', numel(savename{1})));
                testCase.verifyEmpty(leftovers);
            end
        end

        function testSeparatorsAndDotNamesTouchNothing(testCase)
            mkdir(fullfile(testCase.Folder, 'sub'));
            before = OutputNameFixture.snapshot(testCase.Root);
            names = {'sub/x.psyrat', ['sub' '\' 'x.psyrat'], '../x.psyrat', '.', '..'};
            for k = 1:numel(names)
                [problems, leftovers] = psyrat_check_output_names(testCase.Folder, names{k});
                testCase.verifyEmpty(problems, sprintf('%s: the existing contract is kept', names{k}));
                testCase.verifyEmpty(leftovers);
            end
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before, ...
                'no folder, inside or outside the save folder, may be touched');
        end

        function testRandomNumberGeneratorIsUntouched(testCase)
            %fit seeds must not depend on whether a name was checked
            long = [repmat('o', 1, testCase.Lim.save - 6) '.psyrat'];
            state = rng;
            psyrat_check_output_names(testCase.Folder, 'plain.psyrat');
            psyrat_check_output_names(testCase.Folder, long);
            psyrat_check_output_names(testCase.Folder, [repmat(char(233), 1, 40) '.psyrat']);
            testCase.verifyEqual(rng, state, 'the check must not draw from the global stream');
            %positive control: one draw is visible in the same comparison
            rand;
            testCase.verifyNotEqual(rng, state, 'control: a single draw must change the compared state');
        end

        function testProbeNamesAreEqualLengthDisjointAndKeepTheUserStem(testCase)
            stem = [repmat(char(233), 1, 10) '_x'];
            savename = [stem '.psyrat'];
            components = {savename, [stem '_runconfig.json'], [stem '_error.txt']};
            names = testCase.recordProbes(@() psyrat_check_output_names(testCase.Folder, savename));
            testCase.assertNumElements(names, 3, 'one probe per component that does not exist yet');
            for k = 1:3
                p = names{k};
                c = components{k};
                testCase.verifyFalse(strcmpi(p, c), 'a probe must never be the target itself');
                testCase.verifyEqual(numel(p), numel(c), 'same length in UTF-16 units');
                testCase.verifyEqual(numel(unicode2native(p, 'UTF-8')), numel(unicode2native(c, 'UTF-8')), ...
                    'same length in UTF-8 bytes');
                testCase.verifyTrue(startsWith(p, stem), 'the user''s characters must be kept verbatim');
                testCase.verifyTrue(all(p(numel(stem) + 1:end) < 128), 'only an ASCII tail is replaced');
            end
            d = dir(testCase.Folder);
            testCase.verifyEmpty(setdiff({d.name}, {'.', '..'}), 'every probe must be gone');
        end

        function testCleanupFailureIsReportedWithTheExactPathAndNothingIsDeleted(testCase)
            %Another process drops a file into each new probe folder, so the
            %non-recursive removal fails. The check must say so with the exact
            %path, keep that file, and still return its answer about the name.
            lastwarn('');
            [~, problems, leftovers] = testCase.recordProbes( ...
                @() psyrat_check_output_names(testCase.Folder, 'busy.psyrat'), true);
            [msg, id] = lastwarn;
            testCase.verifyEmpty(problems, 'the name itself is creatable');
            testCase.assertNotEmpty(leftovers, 'the failed removals must be returned');
            testCase.verifyEqual(id, 'psyrat:outputNameProbeLeft', 'the failure must be a visible warning');
            testCase.verifyTrue(contains(msg, leftovers{end}), 'the warning must give the exact path');
            testCase.verifyFalse(contains(lower(msg), 'removed the'), 'removal must not be claimed');
            %wording (review L2): the cause is unknown, so the advice is to
            %look before deleting, never to delete unseen content
            testCase.verifyTrue(contains(msg, 'check its contents before deleting it'), ...
                'the warning must ask for a look before any deletion');
            testCase.verifyFalse(contains(msg, 'can be deleted'), ...
                'the warning must not advise deleting unknown content');
            for k = 1:numel(leftovers)
                testCase.verifyTrue(isfolder(leftovers{k}), 'a reported leftover must really exist');
                testCase.verifyTrue(isfile(fullfile(leftovers{k}, 'intruder.txt')), ...
                    'content the check did not create must never be deleted');
                testCase.verifyTrue(startsWith(leftovers{k}, [testCase.Folder filesep]), ...
                    'leftovers can only be in the save folder');
            end
        end

        % ---- a refused test folder is confirmed with a test file (phase 5) ---

        function testDeepUnicodeFolderRefusalIsOverruledByTheWritersFile(testCase)
            %The observed false refusal (T4; root adjudication H1). At this
            %depth MATLAB's mkdir refuses an equal-length folder for a long
            %multibyte name (MATLAB:FILE:PATHMAX), while save and fopen write
            %the files. The check must accept the name, because a file of it
            %can be created, and leave nothing behind.
            %Portability (review F6, F7): the stem that puts the sidecar at the
            %ASCII-measured limit is kept wherever the writers confirm it, as
            %on the measured APFS volume. Where they do not, the longest stem
            %they write at this depth is used instead; no counting rule is
            %assumed for e-acute. Where MATLAB's mkdir makes the equal-length
            %folder, the counterexample is absent: that is logged, and the
            %product property is still checked (the forced-fallback seam test
            %covers the file test on any platform).
            if testCase.limitsFound()
                k = testCase.Lim.sidecar - numel('_runconfig.json') - 2;
            else
                k = 238;
            end
            makeName = @(n) [repmat(char(233), 1, n) '_d.psyrat'];
            %independent evidence at exactly the save folder's depth: a sibling
            %whose name has the same length as 'save'
            twin = fullfile(testCase.Root, 'twin');
            mkdir(twin);
            [writable, detail] = OutputNameFixture.oracle(twin, makeName(k));
            branch = 'ASCII-limit stem';
            if ~writable
                branch = 'stem measured by the writers';
                k = OutputNameFixture.largestValid(twin, makeName, 1, k);
                [writable, detail] = OutputNameFixture.oracle(twin, makeName(max(k, 1)));
            end
            stem = [repmat(char(233), 1, k) '_d'];
            savename = [stem '.psyrat'];
            probeLike = [stem '.abcdef'];            %the output's length
            [st, ~, id] = mkdir(twin, probeLike);
            folderMade = st && isempty(id);
            if folderMade; rmdir(fullfile(twin, probeLike)); end
            rmdir(twin);
            testCase.log(1, sprintf(['deep-Unicode fixture (%s): %d e-acute; output %d UTF-16 ' ...
                'units, %d UTF-8 bytes; equal-length entry full path %d bytes; writers %d ' ...
                '(output %d, sidecar %d, report %d); MATLAB mkdir made the equal-length ' ...
                'folder: %d'], branch, k, numel(savename), ...
                numel(unicode2native(savename, 'UTF-8')), ...
                numel(unicode2native(fullfile(testCase.Folder, probeLike), 'UTF-8')), ...
                writable, detail.output, detail.sidecar, detail.report, folderMade));
            testCase.assertTrue(k >= 1 && writable, sprintf(['premise (oracle at the save ' ...
                'folder''s depth): the writers create every file of this name (output %d, ' ...
                'sidecar %d, report %d)'], detail.output, detail.sidecar, detail.report));
            if folderMade
                testCase.log(1, ['The counterexample is absent here: MATLAB''s mkdir made the ' ...
                    'equal-length folder, so the ordinary folder test decides. The product ' ...
                    'property is still checked below.']);
            end
            before = OutputNameFixture.snapshot(testCase.Root);
            [problems, leftovers, notes] = psyrat_check_output_names(testCase.Folder, savename);
            testCase.verifyEmpty(problems, 'a name the writers can create must not be refused');
            testCase.verifyEmpty(leftovers, 'no test entry may be left');
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before, ...
                'the check must leave the folder exactly as it found it');
            testCase.verifyEqual(unique({notes.outcome}), {'created'}, ...
                'every component must rest on a test entry the check created');
            if ~folderMade
                testCase.verifyTrue(any(strcmp({notes.method}, 'file') & strcmp({notes.outcome}, 'created')), ...
                    'the decision must rest on a test file created after the test folder was refused');
            end
        end

        function testRefusedTestFolderFallsBackToAnExclusiveTestFile(testCase)
            %A seam refuses every test folder, so the file test decides on any
            %platform: an ordinary name is accepted on the evidence of one test
            %file per component, and nothing is left.
            before = OutputNameFixture.snapshot(testCase.Root);
            [seam, problems, leftovers, notes] = testCase.runWithSeams( ...
                @() psyrat_check_output_names(testCase.Folder, 'plain.psyrat'), ...
                struct('mkdirMode', 'fail'));
            testCase.assertNumElements(seam.mkdirNames, 3, 'premise: one refused test folder per component');
            testCase.verifyEmpty(problems, 'the writers can create these files');
            testCase.verifyEmpty(leftovers, 'every test file must be removed');
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before);
            testCase.verifyEqual({notes.method}, {'file', 'file', 'file'});
            testCase.verifyEqual({notes.outcome}, {'created', 'created', 'created'});
        end

        function testRefusedTestFolderAndRefusedTestFileStillRefuseTheName(testCase)
            %Failed creation: with every test folder refused, an output one past
            %the measured limit is still refused, because the test FILE cannot
            %be created either, while the short control test file can.
            savename = [repmat('o', 1, testCase.Lim.save - 6) '.psyrat'];
            [writable, detail] = OutputNameFixture.oracle(testCase.OracleDir, savename);
            before = OutputNameFixture.snapshot(testCase.Root);
            [~, problems, leftovers, notes] = testCase.runWithSeams( ...
                @() psyrat_check_output_names(testCase.Folder, savename), struct('mkdirMode', 'fail'));
            testCase.verifyEmpty(leftovers);
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before);
            if writable
                testCase.verifyEmpty(problems, 'no limit on this volume: nothing to refuse');
                return;
            end
            testCase.assertFalse(detail.output, 'premise (oracle): the output name cannot be written');
            testCase.assertNumElements(problems, 1, 'the name must still be refused before the fit');
            testCase.verifyEqual(problems.component, savename);
            testCase.verifyEqual(notes(1).method, 'file', 'the refusal rests on the test file');
            testCase.verifyEqual(notes(1).outcome, 'refused');
        end

        function testATestNameTakenInTheMeantimeIsNeverClaimed(testCase)
            %Collision: "another process" creates each test name right after
            %the check saw it free and before the file test. The atomic create
            %must answer that it exists; the entry is never treated as ours,
            %each retry uses a new tail of the same length (3 per component),
            %and every planted entry survives byte for byte. With no evidence
            %either way, nothing is refused.
            components = {'plain.psyrat', 'plain_runconfig.json', 'plain_error.txt'};
            [seam, problems, leftovers, notes] = testCase.runWithSeams( ...
                @() psyrat_check_output_names(testCase.Folder, 'plain.psyrat'), ...
                struct('mkdirMode', 'fail', 'isfileMode', 'plant', 'skip', {components}));
            testCase.verifyEmpty(problems, 'no evidence either way: nothing may be refused');
            testCase.verifyEmpty(leftovers, 'a planted entry is not the check''s to report or remove');
            testCase.verifyEqual({notes.outcome}, {'indeterminate', 'indeterminate', 'indeterminate'});
            testCase.assertNumElements(seam.acted, 9, 'premise: three planted names per component');
            lens = cellfun(@(p) numel(localName(p)), seam.acted);
            testCase.verifyEqual(sort(lens), sort(repelem(cellfun(@numel, components), 3)), ...
                'every retry keeps the component''s length (no suffix growth)');
            for k = 1:numel(seam.acted)
                testCase.verifyEqual(fileread(seam.acted{k}), 'planted by the test', ...
                    sprintf('a planted entry must survive byte for byte: %s', seam.acted{k}));
            end
        end

        function testATestFileSomeoneElseWroteIntoIsNotDeleted(testCase)
            %Refused delete: "another process" writes into each test file after
            %the check created it. Content the check did not write is never
            %deleted; each such file is reported with its exact path, and the
            %name itself is accepted (a file of it could be created).
            components = {'busy.psyrat', 'busy_runconfig.json', 'busy_error.txt'};
            [seam, problems, leftovers] = testCase.runWithSeams( ...
                @() psyrat_check_output_names(testCase.Folder, 'busy.psyrat'), ...
                struct('mkdirMode', 'fail', 'isfileMode', 'intrude', 'skip', {components}));
            testCase.verifyEmpty(problems, 'the name itself is creatable');
            testCase.assertNotEmpty(seam.acted, 'premise: foreign data was written into a test file');
            testCase.verifyEqual(sort(leftovers), sort(seam.acted), ...
                'exactly the files holding foreign data are left, and reported');
            testCase.verifyEqual(seam.lastWarningId, 'psyrat:outputNameProbeLeft', ...
                'the leftover must be a visible warning');
            testCase.verifyTrue(contains(seam.lastWarning, leftovers{end}), 'the warning must give the exact path');
            testCase.verifyTrue(contains(seam.lastWarning, 'check its contents before deleting it'));
            for k = 1:numel(leftovers)
                testCase.verifyEqual(fileread(leftovers{k}), 'foreign data', ...
                    'the foreign content must be kept exactly');
            end
        end

        function testATestFileThatCannotBeRemovedIsReported(testCase)
            %Refused removal: the save folder stops accepting changes right
            %after the check created its first test file, so removing it fails.
            %The file is reported with its exact path and never claimed as
            %removed. The locked folder then refuses everything, which is not a
            %name problem, so nothing is refused. Needs a non-root user: root
            %could remove the file anyway, and the leftover premise then fails
            %explicitly rather than skipping.
            testCase.addTeardown(@() fileattrib(testCase.Folder, '+w', 'u'));
            components = {'lock.psyrat', 'lock_runconfig.json', 'lock_error.txt'};
            [seam, problems, leftovers, notes] = testCase.runWithSeams( ...
                @() psyrat_check_output_names(testCase.Folder, 'lock.psyrat'), ...
                struct('mkdirMode', 'fail', 'isfileMode', 'lock', 'skip', {components}));
            testCase.assertNumElements(seam.acted, 1, 'premise: the folder was locked after the first test file');
            testCase.verifyEmpty(problems, 'a folder that refuses everything is not a name problem');
            testCase.assertEqual(leftovers, seam.acted, ...
                'the test file that could not be removed must be reported');
            testCase.verifyTrue(isfile(leftovers{1}), 'the reported file must really exist');
            testCase.verifyEqual(seam.lastWarningId, 'psyrat:outputNameProbeLeft');
            testCase.verifyTrue(contains(seam.lastWarning, leftovers{1}), 'the warning must give the exact path');
            testCase.verifyEqual(notes(1).outcome, 'created', 'the name itself could be created');
        end

        function testAFileTestThatMissesTheSaveFolderRefusesNothing(testCase)
            %Counterfactual for the control: only the first test folder is
            %refused, and the Java route is pointed at a missing folder, so the
            %test file fails for a reason unrelated to the name. A short test
            %FOLDER would still be accepted and turn that failure into a
            %refusal of an ordinary name; the short control FILE, made the same
            %way, fails too, so nothing is refused and nothing is created.
            before = OutputNameFixture.snapshot(testCase.Root);
            [seam, problems, leftovers, notes] = testCase.runWithSeams( ...
                @() psyrat_check_output_names(testCase.Folder, 'plain.psyrat'), ...
                struct('mkdirMode', 'failFirst', 'dirMode', 'elsewhere'));
            testCase.assertGreaterThanOrEqual(seam.redirected, 1, ...
                'premise: the Java route was pointed at a missing folder');
            testCase.verifyEmpty(problems, 'a file test that cannot reach the folder must never refuse a name');
            testCase.verifyEmpty(leftovers);
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before, 'nothing may be created');
            testCase.verifyNumElements(seam.mkdirNames, 1, 'the control is a test file, never a test folder');
            testCase.assertNotEmpty(notes);
            testCase.verifyEqual(notes(end).outcome, 'indeterminate');
            testCase.verifyTrue(contains(notes(end).reason, 'short test file'), ...
                'the reason must say the control failed too');
        end

        function testMissingJavaLeavesTheNameUnrefusedAndSaysWhy(testCase)
            %No Java runtime: a refused test folder cannot be confirmed with a
            %test file, so the name is not refused (the run keeps its existing
            %late handling) and NOTES say why. Taken with a name the volume
            %really cannot hold, so a refusal here would be a missing runtime
            %turned into an invalid name.
            savename = [repmat('o', 1, testCase.Lim.save - 6) '.psyrat'];
            before = OutputNameFixture.snapshot(testCase.Root);
            [~, problems, leftovers, notes] = testCase.runWithSeams( ...
                @() psyrat_check_output_names(testCase.Folder, savename), ...
                struct('mkdirMode', 'fail', 'nojvm', true));
            testCase.verifyEmpty(problems, 'a missing Java runtime must never become an invalid name');
            testCase.verifyEmpty(leftovers);
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before, 'nothing may be created');
            testCase.assertNotEmpty(notes);
            testCase.verifyEqual(notes(1).outcome, 'indeterminate');
            testCase.verifyTrue(contains(notes(1).reason, 'Java'), 'the reason must name the missing runtime');
        end

        function testTheWarningRecorderKeepsOnlyWarningsReallyIssued(testCase)
            %Controls for the warning recorder of runWithSeams (phase 6), on
            %which the two leftover tests rely. Phase 5 read the global last
            %warning after the call; in the lock test that was MATLAB's
            %switched-off name-conflict warning, raised by a later isfile
            %seam call re-joining the path, not the check's warning (which
            %was issued and shown). The recorder instead takes each
            %psyrat:outputNameProbeLeft warning when it is issued.
            target = 'psyrat:outputNameProbeLeft';
            %Positive and filter control: one call issues an unrelated
            %warning, queries the target's state and issues one real target
            %warning. Only the target is recorded, with the identifier and
            %text MATLAB reports, and the query still returns its state.
            [seam, queried] = testCase.runWithSeams(@() localIssueWarnings(target), struct());
            testCase.assertNumElements(seam.warnings, 1, ...
                'exactly the one target warning must be recorded');
            testCase.verifyEqual(seam.warnings(1).id, target);
            testCase.verifyEqual(seam.warnings(1).message, 'recorder control: issued by the test');
            testCase.verifyEqual({seam.lastWarningId, seam.lastWarning}, ...
                {target, 'recorder control: issued by the test'});
            testCase.verifyTrue(isstruct(queried) && isfield(queried, 'identifier') && ...
                strcmp(queried.identifier, target), 'a state query must be forwarded with its output');
            %Negative control on the production call, with the contamination
            %that broke phase 5 present: test names are taken in the meantime
            %by the path-changing isfile seam, the check leaves nothing behind
            %and so issues no leftover warning, and nothing may be recorded.
            components = {'quiet.psyrat', 'quiet_runconfig.json', 'quiet_error.txt'};
            [seam, ~, leftovers] = testCase.runWithSeams( ...
                @() psyrat_check_output_names(testCase.Folder, 'quiet.psyrat'), ...
                struct('mkdirMode', 'fail', 'isfileMode', 'plant', 'skip', {components}));
            testCase.assertEmpty(leftovers, 'premise: the check left nothing behind');
            testCase.assertNotEmpty(seam.acted, 'premise: the path-changing isfile seam ran during the call');
            testCase.assertNotEmpty(seam.globalLastWarningId, ['premise: a warning was raised during ' ...
                'the call (the seams'' own path changes), so a global capture would have seen one']);
            testCase.verifyEmpty(seam.warnings, ...
                'the recorder must not manufacture a warning the check did not issue');
            testCase.verifyEqual({seam.lastWarningId, seam.lastWarning}, {'', ''});
        end

        % ---- review coverage (phase 7) ----------------------------------------

        function testAFileTestMadeInAnotherFolderIsRemovedAndRefusesNothing(testCase)
            %Review F2: Java creates the test file, but MATLAB does not see it
            %in the save folder. Every test folder is refused, and Java is
            %pointed at an existing sibling the test owns.
            %  - An ordinary name: each component's test file is made in the
            %    sibling, is not found in the save folder, is removed, and the
            %    component is left undecided.
            %  - A name past the limit (Java refuses it in any folder): the
            %    short control file is then made in the sibling, fails the same
            %    visibility check and is removed, so nothing is refused.
            %Both folders must end exactly as they began.
            other = fullfile(testCase.Root, 'other');
            mkdir(other);
            seams = struct('mkdirMode', 'fail', 'dirMode', 'sibling', 'dirTarget', other);
            beforeSave = OutputNameFixture.snapshot(testCase.Folder);
            beforeOther = OutputNameFixture.snapshot(other);
            [rec, problems, leftovers, notes] = testCase.runWithSeams( ...
                @() psyrat_check_output_names(testCase.Folder, 'plain.psyrat'), seams);
            testCase.assertEqual(rec.redirected, 3, ...
                'premise: each component''s test file was addressed to the sibling');
            testCase.verifyEmpty(problems, 'a test file MATLAB cannot see decides nothing');
            testCase.verifyEmpty(leftovers, 'each stray test file must be removed');
            testCase.verifyEmpty(rec.warnings, 'nothing was left, so nothing may be announced');
            testCase.verifyEqual({notes.outcome}, {'indeterminate', 'indeterminate', 'indeterminate'});
            testCase.verifyTrue(all(contains({notes.reason}, 'could not be matched')), ...
                'each test file was made (in the sibling) and not seen in the save folder');
            testCase.verifyEqual(OutputNameFixture.snapshot(other), beforeOther, ...
                'the sibling must be exactly as it was');
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Folder), beforeSave, ...
                'the save folder must be exactly as it was');
            %past the limit: the refused test file, then the control in the sibling
            savename = [repmat('o', 1, testCase.Lim.save - 6) '.psyrat'];
            [~, detail] = OutputNameFixture.oracle(testCase.OracleDir, savename);
            testCase.assertFalse(detail.output, ...
                'premise (oracle): the output name is past this volume''s limit');
            [rec, problems, leftovers, notes] = testCase.runWithSeams( ...
                @() psyrat_check_output_names(testCase.Folder, savename), seams);
            testCase.assertEqual(rec.redirected, 2, ...
                'premise: the refused test file and the control were both addressed to the sibling');
            testCase.verifyEmpty(problems, ...
                'a control MATLAB cannot see must not turn a failed test file into a refusal');
            testCase.verifyEmpty(leftovers);
            testCase.verifyEmpty(rec.warnings);
            testCase.assertNumElements(notes, 1);
            testCase.verifyEqual(notes(1).outcome, 'indeterminate');
            testCase.verifyTrue(contains(notes(1).reason, 'short test file'), ...
                'the reason must say the control failed too');
            testCase.verifyEqual(OutputNameFixture.snapshot(other), beforeOther);
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Folder), beforeSave);
        end

        function testARelativeSaveFolderIsTestedInThatFolder(testCase)
            %Review F3: psyrat_run accepts a relative save folder. MATLAB's
            %current folder is set to the test's root and the check receives
            %'save'; only the test folders are refused, so every decision rests
            %on Java test files. This validates a relative MATLAB save folder
            %end to end through the Java route. It does NOT claim to detect a
            %helper that hands Java a relative path: java.io.File documents
            %resolution against user.dir, but root observed on this runtime
            %(2026-09-23) that createNewFile follows MATLAB's current folder,
            %so a relative and an absolute path would both land in Root/save.
            %The writers' verdicts in a same-depth twin fix the truth: an ASCII
            %name whose sidecar is exactly at the measured limit is writable,
            %one character more is not (its output still fits, so the sidecar
            %is refused).
            h = testCase.Lim.sidecar - numel('_runconfig.json');   %stem length
            validName = [repmat('r', 1, h) '.psyrat'];
            invalidName = [repmat('r', 1, h + 1) '.psyrat'];
            twin = fullfile(testCase.Root, 'twin');
            mkdir(twin);
            [okValid, dValid] = OutputNameFixture.oracle(twin, validName);
            [~, dInvalid] = OutputNameFixture.oracle(twin, invalidName);
            rmdir(twin);
            testCase.assertTrue(okValid, sprintf(['premise (writers): the name at the limit ' ...
                'is writable (output %d, sidecar %d, report %d)'], ...
                dValid.output, dValid.sidecar, dValid.report));
            testCase.assertTrue(dInvalid.output && ~dInvalid.sidecar, ['premise (writers): one ' ...
                'more character keeps the output writable and puts the sidecar past the limit']);
            %MATLAB's current folder comes back whatever happens
            startDir = pwd;
            testCase.addTeardown(@() cd(startDir));
            cd(testCase.Root);
            here = pwd;
            javaDir = char(java.lang.System.getProperty('user.dir'));
            %Ownership guard (read-only, before any helper call). Under the
            %documented user.dir resolution, a helper regressed to a relative
            %Java path would address javaDir/save. Where user.dir lies outside
            %the owned root, that folder must not exist, so such a write could
            %only fail, never land outside the root. Where user.dir is the
            %owned root or inside it (a JVM first started after a cd there),
            %any such write stays owned and nothing is required. Nothing is
            %asserted about how this runtime resolves a relative Java path.
            javaCanon = char(java.io.File(javaDir).getCanonicalPath());
            if ~(strcmp(javaCanon, testCase.Root) || startsWith(javaCanon, [testCase.Root filesep]))
                testCase.assertFalse(isfolder(fullfile(javaDir, 'save')), sprintf(['ownership ' ...
                    'premise: %s must not exist, so no relative Java write could land outside ' ...
                    'the test''s root'], fullfile(javaDir, 'save')));
            end
            seams = struct('mkdirMode', 'fail', 'watch', 'save');
            before = OutputNameFixture.snapshot(testCase.Root);
            [rec, problems, leftovers, notes] = testCase.runWithSeams( ...
                @() psyrat_check_output_names('save', validName), seams);
            testCase.assertNumElements(rec.mkdirNames, 3, ...
                'premise: every test folder was refused, so Java decided');
            testCase.verifyEmpty(problems, 'a writable name must not be refused');
            testCase.verifyEmpty(leftovers);
            testCase.verifyEqual({notes.method}, {'file', 'file', 'file'});
            testCase.verifyEqual({notes.outcome}, {'created', 'created', 'created'}, ...
                'each test file must be made where MATLAB sees it, in the relative save folder');
            for k = 1:numel(rec.mkdirNames)
                testCase.verifyFalse(isfile(fullfile(javaDir, 'save', rec.mkdirNames{k})), ...
                    'nothing may be made under Java''s own working folder');
            end
            testCase.verifyEqual(pwd, here, 'the check must not change MATLAB''s current folder');
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before);
            [~, problems, leftovers, notes] = testCase.runWithSeams( ...
                @() psyrat_check_output_names('save', invalidName), seams);
            testCase.assertNumElements(problems, 1, 'the name past the limit must still be refused');
            testCase.verifyEqual(problems.component, [repmat('r', 1, h + 1) '_runconfig.json']);
            testCase.verifyEqual({notes(end).method, notes(end).outcome}, {'file', 'refused'});
            testCase.verifyEmpty(leftovers);
            testCase.verifyEqual(pwd, here);
            testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before);
        end

        function testValidThreeAndFourByteNamesPassTheFileTest(testCase)
            %Review F5: names of three-byte (CJK) and four-byte characters
            %(U+1F600, a surrogate pair in MATLAB), the characters of the
            %2026-09-23 probe, through the Java file test with every test
            %folder refused. Lengths are modest, because no counting rule
            %(bytes, UTF-16 units or code points) is assumed: each name's
            %validity comes from the writers in a same-depth twin. A wrong
            %encoding on the Java route would show as a refusal (the file
            %cannot be made) or as an undecided component (it is made under a
            %name MATLAB does not see).
            astral = native2unicode(uint8([240 159 152 128]), 'UTF-8');
            testCase.assertNumElements(astral, 2, 'premise: MATLAB holds U+1F600 as a surrogate pair');
            cases = { ...
                'three-byte (U+4E2D)', [repmat(char(20013), 1, 60) '_c.psyrat']; ...
                'four-byte (U+1F600)', [repmat(astral, 1, 40) '_e.psyrat']};
            twin = fullfile(testCase.Root, 'twin');
            mkdir(twin);
            state = rng;
            for c = 1:size(cases, 1)
                [label, savename] = cases{c, :};
                [writable, detail] = OutputNameFixture.oracle(twin, savename);
                testCase.log(1, sprintf(['%s: output %d UTF-16 units, %d UTF-8 bytes; ' ...
                    'writers: output %d, sidecar %d, report %d'], label, numel(savename), ...
                    numel(unicode2native(savename, 'UTF-8')), detail.output, detail.sidecar, detail.report));
                testCase.assertTrue(writable, sprintf(['premise (writers at the save folder''s ' ...
                    'depth): the %s name is writable here'], label));
                before = OutputNameFixture.snapshot(testCase.Root);
                [rec, problems, leftovers, notes] = testCase.runWithSeams( ...
                    @() psyrat_check_output_names(testCase.Folder, savename), struct('mkdirMode', 'fail'));
                testCase.assertNumElements(rec.mkdirNames, 3, sprintf( ...
                    'premise (%s): every test folder was refused, so Java decided', label));
                testCase.verifyEmpty(problems, sprintf('%s: a writable name must not be refused', label));
                testCase.verifyEmpty(leftovers, sprintf('%s: no test file may be left', label));
                testCase.verifyEqual({notes.method}, {'file', 'file', 'file'}, label);
                testCase.verifyEqual({notes.outcome}, {'created', 'created', 'created'}, sprintf( ...
                    '%s: MATLAB must see each Java test file under its own name', label));
                testCase.verifyEqual(OutputNameFixture.snapshot(testCase.Root), before, label);
            end
            testCase.verifyEqual(rng, state, 'the check must not draw from the global stream');
            rand;   %positive control: one draw is visible in the same comparison
            testCase.verifyNotEqual(rng, state, 'control: a single draw must change the compared state');
        end

    end

    methods (Access = private)

        function varargout = recordProbes(testCase, fcn, intrude)
            %Run FCN with a recording mkdir shadow first on the path. The
            %first output is the list of names FCN asked mkdir to create in
            %the save folder; further outputs are FCN's own, in order. With
            %INTRUDE, the shadow drops intruder.txt into each folder it has
            %just created there.
            if nargin < 3; intrude = false; end
            global PSYRAT_CON_MKDIR %#ok<GVMIS>
            PSYRAT_CON_MKDIR = struct('folder', testCase.Folder, 'names', {{}}, 'intrude', intrude);
            shadow = fullfile(testCase.ShadowDir, 'mkdir.m');
            fid = fopen(shadow, 'w');
            fprintf(fid, '%s\n', ...
                'function varargout = mkdir(varargin)', ...
                '% TEST SHADOW (TestCheckOutputNames): records mkdir calls in the watched folder, then forwards.', ...
                'global PSYRAT_CON_MKDIR', ...
                'watched = nargin >= 2 && isstruct(PSYRAT_CON_MKDIR) && strcmp(char(varargin{1}), PSYRAT_CON_MKDIR.folder);', ...
                'if watched; PSYRAT_CON_MKDIR.names{end + 1} = char(varargin{2}); end', ...
                'here = fileparts(mfilename(''fullpath''));', ...
                'rmpath(here);', ...
                'restore = onCleanup(@() addpath(here, ''-begin'')); %#ok<NASGU>', ...
                'if nargout == 0', ...
                '    mkdir(varargin{:});', ...
                'else', ...
                '    [varargout{1:nargout}] = mkdir(varargin{:});', ...
                'end', ...
                'if watched && PSYRAT_CON_MKDIR.intrude', ...
                '    made = fullfile(char(varargin{1}), char(varargin{2}));', ...
                '    if isfolder(made)', ...
                '        f = fopen(fullfile(made, ''intruder.txt''), ''w''); fprintf(f, ''not the check''''s''); fclose(f);', ...
                '    end', ...
                'end', ...
                'end');
            fclose(fid);
            cleaner = onCleanup(@() testCase.removeShadow()); %#ok<NASGU>
            addpath(testCase.ShadowDir, '-begin');
            rehash;
            testCase.assertEqual(which('mkdir'), shadow, 'the mkdir shadow is not first on the path');
            nout = max(nargout - 1, 0);
            outs = cell(1, nout);
            if nout == 0
                fcn();
            else
                [outs{:}] = fcn();
            end
            testCase.removeShadow();
            names = PSYRAT_CON_MKDIR.names;
            clear('global', 'PSYRAT_CON_MKDIR');
            varargout = [{names}, outs];
        end

        function removeShadow(testCase)
            if any(strcmp(strsplit(path, pathsep), testCase.ShadowDir))
                rmpath(testCase.ShadowDir);
            end
            rehash;
        end

        function tf = limitsFound(testCase)
            tf = testCase.Lim.saveFound && testCase.Lim.sidecarFound;
        end

        function varargout = runWithSeams(testCase, fcn, seam)
            %Run FCN with test seams first on the path, for this one call only.
            %SEAM fields, all optional:
            %  mkdirMode  'forward' (default), or 'fail': a test folder in the
            %             save folder is refused with MATLAB:FILE:PATHMAX and
            %             nothing is created, so the file test must decide;
            %             'failFirst': only the first one is refused;
            %  dirMode    'forward' (default), or 'elsewhere': dir of the save
            %             folder reports a missing absolute folder, so the Java
            %             route misses the save folder; 'sibling': it reports
            %             DIRTARGET, an existing folder the test owns inside its
            %             root (asserted), so Java creates its test files there;
            %  watch      the save-folder string the seams match, exactly as the
            %             check receives it (default: the absolute save folder;
            %             the relative-folder test passes 'save');
            %  isfileMode 'none' (default); 'plant': a name seen free is created
            %             at once by "another process"; 'intrude': foreign data
            %             is written into a test file once it exists; 'lock':
            %             the save folder stops accepting changes once a test
            %             file exists. Each acts once per path, never for the
            %             names in SKIP;
            %  nojvm      true: usejava('jvm') answers false.
            %A warning recorder is always installed. It forwards every call to
            %the built-in warning unchanged (queries and outputs included, and
            %without touching the path). For psyrat:outputNameProbeLeft it also
            %records, right after the real warning was issued, the identifier
            %and text MATLAB itself reports as the last warning.
            %The first output is the seam record (mkdirNames, acted, redirected,
            %warnings, and lastWarning/lastWarningId, the last RECORDED
            %warning); further outputs are FCN's own. The global last warning
            %at the end of the call is NOT the check's: every call of a seam
            %that forwards by leaving and re-joining the path raises MATLAB's
            %name-conflict warning, which is switched off here for display but
            %still becomes the last warning (observed in phase 5). It is kept as
            %globalLastWarning/globalLastWarningId for diagnosis only. The
            %caller's path, current folder, warning states and last warning are
            %restored on exit, failure included.
            defaults = struct('mkdirMode', 'forward', 'dirMode', 'forward', 'dirTarget', '', ...
                'watch', testCase.Folder, 'isfileMode', 'none', 'skip', {{}}, 'nojvm', false);
            for f = fieldnames(defaults)'
                if ~isfield(seam, f{1}); seam.(f{1}) = defaults.(f{1}); end
            end
            %the caller's last warning and current folder come back on exit,
            %failure included
            [prevMsg, prevId] = lastwarn;
            lastRestore = onCleanup(@() lastwarn(prevMsg, prevId)); %#ok<NASGU>
            startDir = pwd;
            dirRestore = onCleanup(@() cd(startDir)); %#ok<NASGU>
            %the dir seam's target: a missing folder, or an existing folder the
            %test owns beside the save folder; never anything outside the root
            switch seam.dirMode
                case 'elsewhere'
                    dirTarget = fullfile(seam.watch, 'missing_by_the_test');
                case 'sibling'
                    dirTarget = char(seam.dirTarget);
                    testCase.assertTrue(isfolder(dirTarget) && ...
                        startsWith(dirTarget, [testCase.Root filesep]) && ...
                        ~strcmp(dirTarget, testCase.Folder), ['the sibling target must be ' ...
                        'an existing folder inside the test''s own root, not the save folder']);
                otherwise
                    dirTarget = '';
            end
            global PSYRAT_CON_SEAM %#ok<GVMIS>
            PSYRAT_CON_SEAM = struct('folder', seam.watch, 'mkdirMode', seam.mkdirMode, ...
                'dirMode', seam.dirMode, 'dirTarget', dirTarget, ...
                'isfileMode', seam.isfileMode, 'skip', {seam.skip}, ...
                'nojvm', seam.nojvm, 'mkdirNames', {{}}, 'acted', {{}}, 'redirected', 0, ...
                'warnId', 'psyrat:outputNameProbeLeft', ...
                'warnings', {struct('id', {}, 'message', {})});
            %the shadow folder is the test's own: clear any earlier shadow
            old = dir(fullfile(testCase.ShadowDir, '*.m'));
            for k = 1:numel(old)
                delete(fullfile(old(k).folder, old(k).name));
            end
            forward = { ...
                'here = fileparts(mfilename(''fullpath''));', ...
                'rmpath(here);', ...
                'restore = onCleanup(@() addpath(here, ''-begin'')); %#ok<NASGU>'};
            written = {'mkdir'};
            localWrite(fullfile(testCase.ShadowDir, 'mkdir.m'), [{ ...
                'function varargout = mkdir(varargin)', ...
                '% TEST SEAM (TestCheckOutputNames): records mkdir calls in the watched folder; in fail mode refuses them (failFirst: only the first) without creating anything.', ...
                'global PSYRAT_CON_SEAM', ...
                'watched = nargin >= 2 && isstruct(PSYRAT_CON_SEAM) && strcmp(char(varargin{1}), PSYRAT_CON_SEAM.folder);', ...
                'if watched', ...
                '    PSYRAT_CON_SEAM.mkdirNames{end + 1} = char(varargin{2});', ...
                '    if strcmp(PSYRAT_CON_SEAM.mkdirMode, ''fail'') || (strcmp(PSYRAT_CON_SEAM.mkdirMode, ''failFirst'') && numel(PSYRAT_CON_SEAM.mkdirNames) == 1)', ...
                '        out = {false, ''Maximum path length exceeded. (test seam)'', ''MATLAB:FILE:PATHMAX''};', ...
                '        varargout = out(1:max(nargout, 1));', ...
                '        return;', ...
                '    end', ...
                'end'}, forward, { ...
                'if nargout == 0', ...
                '    mkdir(varargin{:});', ...
                'else', ...
                '    [varargout{1:nargout}] = mkdir(varargin{:});', ...
                'end', ...
                'end'}]);
            %The warning recorder forwards through builtin, NOT the path
            %pattern above: leaving and re-joining the path is what raises
            %the name-conflict warning that overwrote the last warning. For
            %the recorded identifier, the last warning is cleared, the real
            %warning is issued, and MATLAB's own last warning is read back,
            %so the recorded identifier and text come from the warning
            %actually issued, never from the recorder.
            written{end + 1} = 'warning';
            localWrite(fullfile(testCase.ShadowDir, 'warning.m'), { ...
                'function varargout = warning(varargin)', ...
                '% TEST SEAM (TestCheckOutputNames): forwards every call to the built-in warning unchanged; for the recorded identifier it also records MATLAB''s own last warning right after the real one was issued.', ...
                'global PSYRAT_CON_SEAM', ...
                'rec = nargin >= 2 && isstruct(PSYRAT_CON_SEAM) && isfield(PSYRAT_CON_SEAM, ''warnId'') && (ischar(varargin{1}) || (isstring(varargin{1}) && isscalar(varargin{1}))) && strcmp(char(varargin{1}), PSYRAT_CON_SEAM.warnId);', ...
                'if rec', ...
                '    lastwarn('''', '''');', ...
                'end', ...
                'if nargout == 0', ...
                '    builtin(''warning'', varargin{:});', ...
                'else', ...
                '    [varargout{1:nargout}] = builtin(''warning'', varargin{:});', ...
                'end', ...
                'if rec', ...
                '    [msg, id] = lastwarn;', ...
                '    PSYRAT_CON_SEAM.warnings(end + 1) = struct(''id'', id, ''message'', msg);', ...
                'end', ...
                'end'});
            if ~isempty(dirTarget)
                written{end + 1} = 'dir';
                localWrite(fullfile(testCase.ShadowDir, 'dir.m'), [{ ...
                    'function varargout = dir(varargin)', ...
                    '% TEST SEAM (TestCheckOutputNames): dir of the watched save folder reports the seam''s target folder (missing, or an owned sibling); forwards everything else.', ...
                    'global PSYRAT_CON_SEAM'}, forward, { ...
                    'if nargout == 0', ...
                    '    dir(varargin{:});', ...
                    '    return;', ...
                    'end', ...
                    '[varargout{1:nargout}] = dir(varargin{:});', ...
                    'if nargin == 1 && isstruct(PSYRAT_CON_SEAM) && strcmp(char(varargin{1}), PSYRAT_CON_SEAM.folder) && ~isempty(varargout{1})', ...
                    '    d = varargout{1};', ...
                    '    [d.folder] = deal(PSYRAT_CON_SEAM.dirTarget);', ...
                    '    varargout{1} = d;', ...
                    '    PSYRAT_CON_SEAM.redirected = PSYRAT_CON_SEAM.redirected + 1;', ...
                    'end', ...
                    'end'}]);
            end
            if ~strcmp(seam.isfileMode, 'none')
                written{end + 1} = 'isfile';
                localWrite(fullfile(testCase.ShadowDir, 'isfile.m'), [{ ...
                    'function tf = isfile(varargin)', ...
                    '% TEST SEAM (TestCheckOutputNames): forwards; for a watched test entry it plants a file, writes foreign data or locks the folder, once per path.', ...
                    'global PSYRAT_CON_SEAM'}, forward, { ...
                    'tf = isfile(varargin{:});', ...
                    'if nargin ~= 1 || ~isstruct(PSYRAT_CON_SEAM); return; end', ...
                    'p = char(varargin{1});', ...
                    '[folder, n, e] = fileparts(p);', ...
                    'if ~strcmp(folder, PSYRAT_CON_SEAM.folder) || any(strcmp([n e], PSYRAT_CON_SEAM.skip)) || any(strcmp(p, PSYRAT_CON_SEAM.acted)); return; end', ...
                    'switch PSYRAT_CON_SEAM.isfileMode', ...
                    '    case ''plant''', ...
                    '        if ~tf; f = fopen(p, ''w''); fprintf(f, ''planted by the test''); fclose(f); PSYRAT_CON_SEAM.acted{end + 1} = p; end', ...
                    '    case ''intrude''', ...
                    '        if tf; f = fopen(p, ''a''); fprintf(f, ''foreign data''); fclose(f); PSYRAT_CON_SEAM.acted{end + 1} = p; end', ...
                    '    case ''lock''', ...
                    '        if tf; fileattrib(folder, ''-w'', ''u''); PSYRAT_CON_SEAM.acted{end + 1} = p; end', ...
                    'end', ...
                    'end'}]);
            end
            if seam.nojvm
                written{end + 1} = 'usejava';
                localWrite(fullfile(testCase.ShadowDir, 'usejava.m'), [{ ...
                    'function tf = usejava(feature)', ...
                    '% TEST SEAM (TestCheckOutputNames): answers usejava(''jvm'') with false; forwards every other feature.', ...
                    'global PSYRAT_CON_SEAM', ...
                    'if isstruct(PSYRAT_CON_SEAM) && PSYRAT_CON_SEAM.nojvm && strcmpi(char(feature), ''jvm'')', ...
                    '    tf = false;', ...
                    '    return;', ...
                    'end'}, forward, { ...
                    'tf = usejava(feature);', ...
                    'end'}]);
            end
            cleaner = onCleanup(@() testCase.removeShadow()); %#ok<NASGU>
            lastwarn('');
            addpath(testCase.ShadowDir, '-begin');
            rehash;
            [~, conflictId] = lastwarn;
            quiet = [];
            if ~isempty(conflictId)
                state = warning('off', conflictId);
                quiet = onCleanup(@() warning(state)); %#ok<NASGU>
            end
            for k = 1:numel(written)
                testCase.assertEqual(which(written{k}), fullfile(testCase.ShadowDir, [written{k} '.m']), ...
                    sprintf('the %s seam is not first on the path', written{k}));
            end
            lastwarn('');
            nout = max(nargout - 1, 0);
            outs = cell(1, nout);
            if nout == 0
                fcn();
            else
                [outs{:}] = fcn();
            end
            [lastmsg, lastid] = lastwarn;
            testCase.removeShadow();
            record = PSYRAT_CON_SEAM;
            %the check's own leftover warning as the recorder took it when it
            %was issued; the global last warning is kept for diagnosis only
            record.lastWarning = '';
            record.lastWarningId = '';
            if ~isempty(record.warnings)
                record.lastWarning = record.warnings(end).message;
                record.lastWarningId = record.warnings(end).id;
            end
            record.globalLastWarning = lastmsg;
            record.globalLastWarningId = lastid;
            clear('global', 'PSYRAT_CON_SEAM');
            varargout = [{record}, outs];
        end

    end
end

function localWrite(file, lines)
fid = fopen(file, 'w');
assert(fid >= 3, 'cannot write the test seam %s', file);
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end

function name = localName(p)
[~, n, e] = fileparts(p);
name = [n e];
end

function queried = localIssueWarnings(target)
%For the recorder control: an unrelated warning, a state query with an
%output, and one real TARGET warning, each through whatever warning is first
%on the path (the recorder, inside runWithSeams)
warning('psyrat:recorderControlOther', 'recorder control: not the target');
queried = warning('query', target);
warning(target, 'recorder control: %s', 'issued by the test');
end
