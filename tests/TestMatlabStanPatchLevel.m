classdef TestMatlabStanPatchLevel < PsyRATTestBase
    % Covers the vendored-MatlabStan patch-level stamp (audit finding B41) and
    % the delete guard that had to go in alongside it.
    %
    % WHY BOTH LIVE IN ONE CLASS. They are one change. Making
    % psyrat_checkversionsofdeps report a real verdict for MatlabStan arms a
    % previously dormant recursive delete in psyrat_installdependents:
    %
    %     if depvercheck.matlabstan == 0
    %         rmdir(fileparts(which('mcmc.m')),'s')
    %
    % which is dead only while that field is hardcoded to -1. The delete target
    % is resolved from the Matlab path, and PsyRAT puts its own vendored bundle
    % on the path under the same folder name the installed copy uses, so on a
    % repository checkout that rmdir resolves to PsyRAT's own source. The stamp
    % and the guard are therefore tested together: the stamp is what arms the
    % wire, and psyrat_safe_rmdepdir is what makes arming it safe.
    %
    % All of this runs without CmdStan, MatlabStan, or any toolbox: the bundle
    % reader takes an explicit directory, so every case is exercised against a
    % temporary fixture rather than the machine's real installation.

    methods (Test)

        % -- psyrat_matlabstan_bundleinfo -------------------------------------

        function testCurrentStampReportsUpToDate(testCase)
            depvers = psyrat_dependentsversions;
            msdir = testCase.makeStampedDir(depvers.matlabstan, ...
                depvers.matlabstanpatch);

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyTrue(info.found);
            testCase.verifyEqual(info.status, 1);
            testCase.verifyEqual(info.version, depvers.matlabstan);
            testCase.verifyEqual(info.patchlevel, depvers.matlabstanpatch);
        end

        function testMissingStampIsReportedOldNotUnknown(testCase)
            % The case B41 is actually about: a directory installed before
            % stamping existed. It cannot contain PsyRAT's local patches, so
            % reporting it as unknown would leave the user with no signal.
            msdir = testCase.makeTempDir();

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyFalse(info.found);
            testCase.verifyEqual(info.status, 0);
            testCase.verifyNotEmpty(info.detail);
        end

        function testOlderPatchLevelIsOld(testCase)
            depvers = psyrat_dependentsversions;
            older = num2str(str2double(depvers.matlabstanpatch) - 1);
            msdir = testCase.makeStampedDir(depvers.matlabstan, older);

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyTrue(info.found);
            testCase.verifyEqual(info.status, 0);
        end

        function testNewerPatchLevelIsUntestedNew(testCase)
            depvers = psyrat_dependentsversions;
            newer = num2str(str2double(depvers.matlabstanpatch) + 1);
            msdir = testCase.makeStampedDir(depvers.matlabstan, newer);

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyEqual(info.status, 2);
        end

        function testUpstreamReleaseDominatesPatchLevel(testCase)
            % A different MatlabStan release is a bigger difference than a patch
            % level within one release, and patch levels across releases are not
            % comparable. An older release with an ABSURDLY high patch level
            % must still read as old, or the dominance rule is not being applied.
            msdir = testCase.makeStampedDir('1.0.0.0', '99999');

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyEqual(info.status, 0);
            testCase.verifyTrue(contains(info.detail, '1.0.0.0'));
        end

        function testFourComponentReleaseCompares(testCase)
            % MatlabStan's release string has FOUR components ('2.15.1.0'). The
            % CmdStan and MatlabProcessManager checks both loop over exactly
            % three, so a three-component comparison would silently ignore the
            % last component. Positive control: a version differing ONLY in the
            % fourth component must be detected.
            msdir = testCase.makeStampedDir('2.15.1.5', '1');

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyEqual(info.status, 2, ...
                ['A release differing only in the fourth component was not ' ...
                'detected, so the comparison is truncating to three parts.']);
        end

        function testTrailingZeroComponentsCompareEqual(testCase)
            % Missing trailing components are zero-filled, so '2.15.1' and
            % '2.15.1.0' are the same release.
            %
            % Phrased as "the declared version with an extra .0 appended" rather
            % than against a hardcoded '2.15.1.0', so that it keeps testing the
            % zero-fill after any version bump. An assumeEqual on the literal
            % would have turned this into a silent skip instead, which is the
            % green-while-skipping failure mode this repo has been bitten by.
            depvers = psyrat_dependentsversions;
            msdir = testCase.makeStampedDir([depvers.matlabstan '.0'], ...
                depvers.matlabstanpatch);

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyEqual(info.status, 1);
        end

        function testUnparseableStampIsUnknownNotOld(testCase)
            % A stamp we cannot read suggests hand editing. Guessing a verdict
            % from a partial read would be worse than admitting ignorance.
            msdir = testCase.makeTempDir();
            fid = fopen(fullfile(msdir, 'PSYRAT_PATCHLEVEL.txt'), 'w');
            fprintf(fid, 'this file has been mangled\n');
            fclose(fid);

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyFalse(info.found);
            testCase.verifyEqual(info.status, -1);
        end

        function testProseMentioningFieldNamesIsNotParsedAsPayload(testCase)
            % The real stamp discusses 'patch_level' and 'upstream_version' in
            % its own explanatory prose. The reader anchors each field to the
            % start of a line, so prose must not be mistaken for the payload.
            msdir = testCase.makeTempDir();
            fid = fopen(fullfile(msdir, 'PSYRAT_PATCHLEVEL.txt'), 'w');
            fprintf(fid, 'Discussion of patch_level = 999 in a sentence.\n');
            fprintf(fid, '    upstream_version = 2.15.1.0\n');
            fprintf(fid, '    patch_level = 1\n');
            fclose(fid);

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyEqual(info.patchlevel, '1');
        end

        function testShippedStampMatchesDeclaredExpectation(testCase)
            % The two halves must be bumped together. If they drift, every user
            % is told to reinstall for nothing, or a real patch ships silently.
            %
            % G56d made the VERDICT insensitive to in-tree drift on purpose
            % (the bundle reports current/pinned even when its stamp is
            % stale, because a dependents update cannot act on it), so this
            % pin compares the RAW stamp fields to the declared expectation
            % rather than going through the status -- a status-1 assert
            % would now pass vacuously on a drifted stamp.
            depvers = psyrat_dependentsversions;
            stamp = fullfile(testCase.projectRoot(), 'bundled_dependents', ...
                ['MatlabStan-' depvers.matlabstan], 'PSYRAT_PATCHLEVEL.txt');
            testCase.assertEqual(exist(stamp, 'file'), 2, ...
                'The shipped bundle has no PSYRAT_PATCHLEVEL.txt.');

            info = psyrat_matlabstan_bundleinfo(fileparts(stamp));

            testCase.assertTrue(info.found, ...
                'The shipped stamp did not parse.');
            testCase.verifyEqual(info.version, depvers.matlabstan, ...
                ['The shipped stamp''s upstream_version and ' ...
                'psyrat_dependentsversions disagree. Bump them in the ' ...
                'SAME commit.']);
            testCase.verifyEqual(info.patchlevel, depvers.matlabstanpatch, ...
                ['The shipped stamp''s patch_level and ' ...
                'depvers.matlabstanpatch disagree. Bump ' ...
                'PSYRAT_PATCHLEVEL.txt and depvers.matlabstanpatch in the ' ...
                'SAME commit.']);
        end

        % -- G56d: the in-tree pinned-bundle verdict --------------------------

        function testStaleStampInTreeReportsPinnedVerdict(testCase)
            % G56d. A would-be-old verdict for a copy INSIDE the toolbox tree
            % reports 3 (bundled, pinned) instead: psyrat_safe_rmdepdir
            % refuses to remove the bundle by design, so an "old" verdict
            % there re-raised the launch update offer every session with an
            % offer that could never succeed (the nag loop). The verdict is a
            % DISTINCT value, not a bare 1, so the launch console can print
            % the pinned-copy line rather than a false "version is current"
            % (wave: info.detail has no production consumer, so the scalar
            % verdict is the only honest channel). The out-of-tree twin of
            % this stamp still reports 0 (testOlderPatchLevelIsOld), the
            % positive control that the override is keyed on containment,
            % not on the stamp.
            depvers = psyrat_dependentsversions;
            older = num2str(str2double(depvers.matlabstanpatch) - 1);
            msdir = testCase.makeInTreeStampedDir(depvers.matlabstan, older);

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyEqual(info.status, 3, ...
                'A stale in-tree stamp must report the pinned verdict, not old.');
            testCase.verifyTrue(contains(info.detail, 'pinned'), ...
                'The detail line must say the copy is the pinned bundle.');
            testCase.verifyTrue(contains(info.detail, 'toolbox'), ...
                'The detail line must point at a toolbox update instead.');
        end

        function testOlderReleaseInTreeReportsPinnedVerdict(testCase)
            % Same override on the release-comparison arm (the upstream
            % release dominates the patch level, so this arm returns before
            % the patch comparison runs). Out-of-tree control:
            % testUpstreamReleaseDominatesPatchLevel.
            msdir = testCase.makeInTreeStampedDir('1.0.0.0', '99999');

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyEqual(info.status, 3);
            testCase.verifyTrue(contains(info.detail, 'pinned'));
        end

        function testMissingStampInTreeReportsPinnedVerdict(testCase)
            % Same override on the no-stamp arm. Out-of-tree control:
            % testMissingStampIsReportedOldNotUnknown.
            msdir = testCase.makeInTreeDir();

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyEqual(info.status, 3);
            testCase.verifyTrue(contains(info.detail, 'pinned'));
        end

        function testPinnedVerdictConsoleAndGateWelds(testCase)
            % G56d wave: the pinned verdict is only honest if psyrat_start
            % actually prints a distinct line for it AND keeps gating the
            % update dialog on == 0 (3 must fall outside every needs-update
            % check). Source welds, since the launch console is GUI-coupled.
            src = fileread(fullfile(testCase.projectRoot(), 'psyrat_start.m'));
            iSwitch = strfind(src, 'switch depvercheck.matlabstan');
            testCase.assertEqual(numel(iSwitch), 1, ...
                'exactly one console switch on the MatlabStan verdict');
            iCase3 = strfind(src, 'case 3');
            testCase.assertNotEmpty(iCase3, ...
                'the console switch must carry a case 3 arm');
            testCase.assertNotEmpty(strfind(src, ...
                'bundled pinned copy is in use '), ...
                ['the case 3 arm must print the pinned-copy line -- the ' ...
                'scalar verdict is the only honest user-facing channel']);
            testCase.assertNotEmpty(strfind(src, ...
                'depvercheck.matlabstan == 0'), ...
                ['the update-dialog gates must stay keyed on == 0, which ' ...
                'verdict 3 deliberately fails']);
        end

        function testContainmentHelperMirrorsTheRemover(testCase)
            % G56d wave: the verdict's canonicalizer is a COPY of
            % psyrat_safe_rmdepdir's local_canonical, and nothing else stops
            % the two from drifting apart -- a divergence would let the
            % verdict and the remover disagree about the same directory.
            % Welded on the executable lines (the docstrings differ on
            % purpose).
            remover = localCanonicalCodeLines(fullfile( ...
                testCase.projectRoot(), 'subroutines', 'installation', ...
                'psyrat_safe_rmdepdir.m'));
            verdict = localCanonicalCodeLines(fullfile( ...
                testCase.projectRoot(), 'subroutines', 'installation', ...
                'psyrat_matlabstan_bundleinfo.m'));
            testCase.assertNotEmpty(remover, ...
                'the remover''s local_canonical must exist');
            testCase.assertNotEmpty(verdict, ...
                'the verdict''s local_canonical must exist');
            testCase.verifyEqual(verdict, remover, ...
                ['local_canonical has diverged between ' ...
                'psyrat_matlabstan_bundleinfo and psyrat_safe_rmdepdir; ' ...
                'change both together or the verdict and the remover can ' ...
                'disagree about the same directory']);
        end

        function testNewerInTreeStampStaysUntestedNew(testCase)
            % The override is scoped to would-be-OLD verdicts only: a NEWER
            % in-tree stamp (a development tree mid-bump) still warns as 2.
            depvers = psyrat_dependentsversions;
            newer = num2str(str2double(depvers.matlabstanpatch) + 1);
            msdir = testCase.makeInTreeStampedDir(depvers.matlabstan, newer);

            info = psyrat_matlabstan_bundleinfo(msdir);

            testCase.verifyEqual(info.status, 2);
        end

        % -- psyrat_safe_rmdepdir ---------------------------------------------

        function testRefusesToDeleteInsideTheToolbox(testCase)
            % The whole reason the guard exists. which('mcmc.m') resolves to the
            % in-repo bundle on a checkout, so an unguarded rmolddeps would
            % delete PsyRAT's own vendored source and every local patch with it.
            target = fullfile(testCase.projectRoot(), 'bundled_dependents');
            testCase.assertEqual(exist(target, 'dir'), 7);

            [removed, reason] = psyrat_safe_rmdepdir(target, 'MatlabStan');

            testCase.verifyFalse(removed);
            testCase.verifyNotEmpty(reason);
            testCase.verifyEqual(exist(target, 'dir'), 7, ...
                'The guard was supposed to refuse, but the directory is gone.');
        end

        function testRefusesTheToolboxRootItself(testCase)
            root = testCase.projectRoot();

            [removed, ~] = psyrat_safe_rmdepdir(root, 'MatlabStan');

            testCase.verifyFalse(removed);
            testCase.verifyEqual(exist(root, 'dir'), 7);
        end

        function testRefusesEmptyTarget(testCase)
            % fileparts(which('missing.m')) is '', and rmdir('','s') must never
            % be reached with it.
            [removed, reason] = psyrat_safe_rmdepdir('', 'MatlabStan');

            testCase.verifyFalse(removed);
            testCase.verifyNotEmpty(reason);
        end

        function testRefusesNonexistentTarget(testCase)
            [removed, reason] = psyrat_safe_rmdepdir( ...
                fullfile(tempdir, 'psyrat_no_such_dir_9c3f1a'), 'CmdStan');

            testCase.verifyFalse(removed);
            testCase.verifyNotEmpty(reason);
        end

        function testRemovesAnInstalledCopyOutsideTheToolbox(testCase)
            % Positive control. Without this the refusal tests above would pass
            % just as well against a guard that refuses EVERYTHING, which would
            % silently break dependency updates rather than protect them.
            target = testCase.makeTempDir();
            nested = fullfile(target, 'nested');
            mkdir(nested);
            fid = fopen(fullfile(nested, 'file.txt'), 'w');
            fprintf(fid, 'x\n');
            fclose(fid);

            [removed, reason] = psyrat_safe_rmdepdir(target, 'MatlabStan');

            testCase.verifyTrue(removed, reason);
            testCase.verifyNotEqual(exist(target, 'dir'), 7);
        end

        function testSiblingDirectoryIsNotTreatedAsInside(testCase)
            % A prefix test without a trailing separator would treat
            % '<root>_other' as living inside '<root>'. That would refuse a
            % legitimate removal, so the boundary must compare path segments.
            root = testCase.projectRoot();
            sibling = [root '_psyrat_guard_probe'];
            testCase.assumeEqual(exist(sibling, 'dir'), 0, ...
                'Probe directory already exists; skipping to avoid touching it.');
            mkdir(sibling);
            c = onCleanup(@() localForceRemove(sibling));

            [removed, reason] = psyrat_safe_rmdepdir(sibling, 'MatlabStan');

            testCase.verifyTrue(removed, reason);
        end

    end

    methods (Access = private)

        function d = makeTempDir(testCase)
            d = tempname;
            mkdir(d);
            testCase.addTeardown(@() localForceRemove(d));
        end

        function d = makeStampedDir(testCase, upstreamVersion, patchLevel)
            d = testCase.makeTempDir();
            testCase.writeStamp(d, upstreamVersion, patchLevel);
        end

        function d = makeInTreeDir(testCase)
            % A scratch directory INSIDE the toolbox tree, for the G56d
            % containment pins. Fixed name (deterministic tests carry no
            % RNG); the suite runs serially, and the teardown removes it.
            %
            % Wave: a leftover from a crashed run is PRE-CLEANED and then
            % ASSERTED gone, never assumed -- runAllUnitTests subtracts
            % assumption-filtered skips from its green signal, so an assume
            % here would let all the in-tree pins skip silently forever
            % (the green-while-skipping failure mode this class's own
            % zero-fill test warns about). The directory holds only a
            % fixture stamp, so removing a leftover destroys nothing.
            d = fullfile(testCase.projectRoot(), 'tests', ...
                'tmp_g56d_intree_fixture');
            localForceRemove(d);
            testCase.assertEqual(exist(d, 'dir'), 0, ...
                ['A leftover in-tree fixture directory could not be ' ...
                'removed; failing loudly rather than skipping.']);
            mkdir(d);
            testCase.addTeardown(@() localForceRemove(d));
        end

        function d = makeInTreeStampedDir(testCase, upstreamVersion, patchLevel)
            d = testCase.makeInTreeDir();
            testCase.writeStamp(d, upstreamVersion, patchLevel);
        end

        function writeStamp(~, d, upstreamVersion, patchLevel)
            fid = fopen(fullfile(d, 'PSYRAT_PATCHLEVEL.txt'), 'w');
            fprintf(fid, 'PsyRAT vendored-bundle stamp\n\n');
            fprintf(fid, '    upstream_version = %s\n', upstreamVersion);
            fprintf(fid, '    patch_level = %s\n', patchLevel);
            fclose(fid);
        end

    end

end

function lines = localCanonicalCodeLines(file)
%LOCALCANONICALCODELINES Executable lines of a file's local_canonical.
%   Slices from the function header to the next 'function ' line (or end of
%   file -- the helper is last in psyrat_matlabstan_bundleinfo), then drops
%   comment-only and blank lines and trailing whitespace, so the weld
%   compares code and tolerates the deliberately different docstrings.
src = fileread(file);
marker = 'function p = local_canonical(p)';
iStart = strfind(src, marker);
lines = {};
if numel(iStart) ~= 1
    return;
end
seg = src(iStart+numel(marker):end);
iNext = strfind(seg, [newline 'function ']);
if ~isempty(iNext)
    seg = seg(1:iNext(1));
end
raw = strsplit(seg, newline, 'CollapseDelimiters', false);
for ii = 1:numel(raw)
    ln = deblank(raw{ii});
    tl = strtrim(ln);
    if isempty(tl) || tl(1) == '%'
        continue;
    end
    lines{end+1} = tl; %#ok<AGROW>
end
end

function localForceRemove(d)
%LOCALFORCEREMOVE Best-effort teardown that never fails a passing test.
if exist(d, 'dir') == 7
    warning('off', 'MATLAB:RMDIR:RemovedFromPath');
    try %#ok<TRYNC>
        rmdir(d, 's');
    end
    warning('on', 'MATLAB:RMDIR:RemovedFromPath');
end
end
