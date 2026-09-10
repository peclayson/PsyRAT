classdef TestCmdstanDiscovery < PsyRATTestBase
    %G43: the repo carried four divergent CmdStan discovery orders (the
    %bundled mstan.stan_home; psyrat_detect_cmdstan_version's path-then-env
    %lookup with reversed env names; psyrat_checkversionsofdeps' path-only
    %read; and the repair template's which()-before-globs emission). The fix
    %centralizes toolbox-side discovery on psyrat_locate_cmdstan, which
    %delegates to mstan.stan_home -- the resolver estimation actually uses --
    %and adds the G43 ruling's newest-version-first ordering to the managed
    %home globs (patch level 2 on the bundled MatlabStan).
    %
    %These tests sandbox discovery three ways: the PSYRAT_TEST_CMDSTAN_PATH
    %environment override (tier 1), a java user.home redirect for the managed
    %glob tier (tier 2), and a fake stanc executable for the version read.
    %The four CmdStan env vars are saved and restored around every test so a
    %machine's real configuration is untouched.
    %
    %G51 extends the coverage to the installer's three per-OS CmdStan gates
    %and the workflow sweep's env seeding, which used to path-gate CmdStan
    %with exist()/which() probes the resolver contradicted; those are pinned
    %as source-text welds below (GUI-coupled local functions and
    %pre-test-run function-body code no headless route executes).
    %
    %NOT covered headlessly: the psyrat_start gate's install dialog (a GUI
    %surface; the env-var-only launch check is queued for the P4
    %click-through), psyrat_checkversionsofdeps' guide.tex legacy branch,
    %and the installer gates' runtime behavior (the welds pin the source,
    %not an execution).

    properties (Access = private)
        SavedEnv
    end

    properties (Constant, Access = private)
        EnvNames = {'PSYRAT_TEST_CMDSTAN_PATH','STAN_HOME','CMDSTAN_HOME','CMDSTAN'};
    end

    methods (TestMethodSetup)
        function saveEnv(testCase)
            testCase.SavedEnv = cellfun(@getenv, testCase.EnvNames, ...
                'UniformOutput', false);
            %psyrat_locate_cmdstan delegates to mstan.stan_home, which lives
            %in bundled_dependents -- a tree PsyRATTestBase does NOT add. On
            %a clean path the resolver would return '' and every discovery
            %test here would fail (or pass vacuously) for the wrong reason,
            %so the prerequisite is declared explicitly; the base class
            %restores the path after each method.
            addpath(genpath(fullfile(testCase.projectRoot(), ...
                'bundled_dependents')));
        end
    end

    methods (TestMethodTeardown)
        function restoreEnv(testCase)
            for k = 1:numel(testCase.EnvNames)
                setenv(testCase.EnvNames{k}, testCase.SavedEnv{k});
            end
        end
    end

    methods (Test)

        function testEnvOverrideWinsAndValidates(testCase)
            % Tier 1: a valid root named by PSYRAT_TEST_CMDSTAN_PATH (the
            % first env name) must win over everything on this machine.
            fake = localMakeFakeCmdstanRoot(testCase, 'cmdstan-env');
            setenv('PSYRAT_TEST_CMDSTAN_PATH', fake);

            home = psyrat_locate_cmdstan();
            testCase.verifyEqual(home, fake, ...
                'a valid env-var root must win the discovery outright');
        end

        function testInvalidEnvRootIsRejected(testCase)
            % Validity is part of discovery: an env var pointing at a
            % directory with no bin/ must NOT be returned. A VALID fake in
            % the next env tier (STAN_HOME) is the anti-vacuity control and
            % makes the outcome machine-independent: the resolver must skip
            % the invalid tier-1 root and land on the tier-2 fake, never on
            % '' and never on whatever real install this machine has (a
            % first draft asserted only ~= fakeBad, which passes vacuously
            % on '' -- adversarial wave, 2026-08-29).
            sandbox = localSandbox(testCase);
            fakeBad = fullfile(sandbox, 'not-cmdstan');
            mkdir(fakeBad);
            fid = fopen(fullfile(fakeBad, 'makefile'), 'w'); fclose(fid);
            % no bin/ directory, so this fails the makefile+bin criterion
            fakeGood = localMakeFakeCmdstanRootAt(sandbox, 'cmdstan-good');

            setenv('PSYRAT_TEST_CMDSTAN_PATH', fakeBad);
            setenv('STAN_HOME', fakeGood);
            setenv('CMDSTAN_HOME', '');
            setenv('CMDSTAN', '');

            home = psyrat_locate_cmdstan();
            testCase.verifyEqual(home, fakeGood, ...
                ['the invalid tier-1 root must be skipped and the valid ' ...
                'tier-2 root found']);
        end

        function testManagedGlobPrefersNewestVersion(testCase)
            % Tier 2 + the G43 ruling: with the env tier cleared and
            % user.home redirected at a sandbox holding cmdstan-2.38.0 and
            % cmdstan-2.40.0 (both valid), discovery must return 2.40.0.
            % 2.38/2.40 is the differential pair on purpose: dir()'s
            % lexicographic order puts 2.38.0 FIRST ('3' < '4'), so the old
            % first-match behavior returns the older install and this test
            % fails against it. (The classic 2.9-vs-2.38 pair is NOT
            % differential here: lexicographically '3' < '9' already put
            % 2.38.0 first.)
            sandbox = localSandbox(testCase);
            older = localMakeFakeCmdstanRootAt(sandbox, 'cmdstan-2.38.0');
            newer = localMakeFakeCmdstanRootAt(sandbox, 'cmdstan-2.40.0');
            testCase.assertTrue(isfolder(older) && isfolder(newer));

            for k = 1:numel(testCase.EnvNames)
                setenv(testCase.EnvNames{k}, '');
            end
            origHome = char(java.lang.System.getProperty('user.home'));
            java.lang.System.setProperty('user.home', sandbox);
            restoreHome = onCleanup(@() ...
                java.lang.System.setProperty('user.home', origHome)); %#ok<NASGU>

            home = psyrat_locate_cmdstan();
            testCase.verifyEqual(home, newer, ...
                'the managed glob tier must resolve the newest parseable version');
        end

        function testDetectVersionSeesAnEnvOnlyInstall(testCase)
            % The G37-deferred trap, closed: an install named only by an env
            % var (nothing on the MATLAB path) must be found AND its version
            % read. The fake stanc prints a sentinel version no real install
            % carries, so the reported string is attributable to the fake.
            testCase.assumeTrue(isunix, ...
                'the fake stanc is a shell script; version read is POSIX-only here');
            fake = localMakeFakeCmdstanRoot(testCase, 'cmdstan-envdetect');
            localAddFakeStanc(fake, '9.99.0');
            setenv('PSYRAT_TEST_CMDSTAN_PATH', fake);

            [verstr, home, ranok] = psyrat_detect_cmdstan_version();
            testCase.verifyEqual(home, fake);
            testCase.verifyTrue(ranok, 'the fake stanc must have executed');
            testCase.verifyEqual(verstr, '9.99.0');
        end

        function testCheckversionsofdepsQuietOnAnEnvOnlyInstall(testCase)
            % End to end through the startup version check: pre-G43 this
            % setup took the path-only branch, found nothing, and warned
            % "Unable to determine CmdStan version" while comparing against
            % a guessed version. Now it must read the env-located stanc and
            % classify 9.99.0 as an untested-newer version (code 2), with no
            % warning. The bundled MatlabStan/MPM are put on the path for the
            % non-CmdStan halves of the check; the base class restores the
            % path afterward.
            testCase.assumeTrue(isunix, ...
                'the fake stanc is a shell script; version read is POSIX-only here');
            root = testCase.projectRoot();
            addpath(genpath(fullfile(root, 'bundled_dependents')));

            fake = localMakeFakeCmdstanRoot(testCase, 'cmdstan-envcheck');
            localAddFakeStanc(fake, '9.99.0');
            setenv('PSYRAT_TEST_CMDSTAN_PATH', fake);

            versions = testCase.verifyWarningFree(@() psyrat_checkversionsofdeps);
            testCase.verifyEqual(versions.cmdstan, 2, ...
                ['the env-located 9.99.0 must be read and classified as ' ...
                'untested-newer, not guessed at with a warning']);
        end

        function testInjectionHelperConsultsTheResolverFirst(testCase)
            % Source-text weld (the TestStartprocPreflightAgreement census
            % pattern). psyrat_find_cmdstan_home -- the FIFTH divergent
            % order, found during the G43 fix itself -- feeds the stan_home
            % argument injected at the stan() call, and it is a local
            % function no headless route reaches, so the wiring is pinned in
            % the source: its body must consult psyrat_locate_cmdstan BEFORE
            % its legacy env tier (which skipped STAN_HOME entirely, so a
            % STAN_HOME-only setup could pass the centralized gate and then
            % fit against a different install). Without this weld the fifth
            % order could quietly return.
            src = fileread(fullfile(testCase.projectRoot(), ...
                'subroutines', 'estimation', 'psyrat_computevarcomp.m'));
            body = extractAfter(src, ...
                'function stan_home = psyrat_find_cmdstan_home');
            testCase.assertFalse(isempty(body), ...
                'the injection helper must still exist under this name');
            iLocate = strfind(body, 'stan_home = psyrat_locate_cmdstan();');
            iLegacy = strfind(body, ...
                'envnames = {''PSYRAT_TEST_CMDSTAN_PATH'',''CMDSTAN_HOME'',''CMDSTAN''}');
            testCase.assertNotEmpty(iLocate, ...
                'the helper must call the central resolver');
            testCase.assertNotEmpty(iLegacy, ...
                'anti-vacuity: the legacy tier must still exist as the fallback');
            testCase.verifyLessThan(iLocate(1), iLegacy(1), ...
                'the resolver consult must precede the legacy tiers');
        end

        function testInstallerAndSweepGatesConsultTheResolver(testCase)
            % G51 source-text weld (same pattern as the injection-helper
            % test above). The installer's three per-OS CmdStan gates live
            % in GUI-coupled local functions and the workflow sweep's env
            % seeding sits in that function's body ahead of any test run,
            % so no headless route executes them; the wiring is pinned in
            % the source instead. Anchored PER FUNCTION (the adversarial
            % wave showed an unanchored global count could pass with a
            % deleted gate and a duplicated one). Each gate must (a)
            % consult psyrat_locate_cmdstan BEFORE its path-only exist()
            % fallback, so an env-var-only install the startup gate accepts
            % (G43) is not reported as "needs to be installed", and (b)
            % record the resolved home as cmdstandir, so the downstream
            % make-build repair and MatlabStan pointing act on the SAME
            % directory the verdict was about instead of falling back to a
            % which() probe that is empty for off-path installs.
            src = fileread(fullfile(testCase.projectRoot(), ...
                'subroutines', 'installation', 'psyrat_installdependents.m'));
            fns = {'psyrat_macdepsinstall', 'psyrat_windepsinstall', ...
                'psyrat_linuxdepsinstall'};
            for k = 1:numel(fns)
                seg = extractAfter(src, ['function ' fns{k}]);
                testCase.assertFalse(isempty(seg), sprintf( ...
                    '%s must still exist under this name', fns{k}));
                nextf = strfind(seg, [newline 'function ']);
                if ~isempty(nextf)
                    seg = seg(1:nextf(1));
                end
                iLocate = strfind(seg, ...
                    'cmdstanhome = psyrat_locate_cmdstan();');
                iUse = strfind(seg, 'cmdstandir = cmdstanhome;');
                iFallback = strfind(seg, ...
                    'elseif exist(''makefile'',''file'') ~= 2');
                testCase.assertNotEmpty(iLocate, sprintf( ...
                    '%s: the gate must consult the resolver', fns{k}));
                testCase.assertNotEmpty(iUse, sprintf( ...
                    ['%s: the resolved home must be recorded as cmdstandir ' ...
                    'for the downstream build/pointing blocks'], fns{k}));
                testCase.assertNotEmpty(iFallback, sprintf( ...
                    '%s: anti-vacuity - the path-only fallback must remain', ...
                    fns{k}));
                testCase.verifyLessThan(iLocate(1), iFallback(1), sprintf( ...
                    '%s: the resolver consult must precede its fallback', ...
                    fns{k}));
            end
            % the one per-OS difference that matters: Windows keys built-ness
            % on the .exe binaries, the POSIX gates on the bare names
            winseg = extractAfter(src, 'function psyrat_windepsinstall');
            testCase.verifyNotEmpty(strfind(winseg, ...
                'fullfile(cmdstanhome,''bin'',''stansummary.exe'')'), ...
                'the Windows gate must check the .exe binaries');
            macseg = extractBefore( ...
                extractAfter(src, 'function psyrat_macdepsinstall'), ...
                'function psyrat_windepsinstall');
            testCase.verifyNotEmpty(strfind(macseg, ...
                'fullfile(cmdstanhome,''bin'',''stansummary'')'), ...
                'the macOS gate must check the bare-name binaries');

            % the sweep's seeding block: resolver in, path-only probe out
            sweep = fileread(fullfile(testCase.projectRoot(), ...
                'runExhaustiveWorkflowBugSweep.m'));
            testCase.verifyNotEmpty(strfind(sweep, ...
                'cmdstanHome = psyrat_locate_cmdstan();'), ...
                'the sweep must seed PSYRAT_TEST_CMDSTAN_PATH from the resolver');
            testCase.verifyEmpty(strfind(sweep, ...
                'cmdstanTestPy = which(''runCmdStanTests.py'')'), ...
                'the sweep''s path-only CmdStan probe must be gone');
        end

        function testUpdateFlowHardeningWelds(testCase)
            % G56 partial hardening, pinned as source welds like the G51
            % test above (the flow is GUI-coupled; no headless route runs
            % it). Facts pinned PER OS ARM (the dispatcher slice between
            % the sys == N anchors) -- the first draft used global counts
            % with pairwise ordering, which cannot prove per-arm
            % containment. (a) Each update arm derives wrkdir, tests
            % isempty(wrkdir) -- the CONDITION is welded, not just the
            % exit, so an inverted guard fails here -- and exits by the
            % B20 convention (relaunch, then the modal dialog) BEFORE
            % rmolddeps, so a doomed update cannot delete working
            % dependents first. (b) Each per-OS installer receives the
            % version verdict and carries the surviving-install notice,
            % condition AND payload (the single-fileparts survivor
            % derivation and the message fragment), between the resolver
            % gate and the download block it explains -- an emptied
            % notice body fails here.
            src = fileread(fullfile(testCase.projectRoot(), ...
                'subroutines', 'installation', 'psyrat_installdependents.m'));

            % (a) dispatcher arms, sliced per OS
            main = extractBefore(src, [newline 'function ']);
            iWin = strfind(main, 'elseif sys == 2');
            iLin = strfind(main, 'elseif sys == 3');
            testCase.assertTrue(isscalar(iWin) && isscalar(iLin) ...
                && iWin < iLin, ...
                'the per-OS dispatcher anchors must be unique and ordered');
            arms = {main(1:iWin-1), main(iWin:iLin-1), main(iLin:end)};
            for k = 1:3
                arm = arms{k};
                iDer = strfind(arm, ...
                    'wrkdir = fileparts(fileparts(which(''runCmdStanTests.py'')));');
                iCond = strfind(arm, 'if isempty(wrkdir)');
                iStart = strfind(arm, 'psyrat_start;');
                iDlg = strfind(arm, ...
                    'errordlg(dlgtext,''Update location not found'',''modal'');');
                iRm = strfind(arm, 'rmolddeps(depvercheck);');
                iCall = strfind(arm, 'depsinstall(wrkdir,urls,depvercheck);');
                testCase.assertEqual(numel(iDer), 1, sprintf( ...
                    'arm %d: exactly one wrkdir derivation', k));
                testCase.assertEqual(numel(iCond), 1, sprintf( ...
                    'arm %d: the empty-wrkdir CONDITION must be welded', k));
                testCase.assertEqual(numel(iDlg), 1, sprintf( ...
                    'arm %d: exactly one update-guard dialog', k));
                testCase.assertEqual(numel(iRm), 1, sprintf( ...
                    'arm %d: rmolddeps must remain, once', k));
                testCase.assertEqual(numel(iCall), 1, sprintf( ...
                    'arm %d: one per-OS installer call with the verdict', k));
                testCase.verifyLessThan(iDer, iCond, sprintf( ...
                    'arm %d: the guard must test the freshly derived wrkdir', k));
                testCase.verifyTrue(any(iStart > iCond & iStart < iDlg), ...
                    sprintf(['arm %d: the relaunch must sit inside the ' ...
                    'guard, before the dialog (B20 creation order)'], k));
                testCase.verifyLessThan(iDlg, iRm, sprintf( ...
                    ['arm %d: the guard must exit before rmolddeps so a ' ...
                    'doomed update cannot delete working dependents first'], k));
                testCase.verifyLessThan(iRm, iCall, sprintf( ...
                    'arm %d: removal precedes the per-OS install', k));
            end

            % (b) per-OS notice: verdict in the signature, notice condition
            % AND payload between the resolver gate and the download block
            fns = {'psyrat_macdepsinstall', 'psyrat_windepsinstall', ...
                'psyrat_linuxdepsinstall'};
            for k = 1:numel(fns)
                seg = extractAfter(src, ['function ' fns{k}]);
                testCase.assertFalse(isempty(seg), sprintf( ...
                    '%s must still exist under this name', fns{k}));
                nextf = strfind(seg, [newline 'function ']);
                if ~isempty(nextf)
                    seg = seg(1:nextf(1));
                end
                testCase.assertTrue(startsWith(seg, '(wrkdir,urls,depvercheck)'), ...
                    sprintf('%s must accept the version verdict', fns{k}));
                iLocate = strfind(seg, ...
                    'cmdstanhome = psyrat_locate_cmdstan();');
                iNote = strfind(seg, ...
                    ['if ~isempty(depvercheck) && depvercheck.cmdstan == 0' ...
                    ' && depcheck.cmdstan == 1']);
                iSurv = strfind(seg, ...
                    'survivor = fileparts(which(''runCmdStanTests.py''));');
                iMsg = strfind(seg, ...
                    'survived the removal step, so no new CmdStan will be downloaded.');
                iDl = strfind(seg, 'if depcheck.cmdstan == 0');
                testCase.assertNotEmpty(iLocate, sprintf( ...
                    '%s: the resolver gate must remain', fns{k}));
                testCase.assertNotEmpty(iNote, sprintf( ...
                    '%s: the verdict-keyed notice condition must exist', fns{k}));
                testCase.assertNotEmpty(iSurv, sprintf( ...
                    ['%s: the notice must derive the survivor with a SINGLE ' ...
                    'fileparts (the double spelling names the parent ' ...
                    'dependents directory in a message inviting deletion)'], ...
                    fns{k}));
                testCase.assertNotEmpty(iMsg, sprintf( ...
                    '%s: the notice payload must state the skipped download', ...
                    fns{k}));
                testCase.assertNotEmpty(iDl, sprintf( ...
                    '%s: anti-vacuity - the download block must remain', ...
                    fns{k}));
                testCase.verifyLessThan(iLocate(1), iNote(1), sprintf( ...
                    ['%s: the notice must sit after the gate whose verdict ' ...
                    'it reports'], fns{k}));
                testCase.verifyTrue(iSurv(1) > iNote(1) && iMsg(1) > iNote(1) ...
                    && iMsg(1) < iDl(1), sprintf( ...
                    ['%s: the payload must live inside the notice, before ' ...
                    'the skip it explains'], fns{k}));
            end
        end

        function testEnvpinUpdateGateWelds(testCase)
            % G56c, pinned as source welds like the G56 test above (the
            % flow is GUI-coupled and the gate is a local function; no
            % headless route runs it). Per OS arm: exactly one gate call,
            % sitting after the empty-wrkdir guard's dialog and BEFORE
            % rmolddeps, so a cancel costs nothing. Gate body: the
            % detection must mirror what makes the download inert -- the
            % resolver's winner matching an env var from the resolver's
            % own env tier -- and the cancel path must relaunch
            % psyrat_start (the update callback already closed every
            % PsyRAT window).
            src = fileread(fullfile(testCase.projectRoot(), ...
                'subroutines', 'installation', 'psyrat_installdependents.m'));

            % per-arm placement, sliced like the G56 weld above
            main = extractBefore(src, [newline 'function ']);
            iWin = strfind(main, 'elseif sys == 2');
            iLin = strfind(main, 'elseif sys == 3');
            testCase.assertTrue(isscalar(iWin) && isscalar(iLin) ...
                && iWin < iLin, ...
                'the per-OS dispatcher anchors must be unique and ordered');
            arms = {main(1:iWin-1), main(iWin:iLin-1), main(iLin:end)};
            for k = 1:3
                arm = arms{k};
                iGate = strfind(arm, 'if ~local_envpin_update_gate(depvercheck)');
                iDlg = strfind(arm, ...
                    'errordlg(dlgtext,''Update location not found'',''modal'');');
                iRm = strfind(arm, 'rmolddeps(depvercheck);');
                testCase.assertEqual(numel(iGate), 1, sprintf( ...
                    'arm %d: exactly one env-pin gate call', k));
                testCase.assertTrue(isscalar(iDlg) && isscalar(iRm), ...
                    sprintf('arm %d: weld anchors must be unique', k));
                testCase.verifyTrue(iGate > iDlg && iGate < iRm, sprintf( ...
                    ['arm %d: the gate must run after the empty-wrkdir ' ...
                    'guard and before rmolddeps'], k));
            end

            % gate body: condition AND payload
            seg = extractAfter(src, ...
                'function proceed = local_envpin_update_gate(depvercheck)');
            testCase.assertFalse(isempty(seg), ...
                'the gate must exist under this name and signature');
            nextf = strfind(seg, [newline 'function ']);
            testCase.assertNotEmpty(nextf, ...
                'the gate must remain a delimited local function');
            seg = seg(1:nextf(1));
            testCase.assertNotEmpty(strfind(seg, ...
                'depvercheck.cmdstan ~= 0'), ...
                ['the gate must key on the SAME verdict rmolddeps reads ' ...
                'before touching CmdStan (wave: unconditional firing made ' ...
                'the dialog noise on MatlabStan/MPM-only updates)']);
            testCase.assertNotEmpty(strfind(seg, ...
                'home = psyrat_locate_cmdstan();'), ...
                'the gate must consult the shared resolver');
            testCase.assertNotEmpty(strfind(seg, ...
                ['{''PSYRAT_TEST_CMDSTAN_PATH'',''STAN_HOME'',' ...
                '''CMDSTAN_HOME'',''CMDSTAN''}']), ...
                ['the gate must test the resolver''s own env tier, in its ' ...
                'precedence order']);
            testCase.assertNotEmpty(strfind(seg, ...
                'val = local_envpin_normalize(getenv(pinvars{ii}));'), ...
                ['the env value must be normalized the way the resolver ' ...
                'normalizes candidates (wave: the raw comparison missed ' ...
                'pins spelled as <home>/makefile or <home>/bin)']);
            testCase.assertNotEmpty(strfind(seg, ...
                'if ~isempty(val) && local_samedir(val,home)'), ...
                ['the CONDITION must be welded: the winner matching a SET ' ...
                'env var, not the mere existence of one']);
            testCase.assertNotEmpty(strfind(seg, ...
                'rmtarget = fileparts(which(''runCmdStanTests.py''));'), ...
                ['the gate must exclude a pin naming rmolddeps'' own ' ...
                'removal target (wave: that pin dies with the removal, so ' ...
                'the update IS effective and the dialog''s premise false)']);
            testCase.assertNotEmpty(strfind(seg, ...
                'local_samedir(rmtarget,home)'), ...
                'the removal-target exclusion must compare canonical homes');
            testCase.assertNotEmpty(strfind(seg, ...
                'so this update will not change which CmdStan estimation uses'), ...
                'the dialog payload must state the inertness plainly');
            testCase.assertNotEmpty(strfind(seg, ...
                '''Update anyway'',''Cancel'',''Cancel'');'), ...
                'the safe choice must be the default button');
            iNo = strfind(seg, 'proceed = false;');
            iRelaunch = strfind(seg, 'psyrat_start;');
            testCase.assertTrue(isscalar(iNo) && isscalar(iRelaunch) ...
                && iNo < iRelaunch, ...
                ['the cancel path must refuse and then relaunch ' ...
                'psyrat_start, so a cancelling user is not left windowless']);
        end

    end
end

%--------------------------------------------------------------------------
function d = localSandbox(testCase)
%Fresh per-test temp dir, removed on teardown.
d = fullfile(tempdir, ['psyrat_g43_' char(java.util.UUID.randomUUID())]);
mkdir(d);
testCase.addTeardown(@() localRmdirIfExists(d));
end

function localRmdirIfExists(d)
if isfolder(d)
    rmdir(d, 's');
end
end

function p = localMakeFakeCmdstanRoot(testCase, name)
p = localMakeFakeCmdstanRootAt(localSandbox(testCase), name);
end

function p = localMakeFakeCmdstanRootAt(parent, name)
%A directory that satisfies stan_home's validity criterion: makefile + bin/.
p = fullfile(parent, name);
mkdir(p);
mkdir(fullfile(p, 'bin'));
fid = fopen(fullfile(p, 'makefile'), 'w');
fclose(fid);
end

function localAddFakeStanc(root, version)
%POSIX shell script standing in for bin/stanc, printing a version banner.
f = fullfile(root, 'bin', 'stanc');
fid = fopen(f, 'w');
fprintf(fid, '#!/bin/sh\necho "stanc version %s"\n', version);
fclose(fid);
fileattrib(f, '+x');
end
