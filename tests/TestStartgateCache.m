classdef TestStartgateCache < PsyRATTestBase
    %TESTSTARTGATECACHE The launcher reuses its dependency and release
    %verdicts within one MATLAB session (G91, owner ruling 2026-09-15).
    %
    % psyrat_start used to spawn `stanc --version` two or three times and ask
    % the GitHub releases API once on EVERY launch, and it is relaunched from
    % 23 bare `psyrat_start;` sites (viewer Back buttons, cancelled dialogs,
    % installer arms). The fix keeps those verdicts in a `persistent` store
    % inside subroutines/installation/psyrat_startgate_cache.m, keyed by a
    % cheap fingerprint of what the checks read (toolbox location and
    % version, the resolved CmdStan home and its stanc file, the MatlabStan
    % and MatlabProcessManager copies on the path), and the launcher skips
    % the three expensive calls on a hit.
    %
    % Three layers, each RED before the helper existed:
    %  1. the helper contract, executed headlessly with fake CmdStan roots and
    %     path shadows (no stanc, no network);
    %  2. the launcher end to end -- the first tests in the suite to EXECUTE
    %     psyrat_start -- with counting shadows on the three expensive calls
    %     (the version check, the minimum-version check, the release check)
    %     and a no-op psyrat_savepath_safe so no test can write pathdef.m;
    %  3. source welds on psyrat_start.m, because the wiring (one lookup
    %     before the version check, one store after the release check, one
    %     clear in each installer/updater callback) is what the runtime tests
    %     cannot distinguish from a coincidence.
    % Every forbid pin is paired with a positive control (an empty file
    % satisfies a prohibition).
    %
    % The leaf functions (psyrat_locate_cmdstan, psyrat_checkversionsofdeps,
    % psyrat_detect_cmdstan_version, psyrat_require_min_cmdstan,
    % psyrat_checkgithubreleases, psyrat_updatecheck) must stay stateless:
    % TestCmdstanDiscovery and TestInstallationAndWorkflowValidation call them
    % live under changing env/path state in the same process and expect fresh
    % answers, and the updater and estimation consult them live too. The last
    % test here pins that.

    properties (Access = private)
        SavedEnv
    end

    properties (Constant, Access = private)
        EnvNames = {'PSYRAT_TEST_CMDSTAN_PATH','STAN_HOME','CMDSTAN_HOME','CMDSTAN'};
        CounterKey = 'psyrat_test_gatecalls';
        HelperName = 'psyrat_startgate_cache';
    end

    methods (TestMethodSetup)
        function saveEnvAndClearStore(testCase)
            testCase.SavedEnv = cellfun(@getenv, testCase.EnvNames, ...
                'UniformOutput', false);
            %the resolver (mstan.stan_home) and the launcher's MatlabStan /
            %MatlabProcessManager exist() gate both need the bundled tree,
            %which PsyRATTestBase does not add; the base class restores the
            %path after each method.
            addpath(genpath(fullfile(testCase.projectRoot(), 'bundled_dependents')));
            %every test starts from an empty store. Guarded so that the RED
            %run (helper absent) fails on the assertions, not in setup.
            if exist([testCase.HelperName '.m'], 'file') == 2
                psyrat_startgate_cache('clear');
            end
            setappdata(0, testCase.CounterKey, struct('deps', 0, 'release', 0, 'detect', 0, 'minver', 0));
        end
    end

    methods (TestMethodTeardown)
        function restoreEnvAndClearStore(testCase)
            for k = 1:numel(testCase.EnvNames)
                setenv(testCase.EnvNames{k}, testCase.SavedEnv{k});
            end
            if exist([testCase.HelperName '.m'], 'file') == 2
                psyrat_startgate_cache('clear');
            end
            if isappdata(0, testCase.CounterKey)
                rmappdata(0, testCase.CounterKey);
            end
        end
    end

    methods (Test)

        % ---------------------------------------------------------------
        % 1. Helper contract
        % ---------------------------------------------------------------

        function testLookupMissesOnAFreshStore(testCase)
            testCase.pointResolverAtFakeRoot('cmdstan-fresh');
            [hit, cached, key, why] = psyrat_startgate_cache('lookup');
            testCase.verifyFalse(hit, 'an empty store must miss');
            testCase.verifyEmpty(cached);
            testCase.verifyTrue(isstruct(key) && ~isempty(fieldnames(key)), ...
                'lookup must still hand back the fingerprint it computed');
            testCase.verifyEqual(why, 'no earlier check this MATLAB session');
            testCase.verifyEmpty(psyrat_startgate_cache('get'));
        end

        function testStoreThenLookupHitsAndIsInspectable(testCase)
            testCase.pointResolverAtFakeRoot('cmdstan-hit');
            dv = struct('cmdstan', 1, 'matlabstan', 3, 'matlabprocessmanager', 1);
            [~, ~, key] = psyrat_startgate_cache('lookup');
            psyrat_startgate_cache('store', key, dv, true, 1, 'current_stable');

            [hit, cached, key2, why] = psyrat_startgate_cache('lookup');
            testCase.verifyTrue(hit, 'the same machine state must hit');
            testCase.verifyEqual(cached.depvercheck, dv);
            testCase.verifyTrue(cached.mincmdstan_ok);
            testCase.verifyEqual(cached.psyratvercheck, 1);
            testCase.verifyEqual(cached.release_status, 'current_stable');
            testCase.verifyEqual(key2, key);
            testCase.verifyEqual(why, '');

            rec = psyrat_startgate_cache('get');
            testCase.verifyTrue(isstruct(rec), 'get must expose the record for inspection');
            testCase.verifyEqual(rec.depvercheck, dv);
            testCase.verifyEqual(rec.key, key);
            testCase.verifyTrue(ischar(rec.checked_at) && ...
                ~isempty(regexp(rec.checked_at, '^\d\d:\d\d:\d\d$', 'once')), ...
                'checked_at is the wall-clock time of the fresh check, HH:mm:ss');
        end

        function testKeyFollowsTheResolvedCmdstanHome(testCase)
            fakeA = testCase.pointResolverAtFakeRoot('cmdstan-a');
            testCase.storeAllCurrent();
            testCase.assertTrue(psyrat_startgate_cache('lookup'), 'control: same root hits');

            fakeB = localMakeFakeCmdstanRootAt(testCase.sandbox(), 'cmdstan-b');
            setenv('PSYRAT_TEST_CMDSTAN_PATH', fakeB);
            [hit, ~, ~, why] = psyrat_startgate_cache('lookup');
            testCase.verifyFalse(hit, 'a different resolved CmdStan home must miss');
            testCase.verifyTrue(contains(why, 'CmdStan location'), ...
                sprintf('why must name the CmdStan location, got "%s"', why));

            %anti-vacuity: the miss is about the root, not a permanent failure
            setenv('PSYRAT_TEST_CMDSTAN_PATH', fakeA);
            testCase.storeAllCurrent();
            testCase.verifyTrue(psyrat_startgate_cache('lookup'), ...
                'after a fresh store under root A the lookup must hit again');
        end

        function testKeyFollowsTheStancBytes(testCase)
            fake = testCase.pointResolverAtFakeRoot('cmdstan-stanc');
            stanc = fullfile(fake, 'bin', 'stanc');
            localWriteText(stanc, 'stanc');
            testCase.storeAllCurrent();
            testCase.assertTrue(psyrat_startgate_cache('lookup'), 'control: unchanged stanc hits');

            %a different byte length, so the outcome does not depend on the
            %file system's mtime resolution
            localWriteText(stanc, 'stanc-rebuilt');
            [hit, ~, ~, why] = psyrat_startgate_cache('lookup');
            testCase.verifyFalse(hit, 'an in-place stanc rebuild must miss');
            testCase.verifyTrue(contains(why, 'CmdStan compiler'), ...
                sprintf('why must name the compiler, got "%s"', why));
        end

        function testKeyFollowsAMatlabStanShadow(testCase)
            testCase.pointResolverAtFakeRoot('cmdstan-ms');
            testCase.storeAllCurrent();
            testCase.assertTrue(psyrat_startgate_cache('lookup'), 'control: hits before the shadow');

            depvers = psyrat_dependentsversions;
            testCase.addShadow({ ...
                'StanModel.m', {'function m = StanModel(varargin)', 'm = struct();', 'end'}; ...
                'PSYRAT_PATCHLEVEL.txt', {'PsyRAT vendored-bundle stamp', '', ...
                    ['    upstream_version = ' depvers.matlabstan], ...
                    ['    patch_level = ' depvers.matlabstanpatch]}});
            [hit, ~, ~, why] = psyrat_startgate_cache('lookup');
            testCase.verifyFalse(hit, 'a different MatlabStan first on the path must miss');
            testCase.verifyTrue(contains(why, 'MatlabStan'), ...
                sprintf('why must name MatlabStan, got "%s"', why));
        end

        function testKeyFollowsAProcessManagerShadow(testCase)
            testCase.pointResolverAtFakeRoot('cmdstan-mpm');
            testCase.storeAllCurrent();
            testCase.assertTrue(psyrat_startgate_cache('lookup'), 'control: hits before the shadow');

            testCase.addShadow({'processManager.m', ...
                {'function p = processManager(varargin)', 'version = ''0.5.1''; %#ok<NASGU>', ...
                 'p = struct();', 'end'}});
            [hit, ~, ~, why] = psyrat_startgate_cache('lookup');
            testCase.verifyFalse(hit, 'a different MatlabProcessManager first on the path must miss');
            testCase.verifyTrue(contains(why, 'MatlabProcessManager'), ...
                sprintf('why must name MatlabProcessManager, got "%s"', why));
        end

        function testMismatchDropsTheStaleEntry(testCase)
            fakeA = testCase.pointResolverAtFakeRoot('cmdstan-drop-a');
            testCase.storeAllCurrent();
            fakeB = localMakeFakeCmdstanRootAt(testCase.sandbox(), 'cmdstan-drop-b');
            setenv('PSYRAT_TEST_CMDSTAN_PATH', fakeB);
            testCase.assertFalse(psyrat_startgate_cache('lookup'), 'control: B misses');

            setenv('PSYRAT_TEST_CMDSTAN_PATH', fakeA);
            [hit, ~, ~, why] = psyrat_startgate_cache('lookup');
            testCase.verifyFalse(hit, ...
                'the entry stored under A must have been DROPPED by the B miss, not shelved');
            testCase.verifyEqual(why, 'no earlier check this MATLAB session');
            testCase.verifyEmpty(psyrat_startgate_cache('get'));
        end

        function testClearForgets(testCase)
            testCase.pointResolverAtFakeRoot('cmdstan-clear');
            testCase.storeAllCurrent();
            testCase.assertTrue(psyrat_startgate_cache('lookup'), 'control: stored');
            psyrat_startgate_cache('clear');
            testCase.verifyFalse(psyrat_startgate_cache('lookup'), 'the clear verb must empty the store');
            testCase.verifyEmpty(psyrat_startgate_cache('get'));
            %clearing an empty store is not an error
            psyrat_startgate_cache('clear');

            %the console line tells users a plain MATLAB clear works too
            testCase.storeAllCurrent();
            testCase.assertTrue(psyrat_startgate_cache('lookup'), 'control: stored again');
            clear psyrat_startgate_cache
            testCase.verifyFalse(psyrat_startgate_cache('lookup'), ...
                'clear <helper> must reset the persistent store');
        end

        function testStoreRefusesAnEmptyKey(testCase)
            testCase.pointResolverAtFakeRoot('cmdstan-nokey');
            dv = struct('cmdstan', 1, 'matlabstan', 1, 'matlabprocessmanager', 1);
            psyrat_startgate_cache('store', [], dv, true, 1, 'current_stable');
            testCase.verifyEmpty(psyrat_startgate_cache('get'), ...
                'an empty key means the fingerprint failed: never store against it');
            testCase.verifyFalse(psyrat_startgate_cache('lookup'));
            %positive control: the same store with a real key is kept
            testCase.storeAllCurrent();
            testCase.verifyTrue(psyrat_startgate_cache('lookup'));
        end

        function testBadActionAndBadInputError(testCase)
            testCase.verifyError(@() psyrat_startgate_cache('bogus'), ...
                'psyrat_startgate_cache:badAction');
            testCase.verifyError(@() psyrat_startgate_cache(), ...
                'psyrat_startgate_cache:badAction');
            testCase.verifyError(@() psyrat_startgate_cache('store', struct('a', 1)), ...
                'psyrat_startgate_cache:badInput');
        end

        function testKeyBuildsWithoutAnyCmdstan(testCase)
            %no env var and a user.home with no managed install: the
            %fingerprint must still be a struct, with the CmdStan fields
            %reporting absence rather than the key throwing or coming back
            %empty (an empty key would disable the cache silently).
            for k = 1:numel(testCase.EnvNames)
                setenv(testCase.EnvNames{k}, '');
            end
            sandbox = testCase.sandbox();
            origHome = char(java.lang.System.getProperty('user.home'));
            java.lang.System.setProperty('user.home', sandbox);
            testCase.addTeardown(@() java.lang.System.setProperty('user.home', origHome));

            [hit, ~, key, why] = psyrat_startgate_cache('lookup');
            testCase.verifyFalse(hit);
            testCase.verifyTrue(isstruct(key), 'the key must be a struct even with no CmdStan');
            testCase.verifyTrue(all(isfield(key, {'toolboxroot', 'psyratver', ...
                'cmdstanhome', 'stanc', 'matlabstan', 'stamp', 'resolver', 'mpm', 'mpmstamp'})), ...
                'the key must carry every named fingerprint field');
            testCase.verifyEqual(why, 'no earlier check this MATLAB session');
        end

        function testKeyBuildingLeavesTheRngUntouched(testCase)
            testCase.pointResolverAtFakeRoot('cmdstan-rng');
            s0 = rng;
            psyrat_startgate_cache('lookup');
            testCase.storeAllCurrent();
            psyrat_startgate_cache('lookup');
            s1 = rng;
            testCase.verifyEqual(s1.Type, s0.Type);
            testCase.verifyEqual(s1.State, s0.State, ...
                'the cache must never draw from the RNG (bit-reproducibility guard)');
        end

        function testKeyBuildingSpawnsNothing(testCase)
            %every stanc spawn in the gate goes through one of these two;
            %the fingerprint must reach neither.
            testCase.pointResolverAtFakeRoot('cmdstan-spawn');
            testCase.addShadow({ ...
                'psyrat_detect_cmdstan_version.m', {'function [verstr, home, ranok] = psyrat_detect_cmdstan_version()', ...
                    'n = getappdata(0, ''psyrat_test_gatecalls''); n.detect = n.detect + 1;', ...
                    'setappdata(0, ''psyrat_test_gatecalls'', n);', ...
                    'verstr = ''''; home = ''''; ranok = false;', 'end'}; ...
                'psyrat_checkversionsofdeps.m', localDepsShadowLines([1 1 1])});
            %positive control: the shadows are live
            psyrat_detect_cmdstan_version();
            psyrat_checkversionsofdeps();
            n = getappdata(0, testCase.CounterKey);
            testCase.assertEqual([n.detect n.deps], [1 1], 'control: shadows must count');

            psyrat_startgate_cache('lookup');
            testCase.storeAllCurrent();
            psyrat_startgate_cache('lookup');
            n = getappdata(0, testCase.CounterKey);
            testCase.verifyEqual([n.detect n.deps], [1 1], ...
                'building the key must not call the version readers (no stanc spawn)');
        end

        % ---------------------------------------------------------------
        % 2. The launcher end to end (no CmdStan, no network)
        % ---------------------------------------------------------------

        function testRelaunchReusesTheChecksWithinASession(testCase)
            testCase.installLauncherShadows([1 1 1], 1, 'current_stable');
            out1 = testCase.launch();
            testCase.assertNotEmpty(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui'), ...
                'control: the launcher must have built the home screen');
            testCase.assertFalse(contains(out1, 'reused from this MATLAB session'), ...
                'the first launch of a session is a fresh check');
            testCase.assertTrue(contains(out1, 'CmdStan version is current'), ...
                'control: the verdict lines print on the fresh launch');

            out2 = testCase.launch();
            n = getappdata(0, testCase.CounterKey);
            testCase.verifyEqual([n.deps n.release n.minver], [1 1 1], ...
                ['the second launch must not call the version check, the minimum-version ' ...
                'check or the release check again']);
            testCase.verifyTrue(contains(out2, 'reused from this MATLAB session'), ...
                'the relaunch must say on the console that it reused the checks');
            testCase.verifyTrue(contains(out2, 'psyrat_startgate_cache(''clear'')'), ...
                'the reuse line must name the recheck command');
            testCase.verifyTrue(contains(out2, 'CmdStan version is current'), ...
                'the three verdict lines must still print from the cached struct');
            testCase.verifyTrue(contains(out2, 'current_stable'), ...
                'the reuse line must name the release status it is reusing');
        end

        function testRelaunchRechecksAfterTheResolverMoves(testCase)
            testCase.installLauncherShadows([1 1 1], 1, 'current_stable');
            testCase.launch();
            fakeB = localMakeFakeCmdstanRootAt(testCase.sandbox(), 'cmdstan-moved');
            setenv('PSYRAT_TEST_CMDSTAN_PATH', fakeB);
            out2 = testCase.launch();
            n = getappdata(0, testCase.CounterKey);
            testCase.verifyEqual([n.deps n.release n.minver], [2 2 2], ...
                'a changed CmdStan location must force a fresh check');
            testCase.verifyTrue(contains(out2, 'checks rerun:') && contains(out2, 'CmdStan location'), ...
                'the relaunch must say why it rechecked');
        end

        function testClearForcesAFreshCheckOnRelaunch(testCase)
            testCase.installLauncherShadows([1 1 1], 1, 'current_stable');
            testCase.launch();
            psyrat_startgate_cache('clear');
            out2 = testCase.launch();
            n = getappdata(0, testCase.CounterKey);
            testCase.verifyEqual([n.deps n.release n.minver], [2 2 2], ...
                'after the clear verb the relaunch must recheck');
            testCase.verifyFalse(contains(out2, 'reused from this MATLAB session'));
        end

        function testGateFailureIsNeverStored(testCase)
            testCase.installLauncherShadows([1 1 1], 1, 'current_stable');
            testCase.addShadow({'psyrat_locate_cmdstan.m', ...
                {'function home = psyrat_locate_cmdstan()', 'home = '''';', 'end'}});
            out = testCase.launch();
            testCase.assertTrue(contains(out, 'Dependencies for PsyRAT were not found'), ...
                'control: the launcher must have taken the not-found path');
            n = getappdata(0, testCase.CounterKey);
            testCase.verifyEqual([n.deps n.release n.minver], [0 0 0], ...
                'the not-found path runs before the version and release checks');
            testCase.verifyEmpty(psyrat_startgate_cache('get'), ...
                'a failed gate must never be stored');
        end

        function testOfflineReleaseVerdictIsReused(testCase)
            %ruling 3 (2026-09-15): an "unable to check" release outcome is
            %reused like any other, so an offline or firewalled user pays the
            %timeout once per session rather than on every Back click.
            testCase.installLauncherShadows([1 1 1], -1, 'unable_to_check');
            testCase.launch();
            out2 = testCase.launch();
            n = getappdata(0, testCase.CounterKey);
            testCase.verifyEqual([n.deps n.release], [1 1], ...
                'a -1 release verdict must be reused, not retried');
            testCase.verifyTrue(contains(out2, 'unable_to_check'), ...
                'the reuse line must show that the reused release check could not complete');
            testCase.verifyEmpty(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui_update'), ...
                'control: no update prompt from a -1 verdict');
        end

        function testCachedVerdictStillRaisesTheUpdatePrompt(testCase)
            %ruling 1 (2026-09-15): the modal update prompts keep firing on
            %every relaunch (B44/G87). The cached verdict must reach the
            %prompt gate exactly as a fresh one does.
            testCase.installLauncherShadows([1 1 1], 0, 'update_available');
            testCase.launch();
            testCase.assertNotEmpty(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui_update'), ...
                'control: a fresh 0 verdict raises the toolbox update prompt');
            close all force
            testCase.launch();
            n = getappdata(0, testCase.CounterKey);
            testCase.verifyEqual([n.deps n.release], [1 1], 'control: the relaunch reused the checks');
            testCase.verifyNotEmpty(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui_update'), ...
                'the cached 0 verdict must raise the prompt again (ruling 1: no silent suppression)');
        end

        function testCachedZeroVerdictStillRaisesTheDepsPrompt(testCase)
            %ruling 1 (2026-09-15), the dependents half: a reused 0 dependent
            %verdict must raise the dependents prompt exactly as a fresh one
            %does. Seen live on 2026-09-15 (Item 10 check D: the prompt fired
            %again from the STORED zeros with no checker call) and pinned here
            %because no suite test fed zero dependent codes through the store.
            %Three zeros, as the live check's stub returned, and a current
            %release, because the toolbox prompt takes the if-branch of the
            %prompt gate and the dependents prompt is its elseif. Both prompts
            %share the Tag psyrat_gui_update (B45), so the prompt TEXT is what
            %tells them apart here.
            %This is a regression pin: it was green when written. Its
            %load-bearing status was proved by a scratch mutation that
            %suppressed the dependents prompt on the reuse branch, which
            %turned this test, and only this test, red (2026-09-19).
            testCase.installLauncherShadows([0 0 0], 1, 'current_stable');
            testCase.launch();
            [txt, nfig] = localDepsPromptText();
            testCase.assertEqual(nfig, 1, 'control: exactly one update prompt figure after the fresh launch');
            testCase.assertTrue(contains(txt, 'update the dependents'), ...
                'control: a fresh 0 dependent verdict raises the DEPENDENTS prompt');
            testCase.assertTrue(contains(txt, 'CmdStan, MatlabStan, MatlabProcessManager'), ...
                'control: the fresh prompt names every dependent flagged 0 (B42)');
            close all force
            testCase.launch();
            n = getappdata(0, testCase.CounterKey);
            testCase.verifyEqual([n.deps n.release], [1 1], 'control: the relaunch reused the checks');
            [txt2, nfig2] = localDepsPromptText();
            testCase.verifyEqual(nfig2, 1, 'exactly one update prompt figure after the relaunch');
            testCase.verifyTrue(contains(txt2, 'update the dependents'), ...
                ['the cached 0 dependent verdict must raise the dependents prompt again ' ...
                '(ruling 1: no silent suppression)']);
            testCase.verifyTrue(contains(txt2, 'CmdStan, MatlabStan, MatlabProcessManager'), ...
                'the reused prompt must name every dependent from the stored verdict');
        end

        % ---------------------------------------------------------------
        % 3. Source welds
        % ---------------------------------------------------------------

        function testLauncherWelds(testCase)
            code = localCode(fullfile(testCase.projectRoot(), 'psyrat_start.m'));

            iLookup = strfind(code, 'psyrat_startgate_cache(''lookup'')');
            iDeps = strfind(code, 'depvercheck = psyrat_checkversionsofdeps;');
            testCase.assertEqual(numel(iDeps), 1, 'control: exactly one live version check');
            testCase.verifyEqual(numel(iLookup), 1, 'exactly one lookup in the launcher');
            testCase.verifyTrue(~isempty(iLookup) && iLookup(1) < iDeps(1), ...
                'the lookup must precede the version check it can skip');

            iStore = strfind(code, 'psyrat_startgate_cache(''store''');
            iRelease = strfind(code, 'psyrat_checkgithubreleases(');
            testCase.assertEqual(numel(iRelease), 1, 'control: exactly one release check call');
            testCase.verifyEqual(numel(iStore), 1, 'exactly one store in the launcher');
            testCase.verifyTrue(~isempty(iStore) && iStore(1) > iRelease(1), ...
                'the store must follow the release check, so only a completed gate is stored');

            %the first figure built in the file is the home screen (the
            %prompts' figures are built by local functions further down)
            iGui = strfind(code, 'psyrat_gui= psyrat_newfigure(');
            testCase.assertGreaterThanOrEqual(numel(iGui), 1, 'control: the home-screen figure is built');
            testCase.verifyTrue(~isempty(iStore) && iStore(1) < iGui(1), ...
                'the store must land before the GUI, so no callback can fire on an unstored launch');

            %one clear as the first statement of each callback that leads
            %into the installer or the updater (the fingerprint is the
            %primary guard; this is the belt and braces for the prompt whose
            %Yes deletes directories)
            blocks = strsplit(code, sprintf('\nfunction '));
            for name = {'updatepsyrat', 'updatedeps', 'installdeps'}
                hits = blocks(startsWith(blocks, [name{1} '(']));
                testCase.assertEqual(numel(hits), 1, ...
                    sprintf('control: exactly one local function %s', name{1}));
                testCase.verifyEqual(numel(strfind(hits{1}, 'psyrat_startgate_cache(''clear'');')), 1, ...
                    sprintf('%s must clear the session store exactly once', name{1}));
            end
            testCase.verifyEqual(numel(strfind(code, 'psyrat_startgate_cache(''clear'');')), 3, ...
                'no clear outside those three callbacks');

            testCase.verifyEqual(numel(strfind(code, 'exist(''psyrat_startgate_cache.m'',''file'') ~= 2')), 1, ...
                'the helper must be in the sentinel chain (a stale saved path must not die on it)');
            testCase.verifyEqual(numel(strfind(code, 'reused from this MATLAB session')), 1, ...
                'the reuse line is printed by exactly one statement');
        end

        function testLeafFunctionsCarryNoMemo(testCase)
            root = testCase.projectRoot();
            leaves = { ...
                fullfile(root, 'subroutines', 'installation', 'psyrat_locate_cmdstan.m'), ...
                fullfile(root, 'subroutines', 'installation', 'psyrat_checkversionsofdeps.m'), ...
                fullfile(root, 'subroutines', 'installation', 'psyrat_detect_cmdstan_version.m'), ...
                fullfile(root, 'subroutines', 'installation', 'psyrat_require_min_cmdstan.m'), ...
                fullfile(root, 'subroutines', 'installation', 'psyrat_checkgithubreleases.m'), ...
                fullfile(root, 'subroutines', 'installation', 'psyrat_updatecheck.m'), ...
                fullfile(root, 'bundled_dependents', 'MatlabStan-2.15.1.0', '+mstan', 'stan_home.m')};
            for k = 1:numel(leaves)
                code = localCode(leaves{k});
                testCase.assertGreaterThan(numel(code), 0, ...
                    sprintf('control: %s was read non-empty', leaves{k}));
                testCase.verifyEqual(numel(regexp(code, '\<persistent\>', 'start')), 0, ...
                    sprintf('%s must stay stateless', leaves{k}));
                testCase.verifyEqual(numel(strfind(code, 'appdata(')), 0, ...
                    sprintf('%s must not keep root appdata', leaves{k}));
                testCase.verifyEqual(numel(strfind(code, testCase.HelperName)), 0, ...
                    sprintf('%s must not consult the session store', leaves{k}));
            end
            %positive controls for the three forbid pins above: the persistent
            %detector fires on the helper, the appdata detector fires on this
            %very file (its counters live in root appdata), the name detector
            %fires on the launcher
            helper = localCode(fullfile(root, 'subroutines', 'installation', [testCase.HelperName '.m']));
            testCase.verifyEqual(numel(regexp(helper, '\<persistent\>', 'start')), 1, ...
                'the helper is where the one persistent store lives');
            self = localCode(fullfile(root, 'tests', 'TestStartgateCache.m'));
            testCase.verifyGreaterThan(numel(strfind(self, 'appdata(')), 0, ...
                'control: the appdata detector must fire on a file that uses it');
            launcher = localCode(fullfile(root, 'psyrat_start.m'));
            testCase.verifyGreaterThan(numel(strfind(launcher, testCase.HelperName)), 0, ...
                'the launcher is what consults the store');

            %estimation never touches the store: a fit must see the live install
            est = dir(fullfile(root, 'subroutines', 'estimation', '**', '*.m'));
            testCase.assertNotEmpty(est, 'control: the estimation tree was found');
            for k = 1:numel(est)
                src = fileread(fullfile(est(k).folder, est(k).name));
                testCase.verifyEqual(numel(strfind(src, testCase.HelperName)), 0, ...
                    sprintf('%s must not consult the session store', est(k).name));
            end
        end
    end

    methods (Access = private)

        function d = sandbox(testCase)
            d = tempname;
            mkdir(d);
            testCase.addTeardown(@() localRmdirIfExists(d));
        end

        function fake = pointResolverAtFakeRoot(testCase, name)
            %tier 1 of mstan.stan_home: a valid root named by the first env
            %name wins outright, which keeps every key machine-independent
            fake = localMakeFakeCmdstanRootAt(testCase.sandbox(), name);
            setenv('PSYRAT_TEST_CMDSTAN_PATH', fake);
        end

        function storeAllCurrent(~)
            [~, ~, key] = psyrat_startgate_cache('lookup');
            dv = struct('cmdstan', 1, 'matlabstan', 1, 'matlabprocessmanager', 1);
            psyrat_startgate_cache('store', key, dv, true, 1, 'current_stable');
        end

        function addShadow(testCase, files)
            %FILES is an n-by-2 cell: {name, cellstr of lines}. The directory
            %goes FIRST on the path so it shadows the toolbox and the bundle;
            %rehash so the function cache cannot serve the previous
            %definition (the repo pairs every shadow addpath with rehash).
            d = testCase.sandbox();
            for k = 1:size(files, 1)
                localWriteLines(fullfile(d, files{k, 1}), files{k, 2});
            end
            addpath(d, '-begin');
            rehash;
            testCase.addTeardown(@() localRmpathQuiet(d));
        end

        function installLauncherShadows(testCase, depcodes, statusCode, statusLabel)
            %the three expensive calls, counted; and psyrat_savepath_safe as a
            %no-op so no test can ever write pathdef.m (it is reached only
            %when a sentinel fails, which a shadow of the release check must
            %not trigger: exist() still finds a file by that name).
            testCase.pointResolverAtFakeRoot('cmdstan-launch');
            testCase.addShadow({ ...
                'psyrat_checkversionsofdeps.m', localDepsShadowLines(depcodes); ...
                'psyrat_require_min_cmdstan.m', { ...
                    'function ok = psyrat_require_min_cmdstan(varargin)', ...
                    'n = getappdata(0, ''psyrat_test_gatecalls''); n.minver = n.minver + 1;', ...
                    'setappdata(0, ''psyrat_test_gatecalls'', n);', ...
                    'ok = true;', 'end'}; ...
                'psyrat_checkgithubreleases.m', { ...
                    'function info = psyrat_checkgithubreleases(varargin)', ...
                    'n = getappdata(0, ''psyrat_test_gatecalls''); n.release = n.release + 1;', ...
                    'setappdata(0, ''psyrat_test_gatecalls'', n);', ...
                    sprintf('info = struct(''statusCode'', %d, ''status'', ''%s'', ''currentVersion'', '''', ''latestReleaseAssetURL'', '''');', statusCode, statusLabel), ...
                    'end'}; ...
                'psyrat_savepath_safe.m', {'function ok = psyrat_savepath_safe(varargin)', 'ok = true;', 'end'}});
            %the launcher must be the repository's own, not a worktree copy
            testCase.assertTrue(startsWith(which('psyrat_start'), testCase.projectRoot()), ...
                sprintf('control: psyrat_start resolves outside the checkout: %s', which('psyrat_start')));
        end

        function out = launch(~)
            %evalc captures the console; the launcher builds its figure(s)
            %and returns (no uiwait), so this is a plain function call
            out = evalc('psyrat_start');
        end
    end
end

function lines = localDepsShadowLines(codes)
lines = { ...
    'function versions = psyrat_checkversionsofdeps', ...
    'n = getappdata(0, ''psyrat_test_gatecalls''); n.deps = n.deps + 1;', ...
    'setappdata(0, ''psyrat_test_gatecalls'', n);', ...
    sprintf('versions = struct(''cmdstan'', %d, ''matlabstan'', %d, ''matlabprocessmanager'', %d);', codes(1), codes(2), codes(3)), ...
    'end'};
end

function localWriteLines(file, lines)
fid = fopen(file, 'w');
assert(fid > 0, 'could not create %s', file);
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end

function localWriteText(file, text)
fid = fopen(file, 'w');
assert(fid > 0, 'could not create %s', file);
fwrite(fid, text);
fclose(fid);
end

function p = localMakeFakeCmdstanRootAt(parent, name)
%A directory that satisfies stan_home's validity criterion: makefile + bin/.
p = fullfile(parent, name);
mkdir(p);
mkdir(fullfile(p, 'bin'));
fid = fopen(fullfile(p, 'makefile'), 'w');
fclose(fid);
end

function localRmdirIfExists(d)
if exist(d, 'dir') == 7
    rmdir(d, 's');
end
end

function localRmpathQuiet(d)
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(d);
warning('on', 'MATLAB:rmpath:DirNotFound');
rehash;
end

function [txt, nfig] = localDepsPromptText()
%the joined text of every update prompt currently up, and how many there
%are. Both launcher prompts share the Tag psyrat_gui_update (B45), so the
%text is what tells the dependents prompt from the toolbox prompt.
figs = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui_update');
nfig = numel(figs);
txt = '';
for k = 1:nfig
    h = findall(figs(k), 'Style', 'text');
    for j = 1:numel(h)
        s = get(h(j), 'String');
        if ischar(s) || isstring(s)
            s = cellstr(s);
        end
        txt = [txt ' ' strjoin(reshape(s, 1, []), ' ')]; %#ok<AGROW>
    end
end
end

function code = localCode(file)
%the file with whole-line comments removed (the launcher's ruling comments
%discuss the helper by name; only code may count)
lines = strsplit(fileread(file), newline);
keep = ~cellfun(@(l) startsWith(strtrim(l), '%'), lines);
code = strjoin(lines(keep), newline);
end
