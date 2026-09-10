classdef TestInstallationAndWorkflowValidation < PsyRATTestBase

    methods (Test)
        function testDependentsVersionStruct(testCase)
            depvers = psyrat_dependentsversions();
            % matlabstanpatch is PsyRAT's own patch level for the vendored
            % MatlabStan bundle, not an upstream release string (B41). It is
            % compared against PSYRAT_PATCHLEVEL.txt inside the bundle; see
            % TestMatlabStanPatchLevel for the two being kept in step.
            testCase.verifyEqual(sort(fieldnames(depvers)), ...
                sort({'cmdstan'; 'matlabstan'; 'matlabstanpatch'; ...
                'matlabprocessmanager'}));
            testCase.verifyClass(depvers.cmdstan, 'char');
            testCase.verifyClass(depvers.matlabstan, 'char');
            testCase.verifyClass(depvers.matlabstanpatch, 'char');
            testCase.verifyClass(depvers.matlabprocessmanager, 'char');
        end

        function testCheckVersionsOfDepsWithMockedFiles(testCase)
            mockDir = localCreateDependencyMockDir();
            c = onCleanup(@() rmpath(mockDir));
            addpath(mockDir, '-begin'); %#ok<NASGU>

            versions = psyrat_checkversionsofdeps();
            testCase.verifyEqual(versions.cmdstan, 1);
            % Was pinned at the -1 'not useful' sentinel until the bundle
            % patch-level stamp made a real MatlabStan check possible.
            testCase.verifyEqual(versions.matlabstan, 1);
            testCase.verifyEqual(versions.matlabprocessmanager, 1);
        end

        function testPrintVersionsWithMockedFiles(testCase)
            mockDir = localCreateDependencyMockDir();
            c = onCleanup(@() rmpath(mockDir));
            addpath(mockDir, '-begin'); %#ok<NASGU>

            txt = evalc('psyrat_printversions');
            testCase.verifyNotEmpty(strfind(txt, 'Matlab version')); %#ok<STREMP>
            testCase.verifyNotEmpty(strfind(txt, 'PsyRAT Toolbox version')); %#ok<STREMP>
            testCase.verifyNotEmpty(strfind(txt, 'CmdStan version')); %#ok<STREMP>
        end

        function testUpdateCheckWithMockedWebread(testCase)
            currentVer = psyrat_defineversion();
            olderVer = localStepVersion(currentVer, -1);
            newerVer = localStepVersion(currentVer, 1);

            mockDir = localCreateWebreadMockDir(currentVer);
            c = onCleanup(@() rmpath(mockDir));
            addpath(mockDir, '-begin'); %#ok<NASGU>

            testCase.verifyEqual(psyrat_updatecheck(currentVer), 1);
            testCase.verifyEqual(psyrat_updatecheck(olderVer), 0);
            testCase.verifyEqual(psyrat_updatecheck(newerVer), 2);
        end

        function testGithubReleaseCheckHandlesNoStableReleases(testCase)
            info = psyrat_checkgithubreleases( ...
                'currentversion', '0.1.0', ...
                'repo', 'peclayson/PsyRAT', ...
                'verbose', false, ...
                'webreadfcn', @localMockNoReleases);

            testCase.verifyEqual(info.statusCode, 3);
            testCase.verifyEqual(info.status, 'no_stable_release');
            testCase.verifyFalse(info.hasStableRelease);
            testCase.verifyTrue(info.isBeta);
            testCase.verifyNotEmpty(strfind(info.message, 'No official stable GitHub releases')); %#ok<STREMP>
        end

        function testGithubReleaseCheckDetectsUpdateAvailable(testCase)
            info = psyrat_checkgithubreleases( ...
                'currentversion', '0.1.0', ...
                'repo', 'peclayson/PsyRAT', ...
                'verbose', false, ...
                'webreadfcn', @localMockStableRelease020);

            testCase.verifyEqual(info.statusCode, 0);
            testCase.verifyEqual(info.status, 'update_available');
            testCase.verifyEqual(info.latestReleaseVersion, '0.2.0');
            testCase.verifyFalse(info.isBeta);
        end

        function testGithubReleaseCheckDetectsBetaBuildAhead(testCase)
            info = psyrat_checkgithubreleases( ...
                'currentversion', '0.3.0', ...
                'repo', 'peclayson/PsyRAT', ...
                'verbose', false, ...
                'webreadfcn', @localMockStableRelease020);

            testCase.verifyEqual(info.statusCode, 2);
            testCase.verifyEqual(info.status, 'beta_ahead_of_latest_release');
            testCase.verifyTrue(info.isBeta);
            testCase.verifyNotEmpty(strfind(info.message, 'potentially unstable beta/development')); %#ok<STREMP>
        end

        function testGithubReleaseCheckHandlesConnectivityErrors(testCase)
            info = psyrat_checkgithubreleases( ...
                'currentversion', '0.1.0', ...
                'repo', 'peclayson/PsyRAT', ...
                'verbose', false, ...
                'webreadfcn', @localMockReleaseError);

            testCase.verifyEqual(info.statusCode, -1);
            testCase.verifyEqual(info.status, 'unable_to_check');
            testCase.verifyNotEmpty(strfind(info.howto, 'internet connection')); %#ok<STREMP>
        end

        function testGithubReleaseCheckTreats404AsBetaContext(testCase)
            info = psyrat_checkgithubreleases( ...
                'currentversion', '0.1.0', ...
                'repo', 'peclayson/PsyRAT', ...
                'verbose', false, ...
                'webreadfcn', @localMockRelease404);

            testCase.verifyEqual(info.statusCode, 3);
            testCase.verifyEqual(info.status, 'no_stable_release');
            testCase.verifyTrue(info.isBeta);
            testCase.verifyNotEmpty(strfind(info.message, 'potentially unstable beta/development')); %#ok<STREMP>
            testCase.verifyNotEmpty(strfind(info.howto, 'private beta repo')); %#ok<STREMP>
        end

        function testGithubReleaseCheckExposesThePsyratZipAssetUrl(testCase)
            %G87 (owner ruling 2026-09-09): the toolbox updater downloads the
            %PsyRAT.zip asset that release/PUBLISHING.md attaches to every
            %release, resolved through this check instead of by scraping the
            %releases web page. A struct-array assets field is what webread
            %returns for a JSON array of same-shaped objects.
            info = psyrat_checkgithubreleases( ...
                'currentversion', '0.1.0', ...
                'repo', 'peclayson/PsyRAT', ...
                'verbose', false, ...
                'webreadfcn', @localMockStableRelease020WithAssets);

            testCase.verifyEqual(info.statusCode, 0);
            testCase.verifyEqual(info.latestReleaseAssetURL, ...
                'https://github.com/peclayson/PsyRAT/releases/download/v0.2.0/PsyRAT.zip');
        end

        function testGithubReleaseCheckFindsTheAssetInACellOfStructs(testCase)
            %webread returns a JSON array as a CELL of structs when the objects'
            %fields differ; the asset lookup tolerates that shape the same way
            %local_fetch_releases tolerates it for the releases themselves.
            info = psyrat_checkgithubreleases( ...
                'currentversion', '0.1.0', ...
                'repo', 'peclayson/PsyRAT', ...
                'verbose', false, ...
                'webreadfcn', @localMockStableRelease020WithCellAssets);

            testCase.verifyEqual(info.statusCode, 0);
            testCase.verifyEqual(info.latestReleaseAssetURL, ...
                'https://github.com/peclayson/PsyRAT/releases/download/v0.2.0/PsyRAT.zip');
        end

        function testGithubReleaseCheckAssetUrlIsEmptyWithoutTheAsset(testCase)
            %no assets field at all (the existing stable-release mock)
            info = psyrat_checkgithubreleases( ...
                'currentversion', '0.1.0', ...
                'repo', 'peclayson/PsyRAT', ...
                'verbose', false, ...
                'webreadfcn', @localMockStableRelease020);
            testCase.verifyEqual(info.statusCode, 0);
            testCase.verifyEqual(info.latestReleaseAssetURL, '');

            %assets present, none of them named PsyRAT.zip
            info = psyrat_checkgithubreleases( ...
                'currentversion', '0.1.0', ...
                'repo', 'peclayson/PsyRAT', ...
                'verbose', false, ...
                'webreadfcn', @localMockStableRelease020WithOtherAssets);
            testCase.verifyEqual(info.statusCode, 0);
            testCase.verifyEqual(info.latestReleaseAssetURL, '');
        end

        function testGithubReleaseCheckNeverExposesAPrereleaseAsset(testCase)
            %a beta published with --prerelease carries PsyRAT.zip too, but the
            %updater installs stable releases only: the asset URL stays empty
            %and the status is the no-stable-release beta context.
            info = psyrat_checkgithubreleases( ...
                'currentversion', '0.1.0', ...
                'repo', 'peclayson/PsyRAT', ...
                'verbose', false, ...
                'webreadfcn', @localMockPrereleaseOnlyWithAsset);

            testCase.verifyEqual(info.statusCode, 3);
            testCase.verifyEqual(info.status, 'no_stable_release');
            testCase.verifyTrue(info.isBeta);
            testCase.verifyEqual(info.latestReleaseAssetURL, '');
        end

        function testComputerelWrapInputValidation(testCase)
            testCase.verifyError(@() psyrat_computevarcompwarp('psyrat_data', struct()), ...
                'psyrat_prefs:notfound');
            testCase.verifyError(@() psyrat_computevarcompwarp('psyrat_prefs', struct()), ...
                'psyrat_data:notfound');

            prefs = psyrat_defaults();
            data = struct();
            data.proc = struct();
            data.proc.data = table({'s1'}, 1.2, 'VariableNames', {'id', 'meas'});
            data.proc.savepath = tempdir;
            data.proc.savename = '';
            testCase.verifyError(@() psyrat_computevarcompwarp('psyrat_prefs', prefs, 'psyrat_data', data), ...
                'psyrat_data:savename');

            data.proc.savename = 'okname';
            data.proc = rmfield(data.proc, 'savepath');
            testCase.verifyError(@() psyrat_computevarcompwarp('psyrat_prefs', prefs, 'psyrat_data', data), ...
                'psyrat_data:savepath');

            data.proc.savepath = tempdir;
            data.proc.savename = 'bad name';
            testCase.verifyError(@() psyrat_computevarcompwarp('psyrat_prefs', prefs, 'psyrat_data', data), ...
                'psyrat_data:savename');
        end

        function testComputeWrapSavesWhenUserDeclinesNonconvergedRerun(testCase)
            mockDir = localCreateNonconvergedRerunMockDir();
            cleanupPath = onCleanup(@() rmpath(mockDir)); %#ok<NASGU>
            addpath(mockDir, '-begin');
            clear psyrat_computevarcomp psyrat_checkconv psyrat_reruncheck psyrat_startproc;

            saveDir = tempname;
            mkdir(saveDir);
            cleanupDir = onCleanup(@() rmdir(saveDir, 's')); %#ok<NASGU>

            prefs = psyrat_defaults();
            prefs.proc.nchains = 3;
            prefs.proc.nwarmup = 50;
            prefs.proc.nsampling = 50;
            prefs.proc.niter = 100;
            prefs.proc.diffest = 1;
            prefs.proc.diffwpcov = 1;
            prefs.proc.diffrescor = 1;
            prefs.proc.ssubjrel = 1;
            prefs.proc.sserrvar = 1;
            prefs.proc.dynrel = 1;
            prefs.proc.traceplots = 1;

            data = struct();
            data.proc = struct();
            data.proc.data = table({'s1'; 's1'}, [1.1; 1.3], ...
                'VariableNames', {'id', 'meas'});
            data.proc.savepath = saveDir;
            data.proc.savename = 'nonconverged.psyrat';
            data.proc.measheader = 'meas';

            out = psyrat_computevarcompwarp('psyrat_prefs', prefs, ...
                'psyrat_data', data);

            testCase.verifyTrue(isfield(out, 'rel'));
            testCase.verifyEqual(out.rel.analysis, 'ic_dynrel');
            testCase.verifyFalse(out.rel.out.conv.converged);
            testCase.verifyEqual(evalin('base', 'PSYRAT_TEST_STARTPROC_CALLED'), 0);
            testCase.verifyTrue(isfile(fullfile(saveDir, 'nonconverged.psyrat')));
        end

        function testComputeWrapTreatsADismissedRerunDialogAsDoNotRerun(testCase)
            %G13: the rerun dialog can be gone before the warp reads its
            %answer (a close handler that deletes, or a user closing MATLAB's
            %figure by other means), so findobj returns empty. The warp must
            %treat an absent dialog as "do not rerun" and save the
            %non-converged result rather than error after the fit finished.
            %Same harness as the sibling above; only the dialog mock differs,
            %creating NO figure at all.
            mockDir = localCreateNonconvergedRerunMockDir();
            fid = fopen(fullfile(mockDir, 'psyrat_reruncheck.m'), 'w');
            fprintf(fid, ['function psyrat_reruncheck\n' ...
                '%%dismissed: no figure survives for the warp to read\n' ...
                'end\n']);
            fclose(fid);
            cleanupPath = onCleanup(@() rmpath(mockDir)); %#ok<NASGU>
            addpath(mockDir, '-begin');
            clear psyrat_computevarcomp psyrat_checkconv psyrat_reruncheck psyrat_startproc;

            saveDir = tempname;
            mkdir(saveDir);
            cleanupDir = onCleanup(@() rmdir(saveDir, 's')); %#ok<NASGU>

            prefs = psyrat_defaults();
            prefs.proc.nchains = 3;
            prefs.proc.nwarmup = 50;
            prefs.proc.nsampling = 50;
            prefs.proc.niter = 100;
            prefs.proc.diffest = 1;
            prefs.proc.diffwpcov = 1;
            prefs.proc.diffrescor = 1;
            prefs.proc.ssubjrel = 1;
            prefs.proc.sserrvar = 1;
            prefs.proc.dynrel = 1;
            prefs.proc.traceplots = 1;

            data = struct();
            data.proc = struct();
            data.proc.data = table({'s1'; 's1'}, [1.1; 1.3], ...
                'VariableNames', {'id', 'meas'});
            data.proc.savepath = saveDir;
            data.proc.savename = 'dismissed.psyrat';
            data.proc.measheader = 'meas';

            out = psyrat_computevarcompwarp('psyrat_prefs', prefs, ...
                'psyrat_data', data);

            testCase.verifyTrue(isfield(out, 'rel'));
            testCase.verifyFalse(out.rel.out.conv.converged);
            testCase.verifyEqual(evalin('base', 'PSYRAT_TEST_STARTPROC_CALLED'), 0);
            testCase.verifyTrue(isfile(fullfile(saveDir, 'dismissed.psyrat')));
        end

        function testComputeWrapPrintsPreflightWarnings(testCase)
            %G8: psyrat_preflight_validate's non-fatal warnings (few
            %participants, low median trial count, incomplete cells) were
            %computed and dropped by the warp, which read only errors(1).
            %The scripted route has no dialog, so the console line is its
            %only channel. One participant with two trials trips the
            %few-participants warning.
            [prefs, data, saveDir, cleanups] = localNonconvergedRerunSetup('warned.psyrat'); %#ok<ASGLU>
            out = evalc(['psyrat_computevarcompwarp(''psyrat_prefs'', prefs, ', ...
                '''psyrat_data'', data);']);
            testCase.verifyTrue(contains(out, 'Preflight warning'), ...
                'the warp must print the preflight warnings before the fit');
            testCase.verifyTrue(contains(out, 'preflight:fewParticipants'), ...
                'a one-participant table must print the few-participants warning');
            testCase.verifyTrue(isfile(fullfile(saveDir, 'warned.psyrat')), ...
                'a warning is non-fatal: the run must still complete and save');
        end

        function testComputeWrapStaysQuietOnAHealthyDataset(testCase)
            %Negative control for the test above: 22 participants with 12
            %trials each trips no preflight warning, so nothing may print.
            [prefs, data, saveDir, cleanups] = localNonconvergedRerunSetup('healthy.psyrat'); %#ok<ASGLU>
            ids = repelem(arrayfun(@(k) sprintf('s%02d', k), 1:22, ...
                'UniformOutput', false)', 12, 1);
            data.proc.data = table(ids, (1:numel(ids))' * 0.01, ...
                'VariableNames', {'id', 'meas'});
            out = evalc(['psyrat_computevarcompwarp(''psyrat_prefs'', prefs, ', ...
                '''psyrat_data'', data);']);
            testCase.verifyFalse(contains(out, 'Preflight warning'), ...
                'a healthy dataset must not print a preflight warning');
            testCase.verifyTrue(isfile(fullfile(saveDir, 'healthy.psyrat')));
        end
    end
end

function mockDir = localCreateDependencyMockDir()
mockDir = tempname;
mkdir(mockDir);

cmdstanGuide = fullfile(mockDir, 'cmdstan-guide.tex');
fid = fopen(cmdstanGuide, 'w');
fprintf(fid, 'a\\b\\c\\d\\version{2.38.0}\n');
fclose(fid);

procMgr = fullfile(mockDir, 'processManager.m');
fid = fopen(procMgr, 'w');
fprintf(fid, [ ...
    'function p = processManager(varargin)\n' ...
    'version = ''0.5.1''; %#ok<NASGU>\n' ...
    'p = struct();\n' ...
    'p.stdout = {''stanc 2.38.0''};\n' ...
    'p.block = @(varargin) [];\n' ...
    'end\n']);
fclose(fid);

% A mock MatlabStan, so the MatlabStan check resolves here rather than falling
% through to the repository's real bundle. psyrat_matlabstan_bundleinfo anchors
% on which('StanModel.m') and then reads the stamp from that same directory, so
% both files are required and both must sit in this directory.
stanModel = fullfile(mockDir, 'StanModel.m');
fid = fopen(stanModel, 'w');
fprintf(fid, [ ...
    'function m = StanModel(varargin)\n' ...
    'm = struct();\n' ...
    'end\n']);
fclose(fid);

depvers = psyrat_dependentsversions;
localWriteStamp(mockDir, depvers.matlabstan, depvers.matlabstanpatch);
end

function localWriteStamp(msdir, upstreamVersion, patchLevel)
%LOCALWRITESTAMP Write a PSYRAT_PATCHLEVEL.txt into a mock MatlabStan directory.
fid = fopen(fullfile(msdir, 'PSYRAT_PATCHLEVEL.txt'), 'w');
fprintf(fid, 'PsyRAT vendored-bundle stamp\n\n');
fprintf(fid, '    upstream_version = %s\n', upstreamVersion);
fprintf(fid, '    patch_level = %s\n', patchLevel);
fclose(fid);
end

function mockDir = localCreateWebreadMockDir(returnVersion)
mockDir = tempname;
mkdir(mockDir);

webreadFile = fullfile(mockDir, 'webread.m');
fid = fopen(webreadFile, 'w');
fprintf(fid, [ ...
    'function out = webread(varargin)\n' ...
    'out = ''%s'';\n' ...
    'end\n'], returnVersion);
fclose(fid);
end

function mockDir = localCreateNonconvergedRerunMockDir()
mockDir = tempname;
mkdir(mockDir);

fid = fopen(fullfile(mockDir, 'psyrat_computevarcomp.m'), 'w');
fprintf(fid, [ ...
    'function REL = psyrat_computevarcomp(varargin)\n' ...
    'REL = struct();\n' ...
    'REL.analysis = ''ic_dynrel'';\n' ...
    'REL.out = struct();\n' ...
    'REL.out.conv = struct(''converged'', false);\n' ...
    'end\n']);
fclose(fid);

fid = fopen(fullfile(mockDir, 'psyrat_checkconv.m'), 'w');
fprintf(fid, [ ...
    'function RELout = psyrat_checkconv(REL)\n' ...
    'RELout = REL;\n' ...
    'end\n']);
fclose(fid);

fid = fopen(fullfile(mockDir, 'psyrat_reruncheck.m'), 'w');
fprintf(fid, [ ...
    'function psyrat_reruncheck\n' ...
    'f = figure(''Visible'', ''off'', ''Tag'', ''psyrat_gui'');\n' ...
    'guidata(f, 0);\n' ...
    'end\n']);
fclose(fid);

fid = fopen(fullfile(mockDir, 'psyrat_startproc.m'), 'w');
fprintf(fid, [ ...
    'function psyrat_startproc(varargin)\n' ...
    'assignin(''base'', ''PSYRAT_TEST_STARTPROC_CALLED'', 1);\n' ...
    'end\n']);
fclose(fid);

assignin('base', 'PSYRAT_TEST_STARTPROC_CALLED', 0);
end

function [prefs, data, saveDir, cleanups] = localNonconvergedRerunSetup(savename)
%The nonconverged-rerun harness above, packaged for reuse: mocks on the path,
%a temporary save directory, the minimal one-participant table, and the
%prefs the sibling tests use. CLEANUPS must stay alive in the caller.
mockDir = localCreateNonconvergedRerunMockDir();
addpath(mockDir, '-begin');
clear psyrat_computevarcomp psyrat_checkconv psyrat_reruncheck psyrat_startproc;
saveDir = tempname;
mkdir(saveDir);
cleanups = {onCleanup(@() rmpath(mockDir)), onCleanup(@() rmdir(saveDir, 's'))};

prefs = psyrat_defaults();
prefs.proc.nchains = 3;
prefs.proc.nwarmup = 50;
prefs.proc.nsampling = 50;
prefs.proc.niter = 100;
prefs.proc.diffest = 1;
prefs.proc.diffwpcov = 1;
prefs.proc.diffrescor = 1;
prefs.proc.ssubjrel = 1;
prefs.proc.sserrvar = 1;
prefs.proc.dynrel = 1;
prefs.proc.traceplots = 1;

data = struct();
data.proc = struct();
data.proc.data = table({'s1'; 's1'}, [1.1; 1.3], ...
    'VariableNames', {'id', 'meas'});
data.proc.savepath = saveDir;
data.proc.savename = savename;
data.proc.measheader = 'meas';
end

function nextVer = localStepVersion(ver, direction)
parts = sscanf(ver, '%d.%d.%d')';
if numel(parts) ~= 3
    error('test:versionparse', 'Version must be semantic (x.y.z).');
end

if direction < 0
    if parts(3) > 0
        parts(3) = parts(3) - 1;
    elseif parts(2) > 0
        parts(2) = parts(2) - 1;
        parts(3) = 9;
    elseif parts(1) > 0
        parts(1) = parts(1) - 1;
        parts(2) = 9;
        parts(3) = 9;
    end
elseif direction > 0
    parts(3) = parts(3) + 1;
end

nextVer = sprintf('%d.%d.%d', parts(1), parts(2), parts(3));
end

function out = localMockNoReleases(varargin) %#ok<INUSD>
out = struct([]);
end

function out = localMockStableRelease020(varargin) %#ok<INUSD>
out = struct( ...
    'tag_name', 'v0.2.0', ...
    'name', 'Version 0.2.0', ...
    'draft', false, ...
    'prerelease', false, ...
    'html_url', 'https://github.com/peclayson/PsyRAT/releases/tag/v0.2.0');
end

function out = localMockReleaseError(varargin) %#ok<INUSD>
error('test:mockReleaseError', ...
    ['Mocked GitHub release lookup failed. ', ...
    'How to fix: this error is intentional for connectivity handling tests.']);
out = struct([]);
end

function out = localMockRelease404(varargin) %#ok<INUSD>
error('MATLAB:webservices:HTTP404StatusCodeError', ...
    ['The server returned the status 404 with message "Not Found" ', ...
    'in response to the request URL.']);
out = struct([]);
end

function out = localMockStableRelease020WithAssets(varargin) %#ok<INUSD>
%the stable release with two assets as a struct array (same-shaped objects)
out = localMockStableRelease020();
out.assets = struct( ...
    'name', {'release-notes.md', 'PsyRAT.zip'}, ...
    'browser_download_url', ...
    {'https://github.com/peclayson/PsyRAT/releases/download/v0.2.0/release-notes.md', ...
    'https://github.com/peclayson/PsyRAT/releases/download/v0.2.0/PsyRAT.zip'});
end

function out = localMockStableRelease020WithCellAssets(varargin) %#ok<INUSD>
%the same two assets as a cell of structs whose fields differ
out = localMockStableRelease020();
out.assets = {struct( ...
    'name', 'release-notes.md', ...
    'browser_download_url', 'https://github.com/peclayson/PsyRAT/releases/download/v0.2.0/release-notes.md', ...
    'size', 12), ...
    struct( ...
    'name', 'PsyRAT.zip', ...
    'browser_download_url', 'https://github.com/peclayson/PsyRAT/releases/download/v0.2.0/PsyRAT.zip')};
end

function out = localMockStableRelease020WithOtherAssets(varargin) %#ok<INUSD>
%assets present, but no PsyRAT.zip among them
out = localMockStableRelease020();
out.assets = struct( ...
    'name', {'release-notes.md'}, ...
    'browser_download_url', ...
    {'https://github.com/peclayson/PsyRAT/releases/download/v0.2.0/release-notes.md'});
end

function out = localMockPrereleaseOnlyWithAsset(varargin) %#ok<INUSD>
%only a beta prerelease exists, and it does carry PsyRAT.zip
out = struct( ...
    'tag_name', 'v0.1.0-beta', ...
    'name', 'PsyRAT v0.1.0-beta', ...
    'draft', false, ...
    'prerelease', true, ...
    'html_url', 'https://github.com/peclayson/PsyRAT/releases/tag/v0.1.0-beta');
out.assets = struct( ...
    'name', {'PsyRAT.zip'}, ...
    'browser_download_url', ...
    {'https://github.com/peclayson/PsyRAT/releases/download/v0.1.0-beta/PsyRAT.zip'});
end
