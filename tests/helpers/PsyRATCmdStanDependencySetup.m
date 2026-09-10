classdef PsyRATCmdStanDependencySetup
    % Resolve and validate CmdStan-related dependencies for integration tests.
    %
    % Environment variable overrides. The three exact-path variables are checked
    % before anything else, one per dependent:
    %   PSYRAT_TEST_MPM_PATH        - exact MatlabProcessManager folder
    %   PSYRAT_TEST_MATLABSTAN_PATH - exact MatlabStan folder
    %   PSYRAT_TEST_CMDSTAN_PATH    - exact CmdStan root
    %
    % PSYRAT_DEPENDENTS is NOT in that tier. It names a PARENT folder holding
    % any/all of the three, searched by name pattern, for a machine whose layout
    % none of the home-relative fallbacks below matches. It is the FIRST search
    % root for MatlabProcessManager and CmdStan, but for MatlabStan the bundled,
    % patched copy in bundled_dependents/ is resolved ahead of it (CLAUDE.md: use
    % the bundled MatlabStan, not upstream), so pointing PSYRAT_DEPENDENTS at a
    % different MatlabStan does not override the bundled one -- set
    % PSYRAT_TEST_MATLABSTAN_PATH for that.
    %
    % Everything after those is home-relative discovery built from getenv('HOME')
    % (see localHomeDir / localDependentRoots), then the generic discovery in
    % localDiscoverCmdStanRoot. No absolute path is hardcoded here, so the helper
    % is not tied to one developer's account.

    methods (Static)
        function cmdstanRoot = configure(testCase)
            projectRoot = localResolveProjectRoot(testCase);
            paths = PsyRATCmdStanDependencySetup.resolvePaths(projectRoot);

            if ~isempty(paths.mpm)
                addpath(genpath(paths.mpm));
            end
            if ~isempty(paths.matlabstan)
                addpath(genpath(paths.matlabstan));
            end

            testCase.assumeFalse(isempty(which('processManager')), ...
                ['processManager is unavailable. Set PSYRAT_TEST_MPM_PATH ', ...
                'or add MatlabProcessManager to MATLAB path.']);
            testCase.assumeFalse(isempty(which('stan')), ...
                ['stan is unavailable. Set PSYRAT_TEST_MATLABSTAN_PATH ', ...
                'or add MatlabStan to MATLAB path.']);

            cmdstanRoot = paths.cmdstan;
            if isempty(cmdstanRoot)
                cmdstanRoot = localDiscoverCmdStanRoot();
            end

            testCase.assumeTrue(~isempty(cmdstanRoot) && isfolder(cmdstanRoot), ...
                ['CmdStan root path not found. Set PSYRAT_TEST_CMDSTAN_PATH ', ...
                'or ensure CmdStan is discoverable in MATLAB path.']);
            testCase.assumeTrue(localHasAnyFile(cmdstanRoot, {'stanc', 'stanc3'}), ...
                'CmdStan compiler binary (stanc/stanc3) is unavailable.');
            testCase.assumeTrue(localHasAnyFile(cmdstanRoot, {'print'}), ...
                'CmdStan print helper binary is unavailable.');
            testCase.assumeTrue(localHasAnyFile(cmdstanRoot, {'stansummary'}), ...
                'CmdStan stansummary helper binary is unavailable.');

            % Keep MatlabStan runtime discovery aligned with test resolution.
            setenv('PSYRAT_TEST_CMDSTAN_PATH', cmdstanRoot);
            setenv('STAN_HOME', cmdstanRoot);
        end

        function paths = resolvePaths(projectRoot)
            paths = struct();
            paths.mpm = getenv('PSYRAT_TEST_MPM_PATH');
            paths.matlabstan = getenv('PSYRAT_TEST_MATLABSTAN_PATH');
            paths.cmdstan = getenv('PSYRAT_TEST_CMDSTAN_PATH');

            % Roots that may hold the three dependents, most specific first.
            % Built once so the three searches below stay in step and no
            % account-specific absolute path appears anywhere in this file.
            roots = localDependentRoots();

            if isempty(paths.mpm)
                paths.mpm = localFirstByPattern( ...
                    localPatternsUnder(roots,'MatlabProcessManager-*'));
            end

            % Bundled MatlabProcessManager, as a LAST-RESORT fallback. Note the
            % ordering is deliberately the OPPOSITE of MatlabStan's below: the
            % bundled MatlabStan is patched and must win outright, whereas this
            % is an unmodified copy, so a machine's real install keeps
            % precedence and nothing changes where discovery already succeeds.
            % Without this tier the FIRST assumption in configure() fails on any
            % machine with no PsyRATDependents folder - notably a CI runner,
            % where the observed skip was 'processManager is unavailable', NOT a
            % missing CmdStan. Installing CmdStan alone would not have helped.
            if isempty(paths.mpm)
                paths.mpm = localFirstExisting({
                    fullfile(projectRoot, 'bundled_dependents', 'MatlabProcessManager-0.5.1')
                    });
            end

            % The bundled, patched MatlabStan wins outright (CLAUDE.md: use it,
            % not upstream), so it stays ahead of any discovered copy.
            if isempty(paths.matlabstan)
                paths.matlabstan = localFirstExisting({
                    fullfile(projectRoot, 'bundled_dependents', 'MatlabStan-2.15.1.0')
                    });
            end
            if isempty(paths.matlabstan)
                paths.matlabstan = localFirstByPattern( ...
                    localPatternsUnder(roots,'MatlabStan-*'));
            end

            if isempty(paths.cmdstan)
                paths.cmdstan = localFirstByPattern( ...
                    [localPatternsUnder(roots,'cmdstan-*'); ...
                    {fullfile(localHomeDir(), 'cmdstan-*')}]);
            end

            if isempty(paths.cmdstan)
                stanHome = getenv('STAN_HOME');
                if ~isempty(stanHome) && isfolder(stanHome)
                    paths.cmdstan = stanHome;
                end
            end

            if isempty(paths.cmdstan) && exist('mstan.stan_home', 'file') == 2
                try
                    candidate = mstan.stan_home();
                    if ischar(candidate) && isfolder(candidate)
                        paths.cmdstan = candidate;
                    end
                catch
                    % ignore discovery errors and keep fallback behavior
                end
            end
        end
    end
end

function projectRoot = localResolveProjectRoot(testCase)
projectRoot = '';

if nargin > 0 && ~isempty(testCase)
    % Some MATLAB versions restrict calling helper methods on TestCase
    % objects from external utility classes. Use best-effort probing and
    % fall back to this file's location if access is restricted.
    if ismethod(testCase,'projectRoot')
        try
            projectRoot = testCase.projectRoot();
        catch
            projectRoot = '';
        end
    end
end

if isempty(projectRoot) || ~isfolder(projectRoot)
    thisFile = mfilename('fullpath');
    % tests/helpers/PsyRATCmdStanDependencySetup.m -> project root
    projectRoot = fileparts(fileparts(fileparts(thisFile)));
end
end

function out = localFirstExisting(candidates)
out = '';
for i = 1:numel(candidates)
    if isfolder(candidates{i})
        out = candidates{i};
        return;
    end
end
end

function out = localFirstByPattern(patterns)
out = '';
for i = 1:numel(patterns)
    hits = dir(patterns{i});
    hits = hits([hits.isdir]);
    hits = hits(~ismember({hits.name}, {'.', '..'}));
    if isempty(hits)
        continue;
    end
    %Reverse LEXICAL order, which is NOT a version sort. sort() accepts only ONE
    %argument for a cell array of char, so the 'descend' flag this line used to
    %pass errored outright ("Only one input argument is supported for cell
    %arrays"). That went unnoticed because the exact-path tier above always
    %answered first on the maintainer's machine and this branch was never
    %reached; making the search home-relative reaches it. Sort ascending, then
    %reverse.
    %
    %The ordering is by character code, so it agrees with version order only
    %while the minor numbers have the same digit count. A folder holding both
    %cmdstan-2.9.0 and cmdstan-2.38.0 resolves 2.9.0 ('9' > '3'), i.e. the OLDER
    %release. That is tolerable here and deliberately not fixed: this is
    %test-helper discovery, it only bites when ONE searched root holds more than
    %one cmdstan-*, and CmdStan minors have been two digits continuously since
    %2.10, so the failure needs a 2.9-or-older install sitting alongside a modern
    %one (the same hazard returns at 2.100). If a real machine ever hits it, set
    %PSYRAT_TEST_CMDSTAN_PATH (an exact path, checked first) rather than teaching
    %this helper to parse versions.
    [~, order] = sort({hits.name});
    order = fliplr(order(:)');
    for j = order
        candidate = fullfile(hits(j).folder, hits(j).name);
        if isfolder(candidate)
            out = candidate;
            return;
        end
    end
end
end

function d = localHomeDir()
%User home directory. HOME is set on macOS/Linux and normally on Windows too;
%when it is not, derive it from userpath, whose factory default is
%<home>/Documents/MATLAB (macOS/Linux) or <home>\Documents\MATLAB (Windows), so
%two fileparts steps recover the home directory. MATLAB does not expand '~' in
%isfolder/dir, so the literal tilde is only a last resort that fails closed.
d = getenv('HOME');
if isempty(d)
    d = getenv('USERPROFILE');
end
if isempty(d)
    up = userpath();
    if ~isempty(up)
        up = strsplit(up, pathsep);
        candidate = fileparts(fileparts(up{1}));
        if ~isempty(candidate) && isfolder(candidate)
            d = candidate;
        end
    end
end
if isempty(d)
    d = '~';
end
end

function roots = localDependentRoots()
%Parent folders that may contain CmdStan / MatlabStan / MatlabProcessManager,
%searched in order. PSYRAT_DEPENDENTS lets a machine whose layout matches none
%of the home-relative conventions point the suite at its own folder without
%editing this file. The remaining entries are the layouts this project has used,
%expressed relative to the home directory rather than to one account.
roots = {};

envRoot = getenv('PSYRAT_DEPENDENTS');
if ~isempty(envRoot)
    roots{end+1,1} = envRoot;
end

h = localHomeDir();
roots{end+1,1} = fullfile(h, 'Documents', 'MATLAB', 'PsyRATDependents');
roots{end+1,1} = fullfile(h, 'Documents', 'Documents_MacStudio', 'MATLAB', ...
    'ERPToolboxes', 'PsyRATDependents');
roots{end+1,1} = fullfile(h, 'Documents', 'Documents_MacStudio', 'MATLAB', ...
    'ERPToolboxes', 'ERADependents');
end

function patterns = localPatternsUnder(roots, leaf)
%Expand a set of root folders into dir() patterns for one dependent, preserving
%the root order so the caller's precedence is the root precedence.
patterns = cell(numel(roots),1);
for i = 1:numel(roots)
    patterns{i} = fullfile(roots{i}, leaf);
end
end

function cmdstanRoot = localDiscoverCmdStanRoot()
cmdstanRoot = '';

% Last resort: infer from MatlabStan defaults.
if exist('mstan.stan_home', 'file') == 2
    try
        candidate = mstan.stan_home();
        if ischar(candidate) && isfolder(candidate)
            cmdstanRoot = candidate;
            return;
        end
    catch
        % ignore and continue
    end
end

% Infer from binaries on MATLAB path.
for name = {'stanc', 'stanc3', 'print', 'stansummary'}
    p = which(name{1});
    if ~isempty(p)
        maybeRoot = fileparts(fileparts(p));
        if isfolder(fullfile(maybeRoot, 'bin'))
            cmdstanRoot = maybeRoot;
            return;
        end
    end
end
end

function tf = localHasAnyFile(cmdstanRoot, names)
tf = false;
for i = 1:numel(names)
    if isfile(fullfile(cmdstanRoot, 'bin', names{i})) || ...
            isfile(fullfile(cmdstanRoot, 'bin', [names{i} '.exe']))
        tf = true;
        return;
    end
end
end
