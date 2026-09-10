classdef TestVerboseQuietRun < PsyRATTestBase
    % Regression test for the verbose==0 CmdStan-path crash.
    %
    % psyrat_computevarcomp maps the user-facing 'verbose' option to the
    % logical 'verbosity' flag it hands to the sampler. The original mapping
    %   if verbose == 1; verbosity = false; elseif verbose == 2; verbosity = true; end
    % had no else branch, so any other value (notably verbose == 0) left
    % 'verbosity' undefined. The default is verbose == 1, so the GUI/CLI default
    % path was safe, but an explicit verbose == 0 crashed with an undefined-
    % variable error the first time 'verbosity' was passed into a Stan call.
    %
    % This test runs a tiny one-facet fit with verbose == 0 and verifies the
    % pipeline reaches the Stan path and returns a populated fit. Before the
    % fix this run errored out before producing any output; after the fix it
    % completes quietly. It does not check recovery (covered by the seeded
    % smoke test); it only locks in that verbose == 0 no longer crashes.
    %
    % Like the other CmdStan integration tests, the dependency setup uses
    % assumptions, so the test is skipped cleanly (Incomplete, not Failed)
    % whenever CmdStan/MatlabStan are unavailable (e.g. on CI runners).

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()));
        end
    end

    methods (Test)
        function testVerboseZeroReachesStanPathWithoutCrash(testCase)
            % Small fixed one-facet dataset; recovery is not the point here, so
            % the RNG seed only needs to make the run deterministic.
            sigU = 2.0; sigE = 1.0; mu = 0;
            nsub = 20; ntrl = 6;
            rng(42);
            ids  = repelem((1:nsub)', ntrl);
            u    = sigU * randn(nsub, 1);
            e    = sigE * randn(nsub * ntrl, 1);
            meas = mu + u(ids) + e;
            dataTbl = table(ids, meas, 'VariableNames', {'id', 'meas'});

            outDir = fullfile(localCmdStanOutputDir(testCase.projectRoot()), ...
                'verbose0_run');
            if ~exist(outDir, 'dir'), mkdir(outDir); end

            % verbose == 0 is the out-of-contract value that previously left
            % 'verbosity' undefined. The call must complete without erroring.
            rel = psyrat_computevarcomp( ...
                'data', dataTbl, ...
                'chains', 1, ...
                'warmup', 200, ...
                'sampling', 200, ...
                'seed', 12345, ...
                'verbose', 0, ...
                'showgui', 1, ...
                'cmdstanoutdir', outDir);

            % Reaching here at all means the verbosity mapping no longer
            % crashed; confirm the fit is the expected populated one-facet run.
            testCase.verifyEqual(rel.analysis, 'ic');
            testCase.verifyTrue(isfield(rel, 'out'));
            for f = {'mu', 'sig_u', 'sig_e'}
                testCase.verifyTrue(isfield(rel.out, f{1}));
                testCase.verifyFalse(isempty(rel.out.(f{1})));
            end
        end
    end
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
% Remove the per-run output directory and any stray CmdStan build/output
% files left in the base artifacts directory.
p = fullfile(baseDir, 'verbose0_run');
if isfolder(p)
    try
        rmdir(p, 's');
    catch
        warning('tests:cleanupDirFailed', ...
            ['Could not remove verbose0 artifact directory ''%s''. ', ...
            'How to fix: close any process using it and delete it manually.'], p);
    end
end

pats = {'cmdstan*.stan', 'cmdstan*.hpp', 'cmdstan18-*', 'output-*.csv', 'temp.data.R'};
for i = 1:numel(pats)
    files = dir(fullfile(baseDir, pats{i}));
    for j = 1:numel(files)
        if ~files(j).isdir
            delete(fullfile(baseDir, files(j).name));
        end
    end
end
end
