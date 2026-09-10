classdef TestSeededEstimationSmoke < PsyRATTestBase
    % Fast seeded one-facet CmdStan smoke test.
    %
    % This test runs the real estimation pipeline on a tiny simulated
    % one-facet dataset with known variance components and checks two things:
    %   1. Recovery: the posterior means of the between-person SD (sig_u) and
    %      the residual SD (sig_e) land near the simulated truth, within a
    %      deliberately loose tolerance.
    %   2. Reproducibility: two runs with the same seed produce bit-identical
    %      posterior draws (maxAbsDiff = 0). This locks in the single-occasion
    %      seed fix that made the basic dependability analysis reproducible.
    %
    % Like the other CmdStan integration tests, the dependency setup uses
    % assumptions, so the test is skipped cleanly (reported Incomplete, not
    % Failed) whenever CmdStan/MatlabStan are unavailable, e.g. on CI runners
    % where CmdStan is not installed.

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
        function testSeededOneFacetRecoveryAndReproducibility(testCase)
            % Known one-facet truth: meas_ij = mu + u_i + e_ij, with
            % u_i ~ N(0, sigU^2) (between person) and e_ij ~ N(0, sigE^2)
            % (within person). Recovery targets are sigU and sigE.
            sigU = 2.0; sigE = 1.0; mu = 0;
            nsub = 25; ntrl = 8;

            % Fix the RNG so the simulated dataset is identical on every run.
            % Combined with the fixed estimation seed below, the whole test is
            % deterministic, which is what the reproducibility check requires.
            rng(42);
            ids  = repelem((1:nsub)', ntrl);
            u    = sigU * randn(nsub, 1);
            e    = sigE * randn(nsub * ntrl, 1);
            meas = mu + u(ids) + e;
            dataTbl = table(ids, meas, 'VariableNames', {'id', 'meas'});

            % Two separate output directories so the two runs never collide on
            % CmdStan build/output files. Teardown removes both.
            baseDir = localCmdStanOutputDir(testCase.projectRoot());
            outA = fullfile(baseDir, 'smoke_runA');
            outB = fullfile(baseDir, 'smoke_runB');
            if ~exist(outA, 'dir'), mkdir(outA); end
            if ~exist(outB, 'dir'), mkdir(outB); end

            seed = 12345;
            nsamp = 250;
            runOnce = @(outDir) psyrat_computevarcomp( ...
                'data', dataTbl, ...
                'chains', 1, ...
                'warmup', 250, ...
                'sampling', nsamp, ...
                'seed', seed, ...
                'verbose', 1, ...
                'showgui', 1, ...
                'cmdstanoutdir', outDir);

            relA = runOnce(outA);
            relB = runOnce(outB);

            % The pipeline produced a single-occasion (internal-consistency)
            % fit with the expected variance-component fields.
            testCase.verifyEqual(relA.analysis, 'ic');
            testCase.verifyTrue(isfield(relA, 'out'));
            for f = {'mu', 'sig_u', 'sig_e'}
                testCase.verifyTrue(isfield(relA.out, f{1}));
                testCase.verifyFalse(isempty(relA.out.(f{1})));
            end
            testCase.verifyEqual(numel(relA.out.sig_u), nsamp);

            % Recovery within a deliberately loose tolerance. The bounds are
            % wide enough to tolerate CmdStan/toolchain version drift while
            % still catching gross pipeline breakage (zeros, NaNs, wildly
            % wrong estimates). Truth: sigU = 2.0, sigE = 1.0.
            testCase.verifyTrue(all(isfinite(relA.out.sig_u)));
            testCase.verifyTrue(all(isfinite(relA.out.sig_e)));
            sigUhat = mean(relA.out.sig_u);
            sigEhat = mean(relA.out.sig_e);
            testCase.verifyGreaterThanOrEqual(sigUhat, 1.0);
            testCase.verifyLessThanOrEqual(sigUhat, 3.0);
            testCase.verifyGreaterThanOrEqual(sigEhat, 0.5);
            testCase.verifyLessThanOrEqual(sigEhat, 1.8);

            % Reproducibility: two identical-seed runs are bit-identical.
            % This is the regression lock for the single-occasion seed fix.
            testCase.verifyEqual(relB.out.sig_u, relA.out.sig_u);
            testCase.verifyEqual(relB.out.sig_e, relA.out.sig_e);
            testCase.verifyEqual(relB.out.mu, relA.out.mu);
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
% Remove the two per-run output directories and any stray CmdStan build or
% output files left in the base artifacts directory.
for sub = {'smoke_runA', 'smoke_runB'}
    p = fullfile(baseDir, sub{1});
    if isfolder(p)
        try
            rmdir(p, 's');
        catch
            warning('tests:cleanupDirFailed', ...
                ['Could not remove smoke artifact directory ''%s''. ', ...
                'How to fix: close any process using it and delete it manually.'], p);
        end
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

dirs = dir(fullfile(baseDir, 'cmdstan_*'));
dirs = dirs([dirs.isdir]);
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
for i = 1:numel(dirs)
    try
        rmdir(fullfile(baseDir, dirs(i).name), 's');
    catch
        % best-effort cleanup
    end
end
end
