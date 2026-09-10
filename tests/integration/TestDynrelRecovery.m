classdef TestDynrelRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the dynamic/conditional reliability
    % model (Rast & Clayson, analysis 11; dynrel = 2).
    %
    % The accuracy oracle and the deterministic calc tests only prove
    % "unchanged", not "recovers truth". Here one-facet location-scale data are
    % simulated with a dimension-conditioned log-residual
    %   log sigma_e = alpha + beta_sigma * z + delta_p
    % a crossed trial main effect (sigma_i), correlated person location/scale
    % effects, and a subject-level standardized dimension z. The model must
    % recover sigma_p, sigma_i, the scale slope beta_sigma, and the dynamic
    % G(z)/D(z) surface within tolerance.
    %
    % Like the other CmdStan integration tests, the dependency setup uses
    % assumptions, so this is skipped cleanly (Incomplete, not Failed) when
    % CmdStan/MatlabStan are unavailable (e.g. on CI runners). It is slow
    % (compiles + samples) and lives in tests/integration so the unit lane
    % (runAllUnitTests) excludes it.

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
        function testOneDimensionRecoversTruth(testCase)
            % ---- known truth ----
            rng(20260622, 'twister');
            nsub = 50; ntrl = 30;
            sigP = 4; sigDP = 0.25; sigI = 2; grand = 15;
            sigLog0 = log(3); bMean = 1.5; bSigma = 0.3;

            zraw = randn(nsub, 1);
            zsub = (zraw - mean(zraw)) / std(zraw);   % standardize (re-std = identity)
            nu_p = sigP * randn(nsub, 1);
            del_p = sigDP * randn(nsub, 1);
            beta_t = sigI * randn(ntrl, 1);

            ids = cell(nsub * ntrl, 1);
            meas = zeros(nsub * ntrl, 1);
            dim1 = zeros(nsub * ntrl, 1);
            r = 0;
            for pp = 1:nsub
                for tt = 1:ntrl
                    r = r + 1;
                    mu = grand + bMean * zsub(pp) + nu_p(pp) + beta_t(tt);
                    sg = exp(sigLog0 + bSigma * zsub(pp) + del_p(pp));
                    ids{r} = sprintf('S%02d', pp);
                    meas(r) = mu + sg * randn;
                    dim1(r) = zsub(pp);
                end
            end
            dataTbl = table(ids, meas, dim1, 'VariableNames', {'id', 'meas', 'dim1'});

            rel = localRunDynrel(testCase, dataTbl, 'dynrecov_1dim');

            % ---- posterior means ----
            gro = cell2mat(rel.out.gro_sds{1});      % draws x 2
            log0 = rel.out.pop_sdlog(:, 1);
            bsig = cell2mat(rel.out.b_sigma{1});     % draws x 1
            sigi = rel.out.sig_trl(:, 1);

            % ---- recovery (50 subjects, 30 trials) ----
            testCase.verifyEqual(mean(gro(:, 1)), sigP, 'AbsTol', 0.9);
            testCase.verifyEqual(mean(sigi), sigI, 'AbsTol', 0.7);
            testCase.verifyEqual(mean(log0), sigLog0, 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsig(:, 1)), bSigma, 'AbsTol', 0.13);
            testCase.verifyEqual(mean(gro(:, 2)), sigDP, 'AbsTol', 0.25);

            % ---- G(z)/D(z) recovered vs truth closed form ----
            zg = [-1.5 0 1.5]; nprime = 20;
            out = psyrat_rel_dynrel('sig_p', gro(:, 1), 'sig_log0', log0, ...
                'b_sigma', bsig(:, 1), 'sig_i', sigi, 'ndim', 1, 'z1', zg, ...
                'obs', nprime, 'CI', .95);
            for a = 1:numel(zg)
                sige = exp(sigLog0 + bSigma * zg(a));
                Gtr = sigP^2 / (sigP^2 + sige^2 / nprime);
                Dtr = sigP^2 / (sigP^2 + (sige^2 + sigI^2) / nprime);
                testCase.verifyEqual(out.G.pt(a), Gtr, 'AbsTol', 0.06);
                testCase.verifyEqual(out.D.pt(a), Dtr, 'AbsTol', 0.06);
            end
        end
    end
end

function rel = localRunDynrel(testCase, dataTbl, tag)
% Run the dynamic-reliability pipeline (dynrel = 2) with a fixed seed.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'dynrel', 2, ...
    'chains', 2, ...
    'warmup', 600, ...
    'sampling', 600, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
testCase.verifyEqual(rel.analysis, 'ic_dynrel');
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
% Remove the per-run output directory and any stray CmdStan build/output files.
p = fullfile(baseDir, 'dynrecov_1dim');
if isfolder(p)
    try
        rmdir(p, 's');
    catch
        warning('tests:cleanupDirFailed', ...
            ['Could not remove recovery artifact directory ''%s''. ', ...
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
