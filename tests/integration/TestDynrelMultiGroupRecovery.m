classdef TestDynrelMultiGroupRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for MULTI-GROUP dynamic reliability
    % (dynrel = 2 with a group column). Fills the coverage gap for the per-group
    % stratification path and confirms the maintainer-decided estimand
    % (2026-06-24): the subject-level dimension is GRAND-mean centered across the
    % whole sample (one common z scale), then each group is fit as its own
    % stratum. Two groups with different dimension means are simulated from the
    % SAME true location-scale parameters; the model must (a) detect two strata,
    % (b) recover sigma_p and the scale slope beta_sigma in each, and (c) carry a
    % grand-centered z (overall mean ~0 with the between-group offset preserved,
    % i.e. NOT centered within group).
    %
    % Skips cleanly (Incomplete, not Failed) when CmdStan/MatlabStan are
    % unavailable. Slow (compiles + samples per group); lives in
    % tests/integration so the unit lane (runAllUnitTests) excludes it.

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
        function testTwoGroupsGrandCenteredRecoverTruth(testCase)
            % ---- known truth (shared across both groups) ----
            rng(20260624, 'twister');
            nper = 40; ntrl = 30;                       % subjects per group
            sigP = 4; sigDP = 0.25; sigI = 2; grand = 15;
            sigLog0 = log(3); bMean = 1.5; bSigma = 0.3;

            % subject-level raw dimension, offset between groups so grand-mean
            % and within-group centering give measurably different z. Both groups
            % keep ample within-group variation so each scale slope is identified.
            rawA = 40 + 12 * randn(nper, 1);
            rawB = 60 + 12 * randn(nper, 1);
            raw  = [rawA; rawB];
            grp  = [repmat({'A'}, nper, 1); repmat({'B'}, nper, 1)];
            nsub = 2 * nper;

            % GRAND-mean standardize across ALL subjects (the behaviour under
            % test) and generate the data from that z; feed the pipeline the RAW
            % dimension so it must reproduce the same grand standardization.
            zsub  = (raw - mean(raw)) / std(raw);
            nu_p  = sigP  * randn(nsub, 1);
            del_p = sigDP * randn(nsub, 1);

            id = cell(nsub * ntrl, 1); meas = zeros(nsub * ntrl, 1);
            group = cell(nsub * ntrl, 1); dim1 = zeros(nsub * ntrl, 1);
            r = 0;
            for s = 1:nsub
                beta_t = sigI * randn(ntrl, 1);         % fresh trial effects
                for t = 1:ntrl
                    r = r + 1;
                    mu = grand + bMean * zsub(s) + nu_p(s) + beta_t(t);
                    sg = exp(sigLog0 + bSigma * zsub(s) + del_p(s));
                    id{r}    = sprintf('S%03d', s);
                    group{r} = grp{s};
                    meas(r)  = mu + sg * randn;
                    dim1(r)  = raw(s);                  % RAW; pipeline standardizes
                end
            end
            T = table(id, meas, group, dim1, ...
                'VariableNames', {'id', 'meas', 'group', 'dim1'});

            rel = localRunDynrelGrouped(testCase, T, 'dynrecov_multigroup');

            % ---- two strata detected ----
            testCase.assertEqual(rel.analysis, 'ic_dynrel');
            testCase.assertNumElements(rel.out.gro_sds, 2, ...
                'Multi-group dynrel must produce one stratum per group.');

            % ---- grand-mean centering carried through the full pipeline ----
            zAll = rel.data.Dim1_z;
            testCase.verifyEqual(mean(zAll), 0, 'AbsTol', 1e-6, ...
                'Pipeline dimension must be grand-mean centered (overall mean 0).');
            gAll = rel.data.group;
            testCase.verifyLessThan(mean(zAll(strcmp(gAll, 'A'))), -0.3, ...
                'Low-mean group must stay below 0 on the common scale.');
            testCase.verifyGreaterThan(mean(zAll(strcmp(gAll, 'B'))), 0.3, ...
                'High-mean group must stay above 0 on the common scale.');

            % ---- per-group recovery of the universe-score SD and scale slope ----
            for s = 1:2
                gro  = cell2mat(rel.out.gro_sds{s});    % draws x 2
                bsig = cell2mat(rel.out.b_sigma{s});    % draws x 1
                lbl  = char(string(rel.out.labels{s}));
                testCase.verifyEqual(mean(gro(:, 1)), sigP, 'AbsTol', 1.3, ...
                    sprintf('sigma_p not recovered in group %s.', lbl));
                testCase.verifyEqual(mean(bsig(:, 1)), bSigma, 'AbsTol', 0.3, ...
                    sprintf('scale slope not recovered in group %s.', lbl));
            end
        end
    end
end

function rel = localRunDynrelGrouped(testCase, dataTbl, tag)
% Run the dynamic-reliability pipeline (dynrel = 2) on grouped data; the group
% column is auto-detected by psyrat_computevarcomp (ngroup > 0 -> per-stratum).
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
p = fullfile(baseDir, 'dynrecov_multigroup');
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
