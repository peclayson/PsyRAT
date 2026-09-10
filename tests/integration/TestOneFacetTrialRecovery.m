classdef TestOneFacetTrialRecovery < PsyRATTestBase
    % One-facet trial-main-effect recovery (S10).
    %
    % The one-facet internal-consistency models estimate a crossed trial random
    % effect so the within-person variance is split into a trial main effect
    % (sig_trl = sigma_i) and a person x trial interaction (sig_e = sigma_pi,e),
    % per Rocha 2026 Table 2. This test simulates one-facet data with a known
    % crossed trial main effect and confirms the two components are recovered
    % SEPARATELY (not lumped into the residual), and that sig_trl collapses
    % toward zero when no trial main effect is present.
    %
    % Like the other CmdStan integration tests, the dependency setup uses
    % assumptions, so the test is skipped cleanly (reported Incomplete, not
    % Failed) when CmdStan/MatlabStan are unavailable (e.g. on CI runners).

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
        function testTrialMainEffectRecoveredSeparately(testCase)
            % Known one-facet truth WITH a crossed trial main effect:
            %   meas_pt = mu + alpha_p + beta_t + e_pt
            % with alpha_p ~ N(0, sigP^2), beta_t ~ N(0, sigI^2) shared across
            % persons by trial index, and e_pt ~ N(0, sigPIE^2). Recovery
            % targets: sig_u ~ sigP, sig_trl ~ sigI, sig_e ~ sigPIE.
            sigP = 3.0; sigI = 1.5; sigPIE = 2.0; mu = 5;
            nsub = 50; ntrl = 40;
            rng(12345);
            alpha = sigP * randn(nsub, 1);
            beta  = sigI * randn(ntrl, 1);
            [ids, meas] = localSimulate(nsub, ntrl, mu, alpha, beta, sigPIE);
            dataTbl = table(ids, meas, 'VariableNames', {'id', 'meas'});

            rel = localRunIc(testCase, dataTbl, 'recov_trl');

            for f = {'sig_u', 'sig_trl', 'sig_e'}
                testCase.verifyTrue(isfield(rel.out, f{1}));
                testCase.verifyTrue(all(isfinite(rel.out.(f{1})(:))));
            end
            sigUhat = mean(rel.out.sig_u(:));
            sigThat = mean(rel.out.sig_trl(:));
            sigEhat = mean(rel.out.sig_e(:));

            % Loose bounds: wide enough to tolerate CmdStan/toolchain drift while
            % catching gross breakage. Truth: sigP=3.0, sigI=1.5, sigPIE=2.0.
            testCase.verifyGreaterThanOrEqual(sigUhat, 2.0);
            testCase.verifyLessThanOrEqual(sigUhat, 4.5);
            testCase.verifyGreaterThanOrEqual(sigThat, 0.8);
            testCase.verifyLessThanOrEqual(sigThat, 2.3);
            testCase.verifyGreaterThanOrEqual(sigEhat, 1.4);
            testCase.verifyLessThanOrEqual(sigEhat, 2.6);

            % Key separation check: the residual must NOT have absorbed the trial
            % main effect. The old (lumped) model gave sig_e ~ sqrt(1.5^2+2^2)=2.5;
            % the split model must leave sig_e near the interaction-only value.
            testCase.verifyLessThan(sigEhat, 2.4);
        end

        function testNoTrialMainEffectCollapsesTowardZero(testCase)
            % No crossed trial main effect (sigI = 0): sig_trl should collapse
            % toward zero and sig_e should capture the full within variance.
            sigP = 3.0; sigPIE = 2.0; mu = 5;
            nsub = 50; ntrl = 40;
            rng(777);
            alpha = sigP * randn(nsub, 1);
            beta  = zeros(ntrl, 1);
            [ids, meas] = localSimulate(nsub, ntrl, mu, alpha, beta, sigPIE);
            dataTbl = table(ids, meas, 'VariableNames', {'id', 'meas'});

            rel = localRunIc(testCase, dataTbl, 'recov_notrl');
            sigThat = mean(rel.out.sig_trl(:));
            sigEhat = mean(rel.out.sig_e(:));

            % With no real trial effect, sig_trl reflects only sampling noise in
            % the column means and stays small; sig_e recovers ~sigPIE.
            testCase.verifyLessThan(sigThat, 1.0);
            testCase.verifyGreaterThanOrEqual(sigEhat, 1.5);
            testCase.verifyLessThanOrEqual(sigEhat, 2.5);
        end
    end
end

function [ids, meas] = localSimulate(nsub, ntrl, mu, alpha, beta, sigPIE)
% Build a balanced one-facet long table from the supplied person (alpha) and
% trial (beta) effects plus i.i.d. residual noise.
ids = cell(nsub * ntrl, 1);
meas = zeros(nsub * ntrl, 1);
r = 0;
for pp = 1:nsub
    for tt = 1:ntrl
        r = r + 1;
        ids{r} = sprintf('S%02d', pp);
        meas(r) = mu + alpha(pp) + beta(tt) + sigPIE * randn;
    end
end
end

function rel = localRunIc(testCase, dataTbl, tag)
% Run the single-occasion (internal-consistency) pipeline with a fixed seed.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'chains', 2, ...
    'warmup', 500, ...
    'sampling', 500, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
testCase.verifyEqual(rel.analysis, 'ic');
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
% Remove the per-run output directories and any stray CmdStan build/output
% files left in the base artifacts directory.
for sub = {'recov_trl', 'recov_notrl'}
    p = fullfile(baseDir, sub{1});
    if isfolder(p)
        try
            rmdir(p, 's');
        catch
            warning('tests:cleanupDirFailed', ...
                ['Could not remove recovery artifact directory ''%s''. ', ...
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
