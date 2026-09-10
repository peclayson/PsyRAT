classdef TestSubjectLevelTrialRecovery < PsyRATTestBase
    % Subject-level (sserr) trial-main-effect recovery (F5).
    %
    % The subject-level error-variance model (sserrvar = 2) now estimates a
    % crossed trial main effect (sig_trl = sigma_i) on the mean, exactly like
    % the population one-facet model, while keeping a PER-SUBJECT residual
    % log-SD. With the trial mean removed, the per-subject residual
    % exp(Intercept_sigma + ind_sd[id]) is the RELATIVE residual
    % sigma_pi,e^2(s), so both a relative (generalizability G_s) and an absolute
    % (dependability phi_s) subject-level coefficient become recoverable.
    %
    % This test simulates one-facet per-subject data with a known crossed trial
    % main effect and heteroscedastic per-subject residuals and confirms:
    %   (a) sig_trl recovers sigma_i and the per-subject residual recovers
    %       sigma_pi,e^2(p) SEPARATELY (the residual does not absorb the trial
    %       main effect);
    %   (b) per-subject phi_s and G_s from psyrat_ssrel obey phi_s <= G_s, with
    %       the gap tracking the sigma_i^2/n gap;
    %   (c) when sigma_i -> 0, phi_s ~ G_s (collapses to the prior behavior).
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
        function testTrialMainEffectAndRelativeAbsoluteSplit(testCase)
            % Known per-subject truth WITH a crossed trial main effect:
            %   meas_pt = mu + alpha_p + beta_t + e_pt
            % alpha_p ~ N(0, sigP^2), beta_t ~ N(0, sigI^2) shared across
            % persons by trial index, and e_pt ~ N(0, sigPIE(p)^2) with a
            % per-subject residual SD. Recovery targets: sig_trl ~ sigI and the
            % per-subject residual ~ sigPIE(p), estimated separately.
            sigP = 3.0; sigI = 1.5; mu = 5;
            nsub = 40; ntrl = 40;
            rng(2468);
            alpha = sigP * randn(nsub, 1);
            beta  = sigI * randn(ntrl, 1);
            % heteroscedastic per-subject residual SDs in [1.2, 3.0]
            sigPIE = 1.2 + 1.8 * rand(nsub, 1);
            [ids, meas] = localSimulate(nsub, ntrl, mu, alpha, beta, sigPIE);
            dataTbl = table(ids, meas, 'VariableNames', {'id', 'meas'});

            rel = localRunSserr(testCase, dataTbl, 'ssrecov_trl');

            % (a) sig_trl is estimated and finite, and recovers sigma_i.
            testCase.verifyTrue(isfield(rel.out, 'sig_trl'));
            testCase.verifyTrue(all(isfinite(rel.out.sig_trl(:))));
            sigThat = mean(rel.out.sig_trl(:));
            testCase.verifyGreaterThanOrEqual(sigThat, 0.8);
            testCase.verifyLessThanOrEqual(sigThat, 2.3);

            % Per-subject relative residual SD = exp(pop_sdlog + ind_sdlog).
            ind_sdlog = cell2mat(rel.out.ind_sdlog{1});
            pop_sdlog = rel.out.pop_sdlog(:, 1);
            resid_sd = mean(exp(pop_sdlog + ind_sdlog), 1)'; % nsub x 1

            % Separation check: the per-subject residual must NOT have absorbed
            % the trial main effect. Lumped truth per subject would be
            % sqrt(sigPIE(p)^2 + sigI^2); the split residual must sit closer to
            % sigPIE(p). Compare mean recovered residual to the two targets.
            relTarget = mean(sigPIE);
            lumpTarget = mean(sqrt(sigPIE.^2 + sigI^2));
            testCase.verifyLessThan(mean(resid_sd), 0.5*(relTarget + lumpTarget));
            % Recovered per-subject residual should track the truth (base
            % corrcoef avoids a Statistics Toolbox dependency).
            cc = corrcoef(resid_sd, sigPIE);
            testCase.verifyGreaterThan(cc(1, 2), 0.4);

            % (b) phi_s (absolute) vs G_s (relative) from psyrat_ssrel.
            [phiT, genT] = localSubjectCoeffs(rel, ntrl);
            % trl_var carries sigma_i^2 (population, shared).
            testCase.verifyGreaterThan(phiT.trl_var(1), 0);
            % Absolute keeps sigma_i^2, so it has MORE error: phi_s <= G_s.
            testCase.verifyTrue(all(phiT.dep_pt <= genT.dep_pt + 1e-9));
            % The gap is real (sigma_i > 0).
            testCase.verifyGreaterThan(mean(genT.dep_pt - phiT.dep_pt), 0);
            % ICC and SEM follow the same direction.
            testCase.verifyTrue(all(phiT.icc_pt <= genT.icc_pt + 1e-9));
            testCase.verifyTrue(all(phiT.sem_pt >= genT.sem_pt - 1e-9));
        end

        function testNoTrialMainEffectCollapsesRelativeToAbsolute(testCase)
            % No crossed trial main effect (sigI = 0): sig_trl collapses toward
            % zero and phi_s ~ G_s (the prior absolute-only behavior).
            sigP = 3.0; mu = 5;
            nsub = 40; ntrl = 40;
            rng(1357);
            alpha = sigP * randn(nsub, 1);
            beta  = zeros(ntrl, 1);
            sigPIE = 1.2 + 1.8 * rand(nsub, 1);
            [ids, meas] = localSimulate(nsub, ntrl, mu, alpha, beta, sigPIE);
            dataTbl = table(ids, meas, 'VariableNames', {'id', 'meas'});

            rel = localRunSserr(testCase, dataTbl, 'ssrecov_notrl');

            sigThat = mean(rel.out.sig_trl(:));
            testCase.verifyLessThan(sigThat, 1.0);

            [phiT, genT] = localSubjectCoeffs(rel, ntrl);
            % With no real trial effect, the absolute and relative coefficients
            % nearly coincide.
            testCase.verifyLessThan(max(abs(phiT.dep_pt - genT.dep_pt)), 0.05);
        end
    end
end

function [phiT, genT] = localSubjectCoeffs(rel, ntrl)
% Build the per-subject absolute (phi_s) and relative (G_s) coefficient tables
% from a fitted sserr result, feeding sig_trl (sigma_i) into psyrat_ssrel.
gcell = rel.out.gro_sds{1};        % [draws x 2] cell
bp = gcell(:, 1);                  % between-person SD draws ([draws x 1] cell)
wp_pop = rel.out.pop_sdlog(:, 1);
wp_ss = rel.out.ind_sdlog{1};      % [draws x nsub] cell
sig_trl = rel.out.sig_trl(:, 1);

nsub = size(cell2mat(wp_ss), 2);
idtable = table((1:nsub)', (1:nsub)', repmat(ntrl, nsub, 1), ...
    'VariableNames', {'id', 'id2', 'trls'});

phiT = psyrat_ssrel('bp', bp, 'wp_pop', wp_pop, 'wp_ss', wp_ss, ...
    'i', sig_trl, 'gcoeff', 1, 'idtable', idtable, 'CI', 0.95);
genT = psyrat_ssrel('bp', bp, 'wp_pop', wp_pop, 'wp_ss', wp_ss, ...
    'i', sig_trl, 'gcoeff', 2, 'idtable', idtable, 'CI', 0.95);
end

function [ids, meas] = localSimulate(nsub, ntrl, mu, alpha, beta, sigPIE)
% Build a balanced one-facet long table from the supplied person (alpha) and
% trial (beta) effects plus per-subject heteroscedastic residual noise.
ids = cell(nsub * ntrl, 1);
meas = zeros(nsub * ntrl, 1);
r = 0;
for pp = 1:nsub
    for tt = 1:ntrl
        r = r + 1;
        ids{r} = sprintf('S%02d', pp);
        meas(r) = mu + alpha(pp) + beta(tt) + sigPIE(pp) * randn;
    end
end
end

function rel = localRunSserr(testCase, dataTbl, tag)
% Run the subject-level error-variance pipeline (sserrvar = 2) with a fixed
% seed.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'sserrvar', 2, ...
    'chains', 2, ...
    'warmup', 500, ...
    'sampling', 500, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
testCase.verifyEqual(rel.analysis, 'ic_sserrvar');
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
for sub = {'ssrecov_trl', 'ssrecov_notrl'}
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
