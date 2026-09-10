classdef TestGammaMarginalRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the MARGINAL observation-variance term
    % (estimand #2) of the Gamma / scaled-chi-square family. Unlike the other gamma
    % recovery tests (which simulate a CONSTANT nu, s_v = 0, where #2 == #1), this
    % test simulates a person-VARYING dispersion, log nu_p = beta + v_p, with v_p
    % correlated with the person mean effect. It therefore exercises the new
    % extractor plumbing (reading the person log-nu SD and the mean<->log-nu
    % correlation from the fit) and confirms the estimator recovers the #2 truth,
    % not the population-nu #1 plug-in.
    %
    % Ground truth is the #2 closed form (validated offline, without CmdStan, in
    % TestGammaDispersionMarginalization). The true generalizability coefficient is
    % computed from the generative parameters via that closed form; the fit must
    % cover it. A large, well-powered design (many persons x many trials) keeps the
    % person-dispersion parameters identified so the recovery is a clean test of
    % the estimand, not of weak identifiability.
    %
    % CmdStan-gated (skips Incomplete, not Failed, when CmdStan/MatlabStan absent).

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
        function testOneFacetMarginalRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            % ---- Known truth (log expected-score scale) ----
            % Strong person-dispersion heterogeneity (s_v = 0.6) so the marginal
            % correction is large (F ~ 1.38: the #2 residual is ~17% above #1), and
            % a well-powered design so s_v / rho are identified and the recovery is
            % a clean estimand test rather than an identifiability stress test.
            alpha = log(5.5);   % log-mean grand intercept
            s_p   = 0.40;       % person log-mean SD
            s_i   = 0.20;       % trial  log-mean SD
            beta  = log(18);    % population log-nu intercept (nu_pop = 18)
            s_v   = 0.60;       % person log-nu SD (dispersion heterogeneity)
            rho   = -0.30;      % person mean <-> log-nu correlation
            nu_pop = exp(beta);
            cov_pv = rho * s_p * s_v;
            F = exp(0.5 * s_v^2 - 2 * cov_pv);   % ~1.38

            nsub = 100; ntrl = 40;
            dataTbl = localSimulateGammaSv(nsub, ntrl, alpha, s_p, s_i, beta, s_v, rho);

            % ---- Closed-form observed-scale truth under #2 and #1 ----
            vc2 = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_pop, s_v, cov_pv);  % marginal (#2)
            vc1 = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_pop);               % population-nu (#1)
            trueG2 = localSingCoeff(sqrt(vc2.sigma_p2), sqrt(vc2.sigma_pi_e2), ...
                sqrt(vc2.sigma_i2), ntrl);
            sigE2 = sqrt(vc2.sigma_pi_e2);
            sigE1 = sqrt(vc1.sigma_pi_e2);
            % At high reliability the estimand difference lives in the RESIDUAL, not
            % the coefficient (both #1 and #2 coefficients here round to ~0.97). So
            % the discriminating checks below are at the residual / dispersion level.
            testCase.verifyGreaterThan(sigE2 / sigE1 - 1, 0.05, ...
                'Chosen truth makes the #2 and #1 residuals nearly identical.');

            % ---- Fit ----
            rel = localRunSingGamma(testCase, dataTbl, 'recov_gamma_marg');

            testCase.verifyEqual(rel.analysis, 'ic');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');
            for f = {'sig_u','sig_trl','sig_e'}
                testCase.verifyTrue(all(isfinite(rel.out.(f{1})(:))), ...
                    sprintf('Non-finite draws in %s', f{1}));
            end

            svHat  = mean(rel.out.gamma.s_v(:));
            rhoHat = mean(rel.out.gamma.rho_pv(:));
            fHat   = mean(rel.out.gamma.disp_factor(:));
            sigEhat = mean(rel.out.sig_e(:));
            [ll, ~, ul] = psyrat_rel_sing('gcoeff', 2, 'metric', 'global', ...
                'bp', rel.out.sig_u(:), ...
                'wp', sqrt(rel.out.sig_e(:).^2 + rel.out.sig_trl(:).^2), ...
                'i',  rel.out.sig_trl(:), 'obs', ntrl, 'CI', 0.95);
            fprintf(['[recovery] s_v hat=%.3f (true %.2f)  rho hat=%.3f (true %.2f)  ', ...
                'F hat=%.3f (true %.3f)\n           sigE hat=%.3f  #2=%.3f  #1=%.3f  ', ...
                'G95%%CI=[%.3f,%.3f]  trueG2=%.3f\n'], ...
                svHat, s_v, rhoHat, rho, fHat, F, sigEhat, sigE2, sigE1, ll, ul, trueG2);

            % (1) IDENTIFIABILITY: the fit recovers the person-dispersion SD (the
            % param the population-nu #1 conversion discards entirely).
            testCase.verifyGreaterThan(svHat, 0.40, sprintf( ...
                'Recovered person log-nu SD %.3f far below the true %.2f.', svHat, s_v));
            testCase.verifyLessThan(svHat, 0.85);
            % rho is a weaker (second-order) parameter; require only that it is not
            % strongly positive (the true value is negative).
            testCase.verifyLessThan(rhoHat, 0.10);

            % (2) ESTIMAND ACTIVE: the extractor propagated (s_v, cov_pv), so the
            % implied correction factor sits materially above the #1 value of 1.0.
            testCase.verifyGreaterThan(fHat, 1.12, sprintf( ...
                'Recovered dispersion factor %.3f not materially above 1 (true %.3f).', ...
                fHat, F));

            % (3) DISCRIMINATOR: the recovered residual tracks the #2 truth, NOT the
            % #1 plug-in truth. If the extractor still used the population nu, sig_e
            % would sit near sigE1.
            testCase.verifyLessThan(abs(sigEhat - sigE2), abs(sigEhat - sigE1), ...
                sprintf(['Recovered residual %.4f closer to the #1 truth %.4f than ', ...
                'the #2 truth %.4f - extractor may not be applying the margin.'], ...
                sigEhat, sigE1, sigE2));

            % (4) RECOVERY SANITY: the #2 true generalizability coefficient falls
            % inside the recovered 95%% CI.
            testCase.verifyGreaterThanOrEqual(trueG2, ll, sprintf( ...
                'True #2 generalizability %.3f below recovered 95%% CI lower %.3f.', ...
                trueG2, ll));
            testCase.verifyLessThanOrEqual(trueG2, ul, sprintf( ...
                'True #2 generalizability %.3f above recovered 95%% CI upper %.3f.', ...
                trueG2, ul));
        end
    end
end

function tbl = localSimulateGammaSv(nsub, ntrl, alpha, s_p, s_i, beta, s_v, rho)
% One-facet Persons x Trials long table from a log-linked mean (person + crossed
% trial main effects) and a person-VARYING scaled-chi-square dispersion. The
% person (mean, log-nu) block is bivariate normal with correlation rho; the
% continuous-df scaled chi-square Y = gamrnd(nu_p/2, 2*mu/nu_p) has E[Y]=mu and
% Var[Y]=2*mu^2/nu_p exactly.
rng(12345);
u_t = s_i * randn(ntrl, 1);                 % crossed trial main effects
z1 = randn(nsub, 1); z2 = randn(nsub, 1);
u_p = s_p * z1;                             % person mean effect
v_p = s_v * (rho * z1 + sqrt(1 - rho^2) * z2);  % person log-nu effect (corr rho)
nu_p = exp(beta + v_p);

nrow = nsub * ntrl;
ids  = cell(nrow, 1);
meas = zeros(nrow, 1);
r = 0;
for p = 1:nsub
    for t = 1:ntrl
        r = r + 1;
        ids{r}  = sprintf('S%03d', p);
        mu = exp(alpha + u_p(p) + u_t(t));
        meas(r) = gamrnd(nu_p(p) / 2, 2 * mu / nu_p(p));
    end
end
tbl = table(ids, meas, 'VariableNames', {'id', 'meas'});
end

function g = localSingCoeff(bp, sig_e, sig_trl, ntrl)
% One-facet generalizability coefficient at scalar (true) SD-equivalents, using
% the same rel_sing call the production summary uses (wp = total within, i =
% trial SD to remove the trial main effect from relative error).
[~, g, ~] = psyrat_rel_sing('gcoeff', 2, 'metric', 'global', ...
    'bp', bp, 'wp', sqrt(sig_e^2 + sig_trl^2), 'i', sig_trl, ...
    'obs', ntrl, 'CI', 0.95);
end

function rel = localRunSingGamma(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 0, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
p = fullfile(baseDir, 'recov_gamma_marg');
if isfolder(p)
    try
        rmdir(p, 's');
    catch
        warning('tests:cleanupDirFailed', ...
            'Could not remove recovery artifact directory ''%s''.', p);
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
