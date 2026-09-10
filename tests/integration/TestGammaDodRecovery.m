classdef TestGammaDodRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % four-event difference-of-differences design (analysis 9, family='gamma',
    % non-concurrent).
    %
    % Four cells (events) are simulated from a one-facet Persons x Trials scaled-
    % chi-square model. Person effects are CORRELATED across cells (a valid 4x4
    % log-scale correlation matrix), so their observed-score cross-cell covariances
    % are non-zero and the difference-of-differences universe-score variance
    % (c = [1 -1 -1 1]) is meaningful; the trial effects are cell-specific
    % (independent). The known log-scale parameters are mapped to observed-scale
    % components with the same public helpers the extractor uses
    % (psyrat_gamma_varcomps + psyrat_gamma_crosscov) and the true GENERALIZABILITY
    % DoD coefficient is computed in closed form via psyrat_dodiffrel. The test fits
    % the model through psyrat_computevarcomp and confirms the true coefficient
    % lands inside the recovered 95% credible interval.
    %
    % The GENERALIZABILITY coefficient is asserted: G = u / (u + relative error)
    % depends on the person universe-score variance and the per-cell residual, but
    % NOT on the trial main-effect term, so it is well identified and independent of
    % the (weakly identified) cross-cell trial covariances the estimator keeps.
    %
    % Skipped cleanly (Incomplete, not Failed) when CmdStan/MatlabStan are
    % unavailable, like the other CmdStan integration tests.

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
        function testDodGeneralizabilityCoefficientRecovered(testCase)
            % ---- Known truth (log expected-score scale), four cells E1..E4 ----
            events = {'E1','E2','E3','E4'};
            alpha  = [log(5.5) log(5.0) log(5.2) log(4.8)];
            sd_p   = [0.40 0.38 0.42 0.39];   % person log-SD per cell
            sd_i   = [0.20 0.18 0.22 0.19];   % trial  log-SD per cell
            nu     = 18;                      % population df (integer for the draw)
            % Person cross-cell log-scale correlation: two difference pairs (1,2) &
            % (3,4) with high within-pair (0.7) and moderate cross-pair (0.4)
            % correlation (PSD; eigenvalues 2.5, 0.9, 0.3, 0.3).
            R = [1.0 0.7 0.4 0.4;
                 0.7 1.0 0.4 0.4;
                 0.4 0.4 1.0 0.7;
                 0.4 0.4 0.7 1.0];

            nsub = 40; ntrl = 12;
            cvec = [1 -1 -1 1];

            dataTbl = localSimulateGammaDod(nsub, ntrl, events, alpha, sd_p, sd_i, nu, R);

            % ---- Fit ----
            rel = localRunDodGamma(testCase, dataTbl, 'recov_gamma_dod');

            testCase.verifyEqual(rel.analysis, 'ic_dodiff');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % Estimator cell order (id_varcov cells follow rel.events).
            evOrder = cellstr(string(rel.events(:)));
            testCase.assertEqual(numel(evOrder), 4);
            idx = zeros(1,4);
            for k = 1:4
                idx(k) = find(strcmp(evOrder{k}, events), 1);
            end
            testCase.assertFalse(any(idx == 0), 'rel.events did not match E1..E4.');

            % ---- Closed-form observed-scale truth, assembled in estimator order --
            [bp_t, bt_t, er_t] = localTruthArrays(alpha, sd_p, sd_i, nu, R, idx);
            trueG = psyrat_dodiffrel('bp', bp_t, 'bt', bt_t, 'er_var', er_t, ...
                'obs', ntrl*ones(1,4), 'est', 'gen', 'ci', .95, 'cvec', cvec);

            % ---- Estimated draws (already in estimator/event order) ----
            est_bp = localCell2Arr(rel.out.id_varcov{1});
            est_bt = localCell2Arr(rel.out.trl_varcov{1});
            est_er = localCell2Arr(rel.out.b_sigma{1});
            testCase.verifyTrue(all(isfinite(est_bp(:))) && all(isfinite(est_bt(:))) ...
                && all(isfinite(est_er(:))), 'Non-finite DoD draws.');

            estG = psyrat_dodiffrel('bp', est_bp, 'bt', est_bt, 'er_var', est_er, ...
                'obs', ntrl*ones(1,4), 'est', 'gen', 'ci', .95, 'cvec', cvec);

            fprintf(['\n[DoD gamma recovery] true G = %.4f | recovered G = ', ...
                '%.4f [%.4f, %.4f]\n'], trueG.pt, estG.pt, estG.ll, estG.ul);

            % Headline recovery: the true generalizability DoD coefficient falls in
            % the recovered 95% CI.
            testCase.verifyGreaterThanOrEqual(trueG.pt, estG.ll, sprintf( ...
                ['True DoD generalizability %.3f below the recovered 95%% CI ', ...
                'lower limit %.3f.'], trueG.pt, estG.ll));
            testCase.verifyLessThanOrEqual(trueG.pt, estG.ul, sprintf( ...
                ['True DoD generalizability %.3f above the recovered 95%% CI ', ...
                'upper limit %.3f.'], trueG.pt, estG.ul));
        end
    end
end

function [bp, bt, er] = localTruthArrays(alpha, sd_p, sd_i, nu, R, idx)
% Assemble the observed-scale truth (person 4x4, trial 4x4, per-cell residual
% log-SD) in the estimator's cell order idx, using the same public helpers the
% extractor uses. Trial cross-covariances are 0 (cell-specific trials); person
% cross-covariances come from R.
bp = zeros(1,4,4); bt = zeros(1,4,4); er = zeros(1,4);
for k = 1:4
    c = idx(k);
    vc = psyrat_gamma_varcomps(alpha(c), sd_p(c), sd_i(c), nu);
    bp(1,k,k) = vc.sigma_p2;
    bt(1,k,k) = vc.sigma_i2;
    er(1,k)   = 0.5 * log(vc.sigma_pi_e2);
end
for a = 1:4
    for b = a+1:4
        ca = idx(a); cb = idx(b);
        cc = psyrat_gamma_crosscov(alpha(ca), sd_p(ca), sd_i(ca), ...
            alpha(cb), sd_p(cb), sd_i(cb), R(ca,cb), 0);
        bp(1,a,b) = cc.cov_p_obs; bp(1,b,a) = cc.cov_p_obs;
        % trial cross-cov is 0 (independent trials); bt off-diagonals stay 0.
    end
end
end

function a = localCell2Arr(c)
% REL.out DoD fields store a raw numeric array wrapped in a cell; unwrap (and
% tolerate the {num2cell(...)} form via cell2mat).
if iscell(c)
    a = cell2mat(c);
else
    a = c;
end
end

function tbl = localSimulateGammaDod(nsub, ntrl, events, alpha, sd_p, sd_i, nu, R)
% Balanced four-cell one-facet Persons x Trials scaled-chi-square long table.
% Person effects are correlated across cells (correlation R, per-cell SD sd_p);
% trial effects are cell-specific (independent). Base-MATLAB randn only;
% chi-square(nu) is a sum of nu (integer) squared standard normals.
rng(12345);

% Correlated person effects: columns with correlation R, then scale per cell.
U = chol(R);                       % R = U'*U
uP = (randn(nsub,4) * U) .* sd_p;  % cov(uP(:,j),uP(:,k)) = R(j,k)*sd_p(j)*sd_p(k)

% Cell-specific (independent) trial effects.
uT = cell(1,4);
for c = 1:4
    uT{c} = sd_i(c) * randn(ntrl,1);
end

nrow = nsub * 4 * ntrl;
ids = cell(nrow,1); evs = cell(nrow,1); meas = zeros(nrow,1);
r = 0;
for c = 1:4
    for p = 1:nsub
        for t = 1:ntrl
            r = r + 1;
            ids{r} = sprintf('S%02d', p);
            evs{r} = events{c};
            logmu = alpha(c) + uP(p,c) + uT{c}(t);
            mu = exp(logmu);
            chi2 = sum(randn(nu,1).^2);
            meas(r) = (mu / nu) * chi2;
        end
    end
end
tbl = table(ids, meas, evs, 'VariableNames', {'id','meas','event'});
end

function rel = localRunDodGamma(testCase, dataTbl, tag)
% Fit the four-event difference-of-differences pipeline under the Gamma family,
% fixed seed. diffest=3 selects DoD; dodmap orders the contrast (E1-E2)-(E3-E4).
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    ... %pin log-nu: this design is movable, and since rulings 33-36 the
    ... %engine default is location-scale. Without this the test would
    ... %silently change estimand and check log-nu targets against an
    ... %LS fit. The location-scale twin has its own *LsRecovery test.
    'gammascale', 1, ...
    'diffest', 3, ...
    'dodmap', {'E1','E2','E3','E4'}, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 1, ...
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
p = fullfile(baseDir, 'recov_gamma_dod');
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
