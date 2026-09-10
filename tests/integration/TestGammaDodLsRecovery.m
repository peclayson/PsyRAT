classdef TestGammaDodLsRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the LOCATION-SCALE Gamma /
    % scaled-chi-square four-event difference-of-differences design (analysis 9,
    % family='gamma', gammascale=2, non-concurrent). Increment 7 of the
    % gammascale rollout; the last static design in scope.
    %
    % Sibling of TestGammaDodRecovery, which covers the log-nu twin. Four cells
    % are simulated from a one-facet Persons x Trials model in which the
    % RESIDUAL SD is modeled directly (log sigma_c,p = b_sigma_c + w_p,c) rather
    % than through a dispersion nu, so Var(Y | mu, sigma) = sigma^2 exactly.
    % Person MEAN effects are correlated across cells, so the observed-score
    % cross-cell covariances are non-zero and the difference-of-differences
    % universe-score variance (c = [1 -1 -1 1]) is meaningful.
    %
    % WHAT IS ASSERTED, and why it is the generalizability coefficient.
    % G = u / (u + relative error) depends on the person universe-score variance
    % and the per-cell residual but NOT on the trial main-effect term, so it is
    % well identified and independent of the weakly identified cross-cell trial
    % covariances the estimator keeps.
    %
    % WHAT IS DELIBERATELY NOT ASSERTED. Increment 6 established, by measurement,
    % that the per-cell residual-SD CONTRAST is far noisier than the coefficient:
    % G depends on the SUM of the per-cell residuals, where estimation errors
    % partly cancel, while a contrast depends on their DIFFERENCE, where the same
    % errors compound (~39% attenuation measured there). The contrast is
    % therefore printed as a diagnostic and NOT gated on. Tightening it to a
    % tolerance that one run happens to pass would be fitting the assertion to
    % that run.
    %
    % The DESIGN POINT IS DELIBERATELY ASYMMETRIC, as in increments 5 and 6:
    % cells 2 and 4 carry lower amplitude and higher residual SD, so the amplitude
    % and precision terms push the implied log-nu contrast the same way and the
    % blend the parameterization exists to break is maximal. Expect that same
    % asymmetry to make those cells' components the weakly identified ones.
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

        function testLocationScaleDodGeneralizabilityRecovered(testCase)
            T = localTruth();
            dataTbl = localSimulateGammaDodLs(T);

            rel = localRunDodGamma(testCase, dataTbl, 'recov_gamma_dod_ls', 2);

            testCase.verifyEqual(rel.analysis, 'ic_dodiff');
            testCase.verifyEqual(rel.family,   'gamma');
            testCase.verifyEqual(rel.engine,   'cmdstan');

            % Convergence first: a coefficient recovered from a non-converged fit
            % says nothing. Use the established idiom - psyrat_conv_maxrhat takes
            % REL.out.conv.data (not the REL), and psyrat_checkconv returns a REL
            % whose out.conv.converged is the flag (not a scalar code).
            [maxRhat, nFinite] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('\n[DoD LS recovery] max R-hat %.4f over %d params (converged = %d)\n', ...
                maxRhat, nFinite, rel.out.conv.converged);
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf( ...
                ['Fit did not converge (max R-hat %.4f); recovery is not ' ...
                'interpretable, so the run is skipped rather than failed.'], maxRhat));

            % Estimator cell order (id_varcov cells follow rel.events).
            evOrder = cellstr(string(rel.events(:)));
            testCase.assertEqual(numel(evOrder), 4);
            idx = zeros(1,4);
            for k = 1:4
                idx(k) = find(strcmp(evOrder{k}, T.events), 1);
            end
            testCase.assertFalse(any(idx == 0), 'rel.events did not match the planted cells.');

            % ---- Closed-form observed-scale truth, in estimator order ----
            [bp_t, bt_t, er_t] = localTruthArrays(T, idx);
            trueG = psyrat_dodiffrel('bp', bp_t, 'bt', bt_t, 'er_var', er_t, ...
                'obs', T.ntrl*ones(1,4), 'est', 'gen', 'ci', .95, 'cvec', T.cvec);

            % ---- Estimated draws (already in estimator order) ----
            est_bp = localCell2Arr(rel.out.id_varcov{1});
            est_bt = localCell2Arr(rel.out.trl_varcov{1});
            est_er = localCell2Arr(rel.out.b_sigma{1});
            testCase.verifyTrue(all(isfinite(est_bp(:))) && all(isfinite(est_bt(:))) ...
                && all(isfinite(est_er(:))), 'Non-finite DoD draws.');

            estG = psyrat_dodiffrel('bp', est_bp, 'bt', est_bt, 'er_var', est_er, ...
                'obs', T.ntrl*ones(1,4), 'est', 'gen', 'ci', .95, 'cvec', T.cvec);

            fprintf(['[DoD LS recovery] true G = %.4f | recovered G = %.4f ' ...
                '[%.4f, %.4f]\n'], trueG.pt, estG.pt, estG.ll, estG.ul);

            % Per-cell observed residual: reported, and the SUM is what G uses.
            resTrue = exp(er_t(:)').^2;
            resEst  = mean(exp(est_er).^2, 1);
            fprintf('[DoD LS recovery] observed residual by cell, true  = [%s]\n', ...
                localFmt(resTrue));
            fprintf('[DoD LS recovery] observed residual by cell, recov = [%s]\n', ...
                localFmt(resEst));
            fprintf(['[DoD LS recovery] residual SUM true %.3f vs recovered ' ...
                '%.3f (%.1f%%)\n'], sum(resTrue), sum(resEst), ...
                100*(sum(resEst)/sum(resTrue) - 1));

            % DIAGNOSTIC ONLY, deliberately not asserted - see the class header.
            % The contrast is this parameterization's interpretive payload and is
            % the quantity increment 6 measured as ~39% attenuated while G landed
            % within 0.014. Printing it keeps that visible without pretending a
            % tolerance was derived rather than fitted.
            cTrue = 0.5*sum(log(resTrue) .* T.cvec);
            cEst  = 0.5*sum(log(resEst)  .* T.cvec);
            fprintf(['[DoD LS recovery] DIAGNOSTIC (not asserted) residual-SD ' ...
                'DoD contrast: true %.4f vs recovered %.4f\n'], cTrue, cEst);

            % Headline recovery: the true generalizability coefficient falls in
            % the recovered 95% credible interval.
            testCase.verifyGreaterThanOrEqual(trueG.pt, estG.ll, sprintf( ...
                ['True DoD generalizability %.3f below the recovered 95%% CI ' ...
                'lower limit %.3f.'], trueG.pt, estG.ll));
            testCase.verifyLessThanOrEqual(trueG.pt, estG.ul, sprintf( ...
                ['True DoD generalizability %.3f above the recovered 95%% CI ' ...
                'upper limit %.3f.'], trueG.pt, estG.ul));
        end

        function testMeanSubmodelAgreesAcrossParameterizations(testCase)
            % Paired identity check: fit the SAME simulated data under both
            % parameterizations, same seed, and confirm the claim the golden diff
            % makes structurally - that the mean submodel is untouched by the
            % scale parameterization - holds on live fits too.
            %
            % Errors are correlated across the pair because the fits share data
            % and seed, which is what makes a tight tolerance informative here and
            % would not be justified across independent runs.
            T = localTruth();
            dataTbl = localSimulateGammaDodLs(T);

            relLs = localRunDodGamma(testCase, dataTbl, 'recov_gamma_dod_ls_pairA', 2);
            relNu = localRunDodGamma(testCase, dataTbl, 'recov_gamma_dod_ls_pairB', 1);

            [rhatLs, nLs] = psyrat_conv_maxrhat(relLs.out.conv.data);
            [rhatNu, nNu] = psyrat_conv_maxrhat(relNu.out.conv.data);
            relLs = psyrat_checkconv(relLs);
            relNu = psyrat_checkconv(relNu);
            fprintf(['\n[DoD LS identity] LS max R-hat %.4f over %d (conv %d) | ' ...
                'log-nu max R-hat %.4f over %d (conv %d)\n'], ...
                rhatLs, nLs, relLs.out.conv.converged, ...
                rhatNu, nNu, relNu.out.conv.converged);
            testCase.assumeEqual(relLs.out.conv.converged, 1, ...
                'Location-scale fit did not converge; identity check skipped.');
            testCase.assumeEqual(relNu.out.conv.converged, 1, ...
                'Log-nu fit did not converge; identity check skipped.');

            % (a) The PERSON and TRIAL observed-scale (co)variance matrices are
            % built from the log-mean submodel alone, so the two fits must agree
            % on them. This is the offline invariance claim (which is exact given
            % identical draws) transported to live fits, where it becomes an
            % estimation-agreement check and needs a real tolerance.
            pLs = mean(localCell2Arr(relLs.out.id_varcov{1}), 1);
            pNu = mean(localCell2Arr(relNu.out.id_varcov{1}), 1);
            pLs = squeeze(pLs); pNu = squeeze(pNu);
            relDiff = abs(pLs - pNu) ./ max(abs(pNu), eps);
            fprintf('[DoD LS identity] person varcov max relative diff %.3f\n', ...
                max(relDiff(:)));
            testCase.verifyLessThan(max(relDiff(:)), 0.35, ...
                ['The person (co)variance matrix is a pure log-mean-submodel ' ...
                'quantity and must agree across parameterizations.']);

            % (b) The structural cost contrast, reproduced a fourth time. The
            % log-nu parameterization needs a substantially LARGER person
            % dispersion SD to describe the same data, because its ind_nu must
            % absorb the amplitude term as well: v_p = const + 2*u_p - 2*w_p.
            % Increments 4/5/6 measured 2.9x, ~2.4x and ~2.3x on different
            % designs, so the direction is a property of the parameterization.
            % Direction is ASSERTED, not merely printed.
            dLs = median(relLs.out.gamma_disp.s_sd, 1);
            dNu = median(relNu.out.gamma_disp.s_v,  1);
            fprintf('[DoD LS identity] s_sd by cell = [%s]\n', localFmt(dLs));
            fprintf('[DoD LS identity] s_v  by cell = [%s]  (ratio [%s])\n', ...
                localFmt(dNu), localFmt(dNu ./ dLs));
            testCase.verifyGreaterThan(mean(dNu ./ dLs), 1.5, ...
                ['The log-nu person dispersion SD must exceed the ' ...
                'location-scale one, because it additionally absorbs amplitude.']);

            % (c) The two stored residuals are NOT the same quantity, and the
            % check says so rather than asserting they agree. The log-nu
            % typical-subject residual carries the trial-averaging factor
            % exp(2*s_i^2) that the location-scale residual has no counterpart
            % for, so the ratio is predicted ABOVE 1, not equal to it.
            rLs = mean(exp(localCell2Arr(relLs.out.b_sigma{1})).^2, 1);
            rNu = mean(exp(localCell2Arr(relNu.out.b_sigma{1})).^2, 1);
            predicted = exp(2 .* T.sd_i.^2);
            fprintf(['[DoD LS identity] residual ratio log-nu/LS by cell = ' ...
                '[%s], predicted exp(2*s_i^2) = [%s]\n'], ...
                localFmt(rNu ./ rLs), localFmt(predicted));
            testCase.verifyGreaterThan(min(rNu ./ rLs), 1.0, ...
                ['The log-nu residual carries a trial-averaging factor the ' ...
                'location-scale residual does not, so it must be LARGER. A ' ...
                'ratio at or below 1 would mean the two are being conflated.']);
        end

    end
end

% =========================================================================
function T = localTruth()
% Planted truth for the location-scale four-cell design. ASYMMETRIC on purpose:
% cells 2 and 4 carry the lower amplitude and the higher residual SD, so the
% amplitude and precision terms push the implied log-nu contrast the same way.
%
% THE ASYMMETRY MUST BE WITHIN-PAIR, NOT ALTERNATING, and that is not obvious.
% The contrast is (E1-E2)-(E3-E4), so a (high, low, high, low) pattern - the
% natural generalization of increment 5's two-event asymmetry - CANCELS: it was
% measured at a residual-SD DoD contrast of -0.063, essentially zero, which would
% have made this design blind to the one quantity the parameterization exists to
% separate. The planted values below instead give pair 1 a large within-pair
% difference in BOTH amplitude and residual SD and pair 2 a small one, so the
% differences do not cancel.
%
% Both terms push the implied log-nu contrast the same way, which is the point:
%   amplitude DoD contrast    = (log 5.5 - log 4.0) - (log 5.3 - log 5.1) = 0.280
%   residual-SD DoD contrast  = (0.606 - 1.006) - (0.640 - 0.660)         = -0.380
%   implied log-nu contrast   = 2*0.280 - 2*(-0.380)                      = 1.320
% so a log-nu fit's cell contrasts blend the two maximally here.
T = struct();
T.events = {'E1','E2','E3','E4'};
T.cvec   = [1 -1 -1 1];
T.alpha  = [log(5.5) log(4.0) log(5.3) log(5.1)];  % per-cell log-mean
T.sd_p   = [0.40 0.38 0.42 0.39];                  % person log-mean SD
T.sd_i   = [0.20 0.18 0.22 0.19];                  % trial  log-mean SD
T.b_sig  = [0.606 1.006 0.640 0.660];              % per-cell log residual-SD
T.s_sd   = [0.50 0.50 0.50 0.50];                  % person log-sigma SD
% Person MEAN cross-cell log-scale correlation: the two difference pairs (1,2)
% and (3,4) correlate 0.7 within pair and 0.4 across (PSD; eigenvalues 2.5, 0.9,
% 0.3, 0.3), matching the log-nu sibling's design so the two are comparable.
T.R = [1.0 0.7 0.4 0.4;
       0.7 1.0 0.4 0.4;
       0.4 0.4 1.0 0.7;
       0.4 0.4 0.7 1.0];
T.nsub = 60;
T.ntrl = 20;
end

% =========================================================================
function [bp, bt, er] = localTruthArrays(T, idx)
% Assemble the observed-scale truth (person 4x4, trial 4x4, per-cell residual
% log-SD) in the estimator's cell order, using the same public helpers the
% extractor uses. Trial cross-covariances are 0 (cell-specific trials); person
% cross-covariances come from T.R. Residual cross-covariance is 0 by design
% (non-concurrent), which is structural in psyrat_dodiffrel's per-cell er_var.
bp = zeros(1,4,4); bt = zeros(1,4,4); er = zeros(1,4);
for k = 1:4
    c = idx(k);
    vc = psyrat_gamma_varcomps_ls(T.alpha(c), T.sd_p(c), T.sd_i(c), ...
        T.b_sig(c), T.s_sd(c));
    bp(1,k,k) = vc.sigma_p2;
    bt(1,k,k) = vc.sigma_i2;
    er(1,k)   = 0.5 * log(vc.sigma_pi_e2);
end
for a = 1:4
    for b = a+1:4
        ca = idx(a); cb = idx(b);
        cc = psyrat_gamma_crosscov(T.alpha(ca), T.sd_p(ca), T.sd_i(ca), ...
            T.alpha(cb), T.sd_p(cb), T.sd_i(cb), T.R(ca,cb), 0);
        bp(1,a,b) = cc.cov_p_obs; bp(1,b,a) = cc.cov_p_obs;
        % trial cross-cov is 0 (independent trials); bt off-diagonals stay 0.
    end
end
end

% =========================================================================
function tbl = localSimulateGammaDodLs(T)
% Balanced four-cell one-facet Persons x Trials scaled-chi-square long table
% under the LOCATION-SCALE generative model: each person gets a per-cell mean
% effect AND a per-cell log residual-SD effect, and the observation's df is
% DERIVED as nu = 2*(mu/sigma)^2 rather than fixed.
%
% Base MATLAB only (randn/rand; no Statistics Toolbox, which CI does not
% install). The log-nu sibling can draw chi2(nu) as a sum of nu squared normals
% because its nu is a fixed INTEGER; here nu = 2*(mu/sigma)^2 is derived and
% non-integer, so that route is unavailable.
%
% The draw is EXACT rather than approximate. y = (mu/nu)*chi2(nu) is exactly
% Gamma(shape = nu/2, scale = 2*mu/nu), which localGammaRand samples by
% Marsaglia-Tsang. A Wilson-Hilferty normal approximation was considered and
% rejected: its error is unquantified in the small-nu tail, and nu here is
% person- and trial-varying rather than fixed, so the tail is genuinely visited.
% Approximating the data-generating process would put an unknown bias in the
% ground truth this whole test compares against.
rng(12345);

% Correlated person MEAN effects; independent person SCALE effects (the scale
% block's cross-cell correlation is not the quantity under test here, and
% leaving it at zero keeps the planted truth readable).
U  = chol(T.R);
uP = (randn(T.nsub,4) * U) .* T.sd_p;
wP = randn(T.nsub,4) .* T.s_sd;

% Cell-specific (independent) trial effects on the mean.
uT = cell(1,4);
for c = 1:4
    uT{c} = T.sd_i(c) * randn(T.ntrl,1);
end

nrow = T.nsub * 4 * T.ntrl;
ids = cell(nrow,1); evs = cell(nrow,1); meas = zeros(nrow,1);
r = 0;
for c = 1:4
    for p = 1:T.nsub
        sigma_cp = exp(T.b_sig(c) + wP(p,c));   % this person's residual SD, cell c
        for t = 1:T.ntrl
            r = r + 1;
            ids{r} = sprintf('S%02d', p);
            evs{r} = T.events{c};
            mu  = exp(T.alpha(c) + uP(p,c) + uT{c}(t));
            nu  = 2 * (mu / sigma_cp)^2;        % derived df, non-integer
            % y = (mu/nu)*chi2(nu) == Gamma(nu/2, 2*mu/nu); E[y] = mu and
            % Var(y) = 2*mu^2/nu = sigma_cp^2 by construction.
            meas(r) = localGammaRand(nu/2) * (2 * mu / nu);
        end
    end
end
tbl = table(ids, meas, evs, 'VariableNames', {'id','meas','event'});
end

% =========================================================================
function x = localGammaRand(shape)
% Exact Gamma(shape, scale = 1) draw by the Marsaglia-Tsang (2000) method, using
% only base-MATLAB randn/rand so the test runs on a CI image installed with
% --products=MATLAB. Statistics Toolbox gamrnd would be simpler and is not
% available there.
%
% The method is defined for shape >= 1; for shape < 1 it uses the standard boost
% Gamma(a) = Gamma(a+1) * U^(1/a). The rejection loop is guaranteed to terminate
% (acceptance probability is bounded well away from zero for every shape), and
% the v <= 0 guard rejects the cube-root branch where the transform is invalid.
if shape < 1
    x = localGammaRand(shape + 1) * rand^(1/shape);
    return;
end
d = shape - 1/3;
c = 1 / sqrt(9*d);
while true
    v = -1;
    while v <= 0
        z = randn;
        v = (1 + c*z)^3;
    end
    u = rand;
    if log(u) < 0.5*z^2 + d - d*v + d*log(v)
        x = d * v;
        return;
    end
end
end

% =========================================================================
function a = localCell2Arr(c)
% REL.out DoD fields store a raw numeric array wrapped in a cell; unwrap (and
% tolerate the {num2cell(...)} form via cell2mat).
if iscell(c)
    a = cell2mat(c);
else
    a = c;
end
end

% =========================================================================
function s = localFmt(v)
s = strjoin(arrayfun(@(x) sprintf('%.3f', x), v(:)', 'UniformOutput', false), ', ');
end

% =========================================================================
function rel = localRunDodGamma(testCase, dataTbl, tag, gammascale)
% Fit the four-event difference-of-differences pipeline under the Gamma family
% at the requested scale parameterization, fixed seed. diffest=3 selects DoD;
% dodmap orders the contrast (E1-E2)-(E3-E4).
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'gammascale', gammascale, ...
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

% =========================================================================
function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

% =========================================================================
function localCleanupArtifacts(baseDir)
tags = {'recov_gamma_dod_ls', 'recov_gamma_dod_ls_pairA', 'recov_gamma_dod_ls_pairB'};
for i = 1:numel(tags)
    p = fullfile(baseDir, tags{i});
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
