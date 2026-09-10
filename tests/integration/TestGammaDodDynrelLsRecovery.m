classdef TestGammaDodDynrelLsRecovery < matlab.unittest.TestCase
    % Live CmdStan recovery for the ONE-FACET person-specific dynamic
    % NONCONCURRENT difference-of-differences (analysis 28, family='gamma',
    % gammascale=2): plant known truth, fit with the real pipeline
    % (psyrat_computevarcomp -> psyrat_build_stan_dodiff_dynrel_margins_chisq
    % -> CmdStan -> extraction), then check that the committed stage-3
    % machinery recovers the planted population coefficient and ranks
    % participants correctly.
    %
    % WHAT IS GATED AND WHAT IS PRINTED. The population G_dod at z = 0 is
    % SUM-based and well identified at this design size, so it is gated. The
    % per-person ordering of the conditional coefficients depends on the
    % planted person scale effects, so a rank-correlation gate checks that
    % participants with genuinely noisier data get genuinely lower
    % coefficients. The per-cell residual-SD CONTRASTS are the
    % attenuation-prone quantities the family's earlier increments measured
    % (~39% attenuation at comparable designs); they are PRINTED as
    % diagnostics and deliberately not gated - tightening them to a tolerance
    % one run happens to pass would be fitting the assertion to that run.
    %
    % PLANTED-TRUTH DOCTRINE (the family's): all-distinct values everywhere -
    % per-cell amplitudes, residual SDs, person location and scale SDs, trial
    % SDs - with the asymmetry WITHIN-PAIR (pair 1 carries the large
    % within-pair differences), because an alternating pattern cancels
    % through the (E1-E2)-(E3-E4) contrast and blinds the design (measured in
    % the static DoD increment). Dimension slopes are nonzero on BOTH
    % submodels so the dynamic machinery is exercised, and the per-cell trial
    % counts are UNEQUAL - the property the long-format layout exists for.
    %
    % Gated on CmdStan availability (skips Incomplete when absent);
    % non-convergence is Incomplete, not Failed.

    % Copyright (C) 2016-2026 Peter E. Clayson
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program (gpl.txt). If not, see
    %     <http://www.gnu.org/licenses/>.
    %

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            outdir = fullfile(localProjectRoot(), 'tests', 'integration', ...
                'cmdstan_artifacts', 'recov_gamma_dod_dynrel_ls');
            if isfolder(outdir)
                try
                    rmdir(outdir, 's');
                catch
                    warning('tests:cleanupDirFailed', ...
                        'Could not remove %s.', outdir);
                end
            end
        end
    end

    methods (Test)

        function testDynrelDodRecovery(testCase)
            T = localTruth();
            [dataTbl, truthFx] = localSimulate(T);

            outDir = fullfile(localProjectRoot(), 'tests', 'integration', ...
                'cmdstan_artifacts', 'recov_gamma_dod_dynrel_ls');
            if ~exist(outDir, 'dir'), mkdir(outDir); end
            rel = psyrat_computevarcomp( ...
                'data', dataTbl, ...
                'family', 'gamma', ...
                'gammascale', 2, ...
                'diffest', 3, ...
                'dynrel', 2, ...
                'sserrvar', 2, ...
                'dodmap', {'E1','E2','E3','E4'}, ...
                'chains', 2, ...
                'warmup', 750, ...
                'sampling', 750, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1, ...
                'cmdstanoutdir', outDir);

            testCase.verifyEqual(rel.analysis, 'ic_dodiff_dynrel_sserrvar');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.gammascale, 2);
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % Convergence first: recovery from a non-converged fit says
            % nothing; skip rather than fail.
            [maxRhat, nFinite] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('\n[DoD dynrel recovery] max R-hat %.4f over %d params (converged = %d)\n', ...
                maxRhat, nFinite, rel.out.conv.converged);
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf( ...
                ['Fit did not converge (max R-hat %.4f); recovery is not ' ...
                'interpretable, so the run is skipped rather than failed.'], maxRhat));

            % ---- Population G_dod at z = 0: closed-form truth vs surface ----
            trueG = localTrueSurfaceG(T);
            summ = psyrat_dynrel_summary(rel, 'ngrid', 3, 'obs', T.ntrl);
            st = summ.strata(1);
            [~, iz0] = min(abs(st.z1));
            estG = st.G.pt(iz0);
            fprintf('[DoD dynrel recovery] true G(0) = %.4f | recovered G(0) = %.4f [%.4f, %.4f]\n', ...
                trueG, estG, st.G.ll(iz0), st.G.ul(iz0));
            testCase.verifyLessThanOrEqual(abs(estG - trueG), 0.05, ...
                ['The population G_dod at z = 0 must recover the planted ' ...
                'truth within 0.05.']);

            % ---- Per-person ordering: planted vs recovered conditional G ----
            ssrel = st.ssrel_table;
            testCase.assertEqual(height(ssrel), T.nsub);
            trueSS = localTruePerPersonG(T, truthFx);
            % align: ssrel.id is the subject label s<k>
            [~, order] = ismember(ssrel.id, truthFx.ids);
            testCase.assertTrue(all(order > 0));
            rho = localSpearman(trueSS(order), ssrel.gen_pt);
            fprintf('[DoD dynrel recovery] per-person G rank correlation (true vs recovered) = %.3f\n', rho);
            testCase.verifyGreaterThanOrEqual(rho, 0.5, ...
                ['Participants with genuinely noisier planted data must ' ...
                'rank lower in the recovered conditional coefficients.']);

            % ---- DIAGNOSTIC ONLY: per-cell residual-SD recovery ----
            bsig = cell2mat(rel.out.b_sigma{1});
            fprintf('[DoD dynrel recovery] DIAGNOSTIC (not asserted) log residual SD, true  = [%s]\n', ...
                sprintf('%.3f ', T.b_sig));
            fprintf('[DoD dynrel recovery] DIAGNOSTIC (not asserted) log residual SD, recov = [%s]\n', ...
                sprintf('%.3f ', mean(bsig, 1)));
        end

    end
end

% ---------------------------------------------------------------------------

function root = localProjectRoot()
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
end

function T = localTruth()
%All-distinct planted truth; asymmetry WITHIN-PAIR (see class header).
T = struct();
T.events = {'E1','E2','E3','E4'};
T.cvec = [1 -1 -1 1];
T.alpha = [log(5.5) log(4.0) log(5.3) log(5.1)];   %per-cell log-mean at z = 0
T.b_dim = [0.14 -0.06 0.10 -0.03];                 %per-cell mean slope on dim1
T.b_sig = [0.606 1.006 0.640 0.660];               %per-cell log residual SD at z = 0
T.b_sig_dim = [-0.08 0.12 -0.05 0.07];             %per-cell scale slope on dim1
T.sd_p = [0.40 0.38 0.42 0.39];                    %person location SDs
T.sd_s = [0.30 0.34 0.28 0.32];                    %person log-sigma SDs
T.sd_i = [0.20 0.18 0.22 0.19];                    %trial SDs
%joint 8x8 person correlation: location block within-pair 0.7 / across 0.4
%(the family's design), scale block 0.3/0.15, location-scale 0.2 within cell.
Rloc = [1.0 0.7 0.4 0.4; 0.7 1.0 0.4 0.4; 0.4 0.4 1.0 0.7; 0.4 0.4 0.7 1.0];
Rscl = [1.0 0.3 0.15 0.15; 0.3 1.0 0.15 0.15; 0.15 0.15 1.0 0.3; 0.15 0.15 0.3 1.0];
T.R8 = [Rloc, 0.2*eye(4); 0.2*eye(4), Rscl];
T.R_i = [1.0 0.25 0.15 0.15; 0.25 1.0 0.15 0.15; ...
    0.15 0.15 1.0 0.25; 0.15 0.15 0.25 1.0];       %cross-cell trial correlation
T.nsub = 60;
T.ntrl = [18 22 15 20];                            %UNEQUAL per-cell counts
end

function [dataTbl, fx] = localSimulate(T)
%Simulate the planted model: per-subject 8-D joint person effects, per-ordinal-
%trial 4-D correlated trial effects (shared crossed levels), scaled chi-square
%observations via the Gamma equivalence, and a standardized dimension column.
rng(20260812, 'twister');
S8 = diag([T.sd_p T.sd_s]) * T.R8 * diag([T.sd_p T.sd_s]);
Si = diag(T.sd_i) * T.R_i * diag(T.sd_i);
u_p = randn(T.nsub, 8) * chol(S8);
maxtrl = max(T.ntrl);
u_i = randn(maxtrl, 4) * chol(Si);
%standardized dimension values without Statistics Toolbox zscore (which
%TestBaseMatlabOnly forbids repo-wide)
dim1 = randn(T.nsub, 1);
dim1 = (dim1 - mean(dim1)) ./ std(dim1);

ids = {}; evs = {}; meas = []; dimcol = [];
for s = 1:T.nsub
    for q = 1:4
        mu = exp(T.alpha(q) + T.b_dim(q)*dim1(s) + u_p(s, q));
        sigma = exp(T.b_sig(q) + T.b_sig_dim(q)*dim1(s) + u_p(s, 4+q));
        for t = 1:T.ntrl(q)
            mu_t = mu * exp(u_i(t, q));
            shape = (mu_t / sigma)^2;
            scale = sigma^2 / mu_t;
            y = localGammaRand(shape) * scale;
            ids{end+1,1} = sprintf('s%03d', s); %#ok<AGROW>
            evs{end+1,1} = T.events{q}; %#ok<AGROW>
            meas(end+1,1) = y; %#ok<AGROW>
            dimcol(end+1,1) = dim1(s); %#ok<AGROW>
        end
    end
end
dataTbl = table(ids, meas, evs, dimcol, ...
    'VariableNames', {'id','meas','event','dim1'});
fx = struct();
fx.ids = arrayfun(@(s) sprintf('s%03d', s), (1:T.nsub)', 'UniformOutput', false);
fx.dim1 = dim1;
fx.u_p = u_p;
end

function G = localTrueSurfaceG(T)
%Closed-form typical-person G_dod at z = 0 from the planted values, through
%the committed kernels (single "draw").
mats = psyrat_gamma_dod_components_ls('facet','trial', ...
    'alpha', T.alpha, 'log_sigma', T.b_sig, ...
    'sd_p', T.sd_p, 'cor_p', reshape(T.R8(1:4,1:4), [1 4 4]), ...
    'sd_trl', T.sd_i, 'cor_trl', reshape(T.R_i, [1 4 4]));
u = psyrat_dod_composite('varcov', mats.p, 'cvec', T.cvec, 'rule', 'unscaled');
r = psyrat_dod_composite('varcov', mats.pi_e, 'cvec', T.cvec, ...
    'rule', 'harmonic', 'n', T.ntrl);
G = u.value / (u.value + r.value);
end

function x = localGammaRand(shape)
%Exact Gamma(shape, scale = 1) draw by the Marsaglia-Tsang (2000) method, base
%MATLAB only (the family recovery tests' shared idiom; Statistics Toolbox
%gamrnd is unavailable on the CI image installed with --products=MATLAB).
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

function rho = localSpearman(a, b)
%Spearman rank correlation, base MATLAB (no Statistics Toolbox corr/tiedrank).
%The inputs here are continuous, so ties have probability zero; plain sort
%ranks suffice.
ra = localRanks(a(:));
rb = localRanks(b(:));
ra = ra - mean(ra); rb = rb - mean(rb);
rho = (ra' * rb) / sqrt((ra' * ra) * (rb' * rb));
end

function r = localRanks(x)
[~, order] = sort(x);
r = zeros(size(x));
r(order) = (1:numel(x))';
end

function g = localTruePerPersonG(T, fx)
%Planted per-person conditional G at each subject's own z with their own
%residual SD (their realized scale effect), through the committed kernels.
n = T.nsub;
g = zeros(n, 1);
for s = 1:n
    z = fx.dim1(s);
    alpha = T.alpha + T.b_dim .* z;
    lsig = T.b_sig + T.b_sig_dim .* z + fx.u_p(s, 5:8);
    mats = psyrat_gamma_dod_components_ls('facet','trial', ...
        'alpha', alpha, 'log_sigma', lsig, ...
        'sd_p', T.sd_p, 'cor_p', reshape(T.R8(1:4,1:4), [1 4 4]), ...
        'sd_trl', T.sd_i, 'cor_trl', reshape(T.R_i, [1 4 4]));
    u = psyrat_dod_composite('varcov', mats.p, 'cvec', T.cvec, 'rule', 'unscaled');
    r = psyrat_dod_composite('varcov', mats.pi_e, 'cvec', T.cvec, ...
        'rule', 'harmonic', 'n', T.ntrl);
    g(s) = u.value / (u.value + r.value);
end
end
