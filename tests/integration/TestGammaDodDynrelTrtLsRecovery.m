classdef TestGammaDodDynrelTrtLsRecovery < matlab.unittest.TestCase
    % Live CmdStan recovery for the TWO-FACET person-specific dynamic
    % NONCONCURRENT difference-of-differences (analysis 29, family='gamma',
    % gammascale=2): trials x occasions, four additional 4-D facet blocks.
    %
    % THIS TEST IS THE LOAD-BEARING VALIDATION GATE FOR ANALYSIS 29 (owner
    % ruling 2026-08-12): the design ships fully supported although the
    % reference bundle's own two-facet runs used a synthetic second occasion
    % (publishable = FALSE there), so the only end-to-end check of the
    % builder -> sampler -> extraction -> stage-3 chain on this analysis is a
    % live planted-truth round trip. It must PASS before merge, not merely
    % exist. Its planted facet components are ALL-DISTINCT across the five
    % facet blocks so a facet-index mis-assignment (the tid/oid/to mapping)
    % fails rather than passing by symmetry - the same closure the
    % modular-cut trt recovery provided for analysis 20.
    %
    % Gating follows the one-facet twin (read its header): the CES-family
    % population G_dod at z = 0 is gated, per-person ordering is gated by
    % rank correlation, per-cell residual contrasts are printed diagnostics.
    % Non-convergence is Incomplete, not Failed.
    %
    % DRAW BUDGET: 4 chains x 1500 warmup / 1000 sampling. The first attempt
    % (2 x 750/750, 2026-08-13) hit max R-hat 1.10 - the two-level occasion
    % facet's scales mix slowly - and was correctly skipped, not passed;
    % this budget is the measured response, not a tuned-to-pass number.

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
                'cmdstan_artifacts', 'recov_gamma_dod_dynrel_trt_ls');
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

        function testDynrelDodTrtRecovery(testCase)
            T = localTruth();
            [dataTbl, truthFx] = localSimulate(T);

            outDir = fullfile(localProjectRoot(), 'tests', 'integration', ...
                'cmdstan_artifacts', 'recov_gamma_dod_dynrel_trt_ls');
            if ~exist(outDir, 'dir'), mkdir(outDir); end
            rel = psyrat_computevarcomp( ...
                'data', dataTbl, ...
                'family', 'gamma', ...
                'gammascale', 2, ...
                'diffest', 3, ...
                'dynrel', 2, ...
                'sserrvar', 2, ...
                'dodmap', {'E1','E2','E3','E4'}, ...
                'chains', 4, ...
                'warmup', 1500, ...
                'sampling', 1000, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1, ...
                'cmdstanoutdir', outDir);

            testCase.verifyEqual(rel.analysis, 'ic_dodiff_dynrel_sserrvar_trt');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.gammascale, 2);
            testCase.verifyEqual(rel.engine, 'cmdstan');

            [maxRhat, nFinite] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('\n[DoD dynrel trt recovery] max R-hat %.4f over %d params (converged = %d)\n', ...
                maxRhat, nFinite, rel.out.conv.converged);
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf( ...
                ['Fit did not converge (max R-hat %.4f); recovery is not ' ...
                'interpretable, so the run is skipped rather than failed.'], maxRhat));

            % ---- Population CES G_dod at z = 0 -----------------------------
            trueG = localTrueSurfaceG(T);
            summ = psyrat_dynrel_summary(rel, 'ngrid', 3, 'obs', T.ntrl, ...
                'nocc', T.nocc, 'reltype', 3);
            st = summ.strata(1);
            [~, iz0] = min(abs(st.z1));
            estG = st.G.pt(iz0);
            fprintf('[DoD dynrel trt recovery] true G(0) = %.4f | recovered G(0) = %.4f [%.4f, %.4f]\n', ...
                trueG, estG, st.G.ll(iz0), st.G.ul(iz0));
            testCase.verifyLessThanOrEqual(abs(estG - trueG), 0.06, ...
                ['The population CES G_dod at z = 0 must recover the ' ...
                'planted truth within 0.06.']);

            % ---- Per-person ordering ---------------------------------------
            ssrel = st.ssrel_table;
            testCase.assertEqual(height(ssrel), T.nsub);
            trueSS = localTruePerPersonG(T, truthFx);
            [~, order] = ismember(ssrel.id, truthFx.ids);
            testCase.assertTrue(all(order > 0));
            rho = localSpearman(trueSS(order), ssrel.gen_pt);
            fprintf('[DoD dynrel trt recovery] per-person G rank correlation = %.3f\n', rho);
            testCase.verifyGreaterThanOrEqual(rho, 0.5);

            % ---- DIAGNOSTIC ONLY: facet SDs (the mapping-sensitive set) ----
            for f = {'sd_trl','sd_occ','sd_tid','sd_oid','sd_to'}
                est = mean(cell2mat(rel.out.(f{1}){1}), 1);
                fprintf('[DoD dynrel trt recovery] DIAGNOSTIC %s recov = [%s]\n', ...
                    f{1}, sprintf('%.3f ', est));
            end
            fprintf('[DoD dynrel trt recovery] DIAGNOSTIC planted trl/occ/tid/oid/to (cell 1) = [%.3f %.3f %.3f %.3f %.3f]\n', ...
                T.sd_i(1), T.sd_o(1), T.sd_pi(1), T.sd_po(1), T.sd_io(1));
        end

    end
end

% ---------------------------------------------------------------------------

function root = localProjectRoot()
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
end

function T = localTruth()
%All-distinct planted truth across cells AND facets (a facet-index swap must
%fail); within-pair asymmetry on the cells as the one-facet twin.
T = struct();
T.events = {'E1','E2','E3','E4'};
T.cvec = [1 -1 -1 1];
T.alpha = [log(5.5) log(4.0) log(5.3) log(5.1)];
T.b_dim = [0.14 -0.06 0.10 -0.03];
T.b_sig = [0.606 1.006 0.640 0.660];
T.b_sig_dim = [-0.08 0.12 -0.05 0.07];
T.sd_p = [0.40 0.38 0.42 0.39];
T.sd_s = [0.30 0.34 0.28 0.32];
%five facet blocks, all-distinct magnitudes so tid/oid/to cannot swap silently
T.sd_i  = [0.20 0.18 0.22 0.19];   %trial
T.sd_o  = [0.12 0.10 0.14 0.11];   %occasion
T.sd_pi = [0.16 0.14 0.18 0.15];   %person x trial
T.sd_po = [0.09 0.07 0.11 0.08];   %person x occasion
T.sd_io = [0.06 0.05 0.08 0.07];   %trial x occasion
Rloc = [1.0 0.7 0.4 0.4; 0.7 1.0 0.4 0.4; 0.4 0.4 1.0 0.7; 0.4 0.4 0.7 1.0];
Rscl = [1.0 0.3 0.15 0.15; 0.3 1.0 0.15 0.15; 0.15 0.15 1.0 0.3; 0.15 0.15 0.3 1.0];
T.R8 = [Rloc, 0.2*eye(4); 0.2*eye(4), Rscl];
R4 = @(r) (1-r)*eye(4) + r*ones(4);
T.R_i = R4(0.25); T.R_o = R4(0.20); T.R_pi = R4(0.15);
T.R_po = R4(0.10); T.R_io = R4(0.12);
T.nsub = 45;
T.ntrl = [12 15 10 13];   %per occasion, per cell (unequal)
T.nocc = 2;
end

function [dataTbl, fx] = localSimulate(T)
rng(20260813, 'twister');
S8 = diag([T.sd_p T.sd_s]) * T.R8 * diag([T.sd_p T.sd_s]);
u_p = randn(T.nsub, 8) * chol(S8);
maxtrl = max(T.ntrl);
mk = @(sd, R) randn(maxtrl, 4) * chol(diag(sd) * R * diag(sd));
u_i = mk(T.sd_i, T.R_i);
u_o = randn(T.nocc, 4) * chol(diag(T.sd_o) * T.R_o * diag(T.sd_o));
u_io = zeros(maxtrl, T.nocc, 4);
Sio = diag(T.sd_io) * T.R_io * diag(T.sd_io);
for o = 1:T.nocc
    u_io(:, o, :) = randn(maxtrl, 4) * chol(Sio);
end
dim1 = randn(T.nsub, 1);
dim1 = (dim1 - mean(dim1)) ./ std(dim1);

Spi = diag(T.sd_pi) * T.R_pi * diag(T.sd_pi);
Spo = diag(T.sd_po) * T.R_po * diag(T.sd_po);

ids = {}; evs = {}; tims = {}; meas = []; dimcol = [];
for s = 1:T.nsub
    u_pi = randn(maxtrl, 4) * chol(Spi);      %this person's trial interactions
    u_po = randn(T.nocc, 4) * chol(Spo);      %this person's occasion interactions
    for o = 1:T.nocc
        for q = 1:4
            sigma = exp(T.b_sig(q) + T.b_sig_dim(q)*dim1(s) + u_p(s, 4+q));
            for t = 1:T.ntrl(q)
                lmu = T.alpha(q) + T.b_dim(q)*dim1(s) + u_p(s, q) + ...
                    u_i(t, q) + u_o(o, q) + u_pi(t, q) + u_po(o, q) + ...
                    u_io(t, o, q);
                mu_t = exp(lmu);
                shape = (mu_t / sigma)^2;
                scale = sigma^2 / mu_t;
                y = localGammaRand(shape) * scale;
                ids{end+1,1} = sprintf('s%03d', s); %#ok<AGROW>
                evs{end+1,1} = T.events{q}; %#ok<AGROW>
                tims{end+1,1} = sprintf('time%d', o); %#ok<AGROW>
                meas(end+1,1) = y; %#ok<AGROW>
                dimcol(end+1,1) = dim1(s); %#ok<AGROW>
            end
        end
    end
end
dataTbl = table(ids, meas, evs, tims, dimcol, ...
    'VariableNames', {'id','meas','event','time','dim1'});
fx = struct();
fx.ids = arrayfun(@(s) sprintf('s%03d', s), (1:T.nsub)', 'UniformOutput', false);
fx.dim1 = dim1;
fx.u_p = u_p;
end

function G = localTrueSurfaceG(T)
%Closed-form typical-person CES G_dod at z = 0 through the committed kernels.
r1 = @(M) reshape(M, [1 4 4]);
mats = psyrat_gamma_dod_components_ls('facet','trial_occasion', ...
    'alpha', T.alpha, 'log_sigma', T.b_sig, ...
    'sd_p', T.sd_p, 'cor_p', r1(T.R8(1:4,1:4)), ...
    'sd_trl', T.sd_i, 'cor_trl', r1(T.R_i), ...
    'sd_occ', T.sd_o, 'cor_occ', r1(T.R_o), ...
    'sd_tid', T.sd_pi, 'cor_tid', r1(T.R_pi), ...
    'sd_oid', T.sd_po, 'cor_oid', r1(T.R_po), ...
    'sd_to', T.sd_io, 'cor_to', r1(T.R_io));
G = localCesG(mats, T);
end

function g = localTruePerPersonG(T, fx)
r1 = @(M) reshape(M, [1 4 4]);
n = T.nsub;
g = zeros(n, 1);
for s = 1:n
    z = fx.dim1(s);
    alpha = T.alpha + T.b_dim .* z;
    lsig = T.b_sig + T.b_sig_dim .* z + fx.u_p(s, 5:8);
    mats = psyrat_gamma_dod_components_ls('facet','trial_occasion', ...
        'alpha', alpha, 'log_sigma', lsig, ...
        'sd_p', T.sd_p, 'cor_p', r1(T.R8(1:4,1:4)), ...
        'sd_trl', T.sd_i, 'cor_trl', r1(T.R_i), ...
        'sd_occ', T.sd_o, 'cor_occ', r1(T.R_o), ...
        'sd_tid', T.sd_pi, 'cor_tid', r1(T.R_pi), ...
        'sd_oid', T.sd_po, 'cor_oid', r1(T.R_po), ...
        'sd_to', T.sd_io, 'cor_to', r1(T.R_io));
    g(s) = localCesG(mats, T);
end
end

function G = localCesG(mats, T)
%CES decomposition (reltype 3) at the planted counts.
ni = T.ntrl; no = repmat(T.nocc, 1, 4);
u = psyrat_dod_composite('varcov', mats.p, 'cvec', T.cvec, 'rule', 'unscaled');
tpi = psyrat_dod_composite('varcov', mats.pi, 'cvec', T.cvec, 'rule', 'harmonic', 'n', ni);
tpo = psyrat_dod_composite('varcov', mats.po, 'cvec', T.cvec, 'rule', 'harmonic', 'n', no);
tpio = psyrat_dod_composite('varcov', mats.pio_e, 'cvec', T.cvec, 'rule', 'harmonic', 'n', ni.*no);
G = u.value / (u.value + tpi.value + tpo.value + tpio.value);
end

function x = localGammaRand(shape)
%Marsaglia-Tsang (2000), base MATLAB only (the family recovery tests' idiom).
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
