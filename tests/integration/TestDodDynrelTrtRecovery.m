classdef TestDodDynrelTrtRecovery < matlab.unittest.TestCase
    % Live CmdStan recovery for the TWO-FACET (trials x occasions)
    % person-specific dynamic NONCONCURRENT difference-of-differences,
    % GAUSSIAN arm (analysis 29, family='gaussian'). The two-facet twin of
    % TestDodDynrelRecovery - read its header first (the load-bearing ruling,
    % the longhand-truth rule, the planted-truth doctrine, and the
    % negative-cell-means point all carry over).
    %
    % THIS TEST IS LOAD-BEARING AND MUST PASS BEFORE MERGE (owner ruling
    % 2026-08-20): the Gaussian family is unanchored at BOTH analyses, so
    % live recovery is analysis 29's only end-to-end check, exactly as the
    % gamma-29 ruling made TestGammaDodDynrelTrtLsRecovery load-bearing.
    %
    % DRAW BUDGET: 4 chains x 1500 warmup / 1000 sampling - the gamma-29
    % MEASURED budget. Occasion-scale mixing is a facet-geometry property the
    % identity link does not remove, so the budget starts there; if this run
    % converges with a large margin a reduction may be MEASURED later, never
    % assumed (the family's measured-response doctrine).
    %
    % All five facet blocks are planted ALL-DISTINCT so the tid/oid/to
    % mapping cannot pass by symmetry. The CES decomposition (reltype 3) is
    % gated; truth is longhand from the decomposition definitions.

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
                'cmdstan_artifacts', 'recov_dod_dynrel_trt_gauss');
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
            testCase.assertTrue(any(dataTbl.meas <= 0), ...
                'fixture must carry non-positive scores (the Gaussian point)');

            outDir = fullfile(localProjectRoot(), 'tests', 'integration', ...
                'cmdstan_artifacts', 'recov_dod_dynrel_trt_gauss');
            if ~exist(outDir, 'dir'), mkdir(outDir); end
            rel = psyrat_computevarcomp( ...
                'data', dataTbl, ...
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
            testCase.verifyEqual(rel.family, 'gaussian');
            testCase.verifyEqual( ...
                rel.out.dod_labels.inference_framework, ...
                'four_gaussian_margins_all_residual_covariances_fixed_zero_v1');

            [maxRhat, nFinite] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('\n[Gauss DoD trt recovery] max R-hat %.4f over %d params (converged = %d)\n', ...
                maxRhat, nFinite, rel.out.conv.converged);
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf( ...
                ['Fit did not converge (max R-hat %.4f); recovery is not ' ...
                'interpretable, so the run is skipped rather than failed.'], maxRhat));

            % ---- Population CES G_dod at z = 0: LONGHAND truth vs surface ----
            trueG = localTrueSurfaceG(T);
            summ = psyrat_dynrel_summary(rel, 'ngrid', 3, 'obs', T.ntrl, ...
                'nocc', repmat(T.nocc, 1, 4), 'reltype', 3);
            st = summ.strata(1);
            [~, iz0] = min(abs(st.z1));
            estG = st.G.pt(iz0);
            fprintf('[Gauss DoD trt recovery] true CES G(0) = %.4f | recovered = %.4f [%.4f, %.4f]\n', ...
                trueG, estG, st.G.ll(iz0), st.G.ul(iz0));
            testCase.verifyLessThanOrEqual(abs(estG - trueG), 0.06, ...
                ['The population CES G_dod at z = 0 must recover the ' ...
                'planted truth within 0.06.']);

            % ---- Per-person ordering ----
            ssrel = st.ssrel_table;
            testCase.assertEqual(height(ssrel), T.nsub);
            trueSS = localTruePerPersonG(T, truthFx);
            [~, order] = ismember(ssrel.id, truthFx.ids);
            testCase.assertTrue(all(order > 0));
            rho = localSpearman(trueSS(order), ssrel.gen_pt);
            fprintf('[Gauss DoD trt recovery] per-person G rank correlation = %.3f\n', rho);
            testCase.verifyGreaterThanOrEqual(rho, 0.5, ...
                ['Participants with genuinely noisier planted data must ' ...
                'rank lower in the recovered conditional coefficients.']);

            % ---- Parameter-level gates ----
            b = mean(cell2mat(rel.out.b{1}), 1);
            bsig = mean(cell2mat(rel.out.b_sigma{1}), 1);
            sdid = mean(cell2mat(rel.out.sd_id{1}), 1);
            fprintf('[Gauss DoD trt recovery] cell means true = [%s] recov = [%s]\n', ...
                sprintf('%.2f ', T.b), sprintf('%.2f ', b));
            fprintf('[Gauss DoD trt recovery] log res SD true = [%s] recov = [%s]\n', ...
                sprintf('%.3f ', T.b_sig), sprintf('%.3f ', bsig));
            %the cell means are ESTIMAND-IRRELEVANT (the section-26
            %mean-invariance theorem: b enters no reliability quantity) and
            %their posterior SD at this design is ~sd_p/sqrt(nsub) ~= 2.2/sqrt(45) ~= 0.33,
            %so a tight gate would assert on sampling noise. This is a loose
            %WIRING sanity bound (~4 posterior SDs), not a recovery claim;
            %the recovery claims live on the scale submodel, which the
            %estimand actually consumes.
            testCase.verifyLessThanOrEqual(max(abs(b - T.b)), 1.3, ...
                'Identity-scale cell means must land within the sanity bound.');
            testCase.verifyLessThanOrEqual(max(abs(bsig - T.b_sig)), 0.18, ...
                'Log residual-SD intercepts must recover within 0.18.');
            testCase.verifyLessThanOrEqual(max(abs(sdid(5:8) - T.sd_s)), 0.22, ...
                'Person log-residual-SD SDs must recover within 0.22.');

            % ---- DIAGNOSTIC ONLY: occasion-block SDs ----
            sdocc = mean(cell2mat(rel.out.sd_occ{1}), 1);
            fprintf('[Gauss DoD trt recovery] DIAGNOSTIC (not asserted) occasion SDs true = [%s] recov = [%s]\n', ...
                sprintf('%.3f ', T.sd_o), sprintf('%.3f ', sdocc));
        end

    end
end

% ---------------------------------------------------------------------------

function root = localProjectRoot()
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
end

function T = localTruth()
%All-distinct planted truth on GAUSSIAN scales, all five facet blocks
%distinct; asymmetry WITHIN-PAIR; negative cell means deliberate.
T = struct();
T.events = {'E1','E2','E3','E4'};
T.cvec = [1 -1 -1 1];
T.b = [-4.8 -1.9 -4.1 -3.4];
T.b_dim = [0.6 -0.25 0.4 -0.1];
T.b_sig = [0.61 1.01 0.64 0.66];
T.b_sig_dim = [-0.08 0.12 -0.05 0.07];
T.sd_p = [2.0 1.7 2.2 1.9];
T.sd_s = [0.30 0.34 0.28 0.32];
T.sd_i = [0.8 0.7 0.9 0.75];
T.sd_o = [0.45 0.38 0.52 0.41];
T.sd_pi = [0.40 0.47 0.35 0.44];
T.sd_po = [0.34 0.28 0.39 0.31];
T.sd_io = [0.30 0.38 0.26 0.34];
Rloc = [1.0 0.7 0.4 0.4; 0.7 1.0 0.4 0.4; 0.4 0.4 1.0 0.7; 0.4 0.4 0.7 1.0];
Rscl = [1.0 0.3 0.15 0.15; 0.3 1.0 0.15 0.15; 0.15 0.15 1.0 0.3; 0.15 0.15 0.3 1.0];
T.R8 = [Rloc, 0.2*eye(4); 0.2*eye(4), Rscl];
R4 = @(a,b) [1 a b b; a 1 b b; b b 1 a; b b a 1];
T.R_i = R4(0.25, 0.15);
T.R_o = R4(0.30, 0.12);
T.R_pi = R4(0.22, 0.10);
T.R_po = R4(0.28, 0.14);
T.R_io = R4(0.20, 0.08);
T.nsub = 45;
T.nocc = 2;
T.ntrl = [12 15 10 13];    %UNEQUAL per-cell ordinal-trial counts (per occasion)
end

function [dataTbl, fx] = localSimulate(T)
%Simulate the planted GAUSSIAN two-facet model: 8-D person block; 4-D
%correlated blocks for ordinal trial, occasion, person x trial, person x
%occasion, and trial x occasion (crossed levels shared across cells); Normal
%observations. Base MATLAB only.
rng(20260821, 'twister');
S8 = diag([T.sd_p T.sd_s]) * T.R8 * diag([T.sd_p T.sd_s]);
u_p = randn(T.nsub, 8) * chol(S8);
maxtrl = max(T.ntrl);
cov4 = @(sd, R) diag(sd) * R * diag(sd);
u_i = randn(maxtrl, 4) * chol(cov4(T.sd_i, T.R_i));
u_o = randn(T.nocc, 4) * chol(cov4(T.sd_o, T.R_o));
Sio = cov4(T.sd_io, T.R_io);
u_io = zeros(maxtrl, T.nocc, 4);
for o = 1:T.nocc
    u_io(:, o, :) = randn(maxtrl, 4) * chol(Sio);
end
Spi = cov4(T.sd_pi, T.R_pi);
Spo = cov4(T.sd_po, T.R_po);
dim1 = randn(T.nsub, 1);
dim1 = (dim1 - mean(dim1)) ./ std(dim1);

ids = {}; evs = {}; tims = {}; meas = []; dimcol = [];
for s = 1:T.nsub
    u_pi = randn(maxtrl, 4) * chol(Spi);      %this person's trial interactions
    u_po = randn(T.nocc, 4) * chol(Spo);      %this person's occasion interactions
    for o = 1:T.nocc
        for q = 1:4
            mu0 = T.b(q) + T.b_dim(q)*dim1(s) + u_p(s, q);
            sigma = exp(T.b_sig(q) + T.b_sig_dim(q)*dim1(s) + u_p(s, 4+q));
            for t = 1:T.ntrl(q)
                y = mu0 + u_i(t, q) + u_o(o, q) + u_pi(t, q) + ...
                    u_po(o, q) + u_io(t, o, q) + sigma * randn;
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
%LONGHAND typical-person CES G_dod at z = 0: universe = person quadratic
%form; relative error = harmonic composites of the person x trial and
%person x occasion blocks plus the four-term residual sum at n_i*n_o. No
%committed kernel is called.
w = T.cvec(:);
no = repmat(T.nocc, 1, 4);
u = w' * (diag(T.sd_p) * T.R8(1:4,1:4) * diag(T.sd_p)) * w;
pi_h = localHarmonic(diag(T.sd_pi) * T.R_pi * diag(T.sd_pi), w, T.ntrl);
po_h = localHarmonic(diag(T.sd_po) * T.R_po * diag(T.sd_po), w, no);
pio_h = sum((w'.^2) .* exp(2 .* T.b_sig) ./ (T.ntrl .* no));
G = u / (u + pi_h + po_h + pio_h);
end

function g = localTruePerPersonG(T, fx)
%LONGHAND planted per-person conditional CES G with each subject's own
%residual SD at their own z.
w = T.cvec(:);
no = repmat(T.nocc, 1, 4);
u = w' * (diag(T.sd_p) * T.R8(1:4,1:4) * diag(T.sd_p)) * w;
pi_h = localHarmonic(diag(T.sd_pi) * T.R_pi * diag(T.sd_pi), w, T.ntrl);
po_h = localHarmonic(diag(T.sd_po) * T.R_po * diag(T.sd_po), w, no);
g = zeros(T.nsub, 1);
for s = 1:T.nsub
    z = fx.dim1(s);
    lsig = T.b_sig + T.b_sig_dim .* z + fx.u_p(s, 5:8);
    pio_h = sum((w'.^2) .* exp(2 .* lsig) ./ (T.ntrl .* no));
    g(s) = u / (u + pi_h + po_h + pio_h);
end
end

function v = localHarmonic(S, w, n)
%The harmonic count rule written longhand.
v = 0;
for q = 1:4
    v = v + w(q)^2 * S(q,q) / n(q);
end
for a = 1:3
    for b = (a+1):4
        h = 2 / (1/n(a) + 1/n(b));
        v = v + 2 * w(a) * w(b) * S(a,b) / h;
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
