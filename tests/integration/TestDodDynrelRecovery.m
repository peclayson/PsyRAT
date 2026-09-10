classdef TestDodDynrelRecovery < matlab.unittest.TestCase
    % Live CmdStan recovery for the ONE-FACET person-specific dynamic
    % NONCONCURRENT difference-of-differences, GAUSSIAN arm (analysis 28,
    % family='gaussian'): plant known truth, fit with the real pipeline
    % (psyrat_computevarcomp -> psyrat_build_stan_dodiff_dynrel_margins ->
    % CmdStan -> extraction), then check that the committed stage-3 machinery
    % recovers the planted population coefficient and ranks participants
    % correctly. The Gaussian class carries no family marker, matching
    % TestSubjDiffDynrelRecovery; the gamma arm's class is
    % TestGammaDodDynrelLsRecovery.
    %
    % THIS TEST IS LOAD-BEARING (owner ruling 2026-08-20, extending the
    % 2026-08-12 analysis-29 precedent): no external reference exists for a
    % Gaussian DoD and none can be minted from the gamma bundle, so live
    % planted-truth recovery is the only end-to-end
    % builder -> sampler -> extraction -> stage-3 check this arm has. It must
    % PASS before the increment merges, not merely run
    % (SCIENTIFIC_FORMULA_AUDIT.md section 26, validation plan of record).
    %
    % WHAT IS GATED AND WHAT IS PRINTED. The population G_dod at z = 0 and
    % the per-person rank correlation are gated as in the gamma class.
    % ADDITIONALLY gated here, which the gamma arm could not: parameter-level
    % recovery of the scale submodel (log residual-SD intercepts at 0.18 and
    % person scale SDs at 0.22, the TestSubjDiffDynrelRecovery precedent's
    % quantities and tolerances), plus a LOOSE sanity bound on the
    % estimand-irrelevant cell means (see the in-line calibration note).
    % Per-cell residual-SD CONTRASTS stay printed diagnostics.
    %
    % TRUTH IS LONGHAND. Unlike the gamma class, every planted-truth quantity
    % is computed from the planted numbers by explicit scalar arithmetic,
    % never through the committed kernels - the Gaussian closed form makes
    % that free, and it removes the partial-mirror hazard the gamma class
    % documents.
    %
    % PLANTED-TRUTH DOCTRINE (the family's): all-distinct values everywhere,
    % asymmetry WITHIN-PAIR, nonzero dimension slopes on BOTH submodels,
    % UNEQUAL per-cell trial counts - and, deliberately, NEGATIVE cell means:
    % they are the point of this arm (representable under Gaussian, rejected
    % under gamma by the preflight positivity check), and the run's family
    % stamp plus the negative data double as proof the Gaussian arm was
    % exercised.
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
                'cmdstan_artifacts', 'recov_dod_dynrel_gauss');
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

            %the planted data must genuinely exercise what gamma cannot
            %represent; with the family stamp below this doubles as proof the
            %Gaussian arm ran
            testCase.assertTrue(any(dataTbl.meas <= 0), ...
                'fixture must carry non-positive scores (the Gaussian point)');

            outDir = fullfile(localProjectRoot(), 'tests', 'integration', ...
                'cmdstan_artifacts', 'recov_dod_dynrel_gauss');
            if ~exist(outDir, 'dir'), mkdir(outDir); end
            rel = psyrat_computevarcomp( ...
                'data', dataTbl, ...
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
            testCase.verifyEqual(rel.family, 'gaussian');
            testCase.verifyEqual(rel.engine, 'cmdstan');
            testCase.verifyEqual( ...
                rel.out.dod_labels.inference_framework, ...
                'four_gaussian_margins_all_residual_covariances_fixed_zero_v1');

            % Convergence first: recovery from a non-converged fit says
            % nothing; skip rather than fail.
            [maxRhat, nFinite] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('\n[Gauss DoD dynrel recovery] max R-hat %.4f over %d params (converged = %d)\n', ...
                maxRhat, nFinite, rel.out.conv.converged);
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf( ...
                ['Fit did not converge (max R-hat %.4f); recovery is not ' ...
                'interpretable, so the run is skipped rather than failed.'], maxRhat));

            % ---- Population G_dod at z = 0: LONGHAND truth vs surface ----
            trueG = localTrueSurfaceG(T);
            summ = psyrat_dynrel_summary(rel, 'ngrid', 3, 'obs', T.ntrl);
            st = summ.strata(1);
            [~, iz0] = min(abs(st.z1));
            estG = st.G.pt(iz0);
            fprintf('[Gauss DoD dynrel recovery] true G(0) = %.4f | recovered G(0) = %.4f [%.4f, %.4f]\n', ...
                trueG, estG, st.G.ll(iz0), st.G.ul(iz0));
            testCase.verifyLessThanOrEqual(abs(estG - trueG), 0.05, ...
                ['The population G_dod at z = 0 must recover the planted ' ...
                'truth within 0.05.']);

            % ---- Per-person ordering: LONGHAND planted vs recovered G ----
            ssrel = st.ssrel_table;
            testCase.assertEqual(height(ssrel), T.nsub);
            trueSS = localTruePerPersonG(T, truthFx);
            [~, order] = ismember(ssrel.id, truthFx.ids);
            testCase.assertTrue(all(order > 0));
            rho = localSpearman(trueSS(order), ssrel.gen_pt);
            fprintf('[Gauss DoD dynrel recovery] per-person G rank correlation (true vs recovered) = %.3f\n', rho);
            testCase.verifyGreaterThanOrEqual(rho, 0.5, ...
                ['Participants with genuinely noisier planted data must ' ...
                'rank lower in the recovered conditional coefficients.']);

            % ---- Parameter-level gates (the Gaussian arm's extra check) ----
            b = mean(cell2mat(rel.out.b{1}), 1);
            bsig = mean(cell2mat(rel.out.b_sigma{1}), 1);
            sdid = mean(cell2mat(rel.out.sd_id{1}), 1);
            fprintf('[Gauss DoD dynrel recovery] cell means true = [%s] recov = [%s]\n', ...
                sprintf('%.2f ', T.b), sprintf('%.2f ', b));
            fprintf('[Gauss DoD dynrel recovery] log res SD true = [%s] recov = [%s]\n', ...
                sprintf('%.3f ', T.b_sig), sprintf('%.3f ', bsig));
            fprintf('[Gauss DoD dynrel recovery] person scale SDs true = [%s] recov = [%s]\n', ...
                sprintf('%.3f ', T.sd_s), sprintf('%.3f ', sdid(5:8)));
            %the cell means are ESTIMAND-IRRELEVANT (the section-26
            %mean-invariance theorem: b enters no reliability quantity) and
            %their posterior SD at this design is ~sd_p/sqrt(nsub) ~= 2.2/sqrt(60) ~= 0.28,
            %so a tight gate would assert on sampling noise. This is a loose
            %WIRING sanity bound (~4 posterior SDs), not a recovery claim;
            %the recovery claims live on the scale submodel, which the
            %estimand actually consumes.
            testCase.verifyLessThanOrEqual(max(abs(b - T.b)), 1.0, ...
                'Identity-scale cell means must land within the sanity bound.');
            testCase.verifyLessThanOrEqual(max(abs(bsig - T.b_sig)), 0.18, ...
                'Log residual-SD intercepts must recover within 0.18.');
            testCase.verifyLessThanOrEqual(max(abs(sdid(5:8) - T.sd_s)), 0.22, ...
                'Person log-residual-SD SDs must recover within 0.22.');

            % ---- DIAGNOSTIC ONLY: per-cell residual-SD contrasts ----
            fprintf('[Gauss DoD dynrel recovery] DIAGNOSTIC (not asserted) SD contrast (1-2), true = %.3f, recov = %.3f\n', ...
                exp(T.b_sig(1)) - exp(T.b_sig(2)), exp(bsig(1)) - exp(bsig(2)));
        end

    end
end

% ---------------------------------------------------------------------------

function root = localProjectRoot()
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
end

function T = localTruth()
%All-distinct planted truth on GAUSSIAN scales; asymmetry WITHIN-PAIR; two
%NEGATIVE cell means deliberately (ERN-like; the point of this arm).
T = struct();
T.events = {'E1','E2','E3','E4'};
T.cvec = [1 -1 -1 1];
T.b = [-4.8 -1.9 -4.1 -3.4];                       %per-cell means at z = 0
T.b_dim = [0.6 -0.25 0.4 -0.1];                    %per-cell mean slope on dim1
T.b_sig = [0.61 1.01 0.64 0.66];                   %per-cell log residual SD at z = 0
T.b_sig_dim = [-0.08 0.12 -0.05 0.07];             %per-cell scale slope on dim1
T.sd_p = [2.0 1.7 2.2 1.9];                        %person location SDs (natural scale)
T.sd_s = [0.30 0.34 0.28 0.32];                    %person log-sigma SDs
T.sd_i = [0.8 0.7 0.9 0.75];                       %trial SDs
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
%Simulate the planted GAUSSIAN model: per-subject 8-D joint person effects,
%per-ordinal-trial 4-D correlated trial effects (shared crossed levels),
%Normal observations, and a standardized dimension column. Base MATLAB only.
rng(20260820, 'twister');
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
        mu = T.b(q) + T.b_dim(q)*dim1(s) + u_p(s, q);
        sigma = exp(T.b_sig(q) + T.b_sig_dim(q)*dim1(s) + u_p(s, 4+q));
        for t = 1:T.ntrl(q)
            y = mu + u_i(t, q) + sigma * randn;
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
%LONGHAND typical-person G_dod at z = 0 from the planted values: identity
%link, so the universe is the person-block quadratic form and the relative
%error is the four-term residual sum. No committed kernel is called.
w = T.cvec(:);
Sp = diag(T.sd_p) * T.R8(1:4,1:4) * diag(T.sd_p);
u = w' * Sp * w;
relerr = sum((w'.^2) .* exp(2 .* T.b_sig) ./ T.ntrl);
G = u / (u + relerr);
end

function g = localTruePerPersonG(T, fx)
%LONGHAND planted per-person conditional G at each subject's own z with their
%own residual SD (their realized scale effect). No committed kernel.
w = T.cvec(:);
Sp = diag(T.sd_p) * T.R8(1:4,1:4) * diag(T.sd_p);
u = w' * Sp * w;
g = zeros(T.nsub, 1);
for s = 1:T.nsub
    z = fx.dim1(s);
    lsig = T.b_sig + T.b_sig_dim .* z + fx.u_p(s, 5:8);
    relerr = sum((w'.^2) .* exp(2 .* lsig) ./ T.ntrl);
    g(s) = u / (u + relerr);
end
end

function rho = localSpearman(a, b)
%Spearman rank correlation, base MATLAB (no Statistics Toolbox corr/tiedrank).
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
