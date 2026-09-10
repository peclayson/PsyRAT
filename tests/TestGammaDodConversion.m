classdef TestGammaDodConversion < matlab.unittest.TestCase
    %TESTGAMMADODCONVERSION Unit tests for the four-event difference-of-differences
    % observed-score assembly of the scaled-chi-square / Gamma family (analysis 9).
    %
    % The gamma DoD path reuses the validated one-facet conversion helpers
    % (psyrat_gamma_varcomps diagonals, psyrat_gamma_crosscov cross-cell covs),
    % assembles the four-cell person / trial (co)variance matrices and the per-cell
    % residual b_sigma, and feeds psyrat_dodiffrel (which applies the c=[1 -1 -1 1]
    % contrast). These checks are INDEPENDENT of the production extractor: they
    % transcribe the supplied reference's extract_dofd_results reliability math
    % (model_fits_dofd_chisq_stancode_standata.R) and confirm the toolbox pipeline
    % reproduces it under the reference's non-concurrent defaults (person + trial
    % cross-covariances kept; residual cross-covariance 0, structural in the
    % per-cell er_var vector psyrat_dodiffrel consumes).
    %
    % SCOPE: these tests isolate the estimand-INDEPENDENT contrast / qform / DoD
    % assembly. Both sides (localAssemble and localRefDod) call the conversion
    % with the POPULATION nu (no dispersion args), matching the reference's
    % population-nu observation term, so they remain a faithful reference-parity
    % check of the assembly. The production extractor (psyrat_gamma_extract_dodiff)
    % now uses the MARGINAL observation term (estimand #2), which multiplies the
    % per-cell e_cond_var by exp(0.5*s_v^2 - 2*cov_pv); that estimand change is a
    % pure per-cell rescale of the residual and is verified separately by
    % TestGammaDispersionMarginalization (offline Monte-Carlo, incl. per-cell
    % independence) and the live DoD recovery test, so it is deliberately NOT
    % re-exercised here.

    methods (Test)

        function testAssemblyMatchesReferenceDodMath(testCase)
            % Random pseudo-posterior draws of the per-cell log-scale parameters;
            % assemble via the toolbox helpers -> psyrat_dodiffrel and compare G, D,
            % both ICCs, and the universe-score variance to an INDEPENDENT
            % transcription of the reference extract_dofd_results math. They must
            % agree to machine precision (the toolbox helpers ARE the reference
            % onefacet functions; psyrat_dodiffrel IS the reference DoD formula).
            rng(24601, 'twister');
            nd = 3000; C = 4;
            pairs = nchoosek(1:C, 2);
            obs = [18 22 15 20];
            cvec = [1 -1 -1 1];

            alpha = 2.5 + 0.30 .* randn(nd, C);
            sd_p  = abs(0.35 + 0.10 .* randn(nd, C));
            sd_i  = abs(0.18 + 0.06 .* randn(nd, C));
            nu    = exp(log(8) + 0.30 .* randn(nd, C));
            cor_p = 0.30 + 0.50 .* rand(nd, size(pairs,1));
            cor_i = 0.00 + 0.40 .* rand(nd, size(pairs,1));

            [idv, tlv, bsig] = localAssemble(alpha, sd_p, sd_i, nu, cor_p, cor_i, pairs);
            gen = psyrat_dodiffrel('bp', idv, 'bt', tlv, 'er_var', bsig, ...
                'obs', obs, 'est', 'gen', 'ci', .95, 'cvec', cvec);
            dep = psyrat_dodiffrel('bp', idv, 'bt', tlv, 'er_var', bsig, ...
                'obs', obs, 'est', 'dep', 'ci', .95, 'cvec', cvec);

            ref = localRefDod(alpha, sd_p, sd_i, nu, cor_p, cor_i, pairs, obs, cvec);

            testCase.verifyEqual(gen.pt,      ref.G_pt,     'AbsTol', 1e-9);
            testCase.verifyEqual(dep.pt,      ref.D_pt,     'AbsTol', 1e-9);
            testCase.verifyEqual(gen.ll,      ref.G_ll,     'AbsTol', 1e-9);
            testCase.verifyEqual(gen.ul,      ref.G_ul,     'AbsTol', 1e-9);
            testCase.verifyEqual(dep.ll,      ref.D_ll,     'AbsTol', 1e-9);
            testCase.verifyEqual(dep.ul,      ref.D_ul,     'AbsTol', 1e-9);
            testCase.verifyEqual(gen.icc_pt,  ref.ICCrel_pt,'AbsTol', 1e-9);
            testCase.verifyEqual(dep.icc_pt,  ref.ICCabs_pt,'AbsTol', 1e-9);
            testCase.verifyEqual(gen.uni_var_pt, ref.uni_pt,'AbsTol', 1e-9);
        end

        function testUniverseVarianceIsWeightedDiagonalSumWhenUncorrelated(testCase)
            % With zero cross-cell person correlations, the assembled person matrix
            % is diagonal, so the DoD universe-score variance reduces to the
            % weighted sum of per-cell person variances, sum_i c_i^2 * sigma_p2_i.
            C = 4; pairs = nchoosek(1:C,2);
            alpha = log([5.0 4.5 5.2 4.8]);
            sd_p  = [0.40 0.38 0.42 0.39];
            sd_i  = [0.20 0.18 0.22 0.19];
            nu    = [8 8 8 8];
            cor_p = zeros(1, size(pairs,1));   % no person cross-cell correlation
            cor_i = zeros(1, size(pairs,1));
            cvec  = [1 -1 -1 1];

            [idv, tlv, bsig] = localAssemble(alpha, sd_p, sd_i, nu, cor_p, cor_i, pairs);
            out = psyrat_dodiffrel('bp', idv, 'bt', tlv, 'er_var', bsig, ...
                'obs', [20 20 20 20], 'est', 'gen', 'ci', .95, 'cvec', cvec);

            sig_p2 = zeros(1,C);
            for c = 1:C
                vc = psyrat_gamma_varcomps(alpha(c), sd_p(c), sd_i(c), nu(c));
                sig_p2(c) = vc.sigma_p2;
            end
            expected = sum((cvec.^2) .* sig_p2);
            testCase.verifyEqual(out.uni_var_pt, expected, 'RelTol', 1e-12);
        end

        function testResidualIsPerCellNonConcurrent(testCase)
            % The residual enters psyrat_dodiffrel as a per-cell vector (er_var),
            % so the unscaled error contrast variance is sum_i c_i^2 * sigma_pi_e2_i
            % (no residual cross-covariance). This pins the b_sigma = 0.5*log(res)
            % encoding round-trip through psyrat_dodiffrel's exp(er_var).^2.
            C = 4; pairs = nchoosek(1:C,2);
            alpha = log([5.0 4.5 5.2 4.8]);
            sd_p  = [0.40 0.38 0.42 0.39];
            sd_i  = [0.20 0.18 0.22 0.19];
            nu    = [7 9 8 10];
            cor_p = 0.5 * ones(1, size(pairs,1));
            cor_i = 0.2 * ones(1, size(pairs,1));
            cvec  = [1 -1 -1 1];

            [idv, tlv, bsig] = localAssemble(alpha, sd_p, sd_i, nu, cor_p, cor_i, pairs);
            out = psyrat_dodiffrel('bp', idv, 'bt', tlv, 'er_var', bsig, ...
                'obs', [20 20 20 20], 'est', 'gen', 'ci', .95, 'cvec', cvec);

            sig_pie2 = zeros(1,C);
            for c = 1:C
                vc = psyrat_gamma_varcomps(alpha(c), sd_p(c), sd_i(c), nu(c));
                sig_pie2(c) = vc.sigma_pi_e2;
            end
            expected = sum((cvec.^2) .* sig_pie2);
            testCase.verifyEqual(out.er_contrast_var_pt, expected, 'RelTol', 1e-12);
        end

    end
end

% =========================================================================
function [id_varcov, trl_varcov, b_sigma] = localAssemble(alpha, sd_p, sd_i, nu, cor_p, cor_i, pairs)
% Replicate the production psyrat_gamma_extract_dodiff assembly using the public
% conversion helpers: per-cell diagonals from psyrat_gamma_varcomps, cross-cell
% off-diagonals from psyrat_gamma_crosscov, per-cell residual encoded as
% b_sigma = 0.5*log(sigma_pi_e2).
nd = size(alpha,1); C = size(alpha,2);
id_varcov  = zeros(nd,C,C);
trl_varcov = zeros(nd,C,C);
b_sigma    = zeros(nd,C);
for c = 1:C
    vc = psyrat_gamma_varcomps(alpha(:,c), sd_p(:,c), sd_i(:,c), nu(:,c));
    id_varcov(:,c,c)  = vc.sigma_p2;
    trl_varcov(:,c,c) = vc.sigma_i2;
    b_sigma(:,c)      = 0.5 .* log(vc.sigma_pi_e2);
end
for pp = 1:size(pairs,1)
    j = pairs(pp,1); k = pairs(pp,2);
    cc = psyrat_gamma_crosscov(alpha(:,j), sd_p(:,j), sd_i(:,j), ...
                               alpha(:,k), sd_p(:,k), sd_i(:,k), ...
                               cor_p(:,pp), cor_i(:,pp));
    id_varcov(:,j,k)  = cc.cov_p_obs; id_varcov(:,k,j)  = cc.cov_p_obs;
    trl_varcov(:,j,k) = cc.cov_i_obs; trl_varcov(:,k,j) = cc.cov_i_obs;
end
end

% =========================================================================
function ref = localRefDod(alpha, sd_p, sd_i, nu, cor_p, cor_i, pairs, obs, cvec)
% Independent transcription of the reference extract_dofd_results reliability
% math (onefacet_components_from_params + onefacet_cov_components_from_params +
% qform / qform_div_n), non-concurrent defaults.
nd = size(alpha,1); C = size(alpha,2); w = cvec(:); n = obs(:);
np = size(pairs,1);

sig_p2 = zeros(nd,C); sig_i2 = zeros(nd,C); sig_pie2 = zeros(nd,C);
for c = 1:C
    vp = sd_p(:,c).^2; vi = sd_i(:,c).^2; vt = vp + vi;
    base = exp(2.*alpha(:,c) + vt);
    sig_p2(:,c)   = base .* (exp(vp) - 1);
    sig_i2(:,c)   = base .* (exp(vi) - 1);
    var_mu        = base .* (exp(vt) - 1);
    e_cond        = (2./nu(:,c)) .* exp(2.*alpha(:,c) + 2.*vt);
    sig_pie2(:,c) = max(var_mu + e_cond - sig_p2(:,c) - sig_i2(:,c), 0);
end

cov_p = zeros(nd,np); cov_i = zeros(nd,np);
for pp = 1:np
    j = pairs(pp,1); k = pairs(pp,2);
    ey = exp(alpha(:,j) + 0.5.*(sd_p(:,j).^2 + sd_i(:,j).^2)) ...
       .* exp(alpha(:,k) + 0.5.*(sd_p(:,k).^2 + sd_i(:,k).^2));
    cov_p(:,pp) = ey .* (exp(cor_p(:,pp).*sd_p(:,j).*sd_p(:,k)) - 1);
    cov_i(:,pp) = ey .* (exp(cor_i(:,pp).*sd_i(:,j).*sd_i(:,k)) - 1);
end

G = zeros(nd,1); D = zeros(nd,1); ICCr = zeros(nd,1); ICCa = zeros(nd,1); U = zeros(nd,1);
for d = 1:nd
    Vp = diag(sig_p2(d,:)); Vi = diag(sig_i2(d,:)); Vpie = diag(sig_pie2(d,:));
    for pp = 1:np
        j = pairs(pp,1); k = pairs(pp,2);
        Vp(j,k)=cov_p(d,pp); Vp(k,j)=cov_p(d,pp);
        Vi(j,k)=cov_i(d,pp); Vi(k,j)=cov_i(d,pp);
    end
    u        = w'*Vp*w;
    trial_ab = localQformDivN(Vi, w, n);
    rel_err  = localQformDivN(Vpie, w, n);
    U(d)   = u;
    G(d)   = u / (u + rel_err);
    D(d)   = u / (u + rel_err + trial_ab);
    ICCr(d)= u / (u + w'*Vpie*w);
    ICCa(d)= u / (u + w'*Vpie*w + w'*Vi*w);
end

qg = quantile(G,[.025 .975]); qd = quantile(D,[.025 .975]);
ref = struct('G_pt',mean(G),'G_ll',qg(1),'G_ul',qg(2), ...
             'D_pt',mean(D),'D_ll',qd(1),'D_ul',qd(2), ...
             'ICCrel_pt',mean(ICCr),'ICCabs_pt',mean(ICCa),'uni_pt',mean(U));
end

% =========================================================================
function v = localQformDivN(V, w, n)
% Reference qform_div_n: diagonal / n_j, off-diagonal / harmonic-mean(n_j,n_k).
C = numel(w); v = 0;
for j = 1:C, v = v + (w(j)^2) * V(j,j) / n(j); end
for j = 1:C-1
    for k = j+1:C
        nharm = 2 / ((1/n(j)) + (1/n(k)));
        v = v + 2*w(j)*w(k)*V(j,k)/nharm;
    end
end
end
