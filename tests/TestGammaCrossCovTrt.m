classdef TestGammaCrossCovTrt < matlab.unittest.TestCase
    %TESTGAMMACROSSCOVTRT Unit tests for the crossed Persons x Occasions x Trials
    % cross-cell observed-score covariance conversion (psyrat_gamma_crosscov_trt),
    % used by the two-facet (test-retest) two-event difference-score reliability of
    % the scaled-chi-square / Gamma family. It is the two-facet analog of
    % psyrat_gamma_crosscov: the two cells' per-factor effects are correlated on
    % the log expected-score scale, and their observed-score cross-covariances come
    % from jointly-lognormal moments.
    %
    % The checks are INDEPENDENT of the production expressions:
    %   - zero cross-cell correlations give zero cross-covariances;
    %   - with only person + trial variation the helper reproduces the validated
    %     one-facet psyrat_gamma_crosscov for the person and trial cross-covs, and
    %     the induced person x trial cross-cov equals its closed form;
    %   - a fixed-seed generative Monte Carlo estimates each cross-covariance from
    %     the covariance-of-cells definition (two cells that share a correlated
    %     factor set have observed cross-covariance equal to the sum of the cross-
    %     covariance components indexed within that set).

    methods (Test)

        function testZeroCorrelationsGiveZeroCovariances(testCase)
            % No cross-cell correlation on any factor => every observed cross-cell
            % covariance is exactly zero, regardless of the marginal SDs.
            aX = log(5); aY = log(4);
            sX = [0.30 0.20 0.15 0.12 0.10 0.08];
            sY = [0.25 0.18 0.14 0.11 0.09 0.07];
            cov = psyrat_gamma_crosscov_trt(aX, sX(1), sX(2), sX(3), sX(4), sX(5), sX(6), ...
                aY, sY(1), sY(2), sY(3), sY(4), sY(5), sY(6), 0, 0, 0, 0, 0, 0);
            for f = {'cov_p','cov_o','cov_t','cov_po','cov_pt','cov_ot', ...
                    'cov_res','cov_mu_total'}
                testCase.verifyEqual(cov.(f{1}), 0, 'AbsTol', 1e-12);
            end
        end

        function testReducesToOneFacetCrossCov(testCase)
            % With occasion and all interaction SDs zero for both cells, the person
            % and trial cross-covariances must match the validated one-facet helper
            % (psyrat_gamma_crosscov), the occasion/interaction-with-occasion cross-
            % covs are zero, and the induced person x trial cross-cov equals its
            % closed form. Vectorized over draws.
            aX = [log(5.5); log(6.0)];  aY = [log(4.8); log(5.2)];
            sPx = [0.40; 0.30];  sTx = [0.12; 0.20];
            sPy = [0.35; 0.28];  sTy = [0.10; 0.18];
            corP = [0.6; 0.4];   corT = [0.3; 0.5];

            onef = psyrat_gamma_crosscov(aX, sPx, sTx, aY, sPy, sTy, corP, corT);
            cov  = psyrat_gamma_crosscov_trt(aX, sPx, 0, sTx, 0, 0, 0, ...
                aY, sPy, 0, sTy, 0, 0, 0, corP, 0, corT, 0, 0, 0);

            testCase.verifyEqual(cov.cov_p, onef.cov_p_obs, 'RelTol', 1e-12);
            testCase.verifyEqual(cov.cov_t, onef.cov_i_obs, 'RelTol', 1e-12);
            testCase.verifyEqual(cov.cov_o,  zeros(2,1), 'AbsTol', 1e-12);
            testCase.verifyEqual(cov.cov_po, zeros(2,1), 'AbsTol', 1e-12);
            testCase.verifyEqual(cov.cov_ot, zeros(2,1), 'AbsTol', 1e-12);

            % Induced person x trial cross-cov (no log-scale interaction): the
            % jointly-lognormal cross second moment minus the two main-effect cross
            % covariances, written independently of the production expression.
            EYx = exp(aX + 0.5.*(sPx.^2 + sTx.^2));
            EYy = exp(aY + 0.5.*(sPy.^2 + sTy.^2));
            covPlog = corP .* sPx .* sPy;
            covTlog = corT .* sTx .* sTy;
            expected_pt = EYx.*EYy .* (exp(covPlog + covTlog) - exp(covPlog) ...
                - exp(covTlog) + 1);
            testCase.verifyEqual(cov.cov_pt, expected_pt, 'RelTol', 1e-12);
        end

        function testMonteCarloCrossCovariances(testCase)
            % Fixed-seed generative Monte Carlo: two cells whose effects for a
            % chosen factor set are cross-correlated (and independent otherwise)
            % have observed cross-covariance equal to the sum of the cross-cov
            % components indexed within that set. Base-MATLAB randn only.
            aX = log(5.0); aY = log(4.4);
            sPx=0.34; sOx=0.22; sTx=0.16; sPOx=0.12; sPTx=0.10; sOTx=0.08;
            sPy=0.30; sOy=0.20; sTy=0.14; sPOy=0.11; sPTy=0.09; sOTy=0.07;
            cP=0.6; cO=0.5; cT=0.4; cPO=0.45; cPT=0.35; cOT=0.30;

            cov = psyrat_gamma_crosscov_trt(aX, sPx, sOx, sTx, sPOx, sPTx, sOTx, ...
                aY, sPy, sOy, sTy, sPOy, sPTy, sOTy, cP, cO, cT, cPO, cPT, cOT);

            rng(2027);
            n = 4e6;
            pcov = @(x,y) mean(x.*y) - mean(x).*mean(y);

            % Correlated bivariate pair for one factor (both cells share it), or an
            % independent pair (each cell draws its own). z-store lets each factor's
            % pair be reused across the "share set" estimators.
            mk = @(sx, sy, c) local_pair(n, sx, sy, c);

            % Precompute both a correlated and an independent pair per factor.
            [pX,pY] = mk(sPx, sPy, cP);   [pXi,pYi] = mk(sPx, sPy, 0);
            [oX,oY] = mk(sOx, sOy, cO);   [oXi,oYi] = mk(sOx, sOy, 0);
            [tX,tY] = mk(sTx, sTy, cT);   [tXi,tYi] = mk(sTx, sTy, 0);
            [poX,poY] = mk(sPOx, sPOy, cPO); [poXi,poYi] = mk(sPOx, sPOy, 0);
            [ptX,ptY] = mk(sPTx, sPTy, cPT); [ptXi,ptYi] = mk(sPTx, sPTy, 0);
            [otX,otY] = mk(sOTx, sOTy, cOT); [otXi,otYi] = mk(sOTx, sOTy, 0);

            % cellX is always the same set of x-effects; cellY shares (correlated)
            % the factors in the set and draws independent y-effects otherwise.
            xsum = aX + pX + oX + tX + poX + ptX + otX;

            covShare = @(yShareSum) pcov(exp(xsum), exp(yShareSum));

            % Share exactly one main-effect factor.
            cov_p_hat = covShare(aY + pY  + oYi + tYi + poYi + ptYi + otYi);
            cov_o_hat = covShare(aY + pYi + oY  + tYi + poYi + ptYi + otYi);
            cov_t_hat = covShare(aY + pYi + oYi + tY  + poYi + ptYi + otYi);
            testCase.verifyEqual(cov_p_hat, cov.cov_p, 'RelTol', 3e-2);
            testCase.verifyEqual(cov_o_hat, cov.cov_o, 'RelTol', 3e-2);
            testCase.verifyEqual(cov_t_hat, cov.cov_t, 'RelTol', 3e-2);

            % Share both parents and the interaction => sum of the three cross-covs.
            cov_po_sum = covShare(aY + pY + oY + tYi + poY + ptYi + otYi);
            testCase.verifyEqual(cov_po_sum - cov_p_hat - cov_o_hat, cov.cov_po, ...
                'RelTol', 5e-2);
        end

    end
end

function [ex, ey] = local_pair(n, sx, sy, c)
% Correlated (c ~= 0) or independent (c == 0) bivariate normal effect pair for
% one factor, one column per cell, using base-MATLAB randn only.
z1 = randn(n, 1);
z2 = randn(n, 1);
ex = sx * z1;
ey = sy * (c * z1 + sqrt(1 - c^2) * z2);
end
