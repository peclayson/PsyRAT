classdef TestGammaCrossCov < matlab.unittest.TestCase
    %TESTGAMMACROSSCOV Unit tests for the Gamma-family cross-cell observed-score
    % covariance conversion (psyrat_gamma_crosscov), used to build the person and
    % trial 2x2 (co)variance arrays for one-facet difference-score reliability.
    %
    % On the log expected-score scale the two cells' person effects are bivariate
    % normal with covariance cor_p*s_p_f*s_p_s (and analogously for trials). The
    % observed-score covariance of the two cells' marginal means is the covariance
    % of two jointly lognormal variables. These tests derive that covariance
    % independently (the general jointly-lognormal covariance identity, and a
    % fixed-seed correlated-normal Monte Carlo) so a green test validates the
    % formula rather than merely pinning it.

    methods (Test)

        function testPersonCovAgainstJointLognormalIdentity(testCase)
            % For X = exp(a + W1), Y = exp(b + W2) with (W1,W2) bivariate normal
            % (variances v1,v2; covariance c12):
            %   Cov(X,Y) = exp(a + b + (v1+v2)/2) * (exp(c12) - 1).
            % The person-marginal mean of cell k is exp(alpha_k + u_pk + si2_k/2),
            % so a = alpha_f + si2_f/2, b = alpha_s + si2_s/2, W = person effect,
            % v = sp2, c12 = cor_p*s_p_f*s_p_s. Written independently of production.
            alpha_f = log(6.3); s_p_f = 0.42; s_i_f = 0.12;
            alpha_s = log(6.0); s_p_s = 0.38; s_i_s = 0.09;
            cor_p = 0.65; cor_i = 0.30;

            cov = psyrat_gamma_crosscov(alpha_f, s_p_f, s_i_f, ...
                alpha_s, s_p_s, s_i_s, cor_p, cor_i);

            sp2f = s_p_f^2; si2f = s_i_f^2; sp2s = s_p_s^2; si2s = s_i_s^2;

            % Person cross covariance via the jointly-lognormal identity.
            a = alpha_f + si2f/2; b = alpha_s + si2s/2;
            c12_p = cor_p * s_p_f * s_p_s;
            expected_p = exp(a + b + (sp2f + sp2s)/2) * (exp(c12_p) - 1);

            % Trial cross covariance: symmetric, marginalizing person instead.
            a_t = alpha_f + sp2f/2; b_t = alpha_s + sp2s/2;
            c12_i = cor_i * s_i_f * s_i_s;
            expected_i = exp(a_t + b_t + (si2f + si2s)/2) * (exp(c12_i) - 1);

            testCase.verifyEqual(cov.cov_p_obs, expected_p, 'RelTol', 1e-12);
            testCase.verifyEqual(cov.cov_i_obs, expected_i, 'RelTol', 1e-12);
            % Non-concurrent difference: residual covariance is defined as 0.
            testCase.verifyEqual(cov.cov_e_obs, 0, 'AbsTol', 1e-12);
        end

        function testMonteCarloCorrelatedPersonEffects(testCase)
            % Fixed-seed correlated-normal Monte Carlo (base MATLAB randn +
            % Cholesky), independent validation of the person cross covariance.
            alpha_f = log(6.3); s_p_f = 0.42; s_i_f = 0.12;
            alpha_s = log(6.0); s_p_s = 0.38; s_i_s = 0.09;
            cor_p = 0.65; cor_i = 0.30;

            cov = psyrat_gamma_crosscov(alpha_f, s_p_f, s_i_f, ...
                alpha_s, s_p_s, s_i_s, cor_p, cor_i);

            rng(7);
            n = 4e6;

            % Correlated person effects with covariance cor_p*s_p_f*s_p_s.
            Sp = [s_p_f^2,               cor_p*s_p_f*s_p_s; ...
                  cor_p*s_p_f*s_p_s,     s_p_s^2];
            Lp = chol(Sp, 'lower');
            Up = (Lp * randn(2,n))';
            Af = exp(alpha_f + Up(:,1) + (s_i_f^2)/2);   % person-marginal mean, cell f
            As = exp(alpha_s + Up(:,2) + (s_i_s^2)/2);   % person-marginal mean, cell s
            cp = cov_(Af, As);
            testCase.verifyEqual(cp, cov.cov_p_obs, 'RelTol', 1.5e-2);

            % Correlated trial effects with covariance cor_i*s_i_f*s_i_s.
            Si = [s_i_f^2,               cor_i*s_i_f*s_i_s; ...
                  cor_i*s_i_f*s_i_s,     s_i_s^2];
            Li = chol(Si, 'lower');
            Ut = (Li * randn(2,n))';
            Bf = exp(alpha_f + Ut(:,1) + (s_p_f^2)/2);   % trial-marginal mean, cell f
            Bs = exp(alpha_s + Ut(:,2) + (s_p_s^2)/2);
            ci = cov_(Bf, Bs);
            testCase.verifyEqual(ci, cov.cov_i_obs, 'RelTol', 1.5e-2);
        end

    end
end

function c = cov_(x, y)
% Population covariance (normalized by N) to match analytic moments.
x = x(:); y = y(:);
c = mean((x - mean(x)) .* (y - mean(y)));
end
