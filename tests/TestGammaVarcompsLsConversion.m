classdef TestGammaVarcompsLsConversion < matlab.unittest.TestCase
    %TESTGAMMAVARCOMPSLSCONVERSION Unit tests for the LOCATION-SCALE
    % scaled-chi-square / Gamma log-scale -> observed-score-scale variance-
    % component conversion (psyrat_gamma_varcomps_ls, gammascale = 2).
    %
    % The location-scale gamma parameterization models the residual SD directly,
    %   log mu    = alpha(z) + u_p + u_t,     u_p ~ N(0, s_p^2), u_t ~ N(0, s_i^2)
    %   log sigma = log_sigma(z) + w_p,       w_p ~ N(0, s_sd^2)
    %   nu        = 2*(mu/sigma)^2,
    % so Var(Y | mu, sigma) = sigma^2 EXACTLY and the residual is decoupled from
    % the mean surface. That decoupling is the whole point of the
    % parameterization and is what gives a clean residual-SD dimension slope, so
    % it is pinned directly below rather than assumed.
    %
    % Contrast with the log-nu sibling (psyrat_gamma_varcomps), where
    % Var(Y | mu, nu) = 2*mu^2/nu is MEAN-COUPLED and its marginal expected
    % conditional variance therefore carries exp(0.5*s_v^2 - 2*cov_pv). The
    % location-scale form has no such term: the person mean <-> scale
    % correlation does not enter e_cond_var at all. Both facts are tested.
    %
    % As in TestGammaVarcompsConversion, the target quantities are derived
    % INDEPENDENTLY of the production expressions (closed forms, an invariance
    % argument, and fixed-seed Monte Carlo) so a green test means the formula is
    % right, not merely unchanged. Base MATLAB only: no Statistics Toolbox.

    methods (Test)

        function testDegenerateZeroRandomEffects(testCase)
            % With no person/trial variation and no person scale heterogeneity
            % the mean surface is constant at exp(alpha) and the only observed
            % variance is the conditional dispersion, which for this
            % parameterization is exactly sigma^2.
            alpha     = log(3);
            sigma     = 1.4;
            log_sigma = log(sigma);

            vc = psyrat_gamma_varcomps_ls(alpha, 0, 0, log_sigma);

            testCase.verifyEqual(vc.sigma_p2,    0,          'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_i2,    0,          'AbsTol', 1e-12);
            testCase.verifyEqual(vc.var_mu,      0,          'AbsTol', 1e-12);
            testCase.verifyEqual(vc.e_cond_var,  sigma^2,    'AbsTol', 1e-12);
            testCase.verifyEqual(vc.total_var,   sigma^2,    'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_pi_e2, sigma^2,    'AbsTol', 1e-12);
            testCase.verifyEqual(vc.mu_bar,      exp(alpha), 'AbsTol', 1e-12);
        end

        function testResidualIsDecoupledFromTheMeanSurface(testCase)
            % THE defining property of gammascale = 2, and the reason the ruling
            % went this way: e_cond_var depends only on the scale submodel. Move
            % alpha, s_p and s_i over a wide range with the scale submodel held
            % fixed and the residual must not budge. Under the log-nu sibling
            % the same sweep changes e_cond_var by orders of magnitude, which the
            % second half of this test asserts so the contrast is pinned, not
            % merely described in a comment.
            log_sigma = log(1.7);
            alpha = [log(2.0); log(5.5); log(40.0)];
            s_p   = [0.10;     0.45;     0.80];
            s_i   = [0.05;     0.20;     0.60];

            vc = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, log_sigma);

            testCase.verifyEqual(vc.e_cond_var, repmat(1.7^2, 3, 1), 'RelTol', 1e-12);

            % The log-nu sibling, at a nu chosen to match the residual for the
            % FIRST row, does not hold it fixed across the other two.
            ssum1   = s_p(1)^2 + s_i(1)^2;
            nu_match = 2 * exp(2*alpha(1) + 2*ssum1) / 1.7^2;
            vc_nu = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_match);
            testCase.verifyEqual(vc_nu.e_cond_var(1), 1.7^2, 'RelTol', 1e-12);
            testCase.verifyGreaterThan(vc_nu.e_cond_var(3) / vc_nu.e_cond_var(1), 100);
        end

        function testPersonScaleHeterogeneityMarginalization(testCase)
            % e_cond_var is the MARGINAL expected conditional variance,
            % E[sigma_p^2] over the person population, so a lognormal second
            % moment: E[exp(2*log_sigma + 2*w_p)] = exp(2*log_sigma + 2*s_sd^2).
            % Checked against fixed-seed Monte Carlo over w_p using base-MATLAB
            % randn, which is independent of the analytic expression.
            alpha = log(5.5); s_p = 0.45; s_i = 0.20;
            log_sigma = log(1.7); s_sd = 0.40;

            vc = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, log_sigma, s_sd);

            rng(4041);
            n  = 4e6;
            wp = s_sd * randn(n,1);
            sigma_p = exp(log_sigma + wp);
            testCase.verifyEqual(mean(sigma_p.^2), vc.e_cond_var, 'RelTol', 1e-2);

            % Total-variance additive identity must hold exactly whenever the
            % residual is not clipped.
            testCase.verifyEqual(vc.sigma_p2 + vc.sigma_i2 + vc.sigma_pi_e2, ...
                vc.total_var, 'RelTol', 1e-12);
        end

        function testPersonMeanScaleCorrelationDoesNotEnterResidual(testCase)
            % Because the residual never touches mu, the person mean <-> scale
            % correlation cannot affect e_cond_var. The function therefore takes
            % no cov_pv argument at all; this test pins that the omission is
            % correct rather than an oversight, by confirming Monte Carlo
            % agreement when the two person effects are strongly correlated.
            alpha = log(5.5); s_p = 0.45; s_i = 0.20;
            log_sigma = log(1.7); s_sd = 0.40; rho = 0.9;

            vc = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, log_sigma, s_sd);

            rng(4042);
            n = 4e6;
            z = randn(n,2);
            up = s_p * z(:,1);
            wp = s_sd * (rho*z(:,1) + sqrt(1-rho^2)*z(:,2));
            sigma_p = exp(log_sigma + wp);

            % Correlated with the person mean effect, yet the marginal second
            % moment is unchanged.
            testCase.verifyGreaterThan(abs(corr_base(up, wp)), 0.85);
            testCase.verifyEqual(mean(sigma_p.^2), vc.e_cond_var, 'RelTol', 1e-2);
        end

        function testMeanSurfaceIdenticalToLogNuSibling(testCase)
            % The mean-surface components depend only on the log-mean submodel,
            % which the two parameterizations share, so they must agree exactly.
            % This pins the plan's claim that only e_cond_var changes.
            alpha = [log(5.5); log(6.0); log(4.2)];
            s_p   = [0.40;     0.25;     0.55];
            s_i   = [0.10;     0.18;     0.05];

            ls = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, log(1.7));
            nu = psyrat_gamma_varcomps(alpha, s_p, s_i, 18);

            testCase.verifyEqual(ls.mu_bar,   nu.mu_bar,   'RelTol', 1e-12);
            testCase.verifyEqual(ls.sigma_p2, nu.sigma_p2, 'RelTol', 1e-12);
            testCase.verifyEqual(ls.sigma_i2, nu.sigma_i2, 'RelTol', 1e-12);
            testCase.verifyEqual(ls.var_mu,   nu.var_mu,   'RelTol', 1e-12);
        end

        function testMatchesLogNuSiblingAtEquivalentDispersion(testCase)
            % The two parameterizations describe the same likelihood family, so
            % at the nu that reproduces a given residual they must return the
            % same components. This defines the correspondence between
            % gammascale 1 and 2 at a single reference point and is the identity
            % any cross-parameterization comparison rests on.
            alpha = log(5.5); s_p = 0.45; s_i = 0.20; sigma = 1.7;
            ssum  = s_p^2 + s_i^2;
            nu_equivalent = 2 * exp(2*alpha + 2*ssum) / sigma^2;

            ls = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, log(sigma));
            nu = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_equivalent);

            testCase.verifyEqual(ls.e_cond_var,  nu.e_cond_var,  'RelTol', 1e-12);
            testCase.verifyEqual(ls.total_var,   nu.total_var,   'RelTol', 1e-12);
            testCase.verifyEqual(ls.sigma_pi_e2, nu.sigma_pi_e2, 'RelTol', 1e-12);
        end

        function testScaledChiSquareMeanSdMomentIdentity(testCase)
            % The parameterization rests on nu = 2*(mu/sigma)^2 giving a scaled
            % chi-square with mean mu and SD sigma. Verified generatively, in
            % base MATLAB only: for INTEGER nu, Y = (mu/nu)*chi2(nu) and chi2(nu)
            % is a sum of nu squared standard normals, so no Statistics Toolbox
            % gamma sampler is needed. This mirrors the reference bundle's own
            % moment check (tests/run_static_validation.R) in the other language.
            mu = 5; sigma = 2.5;
            nu = 2 * (mu/sigma)^2;               % = 8, integer by construction
            testCase.verifyEqual(nu, 8, 'AbsTol', 1e-12);

            rng(4043);
            n = 2e6;
            chi2 = sum(randn(n, nu).^2, 2);
            y = (mu/nu) * chi2;

            testCase.verifyEqual(mean(y), mu,    'RelTol', 5e-3);
            testCase.verifyEqual(std(y),  sigma, 'RelTol', 5e-3);
        end

    end
end

function r = corr_base(x, y)
%Pearson correlation without the Statistics Toolbox (corr lives there; CI
%installs base MATLAB only).
x = x(:) - mean(x(:));
y = y(:) - mean(y(:));
r = (x' * y) / sqrt((x' * x) * (y' * y));
end
