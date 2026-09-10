classdef TestGammaVarcompsConversion < matlab.unittest.TestCase
    %TESTGAMMAVARCOMPSCONVERSION Unit tests for the scaled-chi-square / Gamma
    % log-scale -> observed-score-scale variance-component conversion
    % (psyrat_gamma_varcomps). The conversion is the scientific heart of the
    % Gamma family: because the mean sub-model uses a log link, the Stan
    % person/trial SDs are on the log expected-score scale and must be mapped to
    % observed-score-scale G-theory variance components via lognormal moments
    % BEFORE any reliability formula consumes them.
    %
    % These tests use INDEPENDENT derivations of the target quantities (not a
    % copy of the production expression) so a green test means the production
    % formula is right, not merely unchanged:
    %   - the degenerate zero-random-effects case has a trivial closed form;
    %   - the person/trial components are checked against the "variance of the
    %     conditional mean" lognormal identity;
    %   - a fixed-seed Monte Carlo check confirms the analytic formula matches
    %     the generative model within sampling error.

    methods (Test)

        function testDegenerateZeroRandomEffects(testCase)
            % With no person/trial variation (s_p = s_i = 0) the mean surface is
            % constant at exp(alpha), so sigma_p2 = sigma_i2 = var_mu = 0 and the
            % only observed variance is the scaled-chi-square conditional
            % dispersion E[Var(Y|mu,nu)] = (2/nu)*exp(2*alpha).
            alpha = log(3);
            nu    = 20;
            vc = psyrat_gamma_varcomps(alpha, 0, 0, nu);

            expected_cond = (2/nu) * exp(2*alpha);   % = 0.9
            testCase.verifyEqual(vc.sigma_p2,    0,             'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_i2,    0,             'AbsTol', 1e-12);
            testCase.verifyEqual(vc.var_mu,      0,             'AbsTol', 1e-12);
            testCase.verifyEqual(vc.e_cond_var,  expected_cond, 'AbsTol', 1e-12);
            testCase.verifyEqual(vc.total_var,   expected_cond, 'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_pi_e2, expected_cond, 'AbsTol', 1e-12);
            testCase.verifyEqual(vc.mu_bar,      exp(alpha),    'AbsTol', 1e-12);
        end

        function testPersonTrialAgainstConditionalMeanIdentity(testCase)
            % Independent check of the person/trial components. The person
            % component equals Var_p(E_t[mu|u_p]); E_t[mu|u_p] = exp(alpha + u_p
            % + si2/2) is lognormal in u_p with log-mean (alpha + si2/2) and
            % log-variance sp2, so its variance is
            %   (exp(sp2) - 1) * exp(2*(alpha + si2/2) + sp2).
            % This is written differently from the production expression, so a
            % match catches transcription errors. Vectorized over several draws.
            alpha = [log(5.5);  log(6.0);  log(4.2)];
            s_p   = [0.40;      0.25;      0.55];
            s_i   = [0.10;      0.18;      0.05];
            nu    = [18;        22;        14];

            vc = psyrat_gamma_varcomps(alpha, s_p, s_i, nu);

            sp2 = s_p.^2; si2 = s_i.^2;
            expected_p = (exp(sp2) - 1) .* exp(2.*(alpha + si2./2) + sp2);
            expected_i = (exp(si2) - 1) .* exp(2.*(alpha + sp2./2) + si2);

            testCase.verifyEqual(vc.sigma_p2, expected_p, 'RelTol', 1e-12);
            testCase.verifyEqual(vc.sigma_i2, expected_i, 'RelTol', 1e-12);

            % Total-variance additive identity must hold exactly whenever the
            % residual is not clipped: sigma_p2 + sigma_i2 + sigma_pi_e2 = total.
            testCase.verifyEqual(vc.sigma_p2 + vc.sigma_i2 + vc.sigma_pi_e2, ...
                vc.total_var, 'RelTol', 1e-12);
        end

        function testMonteCarloMeanSurface(testCase)
            % Fixed-seed Monte Carlo validation of the mean-surface components
            % against the generative model, using only base-MATLAB randn (no
            % Statistics Toolbox). This is genuinely independent of the analytic
            % derivation: it would catch a wrong formula, not just a changed one.
            alpha = log(5.5); s_p = 0.45; s_i = 0.20; nu = 18;
            vc = psyrat_gamma_varcomps(alpha, s_p, s_i, nu);

            rng(2024);
            n = 4e6;
            up = s_p * randn(n,1);
            ut = s_i * randn(n,1);
            mu = exp(alpha + up + ut);                 % lognormal mean surface

            % E[mu] and Var(mu) of the mean surface.
            testCase.verifyEqual(mean(mu), vc.mu_bar, 'RelTol', 5e-3);
            testCase.verifyEqual(var(mu,1), vc.var_mu, 'RelTol', 1e-2);

            % Person component = Var over persons of the trial-marginal mean.
            % E_t[mu|u_p] = exp(alpha + u_p + si2/2); MC over persons.
            up2 = s_p * randn(n,1);
            cond_mean_p = exp(alpha + up2 + (s_i^2)/2);
            testCase.verifyEqual(var(cond_mean_p,1), vc.sigma_p2, 'RelTol', 1e-2);

            % Trial component = Var over trials of the person-marginal mean.
            ut2 = s_i * randn(n,1);
            cond_mean_t = exp(alpha + ut2 + (s_p^2)/2);
            testCase.verifyEqual(var(cond_mean_t,1), vc.sigma_i2, 'RelTol', 1e-2);
        end

    end
end
