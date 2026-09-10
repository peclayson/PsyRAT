classdef TestGammaSserrConversion < matlab.unittest.TestCase
    %TESTGAMMASSERRCONVERSION Offline proof of the Gamma (scaled-chi-square)
    % SUBJECT-LEVEL (sserr, analysis 6) observed-scale conversion. No CmdStan.
    %
    % The gamma subject-level design reuses the one-facet gamma model
    % (psyrat_build_stan_sserr_chisq == psyrat_build_stan_single_chisq): the
    % person-varying, log-linked dispersion log nu_p = beta + v_p already gives
    % each participant their own scaled-chi-square residual. The new work is the
    % EXTRACTION (psyrat_gamma_extract_sserr), which CONDITIONS on each person's
    % own log-mean (u_p) and log-nu (v_p) to form the per-subject residual
    %       sigma_pi_e2(p) = 2*exp(2*(alpha + u_p) + 2*s_i^2 - (beta + v_p)),
    % the expected within-person observation variance E_t[2*mu_pt^2/nu_p]. This is
    % EXACT under the ID-only nu structure (nu is subject-constant, so there is no
    % lognormal-approximation error). It is the subject-level COMPLEMENT of the
    % group estimand #2 (psyrat_gamma_varcomps.e_cond_var), which MARGINALIZES the
    % same quantity over the person dispersion population.
    %
    % All draws use base-MATLAB randn only (no Statistics Toolbox), matching
    % TestGammaDispersionMarginalization: E_t[2*mu^2/nu_p] is a function of
    % (mu, nu_p) only, so mean(2*mu.^2 ./ nu_p) is its exact ground truth with no
    % Y draw and no observation sampling noise.

    methods (Static, Access = private)

        function [up, vp] = drawPersonBlock(s_p, s_v, rho, n)
            % Correlated person (mean, log-nu) block via base randn:
            %   Var(up)=s_p^2, Var(vp)=s_v^2, Cov(up,vp)=rho*s_p*s_v.
            z1 = randn(n,1);
            z2 = randn(n,1);
            up = s_p .* z1;
            vp = s_v .* (rho .* z1 + sqrt(1 - rho^2) .* z2);
        end

    end

    methods (Test)

        function testPerSubjectResidualMatchesConditionalMonteCarlo(testCase)
            % DECISIVE. For fixed persons (u_p, v_p) the per-subject residual
            % formula must equal the Monte-Carlo within-person expected
            % observation variance E_t[2*mu_pt^2/nu_p] over the trial effects.
            alpha = log(5.5); s_i = 0.30; beta = log(18);
            % persons spanning the (u_p, v_p) space
            persons = [ 0.60 -0.50;  -0.40 0.70;  0.00 0.00;  0.25 0.35 ];
            rng(6006);
            T = 4e6;
            for k = 1:size(persons,1)
                u_p = persons(k,1); v_p = persons(k,2);
                nu_p = exp(beta + v_p);
                u_t  = s_i * randn(T,1);
                mu_pt = exp(alpha + u_p + u_t);
                mc  = mean(2 .* mu_pt.^2 ./ nu_p);                        % ground truth
                fml = 2 * exp(2*(alpha + u_p) + 2*s_i^2 - (beta + v_p));  % formula
                testCase.verifyEqual(fml, mc, 'RelTol', 5e-3);
            end
        end

        function testPopulationAverageMatchesMarginalEcondVar(testCase)
            % DECISIVE CONSISTENCY. Averaging the per-subject CONDITIONAL residual
            % over the person population must recover the group MARGINAL
            % observation term (estimand #2) from psyrat_gamma_varcomps, by the law
            % of total expectation:
            %   E_p[ sigma_pi_e2(p) ] = e_cond_var(#2).
            % This ties the subject-level (conditional) extraction to the existing
            % group (marginal) helper.
            alpha = log(5.0); s_p = 0.45; s_i = 0.20;
            beta  = log(16);  s_v = 0.50; rho = 0.55;
            nu_pop = exp(beta);
            cov_pv = rho * s_p * s_v;

            rng(6007);
            n = 5e6;
            [up, vp] = TestGammaSserrConversion.drawPersonBlock(s_p, s_v, rho, n);
            % per-subject conditional residual for each drawn person
            sig_pi_e2 = 2 .* exp(2*(alpha + up) + 2*s_i^2 - (beta + vp));
            popAvg = mean(sig_pi_e2);

            vc = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_pop, s_v, cov_pv);
            testCase.verifyEqual(popAvg, vc.e_cond_var, 'RelTol', 1e-2);
        end

        function testReliabilityAssemblyMatchesSsrelFormula(testCase)
            % The Rocha Table 3 one-facet per-subject reliability (numerator = the
            % GROUP person universe variance sigma_p2; denominator = that + the
            % per-subject residual, projected over n_i) must equal the psyrat_ssrel
            % arithmetic bp/(bp + errwp/n_i) when fed the observed-scale
            % SD-equivalents the extractor emits. Pure arithmetic, machine tol.
            alpha = log(5.5); s_p = 0.40; s_i = 0.30;
            beta  = log(18);  s_v = 0.50; rho = -0.30;
            nu_pop = exp(beta); cov_pv = rho * s_p * s_v;
            n_i = 30;

            vc = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_pop, s_v, cov_pv);
            sigma_p2 = vc.sigma_p2; sigma_i2 = vc.sigma_i2;

            persons = [ 0.60 -0.50;  -0.40 0.70;  0.00 0.00 ];
            for k = 1:size(persons,1)
                u_p = persons(k,1); v_p = persons(k,2);
                sig_pi_e2 = 2 * exp(2*(alpha + u_p) + 2*s_i^2 - (beta + v_p));
                % reference subject_level_coefficients (one-facet: n_o = 1)
                rel_err = sig_pi_e2 / n_i;
                abs_err = rel_err + sigma_i2 / n_i;
                G_ref = sigma_p2 / (sigma_p2 + rel_err);   % generalizability
                D_ref = sigma_p2 / (sigma_p2 + abs_err);   % dependability
                % psyrat_ssrel arithmetic with bp = sigma_p2, i2 = sigma_i2
                bp = sigma_p2;
                G_ssrel = bp / (bp + sig_pi_e2 / n_i);
                D_ssrel = bp / (bp + (sig_pi_e2 + sigma_i2) / n_i);
                testCase.verifyEqual(G_ssrel, G_ref, 'RelTol', 1e-12);
                testCase.verifyEqual(D_ssrel, D_ref, 'RelTol', 1e-12);
            end
        end

        function testSsrelLogSdEncodingRoundTrips(testCase)
            % The extractor encodes the per-subject residual as log-SD-equivalents
            % (pop_sdlog + ind_sdlog) so that psyrat_ssrel's wp = exp(wp_pop +
            % wp_ss).^2 reproduces sigma_pi_e2(p) exactly.
            alpha = log(5.5); s_i = 0.30; beta = log(18);
            persons = [ 0.60 -0.50;  -0.40 0.70;  0.00 0.00;  -0.15 -0.25 ];
            sig_pi_e2_pop = 2 * exp(2*alpha + 2*s_i^2 - beta);   % typical subject u=v=0
            pop_sdlog = 0.5 * log(sig_pi_e2_pop);
            for k = 1:size(persons,1)
                u_p = persons(k,1); v_p = persons(k,2);
                sig_pi_e2 = 2 * exp(2*(alpha + u_p) + 2*s_i^2 - (beta + v_p));
                ind_sdlog = 0.5 * log(sig_pi_e2) - pop_sdlog;
                roundtrip = exp(pop_sdlog + ind_sdlog)^2;
                testCase.verifyEqual(roundtrip, sig_pi_e2, 'RelTol', 1e-13);
            end
        end

    end
end
