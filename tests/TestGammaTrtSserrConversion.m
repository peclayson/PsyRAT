classdef TestGammaTrtSserrConversion < matlab.unittest.TestCase
    %TESTGAMMATRTSSERRCONVERSION Offline proof of the Gamma (scaled-chi-square)
    % SUBJECT-LEVEL crossed TEST-RETEST (trt_sserr, analysis 25) observed-scale
    % conversion. No CmdStan.
    %
    % The gamma subject-level test-retest design reuses the crossed gamma TRT model
    % (psyrat_build_stan_trt_sserr_chisq == psyrat_build_stan_trt_chisq): the
    % person-varying, log-linked dispersion log nu_p = beta + v_p already gives each
    % participant their own scaled-chi-square residual, so no location-scale
    % residual submodel is needed. The new work is the EXTRACTION
    % (psyrat_gamma_extract_trt_sserr), which CONDITIONS on each person's own
    % log-mean (u_p) and log-nu (v_p) to form the per-subject residual
    %
    %   sigma_pot_e2(p) = 2*exp(2*(alpha + u_p) + 2*V_w - (beta + v_p)),
    %   V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2,
    %
    % the expected WITHIN-person observation variance E_{o,t,po,pt,ot}[2*mu^2/nu_p].
    % V_w is the sum of every log-mean variance EXCEPT the person main effect s_p^2,
    % which is dropped because u_p has been conditioned on rather than integrated
    % out. This is the case-6 one-facet formula with 2*s_i^2 generalized to 2*V_w,
    % and it is the genuinely new mathematics in analysis 25 - hence
    % testPerSubjectResidualMatchesConditionalMonteCarlo is the decisive test here.
    %
    % Exact under the ID-only nu structure (nu is subject-constant, so there is no
    % lognormal-approximation error). It is the subject-level COMPLEMENT of the
    % group estimand #2 (psyrat_gamma_varcomps_trt.e_cond_var), which MARGINALIZES
    % the same quantity over the person dispersion population.
    %
    % All draws use base-MATLAB randn only (no Statistics Toolbox), matching
    % TestGammaSserrConversion: E[2*mu^2/nu_p] is a function of (mu, nu_p) only, so
    % mean(2*mu.^2 ./ nu_p) is its exact ground truth with no Y draw and no
    % observation sampling noise.

    methods (Static, Access = private)

        function [up, vp] = drawPersonBlock(s_p, s_v, rho, n)
            % Correlated person (mean, log-nu) block via base randn:
            %   Var(up)=s_p^2, Var(vp)=s_v^2, Cov(up,vp)=rho*s_p*s_v.
            z1 = randn(n,1);
            z2 = randn(n,1);
            up = s_p .* z1;
            vp = s_v .* (rho .* z1 + sqrt(1 - rho^2) .* z2);
        end

        function r = condResid(alpha, u_p, V_w, beta, v_p)
            % A RE-IMPLEMENTATION of the formula psyrat_gamma_extract_trt_sserr
            % implements. It is NOT the production function: that is a local
            % function of psyrat_computevarcomp and is unreachable from a test.
            %
            % So be clear about what this file does and does not prove. It proves
            % the MATHEMATICS is right (against Monte-Carlo ground truth, and
            % against the shared psyrat_gamma_varcomps_trt helper). It cannot
            % prove the production extractor implements that mathematics. What
            % covers the IMPLEMENTATION is tests/integration/TestGammaTrtSserrRecovery,
            % which runs the real extractor under CmdStan and asserts the recovered
            % per-subject residual against the planted truth on a scale that
            % discriminates a wrong V_w.
            r = 2 .* exp(2.*(alpha + u_p) + 2.*V_w - (beta + v_p));
        end

    end

    methods (Test)

        function testPerSubjectResidualMatchesConditionalMonteCarlo(testCase)
            % DECISIVE. For a fixed person (u_p, v_p) the per-subject residual
            % formula must equal the Monte-Carlo WITHIN-person expected observation
            % variance E[2*mu^2/nu_p] taken over ALL FIVE within-person log-mean
            % facets (occasion, trial, person x occasion, person x trial,
            % trial x occasion). This is what proves V_w is the right sum, and that
            % the person main effect s_p^2 is correctly EXCLUDED from it.
            alpha = log(5.5); beta = log(18);
            s_o = 0.22; s_t = 0.30; s_po = 0.18; s_pt = 0.26; s_ot = 0.12;
            V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2;

            % persons spanning the (u_p, v_p) space
            persons = [ 0.60 -0.50;  -0.40 0.70;  0.00 0.00;  0.25 0.35 ];
            rng(2501);
            T = 4e6;
            for k = 1:size(persons,1)
                u_p = persons(k,1); v_p = persons(k,2);
                nu_p = exp(beta + v_p);
                % independent draws of every within-person facet
                mu = exp(alpha + u_p ...
                    + s_o  * randn(T,1) ...
                    + s_t  * randn(T,1) ...
                    + s_po * randn(T,1) ...
                    + s_pt * randn(T,1) ...
                    + s_ot * randn(T,1));
                mc  = mean(2 .* mu.^2 ./ nu_p);                      % ground truth
                fml = TestGammaTrtSserrConversion.condResid(alpha, u_p, V_w, beta, v_p);
                testCase.verifyEqual(fml, mc, 'RelTol', 5e-3);
            end
        end

        function testReducesToOneFacetWhenOnlyTrialFacetIsPresent(testCase)
            % REDUCTION. With every crossed facet except the trial main effect set
            % to zero, V_w collapses to s_t^2 and the two-facet conditional residual
            % equals the case-6 one-facet formula (psyrat_gamma_extract_sserr)
            % exactly, with s_i = s_t.
            %
            % LIKE the round-trip test below, this is an algebraic identity given
            % the shared helper - it documents that the generalization is a strict
            % superset of case 6, and would fail only if the two formulas were
            % restated inconsistently HERE. It is not evidence about production.
            alpha = log(5.0); beta = log(16); s_i = 0.30;
            V_w = 0^2 + s_i^2 + 0^2 + 0^2 + 0^2;
            persons = [ 0.60 -0.50;  -0.40 0.70;  0.00 0.00;  -0.15 -0.25 ];
            for k = 1:size(persons,1)
                u_p = persons(k,1); v_p = persons(k,2);
                twoFacet = TestGammaTrtSserrConversion.condResid(alpha, u_p, V_w, beta, v_p);
                oneFacet = 2 * exp(2*(alpha + u_p) + 2*s_i^2 - (beta + v_p));
                testCase.verifyEqual(twoFacet, oneFacet, 'RelTol', 1e-14);
            end
        end

        function testPopulationAverageMatchesMarginalEcondVar(testCase)
            % DECISIVE CONSISTENCY. Averaging the per-subject CONDITIONAL residual
            % over the person population must recover the group MARGINAL observation
            % term (estimand #2) that psyrat_gamma_varcomps_trt computes, by the law
            % of total expectation:  E_p[ sigma_pot_e2(p) ] = e_cond_var(#2).
            % This ties the new subject-level extraction to the existing group-level
            % helper the same way TestGammaSserrConversion does for case 6.
            alpha = log(5.0); s_p = 0.45;
            s_o = 0.20; s_t = 0.25; s_po = 0.15; s_pt = 0.22; s_ot = 0.10;
            beta = log(16); s_v = 0.50; rho = 0.55;
            nu_pop = exp(beta);
            cov_pv = rho * s_p * s_v;
            V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2;

            rng(2502);
            n = 5e6;
            [up, vp] = TestGammaTrtSserrConversion.drawPersonBlock(s_p, s_v, rho, n);
            popAvg = mean(TestGammaTrtSserrConversion.condResid(alpha, up, V_w, beta, vp));

            vc = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, ...
                nu_pop, s_v, cov_pv);
            testCase.verifyEqual(popAvg, vc.e_cond_var, 'RelTol', 1e-2);
        end

        function testConditionalResidualDiffersFromMarginalSigmaRes2(testCase)
            % GUARD against a tempting "reconciliation". psyrat_gamma_varcomps_trt
            % recovers sigma_res2 BY SUBTRACTION from total_var, so it bundles the
            % higher-order lognormal leftover of the mean surface together with the
            % dispersion. The conditional residual is the PURE per-person dispersion.
            % The two are different estimands on purpose (identical, maintainer-
            % sanctioned mismatch to case 6 vs case 1) - assert they differ so a
            % future change that "fixes" one into the other fails loudly.
            alpha = log(5.0); s_p = 0.45;
            s_o = 0.20; s_t = 0.25; s_po = 0.15; s_pt = 0.22; s_ot = 0.10;
            beta = log(16);
            nu_pop = exp(beta);
            V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2;

            % typical subject (u_p = v_p = 0), no dispersion heterogeneity
            vc = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, ...
                nu_pop, 0, 0);
            typical = TestGammaTrtSserrConversion.condResid(alpha, 0, V_w, beta, 0);

            testCase.verifyGreaterThan(vc.sigma_res2, typical, ...
                ['sigma_res2 is subtraction-derived and must exceed the pure ' ...
                'conditional dispersion; if these ever coincide the derivation ' ...
                'has changed - see psyrat_gamma_extract_trt_sserr''s header.']);
        end

        function testSsrelTrtLogSdEncodingRoundTrips(testCase)
            % The extractor encodes the per-subject residual as log-SD-equivalents
            % (pop_sdlog + ind_sdlog) so that psyrat_ssrel_trt's
            % err_ss = exp(pop_sdlog + ind_sdlog) reproduces sqrt(sigma_pot_e2(p)).
            %
            % NOTE this is an algebraic identity by construction, not a test of the
            % extractor: it documents the encoding contract and would fail only if
            % the encoding were restated wrongly HERE. The contract is exercised
            % for real by TestGammaTrtSserrRecovery, which feeds the extractor's
            % actual pop_sdlog/ind_sdlog through psyrat_relsummary.
            alpha = log(5.5); beta = log(18);
            s_o = 0.22; s_t = 0.30; s_po = 0.18; s_pt = 0.26; s_ot = 0.12;
            V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2;

            persons = [ 0.60 -0.50;  -0.40 0.70;  0.00 0.00;  -0.15 -0.25 ];
            sig_pop = TestGammaTrtSserrConversion.condResid(alpha, 0, V_w, beta, 0);
            pop_sdlog = 0.5 * log(sig_pop);
            for k = 1:size(persons,1)
                u_p = persons(k,1); v_p = persons(k,2);
                sig_p_e2 = TestGammaTrtSserrConversion.condResid(alpha, u_p, V_w, beta, v_p);
                ind_sdlog = 0.5 * log(sig_p_e2) - pop_sdlog;
                % psyrat_relsummary forms err_ss = exp(pop_sdlog + ind_sdlog)
                err_ss = exp(pop_sdlog + ind_sdlog);
                testCase.verifyEqual(err_ss^2, sig_p_e2, 'RelTol', 1e-13);
            end
        end

        function testSdVsVarianceConfusionWouldBeVisible(testCase)
            % RC-02 (commit 3b51b23) requires the extractor to emit gro_sds(:,1) as
            % a STANDARD DEVIATION: psyrat_variancet squares it and psyrat_ssrel_trt
            % squares its bp slot, so emitting a variance would report sigma_p^4.
            %
            % This test canNOT pin that contract on the extractor (a local function
            % of psyrat_computevarcomp, unreachable from here) - the live recovery
            % test does, by checking the recovered person component against its
            % planted truth. What this test DOES establish is that the two are far
            % enough apart to be caught: it asserts sqrt(sigma_p2) and sigma_p2
            % differ by a wide margin at realistic parameter values, so the
            % recovery test's bounds are meaningful rather than coincidentally
            % satisfied. (An earlier version asserted sqrt(x)^2 == x, which is an
            % identity and could never fail.)
            alpha = log(5.0); s_p = 0.45;
            s_o = 0.20; s_t = 0.25; s_po = 0.15; s_pt = 0.22; s_ot = 0.10;
            nu_pop = exp(log(16));
            vc = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, ...
                nu_pop, 0, 0);
            sdVal  = sqrt(vc.sigma_p2);   % the correct emission
            varVal = vc.sigma_p2;         % the RC-02 error
            testCase.verifyGreaterThan(abs(varVal - sdVal) / sdVal, 0.5, ...
                ['SD and variance are too close at these parameters for the ', ...
                'recovery test to discriminate them; choose different values.']);
        end

        function testSubjectLevelAssemblyMatchesRelTrtPartition(testCase)
            % The per-subject coefficient must equal psyrat_rel_trt's crossed
            % partition when fed the observed-scale SD-equivalents the extractor
            % emits, with the person's own residual substituted. Exercises the real
            % downstream kernel rather than re-deriving it in the test.
            alpha = log(5.0); s_p = 0.45;
            s_o = 0.20; s_t = 0.25; s_po = 0.15; s_pt = 0.22; s_ot = 0.10;
            beta = log(16); nu_pop = exp(beta);
            V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2;
            n_i = 20; n_o = 2;

            vc = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, ...
                nu_pop, 0, 0);

            u_p = 0.35; v_p = -0.20;
            err_ss = sqrt(TestGammaTrtSserrConversion.condResid(alpha, u_p, V_w, beta, v_p));

            % SD-equivalents exactly as psyrat_gamma_extract_trt_sserr emits them,
            % including the Stan-name/varcomps-name transposition:
            %   oxp <- sigma_po2 (sig_occxid), txp <- sigma_pt2 (sig_trlxid)
            args = {'gcoeff',2,'reltype',3, ...
                'bp',sqrt(vc.sigma_p2),'bo',sqrt(vc.sigma_o2),'bt',sqrt(vc.sigma_t2), ...
                'txp',sqrt(vc.sigma_pt2),'oxp',sqrt(vc.sigma_po2), ...
                'txo',sqrt(vc.sigma_ot2),'err',err_ss, ...
                'obs',n_i,'nocc',n_o,'CI',.95};
            [~, G_kernel, ~] = psyrat_rel_trt(args{:});

            % reference: crossed-events-and-strata generalizability (Rocha Table 3)
            bp = vc.sigma_p2; txp = vc.sigma_pt2; oxp = vc.sigma_po2;
            rel_err = txp/n_i + oxp/n_o + err_ss^2/(n_i*n_o);
            G_ref = bp / (bp + rel_err);

            testCase.verifyEqual(G_kernel, G_ref, 'RelTol', 1e-10);
        end

    end
end
