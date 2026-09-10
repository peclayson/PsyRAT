classdef TestGammaDispersionMarginalization < matlab.unittest.TestCase
    %TESTGAMMADISPERSIONMARGINALIZATION Monte-Carlo proof that the Gamma
    % observed-scale conversion uses the MARGINAL observation-variance term
    % (estimand #2), averaging over the person population's dispersion
    % heterogeneity, rather than plugging in the population dispersion nu_pop
    % (estimand #1).
    %
    % The fitted model has a person-varying log-dispersion, log nu_p = beta + v_p,
    % with v_p correlated with the person MEAN effect at correlation rho. The
    % observation-level conditional variance is Var(Y|mu,nu_p) = 2*mu^2/nu_p, so
    % the group-level observation term is the population average
    %       e_cond_var(#2) = E[ 2*mu^2 / nu_p ],
    % taken over the joint distribution of (mu, nu_p). Because this term is a
    % function of (mu, nu_p) ONLY, its Monte-Carlo ground truth is
    % mean(2*mu.^2 ./ nu_p) from fixed-seed randn draws - no Y draw, no
    % Statistics Toolbox, and no sampling noise from the observation itself.
    %
    % Estimand #1 (population nu) omits a multiplicative factor
    %       F = exp( 0.5*s_v^2 - 2*cov_pv ),   cov_pv = Cov(u_p, v_p),
    % where s_v is the person log-nu SD. These tests confirm #2 matches the
    % ground truth while #1 is off by exactly F, and that F is design-invariant
    % (it depends only on the person mean<->log-nu coupling, never on the
    % occasion/trial/interaction facets or on cross-cell correlations).
    %
    % All draws use base-MATLAB randn only (no Statistics Toolbox), matching
    % TestGammaVarcompsConversion.

    methods (Static, Access = private)

        function [up, vp] = drawPersonBlock(s_p, s_v, rho, n)
            % Draw a correlated person (mean, log-nu) block with base randn:
            %   Var(up)=s_p^2, Var(vp)=s_v^2, Cov(up,vp)=rho*s_p*s_v.
            z1 = randn(n,1);
            z2 = randn(n,1);
            up = s_p .* z1;
            vp = s_v .* (rho .* z1 + sqrt(1 - rho^2) .* z2);
        end

    end

    methods (Test)

        function testOneFacetMarginalMatchesMonteCarlo(testCase)
            % DECISIVE TEST. The marginal observation term from the production
            % helper (with the person dispersion SD and covariance passed in)
            % must equal the Monte-Carlo E[2*mu^2/nu_p] of the generative model.
            alpha = log(5.5); s_p = 0.45; s_i = 0.20;
            beta  = log(18);  s_v = 0.50; rho = 0.60;
            nu_pop = exp(beta);
            cov_pv = rho * s_p * s_v;

            rng(90210);
            n = 4e6;
            [up, vp] = TestGammaDispersionMarginalization.drawPersonBlock(s_p, s_v, rho, n);
            ut = s_i * randn(n,1);
            mu   = exp(alpha + up + ut);
            nu_p = exp(beta + vp);
            mc_e_cond = mean(2 .* mu.^2 ./ nu_p);     % exact ground truth for E[Var(Y|.)]

            % #2: marginal conversion (new trailing args s_v, cov_pv).
            vc2 = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_pop, s_v, cov_pv);
            testCase.verifyEqual(vc2.e_cond_var, mc_e_cond, 'RelTol', 1e-2);

            % The mean-surface variance is unchanged by #2, so total_var(#2)
            % must equal Var(mu) + E[2*mu^2/nu_p] = the true marginal Var(Y).
            mc_var_mu = var(mu, 1);
            testCase.verifyEqual(vc2.var_mu, mc_var_mu, 'RelTol', 1e-2);
            testCase.verifyEqual(vc2.total_var, mc_var_mu + mc_e_cond, 'RelTol', 1e-2);
        end

        function testPopulationNuIsBiasedByExactFactor(testCase)
            % Teeth: the OLD population-nu conversion (#1, no dispersion args) is
            % off from the ground truth by exactly F, and is NOT within 1%.
            alpha = log(5.5); s_p = 0.45; s_i = 0.20;
            beta  = log(18);  s_v = 0.50; rho = 0.60;
            nu_pop = exp(beta);
            cov_pv = rho * s_p * s_v;
            F = exp(0.5 * s_v^2 - 2 * cov_pv);

            rng(90211);
            n = 4e6;
            [up, vp] = TestGammaDispersionMarginalization.drawPersonBlock(s_p, s_v, rho, n);
            ut = s_i * randn(n,1);
            mu   = exp(alpha + up + ut);
            nu_p = exp(beta + vp);
            mc_e_cond = mean(2 .* mu.^2 ./ nu_p);

            vc1 = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_pop);   % #1 (defaults)

            % #1 * F recovers the truth ...
            testCase.verifyEqual(vc1.e_cond_var * F, mc_e_cond, 'RelTol', 1e-2);
            % ... and #1 alone is materially biased (here F ~ 1.14, ~14% low).
            relErr = abs(vc1.e_cond_var - mc_e_cond) / mc_e_cond;
            testCase.verifyGreaterThan(relErr, 0.10);
        end

        function testDegenerateReducesToPopulationNu(testCase)
            % Continuity / backward-compat: with no dispersion heterogeneity
            % (s_v = 0), #2 must equal #1 bit-for-bit, and passing cov_pv is
            % irrelevant because the factor is exp(0) = 1.
            alpha = [log(5.5); log(6.0); log(4.2)];
            s_p   = [0.40; 0.25; 0.55];
            s_i   = [0.10; 0.18; 0.05];
            nu    = [18; 22; 14];

            vc1 = psyrat_gamma_varcomps(alpha, s_p, s_i, nu);
            vc2 = psyrat_gamma_varcomps(alpha, s_p, s_i, nu, 0, 0);
            testCase.verifyEqual(vc2.e_cond_var, vc1.e_cond_var, 'RelTol', 1e-14);
            testCase.verifyEqual(vc2.total_var,  vc1.total_var,  'RelTol', 1e-14);
            testCase.verifyEqual(vc2.sigma_pi_e2, vc1.sigma_pi_e2, 'RelTol', 1e-14);
        end

        function testDesignInvarianceOfFactor(testCase)
            % The correction factor depends only on the person mean<->log-nu
            % coupling (s_p, s_v, rho), never on the trial SD. Two runs with the
            % same (s_p, s_v, rho) but different s_i must show the SAME measured
            % ratio mc_e_cond / e_cond_var(#1) = F.
            alpha = log(5.0); s_p = 0.40; s_v = 0.45; rho = 0.50;
            cov_pv = rho * s_p * s_v;
            F = exp(0.5 * s_v^2 - 2 * cov_pv);
            nu_pop = exp(log(16));

            n = 3e6;
            for s_i = [0.10, 0.35]
                rng(4242);
                [up, vp] = TestGammaDispersionMarginalization.drawPersonBlock(s_p, s_v, rho, n);
                ut = s_i * randn(n,1);
                mu   = exp(alpha + up + ut);
                nu_p = exp(log(16) + vp);
                mc_e_cond = mean(2 .* mu.^2 ./ nu_p);
                vc1 = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_pop);
                testCase.verifyEqual(mc_e_cond / vc1.e_cond_var, F, 'RelTol', 1.5e-2);
                % And the marginal helper matches the truth at this s_i too.
                vc2 = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_pop, s_v, cov_pv);
                testCase.verifyEqual(vc2.e_cond_var, mc_e_cond, 'RelTol', 1.5e-2);
            end
        end

        function testPerCellIndependence(testCase)
            % In the difference / DoD designs each cell is converted with the
            % one-facet helper, so a cell's e_cond_var must depend ONLY on that
            % cell's own (s_v_c, cov(mean_c,lognu_c)) - never on cross-cell
            % person correlations. Simulate a 4-D person block
            % [u_p1, u_p2, v_p1, v_p2] with strong cross-event correlations and
            % confirm each cell matches its own one-facet marginal conversion.
            alpha = [log(5.0), log(4.0)];       % per-event log intercepts
            s_p   = [0.40, 0.30];               % per-event person mean SD
            s_i   = [0.20, 0.15];               % per-event trial SD
            s_v   = [0.45, 0.55];               % per-event person log-nu SD
            beta  = [log(15), log(20)];         % per-event log-nu intercept
            rho_c = [0.50, 0.35];               % within-cell mean<->lognu corr

            % Person correlation matrix over [u_p1, u_p2, v_p1, v_p2], with
            % nonzero cross-event correlations that must NOT affect a cell's own
            % e_cond_var.
            R = eye(4);
            R(1,2) = 0.60; R(2,1) = 0.60;       % cross-event person-mean corr
            R(3,4) = 0.40; R(4,3) = 0.40;       % cross-event log-nu corr
            R(1,3) = rho_c(1); R(3,1) = rho_c(1);   % cell-1 mean<->lognu
            R(2,4) = rho_c(2); R(4,2) = rho_c(2);   % cell-2 mean<->lognu
            R(1,4) = 0.20; R(4,1) = 0.20;       % cross mean1<->lognu2
            R(2,3) = 0.15; R(3,2) = 0.15;       % cross mean2<->lognu1
            sd = [s_p(1), s_p(2), s_v(1), s_v(2)];
            Sigma = (sd' * sd) .* R;            % person covariance
            L = chol(Sigma, 'lower');           % base-MATLAB Cholesky (SPD)

            rng(1234);
            n = 3e6;
            Z = randn(4, n);
            eff = L * Z;                        % 4 x n person effects
            for c = 1:2
                upc = eff(c, :)';               % person mean effect, cell c
                vpc = eff(c + 2, :)';           % person log-nu effect, cell c
                utc = s_i(c) * randn(n,1);      % independent trial effect
                mu   = exp(alpha(c) + upc + utc);
                nu_p = exp(beta(c) + vpc);
                mc_e_cond = mean(2 .* mu.^2 ./ nu_p);

                cov_pv = Sigma(c, c + 2);       % Cov(mean_c, lognu_c)
                vc = psyrat_gamma_varcomps(alpha(c), s_p(c), s_i(c), ...
                    exp(beta(c)), s_v(c), cov_pv);
                testCase.verifyEqual(vc.e_cond_var, mc_e_cond, 'RelTol', 1.5e-2);
            end
        end

        function testTrtMarginalMatchesMonteCarlo(testCase)
            % Two-facet (crossed test-retest) design: six independent mean
            % effects plus a person log-nu effect correlated ONLY with the
            % person mean. The marginal observation term must match
            % E[2*mu^2/nu_p], with the SAME factor form as the one-facet case.
            alpha = log(5.0);
            s_p = 0.35; s_o = 0.25; s_t = 0.20;
            s_po = 0.18; s_pt = 0.15; s_ot = 0.12;
            beta = log(16); s_v = 0.50; rho = 0.55;
            nu_pop = exp(beta);
            cov_pv = rho * s_p * s_v;

            rng(2027);
            n = 4e6;
            [up, vp] = TestGammaDispersionMarginalization.drawPersonBlock(s_p, s_v, rho, n);
            uo  = s_o  * randn(n,1);
            ut  = s_t  * randn(n,1);
            upo = s_po * randn(n,1);
            upt = s_pt * randn(n,1);
            uot = s_ot * randn(n,1);
            mu   = exp(alpha + up + uo + ut + upo + upt + uot);
            nu_p = exp(beta + vp);
            mc_e_cond = mean(2 .* mu.^2 ./ nu_p);

            vc2 = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, ...
                s_ot, nu_pop, s_v, cov_pv);
            testCase.verifyEqual(vc2.e_cond_var, mc_e_cond, 'RelTol', 1e-2);

            % Population-nu (#1) is biased by the identical factor F.
            F = exp(0.5 * s_v^2 - 2 * cov_pv);
            vc1 = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, ...
                s_ot, nu_pop);
            testCase.verifyEqual(vc1.e_cond_var * F, mc_e_cond, 'RelTol', 1e-2);
        end

    end
end
