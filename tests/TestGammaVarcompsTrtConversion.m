classdef TestGammaVarcompsTrtConversion < matlab.unittest.TestCase
    %TESTGAMMAVARCOMPSTRTCONVERSION Unit tests for the crossed Persons x
    % Occasions x Trials log-scale -> observed-score-scale variance-component
    % conversion (psyrat_gamma_varcomps_trt) used by the test-retest (analysis 5)
    % scaled-chi-square / Gamma family. Because the mean sub-model uses a log
    % link, the Stan person/occasion/trial/interaction SDs are on the log
    % expected-score scale and must be mapped to observed-score-scale G-theory
    % variance components via lognormal moments BEFORE psyrat_rel_trt consumes
    % them.
    %
    % The checks are INDEPENDENT of the production expressions so a green test
    % means the formula is right, not merely unchanged:
    %   - the degenerate zero-random-effects case has a trivial closed form;
    %   - the crossed helper MUST reduce exactly to the already-validated
    %     one-facet helper (psyrat_gamma_varcomps) when the occasion and all
    %     three interaction SDs are zero -- a strong exact cross-check of the
    %     main-effect, total-variance, and residual arithmetic;
    %   - the six components plus the residual must sum to the total variance
    %     exactly whenever the residual is not clipped;
    %   - a fixed-seed, fully generative Monte Carlo check estimates each
    %     component directly from the covariance-of-cells definition of a
    %     random-effects ANOVA (two cells that share a factor set have covariance
    %     equal to the sum of the variance components indexed within that set),
    %     which exercises the interaction terms the one-facet reduction cannot.

    methods (Test)

        function testDegenerateZeroRandomEffects(testCase)
            % With no random effects (all six log-scale SDs = 0) the mean surface
            % is constant at exp(alpha), so every component and var_mu are 0 and
            % the only observed variance is the scaled-chi-square conditional
            % dispersion E[Var(Y|mu,nu)] = (2/nu)*exp(2*alpha).
            alpha = log(4);
            nu    = 16;
            vc = psyrat_gamma_varcomps_trt(alpha, 0, 0, 0, 0, 0, 0, nu);

            expected_cond = (2/nu) * exp(2*alpha);
            testCase.verifyEqual(vc.sigma_p2,   0,             'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_o2,   0,             'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_t2,   0,             'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_po2,  0,             'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_pt2,  0,             'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_ot2,  0,             'AbsTol', 1e-12);
            testCase.verifyEqual(vc.var_mu,     0,             'AbsTol', 1e-12);
            testCase.verifyEqual(vc.e_cond_var, expected_cond, 'AbsTol', 1e-12);
            testCase.verifyEqual(vc.total_var,  expected_cond, 'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_res2, expected_cond, 'AbsTol', 1e-12);
            testCase.verifyEqual(vc.mu_bar,     exp(alpha),    'AbsTol', 1e-12);
        end

        function testPartialReductionToOneFacet(testCase)
            % With the occasion main effect and ALL THREE log-scale interaction
            % SDs set to zero, the design has only person and trial log-scale
            % effects -- but it does NOT collapse to the one-facet helper term for
            % term. Under the log link the mean is multiplicative
            % (mu = exp(alpha)*exp(u_p)*exp(u_t)), so persons and trials INTERACT
            % on the observed score even with no log-scale interaction. The
            % crossed helper therefore attributes that induced interaction to an
            % explicit observed-score person x trial component
            %   sigma_pt2 = B*(exp(v_p)-1)*(exp(v_t)-1),   B = exp(2*alpha + v_p + v_t),
            % which the one-facet helper instead folds into its residual (a
            % single-occasion design cannot separate P x T from error). The
            % checks below encode that split, which is the scientifically
            % consequential difference: sigma_pt2 loads the universe score for the
            % stability coefficient, whereas residual divides by trials*occasions.
            alpha = [log(5.5);  log(6.0);  log(4.2)];
            s_p   = [0.40;      0.25;      0.55];
            s_t   = [0.10;      0.18;      0.05];
            nu    = [18;        22;        14];

            onef = psyrat_gamma_varcomps(alpha, s_p, s_t, nu);
            vc   = psyrat_gamma_varcomps_trt(alpha, s_p, 0, s_t, 0, 0, 0, nu);

            % Marginal quantities and the two present main effects match exactly.
            testCase.verifyEqual(vc.mu_bar,     onef.mu_bar,     'RelTol', 1e-12);
            testCase.verifyEqual(vc.sigma_p2,   onef.sigma_p2,   'RelTol', 1e-12);
            testCase.verifyEqual(vc.sigma_t2,   onef.sigma_i2,   'RelTol', 1e-12);
            testCase.verifyEqual(vc.var_mu,     onef.var_mu,     'RelTol', 1e-12);
            testCase.verifyEqual(vc.e_cond_var, onef.e_cond_var, 'RelTol', 1e-12);
            testCase.verifyEqual(vc.total_var,  onef.total_var,  'RelTol', 1e-12);

            % Occasion and the occasion-bearing interactions are exactly zero.
            testCase.verifyEqual(vc.sigma_o2,  zeros(3,1), 'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_po2, zeros(3,1), 'AbsTol', 1e-12);
            testCase.verifyEqual(vc.sigma_ot2, zeros(3,1), 'AbsTol', 1e-12);

            % The induced person x trial component, independently written.
            vp = s_p.^2; vt = s_t.^2;
            B  = exp(2.*alpha + vp + vt);
            expected_pt = B .* (exp(vp) - 1) .* (exp(vt) - 1);
            testCase.verifyEqual(vc.sigma_pt2, expected_pt, 'RelTol', 1e-12);

            % The crossed residual plus that induced interaction reconstruct the
            % one-facet residual, which confounds P x T with dispersion.
            testCase.verifyEqual(vc.sigma_res2 + vc.sigma_pt2, ...
                onef.sigma_pi_e2, 'RelTol', 1e-12);
        end

        function testTotalVarianceAdditiveIdentity(testCase)
            % The six components plus the residual must reconstruct the total
            % observed variance exactly (residual is defined by subtraction and
            % is not clipped for these non-degenerate inputs).
            alpha = log(5);
            s_p = 0.35; s_o = 0.25; s_t = 0.20;
            s_po = 0.18; s_pt = 0.15; s_ot = 0.12; nu = 20;
            vc = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, nu);

            recon = vc.sigma_p2 + vc.sigma_o2 + vc.sigma_t2 ...
                + vc.sigma_po2 + vc.sigma_pt2 + vc.sigma_ot2 + vc.sigma_res2;
            testCase.verifyEqual(recon, vc.total_var, 'RelTol', 1e-12);

            % var_mu and e_cond_var must also add to the total (law of total var).
            testCase.verifyEqual(vc.var_mu + vc.e_cond_var, vc.total_var, ...
                'RelTol', 1e-12);
        end

        function testMonteCarloComponents(testCase)
            % Fixed-seed, fully generative Monte Carlo validation using only
            % base-MATLAB randn (no Statistics Toolbox). Each observed-score
            % component is estimated directly from the random-effects covariance
            % identity: two cells that share exactly a factor set A have
            % covariance equal to the sum of the variance components indexed
            % within A. This is independent of the analytic derivation and would
            % catch a wrong inclusion-exclusion sign, not just a changed formula.
            alpha = log(5);
            s_p = 0.35; s_o = 0.25; s_t = 0.20;
            s_po = 0.18; s_pt = 0.15; s_ot = 0.12; nu = 20;
            vc = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, nu);

            rng(2026);
            n = 4e6;

            % Two independent full sets of the six effects. A "cell" is the
            % exponentiated sum of one effect drawn from each factor; sharing a
            % factor set A between two cells means those effects are taken from
            % the SAME set.
            up1 = s_p*randn(n,1);  uo1 = s_o*randn(n,1);  ut1 = s_t*randn(n,1);
            po1 = s_po*randn(n,1); pt1 = s_pt*randn(n,1);  ot1 = s_ot*randn(n,1);
            up2 = s_p*randn(n,1);  uo2 = s_o*randn(n,1);  ut2 = s_t*randn(n,1);
            po2 = s_po*randn(n,1); pt2 = s_pt*randn(n,1);  ot2 = s_ot*randn(n,1);

            % Reference cell uses set 1 throughout; every "share S" partner reuses
            % set-1 effects for the factors in S and set-2 effects otherwise, so
            % the reference cell is common across all estimators (variance
            % reduction for the interaction differences).
            mu1 = exp(alpha + up1 + uo1 + ut1 + po1 + pt1 + ot1);

            % Population covariance estimator (divide by n, matching var(.,1)).
            pcov = @(x,y) mean(x.*y) - mean(x).*mean(y);

            % Full-surface mean and variance.
            testCase.verifyEqual(mean(mu1),  vc.mu_bar, 'RelTol', 5e-3);
            testCase.verifyEqual(var(mu1,1), vc.var_mu, 'RelTol', 2e-2);

            % Main effects: share exactly one factor.
            b_p = exp(alpha + up1 + uo2 + ut2 + po2 + pt2 + ot2); % share {p}
            b_o = exp(alpha + up2 + uo1 + ut2 + po2 + pt2 + ot2); % share {o}
            b_t = exp(alpha + up2 + uo2 + ut1 + po2 + pt2 + ot2); % share {t}
            cov_p = pcov(mu1, b_p);
            cov_o = pcov(mu1, b_o);
            cov_t = pcov(mu1, b_t);
            testCase.verifyEqual(cov_p, vc.sigma_p2, 'RelTol', 2e-2);
            testCase.verifyEqual(cov_o, vc.sigma_o2, 'RelTol', 2e-2);
            testCase.verifyEqual(cov_t, vc.sigma_t2, 'RelTol', 2e-2);

            % Two-way interactions: share both parent factors AND their
            % interaction. Cov(share {x,y,xy}) = sigma_x + sigma_y + sigma_xy, so
            % the interaction component is recovered by inclusion-exclusion using
            % the SAME common-random-number covariances above.
            b_po = exp(alpha + up1 + uo1 + ut2 + po1 + pt2 + ot2); % share {p,o,po}
            b_pt = exp(alpha + up1 + uo2 + ut1 + po2 + pt1 + ot2); % share {p,t,pt}
            b_ot = exp(alpha + up2 + uo1 + ut1 + po2 + pt2 + ot1); % share {o,t,ot}
            po_hat = pcov(mu1, b_po) - cov_p - cov_o;
            pt_hat = pcov(mu1, b_pt) - cov_p - cov_t;
            ot_hat = pcov(mu1, b_ot) - cov_o - cov_t;
            testCase.verifyEqual(po_hat, vc.sigma_po2, 'RelTol', 4e-2);
            testCase.verifyEqual(pt_hat, vc.sigma_pt2, 'RelTol', 4e-2);
            testCase.verifyEqual(ot_hat, vc.sigma_ot2, 'RelTol', 4e-2);
        end

    end
end
