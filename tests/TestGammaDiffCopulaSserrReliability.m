classdef TestGammaDiffCopulaSserrReliability < matlab.unittest.TestCase
    %TESTGAMMADIFFCOPULASSERRRELIABILITY Offline proof of the Gamma (scaled-chi-square)
    % SUBJECT-LEVEL CONCURRENT two-event DIFFERENCE-score math (analysis 8,
    % ic_diff_sserrvar, family='gamma', diffrescor=2). No CmdStan.
    %
    % The concurrent gamma case-8 model is the case-7 Gaussian-copula model
    % (psyrat_build_stan_diff_copula_chisq; byte-identical golden proven by
    % TestStanModelSnapshots/ic_gamma_diff_sserr_copula). The SUBJECT-LEVEL residual is
    % a downstream extraction (psyrat_gamma_extract_diff_copula_sserr, a local function
    % of psyrat_computevarcomp): per participant it conditions each event's residual on
    % that person's own r_id effects (as non-concurrent case 8) AND adds the per-person
    % copula residual cross-covariance wp_cov_ss, then feeds the SAME subject-level
    % difference machinery (psyrat_ssrel_diff) the Gaussian case-8 path uses.
    %
    % The central scientific point these tests LOCK (owner decision 2026-07-14): the
    % PER-SUBJECT conditional residual cross-covariance
    %       wp_cov_ss(p) = E_i[ Cov(Y1, Y2 | person p, trial i) ]
    % is 0 at rho_e = 0. This is the OPPOSITE of the GROUP-level residual cross-cov
    % (TestGammaDiffCopulaReliability/testMeanInteractionCrossCovNonzeroAtRhoZero),
    % which is NONZERO at rho_e = 0 because it also carries the population person x
    % trial mean-surface interaction. That interaction is a between-person quantity
    % that vanishes once we condition on ONE person, so at the subject level the
    % concurrent path reduces EXACTLY to the non-concurrent case-8 path at rho_e = 0.
    % These are different estimands, consistent with the subject-level DIAGONAL residual
    % being the pure per-person conditional dispersion (which excludes the trial-mean
    % surface; that lives in the trial facet).
    %
    % The extractor computes wp_cov_ss cheaply via the factorization Y_e = mu_e*W_e
    % (W_e a mean-1 Gamma depending only on nu_e): Cov(Y1,Y2|p,i) = mu1_i*mu2_i*
    % Cov(W1,W2), so only the dimensionless per-person copula correlation r_p =
    % corr(W1,W2) needs Monte Carlo; the trial/mean-surface part is analytic. These
    % tests prove that factorization against a direct brute-force simulation of Y, and
    % the fit->extraction plumbing is validated by the live CmdStan recovery in
    % tests/integration/TestGammaDiffCopulaSserrRecovery.
    %
    % Requires the Statistics Toolbox (gaminv/normcdf); skips gracefully if absent.

    methods (TestMethodSetup)
        function requireStats(testCase)
            testCase.assumeTrue(license('test','statistics_toolbox') == 1 && ...
                exist('gaminv','file') == 2 && exist('normcdf','file') == 2, ...
                'Statistics Toolbox (gaminv/normcdf) unavailable; skipping.');
        end
    end

    methods (Test)

        function testPerSubjectCrossCovZeroAtRhoZero(testCase)
            % DECISIVE (the estimand distinction). Conditioning on ONE person (fixed
            % u_p, nu_p), the residual cross-covariance E_i[Cov(Y1,Y2|p,i)] simulated
            % DIRECTLY from the concurrent observations must be ~0 at rho_e = 0 and
            % increase monotonically with rho_e. This is what makes the SUBJECT-level
            % concurrent difference reduce exactly to the non-concurrent case at
            % rho_e = 0 - unlike the GROUP-level residual cross-cov, which is nonzero
            % at rho_e = 0 (that test lives in TestGammaDiffCopulaReliability).
            alpha = [log(5.5) log(4.0)]; beta = [log(18) log(12)];
            s_i   = [0.35 0.30]; cor_t = 0.4;
            u_p   = [0.30 -0.20];               % one person's log-mean offsets
            ni = 50; nrep = 4000;
            prev = -inf;
            for rho = [0.0 0.3 0.6 0.9]
                rng(7700 + round(100*rho));
                cc = localCondCrossCov(alpha,beta,s_i,cor_t,u_p,rho,ni,nrep);
                if rho == 0
                    testCase.verifyLessThan(abs(cc), 0.05, ...
                        ['At rho_e=0 the PER-SUBJECT conditional residual cross-cov ' ...
                        'must be ~0 (conditionally independent dispersions).']);
                else
                    testCase.verifyGreaterThan(cc, prev, ...
                        'Per-subject conditional residual cross-cov must rise with rho_e.');
                end
                prev = cc;
            end
        end

        function testFactorizationMatchesBruteForceMC(testCase)
            % DECISIVE (the extractor's math). The extractor computes wp_cov_ss via
            %   wp_cov_ss(p) = corr(W1,W2) * sqrt((2/nu1)(2/nu2)) * E_i[mu1_i mu2_i],
            % exploiting Y_e = mu_e*W_e (W_e mean-1, depends on nu_e only). This must
            % equal a DIRECT brute-force simulation of E_i[Cov(Y1,Y2|p,i)] that never
            % uses the factorization. Matching validates Cov(Y1,Y2|cell) =
            % mu1*mu2*Cov(W1,W2) (the cell-independence of the copula coupling) and the
            % analytic Var(W_e) = 2/nu_e.
            alpha = [log(6.0) log(4.5)]; beta = [log(16) log(20)];
            s_i   = [0.30 0.25]; cor_t = 0.5;
            for tc = 1:3
                rng(8800 + tc);
                u_p = 0.5*randn(1,2);
                rho = 0.2 + 0.25*tc;                 % 0.45, 0.70, 0.95
                ni  = 40;
                % shared trial effects so BF and the analytic cross-moment agree
                Lt = chol((s_i'*s_i).*[1 cor_t; cor_t 1],'lower');
                u_i = randn(ni,2)*Lt';
                mu1 = exp(alpha(1)+u_p(1)+u_i(:,1));
                mu2 = exp(alpha(2)+u_p(2)+u_i(:,2));
                nu1 = exp(beta(1)+0);                % nu is subject-constant (u_p on mean only)
                nu2 = exp(beta(2)+0);

                % brute force: per cell, cov of nrep concurrent Y-pairs; average over trials
                nrep = 3000; percell = zeros(ni,1);
                for k = 1:ni
                    Y = localBivGammaCopula(mu1(k),mu2(k),nu1,nu2,rho,nrep);
                    c = cov(Y(:,1),Y(:,2));
                    percell(k) = c(1,2);
                end
                BF = mean(percell);

                % factorization: r_p from a big W-only copula MC, analytic pieces
                M = 4e4;
                z1 = randn(M,1); z2 = rho.*z1 + sqrt(1-rho.^2).*randn(M,1);
                W1 = gaminv(normcdf(z1), 0.5*nu1, 2/nu1);   % mean-1 dispersion factors
                W2 = gaminv(normcdf(z2), 0.5*nu2, 2/nu2);
                rr = corr(W1,W2);
                FC = rr * sqrt((2/nu1)*(2/nu2)) * mean(mu1.*mu2);

                testCase.verifyEqual(FC, BF, 'RelTol', 0.06, ...
                    sprintf('mu*W factorization must match brute-force MC (rho=%.2f).', rho));
            end
        end

        function testConcurrentRaisesSubjectDifferenceReliabilityViaSsrelDiff(testCase)
            % INTEGRATION. Feeding the extractor's outputs (per-event person var-cov,
            % per-subject residual log-SDs, per-subject residual cross-cov wp_cov_ss)
            % into the REAL subject-level difference path (psyrat_ssrel_diff) must:
            %  (a) at wp_cov_ss = 0 reproduce the non-concurrent case-8 reliability
            %      (the exact reduction at rho_e = 0), and
            %  (b) with a positive wp_cov_ss RAISE each subject's difference
            %      reliability (a concurrent design shrinks the difference error).
            nd = 200; NSUB = 3;
            rng(9100);
            sp1 = 4.0 + 0.2*randn(nd,1); sp2 = 5.0 + 0.2*randn(nd,1);
            cov_p = -0.6 + 0.1*randn(nd,1);      % cross-event person covariance
            si1 = 0.15 + 0.02*randn(nd,1); si2 = 0.12 + 0.02*randn(nd,1);
            bp = zeros(nd,2,2);
            bp(:,1,1) = sp1; bp(:,2,2) = sp2; bp(:,1,2) = cov_p; bp(:,2,1) = cov_p;
            bt = zeros(nd,2,2); bt(:,1,1) = si1; bt(:,2,2) = si2;

            wp1 = 1.2 + 0.05*randn(nd,NSUB);     % sigma_pi_e2 per subject, event 1
            wp2 = 1.6 + 0.05*randn(nd,NSUB);     % event 2
            er_var_ss1 = 0.5*log(wp1);
            er_var_ss2 = 0.5*log(wp2);
            b_sigma = [0.5*log(mean(wp1,2)), 0.5*log(mean(wp2,2))];
            % positive per-subject residual cross-cov, PSD (|r|<1): wp_cov = r*sqrt(wp1 wp2)
            rho_ss = 0.5;
            wp_cov_ss = rho_ss .* sqrt(wp1 .* wp2);

            idtable0 = table((1:NSUB)', (1:NSUB)', 'VariableNames', {'id','id2'});
            idtable0.trls  = repmat(20, NSUB, 1);
            idtable0.trls1 = repmat(20, NSUB, 1);
            idtable0.trls2 = repmat(20, NSUB, 1);

            args = {'bp',bp,'bt',bt,'er_var',b_sigma,'idtable',idtable0,...
                'CI',0.95,'est','gen','er_var_ss1',er_var_ss1,'er_var_ss2',er_var_ss2};
            outNonc = psyrat_ssrel_diff(args{:},'wp_cov',zeros(nd,1),...
                'wp_cov_ss',zeros(nd,NSUB));
            outConc = psyrat_ssrel_diff(args{:},'wp_cov',zeros(nd,1),...
                'wp_cov_ss',wp_cov_ss);

            % (a) non-concurrent reproduces the hand reference (wp_cov = 0)
            for r = 1:NSUB
                o1 = idtable0.trls1(r); o2 = idtable0.trls2(r);
                var_pdiff = sp1 + sp2 - 2*cov_p;
                rel_err   = wp1(:,r)./o1 + wp2(:,r)./o2;
                Gdraw = var_pdiff ./ (var_pdiff + rel_err);
                testCase.verifyEqual(outNonc.rel_pt(r), mean(Gdraw), 'RelTol', 1e-10, ...
                    sprintf('Subject %d non-concurrent (wp_cov_ss=0) mismatch.', r));
            end

            % (b) a positive residual cross-cov RAISES subject-level reliability
            for r = 1:NSUB
                testCase.verifyGreaterThan(outConc.rel_pt(r), outNonc.rel_pt(r), ...
                    sprintf('Subject %d: positive wp_cov_ss must raise reliability.', r));
            end
        end

    end
end

function cc = localCondCrossCov(alpha,beta,s_i,cor_t,u_p,rho,ni,nrep)
%Brute-force E_i[Cov(Y1,Y2|person p, trial i)] for one fixed person: draw ni
%correlated trial effects, and at each cell simulate nrep concurrent copula pairs and
%take the per-cell covariance; average over trials. nu is subject-constant (the person
%offset u_p loads the MEAN only here), matching the ID-only nu structure.
Lt = chol((s_i'*s_i).*[1 cor_t; cor_t 1],'lower');
u_i = randn(ni,2)*Lt';
mu1 = exp(alpha(1)+u_p(1)+u_i(:,1));
mu2 = exp(alpha(2)+u_p(2)+u_i(:,2));
nu1 = exp(beta(1)); nu2 = exp(beta(2));
percell = zeros(ni,1);
for k = 1:ni
    Y = localBivGammaCopula(mu1(k),mu2(k),nu1,nu2,rho,nrep);
    c = cov(Y(:,1),Y(:,2));
    percell(k) = c(1,2);
end
cc = mean(percell);
end

function Y = localBivGammaCopula(mu1,mu2,nu1,nu2,rho,n)
%Draw n concurrent scaled-chi-square pairs coupled by a Gaussian copula, exactly as
%psyrat_gamma_extract_diff_copula: z1,z2 Gaussian with correlation rho, mapped to
%uniforms, then to Gamma(shape=0.5*nu, scale=2*mu/nu) marginals.
z1 = randn(n,1); z2 = rho.*z1 + sqrt(max(1-rho.^2,0)).*randn(n,1);
y1 = gaminv(normcdf(z1), 0.5*nu1, (2*mu1)/nu1);
y2 = gaminv(normcdf(z2), 0.5*nu2, (2*mu2)/nu2);
Y = [y1 y2];
end
