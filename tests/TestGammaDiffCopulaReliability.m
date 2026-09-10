classdef TestGammaDiffCopulaReliability < matlab.unittest.TestCase
    %TESTGAMMADIFFCOPULARELIABILITY Offline proof of the Gamma (scaled-chi-square)
    % one-facet CONCURRENT two-event DIFFERENCE-score math (analysis 7,
    % family='gamma', diffrescor=2). No CmdStan.
    %
    % The concurrent gamma difference couples the two events' scaled-chi-square
    % residuals with a Gaussian copula correlation rho_e
    % (psyrat_build_stan_diff_copula_chisq; byte-identical golden proven by
    % TestStanModelSnapshots/ic_gamma_diff_copula). The observed-scale variance
    % components AND the residual cross-covariance are recovered by the reference's
    % FULL Monte-Carlo projection (psyrat_gamma_extract_diff_copula, a local function
    % of psyrat_computevarcomp): per posterior draw it simulates the copula design
    % and ANOVA-decomposes theta, alpha, delta.
    %
    % These tests prove the gamma-specific math the extraction implements, mirroring
    % TestGammaSserrConversion / TestGammaDiffSserrReliability; the fit->extraction
    % plumbing (field wiring, subsampling) is validated by the live CmdStan recovery
    % in tests/integration/TestGammaDiffCopulaRecovery.
    %
    % They also LOCK the central scientific finding that distinguishes the gamma
    % concurrent difference from the Gaussian one: because the mean is LOG-linked,
    % the observed-score residual confounds a lognormal mean-surface interaction with
    % the observation dispersion, so the residual cross-covariance is NONZERO even at
    % rho_e = 0 (the Gaussian identity-link residual cross-cov is exactly 0 there).
    %
    % Requires the Statistics Toolbox (gaminv/normcdf), already a PsyRAT dependency
    % (psyrat_diffrel uses quantile); the tests assume it and skip gracefully if
    % absent, keeping the unit lane green.

    methods (TestMethodSetup)
        function requireStats(testCase)
            testCase.assumeTrue(license('test','statistics_toolbox') == 1 && ...
                exist('gaminv','file') == 2 && exist('normcdf','file') == 2, ...
                'Statistics Toolbox (gaminv/normcdf) unavailable; skipping.');
        end
    end

    methods (Test)

        function testCopulaInducesResidualCorrelation(testCase)
            % DECISIVE (copula draw). At a FIXED cell (mu, nu fixed, so no design
            % variance), the only dependence between the two events' observations is
            % the copula. The observed-score residual Pearson correlation must be ~0
            % at rho_e = 0 and increase monotonically with rho_e, and be reproducible
            % across independent large-N draws. This validates the copula sampling
            % mechanic (z2 = rho*z1 + sqrt(1-rho^2)*z_indep -> gamma marginals) the
            % extraction relies on.
            mu1 = 5.5; mu2 = 4.0; nu1 = 18; nu2 = 12;
            N = 4e5;
            prev = 0;
            for rho = [0.0 0.3 0.6 0.9]
                rng(3300 + round(100*rho));
                yA = localBivGammaCopula(mu1,mu2,nu1,nu2,rho,N);
                rng(9900 + round(100*rho));      % independent replicate
                yB = localBivGammaCopula(mu1,mu2,nu1,nu2,rho,N);
                rA = corr(yA(:,1),yA(:,2));
                rB = corr(yB(:,1),yB(:,2));
                testCase.verifyEqual(rA, rB, 'AbsTol', 5e-3, ...
                    sprintf('Copula correlation not reproducible at rho=%.1f.', rho));
                if rho == 0
                    testCase.verifyLessThan(abs(rA), 5e-3, ...
                        'At rho_e=0 the fixed-cell residual correlation must be ~0.');
                else
                    testCase.verifyGreaterThan(rA, prev, ...
                        'Observed residual correlation must increase with rho_e.');
                end
                prev = rA;
            end
        end

        function testMeanInteractionCrossCovNonzeroAtRhoZero(testCase)
            % DECISIVE (the log-link finding). With the full one-facet design (person
            % and trial mean effects correlated across events) and rho_e = 0 (copula
            % OFF), the observed-score RESIDUAL cross-covariance must be POSITIVE and
            % bounded away from zero: the lognormal mean surface exp(alpha+u_p+u_i) has
            % a person x trial interaction whose cross-event covariance survives into
            % the residual (pi_e) component. The Gaussian identity-link residual
            % cross-cov would be exactly 0 here, so this is what makes the gamma
            % concurrent difference genuinely different, and why the residual cross-cov
            % needs the Monte-Carlo projection rather than a rho_e*sigma*sigma formula.
            b = [log(5.5) log(4.0)]; b_nu = [log(18) log(12)];
            sd_id = [0.35 0.30 0.50 0.40]; sd_trl = [0.35 0.30];
            R = eye(4); R(1,2)=.6; R(2,1)=.6; R(3,4)=.5; R(4,3)=.5;
            R(1,3)=-.2; R(3,1)=-.2; R(2,4)=-.2; R(4,2)=-.2;
            Rt = [1 .4; .4 1];
            piecov = zeros(1,2);
            for s = 1:2                              % two seeds: check stability
                rng(4400 + s);
                v = localSimResidualCrossCov(b,b_nu,sd_id,sd_trl,R,Rt,0.0,3000,120);
                piecov(s) = v;
            end
            testCase.verifyGreaterThan(min(piecov), 0.02, ...
                ['At rho_e=0 the log-linked residual cross-cov must be POSITIVE ' ...
                '(mean-surface interaction), not 0.']);
            testCase.verifyEqual(piecov(1), piecov(2), 'RelTol', 0.35, ...
                'Residual cross-cov at rho_e=0 should be stable across seeds.');
        end

        function testConcurrentReducesDifferenceErrorViaDiffrel(testCase)
            % INTEGRATION. Feeding the observed-scale components the extractor emits
            % (per-event person var-cov = bp, trial var-cov = bt, per-event residual
            % log-SD = er_var, residual cross-cov = wp_cov) into the REAL difference
            % reliability path (psyrat_diffrel) must (a) reproduce the difference
            % coefficient computed directly from the delta components, and (b) show a
            % positive residual cross-cov RAISING reliability (a concurrent design
            % shrinks the difference error variance vs the non-concurrent wp_cov=0).
            rng(5500);   % deterministic draws (reproducibility-first; the other
                         % methods in this file seed their randn the same way)
            nd = 200;
            % per-event observed-scale components (draw vectors so the CI path runs)
            bp = zeros(nd,2,2); bt = zeros(nd,2,2);
            sp1 = 5.0 + 0.2*randn(nd,1); sp2 = 4.0 + 0.2*randn(nd,1);
            cp  = 1.8 + 0.1*randn(nd,1);              % positive person cross-cov
            bp(:,1,1)=sp1; bp(:,2,2)=sp2; bp(:,1,2)=cp; bp(:,2,1)=cp;
            si1 = 0.9 + 0.05*randn(nd,1); si2 = 0.8 + 0.05*randn(nd,1);
            ci  = 0.3 + 0.03*randn(nd,1);             % positive trial cross-cov
            bt(:,1,1)=si1; bt(:,2,2)=si2; bt(:,1,2)=ci; bt(:,2,1)=ci;
            wp1 = 7.0 + 0.2*randn(nd,1); wp2 = 4.5 + 0.2*randn(nd,1);
            er_var = [0.5*log(wp1), 0.5*log(wp2)];    % residual log-SD-equivalents
            wp_cov = 2.3 + 0.1*randn(nd,1);           % positive residual cross-cov
            obs = [20 20];

            gConc = psyrat_diffrel('bp',bp,'bt',bt,'er_var',er_var,'obs',obs,...
                'CI',0.95,'est','gen','wp_cov',wp_cov);
            gNonc = psyrat_diffrel('bp',bp,'bt',bt,'er_var',er_var,'obs',obs,...
                'CI',0.95,'est','gen','wp_cov',zeros(nd,1));

            % hand reference (generalizability), obs1=obs2=n' so harmmean=n':
            %   uni = sp1+sp2-2cp; rel_err = (wp1+wp2-2wp_cov)/n'; G = uni/(uni+rel_err)
            uni = sp1 + sp2 - 2*cp;
            relC = (wp1 + wp2 - 2*wp_cov) / obs(1);
            relN = (wp1 + wp2) / obs(1);
            testCase.verifyEqual(gConc.pt, mean(uni ./ (uni + relC)), 'RelTol', 1e-10, ...
                'Concurrent generalizability must match reliability-from-vc_delta.');
            testCase.verifyEqual(gNonc.pt, mean(uni ./ (uni + relN)), 'RelTol', 1e-10, ...
                'Non-concurrent generalizability must match the wp_cov=0 reference.');
            testCase.verifyGreaterThan(gConc.pt, gNonc.pt, ...
                'A positive residual cross-cov (concurrent) must raise reliability.');
        end

        function testTwoFacetConcurrentReducesErrorViaDiffrelTrt(testCase)
            % INTEGRATION (two-facet). The two-facet copula extraction fills every
            % crossed component varcov (person, trial, occasion, and the three two-
            % way interactions) WITH its cross-cov, plus the residual cross-cov
            % (er_cov). Feeding these into psyrat_diffrel_trt must reproduce the
            % Table-6 CES coefficient computed directly from the difference
            % components, and the concurrent design (nonzero trial-indexed +
            % residual cross-covs) must RAISE reliability vs the non-concurrent form
            % (those cross-covs zeroed). This proves the two-facet integration
            % contract the case-10 dispatch relies on.
            rng(6600);   % deterministic draws (reproducibility-first; matches the
                         % seeded convention of the other methods in this file)
            nd = 150; obs = [20 20]; nocc = [2 2]; n = obs(1); no = nocc(1);
            mk = @(v1,v2,cov) local2x2(nd,v1,v2,cov);
            bp  = mk(6.0,5.0,1.8);    % person (cross-cov present either way)
            bpi = mk(5.0,4.5,0.7);    % person x trial (concurrent cross-cov)
            bpo = mk(2.6,2.4,0.5);    % person x occasion
            bt  = mk(4.0,3.6,0.9);    % trial (concurrent cross-cov)
            bo  = mk(2.8,2.6,0.8);    % occasion
            boi = mk(2.0,1.9,0.25);   % trial x occasion (concurrent cross-cov)
            wp1 = 10.5 + 0.2*randn(nd,1); wp2 = 9.0 + 0.2*randn(nd,1);
            er_var = [0.5*log(wp1), 0.5*log(wp2)];
            er_cov = 2.7 + 0.1*randn(nd,1);      % residual cross-cov (copula)

            % non-concurrent counterparts: trial-indexed + residual cross-covs zeroed
            bpiN = mk(5.0,4.5,0);  btN = mk(4.0,3.6,0);  boiN = mk(2.0,1.9,0);

            args = {'bp',bp,'bpo',bpo,'bo',bo,'er_var',er_var,'obs',obs,...
                'nocc',nocc,'reltype',3,'est','gen','CI',0.95};
            gConc = psyrat_diffrel_trt('bpi',bpi,'bt',bt,'boi',boi,'er_cov',er_cov,args{:});
            gNonc = psyrat_diffrel_trt('bpi',bpiN,'bt',btN,'boi',boiN,'er_cov',zeros(nd,1),args{:});

            % hand reference (CES generalizability, obs/nocc equal so harmmeans = n/no)
            uni = bp(:,1,1)+bp(:,2,2)-2*bp(:,1,2);
            bpiD = (bpi(:,1,1)+bpi(:,2,2)-2*bpi(:,1,2))/n;
            bpoD = (bpo(:,1,1)+bpo(:,2,2)-2*bpo(:,1,2))/no;
            bpoiD = (wp1+wp2-2*er_cov)/(n*no);
            testCase.verifyEqual(gConc.pt, mean(uni./(uni+bpiD+bpoD+bpoiD)), 'RelTol',1e-10, ...
                'Two-facet concurrent CES generalizability must match reliability-from-vc_delta.');
            testCase.verifyGreaterThan(gConc.pt, gNonc.pt, ...
                'Concurrent (nonzero trial-indexed + residual cross-covs) must raise two-facet reliability.');
        end

    end
end

function vc = local2x2(nd, v1, v2, cov)
vc = zeros(nd,2,2);
vc(:,1,1) = v1; vc(:,2,2) = v2; vc(:,1,2) = cov; vc(:,2,1) = cov;
end

function Y = localBivGammaCopula(mu1,mu2,nu1,nu2,rho,n)
%Draw n concurrent scaled-chi-square pairs coupled by a Gaussian copula, exactly
%as psyrat_gamma_extract_diff_copula: z1,z2 Gaussian with correlation rho, mapped
%to uniforms, then to Gamma(shape=0.5*nu, scale=2*mu/nu) marginals.
z1 = randn(n,1); z2 = rho.*z1 + sqrt(max(1-rho.^2,0)).*randn(n,1);
y1 = gaminv(normcdf(z1), 0.5*nu1, (2*mu1)/nu1);
y2 = gaminv(normcdf(z2), 0.5*nu2, (2*mu2)/nu2);
Y = [y1 y2];
end

function piecov = localSimResidualCrossCov(b,b_nu,sd_id,sd_trl,R,Rt,rho,np,ni)
%One-facet copula simulation + ANOVA residual cross-covariance, matching the
%extraction's psyrat_vc_onefacet path.
Lp = localCholLower((sd_id'*sd_id).*R);
Lt = localCholLower((sd_trl'*sd_trl).*Rt);
up = randn(np,4)*Lp'; ui = randn(ni,2)*Lt';
mu1 = exp(b(1)+up(:,1)+ui(:,1)'); mu2 = exp(b(2)+up(:,2)+ui(:,2)');
nu1 = exp(b_nu(1)+up(:,3));       nu2 = exp(b_nu(2)+up(:,4));
za = randn(np,ni); zb = randn(np,ni);
ua = normcdf(za); ub = normcdf(rho.*za + sqrt(max(1-rho.^2,0)).*zb);
th = gaminv(ua, (0.5*nu1)+zeros(1,ni), (2*mu1)./nu1);
al = gaminv(ub, (0.5*nu2)+zeros(1,ni), (2*mu2)./nu2);
et = localVcResid(th); ea = localVcResid(al); ed = localVcResid(th-al);
piecov = (et + ea - ed) / 2;
end

function e = localVcResid(Y)
[P,I] = size(Y);
g = mean(Y(:)); pm = mean(Y,2); im = mean(Y,1);
r = Y - pm - im + g;
e = sum(r(:).^2) / ((P-1)*(I-1));
end

function L = localCholLower(C)
C = (C+C')./2; [L,p] = chol(C,'lower');
if p > 0, [V,D] = eig(C); L = V*diag(sqrt(max(real(diag(D)),0))); end
end
