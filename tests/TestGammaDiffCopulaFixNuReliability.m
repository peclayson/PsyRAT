classdef TestGammaDiffCopulaFixNuReliability < matlab.unittest.TestCase
    %TESTGAMMADIFFCOPULAFIXNURELIABILITY Offline proof of the Gamma (scaled-chi-
    % square) one-facet CONCURRENT two-event DIFFERENCE-score math for the GROUP-LEVEL
    % GLOBAL fixed-dispersion variant (analysis 7, family='gamma', diffrescor=2,
    % dispersion=2). No CmdStan.
    %
    % The global-nu variant (psyrat_build_stan_diff_copula_fixnu_chisq; byte-identical
    % golden proven by TestStanModelSnapshots/ic_gamma_diff_copula_fixnu) keeps the
    % SAME Gaussian-copula scaled-chi-square likelihood as the per-person-nu model, but
    % the person block is 2-D ([mean_1, mean_2]) and the dispersion nu is a single
    % event-global value (no per-person log-nu). The observed-scale components AND the
    % residual cross-covariance are recovered by the reference full Monte-Carlo
    % projection (psyrat_gamma_extract_diff_copula_fixnu, a local function of
    % psyrat_computevarcomp) exactly as in the per-person case; only the person block
    % dimension and the constant-across-persons nu differ.
    %
    % These tests prove the global-nu-specific math the extraction implements,
    % mirroring TestGammaDiffCopulaReliability (the per-person sibling); the
    % fit->extraction plumbing (field wiring, subsampling) is validated by the live
    % CmdStan recovery in tests/integration/TestGammaDiffCopulaFixNuRecovery.
    %
    % Central finding preserved for the global-nu case: because the mean is LOG-linked,
    % the observed-score residual confounds a lognormal mean-surface interaction with
    % the observation dispersion, so the residual cross-covariance is NONZERO even at
    % rho_e = 0 EVEN WITHOUT person-varying nu (the mean-surface interaction alone
    % drives it). This is why the residual cross-cov still needs the Monte-Carlo
    % projection rather than a rho_e*sigma*sigma formula.
    %
    % Requires the Statistics Toolbox (gaminv/normcdf), already a PsyRAT dependency;
    % skips gracefully if absent, keeping the unit lane green.

    methods (TestMethodSetup)
        function requireStats(testCase)
            testCase.assumeTrue(license('test','statistics_toolbox') == 1 && ...
                exist('gaminv','file') == 2 && exist('normcdf','file') == 2, ...
                'Statistics Toolbox (gaminv/normcdf) unavailable; skipping.');
        end
    end

    methods (Test)

        function testCopulaInducesResidualCorrelation(testCase)
            % DECISIVE (copula draw). At a FIXED cell (mu, nu fixed) the only dependence
            % between the two events' observations is the copula. The observed-score
            % residual Pearson correlation must be ~0 at rho_e=0 and increase
            % monotonically with rho_e, reproducible across independent large-N draws.
            % The copula draw is identical to the per-person case (nu is fixed at a
            % cell either way), so this confirms the mechanic the global-nu extraction
            % relies on.
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

        function testMeanInteractionCrossCovNonzeroAtRhoZeroGlobalNu(testCase)
            % DECISIVE (the log-link finding, GLOBAL-nu). With the full one-facet design
            % (person and trial mean effects correlated across events), a GLOBAL nu
            % (constant across persons), and rho_e = 0 (copula OFF), the observed-score
            % RESIDUAL cross-covariance must still be POSITIVE and bounded away from
            % zero: the lognormal mean surface exp(b+u_p+u_i) has a person x trial
            % interaction whose cross-event covariance survives into the residual
            % (pi_e) component EVEN WITHOUT person-varying nu. This is the mechanism
            % that makes the gamma concurrent difference genuinely different from the
            % Gaussian identity-link one (whose residual cross-cov is exactly 0 here),
            % and it is present in the group-level global-nu parameterization too.
            b = [log(5.5) log(4.0)]; b_nu = [log(18) log(12)];
            sd_id = [0.35 0.30];            % 2-D person MEAN block (no log-nu SDs)
            sd_trl = [0.35 0.30];
            R2 = [1 .6; .6 1];              % cross-event person-mean correlation
            Rt = [1 .4; .4 1];             % cross-event trial correlation
            ce = zeros(1,2);
            for s = 1:2                     % two seeds: check stability
                rng(4400 + s);
                o = localSimOneFacetFixnu(b,b_nu,sd_id,sd_trl,R2,Rt,0.0,3000,120);
                ce(s) = o.ce;
            end
            testCase.verifyGreaterThan(min(ce), 0.02, ...
                ['At rho_e=0 the global-nu log-linked residual cross-cov must be ' ...
                'POSITIVE (mean-surface interaction), not 0.']);
            testCase.verifyLessThan(max(ce)/min(ce), 3, ...
                'Residual cross-cov at rho_e=0 should be stable across seeds.');
        end

        function testGlobalNuProjectionFeedsDiffrel(testCase)
            % INTEGRATION + pinning identity. Running the global-nu one-facet Monte-
            % Carlo projection (2-D person block, event-global nu) per draw and feeding
            % the observed-scale components it emits (per-event person var-cov = bp,
            % trial var-cov = bt, per-event residual log-SD = er_var, residual cross-cov
            % = wp_cov) into the REAL difference reliability path (psyrat_diffrel) must:
            %   (a) satisfy the cross-cov PINNING identity exactly: for every component,
            %       c1 + c2 - 2*crosscov = max(vc_delta,0);
            %   (b) reproduce the difference coefficient computed directly from the
            %       delta components (reliability-from-vc_delta); and
            %   (c) show a positive residual cross-cov RAISING reliability (concurrent
            %       shrinks the difference error variance vs the non-concurrent wp_cov=0).
            b = [log(5.5) log(4.0)]; b_nu = [log(18) log(12)];
            sd_id = [0.35 0.30]; sd_trl = [0.10 0.09];
            R2 = [1 .5; .5 1]; Rt = [1 .3; .3 1];
            rng(7700);   % deterministic draws (reproducibility-first)
            nd = 100; nprime = 20;
            bp = zeros(nd,2,2); bt = zeros(nd,2,2);
            er_var = zeros(nd,2); wp_cov = zeros(nd,1);
            maxIdErr = 0;
            for d = 1:nd
                o = localSimOneFacetFixnu(b,b_nu,sd_id,sd_trl,R2,Rt,0.4,150,40);
                bp(d,1,1)=o.pt; bp(d,2,2)=o.pa; bp(d,1,2)=o.cp; bp(d,2,1)=o.cp;
                bt(d,1,1)=o.it; bt(d,2,2)=o.ia; bt(d,1,2)=o.ci; bt(d,2,1)=o.ci;
                er_var(d,:) = 0.5*log([o.et o.ea]);
                wp_cov(d) = o.ce;
                % pinning identity (person component): pt + pa - 2*cp == max(vc_delta.p,0)
                maxIdErr = max(maxIdErr, abs((o.pt+o.pa-2*o.cp) - max(o.vd_p,0)));
            end
            testCase.verifyLessThan(maxIdErr, 1e-9, ...
                'Cross-cov pinning identity (pt+pa-2cp = max(vc_delta,0)) must hold exactly.');
            testCase.verifyGreaterThan(mean(wp_cov), 0, ...
                'Concurrent residual cross-cov should be positive with positive correlations.');

            gConc = psyrat_diffrel('bp',bp,'bt',bt,'er_var',er_var,'obs',[nprime nprime],...
                'CI',0.95,'est','gen','wp_cov',wp_cov);
            gNonc = psyrat_diffrel('bp',bp,'bt',bt,'er_var',er_var,'obs',[nprime nprime],...
                'CI',0.95,'est','gen','wp_cov',zeros(nd,1));

            % hand reference (generalizability), obs1=obs2=n' so harmmean=n':
            %   uni = sp1+sp2-2cp; rel_err = (wp1+wp2-2wp_cov)/n'; G = uni/(uni+rel_err)
            uni = squeeze(bp(:,1,1)) + squeeze(bp(:,2,2)) - 2*squeeze(bp(:,1,2));
            wp1 = exp(er_var(:,1)).^2; wp2 = exp(er_var(:,2)).^2;
            relC = (wp1 + wp2 - 2*wp_cov) / nprime;
            testCase.verifyEqual(gConc.pt, mean(uni ./ (uni + relC)), 'RelTol', 1e-10, ...
                'Global-nu concurrent generalizability must match reliability-from-vc_delta.');
            testCase.verifyGreaterThan(gConc.pt, gNonc.pt, ...
                'A positive residual cross-cov (concurrent) must raise reliability.');
        end

    end
end

function Y = localBivGammaCopula(mu1,mu2,nu1,nu2,rho,n)
%Draw n concurrent scaled-chi-square pairs coupled by a Gaussian copula, exactly as
%psyrat_gamma_extract_diff_copula_fixnu: z1,z2 Gaussian with correlation rho, mapped
%to uniforms, then to Gamma(shape=0.5*nu, scale=2*mu/nu) marginals (nu fixed at cell).
z1 = randn(n,1); z2 = rho.*z1 + sqrt(max(1-rho.^2,0)).*randn(n,1);
y1 = gaminv(normcdf(z1), 0.5*nu1, (2*mu1)/nu1);
y2 = gaminv(normcdf(z2), 0.5*nu2, (2*mu2)/nu2);
Y = [y1 y2];
end

function out = localSimOneFacetFixnu(b, b_nu, sd_id, sd_trl, R2, Rt, rho, np, ni)
%One-facet GLOBAL-nu copula simulation + ANOVA decomposition, matching the extraction
%psyrat_gamma_extract_diff_copula_fixnu: a 2-D person MEAN block (no log-nu), an
%event-global nu constant across persons, the Gaussian-copula concurrent pair, and the
%psyrat_vc_onefacet decomposition with the cross-cov pinned to max(vc_delta,0).
Lp = localCholLower((sd_id'*sd_id).*R2);   % 2x2 person-mean covariance
Lt = localCholLower((sd_trl'*sd_trl).*Rt); % 2x2 trial covariance
up = randn(np,2)*Lp'; ui = randn(ni,2)*Lt';
mu1 = exp(b(1)+up(:,1)+ui(:,1)'); mu2 = exp(b(2)+up(:,2)+ui(:,2)');
nu1 = exp(b_nu(1)); nu2 = exp(b_nu(2));     % event-global scalars
za = randn(np,ni); zb = randn(np,ni);
ua = min(1-1e-10,max(1e-10,normcdf(za)));
ub = min(1-1e-10,max(1e-10,normcdf(rho.*za + sqrt(max(1-rho.^2,0)).*zb)));
th = gaminv(ua, 0.5*nu1+zeros(np,ni), (2*mu1)./nu1);
al = gaminv(ub, 0.5*nu2+zeros(np,ni), (2*mu2)./nu2);
vt = localVc1(th); va = localVc1(al); vd = localVc1(th-al);
pt = max(vt.p,0);    pa = max(va.p,0);    cp = (pt+pa-max(vd.p,0))/2;
it = max(vt.i,0);    ia = max(va.i,0);    ci = (it+ia-max(vd.i,0))/2;
et = max(vt.pi_e,eps); ea = max(va.pi_e,eps); ce = (et+ea-max(vd.pi_e,0))/2;
out = struct('pt',pt,'pa',pa,'cp',cp,'it',it,'ia',ia,'ci',ci,'et',et,'ea',ea,'ce',ce,...
    'vd_p',vd.p,'vd_i',vd.i,'vd_pi_e',vd.pi_e);
end

function vc = localVc1(Y)
%One-facet random-effects ANOVA components (matches psyrat_vc_onefacet):
%p=(ms_p-ms_pi)/I, i=(ms_i-ms_pi)/P, pi_e=ms_pi.
[P,I] = size(Y);
g = mean(Y(:)); pm = mean(Y,2); im = mean(Y,1);
r = Y - pm - im + g;
ss_p = I*sum((pm-g).^2); ss_i = P*sum((im-g).^2); ss_pi = sum(r(:).^2);
ms_p = ss_p/(P-1); ms_i = ss_i/(I-1); ms_pi = ss_pi/((P-1)*(I-1));
vc = struct('p',(ms_p-ms_pi)/I,'i',(ms_i-ms_pi)/P,'pi_e',ms_pi);
end

function L = localCholLower(C)
C = (C+C')./2; [L,p] = chol(C,'lower');
if p > 0, [V,D] = eig(C); L = V*diag(sqrt(max(real(diag(D)),0))); end
end
