classdef TestGammaDynrelReliability < PsyRATTestBase
    % Offline proof of the Gamma (scaled-chi-square) dynamic/conditional
    % reliability SURFACE layer (psyrat_rel_dynrel_gamma, one-facet analysis 11;
    % psyrat_rel_dynrel_trt_gamma, trial+occasion two-facet analysis 14). No
    % CmdStan: the surfaces are deterministic functions of the posterior draws, so
    % synthetic draws exercise the full log->observed conversion + coefficient math.
    %
    % The DECISIVE check is the zero-slope reduction: with the dimension slopes
    % b = b_nu = 0 the mean mu(z) and dispersion nu(z) do not depend on z, so the
    % gamma dynrel surface must be (a) CONSTANT across z and (b) EQUAL to the static
    % gamma reliability computed from the same building blocks the static family
    % uses (psyrat_gamma_varcomps[_trt] + psyrat_rel_sing / psyrat_rel_trt). This
    % ties the new dimensional path to the already-validated static gamma family.
    %
    % All draws use base-MATLAB randn only (no Statistics Toolbox).

    methods (Static, Access = private)
        function chol = makeChol(rho, n)
            % 2x2 cholesky_factor_corr draws with R(2,1) = rho: L = [1 0; rho s].
            chol = zeros(n,2,2);
            chol(:,1,1) = 1;
            chol(:,2,1) = rho;
            chol(:,2,2) = sqrt(1 - rho.^2);
        end
    end

    methods (Test)

        function testOneFacetZeroSlopeReducesToStaticGamma(testCase)
            % DECISIVE (one-facet). b = b_nu = 0 -> surface constant across z and
            % equal to the static one-facet gamma reliability.
            rng(1101);
            n = 4000;
            alpha0 = log(5.5) + 0.1*randn(n,1);
            lognu0 = log(18)  + 0.1*randn(n,1);
            s_p    = abs(0.35 + 0.05*randn(n,1));
            s_v    = abs(0.30 + 0.05*randn(n,1));
            sig_i  = abs(0.35 + 0.05*randn(n,1));
            chol   = testCase.makeChol(-0.30 + 0.05*randn(n,1), n);
            rho_pv = chol(:,2,1);
            cov_pv = rho_pv .* s_p .* s_v;
            zgrid  = [-1 -0.5 0 0.5 1];
            obs = 8; ci = 0.95;

            b0 = zeros(n,1);
            rel = psyrat_rel_dynrel_gamma('alpha0',alpha0,'b',b0,'lognu0',lognu0,...
                'b_nu',b0,'sig_p',s_p,'sig_v',s_v,'rho_pv',rho_pv,'sig_i',sig_i,...
                'ndim',1,'z1',zgrid,'obs',obs,'CI',ci);

            % static reference from the same building blocks
            vc = psyrat_gamma_varcomps(alpha0, s_p, sig_i, exp(lognu0), s_v, cov_pv);
            bp = sqrt(vc.sigma_p2); iz = sqrt(vc.sigma_i2);
            wp = sqrt(vc.sigma_pi_e2 + vc.sigma_i2);
            [~,Gref] = psyrat_rel_sing('gcoeff',2,'metric','global','bp',bp,'wp',wp,'i',iz,'obs',obs,'CI',ci);
            [~,Dref] = psyrat_rel_sing('gcoeff',1,'metric','global','bp',bp,'wp',wp,'i',iz,'obs',obs,'CI',ci);

            testCase.verifyEqual(rel.G.pt(:), repmat(rel.G.pt(1),numel(zgrid),1), ...
                'AbsTol', 1e-12, 'One-facet G(z) not constant at zero slope.');
            testCase.verifyEqual(rel.G.pt(1), Gref, 'AbsTol', 1e-12, ...
                'One-facet G(z) does not reduce to the static gamma G.');
            testCase.verifyEqual(rel.D.pt(1), Dref, 'AbsTol', 1e-12, ...
                'One-facet D(z) does not reduce to the static gamma D.');
        end

        function testTwoFacetZeroSlopeReducesToStaticGamma(testCase)
            % DECISIVE (two-facet). b = b_nu = 0 -> surface constant across z and
            % equal to the static crossed test-retest gamma reliability.
            rng(1114);
            n = 4000;
            alpha0 = log(5.5) + 0.1*randn(n,1);
            lognu0 = log(18)  + 0.1*randn(n,1);
            s_p    = abs(0.35 + 0.05*randn(n,1));
            s_v    = abs(0.30 + 0.05*randn(n,1));
            s_o    = abs(0.20 + 0.05*randn(n,1));
            s_t    = abs(0.35 + 0.05*randn(n,1));
            s_pt   = abs(0.20 + 0.05*randn(n,1));  % trial x person   (sig_trlxid)
            s_po   = abs(0.20 + 0.05*randn(n,1));  % occasion x person(sig_occxid)
            s_ot   = abs(0.10 + 0.03*randn(n,1));  % trial x occasion (sig_trlxocc)
            chol   = testCase.makeChol(-0.30 + 0.05*randn(n,1), n);
            rho_pv = chol(:,2,1);
            cov_pv = rho_pv .* s_p .* s_v;
            zgrid  = [-1 0 1]; obs = 8; ci = 0.95;

            b0 = zeros(n,1);
            rel = psyrat_rel_dynrel_trt_gamma('alpha0',alpha0,'b',b0,'lognu0',lognu0,...
                'b_nu',b0,'sig_p',s_p,'sig_v',s_v,'rho_pv',rho_pv,...
                'sig_occ',s_o,'sig_trl',s_t,'sig_trlxid',s_pt,'sig_occxid',s_po,...
                'sig_trlxocc',s_ot,'ndim',1,'z1',zgrid,'obs',obs,'nocc',1,...
                'reltype',3,'CI',ci);

            % static reference: varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, nu, ...)
            vct = psyrat_gamma_varcomps_trt(alpha0, s_p, s_o, s_t, s_po, s_pt, s_ot, ...
                exp(lognu0), s_v, cov_pv);
            fa = {'bo',sqrt(vct.sigma_o2),'bt',sqrt(vct.sigma_t2),...
                  'txp',sqrt(vct.sigma_pt2),'oxp',sqrt(vct.sigma_po2),...
                  'txo',sqrt(vct.sigma_ot2)};
            [~,GrefT] = psyrat_rel_trt('gcoeff',2,'reltype',3,'bp',sqrt(vct.sigma_p2),...
                fa{:},'err',sqrt(vct.sigma_res2),'obs',obs,'nocc',1,'CI',ci);
            [~,DrefT] = psyrat_rel_trt('gcoeff',1,'reltype',3,'bp',sqrt(vct.sigma_p2),...
                fa{:},'err',sqrt(vct.sigma_res2),'obs',obs,'nocc',1,'CI',ci);

            testCase.verifyEqual(rel.G.pt(:), repmat(rel.G.pt(1),numel(zgrid),1), ...
                'AbsTol', 1e-12, 'Two-facet G(z) not constant at zero slope.');
            testCase.verifyEqual(rel.G.pt(1), GrefT, 'AbsTol', 1e-12, ...
                'Two-facet G(z) does not reduce to the static gamma G.');
            testCase.verifyEqual(rel.D.pt(1), DrefT, 'AbsTol', 1e-12, ...
                'Two-facet D(z) does not reduce to the static gamma D.');
        end

        function testDispersionSlopeMovesReliability(testCase)
            % A positive log-nu dimension slope raises nu(z) with z -> lower
            % dispersion -> higher reliability. G(z) must increase strictly and
            % stay within (0,1).
            rng(1120);
            n = 3000;
            alpha0 = log(5.5) + 0.1*randn(n,1);
            lognu0 = log(18)  + 0.1*randn(n,1);
            s_p    = abs(0.35 + 0.05*randn(n,1));
            s_v    = abs(0.30 + 0.05*randn(n,1));
            sig_i  = abs(0.35 + 0.05*randn(n,1));
            chol   = testCase.makeChol(-0.30 + 0.05*randn(n,1), n);
            rho_pv = chol(:,2,1);
            zgrid  = [-1 -0.5 0 0.5 1];

            b0 = zeros(n,1); bN = 0.4*ones(n,1);
            rel = psyrat_rel_dynrel_gamma('alpha0',alpha0,'b',b0,'lognu0',lognu0,...
                'b_nu',bN,'sig_p',s_p,'sig_v',s_v,'rho_pv',rho_pv,'sig_i',sig_i,...
                'ndim',1,'z1',zgrid,'obs',8,'CI',0.95);
            G = rel.G.pt(:);
            testCase.verifyTrue(all(diff(G) > 0), ...
                'G(z) should increase with z for a positive log-nu slope.');
            testCase.verifyTrue(all(G > 0 & G < 1), 'G(z) out of (0,1).');
        end

        function testTwoDimensionSurfaceShape(testCase)
            % Two-dimension (KDIM=3) surface returns an nz1 x nz2 grid.
            rng(1130);
            n = 1500;
            alpha0 = log(5.5) + 0.1*randn(n,1);
            lognu0 = log(18)  + 0.1*randn(n,1);
            s_p    = abs(0.35 + 0.05*randn(n,1));
            s_v    = abs(0.30 + 0.05*randn(n,1));
            sig_i  = abs(0.35 + 0.05*randn(n,1));
            rho_pv = -0.30*ones(n,1);
            rel = psyrat_rel_dynrel_gamma('alpha0',alpha0,'b',0.2*ones(n,3),...
                'lognu0',lognu0,'b_nu',0.1*ones(n,3),'sig_p',s_p,'sig_v',s_v,...
                'rho_pv',rho_pv,'sig_i',sig_i,'ndim',2,'z1',[-1 0 1],'z2',[-1 0 1],...
                'obs',8,'CI',0.95);
            testCase.verifyEqual(size(rel.G.pt), [3 3], ...
                'Two-dimension G surface has the wrong size.');
        end

    end
end
