classdef TestGammaDiffDynrelReliability < PsyRATTestBase
    % Offline proof of the Gamma (scaled-chi-square) DYNAMIC DIFFERENCE-score
    % reliability surface (psyrat_rel_diffdynrel_gamma, analysis 12,
    % family='gamma', group-level, non-concurrent). No CmdStan: the surface is a
    % deterministic function of the posterior draws, so synthetic draws exercise
    % the full log->observed conversion + difference-coefficient math.
    %
    % The DECISIVE check is the zero-slope reduction: with b_dim = b_nu_dim = 0
    % neither mu nor nu depends on z, so the surface must be (a) CONSTANT across z
    % and (b) EQUAL to the STATIC gamma difference score computed from the same
    % building blocks the static case-7 path uses (psyrat_gamma_varcomps +
    % psyrat_gamma_crosscov -> psyrat_diffrel). That ties the new dimensional
    % difference path to the already-validated static gamma difference family.
    %
    % Also pinned here:
    %   - the two-dimension (KDIM=6) event-crossed slope layout,
    %   - per-event n' (difference designs are routinely unbalanced),
    %   - the gamma-vs-Gaussian structural point that a MEAN slope alone moves the
    %     observed-scale variance components (impossible under an identity link),
    %   - the psyrat_dynrel_summary family dispatch, which must NOT route a gamma
    %     case-12 result into the Gaussian difference branch.
    %
    % All draws use base-MATLAB randn only (no Statistics Toolbox).

    methods (Static, Access = private)

        function C = makeCorId(n, rho_p, rho_pv1, rho_pv2)
            % draws x 4 x 4 person correlation array for the block
            % [mean_1, mean_2, lognu_1, lognu_2]. The converters read only
            % (1,2) cross-event mean, (1,3) and (2,4) within-event mean<->log-nu,
            % so the remaining entries are filled symmetrically for realism but
            % are never consumed; positive-definiteness is therefore not exercised
            % and not asserted here.
            C = zeros(n,4,4);
            for k = 1:4
                C(:,k,k) = 1;
            end
            C(:,1,2) = rho_p;   C(:,2,1) = rho_p;
            C(:,1,3) = rho_pv1; C(:,3,1) = rho_pv1;
            C(:,2,4) = rho_pv2; C(:,4,2) = rho_pv2;
        end

        function p = makeParams(n, seed)
            % A common, plausible synthetic posterior for the two-event gamma
            % difference model. Log-scale throughout, matching the Stan model.
            rng(seed);
            p.b      = [log(5.5) + 0.10*randn(n,1), log(4.5) + 0.10*randn(n,1)];
            p.b_nu   = [log(18)  + 0.10*randn(n,1), log(20)  + 0.10*randn(n,1)];
            p.sd_id  = [abs(0.35 + 0.05*randn(n,1)), abs(0.32 + 0.05*randn(n,1)), ...
                        abs(0.30 + 0.05*randn(n,1)), abs(0.28 + 0.05*randn(n,1))];
            p.sd_trl = [abs(0.20 + 0.04*randn(n,1)), abs(0.18 + 0.04*randn(n,1))];
            p.cor_id = TestGammaDiffDynrelReliability.makeCorId(n, ...
                0.60 + 0.05*randn(n,1), -0.30 + 0.05*randn(n,1), ...
                -0.25 + 0.05*randn(n,1));
            p.cor_i  = 0.15 + 0.05*randn(n,1);
        end

        function [G,D,ICCg,ICCd] = staticRef(p, obs, ci)
            % The STATIC gamma difference reference, assembled exactly as
            % psyrat_gamma_extract_diff does before handing off to psyrat_diffrel.
            sp_a = p.sd_id(:,1); sp_b = p.sd_id(:,2);
            sv_a = p.sd_id(:,3); sv_b = p.sd_id(:,4);
            si_a = p.sd_trl(:,1); si_b = p.sd_trl(:,2);
            cov_pv_a = p.cor_id(:,1,3) .* sp_a .* sv_a;
            cov_pv_b = p.cor_id(:,2,4) .* sp_b .* sv_b;

            vc_a = psyrat_gamma_varcomps(p.b(:,1), sp_a, si_a, exp(p.b_nu(:,1)), ...
                sv_a, cov_pv_a);
            vc_b = psyrat_gamma_varcomps(p.b(:,2), sp_b, si_b, exp(p.b_nu(:,2)), ...
                sv_b, cov_pv_b);
            cc = psyrat_gamma_crosscov(p.b(:,1), sp_a, si_a, ...
                p.b(:,2), sp_b, si_b, p.cor_id(:,1,2), p.cor_i);

            nd = size(p.b,1);
            idvc = zeros(nd,2,2); trvc = zeros(nd,2,2);
            idvc(:,1,1) = vc_a.sigma_p2;  idvc(:,2,2) = vc_b.sigma_p2;
            idvc(:,1,2) = cc.cov_p_obs;   idvc(:,2,1) = cc.cov_p_obs;
            trvc(:,1,1) = vc_a.sigma_i2;  trvc(:,2,2) = vc_b.sigma_i2;
            trvc(:,1,2) = cc.cov_i_obs;   trvc(:,2,1) = cc.cov_i_obs;
            er = [0.5*log(vc_a.sigma_pi_e2(:)), 0.5*log(vc_b.sigma_pi_e2(:))];

            g = psyrat_diffrel('bp',idvc,'bt',trvc,'er_var',er,'obs',obs,...
                'CI',ci,'est','gen','wp_cov',0);
            d = psyrat_diffrel('bp',idvc,'bt',trvc,'er_var',er,'obs',obs,...
                'CI',ci,'est','dep','wp_cov',0);
            G = g.pt; D = d.pt; ICCg = g.icc_pt; ICCd = d.icc_pt;
        end
    end

    methods (Test)

        function testZeroSlopeReducesToStaticGammaDifference(testCase)
            % DECISIVE. b_dim = b_nu_dim = 0 -> the surface is constant across z
            % and equals the static gamma difference from the same components.
            n = 3000;
            p = testCase.makeParams(n, 2201);
            zgrid = [-1 -0.5 0 0.5 1];
            obs = [12 30]; ci = 0.95;   % deliberately unbalanced
            z0 = zeros(n,2);            % KDIM = 2 for one dimension

            rel = psyrat_rel_diffdynrel_gamma('b',p.b,'b_nu',p.b_nu,...
                'b_dim',z0,'b_nu_dim',z0,'sd_id',p.sd_id,'sd_trl',p.sd_trl,...
                'cor_id',p.cor_id,'cor_i',p.cor_i,'ndim',1,'z1',zgrid,...
                'obs',obs,'CI',ci);

            [Gref,Dref,ICCgref,ICCdref] = testCase.staticRef(p, obs, ci);

            testCase.verifyEqual(rel.G.pt(:), repmat(rel.G.pt(1),numel(zgrid),1), ...
                'AbsTol', 1e-12, 'G(z) is not constant at zero slope.');
            testCase.verifyEqual(rel.D.pt(:), repmat(rel.D.pt(1),numel(zgrid),1), ...
                'AbsTol', 1e-12, 'D(z) is not constant at zero slope.');
            testCase.verifyEqual(rel.G.pt(1), Gref, 'AbsTol', 1e-12, ...
                'G(z) does not reduce to the static gamma difference G.');
            testCase.verifyEqual(rel.D.pt(1), Dref, 'AbsTol', 1e-12, ...
                'D(z) does not reduce to the static gamma difference D.');
            testCase.verifyEqual(rel.ICCg.pt(1), ICCgref, 'AbsTol', 1e-12, ...
                'ICCg(z) does not reduce to the static gamma difference ICC.');
            testCase.verifyEqual(rel.ICCd.pt(1), ICCdref, 'AbsTol', 1e-12, ...
                'ICCd(z) does not reduce to the static gamma difference ICC.');
        end

        function testTwoDimensionZeroSlopeReduces(testCase)
            % Same reduction for ndim = 2, which pins the KDIM = 6 event-crossed
            % slope layout and the nz1 x nz2 grid shape.
            n = 2000;
            p = testCase.makeParams(n, 2202);
            z1 = [-1 0 1]; z2 = [-1 0 1];
            obs = [15 15]; ci = 0.95;
            z0 = zeros(n,6);

            rel = psyrat_rel_diffdynrel_gamma('b',p.b,'b_nu',p.b_nu,...
                'b_dim',z0,'b_nu_dim',z0,'sd_id',p.sd_id,'sd_trl',p.sd_trl,...
                'cor_id',p.cor_id,'cor_i',p.cor_i,'ndim',2,'z1',z1,'z2',z2,...
                'obs',obs,'CI',ci);

            Gref = testCase.staticRef(p, obs, ci);
            testCase.verifySize(rel.G.pt, [numel(z1) numel(z2)], ...
                'Two-dimension surface has the wrong grid shape.');
            testCase.verifyEqual(rel.G.pt, repmat(Gref,numel(z1),numel(z2)), ...
                'AbsTol', 1e-12, ...
                'Two-dimension G(z1,z2) does not reduce to the static gamma G.');
        end

        function testWrongSlopeWidthErrors(testCase)
            % KDIM must be 2 (one dimension) or 6 (two); a mis-shaped slope block
            % is a wiring error and must fail loudly rather than broadcast.
            n = 50;
            p = testCase.makeParams(n, 2203);
            args = {'b',p.b,'b_nu',p.b_nu,'sd_id',p.sd_id,'sd_trl',p.sd_trl,...
                'cor_id',p.cor_id,'cor_i',p.cor_i,'ndim',1,'z1',[-1 0 1],...
                'obs',[10 10],'CI',0.95};
            testCase.verifyError(@() psyrat_rel_diffdynrel_gamma(args{:},...
                'b_dim',zeros(n,3),'b_nu_dim',zeros(n,2)), 'varargin:bdim');
            testCase.verifyError(@() psyrat_rel_diffdynrel_gamma(args{:},...
                'b_dim',zeros(n,2),'b_nu_dim',zeros(n,3)), 'varargin:bnudim');
        end

        function testPositiveNuSlopeRaisesG(testCase)
            % Directional. nu is the scaled-chi-square degrees of freedom and the
            % conditional variance is 2*mu^2/nu, so a POSITIVE log-nu slope means
            % less dispersion as z rises -> smaller error -> higher G. This is the
            % mechanism the owner-supplied reference data actually exhibit (its
            % largest dimension coefficient is a nu-side slope).
            n = 2000;
            p = testCase.makeParams(n, 2204);
            zgrid = [-1 0 1]; obs = [20 20]; ci = 0.95;
            z0 = zeros(n,2);
            bnud = repmat(0.35, n, 2);   % same positive slope on both events

            rel = psyrat_rel_diffdynrel_gamma('b',p.b,'b_nu',p.b_nu,...
                'b_dim',z0,'b_nu_dim',bnud,'sd_id',p.sd_id,'sd_trl',p.sd_trl,...
                'cor_id',p.cor_id,'cor_i',p.cor_i,'ndim',1,'z1',zgrid,...
                'obs',obs,'CI',ci);

            testCase.verifyGreaterThan(rel.G.pt(3), rel.G.pt(2), ...
                'A positive log-nu slope should raise G with z.');
            testCase.verifyGreaterThan(rel.G.pt(2), rel.G.pt(1), ...
                'A positive log-nu slope should raise G with z.');
            % the residual contrast SD must fall as nu rises
            testCase.verifyLessThan(rel.sige_pt(3), rel.sige_pt(1), ...
                'Residual SD should fall as the dispersion df rises.');
        end

        function testMeanSlopeAloneMovesComponents(testCase)
            % STRUCTURAL: the gamma-vs-Gaussian distinction. Under the log link a
            % pure MEAN slope (b_dim ~= 0, b_nu_dim = 0) rescales every
            % observed-scale component by exp(2*alpha(z)+V), so the surface still
            % varies with z. Under the Gaussian identity link it could not - which
            % is exactly why psyrat_rel_diffdynrel's fixed id_varcov/trl_varcov
            % cannot be reused for this family.
            n = 2000;
            p = testCase.makeParams(n, 2205);
            zgrid = [-1 0 1]; obs = [20 20]; ci = 0.95;
            z0 = zeros(n,2);
            bdim = repmat(0.40, n, 2);   % mean slope only

            rel = psyrat_rel_diffdynrel_gamma('b',p.b,'b_nu',p.b_nu,...
                'b_dim',bdim,'b_nu_dim',z0,'sd_id',p.sd_id,'sd_trl',p.sd_trl,...
                'cor_id',p.cor_id,'cor_i',p.cor_i,'ndim',1,'z1',zgrid,...
                'obs',obs,'CI',ci);

            testCase.verifyGreaterThan(abs(rel.sige_pt(3) - rel.sige_pt(1)), 1e-6, ...
                ['A mean-only slope must still move the observed-scale residual '...
                 'under the log link (mean-coupled dispersion).']);
        end

        function testPerEventTrialCountsAreHonored(testCase)
            % Difference designs are routinely unbalanced (the reference ERP design
            % is 49 error vs 341 correct trials). A per-event obs must differ from
            % the same-mean scalar, and must match psyrat_diffrel driven directly
            % with that per-event count.
            n = 1500;
            p = testCase.makeParams(n, 2206);
            ci = 0.95; z0 = zeros(n,2);
            common = {'b',p.b,'b_nu',p.b_nu,'b_dim',z0,'b_nu_dim',z0,...
                'sd_id',p.sd_id,'sd_trl',p.sd_trl,'cor_id',p.cor_id,...
                'cor_i',p.cor_i,'ndim',1,'z1',0,'CI',ci};

            relUneq = psyrat_rel_diffdynrel_gamma(common{:},'obs',[49 341]);
            relScal = psyrat_rel_diffdynrel_gamma(common{:},'obs',195);

            testCase.verifyGreaterThan(abs(relUneq.G.pt(1) - relScal.G.pt(1)), 1e-8, ...
                ['Per-event n'' must not collapse to the scalar mean; the '...
                 'unbalanced design should give a different coefficient.']);

            Gref = testCase.staticRef(p, [49 341], ci);
            testCase.verifyEqual(relUneq.G.pt(1), Gref, 'AbsTol', 1e-12, ...
                'Per-event n'' surface does not match psyrat_diffrel at the same n''.');

            % a scalar obs is broadcast to both events
            relBc = psyrat_rel_diffdynrel_gamma(common{:},'obs',[195 195]);
            testCase.verifyEqual(relScal.G.pt(1), relBc.G.pt(1), 'AbsTol', 1e-12, ...
                'A scalar obs should broadcast to both events.');
        end

        function testSummaryRoutesGammaAwayFromGaussianDiffBranch(testCase)
            % REGRESSION GUARD for the highest-risk wiring in this feature.
            % psyrat_dynrel_summary decides the difference branch on
            % REL.analysis alone ('ic_diff_dynrel'), which is the SAME label the
            % gamma family emits. Without a family test it would route a gamma
            % result into the Gaussian branch and immediately read
            % REL.out.b_sigma / b_sigma_dim -- parameters a scaled-chi-square model
            % never emits. This builds a gamma case-12 REL and asserts the summary
            % produces a real surface (which is only reachable via the gamma
            % branch, since the gamma REL.out has no b_sigma at all).
            n = 400;
            p = testCase.makeParams(n, 2207);
            nsub = 30;

            REL = struct();
            REL.analysis = 'ic_diff_dynrel';
            REL.family = 'gamma';
            REL.ndim = 1;
            REL.dim_names = {'dim1'};
            REL.out = struct();
            REL.out.labels = {'none'};
            REL.out.b        = {num2cell(p.b)};
            REL.out.b_nu     = {num2cell(p.b_nu)};
            REL.out.b_dim    = {num2cell(zeros(n,2))};
            REL.out.b_nu_dim = {num2cell(zeros(n,2))};
            REL.out.sd_id    = {num2cell(p.sd_id)};
            REL.out.sd_trl   = {num2cell(p.sd_trl)};
            REL.out.cor_id   = {num2cell(p.cor_id)};
            REL.out.cor_i    = {num2cell(p.cor_i)};
            REL.out.dimz     = {linspace(-1.5,1.5,nsub)'};
            REL.out.ntrials    = 20;
            REL.out.ntrials_ev = [12; 30];
            REL.out.elabels  = {'e1';'e2'};
            REL.out.glabels  = {'none'};

            s = psyrat_dynrel_summary(REL,'CI',0.95,'ngrid',5);

            testCase.verifyEqual(numel(s.strata), 1);
            testCase.verifySize(s.strata(1).G.pt, [1 5], ...
                'Gamma difference surface has the wrong grid shape.');
            testCase.verifyTrue(all(isfinite(s.strata(1).G.pt)), ...
                'Gamma difference surface contains non-finite values.');
            testCase.verifyTrue(all(s.strata(1).G.pt > 0 & s.strata(1).G.pt < 1), ...
                'Gamma difference G(z) outside (0,1).');
            % the per-event counts must survive into the summary
            testCase.verifyEqual(s.strata(1).obs_ev, [12 30], ...
                'Per-event n'' did not reach the summary.');
            testCase.verifyNotEmpty(s.strata(1).varcomp, ...
                'Gamma difference variance-component table is empty.');
            % zero slopes -> flat surface, tying the summary path to the reduction
            testCase.verifyEqual(s.strata(1).G.pt(:), ...
                repmat(s.strata(1).G.pt(1),5,1), 'AbsTol', 1e-12, ...
                'Zero-slope gamma difference surface is not flat via the summary.');
        end
    end
end
