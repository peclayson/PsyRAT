classdef TestGammaDynrelLsReliability < PsyRATTestBase
    % Offline proof of the LOCATION-SCALE Gamma (scaled-chi-square) dynamic/
    % conditional reliability SURFACE (psyrat_rel_dynrel_gamma_ls; analysis 11
    % with family = 'gamma' and gammascale = 2). No CmdStan: the surface is a
    % deterministic function of the posterior draws, so synthetic draws exercise
    % the full log -> observed conversion plus the coefficient math.
    %
    % The anchoring check is the zero-slope reduction: with b = b_sigma = 0 the
    % surface must be constant across z and equal to psyrat_gamma_varcomps_ls fed
    % straight into psyrat_rel_sing, which is the same already-validated
    % coefficient the static one-facet gamma path uses.
    %
    % The DISCRIMINATING check is testMeanSlopeMovesReliability. Under the log-nu
    % sibling (gammascale = 1) every observed-scale component - the residual
    % included - carries a common exp(2*alpha(z)) factor, so a mean-level slope
    % cancels exactly and only b_nu moves the coefficient; TestGammaDynrelReliability
    % relies on that. Here the residual is exp(2*log_sigma + 2*s_sd^2) and carries
    % no alpha at all, so the cancellation is GONE and a mean slope moves
    % reliability on its own. A test that merely re-ran the log-nu expectations
    % would pass on either parameterization and prove nothing; this one fails on
    % the wrong model.
    %
    % All draws use base-MATLAB randn only (no Statistics Toolbox).

    methods (Test)

        function testZeroSlopeReducesToTheLocationScaleConverter(testCase)
            % ANCHOR. b = b_sigma = 0 -> mu and sigma do not depend on z, so the
            % surface is flat and equal to the converter + psyrat_rel_sing.
            rng(2101);
            n = 4000;
            alpha0  = log(5.5) + 0.1*randn(n,1);
            logsig0 = log(1.3) + 0.1*randn(n,1);
            s_p     = abs(0.35 + 0.05*randn(n,1));
            s_sd    = abs(0.30 + 0.05*randn(n,1));
            sig_i   = abs(0.35 + 0.05*randn(n,1));
            zgrid   = [-1 -0.5 0 0.5 1];
            obs = 8; ci = 0.95;

            b0 = zeros(n,1);
            rel = psyrat_rel_dynrel_gamma_ls('alpha0',alpha0,'b',b0,...
                'logsig0',logsig0,'b_sigma',b0,'sig_p',s_p,'sig_sd',s_sd,...
                'sig_i',sig_i,'ndim',1,'z1',zgrid,'obs',obs,'CI',ci);

            vc = psyrat_gamma_varcomps_ls(alpha0, s_p, sig_i, logsig0, s_sd);
            bp = sqrt(vc.sigma_p2);
            iz = sqrt(vc.sigma_i2);
            wp = sqrt(vc.sigma_pi_e2 + vc.sigma_i2);
            [~,Gref] = psyrat_rel_sing('gcoeff',2,'metric','global','bp',bp,...
                'wp',wp,'i',iz,'obs',obs,'CI',ci);
            [~,Dref] = psyrat_rel_sing('gcoeff',1,'metric','global','bp',bp,...
                'wp',wp,'i',iz,'obs',obs,'CI',ci);

            testCase.verifyEqual(rel.G.pt(:), repmat(rel.G.pt(1),numel(zgrid),1), ...
                'AbsTol', 1e-12, 'G(z) not constant at zero slopes.');
            testCase.verifyEqual(rel.G.pt(1), Gref, 'AbsTol', 1e-12, ...
                'G(z) does not reduce to the location-scale converter.');
            testCase.verifyEqual(rel.D.pt(1), Dref, 'AbsTol', 1e-12, ...
                'D(z) does not reduce to the location-scale converter.');
        end

        function testResidualSdSlopeLowersReliability(testCase)
            % A positive log residual-SD slope raises sigma(z) with z -> more
            % observation noise -> lower reliability. Strictly decreasing, in
            % (0,1). Note the sign is OPPOSITE the log-nu sibling's b_nu, where a
            % positive slope raises nu and therefore lowers dispersion.
            rng(2110);
            n = 3000;
            alpha0  = log(5.5) + 0.1*randn(n,1);
            logsig0 = log(1.3) + 0.1*randn(n,1);
            s_p     = abs(0.35 + 0.05*randn(n,1));
            s_sd    = abs(0.30 + 0.05*randn(n,1));
            sig_i   = abs(0.35 + 0.05*randn(n,1));
            zgrid   = [-1 -0.5 0 0.5 1];

            rel = psyrat_rel_dynrel_gamma_ls('alpha0',alpha0,'b',zeros(n,1),...
                'logsig0',logsig0,'b_sigma',0.4*ones(n,1),'sig_p',s_p,...
                'sig_sd',s_sd,'sig_i',sig_i,'ndim',1,'z1',zgrid,'obs',8,'CI',0.95);
            G = rel.G.pt(:);
            testCase.verifyTrue(all(diff(G) < 0), ...
                'G(z) should decrease with z for a positive residual-SD slope.');
            testCase.verifyTrue(all(G > 0 & G < 1), 'G(z) out of (0,1).');
        end

        function testMeanSlopeMovesReliability(testCase)
            % DISCRIMINATING. With b_sigma = 0 and a positive mean slope, the
            % log-nu parameterization would return a FLAT surface (its residual
            % carries the same exp(2*alpha(z)) factor as the numerator). Here the
            % residual is alpha-free, so
            %   G(z) = A / ( A + (V-A-B)/n' + E/(n'*exp(2*alpha(z))) )
            % and the last term shrinks as alpha(z) grows: a fixed absolute
            % residual SD is a smaller share of a larger expected score. G must
            % therefore INCREASE strictly with z, and by a margin far above
            % numerical noise.
            rng(2120);
            n = 3000;
            alpha0  = log(5.5) + 0.1*randn(n,1);
            logsig0 = log(1.3) + 0.1*randn(n,1);
            s_p     = abs(0.35 + 0.05*randn(n,1));
            s_sd    = abs(0.30 + 0.05*randn(n,1));
            sig_i   = abs(0.35 + 0.05*randn(n,1));
            zgrid   = [-1 -0.5 0 0.5 1];

            rel = psyrat_rel_dynrel_gamma_ls('alpha0',alpha0,'b',0.4*ones(n,1),...
                'logsig0',logsig0,'b_sigma',zeros(n,1),'sig_p',s_p,...
                'sig_sd',s_sd,'sig_i',sig_i,'ndim',1,'z1',zgrid,'obs',8,'CI',0.95);
            G = rel.G.pt(:);

            testCase.verifyTrue(all(diff(G) > 0), ...
                ['G(z) must increase with z for a positive MEAN slope under ' ...
                'gammascale = 2. A flat surface here means the residual has ' ...
                'picked up an exp(2*alpha) factor it must not have.']);
            % Guard against a pass on floating-point drift: the movement across
            % the grid has to be a real effect, not the ~1e-15 the log-nu
            % cancellation leaves behind.
            testCase.verifyGreaterThan(G(end) - G(1), 1e-3, ...
                'Mean-slope effect is too small to distinguish from round-off.');
            testCase.verifyTrue(all(G > 0 & G < 1), 'G(z) out of (0,1).');
        end

        function testPersonScaleSdRaisesTheMarginalResidual(testCase)
            % s_sd enters ONLY through E[sigma_p^2] = exp(2*log_sigma + 2*s_sd^2),
            % so heterogeneity in person precision inflates the marginal residual
            % and lowers reliability. Pins that the person scale SD is actually
            % propagated rather than silently dropped - the failure mode that
            % would follow from passing gro_sds column 2 to the wrong argument.
            rng(2130);
            n = 2000;
            alpha0  = log(5.5) + 0.1*randn(n,1);
            logsig0 = log(1.3) + 0.1*randn(n,1);
            s_p     = abs(0.35 + 0.05*randn(n,1));
            sig_i   = abs(0.35 + 0.05*randn(n,1));
            common = {'alpha0',alpha0,'b',zeros(n,1),'logsig0',logsig0,...
                'b_sigma',zeros(n,1),'sig_p',s_p,'sig_i',sig_i,...
                'ndim',1,'z1',0,'obs',8,'CI',0.95};

            relFlat = psyrat_rel_dynrel_gamma_ls(common{:},'sig_sd',zeros(n,1));
            relHet  = psyrat_rel_dynrel_gamma_ls(common{:},'sig_sd',0.5*ones(n,1));

            testCase.verifyLessThan(relHet.G.pt(1), relFlat.G.pt(1), ...
                'Person residual-SD heterogeneity must lower the marginal G.');
        end

        function testTwoDimensionSurfaceShape(testCase)
            % Two-dimension (KDIM = 3) surface returns an nz1 x nz2 grid.
            rng(2140);
            n = 1500;
            alpha0  = log(5.5) + 0.1*randn(n,1);
            logsig0 = log(1.3) + 0.1*randn(n,1);
            s_p     = abs(0.35 + 0.05*randn(n,1));
            s_sd    = abs(0.30 + 0.05*randn(n,1));
            sig_i   = abs(0.35 + 0.05*randn(n,1));
            rel = psyrat_rel_dynrel_gamma_ls('alpha0',alpha0,'b',0.2*ones(n,3),...
                'logsig0',logsig0,'b_sigma',0.1*ones(n,3),'sig_p',s_p,...
                'sig_sd',s_sd,'sig_i',sig_i,'ndim',2,'z1',[-1 0 1],'z2',[-1 0 1],...
                'obs',8,'CI',0.95);
            testCase.verifyEqual(size(rel.G.pt), [3 3], ...
                'Two-dimension G surface has the wrong size.');
        end

        function testSlopeBlockShapeIsChecked(testCase)
            % A transposed or mis-shaped slope block must fail loudly rather than
            % broadcast into a wrong answer. Both blocks are checked separately,
            % because b_sigma is the one this parameterization adds.
            rng(2150);
            n = 200;
            base = {'alpha0',randn(n,1),'logsig0',randn(n,1),...
                'sig_p',abs(randn(n,1)),'sig_sd',abs(randn(n,1)),...
                'sig_i',abs(randn(n,1)),'ndim',1,'z1',[-1 0 1],'obs',8,'CI',0.95};
            testCase.verifyError(@() psyrat_rel_dynrel_gamma_ls(...
                'b',zeros(n,3),'b_sigma',zeros(n,1),base{:}), 'varargin:b');
            testCase.verifyError(@() psyrat_rel_dynrel_gamma_ls(...
                'b',zeros(n,1),'b_sigma',zeros(n,3),base{:}), 'varargin:bsigma');
        end

    end
end
