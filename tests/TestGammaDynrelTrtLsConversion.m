classdef TestGammaDynrelTrtLsConversion < PsyRATTestBase
    % Offline proof of the LOCATION-SCALE Gamma (scaled-chi-square) TWO-FACET
    % dynamic/conditional reliability calculators: psyrat_rel_dynrel_trt_gamma_ls
    % (analysis 14, the group-level surface) and psyrat_ssrel_dynrel_trt_gamma_ls
    % (analysis 27, the per-participant conditional read-out), both with
    % family = 'gamma' and gammascale = 2. No CmdStan: the surface is a
    % deterministic function of the posterior draws, so synthetic draws exercise
    % the full log -> observed conversion plus the coefficient math.
    %
    % THREE CHECKS CARRY THE INCREMENT, and each fails on a different wrong model.
    %
    % 1. testSignalComponentsAreParameterizationInvariant is the STRUCTURAL claim
    %    this increment rests on: the six observed-scale signal components are
    %    built from the log-mean submodel alone, so psyrat_gamma_varcomps_trt_ls
    %    and psyrat_gamma_varcomps_trt must agree on all six - including the five
    %    crossed facets - at the same alpha, no matter what the dispersion
    %    arguments are. Only the residual may differ. This is what licenses
    %    reusing the mean-side machinery unchanged, and it is asserted rather than
    %    cited.
    %
    % 2. testMeanSlopeMovesReliability is the DISCRIMINATING check. Under the
    %    log-nu sibling every observed-scale component - the residual included -
    %    carries a common exp(2*alpha(z)) factor, so a mean-level slope cancels
    %    exactly and only b_nu moves the coefficient. Here e_cond_var is
    %    exp(2*log_sigma + 2*s_sd^2) and carries no alpha, so the cancellation is
    %    GONE. A test that merely re-ran the log-nu expectations would pass on
    %    either parameterization and prove nothing; this one fails on the wrong
    %    model.
    %
    % 3. testSubjectResidualCarriesNoWithinPersonFactor pins the two-facet form of
    %    the rule increment 5 established for one facet. The log-nu twin's
    %    per-participant residual carries exp(2*V_w) with
    %    V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2, because its conditional
    %    variance 2*mu^2/nu is an expectation over the within-person facets. Under
    %    location-scale sigma_s is a single per-participant quantity carrying no
    %    facet terms, so there is no expectation to take and NO such factor. Adding
    %    one back "for symmetry" would silently re-couple the residual to the mean
    %    surface.
    %
    % Constant posterior draws are used wherever a closed form is checked:
    % quantiles of a constant collapse to that constant, so the point estimate can
    % be compared at machine tolerance. All draws use base-MATLAB randn/rand only
    % (no Statistics Toolbox), so this runs on CI's --products=MATLAB image.

    methods (Static, Access = private)

        function y = rep(x, nd)
            % nd identical draws of a scalar parameter
            y = repmat(x, nd, 1);
        end

        function t = ssinfo2(z, trls, nocc)
            % Per-subject lookup in the two-facet schema local_dynrel_ssinfo_trt
            % builds: id, id2, z1, trls, occs. id2 = 1:NSUB indexes the ind_sd
            % draw columns.
            n = numel(z);
            t = table((1:n)', (1:n)', z(:), trls(:), repmat(nocc, n, 1), ...
                'VariableNames', {'id','id2','z1','trls','occs'});
        end

        function args = surfaceArgs(n)
            % A shared, well-conditioned draw set for the surface tests. Returned
            % as a name-value list missing only b, b_sigma and the grid.
            args = {'alpha0', log(5.5) + 0.1*randn(n,1), ...
                'logsig0', log(1.3) + 0.1*randn(n,1), ...
                'sig_p',   abs(0.40 + 0.05*randn(n,1)), ...
                'sig_sd',  abs(0.30 + 0.05*randn(n,1)), ...
                'sig_occ', abs(0.22 + 0.03*randn(n,1)), ...
                'sig_trl', abs(0.18 + 0.03*randn(n,1)), ...
                'sig_trlxid',  abs(0.14 + 0.03*randn(n,1)), ...
                'sig_occxid',  abs(0.12 + 0.03*randn(n,1)), ...
                'sig_trlxocc', abs(0.09 + 0.02*randn(n,1))};
        end

    end

    methods (Test)

        % ---------- analysis 14: the group-level surface ----------

        function testSignalComponentsAreParameterizationInvariant(testCase)
            % THE STRUCTURAL CLAIM. All six signal components, plus mu_bar and
            % var_mu, must be IDENTICAL between the two converters at the same
            % alpha and the same five facet SDs, because every one of them is a
            % function of the log-MEAN submodel only. The dispersion arguments are
            % deliberately chosen to be wildly different between the two calls: if
            % any signal component leaked a dependence on nu, s_v, cov_pv,
            % log_sigma or s_sd, this fails.
            %
            % Only e_cond_var (and hence total_var / sigma_res2, which are derived
            % from it) is allowed to differ. That is asserted too, so the test
            % cannot pass by the two converters being accidentally identical.
            rng(3101);
            nd = 500;
            alpha = log(5.5) + 0.1*randn(nd,1);
            s_p  = abs(0.40 + 0.05*randn(nd,1));
            s_o  = abs(0.22 + 0.03*randn(nd,1));
            s_t  = abs(0.18 + 0.03*randn(nd,1));
            s_po = abs(0.12 + 0.03*randn(nd,1));
            s_pt = abs(0.14 + 0.03*randn(nd,1));
            s_ot = abs(0.09 + 0.02*randn(nd,1));

            vcNu = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, ...
                s_ot, exp(log(18) + 0.2*randn(nd,1)), 1.2*ones(nd,1), ...
                0.3*ones(nd,1));
            vcLs = psyrat_gamma_varcomps_trt_ls(alpha, s_p, s_o, s_t, s_po, ...
                s_pt, s_ot, log(1.3) + 0.1*randn(nd,1), 0.5*ones(nd,1));

            shared = {'mu_bar','sigma_p2','sigma_o2','sigma_t2','sigma_po2', ...
                'sigma_pt2','sigma_ot2','var_mu'};
            for k = 1:numel(shared)
                testCase.verifyEqual(vcLs.(shared{k}), vcNu.(shared{k}), ...
                    'AbsTol', 0, sprintf(['%s must be IDENTICAL across the two ' ...
                    'parameterizations - it is a log-mean submodel quantity. A ' ...
                    'difference here means the scale submodel has leaked into a ' ...
                    'signal component.'], shared{k}));
            end

            % ...and the residual MUST differ, or the test above is vacuous.
            testCase.verifyNotEqual(vcLs.e_cond_var, vcNu.e_cond_var, ...
                ['e_cond_var is identical across the parameterizations, so the ' ...
                'invariance check above proves nothing. Check the draw setup.']);
        end

        function testZeroSlopeReducesToTheTwoFacetLocationScaleConverter(testCase)
            % ANCHOR. b = b_sigma = 0 -> mu and sigma do not depend on z, so the
            % surface is flat and equal to psyrat_gamma_varcomps_trt_ls fed
            % straight into the already-validated psyrat_rel_trt.
            rng(3110);
            n = 3000;
            a = TestGammaDynrelTrtLsConversion.surfaceArgs(n);
            zgrid = [-1 -0.5 0 0.5 1];
            obs = 12; nocc = 2; reltype = 3; ci = 0.95;
            b0 = zeros(n,1);

            rel = psyrat_rel_dynrel_trt_gamma_ls(a{:}, 'b', b0, 'b_sigma', b0, ...
                'ndim', 1, 'z1', zgrid, 'obs', obs, 'nocc', nocc, ...
                'reltype', reltype, 'CI', ci);

            get = @(name) a{find(strcmp(a, name), 1) + 1};
            vc = psyrat_gamma_varcomps_trt_ls(get('alpha0'), get('sig_p'), ...
                get('sig_occ'), get('sig_trl'), get('sig_occxid'), ...
                get('sig_trlxid'), get('sig_trlxocc'), get('logsig0'), ...
                get('sig_sd'));
            ref = {'reltype', reltype, 'bp', sqrt(vc.sigma_p2), ...
                'bo', sqrt(vc.sigma_o2), 'bt', sqrt(vc.sigma_t2), ...
                'txp', sqrt(vc.sigma_pt2), 'oxp', sqrt(vc.sigma_po2), ...
                'txo', sqrt(vc.sigma_ot2), 'err', sqrt(vc.sigma_res2), ...
                'obs', obs, 'nocc', nocc, 'CI', ci};
            [~, Gref] = psyrat_rel_trt('gcoeff', 2, ref{:});
            [~, Dref] = psyrat_rel_trt('gcoeff', 1, ref{:});

            testCase.verifyEqual(rel.G.pt(:), repmat(rel.G.pt(1), numel(zgrid), 1), ...
                'AbsTol', 1e-12, 'G(z) not constant at zero slopes.');
            testCase.verifyEqual(rel.G.pt(1), Gref, 'AbsTol', 1e-12, ...
                'G(z) does not reduce to the two-facet location-scale converter.');
            testCase.verifyEqual(rel.D.pt(1), Dref, 'AbsTol', 1e-12, ...
                'D(z) does not reduce to the two-facet location-scale converter.');
        end

        function testResidualSdSlopeLowersReliability(testCase)
            % A positive log residual-SD slope raises sigma(z) with z -> more
            % observation noise -> lower reliability. Strictly decreasing, in
            % (0,1). Note the sign is OPPOSITE the log-nu sibling's b_nu, where a
            % positive slope raises nu and therefore LOWERS dispersion.
            rng(3120);
            n = 2500;
            a = TestGammaDynrelTrtLsConversion.surfaceArgs(n);
            rel = psyrat_rel_dynrel_trt_gamma_ls(a{:}, 'b', zeros(n,1), ...
                'b_sigma', 0.4*ones(n,1), 'ndim', 1, 'z1', [-1 -0.5 0 0.5 1], ...
                'obs', 12, 'nocc', 2, 'reltype', 3, 'CI', 0.95);
            G = rel.G.pt(:);
            testCase.verifyTrue(all(diff(G) < 0), ...
                'G(z) should decrease with z for a positive residual-SD slope.');
            testCase.verifyTrue(all(G > 0 & G < 1), 'G(z) out of (0,1).');
        end

        function testMeanSlopeMovesReliability(testCase)
            % DISCRIMINATING. With b_sigma = 0 and a positive mean slope, the
            % log-nu parameterization would return a FLAT surface (its residual
            % carries the same exp(2*alpha(z)) factor as every signal component).
            % Here the residual is alpha-free apart from the lognormal mean-surface
            % part, so the alpha-free e_cond_var becomes a smaller share of a
            % larger expected score and G must INCREASE strictly with z, by a
            % margin far above numerical noise.
            rng(3130);
            n = 2500;
            a = TestGammaDynrelTrtLsConversion.surfaceArgs(n);
            rel = psyrat_rel_dynrel_trt_gamma_ls(a{:}, 'b', 0.4*ones(n,1), ...
                'b_sigma', zeros(n,1), 'ndim', 1, 'z1', [-1 -0.5 0 0.5 1], ...
                'obs', 12, 'nocc', 2, 'reltype', 3, 'CI', 0.95);
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

        function testPersonScaleSdLowersReliability(testCase)
            % s_sd enters ONLY through E[sigma_p^2] = exp(2*log_sigma + 2*s_sd^2),
            % so heterogeneity in person precision inflates the marginal residual
            % and lowers reliability. Pins that the person scale SD is actually
            % propagated rather than silently dropped - the failure mode that would
            % follow from reading gro_sds column 2 into the wrong argument.
            rng(3140);
            n = 2000;
            a = TestGammaDynrelTrtLsConversion.surfaceArgs(n);
            a(find(strcmp(a,'sig_sd'),1)+1) = [];  %#ok<FNDSB> drop the value
            a(strcmp(a,'sig_sd')) = [];            % ...and its name
            common = [a, {'b', zeros(n,1), 'b_sigma', zeros(n,1), 'ndim', 1, ...
                'z1', 0, 'obs', 12, 'nocc', 2, 'reltype', 3, 'CI', 0.95}];

            relFlat = psyrat_rel_dynrel_trt_gamma_ls(common{:}, 'sig_sd', zeros(n,1));
            relHet  = psyrat_rel_dynrel_trt_gamma_ls(common{:}, 'sig_sd', 0.5*ones(n,1));

            testCase.verifyLessThan(relHet.G.pt(1), relFlat.G.pt(1), ...
                'Person residual-SD heterogeneity must lower the marginal G.');
        end

        function testCrossedFacetsReachTheCoefficient(testCase)
            % Each of the five crossed facet SDs must actually be passed through to
            % the converter. Inflating any one of them alone changes the
            % coefficient; a dropped or mis-slotted argument would leave one of
            % these five identical to the baseline.
            rng(3150);
            n = 1200;
            names = {'sig_occ','sig_trl','sig_trlxid','sig_occxid','sig_trlxocc'};
            base = TestGammaDynrelTrtLsConversion.surfaceArgs(n);
            tail = {'b', zeros(n,1), 'b_sigma', zeros(n,1), 'ndim', 1, ...
                'z1', 0, 'obs', 12, 'nocc', 2, 'reltype', 3, 'CI', 0.95};
            G0 = psyrat_rel_dynrel_trt_gamma_ls(base{:}, tail{:}).G.pt(1);

            for k = 1:numel(names)
                bumped = base;
                idx = find(strcmp(bumped, names{k}), 1) + 1;
                bumped{idx} = bumped{idx} + 0.35;
                Gk = psyrat_rel_dynrel_trt_gamma_ls(bumped{:}, tail{:}).G.pt(1);
                testCase.verifyNotEqual(Gk, G0, sprintf( ...
                    ['Inflating %s left the coefficient unchanged, so that facet ' ...
                    'is not reaching psyrat_gamma_varcomps_trt_ls.'], names{k}));
            end
        end

        function testTwoDimensionSurfaceShape(testCase)
            % Two-dimension (KDIM = 3) surface returns an nz1 x nz2 grid.
            rng(3160);
            n = 900;
            a = TestGammaDynrelTrtLsConversion.surfaceArgs(n);
            rel = psyrat_rel_dynrel_trt_gamma_ls(a{:}, 'b', 0.2*ones(n,3), ...
                'b_sigma', 0.1*ones(n,3), 'ndim', 2, 'z1', [-1 0 1], ...
                'z2', [-1 0 1], 'obs', 12, 'nocc', 2, 'reltype', 3, 'CI', 0.95);
            testCase.verifyEqual(size(rel.G.pt), [3 3], ...
                'Two-dimension G surface has the wrong size.');
        end

        function testSlopeBlockShapeIsChecked(testCase)
            % A transposed or mis-shaped slope block must fail loudly rather than
            % broadcast into a wrong answer. Both blocks are checked separately,
            % because b_sigma is the one this parameterization adds.
            rng(3170);
            n = 200;
            a = TestGammaDynrelTrtLsConversion.surfaceArgs(n);
            tail = {'ndim', 1, 'z1', [-1 0 1], 'obs', 12, 'CI', 0.95};
            testCase.verifyError(@() psyrat_rel_dynrel_trt_gamma_ls(a{:}, ...
                'b', zeros(n,3), 'b_sigma', zeros(n,1), tail{:}), 'varargin:b');
            testCase.verifyError(@() psyrat_rel_dynrel_trt_gamma_ls(a{:}, ...
                'b', zeros(n,1), 'b_sigma', zeros(n,3), tail{:}), 'varargin:bsigma');
        end

        % ---------- analysis 27: the per-participant read-out ----------

        function testSubjectResidualCarriesNoWithinPersonFactor(testCase)
            % DECISIVE for analysis 27, and the two-facet form of increment 5's
            % rule. Each participant's coefficient must equal a direct
            % psyrat_rel_trt call built from the six population components at THAT
            % participant's z, with err = exp(log sigma(z_s) + w_s) EXACTLY - no
            % exp(2*V_w) within-person factor, no leading 2*, and no dependence on
            % the participant's own log-mean.
            %
            % Constant draws, so the point estimate is a closed form at machine
            % tolerance. The reference below deliberately does NOT contain V_w; if
            % the production code grew one, this fails.
            nd = 30;
            alpha0 = log(5.0); b = 0.20; logsig0 = log(1.4); b_sigma = -0.25;
            sig_p = 0.40; s_o = 0.22; s_t = 0.18; s_po = 0.12; s_pt = 0.14; s_ot = 0.09;
            reltype = 3; nocc = 2;

            z    = [-0.8;  0.0;  1.2];
            w    = [ 0.25; -0.15; 0.05];   % each person's own log residual-SD offset
            trls = [ 18;    24;   30];

            ssinfo = TestGammaDynrelTrtLsConversion.ssinfo2(z, trls, nocc);
            rp = @(x) TestGammaDynrelTrtLsConversion.rep(x, nd);

            [gen, dep] = psyrat_ssrel_dynrel_trt_gamma_ls( ...
                'alpha0', rp(alpha0), 'b', rp(b), ...
                'logsig0', rp(logsig0), 'b_sigma', rp(b_sigma), ...
                'sig_p', rp(sig_p), 'sig_occ', rp(s_o), 'sig_trl', rp(s_t), ...
                'sig_occxid', rp(s_po), 'sig_trlxid', rp(s_pt), ...
                'sig_trlxocc', rp(s_ot), ...
                'ind_sd', repmat(w(:)', nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'reltype', reltype, ...
                'nocc', nocc, 'CI', .95);

            for k = 1:numel(z)
                alpha_z  = alpha0  + b       * z(k);
                logsig_z = logsig0 + b_sigma * z(k);
                vc = psyrat_gamma_varcomps_trt_ls(alpha_z, sig_p, s_o, s_t, ...
                    s_po, s_pt, s_ot, logsig_z, 0);
                errSs = exp(logsig_z + w(k));   % NO V_w, NO 2*, NO mu

                args = {'reltype', reltype, ...
                    'bp', sqrt(vc.sigma_p2), 'bo', sqrt(vc.sigma_o2), ...
                    'bt', sqrt(vc.sigma_t2), 'txp', sqrt(vc.sigma_pt2), ...
                    'oxp', sqrt(vc.sigma_po2), 'txo', sqrt(vc.sigma_ot2), ...
                    'err', errSs, 'nocc', nocc, 'CI', .95};

                [~, G_ref] = psyrat_rel_trt('gcoeff', 2, args{:}, 'obs', trls(k));
                [~, D_ref] = psyrat_rel_trt('gcoeff', 1, args{:}, 'obs', trls(k));
                testCase.verifyEqual(gen.dep_pt(k), G_ref, 'RelTol', 1e-12, ...
                    ['Per-participant G does not match the conditional closed ' ...
                    'form. The most likely cause is a within-person exp(2*V_w) ' ...
                    'factor copied over from the log-nu twin, which has no ' ...
                    'counterpart under location-scale.']);
                testCase.verifyEqual(dep.dep_pt(k), D_ref, 'RelTol', 1e-12);
            end
        end

        function testSubjectComponentsMoveWithZ(testCase)
            % The own-z estimand (owner ruling, 2026-07-31), pinned as its log-nu
            % twin pins it. Two participants with the SAME residual-SD offset and a
            % z-invariant residual (b_sigma = 0) must still differ, because all six
            % signal components are evaluated at their own z. It fails if a future
            % change reverts to holding the components at a shared reference z.
            nd = 20;
            alpha0 = log(5.0); b = 0.30; logsig0 = log(1.4);
            sig_p = 0.40; s_o = 0.22; s_t = 0.18; s_po = 0.12; s_pt = 0.14; s_ot = 0.09;
            nocc = 2;

            z    = [-1.0; 1.0];
            w    = [ 0.10; 0.10];   % identical person effects
            trls = [ 24;   24];     % identical trial counts

            ssinfo = TestGammaDynrelTrtLsConversion.ssinfo2(z, trls, nocc);
            rp = @(x) TestGammaDynrelTrtLsConversion.rep(x, nd);

            gen = psyrat_ssrel_dynrel_trt_gamma_ls( ...
                'alpha0', rp(alpha0), 'b', rp(b), ...
                'logsig0', rp(logsig0), 'b_sigma', rp(0), ...
                'sig_p', rp(sig_p), 'sig_occ', rp(s_o), 'sig_trl', rp(s_t), ...
                'sig_occxid', rp(s_po), 'sig_trlxid', rp(s_pt), ...
                'sig_trlxocc', rp(s_ot), ...
                'ind_sd', repmat(w(:)', nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'reltype', 3, ...
                'nocc', nocc, 'CI', .95);

            testCase.verifyNotEqual(gen.dep_pt(1), gen.dep_pt(2), ...
                ['Two participants at different z with identical residuals got ' ...
                'identical coefficients, so the six signal components are not ' ...
                'being evaluated at each participant''s own z.']);
            % Direction: a positive mean slope raises the signal while the residual
            % is held fixed, so the higher-z participant must score higher.
            testCase.verifyGreaterThan(gen.dep_pt(2), gen.dep_pt(1));
        end

        function testSubjectShapeAndIndexGuards(testCase)
            % A transposed slope block or an id2 outside the per-subject draw
            % columns must fail loudly rather than silently assign one
            % participant's residual to another.
            nd = 10;
            ssinfo = TestGammaDynrelTrtLsConversion.ssinfo2([0;1], [20;20], 2);
            rp = @(x) TestGammaDynrelTrtLsConversion.rep(x, nd);
            common = {'alpha0', rp(log(5)), 'logsig0', rp(log(1.4)), ...
                'sig_p', rp(0.4), 'sig_occ', rp(0.2), 'sig_trl', rp(0.2), ...
                'sig_occxid', rp(0.1), 'sig_trlxid', rp(0.1), ...
                'sig_trlxocc', rp(0.1), 'idtable', ssinfo, 'ndim', 1, ...
                'reltype', 3, 'nocc', 2, 'CI', .95};

            testCase.verifyError(@() psyrat_ssrel_dynrel_trt_gamma_ls(common{:}, ...
                'b', repmat([0.1 0.2], nd, 1), 'b_sigma', rp(0), ...
                'ind_sd', zeros(nd,2)), 'varargin:b');
            testCase.verifyError(@() psyrat_ssrel_dynrel_trt_gamma_ls(common{:}, ...
                'b', rp(0), 'b_sigma', repmat([0.1 0.2], nd, 1), ...
                'ind_sd', zeros(nd,2)), 'varargin:bsigma');
            % only one draw column for two participants
            testCase.verifyError(@() psyrat_ssrel_dynrel_trt_gamma_ls(common{:}, ...
                'b', rp(0), 'b_sigma', rp(0), ...
                'ind_sd', zeros(nd,1)), 'varargin:id2');
        end

    end
end
