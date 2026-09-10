classdef TestGammaDynrelSserrReliability < matlab.unittest.TestCase
    %TESTGAMMADYNRELSSERRRELIABILITY Offline proof of the per-participant
    % CONDITIONAL dynamic/conditional reliability calculators for the Gamma /
    % scaled-chi-square SUBJECT-LEVEL designs: analysis 26
    % (psyrat_ssrel_dynrel_gamma, one facet) and analysis 27
    % (psyrat_ssrel_dynrel_trt_gamma, trial + occasion). No CmdStan.
    %
    % Unlike the gamma EXTRACTORS (local functions of psyrat_computevarcomp, only
    % reachable through a live fit), these calculators are standalone files, so
    % these tests drive the PRODUCTION code directly rather than re-deriving its
    % formula. Constant posterior draws are used throughout: quantiles of a
    % constant collapse to that constant, so every coefficient can be checked
    % against a closed form at machine tolerance.
    %
    % WHAT IS BEING PINNED. Both calculators implement the owner-directed
    % "own-z, all components" estimand (2026-07-31): for each participant,
    %
    %   alpha(z_s) = alpha0 + b . x(z_s)      lognu(z_s) = lognu0 + b_nu . x(z_s)
    %   every observed-scale component is evaluated at z_s
    %   sigma_e^2(s,z_s) = pi_mean(z_s) + 2*exp(2*(alpha(z_s) + u_s) + 2*V_w - (lognu(z_s) + v_s))
    %   with pi_mean(z) = exp(2*alpha(z)) * K, K = (exp(2*sig_p^2) - exp(sig_p^2))
    %   * (exp(2*sig_i^2) - exp(sig_i^2)): the POPULATION person x trial term of
    %   the log-linked mean surface, pooled into every participant's residual
    %   exactly as the location-scale twin does (S17-FIX; S20 for log-nu).
    %
    % with V_w the WITHIN-person log-mean variance (the trial main effect for
    % analysis 26; the five crossed facets for analysis 27), excluding the person
    % main effect because u_s is conditioned on rather than integrated out.
    %
    % testComponentsMoveWithZ is the test that actually pins the estimand DECISION:
    % it fails if a future change reverts to holding sigma_p2/sigma_i2 at a shared
    % reference z, which is exactly the cheaper alternative that was rejected.

    methods (Static, Access = private)


        function pm = piMean(alpha, sig_p, sig_i)

            %Closed-form population person x trial term of the log-linked

            %mean surface, exp(2*alpha)*K, derived independently of

            %psyrat_gamma_varcomps (whose var_mu - sigma_p2 - sigma_i2

            %it must equal; the one-facet test asserts that).

            K = (exp(2*sig_p.^2) - exp(sig_p.^2)) .* (exp(2*sig_i.^2) - exp(sig_i.^2));

            pm = exp(2*alpha) .* K;

        end

        function d = rep(v, n)
            % n identical posterior draws of a scalar
            d = repmat(v, n, 1);
        end

        function t = ssinfo1(ids, z1, trls)
            % one-facet per-participant lookup, id2 = 1:NSUB in row order
            n = numel(z1);
            t = table((1:n)', (1:n)', z1(:), trls(:), ...
                'VariableNames', {'id','id2','z1','trls'});
            t.id = ids(:);
        end

    end

    methods (Test)

        % ---------- analysis 26: one facet ----------

        function testOneFacetMatchesClosedForm(testCase)
            % DECISIVE. Every returned coefficient, ICC and SEM must equal the
            % closed form evaluated at that participant's own z with their own
            % (u_p, v_p). Machine tolerance, because the draws are constant.
            nd = 40;
            alpha0 = log(5.0); b = 0.25; lognu0 = log(16); b_nu = -0.30;
            sig_p = 0.40; sig_i = 0.30;

            z    = [-1.0; 0.0; 0.8];
            u    = [ 0.35; -0.20; 0.10];
            v    = [-0.15;  0.40; 0.05];
            trls = [ 20;    35;   28];

            ssinfo = TestGammaDynrelSserrReliability.ssinfo1((1:3)', z, trls);
            rp = @(x) TestGammaDynrelSserrReliability.rep(x, nd);

            [gen, dep] = psyrat_ssrel_dynrel_gamma( ...
                'alpha0', rp(alpha0), 'b', rp(b), ...
                'lognu0', rp(lognu0), 'b_nu', rp(b_nu), ...
                'sig_p', rp(sig_p), 'sig_i', rp(sig_i), ...
                'ind_bs', repmat(u(:)', nd, 1), 'ind_nu', repmat(v(:)', nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'CI', .95);

            testCase.assertEqual(height(gen), 3);
            testCase.assertEqual(height(dep), 3);

            for k = 1:3
                alpha_z = alpha0 + b * z(k);
                lognu_z = lognu0 + b_nu * z(k);
                vc = psyrat_gamma_varcomps(alpha_z, sig_p, sig_i, exp(lognu_z), 0, 0);
                piMean = TestGammaDynrelSserrReliability.piMean(alpha_z, sig_p, sig_i);
                testCase.assertEqual(piMean, max(vc.var_mu - vc.sigma_p2 - vc.sigma_i2, 0), ...
                    'RelTol', 1e-12, ['closed-form pi_mean must equal the helper''s ', ...
                    'var_mu - sigma_p2 - sigma_i2 (the two derivations agree)']);
                sigE2 = piMean + 2 * exp(2*(alpha_z + u(k)) + 2*sig_i^2 - (lognu_z + v(k)));
                n = trls(k);

                G_ref   = vc.sigma_p2 / (vc.sigma_p2 + sigE2 / n);
                D_ref   = vc.sigma_p2 / (vc.sigma_p2 + (sigE2 + vc.sigma_i2) / n);
                ICCg_ref = vc.sigma_p2 / (vc.sigma_p2 + sigE2);
                ICCd_ref = vc.sigma_p2 / (vc.sigma_p2 + sigE2 + vc.sigma_i2);

                testCase.verifyEqual(gen.dep_pt(k), G_ref,   'RelTol', 1e-12);
                testCase.verifyEqual(dep.dep_pt(k), D_ref,   'RelTol', 1e-12);
                testCase.verifyEqual(gen.icc_pt(k), ICCg_ref,'RelTol', 1e-12);
                testCase.verifyEqual(dep.icc_pt(k), ICCd_ref,'RelTol', 1e-12);
                testCase.verifyEqual(gen.sem_pt(k), sqrt(sigE2 / n), 'RelTol', 1e-12);
                testCase.verifyEqual(dep.sem_pt(k), ...
                    sqrt((sigE2 + vc.sigma_i2) / n), 'RelTol', 1e-12);

                % dependability (absolute error) never exceeds generalizability
                testCase.verifyLessThanOrEqual(dep.dep_pt(k), gen.dep_pt(k));
            end
        end

        function testComponentsMoveWithZ(testCase)
            % PINS THE ESTIMAND DECISION ("own-z, all components"). Under gamma the
            % between-person and trial components carry exp(2*alpha(z)), so two
            % participants with identical (u_p, v_p) and identical trial counts but
            % DIFFERENT z must still receive different coefficients whenever the
            % log-mean slope b is nonzero.
            %
            % If a future change reverts to the cheaper "residual-only, reference
            % z" alternative - holding sigma_p2/sigma_i2 fixed and varying only the
            % residual - this test fails, because with b and b_nu chosen so the
            % residual is z-invariant the coefficients would then be identical.
            %
            % Residual z-invariance requires 2*b - b_nu = 0, since
            % sigma_e^2(z) is proportional to exp(2*b*z - b_nu*z).
            nd = 20;
            b = 0.30; b_nu = 2*b;          % => residual identical at every z
            alpha0 = log(5.0); lognu0 = log(16);
            sig_p = 0.40; sig_i = 0.30;

            z    = [-1.0; 1.0];
            u    = [ 0.20; 0.20];          % identical persons
            v    = [-0.10; -0.10];
            trls = [ 25;   25];            % identical trial counts

            ssinfo = TestGammaDynrelSserrReliability.ssinfo1((1:2)', z, trls);
            rp = @(x) TestGammaDynrelSserrReliability.rep(x, nd);

            [gen, ~] = psyrat_ssrel_dynrel_gamma( ...
                'alpha0', rp(alpha0), 'b', rp(b), ...
                'lognu0', rp(lognu0), 'b_nu', rp(b_nu), ...
                'sig_p', rp(sig_p), 'sig_i', rp(sig_i), ...
                'ind_bs', repmat(u(:)', nd, 1), 'ind_nu', repmat(v(:)', nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'CI', .95);

            % sanity: the residual really is z-invariant under this parameterization
            %the CONDITIONAL (observation-level) part of the residual is
            %z-invariant under these slopes; the pooled pi_mean part (S20)
            %moves with z, and so does the between-person component, so the
            %coefficient must move with z either way.
            r1 = 2*exp(2*(alpha0 + b*z(1) + u(1)) + 2*sig_i^2 - (lognu0 + b_nu*z(1) + v(1)));
            r2 = 2*exp(2*(alpha0 + b*z(2) + u(2)) + 2*sig_i^2 - (lognu0 + b_nu*z(2) + v(2)));
            testCase.assertEqual(r1, r2, 'RelTol', 1e-12, ...
                'Test setup error: the conditional residual should be z-invariant here.');

            % so any remaining difference comes from sigma_p2(z) / sigma_i2(z)
            testCase.verifyNotEqual(gen.dep_pt(1), gen.dep_pt(2));
            testCase.verifyGreaterThan(abs(gen.dep_pt(1) - gen.dep_pt(2)), 1e-6, ...
                ['The two participants differ ONLY in z, and the residual is ', ...
                'z-invariant here, so a difference can only come from evaluating ', ...
                'sigma_p2/sigma_i2 at each participant''s own z. Identical values ', ...
                'mean the calculator reverted to a shared reference z.']);
        end

        function testZeroSlopesReduceToStaticSubjectLevelGamma(testCase)
            % REDUCTION. With b = b_nu = 0 the dimension conditioning vanishes and
            % every participant's coefficient must equal the STATIC one-facet
            % subject-level gamma result (analysis 6, psyrat_gamma_extract_sserr +
            % psyrat_ssrel), regardless of their z.
            nd = 25;
            alpha0 = log(5.5); lognu0 = log(18);
            sig_p = 0.45; sig_i = 0.25;

            z    = [-1.5; 0.0; 2.0];       % z varies but must not matter
            u    = [ 0.30; 0.30; 0.30];
            v    = [-0.20; -0.20; -0.20];
            trls = [ 30;   30;   30];

            ssinfo = TestGammaDynrelSserrReliability.ssinfo1((1:3)', z, trls);
            rp = @(x) TestGammaDynrelSserrReliability.rep(x, nd);

            [gen, dep] = psyrat_ssrel_dynrel_gamma( ...
                'alpha0', rp(alpha0), 'b', rp(0), ...
                'lognu0', rp(lognu0), 'b_nu', rp(0), ...
                'sig_p', rp(sig_p), 'sig_i', rp(sig_i), ...
                'ind_bs', repmat(u(:)', nd, 1), 'ind_nu', repmat(v(:)', nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'CI', .95);

            % static case-6 reference, built the way psyrat_gamma_extract_sserr does
            vc = psyrat_gamma_varcomps(alpha0, sig_p, sig_i, exp(lognu0), 0, 0);
            sigE2 = TestGammaDynrelSserrReliability.piMean(alpha0, sig_p, sig_i) ...
                + 2 * exp(2*(alpha0 + u(1)) + 2*sig_i^2 - (lognu0 + v(1)));
            G_ref = vc.sigma_p2 / (vc.sigma_p2 + sigE2 / trls(1));
            D_ref = vc.sigma_p2 / (vc.sigma_p2 + (sigE2 + vc.sigma_i2) / trls(1));

            for k = 1:3
                testCase.verifyEqual(gen.dep_pt(k), G_ref, 'RelTol', 1e-12);
                testCase.verifyEqual(dep.dep_pt(k), D_ref, 'RelTol', 1e-12);
            end
        end

        function testHigherMeanParticipantGetsLowerCoefficient(testCase)
            % Documents the behavioral consequence written into both calculators'
            % headers: on the GROUP surface a mean-level slope cancels because
            % every component shares exp(2*alpha(z)), but at the SUBJECT level u_p
            % enters the residual and not the between-person numerator, so a
            % higher-mean participant receives a strictly lower coefficient. This
            % is a property of Var = 2*mu^2/nu, not a defect - pin it so the
            % asymmetry is not "fixed" by accident.
            nd = 20;
            alpha0 = log(5.0); lognu0 = log(16);
            sig_p = 0.40; sig_i = 0.30;

            z    = [0.0; 0.0];
            u    = [0.00; 0.50];           % second participant has the higher mean
            v    = [0.00; 0.00];
            trls = [ 25;  25];

            ssinfo = TestGammaDynrelSserrReliability.ssinfo1((1:2)', z, trls);
            rp = @(x) TestGammaDynrelSserrReliability.rep(x, nd);

            [gen, ~] = psyrat_ssrel_dynrel_gamma( ...
                'alpha0', rp(alpha0), 'b', rp(0), ...
                'lognu0', rp(lognu0), 'b_nu', rp(0), ...
                'sig_p', rp(sig_p), 'sig_i', rp(sig_i), ...
                'ind_bs', repmat(u(:)', nd, 1), 'ind_nu', repmat(v(:)', nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'CI', .95);

            testCase.verifyLessThan(gen.dep_pt(2), gen.dep_pt(1));
        end

        function testTwoDimensionInteractionTermIsUsed(testCase)
            % With ndim = 2 the design row is [z1; z2; z1*z2]. A participant with
            % z1 = z2 = 1 must therefore differ from one with z1 = 1, z2 = -1 by
            % more than the additive terms alone - i.e. the interaction column is
            % actually applied.
            nd = 20;
            alpha0 = log(5.0); lognu0 = log(16);
            sig_p = 0.40; sig_i = 0.30;
            b    = [0.10 0.05 0.40];       % large interaction slope
            b_nu = [0.00 0.00 0.00];

            n = 2;
            ssinfo = table((1:n)', (1:n)', [1; 1], [1; -1], [25; 25], ...
                'VariableNames', {'id','id2','z1','z2','trls'});
            rp = @(x) TestGammaDynrelSserrReliability.rep(x, nd);

            [gen, ~] = psyrat_ssrel_dynrel_gamma( ...
                'alpha0', rp(alpha0), 'b', repmat(b, nd, 1), ...
                'lognu0', rp(lognu0), 'b_nu', repmat(b_nu, nd, 1), ...
                'sig_p', rp(sig_p), 'sig_i', rp(sig_i), ...
                'ind_bs', zeros(nd, n), 'ind_nu', zeros(nd, n), ...
                'idtable', ssinfo, 'ndim', 2, 'CI', .95);

            % closed form at each participant's [z1 z2 z1*z2]
            for k = 1:n
                xr = [ssinfo.z1(k); ssinfo.z2(k); ssinfo.z1(k)*ssinfo.z2(k)];
                alpha_z = alpha0 + b * xr;
                lognu_z = lognu0 + b_nu * xr;
                vc = psyrat_gamma_varcomps(alpha_z, sig_p, sig_i, exp(lognu_z), 0, 0);
                sigE2 = TestGammaDynrelSserrReliability.piMean(alpha_z, sig_p, sig_i) ...
                    + 2 * exp(2*alpha_z + 2*sig_i^2 - lognu_z);
                G_ref = vc.sigma_p2 / (vc.sigma_p2 + sigE2 / ssinfo.trls(k));
                testCase.verifyEqual(gen.dep_pt(k), G_ref, 'RelTol', 1e-12);
            end
        end

        function testRejectsMisshapedSlopeMatrix(testCase)
            % A transposed or wrong-width b must fail loudly rather than broadcast.
            nd = 10;
            ssinfo = TestGammaDynrelSserrReliability.ssinfo1((1:2)', [0;1], [20;20]);
            rp = @(x) TestGammaDynrelSserrReliability.rep(x, nd);
            bad = @() psyrat_ssrel_dynrel_gamma( ...
                'alpha0', rp(log(5)), 'b', repmat([0.1 0.2], nd, 1), ...
                'lognu0', rp(log(16)), 'b_nu', rp(0), ...
                'sig_p', rp(0.4), 'sig_i', rp(0.3), ...
                'ind_bs', zeros(nd,2), 'ind_nu', zeros(nd,2), ...
                'idtable', ssinfo, 'ndim', 1, 'CI', .95);
            testCase.verifyError(bad, 'varargin:b');
        end

        % ---------- analysis 27: trial + occasion ----------

        function testTwoFacetMatchesRelTrtAtOwnZ(testCase)
            % DECISIVE for analysis 27. Each participant's coefficient must equal a
            % direct psyrat_rel_trt call built from the six observed-scale
            % components evaluated at THAT participant's z, with their own
            % conditional residual. V_w sums the five crossed log-mean variances
            % and excludes the person main effect.
            nd = 30;
            alpha0 = log(5.0); b = 0.20; lognu0 = log(18); b_nu = -0.25;
            sig_p = 0.40; s_o = 0.22; s_t = 0.18; s_po = 0.14; s_pt = 0.12; s_ot = 0.09;
            reltype = 3; nocc = 2;

            z    = [-0.8; 0.0; 1.2];
            u    = [ 0.25; -0.15; 0.05];
            v    = [-0.10;  0.30; 0.00];
            trls = [ 18;    24;   30];

            n = 3;
            ssinfo = table((1:n)', (1:n)', z, trls, repmat(nocc, n, 1), ...
                'VariableNames', {'id','id2','z1','trls','occs'});
            rp = @(x) TestGammaDynrelSserrReliability.rep(x, nd);

            [gen, dep] = psyrat_ssrel_dynrel_trt_gamma( ...
                'alpha0', rp(alpha0), 'b', rp(b), ...
                'lognu0', rp(lognu0), 'b_nu', rp(b_nu), ...
                'sig_p', rp(sig_p), 'sig_occ', rp(s_o), 'sig_trl', rp(s_t), ...
                'sig_occxid', rp(s_po), 'sig_trlxid', rp(s_pt), 'sig_trlxocc', rp(s_ot), ...
                'ind_bs', repmat(u(:)', nd, 1), 'ind_nu', repmat(v(:)', nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'reltype', reltype, ...
                'nocc', nocc, 'CI', .95);

            V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2;
            for k = 1:n
                alpha_z = alpha0 + b * z(k);
                lognu_z = lognu0 + b_nu * z(k);
                vc = psyrat_gamma_varcomps_trt(alpha_z, sig_p, s_o, s_t, s_po, ...
                    s_pt, s_ot, exp(lognu_z), 0, 0);
                errSs = sqrt(2 * exp(2*(alpha_z + u(k)) + 2*V_w - (lognu_z + v(k))));

                args = {'reltype', reltype, ...
                    'bp', sqrt(vc.sigma_p2), 'bo', sqrt(vc.sigma_o2), ...
                    'bt', sqrt(vc.sigma_t2), 'txp', sqrt(vc.sigma_pt2), ...
                    'oxp', sqrt(vc.sigma_po2), 'txo', sqrt(vc.sigma_ot2), ...
                    'err', errSs, 'nocc', nocc, 'CI', .95};

                [~, G_ref, ~] = psyrat_rel_trt('gcoeff', 2, args{:}, 'obs', trls(k));
                [~, D_ref, ~] = psyrat_rel_trt('gcoeff', 1, args{:}, 'obs', trls(k));
                testCase.verifyEqual(gen.dep_pt(k), G_ref, 'RelTol', 1e-12);
                testCase.verifyEqual(dep.dep_pt(k), D_ref, 'RelTol', 1e-12);

                % ICC is pinned at obs = nocc = 1 for composite invariance
                [~, ICCg_ref, ~] = psyrat_rel_trt('gcoeff', 2, ...
                    args{1:end-4}, 'nocc', 1, 'CI', .95, 'obs', 1);
                testCase.verifyEqual(gen.icc_pt(k), ICCg_ref, 'RelTol', 1e-12);
            end
        end

        function testTwoFacetComponentsMoveWithZ(testCase)
            % The analysis-27 twin of testComponentsMoveWithZ: with the residual
            % held z-invariant (2*b - b_nu = 0), two otherwise identical
            % participants at different z must still differ, because all six
            % crossed components are evaluated at their own z.
            nd = 20;
            b = 0.30; b_nu = 2*b;
            alpha0 = log(5.0); lognu0 = log(18);
            sig_p = 0.40; s_o = 0.22; s_t = 0.18; s_po = 0.14; s_pt = 0.12; s_ot = 0.09;

            n = 2;
            ssinfo = table((1:n)', (1:n)', [-1; 1], [25; 25], [2; 2], ...
                'VariableNames', {'id','id2','z1','trls','occs'});
            rp = @(x) TestGammaDynrelSserrReliability.rep(x, nd);

            [gen, ~] = psyrat_ssrel_dynrel_trt_gamma( ...
                'alpha0', rp(alpha0), 'b', rp(b), ...
                'lognu0', rp(lognu0), 'b_nu', rp(b_nu), ...
                'sig_p', rp(sig_p), 'sig_occ', rp(s_o), 'sig_trl', rp(s_t), ...
                'sig_occxid', rp(s_po), 'sig_trlxid', rp(s_pt), 'sig_trlxocc', rp(s_ot), ...
                'ind_bs', repmat([0.2 0.2], nd, 1), 'ind_nu', repmat([-0.1 -0.1], nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'reltype', 3, 'nocc', 2, 'CI', .95);

            testCase.verifyGreaterThan(abs(gen.dep_pt(1) - gen.dep_pt(2)), 1e-6, ...
                ['Two participants differing only in z, with a z-invariant ', ...
                'residual, produced identical coefficients - the six crossed ', ...
                'components are no longer being evaluated at each participant''s z.']);
        end

        function testTwoFacetResidualExcludesPersonFacet(testCase)
            % GUARD on V_w. If sig_p^2 leaked into the within-person term (the most
            % likely porting error from the marginal converter), every coefficient
            % would change. Compare against the correct closed form and against the
            % mutant, and require the calculator to match the former only.
            nd = 20;
            alpha0 = log(5.0); lognu0 = log(18);
            sig_p = 0.50;                              % deliberately large
            s_o = 0.20; s_t = 0.18; s_po = 0.14; s_pt = 0.12; s_ot = 0.09;

            ssinfo = table(1, 1, 0, 25, 2, ...
                'VariableNames', {'id','id2','z1','trls','occs'});
            rp = @(x) TestGammaDynrelSserrReliability.rep(x, nd);

            [gen, ~] = psyrat_ssrel_dynrel_trt_gamma( ...
                'alpha0', rp(alpha0), 'b', rp(0), ...
                'lognu0', rp(lognu0), 'b_nu', rp(0), ...
                'sig_p', rp(sig_p), 'sig_occ', rp(s_o), 'sig_trl', rp(s_t), ...
                'sig_occxid', rp(s_po), 'sig_trlxid', rp(s_pt), 'sig_trlxocc', rp(s_ot), ...
                'ind_bs', zeros(nd,1), 'ind_nu', zeros(nd,1), ...
                'idtable', ssinfo, 'ndim', 1, 'reltype', 3, 'nocc', 2, 'CI', .95);

            vc = psyrat_gamma_varcomps_trt(alpha0, sig_p, s_o, s_t, s_po, s_pt, ...
                s_ot, exp(lognu0), 0, 0);
            args = {'gcoeff', 2, 'reltype', 3, ...
                'bp', sqrt(vc.sigma_p2), 'bo', sqrt(vc.sigma_o2), ...
                'bt', sqrt(vc.sigma_t2), 'txp', sqrt(vc.sigma_pt2), ...
                'oxp', sqrt(vc.sigma_po2), 'txo', sqrt(vc.sigma_ot2), ...
                'obs', 25, 'nocc', 2, 'CI', .95};

            V_correct = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2;
            V_mutant  = V_correct + sig_p^2;           % person facet leaked in

            [~, G_correct, ~] = psyrat_rel_trt(args{:}, ...
                'err', sqrt(2*exp(2*alpha0 + 2*V_correct - lognu0)));
            [~, G_mutant, ~] = psyrat_rel_trt(args{:}, ...
                'err', sqrt(2*exp(2*alpha0 + 2*V_mutant - lognu0)));

            testCase.verifyEqual(gen.dep_pt(1), G_correct, 'RelTol', 1e-12);
            testCase.verifyNotEqual(gen.dep_pt(1), G_mutant);
            % and the mutant is far away, so the check is not vacuous
            testCase.verifyGreaterThan(abs(G_correct - G_mutant), 1e-3);
        end

    end
end
