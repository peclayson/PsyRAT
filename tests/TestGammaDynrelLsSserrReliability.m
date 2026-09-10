classdef TestGammaDynrelLsSserrReliability < matlab.unittest.TestCase
    %TESTGAMMADYNRELLSSSERRRELIABILITY Offline proof of the per-participant
    % CONDITIONAL reliability calculator for the LOCATION-SCALE Gamma /
    % scaled-chi-square subject-level design (analysis 26, family = 'gamma',
    % gammascale = 2): psyrat_ssrel_dynrel_gamma_ls. No CmdStan.
    %
    % The calculator is a standalone file, so these tests drive the PRODUCTION
    % code rather than re-deriving its formula. Constant posterior draws are used
    % throughout: quantiles of a constant collapse to that constant, so every
    % coefficient can be checked against a closed form at machine tolerance.
    %
    % WHAT IS BEING PINNED. For each participant s at their own z_s (the POOLED
    % per-participant residual, S17-FIX):
    %
    %   alpha(z_s)   = alpha0  + b       . x(z_s)
    %   logsig(z_s)  = logsig0 + b_sigma . x(z_s)
    %   sigma_p^2, sigma_i^2  <- psyrat_gamma_varcomps_ls at alpha(z_s)
    %   K            = exp(s_p^2+s_i^2)*(exp(s_p^2)-1)*(exp(s_i^2)-1)
    %   sigma_e^2(s) = exp(2*alpha(z_s))*K + exp(2*(logsig(z_s) + w_s))
    %
    % The closed forms below are written out INDEPENDENTLY (the K identity is
    % section-J algebra, not a call into the converter), so they are an oracle
    % for the estimand rather than a restatement of the implementation.
    %
    % Three things are tested here rather than assumed, because a formula copied
    % from the log-nu twin (TestGammaDynrelSserrReliability) or from the pre-fix
    % conditional-only behavior would still return plausible numbers:
    %
    %  (a) NO mu FACTOR AND NO TRIAL-AVERAGING FACTOR on the CONDITIONAL term.
    %      The twin forms 2*exp(2*(alpha+u_s) + 2*s_i^2 - (lognu+v_s)); here the
    %      conditional variance is read off the scale submodel directly, with no
    %      leading 2, no exp(2*s_i^2) and no alpha. testMatchesClosedForm uses a
    %      nonzero s_i and a nonzero mean slope precisely so each of those stray
    %      factors would break it.
    %  (b) THE POOLED PERSON x TRIAL TERM IS PRESENT AND POPULATION-BUILT
    %      (S17-FIX): exp(2*alpha(z_s))*K rides in every participant's residual,
    %      is common to same-z participants, and vanishes when s_i = 0.
    %      testResidualPoolsThePersonXTrialTerm pins it directly.
    %  (c) THE PARTICIPANT'S OWN LOG-MEAN u_s DOES NOT ENTER AT ALL, so ind_bs
    %      is not an input. testMeanSlopeStillMovesTheCoefficient pins the
    %      related and easily-confused fact that the population mean SLOPE does
    %      still move the coefficient - sigma_p^2 and the pooled term carry
    %      exp(2*alpha(z)) while the conditional term does not, so the
    %      cancellation that holds under gammascale = 1 is gone here.

    methods (Static, Access = private)

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

        function testMatchesClosedForm(testCase)
            % DECISIVE. Every returned coefficient, ICC and SEM must equal the
            % closed form evaluated at that participant's own z with their own
            % w_s. Machine tolerance, because the draws are constant.
            nd = 40;
            alpha0 = log(5.0); b = 0.25;
            logsig0 = log(1.4); b_sigma = -0.30;
            sig_p = 0.40; sig_i = 0.30;   % nonzero: a stray exp(2*s_i^2) would show

            z    = [-1.0; 0.0; 0.8];
            w    = [ 0.35; -0.20; 0.10];
            trls = [ 20;    35;   28];

            ssinfo = TestGammaDynrelLsSserrReliability.ssinfo1((1:3)', z, trls);
            rp = @(x) TestGammaDynrelLsSserrReliability.rep(x, nd);

            [gen, dep] = psyrat_ssrel_dynrel_gamma_ls( ...
                'alpha0', rp(alpha0), 'b', rp(b), ...
                'logsig0', rp(logsig0), 'b_sigma', rp(b_sigma), ...
                'sig_p', rp(sig_p), 'sig_i', rp(sig_i), ...
                'ind_sd', repmat(w(:)', nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'CI', .95);

            testCase.assertEqual(height(gen), 3);
            testCase.assertEqual(height(dep), 3);

            for k = 1:3
                alpha_z  = alpha0  + b       * z(k);
                logsig_z = logsig0 + b_sigma * z(k);
                vc = psyrat_gamma_varcomps_ls(alpha_z, sig_p, sig_i, logsig_z, 0);
                % POOLED residual (S17-FIX): the mean-surface person x trial
                % term exp(2*alpha(z))*K plus the participant's own conditional
                % variance. K is written out from the section-J identity, NOT
                % read off the converter, so this line is an independent oracle.
                K = exp(sig_p^2 + sig_i^2) * (exp(sig_p^2) - 1) * (exp(sig_i^2) - 1);
                sigE2 = exp(2*alpha_z) * K + exp(2*(logsig_z + w(k)));
                n = trls(k);

                G_ref    = vc.sigma_p2 / (vc.sigma_p2 + sigE2 / n);
                D_ref    = vc.sigma_p2 / (vc.sigma_p2 + (sigE2 + vc.sigma_i2) / n);
                ICCg_ref = vc.sigma_p2 / (vc.sigma_p2 + sigE2);
                ICCd_ref = vc.sigma_p2 / (vc.sigma_p2 + sigE2 + vc.sigma_i2);

                testCase.verifyEqual(gen.dep_pt(k), G_ref,    'RelTol', 1e-12);
                testCase.verifyEqual(dep.dep_pt(k), D_ref,    'RelTol', 1e-12);
                testCase.verifyEqual(gen.icc_pt(k), ICCg_ref, 'RelTol', 1e-12);
                testCase.verifyEqual(dep.icc_pt(k), ICCd_ref, 'RelTol', 1e-12);
                testCase.verifyEqual(gen.sem_pt(k), sqrt(sigE2 / n), 'RelTol', 1e-12);
                testCase.verifyEqual(dep.sem_pt(k), ...
                    sqrt((sigE2 + vc.sigma_i2) / n), 'RelTol', 1e-12);

                % dependability (absolute error) never exceeds generalizability
                testCase.verifyLessThanOrEqual(dep.dep_pt(k), gen.dep_pt(k));
            end
        end

        function testConditionalTermComesOnlyFromTheScaleSubmodel(testCase)
            % The participant-specific PART of the residual must be
            % exp(2*(logsig(z) + w_s)) and nothing else. Two participants at the
            % SAME z share the pooled person x trial term exactly (it is built
            % from the POPULATION alpha(z), s_p and s_i), so it cancels from the
            % DIFFERENCE of their residual variances - recoverable from the
            % reported SEMs, since sem^2 * n is the residual variance on the
            % generalizability side. That difference must equal
            % exp(2*logsig(z)) * (exp(2*w_2) - exp(2*w_1)) exactly.
            %
            % SCOPE. Being a difference between two participants, this is blind
            % to any term COMMON to both - which is precisely what makes it the
            % right pin post-S17-FIX: it stays green because the pooled term is
            % participant-invariant at fixed z, and it turns red if the
            % conditional term picks up any mu- or trial-dependent factor, or if
            % the pooled term is ever made to depend on the participant's own
            % effects (the difference would then no longer isolate the
            % conditional part). testResidualPoolsThePersonXTrialTerm pins the
            % common term itself.
            nd = 20;
            d = 0.5;
            z    = [0.4; 0.4];
            w    = [0.0; d];
            trls = [30;  30];
            logsig0 = log(1.4); b_sigma = -0.30;

            ssinfo = TestGammaDynrelLsSserrReliability.ssinfo1((1:2)', z, trls);
            rp = @(x) TestGammaDynrelLsSserrReliability.rep(x, nd);

            [gen, ~] = psyrat_ssrel_dynrel_gamma_ls( ...
                'alpha0', rp(log(5.0)), 'b', rp(0.25), ...
                'logsig0', rp(logsig0), 'b_sigma', rp(b_sigma), ...
                'sig_p', rp(0.40), 'sig_i', rp(0.30), ...
                'ind_sd', repmat(w(:)', nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'CI', .95);

            logsig_z = logsig0 + b_sigma * z(1);
            diffref = exp(2*logsig_z) * (exp(2*d) - 1);
            diffobs = (gen.sem_pt(2)^2 - gen.sem_pt(1)^2) * trls(1);
            testCase.verifyEqual(diffobs, diffref, 'RelTol', 1e-12, ...
                ['The conditional part of the residual must scale as ' ...
                'exp(2*w_s), and the pooled person x trial term must cancel ' ...
                'from a same-z between-participant difference.']);
            % and the less-precise participant must get the lower coefficient
            testCase.verifyLessThan(gen.dep_pt(2), gen.dep_pt(1));
        end

        function testResidualPoolsThePersonXTrialTerm(testCase)
            % THE S17-FIX PIN. The residual recovered from the reported SEM
            % (sem^2 * n on the generalizability side) minus the participant's
            % conditional variance must equal the mean-surface person x trial
            % term exp(2*alpha(z))*K exactly - written out independently, not
            % read off the converter. And the term must VANISH at s_i = 0, where
            % the mean surface has no trial facet to interact with, collapsing
            % the residual to the conditional variance alone (the pre-fix
            % behavior, which is only correct in that degenerate case).
            nd = 20;
            alpha0 = log(5.0); b = 0.25;
            logsig0 = log(1.4); b_sigma = -0.30;
            sig_p = 0.40; sig_i = 0.30;
            z = 0.6; w = 0.15; n = 30;

            ssinfo = TestGammaDynrelLsSserrReliability.ssinfo1(1, z, n);
            rp = @(x) TestGammaDynrelLsSserrReliability.rep(x, nd);
            runcalc = @(si) psyrat_ssrel_dynrel_gamma_ls( ...
                'alpha0', rp(alpha0), 'b', rp(b), ...
                'logsig0', rp(logsig0), 'b_sigma', rp(b_sigma), ...
                'sig_p', rp(sig_p), 'sig_i', rp(si), ...
                'ind_sd', repmat(w, nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'CI', .95);

            alpha_z  = alpha0 + b * z;
            logsig_z = logsig0 + b_sigma * z;
            condvar  = exp(2*(logsig_z + w));

            gen = runcalc(sig_i);
            K = exp(sig_p^2 + sig_i^2) * (exp(sig_p^2) - 1) * (exp(sig_i^2) - 1);
            pooled = gen.sem_pt(1)^2 * n - condvar;
            testCase.verifyEqual(pooled, exp(2*alpha_z) * K, 'RelTol', 1e-10, ...
                ['The residual must pool the population person x trial term ' ...
                'exp(2*alpha(z))*K with the conditional variance (S17-FIX).']);

            gen0 = runcalc(0);
            testCase.verifyEqual(gen0.sem_pt(1)^2 * n, condvar, 'RelTol', 1e-12, ...
                ['At s_i = 0 the person x trial term must vanish and the ' ...
                'residual collapse to the conditional variance alone.']);
        end

        function testMeanSlopeStillMovesTheCoefficient(testCase)
            % PINS BOTH the "own-z, all components" estimand and the fact that the
            % log-mean slope does NOT cancel under gammascale = 2.
            %
            % With b_sigma = 0 the CONDITIONAL variance is z-invariant, so two
            % participants with identical w_s and identical trial counts but
            % different z differ only through exp(2*alpha(z)): in sigma_p^2, in
            % sigma_i^2 and (post-S17-FIX) in the pooled person x trial term.
            % Writing the generalizability coefficient as
            %   G = A*e^{2a} / (A*e^{2a} + (K*e^{2a} + c)/n),   a = alpha(z),
            % dividing through by e^{2a} gives A / (A + K/n + c*e^{-2a}/n): only
            % the conditional term c is deflated as z rises, so the participant
            % at the higher z - larger expected score against the same absolute
            % noise - must still score higher.
            %
            % This fails two ways a future change could go wrong: reverting to a
            % shared reference z for the mean-surface components (the cheaper
            % alternative rejected for the log-nu twin), or re-coupling the
            % conditional term to the mean so the exp(2*alpha) factor cancels
            % out of the coefficient entirely.
            nd = 20;
            b = 0.30; b_sigma = 0;   % => conditional variance identical at every z
            alpha0 = log(5.0); logsig0 = log(1.4);
            sig_p = 0.40; sig_i = 0.30;

            z    = [-1.0; 1.0];
            w    = [ 0.20; 0.20];          % identical persons
            trls = [ 25;   25];            % identical trial counts

            ssinfo = TestGammaDynrelLsSserrReliability.ssinfo1((1:2)', z, trls);
            rp = @(x) TestGammaDynrelLsSserrReliability.rep(x, nd);

            [gen, ~] = psyrat_ssrel_dynrel_gamma_ls( ...
                'alpha0', rp(alpha0), 'b', rp(b), ...
                'logsig0', rp(logsig0), 'b_sigma', rp(b_sigma), ...
                'sig_p', rp(sig_p), 'sig_i', rp(sig_i), ...
                'ind_sd', repmat(w(:)', nd, 1), ...
                'idtable', ssinfo, 'ndim', 1, 'CI', .95);

            % sanity: the residuals differ ONLY by the pooled person x trial
            % term, which carries exp(2*alpha(z)) - the conditional parts are
            % identical by construction (b_sigma = 0, equal w)
            K = exp(sig_p^2 + sig_i^2) * (exp(sig_p^2) - 1) * (exp(sig_i^2) - 1);
            semdiff = (gen.sem_pt(2)^2 - gen.sem_pt(1)^2) * trls(1);
            pooldiff = (exp(2*(alpha0 + b)) - exp(2*(alpha0 - b))) * K;
            testCase.assertEqual(semdiff, pooldiff, 'RelTol', 1e-10, ...
                ['Test setup error: with b_sigma = 0 and equal w, the ' ...
                'residuals should differ only by the pooled term.']);

            testCase.verifyGreaterThan(gen.dep_pt(2), gen.dep_pt(1), ...
                ['The participant at the higher z must score higher: only the ' ...
                'conditional term fails to carry exp(2*alpha(z)), so the mean ' ...
                'slope does not cancel out of the coefficient.']);
            testCase.verifyGreaterThan(gen.dep_pt(2) - gen.dep_pt(1), 1e-3, ...
                'The z effect is too small to distinguish from round-off.');
        end

        function testTwoDimensionUsesTheInteractionTerm(testCase)
            % ndim = 2 must build the design row [z1; z2; z1*z2]. A participant at
            % (1,1) and one at (1,-1) share z1 but differ in the z2 and
            % interaction terms, so with a nonzero third slope they must differ.
            nd = 20;
            bvec = [0.20 0.10 0.40];       % third entry is the interaction
            z1 = [1.0; 1.0]; z2 = [1.0; -1.0];

            t = table((1:2)', (1:2)', z1, z2, [25;25], ...
                'VariableNames', {'id','id2','z1','z2','trls'});
            rp = @(x) TestGammaDynrelLsSserrReliability.rep(x, nd);

            [gen, ~] = psyrat_ssrel_dynrel_gamma_ls( ...
                'alpha0', rp(log(5.0)), 'b', repmat(bvec, nd, 1), ...
                'logsig0', rp(log(1.4)), 'b_sigma', zeros(nd,3), ...
                'sig_p', rp(0.40), 'sig_i', rp(0.30), ...
                'ind_sd', repmat([0.2 0.2], nd, 1), ...
                'idtable', t, 'ndim', 2, 'CI', .95);

            alphaA = log(5.0) + bvec * [1; 1; 1];
            alphaB = log(5.0) + bvec * [1; -1; -1];
            testCase.assertNotEqual(alphaA, alphaB, ...
                'Test setup error: the two design rows should differ.');
            testCase.verifyNotEqual(gen.dep_pt(1), gen.dep_pt(2), ...
                'The two-dimension design row must include the interaction term.');
        end

        function testGuardsRejectMisshapedInput(testCase)
            % Loud failure rather than a silently plausible number.
            nd = 20;
            ssinfo = TestGammaDynrelLsSserrReliability.ssinfo1((1:2)', [0;1], [25;25]);
            rp = @(x) TestGammaDynrelLsSserrReliability.rep(x, nd);
            base = {'alpha0', rp(log(5.0)), 'logsig0', rp(log(1.4)), ...
                'sig_p', rp(0.40), 'sig_i', rp(0.30), 'ndim', 1, 'CI', .95};

            % slope block with the wrong column count
            testCase.verifyError(@() psyrat_ssrel_dynrel_gamma_ls(base{:}, ...
                'b', zeros(nd,3), 'b_sigma', zeros(nd,1), ...
                'ind_sd', zeros(nd,2), 'idtable', ssinfo), 'varargin:b');
            testCase.verifyError(@() psyrat_ssrel_dynrel_gamma_ls(base{:}, ...
                'b', zeros(nd,1), 'b_sigma', zeros(nd,3), ...
                'ind_sd', zeros(nd,2), 'idtable', ssinfo), 'varargin:bsigma');

            % id2 pointing past the available per-subject draw columns: without
            % this guard the participant would be silently assigned another
            % person's residual, or the index would fail far downstream
            testCase.verifyError(@() psyrat_ssrel_dynrel_gamma_ls(base{:}, ...
                'b', zeros(nd,1), 'b_sigma', zeros(nd,1), ...
                'ind_sd', zeros(nd,1), 'idtable', ssinfo), 'varargin:id2');
        end

    end
end
