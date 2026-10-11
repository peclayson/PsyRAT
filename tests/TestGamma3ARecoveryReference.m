classdef TestGamma3ARecoveryReference < PsyRATTestBase
    % Deterministic pins for tests/helpers/PsyRAT3ARecoveryReference.m, the
    % independent planted-truth reference used by the bounded 3A static
    % subject-level gamma recovery cells (TestGammaSserrRecovery,
    % TestGammaDiffSserrRecovery, TestGammaDiffCopulaSserrRecovery). No CmdStan,
    % no random numbers, no production routine.
    %
    % SCORE AND ESTIMAND UNDER TEST. Single-trial gamma observations with a
    % log-linked mean log mu_pi = alpha + U_p + U_i and person dispersion
    % nu_p = exp(beta + V_p). The reference gives, per person, the pooled residual
    %   E_p = M + C_p,   M = B expm1(sP^2) expm1(sI^2) (POPULATION, at alpha),
    %                    C_p = 2 exp(2 (alpha + u_p) + 2 sI^2 - beta - v_p) (OWN),
    % and G_p = P / (P + E_p/n) for the n-trial mean (relative error); for two
    % events the difference signal S = P_1 + P_2 - 2 CovP12 with the
    % non-concurrent error E_1/n + E_2/n or the equal-count concurrent error
    % (E_1 + E_2 - 2K)/n. Counts are explicit arguments (n = 40 in the anchors).
    %
    % LITERAL ANCHORS. The numeric anchors below are copied from
    % the frozen 2C Q8/Q9 expectation set (local development evidence)
    % (S20_model_truth_anchors), which the S20 reference memo records as produced
    % with high-precision Decimal arithmetic and without psyrat_gamma_varcomps
    % or any other production routine. They were not derived by this test's
    % author and are truncated here to binary64. The static typical-person D
    % anchor comes from root's independent Decimal reference
    % (the 3A independent Decimal expectation set, local development evidence).
    % Agreement with the anchors here shows the helper evaluates the stated
    % algebra, and is not by itself an independent scientific proof of the target.
    %
    % Every check either pins a literal, exercises an exact algebraic boundary
    % (expm1(0) = 0, zero factor variance), or shows the helper refuses an
    % invalid reference instead of clipping it into validity.

    properties (Constant)
        % Event 1 fixture of the static/difference tests: amplitude 5.5, nu 15,
        % sP .40, sI .20, n 40. Anchors from S20_model_truth_anchors, event 1.
        E1 = struct('alpha', log(5.5), 'beta', log(15), 'sP', 0.40, 'sI', 0.20, 'n', 40, ...
            'P', 6.410781356105155, 'I', 1.507853362877354, 'M', 0.26162895032077987, ...
            'C00', 4.369257839622333, 'E00', 4.630886789943113, ...
            'G00', 0.9822613620971669, 'G00_stale', 0.9832467436328579, ...
            'C_a', 9.72396214921608,   'E_a', 9.98559109953686,  'G_a', 0.9625189199780029, ...
            'C_b', 1.6073601323742072, 'E_b', 1.868989082694987, 'G_b', 0.9927642767245495);
        % Event 2 fixture of the difference tests: amplitude 6.0, nu 13, sP .42, sI .22.
        E2 = struct('alpha', log(6.0), 'beta', log(13), 'sP', 0.42, 'sI', 0.22, 'n', 40, ...
            'P', 8.695572664629095, 'I', 2.2352678805618144, 'M', 0.43121699215824995, ...
            'C00', 6.101390899838907, 'E00', 6.532607891997157, 'G00', 0.9815648194394949, ...
            'C_a', 13.578895168322973, 'E_a', 14.010112160481223, 'G_a', 0.9612801722489295, ...
            'C_b', 2.2445762746012606, 'E_b', 2.6757932667595105, 'G_b', 0.9923657540872763);
        % person-effect cases shared by both events: (u, v) = (.25, -.3) and (-.4, .2)
        U_A = 0.25;  V_A = -0.30;
        U_B = -0.40; V_B = 0.20;
        % two-event anchors (rho_mean .5, independent trial effects, n 40)
        RHO_MEAN = 0.5;
        COVP12 = 3.576065107863202;
        SIGNAL = 7.954223805007846;
        NONC_ZERO_ERROR = 0.27908736704850673;
        NONC_ZERO_G = 0.9661026577015913;
        COPULA_M_INCREMENT = 0.6928459424790298;   % M_1 + M_2 before division by n
        TOL = 1e-12;
    end

    methods (Test)

        function testEvent1PopulationComponentsMatchFrozenAnchors(testCase)
            % B, P, I, M for the static fixture; exact-form literals.
            f = TestGamma3ARecoveryReference.E1;
            c = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f.alpha, f.sP, f.sI);
            testCase.verifyEqual(c.P, f.P, 'RelTol', TestGamma3ARecoveryReference.TOL);
            testCase.verifyEqual(c.I, f.I, 'RelTol', TestGamma3ARecoveryReference.TOL);
            testCase.verifyEqual(c.M, f.M, 'RelTol', TestGamma3ARecoveryReference.TOL);
            % B is the common factor: P = B expm1(sP^2) and M = P expm1(sI^2)
            testCase.verifyEqual(c.B, exp(2*f.alpha + f.sP^2 + f.sI^2), 'RelTol', 1e-14);
            testCase.verifyEqual(c.M, c.P * expm1(f.sI^2), 'RelTol', 1e-14);
        end

        function testEvent2PopulationComponentsMatchFrozenAnchors(testCase)
            f = TestGamma3ARecoveryReference.E2;
            c = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f.alpha, f.sP, f.sI);
            testCase.verifyEqual(c.P, f.P, 'RelTol', TestGamma3ARecoveryReference.TOL);
            testCase.verifyEqual(c.I, f.I, 'RelTol', TestGamma3ARecoveryReference.TOL);
            testCase.verifyEqual(c.M, f.M, 'RelTol', TestGamma3ARecoveryReference.TOL);
        end

        function testPooledResidualAndCoefficientMatchFrozenAnchors(testCase)
            % Three planted persons per event: (0,0), (.25,-.3), (-.4,.2). C is
            % OWN, M is POPULATION (the same M for every person of the event),
            % E = M + C exactly once, G = P/(P + E/n).
            for fx = {TestGamma3ARecoveryReference.E1, TestGamma3ARecoveryReference.E2}
                f = fx{1};   % loop over the two event fixtures
                c = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f.alpha, f.sP, f.sI);
                u = [0; TestGamma3ARecoveryReference.U_A; TestGamma3ARecoveryReference.U_B];
                v = [0; TestGamma3ARecoveryReference.V_A; TestGamma3ARecoveryReference.V_B];
                C = PsyRAT3ARecoveryReference.ownObservationVariance(f.alpha, u, f.sI, f.beta, v);
                E = PsyRAT3ARecoveryReference.pooledOwnResidual(c.M, C);
                G = PsyRAT3ARecoveryReference.oneFacetCoefficients(c.P, c.I, E, f.n);
                testCase.verifyEqual(C, [f.C00; f.C_a; f.C_b], 'RelTol', TestGamma3ARecoveryReference.TOL);
                testCase.verifyEqual(E, [f.E00; f.E_a; f.E_b], 'RelTol', TestGamma3ARecoveryReference.TOL);
                testCase.verifyEqual(G, [f.G00; f.G_a; f.G_b], 'RelTol', TestGamma3ARecoveryReference.TOL);
                % M is not recomputed at alpha + u: the increment E - C is the
                % same population M for every person of the event
                testCase.verifyEqual(E - C, repmat(c.M, 3, 1), 'AbsTol', 1e-13);
            end
        end

        function testPositiveMDiscriminatesOmittedAndDoubledM(testCase)
            % The anchors carry a strictly positive M. Omitting it reproduces
            % the historical stale coefficient; doubling it moves the coefficient
            % the other way. Both differ from the corrected value by far more
            % than the pin tolerance, so the pins above cannot pass with M
            % omitted, doubled, or recomputed per person.
            f = TestGamma3ARecoveryReference.E1;
            c = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f.alpha, f.sP, f.sI);
            testCase.verifyGreaterThan(c.M, 0);
            C00 = PsyRAT3ARecoveryReference.ownObservationVariance(f.alpha, 0, f.sI, f.beta, 0);
            Gcorr = PsyRAT3ARecoveryReference.oneFacetCoefficients(c.P, c.I, c.M + C00, f.n);
            Gomit = PsyRAT3ARecoveryReference.oneFacetCoefficients(c.P, c.I, C00, f.n);
            Gdbl  = PsyRAT3ARecoveryReference.oneFacetCoefficients(c.P, c.I, 2*c.M + C00, f.n);
            testCase.verifyEqual(Gcorr, f.G00, 'RelTol', TestGamma3ARecoveryReference.TOL);
            testCase.verifyEqual(Gomit, f.G00_stale, 'RelTol', TestGamma3ARecoveryReference.TOL);
            testCase.verifyGreaterThan(abs(Gomit - Gcorr), 1e-4, ...
                'Omitting M must move the coefficient well beyond the pin tolerance.');
            testCase.verifyGreaterThan(abs(Gdbl - Gcorr), 1e-4, ...
                'Doubling M must move the coefficient well beyond the pin tolerance.');
            testCase.verifyGreaterThan(Gomit, Gcorr);   % omitted M is optimistic
            testCase.verifyLessThan(Gdbl, Gcorr);       % doubled M is pessimistic
            % a per-person M recomputed at alpha + u differs from the population M
            cShift = PsyRAT3ARecoveryReference.oneFacetMeanComponents( ...
                f.alpha + TestGamma3ARecoveryReference.U_A, f.sP, f.sI);
            testCase.verifyGreaterThan(abs(cShift.M - c.M) / c.M, 0.1);
        end

        function testDependabilityAbsoluteErrorAddsTrialMainEffect(testCase)
            % D = P/(P + E/n + I/n): 1/D - 1/G = I/(n P) exactly, D < G when I > 0,
            % and the hand-derived event-1 typical-person D (pinned loosely).
            f = TestGamma3ARecoveryReference.E1;
            c = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f.alpha, f.sP, f.sI);
            [G, D] = PsyRAT3ARecoveryReference.oneFacetCoefficients(c.P, c.I, f.E00, f.n);
            testCase.verifyEqual(G, f.G00, 'RelTol', TestGamma3ARecoveryReference.TOL);
            testCase.verifyLessThan(D, G);
            testCase.verifyEqual(1/D - 1/G, c.I / (f.n * c.P), 'RelTol', 1e-10);
            % root's independent Decimal reference (no production oracle):
            % the 3A independent Decimal expectation set (local development evidence),
            % static typical-person D = 0.9766205572463919673
            testCase.verifyEqual(D, 0.97662055724639197, 'RelTol', TestGamma3ARecoveryReference.TOL);
        end

        function testZeroTrialVarianceReduction(testCase)
            % sI = 0: no trial facet, so I = M = 0 EXACTLY (expm1(0) = 0, no
            % clipping), E = C, D = G, and P = exp(2 alpha + sP^2) expm1(sP^2).
            f = TestGamma3ARecoveryReference.E1;
            c = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f.alpha, f.sP, 0);
            testCase.verifyEqual(c.I, 0);
            testCase.verifyEqual(c.M, 0);
            testCase.verifyEqual(c.P, exp(2*f.alpha + f.sP^2) * expm1(f.sP^2), 'RelTol', 1e-14);
            C = PsyRAT3ARecoveryReference.ownObservationVariance(f.alpha, 0.1, 0, f.beta, -0.2);
            testCase.verifyEqual(C, 2 * exp(2*(f.alpha + 0.1) - f.beta + 0.2), 'RelTol', 1e-14);
            E = PsyRAT3ARecoveryReference.pooledOwnResidual(c.M, C);
            testCase.verifyEqual(E, C);
            [G, D] = PsyRAT3ARecoveryReference.oneFacetCoefficients(c.P, c.I, E, f.n);
            testCase.verifyEqual(D, G);
            testCase.verifyEqual(G, c.P / (c.P + C / f.n), 'RelTol', 1e-14);
        end

        function testZeroPersonVarianceBoundary(testCase)
            % sP = 0: P = M = 0 exactly, I > 0, and G = D = 0 exactly (0 over a
            % positive error), with no clipping anywhere.
            f = TestGamma3ARecoveryReference.E1;
            c = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f.alpha, 0, f.sI);
            testCase.verifyEqual(c.P, 0);
            testCase.verifyEqual(c.M, 0);
            testCase.verifyGreaterThan(c.I, 0);
            C = PsyRAT3ARecoveryReference.ownObservationVariance(f.alpha, 0, f.sI, f.beta, 0);
            testCase.verifyEqual(C, f.C00, 'RelTol', TestGamma3ARecoveryReference.TOL);
            [G, D] = PsyRAT3ARecoveryReference.oneFacetCoefficients(c.P, c.I, c.M + C, f.n);
            testCase.verifyEqual(G, 0);
            testCase.verifyEqual(D, 0);
        end

        function testPersonPermutationAndIdAlignment(testCase)
            % The own-person quantities are elementwise in (u_p, v_p): permuting
            % the persons permutes the outputs bit for bit. The live tests rely
            % on this when they align recovered rows to the truth by parsing the
            % 'S###' id. C_p depends on the person ONLY through the combination
            % 2 u_p - v_p, so distinct persons with equal 2u - v deliberately
            % share C and G (persons 2 and 4 below: -.4 both); persons with
            % distinct combinations must differ. A first candidate of this test
            % asserted all-distinct G and failed on exactly that pair.
            f = TestGamma3ARecoveryReference.E1;
            c = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f.alpha, f.sP, f.sI);
            u = [0.10; -0.20; 0.30; 0.00; 0.25];
            v = [0.05;  0.00; -0.10; 0.40; -0.30];
            combo = 2*u - v;                          % [.15; -.4; .7; -.4; .8]
            perm = [3 1 5 2 4];
            C  = PsyRAT3ARecoveryReference.ownObservationVariance(f.alpha, u, f.sI, f.beta, v);
            Cp = PsyRAT3ARecoveryReference.ownObservationVariance(f.alpha, u(perm), f.sI, f.beta, v(perm));
            testCase.verifyEqual(Cp, C(perm));
            G  = PsyRAT3ARecoveryReference.oneFacetCoefficients(c.P, c.I, c.M + C, f.n);
            Gp = PsyRAT3ARecoveryReference.oneFacetCoefficients(c.P, c.I, c.M + Cp, f.n);
            testCase.verifyEqual(Gp, G(perm));
            % equal-combination reduction: persons 2 and 4 (different u, v; same
            % 2u - v) get the same C and G up to rounding of the exponent
            testCase.verifyEqual(combo(2), combo(4));
            testCase.verifyEqual(C(2), C(4), 'RelTol', 1e-14);
            testCase.verifyEqual(G(2), G(4), 'RelTol', 1e-14);
            % distinguishable control: the four distinct combinations give four
            % distinct C and G values, and ordering follows 2u - v
            distinct = [1 2 3 5];
            testCase.verifyEqual(numel(unique(combo(distinct))), 4);
            testCase.verifyEqual(numel(unique(C(distinct))), 4, ...
                'Persons with distinct 2u - v must get distinct C.');
            testCase.verifyEqual(numel(unique(G(distinct))), 4, ...
                'Persons with distinct 2u - v must get distinct G.');
            [~, orderCombo] = sort(combo(distinct));
            [~, orderC] = sort(C(distinct));
            testCase.verifyEqual(orderC, orderCombo, 'C must increase with 2u - v.');
            % a larger own mean or a smaller own log-nu raises C (and lowers G)
            testCase.verifyGreaterThan(C(3), C(1));   % u .30 > .10, v -.10 < .05
            testCase.verifyLessThan(G(3), G(1));
        end

        function testTwoEventSignalMatchesFrozenAnchors(testCase)
            % CovP12 = m1 m2 expm1(rho sP1 sP2), S = P1 + P2 - 2 CovP12, and the
            % typical-person (u = v = 0) non-concurrent error and coefficient.
            f1 = TestGamma3ARecoveryReference.E1; f2 = TestGamma3ARecoveryReference.E2;
            c1 = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f1.alpha, f1.sP, f1.sI);
            c2 = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f2.alpha, f2.sP, f2.sI);
            m1 = PsyRAT3ARecoveryReference.eventMean(f1.alpha, f1.sP, f1.sI);
            m2 = PsyRAT3ARecoveryReference.eventMean(f2.alpha, f2.sP, f2.sI);
            testCase.verifyEqual(m1, exp(f1.alpha + (f1.sP^2 + f1.sI^2)/2), 'RelTol', 1e-14);
            covP = PsyRAT3ARecoveryReference.personCrossCovariance(m1, m2, ...
                TestGamma3ARecoveryReference.RHO_MEAN * f1.sP * f2.sP);
            testCase.verifyEqual(covP, TestGamma3ARecoveryReference.COVP12, ...
                'RelTol', TestGamma3ARecoveryReference.TOL);
            S = c1.P + c2.P - 2*covP;
            testCase.verifyEqual(S, TestGamma3ARecoveryReference.SIGNAL, ...
                'RelTol', TestGamma3ARecoveryReference.TOL);
            err0 = PsyRAT3ARecoveryReference.differenceError(f1.E00, f1.n, f2.E00, f2.n, 0);
            testCase.verifyEqual(err0, TestGamma3ARecoveryReference.NONC_ZERO_ERROR, ...
                'RelTol', TestGamma3ARecoveryReference.TOL);
            G0 = PsyRAT3ARecoveryReference.differenceCoefficient(S, err0);
            testCase.verifyEqual(G0, TestGamma3ARecoveryReference.NONC_ZERO_G, ...
                'RelTol', TestGamma3ARecoveryReference.TOL);
            % the copula diagonal pooling increment before division is M1 + M2
            testCase.verifyEqual(c1.M + c2.M, TestGamma3ARecoveryReference.COPULA_M_INCREMENT, ...
                'RelTol', TestGamma3ARecoveryReference.TOL);
            % the person covariance is SIGNAL: subtracting it from the error
            % instead would give a different coefficient
            testCase.verifyGreaterThan(abs(G0 - (c1.P + c2.P) / (c1.P + c2.P + err0 - 2*covP/f1.n)), 1e-4);
        end

        function testCrossEventCovarianceSignAndZeroControls(testCase)
            % Person covariance: exactly zero at cP12 = 0, sign of cP12 otherwise.
            % Mean-interaction cross-covariance: exactly zero whenever cI12 = 0
            % (the fixtures' independent trial effects), nonzero when both
            % log-scale covariances are nonzero, sign = sign(cP12) sign(cI12).
            m1 = 5.7; m2 = 6.3;
            testCase.verifyEqual(PsyRAT3ARecoveryReference.personCrossCovariance(m1, m2, 0), 0);
            testCase.verifyGreaterThan(PsyRAT3ARecoveryReference.personCrossCovariance(m1, m2, 0.084), 0);
            testCase.verifyLessThan(PsyRAT3ARecoveryReference.personCrossCovariance(m1, m2, -0.084), 0);
            testCase.verifyEqual(PsyRAT3ARecoveryReference.personCrossCovariance(m1, m2, 0.084), ...
                m1 * m2 * expm1(0.084), 'RelTol', 1e-14);
            testCase.verifyEqual(PsyRAT3ARecoveryReference.trialMeanCrossCovariance(m1, m2, 0.084, 0), 0);
            testCase.verifyEqual(PsyRAT3ARecoveryReference.trialMeanCrossCovariance(m1, m2, 0, 0.02), 0);
            testCase.verifyGreaterThan(PsyRAT3ARecoveryReference.trialMeanCrossCovariance(m1, m2, 0.084, 0.02), 0);
            testCase.verifyLessThan(PsyRAT3ARecoveryReference.trialMeanCrossCovariance(m1, m2, 0.084, -0.02), 0);
        end

        function testConcurrentErrorEqualCountsOnly(testCase)
            % (E1 + E2 - 2K)/n at equal counts; K ~= 0 with unequal counts is
            % refused (not generalized); K = 0 allows unequal counts.
            err = PsyRAT3ARecoveryReference.differenceError([4; 6], 40, [5; 7], 40, [1; 2]);
            testCase.verifyEqual(err, [(4 + 5 - 2)/40; (6 + 7 - 4)/40], 'RelTol', 1e-14);
            errU = PsyRAT3ARecoveryReference.differenceError(4, 40, 5, 20, 0);
            testCase.verifyEqual(errU, 4/40 + 5/20, 'RelTol', 1e-14);
            testCase.verifyError(@() PsyRAT3ARecoveryReference.differenceError(4, 40, 5, 20, 1), ...
                'PsyRAT3ARecoveryReference:unequalCounts');
            % concurrency with K > 0 raises the coefficient relative to K = 0
            Gc = PsyRAT3ARecoveryReference.differenceCoefficient(8, err(1));
            Gn = PsyRAT3ARecoveryReference.differenceCoefficient(8, (4 + 5)/40);
            testCase.verifyGreaterThan(Gc, Gn);
        end

        function testCovarianceMcseHandComputed(testCase)
            % x = [1 2 3 4], y = [1 3 2 4]: centered cross-products [2.25 -.25 -.25 2.25],
            % K = 4/3, var(p) = 6.25/3, mcse = sqrt(var(p)/4) = sqrt(25/48).
            x = [1; 2; 3; 4]; y = [1; 3; 2; 4];
            [mcse, ifvar, K] = PsyRAT3ARecoveryReference.covarianceMcse(x, y);
            testCase.verifyEqual(K, 4/3, 'RelTol', 1e-14);
            c = cov(x, y);
            testCase.verifyEqual(K, c(1,2), 'RelTol', 1e-14, 'K must be the cov(x,y) off-diagonal.');
            testCase.verifyEqual(ifvar, 6.25/3, 'RelTol', 1e-14);
            testCase.verifyEqual(mcse, sqrt(25/48), 'RelTol', 1e-14);
            % constant cross-products: zero Monte-Carlo error, nonzero K
            z = [1; -1; 1; -1];
            [mcse0, ~, K0] = PsyRAT3ARecoveryReference.covarianceMcse(z, z);
            testCase.verifyEqual(mcse0, 0);
            testCase.verifyEqual(K0, 4/3, 'RelTol', 1e-14);
            % the row/column orientation of the inputs does not matter
            testCase.verifyEqual(PsyRAT3ARecoveryReference.covarianceMcse(x', y'), mcse);
            testCase.verifyError(@() PsyRAT3ARecoveryReference.covarianceMcse([1 2], [1 2 3]), ...
                'PsyRAT3ARecoveryReference:badMcseInput');
        end

        function testDeltaMethodDerivativeMatchesFiniteDifference(testCase)
            % S = 2, E1 = E2 = 1, n = 4, K = 0: err = .5, G = .8, dG/dK = 2*2/(4*2.5^2) = .16.
            S = 2; E1 = 1; E2 = 1; n = 4;
            err0 = PsyRAT3ARecoveryReference.differenceError(E1, n, E2, n, 0);
            testCase.verifyEqual(PsyRAT3ARecoveryReference.differenceCoefficient(S, err0), 0.8, 'RelTol', 1e-14);
            d = PsyRAT3ARecoveryReference.dGdK(S, err0, n);
            testCase.verifyEqual(d, 0.16, 'RelTol', 1e-12);
            h = 1e-4;
            Gp = PsyRAT3ARecoveryReference.differenceCoefficient(S, ...
                PsyRAT3ARecoveryReference.differenceError(E1, n, E2, n, h));
            Gm = PsyRAT3ARecoveryReference.differenceCoefficient(S, ...
                PsyRAT3ARecoveryReference.differenceError(E1, n, E2, n, -h));
            testCase.verifyEqual((Gp - Gm) / (2*h), d, 'AbsTol', 1e-6);
            % propagation is linear in the K standard error
            testCase.verifyEqual(d * 0.01, 0.0016, 'RelTol', 1e-12);
        end

        function testInvalidReferencesAreRefusedNotClipped(testCase)
            % non-finite inputs, a zero-over-zero coefficient, a negative
            % concurrent error and a negative signal all error out; nothing is
            % floored, clipped or NaN-repaired.
            f = TestGamma3ARecoveryReference.E1;
            testCase.verifyError(@() PsyRAT3ARecoveryReference.oneFacetMeanComponents(Inf, f.sP, f.sI), ...
                'PsyRAT3ARecoveryReference:nonFinite');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.ownObservationVariance(f.alpha, NaN, f.sI, f.beta, 0), ...
                'PsyRAT3ARecoveryReference:nonFinite');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.oneFacetCoefficients(0, 0, 0, f.n), ...
                'PsyRAT3ARecoveryReference:nonFinite');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.oneFacetCoefficients(1, 0, -0.5, f.n), ...
                'PsyRAT3ARecoveryReference:negativeValue');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.differenceError(1, 4, 1, 4, 5), ...
                'PsyRAT3ARecoveryReference:negativeValue');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.differenceCoefficient(-1, 0.5), ...
                'PsyRAT3ARecoveryReference:negativeValue');
            % negative SDs are refused (squaring would launder them into valid variances)
            testCase.verifyError(@() PsyRAT3ARecoveryReference.oneFacetMeanComponents(f.alpha, -f.sP, f.sI), ...
                'PsyRAT3ARecoveryReference:negativeValue');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.eventMean(f.alpha, f.sP, -f.sI), ...
                'PsyRAT3ARecoveryReference:negativeValue');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.ownObservationVariance(f.alpha, 0, -f.sI, f.beta, 0), ...
                'PsyRAT3ARecoveryReference:negativeValue');
            % pooled residual: negative or non-finite M / C are refused
            testCase.verifyError(@() PsyRAT3ARecoveryReference.pooledOwnResidual(-0.1, 4), ...
                'PsyRAT3ARecoveryReference:negativeValue');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.pooledOwnResidual(0.2, [4; -1]), ...
                'PsyRAT3ARecoveryReference:negativeValue');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.pooledOwnResidual(NaN, 4), ...
                'PsyRAT3ARecoveryReference:nonFinite');
            % counts must be finite and strictly positive
            testCase.verifyError(@() PsyRAT3ARecoveryReference.oneFacetCoefficients(1, 0, 1, 0), ...
                'PsyRAT3ARecoveryReference:badCount');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.oneFacetCoefficients(1, 0, 1, Inf), ...
                'PsyRAT3ARecoveryReference:badCount');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.differenceError(1, 0, 1, 4, 0), ...
                'PsyRAT3ARecoveryReference:badCount');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.differenceError(1, 4, 1, -4, 0), ...
                'PsyRAT3ARecoveryReference:badCount');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.dGdK(2, 0.5, NaN), ...
                'PsyRAT3ARecoveryReference:badCount');
            testCase.verifyError(@() PsyRAT3ARecoveryReference.dGdK(2, -0.5, 4), ...
                'PsyRAT3ARecoveryReference:negativeValue');
            % the valid fixture still evaluates (the guards do not touch valid numerics)
            c = PsyRAT3ARecoveryReference.oneFacetMeanComponents(f.alpha, f.sP, f.sI);
            testCase.verifyEqual(c.M, f.M, 'RelTol', TestGamma3ARecoveryReference.TOL);
        end

    end
end
