classdef TestModularCutCopulaKernels < PsyRATTestBase
    %TESTMODULARCUTCOPULAKERNELS Offline checks for the modular cut-copula
    % numerics: Gauss-Hermite quadrature, the stable Gamma probability integral
    % transform, the Gaussian-copula likelihood and its conditional grid
    % posterior, the copula-to-covariance converter, and the stage-2 orchestrator.
    %
    % No CmdStan required: every test drives the kernels directly, either from
    % closed-form identities or from FROZEN REFERENCE VALUES.
    %
    % WHY THE FIXTURES MATTER. The accuracy suite is otherwise regression-only -
    % its oracle duplicates the production formulas, so a green run proves
    % "unchanged", not "correct". The CSVs under
    % tests/baselines/modular_cut_copula/ are different in kind: they were
    % produced by the owner's independent R + CmdStan reference implementation
    % (person_specific_dynamic_concurrent_modular), not by this code. A test that
    % matches them is externally grounded. Regenerate them only from that bundle;
    % see the provenance CSV alongside them for the source file checksums.
    %
    % Organized around the claims each kernel makes:
    %   1 Gauss-Hermite integrates a standard normal exactly to degree 2n-1
    %   2 the mean-SD Gamma parameterization is exact, not approximate
    %   3 the probability transform stays accurate in BOTH tails, including
    %     where the lower-tail form saturates at 1 and would return Inf
    %   4 the sufficient-statistic likelihood equals the rowwise sum
    %   5 the grid posterior recovers a known copula correlation
    %   6 the quadrature reproduces the reference covariance map and respects
    %     the Cauchy-Schwarz bound
    %   7 the stage-2 orchestrator reproduces the reference end to end, for both
    %     the one-facet and the two-facet design

    methods (Test)

        %% -- Claim 1: Gauss-Hermite ------------------------------------------

        function testGaussHermiteIntegratesStandardNormalMoments(testCase)
            %The defining property: sum(w.*x.^k) must equal the standard normal
            %moments 1, 0, 1, 0, 3, 0, 15 for k = 0..6. A rule that was
            %mis-normalized (physicists' weights summing to sqrt(pi)) or
            %mis-scaled (nodes not multiplied by sqrt(2)) fails immediately.
            expected = [1 0 1 0 3 0 15];
            for n = [2 5 24 32 48 64]
                [x, w] = psyrat_gauss_hermite(n);
                testCase.verifyEqual(numel(x), n);
                testCase.verifyEqual(numel(w), n);
                %Nodes ascending, and symmetric about zero for this weight.
                testCase.verifyEqual(x, sort(x), 'AbsTol', 0);
                testCase.verifyEqual(x, -flipud(x), 'AbsTol', 1e-12);
                for k = 0:6
                    if 2*n-1 < k, continue, end
                    testCase.verifyEqual(sum(w .* x.^k), expected(k+1), ...
                        'AbsTol', 1e-11, ...
                        sprintf('n = %d failed the degree-%d moment.', n, k));
                end
            end
        end

        function testGaussHermiteCachingReturnsIdenticalValues(testCase)
            %The cache must return the same doubles, not merely close ones -
            %a cache that rebuilt or reordered would perturb every downstream
            %quadrature by a different amount on the second call.
            [x1, w1] = psyrat_gauss_hermite(32);
            [x2, w2] = psyrat_gauss_hermite(32);
            testCase.verifyEqual(x2, x1, 'AbsTol', 0);
            testCase.verifyEqual(w2, w1, 'AbsTol', 0);
        end

        function testGaussHermiteRejectsDegenerateNodeCounts(testCase)
            testCase.verifyError(@() psyrat_gauss_hermite(1), ...
                'psyrat_gauss_hermite:nodes');
            testCase.verifyError(@() psyrat_gauss_hermite(8.5), ...
                'psyrat_gauss_hermite:nodes');
        end

        %% -- Claim 2: the mean-SD parameterization ---------------------------

        function testGammaMeanSdParameterizationIsExact(testCase)
            %Port of static check 2 in the reference bundle. shape = (mu/sigma)^2
            %and rate = mu/sigma^2 must give E(Y) = mu and Var(Y) = sigma^2
            %EXACTLY, which is what lets the scale submodel be read as a residual
            %SD rather than a dispersion proxy.
            for mu = [0.4 5.5 6.2 140]
                for sigma = [0.1 1.4 3.0]
                    shape = (mu/sigma)^2;
                    rate = mu/sigma^2;
                    testCase.verifyEqual(shape/rate, mu, 'RelTol', 1e-12);
                    testCase.verifyEqual(shape/rate^2, sigma^2, 'RelTol', 1e-12);
                end
            end
        end

        %% -- Claim 3: the probability transform ------------------------------

        function testGammaPitInvertsKnownNormalScores(testCase)
            %Port of static check 5: generate scores from known margins through
            %known normal scores, then transform back. Round-trip error must be
            %negligible, which is what makes the copula stage's normal scores
            %trustworthy.
            n = 12000;
            rng(9102, 'twister');
            latent = randn(n, 1);
            mu = 5.5; sigma = 1.7;
            shape = (mu/sigma)^2; rate = mu/sigma^2;
            y = gammaincinv(0.5*erfc(-latent/sqrt(2)), shape) ./ rate;

            pit = psyrat_gamma_pit(y, mu, sigma);

            testCase.verifyLessThan(max(abs(pit.z - latent)), 1e-6);
            testCase.verifyEqual(pit.lower_clipped, 0);
            testCase.verifyEqual(pit.upper_clipped, 0);
        end

        function testGammaPitSurvivesUpperTailSaturation(testCase)
            %THE FAILURE MODE THIS DESIGN EXISTS TO PREVENT. Past about z = +8.3
            %the lower-tail CDF Phi(z) rounds to exactly 1 in double precision,
            %and inverting a saturated 1 returns Inf. The two-tail form must
            %return a large FINITE positive score instead. Without this the
            %48- and 64-node quadrature rules, which reach |z| ~ 14.9, would
            %poison every covariance they touch.
            mu = 5.5; sigma = 1.7;
            shape = (mu/sigma)^2; rate = mu/sigma^2;

            %A score whose upper tail is ~1e-20: far past double saturation of
            %the lower tail, but well inside the clip floor of 1e-12... so it
            %WILL clip. Use a milder one for the finiteness claim and confirm
            %the naive form would have saturated.
            y_far = gammaincinv(1e-20, shape, 'upper') ./ rate;
            naive = gammainc(rate*y_far, shape);
            testCase.verifyEqual(naive, 1, 'AbsTol', 0, ...
                'Fixture is too mild: the naive lower tail did not saturate.');

            pit = psyrat_gamma_pit(y_far, mu, sigma, 1e-30);
            testCase.verifyTrue(isfinite(pit.z));
            testCase.verifyGreaterThan(pit.z, 9);
            testCase.verifyEqual(pit.upper_clipped, 0);
        end

        function testGammaPitMatchesReferenceNormalScores(testCase)
            %Frozen reference values, including deliberately planted extreme
            %scores in both tails.
            %
            %WHERE THIS PORT IS DELIBERATELY MORE ACCURATE THAN THE REFERENCE.
            %An observation whose UPPER tail is clipped has CDF 1 - 1e-12, and
            %the reference forms that number as a double before inverting it.
            %Doubles are spaced 1.1e-16 apart next to 1, so representing
            %1 - 1e-12 discards most of the tail's significant digits and the
            %resulting score is off by ~3e-6. This implementation never forms
            %1 - tiny: it inverts the small tail directly and reflects the sign,
            %giving exactly -Phi^-1(clip). The clipped rows are therefore
            %checked against that exact value rather than against the reference,
            %and every other row is required to agree tightly.
            fx = localFixture('modular_cut_copula_pit_expected.csv');
            counts = localFixture('modular_cut_copula_pit_clipcounts.csv');
            clip = counts.clip(1);

            pit = psyrat_gamma_pit(fx.y, fx.mu, fx.sigma, clip);

            %Recover which rows hit the floor, from the data rather than from
            %the implementation under test.
            shape = (fx.mu ./ fx.sigma).^2;
            scaled = (fx.mu ./ fx.sigma.^2) .* fx.y;
            smallest_tail = min(gammainc(scaled, shape), ...
                gammainc(scaled, shape, 'upper'));
            clipped = smallest_tail < clip;
            upper_clipped = clipped & fx.z > 0;

            testCase.verifyEqual(pit.lower_clipped, counts.lower_clipped(1));
            testCase.verifyEqual(pit.upper_clipped, counts.upper_clipped(1));
            %Both clip branches must actually be exercised, or the counts above
            %would agree vacuously.
            testCase.verifyGreaterThan(counts.lower_clipped(1), 0);
            testCase.verifyGreaterThan(counts.upper_clipped(1), 0);

            %Unclipped rows: agreement between two independent implementations.
            testCase.verifyEqual(pit.z(~clipped), fx.z(~clipped), ...
                'AbsTol', 1e-9, ...
                'Unclipped normal scores must match the reference.');

            %Clipped rows: the exact reflected quantile at the floor.
            exact = -sqrt(2)*erfcinv(2*clip);
            testCase.verifyEqual(pit.z(upper_clipped), ...
                repmat(-exact, nnz(upper_clipped), 1), 'AbsTol', 1e-13);
            testCase.verifyEqual(pit.z(clipped & ~upper_clipped), ...
                repmat(exact, nnz(clipped & ~upper_clipped), 1), 'AbsTol', 1e-13);
            %The lower-tail floor IS exactly representable, so those rows agree
            %with the reference; only the upper ones diverge, and by the amount
            %the double spacing next to 1 predicts.
            testCase.verifyEqual(pit.z(clipped & ~upper_clipped), ...
                fx.z(clipped & ~upper_clipped), 'AbsTol', 1e-12);
            testCase.verifyLessThan( ...
                max(abs(pit.z(upper_clipped) - fx.z(upper_clipped))), 1e-5);
        end

        function testGammaPitRejectsInvalidSupportAndParameters(testCase)
            testCase.verifyError(@() psyrat_gamma_pit([1;-2], 5, 1), ...
                'psyrat_gamma_pit:support');
            testCase.verifyError(@() psyrat_gamma_pit([1;2], [5;-1], 1), ...
                'psyrat_gamma_pit:parameters');
            testCase.verifyError(@() psyrat_gamma_pit([1;2], [5;5;5], 1), ...
                'psyrat_gamma_pit:size');
        end

        function testGammaPitRejectsAnOutOfRangeClipFloor(testCase)
            %The floor applies to the SMALLER tail before the sign reflection,
            %so clip >= 0.5 would silently reverse the transform's ordering
            %(every lower-tail score becomes +Phi^-1(clip), every upper-tail
            %one its negative) - a wrong-sign rho_e with no error. Zero,
            %negative and non-finite floors previously died later with a
            %message blaming the mean/SD inputs.
            testCase.verifyError(@() psyrat_gamma_pit([1;2], 5, 1, 0.6), ...
                'psyrat_gamma_pit:clip');
            testCase.verifyError(@() psyrat_gamma_pit([1;2], 5, 1, 0.5), ...
                'psyrat_gamma_pit:clip');
            testCase.verifyError(@() psyrat_gamma_pit([1;2], 5, 1, 0), ...
                'psyrat_gamma_pit:clip');
            testCase.verifyError(@() psyrat_gamma_pit([1;2], 5, 1, -1), ...
                'psyrat_gamma_pit:clip');
            testCase.verifyError(@() psyrat_gamma_pit([1;2], 5, 1, NaN), ...
                'psyrat_gamma_pit:clip');
            %the documented default stays valid
            pit = psyrat_gamma_pit([1;2], 5, 1, 1e-12);
            testCase.verifyTrue(all(isfinite(pit.z)));
        end

        %% -- Claim 4: the copula likelihood ----------------------------------

        function testCopulaLoglikEqualsRowwiseSum(testCase)
            %Port of static check 3. The sufficient-statistic form must equal a
            %direct per-observation sum; that equality is what makes a
            %20,001-point grid affordable at 50,000 pairs.
            rng(9101, 'twister');
            n = 2000;
            z1 = randn(n, 1);
            z2 = 0.37*z1 + sqrt(1 - 0.37^2)*randn(n, 1);
            sufficient = psyrat_copula_suffstats(z1, z2);

            for rho = [-0.8 -0.29 0 0.29 0.8]
                direct = sum( ...
                    -0.5*log(1 - rho^2) ...
                    - 0.5*(rho^2*(z1.^2 + z2.^2) - 2*rho*z1.*z2)/(1 - rho^2));
                testCase.verifyEqual(psyrat_copula_loglik(rho, sufficient), ...
                    direct, 'RelTol', 1e-9);
            end
        end

        function testCopulaLoglikIsNegativeInfiniteOutsideTheUnitInterval(testCase)
            sufficient = psyrat_copula_suffstats(randn(50,1), randn(50,1));
            ll = psyrat_copula_loglik([-1 1 1.5 NaN 0.5], sufficient);
            testCase.verifyEqual(ll(1:4), [-Inf -Inf -Inf -Inf]);
            testCase.verifyTrue(isfinite(ll(5)));
            %Shape of the input is preserved, so grid callers can rely on it.
            testCase.verifyEqual(size(ll), [1 5]);
        end

        %% -- Claim 5: the conditional grid posterior -------------------------

        function testRhoGridRecoversKnownGaussianDependence(testCase)
            %Port of static check 4: with normal scores generated at a known
            %correlation, the conditional posterior must concentrate there.
            %
            %Two claims, kept separate because they can fail for different
            %reasons. The SHARP claim is that the posterior tracks the
            %dependence actually present in this sample - that is a statement
            %about the estimator and holds tightly. The looser claim is that the
            %sample dependence tracks the generating parameter, which is
            %ordinary sampling variability and sets its own scale: at n = 20000
            %the standard error of a correlation near 0.37 is about 0.006.
            rng(9101, 'twister');
            n = 20000;
            rho_true = 0.37;
            z1 = randn(n, 1);
            z2 = rho_true*z1 + sqrt(1 - rho_true^2)*randn(n, 1);
            sufficient = psyrat_copula_suffstats(z1, z2);

            post = psyrat_copula_rho_grid(sufficient, 0.5);

            a = z1 - mean(z1);
            b = z2 - mean(z2);
            sample_r = sum(a.*b) / sqrt(sum(a.^2)*sum(b.^2));

            testCase.verifyEqual(post.mean, sample_r, 'AbsTol', 0.01, ...
                'The posterior must track the dependence present in the data.');
            testCase.verifyEqual(post.mean, rho_true, 'AbsTol', 0.03);
            %A proper posterior: mean, median and mode agree closely and the
            %interval brackets them.
            testCase.verifyEqual(post.median, post.mean, 'AbsTol', 0.01);
            testCase.verifyEqual(post.mode, post.mean, 'AbsTol', 0.01);
            testCase.verifyLessThan(post.q025, post.mean);
            testCase.verifyGreaterThan(post.q975, post.mean);
            testCase.verifyLessThan(post.boundary_mass, 1e-6);
        end

        function testRhoGridMatchesReferenceSummaries(testCase)
            %Frozen reference values across five sufficient-statistic regimes,
            %including a near-null case and one at the scale of the full data.
            fx = localFixture('modular_cut_copula_rhogrid_expected.csv');
            spacing = 1.9/20000;   % one grid step

            for k = 1:height(fx)
                sufficient = struct('n', fx.N(k), ...
                    'sum_sq_1', fx.sum_z_theta_sq(k), ...
                    'sum_sq_2', fx.sum_z_alpha_sq(k), ...
                    'sum_cross', fx.sum_cross(k));

                %The likelihood itself, at two probe values.
                testCase.verifyEqual(psyrat_copula_loglik(0.10, sufficient), ...
                    fx.loglik_at_0_10(k), 'RelTol', 1e-12);
                testCase.verifyEqual(psyrat_copula_loglik(-0.40, sufficient), ...
                    fx.loglik_at_neg_0_40(k), 'RelTol', 1e-12);

                post = psyrat_copula_rho_grid(sufficient, fx.uniform_draw(k));

                %Smooth functionals of the whole grid: tight tolerance.
                testCase.verifyEqual(post.mean, fx.rho_conditional_mean(k), ...
                    'RelTol', 1e-10);
                testCase.verifyEqual(post.sd, fx.rho_conditional_sd(k), ...
                    'RelTol', 1e-10);
                testCase.verifyEqual(post.boundary_mass, ...
                    fx.rho_boundary_mass(k), 'AbsTol', 1e-12);

                %Quantile-type outputs are STEP functions of the grid, so they
                %are pinned to the grid spacing rather than to machine epsilon.
                %rho_e is included here because the uniform is supplied, making
                %the inverse-CDF draw deterministic and comparable.
                testCase.verifyEqual(post.rho_e, fx.rho_e(k), 'AbsTol', spacing);
                testCase.verifyEqual(post.median, ...
                    fx.rho_conditional_median(k), 'AbsTol', spacing);
                testCase.verifyEqual(post.q025, fx.rho_conditional_q025(k), ...
                    'AbsTol', spacing);
                testCase.verifyEqual(post.q975, fx.rho_conditional_q975(k), ...
                    'AbsTol', spacing);
                testCase.verifyEqual(post.mode, fx.rho_conditional_mode(k), ...
                    'AbsTol', spacing);
            end
        end

        function testRhoGridIntegratesCorrectlyWhenMassSitsOnTheBound(testCase)
            %The trapezoid rule differs from a plain sum ONLY in its halved end
            %weights, and those weights are negligible whenever the posterior
            %concentrates inside the support - which every other fixture here
            %does. Cases 6 and 7 pile the posterior against the imposed +/-0.95
            %bound, so the endpoint carries real mass and the integration rule
            %becomes observable. This is also the only regime that exercises
            %rho_boundary_mass in the range that raises the numerical-invalidity
            %flag, so a pile-up is reported rather than silently summarized.
            fx = localFixture('modular_cut_copula_rhogrid_expected.csv');
            boundary = find(fx.rho_boundary_mass > 0.5);
            testCase.assertNotEmpty(boundary, ...
                'Fixture no longer contains a boundary-dominated case.');

            for k = boundary'
                sufficient = struct('n', fx.N(k), ...
                    'sum_sq_1', fx.sum_z_theta_sq(k), ...
                    'sum_sq_2', fx.sum_z_alpha_sq(k), ...
                    'sum_cross', fx.sum_cross(k));
                post = psyrat_copula_rho_grid(sufficient, fx.uniform_draw(k));

                testCase.verifyEqual(post.mean, fx.rho_conditional_mean(k), ...
                    'RelTol', 1e-10);
                testCase.verifyEqual(post.sd, fx.rho_conditional_sd(k), ...
                    'RelTol', 1e-8);
                testCase.verifyEqual(post.boundary_mass, ...
                    fx.rho_boundary_mass(k), 'AbsTol', 1e-12);
                %The bound, not the data, is determining the answer here.
                testCase.verifyGreaterThan(post.boundary_mass, 0.01);
                testCase.verifyLessThanOrEqual(abs(post.mean), 0.95);
            end
        end

        function testRhoGridRejectsDegenerateInputs(testCase)
            sufficient = psyrat_copula_suffstats(randn(50,1), randn(50,1));
            testCase.verifyError(@() psyrat_copula_rho_grid(sufficient, 0), ...
                'psyrat_copula_rho_grid:uniform');
            testCase.verifyError(@() psyrat_copula_rho_grid(sufficient, 1), ...
                'psyrat_copula_rho_grid:uniform');
            %Even grid counts would straddle rho = 0 rather than represent it.
            testCase.verifyError( ...
                @() psyrat_copula_rho_grid(sufficient, 0.5, 20000), ...
                'psyrat_copula_rho_grid:grid');
            %Non-integer counts previously routed AROUND the parity test
            %(mod(2002.5,2) is 0.5, not 0) and linspace silently truncated to
            %an even grid; Inf passed too (mod(Inf,2) is NaN) and died later.
            testCase.verifyError( ...
                @() psyrat_copula_rho_grid(sufficient, 0.5, 20002.5), ...
                'psyrat_copula_rho_grid:grid');
            testCase.verifyError( ...
                @() psyrat_copula_rho_grid(sufficient, 0.5, Inf), ...
                'psyrat_copula_rho_grid:grid');
        end

        function testRhoGridSurvivesAUniformDrawAtTheCdfEndpoint(testCase)
            %cumsum over 20k normalized weights can land ~1e-14 BELOW 1, and a
            %uniform draw in that gap previously made the inverse-CDF find
            %return empty - a bare mid-loop failure in the stage, after the
            %fit. The clamp makes the largest representable uniform draw
            %return the top grid point instead (2026-08-11 code review).
            sufficient = psyrat_copula_suffstats(randn(50,1), randn(50,1));
            out = psyrat_copula_rho_grid(sufficient, 1 - 2^-53);
            testCase.verifyTrue(isscalar(out.rho_e) && isfinite(out.rho_e), ...
                'The endpoint draw must return a scalar rho_e, never empty.');
        end

        %% -- Claim 6: the copula-to-covariance quadrature --------------------

        function testQuadratureMatchesReferenceCovarianceMap(testCase)
            %The externally grounded check on the piece of science that is new
            %to the toolbox: the map from a copula correlation to an
            %observed-scale residual covariance. 54 (nu, nu, rho) combinations.
            fx = localFixture('modular_cut_copula_quadrature_expected.csv');

            for k = 1:height(fx)
                out = psyrat_gamma_copula_rescov(fx.nu_theta(k), ...
                    fx.nu_alpha(k), fx.rho(k));
                testCase.verifyEqual(out.value, fx.value(k), ...
                    'AbsTol', 1e-9, 'RelTol', 1e-7, sprintf( ...
                    'nu = (%g, %g), rho = %g', fx.nu_theta(k), ...
                    fx.nu_alpha(k), fx.rho(k)));
                testCase.verifyEqual(out.nodes, fx.nodes(k), ...
                    'Node escalation must terminate where the reference does.');
            end
        end

        function testQuadratureSignAndBoundBehaviour(testCase)
            %Structural claims that hold for every argument, checked
            %independently of the fixture: the covariance takes the sign of rho,
            %is symmetric in its two components, and respects Cauchy-Schwarz.
            for rho = [-0.85 -0.3 0.087 0.55 0.93]
                a = psyrat_gamma_copula_rescov(18, 30, rho);
                b = psyrat_gamma_copula_rescov(30, 18, rho);
                testCase.verifyEqual(sign(a.value), sign(rho));
                testCase.verifyEqual(a.value, b.value, 'RelTol', 1e-7, ...
                    'The covariance must be symmetric in the two components.');
                testCase.verifyLessThanOrEqual(abs(a.value), ...
                    sqrt((2/18)*(2/30)) + 1e-8);
            end
        end

        function testQuadratureShortCircuitsAtZeroCorrelation(testCase)
            %Independent residuals: exactly zero covariance, and the reported
            %node count of 0 records that no integration was performed.
            out = psyrat_gamma_copula_rescov(18, 30, 0);
            testCase.verifyEqual(out.value, 0);
            testCase.verifyEqual(out.nodes, 0);
            testCase.verifyEqual(out.error, 0);
        end

        function testQuadratureEscalatesAndThenFailsLoudly(testCase)
            %The escalation must be real, not decorative. At a tolerance no
            %finite node pair can meet, the function must ERROR rather than
            %quietly return an unconverged value - an unconverged covariance
            %would propagate into every person's reliability estimate.
            testCase.verifyError( ...
                @() psyrat_gamma_copula_rescov(2.5, 4, 0.93, 1e-18), ...
                'psyrat_gamma_copula_rescov:convergence');

            %A tolerance the first pair cannot meet but a later one can must
            %report the larger node count it actually needed.
            loose = psyrat_gamma_copula_rescov(2.5, 4, 0.93);
            tight = psyrat_gamma_copula_rescov(2.5, 4, 0.93, 1e-13);
            testCase.verifyGreaterThanOrEqual(tight.nodes, loose.nodes);
            testCase.verifyEqual(tight.value, loose.value, 'AbsTol', 1e-9);
        end

        function testQuadratureRejectsInvalidArguments(testCase)
            testCase.verifyError(@() psyrat_gamma_copula_rescov(-1, 4, 0.3), ...
                'psyrat_gamma_copula_rescov:nu');
            testCase.verifyError(@() psyrat_gamma_copula_rescov(18, 4, 1), ...
                'psyrat_gamma_copula_rescov:rho');
        end

        %% -- Claim 7: the stage-2 orchestrator -------------------------------

        function testStageTwoMatchesReferenceOneFacet(testCase)
            localVerifyStageAgainstReference(testCase, 'onefacet');
        end

        function testStageTwoMatchesReferenceTwoFacet(testCase)
            %The two-facet path adds occasion, person x trial, person x occasion
            %and trial x occasion effects to the MEAN only. An indexing error in
            %any of the four would move the reconstructed margins and show up
            %here; the fixture uses distinct non-zero values at every level so
            %it cannot pass by symmetry.
            localVerifyStageAgainstReference(testCase, 'twofacet');
        end

        function testStageTwoReconstructsMarginsToMachinePrecision(testCase)
            %THE SHARPEST CHECK IN THIS FILE, and the one that pins the slot
            %mapping. Every conditional mean and residual SD is compared against
            %the reference's own reconstruction of the SAME draw, elementwise.
            %Unlike the normal scores downstream, these involve no incomplete
            %gamma function, so the two implementations must agree to machine
            %precision; a wrong person-block column, a transposed facet index or
            %a dropped offset would move them far more than that.
            %
            %The two-facet design is used because it exercises all four extra
            %facet indices (occasion, person x trial, person x occasion,
            %trial x occasion) on top of person and trial.
            [draws, standata] = localLoadStageFixture('twofacet');
            ref = localFixture( ...
                'modular_cut_copula_twofacet_draw1_margins.csv');

            for m = 1:2
                [mu, sigma] = psyrat_gamma_margin_predictors(draws, standata, 1, m);
                testCase.verifyEqual(mu, ref.(sprintf('mu_%d', m)), ...
                    'RelTol', 1e-12, sprintf( ...
                    'Component %d conditional means must match the reference.', m));
                testCase.verifyEqual(sigma, ref.(sprintf('sigma_%d', m)), ...
                    'RelTol', 1e-12, sprintf( ...
                    'Component %d residual SDs must match the reference.', m));
            end
        end

        function testMarginPredictorsUseTheDocumentedPersonBlockColumns(testCase)
            %The person block's column order is a convention, and reading it
            %wrongly yields a model that still runs. Perturbing ONE column at a
            %time must move exactly the quantity that column is documented to
            %drive: columns 1-2 the two components' means, columns 3-4 their
            %residual SDs.
            [draws, standata] = localLoadStageFixture('onefacet');
            [mu0, sigma0] = psyrat_gamma_margin_predictors(draws, standata, 1, 1);

            bumped = draws;
            bumped.u_p(1, :, 1) = bumped.u_p(1, :, 1) + 0.5;    % component 1 mean
            [mu1, sigma1] = psyrat_gamma_margin_predictors(bumped, standata, 1, 1);
            testCase.verifyEqual(mu1, mu0 * exp(0.5), 'RelTol', 1e-12);
            testCase.verifyEqual(sigma1, sigma0, 'AbsTol', 0);

            bumped = draws;
            bumped.u_p(1, :, 3) = bumped.u_p(1, :, 3) + 0.5;    % component 1 scale
            [mu2, sigma2] = psyrat_gamma_margin_predictors(bumped, standata, 1, 1);
            testCase.verifyEqual(mu2, mu0, 'AbsTol', 0);
            testCase.verifyEqual(sigma2, sigma0 * exp(0.5), 'RelTol', 1e-12);

            bumped = draws;
            bumped.u_p(1, :, 2) = bumped.u_p(1, :, 2) + 0.5;    % the OTHER component
            [mu3, sigma3] = psyrat_gamma_margin_predictors(bumped, standata, 1, 1);
            testCase.verifyEqual(mu3, mu0, 'AbsTol', 0, ...
                'Component 2''s person effect must not touch component 1.');
            testCase.verifyEqual(sigma3, sigma0, 'AbsTol', 0);
        end

        function testMarginPredictorsApplyTheDataDerivedOffsets(testCase)
            %The offsets are centering constants, not parameters: the fitted
            %intercepts are deviations from them. Dropping one silently rescales
            %every mean and residual SD, which is the trap recorded for the
            %earlier external-validation harness. Shifting an offset must move
            %the reconstruction by exactly the corresponding factor, and must
            %move only the component it belongs to.
            [draws, standata] = localLoadStageFixture('onefacet');
            [mu0, sigma0] = psyrat_gamma_margin_predictors(draws, standata, 1, 1);
            [mu0b, sigma0b] = psyrat_gamma_margin_predictors(draws, standata, 1, 2);

            shifted = standata;
            shifted.mu_offset(1) = shifted.mu_offset(1) + 0.25;
            shifted.log_sigma_offset(1) = shifted.log_sigma_offset(1) - 0.10;

            [mu1, sigma1] = psyrat_gamma_margin_predictors(draws, shifted, 1, 1);
            testCase.verifyEqual(mu1, mu0 * exp(0.25), 'RelTol', 1e-12);
            testCase.verifyEqual(sigma1, sigma0 * exp(-0.10), 'RelTol', 1e-12);

            [mu2, sigma2] = psyrat_gamma_margin_predictors(draws, shifted, 1, 2);
            testCase.verifyEqual(mu2, mu0b, 'AbsTol', 0, ...
                'Component 1''s offset must not move component 2.');
            testCase.verifyEqual(sigma2, sigma0b, 'AbsTol', 0);
        end

        function testStageTwoDrawSelectionIsDeterministicAndEvenlySpaced(testCase)
            %Selection must not consume the random stream and must be
            %reproducible from the draw count alone.
            [draws, standata] = localLoadStageFixture('onefacet');
            big = localReplicateDraws(draws, 40);       % 3 -> 120 draws

            a = psyrat_gamma_cut_copula_stage(big, standata, 12345, ...
                'max_draws', 10, 'grid_points', 2001);
            b = psyrat_gamma_cut_copula_stage(big, standata, 999, ...
                'max_draws', 10, 'grid_points', 2001);

            testCase.verifyEqual(a.sel, b.sel, ...
                'Draw selection must not depend on the seed.');
            testCase.verifyEqual(numel(a.sel), 10);
            testCase.verifyEqual(a.sel(1), 1);
            testCase.verifyEqual(a.sel(end), 120);
            %Evenly spaced to within one index.
            testCase.verifyLessThanOrEqual(max(abs(diff(diff(a.sel)))), 1);
            %A different seed must still move the rho draw, or the seed would be
            %doing nothing and a seed-sensitivity check would be meaningless.
            testCase.verifyNotEqual(a.rho_e, b.rho_e);
            %...but never the conditional summaries, which are deterministic.
            testCase.verifyEqual(a.rho_mean, b.rho_mean, 'AbsTol', 0);
        end

        function testPerStratumSeedsBreakTheComonotoneRhoCoupling(testCase)
            %The estimation call sites derive each group stratum's seed as
            %seed + (stratum - 1), taken from the "for i = 1:ndchunks"
            %stratum fit loop (the 2026-08-11 seed-policy increment). With the
            %run seed reused verbatim - the policy before it - every stratum
            %drew its rho_e at IDENTICAL quantiles: the per-draw sequences
            %were COMONOTONE across groups, so a per-draw between-group
            %contrast rho_A(k) - rho_B(k) had systematically deflated spread
            %and an overconfident interval. This pins both halves of the
            %policy at the stage level: reusing a seed still reproduces
            %bit-identically (and collapses every per-draw contrast to zero,
            %the degenerate comonotone case this fixture's identical margins
            %produce), while the +1 offset the call sites now apply decouples
            %the sequences and gives the contrast real spread.
            [draws, standata] = localLoadStageFixture('onefacet');
            big = localReplicateDraws(draws, 40);   % 3 -> 120 draws

            seed = 12345;
            a  = psyrat_gamma_cut_copula_stage(big, standata, seed, ...
                'max_draws', 10, 'grid_points', 2001);
            a2 = psyrat_gamma_cut_copula_stage(big, standata, seed, ...
                'max_draws', 10, 'grid_points', 2001);
            b  = psyrat_gamma_cut_copula_stage(big, standata, seed + 1, ...
                'max_draws', 10, 'grid_points', 2001);

            %same derived seed: bit-identical, so the reproducibility
            %contract survives the policy change...
            testCase.verifyEqual(a2.rho_e, a.rho_e, 'AbsTol', 0);
            %...which is exactly the old cross-strata failure: zero per-draw
            %contrast everywhere, no spread at all
            testCase.verifyEqual(std(a2.rho_e - a.rho_e), 0, 'AbsTol', 0);

            %the ADJACENT derived seed (the smallest offset the call sites
            %produce) must decouple the sequences and give the contrast
            %real spread
            testCase.verifyNotEqual(b.rho_e, a.rho_e);
            testCase.verifyGreaterThan(std(b.rho_e - a.rho_e), 0);
            %and the derived seed is archived, so a stored stratum can be
            %audited for which policy produced it
            testCase.verifyEqual(a.settings.seed, seed);
            testCase.verifyEqual(b.settings.seed, seed + 1);
        end

        function testStageTwoRecordsFlagsWithoutBlocking(testCase)
            %Both flags are diagnostics. Computation must complete and return
            %usable numbers even when they fire, so the problem can be seen.
            [draws, standata] = localLoadStageFixture('onefacet');
            out = psyrat_gamma_cut_copula_stage(draws, standata, 12345, ...
                'grid_points', 2001);

            testCase.verifyTrue(islogical(out.invalidity_flag));
            testCase.verifyTrue(islogical(out.calibration_flag));
            testCase.verifyTrue(all(isfinite(out.rho_e)));
            %This fixture's synthetic draws do not match its data-generating
            %truth, so the PIT scores are over-dispersed and the DESCRIPTIVE
            %calibration flag fires while the NUMERICAL one does not. That
            %separation is the point of having two flags.
            testCase.verifyTrue(all(out.calibration_flag));
            testCase.verifyFalse(any(out.invalidity_flag));
        end

        function testStageTwoCarriesTheModularProvenanceLabels(testCase)
            %These labels state that the numbers came from a cut posterior
            %rather than a full joint one. Downstream surfaces must be able to
            %read them off the result; losing them would misrepresent the
            %inference.
            [draws, standata] = localLoadStageFixture('onefacet');
            out = psyrat_gamma_cut_copula_stage(draws, standata, 12345, ...
                'grid_points', 2001);

            testCase.verifyEqual(out.labels.inference_framework, ...
                'modular_cut_two_stage_gamma_margins_gaussian_copula_v1');
            testCase.verifyFalse(out.labels.copula_feedback_to_margins);
            testCase.verifyEqual(out.labels.joint_posterior_status, ...
                'not_a_full_joint_bayesian_posterior');
            %Every constant that could move the estimand is archived with it.
            for f = {'max_draws','grid_points','prior_sd','clip','seed', ...
                    'n_draws_used','draw_selection','rho_support'}
                testCase.verifyTrue(isfield(out.settings, f{1}), ...
                    sprintf('settings must record %s', f{1}));
            end
        end

        function testStageTwoRejectsShapeMismatches(testCase)
            [draws, standata] = localLoadStageFixture('onefacet');

            bad = draws;
            bad.u_p = bad.u_p(:, :, 1:3);          % three person columns, not four
            testCase.verifyError( ...
                @() psyrat_gamma_cut_copula_stage(bad, standata, 1), ...
                'psyrat_gamma_cut_copula_stage:draws');

            bad = rmfield(draws, 'beta_sigma');
            testCase.verifyError( ...
                @() psyrat_gamma_cut_copula_stage(bad, standata, 1), ...
                'psyrat_gamma_cut_copula_stage:draws');

            bad = standata;
            bad.y_2 = bad.y_2(1:end-1);            % breaks concurrent pairing
            testCase.verifyError( ...
                @() psyrat_gamma_cut_copula_stage(draws, bad, 1), ...
                'psyrat_gamma_cut_copula_stage:standata');

            testCase.verifyError( ...
                @() psyrat_gamma_cut_copula_stage(draws, standata, 1, 'nope', 1), ...
                'psyrat_gamma_cut_copula_stage:options');

            %The documented public 'clip' option is validated at the stage
            %boundary too, so a bad value names the option instead of
            %surfacing per draw inside the PIT kernel.
            testCase.verifyError( ...
                @() psyrat_gamma_cut_copula_stage(draws, standata, 1, 'clip', 0.6), ...
                'psyrat_gamma_cut_copula_stage:options');
            testCase.verifyError( ...
                @() psyrat_gamma_cut_copula_stage(draws, standata, 1, 'clip', NaN), ...
                'psyrat_gamma_cut_copula_stage:options');
        end

    end
end

% ---------------------------------------------------------------------------
% Fixture helpers
% ---------------------------------------------------------------------------

function T = localFixture(name)
%LOCALFIXTURE Reads one frozen reference CSV.
here = fileparts(mfilename('fullpath'));
T = readtable(fullfile(here, 'baselines', 'modular_cut_copula', name));
end

function [draws, standata] = localLoadStageFixture(tag)
%LOCALLOADSTAGEFIXTURE Rebuilds the stage-2 inputs from the frozen CSVs.
%   The draw parameters are stored in LONG format (draw, param, idx1, idx2,
%   value) precisely so that reassembling them here cannot silently depend on a
%   column ordering convention.

meta = localFixture(sprintf('modular_cut_copula_%s_meta.csv', tag));
obs = localFixture(sprintf('modular_cut_copula_%s_standata.csv', tag));
long = localFixture(sprintf('modular_cut_copula_%s_draws.csv', tag));

nd = meta.n_draws(1);
K = meta.K(1);

standata = struct( ...
    'N', meta.N(1), 'K', K, 'P', meta.P(1), 'I', meta.I(1), ...
    'y_1', obs.y_theta, 'y_2', obs.y_alpha, ...
    'p_id', obs.p_id, 'i_id', obs.i_id, ...
    'mu_offset', [meta.mu_offset_1(1) meta.mu_offset_2(1)], ...
    'log_sigma_offset', ...
        [meta.log_sigma_offset_1(1) meta.log_sigma_offset_2(1)]);

Z = zeros(standata.N, K);
for k = 1:K
    Z(:, k) = obs.(sprintf('Z%d', k));
end
standata.Z = Z;

draws = struct( ...
    'alpha_mu', localGather(long, 'alpha_mu', nd, [2 1]), ...
    'alpha_sigma', localGather(long, 'alpha_sigma', nd, [2 1]), ...
    'beta_mu', localGather(long, 'beta_mu', nd, [2 K]), ...
    'beta_sigma', localGather(long, 'beta_sigma', nd, [2 K]), ...
    'u_p', localGather(long, 'u_p_joint', nd, [meta.P(1) 4]), ...
    'u_i', localGather(long, 'u_i', nd, [meta.I(1) 2]));

if ~isnan(meta.O(1))
    standata.O = meta.O(1);
    standata.PI = meta.PI(1);
    standata.PO = meta.PO(1);
    standata.IO = meta.IO(1);
    standata.o_id = obs.o_id;
    standata.pi_id = obs.pi_id;
    standata.po_id = obs.po_id;
    standata.io_id = obs.io_id;
    draws.u_o = localGather(long, 'u_o', nd, [meta.O(1) 2]);
    draws.u_pi = localGather(long, 'u_pi', nd, [meta.PI(1) 2]);
    draws.u_po = localGather(long, 'u_po', nd, [meta.PO(1) 2]);
    draws.u_io = localGather(long, 'u_io', nd, [meta.IO(1) 2]);
end
end

function A = localGather(long, name, nd, dims)
%LOCALGATHER Long-format rows for one parameter into an [nd x dims] array.
%   A trailing singleton second dimension is dropped so vector parameters come
%   back as [nd x n] rather than [nd x n x 1].
param = string(long.param);
rows = find(param == name);
if isempty(rows)
    error('TestModularCutCopulaKernels:fixture', ...
        'Fixture is missing parameter %s.', name);
end

A = nan([nd dims]);
for r = rows'
    i1 = long.idx1(r);
    if isnan(long.idx2(r))
        A(long.draw(r), i1) = long.value(r);
    else
        A(long.draw(r), i1, long.idx2(r)) = long.value(r);
    end
end
if any(isnan(A(:)))
    error('TestModularCutCopulaKernels:fixture', ...
        'Fixture left holes in parameter %s.', name);
end
end

function big = localReplicateDraws(draws, reps)
%LOCALREPLICATEDRAWS Tiles the draw dimension so thinning has something to thin.
big = draws;
for f = fieldnames(draws)'
    A = draws.(f{1});
    rep = ones(1, ndims(A));
    rep(1) = reps;
    big.(f{1}) = repmat(A, rep);
end
end

function localVerifyStageAgainstReference(testCase, tag)
%LOCALVERIFYSTAGEAGAINSTREFERENCE End-to-end stage-2 comparison.
%   rho_e itself is deliberately NOT compared: it is an inverse-CDF draw from a
%   uniform, and MATLAB's random stream is not R's, so the two draw different
%   points from the SAME posterior. Every uniform-independent quantity is
%   compared instead, and the draw path is pinned separately in
%   testRhoGridMatchesReferenceSummaries, where the uniform is supplied.

[draws, standata] = localLoadStageFixture(tag);
expected = localFixture(sprintf('modular_cut_copula_%s_stage2_expected.csv', tag));
meta = localFixture(sprintf('modular_cut_copula_%s_meta.csv', tag));

out = psyrat_gamma_cut_copula_stage(draws, standata, 12345, ...
    'grid_points', meta.grid_points(1), 'prior_sd', meta.prior_sd(1), ...
    'clip', meta.pit_clip(1));

spacing = 1.9/(meta.grid_points(1) - 1);

% TOLERANCE, AND WHY IT IS NOT MACHINE EPSILON. Every quantity below is a
% function of the normal scores, and the two implementations' incomplete-gamma
% routines differ in the last few digits for observations far out in the tails.
% The largest observed effect is on an observation whose tail probability is
% ~5e-12, where the reference's near-1 CDF representation costs it accuracy
% (see testGammaPitMatchesReferenceNormalScores). The reconstruction itself is
% pinned to machine precision separately, by
% testStageTwoReconstructsMarginsToMachinePrecision, so a genuine slot-mapping
% or indexing error cannot hide inside this tolerance.
pit_tol = 1e-7;

testCase.verifyEqual(numel(out.sel), height(expected));

% Conditional posterior summaries: smooth functionals of the whole grid.
testCase.verifyEqual(out.rho_mean, expected.rho_conditional_mean, 'RelTol', pit_tol);
testCase.verifyEqual(out.rho_sd, expected.rho_conditional_sd, 'RelTol', pit_tol);
testCase.verifyEqual(out.rho_boundary_mass, expected.rho_boundary_mass, ...
    'AbsTol', 1e-12);

% Quantile-type outputs: step functions of the grid.
testCase.verifyEqual(out.rho_median, expected.rho_conditional_median, ...
    'AbsTol', spacing);
testCase.verifyEqual(out.rho_q025, expected.rho_conditional_q025, 'AbsTol', spacing);
testCase.verifyEqual(out.rho_q975, expected.rho_conditional_q975, 'AbsTol', spacing);
testCase.verifyEqual(out.rho_mode, expected.rho_conditional_mode, 'AbsTol', spacing);

% Normal-score diagnostics: these depend on every reconstructed mu and sigma,
% so they are a broad check that the linear predictors were rebuilt correctly
% from the draw slots.
testCase.verifyEqual(out.pit_z_1_mean, expected.pit_z_theta_mean, 'RelTol', pit_tol);
testCase.verifyEqual(out.pit_z_1_sd, expected.pit_z_theta_sd, 'RelTol', pit_tol);
testCase.verifyEqual(out.pit_z_2_mean, expected.pit_z_alpha_mean, 'RelTol', pit_tol);
testCase.verifyEqual(out.pit_z_2_sd, expected.pit_z_alpha_sd, 'RelTol', pit_tol);
testCase.verifyEqual(out.pit_normal_score_correlation, ...
    expected.pit_normal_score_correlation, 'RelTol', pit_tol);
testCase.verifyEqual(out.pit_1_abs_gt_1_96, expected.pit_theta_abs_gt_1_96, ...
    'AbsTol', 1e-12);
testCase.verifyEqual(out.pit_2_abs_gt_1_96, expected.pit_alpha_abs_gt_1_96, ...
    'AbsTol', 1e-12);
testCase.verifyEqual(out.pit_joint_upper_1_96, expected.pit_joint_upper_1_96, ...
    'AbsTol', 1e-12);
testCase.verifyEqual(out.pit_joint_lower_1_96, expected.pit_joint_lower_1_96, ...
    'AbsTol', 1e-12);

% Clip accounting.
testCase.verifyEqual(out.pit_1_lower_clipped, expected.pit_theta_lower_clipped);
testCase.verifyEqual(out.pit_1_upper_clipped, expected.pit_theta_upper_clipped);
testCase.verifyEqual(out.pit_2_lower_clipped, expected.pit_alpha_lower_clipped);
testCase.verifyEqual(out.pit_2_upper_clipped, expected.pit_alpha_upper_clipped);
testCase.verifyEqual(out.pit_any_clip_fraction, expected.pit_any_clip_fraction, ...
    'AbsTol', 1e-15);

% Flags must agree with the reference's, since they are functions of the above.
testCase.verifyEqual(out.invalidity_flag, ...
    localLogical(expected.copula_stage_invalidity_flag));
testCase.verifyEqual(out.calibration_flag, ...
    localLogical(expected.pit_calibration_review_flag));

% The drawn rho must at least live on the grid and inside the support.
testCase.verifyTrue(all(abs(out.rho_e) <= 0.95));
end

function tf = localLogical(column)
%LOCALLOGICAL Reference CSVs carry R's TRUE/FALSE, which readtable delivers as
%text; accept that, numeric 0/1, or an already-logical column.
if islogical(column)
    tf = column(:);
elseif isnumeric(column)
    tf = column(:) ~= 0;
else
    tf = strcmpi(strtrim(string(column(:))), "TRUE");
end
end
