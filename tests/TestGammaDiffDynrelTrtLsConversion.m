classdef TestGammaDiffDynrelTrtLsConversion < PsyRATTestBase
    %TESTGAMMADIFFDYNRELTRTLSCONVERSION The TWO-FACET half of the modular
    % cut-copula workflow's observed-scale conversion: the cross-covariance
    % wrapper (psyrat_gamma_crosscov_trt_rescor) and the typical-person surface
    % (psyrat_rel_diffdynrel_trt_gamma_ls) that analysis 20 reports.
    %
    % WHY THE FIXTURES MATTER MORE HERE THAN ANYWHERE ELSE IN THE TOOLBOX. The
    % accuracy suite is otherwise regression-only - its oracle duplicates the
    % production formulas, so green proves "unchanged", not "correct". The
    % one-facet modular design can additionally be replayed against the owner's
    % completed 151-person reference fits. The two-facet design CANNOT: the
    % reference bundle's two-facet runs were prepared against a simulated second
    % occasion and were never sampled (their chain directories are empty and no
    % completion manifest exists). The frozen CSVs under
    % tests/baselines/modular_cut_copula/ are therefore the ONLY externally
    % grounded check this code path has, and they come from the reference's own
    % R implementation rather than from MATLAB.
    %
    % WHAT THEY PIN AND WHAT THEY DO NOT. They pin the formulas - six
    % observed-scale components, their cross-event covariances, the pooled
    % pio_e residual channel, the mean-matched quadrature. They cannot pin the
    % Stan builder's facet-index assignment, because stage 3 consumes SDs by
    % name and never sees an index. That gap stays open.
    %
    % FACET NAMES CROSS OVER between the reference and the toolbox, and the
    % fixture values are deliberately all-distinct so that a mis-mapping fails
    % rather than passing by symmetry:
    %   reference i  (trial)            -> toolbox trl, varcomps arg s_t
    %   reference o  (occasion)         -> toolbox occ, varcomps arg s_o
    %   reference pi (person x trial)   -> toolbox tid, varcomps arg s_pt
    %   reference po (person x occasion)-> toolbox oid, varcomps arg s_po
    %   reference io (trial x occasion) -> toolbox to,  varcomps arg s_ot

    methods (Test)

        %% -- the cross-covariance wrapper ------------------------------------

        function testWrapperReproducesTheShippedTermsExactly(testCase)
            %All seven terms the shipped two-facet helper computes must be
            %BIT-identical, or the concurrent and non-concurrent two-facet paths
            %would silently disagree on components that are supposed to be shared.
            a = localCrossArgs();
            base = psyrat_gamma_crosscov_trt(a{:});
            got  = psyrat_gamma_crosscov_trt_rescor(a{:}, 0);

            for f = {'cov_p','cov_o','cov_t','cov_po','cov_pt','cov_ot', ...
                    'cov_res','cov_mu_total'}
                testCase.verifyEqual(got.(f{1}), base.(f{1}), 'AbsTol', 0, ...
                    sprintf('%s must pass through unchanged.', f{1}));
            end
        end

        function testPooledResidualChannelIsTheThreeWayTermPlusTheSuppliedCovariance(testCase)
            %cov_res_e_obs is the reference's covariance$pio_e =
            %covariance$pio_mean + residual_covariance (R/reliability.R:266-268).
            a = localCrossArgs();
            base = psyrat_gamma_crosscov_trt(a{:});
            cov_e = 0.017 * (1:numel(base.cov_res))';
            got = psyrat_gamma_crosscov_trt_rescor(a{:}, cov_e);

            testCase.verifyEqual(got.cov_res_e_obs, base.cov_res + cov_e, ...
                'AbsTol', 0);
            testCase.verifyEqual(got.cov_e_obs, cov_e, 'AbsTol', 0);
        end

        function testPersonByTrialIsNotFoldedIntoTheResidual(testCase)
            %THE CENTRAL SCIENTIFIC CLAIM OF THIS FILE. In the one-facet design
            %person x trial is confounded with error and belongs in the residual
            %slot; here it is separately estimable, carries its own D-study
            %divisor n_i (against the residual's n_i*n_o), and must stay out.
            %A future "reconciliation" with the one-facet sibling would show up
            %here.
            a = localCrossArgs();
            got = psyrat_gamma_crosscov_trt_rescor(a{:}, 0);

            testCase.verifyGreaterThan(max(abs(got.cov_pt)), 1e-8, ...
                'The person x trial cross term must be nonzero in this fixture.');
            testCase.verifyEqual(got.cov_res_e_obs, got.cov_res, 'AbsTol', 0, ...
                ['With no supplied residual covariance the pooled channel is '...
                'the three-way term ALONE - person x trial must not appear.']);
        end

        function testWrapperRejectsAMismatchedResidualLength(testCase)
            a = localCrossArgs();
            testCase.verifyError( ...
                @() psyrat_gamma_crosscov_trt_rescor(a{:}, [1;2]), ...
                'psyrat_gamma_crosscov_trt_rescor:size');
        end

        function testResidualCrossCovCancellation(testCase)
            %The shipped cov_res is an eight-term alternating inclusion-exclusion
            %in which EVERY term tends to +/-1 as the log-scale covariances
            %shrink while the sum is second order in them - a cancellation of
            %about 1.7e5 at a common log-covariance of 1e-3. This MEASURES the
            %resulting loss against an algebraically identical expm1 regrouping
            %(itself only better conditioned, not exact), so the decision to keep
            %the shipped form for bit-identity with analysis 10 rests on a number
            %rather than an assumption. See the wrapper's header.
            %
            %The ceilings sit just above the observed values (3.8e-12, 2.8e-9,
            %1.3e-7) so that those numbers are genuinely pinned and a regression
            %in either form fails here rather than passing a loose bound.
            %Realistic ERP designs sit near a log-scale covariance of 1e-3.
            scales = [1e-3 1e-4 1e-5];
            tols   = [5e-12 4e-9 2e-7];
            for k = 1:numel(scales)
                sc = scales(k);
                [shipped, exact] = localResidualBothWays(sc);
                relerr = abs(shipped - exact) / abs(exact);
                testCase.verifyLessThan(relerr, tols(k), sprintf( ...
                    ['At log-covariance %g the shipped and exact residual '...
                    'cross-covariances diverged by %g relative. If this grew, '...
                    'the eight-term form is no longer adequate.'], sc, relerr));
            end

            %And the regrouping really is the same quantity: at an ordinary
            %magnitude, where nothing cancels, the two agree to machine precision.
            [shipped, exact] = localResidualBothWays(1e-2);
            testCase.verifyEqual(shipped, exact, 'RelTol', 1e-11, ...
                'The expm1 regrouping must be algebraically identical.');
        end

        %% -- the reference fixtures: components and covariances ---------------

        function testComponentsMatchTheReferenceOracle(testCase)
            %Every observed-scale component, per event, against the reference's
            %own twofacet_mean_components. This is the externally grounded check.
            fx = localLoadStage3();

            for r = 1:numel(fx.rows)
                s = fx.rows(r);
                vcA = psyrat_gamma_varcomps_trt_ls(s.alpha_a, s.sp_a, s.so_a, ...
                    s.st_a, s.spo_a, s.spt_a, s.sot_a, log(s.sigma_a), 0);
                vcB = psyrat_gamma_varcomps_trt_ls(s.alpha_b, s.sp_b, s.so_b, ...
                    s.st_b, s.spo_b, s.spt_b, s.sot_b, log(s.sigma_b), 0);

                localCheck(testCase, vcA.sigma_p2,   s.var_theta_p,   'var_theta_p', r);
                localCheck(testCase, vcA.sigma_t2,   s.var_theta_i,   'var_theta_i', r);
                localCheck(testCase, vcA.sigma_o2,   s.var_theta_o,   'var_theta_o', r);
                localCheck(testCase, vcA.sigma_pt2,  s.var_theta_pi,  'var_theta_pi', r);
                localCheck(testCase, vcA.sigma_po2,  s.var_theta_po,  'var_theta_po', r);
                localCheck(testCase, vcA.sigma_ot2,  s.var_theta_io,  'var_theta_io', r);
                %the POOLED residual: three-way mean-surface remainder + the
                %participant's own sigma^2, which is what passing their own log
                %residual SD with s_sd = 0 produces.
                localCheck(testCase, vcA.sigma_res2, s.var_theta_pio_e, 'var_theta_pio_e', r);

                localCheck(testCase, vcB.sigma_p2,   s.var_alpha_p,   'var_alpha_p', r);
                localCheck(testCase, vcB.sigma_res2, s.var_alpha_pio_e, 'var_alpha_pio_e', r);
            end
        end

        function testCrossCovariancesMatchTheReferenceOracle(testCase)
            %The wrapper's seven cross terms and its pooled channel, against the
            %reference's twofacet_cross_components plus its residual_covariance.
            fx = localLoadStage3();

            for r = 1:numel(fx.rows)
                s = fx.rows(r);
                cc = psyrat_gamma_crosscov_trt_rescor( ...
                    s.alpha_a, s.sp_a, s.so_a, s.st_a, s.spo_a, s.spt_a, s.sot_a, ...
                    s.alpha_b, s.sp_b, s.so_b, s.st_b, s.spo_b, s.spt_b, s.sot_b, ...
                    s.cor_p, s.cor_o, s.cor_t, s.cor_po, s.cor_pt, s.cor_ot, ...
                    s.residual_covariance_person);

                localCheck(testCase, cc.cov_p,  s.cov_p,  'cov_p', r);
                localCheck(testCase, cc.cov_o,  s.cov_o,  'cov_o', r);
                localCheck(testCase, cc.cov_t,  s.cov_i,  'cov_t <- reference cov_i', r);
                localCheck(testCase, cc.cov_pt, s.cov_pi, 'cov_pt <- reference cov_pi', r);
                localCheck(testCase, cc.cov_po, s.cov_po, 'cov_po', r);
                localCheck(testCase, cc.cov_ot, s.cov_io, 'cov_ot <- reference cov_io', r);
                localCheck(testCase, cc.cov_res_e_obs, s.cov_pio_e, 'cov_res_e_obs', r);
            end
        end

        function testMeanMatchedQuadraturePointMatchesTheReference(testCase)
            %The participant's expected mean marginalizes over ALL FIVE non-person
            %facets - not the trial facet alone, which is the one-facet rule - and
            %carries their realized location effect. nu follows from it, and the
            %residual covariance is the quadrature rescaled by the two means.
            fx = localLoadStage3();

            for r = 1:numel(fx.rows)
                s = fx.rows(r);
                localCheck(testCase, s.E_a_recomputed, s.expected_theta_person, ...
                    'expected_theta_person', r);
                localCheck(testCase, s.E_b_recomputed, s.expected_alpha_person, ...
                    'expected_alpha_person', r);

                nu_a = 2 * (s.expected_theta_person / s.sigma_a)^2;
                nu_b = 2 * (s.expected_alpha_person / s.sigma_b)^2;
                localCheck(testCase, nu_a, s.nu_theta, 'nu_theta', r);
                localCheck(testCase, nu_b, s.nu_alpha, 'nu_alpha', r);

                q = psyrat_gamma_copula_rescov(nu_a, nu_b, s.rho_e);
                got = s.expected_theta_person * s.expected_alpha_person * q.value;
                testCase.verifyEqual(got, s.residual_covariance_person, ...
                    'RelTol', 1e-7, sprintf( ...
                    'row %d: quadrature residual covariance', r));
            end
        end

        %% -- the surface -------------------------------------------------------

        function testSurfaceReducesToTheOneFacetSurfaceWithTheCopulaChannelOn(testCase)
            %THE END-TO-END REDUCTION, on the path analysis 20 actually takes.
            %With the occasion facets AND the person x trial facet zeroed, and
            %nocc = 1, the two-facet coefficient IS the one-facet coefficient -
            %but only because the pieces recombine, not because the residual
            %slots match. They do not: sigma_res2 here is dispersion alone while
            %the one-facet sigma_pi_e2 additionally carries the mean-surface
            %person x trial term, and psyrat_gamma_varcomps_trt_ls's header
            %records that s_pt = 0 is REQUIRED and is not implied by zeroing the
            %occasion facets. What makes the coefficients agree is
            %   sigma_res2 + sigma_pt2 = sigma_pi_e2   and   cov_pt + cov_res = cov_pimean,
            %with the two-facet's separate n_i and n_i*n_o divisors coinciding at
            %nocc = 1. Break any one of those and this fails.
            [args1, args2] = localMatchedOneAndTwoFacetArgs();
            rho = 0.35 * ones(localNDraws(), 1);

            one = psyrat_rel_diffdynrel_gamma_ls(args1{:}, 'rho', rho);
            two = psyrat_rel_diffdynrel_trt_gamma_ls(args2{:}, 'rho', rho);

            testCase.verifyEqual(two.G.pt, one.G.pt, 'RelTol', 1e-10, ...
                'The reduced two-facet generalizability surface must match the one-facet.');
            testCase.verifyEqual(two.D.pt, one.D.pt, 'RelTol', 1e-10);
            testCase.verifyEqual(two.ICCg.pt, one.ICCg.pt, 'RelTol', 1e-10);
            testCase.verifyEqual(two.ICCd.pt, one.ICCd.pt, 'RelTol', 1e-10);
        end

        function testTheReportedContrastSDsDifferByThePersonByTrialTerm(testCase)
            %AND THEY ARE SUPPOSED TO. sige_pt is the residual CONTRAST SD, and
            %"residual" does not name the same set of effects in the two designs:
            %the one-facet slot carries person x trial, this one breaks it out.
            %So the reduction above holds for the COEFFICIENTS while the reported
            %contrast SDs legitimately differ, by exactly the person x trial
            %contrast variance:
            %
            %   sige_1f^2 - sige_2f^2 = sigma_pt2_A + sigma_pt2_B - 2*cov_pt
            %
            %Asserting equality here instead would look like a stronger test and
            %would in fact be asserting that the two designs decompose the same
            %way, which is the one thing this file exists to deny.
            [args1, args2] = localMatchedOneAndTwoFacetArgs();
            rho = 0.35 * ones(localNDraws(), 1);

            one = psyrat_rel_diffdynrel_gamma_ls(args1{:}, 'rho', rho);
            two = psyrat_rel_diffdynrel_trt_gamma_ls(args2{:}, 'rho', rho);

            testCase.verifyLessThan(two.sige_pt, one.sige_pt, ...
                'Breaking person x trial out must SHRINK the residual contrast.');

            %Rebuild both contrast VARIANCES per draw from the components, and
            %check each surface reports the mean of its own square root. Both
            %sige_pt values average a square root over draws, so there is no
            %algebraic identity between them to assert - what is pinned instead
            %is that each names the set of effects it claims to, and that the
            %per-draw variances differ by exactly the person x trial contrast.
            p = localContrastPieces(args2, rho);

            testCase.verifyEqual(one.sige_pt, mean(sqrt(p.varOne),1), ...
                'RelTol', 1e-9, 'one-facet sige_pt must be the pi_e contrast.');
            testCase.verifyEqual(two.sige_pt, mean(sqrt(p.varTwo),1), ...
                'RelTol', 1e-9, 'two-facet sige_pt must be the pio_e contrast.');
            testCase.verifyEqual(p.varOne - p.varTwo, p.varPersonByTrial, ...
                'RelTol', 1e-9, ...
                'The per-draw gap must be the person x trial contrast variance.');
        end

        function testTheNonConcurrentPathsDisagreeInTheDocumentedDirection(testCase)
            %Without rho the two surfaces do NOT agree, and that is a property of
            %the one-facet path rather than a defect here: it sets wp_cov = 0,
            %which drops not only the residual covariance but also the
            %model-implied mean-surface cross term EY_A*EY_B*expm1(c_p)*expm1(c_i)
            %- a measured approximation the owner accepted on 2026-08-05
            %(gamma_scale_submodel_decision.md section I.4). The two-facet path
            %cannot make that omission, because the same quantity arrives as the
            %person x trial covariance BLOCK rather than through wp_cov.
            %
            %Dropping a positive cross term inflates the error, so the one-facet
            %surface must come out LOWER. Pinned so that if the one-facet
            %approximation is ever revisited, this fails and says why.
            [args1, args2] = localMatchedOneAndTwoFacetArgs();

            one = psyrat_rel_diffdynrel_gamma_ls(args1{:});
            two = psyrat_rel_diffdynrel_trt_gamma_ls(args2{:});

            testCase.verifyGreaterThan(two.G.pt, one.G.pt, ...
                ['Without rho the two-facet surface retains the mean-surface '...
                'cross term the one-facet path drops, so it must be higher.']);
        end

        function testZeroCopulaCorrelationIsTheNonConcurrentReduction(testCase)
            %rho = 0 must give exactly the omitted-rho surface: the quadrature
            %short-circuits and the pooled channel collapses to the mean-surface
            %three-way term.
            args = localSurfaceArgs();
            withZero = psyrat_rel_diffdynrel_trt_gamma_ls(args{:}, ...
                'rho', zeros(localNDraws(),1));
            without = psyrat_rel_diffdynrel_trt_gamma_ls(args{:});

            testCase.verifyEqual(withZero.G.pt, without.G.pt, 'AbsTol', 0);
            testCase.verifyEqual(withZero.sige_pt, without.sige_pt, 'AbsTol', 0);
        end

        function testPositiveResidualCouplingShrinksTheDifferenceError(testCase)
            %Direction check, and the reason the copula channel exists: events
            %measured on the same trials share residual noise, which cancels in
            %the difference. A positive correlation must therefore REDUCE the
            %residual contrast SD and raise reliability. A sign error anywhere in
            %the chain flips this.
            args = localSurfaceArgs();
            nd = localNDraws();
            none = psyrat_rel_diffdynrel_trt_gamma_ls(args{:}, 'rho', zeros(nd,1));
            pos  = psyrat_rel_diffdynrel_trt_gamma_ls(args{:}, 'rho', 0.5*ones(nd,1));
            neg  = psyrat_rel_diffdynrel_trt_gamma_ls(args{:}, 'rho', -0.5*ones(nd,1));

            testCase.verifyLessThan(pos.sige_pt, none.sige_pt);
            testCase.verifyGreaterThan(neg.sige_pt, none.sige_pt);
            testCase.verifyGreaterThan(pos.G.pt, none.G.pt);
            testCase.verifyLessThan(neg.G.pt, none.G.pt);
        end

        function testSurfaceRejectsAMisalignedOrInadmissibleRho(testCase)
            %Draw alignment is the caller's job: the copula correlation exists
            %only for the retained stage-1 draws, so a short rho must be rejected
            %rather than recycled against unrelated margins.
            args = localSurfaceArgs();
            nd = localNDraws();

            testCase.verifyError( ...
                @() psyrat_rel_diffdynrel_trt_gamma_ls(args{:}, 'rho', ones(nd-1,1)*0.2), ...
                'varargin:rho');
            testCase.verifyError( ...
                @() psyrat_rel_diffdynrel_trt_gamma_ls(args{:}, 'rho', ones(nd,1)), ...
                'varargin:rho');
            testCase.verifyError( ...
                @() psyrat_rel_diffdynrel_trt_gamma_ls(args{:}, 'rho', nan(nd,1)), ...
                'varargin:rho');
        end

        function testOneFacetSurfaceRejectsAnInadmissibleRhoToo(testCase)
            %The one-facet surface was the single entry point in this
            %workstream without the |rho| < 1 / finiteness rejection its
            %two-facet twin and both ssrel calculators carry; a bad rho died
            %mid-grid inside the quadrature kernel instead (pre-merge review
            %finding, closed here).
            [args1, ~] = localMatchedOneAndTwoFacetArgs();
            nd = localNDraws();

            testCase.verifyError( ...
                @() psyrat_rel_diffdynrel_gamma_ls(args1{:}, 'rho', ones(nd,1)), ...
                'varargin:rho');
            testCase.verifyError( ...
                @() psyrat_rel_diffdynrel_gamma_ls(args1{:}, 'rho', nan(nd,1)), ...
                'varargin:rho');
            testCase.verifyError( ...
                @() psyrat_rel_diffdynrel_gamma_ls(args1{:}, 'rho', ones(nd-1,1)*0.2), ...
                'varargin:rho');
        end

        function testSurfaceRejectsNonFiniteDesignCounts(testCase)
            %NaN < 1 is FALSE in MATLAB, so a bare magnitude check would let a
            %NaN n' through and produce an all-NaN surface with no diagnostic.
            badObs = localReplace(localSurfaceArgs(), 'obs', NaN);
            testCase.verifyError( ...
                @() psyrat_rel_diffdynrel_trt_gamma_ls(badObs{:}), 'varargin:obs');

            badOcc = localReplace(localSurfaceArgs(), 'nocc', NaN);
            testCase.verifyError( ...
                @() psyrat_rel_diffdynrel_trt_gamma_ls(badOcc{:}), 'varargin:nocc');
        end

        function testSurfaceRejectsAPerEventFacetGivenOnce(testCase)
            %Every crossed facet SD is per-event. A single column would broadcast
            %silently and give both events the same facet SD, which is a
            %different model rather than an input error.
            args = localSurfaceArgs();
            bad = localReplace(args, 'sd_oid', ones(localNDraws(),1)*0.07);
            testCase.verifyError( ...
                @() psyrat_rel_diffdynrel_trt_gamma_ls(bad{:}), 'varargin:facetsd');
        end

    end
end

%% ------------------------------------------------------------------------
%% helpers
%% ------------------------------------------------------------------------

function a = localCrossArgs()
%LOCALCROSSARGS Twenty draw-vector arguments for the two-facet cross-covariance
%helpers. Every SD and correlation is DISTINCT so an argument-order slip cannot
%pass by symmetry.
n = 6;
step = @(v) v + (0:n-1)' * (v/50);
a = {step(1.62), step(0.301), step(0.104), step(0.121), step(0.082), ...
     step(0.093), step(0.061), ...
     step(1.44), step(0.283), step(0.096), step(0.113), step(0.074), ...
     step(0.086), step(0.052), ...
     step(0.42), step(0.24), step(0.31), step(-0.14), step(0.19), step(0.27)};
end

function [shipped, exact] = localResidualBothWays(sc)
%LOCALRESIDUALBOTHWAYS The shipped eight-term residual cross-covariance and the
%algebraically identical expm1 regrouping, at a target log-scale covariance sc.
%
%   P^-1*cov_res = ABCDEF - ABD - ACE - BCF + A + B + C - 1
%                = A*B*D*expm1(c+e+f) - A*expm1(c+e) - B*expm1(c+f) + expm1(c)
%
%with A = exp(a) for the person log-covariance a, B = exp(b) occasion, C = exp(c)
%trial, D/E/F the three two-way interactions.
al1 = 1.6; al2 = 1.4;
sp1 = .30; sp2 = .28; so1 = .10; so2 = .09; st1 = .12; st2 = .11;
spo1 = .08; spo2 = .07; spt1 = .09; spt2 = .08; sot1 = .06; sot2 = .05;
cp = sc/(sp1*sp2); co = sc/(so1*so2); ct = sc/(st1*st2);
cpo = sc/(spo1*spo2); cpt = sc/(spt1*spt2); cot = sc/(sot1*sot2);

b = psyrat_gamma_crosscov_trt(al1,sp1,so1,st1,spo1,spt1,sot1, ...
                              al2,sp2,so2,st2,spo2,spt2,sot2, ...
                              cp,co,ct,cpo,cpt,cot);
shipped = b.cov_res;

a = cp*sp1*sp2; bb = co*so1*so2; c = ct*st1*st2;
d = cpo*spo1*spo2; e = cpt*spt1*spt2; f = cot*sot1*sot2;
V1 = sp1^2+so1^2+st1^2+spo1^2+spt1^2+sot1^2;
V2 = sp2^2+so2^2+st2^2+spo2^2+spt2^2+sot2^2;
P = exp(al1 + 0.5*V1) * exp(al2 + 0.5*V2);
exact = P * ( exp(a)*exp(bb)*exp(d)*expm1(c+e+f) ...
            - exp(a)*expm1(c+e) - exp(bb)*expm1(c+f) + expm1(c) );
end

function localCheck(testCase, got, want, name, row)
%LOCALCHECK Fixture comparison with a message that names the row and quantity.
%RelTol 1e-9 throughout: the reference writes its CSVs at full double precision,
%so the only gap is accumulated floating-point ordering, not method.
testCase.verifyEqual(got, want, 'RelTol', 1e-9, ...
    sprintf('fixture row %d: %s', row, name));
end

function n = localNDraws()
n = 8;
end

function args = localSurfaceArgs()
%LOCALSURFACEARGS A small but fully general two-facet surface call: one
%dimension, unequal per-event facet SDs, nonzero correlations on every facet.
n = localNDraws();
step = @(v) v + (0:n-1)' * (v/40);
corid = zeros(n,4,4);
for q = 1:n
    corid(q,:,:) = eye(4);
    corid(q,1,2) = 0.45; corid(q,2,1) = 0.45;
end
args = {'b',[step(1.62) step(1.44)], 'b_sigma',[step(-0.5) step(-0.62)], ...
    'b_dim',[step(0.10) step(0.08)], 'b_sigma_dim',[step(-0.05) step(-0.04)], ...
    'sd_id',[step(0.301) step(0.283) step(0.201) step(0.184)], ...
    'sd_trl',[step(0.121) step(0.113)], 'sd_occ',[step(0.104) step(0.096)], ...
    'sd_tid',[step(0.093) step(0.086)], 'sd_oid',[step(0.082) step(0.074)], ...
    'sd_to',[step(0.061) step(0.052)], 'cor_id',corid, ...
    'cor_t',step(0.31), 'cor_o',step(0.24), 'cor_pt',step(0.19), ...
    'cor_po',step(-0.14), 'cor_ot',step(0.27), ...
    'ndim',1, 'z1',linspace(-1,1,4), 'obs',25, 'nocc',2, 'reltype',3, 'CI',.95};
end

function [args1, args2] = localMatchedOneAndTwoFacetArgs()
%LOCALMATCHEDONEANDTWOFACETARGS One-facet arguments and the two-facet arguments
%that must reproduce them: occasion, person x occasion, trial x occasion AND
%person x trial all zeroed, nocc = 1, reltype 3 (both facets random).
n = localNDraws();
step = @(v) v + (0:n-1)' * (v/40);
z = zeros(n,1);
corid = zeros(n,4,4);
for q = 1:n
    corid(q,:,:) = eye(4);
    corid(q,1,2) = 0.45; corid(q,2,1) = 0.45;
end
common = {'b',[step(1.62) step(1.44)], 'b_sigma',[step(-0.5) step(-0.62)], ...
    'b_dim',[step(0.10) step(0.08)], 'b_sigma_dim',[step(-0.05) step(-0.04)], ...
    'sd_id',[step(0.301) step(0.283) step(0.201) step(0.184)], ...
    'sd_trl',[step(0.121) step(0.113)], 'cor_id',corid, ...
    'ndim',1, 'z1',linspace(-1,1,4), 'CI',.95};

args1 = [common, {'cor_i',step(0.31), 'obs',[25 25]}];
args2 = [common, {'sd_occ',[z z], 'sd_tid',[z z], 'sd_oid',[z z], ...
    'sd_to',[z z], 'cor_t',step(0.31), 'cor_o',z, 'cor_pt',z, ...
    'cor_po',z, 'cor_ot',z, 'obs',25, 'nocc',1, 'reltype',3}];
end

function p = localContrastPieces(args2, rho)
%LOCALCONTRASTPIECES Rebuilds, per grid point and per draw, the two residual
%CONTRAST variances the two designs report and the person x trial contrast that
%separates them. Independent of both surface functions on purpose - it is built
%from the component converters directly, so it checks what each sige_pt names
%rather than restating either implementation.
b     = localArg(args2,'b');      bsg = localArg(args2,'b_sigma');
bd    = localArg(args2,'b_dim');  bsd = localArg(args2,'b_sigma_dim');
sdid  = localArg(args2,'sd_id');  sdt = localArg(args2,'sd_trl');
cor_t = localArg(args2,'cor_t');  z1  = localArg(args2,'z1');
corid = localArg(args2,'cor_id');
cor_p = corid(:,1,2);
n = size(b,1);
zc = zeros(n,1);

varOne = zeros(n,numel(z1)); varTwo = varOne; varPt = varOne;
for a = 1:numel(z1)
    alpha_a = b(:,1) + bd(:,1)*z1(a);   alpha_b = b(:,2) + bd(:,2)*z1(a);
    lsig_a  = bsg(:,1) + bsd(:,1)*z1(a); lsig_b = bsg(:,2) + bsd(:,2)*z1(a);
    sp_a = sdid(:,1); sp_b = sdid(:,2); ssd_a = sdid(:,3); ssd_b = sdid(:,4);
    st_a = sdt(:,1);  st_b = sdt(:,2);

    E_a = exp(alpha_a + 0.5*st_a.^2);  E_b = exp(alpha_b + 0.5*st_b.^2);
    cov_e = zeros(n,1);
    for q = 1:n
        qq = psyrat_gamma_copula_rescov(2*(E_a(q)/exp(lsig_a(q)))^2, ...
                                        2*(E_b(q)/exp(lsig_b(q)))^2, rho(q));
        cov_e(q) = E_a(q)*E_b(q)*qq.value;
    end

    v1a = psyrat_gamma_varcomps_ls(alpha_a, sp_a, st_a, lsig_a, ssd_a);
    v1b = psyrat_gamma_varcomps_ls(alpha_b, sp_b, st_b, lsig_b, ssd_b);
    c1  = psyrat_gamma_crosscov_rescor(alpha_a,sp_a,st_a, ...
                                       alpha_b,sp_b,st_b, cor_p, cor_t, cov_e);
    v2a = psyrat_gamma_varcomps_trt_ls(alpha_a,sp_a,zc,st_a,zc,zc,zc,lsig_a,ssd_a);
    v2b = psyrat_gamma_varcomps_trt_ls(alpha_b,sp_b,zc,st_b,zc,zc,zc,lsig_b,ssd_b);
    c2  = psyrat_gamma_crosscov_trt_rescor(alpha_a,sp_a,zc,st_a,zc,zc,zc, ...
                                           alpha_b,sp_b,zc,st_b,zc,zc,zc, ...
                                           cor_p,zc,cor_t,zc,zc,zc, cov_e);

    varOne(:,a) = max(v1a.sigma_pi_e2 + v1b.sigma_pi_e2 - 2*c1.cov_pi_e_obs, 0);
    varTwo(:,a) = max(v2a.sigma_res2 + v2b.sigma_res2 - 2*c2.cov_res_e_obs, 0);
    varPt(:,a)  = v2a.sigma_pt2 + v2b.sigma_pt2 - 2*c2.cov_pt;
end
p = struct('varOne',varOne,'varTwo',varTwo,'varPersonByTrial',varPt);
end

function v = localArg(args, name)
%LOCALARG Reads one value out of a name-value argument cell array.
idx = find(strcmp(args, name), 1);
assert(~isempty(idx), 'localArg: %s not present', name);
v = args{idx+1};
end

function args = localReplace(args, name, value)
%LOCALREPLACE Swaps one name-value pair in an argument cell array.
idx = find(strcmp(args, name), 1);
assert(~isempty(idx), 'localReplace: %s not present', name);
args{idx+1} = value;
end

function fx = localLoadStage3()
%LOCALLOADSTAGE3 Rebuilds every stage-3 input from the frozen reference CSVs.
%
%The reference stores its draw parameters in LONG format precisely so that
%reassembling them here cannot depend on a column-ordering convention, and the
%stage-3-only population parameters are in their own file because adding them to
%the stage-2 draws would have rewritten a fixture this test has no business
%changing.

here = fileparts(mfilename('fullpath'));
rd = @(f) readtable(fullfile(here,'baselines','modular_cut_copula',f), ...
    'VariableNamingRule','preserve');

%The person table (trial and occasion counts per participant) is emitted
%alongside these but is not read here: this file checks the conversion, which
%carries no D-study counts. The per-participant calculator's tests consume it.
meta   = rd('modular_cut_copula_twofacet_meta.csv');
long   = rd('modular_cut_copula_twofacet_draws.csv');
pars   = rd('modular_cut_copula_twofacet_stage3_params.csv');
exp3   = rd('modular_cut_copula_twofacet_stage3_expected.csv');

%component columns are identical across the two D-study designs (design changes
%only the divisors), so one design's rows carry every quantity checked here.
exp3 = exp3(strcmp(exp3.design,'actual'), :);

K = meta.K(1);
mu_off  = [meta.mu_offset_1(1) meta.mu_offset_2(1)];
sig_off = [meta.log_sigma_offset_1(1) meta.log_sigma_offset_2(1)];

%NOTE these are plain local functions taking the table explicitly, NOT nested
%functions. A nested function shares the parent's workspace, so an internal mask
%named `m` silently overwrote the `m` used as the event index in the loop below
%and the lookups started matching on a logical vector.
draw_par = @(name,i1,i2,d) localLookup(long, name, i1, i2, d);
pop_par  = @(name,ev,d)    localLookupPar(pars, name, ev, d);

rows = struct([]);
for r = 1:height(exp3)
    d = exp3.posterior_row(r);
    p = exp3.p_index(r);
    zz = [exp3.Dim1_z(r); exp3.Dim2_z(r); exp3.Dim1_z(r)*exp3.Dim2_z(r)];
    zz = zz(1:K);

    s = struct();
    for m = 1:2
        alpha = mu_off(m) + draw_par('alpha_mu',m,NaN,d);
        lsig  = sig_off(m) + draw_par('alpha_sigma',m,NaN,d);
        for k = 1:K
            alpha = alpha + draw_par('beta_mu',m,k,d) * zz(k);
            lsig  = lsig  + draw_par('beta_sigma',m,k,d) * zz(k);
        end
        %the person's realized effects: slots 1-2 location, 3-4 log residual SD
        mu_re  = draw_par('u_p_joint',p,m,d);
        sig_re = draw_par('u_p_joint',p,2+m,d);

        %reference facet name -> toolbox varcomps argument (see the class header)
        sp  = pop_par('sd_p_joint',m,d);
        st  = pop_par('sd_i',m,d);     % trial
        so  = pop_par('sd_o',m,d);     % occasion
        spt = pop_par('sd_pi',m,d);    % person x trial
        spo = pop_par('sd_po',m,d);    % person x occasion
        sot = pop_par('sd_io',m,d);    % trial x occasion
        v_other = st^2 + so^2 + spt^2 + spo^2 + sot^2;

        tag = {'a','b'}; t = tag{m};
        s.(['alpha_' t]) = alpha;
        s.(['sp_'  t]) = sp;  s.(['st_'  t]) = st;  s.(['so_'  t]) = so;
        s.(['spt_' t]) = spt; s.(['spo_' t]) = spo; s.(['sot_' t]) = sot;
        s.(['sigma_' t]) = exp(lsig + sig_re);
        s.(['E_' t '_recomputed']) = exp(alpha + mu_re + 0.5*v_other);
    end

    s.cor_p  = pop_par('cor_p_mu',NaN,d);
    s.cor_t  = pop_par('cor_i',NaN,d);
    s.cor_o  = pop_par('cor_o',NaN,d);
    s.cor_pt = pop_par('cor_pi',NaN,d);
    s.cor_po = pop_par('cor_po',NaN,d);
    s.cor_ot = pop_par('cor_io',NaN,d);
    s.rho_e  = pop_par('rho_e',NaN,d);

    for f = {'p','i','o','pi','po','io','pio_e'}
        s.(['var_theta_' f{1}]) = exp3.(['var_theta_' f{1}])(r);
        s.(['var_alpha_' f{1}]) = exp3.(['var_alpha_' f{1}])(r);
        s.(['cov_' f{1}]) = exp3.(['cov_theta_alpha_' f{1}])(r);
    end
    s.expected_theta_person = exp3.expected_theta_person(r);
    s.expected_alpha_person = exp3.expected_alpha_person(r);
    s.nu_theta = exp3.nu_theta_at_person_expected_mean(r);
    s.nu_alpha = exp3.nu_alpha_at_person_expected_mean(r);
    s.residual_covariance_person = exp3.residual_covariance_person(r);

    if isempty(rows); rows = s; else; rows(end+1) = s; end %#ok<AGROW>
end
fx.rows = rows;
end

function v = localLookup(long, name, i1, i2, d)
%LOCALLOOKUP One value out of the long-format draw table.
sel = strcmp(long.param,name) & long.draw == d & long.idx1 == i1;
if ~isnan(i2); sel = sel & long.idx2 == i2; end
v = long.value(sel);
assert(isscalar(v), 'localLookup: %s[%g,%g] draw %g matched %d rows', ...
    name, i1, i2, d, numel(v));
end

function v = localLookupPar(pars, name, ev, d)
%LOCALLOOKUPPAR One value out of the stage-3 population-parameter table.
sel = strcmp(pars.param,name) & pars.draw == d;
if ~isnan(ev); sel = sel & pars.event == ev; end
v = pars.value(sel);
assert(isscalar(v), 'localLookupPar: %s (event %g) draw %g matched %d rows', ...
    name, ev, d, numel(v));
end
