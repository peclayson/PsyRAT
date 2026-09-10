classdef TestGammaDodDynrelLsConversion < matlab.unittest.TestCase
    %TESTGAMMADODDYNRELLSCONVERSION Offline formula tests for the person-
    % specific dynamic NONCONCURRENT difference-of-differences kernels
    % (analyses 28/29): PSYRAT_DOD_COMPOSITE and
    % PSYRAT_GAMMA_DOD_COMPONENTS_LS.
    %
    % These are the MATLAB ports of the reference bundle's deterministic
    % checks (person_specific_dynamic_nonconcurrent_dod_modular,
    % tests/run_static_validation.R; provenance hashes in
    % tests/baselines/nonconcurrent_dod): check 1 (Gamma mean/SD mapping),
    % check 4 (the six-pair zero residual-covariance boundary), check 5 (the
    % signed six-covariance expansion equals the quadratic form), and the
    % HARMONIC half of check 6 (the exact-shared half is frozen in the
    % fixtures but out of toolbox scope - harmonic-only owner ruling,
    % 2026-08-12).
    %
    % Two structural claims ride along:
    %  - VETTED-KERNEL AGREEMENT: on the shared subspace (harmonic rule, a
    %    diagonal-only residual channel) the new composite kernel must
    %    reproduce the shipped static-DoD kernel PSYRAT_DODIFFREL draw for
    %    draw. The tie is this TEST, not a call - the shipped kernel stays
    %    untouched (its admissibility commentary is load-bearing for a
    %    channel structure this design does not have).
    %  - THE NONCONCURRENT BOUNDARY IS STRUCTURAL: the component assembler
    %    hard-rejects any residual-coupling vocabulary on its signature, and
    %    its pooled residual cross terms are BIT-identical to the
    %    mean-surface-only cross terms (cov_e_obs = 0), with a positive
    %    control proving that assertion could fail.
    %
    % Tolerance grammar (repo convention): AbsTol 1e-12 for same-arithmetic
    % routes evaluated in a different order; RelTol 1e-9 (with a 1e-12
    % absolute floor for near-zero entries) for assembled quantities checked
    % against an independent transcription of the reference's subtraction
    % forms; AbsTol 0 for deliberate bit-identity claims, each with a comment
    % saying why bit-identity is the correct expectation.

    methods (Test)

        function testGammaMeanSdMappingIdentity(testCase)
            %Reference check 1 at the bundle's own numbers: the scaled
            %chi-square / Gamma(shape,rate) mapping reproduces the mean and
            %SD exactly, and the derived degrees of freedom obey
            %nu = 2*(mu/sigma)^2 (so sigma = mu*sqrt(2/nu)).
            mu = 6.3; sigma = 1.7;
            shape = (mu / sigma)^2;
            rate = mu / sigma^2;
            testCase.verifyEqual(shape / rate, mu, 'AbsTol', 1e-12, ...
                'Gamma mean mapping failed.');
            testCase.verifyEqual(shape / rate^2, sigma^2, 'AbsTol', 1e-12, ...
                'Gamma variance mapping failed.');
            nu = 2 * (mu / sigma)^2;
            testCase.verifyEqual(mu * sqrt(2 / nu), sigma, 'AbsTol', 1e-12, ...
                'The derived-nu identity sigma = mu*sqrt(2/nu) failed.');
        end

        function testSignedExpansionMatchesQuadraticForm(testCase)
            %Reference check 5: the explicit six-covariance signed expansion,
            %the matrix quadratic form, and the composite kernel agree - first
            %on the bundle's exact S matrix, then on seeded random symmetric
            %draws. tol = 0 disables the tiny-negative floor so the raw
            %arithmetic is compared.
            w = [1 -1 -1 1]';
            S = [4 .3 .2 .1; .3 3 .1 .2; .2 .1 2 .25; .1 .2 .25 2.5];
            manual = sum(diag(S));
            for a = 1:3
                for b = (a+1):4
                    manual = manual + 2 * w(a) * w(b) * S(a,b);
                end
            end
            quad = w' * S * w;
            testCase.verifyEqual(manual, quad, 'AbsTol', 1e-12, ...
                'The signed six-covariance expansion must equal w''*S*w.');

            comp = psyrat_dod_composite('varcov', reshape(S, [1 4 4]), ...
                'rule', 'unscaled', 'tol', 0);
            testCase.verifyEqual(comp.value, quad, 'AbsTol', 1e-12, ...
                'The composite kernel must equal the quadratic form.');

            rng(24601, 'twister');
            nd = 200;
            V = zeros(nd, 4, 4);
            expect = zeros(nd, 1);
            for d = 1:nd
                A = randn(4);
                Sd = A + A';   %symmetric, not necessarily PSD - the identity
                               %holds for ANY symmetric matrix
                V(d, :, :) = Sd;
                expect(d) = w' * Sd * w;
            end
            comp = psyrat_dod_composite('varcov', V, 'rule', 'unscaled', ...
                'tol', 0);
            testCase.verifyEqual(comp.value, expect, 'AbsTol', 1e-12, ...
                'Draw-wise agreement with the quadratic form failed.');
        end

        function testHarmonicCountRuleMatchesExplicitFormula(testCase)
            %The harmonic half of reference check 6, at its exact numbers:
            %variance terms carry their own cell's count; each covariance pair
            %carries the PAIRWISE harmonic mean - including a same-count pair,
            %where the harmonic mean degenerates to the shared count, so a
            %global-harmonic misimplementation fails here.
            w = [1 -1 -1 1];
            n = [5 9 5 9];
            V = diag([4 3 2 2.5]);
            V(1,2) = .3; V(2,1) = .3;    %cross-count pair: h(5,9)
            V(1,3) = .2; V(3,1) = .2;    %same-count pair:  h(5,5) = 5
            h12 = 2 / (1/5 + 1/9);
            expected = 4/5 + 3/9 + 2/5 + 2.5/9 ...
                + 2 * w(1) * w(2) * .3 / h12 ...
                + 2 * w(1) * w(3) * .2 / 5;
            comp = psyrat_dod_composite('varcov', reshape(V, [1 4 4]), ...
                'rule', 'harmonic', 'n', n, 'tol', 0);
            testCase.verifyEqual(comp.value, expected, 'AbsTol', 1e-12, ...
                'Harmonic count scaling failed against the explicit formula.');

            %the harmonic rule without counts is an error, not a default
            testCase.verifyError(@() psyrat_dod_composite('varcov', ...
                reshape(V, [1 4 4]), 'rule', 'harmonic'), 'varargin:n');
        end

        function testVettedKernelAgreementOnSharedSubspace(testCase)
            %On the static kernel's own subspace - harmonic rule, residual
            %channel diagonal-only - the new composite must reproduce
            %PSYRAT_DODIFFREL draw for draw: the reliability and ICC draws
            %elementwise, and the summary points of every recorded quantity.
            %AbsTol 1e-12: same arithmetic assembled in a different order.
            %The fixture is non-vacuous by construction: all four counts
            %differ, and the between-trial block carries nonzero
            %off-diagonals, so the pairwise harmonic divisors genuinely vary
            %across pairs.
            rng(24601, 'twister');
            nd = 500;
            obs = [18 22 15 20];
            Lp = [0.55 0.30; 0.48 -0.35; -0.42 0.28; 0.36 0.22];
            Lt = [0.28 0.12; 0.24 -0.15; -0.20 0.14; 0.18 0.11];
            base_bp = Lp * Lp' + diag([0.35 0.42 0.30 0.38]);
            base_bt = Lt * Lt' + diag([0.06 0.075 0.055 0.07]);
            testCase.assertTrue(all(abs([base_bt(1,2) base_bt(1,3)]) > 1e-3), ...
                'Fixture must carry nonzero between-trial off-diagonals.');
            bp = zeros(nd, 4, 4); bt = zeros(nd, 4, 4);
            for d = 1:nd
                s = 1 + 0.3 * rand;   %PSD-preserving per-draw scaling
                bp(d, :, :) = base_bp * s;
                bt(d, :, :) = base_bt * (2 - s);
            end
            er = [log(1.5) log(1.2) log(1.8) log(1.35)] + 0.2 * randn(nd, 4);

            dep = psyrat_dodiffrel('bp', bp, 'bt', bt, 'er_var', er, ...
                'obs', obs, 'est', 'dep', 'ci', .95);
            gen = psyrat_dodiffrel('bp', bp, 'bt', bt, 'er_var', er, ...
                'obs', obs, 'est', 'gen', 'ci', .95);

            u   = psyrat_dod_composite('varcov', bp, 'rule', 'unscaled', 'tol', 0);
            btu = psyrat_dod_composite('varcov', bt, 'rule', 'unscaled', 'tol', 0);
            bts = psyrat_dod_composite('varcov', bt, 'rule', 'harmonic', ...
                'n', obs, 'tol', 0);
            ermat = zeros(nd, 4, 4);
            for q = 1:4
                ermat(:, q, q) = exp(er(:, q)).^2;
            end
            eru = psyrat_dod_composite('varcov', ermat, 'rule', 'unscaled', 'tol', 0);
            ers = psyrat_dod_composite('varcov', ermat, 'rule', 'harmonic', ...
                'n', obs, 'tol', 0);

            testCase.verifyEqual(u.value ./ (u.value + bts.value + ers.value), ...
                dep.rel_draws, 'AbsTol', 1e-12, ...
                'Dependability draws must match the vetted static kernel.');
            testCase.verifyEqual(u.value ./ (u.value + ers.value), ...
                gen.rel_draws, 'AbsTol', 1e-12, ...
                'Generalizability draws must match the vetted static kernel.');
            testCase.verifyEqual(u.value ./ (u.value + btu.value + eru.value), ...
                dep.icc_draws, 'AbsTol', 1e-12, ...
                'Dependability ICC draws must match the vetted static kernel.');
            testCase.verifyEqual(mean(u.value), dep.uni_var_pt, 'AbsTol', 1e-12, ...
                'Universe-score point estimate must match.');
            testCase.verifyEqual(mean(bts.value), dep.bt_scaled_var_pt, ...
                'AbsTol', 1e-12, 'Scaled trial contrast must match.');
            testCase.verifyEqual(mean(ers.value), dep.er_scaled_var_pt, ...
                'AbsTol', 1e-12, 'Scaled residual contrast must match.');
        end

        function testZeroResidualCovarianceBitIdentity(testCase)
            %Reference check 4, as structure. The assembled residual channel's
            %off-diagonals must be BIT-identical (AbsTol 0) to the
            %mean-surface-only cross terms: bit-identity is the correct
            %expectation because the assembler reaches them through the same
            %concurrent helper with cov_e_obs wired to literal zero, and this
            %test pins that wiring. The positive control shows a nonzero
            %coupling WOULD move the value, so the bit-identity assertion can
            %fail.
            fx = localOnefacetFixture(50);
            mats = psyrat_gamma_dod_components_ls('facet', 'trial', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);
            for a = 1:3
                for b = (a+1):4
                    cv0 = psyrat_gamma_crosscov_rescor( ...
                        fx.alpha(:,a), fx.sd_p(:,a), fx.sd_trl(:,a), ...
                        fx.alpha(:,b), fx.sd_p(:,b), fx.sd_trl(:,b), ...
                        fx.cor_p(:,a,b), fx.cor_trl(:,a,b), 0);
                    testCase.verifyEqual(mats.pi_e(:,a,b), cv0.cov_pi_e_obs, ...
                        'AbsTol', 0, sprintf(['Residual channel (%d,%d) must ', ...
                        'be bit-identical to the zero-coupling cross term.'], a, b));
                    %positive control: a coupled call differs, so the
                    %bit-identity above is falsifiable
                    cv1 = psyrat_gamma_crosscov_rescor( ...
                        fx.alpha(:,a), fx.sd_p(:,a), fx.sd_trl(:,a), ...
                        fx.alpha(:,b), fx.sd_p(:,b), fx.sd_trl(:,b), ...
                        fx.cor_p(:,a,b), fx.cor_trl(:,a,b), 0.37);
                    testCase.verifyTrue(any(mats.pi_e(:,a,b) ~= cv1.cov_pi_e_obs), ...
                        'The positive control must differ, or the identity is vacuous.');
                end
            end

            %the coupling vocabulary is rejected at the signature, never
            %silently swallowed
            for coupled = {'rho', 'rho_e', 'cov_e_obs', 'cov_res_e_obs', ...
                    'er_cov', 'rescor', 'copula'}
                testCase.verifyError(@() psyrat_gamma_dod_components_ls( ...
                    coupled{1}, 0.3, 'facet', 'trial'), ...
                    'psyrat_gamma_dod_components_ls:noresidualcoupling');
            end
        end

        function testComponentAssemblyParityOnefacet(testCase)
            %The assembler against an independent transcription of the
            %reference's SUBTRACTION forms (reliability.R:93-121; the
            %toolbox converters use algebraically identical product forms, so
            %agreement here is a genuine two-route check, not a mirror).
            fx = localOnefacetFixture(60);
            mats = psyrat_gamma_dod_components_ls('facet', 'trial', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);
            ref = localRefOnefacet(fx);
            for key = {'p', 'i', 'pi_e'}
                testCase.verifyEqual(mats.(key{1}), ref.(key{1}), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, sprintf( ...
                    'One-facet %s matrix disagrees with the reference forms.', ...
                    key{1}));
            end
        end

        function testComponentAssemblyParityTwofacet(testCase)
            %The two-facet assembler against the reference's subtraction forms
            %(reliability.R:123-153), with the toolbox facet-name mapping
            %(trl->i, occ->o, tid->pi, oid->po, to->io) exercised on
            %all-distinct values so a facet swap fails rather than passing by
            %symmetry.
            fx = localTwofacetFixture(40);
            mats = psyrat_gamma_dod_components_ls('facet', 'trial_occasion', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl, ...
                'sd_occ', fx.sd_occ, 'cor_occ', fx.cor_occ, ...
                'sd_tid', fx.sd_tid, 'cor_tid', fx.cor_tid, ...
                'sd_oid', fx.sd_oid, 'cor_oid', fx.cor_oid, ...
                'sd_to', fx.sd_to, 'cor_to', fx.cor_to);
            ref = localRefTwofacet(fx);
            for key = {'p', 'i', 'o', 'pi', 'po', 'io', 'pio_e'}
                testCase.verifyEqual(mats.(key{1}), ref.(key{1}), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, sprintf( ...
                    'Two-facet %s matrix disagrees with the reference forms.', ...
                    key{1}));
            end
        end

        function testCompositeNegativeFlagging(testCase)
            %The record-never-repair convention: tiny floating-point negatives
            %(within tol at the draw's own magnitude) become exactly zero with
            %no flag; material negatives keep their value and raise the flag.
            w = [1 -1 -1 1];
            V = zeros(3, 4, 4);
            %draw 1: comfortably positive
            V(1, :, :) = diag([1 1 1 1]);
            %draw 2: tiny negative composite (-4.6e-15, scale ~ 1)
            V(2, :, :) = diag(repmat(1e-16, 1, 4));
            V(2, 1, 2) = 2.5e-15; V(2, 2, 1) = 2.5e-15;
            %draw 3: materially negative composite (-0.5)
            V(3, :, :) = diag([.1 .1 .1 .1]);
            V(3, 1, 2) = 0.45; V(3, 2, 1) = 0.45;

            comp = psyrat_dod_composite('varcov', V, 'cvec', w, ...
                'rule', 'unscaled');
            testCase.verifyEqual(comp.value(1), 4, 'AbsTol', 1e-12);
            testCase.verifyEqual(comp.invalid_negative(1), 0);
            testCase.verifyEqual(comp.value(2), 0, 'AbsTol', 0, ...
                'A tiny negative must become exactly zero.');
            testCase.verifyEqual(comp.invalid_negative(2), 0, ...
                'A tiny negative is not an invalidity.');
            testCase.verifyEqual(comp.value(3), 0.4 - 0.9, 'AbsTol', 1e-12, ...
                'A material negative must be reported UNCHANGED.');
            testCase.verifyEqual(comp.invalid_negative(3), 1, ...
                'A material negative must raise the flag.');

            %non-finite input is an error, never a flagged pass-through
            Vbad = V; Vbad(2, 3, 3) = NaN;
            testCase.verifyError(@() psyrat_dod_composite('varcov', Vbad, ...
                'rule', 'unscaled'), 'psyrat_dod_composite:nonfinite');

            %THE SCALE RULE ITSELF: the magnitude scale uses the UNSCALED
            %absolute covariances (the reference's transcribed choice), not
            %the harmonic-scaled ones. This fixture separates the two: with
            %the correct scale (~5.0) the -0.05 composite sits exactly at the
            %tolerance boundary and is floored unflagged; under the
            %wrong-scale variant (max(1, 0.05) = 1) it would be flagged
            %material. Inspection-verified before; suite-pinned now.
            Vs = zeros(1, 4, 4);
            for q = 1:4
                Vs(1, q, q) = 1e-6;
            end
            Vs(1, 1, 2) = 2.5; Vs(1, 2, 1) = 2.5;
            comp = psyrat_dod_composite('varcov', Vs, 'cvec', w, ...
                'rule', 'harmonic', 'n', [100 100 100 100], 'tol', 0.01);
            testCase.verifyEqual(comp.value, 0, 'AbsTol', 0, ...
                ['The unscaled-covariance magnitude scale must floor this ' ...
                'boundary composite to exactly zero.']);
            testCase.verifyEqual(comp.invalid_negative, 0, ...
                ['A scaled-covariance magnitude scale would have flagged ' ...
                'this draw; the reference''s rule must not.']);
        end

        function testCorrelationDiagonalGuard(testCase)
            %Passing a covariance block where a correlation block belongs is
            %the likeliest scale error on the assembler signature.
            fx = localOnefacetFixture(5);
            bad = fx.cor_p .* 1.7;   %no longer unit-diagonal
            testCase.verifyError(@() psyrat_gamma_dod_components_ls( ...
                'facet', 'trial', 'alpha', fx.alpha, ...
                'log_sigma', fx.log_sigma, 'sd_p', fx.sd_p, 'cor_p', bad, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl), ...
                'psyrat_gamma_dod_components_ls:notcorrelation');
        end

    end
end

% ---------------------------------------------------------------------------
% Local fixtures: seeded synthetic pseudo-posterior draws with all-distinct
% per-cell values and at least one negative correlation per block, so index,
% sign, and facet swaps fail rather than passing by symmetry.
% ---------------------------------------------------------------------------

function fx = localOnefacetFixture(nd)
rng(4242, 'twister');
fx.alpha = [1.65 1.42 1.28 1.20] + 0.05 * randn(nd, 4);
fx.log_sigma = [0.48 0.31 0.22 0.15] + 0.04 * randn(nd, 4);
fx.sd_p = abs([0.26 0.235 0.21 0.185] + 0.02 * randn(nd, 4));
fx.sd_trl = abs([0.112 0.128 0.098 0.121] + 0.01 * randn(nd, 4));
fx.cor_p = localCorrDraws(nd, ...
    [0.55 0.30; 0.48 -0.35; -0.42 0.28; 0.36 0.22], [0.60 0.68 0.75 0.83]);
fx.cor_trl = localCorrDraws(nd, ...
    [0.50 0.25; -0.38 0.30; 0.27 -0.33; 0.36 0.22], [0.62 0.70 0.78 0.85]);
end

function fx = localTwofacetFixture(nd)
fx = localOnefacetFixture(nd);
rng(31337, 'twister');
fx.sd_occ = abs([0.075 0.062 0.088 0.070] + 0.008 * randn(nd, 4));
fx.sd_tid = abs([0.068 0.081 0.059 0.074] + 0.008 * randn(nd, 4));
fx.sd_oid = abs([0.059 0.047 0.066 0.052] + 0.008 * randn(nd, 4));
fx.sd_to = abs([0.052 0.066 0.044 0.058] + 0.008 * randn(nd, 4));
fx.cor_occ = localCorrDraws(nd, ...
    [0.44 -0.29; 0.35 0.26; -0.31 0.28; 0.24 0.37], [0.64 0.71 0.79 0.86]);
fx.cor_tid = localCorrDraws(nd, ...
    [0.41 0.33; 0.29 -0.36; 0.30 0.26; -0.28 0.24], [0.63 0.72 0.80 0.87]);
fx.cor_oid = localCorrDraws(nd, ...
    [-0.39 0.31; 0.27 0.34; 0.25 -0.29; 0.32 0.21], [0.65 0.73 0.81 0.88]);
fx.cor_to = localCorrDraws(nd, ...
    [0.37 0.28; -0.34 0.25; 0.29 0.35; 0.23 -0.27], [0.66 0.74 0.82 0.89]);
end

function C = localCorrDraws(nd, L, dvec)
%LOCALCORRDRAWS PSD correlation draws from fixed loadings with a small
%per-draw diagonal inflation, mirroring the frozen-fixture generator: valid
%correlation matrices at every draw, distinct across draws.
C = zeros(nd, 4, 4);
for d = 1:nd
    S = L * L' + diag(dvec + 0.002 * (d - 1));
    R = diag(1 ./ sqrt(diag(S))) * S * diag(1 ./ sqrt(diag(S)));
    C(d, :, :) = R;
end
end

% ---------------------------------------------------------------------------
% Independent transcription of the reference's SUBTRACTION forms
% (reliability.R:93-153). Deliberately not the toolbox's product-identity
% forms: two algebraically identical routes evaluated differently.
% ---------------------------------------------------------------------------

function ref = localRefOnefacet(fx)
nd = size(fx.alpha, 1);
z44 = zeros(nd, 4, 4);
ref = struct('p', z44, 'i', z44, 'pi_e', z44);
vt = fx.sd_p.^2 + fx.sd_trl.^2;
for q = 1:4
    B = exp(2 .* fx.alpha(:,q) + vt(:,q));
    p = B .* expm1(fx.sd_p(:,q).^2);
    i = B .* expm1(fx.sd_trl(:,q).^2);
    total = B .* expm1(vt(:,q));
    ref.p(:,q,q) = p;
    ref.i(:,q,q) = i;
    ref.pi_e(:,q,q) = (total - p - i) + exp(2 .* fx.log_sigma(:,q));
end
for a = 1:3
    for b = (a+1):4
        pref = exp(fx.alpha(:,a) + fx.alpha(:,b) + 0.5 .* (vt(:,a) + vt(:,b)));
        cp = fx.cor_p(:,a,b) .* fx.sd_p(:,a) .* fx.sd_p(:,b);
        ci = fx.cor_trl(:,a,b) .* fx.sd_trl(:,a) .* fx.sd_trl(:,b);
        p = pref .* expm1(cp);
        i = pref .* expm1(ci);
        total = pref .* expm1(cp + ci);
        vals = {p, i, total - p - i};
        keys = {'p', 'i', 'pi_e'};
        for k = 1:3
            ref.(keys{k})(:,a,b) = vals{k};
            ref.(keys{k})(:,b,a) = vals{k};
        end
    end
end
end

function ref = localRefTwofacet(fx)
%Reference component keys and the toolbox input mapping: i=trl, o=occ,
%pi=tid, po=oid, io=to.
nd = size(fx.alpha, 1);
z44 = zeros(nd, 4, 4);
ref = struct('p', z44, 'i', z44, 'o', z44, 'pi', z44, 'po', z44, ...
    'io', z44, 'pio_e', z44);
vp = fx.sd_p.^2; vi = fx.sd_trl.^2; vo = fx.sd_occ.^2;
vpi = fx.sd_tid.^2; vpo = fx.sd_oid.^2; vio = fx.sd_to.^2;
vt = vp + vi + vo + vpi + vpo + vio;
for q = 1:4
    B = exp(2 .* fx.alpha(:,q) + vt(:,q));
    comp = @(v) B .* expm1(v);
    p = comp(vp(:,q)); i = comp(vi(:,q)); o = comp(vo(:,q));
    pi_ = comp(vp(:,q) + vi(:,q) + vpi(:,q)) - p - i;
    po_ = comp(vp(:,q) + vo(:,q) + vpo(:,q)) - p - o;
    io_ = comp(vi(:,q) + vo(:,q) + vio(:,q)) - i - o;
    total = comp(vt(:,q));
    pio_mean = total - p - i - o - pi_ - po_ - io_;
    ref.p(:,q,q) = p; ref.i(:,q,q) = i; ref.o(:,q,q) = o;
    ref.pi(:,q,q) = pi_; ref.po(:,q,q) = po_; ref.io(:,q,q) = io_;
    ref.pio_e(:,q,q) = pio_mean + exp(2 .* fx.log_sigma(:,q));
end
for a = 1:3
    for b = (a+1):4
        pref = exp(fx.alpha(:,a) + fx.alpha(:,b) + 0.5 .* (vt(:,a) + vt(:,b)));
        cl = @(cor, sd) cor(:,a,b) .* sd(:,a) .* sd(:,b);
        cp = cl(fx.cor_p, fx.sd_p); ci = cl(fx.cor_trl, fx.sd_trl);
        co = cl(fx.cor_occ, fx.sd_occ); cpi = cl(fx.cor_tid, fx.sd_tid);
        cpo = cl(fx.cor_oid, fx.sd_oid); cio = cl(fx.cor_to, fx.sd_to);
        comp = @(c) pref .* expm1(c);
        p = comp(cp); i = comp(ci); o = comp(co);
        pi_ = comp(cp + ci + cpi) - p - i;
        po_ = comp(cp + co + cpo) - p - o;
        io_ = comp(ci + co + cio) - i - o;
        total = comp(cp + ci + co + cpi + cpo + cio);
        pio_mean = total - p - i - o - pi_ - po_ - io_;
        vals = {p, i, o, pi_, po_, io_, pio_mean};
        keys = {'p', 'i', 'o', 'pi', 'po', 'io', 'pio_e'};
        for k = 1:7
            ref.(keys{k})(:,a,b) = vals{k};
            ref.(keys{k})(:,b,a) = vals{k};
        end
    end
end
end
