classdef TestGaussDodComponentsConversion < matlab.unittest.TestCase
    %TESTGAUSSDODCOMPONENTSCONVERSION Offline formula tests for the GAUSSIAN
    % component assembler of the person-specific dynamic NONCONCURRENT
    % difference-of-differences designs (analyses 28/29, family='gaussian'):
    % PSYRAT_DOD_COMPONENTS. The shared composite kernel PSYRAT_DOD_COMPOSITE
    % is pinned family-agnostically by TestGammaDodDynrelLsConversion and is
    % not re-pinned here.
    %
    % These tests pin the two IDENTITY-LINK THEOREMS the Gaussian arm rests on
    % (SCIENTIFIC_FORMULA_AUDIT.md section 26):
    %  - THE RESIDUAL CHANNEL IS EXACTLY DIAGONAL. Under the identity link
    %    there is no induced mean-surface interaction, and the six
    %    observation-level residual covariances are zero by estimand, so the
    %    pooled residual cross channel vanishes IDENTICALLY - unlike the gamma
    %    assembler, whose pi_e off-diagonals stay nonzero at zero coupling
    %    (the log-link mean-surface term). The gamma assembler on the same
    %    numbers is the positive control proving the assertion can fail.
    %  - SIGNAL CHANNELS ARE MEAN- AND Z-INVARIANT. Location parameters enter
    %    no second moment: perturbing alpha moves only mats.expected, and
    %    perturbing log_sigma moves only the residual diagonal.
    %
    % Tolerance grammar (repo convention): AbsTol 0 for bit-identity claims,
    % each with a comment saying why bit-identity is the correct expectation;
    % AbsTol 1e-12 for same-arithmetic routes evaluated in a different order.
    %
    % The vetted-kernel tie (the Gaussian stage-3 against the shipped static
    % kernel PSYRAT_DODIFFREL, full-space for analysis 28) lands with the
    % stage-3 calculators; it needs them, not only this assembler.

    methods (Test)

        function testResidualChannelExactlyDiagonal(testCase)
            %The defining Gaussian theorem, as a bit-identity. AbsTol 0 is the
            %correct expectation because the assembler must write LITERAL
            %zeros into the residual off-diagonals - never compute them as a
            %difference that happens to cancel.
            fx = localOnefacetFixture(50);
            mats = psyrat_dod_components('facet', 'trial', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);
            for a = 1:3
                for b = (a+1):4
                    testCase.verifyEqual(mats.pi_e(:,a,b), ...
                        zeros(size(fx.alpha,1),1), 'AbsTol', 0, sprintf( ...
                        ['Residual channel (%d,%d) must be EXACTLY zero ', ...
                        'under the identity link.'], a, b));
                end
            end

            %positive control: the GAMMA assembler on the same numbers keeps
            %nonzero pi_e off-diagonals (the log-link mean-surface term), so
            %the bit-zero assertion above is falsifiable and the theorem is a
            %genuine family difference, not a vacuous claim. The gamma
            %assembler reads alpha/sd_p as log-scale quantities; for this
            %control they are just numbers, and positive alpha keeps them in
            %its ordinary range.
            gx = localOnefacetFixture(50);
            gx.alpha = abs(gx.alpha) ./ 4;   %log-scale-plausible magnitudes
            gmats = psyrat_gamma_dod_components_ls('facet', 'trial', ...
                'alpha', gx.alpha, 'log_sigma', gx.log_sigma, ...
                'sd_p', gx.sd_p ./ 10, 'cor_p', gx.cor_p, ...
                'sd_trl', gx.sd_trl ./ 10, 'cor_trl', gx.cor_trl);
            offmax = 0;
            for a = 1:3
                for b = (a+1):4
                    offmax = max(offmax, max(abs(gmats.pi_e(:,a,b))));
                end
            end
            testCase.verifyTrue(offmax > 0, ...
                ['The gamma positive control must carry nonzero residual ', ...
                'off-diagonals, or the diagonal theorem test is vacuous.']);
        end

        function testRelativeErrorReducesToFourTermSum(testCase)
            %Corollary of the diagonal theorem, pinned end to end through the
            %shared composite kernel: with a diagonal residual channel the
            %harmonic covariance divisors have nothing to divide, and the
            %relative-error composite is exactly sum_q c_q^2*sigma_q^2/n_q,
            %written longhand here.
            fx = localOnefacetFixture(40);
            mats = psyrat_dod_components('facet', 'trial', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);
            cvec = [1 -1 -1 1];
            n = [18 22 15 20];
            comp = psyrat_dod_composite('varcov', mats.pi_e, 'cvec', cvec, ...
                'rule', 'harmonic', 'n', n, 'tol', 0);
            longhand = zeros(size(fx.alpha,1),1);
            for q = 1:4
                longhand = longhand + ...
                    cvec(q)^2 .* exp(2 .* fx.log_sigma(:,q)) ./ n(q);
            end
            testCase.verifyEqual(comp.value, longhand, 'AbsTol', 1e-12, ...
                ['The relative-error composite must reduce to the four-term ', ...
                'sum under the identity link.']);
        end

        function testSignalChannelsAreDirectCovarianceBlocks(testCase)
            %Under the identity link the observed-scale components are the
            %model covariance blocks DIRECTLY - no lognormal conversion. The
            %transcription below is written longhand, entry by entry.
            fx = localOnefacetFixture(60);
            mats = psyrat_dod_components('facet', 'trial', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);
            nd = size(fx.alpha,1);
            for q = 1:4
                testCase.verifyEqual(mats.p(:,q,q), fx.sd_p(:,q).^2, ...
                    'AbsTol', 1e-12, sprintf('p diagonal %d', q));
                testCase.verifyEqual(mats.i(:,q,q), fx.sd_trl(:,q).^2, ...
                    'AbsTol', 1e-12, sprintf('i diagonal %d', q));
                testCase.verifyEqual(mats.pi_e(:,q,q), ...
                    exp(2 .* fx.log_sigma(:,q)), 'AbsTol', 1e-12, ...
                    sprintf('residual diagonal %d', q));
            end
            for a = 1:3
                for b = (a+1):4
                    testCase.verifyEqual(mats.p(:,a,b), ...
                        fx.sd_p(:,a) .* fx.sd_p(:,b) .* fx.cor_p(:,a,b), ...
                        'AbsTol', 1e-12, sprintf('p cross (%d,%d)', a, b));
                    testCase.verifyEqual(mats.i(:,a,b), ...
                        fx.sd_trl(:,a) .* fx.sd_trl(:,b) .* fx.cor_trl(:,a,b), ...
                        'AbsTol', 1e-12, sprintf('i cross (%d,%d)', a, b));
                    testCase.verifyEqual(mats.p(:,a,b), mats.p(:,b,a), ...
                        'AbsTol', 0, 'p must be symmetric.');
                    testCase.verifyEqual(mats.i(:,a,b), mats.i(:,b,a), ...
                        'AbsTol', 0, 'i must be symmetric.');
                end
            end
            %contract fields: expected is alpha itself (identity link;
            %bit-identity because no arithmetic may touch it), v_total is the
            %summed random-effect variance excluding the residual (the same
            %structural meaning as the gamma assembler's log-scale v_total)
            testCase.verifyEqual(mats.expected, fx.alpha, 'AbsTol', 0, ...
                'expected must be alpha itself under the identity link.');
            testCase.verifyEqual(mats.v_total, fx.sd_p.^2 + fx.sd_trl.^2, ...
                'AbsTol', 1e-12, 'v_total must sum the random-effect variances.');
            testCase.verifyEqual(size(mats.pi_e), [nd 4 4]);
        end

        function testTwofacetBlocksMapDistinctly(testCase)
            %The two-facet assembler on all-distinct values, so a facet swap
            %(occ vs tid vs oid vs to) fails rather than passing by symmetry.
            fx = localTwofacetFixture(40);
            mats = psyrat_dod_components('facet', 'trial_occasion', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl, ...
                'sd_occ', fx.sd_occ, 'cor_occ', fx.cor_occ, ...
                'sd_tid', fx.sd_tid, 'cor_tid', fx.cor_tid, ...
                'sd_oid', fx.sd_oid, 'cor_oid', fx.cor_oid, ...
                'sd_to', fx.sd_to, 'cor_to', fx.cor_to);
            pairs = { ...
                'p',  fx.sd_p,   fx.cor_p; ...
                'i',  fx.sd_trl, fx.cor_trl; ...
                'o',  fx.sd_occ, fx.cor_occ; ...
                'pi', fx.sd_tid, fx.cor_tid; ...
                'po', fx.sd_oid, fx.cor_oid; ...
                'io', fx.sd_to,  fx.cor_to};
            for k = 1:size(pairs,1)
                key = pairs{k,1}; sd = pairs{k,2}; cor = pairs{k,3};
                for q = 1:4
                    testCase.verifyEqual(mats.(key)(:,q,q), sd(:,q).^2, ...
                        'AbsTol', 1e-12, sprintf('%s diagonal %d', key, q));
                end
                for a = 1:3
                    for b = (a+1):4
                        testCase.verifyEqual(mats.(key)(:,a,b), ...
                            sd(:,a) .* sd(:,b) .* cor(:,a,b), ...
                            'AbsTol', 1e-12, ...
                            sprintf('%s cross (%d,%d)', key, a, b));
                    end
                end
            end
            %two-facet residual: pio_e diagonal = sigma^2, off-diagonals
            %exactly zero (same theorem, same bit-identity expectation)
            for q = 1:4
                testCase.verifyEqual(mats.pio_e(:,q,q), ...
                    exp(2 .* fx.log_sigma(:,q)), 'AbsTol', 1e-12, ...
                    sprintf('pio_e diagonal %d', q));
            end
            for a = 1:3
                for b = (a+1):4
                    testCase.verifyEqual(mats.pio_e(:,a,b), ...
                        zeros(size(fx.alpha,1),1), 'AbsTol', 0, sprintf( ...
                        'pio_e cross (%d,%d) must be exactly zero.', a, b));
                end
            end
            testCase.verifyEqual(mats.v_total, ...
                fx.sd_p.^2 + fx.sd_trl.^2 + fx.sd_occ.^2 + ...
                fx.sd_tid.^2 + fx.sd_oid.^2 + fx.sd_to.^2, ...
                'AbsTol', 1e-12, 'two-facet v_total must sum all six blocks.');
        end

        function testMeanAndZInvariance(testCase)
            %The second identity-link theorem, pinned by perturbation. Moving
            %alpha (which is where b, b_dim, and the evaluation point z enter
            %this assembler) must move ONLY mats.expected; moving log_sigma
            %must move ONLY the residual diagonal. AbsTol 0: the unchanged
            %matrices must be BIT-identical, because alpha and log_sigma must
            %not appear in their assembly at all.
            fx = localOnefacetFixture(30);
            base = psyrat_dod_components('facet', 'trial', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);
            shifted = psyrat_dod_components('facet', 'trial', ...
                'alpha', fx.alpha + 7.3, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);
            for key = {'p', 'i', 'pi_e', 'v_total'}
                testCase.verifyEqual(shifted.(key{1}), base.(key{1}), ...
                    'AbsTol', 0, sprintf( ...
                    'A location shift must not move %s.', key{1}));
            end
            testCase.verifyEqual(shifted.expected, fx.alpha + 7.3, ...
                'AbsTol', 0, 'The location shift must reach expected.');

            scaled = psyrat_dod_components('facet', 'trial', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma + 0.4, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);
            for key = {'p', 'i', 'expected', 'v_total'}
                testCase.verifyEqual(scaled.(key{1}), base.(key{1}), ...
                    'AbsTol', 0, sprintf( ...
                    'A scale shift must not move %s.', key{1}));
            end
            for q = 1:4
                testCase.verifyEqual(scaled.pi_e(:,q,q), ...
                    exp(2 .* (fx.log_sigma(:,q) + 0.4)), 'AbsTol', 1e-12, ...
                    'The scale shift must reach the residual diagonal.');
            end
        end

        function testCouplingVocabularyRejected(testCase)
            %The nonconcurrent boundary at the signature, identical in kind to
            %the gamma assembler: coupling vocabulary is rejected by name,
            %never silently swallowed by the name-value parser.
            for coupled = {'rho', 'rho_e', 'cov_e_obs', 'cov_res_e_obs', ...
                    'er_cov', 'rescor', 'copula'}
                testCase.verifyError(@() psyrat_dod_components( ...
                    coupled{1}, 0.3, 'facet', 'trial'), ...
                    'psyrat_dod_components:noresidualcoupling');
            end
        end

        function testCorrelationDiagonalGuard(testCase)
            %Passing a covariance block where a correlation block belongs is
            %the likeliest scale error on the assembler signature.
            fx = localOnefacetFixture(5);
            bad = fx.cor_p .* 1.7;   %no longer unit-diagonal
            testCase.verifyError(@() psyrat_dod_components( ...
                'facet', 'trial', 'alpha', fx.alpha, ...
                'log_sigma', fx.log_sigma, 'sd_p', fx.sd_p, 'cor_p', bad, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl), ...
                'psyrat_dod_components:notcorrelation');
        end

    end
end

% ---------------------------------------------------------------------------
% Local fixtures: seeded synthetic pseudo-posterior draws on GAUSSIAN
% (identity-link) scales - negative cell means included deliberately, because
% they are representable under this family and not under gamma - with
% all-distinct per-cell values and at least one negative correlation per
% block, so index, sign, and facet swaps fail rather than passing by symmetry.
% Seeds differ from TestGammaDodDynrelLsConversion's so the two classes cannot
% accidentally share draws.
% ---------------------------------------------------------------------------

function fx = localOnefacetFixture(nd)
rng(5252, 'twister');
fx.alpha = [-4.8 -1.9 -4.1 -3.4] + 0.3 * randn(nd, 4);
fx.log_sigma = [0.61 1.01 0.64 0.66] + 0.05 * randn(nd, 4);
fx.sd_p = abs([2.0 1.7 2.2 1.9] + 0.15 * randn(nd, 4));
fx.sd_trl = abs([0.80 0.70 0.90 0.75] + 0.06 * randn(nd, 4));
fx.cor_p = localCorrDraws(nd, ...
    [0.55 0.30; 0.48 -0.35; -0.42 0.28; 0.36 0.22], [0.60 0.68 0.75 0.83]);
fx.cor_trl = localCorrDraws(nd, ...
    [0.50 0.25; -0.38 0.30; 0.27 -0.33; 0.36 0.22], [0.62 0.70 0.78 0.85]);
end

function fx = localTwofacetFixture(nd)
fx = localOnefacetFixture(nd);
rng(41414, 'twister');
fx.sd_occ = abs([0.45 0.38 0.52 0.41] + 0.04 * randn(nd, 4));
fx.sd_tid = abs([0.40 0.47 0.35 0.44] + 0.04 * randn(nd, 4));
fx.sd_oid = abs([0.34 0.28 0.39 0.31] + 0.04 * randn(nd, 4));
fx.sd_to = abs([0.30 0.38 0.26 0.34] + 0.04 * randn(nd, 4));
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
