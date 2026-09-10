classdef TestCorrCholLkj < matlab.unittest.TestCase
    % Unit tests for psyrat_corr_chol_lkj, the general K x K LKJ correlation-matrix
    % Cholesky transform used by the native HMC engine for the > 2-D correlation
    % blocks (e.g. the case-13 joint 4-D person block). These are fast,
    % deterministic checks (no sampling): the implied constrained-correlation
    % marginals are validated against the known LKJ(eta) marginals by a separate
    % CmdStan-free prior-marginal Monte Carlo check during development, and by the
    % case-13 recovery integration test downstream.

    methods (Test)
        function testK1Trivial(testCase)
            % K = 1 is a degenerate 1x1 correlation "matrix": L = 1, lp = 0.
            [L, lp] = psyrat_corr_chol_lkj(zeros(0,1), 1, 2);
            testCase.verifyEqual(L, 1);
            testCase.verifyEqual(lp, 0);
        end

        function testK2ReducesToValidatedTwoByTwoForm(testCase)
            % For K = 2, eta = 2 the log-density must equal the validated 2x2 form
            % used by the diff family: 2*log(1 - rho^2) with rho = tanh(y).
            for y = [-1.3, -0.2, 0.0, 0.6, 1.1]
                [L, lp] = psyrat_corr_chol_lkj(y, 2, 2);
                rho = tanh(y);
                testCase.verifyEqual(lp, 2*log(1 - rho^2), 'AbsTol', 1e-12);
                % L is the 2x2 correlation Cholesky [1 0; rho sqrt(1-rho^2)]
                testCase.verifyEqual(L(1,1), 1, 'AbsTol', 1e-12);
                testCase.verifyEqual(L(2,1), rho, 'AbsTol', 1e-12);
                testCase.verifyEqual(L(2,2), sqrt(1-rho^2), 'AbsTol', 1e-12);
            end
        end

        function testK2EtaOneReducesToUniform(testCase)
            % eta = 1 (uniform over correlation matrices): the 2x2 unconstrained
            % density is log(1 - rho^2) (LKJ(1) density is flat; only the tanh
            % Jacobian remains).
            y = 0.4; rho = tanh(y);
            [~, lp] = psyrat_corr_chol_lkj(y, 2, 1);
            testCase.verifyEqual(lp, log(1 - rho^2), 'AbsTol', 1e-12);
        end

        function testProducesValidCorrelationCholesky(testCase)
            % For K = 2..6 and many random unconstrained vectors, L must be lower-
            % triangular with R = L*L' a valid correlation matrix (unit diagonal,
            % positive semidefinite).
            rng(20260627, 'twister');
            for K = [2 3 4 5 6 8]   % K = 8 is the DoD (case 9) person/trial block size
                D = K*(K-1)/2;
                for rep = 1:50
                    y = randn(D,1) * 1.0;
                    L = psyrat_corr_chol_lkj(y, K, 2);
                    testCase.verifySize(L, [K K]);
                    % lower-triangular
                    testCase.verifyEqual(max(max(abs(triu(L,1)))), 0, 'AbsTol', 1e-12);
                    R = L*L.';
                    % unit diagonal
                    testCase.verifyEqual(diag(R), ones(K,1), 'AbsTol', 1e-9);
                    % positive semidefinite
                    testCase.verifyGreaterThan(min(eig((R+R.')/2)), -1e-9);
                    % valid correlations
                    testCase.verifyLessThanOrEqual(max(abs(R(:))), 1 + 1e-9);
                end
            end
        end

        function testADGradientMatchesFiniteDifference(testCase)
            % The AD (dlarray) gradient of the log-density must match central
            % finite differences. Gated on the Deep Learning Toolbox (dlarray);
            % skips cleanly otherwise.
            % NB: dlgradient is callable only inside dlfeval, so
            % exist('dlgradient','file') is always 0; probe dlfeval instead.
            testCase.assumeTrue((exist('dlarray','file') > 0) && ...
                (exist('dlfeval','file') > 0), ...
                'Deep Learning Toolbox (dlarray AD) not available.');
            rng(7, 'twister');
            for K = [2 3 4 5 8]   % K = 8 is the DoD (case 9) person/trial block size
                D = K*(K-1)/2;
                y0 = randn(D,1) * 0.6;
                [~, gad] = dlfeval(@localGrad, dlarray(y0), K);
                gad = double(extractdata(gad));
                h = 1e-5; gfd = zeros(D,1);
                for k = 1:D
                    e = zeros(D,1); e(k) = h;
                    [~, lpp] = psyrat_corr_chol_lkj(y0+e, K, 2);
                    [~, lpm] = psyrat_corr_chol_lkj(y0-e, K, 2);
                    gfd(k) = (lpp - lpm) / (2*h);
                end
                relerr = norm(gad - gfd) / max(1, norm(gfd));
                testCase.verifyLessThan(relerr, 1e-6, ...
                    sprintf('AD vs FD gradient mismatch at K=%d (relerr %.2e).', K, relerr));
            end
        end
    end
end

function [lp, g] = localGrad(yd, K)
% Traced inside dlfeval: log-density and its reverse-mode gradient. yd must be
% the dlarray passed THROUGH dlfeval so the operations are traced.
[~, lp] = psyrat_corr_chol_lkj(yd, K, 2);
g = dlgradient(lp, yd);
end
