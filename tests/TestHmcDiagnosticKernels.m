classdef TestHmcDiagnosticKernels < matlab.unittest.TestCase
    %TESTHMCDIAGNOSTICKERNELS The native HMC engine's R-hat and effective sample
    %size must FAIL on frozen chains and on chains that disagree (work queue Q5;
    %post-beta review PBA-03).
    %
    % Before Q5 the kernels in psyrat_hmc_sample returned R-hat 1 when the
    % within-chain variance was zero and ESS = N when the autocovariance was
    % zero, so four chains stuck at four different constants reported R-hat 1
    % and ESS 400 and passed psyrat_checkconv; the ESS also centered each
    % chain on its own mean, so chains sitting at different values reported a
    % full ESS. The acceptance ratio drawSamples returns was discarded.
    %
    % The kernels are local functions, so they are exercised through the
    % documented diagnostics seam of psyrat_hmc_sample ([~,conv] = ...), which builds the
    % convergence cell from a fixed draw array through the same code the
    % sampling path uses. The sampling path itself keeps its positive control
    % in TestNativeEngine.testHmcHarnessRecoversAndConverges.

    methods (Test)

        function testChainsFrozenAtDifferentConstantsFail(testCase)
            x = repmat([0 1 2 3], 100, 1);      % 100 draws x 4 chains, each constant
            [~, conv] = psyrat_hmc_sample('diagnostics', x, {'sig_u'}, [0.7 0.7 0.7 0.7]);
            testCase.verifyEqual(conv{2, 1}, 'sig_u');
            testCase.verifyEqual(conv{2, 3}, Inf, 'R-hat must be Inf for frozen, disagreeing chains');
            testCase.verifyLessThan(conv{2, 2}, 40, 'ESS must fail the 10*nchains floor');
        end

        function testChainsFrozenAtTheSameConstantFail(testCase)
            x = ones(100, 4);
            [~, conv] = psyrat_hmc_sample('diagnostics', x, {'sig_u'}, [0.7 0.7 0.7 0.7]);
            testCase.verifyEqual(conv{2, 3}, Inf);
            testCase.verifyEqual(conv{2, 2}, 0);
        end

        function testSeparatedMeansGiveTinyEss(testCase)
            rng(123);
            x = randn(2000, 4) + [0 10 20 30];
            [~, conv] = psyrat_hmc_sample('diagnostics', x, {'mu'}, [0.7 0.7 0.7 0.7]);
            testCase.verifyGreaterThan(conv{2, 3}, 1.1, 'R-hat must flag separated chains');
            testCase.verifyLessThan(conv{2, 2}, 40, ...
                'ESS must account for between-chain disagreement (it reported 8000 before Q5)');
        end

        function testIidChainsPass(testCase)
            rng(7);
            x = randn(500, 4, 2);               % two monitored parameters
            [~, conv] = psyrat_hmc_sample('diagnostics', x, {'a', 'b'}, [0.65 0.7 0.6 0.66]);
            testCase.verifySize(conv, [3, 10]);  % header + 2 rows; 3 cols + cell 4 + label + 4 ratios + metadata
            for r = 2:3
                testCase.verifyLessThan(conv{r, 3}, 1.02);
                testCase.verifyGreaterThan(conv{r, 2}, 0.6 * 2000);
                testCase.verifyLessThan(conv{r, 2}, 1.5 * 2000);
            end
        end

        function testAutocorrelatedChainsHaveReducedEss(testCase)
            % an AR(1) chain with phi = 0.9 has ESS about N*(1-phi)/(1+phi) = N/19
            rng(11);
            n = 4000; m = 4; phi = 0.9;
            x = zeros(n, m);
            for c = 1:m
                e = randn(n, 1);
                for t = 2:n
                    x(t, c) = phi * x(t-1, c) + e(t);
                end
            end
            [~, conv] = psyrat_hmc_sample('diagnostics', x, {'a'}, [0.7 0.7 0.7 0.7]);
            expected = n * m * (1 - phi) / (1 + phi);
            testCase.verifyGreaterThan(conv{2, 2}, 0.5 * expected);
            testCase.verifyLessThan(conv{2, 2}, 1.6 * expected);
            testCase.verifyLessThan(conv{2, 3}, 1.05);
        end

        function testAcceptanceFloorFailsEveryRow(testCase)
            rng(3);
            x = randn(500, 4, 2);
            ratios = [0.7 0.05 0.7 0.7];
            [~, conv] = testCase.verifyWarning( ...
                @() psyrat_hmc_sample('diagnostics', x, {'a', 'b'}, ratios), ...
                'PsyRAT:hmcAcceptanceFloor');
            testCase.verifyEqual(conv{2, 3}, Inf);
            testCase.verifyEqual(conv{3, 3}, Inf);
            testCase.verifyEqual(conv{1, 5}, 'accept_ratio');
            testCase.verifyEqual([conv{1, 6:9}], ratios, 'AbsTol', 1e-12);
        end

        function testHealthyAcceptanceIsRecordedWithoutWarning(testCase)
            rng(5);
            x = randn(300, 4);
            ratios = [0.65 0.71 0.62 0.68];
            [~, conv] = testCase.verifyWarningFree( ...
                @() psyrat_hmc_sample('diagnostics', x, {'a'}, ratios));
            testCase.verifyEqual(conv{1, 5}, 'accept_ratio');
            testCase.verifyEqual([conv{1, 6:9}], ratios, 'AbsTol', 1e-12);
            testCase.verifyLessThan(conv{2, 3}, 1.05);
        end

    end
end
