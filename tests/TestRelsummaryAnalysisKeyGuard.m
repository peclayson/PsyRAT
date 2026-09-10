classdef TestRelsummaryAnalysisKeyGuard < PsyRATTestBase
    %Audit 2026-09-05: psyrat_relsummary's 'analysis' argument is a dispatch
    %key ('sing', 'trt', ...), not the stored rel.analysis string ('ic',
    %'trt_diff', ...), and neither dispatch switch has an otherwise arm. An
    %unknown key, including the stored value the manual calls the most common
    %scripted mistake, fell through both switches: relerr was never assigned,
    %so the documented two-output call failed with MATLAB's generic
    %output-argument error. The key is now validated with varargin:analysis.
    %
    %RED against the pre-fix code: both rejection tests (the pre-fix failure
    %is a different, unnamed error). The accepted-key test is the positive
    %control and passes on both builds.

    methods (Test)

        function testStoredRelAnalysisValueIsRejectedWithTheMapping(testCase)
            pd = PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();
            pd.proc = struct('interactive', false);
            testCase.verifyError(@() psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'ic', 'depcutoff', 0.80, 'meascutoff', 2, 'depcentmeas', 1), ...
                'varargin:analysis');
            %the message check cannot pass vacuously: a call that does not
            %throw fails the assertion below
            threw = false;
            try
                psyrat_relsummary('psyrat_data', pd, 'analysis', 'ic', ...
                    'depcutoff', 0.80, 'meascutoff', 2, 'depcentmeas', 1);
            catch err
                threw = true;
                testCase.verifyEqual(err.identifier, 'varargin:analysis');
                testCase.verifyTrue(contains(err.message, '''sing'''), ...
                    'the message must name the dispatch key to use');
            end
            testCase.assertTrue(threw, 'the stored rel.analysis value must be rejected');
        end

        function testNonTextKeyIsRejected(testCase)
            pd = PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();
            pd.proc = struct('interactive', false);
            testCase.verifyError(@() psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 1, 'depcutoff', 0.80, 'meascutoff', 2, 'depcentmeas', 1), ...
                'varargin:analysis');
        end

        function testAcceptedKeyStillRuns(testCase)
            pd = PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();
            pd.proc = struct('interactive', false);
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, 'analysis', 'sing', ...
                'depcutoff', 0.30, 'meascutoff', 2, 'depcentmeas', 1, 'CI', 0.95);
            testCase.verifyTrue(isfield(out, 'relsummary'));
            testCase.verifyTrue(isfield(relerr, 'nogooddata'));
        end

    end
end
