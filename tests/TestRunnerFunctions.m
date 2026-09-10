classdef TestRunnerFunctions < PsyRATTestBase

    methods (Test)
        function testRunAccuracyValidationSuiteSelectsOracleTests(testCase)
            %The accuracy lane deliberately spans THREE classes with different
            %evidentiary weight, so this asserts membership rather than a single class:
            %  TestCalculationAccuracyOracle    regression only (the oracle mirrors the
            %                                   production formulas)
            %  TestIndependentOracleAgreement   external corroboration (paper-sourced)
            %  TestExternalBenchmark            external numeric anchor (published tables)
            %The last two were previously absent, which meant the lane named "accuracy"
            %contained nothing that could establish accuracy. Every selected test must
            %still come from one of the three -- an unrelated class leaking in would mean
            %the substring selection had gone wrong.
            expected = {'TestCalculationAccuracyOracle', ...
                        'TestIndependentOracleAgreement', ...
                        'TestExternalBenchmark'};

            [allPassed, results] = runAccuracyValidationSuite('Verbose', false);

            testCase.verifyGreaterThan(numel(results), 0);
            inLane = false(1, numel(results));
            for k = 1:numel(expected)
                inLane = inLane | contains({results.Name}, expected{k});
            end
            testCase.verifyTrue(all(inLane), ...
                'runAccuracyValidationSuite selected a class outside the accuracy lane.');

            %Each named class must actually be present, so silently dropping one (e.g. by
            %renaming a file) fails loudly instead of shrinking the lane unnoticed.
            for k = 1:numel(expected)
                testCase.verifyTrue(any(contains({results.Name}, expected{k})), ...
                    sprintf('accuracy lane is missing %s', expected{k}));
            end

            testCase.verifyEqual(allPassed, all([results.Passed]));
        end

        function testRunAllUnitTestsSelectionForPreflightClass(testCase)
            [allPassed, results] = runAllUnitTests('Verbose', false, 'Select', 'TestPreflightValidation');

            testCase.verifyGreaterThan(numel(results), 0);
            testCase.verifyTrue(all(contains({results.Name}, 'TestPreflightValidation')));
            testCase.verifyEqual(allPassed, all([results.Passed]));
        end

        function testRunAccuracyValidationSuiteCanWriteJUnit(testCase)
            junitFile = [tempname '.xml'];
            c = onCleanup(@() localDeleteIfExists(junitFile)); %#ok<NASGU>

            [~, ~] = runAccuracyValidationSuite( ...
                'Verbose', false, ...
                'ReportJUnit', true, ...
                'JUnitFile', junitFile);

            testCase.verifyTrue(isfile(junitFile));
        end
    end
end

function localDeleteIfExists(path)
if isfile(path)
    delete(path);
end
end
