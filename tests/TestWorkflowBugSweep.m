classdef TestWorkflowBugSweep < PsyRATTestBase

    methods (Test)
        function testWorkflowRunMatrixIncludesDiffDecisionModes(testCase)
            runs = psyrat_workflow_runs('ic_diff');
            testCase.verifyEqual(numel(runs), 1);
            testCase.verifyEqual(runs{1}.analysis, 'sing_diff');
            testCase.verifyEqual(runs{1}.diffgcoeff, [1 2]);
            
            %One run, not two. psyrat_startview_sing maps rel.analysis to a
            %builder label 1:1, so ic_diff_sserrvar reaches ONLY the
            %subject-level builder on every product surface. The matrix used
            %to sweep it through 'sing_diff' as well, which built a relsummary
            %with no diffscore.ssrel_table while the output functions still
            %dispatched on rel.analysis and took the subject-level branch --
            %48 sweep failures in a configuration no user can create (RC-38).
            %Assert the analysis label, not just the count: the old test
            %checked numel and diffgcoeff only, so it mirrored the defect
            %instead of catching it.
            runs2 = psyrat_workflow_runs('ic_diff_sserrvar');
            testCase.verifyEqual(numel(runs2), 1);
            testCase.verifyEqual(runs2{1}.analysis, 'sing_diff_sserr');
            testCase.verifyEqual(runs2{1}.diffgcoeff, [1 2]);
            testCase.verifyFalse( ...
                any(cellfun(@(r) strcmp(r.analysis,'sing_diff'), runs2)), ...
                'ic_diff_sserrvar must not be swept through the sing_diff builder.');
        end
        
        function testWorkflowRunMatrixCoversKnownFamilies(testCase)
            runsIC = psyrat_workflow_runs('ic');
            runsTRT = psyrat_workflow_runs('trt');
            runsUNK = psyrat_workflow_runs('bad_analysis_name');
            
            testCase.verifyEqual(runsIC{1}.analysis, 'sing');
            testCase.verifyEqual(runsTRT{1}.analysis, 'trt');
            testCase.verifyEqual(runsUNK{1}.analysis, 'unknown');
        end
        
        function testRunExhaustiveWorkflowBugSweepNoScenariosValidation(testCase)
            testCase.verifyError(@() runExhaustiveWorkflowBugSweep(...
                'ScenarioFilter', {'__no_match__'}, ...
                'Verbose', false, ...
                'SaveReport', false), ...
                'workflow:noscenarios');
        end

        function testRunExhaustiveWorkflowBugSweepAcceptsDiffWpCovOption(testCase)
            testCase.verifyError(@() runExhaustiveWorkflowBugSweep(...
                'ScenarioFilter', {'__no_match__'}, ...
                'Verbose', false, ...
                'SaveReport', false, ...
                'DiffWpCov', 1), ...
                'workflow:noscenarios');
        end
    end
end
