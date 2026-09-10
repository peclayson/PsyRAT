classdef TestCmdStanIntegration < PsyRATTestBase
    % Integration tests that execute real CmdStan models through psyrat_computevarcomp.

    properties (TestParameter)
        Scenario = PsyRATCmdStanIntegrationFactory.scenarios();
    end

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            outDir = localCmdStanOutputDir(testCase.projectRoot());
            pats = {'cmdstan*.stan', 'cmdstan*.hpp', 'cmdstan18-*', 'output-*.csv', 'temp.data.R'};
            for i = 1:numel(pats)
                files = dir(fullfile(outDir, pats{i}));
                for j = 1:numel(files)
                    if ~files(j).isdir
                        delete(fullfile(outDir, files(j).name));
                    end
                end
            end

            dirs = dir(fullfile(outDir, 'cmdstan_*'));
            dirs = dirs([dirs.isdir]);
            dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
            for i = 1:numel(dirs)
                try
                    rmdir(fullfile(outDir, dirs(i).name), 's');
                catch
                    warning('tests:cleanupDirFailed', ...
                        ['Could not remove CmdStan artifact directory ''%s''. ', ...
                        'How to fix: close any process that may still be using it and delete it manually.'], ...
                        fullfile(outDir, dirs(i).name));
                end
            end
        end
    end

    methods (Test)
        function testCmdStanModelsAcrossInputCombinations(testCase, Scenario)
            psyrat_assume_dev_testdata(testCase);
            dataTbl = PsyRATCmdStanIntegrationFactory.makeScenarioTable(testCase.projectRoot(), Scenario);
            cmdstanOutDir = localCmdStanOutputDir(testCase.projectRoot());

            runScenario = @() psyrat_computevarcomp( ...
                'data', dataTbl, ...
                'chains', 1, ...
                'iter', 20, ...
                'verbose', 1, ...
                'showgui', 1, ...
                'cmdstanoutdir', cmdstanOutDir, ...
                'sserrvar', Scenario.sserrvar, ...
                'diffest', Scenario.diffest, ...
                'diffwpcov', Scenario.diffwpcov, ...
                'diffrescor', Scenario.diffrescor);

            if isfield(Scenario, 'expectedErrorId') && ~isempty(Scenario.expectedErrorId)
                testCase.verifyError(runScenario, Scenario.expectedErrorId);
                return;
            end

            rel = runScenario();

            testCase.verifyTrue(isfield(rel, 'analysis'));
            testCase.verifyEqual(rel.analysis, Scenario.expectedAnalysis);
            testCase.verifyTrue(isfield(rel, 'out'));

            switch Scenario.expectedAnalysis
                case 'ic'
                    required = {'mu', 'sig_u', 'sig_e', 'labels', 'conv'};
                case 'trt'
                    required = {'mu', 'sig_id', 'sig_occ', 'sig_trl', ...
                        'sig_trlxid', 'sig_occxid', 'sig_trlxocc', 'sig_err', ...
                        'labels', 'conv'};
                case 'ic_sserrvar'
                    required = {'mu', 'ind_bs', 'gro_sds', 'pop_sdlog', ...
                        'ind_sdlog', 'sig_trl', 'labels', 'id_matches', 'conv'};
                case 'ic_diff'
                    required = {'b', 'b_sigma', 'sd_id', 'sd_trl', ...
                        'id_varcov', 'trl_varcov', 'rescor', 'err_varcov', ...
                        'elabels', 'glabels', 'conv'};
                case 'ic_diff_sserrvar'
                    required = {'b', 'b_sigma', 'sd_id', 'sd_trl', ...
                        'id_varcov', 'trl_varcov', 'rescor', 'err_varcov', ...
                        'er_var_ss1', 'er_var_ss2', 'wp_cov_ss', ...
                        'id_matches', 'elabels', 'glabels', 'conv'};
                otherwise
                    required = {};
            end

            for i = 1:numel(required)
                testCase.verifyTrue(isfield(rel.out, required{i}));
                testCase.verifyFalse(isempty(rel.out.(required{i})));
            end

            if isfield(rel.out, 'labels')
                expectedLabels = localExpectedLabelCount(dataTbl);
                if ischar(rel.out.labels)
                    actualLabels = 1;
                else
                    actualLabels = numel(rel.out.labels);
                end
                testCase.verifyEqual(actualLabels, expectedLabels);
            end
        end
    end
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function n = localExpectedLabelCount(tbl)
hasGroup = any(strcmp(tbl.Properties.VariableNames, 'group'));
hasEvent = any(strcmp(tbl.Properties.VariableNames, 'event'));
hasTime = any(strcmp(tbl.Properties.VariableNames, 'time'));

if hasTime
    nGroup = 1;
    nEvent = 1;
    if hasGroup
        nGroup = numel(unique(tbl.group));
    end
    if hasEvent
        nEvent = numel(unique(tbl.event));
    end
    n = nGroup * nEvent;
else
    if hasGroup && hasEvent
        n = numel(unique(tbl.group)) * numel(unique(tbl.event));
    elseif hasGroup
        n = numel(unique(tbl.group));
    elseif hasEvent
        n = numel(unique(tbl.event));
    else
        n = 1;
    end
end
end
