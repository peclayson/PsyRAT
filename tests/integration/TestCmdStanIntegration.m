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

            % Q4 (2026-09-11; post-beta review PBA-01): the CmdStan diagnostics
            % come from the stansummary CSV by column name. Every convergence
            % table carries bulk ESS, tail ESS and MCSE from column 4 of its
            % data rows, labels them from header cell 5 (cell 4 stays the
            % numeric divergence count), and gates on the smaller ESS.
            conv = rel.out.conv.data;
            if iscell(conv) && ~isempty(conv) && iscell(conv{1})
                conv = conv{1};             % nested layout: the first fit
            end
            testCase.assertTrue(iscell(conv) && size(conv,1) >= 2, ...
                'the convergence table must carry at least one parameter row');
            testCase.assertFalse(strcmp(conv{2,1}, 'convergence_unavailable'), ...
                ['diagnostics were not extracted for this fit (print(fit) and the ' ...
                'direct stansummary call both returned nothing)']);
            testCase.verifyGreaterThanOrEqual(size(conv,2), 7, ...
                'the table must carry the ess_bulk/ess_tail/mcse columns (CSV route)');
            testCase.verifyEqual(conv{1,5}, 'ess_bulk');
            metadata=psyrat_diagnostic_metadata(conv);
            testCase.verifyEqual(metadata.format,'modern_csv');
            testCase.verifyEqual(metadata.precision_significant_digits,8);
            testCase.verifyFalse(metadata.degraded);
            testCase.verifyEqual(metadata.completeness_source,'chain_csv_headers');
            testCase.verifyTrue(all(ismember(metadata.expected_names,conv(2:end,1))));
            testCase.verifyNotEmpty(metadata.model_source_hash);

            testCase.verifyTrue(isnumeric(conv{1,4}) && isscalar(conv{1,4}), ...
                'the divergence count must stay in header cell 4');
            for rr = 2:size(conv,1)
                b = conv{rr,4}; t = conv{rr,5}; n = conv{rr,2};
                if isnan(b) || isnan(t)
                    testCase.verifyTrue(isnan(n), 'a constant row keeps nan as n_eff');
                else
                    testCase.verifyEqual(n, min(b,t), 'AbsTol', 1e-12, ...
                        'n_eff must be the smaller of bulk and tail ESS');
                end
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
