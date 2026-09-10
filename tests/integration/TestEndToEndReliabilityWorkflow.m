classdef TestEndToEndReliabilityWorkflow < PsyRATTestBase
    % End-to-end checks from CmdStan variance estimation to reliability summaries.

    properties (Constant)
        Chains = 1;
        Iter = 30;
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
        function testSingleSessionInternalConsistencyEndToEnd(testCase)
            scenario = PsyRATCmdStanIntegrationFactory.scenarios().ic_group_event;
            psyrat_data = localRunScenario(testCase, scenario);

            [psyrat_data, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing', ...
                'depcutoff', 0.5, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'CI', 0.95);

            testCase.verifyEqual(relerr.nogooddata, 0);
            ev = psyrat_data.relsummary.group(1).event(1);
            testCase.verifyGreaterThanOrEqual(ev.dep.m, 0);
            testCase.verifyLessThanOrEqual(ev.dep.m, 1);
            testCase.verifyGreaterThanOrEqual(ev.icc.m, 0);
            testCase.verifyLessThanOrEqual(ev.icc.m, 1);

            depCut = psyrat_depcutofft('psyrat_data', psyrat_data, 'gui', 0);
            depOverall = psyrat_depoverallt('psyrat_data', psyrat_data, 'gui', 0);
            varTbl = psyrat_variancet('psyrat_data', psyrat_data, 'gui', 0);
            testCase.verifyGreaterThanOrEqual(height(depCut), 1);
            testCase.verifyGreaterThanOrEqual(height(depOverall), 1);
            testCase.verifyGreaterThanOrEqual(height(varTbl), 1);
        end

        function testSubjectLevelInternalConsistencyEndToEnd(testCase)
            scenario = PsyRATCmdStanIntegrationFactory.scenarios().ic_sserr_group_event;
            psyrat_data = localRunScenario(testCase, scenario);

            [psyrat_data, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing_sserr', ...
                'depcutoff', 0.5, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'CI', 0.95);

            testCase.verifyEqual(relerr.nogooddata, 0);
            ss = psyrat_data.relsummary.group(1).event(1).ssrel_table;
            testCase.verifyTrue(any(strcmp(ss.Properties.VariableNames, 'dep_pt')));
            testCase.verifyTrue(any(strcmp(ss.Properties.VariableNames, 'icc_pt')));
            testCase.verifyTrue(any(strcmp(ss.Properties.VariableNames, 'sem_pt')));
            testCase.verifyTrue(any(strcmp(ss.Properties.VariableNames, 'ind2include')));
            testCase.verifyTrue(any(strcmp(ss.Properties.VariableNames, 'ind2exclude')));

            % S12 regression: this scenario deliberately reuses subject ids
            % across groups (PsyRATCmdStanIntegrationFactory preserves the
            % overlap for includeGroup=true scenarios). Independently derive
            % the ground-truth per-(group,event) trial count directly from
            % the same raw scenario table psyrat_computevarcomp consumed, and
            % confirm relsummary's ssrel_table.trls matches it exactly for
            % EVERY group/event cell -- not the pooled-across-groups value a
            % regression of the S12 bug would produce.
            dataTbl = PsyRATCmdStanIntegrationFactory.makeScenarioTable( ...
                testCase.projectRoot(), scenario);
            groupNames = psyrat_data.rel.groups(:);
            eventNames = psyrat_data.rel.events(:);
            for g = 1:numel(groupNames)
                for e = 1:numel(eventNames)
                    ss = psyrat_data.relsummary.group(g).event(e).ssrel_table;
                    cell = dataTbl(strcmp(string(dataTbl.group), string(groupNames{g})) & ...
                        strcmp(string(dataTbl.event), string(eventNames{e})), :);
                    trueCounts = varfun(@length, cell, 'GroupingVariables', {'id'}, ...
                        'InputVariables', 'meas');
                    trueById = containers.Map( ...
                        cellstr(string(trueCounts.id)), num2cell(trueCounts.GroupCount));
                    obsById = containers.Map( ...
                        cellstr(string(ss.id)), num2cell(ss.trls));
                    ks = keys(obsById);
                    for k = 1:numel(ks)
                        testCase.verifyEqual(obsById(ks{k}), trueById(ks{k}), ...
                            sprintf('group %s, event %s, id %s', ...
                            groupNames{g}, eventNames{e}, ks{k}));
                    end
                end
            end
        end

        function testDifferenceScoreReliabilityEndToEnd(testCase)
            scenario = PsyRATCmdStanIntegrationFactory.scenarios().ic_diff_group_event;
            psyrat_data = localRunScenario(testCase, scenario);

            [psyrat_data, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing_diff', ...
                'depcutoff', 0.5, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'CI', 0.95);

            testCase.verifyEqual(relerr.nogooddata, 0);
            grp = psyrat_data.relsummary.group(1);
            testCase.verifyEqual(numel(grp.event), 2);
            testCase.verifyTrue(isfield(grp, 'diffscore'));
            testCase.verifyTrue(isfield(grp.diffscore, 'pt'));
            testCase.verifyGreaterThanOrEqual(grp.diffscore.pt, 0);
            testCase.verifyLessThanOrEqual(grp.diffscore.pt, 1);
        end

        function testSubjectLevelTrtAndDiffTagsAndTables(testCase)
            trtScenario = PsyRATCmdStanIntegrationFactory.scenarios().trt_group_event_time_sserr;
            trtData = localRunScenario(testCase, trtScenario);
            testCase.verifyEqual(trtData.rel.analysis, 'trt_sserrvar');

            diffScenario = PsyRATCmdStanIntegrationFactory.scenarios().ic_diff_group_event_sserr;
            diffData = localRunScenario(testCase, diffScenario);
            testCase.verifyEqual(diffData.rel.analysis, 'ic_diff_sserrvar');
        end

        function testTestRetestReliabilityCombinationsEndToEnd(testCase)
            scenario = PsyRATCmdStanIntegrationFactory.scenarios().trt_group_event_time;
            baseData = localRunScenario(testCase, scenario);

            combos = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3];
            validCombos = 0;
            tableData = [];
            for i = 1:size(combos, 1)
                gcoeff = combos(i, 1);
                reltype = combos(i, 2);

                [outData, relerr] = psyrat_relsummary( ...
                    'psyrat_data', baseData, ...
                    'analysis', 'trt', ...
                    'gcoeff', gcoeff, ...
                    'reltype', reltype, ...
                    'relcutoff', 0.2, ...
                    'meascutoff', 2, ...
                    'relcentmeas', 1, ...
                    'CI', 0.95);

                if relerr.nogooddata ~= 0
                    continue;
                end
                validCombos = validCombos + 1;
                if isempty(tableData)
                    tableData = outData;
                end

                ev = outData.relsummary.group(1).event(1);
                testCase.verifyGreaterThanOrEqual(ev.rel.m, 0);
                testCase.verifyLessThanOrEqual(ev.rel.m, 1);
                testCase.verifyGreaterThanOrEqual(ev.icc.m, 0);
                testCase.verifyLessThanOrEqual(ev.icc.m, 1);

                if gcoeff == 1
                    testCase.verifyEqual(outData.relsummary.gcoeff_name, 'dep');
                else
                    testCase.verifyEqual(outData.relsummary.gcoeff_name, 'gen');
                end

                if reltype == 1
                    testCase.verifyEqual(outData.relsummary.reltype_name, 'ic');
                elseif reltype == 2
                    testCase.verifyEqual(outData.relsummary.reltype_name, 'trt');
                else
                    testCase.verifyEqual(outData.relsummary.reltype_name, 'ic_trt');
                end
            end

            testCase.verifyGreaterThan(validCombos, 0, ...
                ['No TRT combinations produced valid summary data with this ', ...
                'integration fixture. Increase iterations or lower relcutoff.']);
            if validCombos == 0
                return;
            end

            % Validate TRT table generators on one representative successful combination.
            relCut = psyrat_trt_relcutofft('psyrat_data', tableData, 'gui', 0);
            relOverall = psyrat_trt_reloverallt('psyrat_data', tableData, 'gui', 0);
            varTbl = psyrat_variancet('psyrat_data', tableData, 'gui', 0);
            testCase.verifyGreaterThanOrEqual(height(relCut), 1);
            testCase.verifyGreaterThanOrEqual(height(relOverall), 1);
            testCase.verifyGreaterThanOrEqual(height(varTbl), 1);
        end
    end
end

function psyrat_data = localRunScenario(testCase, scenario)
psyrat_assume_dev_testdata(testCase);
dataTbl = PsyRATCmdStanIntegrationFactory.makeScenarioTable(testCase.projectRoot(), scenario);
outDir = localCmdStanOutputDir(testCase.projectRoot());

rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'chains', testCase.Chains, ...
    'iter', testCase.Iter, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir, ...
    'sserrvar', scenario.sserrvar, ...
    'diffest', scenario.diffest, ...
    'diffwpcov', scenario.diffwpcov, ...
    'diffrescor', scenario.diffrescor);

psyrat_data = struct();
psyrat_data.rel = rel;
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end
