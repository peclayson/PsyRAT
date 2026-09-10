classdef TestSummaryAndPresentationFunctions < PsyRATTestBase

    methods (Test)
        function testRelSummarySingBuildsSummaryAndErrorsStruct(testCase)
            psyrat_data = PsyRATTestDataFactory.makeRawSingRelData();

            [out, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing', ...
                'depcutoff', 0.70, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'CI', 0.95);

            testCase.verifyTrue(isfield(out, 'relsummary'));
            testCase.verifyTrue(isfield(relerr, 'trlcutoff'));
            testCase.verifyTrue(isfield(relerr, 'trlmax'));
            testCase.verifyTrue(isfield(relerr, 'nogooddata'));

            ev = out.relsummary.group(1).event(1);
            testCase.verifyTrue(isfield(ev, 'dep'));
            testCase.verifyTrue(isfield(ev, 'icc'));
            testCase.verifyTrue(isfield(ev, 'trlinfo'));
        end

        function testRelSummaryValidation(testCase)
            psyrat_data = PsyRATTestDataFactory.makeRawSingRelData();
            testCase.verifyError(@() psyrat_relsummary('psyrat_data', psyrat_data), ...
                'varargin:noanalysistype');
            testCase.verifyError(@() psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing', ...
                'depcutoff', 0.8, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'CI', 1.1), ...
                'varargin:ci');
        end

        function testDepCutoffAndOverallTables(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();

            cutoffTbl = psyrat_depcutofft('psyrat_data', psyrat_data, 'gui', 0);
            overallTbl = psyrat_depoverallt('psyrat_data', psyrat_data, 'gui', 0);

            testCase.verifyEqual(cutoffTbl.Properties.VariableNames, ...
                {'Label', 'Trial_Cutoff', 'Dependability'});
            testCase.verifyEqual(overallTbl.Properties.VariableNames, ...
                {'Label', 'n_Included', 'n_Excluded', 'Dependability', ...
                'Mean_Num_Trials', 'Med_Num_Trials', 'Std_Num_Trials', ...
                'Min_Num_Trials', 'Max_Num_Trials'});
            testCase.verifyEqual(height(cutoffTbl), 1);
            testCase.verifyEqual(height(overallTbl), 1);
        end
        
        function testDepCutoffAndOverallTablesGeneralizabilityLabels(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            psyrat_data.relsummary.gcoeff_name = 'gen';
            
            cutoffTbl = psyrat_depcutofft('psyrat_data', psyrat_data, 'gui', 0);
            overallTbl = psyrat_depoverallt('psyrat_data', psyrat_data, 'gui', 0);
            
            testCase.verifyEqual(cutoffTbl.Properties.VariableNames, ...
                {'Label', 'Trial_Cutoff', 'Generalizability'});
            testCase.verifyEqual(overallTbl.Properties.VariableNames, ...
                {'Label', 'n_Included', 'n_Excluded', 'Generalizability', ...
                'Mean_Num_Trials', 'Med_Num_Trials', 'Std_Num_Trials', ...
                'Min_Num_Trials', 'Max_Num_Trials'});
        end

        function testTrtCutoffAndOverallTables(testCase)
            psyrat_data = PsyRATTestDataFactory.makeTRTSummaryData();

            cutoffTbl = psyrat_trt_relcutofft('psyrat_data', psyrat_data, 'gui', 0);
            overallTbl = psyrat_trt_reloverallt('psyrat_data', psyrat_data, 'gui', 0);

            testCase.verifyEqual(cutoffTbl.Properties.VariableNames, ...
                {'Label', 'Trial_Cutoff', 'Dependability'});
            testCase.verifyEqual(overallTbl.Properties.VariableNames, ...
                {'Label', 'n_Included', 'n_Excluded', 'Dependability', ...
                'Mean_Num_Trials', 'Med_Num_Trials', 'Std_Num_Trials', ...
                'Min_Num_Trials', 'Max_Num_Trials'});
            testCase.verifyEqual(height(cutoffTbl), 1);
            testCase.verifyEqual(height(overallTbl), 1);
        end

        function testVarianceTableCreation(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            tbl = psyrat_variancet('psyrat_data', psyrat_data, 'gui', 0);

            %The one-facet variance table now reports the trial main-effect SD
            %(sigma_i) as its own column between the between-person and the
            %person x trial residual SDs (Rocha Table 2 / S10). ICC and SEM
            %carry the decision type (RC-11); this fixture is gcoeff 1 (dep).
            testCase.verifyEqual(tbl.Properties.VariableNames, ...
                {'Label', 'Between_StdDev', 'Trial_StdDev', 'Within_StdDev', ...
                'SEM_Dep', 'ICC_Dep'});
            testCase.verifyEqual(height(tbl), 1);
        end

        function testVarianceTableHeaderTracksEstimand(testCase)
            %RC-11 regression: the ICC/SEM VALUE changes with the coefficient
            %selector, so the HEADER must change with it too. Previously both
            %estimands were printed under a bare 'ICC'/'SEM', which made an
            %archived export non-reconstructible.
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();

            psyrat_data.relsummary.gcoeff = 1;
            psyrat_data.relsummary.gcoeff_name = 'dep';
            tblDep = psyrat_variancet('psyrat_data', psyrat_data, 'gui', 0);
            testCase.verifyEqual(tblDep.Properties.VariableNames{5}, 'SEM_Dep');
            testCase.verifyEqual(tblDep.Properties.VariableNames{6}, 'ICC_Dep');

            psyrat_data.relsummary.gcoeff = 2;
            psyrat_data.relsummary.gcoeff_name = 'gen';
            tblGen = psyrat_variancet('psyrat_data', psyrat_data, 'gui', 0);
            testCase.verifyEqual(tblGen.Properties.VariableNames{5}, 'SEM_Gen');
            testCase.verifyEqual(tblGen.Properties.VariableNames{6}, 'ICC_Gen');

            %structs saved before gcoeff_name existed still carry the numeric
            %selector, and the label must be recovered from it rather than
            %silently dropping back to an unqualified header
            legacy = rmfield(psyrat_data.relsummary, 'gcoeff_name');
            legacy.gcoeff = 2;
            psyrat_data.relsummary = legacy;
            tblLegacy = psyrat_variancet('psyrat_data', psyrat_data, 'gui', 0);
            testCase.verifyEqual(tblLegacy.Properties.VariableNames{5}, 'SEM_Gen');
            testCase.verifyEqual(tblLegacy.Properties.VariableNames{6}, 'ICC_Gen');
        end

        function testVarianceTableFallbackForICDiffSSErrVar(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            psyrat_data.rel.analysis = 'ic_diff_sserrvar';

            tbl = psyrat_variancet('psyrat_data', psyrat_data, 'gui', 0);

            testCase.verifyEqual(tbl.Properties.VariableNames, ...
                {'Label', 'Group_Between_StdDev', 'Group_Within_StdDev', ...
                'Group_SEM_Dep', 'Group_ICC_Dep'});
            testCase.verifyEqual(height(tbl), 1);
        end

        function testVarianceTableSemFallsOpenOnMismatchedDataStruct(testCase)
            %RC-20 guard. The new per-analysis SEM cases read draws by field
            %name, so a struct carrying one of those analysis labels while its
            %relsummary.data was built by a DIFFERENT branch must fall back to
            %psyrat_sem_string's residual SEM, not throw. That shape is not
            %hypothetical: rel.analysis is the estimation identity the viewers
            %preserve, and mapping it to the wrong summary builder is the exact
            %defect behind the earlier "Unrecognized field name ssrel_table"
            %crash (see testICDiffSSErrVarRequiresSubjectLevelSummaryLabel).
            %
            %Both relabeled shapes are covered: the one-facet fixture has no
            %b_sigma (difference case) and no pop_sdlog (sserr case). Without
            %the field guards this errors with MATLAB:nonExistentField and
            %takes the whole variance table down over one column.
            for a = {'ic_diff', 'ic_diff_sserrvar', 'trt_sserrvar'}
                psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
                psyrat_data.rel.analysis = a{1};

                tbl = psyrat_variancet('psyrat_data', psyrat_data, 'gui', 0);

                testCase.verifyEqual(height(tbl), 1, ...
                    sprintf('Relabeled %s must still render one row.', a{1}));
                %the SEM cell is populated (the fallback ran) rather than the
                %call having errored out
                semcol = tbl.Properties.VariableNames{4};
                testCase.verifyTrue(startsWith(semcol, {'SEM', 'Group_SEM'}), ...
                    sprintf('Column 4 should be the SEM for %s, got %s.', ...
                    a{1}, semcol));
                testCase.verifyNotEmpty(tbl.(semcol){1}, ...
                    sprintf('Relabeled %s must still print a SEM.', a{1}));
            end
        end

        function testVarianceTableSserrSquaresBetweenPersonSD(testCase)
            %RC-09 regression, VALUE test. gro_sds holds sigma_p as a standard
            %deviation, so the sserr branch of psyrat_variancet must square it
            %before forming the ICC (a ratio of variances) and must summarize it
            %on the SD scale for Group_Between_StdDev.
            %
            %Anchored at sigma_p = 4.5 with a residual SD of exactly 2 and a
            %trial main-effect SD of exactly 1, so both columns have closed
            %forms with no reliance on the production code:
            %   Group_Between_StdDev = 4.50
            %   Group_ICC (gcoeff 2, relative) = 4.5^2/(4.5^2 + 2^2)     = 0.8351
            %   Group_ICC (gcoeff 1, absolute) = 4.5^2/(4.5^2 + 2^2 + 1) = 0.8020
            %
            %sigma_p = 1 is deliberately avoided: it is the exact fixed point
            %where squaring and not squaring agree, so the defect is invisible
            %there. Pre-fix this branch reported 2.12 (= sqrt(4.5)) and 0.53.
            %Draws are constant, so every quantile equals the point estimate.
            nd = 6;
            sigP = 4.5; residSD = 2; trlSD = 1;

            psyrat_data = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
            [pd, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing_sserr', ...
                'depcutoff', 0.01, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'gcoeff', 1, ...
                'CI', 0.95);
            testCase.verifyEqual(relerr.nogooddata, 0);

            %overwrite the person block with known constants in both groups
            for g = 1:2
                pd.relsummary.data.g(g).e(1).gro_sds = num2cell(sigP * ones(nd, 1));
                pd.relsummary.data.g(g).e(1).pop_sdlog = log(residSD) * ones(nd, 1);
                pd.relsummary.data.g(g).e(1).sig_trl = trlSD * ones(nd, 1);
            end

            %relative (gcoeff 2): the trial main effect is excluded
            pd.relsummary.gcoeff = 2;
            pd.relsummary.gcoeff_name = 'gen';
            tblRel = psyrat_variancet('psyrat_data', pd, 'gui', 0);
            testCase.verifyEqual(strtrim(tblRel.Group_Between_StdDev{1}), ...
                '4.50 CI [4.50 4.50]');
            testCase.verifyEqual(strtrim(tblRel.Group_ICC_Gen{1}), ...
                '0.84 CI [0.84 0.84]');

            %absolute (gcoeff 1): sigma_i^2 joins the error, so the ICC drops
            pd.relsummary.gcoeff = 1;
            pd.relsummary.gcoeff_name = 'dep';
            tblAbs = psyrat_variancet('psyrat_data', pd, 'gui', 0);
            testCase.verifyEqual(strtrim(tblAbs.Group_Between_StdDev{1}), ...
                '4.50 CI [4.50 4.50]');
            testCase.verifyEqual(strtrim(tblAbs.Group_ICC_Dep{1}), ...
                '0.80 CI [0.80 0.80]');
        end

        function testVarianceTableTrtSserrSquaresBetweenPersonSD(testCase)
            %RC-09 twin for the crossed test-retest subject-level design
            %(trt_sserrvar, analysis 25), which reaches the same branch of
            %psyrat_variancet through the guard at :170. Same anchor as the
            %one-facet test above; gro_sds is draws x 2 here (person location
            %SD, person residual log-SD spread) so only column 1 is pinned.
            %
            %RC-06 group twin (2026-07-30) MOVED the expected ICC, and the value
            %is re-pinned rather than the test retired: its purpose is to catch
            %a bp that reaches the coefficient already squared, and that purpose
            %survives the estimand change intact. The old one-facet expectation
            %was sigma_p^2/(sigma_p^2 + sigma_pi,e^2) = 0.8351 -> '0.84'. The
            %group ICC is now the two-facet coefficient, whose denominator also
            %carries the occasion facet, giving 0.8257 -> '0.83'. Feeding
            %sigma_p^2 in where an SD belongs still yields 0.99, so the
            %regression this test exists for remains far outside tolerance.
            %
            %The CI is no longer degenerate ([0.82 0.83] rather than
            %[0.84 0.84]) because only sigma_p and the residual are pinned to
            %constants here; the crossed components come from the fixture and
            %vary across draws. That is expected, not a loss of control.
            nd = 8;
            sigP = 4.5; residSD = 2;

            psyrat_data = PsyRATTestDataFactory.makeTRTSSErrRelData();
            [pd, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'trt_sserr', ...
                'gcoeff', 2, ...
                'reltype', 3, ...
                'depcutoff', 0.01, ...
                'meascutoff', 2, ...
                'CI', 0.95);
            testCase.verifyEqual(relerr.nogooddata, 0);

            pd.relsummary.data.g(1).e(1).gro_sds = ...
                num2cell([sigP * ones(nd, 1), 0.4 * ones(nd, 1)]);
            pd.relsummary.data.g(1).e(1).pop_sdlog = log(residSD) * ones(nd, 1);
            pd.relsummary.gcoeff = 2;
            pd.relsummary.gcoeff_name = 'gen';

            tbl = psyrat_variancet('psyrat_data', pd, 'gui', 0);
            testCase.verifyEqual(strtrim(tbl.Group_Between_StdDev{1}), ...
                '4.50 CI [4.50 4.50]');
            testCase.verifyEqual(strtrim(tbl.Group_ICC_Gen{1}), ...
                '0.83 CI [0.82 0.83]');
        end

        function testVarianceTableIcDiffSemHonorsEstimand(testCase)
            %RC-20. 'ic_diff' was absent from the psyrat_sem_stats switch, so
            %its SEM fell through to psyrat_sem_string's witsd ./ sqrt(obs)
            %fallback: the person x trial,e residual alone, with no gcoeff
            %branch and no sigma_i, printed under the _Dep/_Gen header RC-11
            %added -- a positive claim the number did not honor. Observed live
            %on the GUI 2026-07-27: flipping the difference-score coefficient
            %moved the per-event ICCs while this column stayed byte-identical
            %and only its header changed.
            %
            %Closed forms, with constant draws so every quantile collapses to
            %the point estimate and n pinned by overwriting dep.meas (the field
            %psyrat_sem_obs reads for this analysis):
            %   SEM_Gen = sqrt(sigma_pi,e^2 / n)  = sqrt(2.0^2 / 4)       = 1.00
            %   SEM_Dep = sqrt((sigma_pi,e^2 + sigma_i^2) / n)
            %                                     = sqrt((4 + 2.25) / 4)  = 1.25
            %Both constants stay away from 1, the fixed point where squaring
            %and not squaring agree (see the comment in
            %testVarianceTableSserrSquaresBetweenPersonSD -- cited by method
            %name, because a line cite here decays on the next insertion). On
            %the fixture as summarized, the pre-fix column printed 0.22/0.23
            %under BOTH selectors.
            %
            %Built from the DoD fixture's difference mode, which emits a
            %genuine ic_diff structure, so this runs the real summary path
            %rather than a hand-assembled struct.
            residSD = 2; trlSD = 1.5; nObs = 4;
            colname = {'SEM_Dep', 'SEM_Gen'};
            expected = {'1.25 CI [1.25 1.25]', '1.00 CI [1.00 1.00]'};

            pd0 = PsyRATTestDataFactory.makeDoDRelData();
            vd = psyrat_dod_build_virtual_rel('psyrat_data', pd0, ...
                'mode', 'diff', 'groupidx', 1, 'eventpair', [1 4]);

            for gc = 1:2
                pd = localSingDiff(vd, gc);
                nd = numel(pd.relsummary.data.g(1).e(1).b_sigma.raw);
                for e = 1:2
                    pd.relsummary.data.g(1).e(e).b_sigma.raw = ...
                        num2cell(log(residSD) * ones(nd, 1));
                    pd.relsummary.data.g(1).e(e).sd_trl.raw = ...
                        num2cell(trlSD * ones(nd, 1));
                    pd.relsummary.group(1).event(e).dep.meas = nObs;
                end

                tbl = psyrat_variancet('psyrat_data', pd, 'gui', 0);
                for e = 1:2
                    testCase.verifyEqual(strtrim(tbl.(colname{gc}){e}), ...
                        expected{gc}, sprintf( ...
                        'ic_diff SEM, diffgcoeff %d, event %d.', gc, e));
                end
            end
        end

        function testVarianceTableIcDiffSserrGroupRowHonorsEstimand(testCase)
            %RC-20, extended. On the 'ic_diff_sserrvar' group row BOTH columns
            %were estimand-blind, not just the SEM. That branch of
            %psyrat_relsummary stores sd_id/b_sigma/trl_varcov rather than the
            %gro_sds/pop_sdlog pair, so it lands in the second display branch
            %of psyrat_variancet, which had no gcoeff test at all -- flipping
            %the coefficient changed the header of both columns and neither
            %number. sigma_i was available the whole time as the diagonal of
            %trl_varcov, stored as a VARIANCE.
            %
            %Closed forms with constant draws, sigma_p = 4.5, sigma_pi,e = 2,
            %sigma_i^2 = 2.25 and n = 4:
            %   ICC_Gen = 4.5^2 / (4.5^2 + 4)        = 0.8351 -> 0.84
            %   ICC_Dep = 4.5^2 / (4.5^2 + 4 + 2.25) = 0.7642 -> 0.76
            %   SEM_Gen = sqrt(4 / 4)                = 1.00
            %   SEM_Dep = sqrt(6.25 / 4)             = 1.25
            %Pre-fix this fixture printed ICC 0.54/0.48 and SEM 0.27/0.32 under
            %both selectors.
            sigP = 4.5; residSD = 2; trlVar = 2.25; nObs = 4;
            iccname = {'Group_ICC_Dep', 'Group_ICC_Gen'};
            semname = {'Group_SEM_Dep', 'Group_SEM_Gen'};
            iccexp = {'0.76 CI [0.76 0.76]', '0.84 CI [0.84 0.84]'};
            semexp = {'1.25 CI [1.25 1.25]', '1.00 CI [1.00 1.00]'};

            pd0 = PsyRATTestDataFactory.makeICDiffSSErrRelData();

            for gc = 1:2
                [pd, relerr] = psyrat_relsummary( ...
                    'psyrat_data', pd0, ...
                    'analysis', 'sing_diff_sserr', ...
                    'depcutoff', 0.5, ...
                    'meascutoff', 2, ...
                    'depcentmeas', 1, ...
                    'diffgcoeff', gc, ...
                    'CI', 0.95);
                testCase.assertEqual(relerr.nogooddata, 0);

                nd = numel(pd.relsummary.data.g(1).e(1).b_sigma.raw);
                for e = 1:2
                    pd.relsummary.data.g(1).e(e).sd_id.raw = ...
                        num2cell(sigP * ones(nd, 1));
                    pd.relsummary.data.g(1).e(e).b_sigma.raw = ...
                        num2cell(log(residSD) * ones(nd, 1));
                    pd.relsummary.data.g(1).e(e).trl_varcov.raw = ...
                        num2cell(trlVar * ones(nd, 1));
                    pd.relsummary.group(1).event(e).dep.meas = nObs;
                end

                tbl = psyrat_variancet('psyrat_data', pd, 'gui', 0);
                for e = 1:2
                    testCase.verifyEqual(strtrim(tbl.(iccname{gc}){e}), ...
                        iccexp{gc}, sprintf( ...
                        'ic_diff_sserrvar ICC, diffgcoeff %d, event %d.', gc, e));
                    testCase.verifyEqual(strtrim(tbl.(semname{gc}){e}), ...
                        semexp{gc}, sprintf( ...
                        'ic_diff_sserrvar SEM, diffgcoeff %d, event %d.', gc, e));
                end
            end
        end

        function testVarianceTableTrtSserrGroupRowIsTwoFacet(testCase)
            %RC-06 group twin, owner ruling 2026-07-30. Replaces
            %testVarianceTableTrtSserrSemMatchesItsOwnIcc, which is REWRITTEN
            %rather than re-pinned because the property it locked no longer
            %holds by design.
            %
            %That test asserted sigma_p^2 / (sigma_p^2 + SEM^2 * n) == ICC, i.e.
            %that SEM and ICC shared an error TERM. Under the ruling they share
            %the variance PARTITION but not n': the ICC is the two-facet
            %coefficient at obs = 1, n'_o = 1, while the SEM uses the trial-count
            %central tendency and relsummary.nocc, which is the convention the
            %plain-'trt' row already ships. Re-pinning the old identity with new
            %constants would have re-imposed the very coupling that was lifted.
            %
            %What is locked instead, all recomputed here rather than read from
            %production:
            %  1. the printed ICC IS the two-facet coefficient at obs=1,n'_o=1;
            %  2. it equals gro_icc, which the SAME run already stores for the
            %     caterpillar plot -- before the fix those two disagreed
            %     (0.6919 in the table vs 0.6288 on the plot at CES);
            %  3. both columns MOVE with reltype. That is the structural half:
            %     the old expression could not read reltype at all, so CE, CS
            %     and CES printed one number.
            pd0 = PsyRATTestDataFactory.makeTRTSSErrRelData();
            iccPrinted = zeros(1,3);
            semPrinted = zeros(1,3);

            for rt = 1:3
                [pd, relerr] = psyrat_relsummary( ...
                    'psyrat_data', pd0, ...
                    'analysis', 'trt_sserr', ...
                    'gcoeff', 1, ...
                    'reltype', rt, ...
                    'depcutoff', 0.01, ...
                    'meascutoff', 2, ...
                    'CI', 0.95);
                testCase.assertEqual(relerr.nogooddata, 0);

                tbl = psyrat_variancet('psyrat_data', pd, 'gui', 0);
                iccPrinted(rt) = localPointEstimate(tbl.Group_ICC_Dep{1});
                semPrinted(rt) = localPointEstimate(tbl.Group_SEM_Dep{1});

                d = pd.relsummary.data.g(1).e(1);
                bp_sd = cell2mat(d.gro_sds(:,1));

                %(1) independent two-facet coefficient at a single trial on a
                %single occasion, built from the field mapping psyrat_ssrel_trt
                %uses, not from psyrat_variancet
                [~, iccExpected, ~] = psyrat_rel_trt( ...
                    'gcoeff', 1, 'reltype', rt, ...
                    'bp', bp_sd, 'bo', d.sig_occ, 'bt', d.sig_trl, ...
                    'txp', d.sig_trlxid, 'oxp', d.sig_occxid, ...
                    'txo', d.sig_trlxocc, 'err', exp(d.pop_sdlog), ...
                    'obs', 1, 'nocc', 1, 'CI', 0.95);
                testCase.verifyEqual(iccPrinted(rt), round(iccExpected,2), ...
                    'AbsTol', 5e-3, sprintf(['trt_sserrvar group ICC must be ' ...
                    'the two-facet coefficient at obs=1, n''_o=1 (reltype %d)'], rt));

                %(2) and it must agree with the number the same run already
                %carries for its own caterpillar plot
                sst = pd.relsummary.group(1).event(1).ssrel_table;
                testCase.verifyEqual(iccPrinted(rt), round(sst.gro_icc(1),2), ...
                    'AbsTol', 5e-3, sprintf(['the variance table and the ' ...
                    'subject-level plot must report ONE group ICC (reltype %d)'], rt));

                %the SEM must sit on the same partition, at the SEM's own n'
                obs = pd.relsummary.group(1).event(1).trlinfo.mean;
                nocc = 1;
                if isfield(pd.relsummary,'nocc') && ~isempty(pd.relsummary.nocc)
                    nocc = pd.relsummary.nocc;
                end
                txp = d.sig_trlxid(:).^2;  oxp = d.sig_occxid(:).^2;
                bt  = d.sig_trl(:).^2;     bo  = d.sig_occ(:).^2;
                txo = d.sig_trlxocc(:).^2; err = exp(d.pop_sdlog(:)).^2;
                switch rt
                    case 1
                        e2 = txp./obs + err./(obs*nocc) + bt./obs + txo./(obs*nocc);
                    case 2
                        e2 = oxp./nocc + err./(obs*nocc) + bo./nocc + txo./(obs*nocc);
                    case 3
                        e2 = txp./obs + oxp./nocc + err./(obs*nocc) + ...
                            bt./obs + bo./nocc + txo./(obs*nocc);
                end
                testCase.verifyEqual(semPrinted(rt), round(mean(sqrt(e2)),2), ...
                    'AbsTol', 5e-3, sprintf(['trt_sserrvar group SEM must be ' ...
                    'the two-facet error term at the reported n'' (reltype %d)'], rt));
            end

            %(3) the structural assertion. Before the fix the ICC expression
            %could not read reltype, so all three printed the same number, and
            %the SEM was on the one-facet branch and likewise could not.
            testCase.verifyNotEqual(iccPrinted(1), iccPrinted(3), ...
                'group ICC must respond to the CE/CES coefficient selection');
            testCase.verifyNotEqual(iccPrinted(2), iccPrinted(3), ...
                'group ICC must respond to the CS/CES coefficient selection');
            testCase.verifyNotEqual(semPrinted(1), semPrinted(3), ...
                'group SEM must respond to the CE/CES coefficient selection');
            testCase.verifyNotEqual(semPrinted(2), semPrinted(3), ...
                'group SEM must respond to the CS/CES coefficient selection');
        end

        function testVarianceTableTrtSserrLabelsTheSemOccasionEstimand(testCase)
            %RC-06 group twin, second half. Found by adversarially reviewing the
            %FIX, not the code it replaced -- the defect did not exist before
            %this commit.
            %
            %Until the twin was closed this table carried no n'_o term anywhere:
            %the SEM came from the one-facet branch, so BOTH occasion settings
            %printed 0.60 and the absent occasion label was inert. Giving the SEM
            %the plain-'trt' convention makes n'_o load-bearing here for the
            %first time -- 0.88 at n'_o = 1 versus 0.66 at n'_o = 2 on this
            %fixture -- while the ICC stays pinned at n'_o = 1. A silent 25%
            %swing under identical on-screen text is the exact fault RC-20 and
            %RC-27/RC-29 exist to prevent.
            %
            %The label names the SEM's n'_o and the ICC's SEPARATELY because they
            %genuinely differ on this row. psyrat_trt_occlabel's wording is about
            %"the reported coefficient", so reusing the {'trt','trt_diff'} clause
            %would attribute n'_o = k to a coefficient computed at 1.
            %
            %This needs gui = 1. Every other variance-table test runs gui = 0,
            %which never constructs the figure, so the window title was
            %untestable by construction and the gap was invisible to the suite.
            pd0 = PsyRATTestDataFactory.makeTRTSSErrRelData();
            sem = cell(1,2);
            figname = cell(1,2);

            for mode = 1:2
                args = {'psyrat_data', pd0, 'analysis', 'trt_sserr', ...
                    'gcoeff', 1, 'reltype', 3, 'depcutoff', 0.01, ...
                    'meascutoff', 2, 'CI', 0.95};
                if mode == 2
                    args = [args, {'noccmode', 2, 'nocc', 2}]; %#ok<AGROW>
                end
                [pd, relerr] = psyrat_relsummary(args{:});
                testCase.assertEqual(relerr.nogooddata, 0);

                tbl = psyrat_variancet('psyrat_data', pd, 'gui', 0);
                sem{mode} = strtrim(tbl.Group_SEM_Dep{1});

                close all force;
                psyrat_variancet('psyrat_data', pd, 'gui', 1);
                f = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
                testCase.assertNotEmpty(f, 'variance table figure was not created');
                figname{mode} = f(1).Name;
                close all force;
            end

            %the premise: the SEM really does move with the occasion estimand
            testCase.verifyNotEqual(sem{1}, sem{2}, ...
                'trt_sserrvar SEM should depend on the selected occasion estimand');

            %so the window must say which n'_o produced it, and must not
            %attribute that n'_o to the ICC printed beside it
            testCase.verifyTrue(contains(string(figname{1}), "SEM at n'_o = 1"), ...
                'single-occasion SEM must be labeled on the window');
            testCase.verifyTrue(contains(string(figname{2}), "SEM at n'_o = 2"), ...
                'multi-occasion composite SEM must be labeled with its own k');
            testCase.verifyTrue(contains(string(figname{2}), "ICC at n'_o = 1"), ...
                'the label must not attribute the SEM n''_o to the pinned ICC');
            testCase.verifyNotEqual(figname{1}, figname{2}, ...
                'two runs printing different SEMs must not share a window title');
        end

        function testSubjectLevelDiffTablesAndPlots(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICDiffSSErrRelData();
            [psyrat_data, ~] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing_diff_sserr', ...
                'depcutoff', 0.5, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'diffgcoeff', 1, ...
                'CI', 0.95);

            overallTbl = psyrat_depoverallt('psyrat_data', psyrat_data, 'gui', 0);
            cutoffTbl = psyrat_depcutofft('psyrat_data', psyrat_data, 'gui', 0);
            varianceTbl = psyrat_variancet('psyrat_data', psyrat_data, 'gui', 0);
            psyrat_ssrelplot('psyrat_data', psyrat_data, 'stat', 'dep');
            psyrat_ssrelplot('psyrat_data', psyrat_data, 'stat', 'icc');

            testCase.verifyGreaterThanOrEqual(height(overallTbl), 1);
            testCase.verifyGreaterThanOrEqual(height(cutoffTbl), 1);
            testCase.verifyGreaterThanOrEqual(height(varianceTbl), 1);
        end

        function testSubjectLevelDiffOverallTableUsesDiffScoreWhenEventTablesMissing(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICDiffSSErrRelData();
            [psyrat_data, ~] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing_diff_sserr', ...
                'depcutoff', 0.5, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'diffgcoeff', 1, ...
                'CI', 0.95);
            psyrat_data.relsummary.group(1).event = ...
                rmfield(psyrat_data.relsummary.group(1).event, 'ssrel_table');

            overallTbl = psyrat_depoverallt('psyrat_data', psyrat_data, 'gui', 0);

            testCase.verifyEqual(height(overallTbl), 1);
            testCase.verifyEqual(overallTbl.Label{1}, 'diff score');
        end

        function testICDiffSSErrVarRequiresSubjectLevelSummaryLabel(testCase)
            %Guards the viewer-dispatch contract behind the "Unrecognized field
            %name ssrel_table" crash. A subject-level difference-score run
            %(rel.analysis == 'ic_diff_sserrvar') must be summarized with the
            %'sing_diff_sserr' builder, which is the only writer of
            %diffscore.ssrel_table (psyrat_relsummary.m). The single-facet
            %viewer psyrat_startview_sing previously mislabeled these runs as
            %'sing_diff' (the non-subject-level builder), which omits
            %ssrel_table and crashes psyrat_depoverallt. This test does not
            %drive the GUI callback itself (GUI callbacks are not headlessly
            %testable); it locks the downstream invariant the callback must
            %honor.

            %Correct label: builds diffscore.ssrel_table and the overall table
            %renders a difference-score row without error.
            psyrat_data = PsyRATTestDataFactory.makeICDiffSSErrRelData();
            [good, ~] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing_diff_sserr', ...
                'depcutoff', 0.5, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'diffgcoeff', 1, ...
                'CI', 0.95);
            testCase.verifyTrue( ...
                isfield(good.relsummary.group(1).diffscore, 'ssrel_table'), ...
                'sing_diff_sserr must populate diffscore.ssrel_table.');
            overallTbl = psyrat_depoverallt('psyrat_data', good, 'gui', 0);
            testCase.verifyTrue(any(strcmp(overallTbl.Label, 'diff score')), ...
                'Overall table should include a difference-score row.');

            %Failure state the wrong viewer label produced: a diffscore for an
            %ic_diff_sserrvar run that lacks ssrel_table (the non-subject-level
            %'sing_diff' builder omits it). rel.analysis stays 'ic_diff_sserrvar'
            %(the estimation identity the viewer preserves), so psyrat_depoverallt
            %takes the subject-level branch and crashes on the missing field.
            %ssrel_table is removed directly here because the 'sing_diff'
            %builder needs raw event scores (REL.out.b) that the subject-level
            %factory does not synthesize; this reproduces the same crash state.
            bad = good;
            bad.relsummary.group(1).diffscore = ...
                rmfield(bad.relsummary.group(1).diffscore, 'ssrel_table');
            testCase.verifyEqual(bad.rel.analysis, 'ic_diff_sserrvar');
            testCase.verifyError( ...
                @() psyrat_depoverallt('psyrat_data', bad, 'gui', 0), ...
                'MATLAB:nonExistentField');
        end

        function testSubjectLevelOverallTableDashesAZeroedGroup(testCase)
            %G44. Post-G39 strict intersection, a group can retain nobody
            %while its ssrel_table still holds informative per-subject rows.
            %The sserr overall table used to render such a group as
            %' NaN SD: NaN' with empty min/max cells; it must dash the six
            %summary cells instead (the diff rows' existing '---' convention)
            %and keep the complementing 0 / N counts numeric. The un-doctored
            %group in the same table is the positive control: its cells must
            %stay numeric, so the dash cannot come from a branch that fires
            %on every row.
            psyrat_data = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
            [psyrat_data, ~] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing_sserr', ...
                'depcutoff', 0.1, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'CI', 0.95);

            %anti-vacuity: both groups genuinely retained someone before the
            %zeroed state is staged, so the dashes below are attributable to
            %the doctoring rather than to the fixture.
            testCase.assertNotEmpty(psyrat_data.relsummary.group(1).goodids);
            testCase.assertNotEmpty(psyrat_data.relsummary.group(2).goodids);

            %stage the newly reachable state: group 1 retains nobody (the G39
            %fold empties the group's retained set through a dead sibling
            %event; only goodids/badids carry that state into this renderer)
            roster = psyrat_data.relsummary.group(1).event(1).ssrel_table.id;
            psyrat_data.relsummary.group(1).goodids = cell(0, 1);
            psyrat_data.relsummary.group(1).badids = roster;

            overallTbl = psyrat_depoverallt('psyrat_data', psyrat_data, 'gui', 0);
            vn = overallTbl.Properties.VariableNames;

            %row 1 = the zeroed group: columns 4-9 (coefficient + the five
            %trial summaries) dash out; the counts stay numeric and complement
            for col = 4:9
                testCase.verifyEqual(overallTbl.(vn{col}){1}, '---', sprintf( ...
                    'zeroed group: column %s must dash out', vn{col}));
            end
            testCase.verifyEqual(overallTbl.(vn{2}){1}, 0, ...
                'zeroed group: n Included must stay the numeric 0');
            testCase.verifyEqual(overallTbl.(vn{3}){1}, numel(roster), ...
                'zeroed group: n Excluded must complement to the roster');

            %row 2 = the surviving group: everything stays numeric. The exact
            %count is also the G53 byte-identity pin for the ic_sserrvar
            %route: the rendered goodn is the group-carrier count, and on
            %this route the carrier IS the event's include set, so the
            %carrier count and the event indicator summed over the carrier
            %rows must agree -- the weld below pins that equivalence, and
            %the render pin then holds under either spelling.
            g2tbl = psyrat_data.relsummary.group(2).event(1).ssrel_table;
            g2good = psyrat_data.relsummary.group(2).goodids;
            testCase.assertEqual( ...
                sum(g2tbl.ind2include(ismember(g2tbl.id, g2good))), ...
                numel(g2good), ...
                'ic_sserrvar weld: the carrier must equal the event include set');
            testCase.verifyEqual(overallTbl.(vn{2}){2}, numel(g2good), ...
                'surviving group: n Included must be the exact carrier count');
            testCase.verifyTrue(contains(overallTbl.(vn{4}){2}, 'SD:'), ...
                'surviving group: the coefficient cell must stay the numeric render');
            for col = 5:9
                testCase.verifyTrue(isnumeric(overallTbl.(vn{col}){2}), sprintf( ...
                    'surviving group: column %s must stay numeric', vn{col}));
            end
        end

        function testSubjectLevelDiffGroupRowConsumesResidualCrossCov(testCase)
            %Locks the case-8 contract the gamma concurrent path was violating:
            %REL.out.err_varcov(:,2,1) is the GROUP-level residual cross-covariance,
            %and psyrat_relsummary must feed it to BOTH the group diffscore
            %(psyrat_diffrel) and the per-subject rows (psyrat_ssrel_diff). The gamma
            %case-8 assembly used to rebuild err_varcov diagonal-only, so a concurrent
            %run reported the NON-concurrent group coefficient while its own
            %per-participant rows used the coupling -- an internal contradiction no
            %test caught. This drives the real summary path with the cross-cov present
            %vs zeroed and pins the direction of both.
            args = {'analysis','sing_diff_sserr','depcutoff',0.5,'meascutoff',2, ...
                'depcentmeas',1,'diffgcoeff',1,'CI',0.95};

            withCov = PsyRATTestDataFactory.makeICDiffSSErrRelData();
            evc = cell2mat(withCov.rel.out.err_varcov{:,1});
            testCase.assertGreaterThan(evc(1,1,2), 0, ...
                'Fixture must carry a positive residual cross-cov for this test.');

            %same fixture with ONLY the residual cross-covariance zeroed
            noCov = withCov;
            evc0 = evc;
            evc0(:,1,2) = 0;
            evc0(:,2,1) = 0;
            noCov.rel.out.err_varcov(:,1) = {num2cell(evc0)};
            noCov.rel.out.wp_cov_ss(:,1) = {zeros(size(withCov.rel.out.wp_cov_ss{:,1}))};

            [a, ~] = psyrat_relsummary('psyrat_data', withCov, args{:});
            [b, ~] = psyrat_relsummary('psyrat_data', noCov,  args{:});

            %(1) GROUP row: a positive residual cross-cov shrinks the difference error
            %    variance (rel_err subtracts 2*wp_cov/harmmean(obs)), so it must RAISE
            %    the coefficient. If err_varcov's off-diagonal were dropped, these two
            %    would be identical.
            testCase.verifyGreaterThan( ...
                a.relsummary.group(1).diffscore.pt, ...
                b.relsummary.group(1).diffscore.pt, ...
                ['A positive residual cross-covariance must raise the GROUP '...
                 'difference coefficient; equality means err_varcov(:,2,1) was ignored.']);

            %(2) PER-SUBJECT rows must move the SAME way. The bug was that the group
            %    row said rho_e = 0 while the subject rows said otherwise.
            ta = a.relsummary.group(1).diffscore.ssrel_table;
            tb = b.relsummary.group(1).diffscore.ssrel_table;
            testCase.verifyEqual(height(ta), height(tb));
            testCase.verifyTrue(all(ta.rel_pt > tb.rel_pt), ...
                ['Per-subject difference coefficients must rise with the residual '...
                 'cross-covariance, in the same direction as the group row.']);
        end

        function testPointIntervalPlotAndValidation(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            h = psyrat_ptintervalplot('psyrat_data', psyrat_data, 'stat', 'ICC', 'CI', 0.95);

            testCase.verifyTrue(all(isgraphics(h)));
            testCase.verifyError(@() psyrat_ptintervalplot('psyrat_data', psyrat_data, 'stat', 'bad'), ...
                'varargin:stat');
        end

        function testDepTrialsPlotAndValidation(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            fig = psyrat_depvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [1 6], ...
                'depcutoff', 0.8, ...
                'depline', 2, ...
                'CI', 0.95);

            testCase.verifyTrue(isgraphics(fig, 'figure'));
            
            figGen = psyrat_depvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [1 6], ...
                'depcutoff', 0.8, ...
                'depline', 2, ...
                'gcoeff', 2, ...
                'CI', 0.95);
            testCase.verifyTrue(contains(figGen.Name, 'Generalizability'));
            
            testCase.verifyError(@() psyrat_depvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [1 6], ...
                'depcutoff', 0.8, ...
                'depline', 7), ...
                'varargin:depline');
            testCase.verifyError(@() psyrat_depvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [1 6], ...
                'depcutoff', 0.8, ...
                'depline', 2, ...
                'gcoeff', 9), ...
                'varargin:gcoeff');
        end

        function testTrtTrialsPlotAndValidation(testCase)
            psyrat_data = PsyRATTestDataFactory.makeTRTSummaryData();
            fig = psyrat_trt_relvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [1 6], ...
                'relcutoff', 0.8, ...
                'relline', 2, ...
                'CI', 0.95);

            testCase.verifyTrue(isgraphics(fig, 'figure'));
            testCase.verifyError(@() psyrat_trt_relvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [1 6], ...
                'relcutoff', 0.8, ...
                'relline', 99), ...
                'varargin:relline');
        end

        function testDepTrialsPlotHonorsRangeStartAboveOne(testCase)
            %RC-30. 'trials' is documented as a [min max] range (help line 12)
            %and psyrat_rel_sing deliberately aligns its return vector to
            %obs(1):obs(2) (psyrat_rel_sing.m:169), but this plot indexed its
            %buffer by ABSOLUTE trial number, so the buffer grew to trials(2)
            %rows while x had only trials(2)-trials(1)+1 elements and
            %plot(x,plotrel) threw MATLAB:samelen for any trials(1) > 1.
            %Separately, the y-axis upper limit read trials(1) instead of 1.
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();

            fig = psyrat_depvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [5 20], ...
                'depcutoff', 0.8, ...
                'depline', 2, ...
                'CI', 0.95);
            testCase.verifyTrue(isgraphics(fig, 'figure'), ...
                'A range starting above 1 must draw, not error.');

            curves = localFindCurveLines(fig, 5:20);
            testCase.verifyNumElements(curves, 1, ...
                'This fixture is one group and one event, so exactly one curve.');
            curve = curves(1);
            testCase.verifyEqual(numel(curve.YData), 16, ...
                ['The curve must carry one point per requested trial count '...
                 '(20-5+1 = 16), not one per trial number up to the maximum.']);

            %the coefficient axis runs 0 to 1; before the fix its upper limit
            %was trials(1), i.e. 5 here
            ylimits = get(ancestor(curve, 'axes'), 'YLim');
            testCase.verifyEqual(ylimits(2), 1, ...
                sprintf(['y-axis upper limit must be 1 (a reliability '...
                 'coefficient), not trials(1) = 5. Observed %g.'], ylimits(2)));

            %Alignment, not merely length: row 1 of the buffer must be
            %trials(1). If the range were misaligned by the offset, the curve
            %would still have 16 points but would report the coefficients for
            %trials 1-16 under x-values 5-20. Compare against the same
            %quantities computed over the full range, which is deterministic
            %for depline = 2 (a posterior mean, no resampling).
            figFull = psyrat_depvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [1 20], ...
                'depcutoff', 0.8, ...
                'depline', 2, ...
                'CI', 0.95);
            full = localFindCurveLines(figFull, 1:20);
            testCase.verifyEqual(curve.YData, full(1).YData(5:20), 'AbsTol', 1e-12, ...
                ['Trials 5-20 of the ranged curve must equal trials 5-20 of '...
                 'the full curve; a mismatch means the buffer is offset.']);

            %Each depline branch writes the buffer through its OWN indexed
            %assignment, so exercising only the point estimate leaves two of the
            %three edited lines in this file unlocked at a ranged start.
            for dl = [1 3]
                figLine = psyrat_depvtrialsplot( ...
                    'psyrat_data', psyrat_data, ...
                    'trials', [5 20], ...
                    'depcutoff', 0.8, ...
                    'depline', dl, ...
                    'CI', 0.95);
                testCase.verifyNumElements(localFindCurveLines(figLine, 5:20), 1, ...
                    sprintf('depline = %d must draw one 16-point curve over [5 20].', dl));
            end
        end

        function testTrtTrialsPlotHonorsRangeStartAboveOne(testCase)
            %RC-30, test-retest twin. Same two defects (at different lines in
            %each file - only the buffer lines coincide); the finding as filed
            %named only this file, but psyrat_depvtrialsplot carried them too.
            %Both are fixed together, so both are locked together.
            %psyrat_rel_trt.m:428 is the alignment convention here.
            psyrat_data = PsyRATTestDataFactory.makeTRTSummaryData();

            fig = psyrat_trt_relvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [5 20], ...
                'relcutoff', 0.8, ...
                'relline', 2, ...
                'CI', 0.95);
            testCase.verifyTrue(isgraphics(fig, 'figure'), ...
                'A range starting above 1 must draw, not error.');

            curves = localFindCurveLines(fig, 5:20);
            testCase.verifyNumElements(curves, 1, ...
                'This fixture is one group and one event, so exactly one curve.');
            curve = curves(1);
            testCase.verifyEqual(numel(curve.YData), 16, ...
                ['The curve must carry one point per requested trial count '...
                 '(20-5+1 = 16), not one per trial number up to the maximum.']);

            ylimits = get(ancestor(curve, 'axes'), 'YLim');
            testCase.verifyEqual(ylimits(2), 1, ...
                sprintf(['y-axis upper limit must be 1 (a reliability '...
                 'coefficient), not trials(1) = 5. Observed %g.'], ylimits(2)));

            figFull = psyrat_trt_relvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [1 20], ...
                'relcutoff', 0.8, ...
                'relline', 2, ...
                'CI', 0.95);
            full = localFindCurveLines(figFull, 1:20);
            testCase.verifyEqual(curve.YData, full(1).YData(5:20), 'AbsTol', 1e-12, ...
                ['Trials 5-20 of the ranged curve must equal trials 5-20 of '...
                 'the full curve; a mismatch means the buffer is offset.']);

            %as above: all three relline branches have their own indexed
            %assignment into the buffer, so all three need a ranged start
            for rl = [1 3]
                figLine = psyrat_trt_relvtrialsplot( ...
                    'psyrat_data', psyrat_data, ...
                    'trials', [5 20], ...
                    'relcutoff', 0.8, ...
                    'relline', rl, ...
                    'CI', 0.95);
                testCase.verifyNumElements(localFindCurveLines(figLine, 5:20), 1, ...
                    sprintf('relline = %d must draw one 16-point curve over [5 20].', rl));
            end
        end

        function testDepTrialsPlotRangeAcrossGroupsAndEvents(testCase)
            %RC-30, multi-column path. The two tests above run a one-group,
            %one-event fixture, so they never exercise the case where plotrel
            %grows a SECOND COLUMN, nor the case where the figure holds several
            %subplots each carrying several curves. Two groups x two events
            %covers both, and is the configuration in which a first-match-wins
            %curve lookup would silently check one curve out of four.
            raw = PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();
            [psyrat_data, relerr] = psyrat_relsummary( ...
                'psyrat_data', raw, ...
                'analysis', 'sing', ...
                'depcutoff', 0.01, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'gcoeff', 1, ...
                'CI', 0.95);
            testCase.assertEqual(relerr.nogooddata, 0, ...
                'Fixture must yield usable data or the plot assertions are vacuous.');

            fig = psyrat_depvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [5 20], ...
                'depcutoff', 0.8, ...
                'depline', 2, ...
                'gcoeff', 1, ...
                'CI', 0.95);

            %2 groups x 2 events = 4 curves, each with one point per trial count
            curves = localFindCurveLines(fig, 5:20);
            testCase.verifyNumElements(curves, 4, ...
                ['Two groups x two events must yield four curves; fewer means '...
                 'a group column or an event subplot was dropped.']);

            %every subplot's y-axis is a coefficient axis, not a trial-count axis
            for k = 1:numel(curves)
                yl = get(ancestor(curves(k), 'axes'), 'YLim');
                testCase.verifyEqual(yl(2), 1, ...
                    sprintf('Curve %d: y-axis upper limit must be 1, observed %g.', ...
                        k, yl(2)));
            end

            %same alignment invariant as the single-group tests, applied to the
            %whole set. Rows are sorted so the assertion does not depend on the
            %order findall happens to return handles in.
            figFull = psyrat_depvtrialsplot( ...
                'psyrat_data', psyrat_data, ...
                'trials', [1 20], ...
                'depcutoff', 0.8, ...
                'depline', 2, ...
                'gcoeff', 1, ...
                'CI', 0.95);
            fullCurves = localFindCurveLines(figFull, 1:20);
            testCase.assertNumElements(fullCurves, 4);

            ranged = sortrows(vertcat(curves.YData));
            fullTail = vertcat(fullCurves.YData);
            testCase.verifyEqual(ranged, sortrows(fullTail(:, 5:20)), ...
                'AbsTol', 1e-12, ...
                ['Every group x event curve over [5 20] must equal trials 5-20 '...
                 'of its own full-range curve.']);

            %What this fixture can and cannot distinguish, stated so the count
            %above is not read as more than it is. sig_u, sig_e and sig_trl vary
            %with GROUP only (PsyRATTestDataFactory.m:488-490: 0.03*(g-1),
            %0.02*(g-1), and a constant); only mu varies with event, and mu does
            %not enter a reliability coefficient. So the two EVENTS are
            %numerically identical by construction and no assertion here can tell
            %them apart. What is provable is that the two GROUPS are distinct --
            %i.e. the second plotrel column carries that group's own curve rather
            %than a copy of the first, which is the column-dimension claim the
            %single-group tests cannot make at all.
            distinct = unique(round(vertcat(curves.YData), 12), 'rows');
            testCase.verifySize(distinct, [2 16], ...
                ['Expected exactly two distinct curves (two groups; the two '...
                 'events are identical by fixture construction). One row means '...
                 'a group column was duplicated; four means the fixture changed '...
                 'and this comment is stale.']);
        end

        function testSubjectLevelPlots(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            psyrat_data.relsummary.group(1).event(1).ssrel_table.dep_pt(:) = 0.50;
            psyrat_data.relsummary.group(1).event(1).ssrel_table.dep_ll(:) = 0.30;
            psyrat_data.relsummary.group(1).event(1).ssrel_table.dep_ul(:) = 0.70;
            psyrat_data.relsummary.group(1).event(1).ssrel_table.icc_pt(:) = 0.45;
            psyrat_data.relsummary.group(1).event(1).ssrel_table.icc_ll(:) = 0.25;
            psyrat_data.relsummary.group(1).event(1).ssrel_table.icc_ul(:) = 0.65;

            %psyrat_depssplot, once smoke-tested here, was deleted in the
            %pre-beta batch (2026-09-05): no production caller, a stale fork
            %of psyrat_depvtrialsplot, and unable to draw two groups.
            psyrat_ssrelplot('psyrat_data', psyrat_data, 'stat', 'dep');
            psyrat_ssrelplot('psyrat_data', psyrat_data, 'stat', 'icc');
            testCase.verifyGreaterThanOrEqual(numel(findall(groot, 'Type', 'figure')), 2);
        end

        function testTwoFacetSubjectLevelPlotsUseTwoFacetReferenceLine(testCase)
            % RC-06 companion. psyrat_ssrelplot's own bp_var/pop_errvar
            % expression is a ONE-FACET reference line. Under a two-facet
            % (test-retest) table it would be a different quantity from the bars
            % it is drawn against, so psyrat_ssrel_trt now carries gro_rel and
            % gro_icc and the plot reads those. Both existing ssrelplot tests use
            % one-facet tables, so without this the new branch is never executed.
            n = 200;
            idtable = table({'s1';'s2';'s3'}, [4;9;16], ...
                'VariableNames', {'id','trls'});
            ss = psyrat_ssrel_trt('gcoeff',1,'reltype',3, ...
                'bp', repmat(sqrt(2),n,1), 'bo', repmat(sqrt(3),n,1), ...
                'bt', repmat(sqrt(5),n,1), 'txp', repmat(sqrt(7),n,1), ...
                'oxp', repmat(sqrt(11),n,1), 'txo', repmat(sqrt(13),n,1), ...
                'err', repmat(sqrt(17),n,1), ...
                'idtable', idtable, 'CI', 0.95, 'nocc', 2);

            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            psyrat_data.relsummary.group(1).event(1).ssrel_table = ss;

            % the one-facet expressions the plot used before, for contrast
            oneFacetDep = ss.bp_var(1) / ...
                (ss.bp_var(1) + ss.pop_errvar(1)/mean(ss.trls));
            oneFacetIcc = ss.bp_var(1) / (ss.bp_var(1) + ss.pop_errvar(1));

            spec = {'dep', ss.gro_rel(1), oneFacetDep; ...
                    'icc', ss.gro_icc(1), oneFacetIcc};
            for k = 1:size(spec,1)
                psyrat_ssrelplot('psyrat_data', psyrat_data, 'stat', spec{k,1});
                fig = gcf;
                % the reference line is the only Line at LineWidth 1.5; the
                % caterpillar bars are errorbar objects, not lines
                h = findobj(fig, 'Type', 'line', 'LineWidth', 1.5);
                testCase.verifyNotEmpty(h, ...
                    sprintf('no reference line found on the %s panel', spec{k,1}));
                y = get(h(1), 'YData');
                testCase.verifyEqual(y(1), spec{k,2}, 'AbsTol', 1e-12, sprintf( ...
                    '%s reference line is not the stored two-facet value', spec{k,1}));
                testCase.verifyNotEqual(y(1), spec{k,3}, sprintf( ...
                    '%s reference line still equals the one-facet expression', spec{k,1}));
                close(fig);
            end
        end

        function testRelFiguresAndTraceplotValidation(testCase)
            testCase.verifyError(@() psyrat_relfigures('analysis', 'sing'), 'varargin:nofile');
            testCase.verifyError(@() psyrat_checktraceplots([]), 'REL:notfound');
        end
        
        function testCriterionFiguresSingAndTrt(testCase)
            prefs = psyrat_defaults();
            prefs.view.plotcriterion = 0;
            prefs.view.tablecriterion = 1;
            prefs.view.ntrials = 6;
            prefs.view.depcentmeas = 1;
            prefs.view.relcentmeas = 1;
            prefs.view.reltype = 3;
            
            singData = PsyRATTestDataFactory.makeICSummaryData();
            prefs.view.gcoeff = 2; %criterion outputs should ignore this and use dependability
            cut1 = 0.1;
            cut2 = 1.4;
            
            close all force;
            prefs.view.criterioncutoff = cut1;
            psyrat_criterionfigures( ...
                'psyrat_prefs', prefs, ...
                'psyrat_data', singData, ...
                'analysis', 'sing');
            [ptSing1, singName1] = localCriterionPointEstimateFromOutputFigure();
            testCase.verifyTrue(contains(string(singName1), 'Dependability'));
            
            dSing = singData.relsummary.data.g(1).e(1);
            obsSing = max(1, round(singData.relsummary.group(1).event(1).trlinfo.mean));
            [~, expSing1, ~] = psyrat_cutscore_sing( ...
                'gcoeff', 1, ...
                'bp', dSing.sig_u.raw, ...
                'wp', dSing.sig_e.raw, ...
                'mu', dSing.mu.raw, ...
                'cut', cut1, ...
                'obs', obsSing, ...
                'CI', singData.relsummary.ciperc);
            testCase.verifyEqual(ptSing1, expSing1, 'AbsTol', 1e-12);
            
            close all force;
            prefs.view.criterioncutoff = cut2;
            psyrat_criterionfigures( ...
                'psyrat_prefs', prefs, ...
                'psyrat_data', singData, ...
                'analysis', 'sing');
            [ptSing2, singName2] = localCriterionPointEstimateFromOutputFigure();
            testCase.verifyTrue(contains(string(singName2), 'Dependability'));
            testCase.verifyNotEqual(ptSing1, ptSing2);
            
            trtData = PsyRATTestDataFactory.makeTRTSummaryData();
            prefs.view.gcoeff = 2; %criterion outputs should ignore this and use dependability
            
            close all force;
            prefs.view.criterioncutoff = cut1;
            psyrat_criterionfigures( ...
                'psyrat_prefs', prefs, ...
                'psyrat_data', trtData, ...
                'analysis', 'trt');
            [ptTrt1, trtName1] = localCriterionPointEstimateFromOutputFigure();
            testCase.verifyTrue(contains(string(trtName1), 'Dependability'));
            
            dTrt = trtData.relsummary.data.g(1).e(1);
            obsTrt = max(1, round(trtData.relsummary.group(1).event(1).trlinfo.mean));
            [~, expTrt1, ~] = psyrat_cutscore_trt( ...
                'gcoeff', 1, ...
                'reltype', prefs.view.reltype, ...
                'bp', dTrt.sig_id.raw, ...
                'bo', dTrt.sig_occ.raw, ...
                'bt', dTrt.sig_trl.raw, ...
                'txp', dTrt.sig_trlxid.raw, ...
                'oxp', dTrt.sig_occxid.raw, ...
                'txo', dTrt.sig_trlxocc.raw, ...
                'err', dTrt.sig_err.raw, ...
                'mu', dTrt.mu.raw, ...
                'cut', cut1, ...
                'obs', obsTrt, ...
                'nocc', 1, ...
                'CI', trtData.relsummary.ciperc);
            testCase.verifyEqual(ptTrt1, expTrt1, 'AbsTol', 1e-12);
            
            close all force;
            prefs.view.criterioncutoff = cut2;
            psyrat_criterionfigures( ...
                'psyrat_prefs', prefs, ...
                'psyrat_data', trtData, ...
                'analysis', 'trt');
            [ptTrt2, trtName2] = localCriterionPointEstimateFromOutputFigure();
            testCase.verifyTrue(contains(string(trtName2), 'Dependability'));
            testCase.verifyNotEqual(ptTrt1, ptTrt2);
        end

        function testCriterionScreensReportTrialCutoffFailure(testCase)
            %RC-33. The criterion screens rebuild the whole D-study themselves,
            %so they recompute trlcutoff and the eventgoodids/eventbadids
            %partition -- and under noccmode = 2 they do it at n'_o = 1 while the
            %parent view screen used the k-occasion composite. They previously
            %reported none of that: psyrat_criterionfigures returned only on
            %relerr.nogooddata and ignored relerr.trlcutoff, unlike
            %psyrat_relfigures.
            %
            %The fixture must be PRE-relsummary. makeTRTSummaryData ships a
            %prebuilt relsummary, which short-circuits the rebuild at the top of
            %psyrat_criterionfigures entirely, so relerr is never populated and
            %this test would pass vacuously against any implementation.
            prefs = localCriterionPrefs();
            trtData = localRawTrtRelData();

            %An unreachable threshold makes the cutoff search come up empty.
            %At noccmode = 2 the parent screen used a different estimand, so the
            %dialog must also say the retained sample can differ.
            close all force;
            prefs.view.relvalue = 0.9999;
            prefs.view.noccmode = 2;
            psyrat_criterionfigures('psyrat_prefs', prefs, ...
                'psyrat_data', trtData, 'analysis', 'trt');
            msg = localErrorDialogText();
            testCase.verifyNotEmpty(msg, ...
                'An unreachable threshold must raise the trial-cutoff dialog');
            testCase.verifyTrue(contains(msg, "could not be calculated"), ...
                'Dialog must mirror the psyrat_relfigures trial-cutoff message');
            %NOTE the quoting: inside a DOUBLE-quoted MATLAB string an apostrophe
            %is a literal, so the pattern is "n'_o", not "n''_o". The doubled
            %form searches for two apostrophes and never matches -- which would
            %make the negative assertions below pass for the wrong reason.
            testCase.verifyTrue(contains(msg, "n'_o = 1"), ...
                'At noccmode = 2 the dialog must name the divergence');
            testCase.verifyTrue(contains(msg, "retained participants"), ...
                'At noccmode = 2 the dialog must say the retained sample can differ');

            %Same failure at noccmode = 1: the parent screen is at n'_o = 1 too,
            %so there is nothing to diverge and the extra line would mislead.
            close all force;
            prefs.view.noccmode = 1;
            psyrat_criterionfigures('psyrat_prefs', prefs, ...
                'psyrat_data', trtData, 'analysis', 'trt');
            msg1 = localErrorDialogText();
            testCase.verifyTrue(contains(msg1, "could not be calculated"), ...
                'The cutoff failure is reported at either occasion mode');
            testCase.verifyFalse(contains(msg1, "retained participants"), ...
                'At noccmode = 1 there is no divergence to report');

            %And a reachable threshold must stay silent, so the branch is not
            %firing unconditionally.
            close all force;
            prefs.view.relvalue = 0.20;
            prefs.view.noccmode = 2;
            psyrat_criterionfigures('psyrat_prefs', prefs, ...
                'psyrat_data', trtData, 'analysis', 'trt');
            testCase.verifyEmpty(localErrorDialogText(), ...
                'A reachable threshold must raise no dialog');
            close all force;
        end

        function testCriterionRouteRefusesTheVirtualDiffRel(testCase)
            %G42, superseding this test's previous role. Before G42 this test
            %drove the DoD virtual-diff REL (rel.analysis = 'ic_diff') through
            %the criterion screens to exercise the RC-33 trlmax dialog via
            %'sing_diff''s derived-count producer. G42 closes that route ON
            %PURPOSE: the ic_diff criterion output came from
            %local_cutscore_sing's RC-05 sd_id/b_sigma stand-in, whose mu is
            %the posterior mean of a between-person SD -- not a location --
            %so the cut-score it produced is not an estimand. A dialog
            %reachable through a non-estimand does not justify the route.
            %
            %Scope of the closure, corrected by the adversarial wave
            %(2026-08-29): this closed ONE trlmax producer, not the last one.
            %The trlmax dialog remains live through this function when a
            %found cutoff excludes every participant -- the relsummary routes
            %then fall back to the bad-id roster for the trial table, so the
            %cutoff genuinely exceeds trlinfo.max (the trt producer is
            %pinned by testTrtExtrapolatedCutoffWithNoRetainedIdsStill-
            %WritesGoodN). The "tautologically dead per-event producers"
            %argument holds only while somebody is retained.
            pd = PsyRATTestDataFactory.makeDoDRelData();
            vd = psyrat_dod_build_virtual_rel('psyrat_data', pd, ...
                'mode', 'diff', 'groupidx', 1, 'eventpair', [1 4]);
            testCase.assertEqual(vd.rel.analysis, 'ic_diff', ...
                'weld: the virtual-diff REL must carry the label the guard keys on');

            prefs = localCriterionPrefs();
            prefs.view.diffgcoeff = 1;
            prefs.view.depvalue = 0.80;

            testCase.verifyError(@() psyrat_criterionfigures( ...
                'psyrat_prefs', prefs, 'psyrat_data', vd, ...
                'analysis', 'sing'), 'criterion:unsupportedDesign', ...
                ['the DoD virtual-diff REL must be refused by the G42 ' ...
                'guard rather than produce the RC-05 non-estimand cut-score']);
        end

        function testCriterionScreensStateTheOccasionEstimand(testCase)
            %RC-27/RC-29. The criterion surfaces pin n'_o = 1 by RULING, not by
            %validity: psyrat_cutscore_trt rejects gcoeff = 2 outright, but a
            %multi-occasion cut-score is a perfectly valid quantity, so pinning
            %it is an estimand choice. A choice has to be visible. Under
            %noccmode = 2 the parent view screen reports a k-occasion composite
            %while these screens report a single administration, and before this
            %the only place that said so was the saved file's header.
            %
            %The one-facet half is the more important assertion of the two: it
            %pins that a design with no occasion facet contributes no note, so
            %every one-facet title stays byte-identical.
            prefs = psyrat_defaults();
            prefs.view.plotcriterion = 1;
            prefs.view.tablecriterion = 1;
            prefs.view.ntrials = 6;
            prefs.view.reltype = 3;
            prefs.view.criterioncutoff = 0.5;
            prefs.view.depcentmeas = 1;
            prefs.view.relcentmeas = 1;

            close all force;
            psyrat_criterionfigures('psyrat_prefs', prefs, ...
                'psyrat_data', PsyRATTestDataFactory.makeTRTSummaryData(), ...
                'analysis', 'trt');
            [trtName, trtSub] = localCriterionTableTitles();
            %NOTE the quoting: inside a DOUBLE-quoted MATLAB string an
            %apostrophe is a literal character, so the pattern is "n'_o", not
            %"n''_o". Writing the doubled form here searches for two apostrophes
            %and never matches -- which makes the positive assertions below fail
            %loudly, but would make the negative ones further down pass for the
            %wrong reason.
            testCase.verifyTrue(contains(string(trtName), "n'_o = 1"), ...
                'Two-facet criterion table window title should state n''_o = 1');
            testCase.verifyTrue(contains(string(trtSub), "n'_o = 1"), ...
                'Two-facet criterion table should carry an on-screen occasion line');
            testCase.verifyTrue(contains(string(localCriterionPlotName()), "n'_o = 1"), ...
                'Two-facet criterion plot title should state n''_o = 1');

            close all force;
            psyrat_criterionfigures('psyrat_prefs', prefs, ...
                'psyrat_data', PsyRATTestDataFactory.makeICSummaryData(), ...
                'analysis', 'sing');
            [singName, singSub] = localCriterionTableTitles();
            testCase.verifyFalse(contains(string(singName), "n'_o"), ...
                'One-facet criterion title must not claim an occasion estimand');
            testCase.verifyEmpty(singSub, ...
                'One-facet criterion table must not gain an occasion line');
            testCase.verifyFalse(contains(string(localCriterionPlotName()), "n'_o"), ...
                'One-facet criterion plot title must not claim an occasion estimand');
            close all force;
        end

        function testCriterionFiguresSingIncludesTrialMainEffect(testCase)
            %Wiring-level regression: the one-facet criterion (cut-score)
            %dependability absolute error must include the trial main effect
            %sigma_i^2 (Rocha 2026 Table 2 denominator
            %sigma_p^2+(mu-C)^2+sigma_pi,e^2/n'+sigma_i^2/n'). This drives
            %psyrat_criterionfigures end-to-end so the local_cutscore_sing
            %wiring is exercised; the calculation-level tests call
            %psyrat_cutscore_sing directly (with i) and miss this class of bug.
            prefs = psyrat_defaults();
            prefs.view.plotcriterion = 0;
            prefs.view.tablecriterion = 1;
            prefs.view.depcentmeas = 1;
            prefs.view.gcoeff = 2; %must be ignored: criterion output is dependability
            cut = 0.1;
            prefs.view.criterioncutoff = cut;

            singData = PsyRATTestDataFactory.makeICSummaryDataWithTrial();
            d = singData.relsummary.data.g(1).e(1);
            obs = max(1, round(singData.relsummary.group(1).event(1).trlinfo.mean));

            close all force;
            psyrat_criterionfigures( ...
                'psyrat_prefs', prefs, ...
                'psyrat_data', singData, ...
                'analysis', 'sing');
            [ptActual, figName] = localCriterionPointEstimateFromOutputFigure();
            testCase.verifyTrue(contains(string(figName), 'Dependability'));

            %Correct cut-score: total within-person variance keeps sigma_i^2
            %(wp^2 = sig_e^2 + sig_trl^2), passing i = sig_trl.
            [~, ptCorrect, ~] = psyrat_cutscore_sing( ...
                'gcoeff', 1, ...
                'bp', d.sig_u.raw, ...
                'wp', sqrt(d.sig_e.raw.^2 + d.sig_trl.raw.^2), ...
                'i', d.sig_trl.raw, ...
                'mu', d.mu.raw, ...
                'cut', cut, ...
                'obs', obs, ...
                'CI', singData.relsummary.ciperc);

            %Value the pre-fix wiring produced: residual-only within-person
            %variance, omitting sigma_i^2.
            [~, ptBuggy, ~] = psyrat_cutscore_sing( ...
                'gcoeff', 1, ...
                'bp', d.sig_u.raw, ...
                'wp', d.sig_e.raw, ...
                'mu', d.mu.raw, ...
                'cut', cut, ...
                'obs', obs, ...
                'CI', singData.relsummary.ciperc);

            testCase.verifyEqual(ptActual, ptCorrect, 'AbsTol', 1e-12, ...
                'Criterion cut-score must include the trial main effect sigma_i^2.');
            testCase.verifyGreaterThan(abs(ptCorrect - ptBuggy), 1e-3, ...
                'Fixture must make the sigma_i^2 contribution observable.');
            testCase.verifyNotEqual(ptActual, ptBuggy);
        end

        function testCriterionPathForcesAbsoluteCutScoreRegardlessOfGcoeff(testCase)
            %F7 wiring regression: cut-scores are absolute-error only. Even
            %when the user selects a relative/generalizability decision
            %(prefs.view.gcoeff = 2), psyrat_criterionfigures must force
            %gcoeff = 1 so the criterion output is the absolute (dependability)
            %cut-score for both the one-facet and two-facet designs. Because the
            %cut-score functions now REJECT gcoeff = 2, the only way the
            %criterion path can succeed with gcoeff = 2 selected is by
            %hardcoding the absolute decision -- so this exercises the wiring,
            %not only the functions (the F1 lesson).
            prefs = psyrat_defaults();
            prefs.view.plotcriterion = 0;
            prefs.view.tablecriterion = 1;
            prefs.view.ntrials = 6;
            prefs.view.depcentmeas = 1;
            prefs.view.relcentmeas = 1;
            prefs.view.reltype = 3;
            prefs.view.gcoeff = 2; %relative selection that must be ignored
            cut = 0.3;
            prefs.view.criterioncutoff = cut;

            %--- one-facet criterion path ---
            singData = PsyRATTestDataFactory.makeICSummaryData();
            close all force;
            psyrat_criterionfigures( ...
                'psyrat_prefs', prefs, 'psyrat_data', singData, 'analysis', 'sing');
            [ptSing, singName] = localCriterionPointEstimateFromOutputFigure();
            testCase.verifyTrue(contains(string(singName), 'Dependability'));

            dSing = singData.relsummary.data.g(1).e(1);
            obsSing = max(1, round(singData.relsummary.group(1).event(1).trlinfo.mean));
            [~, expSing, ~] = psyrat_cutscore_sing('gcoeff', 1, ...
                'bp', dSing.sig_u.raw, 'wp', dSing.sig_e.raw, 'mu', dSing.mu.raw, ...
                'cut', cut, 'obs', obsSing, 'CI', singData.relsummary.ciperc);
            testCase.verifyEqual(ptSing, expSing, 'AbsTol', 1e-12, ...
                'One-facet criterion output must be the absolute (dependability) cut-score.');

            %--- two-facet (test-retest) criterion path ---
            trtData = PsyRATTestDataFactory.makeTRTSummaryData();
            close all force;
            psyrat_criterionfigures( ...
                'psyrat_prefs', prefs, 'psyrat_data', trtData, 'analysis', 'trt');
            [ptTrt, trtName] = localCriterionPointEstimateFromOutputFigure();
            testCase.verifyTrue(contains(string(trtName), 'Dependability'));

            dTrt = trtData.relsummary.data.g(1).e(1);
            obsTrt = max(1, round(trtData.relsummary.group(1).event(1).trlinfo.mean));
            [~, expTrt, ~] = psyrat_cutscore_trt('gcoeff', 1, 'reltype', prefs.view.reltype, ...
                'bp', dTrt.sig_id.raw, 'bo', dTrt.sig_occ.raw, 'bt', dTrt.sig_trl.raw, ...
                'txp', dTrt.sig_trlxid.raw, 'oxp', dTrt.sig_occxid.raw, ...
                'txo', dTrt.sig_trlxocc.raw, 'err', dTrt.sig_err.raw, 'mu', dTrt.mu.raw, ...
                'cut', cut, 'obs', obsTrt, 'nocc', 1, 'CI', trtData.relsummary.ciperc);
            testCase.verifyEqual(ptTrt, expTrt, 'AbsTol', 1e-12, ...
                'Two-facet criterion output must be the absolute (dependability) cut-score.');

            %The criterion path could only succeed above because it forced
            %gcoeff = 1: the cut-score functions reject the relative decision.
            testCase.verifyError(@() psyrat_cutscore_sing('gcoeff', 2, ...
                'bp', dSing.sig_u.raw, 'wp', dSing.sig_e.raw, 'mu', dSing.mu.raw, ...
                'cut', cut, 'obs', obsSing, 'CI', singData.relsummary.ciperc), ...
                'varargin:gcoeff');
            testCase.verifyError(@() psyrat_cutscore_trt('gcoeff', 2, 'reltype', prefs.view.reltype, ...
                'bp', dTrt.sig_id.raw, 'bo', dTrt.sig_occ.raw, 'bt', dTrt.sig_trl.raw, ...
                'txp', dTrt.sig_trlxid.raw, 'oxp', dTrt.sig_occxid.raw, ...
                'txo', dTrt.sig_trlxocc.raw, 'err', dTrt.sig_err.raw, 'mu', dTrt.mu.raw, ...
                'cut', cut, 'obs', obsTrt, 'nocc', 1, 'CI', trtData.relsummary.ciperc), ...
                'varargin:gcoeff');
        end

        function testFormulaCatalogAvailabilitySummary(testCase)
            catalog = psyrat_formula_catalog;

            requiredCols = {'formula_id','paper_table','design','coefficient','decision_type',...
                'formula_ascii','psyrat_function','gui_status','implementation_status','notes'};
            testCase.verifyTrue(all(ismember(requiredCols, catalog.Properties.VariableNames)));
            testCase.verifyGreaterThan(height(catalog), 0);
            testCase.verifyTrue(all(~strcmp(catalog.implementation_status,'Needs implementation')));

            singCatalog = psyrat_formula_catalog('analysis','sing');
            testCase.verifyGreaterThan(height(singCatalog), 0);
            testCase.verifyTrue(all(ismember(singCatalog.paper_table, {'2','6'})));
            
            singDiffCatalog = psyrat_formula_catalog('analysis','sing_diff');
            testCase.verifyFalse(any(ismember(singDiffCatalog.formula_id, {'t6_diff2_d','t6_diff2_g'})));
            
            trtCatalog = psyrat_formula_catalog('analysis','trt');
            testCase.verifyFalse(any(ismember(trtCatalog.formula_id, {'t6_diff1_d','t6_diff1_g'})));

            % dynamic/conditional reliability (Rast & Clayson) rows are tagged
            % paper_table 'R&C'; the dynrel filter returns them and the Rocha
            % sing/trt filters must exclude them.
            dynrelCatalog = psyrat_formula_catalog('analysis','dynrel');
            testCase.verifyGreaterThan(height(dynrelCatalog), 0);
            testCase.verifyTrue(all(strcmp(dynrelCatalog.paper_table, 'R&C')));
            testCase.verifyTrue(any(strcmp(dynrelCatalog.formula_id, 'dyn_g')));
            testCase.verifyFalse(any(strcmp(singCatalog.paper_table, 'R&C')));
            testCase.verifyFalse(any(strcmp(trtCatalog.paper_table, 'R&C')));

            % nonparallel splits (Rocha Tables 4/5) rows are tagged paper_table
            % '4'/'5'; the splits filter returns them and the Rocha sing/trt
            % filters must exclude them.
            splitsCatalog = psyrat_formula_catalog('analysis','splits');
            testCase.verifyGreaterThan(height(splitsCatalog), 0);
            testCase.verifyTrue(all(ismember(splitsCatalog.paper_table, {'4','5'})));
            testCase.verifyTrue(any(strcmp(splitsCatalog.formula_id, 't4_d')));
            testCase.verifyTrue(any(strcmp(splitsCatalog.formula_id, 't5_d')));
            testCase.verifyFalse(any(ismember(singCatalog.paper_table, {'4','5'})));
            testCase.verifyFalse(any(ismember(trtCatalog.paper_table, {'4','5'})));

            % RC-14: a formula string must determine its own decision type. Two
            % rows may legitimately share a formula_ascii when they describe the
            % same decision in different designs -- t3_rel_d/t5_d and
            % t3_rel_g/t5_g do, Rocha Table 3 two-facet vs Table 5 splits
            % p x (i:s) x o -- so plain uniqueness is the wrong invariant and
            % would fail on correct rows. What can never be right is one string
            % serving BOTH an Absolute and a Relative row, which is the defect
            % t3_icc_abs/t3_icc_rel carried: byte-identical strings evaluating
            % to the absolute quantity under opposite decision_type labels.
            [~, ~, grp] = unique(catalog.formula_ascii);
            for g = 1:max(grp)
                rows = find(grp == g);
                dtypes = unique(catalog.decision_type(rows));
                testCase.verifyEqual(numel(dtypes), 1, sprintf( ...
                    ['formula_ascii "%s" is shared by rows with different ' ...
                    'decision_type (%s); a formula string cannot be correct ' ...
                    'for both. Offending formula_id: %s'], ...
                    catalog.formula_ascii{rows(1)}, ...
                    strjoin(dtypes(:)', ', '), ...
                    strjoin(catalog.formula_id(rows)', ', ')));
            end

            for i = 1:height(catalog)
                [tmpl, ttl] = psyrat_formula_template('formula_id', catalog.formula_id{i});
                testCase.verifyNotEmpty(tmpl);
                testCase.verifyNotEmpty(ttl);
            end
        end

        function testStartviewSingShowsCellEventLabelForSingleEvent(testCase)
            prefs = psyrat_defaults();
            data = PsyRATTestDataFactory.makeSingleEventViewData();

            close all force;
            psyrat_startview_sing('psyrat_prefs', prefs, 'psyrat_data', data);

            guiFig = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui');
            testCase.verifyFalse(isempty(guiFig));

            textControls = findall(guiFig(1), 'Type', 'uicontrol', 'Style', 'text');
            labels = string(get(textControls, 'String'));
            testCase.verifyTrue(any(contains(labels, "Cell/Event:  Reward")));
        end

        function testStartviewSingDoesNotShowCellEventLabelForDiffMode(testCase)
            prefs = psyrat_defaults();
            data = PsyRATTestDataFactory.makeSingleEventViewData();
            data.rel.analysis = 'ic_diff';
            data.rel.events = {'Reward'; 'Loss'};

            close all force;
            psyrat_startview_sing('psyrat_prefs', prefs, 'psyrat_data', data);

            guiFig = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui');
            testCase.verifyFalse(isempty(guiFig));

            textControls = findall(guiFig(1), 'Type', 'uicontrol', 'Style', 'text');
            labels = string(get(textControls, 'String'));
            testCase.verifyFalse(any(contains(labels, "Cell/Event:")));
        end

        function testSingDiffPerEventCoefficientsIncludeTrialMainEffect(testCase)
            %RC-02. The six per-event psyrat_rel_sing calls in the 'sing_diff'
            %branch passed the person x trial,e residual as the within-person SD
            %and omitted the trial main effect entirely, so they returned the
            %RELATIVE coefficient no matter which estimand the user picked -- and
            %the output labeled it "Dependability". Because the two coefficients
            %were bit-identical, nothing in the numbers could reveal it.
            %
            %The route is built from the DoD fixture's difference mode, which
            %emits a genuine ic_diff structure (b/b_sigma/sd_id/sd_trl/...), so
            %this exercises the real summary path rather than a hand-assembled
            %struct.
            pd = PsyRATTestDataFactory.makeDoDRelData();
            vd = psyrat_dod_build_virtual_rel('psyrat_data', pd, ...
                'mode', 'diff', 'groupidx', 1, 'eventpair', [1 4]);

            dep = localSingDiff(vd, 1);
            gen = localSingDiff(vd, 2);

            for e = 1:2
                d = dep.relsummary.group(1).event(e);
                g = gen.relsummary.group(1).event(e);

                %Absolute error carries sigma_i, relative error does not, so
                %dependability must sit strictly below generalizability. Equality
                %is the exact signature of the defect.
                testCase.verifyLessThan(d.dep.m, g.dep.m, sprintf( ...
                    'event %d: dependability must be below generalizability', e));
                testCase.verifyGreaterThan(g.dep.m - d.dep.m, 1e-6, sprintf( ...
                    'event %d: coefficients are effectively identical, so sigma_i is not reaching the calc', e));

                %The ICC surface takes the same components and had the same defect.
                testCase.verifyLessThan(d.icc.m, g.icc.m, sprintf( ...
                    'event %d: absolute ICC must be below relative ICC', e));
            end

            %The consequence that reaches a user decision: the per-event curve
            %drives the trial cutoff, and the cutoff drives participant
            %inclusion. An optimistic coefficient retains participants a correct
            %one would drop, so the two estimands must not agree here either.
            testCase.verifyGreaterThanOrEqual( ...
                dep.relsummary.group(1).event(1).trlcutoff, ...
                gen.relsummary.group(1).event(1).trlcutoff, ...
                'The absolute-error estimand must require at least as many trials.');
        end

        function testSingDiffPerEventFollowsDiffGcoeffNotStaleGcoeff(testCase)
            %RC-02b. The ic_diff viewer builds only a "Difference-Score
            %Coefficient" popup; it never builds a G-Theory Coefficient popup. So
            %the population gcoeff reaching the per-event sites came from a stored
            %preference the user cannot see on that screen, and it can carry over
            %from a previously viewed one-facet run. Passing the two in conflict
            %pins the resolution: the visible control wins.
            pd = PsyRATTestDataFactory.makeDoDRelData();
            vd = psyrat_dod_build_virtual_rel('psyrat_data', pd, ...
                'mode', 'diff', 'groupidx', 1, 'eventpair', [1 4]);

            %Stale pref says generalizability; the user's visible choice says
            %dependability. The output must follow the user.
            conflict = localSingDiff(vd, 1, 2);
            testCase.verifyEqual(conflict.relsummary.gcoeff, 1);
            testCase.verifyEqual(conflict.relsummary.gcoeff_name, 'dep');

            %And it must agree numerically with the run where nothing conflicts,
            %not with the run matching the stale pref.
            agree = localSingDiff(vd, 1);
            other = localSingDiff(vd, 2);
            testCase.verifyEqual( ...
                conflict.relsummary.group(1).event(1).dep.m, ...
                agree.relsummary.group(1).event(1).dep.m, 'AbsTol', 1e-12);
            testCase.verifyNotEqual( ...
                conflict.relsummary.group(1).event(1).dep.m, ...
                other.relsummary.group(1).event(1).dep.m);

            %Mirror case, so the test cannot pass by hardcoding dependability.
            conflict2 = localSingDiff(vd, 2, 1);
            testCase.verifyEqual(conflict2.relsummary.gcoeff, 2);
            testCase.verifyEqual(conflict2.relsummary.gcoeff_name, 'gen');
        end

        function testDoDObservedTableAndVarComponentExport(testCase)
            data = PsyRATTestDataFactory.makeDoDRelData();

            observed = psyrat_dod_observedt( ...
                'psyrat_data', data, ...
                'groupidx', 1, ...
                'map', [1 2 3 4], ...
                'gui', 0);

            testCase.verifyEqual(observed.Properties.VariableNames, ...
                {'Label', 'n_Included', 'n_Excluded', 'Dependability', ...
                'Generalizability', 'ICC_Dep', 'ICC_Gen', ...
                'MeanTrials_ERP_1', 'MeanTrials_ERP_2', ...
                'MeanTrials_ERP_3', 'MeanTrials_ERP_4'});
            testCase.verifyEqual(height(observed), 1);
            testCase.verifyEqual(observed.n_Included(1), 3);
            testCase.verifyEqual(observed.n_Excluded(1), 0);
            testCase.verifyEqual(observed.MeanTrials_ERP_1(1), 31/3, 'AbsTol', 1e-12);
            testCase.verifyEqual(observed.MeanTrials_ERP_2(1), 38/3, 'AbsTol', 1e-12);
            testCase.verifyEqual(observed.MeanTrials_ERP_3(1), 28/3, 'AbsTol', 1e-12);
            testCase.verifyEqual(observed.MeanTrials_ERP_4(1), 35/3, 'AbsTol', 1e-12);

            varcomp = psyrat_dod_varcompt( ...
                'psyrat_data', data, ...
                'groupidx', 1, ...
                'map', [1 2 3 4], ...
                'savefile', 0, ...
                'gui', 0);

            testCase.verifyEqual(varcomp.Properties.VariableNames, ...
                {'Component', 'Point', 'Low_CrI', 'High_CrI'});
            testCase.verifyEqual(height(varcomp), 9);
            testCase.verifyEqual(varcomp.Component{1}, 'Universe variance of DoD contrast');
            testCase.verifyEqual(varcomp.Component{6}, 'Dependability');
            testCase.verifyEqual(varcomp.Component{7}, 'Generalizability');
        end

        % --- engine provenance notice (viewer labeling) --------------------
        function testEngineNoticeCmdstanAndLegacyAreSilent(testCase)
            % CmdStan runs and legacy files (no engine field) produce no notice,
            % so the default workflow is visually unchanged.
            testCase.verifyEqual(psyrat_engine_notice(struct('engine','cmdstan')), '');
            testCase.verifyEqual(psyrat_engine_notice(struct('engine','')), '');
            testCase.verifyEqual(psyrat_engine_notice(struct('analysis','ic')), '');
            testCase.verifyEqual(psyrat_engine_notice([]), '');
        end

        function testEngineNoticeHmcFlagsNonBitIdentical(testCase)
            msg = psyrat_engine_notice(struct('engine','hmc'));
            testCase.verifyNotEmpty(msg);
            testCase.verifyTrue(contains(msg, 'HMC'));
            testCase.verifyTrue(contains(msg, 'not bit-identical', 'IgnoreCase', true) || ...
                contains(msg, 'NOT bit-identical', 'IgnoreCase', true));
        end

        function testEngineNoticeFitlmeFlagsApproximate(testCase)
            msg = psyrat_engine_notice(struct('engine','fitlme'));
            testCase.verifyNotEmpty(msg);
            testCase.verifyTrue(contains(msg, 'fitlme'));
            testCase.verifyTrue(contains(msg, 'APPROXIMATE') && contains(msg, 'frequentist'));
        end

        function testEngineNoticeUnknownEngineStillFlags(testCase)
            msg = psyrat_engine_notice(struct('engine','someother'));
            testCase.verifyNotEmpty(msg);
            testCase.verifyTrue(contains(msg, 'someother'));
        end

        function testTrtExtrapolatedCutoffWithNoRetainedIdsStillWritesGoodN(testCase)
            %RC-38 Cause A. In the test-retest group branch (psyrat_relsummary
            %case 'trt', analysis 2) goodids/badids are a CHAR 'none' when no
            %trial cutoff was found, but an EMPTY CELL when a cutoff WAS found
            %and then excluded every participant -- an extrapolated cutoff above
            %everyone's trial count. strcmp on an empty cell returns [], so the
            %two "elseif strcmp(...)" guards were non-exhaustive: neither arm
            %fired, leaving 'goodids' unassigned (a hard crash on the isempty
            %check) and 'goodn' unwritten.
            %
            %The fixture makes the cutoff land at ~14 trials against 5 observed,
            %which is the only route to the empty-cell state. Asserting goodn is
            %what separates fixing BOTH guards from fixing only the first: with
            %just the first repaired the function stops throwing but silently
            %omits the retained-participant count.
            psyrat_data = localRawTrtRelDataMultiGroup();

            [out, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'trt', ...
                'gcoeff', 1, ...
                'reltype', 1, ...
                'relcutoff', 0.70, ...
                'meascutoff', 2, ...
                'relcentmeas', 1, ...
                'CI', 0.95);

            for gloc = 1:2
                ev = out.relsummary.group(gloc).event(1);

                %the empty-cell state is genuinely reached (guards the fixture
                %itself: if a later change made the cutoff attainable this test
                %would otherwise pass without exercising the defect at all)
                testCase.verifyTrue(iscell(ev.eventgoodids), ...
                    sprintf(['Group %d must reach the EMPTY CELL path, not ' ...
                    'the char ''none'' path.'], gloc));
                testCase.verifyEmpty(ev.eventgoodids, ...
                    sprintf('Group %d should retain no participants.', gloc));
                testCase.verifyGreaterThan(ev.trlcutoff, 5, ...
                    sprintf(['Group %d cutoff must be EXTRAPOLATED beyond the ' ...
                    '5 observed trials.'], gloc));

                %the second guard: goodn must exist AND be zero
                testCase.verifyTrue(isfield(ev, 'goodn'), ...
                    sprintf('Group %d must have goodn written.', gloc));
                testCase.verifyEqual(ev.goodn, 0, ...
                    sprintf(['Group %d retained no participants, so goodn ' ...
                    'must be 0.'], gloc));
            end

            %the extrapolation is reported through relerr rather than silently
            testCase.verifyEqual(relerr.trlmax, 1);
        end

        function testSubjectLevelDiffWithNoRetainedIdsRaisesTheNotice(testCase)
            %RC-40, superseded on this route by G48 (owner ruling
            %2026-08-29). History: the sing_diff_sserr helper used to set
            %relerr.nogooddata silently (a user who asked for results got an
            %empty screen with no reason given); RC-40 added the one-per-run
            %modal; G48 then removed the flag-and-modal entirely, aligning
            %the route with the G39/G45 strict-fold conventions. The notice
            %is now a per-group CONSOLE line, the flag stays 0, and the
            %summary returns in full -- so the modal threshold dialog must
            %NOT appear even on an interactive run. The full alignment pins
            %(complementing tables, relfigures drawing output, one group
            %not blanking another) live in TestSingDiffSserrStrictAlignment.
            %
            %The RC-40 reachability record for the TEN guarded sites (six
            %live, four dead, adjudicated 2026-08-18 by profiler) is
            %unaffected by G48 and lives above the helper in the source.
            psyrat_data = PsyRATTestDataFactory.makeICDiffSSErrRelData();
            psyrat_data.proc.interactive = true;

            dlgname = 'Data do not meet reliability threshold';
            delete(findall(0, 'Type', 'figure', 'Name', dlgname));

            txt = evalc(['[out, relerr] = psyrat_relsummary(' ...
                '''psyrat_data'', psyrat_data, ' ...
                '''analysis'', ''sing_diff_sserr'', ' ...
                '''depcutoff'', 0.9999, ''meascutoff'', 2, ' ...
                '''depcentmeas'', 1, ''diffgcoeff'', 1, ''CI'', 0.95);']);

            %guards the fixture itself: if a later change made this cutoff
            %attainable the test would pass without exercising the state
            testCase.assertTrue( ...
                isempty(out.relsummary.group(1).goodids), ...
                'The fixture must actually empty goodids, or this proves nothing.');

            testCase.verifyEqual(relerr.nogooddata, 0, ...
                'G48: an empty inclusion must not set the abort flag');
            testCase.verifyTrue(contains(txt, 'difference score'), ...
                'a run that retains nobody must still say so, on the console');
            testCase.verifyEmpty( ...
                findall(0, 'Type', 'figure', 'Name', dlgname), ...
                'G48: the modal threshold dialog must not appear on this route');
        end

        function testHeadlessEmptyInclusionStaysDialogFree(testCase)
            %The headless half of the same contract, kept because a batch
            %caller must never block on a modal. Post-G48 there is no dialog
            %on this route in EITHER mode; the console line carries the
            %notice and the summary is returned in full, so a batch caller
            %reads the zeroed group straight from the output instead of
            %detecting an aborted run through relerr.
            %
            %Honesty note (adversarial wave, 2026-08-29): renamed from
            %testTheNogooddataNoticeStaysOffHeadlessRuns -- the notice that
            %name referred to no longer exists on this route. The no-dialog
            %assert could not fail pre-G48 either (the modal was already
            %gated on interactive); the discriminating asserts are
            %nogooddata == 0 and the console-notice text, and the test's
            %residual value is pinning that no dialog regresses back in.
            psyrat_data = PsyRATTestDataFactory.makeICDiffSSErrRelData();
            psyrat_data.proc.interactive = false;

            dlgname = 'Data do not meet reliability threshold';
            delete(findall(0, 'Type', 'figure', 'Name', dlgname));

            txt = evalc(['[out, relerr] = psyrat_relsummary(' ...
                '''psyrat_data'', psyrat_data, ' ...
                '''analysis'', ''sing_diff_sserr'', ' ...
                '''depcutoff'', 0.9999, ''meascutoff'', 2, ' ...
                '''depcentmeas'', 1, ''diffgcoeff'', 1, ''CI'', 0.95);']);

            testCase.assertTrue( ...
                isempty(out.relsummary.group(1).goodids), ...
                'The fixture must actually empty goodids, or this proves nothing.');
            testCase.verifyEqual(relerr.nogooddata, 0);
            testCase.verifyTrue(contains(txt, 'difference score'), ...
                'the console notice must fire on headless runs too');
            testCase.verifyEmpty( ...
                findall(0, 'Type', 'figure', 'Name', dlgname), ...
                'A headless run must not raise a dialog.');
        end

    end
end

function psyrat_data = localRawTrtRelDataMultiGroup()
%RC-38 helper. Pre-relsummary test-retest fixture with TWO groups and no
%events, which is the only shape that routes into psyrat_relsummary's
%case 'trt' / analysis 2 branch (ngroups > 1 && nevents == 1).
%
%Three participants per group x two occasions x five trials. The variance
%components are deliberately noisy (sigma_id 2.0 against sigma_err 4.9) so the
%dependability curve crosses 0.70 only at ~14 trials. That matters: the cutoff
%must be FOUND (otherwise the code takes the char 'none' path) but must exceed
%every participant's five observed trials, which is what drives the per-id
%sweep to exclude everyone and leaves eventgoodids as an EMPTY CELL.
subs = {'s1','s2','s3'};
occ = {'t1','t2'};
gnames = {'GroupA'; 'GroupB'};
ntrl = 5;

ids = {}; times = {}; grps = {}; meas = [];
for g = 1:numel(gnames)
    for s = 1:numel(subs)
        for o = 1:numel(occ)
            for t = 1:ntrl
                %ids are namespaced per group so the two groups do not pool
                ids{end+1,1} = sprintf('%s_%s', gnames{g}, subs{s}); %#ok<AGROW>
                times{end+1,1} = occ{o}; %#ok<AGROW>
                grps{end+1,1} = gnames{g}; %#ok<AGROW>
                meas(end+1,1) = numel(meas) + 1; %#ok<AGROW>
            end
        end
    end
end

ndraw = 5;
ng = numel(gnames);

rel = struct();
rel.analysis = 'trt';
rel.groups = gnames;
rel.events = 'none';
rel.filename = 'rc38_trt_multigroup_fixture.psyrat';
rel.nchains = 2;
rel.niter = 100;
rel.time = occ(:);
rel.data = table(ids, grps, meas, times, ...
    'VariableNames', {'id','group','meas','time'});

%analysis 2 reads rel.out.* as one COLUMN per group
rel.out = struct();
rel.out.labels      = gnames;
rel.out.mu          = repmat(0.10, ndraw, ng);
rel.out.sig_id      = repmat(2.00, ndraw, ng);
rel.out.sig_occ     = repmat(0.30, ndraw, ng);
rel.out.sig_trl     = repmat(0.20, ndraw, ng);
rel.out.sig_trlxid  = repmat(0.20, ndraw, ng);
rel.out.sig_occxid  = repmat(0.40, ndraw, ng);
rel.out.sig_trlxocc = repmat(0.09, ndraw, ng);
rel.out.sig_err     = repmat(4.90, ndraw, ng);

psyrat_data = struct('rel', rel);
%headless: suppress the "does not reach threshold" modal
psyrat_data.proc.interactive = false;
end

function val = localPointEstimate(cellstr_entry)
%Pull the point estimate out of a rendered ' %0.2f CI [%0.2f %0.2f]' cell, so a
%cross-column invariant can be checked against what the table actually PRINTS
%rather than against a re-implementation of how it was computed.
val = sscanf(strtrim(cellstr_entry), '%f', 1);
if isempty(val)
    error('tests:pointestimate', ...
        'Could not parse a point estimate from "%s".', cellstr_entry);
end
end

function out = localSingDiff(psyrat_data, diffgcoeff, gcoeff)
%Summarize a single-session difference-score result. diffgcoeff is the user's
%visible "Difference-Score Coefficient" choice; gcoeff is the population
%preference, passed only when a test needs the two to disagree so the resolution
%is pinned rather than assumed.
args = {'psyrat_data', psyrat_data, ...
    'analysis', 'sing_diff', ...
    'depcutoff', 0.5, ...
    'meascutoff', 2, ...
    'depcentmeas', 1, ...
    'diffgcoeff', diffgcoeff, ...
    'CI', 0.95};
if nargin >= 3
    args = [args, {'gcoeff', gcoeff}];
end
[out, relerr] = psyrat_relsummary(args{:});
if relerr.nogooddata ~= 0
    error('tests:singdiff', ...
        'sing_diff fixture produced no usable data (nogooddata = %d).', ...
        relerr.nogooddata);
end
end

function [ptEstimate, figName] = localCriterionPointEstimateFromOutputFigure()
figs = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output');

for i = 1:numel(figs)
    nm = figs(i).Name;
    if contains(string(nm), 'Overall Criterion')
        figName = nm;
        t = findall(figs(i), 'Type', 'uitable');
        if isempty(t)
            error('tests:criteriontable', ...
                'Criterion summary figure did not contain a uitable.');
        end
        data = t(1).Data;
        if isempty(data)
            error('tests:criteriontable', ...
                'Criterion summary table was empty.');
        end
        %Locate the point-estimate column by NAME, not by position. This used to
        %read data{1,5}, which silently encoded the column ORDER of
        %local_tablecriterion into a helper three tests share -- so RC-36's
        %n_Included column (inserted at position 2 to match the parent overall
        %tables) shifted the estimate to column 6 and failed all three at once
        %with no hint that a column had moved. Re-hardcoding 6 would just reset
        %the same trap for the next column.
        cn = cellstr(t(1).ColumnName);
        col = find(contains(cn, 'Point Estimate'), 1);
        if isempty(col)
            error('tests:criteriontable', ...
                'Criterion summary table has no point-estimate column.');
        end
        ptEstimate = data{1,col};
        return;
    end
end

error('tests:criteriontable', 'Could not find overall criterion summary figure.');
end

function msg = localErrorDialogText()
%RC-33 helper. Return the text of the errordlg the criterion screens raise, or
%'' when there is none. errordlg builds a msgbox figure tagged
%'Msgbox_<title>'; the message lines live in its text object as a cellstr, so
%they are joined with newlines and matched with contains(). findall is required
%because the dialog's handle visibility is off.
%
%The empty return is what the negative assertions rely on, so it must mean "no
%dialog" and not "dialog I failed to locate" -- hence keying off the tag rather
%than off the message text.
%
%The 2026-08-17 notification ruling gave the two notices distinct purpose
%titles (still non-modal, still raised last), and errordlg tags follow the
%title, so the sweep covers both purpose tags plus the legacy untitled tag --
%an empty return still means "no notice dialog at all".
msg = '';
tags = {'Msgbox_Cutoff not calculable', ...
    'Msgbox_Extrapolation beyond data', ...
    'Msgbox_Error Dialog'};
dlgs = [];
for k = 1:numel(tags)
    dlgs = [dlgs; findall(groot, 'Type', 'figure', 'Tag', tags{k})]; %#ok<AGROW>
end
if isempty(dlgs)
    return;
end
lines = {};
for k = 1:numel(dlgs)
    txt = findall(dlgs(k), 'Type', 'text');
    for j = 1:numel(txt)
        lines = [lines; cellstr(txt(j).String)]; %#ok<AGROW>
    end
end
msg = strjoin(lines, newline);
end

function prefs = localCriterionPrefs()
%RC-33 helper. Shared criterion-screen preferences for the dialog tests. Both
%output surfaces stay on so the dialogs are raised after real figures exist,
%which is the ordering psyrat_criterionfigures relies on.
prefs = psyrat_defaults();
prefs.view.plotcriterion = 1;
prefs.view.tablecriterion = 1;
prefs.view.ntrials = 6;
prefs.view.reltype = 3;
prefs.view.criterioncutoff = 0.5;
prefs.view.depcentmeas = 1;
prefs.view.relcentmeas = 1;
prefs.view.meascutoff = 2;
end

function psyrat_data = localRawTrtRelData()
%RC-33 helper. Minimal PRE-relsummary test-retest struct (rel.out draws plus a
%rel.data trial table) -- the input psyrat_relsummary consumes, and therefore
%the only kind of fixture that makes psyrat_criterionfigures actually rebuild.
%Three participants x two occasions x five trials, so an unreachable threshold
%reports the cutoff failure through relerr rather than aborting.
%
%The stated reason used to be "leaves the cutoff search empty without emptying
%the bad-id list, which would trip the nogooddata guard" -- measured false on
%2026-08-18. This fixture FIRES no nogooddata guard: an unreachable threshold
%takes the trlcutoff = -1 arm, and since G54 respelled the second guard on the
%sentinel value, that arm (everyone bad) now EVALUATES its isempty(badids)
%guard -- but badids holds the full three-subject roster there, so the guard
%cannot fire. The pre-G54 record (all ten guarded sites at profiler hits = 0
%across this test's three calls) described the era when the second
%"if isempty(trlcutoff)" could never be entered at all; the setter bodies
%remain unexecuted either way.
subs = {'s1','s2','s3'};
occ = {'t1','t2'};
ntrl = 5;
ids = {}; times = {}; meas = [];
for s = 1:numel(subs)
    for o = 1:numel(occ)
        for t = 1:ntrl
            ids{end+1,1} = subs{s}; %#ok<AGROW>
            times{end+1,1} = occ{o}; %#ok<AGROW>
            meas(end+1,1) = numel(meas) + 1; %#ok<AGROW>
        end
    end
end

rel = struct();
rel.analysis = 'trt';
rel.groups = 'none';
rel.events = 'none';
rel.filename = 'rc33_criterion_fixture.psyrat';
rel.nchains = 2;
rel.niter = 100;
rel.time = occ(:);
rel.data = table(ids, meas, times, 'VariableNames', {'id','meas','time'});

rel.out = struct();
rel.out.labels = {'measure'};
rel.out.mu = [0.10; 0.12; 0.11; 0.09; 0.10];
rel.out.sig_id      = [2.00; 2.10; 1.90; 2.05; 1.95];
rel.out.sig_occ     = [0.30; 0.35; 0.28; 0.32; 0.31];
rel.out.sig_trl     = [0.20; 0.25; 0.23; 0.21; 0.22];
rel.out.sig_trlxid  = [0.20; 0.22; 0.18; 0.21; 0.19];
rel.out.sig_occxid  = [0.40; 0.42; 0.38; 0.41; 0.39];
rel.out.sig_trlxocc = [0.09; 0.10; 0.11; 0.09; 0.10];
rel.out.sig_err     = [4.90; 5.00; 5.10; 4.80; 5.20];

psyrat_data = struct('rel', rel);
end

function [figName, subtitle] = localCriterionTableTitles()
%RC-27/RC-29 helper. Return the criterion SUMMARY TABLE figure's window title
%and its on-screen occasion line. The occasion line is a separate text
%uicontrol rather than an addition to the heading, because the heading is
%already ~47 characters at fontsize 16 in a 700 px control; it is identified by
%its leading text, and returns '' when the figure has none (the one-facet case,
%which is the point of the negative assertion).
figs = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
for i = 1:numel(figs)
    if contains(string(figs(i).Name), 'Overall Criterion')
        figName = figs(i).Name;
        subtitle = '';
        txt = findall(figs(i), 'Style', 'text');
        for k = 1:numel(txt)
            s = txt(k).String;
            if iscell(s); s = strjoin(s, ' '); end
            if contains(string(s), 'Cut-scores computed at')
                subtitle = char(s);
                return;
            end
        end
        return;
    end
end
error('tests:criteriontable', 'Could not find overall criterion summary figure.');
end

function figName = localCriterionPlotName()
%RC-27/RC-29 helper. Return the criterion CURVE figure's window title. The two
%criterion figures share the 'psyrat_output' tag, so they are told apart by
%name: the table's begins 'Overall Criterion', the plot's is built from
%L.NumberOf and so contains ' v Number of'.
figs = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
for i = 1:numel(figs)
    if contains(string(figs(i).Name), ' v Number of')
        figName = figs(i).Name;
        return;
    end
end
error('tests:criterionplot', 'Could not find criterion curve figure.');
end

function lns = localFindCurveLines(fig, expectedX)
%RC-30 helper. A reliability-vs-trials figure holds two kinds of line: the
%coefficient curves, one per group with one point per requested trial count,
%and the cutoff reference line, which psyrat_addhline draws with exactly two
%points (psyrat_addhline.m:34). Select the curves by matching XData to the
%requested range, so a buffer of the wrong length or offset fails here with a
%named diagnostic rather than passing a bare "the figure exists" check.
%
%Returns ALL matches, not the first. A ngroups x nevents figure holds one
%curve per group per subplot, all sharing the same XData; returning the first
%would silently check one of four and report success for the rest.
lines = findall(fig, 'Type', 'line');
keep = false(size(lines));

for i = 1:numel(lines)
    xd = lines(i).XData;
    keep(i) = numel(xd) == numel(expectedX) && ...
        isequal(xd(:).', double(expectedX(:).'));
end

if ~any(keep)
    error('tests:rc30curve', ...
        ['No line in the figure had XData equal to the requested trial range '...
         '%d:%d. Lines found had %s points.'], ...
        expectedX(1), expectedX(end), ...
        strjoin(arrayfun(@(h) num2str(numel(h.XData)), lines, ...
            'UniformOutput', false), ', '));
end

lns = lines(keep);
end
