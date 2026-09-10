classdef TestPreflightValidation < PsyRATTestBase

    methods (Test)
        function testPreflightPassesValidData(testCase)
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.2; 0.9; 1.0], ...
                {'A';'B';'A';'B'}, 'VariableNames', {'id','meas','event'});

            report = psyrat_preflight_validate( ...
                'datatable', tbl, ...
                'sserrvar', 1, ...
                'diffest', 2, ...
                'savename', 'output.psyrat', ...
                'savepath', tempdir, ...
                'requireSaveName', true, ...
                'requireSavePath', true);

            testCase.verifyTrue(report.ok);
            testCase.verifyEmpty(report.errors);
        end

        function testPreflightDetectsMissingColumns(testCase)
            tbl = table({'s1';'s2'}, [1.1; 2.2], ...
                'VariableNames', {'subject','measurement'});

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyFalse(report.ok);
            testCase.verifyEqual(report.errors(1).id, 'varargin:colheaders');
            testCase.verifyNotEmpty(strfind(report.errors(1).message, 'How to fix:')); %#ok<STREMP>
        end

        function testPreflightRejectsZeroRowTable(testCase)
            %A table with headers but no data rows must be rejected here. It
            %was not: a 0-row table has size [0 nvar], so isempty() is true for
            %it, and the old "if ~isempty(opts.datatable)" gate treated it as
            %"no table supplied" and skipped every table validation. Measured
            %against the pre-guard source, all four cases below returned
            %report.ok = 1 with zero errors -- including the wrong-headers one,
            %because the missing-header check never ran either.
            tbl = table(cell(0,1), zeros(0,1), cell(0,1), ...
                'VariableNames', {'id','meas','event'});

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyFalse(report.ok, ...
                'A dataset with no rows must not pass preflight.');
            testCase.verifyEqual(report.errors(1).id, 'preflight:emptyTable');
            testCase.verifyTrue(contains(report.errors(1).message, 'How to fix:'));

            %guards the fixture: if the table ever gained a row this would pass
            %for the wrong reason
            testCase.verifyEqual(height(tbl), 0, ...
                'The fixture must actually have zero rows, or this proves nothing.');
        end

        function testPreflightZeroRowRejectionPreemptsTheEventChecks(testCase)
            %The emptiness error must be the ONLY one raised, and must not be
            %replaced by the diff/DoD event-count errors, whose advice ("select
            %two event conditions", "limit selected events to exactly 4")
            %cannot be acted on for a file that has no rows at all.
            tbl = table(cell(0,1), zeros(0,1), cell(0,1), ...
                'VariableNames', {'id','meas','event'});

            for diffest = [2 3]
                report = psyrat_preflight_validate('datatable', tbl, ...
                    'diffest', diffest);

                testCase.verifyFalse(report.ok);
                testCase.verifyNumElements(report.errors, 1, ...
                    sprintf('diffest=%d must raise exactly one error', diffest));
                testCase.verifyEqual(report.errors(1).id, 'preflight:emptyTable', ...
                    sprintf('diffest=%d must blame the empty file, not the event count', diffest));
            end
        end

        function testPreflightZeroRowRejectionOutranksBadHeaders(testCase)
            %Same table, entirely wrong headers. The file is the problem, so the
            %emptiness error is what the user gets -- not colheaders advice they
            %cannot act on. Pre-guard this returned ok = 1 with no errors at all.
            tbl = table(cell(0,1), zeros(0,1), ...
                'VariableNames', {'wrongname','alsowrong'});

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyFalse(report.ok);
            testCase.verifyEqual(report.errors(1).id, 'preflight:emptyTable');
        end

        function testPreflightStillIgnoresAnAbsentDatatable(testCase)
            %Positive control for the gate restructure. opts.datatable defaults
            %to [], which is ALSO isempty, and that case must keep meaning "no
            %table was supplied" -- no table error at all. Without this pin the
            %new height() branch could be reached by the default and every
            %caller that omits the datatable would start failing.
            report = psyrat_preflight_validate('savename', 'output.psyrat');

            ids = {report.errors.id};
            testCase.verifyFalse(any(strcmp(ids, 'preflight:emptyTable')), ...
                'Omitting the datatable must not raise the empty-table error.');
            testCase.verifyFalse(any(strcmp(ids, 'varargin:datatabletype')), ...
                'Omitting the datatable must not raise the wrong-type error.');
        end

        function testPreflightStillRejectsANonTableDatatable(testCase)
            %The other half of the same control: a non-empty non-table must keep
            %producing varargin:datatabletype. The restructure moved that arm
            %from an inner "if ~istable" to an elseif, so it needs its own pin.
            report = psyrat_preflight_validate('datatable', magic(3));

            testCase.verifyFalse(report.ok);
            testCase.verifyEqual(report.errors(1).id, 'varargin:datatabletype');
        end

        function testPreflightRejectsLogicalGroupColumn(testCase)
            %B27 (owner ruling 2026-08-17: reject, do not convert): a LOGICAL
            %label column bypasses psyrat_normalize_label_columns (whose own
            %comment deliberately declines to convert logical -- inventing a
            %third label spelling for a type no shipped path converts was
            %outside that fix's scope) and reaches the engine raw, where the
            %downstream behavior is unpredictable and case-dependent: some
            %chunkers die on control-character labels, the gamma
            %string()-keyed cases complete while INVENTING 'false'/'true'
            %labels, and a logical id can silently no-op the test-retest
            %participant-wide exclusion. Preflight runs on every entry route
            %including the direct-table one, so the TYPE is rejected here
            %with a convert-it-yourself message instead.
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.2; 0.9; 1.0], ...
                logical([0; 1; 0; 1]), 'VariableNames', {'id','meas','group'});

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'labels:logicaltype')), ...
                'a logical group column must be rejected at preflight (B27)');
        end

        function testPreflightRejectsLogicalEventColumn(testCase)
            %B27, second declared label column: same rejection for event
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.2; 0.9; 1.0], ...
                logical([0; 1; 0; 1]), 'VariableNames', {'id','meas','event'});

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'labels:logicaltype')), ...
                'a logical event column must be rejected at preflight (B27)');
        end

        function testPreflightAcceptsNumericLabelColumns(testCase)
            %B27 negative control: a numeric 0/1 group column is exactly what
            %the rejection message tells the user to convert TO (the
            %normalizer spells it '0'/'1'), so the type guard must not fire
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.2; 0.9; 1.0], ...
                [0; 1; 0; 1], 'VariableNames', {'id','meas','group'});

            report = psyrat_preflight_validate('datatable', tbl);

            ids = {};
            if ~isempty(report.errors)
                ids = {report.errors.id};
            end
            testCase.verifyFalse(any(strcmp(ids, 'labels:logicaltype')), ...
                'a numeric 0/1 group column must pass the B27 type guard');
        end

        function testPreflightLogicalGuardResolvesExactNames(testCase)
            %B27 resolution pin (refuter-driven): the guard must inspect the
            %EXACT-named column -- the one the normalizer's ismember and the
            %engine's case-sensitive datatable.group dot access consume --
            %never a case-variant decoy. A non-logical 'Group' decoy ordered
            %before a logical 'group' must not hide the logical column the
            %engine actually reads.
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.2; 0.9; 1.0], ...
                {'a';'a';'b';'b'}, logical([0; 1; 0; 1]), ...
                'VariableNames', {'id','meas','Group','group'});

            report = psyrat_preflight_validate('datatable', tbl);

            ids = {};
            if ~isempty(report.errors)
                ids = {report.errors.id};
            end
            testCase.verifyTrue(any(strcmp(ids, 'labels:logicaltype')), ...
                'a case-variant decoy must not hide the exact-named logical column');
        end

        function testPreflightLogicalGuardIgnoresCaseVariantDecoys(testCase)
            %B27 resolution pin, mirror image: a LOGICAL case-variant decoy
            %('Group') beside a healthy exact-named 'group' must not trigger
            %the rejection -- the engine reads only the exact name, so the
            %table runs fine.
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.2; 0.9; 1.0], ...
                logical([0; 1; 0; 1]), {'a';'a';'b';'b'}, ...
                'VariableNames', {'id','meas','Group','group'});

            report = psyrat_preflight_validate('datatable', tbl);

            ids = {};
            if ~isempty(report.errors)
                ids = {report.errors.id};
            end
            testCase.verifyFalse(any(strcmp(ids, 'labels:logicaltype')), ...
                'a logical case-variant decoy must not reject a healthy exact-named column');
        end

        function testPreflightDetectsNonNumericMeasurement(testCase)
            tbl = table({'s1';'s2'}, {'x';'y'}, ...
                'VariableNames', {'id','meas'});

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'meas:notnumeric')));
        end

        function testPreflightDetectsMissingCells(testCase)
            tbl = table({'s1';'s2'}, [1.2; NaN], ...
                'VariableNames', {'id','meas'});

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'varargin:missingcells')));
        end

        function testPreflightGammaRejectsNonPositiveMeas(testCase)
            % The Gamma / scaled-chi-square family models strictly positive,
            % mean-linked data (e.g. single-trial time-frequency power). A
            % non-positive measurement is outside its support and must be a
            % blocking error before any expensive estimation.
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 0.0; 0.9; 1.0], ...
                'VariableNames', {'id','meas'});

            report = psyrat_preflight_validate('datatable', tbl, 'family', 'gamma');

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'meas:notpositive')));
        end

        function testPreflightGammaAllowsPositiveMeas(testCase)
            % All strictly positive scores are valid for the Gamma family.
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.2; 0.9; 1.0], ...
                'VariableNames', {'id','meas'});

            report = psyrat_preflight_validate('datatable', tbl, 'family', 'gamma');

            ids = {report.errors.id};
            testCase.verifyFalse(any(strcmp(ids, 'meas:notpositive')));
        end

        function testPreflightGaussianAllowsNonPositiveMeas(testCase)
            % The default Gaussian family must still accept non-positive values
            % (ERP amplitudes are routinely negative). The positivity guard is
            % gamma-only and must not fire here.
            tbl = table({'s1';'s1';'s2';'s2'}, [-1.1; 0.0; 0.9; 1.0], ...
                'VariableNames', {'id','meas'});

            report = psyrat_preflight_validate('datatable', tbl, 'family', 'gaussian');

            ids = {report.errors.id};
            testCase.verifyFalse(any(strcmp(ids, 'meas:notpositive')));
        end

        function testDefaultsFamilyIsGaussian(testCase)
            % The observation-likelihood family defaults to Gaussian, so existing
            % workflows are unchanged until a user opts into Gamma.
            prefs = psyrat_defaults();
            testCase.verifyEqual(prefs.proc.family, 'gaussian');
        end

        function testComputevarcompForwardsFamilyToPreflight(testCase)
            % Regression for the family wiring: psyrat_computevarcomp must forward
            % the family to the shared preflight so the gamma positivity guard
            % actually fires. Non-positive data under 'family','gamma' therefore
            % errors at the preflight (identifier meas:notpositive) before any
            % estimation. If the family were not forwarded the preflight would
            % default to Gaussian and this call would sail past the guard.
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 0.0; 0.9; 1.0], ...
                'VariableNames', {'id','meas'});
            % Minimal MCMC args so the call reaches the preflight (validated
            % before them); tiny values, since it errors before any estimation.
            testCase.verifyError( ...
                @() psyrat_computevarcomp('data', tbl, 'family', 'gamma', ...
                'chains', 2, 'warmup', 10, 'sampling', 10), ...
                'meas:notpositive');
        end

        function testPreflightAllowsTrtSserrWithTime(testCase)
            % Unconditional per-subject test-retest reliability (trt_sserrvar,
            % analysis 25): subject-level error variance WITH a time/occasion
            % facet is now supported when non-difference (diffest=1) and
            % non-dynamic (dynrel~=2). It must NOT trigger the retest-mode
            % subject-reliability block.
            tbl = table({'s1';'s1'}, [1.0; 1.2], {'t1';'t2'}, ...
                'VariableNames', {'id','meas','time'});

            report = psyrat_preflight_validate('datatable', tbl, ...
                'sserrvar', 2, 'diffest', 1, 'dynrel', 1);

            ids = {report.errors.id};
            testCase.verifyFalse(any(strcmp(ids, ...
                'preflight:timeWithSubjectReliabilityUnsupported')));
        end

        function testPreflightAllowsDynrelSubjNondiffTrtWithTime(testCase)
            % Now supported: NON-difference subject-level dynamic reliability WITH
            % a time/occasion facet (dynrel=2, diffest=1, sserrvar=2) is the
            % two-facet subject-level dynamic variant (ic_dynrel_sserrvar_trt,
            % analysis 27), reporting each participant's conditional reliability
            % phi_s(z) from the per-subject residual the location-scale model
            % already estimates. It must NOT trigger the retest-mode
            % subject-reliability block.
            tbl = table({'s1';'s1'}, [1.0; 1.2], {'t1';'t2'}, ...
                'VariableNames', {'id','meas','time'});

            report = psyrat_preflight_validate('datatable', tbl, ...
                'sserrvar', 2, 'dynrel', 2, 'diffest', 1);

            ids = {report.errors.id};
            testCase.verifyFalse(any(strcmp(ids, ...
                'preflight:timeWithSubjectReliabilityUnsupported')));
        end

        function testPreflightAllowsTwoFacetDiffWithTime(testCase)
            % Two-event difference scores WITH a time/occasion facet is the
            % supported two-facet (test-retest) difference pipeline (analysis
            % 10); it must NOT trigger timeWithDiffUnsupported.
            tbl = table({'s1';'s1';'s1';'s1'}, [1.0; 1.2; 0.9; 1.1], ...
                {'A';'B';'A';'B'}, {'t1';'t1';'t2';'t2'}, ...
                'VariableNames', {'id','meas','event','time'});

            report = psyrat_preflight_validate('datatable', tbl, 'diffest', 2);

            ids = {report.errors.id};
            testCase.verifyFalse(any(strcmp(ids, 'preflight:timeWithDiffUnsupported')));
        end

        function testPreflightDetectsUnsupportedTimeWithDiff(testCase)
            % Still unsupported with a time/occasion facet: difference-of-
            % differences (diffest=3) and subject-level difference scores
            % (diffest=2 with sserrvar=2).
            tblDoD = table({'s1';'s1';'s1';'s1'}, [1.0; 1.2; 0.9; 1.1], ...
                {'A';'B';'C';'D'}, {'t1';'t1';'t2';'t2'}, ...
                'VariableNames', {'id','meas','event','time'});
            reportDoD = psyrat_preflight_validate('datatable', tblDoD, 'diffest', 3);
            testCase.verifyFalse(reportDoD.ok);
            testCase.verifyTrue(any(strcmp({reportDoD.errors.id}, ...
                'preflight:timeWithDiffUnsupported')));

            tblSS = table({'s1';'s1';'s1';'s1'}, [1.0; 1.2; 0.9; 1.1], ...
                {'A';'B';'A';'B'}, {'t1';'t1';'t2';'t2'}, ...
                'VariableNames', {'id','meas','event','time'});
            reportSS = psyrat_preflight_validate('datatable', tblSS, ...
                'diffest', 2, 'sserrvar', 2);
            testCase.verifyFalse(reportSS.ok);
            testCase.verifyTrue(any(strcmp({reportSS.errors.id}, ...
                'preflight:timeWithDiffUnsupported')));
        end

        function testPreflightAllowsDiffWithSubjectReliability(testCase)
            tbl = table({'s1';'s1';'s2';'s2'}, [1.0; 1.2; 0.9; 1.1], ...
                {'A';'B';'A';'B'}, 'VariableNames', {'id','meas','event'});

            report = psyrat_preflight_validate('datatable', tbl, 'diffest', 2, 'sserrvar', 2);

            testCase.verifyTrue(report.ok);
            ids = {report.errors.id};
            testCase.verifyFalse(any(strcmp(ids, 'preflight:diffWithSubjectReliabilityUnsupported')));
        end

        function testPreflightBlocksDoDWithSubjectReliability(testCase)
            tbl = table({'s1';'s1';'s1';'s1'}, [1.0; 1.2; 0.9; 1.1], ...
                {'A';'B';'C';'D'}, 'VariableNames', {'id','meas','event'});

            report = psyrat_preflight_validate('datatable', tbl, 'diffest', 3, 'sserrvar', 2);

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'preflight:dodWithSubjectReliabilityUnsupported')));
        end

        function testPreflightDetectsDiffCovarianceWithoutPairs(testCase)
            tbl = table({'s1';'s2'}, [1.0; 1.2], {'A';'B'}, ...
                'VariableNames', {'id','meas','event'});

            report = psyrat_preflight_validate('datatable', tbl, ...
                'diffest', 2, 'sserrvar', 2, 'diffwpcov', 2);

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'preflight:diffCovNoPairs')));
        end

        function testPreflightDetectsDiffNeedsTwoEvents(testCase)
            tbl = table({'s1';'s1';'s2'}, [1.0; 1.2; 0.9], {'A';'A';'A'}, ...
                'VariableNames', {'id','meas','event'});

            report = psyrat_preflight_validate('datatable', tbl, 'diffest', 2);

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'preflight:diffNeedsTwoEvents')));
        end

        function testPreflightAllowsUnbalancedDiffEvents(testCase)
            tbl = table({'doors';'doors';'doors';'doors';'doors'}, ...
                {'s1';'s1';'s1';'s2';'s2'}, ...
                [1.1; 1.2; 1.3; 0.9; 0.8], ...
                {'loss';'loss';'gain';'loss';'gain'}, ...
                'VariableNames', {'group','id','meas','event'});

            report = psyrat_preflight_validate('datatable', tbl, 'diffest', 2);

            testCase.verifyTrue(report.ok);
            ids = {report.errors.id};
            testCase.verifyFalse(any(strcmp(ids, 'varargin:diffpairing')));
        end

        function testPreflightDetectsWhitespaceInOutputTargets(testCase)
            report = psyrat_preflight_validate( ...
                'savename', 'my output.psyrat', ...
                'savepath', fullfile(tempdir, 'path with spaces'), ...
                'requireSaveName', true, ...
                'requireSavePath', true);

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'psyrat_data:savename')));
            testCase.verifyTrue(any(strcmp(ids, 'psyrat_data:savepath')));
        end

        function testPreflightWarnsFewParticipants(testCase)
            % 5 participants (< 20) with 12 trials each: only the
            % few-participants warning should fire, not the trials warning.
            tbl = buildBalancedTable(5, 12);

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyTrue(report.ok);
            testCase.verifyEmpty(report.errors);
            warnIds = {report.warnings.id};
            testCase.verifyTrue(any(strcmp(warnIds, 'preflight:fewParticipants')));
            testCase.verifyFalse(any(strcmp(warnIds, 'preflight:lowMedianTrials')));
        end

        function testPreflightWarnsLowMedianTrials(testCase)
            % 25 participants (>= 20) with 4 trials each: only the
            % low-trials warning should fire, not the participant warning.
            tbl = buildBalancedTable(25, 4);

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyTrue(report.ok);
            warnIds = {report.warnings.id};
            testCase.verifyTrue(any(strcmp(warnIds, 'preflight:lowMedianTrials')));
            testCase.verifyFalse(any(strcmp(warnIds, 'preflight:fewParticipants')));
        end

        function testPreflightWarnsPerGroupFewParticipants(testCase)
            % Two groups, 12 trials each so no trials warning: group A has
            % 25 participants (enough), group B has 5 (too few). Exactly one
            % few-participants warning should fire, naming group B.
            ntrials = 12;
            idsA = repelem((1:25)', ntrials);
            idsB = repelem((100 + (1:5))', ntrials);
            ids = [idsA; idsB];
            grp = [repmat("A", numel(idsA), 1); repmat("B", numel(idsB), 1)];
            meas = 0.01 * (1:numel(ids))';
            tbl = table(ids, meas, grp, 'VariableNames', {'id','meas','group'});

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyTrue(report.ok);
            warnIds = {report.warnings.id};
            fewMask = strcmp(warnIds, 'preflight:fewParticipants');
            testCase.verifyEqual(sum(fewMask), 1);
            whats = {report.warnings.what};
            testCase.verifyNotEmpty(strfind(whats{find(fewMask, 1)}, 'B')); %#ok<STREMP>
        end

        function testPreflightNoWarningsForHealthyDataset(testCase)
            % 22 participants (>= 20) with 12 trials each: no warnings.
            tbl = buildBalancedTable(22, 12);

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyTrue(report.ok);
            testCase.verifyEmpty(report.warnings);
        end

        function testPreflightWarningsDoNotBlock(testCase)
            % A warning-triggering dataset is still a valid run: warnings
            % never set report.ok to false or add errors.
            tbl = buildBalancedTable(3, 3);

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyTrue(report.ok);
            testCase.verifyEmpty(report.errors);
            testCase.verifyNotEmpty(report.warnings);
        end

        function testPreflightSplitsRequiresWeightColumn(testCase)
            % Reliability with data splits requires an items-per-split (n_i)
            % 'weight' column; its absence is a blocking error.
            tbl = table((1:24)', (0.01*(1:24))', 'VariableNames', {'id','meas'});

            report = psyrat_preflight_validate('datatable', tbl, 'splits', 3);

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'preflight:splitsNeedsWeight')));
        end

        function testPreflightSplitsBlocksEqualWeights(testCase)
            % Nonparallel splits (splits=3) with identical n_i across splits
            % leave the per-item residual and person-by-split interaction
            % confounded, so it is a blocking error.
            tbl = buildSplitsTable(24, [20 20 20 20]);

            report = psyrat_preflight_validate('datatable', tbl, 'splits', 3);

            testCase.verifyFalse(report.ok);
            ids = {report.errors.id};
            testCase.verifyTrue(any(strcmp(ids, 'preflight:splitsNoVariation')));
        end

        function testPreflightSplitsWarnsLowVariation(testCase)
            % splits=3 with low but nonzero n_i variation (CV ~ 0.05) is
            % permitted but warns that components may be weakly identified.
            tbl = buildSplitsTable(24, [26 27 28 29 30]);

            report = psyrat_preflight_validate('datatable', tbl, 'splits', 3);

            testCase.verifyTrue(report.ok);
            ids = {report.warnings.id};
            testCase.verifyTrue(any(strcmp(ids, 'preflight:splitsLowVariation')));
        end

        function testPreflightSplitsAllowsWideVariation(testCase)
            % splits=3 with well-spread n_i (CV ~ 0.7) raises no splits warning.
            tbl = buildSplitsTable(24, [8 30 60 100]);

            report = psyrat_preflight_validate('datatable', tbl, 'splits', 3);

            testCase.verifyTrue(report.ok);
            ids = {report.warnings.id};
            testCase.verifyFalse(any(strcmp(ids, 'preflight:splitsLowVariation')));
            testCase.verifyFalse(any(strcmp(ids, 'preflight:splitsNoVariation')));
        end

        function testPreflightSplitsWarnsParallelUnequalWeights(testCase)
            % Parallel splits (splits=2) declared with unequal n_i is a
            % non-fatal warning (the nonparallel model may fit better).
            tbl = buildSplitsTable(24, [18 20 22 24]);

            report = psyrat_preflight_validate('datatable', tbl, 'splits', 2);

            testCase.verifyTrue(report.ok);
            ids = {report.warnings.id};
            testCase.verifyTrue(any(strcmp(ids, 'preflight:splitsUnequalWeight')));
        end

        function testPreflightSplitsParallelEqualWeightsClean(testCase)
            % Parallel splits with equal n_i raises no splits-specific warning.
            tbl = buildSplitsTable(24, [20 20 20 20]);

            report = psyrat_preflight_validate('datatable', tbl, 'splits', 2);

            testCase.verifyTrue(report.ok);
            ids = {report.warnings.id};
            testCase.verifyFalse(any(strcmp(ids, 'preflight:splitsUnequalWeight')));
        end

        function testPreflightSplitsExemptFromLowMedianTrials(testCase)
            % A split-mean row is a chunk of trials, not one trial, so the
            % splits-per-person row count must NOT be read as trials-per-cell.
            % This thin table (5 split rows per participant) would otherwise
            % trip preflight:lowMedianTrials; splits runs are exempt for both
            % parallel (splits=2) and nonparallel (splits=3).
            tbl = buildSplitsTable(24, [26 27 28 29 30]);

            for splitsmode = [2 3]
                report = psyrat_preflight_validate('datatable', tbl, ...
                    'splits', splitsmode);
                ids = {report.warnings.id};
                testCase.verifyFalse(...
                    any(strcmp(ids, 'preflight:lowMedianTrials')), ...
                    sprintf('lowMedianTrials should be suppressed for splits=%d', ...
                    splitsmode));
            end
        end

        function testPreflightWarnsIncompleteEventTimeCells(testCase)
            % B24's reachability lesson: the loader's coverage guards enforce
            % UNIFORMITY of per-participant counts, not COMPLETENESS of the
            % (event x occasion) grid, so uniform-but-incomplete data load and
            % estimate silently. The ruled fix is a preflight WARNING (never an
            % error): the user is told which cells are empty; the run proceeds.
            %
            % Fixture shape mirrors the B24 measurement fixture: every
            % participant holds 3 of the 4 (event x time) cells -- participant
            % 1 misses (correct, t1), everyone else misses (error, t2) -- so
            % BOTH loader uniformities hold (2 events and 3 cells per
            % participant). 22 participants and 12 trials per cell keep the
            % fewParticipants and lowMedianTrials warnings out of the way, so
            % the new warning is asserted in isolation.
            tbl = buildEventTimeTable(22, 12, false);

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyTrue(report.ok, ...
                'incomplete coverage must WARN, never block (ruled 2026-08-16)');
            testCase.verifyEmpty(report.errors);
            warnIds = {report.warnings.id};
            testCase.assertTrue( ...
                any(strcmp(warnIds, 'preflight:incompleteEventTimeCells')), ...
                'the incomplete-coverage warning did not fire on empty cells');
            testCase.verifyFalse(any(strcmp(warnIds, 'preflight:fewParticipants')));
            testCase.verifyFalse(any(strcmp(warnIds, 'preflight:lowMedianTrials')));

            % The message should quantify the gap (22 empty cells of the
            % 22 x 2 x 2 = 88-cell grid, every participant affected).
            mask = strcmp(warnIds, 'preflight:incompleteEventTimeCells');
            what = report.warnings(find(mask, 1)).what;
            testCase.verifyNotEmpty(strfind(what, '22 of 88')); %#ok<STREMP>
        end

        function testPreflightIncompleteCellsAbsentForCompleteData(testCase)
            % Negative control: the same builder with complete coverage must
            % not fire the warning (guards against over-firing).
            tbl = buildEventTimeTable(22, 12, true);

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyTrue(report.ok);
            warnIds = {report.warnings.id};
            testCase.verifyFalse( ...
                any(strcmp(warnIds, 'preflight:incompleteEventTimeCells')), ...
                'the coverage warning fired on a complete grid');
        end

        function testPreflightIncompleteCellsImmuneToDelimiterLabels(testCase)
            % Regression for a refuter-caught defect in the first cut of this
            % warning: present cells were keyed by joining the three labels
            % with '|', so labels CONTAINING '|' could make two DIFFERENT
            % real cells produce one key -- events {a, a|b} with times
            % {b|c, c} give (a, b|c) and (a|b, c) the identical joined key --
            % and the cardinality arithmetic fabricated missing cells out of
            % a COMPLETE grid (and the key lookup could hide real ones).
            % Membership is now row-wise; a complete grid must stay silent
            % whatever the labels. (A '|'-free label pair does NOT trigger
            % the old defect -- the collision needs both halves of the
            % ambiguous parse to be real labels -- so this fixture is the
            % colliding pair, watched to fail against the joined-key code.)
            tbl = buildEventTimeTable(22, 12, true);
            tbl.event = strrep(tbl.event, 'correct', 'a');
            tbl.event = strrep(tbl.event, 'error', 'a|b');
            tbl.time = strrep(tbl.time, 't1', 'b|c');
            tbl.time = strrep(tbl.time, 't2', 'c');

            report = psyrat_preflight_validate('datatable', tbl);

            testCase.verifyTrue(report.ok);
            warnIds = {report.warnings.id};
            testCase.verifyFalse( ...
                any(strcmp(warnIds, 'preflight:incompleteEventTimeCells')), ...
                ['the coverage warning fired on a COMPLETE grid whose labels ' ...
                'contain the old join delimiter']);
        end

        function testPreflightIncompleteCellsRequiresEventAndTime(testCase)
            % Scope control: the warning is defined on the event x occasion
            % grid, so it must be silent when either facet column is absent.
            % (An id x event gap without a time facet is the loader's own
            % mismatchedevents ERROR territory, not this warning's.)
            base = buildEventTimeTable(22, 12, false);

            report = psyrat_preflight_validate( ...
                'datatable', removevars(base, 'time'));
            warnIds = {report.warnings.id};
            testCase.verifyFalse( ...
                any(strcmp(warnIds, 'preflight:incompleteEventTimeCells')), ...
                'the coverage warning fired without a time column');

            report = psyrat_preflight_validate( ...
                'datatable', removevars(base, 'event'));
            warnIds = {report.warnings.id};
            testCase.verifyFalse( ...
                any(strcmp(warnIds, 'preflight:incompleteEventTimeCells')), ...
                'the coverage warning fired without an event column');
        end
    end
end

function tbl = buildSplitsTable(nid, weights)
% Build a long-format splits table: nid participants, each with numel(weights)
% split-mean rows carrying the given items-per-split (n_i) 'weight' pattern.
% Deterministic numeric measurements; no missing cells.
weights = weights(:)';
k = numel(weights);
ids = repelem((1:nid)', k);
w = repmat(weights', nid, 1);
meas = 0.01 * (1:nid*k)';
tbl = table(ids, meas, w, 'VariableNames', {'id','meas','weight'});
end

function tbl = buildBalancedTable(nid, ntrials)
% Build a balanced long-format id/meas table with nid participants and
% ntrials rows each. Deterministic numeric measurements (no missing cells).
ids = repelem((1:nid)', ntrials);
meas = 0.01 * (1:nid * ntrials)';
tbl = table(ids, meas, 'VariableNames', {'id','meas'});
end

function tbl = buildEventTimeTable(nid, ntrials, complete)
% Long-format id/meas/event/time table over a 2-event ('correct','error') x
% 2-occasion ('t1','t2') grid, ntrials rows per present cell. complete=true
% fills all 4 cells per participant. complete=false reproduces the B24
% uniform-but-incomplete shape: participant 1 misses (correct, t1), every
% other participant misses (error, t2), so each participant carries exactly
% 2 distinct events and 3 cells (both loader uniformity guards pass) while
% the grid has one empty cell per participant. Deterministic measurements.
ids = {}; meas = []; event = {}; time = {};
k = 0;
for s = 1:nid
    for cellspec = {{'correct','t1'}, {'correct','t2'}, {'error','t1'}, {'error','t2'}}
        ev = cellspec{1}{1}; tm = cellspec{1}{2};
        if ~complete
            if s == 1 && strcmp(ev, 'correct') && strcmp(tm, 't1')
                continue;
            elseif s > 1 && strcmp(ev, 'error') && strcmp(tm, 't2')
                continue;
            end
        end
        for t = 1:ntrials %#ok<*AGROW>
            k = k + 1;
            ids{end+1,1} = sprintf('s%02d', s);
            event{end+1,1} = ev;
            time{end+1,1} = tm;
            meas(end+1,1) = 0.01 * k;
        end
    end
end
tbl = table(ids, meas, event, time, 'VariableNames', {'id','meas','event','time'});
end
