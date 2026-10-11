classdef TestRestrictedDisclosureWording < PsyRATTestBase
    %TESTRESTRICTEDDISCLOSUREWORDING Wording contract for the RESTRICTED
    %disclosure of the three gamma two-facet designs (analyses 14, 25 and 27,
    %location-scale and log-nu; the LS27 wording residual) and for the
    %dynamic-migration wording of psyrat_reprocess_gamma_result (the Q8
    %wording residual). Row 27 of the 0.1.1-beta release queue.
    %
    %WHAT IS PINNED. The scope of the label stays exactly the three designs
    %under gamma (CONTROL). The five wording corrections, each checked on
    %every surface that carries the finding: (1) the three designs are named
    %and no statement covers every in-repo check; (2) analysis 27's D caveat
    %is stated beside its G caveat, with dataset-specific offsets on the long
    %surfaces; (3) the component ratios are on the variance scale, the ranges
    %belong to the five scope captures, and the retained control's 1.42, 1.66
    %and 1.02 are given without being called a pass; (4) reading the
    %components against their intervals does not remedy the evidence, and the
    %interval-miss statement stays; (5) the 2e-12 agreement is implementation
    %consistency, not recovery. No universal maximum is claimed (CONTROL). The
    %"within about .03 and .06" observation stays in a clause of its own that
    %names analyses 14 and 25, and analysis 27's offsets are named as medians
    %across participants of the posterior-mean coefficient's signed error, with
    %B's positive D offset only beside its non-convergence. For Q8, the
    %migration tool's help and record (including a first migration with an
    %explicit surface target), and Appendix C, say that the raw draws are
    %unchanged while derived values can change; three draw-backed CONTROLS pin
    %the facts that wording states. Every rule is also shown to catch a seeded
    %bad sentence and to accept its correction, so no rule is green only
    %because its subject is absent.
    %
    %HOW. LS27WordingRules holds the surfaces and the content rules; no rule
    %compares a paragraph with a stored copy. LS27WordingFixture holds the
    %draw fixtures (the factory's deterministic analysis-27, -26 and -14
    %draws) and an independent recomputation of analysis 27's per-participant
    %coefficients from the raw draws under the old and the new residual
    %target. No test fits a model, writes a file or skips.
    %
    %NOT PINNED HERE. Manual chapter 16's status cards (read, unchanged; manual
    %chapters 9 and 15 are pinned for the scope rules only), the saved-figure
    %carriers (TestRestrictedFigureDisclosure), the scripted route (manual
    %chapter 14, another increment), the rendered layout (a native look), and
    %the LS27 recovery test's assertion and recorded outcome, which stay
    %unchanged.

    methods (Test)

        % ------------------------------------------------------------ scope

        function testScopeStaysAnalyses14And25And27UnderGammaOnly(testCase)
            % CONTROL. The label covers the six arms and nothing else: the
            % three exact analysis labels under gamma, both parameterizations.
            % The gamma occasion designs outside the ruling (5, 10 and 29), the
            % one-facet dynamic designs (11, 26), the prefix cousin 15 and the
            % Gaussian arms stay silent.
            R = psyrat_gamma_twofacet_restriction();
            testCase.verifyEqual(sort(R.designs(:)), ...
                sort({'ic_dynrel_trt'; 'ic_dynrel_sserrvar_trt'; 'trt_sserrvar'}), ...
                'the label must cover exactly the three analysis labels');
            testCase.verifyEqual(R.label, 'RESTRICTED');
            key = 'Design maturity: RESTRICTED';
            for k = 1:numel(R.designs)
                for gs = [1 2]
                    L = psyrat_provenance_lines(struct('rel', struct( ...
                        'family', 'gamma', 'analysis', R.designs{k}, 'gammascale', gs)));
                    testCase.verifyEqual(nnz(contains(L, key)), 1, sprintf( ...
                        '%s, gammascale %d: expected one RESTRICTED line', R.designs{k}, gs));
                end
                L = psyrat_provenance_lines(struct('rel', struct( ...
                    'family', 'gaussian', 'analysis', R.designs{k})));
                testCase.verifyFalse(any(contains(L, key)), sprintf( ...
                    'Gaussian %s must stay Supported', R.designs{k}));
            end
            outside = {'trt', 'analysis 5'; 'trt_diff', 'analysis 10'; ...
                'ic_dodiff_dynrel_sserrvar_trt', 'analysis 29'; ...
                'ic_dynrel', 'analysis 11'; 'ic_dynrel_sserrvar', 'analysis 26'; ...
                'ic_diff_dynrel_trt', 'analysis 15 (label contains dynrel_trt)'};
            for k = 1:size(outside, 1)
                L = psyrat_provenance_lines(struct('rel', struct( ...
                    'family', 'gamma', 'analysis', outside{k, 1})));
                testCase.verifyFalse(any(contains(L, key)), sprintf( ...
                    'gamma %s is outside the ruled scope and must not be labeled', outside{k, 2}));
            end
            for n = [14 25 27]
                Rn = psyrat_gamma_twofacet_restriction(n);
                testCase.verifySubstring(Rn.dialog{1}, sprintf('analysis %d', n));
            end
        end

        % ------------------------------------------------ the five corrections

        function testScopeWordingNamesThreeDesignsWithoutUniversals(testCase)
            % Correction (1). The export header names "three gamma two-facet
            % designs (14, 25 and 27)"; no surface says "the gamma two-facet
            % designs" (analyses 5, 10 and 29 are gamma occasion designs too),
            % and no surface makes a statement about every in-repo check:
            % analysis 5 is fitted at five occasions in the repository, so "no
            % in-repo check covers more than four occasion levels" holds only
            % for the reviewed captures of these designs.
            testCase.verifyRule(@LS27WordingRules.scopeWording);
        end

        function testAnalysis27DependabilityCaveatIsStated(testCase)
            % Correction (2). Every surface that states the coefficients names
            % analysis 27's D bias beside its G bias IN ONE CLAUSE. A check
            % that a surface merely contains "dependab" is not used: the
            % header's clause on analyses 14 and 25 already says
            % "dependability", so such a check passes on the text that omits
            % the caveat. The export header, the help text and the manual also
            % give the dataset-specific D offsets (-.13 in the converged
            % prospective fit, -.02 in the retained control); B's +.0766 may
            % appear only beside its non-convergence, and no clause puts 27
            % inside the "within about .03 and .06" statement.
            testCase.verifyRule(@LS27WordingRules.analysis27Caveat);
        end

        function testRatiosAreVarianceScaleAndAttributedToFiveCaptures(testCase)
            % Correction (3), with the unit clarification of (2). A surface
            % that prints component ratios says they are on the variance
            % scale (viewer tables show SD summaries); every sentence with the
            % 2.4-12.1 range attributes it to the five scope captures; the long
            % surfaces also give the retained control's 1.42, 1.66 and 1.02;
            % no sentence with 1.02 reads as a pass.
            testCase.verifyRule(@LS27WordingRules.ratioWording);
        end

        function testIntervalsAreNotPresentedAsARemedy(testCase)
            % Correction (4). The statement that some intervals miss the truth
            % stays wherever it is today; the surfaces that advise reading the
            % components against their intervals say that this does not
            % remedy the weak evidence.
            testCase.verifyRule(@LS27WordingRules.intervalWording);
        end

        function testEnumerationAgreementIsImplementationConsistencyOnly(testCase)
            % Correction (5). The 2e-12 agreement between production and an
            % independent enumeration is stated as implementation consistency,
            % not recovery, identifiability or sampling validation; "not a
            % formula error", "not at fault" and an unhedged causal account of
            % the mechanism are gone.
            testCase.verifyRule(@LS27WordingRules.consistencyWording);
        end

        function testNoUniversalMaximumIsStated(testCase)
            % CONTROL (guard for the rewrite). "Up to .13" and its variants are
            % not bounds: the offsets are dataset-specific.
            testCase.verifyRule(@LS27WordingRules.maximumWording);
        end

        function testRulesCatchSeededViolationsAndAcceptTheirCorrections(testCase)
            % A forbid-only rule is green when its subject is absent, so each
            % rule must catch a seeded bad sentence and accept the seeded
            % correction. Rows: rule, bad text, corrected text, surface flags.
            coef = {'coefficients', true};
            rules = struct('maximum', @LS27WordingRules.maximumWording, ...
                'caveat', @LS27WordingRules.analysis27Caveat, ...
                'ratio', @LS27WordingRules.ratioWording, ...
                'scope', @LS27WordingRules.scopeWording, ...
                'consistency', @LS27WordingRules.consistencyWording, ...
                'interval', @LS27WordingRules.intervalWording);
            seeds = { ...
                'maximum', 'D was up to .13 low for analysis 27.', ...
                'The median D error in this dataset was -.13.', {}; ...
                'caveat', 'For analysis 27, G and D are biased downward. The second fit put D .08 high.', ...
                ['For analysis 27, G and D are biased downward. The second fit put D .08 high, ' ...
                'but that fit did not converge.'], coef; ...
                'ratio', 'The control read 1.42, 1.66 and 1.02, which is a new pass.', ...
                'The control read 1.42, 1.66 and 1.02, which is not a new pass.', {}; ...
                'caveat', 'G and D within about .03/.06 for analyses 14 and 25, biased low for analysis 27.', ...
                'G and D within about .03/.06 for analyses 14 and 25; G and D biased low for analysis 27.', coef; ...
                'caveat', 'G and D stay within about .03 and .06 of the truth; G and D are biased low for analysis 27.', ...
                ['G and D stay within about .03 and .06 of the truth for analyses 14 and 25; G and D ' ...
                'are biased low for analysis 27.'], coef; ...
                'caveat', 'G is biased low for analysis 27.', 'G and D are biased low for analysis 27.', coef; ...
                'caveat', 'For analysis 27, G and D are both biased downward; D was -.13 in the converged fit.', ...
                ['For analysis 27, G and D are both biased downward; as medians across participants of ' ...
                'the posterior-mean coefficient minus its planted value, D was -.13 in the converged fit.'], coef; ...
                'scope', 'The gamma two-facet designs are Restricted.', ...
                'Three gamma two-facet designs (analyses 14, 25 and 27) are Restricted.', {}; ...
                'scope', 'No in-repo check covers more than four occasion levels.', ...
                'No reviewed capture of these designs had more than four occasion levels.', {}; ...
                'ratio', 'In the captures the occasion component came out 2.4 to 12.1 times its planted value.', ...
                ['On the variance scale, in the five scope captures the occasion component came out ' ...
                '2.4 to 12.1 times its planted value.'], {}; ...
                'consistency', 'Production agreed with an independent enumeration to 2e-12, so this is not a formula error.', ...
                ['Production agreed with an independent enumeration to 2e-12, which shows ' ...
                'implementation consistency, not recovery.'], {}; ...
                'interval', 'Read the component estimates against their intervals.', ...
                'Reading the component estimates against their intervals does not remedy the weak evidence.', ...
                {'remedy', true}};
            for k = 1:size(seeds, 1)
                rule = rules.(seeds{k, 1});
                flags = seeds{k, 4};
                bad = rule(LS27WordingRules.surface(sprintf('seed %d (bad)', k), seeds{k, 2}, flags{:}));
                good = rule(LS27WordingRules.surface(sprintf('seed %d (good)', k), seeds{k, 3}, flags{:}));
                testCase.verifyNotEmpty(bad, sprintf('%s rule missed seeded violation %d: "%s"', ...
                    seeds{k, 1}, k, seeds{k, 2}));
                testCase.verifyEmpty(good, strjoin([{sprintf('%s rule rejected seeded correction %d:', ...
                    seeds{k, 1}, k)}, good], newline));
            end
            % Q8: ambiguous "no numeric change"; the surface-change and
            % analysis-27 statements a record or paragraph must carry
            q8 = { ...
                {'raw draws unchanged; no numeric change'}, {'raw draws unchanged'}, [false false false]; ...
                {'raw draws unchanged; target stamped'}, ...
                {'raw draws unchanged; a surface_target other than the default can change the dynamic surface'}, ...
                [true false true]; ...
                {'raw draws unchanged'}, ...
                {'raw draws unchanged; analysis-27 per-participant coefficients can change at the next summary'}, ...
                [true true false]};
            for k = 1:size(q8, 1)
                f = num2cell(q8{k, 3});
                bad = LS27WordingRules.q8Wording(sprintf('Q8 seed %d (bad)', k), q8{k, 1}, f{:});
                good = LS27WordingRules.q8Wording(sprintf('Q8 seed %d (good)', k), q8{k, 2}, f{:});
                testCase.verifyNotEmpty(bad, sprintf('Q8 rule missed seeded violation %d', k));
                testCase.verifyEmpty(good, strjoin([{sprintf('Q8 rule rejected seeded correction %d:', k)}, good], newline));
            end
        end

        % ------------------------------------------------------- Q8 wording

        function testMigrationHelpAndRecordSeparateRawDrawsFromDerivedValues(testCase)
            % The tool's help and the record it returns (and stores in the
            % result) must not say "no numeric change" without saying WHAT is
            % unchanged: the raw draws are, but a stamped analysis-27 result
            % is summarized under the version-2 remainder scheme, and an
            % explicit surface target changes a dynamic surface.
            v = LS27WordingRules.q8Wording('psyrat_reprocess_gamma_result help', ...
                LS27WordingRules.sentences(LS27WordingRules.normalize( ...
                help('psyrat_reprocess_gamma_result'))), true, true, false);
            for param = {'lognu', 'ls'}
                [~, r27] = psyrat_reprocess_gamma_result(LS27WordingFixture.dyn27(param{1}));
                v = [v, LS27WordingRules.q8Wording(sprintf('record, analysis 27 %s migrated', param{1}), ...
                    LS27WordingRules.recordUnits(r27), true, true, false)]; %#ok<AGROW>
                [u14, r14] = psyrat_reprocess_gamma_result(LS27WordingFixture.dyn14(param{1}));
                v = [v, LS27WordingRules.q8Wording(sprintf('record, analysis 14 %s migrated', param{1}), ...
                    LS27WordingRules.recordUnits(r14), true, false, false)]; %#ok<AGROW>
                [~, rT] = psyrat_reprocess_gamma_result(u14, 'surface_target', 'typical_person');
                testCase.assertEqual(rT.status, 'target_restamped', 'fixture: the restamp route was not reached');
                v = [v, LS27WordingRules.q8Wording(sprintf('record, analysis 14 %s restamped', param{1}), ...
                    LS27WordingRules.recordUnits(rT), true, false, true)]; %#ok<AGROW>
                % a FIRST migration with an explicit target (Appendix C.2's
                % example) also moves the surface, so its record must say so
                [~, rF14] = psyrat_reprocess_gamma_result(LS27WordingFixture.dyn14(param{1}), ...
                    'surface_target', 'typical_person');
                testCase.assertEqual(rF14.status, 'migrated', 'fixture: the first-migration route was not reached');
                v = [v, LS27WordingRules.q8Wording(sprintf('record, analysis 14 %s first migration to typical_person', ...
                    param{1}), LS27WordingRules.recordUnits(rF14), true, false, true)]; %#ok<AGROW>
                [~, rF27] = psyrat_reprocess_gamma_result(LS27WordingFixture.dyn27(param{1}), ...
                    'surface_target', 'typical_person');
                testCase.assertEqual(rF27.status, 'migrated', 'fixture: the first-migration route was not reached');
                v = [v, LS27WordingRules.q8Wording(sprintf('record, analysis 27 %s first migration to typical_person', ...
                    param{1}), LS27WordingRules.recordUnits(rF27), true, true, true)]; %#ok<AGROW>
            end
            testCase.verifyEmpty(v, strjoin([{'Q8 wording, migration tool'}, v], newline));
        end

        function testAppendixCSeparatesRawDrawsFromDerivedValues(testCase)
            % Appendix C's paragraph on the dynamic designs says the raw draws
            % are unchanged and that analysis 27's per-participant
            % coefficients change at the next summary; the appendix says an
            % explicit surface target can change a surface.
            root = testCase.projectRoot();
            para = LS27WordingRules.appendixDynamicParagraph(root);
            testCase.assertNotEmpty(para, ['Appendix C''s dynamic-designs paragraph was not ' ...
                'found (anchor: a block starting "**Dynamic" that mentions 27)']);
            v = LS27WordingRules.q8Wording('Appendix C, dynamic designs', ...
                LS27WordingRules.sentences(LS27WordingRules.normalize(para)), true, true, false);
            whole = fileread(fullfile(root, 'documentation', 'manual', 'appendix-c_reprocessing.md'));
            v = [v, LS27WordingRules.q8Wording('Appendix C', ...
                LS27WordingRules.sentences(LS27WordingRules.normalize(whole)), false, false, true)];
            testCase.verifyEmpty(v, strjoin([{'Q8 wording, Appendix C'}, v], newline));
        end

        function testMigrationKeepsRawDrawsWhileAnalysis27SubjectCoefficientsChange(testCase)
            % CONTROL (the facts the Q8 wording states). For an unstamped
            % analysis-27 result, both parameterizations, constant and varied
            % draws: the migration adds the stamp and changes no raw array; the
            % summary before it matches an independent recomputation under the
            % old target (E = C), the summary after it one under the new
            % target (E = R + C); every participant's G and D fall because R
            % is positive; the surface does not move; and the old values stay
            % reproducible from the same draws on request.
            for param = {'lognu', 'ls'}
                for variant = {'constant', 'varied'}
                    label = sprintf('analysis 27 %s %s', param{1}, variant{1});
                    pd0 = LS27WordingFixture.dyn27(param{1}, variant{1});
                    s0 = psyrat_dynrel_summary(pd0.rel, 'CI', .95);
                    [pd1, rec] = psyrat_reprocess_gamma_result(pd0);
                    s1 = psyrat_dynrel_summary(pd1.rel, 'CI', .95);
                    testCase.verifyEqual(rec.status, 'migrated', label);

                    % raw arrays and result metadata: only the stamp was added
                    testCase.verifyEqual(rmfield(pd1.rel.out, 'gamma_target'), pd0.rel.out, ...
                        [label ': the migration changed a raw array']);
                    testCase.verifyEqual(rmfield(pd1.rel, 'out'), rmfield(pd0.rel, 'out'), ...
                        [label ': the migration changed result metadata']);
                    testCase.verifyEqual(s0.gamma_target.remainder_scheme, 'legacy_observation_only_v1', label);
                    testCase.verifyEqual(s1.gamma_target.remainder_scheme, ...
                        'population_highest_order_plus_person_observation_v2', label);

                    % independent recomputation under the two targets
                    E0 = LS27WordingFixture.expected27(pd0.rel, 'legacy');
                    E1 = LS27WordingFixture.expected27(pd0.rel, 'v2');
                    testCase.assertGreaterThan(min(E1.R_mean), 0, ...
                        [label ': fixture design: the remainder R must be positive']);
                    T0 = s0.strata(1).ssrel_table;
                    T1 = s1.strata(1).ssrel_table;
                    testCase.assertEqual(T0.id, E0.id, [label ': row order']);
                    for c = {'gen_pt', 'gen_ll', 'gen_ul', 'dep_pt', 'dep_ll', 'dep_ul'}
                        testCase.verifyEqual(T0.(c{1}), E0.(c{1}), 'AbsTol', 1e-12, ...
                            sprintf('%s before migration: %s', label, c{1}));
                        testCase.verifyEqual(T1.(c{1}), E1.(c{1}), 'AbsTol', 1e-12, ...
                            sprintf('%s after migration: %s', label, c{1}));
                    end
                    testCase.verifyTrue(all(E1.gen_pt < E0.gen_pt) && all(E1.dep_pt < E0.dep_pt), ...
                        [label ': adding R >= 0 to the residual must lower every G and D']);

                    % the surface is the same raw-draw computation either way
                    testCase.verifyEqual(s1.strata(1).G, s0.strata(1).G, [label ': surface G moved']);
                    testCase.verifyEqual(s1.strata(1).D, s0.strata(1).D, [label ': surface D moved']);

                    % the pre-migration values stay reproducible from the draws
                    sX = psyrat_dynrel_summary(pd1.rel, 'CI', .95, ...
                        'remainder_scheme', 'legacy_observation_only_v1');
                    testCase.verifyEqual(sX.strata(1).ssrel_table.gen_pt, T0.gen_pt, label);
                    testCase.verifyEqual(sX.strata(1).ssrel_table.dep_pt, T0.dep_pt, label);
                end
            end
        end

        function testFirstMigrationWithExplicitTargetMovesOnlyTheDerivedSurface(testCase)
            % CONTROL (the fact the corrected Q8 record states). A first migration
            % of an unstamped result with an explicit typical_person target keeps
            % every raw array. The next summary's surface is exactly the
            % typical-person surface of the unstamped draws, so the change is the
            % target and not the draws, and it differs from the historical
            % population-marginal surface. For analysis 27 the per-participant
            % coefficients follow the version-2 scheme (independent
            % recomputation from the raw draws). No refit.
            cases = {'14', 'lognu'; '14', 'ls'; '27', 'lognu'; '27', 'ls'};
            for k = 1:size(cases, 1)
                label = sprintf('analysis %s %s', cases{k, 1}, cases{k, 2});
                if strcmp(cases{k, 1}, '14')
                    a = LS27WordingFixture.dyn14(cases{k, 2});
                else
                    a = LS27WordingFixture.dyn27(cases{k, 2}, 'constant');
                end
                s0 = psyrat_dynrel_summary(a.rel, 'CI', .95);
                sT = psyrat_dynrel_summary(a.rel, 'CI', .95, 'surface_target', 'typical_person');
                [b, rec] = psyrat_reprocess_gamma_result(a, 'surface_target', 'typical_person');
                s1 = psyrat_dynrel_summary(b.rel, 'CI', .95);
                testCase.verifyEqual(rec.status, 'migrated', label);
                testCase.verifyEqual(b.rel.out.gamma_target.surface_target, 'typical_person', label);
                testCase.verifyEqual(rmfield(b.rel.out, 'gamma_target'), a.rel.out, ...
                    [label ': the first migration changed a raw array']);
                testCase.verifyEqual(rmfield(b.rel, 'out'), rmfield(a.rel, 'out'), ...
                    [label ': the first migration changed result metadata']);
                testCase.verifyEqual(s1.strata(1).G, sT.strata(1).G, ...
                    [label ': the migrated surface must be the typical-person surface of the same draws']);
                testCase.verifyEqual(s1.strata(1).D, sT.strata(1).D, label);
                testCase.verifyGreaterThan(max(abs(s1.strata(1).G.pt(:) - s0.strata(1).G.pt(:))), 1e-6, ...
                    [label ': the explicit target must move the surface from the historical default']);
                if strcmp(cases{k, 1}, '27')
                    E1 = LS27WordingFixture.expected27(a.rel, 'v2');
                    testCase.verifyEqual(s1.strata(1).ssrel_table.gen_pt, E1.gen_pt, 'AbsTol', 1e-12, label);
                    testCase.verifyEqual(s1.strata(1).ssrel_table.dep_pt, E1.dep_pt, 'AbsTol', 1e-12, label);
                end
            end
        end

        function testOtherDynamicDesignsChangeOnlyThroughAnExplicitTarget(testCase)
            % CONTROL (no false assurance in either direction). A default
            % migration changes no derived value of analysis 14 (surface) or
            % analysis 26 (surface and per-participant table); an explicit
            % typical-person target moves the analysis-14 surface. The raw
            % arrays never change.
            for param = {'lognu', 'ls'}
                pd0 = LS27WordingFixture.dyn14(param{1});
                s0 = psyrat_dynrel_summary(pd0.rel, 'CI', .95);
                [pd1, r1] = psyrat_reprocess_gamma_result(pd0);
                [pd2, r2] = psyrat_reprocess_gamma_result(pd1, 'surface_target', 'typical_person');
                s1 = psyrat_dynrel_summary(pd1.rel, 'CI', .95);
                s2 = psyrat_dynrel_summary(pd2.rel, 'CI', .95);
                label = ['analysis 14 ' param{1}];
                testCase.verifyEqual({r1.status, r2.status}, {'migrated', 'target_restamped'}, label);
                testCase.verifyEqual(rmfield(pd1.rel.out, 'gamma_target'), pd0.rel.out, label);
                testCase.verifyEqual(rmfield(pd2.rel.out, 'gamma_target'), pd0.rel.out, label);
                testCase.verifyEqual(s1.strata(1).G, s0.strata(1).G, [label ': default migration moved G']);
                testCase.verifyEqual(s1.strata(1).D, s0.strata(1).D, [label ': default migration moved D']);
                testCase.verifyGreaterThan(max(abs(s2.strata(1).G.pt(:) - s1.strata(1).G.pt(:))), 1e-6, ...
                    [label ': the typical-person target must move the surface']);

                q0 = LS27WordingFixture.dyn26(param{1});
                t0 = psyrat_dynrel_summary(q0.rel, 'CI', .95);
                [q1, ~] = psyrat_reprocess_gamma_result(q0);
                t1 = psyrat_dynrel_summary(q1.rel, 'CI', .95);
                label = ['analysis 26 ' param{1}];
                testCase.verifyEqual(rmfield(q1.rel.out, 'gamma_target'), q0.rel.out, label);
                testCase.verifyEqual(t1.strata(1).ssrel_table, t0.strata(1).ssrel_table, ...
                    [label ': default migration moved a per-participant value']);
                testCase.verifyEqual(t1.strata(1).G, t0.strata(1).G, [label ': default migration moved G']);
            end
        end

    end

    methods (Access = private)

        function verifyRule(testCase, rule)
            %Apply RULE to every surface; one verification per surface, so the
            %failure lists each surface's violations under its own name.
            S = LS27WordingRules.surfaces(testCase.projectRoot());
            for k = 1:numel(S)
                v = rule(S(k));
                testCase.verifyEmpty(v, strjoin([{S(k).name}, v], newline));
            end
        end

    end
end
