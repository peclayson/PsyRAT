classdef TestCriterionUnsupportedDesignGuard < PsyRATTestBase
    % G42: psyrat_criterionfigures must refuse the designs whose result
    % structs lack the fields its cut-score consumers dereference, with a
    % named error where the route resolves -- not a raw
    % MATLAB:nonExistentField from deep inside local_build_relsummary's
    % downstream consumers.
    %
    % The GUI half of this is already gated (G38: the trt viewer hides the
    % criterion button for trt_sserrvar/trt_diff -- a blacklist, with the
    % startview router keeping the other families out of that viewer;
    % psyrat_startview_sing gates on ~sserrvar && ~isdiff &&
    % hasEstimatedMu). These tests cover the scripted-caller half those
    % gates cannot reach.
    %
    % The positive controls live in TestSummaryAndPresentationFunctions,
    % which drives psyrat_criterionfigures end to end on plain 'ic' and
    % plain 'trt' data and asserts the plotted point estimates against the
    % cut-score kernels. testPlainDesignsPassTheGuard below is a cheap
    % in-class sentinel of the same fact (no figures drawn), so this class
    % cannot pass vacuously if the guard ever starts rejecting everything.
    %
    % G50 widened the coverage past the five-label blacklist with a
    % CATEGORICAL family refusal ('dod'/'dynrel' plus the exact label
    % ic_splits); G55 (owner ruling 2026-08-30, reversing the G50-era
    % constraint) widened it to the whole 'splits' family. The refused
    % structs either do not carry the components this screen consumes and
    % used to die on raw errors deep inside the rebuild, or (trt_splits)
    % rebuild to numbers that are not the splits estimand. The adversarial wave
    % rejected an earlier error-translation design (catching the rebuild's
    % nonExistentField) because it mislabeled SUPPORTED designs on a
    % mismatched route as unsupported; the family check cannot capture
    % those -- they fall through and fail with the raw, factually true
    % error, pinned by the negative control below. Legacy label-free
    % structs (family 'sing' by default) are deliberately not refused.

    methods (Test)

        function testBlockedDesignsErrorWithTheNamedId(testCase)
            % Each unsupported design, paired with the route a scripted
            % caller would plausibly use. The minimal struct is enough:
            % the guard fires on rel.analysis before any data field is
            % touched, which is the point -- fail where the route resolves.
            spec = { ...
                'ic_sserrvar',      'sing'; ...
                'ic_diff',          'sing'; ...
                'ic_diff_sserrvar', 'sing'; ...
                'trt_sserrvar',     'trt'; ...
                'trt_diff',         'trt'};

            for k = 1:size(spec, 1)
                pd = struct('rel', struct('analysis', spec{k, 1}));
                testCase.verifyError(@() psyrat_criterionfigures( ...
                    'psyrat_data', pd, 'analysis', spec{k, 2}), ...
                    'criterion:unsupportedDesign', sprintf( ...
                    '%s via ''%s'' must be refused by the G42 guard', ...
                    spec{k, 1}, spec{k, 2}));
            end
        end

        function testGuardKeysOnTheDataNotTheRouteArgument(testCase)
            % The guard checks rel.analysis, not the 'analysis' argument, so
            % a mislabeled route cannot slip a blocked design past it. Pin
            % one cross-pairing: trt_sserrvar data sent down the 'sing'
            % route must still be refused.
            pd = struct('rel', struct('analysis', 'trt_sserrvar'));
            testCase.verifyError(@() psyrat_criterionfigures( ...
                'psyrat_data', pd, 'analysis', 'sing'), ...
                'criterion:unsupportedDesign');
        end

        function testExoticFamiliesRefusedCategorically(testCase)
            % G50: the DoD and dynrel families plus ic_splits must surface
            % the named refusal at the guard, whether or not a relsummary is
            % present -- their structs categorically lack the components this
            % screen consumes. Real fixtures rather than minimal structs, so
            % the label/family classification under test is the shipped one.
            prefs = struct('view', struct( ...
                'criterioncutoff', 0.5, ...
                'plotcriterion', 0, ...
                'tablecriterion', 0), 'ver', '0-test');

            fixtures = { ...
                'ic_dodiff', @() PsyRATTestDataFactory.makeDoDRelData(); ...
                'ic_splits', @() PsyRATTestDataFactory.makeSplitsSummaryData(); ...
                'ic_diff_dynrel', @() PsyRATTestDataFactory.makeDynrelDiffSummaryData()};

            for k = 1:size(fixtures, 1)
                pd = fixtures{k, 2}();
                testCase.assertEqual(pd.rel.analysis, fixtures{k, 1}, ...
                    'fixture precondition broken: unexpected analysis label');
                testCase.verifyError(@() psyrat_criterionfigures( ...
                    'psyrat_prefs', prefs, 'psyrat_data', pd, ...
                    'analysis', 'sing'), 'criterion:unsupportedDesign', ...
                    sprintf('%s must be refused with the named id, not a raw error', ...
                    fixtures{k, 1}));
                % the refusal is categorical, so the trt route must refuse
                % the same struct too (pre-G50 the trt route died raw on
                % MATLAB:cellRefFromNonCell for the splits fixtures)
                testCase.verifyError(@() psyrat_criterionfigures( ...
                    'psyrat_prefs', prefs, 'psyrat_data', pd, ...
                    'analysis', 'trt'), 'criterion:unsupportedDesign', ...
                    sprintf('%s must be refused on the trt route too', ...
                    fixtures{k, 1}));
            end
        end

        function testMismatchedSupportedDesignIsNotMislabeled(testCase)
            % Negative control the adversarial wave demanded: a SUPPORTED
            % design sent down the wrong route (plain trt data on the
            % default 'sing' route -- the route psyrat_criterionfigures
            % assumes when 'analysis' is omitted) must NOT be captured by
            % the G50 refusal. It dies on the raw, factually true
            % missing-field error at the rebuild's first sing-component
            % read, exactly as before G50; reporting it as an unsupported
            % design would be a false diagnosis (the design is supported,
            % the route was wrong) and would mask genuine field problems.
            prefs = struct('view', struct( ...
                'criterioncutoff', 0.5, ...
                'plotcriterion', 0, ...
                'tablecriterion', 0), 'ver', '0-test');
            % a REL-only trt struct (no relsummary, so the rebuild runs):
            % trt components present, sing components (sig_e) absent
            pd = PsyRATTestDataFactory.makeSplitsTRTSummaryData();
            pd.rel.analysis = 'trt';
            pd.rel = rmfield(pd.rel, 'splits');
            testCase.assertFalse(isfield(pd.rel.out, 'sig_e'), ...
                'fixture precondition broken: a sing component is present');
            testCase.verifyError(@() psyrat_criterionfigures( ...
                'psyrat_prefs', prefs, 'psyrat_data', pd, ...
                'analysis', 'sing'), 'MATLAB:nonExistentField', ...
                ['a supported design on a mismatched route must fail with ' ...
                'the raw field error, not the unsupported-design refusal']);
        end

        function testTrtSplitsRefusedSinceG55(testCase)
            % G55 (owner ruling 2026-08-30, reversing the G50-era
            % constraint): trt_splits must surface the named refusal on
            % BOTH routes. Under G50 this test's predecessor pinned the
            % opposite -- that the trt-route rebuild completed past the
            % refusal -- while its own docstring conceded the numbers were
            % not established as the splits estimand (the route drops the
            % sig_posxid term and the items-per-split weighting the splits
            % error folds in, and runs a D-study projection on a design
            % case 24 documents as observed-design-only; the sig_trl ->
            % trial-slot binding is NOT a defect -- psyrat_splits_summary
            % makes the same mapping on purpose). The estimand question
            % was adjudicated by blocking the label; a splits-aware
            % criterion path at the observed design would be a wiring
            % effort over psyrat_splits_summary's folded error, and the
            % projection machinery is the part that needs real design.
            prefs = struct('view', struct( ...
                'criterioncutoff', 0.5, ...
                'plotcriterion', 0, ...
                'tablecriterion', 0), 'ver', '0-test');
            pd = PsyRATTestDataFactory.makeSplitsTRTSummaryData();
            testCase.assertEqual(pd.rel.analysis, 'trt_splits', ...
                'fixture precondition broken: unexpected analysis label');
            testCase.assertFalse(isfield(pd, 'relsummary'), ...
                'fixture precondition broken: fixture already summarized');

            % try/catch rather than verifyError: the identifier is shared
            % by the G42 five-label guard in the same file, so the message
            % pin below is what proves the FAMILY guard fired (the G42
            % message carries no splits sentence -- adding trt_splits to
            % the G42 blacklist instead would fail this pin). It is also
            % the only coverage the user-facing refusal text has.
            caught = [];
            try
                psyrat_criterionfigures('psyrat_prefs', prefs, ...
                    'psyrat_data', pd, 'analysis', 'trt');
            catch err
                caught = err;
            end
            testCase.assertNotEmpty(caught, ...
                'trt_splits must be refused on the trt route');
            testCase.verifyEqual(caught.identifier, ...
                'criterion:unsupportedDesign');
            testCase.verifyTrue(contains(caught.message, ...
                'items-per-split weighting'), ...
                'the refusal must state the splits-specific reason');
            testCase.verifyTrue(contains(caught.message, 'trt_splits'), ...
                'the refusal must name the refused label');

            testCase.verifyError(@() psyrat_criterionfigures( ...
                'psyrat_prefs', prefs, 'psyrat_data', pd, ...
                'analysis', 'sing'), 'criterion:unsupportedDesign', ...
                'the refusal is categorical, so the sing route refuses too');
        end

        function testPlainDesignsPassTheGuard(testCase)
            % Positive control: plain 'ic' and plain 'trt' summary data must
            % clear the guard and complete. Plot and table outputs are
            % switched off so this exercises the guard plus the settings and
            % rebuild plumbing without drawing figures; the full end-to-end
            % assertions (figures drawn, values matching the kernels) are
            % TestSummaryAndPresentationFunctions'.
            prefs = struct('view', struct( ...
                'criterioncutoff', 0.5, ...
                'plotcriterion', 0, ...
                'tablecriterion', 0), 'ver', '0-test');

            singData = PsyRATTestDataFactory.makeICSummaryData();
            psyrat_criterionfigures('psyrat_prefs', prefs, ...
                'psyrat_data', singData, 'analysis', 'sing');

            trtData = PsyRATTestDataFactory.makeTRTSummaryData();
            psyrat_criterionfigures('psyrat_prefs', prefs, ...
                'psyrat_data', trtData, 'analysis', 'trt');

            % Reaching this line is the assertion: neither call errored.
            testCase.verifyTrue(true);
        end

    end
end
