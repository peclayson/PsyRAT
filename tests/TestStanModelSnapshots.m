classdef TestStanModelSnapshots < PsyRATTestBase
    % Byte-identical golden snapshots of the generated Stan model code.
    %
    % This is the safety net for the W6 Stan-model-builder extraction: it
    % captures the exact Stan source that psyrat_computevarcomp generates for
    % each design/sserrvar/difference scenario, in both the non-threaded and
    % within-chain threaded (reduce_sum) forms, and compares it byte-for-byte
    % to committed golden files. Any change to the generated Stan (intended or
    % not) makes the corresponding golden fail, so the planned extraction can be
    % proven behavior-preserving at the source level.
    %
    % Capture works through the PSYRAT_EMIT_MODEL_FILE hook in psyrat_run_stan,
    % which writes the assembled model code and aborts before any CmdStan call.
    % The threaded form is forced deterministically with the
    % PSYRAT_FORCE_AVAILABLE_CORES override (so within-chain threading turns on
    % regardless of the host core count). The test therefore needs NO
    % CmdStan/MatlabStan and runs in the standard unit lane on every push.
    %
    % Goldens live in tests/golden/stan/<scenario>.stan (non-threaded) and
    % tests/golden/stan/<scenario>__threaded.stan (threaded). To (re)generate
    % them after an intended change, run with PSYRAT_REFRESH_STAN_GOLDENS=1,
    % review the diff, and commit the updated files.

    properties (TestParameter)
        Scenario = psyrat_snapshot_scenarios();
    end

    methods (Test)
        function testGeneratedStanMatchesGolden(testCase, Scenario)
            % Scenarios that intentionally fail preflight never generate a
            % model. Assert the preflight guard still fires (asserting rather
            % than skipping keeps the unit lane free of Incomplete results).
            if ~isempty(Scenario.expectedErrorId)
                runScenario = localRunner(testCase, Scenario, false);
                testCase.verifyError(runScenario, Scenario.expectedErrorId);
                return;
            end

            % Non-threaded and within-chain threaded forms are separate twins
            % that the builder extraction will collapse, so both are pinned.
            localCheckGolden(testCase, Scenario, false, '');
            localCheckGolden(testCase, Scenario, true, '__threaded');
        end

        function testNativeKeyMatchesCaseSite(testCase, Scenario)
            % Pin the native design key each analysis's case site actually
            % passes to psyrat_run_stan.
            %
            % WHY THIS EXISTS. TestNativeEngine.testRegistryDispatchConsistency
            % documents a THREE-way contract - the key is spelled at the
            % psyrat_run_stan case site, in psyrat_native_supported, and in the
            % psyrat_native_hmc / psyrat_native_fitlme switches - but enforces
            % only the last two. Nothing read the case site. A mistyped or
            % dropped key there resolves through psyrat_native_supported's
            % "otherwise -> engines = {}" arm, which is INDISTINGUISHABLE from
            % "CmdStan-only by design" (analyses 25-29 and the gamma arms of
            % 19/20 pass '' on purpose). The design silently loses its native
            % engine and every suite stays green. This is the missing leg.
            %
            % WHY RUNTIME CAPTURE, NOT A SOURCE SCAN. Eight of the 29 cases
            % choose the key conditionally - on family (19/20), on the analysis
            % number inside a shared "case {11,26}" / "case {14,27}" arm, on
            % concurrency (7/10), and through two levels of function indirection
            % (8, whose key lives past the switch's own end). A parser cannot
            % resolve any of those. The PSYRAT_EMIT_NATIVE_FILE sidecar reads the
            % key off the fully-resolved spec at run time instead, so what is
            % pinned is what the case site actually chose.
            %
            % Analyses 26 and 27 are reached only by gamma scenarios here. Their
            % key assignment keys on the analysis NUMBER, not the family, so that
            % still exercises the branch.
            if ~isempty(Scenario.expectedErrorId)
                % Never reaches psyrat_run_stan; assert the preflight guard
                % instead of skipping, so the lane stays free of Incomplete.
                testCase.verifyError(localRunner(testCase, Scenario, false), ...
                    Scenario.expectedErrorId);
                return;
            end

            golden = localNativeKeyGolden();
            testCase.assertTrue(golden.isKey(Scenario.name), sprintf( ...
                ['Scenario ''%s'' has no entry in localNativeKeyGolden. Every ' ...
                'snapshot scenario must pin the native key its design routes ' ...
                'to; add a row (with its analysis number) rather than ' ...
                'exempting the scenario.'], Scenario.name));
            expected = golden(Scenario.name);

            % (1) The GOLDEN itself may only claim an empty key where an empty
            % key is legitimate. Without this, a failing row could be "fixed" by
            % blanking its expectation, which is precisely the silent
            % CmdStan-only degradation this test exists to catch.
            if isempty(expected.key)
                testCase.verifyTrue( ...
                    ismember(expected.analysis, localCmdStanOnlyAnalyses()), ...
                    sprintf(['localNativeKeyGolden claims an EMPTY native key ' ...
                    'for scenario ''%s'' (analysis %d), but only analyses %s ' ...
                    'refuse a native fallback by design. An empty key elsewhere ' ...
                    'is a dropped key, not a decision.'], Scenario.name, ...
                    expected.analysis, mat2str(localCmdStanOnlyAnalyses())));
            end

            actual = localCaptureNativeKey(testCase, Scenario);

            % (2) The case site resolved to the pinned key.
            testCase.verifyEqual(actual, expected.key, sprintf( ...
                ['Scenario ''%s'' (analysis %d) passed native key ''%s'' to ' ...
                'psyrat_run_stan; localNativeKeyGolden pins ''%s''. If the ' ...
                'routing changed deliberately, update the golden row; if not, ' ...
                'the case site lost or mistyped its key and this design has ' ...
                'silently become CmdStan-only.'], Scenario.name, ...
                expected.analysis, actual, expected.key));

            % (3) The key the case site passed is one the registry recognizes.
            % This is the typo-catcher: an unregistered key is not an error at
            % run time, it just yields no engine.
            if ~isempty(actual)
                testCase.verifyTrue(ismember(actual, psyrat_native_supported()), ...
                    sprintf(['Scenario ''%s'' (analysis %d) passes native key ' ...
                    '''%s'', which psyrat_native_supported does not enumerate. ' ...
                    'An unrecognized key returns no engine, so this design would ' ...
                    'fall back to CmdStan with no error anywhere.'], ...
                    Scenario.name, expected.analysis, actual));
            end
        end

        function testNativeKeySidecarIsOptIn(testCase)
            % NEGATIVE CONTROL for the sidecar's opt-in guard, run as a
            % SAME-PATH differential.
            %
            % psyrat_write_emitted_model must write the native-key file ONLY when
            % PSYRAT_EMIT_NATIVE_FILE is set. Four other test classes drive the
            % PSYRAT_EMIT_MODEL_FILE hook (TestSplitsParallel, TestGammaScaleOption,
            % TestDispersionDefault, TestDodConcurrentGuard) and must be left
            % untouched, with no stray file and no new error.
            %
            % WHY THE SAME PATH, TWICE. An earlier version of this test watched a
            % fresh tempname while the variable was unset and asserted no file
            % appeared there. That cannot fail: nothing would ever write to a path
            % nobody was given. Verified by deleting the guard - the test still
            % passed. The assertion only means something if the SAME path is
            % written in one arm and then, after deletion, stays absent in the
            % other. Arm 1 is the positive control that the path and the run are
            % capable of producing the file at all.
            cfg = psyrat_snapshot_scenarios().ic_base;
            modelFile  = [tempname '.stan'];
            nativeFile = [tempname '.nativekey'];
            restoreModel  = onCleanup(@() setenv('PSYRAT_EMIT_MODEL_FILE', ''));  %#ok<NASGU>
            restoreNative = onCleanup(@() setenv('PSYRAT_EMIT_NATIVE_FILE', '')); %#ok<NASGU>
            cleanupModel  = onCleanup(@() localDeleteIfExists(modelFile));  %#ok<NASGU>
            cleanupNative = onCleanup(@() localDeleteIfExists(nativeFile)); %#ok<NASGU>

            % --- arm 1: variable SET -> the sidecar appears at this path
            setenv('PSYRAT_EMIT_MODEL_FILE', modelFile);
            setenv('PSYRAT_EMIT_NATIVE_FILE', nativeFile);
            testCase.verifyError(localRunner(testCase, cfg, false), ...
                'PsyRAT:emitModelOnly');
            testCase.assertTrue(isfile(nativeFile), ...
                ['Sidecar was not written even with PSYRAT_EMIT_NATIVE_FILE set, ' ...
                'so the negative arm below would prove nothing.']);

            % --- arm 2: same path, variable UNSET -> it must not come back
            delete(nativeFile);
            setenv('PSYRAT_EMIT_NATIVE_FILE', '');
            testCase.verifyError(localRunner(testCase, cfg, false), ...
                'PsyRAT:emitModelOnly');
            testCase.verifyTrue(isfile(modelFile), ...
                'Model emit did not run in arm 2, so its absence proves nothing.');
            testCase.verifyFalse(isfile(nativeFile), ...
                ['The native-key sidecar reappeared with PSYRAT_EMIT_NATIVE_FILE ' ...
                'unset. The opt-in guard in psyrat_write_emitted_model is gone, ' ...
                'and the other PSYRAT_EMIT_MODEL_FILE users now leave stray files.']);
        end

        function testDelegatingGammaBuildersMatchTheirGroupLevelTwin(testCase)
            % Several SUBJECT-LEVEL gamma designs reuse their group-level twin's
            % Stan model verbatim: the person-varying log-nu already gives each
            % participant their own dispersion, so no location-scale residual
            % submodel is needed and the subject-level distinction is purely a
            % downstream extraction. Their builders are one-line delegations
            % (psyrat_build_stan_sserr_chisq, _diff_sserr_chisq,
            % _trt_sserr_chisq), which exist only to keep a distinct golden and
            % clean per-analysis provenance.
            %
            % Each design's own golden is pinned above; this test additionally
            % pins the DELEGATION itself, so that a future edit which quietly
            % gives one of these designs a divergent model fails here with a clear
            % reason rather than silently rewriting a golden on the next refresh.
            % A location-scale delegation is pinned wherever BOTH members of the
            % pair have a gammascale=2 variant. ic_gamma_ls_sserr (case 6) and
            % ic_gamma_ls_trt_sserr (case 25) are deliberately absent, because
            % their group-level twins (analyses 1 and 5) have no location-scale
            % variant to match - those are single-mean designs whose coefficients
            % the parameterization leaves unchanged. Cases 7/8 differ: the
            % subject-level builder forwards gammascale to the group-level one, so
            % the delegation identity survives under BOTH parameterizations. The
            % dynamic designs (26<-11, 27<-14) do too, for a different reason -
            % they share a case block and a builder outright, under either
            % parameterization.
            %
            % The 26<-11 location-scale row was MISSING until 2026-08-05: the two
            % goldens were byte-identical from the day increment 4 landed them, but
            % nothing asserted it, so a future edit could have silently given one
            % of them a divergent model. Added alongside its two-facet counterpart.
            pairs = { ...
                'ic_gamma_sserr',                    'ic_gamma_base'; ...               % case 6  <- case 1
                'ic_gamma_diff_sserr',               'ic_gamma_diff_event'; ...         % case 8  <- case 7
                'ic_gamma_ls_diff_sserr',            'ic_gamma_ls_diff_event'; ...      % case 8 LS <- case 7 LS
                'ic_gamma_diff_sserr_copula',        'ic_gamma_diff_copula'; ...        % case 8 concurrent <- case 7 concurrent
                'ic_gamma_trt_sserr',                'ic_gamma_trt'; ...                % case 25 <- case 5
                'dynrel_gamma_sserr_onefacet',       'dynrel_gamma_onefacet'; ...       % case 26 <- case 11
                'dynrel_gamma_ls_sserr_onefacet',    'dynrel_gamma_ls_onefacet'; ...    % case 26 LS <- case 11 LS
                'dynrel_gamma_sserr_trt_onefacet',   'dynrel_gamma_trt_onefacet'; ...   % case 27 <- case 14
                'dynrel_gamma_ls_sserr_trt_onefacet','dynrel_gamma_ls_trt_onefacet'; ...% case 27 LS <- case 14 LS
                'dynrel_gamma_ls_diff_sserr_onefacet','dynrel_gamma_ls_diff_onefacet'};% case 13 LS <- case 12 LS

            goldenDir = fullfile(testCase.projectRoot(), 'tests', 'golden', 'stan');
            for k = 1:size(pairs, 1)
                for suffix = {'', '__threaded'}
                    sub = fullfile(goldenDir, [pairs{k,1} suffix{1} '.stan']);
                    grp = fullfile(goldenDir, [pairs{k,2} suffix{1} '.stan']);
                    testCase.assertEqual(exist(sub, 'file'), 2, ...
                        sprintf('Missing golden: %s', sub));
                    testCase.assertEqual(exist(grp, 'file'), 2, ...
                        sprintf('Missing golden: %s', grp));
                    testCase.verifyEqual(fileread(sub), fileread(grp), ...
                        sprintf(['Golden %s%s must be BYTE-IDENTICAL to its ' ...
                        'group-level twin %s%s. These pairs share their emitted ' ...
                        'model for one of THREE reasons, and which one decides ' ...
                        'where to look: the STATIC rows (cases 6, 8, 25) have a ' ...
                        'subject-level builder that is a one-line delegation, so ' ...
                        'annotate that builder''s header; the NON-DIFFERENCE ' ...
                        'DYNAMIC rows (cases 26, 27, in both parameterizations) ' ...
                        'have no subject-level builder at all - they share a case ' ...
                        'block AND a builder outright with cases 11/14, so ' ...
                        'annotate the group-level builder; and case 13 LS is a ' ...
                        'THIRD kind - a separate case block with its own ' ...
                        'delegating builder (psyrat_build_stan_diffdynrel_sserr_' ...
                        'chisq), because unlike the Gaussian pair the gamma ' ...
                        'case-12 model is ALREADY subject-level-capable (its ' ...
                        'person block is already 4-wide with the two events'' ' ...
                        'scale effects), so there is nothing to widen. If the ' ...
                        'model genuinely ' ...
                        'diverged, update this pair list and say why there.'], ...
                        pairs{k,1}, suffix{1}, pairs{k,2}, suffix{1}));
                end
            end
        end
    end
end

function localCheckGolden(testCase, cfg, threaded, suffix)
actualStan = localCaptureGeneratedStan(testCase, cfg, threaded);

goldenPath = fullfile(testCase.projectRoot(), 'tests', 'golden', 'stan', ...
    [cfg.name suffix '.stan']);

if localRefreshRequested()
    localEnsureDir(fileparts(goldenPath));
    fid = fopen(goldenPath, 'w');
    testCase.assertGreaterThanOrEqual(fid, 0, ...
        sprintf('Could not open golden file for writing: %s', goldenPath));
    closeFid = onCleanup(@() fclose(fid));
    fwrite(fid, actualStan, 'char');
    clear closeFid;
    fprintf('Refreshed Stan golden: %s\n', goldenPath);
    return;
end

testCase.assertTrue(isfile(goldenPath), ...
    sprintf(['Missing Stan golden for scenario ''%s''%s (%s). How to fix: ', ...
    'run with PSYRAT_REFRESH_STAN_GOLDENS=1 to capture goldens, review, commit.'], ...
    cfg.name, suffix, goldenPath));

expectedStan = fileread(goldenPath);
testCase.verifyEqual(actualStan, expectedStan, ...
    sprintf(['Generated Stan for scenario ''%s''%s differs from its committed ', ...
    'golden. If the change is intended, refresh with ', ...
    'PSYRAT_REFRESH_STAN_GOLDENS=1 and review the diff.'], cfg.name, suffix));
end

function stanText = localCaptureGeneratedStan(testCase, cfg, threaded)
% Drive psyrat_computevarcomp for one scenario with the emit hook active and
% return the generated Stan source. The run aborts before CmdStan via the
% PSYRAT_EMIT_MODEL_FILE hook (PsyRAT:emitModelOnly). When threaded is true,
% within-chain threading is forced on with PSYRAT_FORCE_AVAILABLE_CORES.
emitFile = [tempname '.stan'];
setenv('PSYRAT_EMIT_MODEL_FILE', emitFile);
restoreEmit = onCleanup(@() setenv('PSYRAT_EMIT_MODEL_FILE', '')); %#ok<NASGU>

if threaded
    setenv('PSYRAT_FORCE_AVAILABLE_CORES', '8');
    restoreCores = onCleanup(@() setenv('PSYRAT_FORCE_AVAILABLE_CORES', '')); %#ok<NASGU>
end
cleanupFile = onCleanup(@() localDeleteIfExists(emitFile)); %#ok<NASGU>

testCase.verifyError(localRunner(testCase, cfg, threaded), 'PsyRAT:emitModelOnly');

testCase.assertTrue(isfile(emitFile), ...
    sprintf('Emit hook did not write a model file for scenario ''%s'' (threaded=%d).', ...
    cfg.name, threaded));
stanText = fileread(emitFile);
end

function nativeKey = localCaptureNativeKey(testCase, cfg)
% Drive psyrat_computevarcomp for one scenario and return the native design key
% its case site passed to psyrat_run_stan.
%
% Same mechanism as localCaptureGeneratedStan: the run aborts before any CmdStan
% interaction via PSYRAT_EMIT_MODEL_FILE (PsyRAT:emitModelOnly). BOTH variables
% are set because the model variable is what fires the hook, while
% PSYRAT_EMIT_NATIVE_FILE is what makes psyrat_write_emitted_model also record
% the key. The model file is a throwaway here.
%
% Non-threaded only. Within-chain threading changes the generated Stan but not
% the native key, so capturing the threaded twin would double the cost of this
% test for no additional coverage.
modelFile  = [tempname '.stan'];
nativeFile = [tempname '.nativekey'];
setenv('PSYRAT_EMIT_MODEL_FILE', modelFile);
setenv('PSYRAT_EMIT_NATIVE_FILE', nativeFile);
restoreModel  = onCleanup(@() setenv('PSYRAT_EMIT_MODEL_FILE', ''));  %#ok<NASGU>
restoreNative = onCleanup(@() setenv('PSYRAT_EMIT_NATIVE_FILE', '')); %#ok<NASGU>
cleanupModel  = onCleanup(@() localDeleteIfExists(modelFile));  %#ok<NASGU>
cleanupNative = onCleanup(@() localDeleteIfExists(nativeFile)); %#ok<NASGU>

testCase.verifyError(localRunner(testCase, cfg, false), 'PsyRAT:emitModelOnly');

testCase.assertTrue(isfile(nativeFile), sprintf( ...
    ['Native-key sidecar was not written for scenario ''%s''. The run must ' ...
    'reach psyrat_run_stan for the key to be observable.'], cfg.name));

% The sidecar is written even when the key is EMPTY, which is a meaningful
% value here, so an empty file means "this case refuses a native fallback",
% not "the hook never ran" - the assert above separates those two.
nativeKey = fileread(nativeFile);
if isempty(nativeKey)
    nativeKey = '';   % normalize 1x0 vs 0x0 char so verifyEqual compares content
end
end

function analyses = localCmdStanOnlyAnalyses()
% Analyses that pass an EMPTY native key deliberately, refusing a native
% fallback rather than silently fitting a different model.
%
% 25, 26, 27 are CmdStan-only "for now" (no native port of their model). 28/29
% null the key in BOTH families: the gamma arm is pinned to CmdStan upstream,
% and the Gaussian arm is CmdStan-only by the 2026-08-20 owner ruling (no
% native port of the four-margin model exists; a native HMC port is a
% separately scoped follow-up). 19/20 are listed because their GAMMA arm
% nulls the key for the same reason; their Gaussian arm passes a real key and is
% pinned to it by the golden map.
analyses = [19 20 25 26 27 28 29];
end

function m = localNativeKeyGolden()
% Golden map: snapshot scenario -> the native design key its case site resolves
% to, with the analysis number that scenario routes to.
%
% Keyed by SCENARIO, not by analysis, because several analyses resolve to
% different keys depending on the run: 7 and 10 switch on concurrency
% (diff/diff_rescor, diff_trt/diff_trt_rescor) and 19/20 null their key under
% the gamma family. Both arms of every such conditional are represented below,
% which is what makes this a pin rather than a sample.
%
% The keys were derived from the case sites independently of this test, then
% compared against what the run actually emitted. To regenerate after an
% intended routing change, read the failure message: it names the scenario, its
% analysis, and both keys.
pairs = { ...
    ...% analysis 1
    'ic_base',                                         1, 'sing';                            ...
    'ic_base_priors',                                  1, 'sing';                            ...
    'ic_gamma_base',                                   1, 'sing';                            ...
    ...% analysis 2
    'ic_gamma_group',                                  2, 'sing_facet';                      ...
    'ic_group',                                        2, 'sing_facet';                      ...
    ...% analysis 3
    'ic_event',                                        3, 'sing_facet';                      ...
    'ic_event_priors',                                 3, 'sing_facet';                      ...
    'ic_gamma_event',                                  3, 'sing_facet';                      ...
    ...% analysis 4
    'ic_group_event',                                  4, 'sing_facet';                      ...
    ...% analysis 5
    'ic_gamma_trt',                                    5, 'trt';                             ...
    'trt_event_time',                                  5, 'trt';                             ...
    'trt_group_event_time',                            5, 'trt';                             ...
    'trt_group_time',                                  5, 'trt';                             ...
    'trt_time',                                        5, 'trt';                             ...
    'trt_time_priors',                                 5, 'trt';                             ...
    ...% analysis 6
    'ic_gamma_ls_sserr',                               6, 'sserr';                           ...
    'ic_gamma_sserr',                                  6, 'sserr';                           ...
    'ic_sserr_base',                                   6, 'sserr';                           ...
    'ic_sserr_event',                                  6, 'sserr';                           ...
    'ic_sserr_group',                                  6, 'sserr';                           ...
    'ic_sserr_group_event',                            6, 'sserr';                           ...
    ...% analysis 7
    'ic_diff_event',                                   7, 'diff';                            ...
    'ic_diff_event_priors',                            7, 'diff';                            ...
    'ic_diff_group_event',                             7, 'diff';                            ...
    'ic_diff_rescor_event',                            7, 'diff_rescor';                     ...
    'ic_gamma_diff_copula',                            7, 'diff_rescor';                     ...
    'ic_gamma_diff_copula_fixnu',                      7, 'diff_rescor';                     ...
    'ic_gamma_diff_event',                             7, 'diff';                            ...
    'ic_gamma_ls_diff_event',                          7, 'diff';                            ...
    ...% analysis 8
    'ic_diff_group_event_sserr',                       8, 'diff_sserr';                      ...
    'ic_diff_group_event_sserr_priors',                8, 'diff_sserr';                      ...
    'ic_gamma_diff_sserr',                             8, 'diff_sserr';                      ...
    'ic_gamma_diff_sserr_copula',                      8, 'diff_sserr';                      ...
    'ic_gamma_ls_diff_sserr',                          8, 'diff_sserr';                      ...
    ...% analysis 9
    'dod_event',                                       9, 'dod';                             ...
    'dod_event_priors',                                9, 'dod';                             ...
    'ic_gamma_dodiff',                                 9, 'dod';                             ...
    'ic_gamma_ls_dodiff',                              9, 'dod';                             ...
    ...% analysis 10
    'ic_gamma_diff_trt',                              10, 'diff_trt';                        ...
    'ic_gamma_diff_trt_copula',                       10, 'diff_trt_rescor';                 ...
    'ic_gamma_diff_trt_copula_fixnu',                 10, 'diff_trt_rescor';                 ...
    'ic_gamma_ls_diff_trt',                           10, 'diff_trt';                        ...
    'trt_diff_event_time',                            10, 'diff_trt';                        ...
    'trt_diff_event_time_priors',                     10, 'diff_trt';                        ...
    'trt_diff_rescor_event_time',                     10, 'diff_trt_rescor';                 ...
    ...% analysis 11
    'dynrel_gamma_ls_onefacet',                       11, 'dynrel';                          ...
    'dynrel_gamma_onefacet',                          11, 'dynrel';                          ...
    'dynrel_onefacet',                                11, 'dynrel';                          ...
    ...% analysis 12
    'dynrel_diff_onefacet',                           12, 'diff_dynrel';                     ...
    'dynrel_gamma_diff_onefacet',                     12, 'diff_dynrel';                     ...
    'dynrel_gamma_ls_diff_onefacet',                  12, 'diff_dynrel';                     ...
    ...% analysis 13
    'dynrel_diff_sserr_onefacet',                     13, 'diff_dynrel_sserr';               ...
    'dynrel_gamma_ls_diff_sserr_onefacet',            13, 'diff_dynrel_sserr';               ...
    ...% analysis 14
    'dynrel_gamma_ls_trt_onefacet',                   14, 'dynrel_trt';                      ...
    'dynrel_gamma_trt_onefacet',                      14, 'dynrel_trt';                      ...
    'dynrel_trt_onefacet',                            14, 'dynrel_trt';                      ...
    ...% analysis 15
    'dynrel_diff_trt_onefacet',                       15, 'diff_dynrel_trt';                 ...
    ...% analysis 16
    'dynrel_diff_sserr_trt_onefacet',                 16, 'diff_dynrel_sserr_trt';           ...
    ...% analysis 17
    'dynrel_diff_rescor_onefacet',                    17, 'diff_dynrel_rescor';              ...
    ...% analysis 18
    'dynrel_diff_rescor_trt_onefacet',                18, 'diff_dynrel_trt_rescor';          ...
    ...% analysis 19
    'dynrel_diff_sserr_rescor_onefacet',              19, 'diff_dynrel_sserr_rescor';        ...
    'dynrel_gamma_ls_diff_sserr_rescor_onefacet',     19, '';                                ...
    ...% analysis 20
    'dynrel_diff_sserr_rescor_trt_onefacet',          20, 'diff_dynrel_sserr_trt_rescor';    ...
    'dynrel_gamma_ls_diff_sserr_rescor_trt_onefacet', 20, '';                                ...
    ...% analysis 21
    'dynrel_diff_sserr_ssrescor_onefacet',            21, 'diff_dynrel_sserr_rho';           ...
    ...% analysis 22
    'dynrel_diff_sserr_ssrescor_trt_onefacet',        22, 'diff_dynrel_sserr_trt_rho';       ...
    ...% analysis 23
    'splits_iins_single',                             23, 'splits';                          ...
    ...% analysis 24
    'splits_iins_trt',                                24, 'splits_trt';                      ...
    ...% analysis 25
    'ic_gamma_ls_trt_sserr',                          25, '';                                ...
    'ic_gamma_trt_sserr',                             25, '';                                ...
    'trt_group_event_time_sserr',                     25, '';                                ...
    ...% analysis 26
    'dynrel_gamma_ls_sserr_onefacet',                 26, '';                                ...
    'dynrel_gamma_sserr_onefacet',                    26, '';                                ...
    ...% analysis 27
    'dynrel_gamma_ls_sserr_trt_onefacet',             27, '';                                ...
    'dynrel_gamma_sserr_trt_onefacet',                27, '';                                ...
    ...% analysis 28
    'dynrel_dodiff_sserr_onefacet',                   28, '';                                ...
    'dynrel_gamma_ls_dodiff_sserr_onefacet',          28, '';                                ...
    ...% analysis 29
    'dynrel_dodiff_sserr_trt_onefacet',               29, '';                                ...
    'dynrel_gamma_ls_dodiff_sserr_trt_onefacet',      29, '';                                ...
    };

m = containers.Map('KeyType','char','ValueType','any');
for k = 1:size(pairs,1)
    m(pairs{k,1}) = struct('analysis',pairs{k,2},'key',pairs{k,3});
end
end

function runner = localRunner(testCase, cfg, threaded)
% Build the data table and the psyrat_computevarcomp argument list for one
% scenario, returning a no-arg function handle that runs it.
if isfield(cfg, 'isDoDDynrelTrt') && cfg.isDoDDynrelTrt
    % two-facet person-specific dynamic DoD (analysis 29): the four-event
    % dimension table with a two-occasion time column. Checked FIRST with its
    % one-facet sibling below: both also set isDynrel/isDoD for the argument
    % assembly, so the generic dynrel/DoD table branches must not shadow them.
    dataTbl = localBuildDoDDynrelTrtTable();
elseif isfield(cfg, 'isDoDDynrel') && cfg.isDoDDynrel
    % person-specific dynamic NONCONCURRENT difference-of-differences (analysis
    % 28): four-event long-format table with a dimension column and deliberately
    % UNEQUAL per-cell trial counts and amplitudes (equal cells would not
    % exercise the cell contrasts, and equal counts would not exercise the
    % long-format property the design exists for).
    dataTbl = localBuildDoDDynrelTable();
elseif isfield(cfg, 'isDynrelDiffSsRhoTrt') && cfg.isDynrelDiffSsRhoTrt
    % two-facet SUBJECT-LEVEL CONCURRENT difference with a PER-SUBJECT residual
    % correlation (analysis 22, Stage 4e): same balanced two-event + occasion +
    % dimension table; sserrvar=2 + diffwpcov=2 + ssrescor=2 routes to 22.
    dataTbl = localBuildDynrelDiffTrtTable();
elseif isfield(cfg, 'isDynrelDiffSsRho') && cfg.isDynrelDiffSsRho
    % one-facet SUBJECT-LEVEL CONCURRENT difference with a PER-SUBJECT residual
    % correlation (analysis 21, Stage 4e): balanced two-event dimension table;
    % sserrvar=2 + diffwpcov=2 + ssrescor=2 routes to 21.
    dataTbl = localBuildDynrelDiffTable();
elseif isfield(cfg, 'isDynrelDiffSserrRescorTrt') && cfg.isDynrelDiffSserrRescorTrt
    % two-facet SUBJECT-LEVEL CONCURRENT difference (analysis 20): balanced
    % two-event + occasion + dimension table; sserrvar=2 + diffwpcov=2 routes to 20.
    dataTbl = localBuildDynrelDiffTrtTable();
elseif isfield(cfg, 'isDynrelDiffSserrRescor') && cfg.isDynrelDiffSserrRescor
    % one-facet SUBJECT-LEVEL CONCURRENT difference (analysis 19): balanced
    % two-event dimension table; sserrvar=2 + diffwpcov=2 routes to 19.
    dataTbl = localBuildDynrelDiffTable();
elseif isfield(cfg, 'isDynrelDiffRescorTrt') && cfg.isDynrelDiffRescorTrt
    % two-facet group-level CONCURRENT difference (analysis 18): same balanced
    % two-event + occasion + dimension table; diffwpcov=2 routes it to 18.
    dataTbl = localBuildDynrelDiffTrtTable();
elseif isfield(cfg, 'isDynrelDiffRescor') && cfg.isDynrelDiffRescor
    % one-facet group-level CONCURRENT difference (analysis 17): the balanced
    % two-event dimension table (equal event counts pair); diffwpcov=2 routes to 17.
    dataTbl = localBuildDynrelDiffTable();
elseif isfield(cfg, 'isDynrelDiffSserrTrt') && cfg.isDynrelDiffSserrTrt
    % subject-level two-facet difference (analysis 16): same two-event + occasion
    % + dimension table as the group-level case 15; sserrvar=2 routes it to 16.
    dataTbl = localBuildDynrelDiffTrtTable();
elseif isfield(cfg, 'isDynrelDiffTrt') && cfg.isDynrelDiffTrt
    dataTbl = localBuildDynrelDiffTrtTable();
elseif isfield(cfg, 'isDynrelTrt') && cfg.isDynrelTrt
    dataTbl = localBuildDynrelTrtTable();
elseif isfield(cfg, 'isDynrelDiff') && cfg.isDynrelDiff
    dataTbl = localBuildDynrelDiffTable();
elseif isfield(cfg, 'isDynrel') && cfg.isDynrel
    dataTbl = localBuildDynrelTable();
elseif isfield(cfg, 'isSplitsTrt') && cfg.isSplitsTrt
    % nonparallel data-splits, test-retest (analysis 24): id/meas/weight + time
    dataTbl = localBuildSplitsTrtTable();
elseif isfield(cfg, 'isSplits') && cfg.isSplits
    % nonparallel data-splits, single occasion (analysis 23): id/meas/weight
    dataTbl = localBuildSplitsTable();
elseif isfield(cfg, 'isDoD') && cfg.isDoD
    dataTbl = localBuildDoDTable();
elseif isfield(cfg, 'isGamma') && cfg.isGamma
    % Gamma / scaled-chi-square one-facet family: reuse the one-facet table but
    % force strictly positive scores (the family's support). The generated Stan
    % text does not depend on the data values, so the shift is inconsequential
    % to the golden while keeping the fixture coherent with the family.
    psyrat_assume_dev_testdata(testCase);
    dataTbl = PsyRATCmdStanIntegrationFactory.makeScenarioTable( ...
        testCase.projectRoot(), cfg);
    dataTbl.meas = dataTbl.meas - min(dataTbl.meas) + 1;
else
    psyrat_assume_dev_testdata(testCase);
    dataTbl = PsyRATCmdStanIntegrationFactory.makeScenarioTable( ...
        testCase.projectRoot(), cfg);
end

% Gamma / scaled-chi-square dynrel scenarios reuse the dimension-bearing dynrel
% tables (isDynrel*) whose deterministic meas can be non-positive; the gamma
% family requires strictly positive scores (its support), so shift them exactly
% as the standalone isGamma branch does. The generated Stan is data-value
% independent, so the shift is inconsequential to the golden.
if isfield(cfg, 'family') && strcmp(cfg.family, 'gamma') && any(dataTbl.meas <= 0)
    dataTbl.meas = dataTbl.meas - min(dataTbl.meas) + 1;
end

args = {'data', dataTbl, 'chains', 1, 'warmup', 1, 'sampling', 1, ...
    'seed', 12345, 'showgui', 0, 'verbose', 1, ...
    'sserrvar', cfg.sserrvar, 'diffest', cfg.diffest, ...
    'diffwpcov', cfg.diffwpcov, 'diffrescor', cfg.diffrescor};

if isfield(cfg, 'family') && ~isempty(cfg.family)
    args = [args, {'family', cfg.family}];
end

if isfield(cfg, 'isDynrel') && cfg.isDynrel
    args = [args, {'dynrel', 2}];
end
if isfield(cfg, 'isSplits') && cfg.isSplits
    args = [args, {'splits', cfg.splits}];
end
if isfield(cfg, 'ssrescor') && ~isempty(cfg.ssrescor)
    args = [args, {'ssrescor', cfg.ssrescor}];
end
if isfield(cfg, 'dispersion') && ~isempty(cfg.dispersion)
    args = [args, {'dispersion', cfg.dispersion}];
end
if isfield(cfg, 'gammascale') && ~isempty(cfg.gammascale)
    args = [args, {'gammascale', cfg.gammascale}];
end
if threaded
    args = [args, {'withinchain', 2}];
end
if isfield(cfg, 'isDoD') && cfg.isDoD
    args = [args, {'dodmap', cfg.dodmap}];
end
if isfield(cfg, 'priors') && ~isempty(cfg.priors)
    args = [args, {'priors', cfg.priors}];
end

runner = @() psyrat_computevarcomp(args{:});
end

function tbl = localBuildDoDTable()
% Deterministic four-event long-format table for the difference-of-differences
% scenario: 4 participants x 4 events x 6 trials, numeric measurements.
events = {'E1', 'E2', 'E3', 'E4'};
ids = {};
evs = {};
meas = [];
for s = 1:4
    for e = 1:4
        for t = 1:6
            ids{end+1, 1} = sprintf('s%d', s); %#ok<AGROW>
            evs{end+1, 1} = events{e}; %#ok<AGROW>
            meas(end+1, 1) = 0.1 * t + 0.01 * e; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, evs, 'VariableNames', {'id', 'meas', 'event'});
end

function tbl = localBuildDoDDynrelTable()
% Deterministic four-event long-format table for the person-specific dynamic
% difference-of-differences scenario (analysis 28): 5 participants x 4 events
% with a subject-level dimension column, STRICTLY POSITIVE measurements (the
% gamma family's support), deliberately UNEQUAL per-cell trial counts
% ([6 8 5 7] - the long format must accept them; the paired layout errors) and
% unequal per-cell amplitudes (equal cell means would not exercise the cell
% contrasts).
events = {'E1', 'E2', 'E3', 'E4'};
ntrl = [6 8 5 7];
amp = [5.0 4.2 5.4 4.6];
nsub = 5;
ids = {};
evs = {};
meas = [];
dim1 = [];
for s = 1:nsub
    dval = -1 + 0.45 * (s - 1);   %distinct standardizable dimension values
    for e = 1:4
        for t = 1:ntrl(e)
            ids{end+1, 1} = sprintf('s%d', s); %#ok<AGROW>
            evs{end+1, 1} = events{e}; %#ok<AGROW>
            meas(end+1, 1) = amp(e) + 0.1 * t + 0.05 * s; %#ok<AGROW>
            dim1(end+1, 1) = dval; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, evs, dim1, ...
    'VariableNames', {'id', 'meas', 'event', 'dim1'});
end

function tbl = localBuildDoDDynrelTrtTable()
% Two-facet twin of localBuildDoDDynrelTable (analysis 29): the same
% four-event + dimension structure repeated over two occasions, with a small
% occasion shift on the measurements so the occasion facet is exercised.
events = {'E1', 'E2', 'E3', 'E4'};
ntrl = [6 8 5 7];
amp = [5.0 4.2 5.4 4.6];
nsub = 5;
ids = {};
evs = {};
tims = {};
meas = [];
dim1 = [];
for s = 1:nsub
    dval = -1 + 0.45 * (s - 1);
    for o = 1:2
        for e = 1:4
            for t = 1:ntrl(e)
                ids{end+1, 1} = sprintf('s%d', s); %#ok<AGROW>
                evs{end+1, 1} = events{e}; %#ok<AGROW>
                tims{end+1, 1} = sprintf('time%d', o); %#ok<AGROW>
                meas(end+1, 1) = amp(e) + 0.1 * t + 0.05 * s + 0.15 * (o - 1); %#ok<AGROW>
                dim1(end+1, 1) = dval; %#ok<AGROW>
            end
        end
    end
end
tbl = table(ids, meas, evs, tims, dim1, ...
    'VariableNames', {'id', 'meas', 'event', 'time', 'dim1'});
end

function tbl = localBuildSplitsTable()
% Deterministic id/meas/weight long-format table for the nonparallel single-
% occasion split-mean scenario (analysis 23): 6 participants x 5 splits, numeric
% split-mean measurements, plus a REQUIRED items-per-split weight column (n_i)
% with real variation. The variation is load-bearing: the splits=3 preflight
% (localValidateSplits) errors on all-equal n_i, and that preflight runs before
% the Stan emit hook, so an all-equal weight would abort capture early. Weights
% are integers >= 1, matching the Stan vector<lower=1> declaration.
nsub = 6; ksplit = 5;
ids = {}; meas = []; weight = []; w = 0;
for s = 1:nsub
    for k = 1:ksplit
        ids{end+1, 1} = sprintf('s%d', s); %#ok<AGROW>
        meas(end+1, 1) = 0.1 * k + 0.5 * s; %#ok<AGROW>
        weight(end+1, 1) = 8 + 4 * mod(w, 7); %#ok<AGROW> varying n_i in [8, 32]
        w = w + 1;
    end
end
tbl = table(ids, meas, weight, 'VariableNames', {'id', 'meas', 'weight'});
end

function tbl = localBuildSplitsTrtTable()
% Deterministic id/meas/weight/time long-format table for the nonparallel
% test-retest split-mean scenario (analysis 24): 6 participants x 3 occasions x 4
% splits, numeric split-mean measurements, the same varying items-per-split
% weight column (n_i; see localBuildSplitsTable), and a time/occasion column
% (ntime>0 routes splits=3 to the test-retest split model, case 24).
nsub = 6; nocc = 3; ksplit = 4;
ids = {}; meas = []; weight = []; tm = {}; w = 0;
for s = 1:nsub
    for o = 1:nocc
        for k = 1:ksplit
            ids{end+1, 1} = sprintf('s%d', s); %#ok<AGROW>
            tm{end+1, 1} = sprintf('t%d', o);  %#ok<AGROW>
            meas(end+1, 1) = 0.1 * k + 0.3 * o + 0.5 * s; %#ok<AGROW>
            weight(end+1, 1) = 8 + 4 * mod(w, 7); %#ok<AGROW> varying n_i in [8, 32]
            w = w + 1;
        end
    end
end
tbl = table(ids, meas, weight, tm, ...
    'VariableNames', {'id', 'meas', 'weight', 'time'});
end

function tbl = localBuildDynrelTable()
% Deterministic long-format table for the dynamic-reliability scenario: 20
% participants x 8 trials, numeric measurements, plus a between-person
% dimension column (one constant value per participant).
ids = {}; meas = []; dim1 = [];
for s = 1:20
    sval = -1 + (s - 1) * (2 / 19); % subject-level covariate spanning [-1, 1]
    for t = 1:8
        ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
        meas(end+1, 1) = 0.1 * t + 0.5 * sval; %#ok<AGROW>
        dim1(end+1, 1) = sval; %#ok<AGROW>
    end
end
tbl = table(ids, meas, dim1, 'VariableNames', {'id', 'meas', 'dim1'});
end

function tbl = localBuildDynrelTrtTable()
% Deterministic long-format table for the trial + occasion two-facet
% dynamic-reliability scenario: 16 participants x 3 occasions x 5 trials, numeric
% measurements, plus a between-person dimension column (one constant value per
% participant) and a time/occasion column (so ntime > 0 routes to analysis 14).
ids = {}; meas = []; tm = {}; dim1 = [];
for s = 1:16
    sval = -1 + (s - 1) * (2 / 15); % subject-level covariate spanning [-1, 1]
    for o = 1:3
        for t = 1:5
            ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
            tm{end+1, 1} = sprintf('t%d', o);    %#ok<AGROW>
            meas(end+1, 1) = 0.1 * t + 0.3 * o + 0.5 * sval; %#ok<AGROW>
            dim1(end+1, 1) = sval; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, tm, dim1, 'VariableNames', {'id', 'meas', 'time', 'dim1'});
end

function tbl = localBuildDynrelDiffTrtTable()
% Deterministic long-format table for the trial + occasion two-facet group-level
% dynamic difference scenario: 12 participants x 2 events x 2 occasions x 4
% trials, numeric measurements, a between-person dimension column (one constant
% value per participant), and a time/occasion column (so ntime > 0 and the two
% events route to analysis 15).
events = {'gain', 'loss'};
ids = {}; meas = []; evs = {}; tm = {}; dim1 = [];
for s = 1:12
    sval = -1 + (s - 1) * (2 / 11); % subject-level covariate spanning [-1, 1]
    for e = 1:2
        for o = 1:2
            for t = 1:4
                ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
                evs{end+1, 1} = events{e};           %#ok<AGROW>
                tm{end+1, 1} = sprintf('t%d', o);    %#ok<AGROW>
                meas(end+1, 1) = 0.1 * t + 0.3 * o + 0.5 * sval + (e == 1); %#ok<AGROW>
                dim1(end+1, 1) = sval; %#ok<AGROW>
            end
        end
    end
end
tbl = table(ids, meas, evs, tm, dim1, ...
    'VariableNames', {'id', 'meas', 'event', 'time', 'dim1'});
end

function tbl = localBuildDynrelDiffTable()
% Deterministic two-event long-format table for the dynamic difference-score
% scenario: 20 participants x 2 events x 8 trials, numeric measurements, plus a
% between-person dimension column (one constant value per participant).
events = {'gain', 'loss'};
ids = {}; meas = []; evs = {}; dim1 = [];
for s = 1:20
    sval = -1 + (s - 1) * (2 / 19); % subject-level covariate spanning [-1, 1]
    for e = 1:2
        for t = 1:8
            ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
            evs{end+1, 1} = events{e}; %#ok<AGROW>
            meas(end+1, 1) = 0.1 * t + 0.5 * sval + 0.3 * (e == 1); %#ok<AGROW>
            dim1(end+1, 1) = sval; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, evs, dim1, 'VariableNames', {'id', 'meas', 'event', 'dim1'});
end

function tf = localRefreshRequested()
tf = strcmp(getenv('PSYRAT_REFRESH_STAN_GOLDENS'), '1');
end

function localEnsureDir(d)
if ~isempty(d) && ~isfolder(d)
    mkdir(d);
end
end

function localDeleteIfExists(f)
if isfile(f)
    delete(f);
end
end
