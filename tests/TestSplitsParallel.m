classdef TestSplitsParallel < PsyRATTestBase
    % Deterministic guard that the PARALLEL data-splits path (splits=2) reuses the
    % validated single-trial machinery byte-for-byte at the Stan-emission level. No
    % CmdStan/MatlabStan required: the model code is captured through the
    % PSYRAT_EMIT_MODEL_FILE hook in psyrat_run_stan (which writes the assembled
    % model and aborts before any CmdStan call), exactly as TestStanModelSnapshots
    % does.
    %
    % Parallel splits declare all items-per-split counts (n_i) equal, so the split
    % MEAN is the score and n' becomes the number of splits; the engine routes a
    % splits=2 run straight through the standard one-facet (case 1) / test-retest
    % (case 5) chain, carrying REL.splits=2 only so the VIEWERS relabel "trial" ->
    % "split". The Stan itself is therefore unchanged (it still reads "trials"), and
    % these tests pin that: splits=2 must emit the committed single-trial goldens
    % (ic_base / trt_time, plain and __threaded) byte-for-byte.
    %
    % The relabeling that splits=2 DOES change lives in the viewers and is pinned
    % by TestUnitLabels; the live routing (REL.analysis is the non-_splits label)
    % and known-truth recovery for splits=2 are covered by the slow CmdStan
    % TestSplitsRecovery. The numeric calc is unchanged by construction:
    % psyrat_relsummary never reads REL.splits, so a parallel run computes the same
    % reliability as the single-trial run it routes to.

    methods (Test)

        % ---- one-facet (splits=2 -> case 1, ic_base) -------------------------

        function testParallelSingleEmitsSingleFacetGolden(testCase)
            % splits=2 with id/meas + equal-weight items-per-split and no facets
            % routes to the one-facet model; the emitted Stan must equal the
            % committed single-trial ic_base golden, plain and threaded.
            tbl = localParallelSingleTable();
            localVerifyMatchesGolden(testCase, tbl, 'ic_base');
        end

        % ---- test-retest (splits=2 + time -> case 5, trt_time) ---------------

        function testParallelTrtEmitsTestRetestGolden(testCase)
            % splits=2 with a time/occasion column routes to the test-retest
            % model; the emitted Stan must equal the committed trt_time golden,
            % plain and threaded.
            tbl = localParallelTrtTable();
            localVerifyMatchesGolden(testCase, tbl, 'trt_time');
        end

        % ---- table-independent equivalence -----------------------------------

        function testParallelMatchesSingleTrialEmission(testCase)
            % Stronger than the golden anchor: on the SAME measurements, splits=2
            % (with the required weight column) must emit exactly what a plain
            % single-trial run (splits=1, no weight column) emits. This proves the
            % parallel path routes to the identical case regardless of the
            % golden's own table, for both the one-facet and test-retest designs.
            single2 = localCaptureStan(testCase, localParallelSingleTable(), 2, false);
            single1 = localCaptureStan(testCase, ...
                removevars(localParallelSingleTable(), 'weight'), 1, false);
            testCase.verifyEqual(single2, single1, ...
                'Parallel splits=2 one-facet Stan differs from the single-trial run.');

            trt2 = localCaptureStan(testCase, localParallelTrtTable(), 2, false);
            trt1 = localCaptureStan(testCase, ...
                removevars(localParallelTrtTable(), 'weight'), 1, false);
            testCase.verifyEqual(trt2, trt1, ...
                'Parallel splits=2 test-retest Stan differs from the single-trial run.');
        end

    end
end

% =============================== local helpers ===============================

function localVerifyMatchesGolden(testCase, tbl, goldenName)
% Capture the plain and threaded splits=2 Stan for tbl and assert each equals the
% committed golden of the same name.
for threaded = [false true]
    if threaded
        suffix = '__threaded';
    else
        suffix = '';
    end
    actual = localCaptureStan(testCase, tbl, 2, threaded);
    goldenPath = fullfile(testCase.projectRoot(), 'tests', 'golden', 'stan', ...
        [goldenName suffix '.stan']);
    testCase.assertTrue(isfile(goldenPath), ...
        sprintf('Missing golden for parallel comparison: %s', goldenPath));
    expected = fileread(goldenPath);
    testCase.verifyEqual(actual, expected, ...
        sprintf(['Parallel splits=2 Stan does not match the single-trial golden ', ...
        '''%s%s''. The parallel path must reuse the standard model byte-for-byte.'], ...
        goldenName, suffix));
end
end

function stanText = localCaptureStan(testCase, tbl, splitsVal, threaded)
% Drive psyrat_computevarcomp for one table with the emit hook active and return
% the generated Stan source. The run aborts before CmdStan via the
% PSYRAT_EMIT_MODEL_FILE hook (PsyRAT:emitModelOnly). When threaded is true,
% within-chain threading is forced on with PSYRAT_FORCE_AVAILABLE_CORES. This
% mirrors localCaptureGeneratedStan in TestStanModelSnapshots.
emitFile = [tempname '.stan'];
setenv('PSYRAT_EMIT_MODEL_FILE', emitFile);
restoreEmit = onCleanup(@() setenv('PSYRAT_EMIT_MODEL_FILE', '')); %#ok<NASGU>

if threaded
    setenv('PSYRAT_FORCE_AVAILABLE_CORES', '8');
    restoreCores = onCleanup(@() setenv('PSYRAT_FORCE_AVAILABLE_CORES', '')); %#ok<NASGU>
end
cleanupFile = onCleanup(@() localDeleteIfExists(emitFile)); %#ok<NASGU>

args = {'data', tbl, 'chains', 1, 'warmup', 1, 'sampling', 1, ...
    'seed', 12345, 'showgui', 0, 'verbose', 1, 'splits', splitsVal};
if threaded
    args = [args, {'withinchain', 2}];
end

testCase.verifyError(@() psyrat_computevarcomp(args{:}), 'PsyRAT:emitModelOnly');

testCase.assertTrue(isfile(emitFile), ...
    sprintf('Emit hook did not write a model file (splits=%d, threaded=%d).', ...
    splitsVal, threaded));
stanText = fileread(emitFile);
end

function tbl = localParallelSingleTable()
% Parallel one-facet split-mean data: id/meas plus an EQUAL items-per-split
% weight column (parallel => all n_i equal, so no preflight unequal-weight
% warning). 6 participants x 5 splits. splits=2 routes this to the one-facet
% model (case 1).
nsub = 6; ksplit = 5; w = 20;
ids = {}; meas = []; weight = [];
for s = 1:nsub
    for k = 1:ksplit
        ids{end+1, 1} = sprintf('s%d', s); %#ok<AGROW>
        meas(end+1, 1) = 0.1 * k + 0.5 * s; %#ok<AGROW>
        weight(end+1, 1) = w; %#ok<AGROW>
    end
end
tbl = table(ids, meas, weight, 'VariableNames', {'id', 'meas', 'weight'});
end

function tbl = localParallelTrtTable()
% Parallel test-retest split-mean data: id/meas/weight (equal n_i) plus a
% time/occasion column. 6 participants x 3 occasions x 4 splits. splits=2 with a
% time facet routes this to the test-retest model (case 5).
nsub = 6; nocc = 3; ksplit = 4; w = 20;
ids = {}; meas = []; weight = []; tm = {};
for s = 1:nsub
    for o = 1:nocc
        for k = 1:ksplit
            ids{end+1, 1} = sprintf('s%d', s); %#ok<AGROW>
            tm{end+1, 1} = sprintf('t%d', o);  %#ok<AGROW>
            meas(end+1, 1) = 0.1 * k + 0.3 * o + 0.5 * s; %#ok<AGROW>
            weight(end+1, 1) = w; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, weight, tm, ...
    'VariableNames', {'id', 'meas', 'weight', 'time'});
end

function localDeleteIfExists(f)
if isfile(f)
    delete(f);
end
end
