classdef TestEngineAliasOptions < PsyRATTestBase
    %Pins for the engine-boundary halves of the alias canonicalization:
    %
    % (a) psyrat_computevarcomp accepts 'ssubjrel' (the canonical public
    %     name) as an alias for its historical 'sserrvar' option. Before
    %     this, 'ssubjrel' passed directly to the engine was silently
    %     dropped by the name/value scan (the engine has no unknown-argument
    %     validation), so the GROUP-LEVEL model was fit with no signal --
    %     the same silent-wrong-model class as B33 -- while the engine's own
    %     error text advertised the name it could not accept.
    % (b) the diffrescor/diffwpcov pair is range-validated to {1,2} at the
    %     engine boundary, mirroring the ssrescor/dispersion siblings.
    %     Before this, an out-of-range value was silently normalized away
    %     while REL stamped the raw value (false provenance), and an
    %     explicit empty [] for diffwpcov counted as "provided" and clamped
    %     a lone diffrescor=2 to 1. Empty now means "unset" (the dispersion
    %     convention) for both names -- a deliberate, unreachable-from-
    %     production micro-change pinned here.
    %
    %All model-content assertions use the PSYRAT_EMIT_MODEL_FILE hook (the
    %TestDodConcurrentGuard pattern): the run aborts with
    %PsyRAT:emitModelOnly after the Stan model is written, before CmdStan.

    methods (Test)

        function testSsubjrelAliasRoutesLikeSserrvar(testCase)
            %the alias must select the SAME model as the legacy name, and
            %the positive control proves the option is load-bearing at all
            %(without it, all three emits being equal would also "pass")
            stanLegacy = localEmit(testCase, localOneFacetTable(), ...
                {'sserrvar', 2});
            stanAlias = localEmit(testCase, localOneFacetTable(), ...
                {'ssubjrel', 2});
            stanPlain = localEmit(testCase, localOneFacetTable(), {});

            testCase.verifyEqual(stanAlias, stanLegacy, ...
                'ssubjrel=2 must emit byte-identically to sserrvar=2');
            testCase.verifyNotEqual(stanPlain, stanLegacy, ...
                'positive control: the subject-level toggle must change the model');
        end

        function testSsubjrelWinsOverDisagreeingSserrvarWithWarning(testCase)
            %both names, disagreeing: canonical wins, with the same warning
            %id the shared canonicalizer uses upstream
            emitFile = [tempname '.stan'];
            setenv('PSYRAT_EMIT_MODEL_FILE', emitFile);
            testCase.addTeardown(@() setenv('PSYRAT_EMIT_MODEL_FILE', ''));
            testCase.addTeardown(@() localDeleteIfExists(emitFile));

            testCase.verifyWarning(@() localRunSwallowEmit(...
                localOneFacetTable(), {'ssubjrel', 2, 'sserrvar', 1}), ...
                'psyrat_prefs:sserrvarMismatch');
            testCase.assertTrue(isfile(emitFile), ...
                'the disagreeing pair must still reach the model build');
            stanBoth = fileread(emitFile);

            stanLegacy = localEmit(testCase, localOneFacetTable(), ...
                {'sserrvar', 2});
            testCase.verifyEqual(stanBoth, stanLegacy, ...
                'ssubjrel must win the disagreement (subject-level model emitted)');
        end

        function testDiffrescorRejectsOutOfRange(testCase)
            %the emit hook is armed so a guard regression fails fast at the
            %model build instead of launching a live CmdStan compile
            localArmEmitHook(testCase);
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDiffTable(), 'diffest', 2, 'diffrescor', 3, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:diffrescor');
        end

        function testDiffwpcovRejectsOutOfRange(testCase)
            %before the guard, diffwpcov=3 silently forced diffrescor=1 (via
            %the provided-flag branch of the normalization) while
            %REL.diffwpcov stamped the raw 3
            localArmEmitHook(testCase);
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDiffTable(), 'diffest', 2, 'diffwpcov', 3, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:diffwpcov');
        end

        function testSserrvarRejectsOutOfRangeEitherSpelling(testCase)
            %completes the range-validation family for the third aliased
            %option: without this, an out-of-range value under the NEW
            %'ssubjrel' spelling would flip from silently-dropped to an
            %untyped undefined-variable crash deep in the analysis selector
            %(and the legacy spelling always crashed that way). One id
            %covers both spellings.
            localArmEmitHook(testCase);
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localOneFacetTable(), 'ssubjrel', 3, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:sserrvar');
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localOneFacetTable(), 'sserrvar', 3, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:sserrvar');

            %and [] means unset for this option too (dispersion convention)
            stanPlain = localEmit(testCase, localOneFacetTable(), {});
            stanEmptyS = localEmit(testCase, localOneFacetTable(), ...
                {'sserrvar', []});
            testCase.verifyEqual(stanEmptyS, stanPlain, ...
                'sserrvar=[] must behave exactly like an omitted sserrvar');
        end

        function testEmptyDiffOptionsTreatedAsUnset(testCase)
            %empty [] = "unset", the dispersion convention. The load-bearing
            %half is diffwpcov: previously 'diffwpcov',[] counted as
            %provided and CLAMPED a lone diffrescor=2 to 1, so this emit
            %came out non-concurrent -- the pinned micro-change is that the
            %concurrent request now survives
            stanConcurrent = localEmit(testCase, localDiffTable(), ...
                {'diffest', 2, 'diffrescor', 2});
            testCase.verifySubstring(stanConcurrent, 'Lrescor', ...
                'baseline: a lone diffrescor=2 emits the concurrent model');
            %positive control: the toggle must actually discriminate, so a
            %dead-always-concurrent regression cannot keep all emits equal
            testCase.verifyNotEqual(stanConcurrent, ...
                localEmit(testCase, localDiffTable(), {'diffest', 2}), ...
                'the concurrent and plain difference models must differ');

            stanEmptyW = localEmit(testCase, localDiffTable(), ...
                {'diffest', 2, 'diffrescor', 2, 'diffwpcov', []});
            testCase.verifyEqual(stanEmptyW, stanConcurrent, ...
                'diffwpcov=[] must behave exactly like an omitted diffwpcov');

            %and the diffrescor half: [] behaves like omitted (plain model)
            stanPlain = localEmit(testCase, localDiffTable(), ...
                {'diffest', 2});
            stanEmptyR = localEmit(testCase, localDiffTable(), ...
                {'diffest', 2, 'diffrescor', []});
            testCase.verifyEqual(stanEmptyR, stanPlain, ...
                'diffrescor=[] must behave exactly like an omitted diffrescor');
        end
    end
end

%--------------------------------------------------------------------------
function stanText = localEmit(testCase, dataTbl, extraArgs)
%Run psyrat_computevarcomp with the emit hook so the generated Stan is written
%and the run aborts with PsyRAT:emitModelOnly before any CmdStan involvement;
%return the model text. Same pattern as TestDodConcurrentGuard's localEmit.
emitFile = [tempname '.stan'];
setenv('PSYRAT_EMIT_MODEL_FILE', emitFile);
restoreEmit = onCleanup(@() setenv('PSYRAT_EMIT_MODEL_FILE', '')); %#ok<NASGU>
cleanupFile = onCleanup(@() localDeleteIfExists(emitFile)); %#ok<NASGU>

args = [{'data', dataTbl, 'chains',1, 'warmup',1, 'sampling',1, ...
    'seed',12345, 'showgui',0, 'verbose',1}, extraArgs];
testCase.verifyError(@() psyrat_computevarcomp(args{:}), 'PsyRAT:emitModelOnly');
testCase.assertTrue(isfile(emitFile), 'Emit hook did not write a model file.');
stanText = fileread(emitFile);
end

function localArmEmitHook(testCase)
%Arm the emit hook without caring about the file: tests that expect a typed
%parse error use this so a guard REGRESSION aborts at the model build
%(PsyRAT:emitModelOnly, a fast clean failure) instead of compiling CmdStan.
emitFile = [tempname '.stan'];
setenv('PSYRAT_EMIT_MODEL_FILE', emitFile);
testCase.addTeardown(@() setenv('PSYRAT_EMIT_MODEL_FILE', ''));
testCase.addTeardown(@() localDeleteIfExists(emitFile));
end

function localRunSwallowEmit(dataTbl, extraArgs)
%Run the engine under the emit hook, swallowing ONLY the expected abort so a
%verifyWarning wrapper can observe the warning without the error failing it.
args = [{'data', dataTbl, 'chains',1, 'warmup',1, 'sampling',1, ...
    'seed',12345, 'showgui',0, 'verbose',1}, extraArgs];
try
    psyrat_computevarcomp(args{:});
catch ME
    if ~strcmp(ME.identifier, 'PsyRAT:emitModelOnly')
        rethrow(ME);
    end
end
end

function localDeleteIfExists(f)
if isfile(f)
    delete(f);
end
end

function tbl = localOneFacetTable()
%Deterministic one-facet long table (no event/time columns): 20 participants
%x 8 trials, enough shape for the one-facet and subject-level routings.
ids = {}; meas = [];
for s = 1:20
    for t = 1:8
        ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
        meas(end+1, 1) = 5 + 0.1 * t + 0.05 * s; %#ok<AGROW>
    end
end
tbl = table(ids, meas, 'VariableNames', {'id', 'meas'});
end

function tbl = localDiffTable()
%Deterministic two-event long table for the static two-event difference
%routing: 20 participants x 2 events x 8 trials (TestDodConcurrentGuard's
%shape); no dim column, so the run routes static.
ids = {}; meas = []; evt = {};
for s = 1:20
    for e = 1:2
        for t = 1:8
            ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
            evt{end+1, 1} = sprintf('e%d', e); %#ok<AGROW>
            meas(end+1, 1) = 4 + e + 0.1 * t + 0.05 * s; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, evt, 'VariableNames', {'id', 'meas', 'event'});
end
