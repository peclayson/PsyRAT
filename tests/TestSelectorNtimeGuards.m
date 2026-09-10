classdef TestSelectorNtimeGuards < PsyRATTestBase
    % G40, CLOSED AS MOOT 2026-08-29 -- this class is the regression weld.
    %
    % The filing claimed the engine selector's arms for analysis 9
    % (diffest=3, sserrvar=1) and analysis 8 (diffest=2, sserrvar=2) carry no
    % ntime test, so a DIRECT psyrat_computevarcomp call with a time column
    % would silently route to the single-session case, with
    % psyrat_preflight_validate protecting "only the GUI/psyrat_run paths".
    % The trace (Golden Rule 3) refuted the reachability half:
    % psyrat_computevarcomp runs psyrat_preflight_validate ITSELF,
    % unconditionally, before the selector ("run consolidated validation
    % before any CmdStan setup", the call above the colnames scan), and
    % errors on ~preflight.ok. A direct call with either combination plus a
    % time column therefore already fails loudly with
    % preflight:timeWithDiffUnsupported -- probed live 2026-08-29. The
    % owner-ruled outcome (fail loud for direct callers, never a silently
    % wrong model) already holds, so no selector guard was added.
    %
    % The selector arms really are ntime-blind, which is why these pins
    % matter: they weld the closure to behavior. If the in-engine preflight
    % call is ever removed, made conditional, or these ids stop covering the
    % two combinations, the negative pins fail and G40 reopens.
    %
    % Positive controls run through the emit hook (PsyRAT:emitModelOnly, the
    % TestDodConcurrentGuard idiom): reaching the hook proves the engine
    % accepted the combination and built a model, so the negative pins
    % cannot pass by the engine rejecting everything.

    methods (Test)

        function testDodWithTimeColumnFailsLoudlyInTheEngine(testCase)
            % diffest=3 + a time column: blocked inside psyrat_computevarcomp
            % by its own preflight call, not only on the GUI/psyrat_run paths.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDodTableWithTime(), 'family', 'gaussian', ...
                'diffest', 3, 'dodmap', {'E1','E2','E3','E4'}, ...
                'chains', 1, 'warmup', 1, 'sampling', 1, 'seed', 12345, ...
                'showgui', 0, 'verbose', 1), ...
                'preflight:timeWithDiffUnsupported');
        end

        function testSubjectLevelDiffWithTimeColumnFailsLoudlyInTheEngine(testCase)
            % diffest=2 + sserrvar=2 + a time column: same in-engine block,
            % via a different preflight id -- the subject-level-with-time
            % check fires before blockDiffWithTime for this combination
            % (probed 2026-08-29). The pin is on the id that actually fires,
            % so a reordering of preflight's checks that silently ADMITTED
            % the combination would fail here rather than pass on a stale id.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDiffTableWithTime(), ...
                'diffest', 2, 'sserrvar', 2, ...
                'chains', 1, 'warmup', 1, 'sampling', 1, 'seed', 12345, ...
                'showgui', 0, 'verbose', 1), ...
                'preflight:timeWithSubjectReliabilityUnsupported');
        end

        function testDodWithoutTimeStillRoutes(testCase)
            % Positive control: the same diffest=3 request without a time
            % column must still reach the model build (analysis 9).
            stanText = localEmit(testCase, localDodTable(), ...
                {'family', 'gaussian', 'diffest', 3, ...
                'dodmap', {'E1','E2','E3','E4'}});
            testCase.verifyNotEmpty(stanText, ...
                'diffest=3 without time must still build its model');
        end

        function testSubjectLevelDiffWithoutTimeStillRoutes(testCase)
            % Positive control: diffest=2 + sserrvar=2 without a time column
            % must still reach the model build (analysis 8).
            stanText = localEmit(testCase, localDiffTable(), ...
                {'diffest', 2, 'sserrvar', 2});
            testCase.verifyNotEmpty(stanText, ...
                'diffest=2/sserrvar=2 without time must still build its model');
        end

        function testTrtDiffRouteStillAcceptsTime(testCase)
            % Positive control for the neighboring supported combination:
            % analysis 10 (diffest=2, sserrvar=1, ntime>0) is the two-facet
            % difference pipeline (the G34 un-gated route) and must keep
            % routing with a time column present.
            stanText = localEmit(testCase, localDiffTableWithTime(), ...
                {'diffest', 2, 'sserrvar', 1});
            testCase.verifyNotEmpty(stanText, ...
                'analysis 10 must keep accepting a time column');
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

args = [{'data', dataTbl, 'chains', 1, 'warmup', 1, 'sampling', 1, ...
    'seed', 12345, 'showgui', 0, 'verbose', 1}, extraArgs];
testCase.verifyError(@() psyrat_computevarcomp(args{:}), 'PsyRAT:emitModelOnly');
testCase.assertTrue(isfile(emitFile), 'Emit hook did not write a model file.');
stanText = fileread(emitFile);
end

function localDeleteIfExists(f)
if isfile(f)
    delete(f);
end
end

function tbl = localDodTable()
%Deterministic FOUR-event long table for difference-of-differences routing
%(analysis 9), same shape as TestDodConcurrentGuard's fixture: 20 participants
%x 4 events x 8 trials, distinct cell amplitudes.
ids = {}; meas = []; evt = {};
amp = [5.0 4.2 5.4 4.6];
for s = 1:20
    for e = 1:4
        for t = 1:8
            ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
            evt{end+1, 1} = sprintf('E%d', e); %#ok<AGROW>
            meas(end+1, 1) = amp(e) + 0.1 * t + 0.05 * s; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, evt, 'VariableNames', {'id','meas','event'});
end

function tbl = localDodTableWithTime()
%The four-event table crossed with two occasions: the state the G40 filing
%said would silently drop the occasion facet (ntime > 0 with diffest = 3).
ids = {}; meas = []; evt = {}; tim = {};
amp = [5.0 4.2 5.4 4.6];
for s = 1:20
    for o = 1:2
        for e = 1:4
            for t = 1:8
                ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
                evt{end+1, 1} = sprintf('E%d', e); %#ok<AGROW>
                tim{end+1, 1} = sprintf('t%d', o); %#ok<AGROW>
                meas(end+1, 1) = amp(e) + 0.1 * t + 0.05 * s + 0.02 * o; %#ok<AGROW>
            end
        end
    end
end
tbl = table(ids, meas, evt, tim, 'VariableNames', {'id','meas','event','time'});
end

function tbl = localDiffTable()
%Deterministic two-event long table for the two-event difference designs
%(analyses 7/8): 20 participants x 2 events x 8 trials.
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
tbl = table(ids, meas, evt, 'VariableNames', {'id','meas','event'});
end

function tbl = localDiffTableWithTime()
%The two-event table crossed with two occasions: analysis 10's supported
%shape, and the in-engine preflight error state when sserrvar = 2 rides
%along.
ids = {}; meas = []; evt = {}; tim = {};
for s = 1:20
    for o = 1:2
        for e = 1:2
            for t = 1:8
                ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
                evt{end+1, 1} = sprintf('e%d', e); %#ok<AGROW>
                tim{end+1, 1} = sprintf('t%d', o); %#ok<AGROW>
                meas(end+1, 1) = 4 + e + 0.1 * t + 0.05 * s + 0.02 * o; %#ok<AGROW>
            end
        end
    end
end
tbl = table(ids, meas, evt, tim, 'VariableNames', {'id','meas','event','time'});
end
