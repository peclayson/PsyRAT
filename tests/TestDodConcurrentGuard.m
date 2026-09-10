classdef TestDodConcurrentGuard < PsyRATTestBase
    % The four-event difference-of-differences (analysis 9) supports only
    % non-concurrent events in EVERY family: the Gaussian builder
    % (psyrat_build_stan_dodiff) has no residual-covariance term, the case-9
    % body never reads diffrescor/diffwpcov, and the gamma copula twin is not
    % derived. Before the case-9 guard was widened it was gamma-only, so a
    % Gaussian diffest=3 run with diffrescor=2 or diffwpcov=2 silently fit the
    % non-concurrent model while REL and the run-config sidecar recorded a
    % concurrent model that was never fit.
    %
    % These tests pin the widened guard on the Gaussian side. The gamma side
    % of the SAME guard (same id) is pinned by TestGammaScaleOption/
    % testConcurrentDifferenceOfDifferencesRejectedUnderLogNuToo, and the
    % plain-DoD Stan text is separately pinned byte-for-byte by
    % TestStanModelSnapshots (dod_event golden).

    methods (Test)

        function testGaussianConcurrentDodRejectedViaDiffrescor(testCase)
            % Gaussian + diffest=3 + diffrescor=2 must error, not silently
            % drop the covariance request: the flags are stamped into REL and
            % the run-config sidecar before the guard, so an ignored request
            % would save provenance for a model that was never fit.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDodTable(), 'family','gaussian', ...
                'diffest',3, 'dodmap',{'E1','E2','E3','E4'}, ...
                'diffrescor',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:familyconcurrent');
        end

        function testGaussianConcurrentDodRejectedViaDiffwpcov(testCase)
            % The same rejection through the other spelling of the setting.
            % The guard keys on (diffrescor == 2 || diffwpcov == 2), and a
            % caller who passes only diffwpcov (the GUI's name) must not slip
            % through the diffrescor half.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDodTable(), 'family','gaussian', ...
                'diffest',3, 'dodmap',{'E1','E2','E3','E4'}, ...
                'diffwpcov',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:familyconcurrent');
        end

        function testPlainGaussianDodStillBuilds(testCase)
            % The widened guard must not touch the plain (non-concurrent) DoD
            % path: with both flags at their default of 1 the run must still
            % REACH the model build (the emit hook aborts it there, before
            % CmdStan). That reaching is the real content of this test -- it is
            % what a mutant dropping the (diffrescor == 2 || diffwpcov == 2)
            % conjunct fails. The Lrescor assertion below is a belt-and-braces
            % check, not the discriminating one: psyrat_build_stan_dodiff takes
            % no covariance flag and has no rescor branch in any family, so it
            % cannot fail for any reachable configuration.
            stanText = localEmit(testCase, localDodTable(), ...
                {'family','gaussian', 'diffest',3, 'dodmap',{'E1','E2','E3','E4'}});
            testCase.verifyEmpty(strfind(stanText, 'Lrescor'), ...
                'A plain DoD run must not emit a residual-covariance model.');
        end

        function testGaussianTwoEventConcurrentDiffUnaffected(testCase)
            % Negative control pinning the analysis == 9 conjunct.
            % The reason the conjunct is load-bearing is simply that
            % diffrescor = 2 is a DELIBERATE and SUPPORTED request on the
            % two-event difference designs (analysis 7 here, and 10): a guard
            % keyed on the flags alone would reject them outright. (An earlier
            % draft of this comment justified the conjunct by preference LEAK
            % across runs, citing the diffest ~= 1 rationale on the gammascale
            % guard. That is a non sequitur -- the gammascale rationale is
            % about protecting NON-difference designs that never read the flag,
            % and a flags-only guard would break analysis 7 whether or not any
            % value leaked.) That design must still reach the model build, and
            % the emitted model must actually be the concurrent one (Lrescor is
            % its residual-correlation Cholesky).
            stanText = localEmit(testCase, localDiffTable(), ...
                {'family','gaussian', 'diffest',2, 'diffrescor',2});
            testCase.verifySubstring(stanText, 'Lrescor', ...
                'The concurrent two-event difference must emit the rescor model.');
        end

    end
end

function stanText = localEmit(testCase, dataTbl, extraArgs)
%Run psyrat_computevarcomp with the emit hook so the generated Stan is written
%and the run aborts with PsyRAT:emitModelOnly before any CmdStan involvement;
%return the model text. Same pattern as TestGammaScaleOption's localEmit.
emitFile = [tempname '.stan'];
setenv('PSYRAT_EMIT_MODEL_FILE', emitFile);
restoreEmit = onCleanup(@() setenv('PSYRAT_EMIT_MODEL_FILE', ''));
cleanupFile = onCleanup(@() localDeleteIfExists(emitFile));

args = [{'data', dataTbl, 'chains',1, 'warmup',1, 'sampling',1, ...
    'seed',12345, 'showgui',0, 'verbose',1}, extraArgs];
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
%Deterministic FOUR-event long table for the difference-of-differences routing
%(analysis 9): 20 participants x 4 events x 8 trials, same shape as
%TestGammaScaleOption's localGammaDodTable. The four event labels must match
%the dodmap the caller passes; diffest==3 plus that map routes to analysis 9.
%The four cells carry DIFFERENT amplitudes on purpose so the cell contrasts
%the design exists to separate are exercised.
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

function tbl = localDiffTable()
%Deterministic two-event long table for the static two-event difference routing
%(analysis 7): 20 participants x 2 events x 8 trials, same shape as
%TestGammaScaleOption's localGammaDiffTable. The two events carry DIFFERENT
%amplitudes on purpose, and there is no dim column, so the run routes static.
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
