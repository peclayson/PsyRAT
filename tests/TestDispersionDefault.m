classdef TestDispersionDefault < matlab.unittest.TestCase
    %TESTDISPERSIONDEFAULT Pins the CONTEXT-AWARE default of the 'dispersion' option
    % (psyrat_computevarcomp), which selects between two structurally different
    % Gaussian-copula scaled-chi-square models for the concurrent Gamma difference
    % designs:
    %   dispersion = 1 -> per-person log-nu   (4-D person block [m1,m2,lognu1,lognu2])
    %   dispersion = 2 -> GLOBAL fixed nu     (2-D person MEAN block)
    %
    % The default is deliberately CONTEXT-AWARE: GLOBAL (2) for the GROUP-LEVEL
    % concurrent Gamma difference design (family='gamma' + diffest=2 + diffrescor=2 +
    % sserrvar=1, analyses 7/10) and per-person (1) everywhere else. Every OTHER test
    % and snapshot scenario passes 'dispersion' explicitly, so without this class the
    % default itself is untested: inverting or deleting it would leave the whole suite
    % green while silently changing the estimand of every default group-level run.
    %
    % No CmdStan required: the PSYRAT_EMIT_MODEL_FILE hook writes the generated Stan
    % and aborts with PsyRAT:emitModelOnly before any sampling, so these tests assert
    % on the emitted model text (the 2-D vs 4-D person block is the discriminator).

    methods (Test)

        function testGroupLevelConcurrentGammaDiffDefaultsToGlobalNu(testCase)
            % DEFAULT (no 'dispersion' passed) on the group-level concurrent gamma
            % difference must emit the GLOBAL-nu model: a 2-D person block and no
            % per-person log-nu SDs.
            stanText = localEmit(testCase, localDiffTable(), ...
                {'family','gamma','diffest',2,'diffwpcov',2});
            testCase.verifySubstring(stanText, 'vector<lower=0>[2] sd_id', ...
                ['Group-level concurrent gamma difference must DEFAULT to the global-nu ' ...
                'model (2-D person mean block).']);
            testCase.verifyEmpty(strfind(stanText, 'vector<lower=0>[4] sd_id'), ...
                'The default group-level model must NOT emit the 4-D per-person-log-nu block.');
        end

        function testExplicitPerPersonStillReachable(testCase)
            % dispersion=1 must still select the per-person log-nu model (the pinned
            % variant the existing goldens and S11 recovery tests depend on).
            stanText = localEmit(testCase, localDiffTable(), ...
                {'family','gamma','diffest',2,'diffwpcov',2,'dispersion',1});
            testCase.verifySubstring(stanText, 'vector<lower=0>[4] sd_id', ...
                'dispersion=1 must select the per-person log-nu (4-D person block) model.');
        end

        function testExplicitGlobalMatchesDefault(testCase)
            % An explicit dispersion=2 and the default must emit byte-identical Stan -
            % i.e. the default really is 2, not merely "something global-looking".
            args = {'family','gamma','diffest',2,'diffwpcov',2};
            defaultText  = localEmit(testCase, localDiffTable(), args);
            explicitText = localEmit(testCase, localDiffTable(), [args {'dispersion',2}]);
            testCase.verifyEqual(defaultText, explicitText, ...
                'The context-aware default must be exactly dispersion=2 for this design.');
        end

        function testNonConcurrentGammaDiffDefaultsToPerPerson(testCase)
            % Outside the CONCURRENT design (diffwpcov=1 -> diffrescor=1) the default
            % must stay per-person: the non-concurrent gamma difference builder keeps
            % its 4-D person block, so the flip must not leak into it.
            %
            % gammascale=1 is passed EXPLICITLY and is load-bearing. This is a
            % log-nu-builder test, and since the 2026-08-07 default flip an omitted
            % gammascale resolves to 2 on analysis 7, which routes to the
            % location-scale builder instead. The assertion below still passed
            % because both builders declare a 4-D sd_id - so the test silently
            % stopped exercising the branch its name and failure message describe.
            % Do not drop this argument to "use the default".
            stanText = localEmit(testCase, localDiffTable(), ...
                {'family','gamma','diffest',2,'diffwpcov',1,'gammascale',1});
            testCase.verifySubstring(stanText, 'vector<lower=0>[4] sd_id', ...
                ['The non-concurrent gamma difference must keep the per-person log-nu ' ...
                'block; the global-nu default is scoped to the CONCURRENT design.']);
        end

        function testExplicitGlobalRejectedOutsideGroupLevelDesign(testCase)
            % dispersion=2 is a GROUP-LEVEL estimand: it must be REJECTED (not silently
            % ignored) for subject-level reliability, which needs per-person dispersion.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDiffTable(), 'family','gamma', 'diffest',2, ...
                'diffwpcov',2, 'sserrvar',2, 'dispersion',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:dispersion', ...
                ['dispersion=2 must error for the subject-level design rather than ' ...
                'silently falling back to per-person dispersion.']);
        end

        function testInvalidDispersionRejected(testCase)
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDiffTable(), 'family','gamma', 'diffest',2, ...
                'diffwpcov',2, 'dispersion',0, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:dispersion', ...
                'An out-of-range dispersion must error rather than fall through.');
        end

    end
end

function stanText = localEmit(testCase, dataTbl, extraArgs)
%Run psyrat_computevarcomp with the emit hook so the generated Stan is written and the
%run aborts before CmdStan; return the model text.
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

function tbl = localDiffTable()
%Deterministic balanced two-event long table with strictly positive scores (the gamma
%family's support): 6 participants x 2 events x 8 trials.
nsub = 6; ntrl = 8; evn = {'E1','E2'};
ids = cell(nsub*2*ntrl,1); evs = cell(nsub*2*ntrl,1); meas = zeros(nsub*2*ntrl,1);
r = 0;
for p = 1:nsub
    for e = 1:2
        for t = 1:ntrl
            r = r + 1;
            ids{r} = sprintf('S%02d',p);
            evs{r} = evn{e};
            % deterministic, strictly positive, varying across p/e/t
            meas(r) = 5 + 0.25*p + 0.5*e + 0.1*t;
        end
    end
end
tbl = table(ids, meas, evs, 'VariableNames', {'id','meas','event'});
end

function localDeleteIfExists(f)
if exist(f,'file') == 2
    try, delete(f); catch, end
end
end
