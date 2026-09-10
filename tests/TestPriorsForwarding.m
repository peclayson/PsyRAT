classdef TestPriorsForwarding < PsyRATTestBase
    %Regression pins for finding G33: priors set through psyrat_run (or a
    %replayed sidecar) were stored in psyrat_prefs.proc.priors and recorded in
    %the run-config sidecar, but psyrat_computevarcompwarp never forwarded a
    %'priors' argument to the engine, so every warp-mediated run fit the
    %defaults while the sidecar claimed the custom set. These are the FIRST
    %tests that pass a non-default prior through the warp at all -- before the
    %fix, the only non-default priors in the suite reached the engine via a
    %DIRECT psyrat_computevarcomp call (TestStanModelSnapshots), which is
    %exactly the path that always worked.
    %
    %Two guarantees are pinned from both sides:
    %  1. A customized proc.priors reaches the generated Stan model, and a
    %     partial override keeps the untouched fields at their defaults
    %     (psyrat_merge_priors semantics, end to end).
    %  2. A pure default run emits the SAME engine argument list as before the
    %     fix (no 'priors' entry) -- proc.priors is never empty (psyrat_defaults
    %     seeds the full struct), so the warp's isequal-against-defaults guard
    %     is what keeps default runs bit-identical.
    %The sidecar half: the warp pins the priors the engine ACTUALLY used
    %(REL.priors) back into the prefs before psyrat_write_runconfig, mirroring
    %the dispersion/gammascale resolved-value pins next to it.

    methods (Test)

        function testRunForwardsCustomPriorsToEmittedStan(testCase)
            % The G33 defect, pinned behaviorally through the documented entry
            % point: psyrat_run('priors',...) must change the generated model.
            % Before the fix this emitted cauchy(0,40) -- the default -- and
            % this assertion failed.
            emitFile = [tempname '.stan'];
            setenv('PSYRAT_EMIT_MODEL_FILE', emitFile);
            testCase.addTeardown(@() setenv('PSYRAT_EMIT_MODEL_FILE', ''));
            testCase.addTeardown(@() localDeleteIfExists(emitFile));
            tdir = localTempDir(testCase);

            % a partial override: only single.sig_u moves, everything else in
            % the single family must stay at its documented default
            pr = struct('single', struct('sig_u', 7.25));
            testCase.verifyError(@() psyrat_run('data', localOneFacetTable(), ...
                'priors', pr, 'savepath', tdir, 'savename', 'x.psyrat'), ...
                'PsyRAT:emitModelOnly');
            testCase.assertTrue(isfile(emitFile), ...
                'Emit hook did not write a model file.');
            stanText = fileread(emitFile);

            testCase.verifySubstring(stanText, 'sig_u ~ cauchy(0,7.25);', ...
                'The custom prior scale must reach the generated Stan model.');
            % partial-override safety: the untouched residual scale keeps its
            % default (psyrat_default_priors single.sig_e = 40)
            testCase.verifySubstring(stanText, 'sig_e ~ cauchy(0,40);', ...
                'A partial override must leave unspecified fields at defaults.');
        end

        function testWarpOmitsPriorsOnDefaultRunAndForwardsCustom(testCase)
            % Bit-identical default runs: proc.priors is NEVER empty (it is
            % seeded with the full default struct), so the warp must key on
            % isequal-with-defaults, not emptiness. A bare non-empty guard
            % would add a 'priors' entry to every run's engine argument list.
            prefs = psyrat_defaults();
            args = localWarpCapture(testCase, prefs);
            testCase.verifyFalse(isfield(args, 'priors'), ...
                'A pure default run must not forward a ''priors'' argument.');

            % ... and a customized struct must arrive verbatim
            prefs2 = psyrat_defaults();
            prefs2.proc.priors.single.sig_u = 7.25;
            args2 = localWarpCapture(testCase, prefs2);
            testCase.assertTrue(isfield(args2, 'priors'), ...
                'A customized proc.priors must be forwarded to the engine.');
            testCase.verifyEqual(args2.priors.single.sig_u, 7.25);
        end

        function testRunconfigRecordsFittedPriors(testCase)
            % The sidecar must record the priors the engine ACTUALLY used
            % (REL.priors, the merged result), not the raw preference. The stub
            % engine here reports a marker value distinct from the preference;
            % the sidecar must carry the marker. Before the fix the sidecar
            % carried the preference of a run that never applied it.
            prefs = psyrat_defaults();
            prefs.proc.priors = struct('single', struct('sig_u', 7.25));
            [~, sidecar] = localWarpCapture(testCase, prefs, ...
                'relPriors', struct('single', struct('sig_u', 99)));
            cfg = jsondecode(fileread(sidecar));
            testCase.verifyEqual(cfg.priors.single.sig_u, 99, ...
                'The sidecar must record the priors the engine reported fitting.');
        end

        function testDefaultPriorsSurviveJsonRoundTrip(testCase)
            % CHARACTERIZATION GUARD, not a regression test of the fix: it
            % never touches the warp. It pins the ASSUMPTION the warp's
            % isequal-against-defaults guard rests on for replays -- a default
            % sidecar's cfg.priors comes back through jsondecode, and if the
            % round trip perturbed the struct, a default replay would be
            % forwarded as "custom" (numerically identical, but a different
            % engine argument list than a fresh default run). If MATLAB's JSON
            % handling or the priors struct's shapes ever change, this fails
            % first and names the assumption.
            rt = jsondecode(jsonencode(psyrat_default_priors));
            testCase.verifyEqual(rt, psyrat_default_priors, ...
                'psyrat_default_priors must round-trip jsonencode/jsondecode.');
        end

    end
end

%--------------------------------------------------------------------------
function [args, sidecar] = localWarpCapture(testCase, prefs, varargin)
%Run the REAL psyrat_computevarcompwarp with a stub engine first on the path
%and return the name/value arguments the engine was called with, plus the path
%of the run-config sidecar the warp wrote. Mirrors TestCanonicalizeProc's
%scaffolding (the same stub TestCliEntryPoint uses). Optional 'relPriors'
%makes the stub report that struct as REL.priors, exercising the warp's
%fitted-value pin-back into the sidecar.
p = inputParser;
addParameter(p, 'relPriors', []);
parse(p, varargin{:});

mockDir = localCreateEngineCaptureMockDir(p.Results.relPriors);
testCase.addTeardown(@() localCleanupMock(mockDir));
testCase.addTeardown(@() evalin('base', 'clear PSYRAT_TEST_ENGINE_ARGS'));
addpath(mockDir, '-begin');
clear psyrat_computevarcomp psyrat_checkconv;

saveDir = localTempDir(testCase);

prefs.proc.nchains = 3;
prefs.proc.nwarmup = 50;
prefs.proc.nsampling = 50;
prefs.proc.niter = 100;
prefs.proc.diffest = 1;
prefs.proc.dynrel = 1;
prefs.proc.traceplots = 1;

data = struct();
data.proc = struct();
data.proc.data = localOneFacetTable();
data.proc.savepath = saveDir;
data.proc.savename = 'priorsfwd.psyrat';
data.proc.measheader = 'meas';

psyrat_computevarcompwarp('psyrat_prefs', prefs, 'psyrat_data', data);
args = evalin('base', 'PSYRAT_TEST_ENGINE_ARGS');
evalin('base', 'clear PSYRAT_TEST_ENGINE_ARGS');

sc = dir(fullfile(saveDir, '*_runconfig.json'));
if isempty(sc)
    sidecar = '';
else
    sidecar = fullfile(saveDir, sc(1).name);
end
end

function localCleanupMock(mockDir)
rmpath(mockDir);
rmdir(mockDir, 's');
end

function mockDir = localCreateEngineCaptureMockDir(relPriors)
%Stub psyrat_computevarcomp that records its arguments and reports a converged
%fit, plus a pass-through psyrat_checkconv. When relPriors is non-empty the
%stub also reports it as REL.priors, the field the warp's pin-back reads.
mockDir = tempname;
mkdir(mockDir);

if isempty(relPriors)
    priorsLine = '';
else
    priorsLine = sprintf('REL.priors = %s;\n', localStructLiteral(relPriors));
end

fid = fopen(fullfile(mockDir, 'psyrat_computevarcomp.m'), 'w');
fprintf(fid, [ ...
    'function REL = psyrat_computevarcomp(varargin)\n' ...
    'args = struct();\n' ...
    'for k = 1:2:numel(varargin)\n' ...
    '    if ischar(varargin{k})\n' ...
    '        args.(varargin{k}) = varargin{k+1};\n' ...
    '    end\n' ...
    'end\n' ...
    'assignin(''base'', ''PSYRAT_TEST_ENGINE_ARGS'', args);\n' ...
    'REL = struct();\n' ...
    'REL.analysis = ''ic'';\n' ...
    'REL.out = struct();\n' ...
    'REL.out.conv = struct(''converged'', true);\n' ...
    priorsLine ...
    'end\n']);
fclose(fid);

fid = fopen(fullfile(mockDir, 'psyrat_checkconv.m'), 'w');
fprintf(fid, [ ...
    'function RELout = psyrat_checkconv(REL)\n' ...
    'RELout = REL;\n' ...
    'end\n']);
fclose(fid);
end

function s = localStructLiteral(v)
%Render a one-level-nested scalar-double struct as MATLAB source for the stub.
%Only the shapes this test file passes are supported; anything else errors so
%a widened use fails loudly rather than generating wrong code.
assert(isstruct(v) && isscalar(v), 'localStructLiteral: scalar struct only');
parts = {};
fams = fieldnames(v);
for i = 1:numel(fams)
    inner = v.(fams{i});
    assert(isstruct(inner) && isscalar(inner), ...
        'localStructLiteral: one nested struct level only');
    keys = fieldnames(inner);
    innerParts = {};
    for j = 1:numel(keys)
        val = inner.(keys{j});
        assert(isnumeric(val) && isscalar(val), ...
            'localStructLiteral: scalar numeric leaves only');
        innerParts{end+1} = sprintf('''%s'', %.17g', keys{j}, val); %#ok<AGROW>
    end
    parts{end+1} = sprintf('''%s'', struct(%s)', fams{i}, ...
        strjoin(innerParts, ', ')); %#ok<AGROW>
end
s = sprintf('struct(%s)', strjoin(parts, ', '));
end

function tbl = localOneFacetTable()
%Deterministic one-facet long table (no event/time columns): 20 participants
%x 8 trials, enough shape for the shared preflight.
ids = {}; meas = [];
for s = 1:20
    for t = 1:8
        ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
        meas(end+1, 1) = 5 + 0.1 * t + 0.05 * s; %#ok<AGROW>
    end
end
tbl = table(ids, meas, 'VariableNames', {'id', 'meas'});
end

function tdir = localTempDir(testCase)
tdir = tempname;
mkdir(tdir);
testCase.addTeardown(@() rmdir(tdir, 's'));
end

function localDeleteIfExists(f)
if isfile(f)
    delete(f);
end
end
