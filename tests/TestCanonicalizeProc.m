classdef TestCanonicalizeProc < PsyRATTestBase
    %Unit pins for psyrat_canonicalize_proc, the single shared canonicalizer
    %for the two alias pairs in the proc prefs struct:
    %  diffrescor/diffwpcov  (one setting under two names; diffwpcov wins)
    %  ssubjrel/sserrvar     (canonical + legacy alias; ssubjrel wins)
    %
    %The function extracted the two reconciliation blocks that lived at
    %psyrat_computevarcompwarp.m:73-122, and these are the FIRST tests the
    %two warning ids (psyrat_prefs:diffrescorClamped and
    %psyrat_prefs:sserrvarMismatch) have ever had -- before the extraction
    %the only place either id appeared in the repo was the warning call
    %itself, so a typo in a moved id would have shipped green.
    %
    %Policy of record: the END STATES here must match the pre-change warp
    %blocks byte-for-byte. testGoldenTableEndStates pins the full input
    %grid; the named methods above it pin the load-bearing rows readably.

    methods (Test)

        %------------------------------------------------------------------
        % diffrescor/diffwpcov pair
        %------------------------------------------------------------------

        function testDiffPairBothMissingDefaultsOff(testCase)
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(struct()));
            testCase.verifyEqual(proc.diffrescor, 1);
            testCase.verifyEqual(proc.diffwpcov, 1);
        end

        function testDiffPairEmptyFieldsTreatedAsMissing(testCase)
            %empty ([]) means "unset" for this pair, exactly as the old warp
            %backfills treated it
            in = struct('diffrescor', [], 'diffwpcov', []);
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(in));
            testCase.verifyEqual(proc.diffrescor, 1);
            testCase.verifyEqual(proc.diffwpcov, 1);
        end

        function testDiffPairLoneDiffrescorOnPromotes(testCase)
            %a lone diffrescor=2 must reach estimation as a consistent pair --
            %this is the contract psyrat_run's resolver comment calls out (a
            %lone diffrescor=2 used to be silently discarded downstream)
            in = struct('diffrescor', 2);
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(in));
            testCase.verifyEqual(proc.diffrescor, 2);
            testCase.verifyEqual(proc.diffwpcov, 2);
        end

        function testDiffPairLoneDiffwpcovOnClampsWithWarning(testCase)
            %SUBTLE, deliberate: the old warp block had no provided-flag; it
            %backfilled a missing diffrescor to 1 FIRST and then compared, so
            %a lone hand-built diffwpcov=2 warned even though nothing the
            %caller supplied was changed. The end state is right either way,
            %and reproducing the old block exactly (warning included) is the
            %contract, so the canonicalizer keeps that behavior.
            in = struct('diffwpcov', 2);
            proc = testCase.verifyWarning(@() psyrat_canonicalize_proc(in), ...
                'psyrat_prefs:diffrescorClamped');
            testCase.verifyEqual(proc.diffrescor, 2);
            testCase.verifyEqual(proc.diffwpcov, 2);
        end

        function testDiffPairBothAgreeWarningFree(testCase)
            in = struct('diffrescor', 2, 'diffwpcov', 2);
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(in));
            testCase.verifyEqual(proc.diffrescor, 2);
            testCase.verifyEqual(proc.diffwpcov, 2);
        end

        function testDiffPairDisagreeWarnsAndFollowsWpcov(testCase)
            %diffwpcov is decisive on a disagreeing in-range pair -- the same
            %direction as the engine's own normalization, which is the policy
            %of record (B33 ruling)
            in = struct('diffrescor', 2, 'diffwpcov', 1);
            proc = testCase.verifyWarning(@() psyrat_canonicalize_proc(in), ...
                'psyrat_prefs:diffrescorClamped');
            testCase.verifyEqual(proc.diffrescor, 1);
            testCase.verifyEqual(proc.diffwpcov, 1);
        end

        function testDiffPairOutOfRangeWpcovSkipsReconciliation(testCase)
            %DELIBERATE non-fix: an out-of-range diffwpcov skips the clamp
            %entirely and the pair is forwarded as-is (old warp guard
            %'any(diffwpcov == [1 2])'). The engine now range-rejects both
            %names at its own boundary, so the hole this skip used to leave
            %open is closed there, not here. Do not "fix" this row without a
            %ruling -- it is pinned as current behavior on purpose.
            in = struct('diffrescor', 2, 'diffwpcov', 3);
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(in));
            testCase.verifyEqual(proc.diffrescor, 2);
            testCase.verifyEqual(proc.diffwpcov, 3);
        end

        %------------------------------------------------------------------
        % ssubjrel/sserrvar pair
        %------------------------------------------------------------------

        function testSsubjrelBothMissingDefaultsOff(testCase)
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(struct()));
            testCase.verifyEqual(proc.ssubjrel, 1);
            testCase.verifyEqual(proc.sserrvar, 1);
        end

        function testSsubjrelLoneLegacyPromotes(testCase)
            in = struct('sserrvar', 2);
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(in));
            testCase.verifyEqual(proc.ssubjrel, 2);
            testCase.verifyEqual(proc.sserrvar, 2);
        end

        function testSsubjrelLoneCanonicalMirrors(testCase)
            in = struct('ssubjrel', 2);
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(in));
            testCase.verifyEqual(proc.ssubjrel, 2);
            testCase.verifyEqual(proc.sserrvar, 2);
        end

        function testSsubjrelBothAgreeWarningFree(testCase)
            in = struct('ssubjrel', 2, 'sserrvar', 2);
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(in));
            testCase.verifyEqual(proc.ssubjrel, 2);
            testCase.verifyEqual(proc.sserrvar, 2);
        end

        function testSsubjrelDisagreeWarnsCanonicalWins(testCase)
            %opposite precedence direction from the diff pair, on purpose:
            %here the CANONICAL name wins over the legacy alias
            in = struct('ssubjrel', 2, 'sserrvar', 1);
            proc = testCase.verifyWarning(@() psyrat_canonicalize_proc(in), ...
                'psyrat_prefs:sserrvarMismatch');
            testCase.verifyEqual(proc.ssubjrel, 2);
            testCase.verifyEqual(proc.sserrvar, 2);
        end

        function testSsubjrelEmptyLegacyTreatedAsMissing(testCase)
            in = struct('ssubjrel', [], 'sserrvar', []);
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(in));
            testCase.verifyEqual(proc.ssubjrel, 1);
            testCase.verifyEqual(proc.sserrvar, 1);
        end

        %------------------------------------------------------------------
        % option flags
        %------------------------------------------------------------------

        function testQuietSuppressesWarningsNotOutcomes(testCase)
            %'quiet' exists so the GUI call sites keep their historical
            %silence (their inline blocks never warned); it must change
            %NOTHING but warning emission
            inDiff = struct('diffrescor', 2, 'diffwpcov', 1);
            loud = testCase.verifyWarning(@() psyrat_canonicalize_proc(inDiff), ...
                'psyrat_prefs:diffrescorClamped');
            quiet = testCase.verifyWarningFree(...
                @() psyrat_canonicalize_proc(inDiff, 'quiet', true));
            testCase.verifyEqual(quiet, loud);

            inSs = struct('ssubjrel', 2, 'sserrvar', 1);
            loud = testCase.verifyWarning(@() psyrat_canonicalize_proc(inSs), ...
                'psyrat_prefs:sserrvarMismatch');
            quiet = testCase.verifyWarningFree(...
                @() psyrat_canonicalize_proc(inSs, 'quiet', true));
            testCase.verifyEqual(quiet, loud);
        end

        function testPairsSsubjrelLeavesDiffPairUntouched(testCase)
            %'pairs','ssubjrel' exists because two GUI call sites
            %(psyrat_prefs_back and psyrat_exec) historically canonicalized
            %ONLY the ssubjrel pair. Canonicalizing the diff pair there too
            %would clamp a hand-built diffrescor=2/diffwpcov=1 struct before
            %preflight and flip localDiffCovRequested's symmetric-OR answer,
            %changing which runs preflight rejects -- so the restriction is
            %load-bearing, not cosmetic.
            in = struct('diffrescor', 2, 'diffwpcov', 1, 'sserrvar', 2);
            proc = testCase.verifyWarningFree(...
                @() psyrat_canonicalize_proc(in, 'quiet', true, 'pairs', 'ssubjrel'));
            testCase.verifyEqual(proc.diffrescor, 2); %untouched, still disagreeing
            testCase.verifyEqual(proc.diffwpcov, 1);
            testCase.verifyEqual(proc.ssubjrel, 2);   %promoted from legacy
            testCase.verifyEqual(proc.sserrvar, 2);
        end

        function testOtherProcFieldsUntouched(testCase)
            %the canonicalizer owns exactly four fields; everything else in
            %proc must pass through untouched (it is called on full prefs
            %structs at every site)
            in = struct('diffrescor', 1, 'diffwpcov', 1, ...
                'ssubjrel', 1, 'sserrvar', 1, ...
                'family', 'gamma', 'seed', 12345, 'gammascale', [], ...
                'dodmap', struct('a', 1), 'nchains', 4);
            proc = testCase.verifyWarningFree(@() psyrat_canonicalize_proc(in));
            testCase.verifyEqual(proc, in); %full struct equality: nothing moved
        end

        function testBadOptionUsageRejected(testCase)
            %internal-seam misuse guard: a typo'd 'pairs' value silently
            %canonicalizing the diff pair at a site that must not is exactly
            %the hazard class this function exists to remove, so misuse
            %fails loudly instead of falling through
            testCase.verifyError(...
                @() psyrat_canonicalize_proc(struct(), 'pairs', 'diff'), ...
                'psyrat_canonicalize_proc:badoption');
            testCase.verifyError(...
                @() psyrat_canonicalize_proc(struct(), 'quiet'), ...
                'psyrat_canonicalize_proc:badoption');
        end

        %------------------------------------------------------------------
        % the executable spec
        %------------------------------------------------------------------

        function testGoldenTableEndStates(testCase)
            %Every row below was GENERATED from the pre-change reconciliation
            %blocks at psyrat_computevarcompwarp.m:73-122 (main @ 5ce5824) by
            %running a verbatim copy of that code over the full input grid and
            %recording the end states plus the warning id raised (2026-08-26;
            %the enumeration harness was throwaway, this table is the
            %committed executable spec). The canonicalizer replaced those
            %blocks, so it must reproduce every row exactly.
            %
            %Encoding: 'missing' = field absent; 'empty' = field present, [].
            %DIFF rows hold the ssubjrel pair at a consistent 1/1;
            %SSUBJ rows hold the diff pair at a consistent 1/1.
            CL = 'psyrat_prefs:diffrescorClamped';
            MM = 'psyrat_prefs:sserrvarMismatch';
            rows = { ...
                %pair    in1        in2        out1 out2 warning
                'diff',  'missing', 'missing', 1,   1,   'none'; ...
                'diff',  'missing', 'empty',   1,   1,   'none'; ...
                'diff',  'missing', 1,         1,   1,   'none'; ...
                'diff',  'missing', 2,         2,   2,   CL; ...
                'diff',  'missing', 3,         1,   3,   'none'; ...
                'diff',  'empty',   'missing', 1,   1,   'none'; ...
                'diff',  'empty',   'empty',   1,   1,   'none'; ...
                'diff',  'empty',   1,         1,   1,   'none'; ...
                'diff',  'empty',   2,         2,   2,   CL; ...
                'diff',  'empty',   3,         1,   3,   'none'; ...
                'diff',  1,         'missing', 1,   1,   'none'; ...
                'diff',  1,         'empty',   1,   1,   'none'; ...
                'diff',  1,         1,         1,   1,   'none'; ...
                'diff',  1,         2,         2,   2,   CL; ...
                'diff',  1,         3,         1,   3,   'none'; ...
                'diff',  2,         'missing', 2,   2,   'none'; ...
                'diff',  2,         'empty',   2,   2,   'none'; ...
                'diff',  2,         1,         1,   1,   CL; ...
                'diff',  2,         2,         2,   2,   'none'; ...
                'diff',  2,         3,         2,   3,   'none'; ...
                'diff',  3,         'missing', 1,   1,   CL; ...
                'diff',  3,         'empty',   1,   1,   CL; ...
                'diff',  3,         1,         1,   1,   CL; ...
                'diff',  3,         2,         2,   2,   CL; ...
                'diff',  3,         3,         3,   3,   'none'; ...
                'ssubj', 'missing', 'missing', 1,   1,   'none'; ...
                'ssubj', 'missing', 'empty',   1,   1,   'none'; ...
                'ssubj', 'missing', 1,         1,   1,   'none'; ...
                'ssubj', 'missing', 2,         2,   2,   'none'; ...
                'ssubj', 'empty',   'missing', 1,   1,   'none'; ...
                'ssubj', 'empty',   'empty',   1,   1,   'none'; ...
                'ssubj', 'empty',   1,         1,   1,   'none'; ...
                'ssubj', 'empty',   2,         2,   2,   'none'; ...
                'ssubj', 1,         'missing', 1,   1,   'none'; ...
                'ssubj', 1,         'empty',   1,   1,   'none'; ...
                'ssubj', 1,         1,         1,   1,   'none'; ...
                'ssubj', 1,         2,         1,   1,   MM; ...
                'ssubj', 2,         'missing', 2,   2,   'none'; ...
                'ssubj', 2,         'empty',   2,   2,   'none'; ...
                'ssubj', 2,         1,         2,   2,   MM; ...
                'ssubj', 2,         2,         2,   2,   'none'; ...
                };

            for r = 1:size(rows, 1)
                [pair, in1, in2, out1, out2, warnid] = rows{r, :};
                label = sprintf('row %d: %s(%s,%s)', r, pair, ...
                    localDescribe(in1), localDescribe(in2));

                %assemble the input struct; the OTHER pair is held consistent
                %so any cross-talk between the two blocks would surface as a
                %failed field check
                if strcmp(pair, 'diff')
                    in = struct('ssubjrel', 1, 'sserrvar', 1);
                    in = localSetField(in, 'diffrescor', in1);
                    in = localSetField(in, 'diffwpcov', in2);
                    f1 = 'diffrescor'; f2 = 'diffwpcov';
                    o1 = 'ssubjrel';   o2 = 'sserrvar';
                else
                    in = struct('diffrescor', 1, 'diffwpcov', 1);
                    in = localSetField(in, 'ssubjrel', in1);
                    in = localSetField(in, 'sserrvar', in2);
                    f1 = 'ssubjrel';   f2 = 'sserrvar';
                    o1 = 'diffrescor'; o2 = 'diffwpcov';
                end

                if strcmp(warnid, 'none')
                    proc = testCase.verifyWarningFree(...
                        @() psyrat_canonicalize_proc(in), label);
                else
                    proc = testCase.verifyWarning(...
                        @() psyrat_canonicalize_proc(in), warnid, label);
                end

                testCase.verifyEqual(proc.(f1), out1, [label ' -> ' f1]);
                testCase.verifyEqual(proc.(f2), out2, [label ' -> ' f2]);
                %the held-consistent other pair must come through unchanged
                testCase.verifyEqual(proc.(o1), 1, [label ' -> ' o1]);
                testCase.verifyEqual(proc.(o2), 1, [label ' -> ' o2]);
            end
        end

        %------------------------------------------------------------------
        % wiring pins: the warp delegates to the canonicalizer
        %------------------------------------------------------------------

        function testWarpDelegatesToCanonicalizer(testCase)
            %the warp must delegate BOTH alias blocks to the shared
            %canonicalizer: exactly one call, with the two warning ids living
            %only in the canonicalizer now. The positive controls for the
            %zero-counts are the canonicalizer-side counts -- a renamed or
            %typo'd id cannot make the warp-side zeros pass vacuously.
            warpSrc = fileread(fullfile(testCase.projectRoot(), ...
                'subroutines', 'estimation', 'psyrat_computevarcompwarp.m'));
            canonSrc = fileread(fullfile(testCase.projectRoot(), ...
                'subroutines', 'filefunctions', 'psyrat_canonicalize_proc.m'));
            testCase.verifyEqual(numel(strfind(warpSrc, 'psyrat_canonicalize_proc(')), 1, ...
                'the warp must call the shared canonicalizer exactly once');
            testCase.verifyEqual(numel(strfind(warpSrc, 'psyrat_prefs:diffrescorClamped')), 0, ...
                'the diff-pair warning id must have moved out of the warp');
            testCase.verifyEqual(numel(strfind(warpSrc, 'psyrat_prefs:sserrvarMismatch')), 0, ...
                'the ssubjrel-pair warning id must have moved out of the warp');
            testCase.verifyEqual(numel(strfind(canonSrc, 'psyrat_prefs:diffrescorClamped')), 1, ...
                'positive control: the diff-pair id lives in the canonicalizer');
            testCase.verifyEqual(numel(strfind(canonSrc, 'psyrat_prefs:sserrvarMismatch')), 1, ...
                'positive control: the ssubjrel-pair id lives in the canonicalizer');
        end

        function testWarpForwardsCanonicalValuesToEngine(testCase)
            %end-to-end through the REAL warp with a stub engine first on the
            %path (the same capture pattern as TestCliEntryPoint): inputs that
            %arrive under only ONE of a pair's two names must reach the engine
            %as a consistent, canonicalized pair. This is the behavioral proof
            %that the extraction changed nothing at the warp's engine hop.

            %lone legacy sserrvar=2 (canonical field absent) -> the engine's
            %'sserrvar' argument (fed from proc.ssubjrel at the call site)
            %must be 2, i.e. the legacy value was promoted into the canonical
            %field before the hop
            prefs = psyrat_defaults();
            prefs.proc = rmfield(prefs.proc, 'ssubjrel');
            prefs.proc.sserrvar = 2;
            args = localWarpCapture(testCase, prefs);
            testCase.verifyEqual(args.sserrvar, 2, ...
                'a lone legacy sserrvar=2 must reach the engine as 2');

            %lone diffrescor=2 (diffwpcov absent) -> BOTH names arrive as 2
            %(the lone-name promotion psyrat_run's resolver comment documents)
            prefs = psyrat_defaults();
            prefs.proc = rmfield(prefs.proc, 'diffwpcov');
            prefs.proc.diffrescor = 2;
            args = localWarpCapture(testCase, prefs);
            testCase.verifyEqual(args.diffrescor, 2, ...
                'a lone diffrescor=2 must survive to the engine');
            testCase.verifyEqual(args.diffwpcov, 2, ...
                'the other name for the same setting must follow');
        end

        function testStartprocDelegatesToCanonicalizer(testCase)
            %the five GUI reconciliation sites must all route through the
            %shared canonicalizer: four calls total (psyrat_procprefs and
            %psyrat_prefs_save reconcile both pairs; psyrat_prefs_back and
            %psyrat_exec reconcile only the ssubjrel pair), and the old
            %inline blocks must be gone. The zero-count needles carry
            %anti-vacuity fixtures below (the B35 lesson: a substring pin
            %whose needle drifts from the text it hunts cannot fail).
            guiSrc = fileread(fullfile(testCase.projectRoot(), ...
                'subroutines', 'guis', 'psyrat_startproc.m'));
            canonSrc = fileread(fullfile(testCase.projectRoot(), ...
                'subroutines', 'filefunctions', 'psyrat_canonicalize_proc.m'));

            testCase.verifyEqual(numel(strfind(guiSrc, 'psyrat_canonicalize_proc(')), 4, ...
                'psyrat_startproc must call the canonicalizer at exactly its four sites');
            testCase.verifyEqual(numel(strfind(guiSrc, '''pairs'',''ssubjrel''')), 2, ...
                'exactly the back and exec sites restrict to the ssubjrel pair');

            %old inline ssubjrel blocks gone; the needle is validated against
            %a verbatim copy of the historical text so a typo'd needle fails
            oldSsubjGuard = 'if ~isfield(psyrat_prefs.proc,''ssubjrel'')';
            legacyFixture = ['if ~isfield(psyrat_prefs.proc,''ssubjrel'') || ' ...
                'isempty(psyrat_prefs.proc.ssubjrel)'];
            testCase.assertEqual(numel(strfind(legacyFixture, oldSsubjGuard)), 1, ...
                'anti-vacuity: the needle must match the historical block text');
            testCase.verifyEqual(numel(strfind(guiSrc, oldSsubjGuard)), 0, ...
                'the inline ssubjrel reconcilers must be gone from the GUI');

            %old inline diff blocks gone; positive control: the block-local
            %name still exists inside the canonicalizer itself
            testCase.assertGreaterThan(numel(strfind(canonSrc, 'diffwpcov_missing')), 0, ...
                'anti-vacuity: the needle matches the canonicalizer''s own block');
            testCase.verifyEqual(numel(strfind(guiSrc, 'diffwpcov_missing')), 0, ...
                'the inline diff-pair reconcilers must be gone from the GUI');
        end
    end
end

%--------------------------------------------------------------------------
function args = localWarpCapture(testCase, prefs)
%Run the REAL psyrat_computevarcompwarp with a stub engine first on the path
%and return the name/value arguments the engine was actually called with.
%Mirrors TestInstallationAndWorkflowValidation's mocked-warp-run scaffolding.
mockDir = localCreateEngineCaptureMockDir();
testCase.addTeardown(@() localCleanupMock(mockDir));
%clear the capture variable even when the warp errors mid-run, so a failure
%here cannot leak state into later tests
testCase.addTeardown(@() evalin('base', 'clear PSYRAT_TEST_ENGINE_ARGS'));
addpath(mockDir, '-begin');
clear psyrat_computevarcomp psyrat_checkconv;

saveDir = tempname;
mkdir(saveDir);
testCase.addTeardown(@() rmdir(saveDir, 's'));

%fields psyrat_defaults does not centralize (dynrel) and small iteration
%counts; the fields under test are set by the caller
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
data.proc.savename = 'canon.psyrat';
data.proc.measheader = 'meas';

psyrat_computevarcompwarp('psyrat_prefs', prefs, 'psyrat_data', data);
args = evalin('base', 'PSYRAT_TEST_ENGINE_ARGS');
evalin('base', 'clear PSYRAT_TEST_ENGINE_ARGS');
end

function localCleanupMock(mockDir)
rmpath(mockDir);
rmdir(mockDir, 's');
end

function mockDir = localCreateEngineCaptureMockDir()
%Stub psyrat_computevarcomp that records its arguments and reports a converged
%fit, plus a pass-through psyrat_checkconv (the same stub TestCliEntryPoint
%uses). Written to a temp folder the caller puts first on the path.
mockDir = tempname;
mkdir(mockDir);

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
    'end\n']);
fclose(fid);

fid = fopen(fullfile(mockDir, 'psyrat_checkconv.m'), 'w');
fprintf(fid, [ ...
    'function RELout = psyrat_checkconv(REL)\n' ...
    'RELout = REL;\n' ...
    'end\n']);
fclose(fid);
end

function tbl = localOneFacetTable()
%Deterministic one-facet long table (no event/time columns): 20 participants
%x 8 trials, enough shape for the shared preflight regardless of which alias
%fields a scenario sets.
ids = {}; meas = [];
for s = 1:20
    for t = 1:8
        ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
        meas(end+1, 1) = 5 + 0.1 * t + 0.05 * s; %#ok<AGROW>
    end
end
tbl = table(ids, meas, 'VariableNames', {'id', 'meas'});
end

function s = localSetField(s, name, value)
%'missing' = leave the field absent; 'empty' = present but []; else numeric
if ischar(value)
    switch value
        case 'missing'
            %leave absent
        case 'empty'
            s.(name) = [];
        otherwise
            error('TestCanonicalizeProc:badRow', ...
                'Unknown golden-table sentinel ''%s''.', value);
    end
else
    s.(name) = value;
end
end

function d = localDescribe(value)
if ischar(value)
    d = value;
else
    d = num2str(value);
end
end
