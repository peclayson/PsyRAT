classdef TestCliEntryPoint < PsyRATTestBase
    %Tests for the scripted entry point (psyrat_run) and the run-config sidecar
    %(psyrat_write_runconfig). These cover input handling and reproducibility
    %plumbing only; they do not launch CmdStan, so they run in the standard
    %unit suite.

    methods (Test)

        function testRunRequiresOneDataSource(testCase)
            % Neither 'file' nor 'data' supplied.
            testCase.verifyError( ...
                @() psyrat_run('savepath', tempdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:datasource');
        end

        function testRunRejectsBothDataSources(testCase)
            % Both 'file' and 'data' supplied.
            testCase.verifyError( ...
                @() psyrat_run('file', 'x.csv', 'data', table(1), ...
                    'savepath', tempdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:datasource');
        end

        function testRunRequiresSavepath(testCase)
            % Data source present but no output location.
            testCase.verifyError( ...
                @() psyrat_run('data', table(1)), ...
                'psyrat_run:savepath');
        end

        function testRunDetectsMissingFile(testCase)
            testCase.verifyError( ...
                @() psyrat_run('file', fullfile(tempdir, 'no_such_psyrat_input.csv'), ...
                    'savepath', tempdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:filenotfound');
        end

        function testRunconfigRoundTrip(testCase)
            tdir = localTempDir(testCase);

            prefs = psyrat_defaults();
            prefs.proc.seed = 777;
            prefs.proc.nchains = 3;
            prefs.proc.nwarmup = 1000;
            prefs.proc.nsampling = 2000;

            data = struct();
            data.proc.savepath = tdir;
            data.proc.savename = 'demo.psyrat';
            data.proc.idheader = 'id';
            data.proc.measheader = 'ern_amp';
            data.proc.eventheader = 'cond';
            data.proc.groupheader = '';
            data.proc.timeheader = '';
            data.proc.whichevents = '';
            data.proc.whichgroups = '';
            data.proc.whichtimes = '';
            data.rel.analysis = 'ic';

            sidecar = psyrat_write_runconfig(data, prefs);

            % the sidecar is named <base>_runconfig.json next to the output
            testCase.verifyEqual(exist(sidecar, 'file'), 2);
            [~, nm, ext] = fileparts(sidecar);
            testCase.verifyEqual([nm ext], 'demo_runconfig.json');

            % the recorded settings round-trip through JSON
            cfg = jsondecode(fileread(sidecar));
            testCase.verifyEqual(cfg.estimation.seed, 777);
            testCase.verifyEqual(cfg.estimation.nchains, 3);
            testCase.verifyEqual(cfg.estimation.nwarmup, 1000);
            testCase.verifyEqual(cfg.estimation.nsampling, 2000);
            testCase.verifyEqual(cfg.analysis, 'ic');
            testCase.verifyEqual(cfg.columns.measheader, 'ern_amp');
            testCase.verifyEqual(cfg.columns.idheader, 'id');
            % schema 1.2.0: the source-revision field exists and is char
            % (value is machine-dependent -- a hash in a git checkout, ''
            % elsewhere -- so only presence and type are pinned)
            testCase.verifyTrue(isfield(cfg, 'source_revision'));
            testCase.verifyTrue(ischar(cfg.source_revision) || ...
                isstring(cfg.source_revision) || isempty(cfg.source_revision));
        end

        function testRunconfigRecordsFamily(testCase)
            % The observation-likelihood family is recorded in the sidecar so a
            % gamma run is not silently replayed as Gaussian.
            tdir = localTempDir(testCase);

            prefs = psyrat_defaults();
            prefs.proc.family = 'gamma';
            data = struct();
            data.proc.savepath = tdir;
            data.proc.savename = 'gam.psyrat';
            data.rel.analysis = 'ic';
            cfg = jsondecode(fileread(psyrat_write_runconfig(data, prefs)));
            testCase.verifyEqual(cfg.design.family, 'gamma');

            % Default (Gaussian) prefs record 'gaussian'.
            prefs2 = psyrat_defaults();
            data.proc.savename = 'gau.psyrat';
            cfg2 = jsondecode(fileread(psyrat_write_runconfig(data, prefs2)));
            testCase.verifyEqual(cfg2.design.family, 'gaussian');
        end

        function testRunAcceptsGammaFamily(testCase)
            % A CLI/scripting user can request Gamma via psyrat_run; it threads
            % to the shared preflight, which rejects non-positive data before any
            % estimation (proving the family reaches the engine path).
            tdir = localTempDir(testCase);
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 0.0; 0.9; 1.0], ...
                'VariableNames', {'id','meas'});
            testCase.verifyError(@() psyrat_run('data', tbl, 'family', 'gamma', ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'meas:notpositive');
        end

        function testRunConfigReplaysFamily(testCase)
            % A gamma run's sidecar replays AS gamma (not silently downgraded to
            % Gaussian): family flows config -> proc -> preflight, which rejects
            % the non-positive data with the gamma positivity error.
            tdir = localTempDir(testCase);
            prefs = psyrat_defaults();
            prefs.proc.family = 'gamma';
            d = struct();
            d.proc.savepath = tdir;
            d.proc.savename = 'gam.psyrat';
            d.proc.idheader = 'id';       % recorded columns so the replay maps
            d.proc.measheader = 'meas';   % id/meas rather than empty headers
            d.rel.analysis = 'ic';
            sidecar = psyrat_write_runconfig(d, prefs);

            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 0.0; 0.9; 1.0], ...
                'VariableNames', {'id','meas'});
            testCase.verifyError(@() psyrat_run('config', sidecar, 'data', tbl, ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'meas:notpositive');
        end

        function testRunconfigRecordsDataProvenance(testCase)
            % The sidecar records input-data provenance: a deterministic content
            % hash of the loaded table, the source file name + file hash (when a
            % file exists), the row count, and the detected CmdStan version.
            tdir = localTempDir(testCase);

            prefs = psyrat_defaults();

            % A small long-format table plus a real source file on disk to hash.
            tbl = table([1;1;2;2], [0.1;0.2;0.3;0.4], ...
                'VariableNames', {'id', 'meas'});
            csvpath = fullfile(tdir, 'src.csv');
            writetable(tbl, csvpath);

            data = struct();
            data.proc.savepath = tdir;
            data.proc.savename = 'demo.psyrat';
            data.proc.idheader = 'id';
            data.proc.measheader = 'meas';
            data.proc.eventheader = '';
            data.proc.groupheader = '';
            data.proc.timeheader = '';
            data.proc.whichevents = '';
            data.proc.whichgroups = '';
            data.proc.whichtimes = '';
            data.rel.analysis = 'ic';
            data.raw.filename = csvpath;
            data.raw.data = tbl;

            sidecar = psyrat_write_runconfig(data, prefs);
            cfg = jsondecode(fileread(sidecar));

            % CmdStan version field is always present and a string (the detected
            % version, or '' when CmdStan is unavailable -- environment-agnostic).
            testCase.verifyTrue(isfield(cfg, 'cmdstan_version'));
            testCase.verifyTrue(ischar(cfg.cmdstan_version));

            % Data-provenance block.
            testCase.verifyTrue(isfield(cfg, 'data'));
            testCase.verifyEqual(cfg.data.nrows, 4);
            testCase.verifyEqual(cfg.data.source_file, csvpath);
            % SHA-256 hex digests are 64 lowercase hex characters.
            testCase.verifyTrue(localIsSha256(cfg.data.table_sha256));
            testCase.verifyTrue(localIsSha256(cfg.data.file_sha256));

            % Determinism: the same table hashes the same on a second write.
            sidecar2 = psyrat_write_runconfig(data, prefs);
            cfg2 = jsondecode(fileread(sidecar2));
            testCase.verifyEqual(cfg2.data.table_sha256, cfg.data.table_sha256);

            % Sensitivity: a changed value changes the table hash.
            data.raw.data.meas(1) = 99;
            sidecar3 = psyrat_write_runconfig(data, prefs);
            cfg3 = jsondecode(fileread(sidecar3));
            testCase.verifyNotEqual(cfg3.data.table_sha256, cfg.data.table_sha256);
        end

        function testRunconfigDataProvenanceDegradesWithoutRaw(testCase)
            % When no .raw is present (e.g. a legacy/edge struct), the provenance
            % block is written with empty defaults rather than erroring.
            tdir = localTempDir(testCase);
            prefs = psyrat_defaults();

            data = struct();
            data.proc.savepath = tdir;
            data.proc.savename = 'demo.psyrat';
            data.proc.idheader = 'id';
            data.proc.measheader = 'meas';
            data.proc.eventheader = '';
            data.proc.groupheader = '';
            data.proc.timeheader = '';
            data.proc.whichevents = '';
            data.proc.whichgroups = '';
            data.proc.whichtimes = '';
            data.rel.analysis = 'ic';

            sidecar = psyrat_write_runconfig(data, prefs);
            cfg = jsondecode(fileread(sidecar));

            testCase.verifyTrue(isfield(cfg, 'data'));
            testCase.verifyEqual(cfg.data.source_file, '');
            testCase.verifyEqual(cfg.data.table_sha256, '');
            testCase.verifyEmpty(cfg.data.nrows);
        end

        %% Residual covariance for co-occurring difference events -------------
        % 'diffrescor' and 'diffwpcov' are ONE setting under two names. The GUI
        % shows a single control and writes both; psyrat_run must do the same,
        % because psyrat_computevarcompwarp re-derives diffrescor FROM diffwpcov.
        % Before the reconciliation, psyrat_run(...,'diffrescor',2) with no
        % 'diffwpcov' was silently reset to 1 and the run fit the UNCORRELATED
        % residual model, with every difference-score coefficient, SEM and
        % cut-score computed under the wrong residual covariance structure.

        function testDiffrescorAloneReachesEstimation(testCase)
            % 'diffrescor',2 with no 'diffwpcov' must arrive at
            % psyrat_computevarcomp as BOTH names set to 2. A stub engine on the
            % path captures what the engine was actually called with, so this
            % pins the whole psyrat_run -> warp -> engine hop without CmdStan.
            args = localCaptureEngineArgs(testCase, {'diffest', 2, 'diffrescor', 2});
            testCase.verifyEqual(args.diffrescor, 2, ...
                'diffrescor=2 must survive to the engine');
            testCase.verifyEqual(args.diffwpcov, 2, ...
                'the other name for the same setting must follow');
        end

        function testDiffwpcovAloneReachesEstimation(testCase)
            % The mirror case (the one the GUI produces) must be unchanged.
            args = localCaptureEngineArgs(testCase, {'diffest', 2, 'diffwpcov', 2});
            testCase.verifyEqual(args.diffrescor, 2);
            testCase.verifyEqual(args.diffwpcov, 2);
        end

        function testDiffCovDefaultIsOff(testCase)
            % Neither name supplied -> the documented default (covariance = 0).
            args = localCaptureEngineArgs(testCase, {'diffest', 2});
            testCase.verifyEqual(args.diffrescor, 1);
            testCase.verifyEqual(args.diffwpcov, 1);
        end

        function testDiffrescorConflictErrors(testCase)
            % Two explicit arguments that disagree select structurally different
            % models, so PsyRAT must refuse rather than silently pick one.
            tdir = localTempDir(testCase);
            tbl = localDiffTable();
            testCase.verifyError(@() psyrat_run('data', tbl, 'eventcol', 'event', ...
                'diffest', 2, 'diffrescor', 2, 'diffwpcov', 1, ...
                'savepath', tdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:diffrescorConflict');
            % ...in either order of disagreement.
            testCase.verifyError(@() psyrat_run('data', tbl, 'eventcol', 'event', ...
                'diffest', 2, 'diffrescor', 1, 'diffwpcov', 2, ...
                'savepath', tdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:diffrescorConflict');
        end

        function testDiffrescorAgreeingArgumentsAreAccepted(testCase)
            % Passing both with the SAME value is not a conflict.
            args = localCaptureEngineArgs(testCase, ...
                {'diffest', 2, 'diffrescor', 2, 'diffwpcov', 2});
            testCase.verifyEqual(args.diffrescor, 2);
            testCase.verifyEqual(args.diffwpcov, 2);
        end

        function testDiffrescorRejectsOutOfRangeValue(testCase)
            tdir = localTempDir(testCase);
            testCase.verifyError(@() psyrat_run('data', localDiffTable(), ...
                'diffrescor', 3, 'savepath', tdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:diffrescor');
        end

        function testDiffrescorConfigConflictErrors(testCase)
            % A sidecar that records the two names disagreeing is internally
            % inconsistent (PsyRAT never writes one), so which model was fit
            % cannot be recovered from it: error rather than guess.
            tdir = localTempDir(testCase);
            prefs = psyrat_defaults();
            prefs.proc.diffest = 2;
            prefs.proc.diffrescor = 2;
            prefs.proc.diffwpcov = 1;   % hand-edited / corrupt combination
            d = localMinimalDataStruct(tdir, 'conflict.psyrat');
            sidecar = psyrat_write_runconfig(d, prefs);

            testCase.verifyError(@() psyrat_run('config', sidecar, ...
                'data', localDiffTable(), 'eventcol', 'event', ...
                'savepath', tdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:diffrescorConfigConflict');
        end

        %% Chain / iteration bounds (aligned with the GUI) --------------------

        function testRunRejectsSingleChain(testCase)
            % R-hat is a between-chain to within-chain variance ratio, so a
            % 1-chain run reports a convergence verdict that means nothing. The
            % GUI refuses fewer than 3 chains; the CLI must refuse them too.
            tdir = localTempDir(testCase);
            tbl = table([1;1;2;2], [0.1;0.2;0.3;0.4], 'VariableNames', {'id','meas'});
            testCase.verifyError(@() psyrat_run('data', tbl, 'chains', 1, ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'psyrat_run:chains');
            testCase.verifyError(@() psyrat_run('data', tbl, 'chains', 2, ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'psyrat_run:chains');
            testCase.verifyError(@() psyrat_run('data', tbl, 'chains', 3.5, ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'psyrat_run:chains');
        end

        function testRunRejectsInvalidIterationsAndSeed(testCase)
            tdir = localTempDir(testCase);
            tbl = table([1;1;2;2], [0.1;0.2;0.3;0.4], 'VariableNames', {'id','meas'});
            testCase.verifyError(@() psyrat_run('data', tbl, 'warmup', -1, ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'psyrat_run:warmup');
            testCase.verifyError(@() psyrat_run('data', tbl, 'sampling', 0, ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'psyrat_run:sampling');
            testCase.verifyError(@() psyrat_run('data', tbl, 'seed', -5, ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'psyrat_run:seed');
            testCase.verifyError(@() psyrat_run('data', tbl, 'seed', 1.5, ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'psyrat_run:seed');
        end

        %% Sampler overrides on the CLI --------------------------------------
        % adapt_delta and max_treedepth change NUTS step-size adaptation during
        % warmup, so they change the draws at a FIXED seed. They were exposed in
        % the GUI and honored downstream but unreachable from psyrat_run and
        % absent from the sidecar, so a GUI run that used them replayed as a
        % different sampler under the same recorded seed.

        function testSamplerOverridesReachEstimation(testCase)
            args = localCaptureEngineArgs(testCase, ...
                {'adapt_delta', 0.97, 'max_treedepth', 13});
            testCase.verifyEqual(args.adapt_delta, 0.97);
            testCase.verifyEqual(args.max_treedepth, 13);
        end

        function testSamplerOverridesDefaultToEmpty(testCase)
            % The default must stay empty, which is what keeps default runs (and
            % the golden Stan snapshots) byte-identical to before.
            args = localCaptureEngineArgs(testCase, {});
            testCase.verifyEmpty(args.adapt_delta);
            testCase.verifyEmpty(args.max_treedepth);
        end

        function testSamplerOverridesRoundTripThroughConfig(testCase)
            % Recorded in the sidecar AND restored by a replay, all the way to
            % the engine call. Without both halves a replay silently resamples.
            tdir = localTempDir(testCase);
            prefs = psyrat_defaults();
            prefs.proc.adapt_delta = 0.97;
            prefs.proc.max_treedepth = 13;
            d = localMinimalDataStruct(tdir, 'sampler.psyrat');
            sidecar = psyrat_write_runconfig(d, prefs);

            cfg = jsondecode(fileread(sidecar));
            testCase.verifyEqual(cfg.estimation.adapt_delta, 0.97);
            testCase.verifyEqual(cfg.estimation.max_treedepth, 13);

            args = localCaptureEngineArgs(testCase, {'config', sidecar});
            testCase.verifyEqual(args.adapt_delta, 0.97);
            testCase.verifyEqual(args.max_treedepth, 13);
        end

        function testRunconfigRecordsEmptySamplerOverridesByDefault(testCase)
            tdir = localTempDir(testCase);
            d = localMinimalDataStruct(tdir, 'samplerdef.psyrat');
            cfg = jsondecode(fileread(psyrat_write_runconfig(d, psyrat_defaults())));
            testCase.verifyTrue(isfield(cfg.estimation, 'adapt_delta'));
            testCase.verifyTrue(isfield(cfg.estimation, 'max_treedepth'));
            testCase.verifyEmpty(cfg.estimation.adapt_delta);
            testCase.verifyEmpty(cfg.estimation.max_treedepth);
        end

        %% Sidecar coefficient selectors (reltype / diffgcoeff) ---------------

        function testRunconfigRecordsCoefficientSelectors(testCase)
            % reltype and diffgcoeff both determine WHICH coefficient the
            % reported numbers are, so a sidecar that omitted them could name a
            % coefficient the user never saw. Both are now recorded.
            tdir = localTempDir(testCase);

            prefs = psyrat_defaults();
            prefs.view.reltype = 3;      % equivalence and stability
            prefs.view.diffgcoeff = 2;   % generalizability on the difference row
            d = localMinimalDataStruct(tdir, 'coef.psyrat');
            cfg = jsondecode(fileread(psyrat_write_runconfig(d, prefs)));

            testCase.verifyEqual(cfg.view.reltype, 3);
            testCase.verifyEqual(cfg.view.diffgcoeff, 2);
            % the pre-existing view fields are untouched
            testCase.verifyEqual(cfg.view.gcoeff, 1);
            testCase.verifyEqual(cfg.view.noccmode, 1);

            % Defaults record the psyrat_defaults values, not nothing.
            d.proc.savename = 'coefdef.psyrat';
            cfg2 = jsondecode(fileread(psyrat_write_runconfig(d, psyrat_defaults())));
            testCase.verifyEqual(cfg2.view.reltype, 1);
            testCase.verifyEqual(cfg2.view.diffgcoeff, 1);
        end

        function testRunconfigCoefficientSelectorsAreRecordOnly(testCase)
            % Honest scope: the view block is record-only by construction
            % (psyrat_run walks a fixed option map and exposes no viewer
            % options). This pins that a replay reads the sidecar without
            % choking on the two new fields, and that they are NOT silently
            % turned into estimation settings.
            tdir = localTempDir(testCase);
            prefs = psyrat_defaults();
            prefs.view.reltype = 2;
            prefs.view.diffgcoeff = 2;
            d = localMinimalDataStruct(tdir, 'ro.psyrat');
            sidecar = psyrat_write_runconfig(d, prefs);

            args = localCaptureEngineArgs(testCase, {'config', sidecar});
            testCase.verifyFalse(isfield(args, 'reltype'));
            testCase.verifyFalse(isfield(args, 'diffgcoeff'));
        end

    end

    methods (Access = private)
        function tdir = localTempDir(testCase)
            % Create a unique temp directory and remove it after the test.
            tdir = tempname;
            mkdir(tdir);
            testCase.addTeardown(@() rmdir(tdir, 's'));
        end

        function args = localCaptureEngineArgs(testCase, extraArgs)
            % Run psyrat_run end to end with a stub psyrat_computevarcomp on the
            % path, and return the name/value arguments the stub was called
            % with. This is the only headless way to prove that a psyrat_run
            % option actually REACHES estimation: the shared warp resolves
            % several settings on the way, so asserting on psyrat_run's inputs
            % alone would prove nothing about what gets fit.
            mockDir = localCreateEngineCaptureMockDir();
            testCase.addTeardown(@() rmdir(mockDir, 's'));
            cleanupPath = onCleanup(@() rmpath(mockDir)); %#ok<NASGU>
            addpath(mockDir, '-begin');
            clear psyrat_computevarcomp psyrat_checkconv;

            % Clear first, so a psyrat_run that never reached the engine cannot
            % be mistaken for a successful capture of a previous test's args.
            evalin('base', 'clear PSYRAT_TEST_ENGINE_ARGS');

            tdir = localTempDir(testCase);
            psyrat_run('data', localDiffTable(), 'eventcol', 'event', ...
                'savepath', tdir, 'savename', 'capture.psyrat', extraArgs{:});

            testCase.assertTrue(evalin('base', ...
                'exist(''PSYRAT_TEST_ENGINE_ARGS'', ''var'') == 1'), ...
                'the stub engine was never called');
            args = evalin('base', 'PSYRAT_TEST_ENGINE_ARGS');
        end
    end

end

function d = localMinimalDataStruct(tdir, savename)
% The smallest psyrat_data struct psyrat_write_runconfig needs: an output
% location, the column mapping it records, and an analysis label.
d = struct();
d.proc.savepath = tdir;
d.proc.savename = savename;
d.proc.idheader = 'id';
d.proc.measheader = 'meas';
d.proc.eventheader = 'event';
d.proc.groupheader = '';
d.proc.timeheader = '';
d.proc.whichevents = '';
d.proc.whichgroups = '';
d.proc.whichtimes = '';
d.rel.analysis = 'ic_diff';
end

function tbl = localDiffTable()
% Long-format two-event table with within-participant pairs, which is what the
% shared preflight requires once residual covariance is requested. Values are
% arbitrary: every test using it stubs out estimation.
ids = repelem((1:4)', 4, 1);
events = repmat({'A';'A';'B';'B'}, 4, 1);
meas = (1:16)' / 10;
tbl = table(ids, meas, events, 'VariableNames', {'id', 'meas', 'event'});
end

function mockDir = localCreateEngineCaptureMockDir()
% Stub psyrat_computevarcomp that records its arguments and reports a converged
% fit, plus a pass-through psyrat_checkconv. Written to a temp folder that the
% caller puts first on the path.
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
    'REL.analysis = ''ic_diff'';\n' ...
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

function tf = localIsSha256(s)
% True when s is a 64-character lowercase hex string (a SHA-256 digest).
tf = (ischar(s) || isstring(s)) && ~isempty(regexp(char(s), '^[0-9a-f]{64}$', 'once'));
end
