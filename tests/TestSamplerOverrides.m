classdef TestSamplerOverrides < PsyRATTestBase
    % Coverage for the user-exposed CmdStan sampler overrides adapt_delta and
    % max_treedepth.
    %
    % Three layers are pinned:
    %   0. psyrat_run - the scripted entry point's own bounds on the two
    %      overrides, which mirror the GUI's preference-window validation
    %      (psyrat_startproc). The rest of the CLI surface for these overrides
    %      (that they reach the engine, and that they round-trip through the
    %      run-config sidecar) lives in TestCliEntryPoint, which owns the stub
    %      engine those checks need.
    %   1. psyrat_apply_sampler_overrides - the PsyRAT-side merge that folds the
    %      bare overrides into the CmdStan 'control' struct inside psyrat_run_stan.
    %      This is a pure name/value transform and runs with no CmdStan/MatlabStan.
    %   2. The bundled MatlabStan control interface - that 'adapt_delta' reaches
    %      the sampler as delta= and that 'max_treedepth' reaches it as max_depth=
    %      in the generated CmdStan command. max_treedepth support is the patch
    %      added to StanModel.set.control; this layer needs CmdStan to construct a
    %      StanModel, so it self-skips (assumeTrue) when CmdStan is unavailable.

    methods (Test)

        %% Layer 0: psyrat_run bounds (pure, no CmdStan) -----------------------

        function testRunRejectsOutOfRangeAdaptDelta(testCase)
            % Same open interval (0,1) the GUI enforces, so a value the GUI
            % refuses cannot be smuggled in through the command line.
            tdir = localTempDirFor(testCase);
            T = localTinyTable();
            testCase.verifyError(@() psyrat_run('data', T, 'adapt_delta', 1.5, ...
                'savepath', tdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:adapt_delta');
            testCase.verifyError(@() psyrat_run('data', T, 'adapt_delta', 0, ...
                'savepath', tdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:adapt_delta');
        end

        function testRunRejectsNonIntegerMaxTreedepth(testCase)
            tdir = localTempDirFor(testCase);
            T = localTinyTable();
            testCase.verifyError(@() psyrat_run('data', T, 'max_treedepth', 12.5, ...
                'savepath', tdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:max_treedepth');
            testCase.verifyError(@() psyrat_run('data', T, 'max_treedepth', 0, ...
                'savepath', tdir, 'savename', 'x.psyrat'), ...
                'psyrat_run:max_treedepth');
        end

        function testRunValidatesOverridesRestoredFromConfig(testCase)
            % The bounds run AFTER the run-config merge, so a sidecar carrying a
            % bad value is caught too. This doubles as proof that the recorded
            % value actually reaches psyrat_run's options rather than being
            % dropped by the config reader.
            tdir = localTempDirFor(testCase);
            prefs = psyrat_defaults();
            prefs.proc.adapt_delta = 1.5;      % out of range on purpose
            d = struct();
            d.proc.savepath = tdir;
            d.proc.savename = 'bad.psyrat';
            d.rel.analysis = 'ic';
            sidecar = psyrat_write_runconfig(d, prefs);

            testCase.verifyError(@() psyrat_run('config', sidecar, ...
                'data', localTinyTable(), 'savepath', tdir, ...
                'savename', 'x.psyrat'), 'psyrat_run:adapt_delta');
        end

        %% Layer 1: psyrat_apply_sampler_overrides (pure, no CmdStan) ----------

        function testNoOverridesLeavesArgsUnchanged(testCase)
            % With no override pairs the argument list must be byte-identical, so
            % default estimation runs are provably unaffected.
            args = {'data', struct('x', 1), 'chains', 4, ...
                'control', struct('delta', 0.9)};
            out = psyrat_apply_sampler_overrides(args);
            testCase.verifyEqual(out, args);
        end

        function testEmptyOverridesAreNoOp(testCase)
            % Empty override values are treated as "not supplied": the bare pairs
            % are stripped and nothing is folded into control.
            args = {'chains', 4, 'adapt_delta', [], 'max_treedepth', [], ...
                'control', struct('delta', 0.9)};
            out = psyrat_apply_sampler_overrides(args);
            [hasAd, ~] = localGet(out, 'adapt_delta');
            [hasMt, ~] = localGet(out, 'max_treedepth');
            [hasCtrl, ctrl] = localGet(out, 'control');
            testCase.verifyFalse(hasAd, 'adapt_delta pair should be removed');
            testCase.verifyFalse(hasMt, 'max_treedepth pair should be removed');
            testCase.verifyTrue(hasCtrl);
            testCase.verifyEqual(ctrl.delta, 0.9);
            testCase.verifyFalse(isfield(ctrl, 'adapt_delta'));
            testCase.verifyFalse(isfield(ctrl, 'max_treedepth'));
        end

        function testAdaptDeltaFoldsIntoNewControl(testCase)
            % adapt_delta with no existing control -> a control struct is created.
            args = {'chains', 4, 'adapt_delta', 0.95};
            out = psyrat_apply_sampler_overrides(args);
            [hasAd, ~] = localGet(out, 'adapt_delta');
            [hasCtrl, ctrl] = localGet(out, 'control');
            testCase.verifyFalse(hasAd, 'bare adapt_delta pair should be removed');
            testCase.verifyTrue(hasCtrl, 'a control struct should be created');
            testCase.verifyEqual(ctrl.adapt_delta, 0.95);
            testCase.verifyEqual(mod(numel(out), 2), 0, 'result must stay name/value');
        end

        function testMaxTreedepthFoldsIntoExistingControl(testCase)
            % max_treedepth with an existing per-case control -> the tuned delta
            % is preserved and max_treedepth is added alongside it.
            args = {'chains', 4, 'control', struct('delta', 0.98), ...
                'max_treedepth', 12};
            out = psyrat_apply_sampler_overrides(args);
            [hasMt, ~] = localGet(out, 'max_treedepth');
            [~, ctrl] = localGet(out, 'control');
            testCase.verifyFalse(hasMt, 'bare max_treedepth pair should be removed');
            testCase.verifyEqual(ctrl.delta, 0.98, 'per-case delta must be preserved');
            testCase.verifyEqual(ctrl.max_treedepth, 12);
        end

        function testBothOverridesFoldTogether(testCase)
            args = {'control', struct('delta', 0.9), ...
                'adapt_delta', 0.99, 'max_treedepth', 15};
            out = psyrat_apply_sampler_overrides(args);
            [~, ctrl] = localGet(out, 'control');
            testCase.verifyEqual(ctrl.adapt_delta, 0.99);
            testCase.verifyEqual(ctrl.max_treedepth, 15);
            % control must appear exactly once
            keys = out(1:2:end);
            testCase.verifyEqual(sum(strcmpi(keys, 'control')), 1);
        end

        %% Layer 1b: psyrat_computevarcomp input validation (pure, no CmdStan) -

        function testComputevarcompRejectsOutOfRangeAdaptDelta(testCase)
            % A supplied adapt_delta must lie in the open interval (0,1); a value
            % >= 1 is rejected before any CmdStan setup. Guards against the
            % override silently reaching CmdStan (which would error later).
            T = localTinyTable();
            bad = @() psyrat_computevarcomp('data', T, 'chains', 4, ...
                'warmup', 1, 'sampling', 1, 'seed', 12345, 'adapt_delta', 1.5);
            testCase.verifyError(bad, 'varargin:adapt_delta');
        end

        function testComputevarcompRejectsNonIntegerMaxTreedepth(testCase)
            % A supplied max_treedepth must be a positive integer; a fractional
            % value is rejected up front.
            T = localTinyTable();
            bad = @() psyrat_computevarcomp('data', T, 'chains', 4, ...
                'warmup', 1, 'sampling', 1, 'seed', 12345, 'max_treedepth', 12.5);
            testCase.verifyError(bad, 'varargin:max_treedepth');
        end

        function testComputevarcompAcceptsEmptyOverrides(testCase)
            % Omitting the overrides must NOT trip the validation (empty is the
            % documented default). Reaching the CmdStan/preflight stage without a
            % 'varargin:*' error confirms the guards let the default path through.
            T = localTinyTable();
            try
                psyrat_computevarcomp('data', T, 'chains', 4, ...
                    'warmup', 1, 'sampling', 1, 'seed', 12345);
            catch err
                testCase.verifyEmpty(regexp(err.identifier, ...
                    '^varargin:(adapt_delta|max_treedepth)$', 'once'), ...
                    sprintf(['Default (no-override) call must not fail sampler ', ...
                    'validation, but errored with %s'], err.identifier));
            end
        end

        %% Layer 2: bundled MatlabStan control interface (needs CmdStan) -------

        function testStanCommandCarriesOverrides(testCase)
            model = localBuildStanModel(testCase);

            model.control = struct('adapt_delta', 0.95, 'max_treedepth', 12);
            cmd = sprintf('%s', model.command{:});
            testCase.verifyTrue(contains(cmd, 'delta=0.95'), ...
                'adapt_delta override should appear as delta=0.95');
            testCase.verifyTrue(contains(cmd, 'max_depth=12'), ...
                'max_treedepth override should appear as max_depth=12');
        end

        function testStanCommandDefaultsUnchanged(testCase)
            model = localBuildStanModel(testCase);
            cmd = sprintf('%s', model.command{:});
            testCase.verifyTrue(contains(cmd, 'max_depth=10'), ...
                'default NUTS max_depth of 10 should be emitted');
            testCase.verifyTrue(contains(cmd, 'delta=0.8'), ...
                'default adapt delta of 0.8 should be emitted');
        end

    end
end

function model = localBuildStanModel(testCase)
% Construct a trivial StanModel using the worktree's bundled MatlabStan. Skips
% (assumeTrue) when CmdStan/MatlabStan is unavailable so CmdStan-less CI lanes do
% not fail. Requires the patched worktree StanModel to be the loaded classdef
% (run the suite with clear classes if a second copy shadows the path).
root = testCase.projectRoot();
addpath(genpath(fullfile(root, 'bundled_dependents', 'MatlabStan-2.15.1.0')));

wd = tempname;
mkdir(wd);
% Build inside the temp dir so any stray MatlabStan artifacts (anon_model.stan,
% output.csv, compiled binary) land there and are removed on teardown rather
% than polluting the tests/ folder.
origDir = cd(wd);
% teardowns run LIFO: register rmdir first so cd-back runs before it
testCase.addTeardown(@() rmdir(wd, 's'));
testCase.addTeardown(@() cd(origDir));

code = {'parameters { real y; } model { y ~ normal(0,1); }'};
model = [];
try
    model = StanModel('model_code', code, 'working_dir', wd, ...
        'file_overwrite', true);
catch err
    testCase.assumeFail(sprintf(...
        ['CmdStan/MatlabStan not available for the StanModel command check ', ...
        '(%s). Skipping this layer.'], err.message));
end
end

function [found, value] = localGet(args, name)
found = false;
value = [];
for k = 1:2:numel(args)
    if strcmpi(args{k}, name)
        found = true;
        value = args{k+1};
        return;
    end
end
end

function tdir = localTempDirFor(testCase)
% Unique temp directory, removed after the test. CmdStan (and therefore the
% shared save-path validator) rejects whitespace, which tempname never produces.
tdir = tempname;
mkdir(tdir);
testCase.addTeardown(@() rmdir(tdir, 's'));
end

function T = localTinyTable()
% Minimal long-format table (id, meas) sufficient to reach the sampler-override
% validation in psyrat_computevarcomp, which runs before any CmdStan/preflight
% work. Kept trivial on purpose - the validation being tested does not inspect
% the data.
T = table([1; 1; 2; 2], [0.1; 0.2; 0.3; 0.4], 'VariableNames', {'id', 'meas'});
T.Properties.Description = 'sampler_override_validation';
end
