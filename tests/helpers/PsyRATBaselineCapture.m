classdef PsyRATBaselineCapture
    % Capture and compare seeded CmdStan estimation baselines (W0).
    %
    % Purpose: freeze the current posterior-draw output of each design path so
    % later changes (especially the W3 prior edits) can be shown to change
    % results only where intended. After the seed fix (commit 3a3aecd), seeded
    % CmdStan runs are bit-reproducible on a given machine/CmdStan version, so
    % a re-run at the same config reproduces the baseline exactly.
    %
    % CAVEAT: bit-reproducibility holds within a machine/CmdStan version. These
    % baselines are local regression references, not cross-platform golden
    % values, so they are not asserted in CI.
    %
    % Fixed capture configuration (documented, not production length): 4 chains,
    % 1000 warmup, 1000 sampling, seed 12345. Reproducibility is exact at any
    % length; the W3 prior-change check re-runs at this same config and diffs.
    %
    % Usage:
    %   root = '/path/to/PsyRAT_test';
    %   PsyRATBaselineCapture.capture(root,'ic_base', ...
    %       fullfile(root,'tests','baselines','ic_base.mat'));
    %   rpt = PsyRATBaselineCapture.compare(root,'ic_base', ...
    %       fullfile(root,'tests','baselines','ic_base.mat'));

    properties (Constant)
        Chains = 4;
        Warmup = 1000;
        Sampling = 1000;
        Seed = 12345;
    end

    methods (Static)
        function meta = capture(root, scenarioName, outFile)
            relout = PsyRATBaselineCapture.runScenario(root, scenarioName);
            meta = struct( ...
                'scenario', scenarioName, ...
                'chains', PsyRATBaselineCapture.Chains, ...
                'warmup', PsyRATBaselineCapture.Warmup, ...
                'sampling', PsyRATBaselineCapture.Sampling, ...
                'seed', PsyRATBaselineCapture.Seed, ...
                'psyrat_version', psyrat_defineversion(), ...
                'captured_utc', datestr(now, 'yyyy-mm-ddTHH:MM:SS')); %#ok<TNOW1,DATST>
            baselineDir = fileparts(outFile);
            if ~isempty(baselineDir) && ~isfolder(baselineDir)
                mkdir(baselineDir);
            end
            save(outFile, 'relout', 'meta', '-v7');
            fprintf('Captured baseline %s -> %s\n', scenarioName, outFile);
        end

        function report = compare(root, scenarioName, baselineFile)
            % Re-run the scenario and compare numeric draw fields of REL.out
            % to the saved baseline. Returns per-field maximum absolute
            % difference and an overall bit-identical flag.
            loaded = load(baselineFile, 'relout');
            fresh = PsyRATBaselineCapture.runScenario(root, scenarioName);
            report = struct('scenario', scenarioName, 'fields', struct(), ...
                'maxAbsDiff', 0, 'bitIdentical', true);
            report = PsyRATBaselineCapture.diffStruct(loaded.relout, fresh, '', report);
        end

        function relout = runScenario(root, scenarioName)
            cfg = PsyRATCmdStanIntegrationFactory.scenarios().(scenarioName);
            tbl = PsyRATCmdStanIntegrationFactory.makeScenarioTable(root, cfg);
            outdir = fullfile(root, 'tests', 'integration', 'cmdstan_artifacts');
            if ~isfolder(outdir)
                mkdir(outdir);
            end
            REL = psyrat_computevarcomp( ...
                'data', tbl, ...
                'chains', PsyRATBaselineCapture.Chains, ...
                'warmup', PsyRATBaselineCapture.Warmup, ...
                'sampling', PsyRATBaselineCapture.Sampling, ...
                'seed', PsyRATBaselineCapture.Seed, ...
                'verbose', 1, ...
                'showgui', 1, ...
                'cmdstanoutdir', outdir, ...
                'sserrvar', cfg.sserrvar, ...
                'diffest', cfg.diffest, ...
                'diffwpcov', cfg.diffwpcov, ...
                'diffrescor', cfg.diffrescor);
            relout = REL.out;
        end
    end

    methods (Static, Access = private)
        function report = diffStruct(a, b, prefix, report)
            % Recursively compare numeric leaves of two structs.
            fn = fieldnames(a);
            for i = 1:numel(fn)
                key = fn{i};
                label = key;
                if ~isempty(prefix)
                    label = [prefix '.' key];
                end
                if ~isfield(b, key)
                    report.fields.(matlab.lang.makeValidName(label)) = NaN;
                    report.bitIdentical = false;
                    continue;
                end
                av = a.(key);
                bv = b.(key);
                if isstruct(av) && isstruct(bv) && isscalar(av) && isscalar(bv)
                    report = PsyRATBaselineCapture.diffStruct(av, bv, label, report);
                elseif isnumeric(av) && isnumeric(bv) && isequal(size(av), size(bv))
                    d = max(abs(double(av(:)) - double(bv(:))), [], 'omitnan');
                    if isempty(d)
                        d = 0;
                    end
                    report.fields.(matlab.lang.makeValidName(label)) = d;
                    report.maxAbsDiff = max(report.maxAbsDiff, d);
                    if d ~= 0
                        report.bitIdentical = false;
                    end
                end
                % Non-numeric / non-struct fields (labels, etc.) are ignored.
            end
        end
    end
end
