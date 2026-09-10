classdef TestModularCutCopulaFixtureIntegrity < PsyRATTestBase
    %TESTMODULARCUTCOPULAFIXTUREINTEGRITY Integrity of the frozen modular
    % cut-copula reference fixtures themselves, independent of any calculator.
    %
    % These fixtures are different in kind from the rest of the accuracy
    % material: they come from the owner's INDEPENDENT R + CmdStan reference
    % implementation, so a MATLAB test matching them is externally grounded
    % rather than self-referential. That standing depends on two properties
    % this class pins:
    %
    %   1 PROVENANCE. The fixture directory hashes every reference-bundle file
    %     the generator sources - all FIVE, plus the R version and RNG kinds -
    %     so a changed bundle is detectable from the repo alone. Where the
    %     bundle is present, the recorded hashes are verified against it
    %     (BUNDLE-GATED, mirroring TestModularCutCopulaExternalReplay: skips
    %     Incomplete, not Failed, when the bundle is absent).
    %
    %   2 DESIGN SEPARATION. The stage-3 'standardized' D-study block must
    %     genuinely differ from the 'actual' block. The original two-facet
    %     fixture set n' = (4, 2), which coincided exactly with the simulated
    %     design's (I, O), so the two blocks were byte-identical and the
    %     design-selector divisors were unpinned - a swapped or mis-assigned
    %     divisor would have passed (2026-08-11 post-merge follow-up). This
    %     guard fails if a regeneration ever collapses them again.

    methods (Test)

        %% -- provenance -------------------------------------------------------

        function testProvenanceHashesAllFiveSourcedReferenceFiles(testCase)
            %No bundle needed: this pins the STRUCTURE of the provenance file,
            %so the repo alone proves the hash list is complete.
            prov = localReadProvenance();
            want = {'R/common.R'; 'R/model_generation.R'; 'R/copula_stage.R'; ...
                'R/reliability.R'; 'tests/simulation_helpers.R'};
            for i = 1:numel(want)
                row = strcmp(prov.file, want{i});
                testCase.verifyEqual(sum(row), 1, ...
                    sprintf('%s must be hashed exactly once.', want{i}));
                md5 = prov.md5{row};
                testCase.verifyTrue( ...
                    ~isempty(regexp(md5, '^[0-9a-f]{32}$', 'once')), ...
                    sprintf('%s: not a well-formed md5: ''%s''', want{i}, md5));
            end
            %the environment rows ride along so a regeneration under a
            %different R shows up in the diff, not just in the numbers
            testCase.verifyEqual(sum(strcmp(prov.file, 'R_version')), 1);
            testCase.verifyEqual(sum(strcmp(prov.file, 'RNG_kinds')), 1);
        end

        function testProvenanceHashesMatchTheReferenceBundle(testCase)
            %BUNDLE-GATED. Where the reference bundle is present, every
            %recorded hash must match the live file: the generator's rule is
            %that an unexpected fixture diff means the reference moved and the
            %MATLAB tolerances need re-examination, and this is the check that
            %notices the move BEFORE a regeneration.
            root = localBundleRoot();
            testCase.assumeTrue( ...
                isfile(fullfile(root, 'R', 'copula_stage.R')), ...
                localSkipMessage());

            prov = localReadProvenance();
            hashed = prov(~cellfun(@isempty, prov.md5), :);
            testCase.verifyEqual(height(hashed), 5, ...
                'Exactly the five sourced files carry hashes.');
            for i = 1:height(hashed)
                parts = strsplit(hashed.file{i}, '/');
                f = fullfile(root, parts{:});
                testCase.verifyTrue(isfile(f), ...
                    sprintf('Bundle present but %s is missing.', hashed.file{i}));
                testCase.verifyEqual(localMd5(f), hashed.md5{i}, ...
                    sprintf(['%s: bundle file no longer matches the hash the ' ...
                    'fixtures were generated from.'], hashed.file{i}));
            end
        end

        %% -- design separation ------------------------------------------------

        function testStandardizedDesignDiffersFromActualTwofacet(testCase)
            two = localReadFixture('modular_cut_copula_twofacet_stage3_expected.csv');
            [act, std_] = localSplitDesigns(testCase, two);

            testCase.verifyNotEqual(std_.n_i_prime(1), act.n_i_prime(1), ...
                'Standardized n_i'' must differ from the actual design.');
            testCase.verifyNotEqual(std_.n_o_prime(1), act.n_o_prime(1), ...
                'Standardized n_o'' must differ from the actual design.');
            testCase.verifyNotEqual(std_.n_i_prime(1), std_.n_o_prime(1), ...
                'n_i'' and n_o'' must differ from EACH OTHER, or a divisor swap passes.');
            testCase.verifyTrue(all(act.G_CES_delta ~= std_.G_CES_delta), ...
                'Every standardized coefficient must differ from its actual twin.');
            %...while the components must NOT differ: the design selects
            %divisors only, so a component that moves with the design would
            %mean the selector leaked into the G-study
            testCase.verifyEqual(std_.var_delta_p, act.var_delta_p, ...
                'AbsTol', 0, 'Components are design-invariant.');
        end

        function testStandardizedDesignDiffersFromActualOnefacet(testCase)
            one = localReadFixture('modular_cut_copula_onefacet_stage3_expected.csv');
            [act, std_] = localSplitDesigns(testCase, one);

            testCase.verifyNotEqual(std_.n_i_prime(1), act.n_i_prime(1), ...
                'Standardized n_i'' must differ from the actual design.');
            testCase.verifyTrue(all(act.G_delta ~= std_.G_delta), ...
                'Every standardized coefficient must differ from its actual twin.');
            testCase.verifyEqual(std_.var_delta_p, act.var_delta_p, ...
                'AbsTol', 0, 'Components are design-invariant.');
        end

    end
end

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function t = localReadFixture(name)
here = fileparts(mfilename('fullpath'));
t = readtable(fullfile(here, 'baselines', 'modular_cut_copula', name), ...
    'VariableNamingRule', 'preserve');
end

function t = localReadProvenance()
%LOCALREADPROVENANCE Force every column to char: an md5 that happened to be all
%digits around a single 'e' would otherwise be parsed as a number.
here = fileparts(mfilename('fullpath'));
f = fullfile(here, 'baselines', 'modular_cut_copula', ...
    'modular_cut_copula_provenance.csv');
opts = detectImportOptions(f, 'VariableNamingRule', 'preserve');
opts = setvartype(opts, 'char');
t = readtable(f, opts);
end

function [act, std_] = localSplitDesigns(testCase, t)
%LOCALSPLITDESIGNS The two D-study blocks, aligned row-for-row so each
%standardized row is compared to its own (draw, participant) actual twin.
act  = sortrows(t(strcmp(t.design, 'actual'), :), {'posterior_row', 'p_index'});
std_ = sortrows(t(strcmp(t.design, 'standardized'), :), {'posterior_row', 'p_index'});
testCase.assertEqual(height(act), height(std_), ...
    'The two design blocks must pair off exactly.');
testCase.assertGreaterThan(height(act), 0, 'No design rows found.');
end

function h = localMd5(f)
%LOCALMD5 MD5 of a file's bytes via the JVM, so no toolbox is required.
fid = fopen(f, 'r');
assert(fid > 0, 'localMd5: cannot open %s', f);
cleaner = onCleanup(@() fclose(fid));
bytes = fread(fid, inf, '*uint8');
md = java.security.MessageDigest.getInstance('MD5');
md.update(typecast(bytes(:), 'int8'));
h = lower(sprintf('%02x', typecast(md.digest(), 'uint8')));
end

function root = localBundleRoot()
%LOCALBUNDLEROOT Reference bundle location, mirroring the fixture generator:
%   the MPSDC_REFERENCE_DIR environment variable, else the Desktop default.
root = getenv('MPSDC_REFERENCE_DIR');
if isempty(root)
    root = fullfile(localHome(), 'Desktop', 'stan_chisquare', ...
        'person_specific_dynamic_concurrent_modular');
end
end

function home = localHome()
%LOCALHOME The user home directory without the Statistics/OS toolboxes.
if ispc
    home = getenv('USERPROFILE');
else
    home = getenv('HOME');
end
end

function msg = localSkipMessage()
%LOCALSKIPMESSAGE One place for the assume text, matching the external replay.
msg = ['The modular cut-copula reference bundle is not present (set ', ...
    'MPSDC_REFERENCE_DIR or place it at ~/Desktop/stan_chisquare/', ...
    'person_specific_dynamic_concurrent_modular). Hash verification skipped.'];
end
