classdef TestNonconcurrentDodFixtureIntegrity < PsyRATTestBase
    %TESTNONCONCURRENTDODFIXTUREINTEGRITY Integrity of the frozen nonconcurrent
    % difference-of-differences reference fixtures themselves, independent of
    % any calculator (PsyRAT analyses 28/29).
    %
    % These fixtures are different in kind from the rest of the accuracy
    % material: they come from the owner's INDEPENDENT R reference
    % implementation (person_specific_dynamic_nonconcurrent_dod_modular), so a
    % MATLAB test matching them is externally grounded rather than
    % self-referential. That standing depends on the properties this class
    % pins, mirroring TestModularCutCopulaFixtureIntegrity:
    %
    %   1 PROVENANCE. The fixture directory hashes every reference-bundle file
    %     the generator sources - all FIVE, plus the R version and RNG kinds -
    %     so a changed bundle is detectable from the repo alone. Where the
    %     bundle is present, the recorded hashes are verified against it
    %     (BUNDLE-GATED: skips Incomplete, not Failed, when absent).
    %
    %   2 DESIGN SEPARATION. The stage-3 'standardized' D-study block must
    %     genuinely differ from the 'actual' block, and the two standardized
    %     trial counts must differ from each other - otherwise the per-cell
    %     divisor slots are unpinned and a swapped divisor passes (the lesson
    %     recorded in the copula fixture set's 2026-08-11 follow-up).
    %
    %   3 THE NONCONCURRENT BOUNDARY. The frozen six-pair residual audit must
    %     enumerate all six pairs at exactly zero, and the frozen per-person
    %     DoD residual variance must equal the sum of the four per-cell
    %     residual variances - the identity that holds BECAUSE the six
    %     residual covariances are zero (reference bundle check 8).
    %
    % SCOPE NOTE: the *_exact_shared columns are frozen but NOT consumed by the
    % toolbox (harmonic-only scope, owner ruling 2026-08-12). One guard below
    % keeps them separated from the harmonic columns so they retain pinning
    % value for a later increment.

    methods (Test)

        %% -- provenance -------------------------------------------------------

        function testProvenanceHashesAllFiveSourcedReferenceFiles(testCase)
            %No bundle needed: this pins the STRUCTURE of the provenance file,
            %so the repo alone proves the hash list is complete.
            prov = localReadProvenance();
            want = {'R/common.R'; 'R/model_generation.R'; ...
                'R/independence_stage.R'; 'R/reliability.R'; ...
                'tests/simulation_helpers.R'};
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
            %recorded hash must match the live file: an unexpected fixture
            %diff means the reference moved and the MATLAB tolerances need
            %re-examination, and this is the check that notices the move
            %BEFORE a regeneration.
            root = localBundleRoot();
            testCase.assumeTrue( ...
                isfile(fullfile(root, 'R', 'independence_stage.R')), ...
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

        function testStandardizedDesignDiffersFromActualOnefacet(testCase)
            one = localReadFixture('nonconcurrent_dod_onefacet_stage3_expected.csv');
            [act, std_] = localSplitDesigns(testCase, one);

            testCase.verifyNotEqual(std_.n_i_error(1), act.n_i_error(1), ...
                'Standardized n_i(error) must differ from the actual design.');
            testCase.verifyNotEqual(std_.n_i_correct(1), act.n_i_correct(1), ...
                'Standardized n_i(correct) must differ from the actual design.');
            testCase.verifyNotEqual(std_.n_i_error(1), std_.n_i_correct(1), ...
                ['The two standardized trial counts must differ from EACH ', ...
                'OTHER, or a per-cell divisor swap passes.']);
            testCase.verifyTrue(all(act.G_dod ~= std_.G_dod), ...
                'Every standardized coefficient must differ from its actual twin.');
            %...while the components must NOT differ: the design selects
            %divisors only, so a component that moves with the design would
            %mean the selector leaked into the G-study
            testCase.verifyEqual(std_.var_dod_p, act.var_dod_p, ...
                'AbsTol', 0, 'Components are design-invariant.');
        end

        function testStandardizedDesignDiffersFromActualTwofacet(testCase)
            two = localReadFixture('nonconcurrent_dod_twofacet_stage3_expected.csv');
            [act, std_] = localSplitDesigns(testCase, two);

            testCase.verifyNotEqual(std_.n_i_error(1), act.n_i_error(1), ...
                'Standardized n_i(error) must differ from the actual design.');
            testCase.verifyNotEqual(std_.n_o_error(1), act.n_o_error(1), ...
                'Standardized n_o(error) must differ from the actual design.');
            testCase.verifyNotEqual(std_.n_o_error(1), std_.n_o_correct(1), ...
                ['The two standardized occasion counts must differ from EACH ', ...
                'OTHER, or a per-cell occasion-divisor swap passes.']);
            testCase.verifyTrue(all(act.G_CES_dod ~= std_.G_CES_dod), ...
                'Every standardized coefficient must differ from its actual twin.');
            testCase.verifyEqual(std_.var_dod_p, act.var_dod_p, ...
                'AbsTol', 0, 'Components are design-invariant.');
        end

        function testExactSharedColumnsDifferFromHarmonic(testCase)
            %The exact-shared sensitivity is frozen but out of toolbox scope
            %(harmonic-only ruling). Keeping the two column families separated
            %preserves their pinning value: n_common was deliberately set
            %BELOW min(n_i) at generation, so identical columns would mean a
            %regeneration collapsed the two count rules.
            one = localReadFixture('nonconcurrent_dod_onefacet_stage3_expected.csv');
            testCase.verifyTrue( ...
                any(one.relative_error_dod ~= one.relative_error_dod_exact_shared), ...
                ['The exact-shared error columns must differ from the ', ...
                'harmonic-primary columns on at least one row.']);
        end

        %% -- the nonconcurrent boundary ---------------------------------------

        function testResidualAuditEnumeratesSixZeroPairs(testCase)
            for tag = {'onefacet', 'twofacet'}
                audit = localReadFixture( ...
                    sprintf('nonconcurrent_dod_%s_residual_audit.csv', tag{1}));
                testCase.verifyEqual(height(audit), 6, ...
                    sprintf('%s: the audit must enumerate all six pairs.', tag{1}));
                testCase.verifyEqual(audit.conditional_residual_covariance, ...
                    zeros(6, 1), 'AbsTol', 0, sprintf( ...
                    '%s: every conditional residual covariance is exactly zero.', ...
                    tag{1}));
                testCase.verifyTrue(all(~localAsLogical(audit.estimated)), sprintf( ...
                    '%s: no residual dependence parameter is estimated.', tag{1}));
            end
        end

        function testDodResidualVarianceIsTheSumOfFour(testCase)
            %Reference bundle check 8: with all six residual covariances at
            %zero the DoD residual variance is EXACTLY the sum of the four
            %per-cell residual variances - no covariance corrections. The
            %tolerance covers only the CSV decimal round-trip (each column is
            %independently printed at ~15 significant digits).
            for tag = {'onefacet', 'twofacet'}
                t = localReadFixture( ...
                    sprintf('nonconcurrent_dod_%s_stage3_expected.csv', tag{1}));
                total = t.residual_var_error_theta_person + ...
                    t.residual_var_correct_theta_person + ...
                    t.residual_var_error_alpha_person + ...
                    t.residual_var_correct_alpha_person;
                testCase.verifyEqual(t.residual_var_dod_person, total, ...
                    'AbsTol', 1e-10, sprintf( ...
                    '%s: residual_var_dod_person must be the sum of four.', ...
                    tag{1}));
            end
        end

    end
end

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function t = localReadFixture(name)
here = fileparts(mfilename('fullpath'));
t = readtable(fullfile(here, 'baselines', 'nonconcurrent_dod', name), ...
    'VariableNamingRule', 'preserve');
end

function t = localReadProvenance()
%LOCALREADPROVENANCE Force every column to char: an md5 that happened to be all
%digits around a single 'e' would otherwise be parsed as a number.
here = fileparts(mfilename('fullpath'));
f = fullfile(here, 'baselines', 'nonconcurrent_dod', ...
    'nonconcurrent_dod_provenance.csv');
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

function tf = localAsLogical(col)
%LOCALASLOGICAL Normalize a CSV boolean column: R's write.csv emits the
%strings TRUE/FALSE, which readtable imports as cellstr; a future import
%change to logical or numeric must keep working too.
if islogical(col)
    tf = col;
elseif isnumeric(col)
    tf = col ~= 0;
else
    tf = strcmpi(cellstr(col), 'TRUE');
end
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
%   the MPSDOD_REFERENCE_DIR environment variable, else the Desktop default
%   (the canonical bundle; working copies elsewhere are temporary).
root = getenv('MPSDOD_REFERENCE_DIR');
if isempty(root)
    root = fullfile(localHome(), 'Desktop', 'stan_chisquare', ...
        'person_specific_dynamic_nonconcurrent_dod_modular');
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
msg = ['The nonconcurrent DoD reference bundle is not present (set ', ...
    'MPSDOD_REFERENCE_DIR or place it at ~/Desktop/stan_chisquare/', ...
    'person_specific_dynamic_nonconcurrent_dod_modular). Hash verification ', ...
    'skipped.'];
end
