classdef TestModularCutCopulaStanBoundary < PsyRATTestBase
    %TESTMODULARCUTCOPULASTANBOUNDARY The stage-1 Stan programs of the modular
    % cut-copula designs (analyses 19/20, family='gamma') must contain the exact
    % Gamma margins and NOTHING of the copula.
    %
    % WHY THIS IS A TEST AND NOT A COMMENT. The defining property of a cut
    % posterior is what is ABSENT from stage 1. If a future change moved the
    % probability transforms or the copula density back inside HMC - as an
    % "optimization", or by copying a line from the Gaussian concurrent sibling -
    % the model would still compile, still sample, and still produce plausible
    % reliabilities. It would simply be estimating a DIFFERENT posterior from the
    % one every output surface claims. No numeric regression test can see that,
    % because both posteriors are defensible; only the boundary itself can.
    %
    % HOW THIS PINS THE LIVE BUILDER, not just a stale file. The assertions below
    % read the committed goldens, and TestStanModelSnapshots separately requires
    % those goldens to equal what the builders emit right now (it regenerates and
    % compares byte for byte). The two together are what pin the live output;
    % neither alone would.
    %
    % Ported from the reference bundle's own static validation
    % (tests/run_static_validation.R, check 1), which asserts the same absences
    % against the same design.

    properties (Constant)
        % Goldens for both designs, in both the plain and threaded forms.
        Goldens = {
            'dynrel_gamma_ls_diff_sserr_rescor_onefacet.stan'
            'dynrel_gamma_ls_diff_sserr_rescor_onefacet__threaded.stan'
            'dynrel_gamma_ls_diff_sserr_rescor_trt_onefacet.stan'
            'dynrel_gamma_ls_diff_sserr_rescor_trt_onefacet__threaded.stan'
            };

        % Tokens whose presence would mean the copula had re-entered HMC. Each is
        % a real construct of the FULL-JOINT concurrent models, so none of them
        % can appear here by coincidence:
        %   gamma_lcdf/gamma_lccdf - the probability integral transform
        %   inv_Phi                - the inverse-normal step after it
        %   rho_e                  - the copula correlation as a sampled parameter
        %   Lrescor/multi_normal   - the Gaussian sibling's residual coupling
        %   biv_                   - the bivariate copula density helpers
        Forbidden = {'gamma_lcdf', 'gamma_lccdf', 'inv_Phi', 'rho_e', ...
            'Lrescor', 'multi_normal', 'biv_'};
    end

    methods (Test)

        function testStageOneContainsNoCopulaConstruct(testCase)
            for k = 1:numel(TestModularCutCopulaStanBoundary.Goldens)
                name = TestModularCutCopulaStanBoundary.Goldens{k};
                code = localReadGolden(testCase, name);
                for f = TestModularCutCopulaStanBoundary.Forbidden
                    testCase.verifyEmpty(strfind(code, f{1}), sprintf( ...
                        ['%s contains "%s". The stage-1 margins model must not '...
                        'carry any copula construct: the coupling is estimated '...
                        'in the separate cut stage, and moving it back inside '...
                        'HMC silently changes which posterior is being '...
                        'sampled.'], name, f{1}));
                end
            end
        end

        function testStageOneContainsTheExactGammaMargins(testCase)
            % The other half of the claim: the margins must be the moment-matched
            % scaled chi-square, not an approximation, and both components must
            % actually enter the likelihood.
            required = {
                'real scaled_chi_square_lpdf(real y, real mu, real nu)'
                'nu1 = 2 * square(mu1 ./ sigma1);'
                'nu2 = 2 * square(mu2 ./ sigma2);'
                'scaled_chi_square_lpdf(meas1[n] | mu1[n], nu1[n])'
                'scaled_chi_square_lpdf(meas2[n] | mu2[n], nu2[n])'
                };
            for k = 1:numel(TestModularCutCopulaStanBoundary.Goldens)
                name = TestModularCutCopulaStanBoundary.Goldens{k};
                code = localReadGolden(testCase, name);
                for r = required'
                    testCase.verifyNotEmpty(strfind(code, r{1}), ...
                        sprintf('%s is missing: %s', name, r{1}));
                end
            end
        end

        function testOffsetsEnterAsDataAndInterceptsAreDeviations(testCase)
            % The offsets are what make the intercepts deviations. If they ever
            % stopped being data, every stored quantity that adds them back
            % downstream would double-count or drop a centring constant.
            for k = 1:numel(TestModularCutCopulaStanBoundary.Goldens)
                name = TestModularCutCopulaStanBoundary.Goldens{k};
                code = localReadGolden(testCase, name);
                testCase.verifyNotEmpty(strfind(code, 'vector[2] mu_offset;'));
                testCase.verifyNotEmpty(strfind(code, 'vector[2] log_sigma_offset;'));
                testCase.verifyNotEmpty(strfind(code, ...
                    'rep_vector(mu_offset[1] + b[1], NOBS)'));
                testCase.verifyNotEmpty(strfind(code, ...
                    'rep_vector(log_sigma_offset[1] + b_sigma[1], NOBS)'));
                % They must be declared in the DATA block, not the parameters
                % block: find the data block and check they fall inside it.
                dataStart = strfind(code, sprintf('data {\n'));
                paramStart = strfind(code, sprintf('parameters {\n'));
                offsetAt = strfind(code, 'vector[2] mu_offset;');
                testCase.assertNotEmpty(dataStart);
                testCase.assertNotEmpty(paramStart);
                testCase.verifyTrue(offsetAt(1) > dataStart(1) && ...
                    offsetAt(1) < paramStart(1), sprintf( ...
                    '%s declares mu_offset outside the data block.', name));
            end
        end

        function testResidualSdCarriesNoFacetEffect(testCase)
            % The person-specific residual depends on this: sigma varies across
            % PEOPLE and with the dimension predictors, never across trials or
            % occasions. A facet term on sigma would make the stored per-person
            % residual an average over that facet rather than a single-trial
            % conditional variance.
            for k = 1:numel(TestModularCutCopulaStanBoundary.Goldens)
                name = TestModularCutCopulaStanBoundary.Goldens{k};
                code = localReadGolden(testCase, name);
                lines = regexp(code, '\r\n|\n|\r', 'split');
                sigmaUpdates = lines(contains(lines, 'sigma1[n] +=') | ...
                    contains(lines, 'sigma2[n] +='));
                testCase.verifyNotEmpty(sigmaUpdates, ...
                    sprintf('%s has no residual accumulation at all.', name));
                for u = sigmaUpdates
                    % Only the person block may touch it.
                    testCase.verifyNotEmpty(strfind(u{1}, 'r_id['), sprintf( ...
                        '%s: residual update does not read the person block: %s', ...
                        name, strtrim(u{1})));
                    for facet = {'r_trl[', 'r_occ[', 'r_tid[', 'r_oid[', 'r_to['}
                        testCase.verifyEmpty(strfind(u{1}, facet{1}), sprintf( ...
                            ['%s: a facet effect (%s) reaches the residual SD. '...
                            'The residual must vary across people and the '...
                            'dimension predictors only.'], name, facet{1}));
                    end
                end
            end
        end

        function testTwoFacetLoadsEveryFacetOnTheMean(testCase)
            % The converse of the test above: all five crossed facets must reach
            % the mean, on BOTH components. A dropped facet would silently move
            % its variance into the residual.
            for name = {'dynrel_gamma_ls_diff_sserr_rescor_trt_onefacet.stan', ...
                    'dynrel_gamma_ls_diff_sserr_rescor_trt_onefacet__threaded.stan'}
                code = localReadGolden(testCase, name{1});
                for facet = {'r_trl', 'r_occ', 'r_tid', 'r_oid', 'r_to'}
                    testCase.verifyNotEmpty(strfind(code, ...
                        sprintf('mu1[n] += %s[', facet{1})), sprintf( ...
                        '%s: %s does not reach component 1''s mean.', ...
                        name{1}, facet{1}));
                    testCase.verifyNotEmpty(strfind(code, ...
                        sprintf('mu2[n] += %s[', facet{1})), sprintf( ...
                        '%s: %s does not reach component 2''s mean.', ...
                        name{1}, facet{1}));
                end
            end
        end

        function testEveryNameTheExtractionReadsIsEmitted(testCase)
            % The estimation arm reads these by name after the fit. A missing one
            % would abort AFTER sampling completes - the failure mode that costs
            % a multi-hour run, and the one this repo has already been bitten by.
            % The two lists are written out separately rather than one extending
            % the other, because the trial correlation is NOT named the same in
            % both: the one-facet model calls it cor_i (matching the gamma
            % one-facet family), while the two-facet model calls it cor_t and
            % reserves the -i suffix for nothing. Deriving one list from the
            % other hides that.
            oneFacet = {'r_id', 'r_trl', 'cor_id', 'cor_i', 'sd_id', 'sd_trl', ...
                'b_dim1', 'b_dim2', 'b_sigma_dim1', 'b_sigma_dim2'};
            twoFacet = {'r_id', 'r_trl', 'r_occ', 'r_tid', 'r_oid', 'r_to', ...
                'cor_id', 'sd_id', ...
                'sd_trl', 'sd_occ', 'sd_tid', 'sd_oid', 'sd_to', ...
                'cor_t', 'cor_o', 'cor_pt', 'cor_po', 'cor_ot', ...
                'b_dim1', 'b_dim2', 'b_sigma_dim1', 'b_sigma_dim2'};

            localVerifyNamesDeclared(testCase, ...
                'dynrel_gamma_ls_diff_sserr_rescor_onefacet.stan', oneFacet);
            localVerifyNamesDeclared(testCase, ...
                'dynrel_gamma_ls_diff_sserr_rescor_trt_onefacet.stan', twoFacet);
        end

    end
end
% Convergence codes 57/58 are NOT checked here. psyrat_convprefix is a local
% function of psyrat_computevarcomp.m and cannot be called from outside it, so
% the repo already solves that problem the only way available - by parsing the
% source text and cross-checking it against these same goldens, in
% TestConvPrefixCoverage. The 57/58 entries live in that file's golden map.

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function code = localReadGolden(testCase, name)
%LOCALREADGOLDEN Reads one committed golden Stan file.
path = fullfile(testCase.projectRoot(), 'tests', 'golden', 'stan', name);
testCase.assertTrue(isfile(path), sprintf('Missing golden: %s', name));
code = fileread(path);
end

function localVerifyNamesDeclared(testCase, golden, names)
%LOCALVERIFYNAMESDECLARED Each name must appear as a declared identifier.
code = localReadGolden(testCase, golden);
for n = names
    pattern = ['(?<![\w])' regexptranslate('escape', n{1}) '(?![\w])'];
    testCase.verifyNotEmpty(regexp(code, pattern, 'once'), sprintf( ...
        '%s never declares %s, but the estimation arm extracts it.', ...
        golden, n{1}));
end
end
