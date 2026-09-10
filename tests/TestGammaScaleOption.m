classdef TestGammaScaleOption < matlab.unittest.TestCase
    %TESTGAMMASCALEOPTION Pins the 'gammascale' option (psyrat_computevarcomp),
    % which selects which quantity the Gamma / scaled-chi-square DYNAMIC
    % (dimension-conditional) models put the linear predictor and the person
    % random effect on:
    %   gammascale = 1 -> log nu    (dispersion; the DEFAULT, unchanged behavior)
    %   gammascale = 2 -> log sigma (residual SD; the location-scale form)
    %
    % WHY THE OPTION EXISTS. Under gammascale = 1 the conditional variance is
    % 2*mu^2/nu, so a dimension slope on log-nu blends a mean-level effect with a
    % dispersion effect and cannot be read as a statement about precision alone.
    % gammascale = 2 makes b_sigma a pure residual-SD slope. See
    % documentation/gamma_scale_submodel_decision.md.
    %
    % The default is pinned here as well as the new branch: without that test,
    % inverting the default would leave the whole suite green while silently
    % changing the estimand of every existing gamma dynrel run.
    %
    % No CmdStan required: the PSYRAT_EMIT_MODEL_FILE hook writes the generated
    % Stan and aborts with PsyRAT:emitModelOnly before any sampling, so these
    % tests assert on the emitted model text.

    methods (Test)

        function testDefaultIsLocationScaleForMovableDesigns(testCase)
            % DEFAULT (no 'gammascale' passed) on a design that SUPPORTS
            % location-scale must now emit the location-scale model.
            %
            % This replaces the former testDefaultIsLogNuDispersion, which pinned
            % the opposite. Rulings 33-36 (2026-08-07, gamma_scale_submodel_
            % decision.md section G) flipped the default for the eleven movable
            % designs once every one of them had been validated. The old test was
            % the tripwire ruling 31 asked for; it fired exactly as intended, and
            % is inverted rather than deleted so the guarantee survives.
            stanText = localEmit(testCase, localDynrelTable(), ...
                {'family','gamma','dynrel',2});
            testCase.verifySubstring(stanText, 'real Intercept_sigma;', ...
                'The gamma dynrel default must now be the location-scale model.');
            testCase.verifyEmpty(strfind(stanText, 'real Intercept_nu;'), ...
                'The default must NOT emit the log-nu scale submodel.'); %#ok<STREMP>
        end

        function testDefaultStaysLogNuForNonMovableDesigns(testCase)
            % The other half of the flip, and the half that makes it CONDITIONAL
            % rather than global. Analyses 1-5 have no location-scale variant at
            % all, so their default must still be log-nu; a wholesale flip would
            % emit a model here for which no extractor exists.
            %
            % Without this test the conditional could be silently widened to a
            % blanket default and nothing in the offline suite would notice until
            % a run failed downstream.
            stanText = localEmit(testCase, localPlainGammaTable(), ...
                {'family','gamma'});
            testCase.verifySubstring(stanText, 'real Intercept_nu;', ...
                'Analysis 1 has no location-scale variant; it must stay log-nu.');
            testCase.verifyEmpty(strfind(stanText, 'real Intercept_sigma;'), ...
                'A non-movable design must not be moved by the new default.'); %#ok<STREMP>
        end

        function testDefaultStaysLogNuForConcurrentDesigns(testCase)
            % The third arm of the condition, and the one guarding the QUIETEST
            % failure mode in this feature.
            %
            % Analysis 7 IS movable, but its concurrent branch calls a copula
            % builder whose signature is (priors, threaded) - it has no
            % gammascale parameter and cannot carry a location-scale request. If
            % the default flip reached a concurrent run, the concurrency guard
            % could NOT save it: that guard runs at option-parse time, before the
            % context-aware resolution, so it never sees a defaulted value. The
            % run would stamp REL.gammascale = 2 while the emitted model was
            % log-nu, and the downstream extraction would read log-nu draws as
            % location-scale ones. Wrong numbers, no error.
            %
            % WHAT THIS TEST DOES AND DOES NOT COVER - measured, not assumed.
            % It pins that a concurrent run EMITS log-nu. It does NOT protect the
            % diffrescor conjunct: with that conjunct deleted this test still
            % passes, because the copula builder emits log-nu either way and only
            % the REL.gammascale stamp diverges. Two earlier drafts of this
            % comment were wrong here - first claiming the failure mode was a
            % hard error, then claiming this test was the tripwire for it.
            % The actual guard is testResolutionBlockKeepsAllThreeGuardConditions,
            % which reads the source because nothing behavioral can see this.
            % NAMING: assert on b_nu, NOT Intercept_nu. In a DIFFERENCE design the
            % log-nu submodel is a vector of per-event intercepts called b_nu;
            % Intercept_nu is the non-difference dynrel spelling and never appears
            % here. An earlier draft of this test asserted Intercept_nu and failed
            % against a model that was correctly log-nu throughout.
            stanText = localEmit(testCase, localGammaDiffTable(), ...
                {'family','gamma','diffest',2,'diffrescor',2});
            testCase.verifySubstring(stanText, 'b_nu', ...
                'A concurrent gamma difference run must still be log-nu.');
            testCase.verifyEmpty(strfind(stanText, 'b_sigma'), ...
                'The default must not move a design that rejects gammascale=2.'); %#ok<STREMP>
        end

        function testGammascaleTwoEmitsResidualSdSubmodel(testCase)
            % gammascale=2 must swap the scale submodel wholesale: the log-nu
            % intercept, slopes and person effect are replaced by their
            % residual-SD counterparts. Assert BOTH directions so a partial
            % rename cannot pass.
            stanText = localEmit(testCase, localDynrelTable(), ...
                {'family','gamma','dynrel',2,'gammascale',2});

            testCase.verifySubstring(stanText, 'real Intercept_sigma;');
            testCase.verifySubstring(stanText, 'vector[KDIM] b_sigma;');
            testCase.verifySubstring(stanText, 'vector[NSUB] ind_sd;');

            testCase.verifyEmpty(strfind(stanText, 'real Intercept_nu;'), ...
                'gammascale=2 must not retain the log-nu intercept.');
            testCase.verifyEmpty(strfind(stanText, 'vector[KDIM] b_nu;'), ...
                'gammascale=2 must not retain the log-nu dimension slopes.');
            testCase.verifyEmpty(strfind(stanText, 'vector[NSUB] ind_nu;'), ...
                'gammascale=2 must not retain the person log-nu effect.');
        end

        function testGammascaleTwoDerivesNuFromMeanAndSd(testCase)
            % The whole parameterization rests on nu = 2*(mu/sigma)^2, which is
            % what makes Var(Y|mu,sigma) exactly sigma^2. Pin the derivation
            % line, and pin that the family still uses the SINGLE existing
            % scaled_chi_square_lpdf density rather than gaining a second one.
            stanText = localEmit(testCase, localDynrelTable(), ...
                {'family','gamma','dynrel',2,'gammascale',2});

            testCase.verifySubstring(stanText, 'nu = 2 * square(mu ./ sigma);', ...
                'gammascale=2 must derive nu from the mean and residual SD.');
            testCase.verifyEqual( ...
                numel(strfind(stanText, 'real scaled_chi_square_lpdf')), 1, ...
                'The gamma family must keep exactly one density function.');
        end

        function testGammascaleTwoRejectedOutsideGammaFamily(testCase)
            % The option names a Gamma-family parameterization. Under the
            % Gaussian family it is meaningless and must error rather than be
            % silently ignored.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDynrelTable(), 'dynrel',2, 'gammascale',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:gammascale');
        end

        function testGammascaleTwoRejectedForNonDynrelGamma(testCase)
            % A plain one-facet gamma run (analysis 1) with gammascale=2 must
            % error rather than silently fit a different estimand than every
            % other case-1 result.
            %
            % RATIONALE CORRECTED 2026-08-04. This comment used to read "Only the
            % dynrel gamma designs carry dimension slopes, so only they have a
            % blended slope to clean up." That is the scoping argument owner
            % ruling H overturned: a static design has no dimension SLOPE, but its
            % event/cell and person log-nu CONTRASTS blend amplitude with absolute
            % residual SD just as a slope does (finding S16). Analyses 6 and 25 are
            % now accepted, and the remaining static difference designs are in
            % scope but unbuilt.
            %
            % Analyses 1-5 remain rejected ON THE MERITS, not as unbuilt work:
            % they are single-mean designs in which every observed-scale component
            % carries a common exp(2*alpha) factor that divides out of every
            % coefficient, so the parameterization genuinely changes nothing
            % reported. That is why this test keeps asserting an error here while
            % the ic_gamma_ls_sserr / ic_gamma_ls_trt_sserr snapshot scenarios
            % assert the opposite for analyses 6 and 25.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDynrelTable(), 'family','gamma', 'gammascale',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:gammascale');
        end

        function testGammascaleTwoRejectedWithGlobalDispersion(testCase)
            % dispersion=2 fixes a single global nu; under gammascale=2 nu is
            % derived, so there is nothing left to fix and the combination must
            % be rejected rather than have one option silently win.
            %
            % For THIS scenario it is rejected by the EXISTING dispersion guard,
            % not a gammascale-side one, and the identifier asserted here is
            % 'varargin:dispersion' for that reason: the run is dynrel (diffest=1),
            % so isGroupConcurrentGammaDiff is false and the
            % 'dispersion==2 && ~isGroupConcurrentGammaDiff' guard fires first.
            %
            % RATIONALE CORRECTED 2026-08-04 (increment 5). This comment used to
            % add that the two options "are mutually exclusive by construction and
            % a gammascale-side guard would be unreachable code in a fragile
            % core." That stopped being true when owner ruling H put the static
            % DIFFERENCE designs (analyses 7/8, diffest=2) in scope: for a
            % concurrent case-7 run isGroupConcurrentGammaDiff is TRUE, so the
            % dispersion guard cannot fire and nothing would have caught
            % gammascale=2 + diffrescor=2. A real gammascale-side guard now exists
            % and is pinned by the two tests below. The scenario HERE is
            % unaffected and still routes through the dispersion guard, which is
            % why this test is unchanged apart from its rationale.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDynrelTable(), 'family','gamma', 'dynrel',2, ...
                'gammascale',2, 'dispersion',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:dispersion');
        end

        function testGammascaleTwoRejectedForConcurrentGroupDifference(testCase)
            % The gap the dispersion guard CANNOT cover (increment 5). For a
            % group-level concurrent gamma difference (analysis 7, diffest=2 +
            % diffrescor=2 + sserrvar=1) isGroupConcurrentGammaDiff is true, so
            % 'dispersion' defaults to 2 and its guard is inert. Before the
            % gammascale-side guard, this combination reached
            % psyrat_build_stan_diff_copula_fixnu_chisq - which has no
            % location-scale branch - and silently emitted a log-nu model for a run
            % the user explicitly asked to be location-scale.
            %
            % Concurrent variants are deferred by ruling E until the observed-scale
            % residual COVARIANCE of the Gaussian-copula scaled-chi-square model is
            % defined under direct sigma. Rejecting is the correct behavior, not a
            % placeholder for missing plumbing.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localGammaDiffTable(), 'family','gamma', ...
                'diffest',2, 'diffrescor',2, 'gammascale',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:gammascale');
        end

        function testGammascaleTwoRejectedForConcurrentSubjectDifference(testCase)
            % The same gap on the subject-level side (analysis 8 concurrent).
            % Here sserrvar=2 makes isGroupConcurrentGammaDiff FALSE, so
            % 'dispersion' defaults to 1 and its guard is inert for the opposite
            % reason. Both concurrent entry points therefore need the
            % gammascale-side guard, and both are pinned.
            %
            % diffwpcov=2 is exercised rather than diffrescor=2 because the option
            % parser normalizes diffwpcov=2 into diffrescor=2; asserting through
            % the alias proves the guard sees the normalized value and cannot be
            % bypassed by choosing the other flag.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localGammaDiffTable(), 'family','gamma', ...
                'diffest',2, 'sserrvar',2, 'diffwpcov',2, 'gammascale',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:gammascale');
        end

        function testGammascaleTwoAcceptedForNonconcurrentDifference(testCase)
            % The positive counterpart, so the two rejection tests above cannot
            % pass by rejecting the whole design family. A NON-concurrent gamma
            % two-event difference (analysis 7) must reach the builder under
            % gammascale=2 and emit the location-scale scale submodel.
            stanText = localEmit(testCase, localGammaDiffTable(), ...
                {'family','gamma','diffest',2,'gammascale',2});
            testCase.verifySubstring(stanText, 'vector[2] b_sigma;', ...
                'The location-scale difference model must declare per-event b_sigma.');
            testCase.verifySubstring(stanText, 'nu = 2 * square(mu ./ sigma);', ...
                'nu must be DERIVED from the mean and residual SD.');
            testCase.verifyEmpty(strfind(stanText, 'b_nu'), ...
                'A location-scale run must not declare or sample b_nu.'); %#ok<STREMP>
        end

        function testGammascaleTwoAcceptedForNonconcurrentTwoFacetDifference(testCase)
            % Increment 6: the TWO-FACET (test-retest) difference, analysis 10.
            % Same acceptance check as the one-facet case above, plus the
            % assertions that make it specifically two-facet - the five crossed
            % factor blocks must still be present, because the location-scale
            % variant is supposed to touch the scale submodel and nothing else.
            stanText = localEmit(testCase, localGammaDiffTrtTable(), ...
                {'family','gamma','diffest',2,'gammascale',2});

            testCase.verifySubstring(stanText, 'vector[2] b_sigma;', ...
                'The location-scale two-facet difference model must declare per-event b_sigma.');
            testCase.verifySubstring(stanText, 'nu = 2 * square(mu ./ sigma);', ...
                'nu must be DERIVED from the mean and residual SD.');
            testCase.verifyEmpty(strfind(stanText, 'b_nu'), ...
                'A location-scale run must not declare or sample b_nu.'); %#ok<STREMP>

            % The mean submodel must be intact: all five crossed factors plus the
            % six cross-event correlations. If a future edit narrows the
            % location-scale branch to the one-facet structure, this fails here
            % rather than in a live fit.
            for name = {'sd_trl','sd_occ','sd_tid','sd_oid','sd_to'}
                testCase.verifySubstring(stanText, name{1}, ...
                    sprintf('The two-facet mean submodel must retain %s.', name{1}));
            end
            for name = {'cor_p','cor_o','cor_t','cor_po','cor_pt','cor_ot'}
                testCase.verifySubstring(stanText, name{1}, ...
                    sprintf('The two-facet model must expose %s.', name{1}));
            end
        end

        function testGammascaleTwoRejectedForConcurrentTwoFacetDifference(testCase)
            % The concurrency guard keys on the two OPTIONS
            % (gammascale==2 && diffrescor==2), not on the analysis number, so it
            % covered analysis 10 the moment that design entered scope. This pins
            % that rather than leaving it to inspection - the increment-5 lesson
            % was that a guard's reachability argument can go stale silently while
            % every existing test stays green.
            %
            % Its concurrent builders (psyrat_build_stan_diff_trt_copula_chisq and
            % the fixnu variant) have no location-scale branch, so without the
            % guard this combination would emit a log-nu model for a run the user
            % explicitly asked to be location-scale. Deferred by ruling E.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localGammaDiffTrtTable(), 'family','gamma', ...
                'diffest',2, 'diffrescor',2, 'gammascale',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:gammascale');
        end

        function testGammascaleTwoAcceptedForDifferenceOfDifferences(testCase)
            % Increment 7: the four-event difference-of-differences, analysis 9 -
            % the last static design in scope. Same acceptance shape as the two
            % cases above, plus the assertions that make it specifically
            % four-cell, because the location-scale variant is supposed to touch
            % the scale submodel and nothing else.
            stanText = localEmit(testCase, localGammaDodTable(), ...
                {'family','gamma','diffest',3,'dodmap',{'E1','E2','E3','E4'}, ...
                'gammascale',2});

            testCase.verifySubstring(stanText, 'vector[4] b_sigma;', ...
                'The location-scale DoD model must declare per-CELL b_sigma (four of them).');
            testCase.verifySubstring(stanText, 'nu = 2 * square(mu ./ sigma);', ...
                'nu must be DERIVED from the mean and residual SD.');
            testCase.verifyEmpty(strfind(stanText, 'b_nu'), ...
                'A location-scale run must not declare or sample b_nu.'); %#ok<STREMP>

            % The mean submodel must be intact: the 8-D person block (four means
            % plus four log-sigmas) and the 4-D trial block. If a future edit
            % narrows the location-scale branch to the two-event structure, this
            % fails here rather than in a live fit.
            testCase.verifySubstring(stanText, 'vector<lower=0>[8] sd_id;', ...
                'The DoD person block must stay 8-D under location-scale.');
            testCase.verifySubstring(stanText, 'cholesky_factor_corr[8] L_id;', ...
                'The DoD person correlation must stay 8-D under location-scale.');
            testCase.verifySubstring(stanText, 'vector<lower=0>[4] sd_trl;', ...
                'The DoD trial block must stay 4-D (per-cell mean only).');
            testCase.verifySubstring(stanText, 'r_id[ID[n], c + 4]', ...
                'The scale effect must read person dims 5-8, keeping means in 1-4.');
            % The scale prior must be the log-sigma one, not the log-nu one: a
            % branch that emitted b_sigma but kept sig_nu would look right above
            % and be wrong in the only place the two priors differ.
            testCase.verifySubstring(stanText, 'sd_id[5:8] ~ student_t', ...
                'The person log-sigma SDs must carry their own prior line.');
            testCase.verifyEmpty(strfind(stanText, 'nu_center'), ...
                'A location-scale run must not use the log-nu prior centre.'); %#ok<STREMP>
        end

        function testGammascaleTwoAcceptedForTwoFacetDynrel(testCase)
            % Increment 8: the TWO-FACET dynamic design, analysis 14 - the first
            % location-scale design carrying dimension slopes on both submodels
            % TOGETHER WITH crossed facets. The assertions that make it
            % specifically two-facet are the five crossed facet SDs: the
            % location-scale variant is supposed to touch the scale submodel and
            % nothing else, so if a future edit narrows this branch to the
            % one-facet structure it fails here rather than in a live fit.
            stanText = localEmit(testCase, localDynrelTrtTable(), ...
                {'family','gamma','dynrel',2,'gammascale',2});

            testCase.verifySubstring(stanText, 'vector[KDIM] b_sigma;', ...
                'The location-scale two-facet dynrel model must declare b_sigma.');
            testCase.verifySubstring(stanText, 'nu = 2 * square(mu ./ sigma);', ...
                'nu must be DERIVED from the mean and residual SD.');
            testCase.verifyEmpty(strfind(stanText, 'b_nu'), ...
                'A location-scale run must not declare or sample b_nu.'); %#ok<STREMP>
            testCase.verifyEmpty(strfind(stanText, 'Intercept_nu'), ...
                'A location-scale run must not declare Intercept_nu.'); %#ok<STREMP>

            % The five crossed facets must remain on the MEAN submodel, with their
            % raw vectors and their transformed terms. This is the structural
            % counterpart of the estimand-invariance argument: the six
            % observed-score signal components are built from the log-mean
            % submodel alone, so they must be parameterization-invariant.
            for name = {'sig_occ','sig_trl','sig_trlxid','sig_occxid','sig_trlxocc'}
                testCase.verifySubstring(stanText, name{1}, ...
                    sprintf('The two-facet mean submodel must retain %s.', name{1}));
            end
            for name = {'occ_terms','trl_terms','trlxid_terms','occxid_terms', ...
                    'trlxocc_terms'}
                testCase.verifySubstring(stanText, name{1}, ...
                    sprintf('The two-facet mean submodel must retain %s.', name{1}));
            end
            % ...and NOT be duplicated onto the scale submodel. The residual is the
            % highest-order p x t x o term; giving log sigma a facet term would be
            % a different model, not a reparameterization of this one.
            testCase.verifySubstring(stanText, ...
                'vector[NOBS] sigma = Intercept_sigma + Xdim * b_sigma;', ...
                'The scale submodel must carry the intercept and slopes only.');
        end

        function testGammascaleTwoAcceptedForTwoFacetDynrelSubjectLevel(testCase)
            % Analysis 27 shares its case block AND its builder with analysis 14,
            % so the emitted model must be identical - the subject-level
            % distinction is a downstream extraction (ind_sd is a transformed
            % parameter of that same model) plus a conditional calculator. Pinning
            % the routing here is what proves sserrvar=2 does not fall out of the
            % location-scale branch.
            groupText = localEmit(testCase, localDynrelTrtTable(), ...
                {'family','gamma','dynrel',2,'gammascale',2});
            subjText = localEmit(testCase, localDynrelTrtTable(), ...
                {'family','gamma','dynrel',2,'sserrvar',2,'gammascale',2});
            testCase.verifyEqual(subjText, groupText, ...
                ['The analysis-27 location-scale model must be byte-identical to ' ...
                'the analysis-14 one; the builder is shared.']);
            testCase.verifySubstring(subjText, 'vector[NSUB] ind_sd;', ...
                'The subject-level read-out needs ind_sd as a transformed parameter.');
        end

        function testTwoFacetDynrelUnderLogNuStillEmitsLogNu(testCase)
            % The NEGATIVE half of the pair above, and the reason it exists: a
            % green location-scale test proves nothing about whether the log-nu
            % path still works. The analysis-14 gamma run must still emit the
            % log-nu submodel untouched when log-nu is asked for.
            %
            % gammascale is now passed EXPLICITLY. Analysis 14 is movable, so
            % since rulings 33-36 omitting it selects location-scale and this test
            % would be asserting log-nu against an LS model. The explicit 1 keeps
            % the test pinning what its name says.
            stanText = localEmit(testCase, localDynrelTrtTable(), ...
                {'family','gamma','dynrel',2,'gammascale',1});
            testCase.verifySubstring(stanText, 'vector[KDIM] b_nu;', ...
                'The default two-facet dynrel gamma run must still declare b_nu.');
            testCase.verifyEmpty(strfind(stanText, 'b_sigma'), ...
                'A log-nu run must not declare b_sigma.'); %#ok<STREMP>
            testCase.verifyEmpty(strfind(stanText, 'nu = 2 * square'), ...
                'A log-nu run must SAMPLE nu, not derive it.'); %#ok<STREMP>
        end

        function testGammascaleTwoAcceptedForDynamicDifference(testCase)
            % Increment 8b: the DYNAMIC DIFFERENCE design, analysis 12 - the first
            % location-scale design with cross-event covariance AND dimension
            % slopes on both submodels. The assertions that make it specifically
            % case 12 are the 4-D person block and the event-crossed Xdim design;
            % if a future edit narrowed this branch to the non-difference dynrel
            % structure it would fail here rather than in a live fit.
            stanText = localEmit(testCase, localGammaDiffDynrelTable(), ...
                {'family','gamma','diffest',2,'dynrel',2,'gammascale',2});

            % NAMING: unlike the non-difference dynrel builders, here b_sigma is
            % the per-EVENT INTERCEPT and b_sigma_dim the slope, following the
            % static difference builder.
            testCase.verifySubstring(stanText, 'vector[2] b_sigma;', ...
                'Case 12 location-scale must declare per-event b_sigma intercepts.');
            testCase.verifySubstring(stanText, 'vector[KDIM] b_sigma_dim;', ...
                'Case 12 location-scale must declare event-crossed b_sigma_dim.');
            testCase.verifySubstring(stanText, 'nu = 2 * square(mu ./ sigma);', ...
                'nu must be DERIVED from the mean and residual SD.');
            testCase.verifyEmpty(strfind(stanText, 'b_nu'), ...
                'A location-scale run must not declare or sample b_nu.'); %#ok<STREMP>

            % The difference structure must survive: the 4-D person block with its
            % correlation, the 2-D cross-event trial block, and the event-crossed
            % dimension design. These are what distinguish case 12 from case 11.
            testCase.verifySubstring(stanText, ...
                'vector<lower=0>[4] sd_id;  // person SDs: [mean_1, mean_2, logsigma_1, logsigma_2]', ...
                'The 4-D person block must be retained with log-sigma in slots 3-4.');
            testCase.verifySubstring(stanText, 'cholesky_factor_corr[4] L_id;', ...
                'The 4x4 person correlation must be retained.');
            testCase.verifySubstring(stanText, 'cholesky_factor_corr[2] L_trl;', ...
                'The cross-event trial correlation must be retained.');
            testCase.verifySubstring(stanText, ...
                'matrix[NOBS, KDIM] Xdim = to_matrix(Xdim_vec, NOBS, KDIM);', ...
                'The event-by-dimension design must be rebuilt from Xdim_vec.');
            testCase.verifySubstring(stanText, 'real cor_p = cor_id[1, 2];', ...
                'The cross-event person-mean correlation must be emitted.');

            % The scale submodel carries the person effect but NO trial term: the
            % residual stays one sigma_e,p per participant per event.
            testCase.verifySubstring(stanText, ...
                'vector[NOBS] sigma = X * b_sigma + Xdim * b_sigma_dim;', ...
                'The scale submodel must carry per-event intercepts and slopes.');
            testCase.verifySubstring(stanText, ...
                'sigma[n] += r_id[ID[n], 3] * meas1[n] + r_id[ID[n], 4] * meas2[n];', ...
                'The scale submodel must take its person effects from slots 3-4.');
            testCase.verifyEmpty(strfind(stanText, 'sigma[n] += r_trl'), ...
                'The scale submodel must carry no trial term.'); %#ok<STREMP>
        end

        function testDynamicDifferenceUnderLogNuStillEmitsLogNu(testCase)
            % The NEGATIVE half of the pair above: an explicit gammascale=1 on
            % case 12 must still produce the log-nu model.
            %
            % This test USED to omit gammascale, and in that form it was the
            % tripwire that would catch a silent default flip (rulings 31-32). The
            % flip has since been made deliberately (rulings 33-36) and case 12 is
            % movable, so the omission now selects location-scale. The tripwire
            % role moved to testExplicitGammascaleTwoMatchesDefault below; this
            % test reverts to pinning the log-nu builder, which is still reachable
            % and still needs coverage.
            stanText = localEmit(testCase, localGammaDiffDynrelTable(), ...
                {'family','gamma','diffest',2,'dynrel',2,'gammascale',1});
            testCase.verifySubstring(stanText, 'vector[2] b_nu;', ...
                'The default case-12 gamma run must still declare b_nu.');
            testCase.verifySubstring(stanText, 'vector[KDIM] b_nu_dim;', ...
                'The default case-12 gamma run must still declare b_nu_dim.');
            testCase.verifyEmpty(strfind(stanText, 'b_sigma'), ...
                'A log-nu run must not declare b_sigma.'); %#ok<STREMP>
            testCase.verifyEmpty(strfind(stanText, 'nu = 2 * square'), ...
                'A log-nu run must SAMPLE nu, not derive it.'); %#ok<STREMP>
        end

        function testDynamicDifferenceMeanSubmodelIsParameterizationInvariant(testCase)
            % The structural form of the estimand-invariance argument, and the
            % reason psyrat_gamma_crosscov is reused unchanged: the MEAN submodel
            % must be identical between the two parameterizations, so the
            % observed-scale signal components and BOTH cross-event covariances
            % cannot depend on the scale parameterization.
            %
            % BOTH parameterizations are named EXPLICITLY. nuText previously
            % relied on the default; once case 12's default became location-scale
            % (rulings 33-36) that made this test compare an LS model against
            % itself, so it would have kept passing while testing nothing. It did
            % not fail - which is why it is called out here rather than left to be
            % noticed.
            lsText = localEmit(testCase, localGammaDiffDynrelTable(), ...
                {'family','gamma','diffest',2,'dynrel',2,'gammascale',2});
            nuText = localEmit(testCase, localGammaDiffDynrelTable(), ...
                {'family','gamma','diffest',2,'dynrel',2,'gammascale',1});

            meanLines = {
                'vector[NOBS] mu = X * b + Xdim * b_dim;'
                '  mu[n] += r_id[ID[n], 1] * meas1[n] + r_id[ID[n], 2] * meas2[n]'
                '        + r_trl[TRL[n], 1] * meas1[n] + r_trl[TRL[n], 2] * meas2[n];'
                'mu = exp(mu);'
                'vector[2] b;  // per-event log-mean intercepts (at z = 0)'
                'vector[KDIM] b_dim;  // event-specific log-mean dimension slopes'
                'vector<lower=0>[2] sd_trl;  // per-event trial main-effect log-SDs'
                };
            for k = 1:numel(meanLines)
                testCase.verifySubstring(lsText, meanLines{k}, ...
                    sprintf('LS model must retain the mean-submodel line: %s', meanLines{k}));
                testCase.verifySubstring(nuText, meanLines{k}, ...
                    sprintf('log-nu model must retain the mean-submodel line: %s', meanLines{k}));
            end
        end

        function testGammascaleTwoAcceptedForSubjectLevelDynamicDifference(testCase)
            % Increment 8c: the SUBJECT-LEVEL dynamic difference, analysis 13
            % (dynrel=2 + diffest=2 + sserrvar=2). The emitted model must be the
            % case-12 location-scale model VERBATIM - the builder is a one-line
            % delegation, because the gamma case-12 person block is already 4-wide
            % with the two events' scale effects, unlike the Gaussian pair where
            % the subject-level builder has to widen a 2-D block.
            %
            % The byte-identity itself is pinned by TestStanModelSnapshots; what
            % this test pins is that the ROUTING reaches that builder at all.
            stanText = localEmit(testCase, localGammaDiffDynrelTable(), ...
                {'family','gamma','diffest',2,'dynrel',2,'sserrvar',2, ...
                'gammascale',2});

            testCase.verifySubstring(stanText, 'vector[2] b_sigma;', ...
                'Case 13 location-scale must declare per-event b_sigma intercepts.');
            testCase.verifySubstring(stanText, 'vector[KDIM] b_sigma_dim;', ...
                'Case 13 location-scale must declare event-crossed b_sigma_dim.');
            testCase.verifySubstring(stanText, 'nu = 2 * square(mu ./ sigma);', ...
                'nu must be DERIVED from the mean and residual SD.');
            testCase.verifyEmpty(strfind(stanText, 'b_nu'), ...
                'A location-scale run must not declare or sample b_nu.'); %#ok<STREMP>
            testCase.verifySubstring(stanText, ...
                'vector<lower=0>[4] sd_id;  // person SDs: [mean_1, mean_2, logsigma_1, logsigma_2]', ...
                'The 4-D person block must carry log-sigma in slots 3-4.');
            testCase.verifySubstring(stanText, ...
                'sigma[n] += r_id[ID[n], 3] * meas1[n] + r_id[ID[n], 4] * meas2[n];', ...
                'The scale submodel must take its person effects from slots 3-4.');

            % NEGATIVE: the gamma path must NOT pick up the Gaussian
            % subject-level builder's er_var_ss generated quantity. Under gamma
            % the per-subject residual is formed in MATLAB from r_id, so a model
            % emitting er_var_ss would mean the Gaussian builder was reached.
            testCase.verifyEmpty(strfind(stanText, 'er_var_ss'), ...
                ['The gamma case-13 model must not emit the Gaussian ' ...
                'er_var_ss generated quantity.']); %#ok<STREMP>
        end

        function testSubjectLevelDynamicDifferenceRejectedUnderLogNu(testCase)
            % The NEGATIVE half of the pair above, and the test that pins the
            % scope decision itself. Analysis 13 is the ONLY gamma analysis
            % supported in a single parameterization, so unlike every other
            % design in this file its log-nu counterpart must ERROR rather than
            % emit a log-nu model.
            %
            % This matters because the rejection is NOT fail-safe by default: the
            % family guard keys on the analysis number alone, so once 13 is in its
            % list a DEFAULT gammascale=1 run passes it, and the gammascale guard
            % only fires when gammascale==2. Without the dedicated guard the run
            % would reach the case-13 block and build the GAUSSIAN model for a
            % gamma request - the silent Gaussian-on-gamma outcome the family
            % guard exists to prevent.
            args = {'data', localGammaDiffDynrelTable(), 'chains',1, ...
                'warmup',1, 'sampling',1, 'seed',12345, 'showgui',0, ...
                'verbose',1, 'family','gamma', 'diffest',2, 'dynrel',2, ...
                'sserrvar',2};

            % explicit gammascale = 1 - still rejected, and this is now the ONLY
            % way to reach the rejection.
            testCase.verifyError(@() psyrat_computevarcomp(args{:}, ...
                'gammascale',1), 'varargin:gammascale');

            % The DEFAULT no longer errors, and that is deliberate (ruling 36).
            % Case 13 is in the set that resolves to gammascale=2, so a user who
            % omits the option now gets the one estimand the design HAS, instead
            % of being told to pass the value that is already the default. The
            % guard above is retained and fires only on an explicit 1, which is
            % what keeps it fail-closed if the default is ever moved back.
            %
            % Assert the default resolves to 2 rather than merely "does not
            % error", so this cannot pass for the wrong reason.
            stanText = localEmit(testCase, localGammaDiffDynrelTable(), ...
                {'family','gamma','diffest',2,'dynrel',2,'sserrvar',2});
            testCase.verifySubstring(stanText, 'b_sigma', ...
                'A default case-13 run must resolve to the location-scale model.');
            testCase.verifyEmpty(strfind(stanText, 'Intercept_nu'), ...
                'A default case-13 run must not emit a log-nu submodel.'); %#ok<STREMP>
        end

        function testGammascaleTwoRejectedForConcurrentDifferenceOfDifferences(testCase)
            % Analysis 9 is rejected for a concurrent request by TWO independent
            % guards, and this pins that the combination cannot get through. The
            % generic gammascale==2 && diffrescor==2 guard fires first (ruling E);
            % a family-level case-9 guard would reject it regardless of
            % gammascale, because case 9 has no copula twin in EITHER
            % parameterization. Neither needed widening for this increment - which
            % is exactly the kind of claim that goes stale silently, so it is
            % pinned rather than left to inspection.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localGammaDodTable(), 'family','gamma', ...
                'diffest',3, 'dodmap',{'E1','E2','E3','E4'}, ...
                'diffrescor',2, 'gammascale',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:gammascale');
        end

        function testConcurrentDifferenceOfDifferencesRejectedUnderLogNuToo(testCase)
            % The second of the two guards, isolated. At gammascale = 1 the
            % generic guard cannot fire, so this reaches the family-level case-9
            % guard and must still be rejected - with a DIFFERENT error id. Without
            % this, the pair above could pass while the case-9 guard had quietly
            % stopped working, since the generic guard would mask it.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localGammaDodTable(), 'family','gamma', ...
                'diffest',3, 'dodmap',{'E1','E2','E3','E4'}, ...
                'diffrescor',2, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:familyconcurrent');
        end

        function testInvalidGammascaleRejected(testCase)
            % An out-of-range value must error rather than fall through to a
            % default, so a typo cannot silently change the estimand.
            testCase.verifyError(@() psyrat_computevarcomp( ...
                'data', localDynrelTable(), 'family','gamma', 'dynrel',2, ...
                'gammascale',0, ...
                'chains',1,'warmup',1,'sampling',1,'seed',12345,'showgui',0,'verbose',1), ...
                'varargin:gammascale');
        end

        function testRunThreadsGammascaleToTheEngine(testCase)
            % The option is useless if it cannot be reached from the documented
            % entry point. psyrat_run must forward it through the warp to the
            % engine; reaching the engine's scope guard (this run is analysis 1,
            % not 11/26) is the proof that it arrived.
            tdir = localTempDir(testCase);
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.4; 0.9; 1.0], ...
                'VariableNames', {'id','meas'});
            testCase.verifyError(@() psyrat_run('data', tbl, 'family','gamma', ...
                'gammascale', 2, 'savepath', tdir, 'savename', 'x.psyrat'), ...
                'varargin:gammascale');
        end

        function testRunconfigRecordsGammascale(testCase)
            % gammascale = 1 and 2 are DIFFERENT ESTIMANDS, so a sidecar that
            % lost it would replay into whichever default the running version
            % happens to have and silently produce different variance
            % components. Recorded as-is, [] when unset.
            tdir = localTempDir(testCase);

            prefs = psyrat_defaults();
            prefs.proc.gammascale = 2;
            data = struct();
            data.proc.savepath = tdir;
            data.proc.savename = 'ls.psyrat';
            data.rel.analysis = 'ic';
            cfg = jsondecode(fileread(psyrat_write_runconfig(data, prefs)));
            testCase.verifyEqual(cfg.design.gammascale, 2);

            % An unset value still records []. What that MEANS changed with
            % rulings 33-36 and is worth stating, because the assertion below
            % survived the flip unchanged while its guarantee did not: re-deriving
            % the default now yields location-scale for a movable design, so a []
            % sidecar no longer replays as log-nu.
            %
            % This is not a live hazard for sidecars the toolbox writes. The warp
            % back-copies the RESOLVED REL.gammascale into prefs before calling
            % psyrat_write_runconfig, so any sidecar written through the normal
            % path pins a concrete 1 or 2. [] arises only when the writer is
            % called directly, as here, or from a version predating that
            % back-copy - which is the migration case the release note covers.
            prefs2 = psyrat_defaults();
            data.proc.savename = 'unset.psyrat';
            cfg2 = jsondecode(fileread(psyrat_write_runconfig(data, prefs2)));
            testCase.verifyEmpty(cfg2.design.gammascale, ...
                'An unset gammascale must record [] rather than a guessed value.');
        end

        function testRunConfigReplaysGammascale(testCase)
            % A location-scale run's sidecar must replay AS location-scale.
            % gammascale flows config -> proc -> warp -> engine; hitting the
            % engine's scope guard proves the whole chain, since a lost value
            % would default to 1 and run without error.
            tdir = localTempDir(testCase);
            prefs = psyrat_defaults();
            prefs.proc.family = 'gamma';
            prefs.proc.gammascale = 2;
            d = struct();
            d.proc.savepath = tdir;
            d.proc.savename = 'ls.psyrat';
            d.proc.idheader = 'id';
            d.proc.measheader = 'meas';
            d.rel.analysis = 'ic';
            sidecar = psyrat_write_runconfig(d, prefs);

            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.4; 0.9; 1.0], ...
                'VariableNames', {'id','meas'});
            testCase.verifyError(@() psyrat_run('config', sidecar, 'data', tbl, ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'varargin:gammascale');
        end

        function testExplicitGammascaleTwoMatchesDefault(testCase)
            % THE TRIPWIRE. It pins the default's VALUE, not merely its shape: an
            % explicit 2 and the default must emit byte-identical Stan on a
            % movable design, so the default really is 2 rather than "something
            % location-scale shaped".
            %
            % This is the inversion of testExplicitGammascaleOneMatchesDefault,
            % which pinned the default at exactly 1 and did its job - it failed
            % the moment rulings 33-36 moved the default, which is what ruling 31
            % asked a tripwire to do. It is inverted rather than deleted, because
            % deleting it would leave nothing asserting the default's value and
            % the next move could then happen silently.
            defaultText  = localEmit(testCase, localDynrelTable(), ...
                {'family','gamma','dynrel',2});
            explicitText = localEmit(testCase, localDynrelTable(), ...
                {'family','gamma','dynrel',2,'gammascale',2});
            testCase.verifyEqual(defaultText, explicitText, ...
                'The gammascale default must be exactly 2 on a movable design.');

            % ... and must NOT match an explicit 1. Without this half the test
            % would still pass if the two parameterizations ever collapsed into
            % the same emitted model.
            logNuText = localEmit(testCase, localDynrelTable(), ...
                {'family','gamma','dynrel',2,'gammascale',1});
            testCase.verifyNotEqual(defaultText, logNuText, ...
                'The default must differ from the log-nu parameterization.');
        end

        function testStaleConcurrentPreferenceDoesNotVetoTheDefault(testCase)
            % REGRESSION for the SILENT half of the diffrescor keying defect.
            %
            % diffrescor/diffwpcov is a PERSISTENT processing preference, and the
            % option parser normalizes diffwpcov=2 into diffrescor=2 without ever
            % consulting diffest. So a value left set during an earlier CONCURRENT
            % DIFFERENCE session travels into later, unrelated runs.
            %
            % The resolution block used to exclude concurrency with a bare
            % `diffrescor ~= 2`. That also vetoed the location-scale flip on the
            % NON-difference designs (6/11/14/25/26/27), which never read diffrescor
            % at all - analysis 11's builder call passes only (priors, threaded,
            % gammascale). The run then fitted LOG-NU while the documented default
            % is location-scale: no error, no warning, a different estimand on disk.
            %
            % Nothing behavioral covered this before: every diffrescor case in this
            % class passed diffest 2 or 3.
            stanText = localEmit(testCase, localDynrelTable(), ...
                {'family','gamma','dynrel',2,'diffrescor',2});

            testCase.verifySubstring(stanText, 'real Intercept_sigma;', ...
                ['A stale concurrent preference must not veto the location-scale ' ...
                'default on a NON-difference design, which never reads diffrescor.']);
            testCase.verifyEmpty(strfind(stanText, 'real Intercept_nu;'), ...
                ['The run fell back to log-nu: diffrescor leaked into a design ' ...
                'that does not use it, silently changing the estimand.']); %#ok<STREMP>
        end

        function testExplicitGammascaleTwoSurvivesStaleConcurrentPreference(testCase)
            % REGRESSION for the LOUD half of the same defect. Same stale
            % preference, but here the user has EXPLICITLY chosen location-scale.
            %
            % The rejection guard keyed on diffrescor alone, so it aborted with
            % 'varargin:gammascale' and advised omitting diffrescor "for the
            % non-concurrent difference score" - on a run with no difference score.
            % The combination ran on main, so the guard introduced a regression.
            %
            % localEmit asserts the run reaches the emit hook, so if the guard fires
            % again this fails on the error identifier (varargin:gammascale instead
            % of PsyRAT:emitModelOnly) rather than on the substring below.
            stanText = localEmit(testCase, localDynrelTable(), ...
                {'family','gamma','dynrel',2,'diffrescor',2,'gammascale',2});

            testCase.verifySubstring(stanText, 'real Intercept_sigma;', ...
                ['An explicit gammascale=2 must still emit the location-scale ' ...
                'model when an unrelated diffrescor=2 is left set.']);
        end

        % NOTE: the protection these two must NOT weaken is pinned elsewhere in
        % this class and was left untouched - concurrent DIFFERENCE runs still
        % reject an explicit gammascale=2 (testGammascaleTwoRejectedForConcurrent-
        % GroupDifference and ...SubjectDifference) and still default to log-nu
        % (testDefaultStaysLogNuForConcurrentDesigns). All three pass diffest=2,
        % which is exactly the term the fix added.

        function testResolutionBlockKeepsAllThreeGuardConditions(testCase)
            % SOURCE-SCANNING GUARD, and the offline lane's only protection for
            % the most dangerous edit in this feature.
            %
            % WHY A BEHAVIORAL TEST CANNOT DO THIS. Measured, not assumed: with
            % the `diffrescor ~= 2` conjunct deleted from the resolution block,
            % ALL 30 tests in this class and ALL 73 Stan snapshots still pass.
            % The reason is that analyses 7/8/10 reach copula builders whose
            % signature is (priors, threaded) - they have no gammascale parameter
            % and emit log-nu unconditionally. So the generated Stan is IDENTICAL
            % whether or not a concurrent run was wrongly defaulted to 2. What
            % actually breaks is invisible to any Stan-text assertion:
            % REL.gammascale is stamped 2 while the fitted model is log-nu, and
            % the downstream extraction then reads log-nu draws as location-scale
            % ones. Wrong numbers, no error, green suite.
            %
            % testDefaultStaysLogNuForConcurrentDesigns does NOT close this - it
            % asserts on emitted Stan and passes either way. An earlier version
            % of this file claimed it was the tripwire; that was wrong, and this
            % test exists because the claim was checked rather than trusted.
            %
            % Reading the source is the same tactic TestRelOutFieldInit uses for
            % the same reason: the failure only manifests in a real CmdStan run,
            % so the offline lane has to inspect the code instead of running it.
            src = fileread(localEngineSourcePath());

            resolveLine = localMatchOneLine(testCase, src, ...
                '^\s*if\s+~gammascale_user_set\s+.*$', ...
                'the context-aware gammascale resolution');

            testCase.verifySubstring(resolveLine, '~gammascale_user_set', ...
                ['The resolution must fire only when the user did NOT set ' ...
                'gammascale, or an explicit 1 would be overwritten.']);
            testCase.verifySubstring(resolveLine, 'strcmp(family,''gamma'')', ...
                ['The resolution must require the gamma family. These analysis ' ...
                'numbers also exist under Gaussian, where gammascale=2 is ' ...
                'rejected outright, so a Gaussian run would default into a ' ...
                'hard error.']);
            testCase.verifySubstring(resolveLine, 'diffrescor ~= 2', ...
                ['The resolution MUST exclude concurrent designs. Without this ' ...
                'a concurrent gamma run is stamped gammascale=2 while its ' ...
                'copula builder emits log-nu - a silent estimand mismatch that ' ...
                'no behavioral test in this repo detects.']);
            % The exclusion must be QUALIFIED by diffest, not keyed on diffrescor
            % alone. diffrescor/diffwpcov is a persistent preference normalized
            % without reference to diffest, so a stale 2 carried in from an earlier
            % concurrent session would otherwise veto the flip on the NON-difference
            % designs (6/11/14/25/26/27), which never read diffrescor, and fit
            % log-nu on a run whose documented default is location-scale.
            testCase.verifySubstring(resolveLine, 'diffest == 1', ...
                ['The concurrent exclusion must be qualified by diffest. Keyed ' ...
                'on diffrescor alone it also vetoes the non-difference designs, ' ...
                'where diffrescor is inert - a silent fallback to log-nu.']);

            % The MODULAR CUT-COPULA escape (analyses 19/20) is the one hole in
            % that exclusion, and its scope is the whole of its safety. Those two
            % designs are concurrent AND location-scale-only, so they must be let
            % through; 7/8/10 must not be, because their copula builders take no
            % gammascale argument and emit log-nu regardless - the exact silent
            % mismatch this exclusion exists to prevent. Pin that the escape is
            % keyed on an analysis-number test and on nothing broader.
            testCase.verifySubstring(resolveLine, 'isModularCutAnalysis', ...
                ['The concurrent exclusion may only be escaped through the ' ...
                'named modular-cut analysis test.']);
            escapeLine = localMatchOneLine(testCase, src, ...
                '^\s*isModularCutAnalysis\s*=.*$', ...
                'the modular-cut escape definition');
            testCase.verifySubstring(escapeLine, 'any(analysis == [19 20])', ...
                ['The escape must name analyses 19 and 20 explicitly. Widening ' ...
                'it - to diffest, to sserrvar, to a family test - would let the ' ...
                'full-joint copula designs default to location-scale while ' ...
                'their builders still emit log-nu.']);
        end

        function testResolutionListMatchesTheSupportedDesignGuard(testCase)
            % The invariant that keeps the two lists honest: the default resolves
            % to location-scale EXACTLY where location-scale is supported. If
            % they drift, either a defaulted run trips the guard (a user is told
            % to remove an option they never passed) or a design silently keeps
            % log-nu after being declared movable.
            src = fileread(localEngineSourcePath());

            resolveList = localAnalysisList(testCase, src, ...
                'if\s+~gammascale_user_set[^\n]*\n[^\n]*any\(analysis\s*==\s*\[([^\]]*)\]', ...
                'the resolution block');
            guardList = localAnalysisList(testCase, src, ...
                'if\s+gammascale\s*==\s*2\s*&&\s*~any\(analysis\s*==\s*\[([^\]]*)\]', ...
                'the supported-design guard');

            testCase.verifyEqual(resolveList, guardList, ...
                ['The resolution list and the supported-design guard list must ' ...
                'stay identical. Update both, or neither.']);
        end

        function testResolutionBlockRestampsRelGammascale(testCase)
            % SOURCE-SCANNING GUARD, third of three, and the only thing standing
            % between a deleted line and a wrong estimand recorded on disk.
            %
            % The resolution block does two things. `gammascale = 2` steers the
            % BUILDER, and every behavioral test in this class sees it. The
            % re-stamp `REL.gammascale = gammascale` records which model was
            % actually fitted, and NOTHING in the offline lane sees it: the emit
            % hook aborts ~17,900 lines further down and REL is never returned,
            % and there is no CmdStan-free path to a completed gamma fit.
            %
            % So deleting the re-stamp leaves the emitted Stan byte-identical and
            % every OTHER test in the UNIT lane green - measured, not assumed: with
            % the line deleted that lane ran 904 tests, 903 passed, and the single
            % failure was this one. (The live recovery lane is no help either: its
            % tests pass gammascale explicitly, which sets gammascale_user_set and
            % skips this block entirely.)
            %
            % WHAT ACTUALLY BREAKS, in two different ways. A defaulted movable run
            % fits location-scale while REL.gammascale says 1.
            %   (i) On the six DYNAMIC movable designs - 11/12/13/14/26/27, the
            %   ones whose results are read through psyrat_dynrel_summary, whose
            %   only callers are the dynrel viewer (psyrat_startview_dynrel:267)
            %   and psyrat_report:491; there is no numeric list in the dispatch to
            %   read this off, so it is derived from the design numbering - the
            %   isequal(REL.gammascale,2) gate at psyrat_dynrel_summary:193 drops
            %   the run into the log-nu branch, which reads REL.out.pop_lognu and
            %   b_nu - fields the location-scale layout deliberately never creates
            %   (psyrat_computevarcomp.m:4736). That is a HARD CRASH on the first
            %   stratum, NOT silent wrong numbers: MATLAB:nonExistentField at
            %   psyrat_dynrel_summary:678. That file's own comment says so at
            %   :452-455, and an earlier draft of THIS comment asserted silent
            %   wrong numbers and so contradicted the file it was describing. The
            %   crash lands at view/report time, after the .psyrat has already been
            %   saved carrying the wrong stamp.
            %   (ii) SILENTLY, on all twelve movable designs - including the six
            %   static ones that never reach that gate, since psyrat_relsummary
            %   reads the field nowhere - psyrat_computevarcompwarp.m:587-588
            %   back-copies the wrong 1 into the run-config sidecar. Replaying that
            %   sidecar passes an explicit 1, which suppresses this very block and
            %   fits LOG-NU: an internally consistent run that quietly reproduces
            %   the WRONG ESTIMAND rather than the model originally fitted.
            % A loud crash on half the designs, and a permanently wrong
            % reproducibility record on all of them.
            %
            % Neither sibling scan covers it: testResolutionBlockKeepsAllThree-
            % GuardConditions reads only the `if` line, and testResolutionList-
            % MatchesTheSupportedDesignGuard reads only the two analysis lists.
            % The consumer-side tests in TestGammaLsDynrelSummaryDispatch set
            % rel.gammascale by hand, so they exercise the gate and never the
            % producer that fills it.
            %
            % WHAT THIS DOES NOT PROVE. It pins that the assignment is PRESENT in
            % the block, not that the returned REL carries 2. That still needs a
            % live CmdStan gamma fit and is deliberately left as an opportunistic
            % check rather than a reason to run one.
            src = fileread(localEngineSourcePath());

            % Capture the block BODY, not just its header - the header is what the
            % sibling scans already cover. Lazy to the first line-initial `end`,
            % which is stricter than a word boundary: it requires `end` to stand
            % alone on its line, which is exactly the block terminator.
            %
            % Do NOT reach for \b here. In MATLAB regexp \b is a literal BACKSPACE
            % (char 8), NOT a word boundary - measured, not assumed: 'cat\b' finds
            % nothing in 'cat dog', while 'end\b' does match ['end' char(8)].
            % MATLAB spells the word boundary \< and \>. The first draft of this
            % scanner ended its pattern with end\b, matched nothing at all, and was
            % caught only by the assert below - which is why that assert checks for
            % exactly one match instead of just indexing block{1}.
            block = regexp(src, ...
                'if\s+~gammascale_user_set[\s\S]*?\n[ \t]*end[ \t]*(\r?\n|$)', ...
                'match');
            testCase.assertEqual(numel(block), 1, ...
                sprintf(['Expected exactly one context-aware gammascale ' ...
                'resolution block, found %d. It was renamed or duplicated - ' ...
                'fix this scanner rather than deleting it.'], numel(block)));

            % Strip FULL-LINE comments before asserting. Two of the block's four
            % comment lines mention REL.gammascale (neither quotes the assignment
            % verbatim today) and a future one could. Know the limit: '^\s*%'
            % catches full-line comments ONLY - a trailing comment and the body of
            % a %{ %} block both survive it - so this narrows the assertion to
            % lines that do not START as comments, not to code as such. It cannot
            % over-strip, because in MATLAB a line whose first non-blank character
            % is % is always a comment.
            lines = regexp(block{1}, '\r\n|\n|\r', 'split');
            lines = lines(cellfun(@isempty, regexp(lines, '^\s*%', 'once')));
            body = strjoin(lines, newline);

            testCase.verifySubstring(body, 'REL.gammascale = gammascale;', ...
                ['The resolution block MUST re-stamp REL.gammascale. Without ' ...
                'it a defaulted run fits the location-scale model while the ' ...
                'saved result and its run-config sidecar both record the ' ...
                'provisional 1 - which crashes psyrat_dynrel_summary on the ' ...
                'dynamic designs and silently pins the wrong estimand in the ' ...
                'sidecar on all of them.']);
            testCase.verifySubstring(body, 'gammascale = 2;', ...
                ['The resolution block must still set gammascale = 2, or the ' ...
                'test above would pass on a block that stopped resolving and ' ...
                're-stamped the provisional value onto itself.']);
        end

        function testNonMovableDefaultStillMatchesExplicitOne(testCase)
            % The tripwire's counterpart on the other side of the condition: for a
            % design with no location-scale variant the default must still be
            % exactly 1, byte for byte. Together the two pin BOTH branches of the
            % context-aware default, so neither can drift unnoticed.
            defaultText  = localEmit(testCase, localPlainGammaTable(), ...
                {'family','gamma'});
            explicitText = localEmit(testCase, localPlainGammaTable(), ...
                {'family','gamma','gammascale',1});
            testCase.verifyEqual(defaultText, explicitText, ...
                'A non-movable design''s default must be exactly 1.');
        end

    end
end

function stanText = localEmit(testCase, dataTbl, extraArgs)
%Run psyrat_computevarcomp with the emit hook so the generated Stan is written
%and the run aborts before CmdStan; return the model text.
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

function tbl = localDynrelTable()
%Deterministic dimension-bearing long table for the one-facet dynrel routing:
%20 participants x 8 trials, one constant dimension value per participant, and
%strictly positive scores (the gamma family's support).
ids = {}; meas = []; dim1 = [];
for s = 1:20
    sval = -1 + (s - 1) * (2 / 19); % subject-level covariate spanning [-1, 1]
    for t = 1:8
        ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
        meas(end+1, 1) = 5 + 0.1 * t + 0.5 * sval; %#ok<AGROW>
        dim1(end+1, 1) = sval; %#ok<AGROW>
    end
end
tbl = table(ids, meas, dim1, 'VariableNames', {'id','meas','dim1'});
end

function p = localEngineSourcePath()
%Absolute path to the estimator core, derived from THIS file's location so the
%scan does not depend on pwd or on which copy of the toolbox is first on the
%MATLAB path (a .claude/worktrees copy shadowing the real one has bitten this
%repo before).
here = fileparts(mfilename('fullpath'));          % <root>/tests
p = fullfile(fileparts(here), 'subroutines', 'estimation', 'psyrat_computevarcomp.m');
end

function line = localMatchOneLine(testCase, src, pattern, what)
%Return the single source line matching pattern; fail loudly on 0 or 2+ matches
%rather than silently scanning the wrong one.
lines = regexp(src, '\r\n|\n|\r', 'split');
hits = ~cellfun(@isempty, regexp(lines, pattern, 'once'));
testCase.assertEqual(nnz(hits), 1, ...
    sprintf(['Expected exactly one source line for %s, found %d. The block ' ...
    'was renamed or duplicated - fix this scanner rather than deleting it.'], ...
    what, nnz(hits)));
line = lines{find(hits, 1)};
end

function nums = localAnalysisList(testCase, src, pattern, what)
%Extract an `any(analysis == [...])` list as a numeric row vector.
tok = regexp(src, pattern, 'tokens', 'once');
testCase.assertNotEmpty(tok, ...
    sprintf('Could not locate the analysis list for %s.', what));
nums = sort(str2double(strsplit(strtrim(tok{1}))));
testCase.assertFalse(any(isnan(nums)), ...
    sprintf('Non-numeric entry in the analysis list for %s.', what));
end

function tbl = localPlainGammaTable()
%Deterministic one-facet long table with NO event, time or dimension column, so
%the run routes to analysis 1 - a design with no location-scale variant. Used to
%pin that the context-aware default leaves the non-movable designs on log-nu.
%20 participants x 8 trials, strictly positive scores (the gamma family's support).
ids = {}; meas = [];
for s = 1:20
    for t = 1:8
        ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
        meas(end+1, 1) = 5 + 0.1 * t + 0.05 * s; %#ok<AGROW>
    end
end
tbl = table(ids, meas, 'VariableNames', {'id','meas'});
end

function tbl = localGammaDiffTable()
%Deterministic two-event long table for the static difference routing (analyses
%7/8): 20 participants x 8 trials x 2 events, strictly positive scores (the gamma
%family's support). The two events carry DIFFERENT amplitudes on purpose - a
%difference design with equal event means would not exercise the event contrast
%this parameterization exists to unblend - and no dim column, so the run routes
%static rather than dynamic.
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

function tbl = localGammaDiffDynrelTable()
%Deterministic two-event AND dimension-bearing long table for the DYNAMIC
%DIFFERENCE routing (analysis 12): 20 participants x 8 trials x 2 events, one
%constant dimension value per participant, strictly positive scores. It needs
%BOTH an event column (so diffest=2 routes to a difference design) and a dim
%column (so dynrel=2 routes dynamic) - neither localGammaDiffTable nor
%localDynrelTable has both, which is why this helper exists.
%
%The two events carry DIFFERENT amplitudes on purpose, for the same reason
%localGammaDiffTable does: a difference design with equal event means would not
%exercise the event contrast this parameterization exists to unblend.
ids = {}; meas = []; evt = {}; dim1 = [];
for s = 1:20
    sval = -1 + (s - 1) * (2 / 19); % subject-level covariate spanning [-1, 1]
    for e = 1:2
        for t = 1:8
            ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
            evt{end+1, 1} = sprintf('e%d', e); %#ok<AGROW>
            meas(end+1, 1) = 4 + e + 0.1 * t + 0.5 * sval; %#ok<AGROW>
            dim1(end+1, 1) = sval; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, evt, dim1, ...
    'VariableNames', {'id','meas','event','dim1'});
end

function tbl = localDynrelTrtTable()
%Deterministic dimension-bearing, two-OCCASION long table for the TWO-FACET
%dynrel routing (analyses 14/27): 20 participants x 2 occasions x 8 trials, one
%constant dimension value per participant, and strictly positive scores (the
%gamma family's support). The time column is what routes this to analysis 14
%rather than the one-facet analysis 11 (dynrel==2 && ntime>0); adding
%sserrvar==2 moves it to analysis 27. There is no event column, so it routes
%non-difference.
ids = {}; meas = []; dim1 = []; tme = {};
for s = 1:20
    sval = -1 + (s - 1) * (2 / 19); % subject-level covariate spanning [-1, 1]
    for o = 1:2
        for t = 1:8
            ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
            tme{end+1, 1} = sprintf('t%d', o); %#ok<AGROW>
            meas(end+1, 1) = 5 + 0.2 * o + 0.1 * t + 0.5 * sval; %#ok<AGROW>
            dim1(end+1, 1) = sval; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, dim1, tme, 'VariableNames', {'id','meas','dim1','time'});
end

function tbl = localGammaDiffTrtTable()
%Deterministic two-event, two-OCCASION long table for the static two-facet
%difference routing (analysis 10): 20 participants x 2 events x 2 occasions x 8
%trials, strictly positive scores (the gamma family's support). The time column
%is what routes this to analysis 10 rather than analysis 7
%(diffest==2 && sserrvar==1 && ntime>0). As in the one-facet table the two
%events carry DIFFERENT amplitudes on purpose, and there is no dim column, so
%the run routes static rather than dynamic.
ids = {}; meas = []; evt = {}; tme = {};
for s = 1:20
    for e = 1:2
        for o = 1:2
            for t = 1:8
                ids{end+1, 1} = sprintf('s%02d', s); %#ok<AGROW>
                evt{end+1, 1} = sprintf('e%d', e); %#ok<AGROW>
                tme{end+1, 1} = sprintf('t%d', o); %#ok<AGROW>
                meas(end+1, 1) = 4 + e + 0.2 * o + 0.1 * t + 0.05 * s; %#ok<AGROW>
            end
        end
    end
end
tbl = table(ids, meas, evt, tme, 'VariableNames', {'id','meas','event','time'});
end

function tbl = localGammaDodTable()
%Deterministic FOUR-event long table for the difference-of-differences routing
%(analysis 9): 20 participants x 4 events x 8 trials, strictly positive scores
%(the gamma family's support). The four event labels must match the dodmap the
%caller passes; diffest==3 plus that map is what routes to analysis 9.
%
%The four cells carry DIFFERENT amplitudes on purpose. A difference-of-
%differences with equal cell means would not exercise the cell contrasts, which
%are the quantity this parameterization exists to separate from amplitude.
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

function d = localTempDir(testCase)
%Self-contained scratch directory for the CLI round-trip tests, removed when the
%test finishes.
d = tempname;
mkdir(d);
testCase.addTeardown(@() localRemoveDir(d));
end

function localRemoveDir(d)
if exist(d, 'dir') == 7
    try, rmdir(d, 's'); catch, end
end
end

function localDeleteIfExists(f)
if exist(f,'file') == 2
    try, delete(f); catch, end
end
end
