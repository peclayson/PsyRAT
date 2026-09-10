classdef TestNativeEngine < PsyRATTestBase
    %Tests for the native (non-CmdStan) estimation engines and the engine
    %dispatch backbone: capability detection (psyrat_estimation_engines_available),
    %engine resolution/cascade (psyrat_resolve_engine), the duck-typed fit shim
    %(PsyRATNativeFit), and the native fitlme one-facet engine
    %(psyrat_native_fitlme via psyrat_run_native).
    %
    %None of these launch CmdStan, so they run in the standard unit suite. They
    %do require the Statistics and Machine Learning Toolbox (fitlme); tests that
    %need it are skipped with assumeTrue when it is unavailable so the suite
    %still passes on a minimal install.

    methods (Test)

        % --- capability detection ------------------------------------------
        function testEnginesAvailableShape(testCase)
            a = psyrat_estimation_engines_available();
            for f = {'cmdstan','hmc','fitlme','ad','native'}
                testCase.verifyTrue(isfield(a,f{1}), ...
                    sprintf('avail missing field %s',f{1}));
                testCase.verifyTrue(islogical(a.(f{1})) && isscalar(a.(f{1})), ...
                    sprintf('avail.%s should be a logical scalar',f{1}));
            end
            testCase.verifyTrue(isfield(a,'details') && isstruct(a.details));
            % cmdstan_status is the structured reason behind the cmdstan flag
            testCase.verifyTrue(isfield(a,'cmdstan_status') && ischar(a.cmdstan_status));
            % native is exactly the OR of the two native engines
            testCase.verifyEqual(a.native, a.fitlme || a.hmc);
        end

        % --- engine resolution / cascade -----------------------------------
        function testResolveForcedCmdstanNoDetection(testCase)
            % engine 1 must return 'cmdstan' without consulting availability,
            % so passing an empty avail (which would error if dereferenced)
            % still works.
            spec = localSpec(1,'sing');
            testCase.verifyEqual(psyrat_resolve_engine(1,spec,[]),'cmdstan');
        end

        function testResolveAutoPrefersCmdstan(testCase)
            spec = localSpec(0,'sing');
            av = localMockAvail(true,true,true);
            testCase.verifyEqual(psyrat_resolve_engine(0,spec,av),'cmdstan');
        end

        function testResolveAutoFallsBackToHmcWhenSupported(testCase)
            % The cascade prefers HMC over fitlme when CmdStan is off and BOTH are
            % available and the design supports HMC. 'trt' (case 5) is the one
            % design with both native engines, so it exercises this genuinely.
            spec = localSpec(0,'trt');
            av = localMockAvail(false,true,true); % no cmdstan, hmc+fitlme present
            testCase.verifyEqual(psyrat_resolve_engine(0,spec,av),'hmc');
        end

        function testResolveAutoFitlmeOnlyDesignPicksFitlmeDespiteHmc(testCase)
            % A fitlme-only design ('sing') picks fitlme even when HMC is also
            % installed, because the registry does not list 'hmc' for it.
            spec = localSpec(0,'sing');
            av = localMockAvail(false,true,true); % no cmdstan, hmc+fitlme present
            testCase.verifyEqual(psyrat_resolve_engine(0,spec,av),'fitlme');
        end

        function testResolveAutoFallsBackToFitlme(testCase)
            spec = localSpec(0,'sing');
            av = localMockAvail(false,false,true); % only fitlme available
            testCase.verifyEqual(psyrat_resolve_engine(0,spec,av),'fitlme');
        end

        function testResolveForcedFitlmeUnsupportedDesignErrors(testCase)
            spec = localSpec(2,'sserr'); % per-subject location-scale: not fitlme-doable
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        function testResolveForcedFitlmeUnavailableErrors(testCase)
            spec = localSpec(2,'sing');
            av = localMockAvail(true,true,false); % fitlme not installed
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnavailable');
        end

        function testResolveForcedHmcUnavailableErrors(testCase)
            spec = localSpec(3,'sing');
            av = localMockAvail(true,false,true); % hmc not installed
            testCase.verifyError(@() psyrat_resolve_engine(3,spec,av), ...
                'PsyRAT:engineUnavailable');
        end

        function testResolveNoEngineAvailableErrors(testCase)
            spec = localSpec(0,'sing');
            av = localMockAvail(false,false,false); % nothing available
            testCase.verifyError(@() psyrat_resolve_engine(0,spec,av), ...
                'PsyRAT:noEngineAvailable');
        end

        function testResolveBadEngineCodeErrors(testCase)
            spec = localSpec(0,'sing');
            testCase.verifyError(@() psyrat_resolve_engine(9,spec, ...
                localMockAvail(true,true,true)),'PsyRAT:badEngineCode');
        end

        % --- native HMC: dynrel (case 11) registry + routing ---------------
        function testNativeSupportedDynrelIsHmc(testCase)
            % The one-facet location-scale dynamic-reliability design is the
            % first HMC-only native design; fitlme cannot fit its per-subject
            % heteroscedastic residual.
            testCase.verifyEqual(psyrat_native_supported('dynrel'),{'hmc'});
        end

        function testResolveDynrelAutoPrefersHmc(testCase)
            % With CmdStan off and both native engines present, a dynrel run
            % must cascade to HMC (the design supports HMC, not fitlme).
            spec = localSpec(0,'dynrel');
            av = localMockAvail(false,true,true);
            testCase.verifyEqual(psyrat_resolve_engine(0,spec,av),'hmc');
        end

        function testResolveDynrelForcedHmc(testCase)
            % Forced HMC (engine 3) resolves to HMC when hmcSampler is present.
            spec = localSpec(3,'dynrel');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDynrelForcedFitlmeUnsupportedErrors(testCase)
            % fitlme cannot fit the location-scale dynrel model, so a forced
            % fitlme request for it must error (not silently misroute).
            spec = localSpec(2,'dynrel');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: sserr (case 6) registry + routing -----------------
        function testNativeSupportedSserrIsHmc(testCase)
            % Subject-level error variance is the base location-scale design
            % (KDIM = 0); a per-subject heteroscedastic residual is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('sserr'),{'hmc'});
        end

        function testResolveSserrForcedHmc(testCase)
            spec = localSpec(3,'sserr');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveSserrForcedFitlmeUnsupportedErrors(testCase)
            % fitlme cannot fit the per-subject heteroscedastic residual.
            spec = localSpec(2,'sserr');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff (case 7) registry + routing ------------------
        function testNativeSupportedDiffIsHmc(testCase)
            % Two-event difference (no rescor) is a multivariate location-scale
            % design fitlme cannot fit; HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff'),{'hmc'});
        end

        function testResolveDiffForcedHmc(testCase)
            spec = localSpec(3,'diff');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_rescor (case 7, diffrescor=2) registry/routing -
        function testNativeSupportedDiffRescorIsHmc(testCase)
            % Two-event difference WITH residual covariance is a bivariate
            % location-scale design fitlme cannot fit; HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_rescor'),{'hmc'});
        end

        function testResolveDiffRescorForcedHmc(testCase)
            spec = localSpec(3,'diff_rescor');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffRescorForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_rescor');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_dynrel (case 12) registry + routing ----------
        function testNativeSupportedDiffDynrelIsHmc(testCase)
            % Group-level dynamic difference (diff + event-by-dimension slopes on
            % the location-scale residual) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_dynrel'),{'hmc'});
        end

        function testResolveDiffDynrelForcedHmc(testCase)
            spec = localSpec(3,'diff_dynrel');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffDynrelForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_dynrel');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_dynrel_sserr (case 13) registry + routing ----
        function testNativeSupportedDiffDynrelSserrIsHmc(testCase)
            % Subject-level dynamic difference (4-D joint person block, per-subject
            % residual) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_dynrel_sserr'),{'hmc'});
        end

        function testResolveDiffDynrelSserrForcedHmc(testCase)
            spec = localSpec(3,'diff_dynrel_sserr');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffDynrelSserrForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_dynrel_sserr');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: dynrel_trt (case 14) registry + routing -----------
        function testNativeSupportedDynrelTrtIsHmc(testCase)
            % Two-facet (trial+occasion) dynamic reliability (locscale + five
            % scalar crossed mean facets) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('dynrel_trt'),{'hmc'});
        end

        function testResolveDynrelTrtForcedHmc(testCase)
            spec = localSpec(3,'dynrel_trt');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDynrelTrtForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'dynrel_trt');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: dod (case 9) registry + routing -------------------
        function testNativeSupportedDodIsHmc(testCase)
            % Four-event difference-of-differences (4-cell location-scale with 8-D
            % person/trial blocks) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('dod'),{'hmc'});
        end

        function testResolveDodForcedHmc(testCase)
            spec = localSpec(3,'dod');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDodForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'dod');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_dynrel_trt (case 15) registry + routing ------
        function testNativeSupportedDiffDynrelTrtIsHmc(testCase)
            % Group-level two-facet dynamic difference (six crossed 2x2 blocks +
            % dimension slopes) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_dynrel_trt'),{'hmc'});
        end

        function testResolveDiffDynrelTrtForcedHmc(testCase)
            spec = localSpec(3,'diff_dynrel_trt');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffDynrelTrtForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_dynrel_trt');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_dynrel_sserr_trt (case 16) registry + routing -
        function testNativeSupportedDiffDynrelSserrTrtIsHmc(testCase)
            % Subject-level two-facet dynamic difference (4-D person block + five
            % crossed 2x2 factors) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_dynrel_sserr_trt'),{'hmc'});
        end

        function testResolveDiffDynrelSserrTrtForcedHmc(testCase)
            spec = localSpec(3,'diff_dynrel_sserr_trt');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffDynrelSserrTrtForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_dynrel_sserr_trt');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_dynrel_rescor (case 17) registry + routing ---
        function testNativeSupportedDiffDynrelRescorIsHmc(testCase)
            % Concurrent one-facet group-level dynamic difference (bivariate
            % residual + dimension slopes) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_dynrel_rescor'),{'hmc'});
        end

        function testResolveDiffDynrelRescorForcedHmc(testCase)
            spec = localSpec(3,'diff_dynrel_rescor');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffDynrelRescorForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_dynrel_rescor');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_dynrel_trt_rescor (case 18) registry/routing -
        function testNativeSupportedDiffDynrelTrtRescorIsHmc(testCase)
            % Concurrent two-facet group-level dynamic difference (six crossed 2x2
            % factors + per-row bivariate residual) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_dynrel_trt_rescor'),{'hmc'});
        end

        function testResolveDiffDynrelTrtRescorForcedHmc(testCase)
            spec = localSpec(3,'diff_dynrel_trt_rescor');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffDynrelTrtRescorForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_dynrel_trt_rescor');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_dynrel_sserr_rescor (case 19) registry/routing
        function testNativeSupportedDiffDynrelSserrRescorIsHmc(testCase)
            % Concurrent one-facet subject-level dynamic difference (4-D person
            % block + per-row bivariate residual) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_dynrel_sserr_rescor'),{'hmc'});
        end

        function testResolveDiffDynrelSserrRescorForcedHmc(testCase)
            spec = localSpec(3,'diff_dynrel_sserr_rescor');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffDynrelSserrRescorForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_dynrel_sserr_rescor');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_dynrel_sserr_trt_rescor (case 20) routing -----
        function testNativeSupportedDiffDynrelSserrTrtRescorIsHmc(testCase)
            % Concurrent two-facet subject-level dynamic difference (4-D person
            % block + five crossed factors + per-row bivariate residual) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_dynrel_sserr_trt_rescor'),{'hmc'});
        end

        function testResolveDiffDynrelSserrTrtRescorForcedHmc(testCase)
            spec = localSpec(3,'diff_dynrel_sserr_trt_rescor');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffDynrelSserrTrtRescorForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_dynrel_sserr_trt_rescor');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_dynrel_sserr_rho (case 21) registry/routing --
        function testNativeSupportedDiffDynrelSserrRhoIsHmc(testCase)
            % Concurrent one-facet subject-level dynamic difference with a
            % per-subject residual correlation is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_dynrel_sserr_rho'),{'hmc'});
        end

        function testResolveDiffDynrelSserrRhoForcedHmc(testCase)
            spec = localSpec(3,'diff_dynrel_sserr_rho');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffDynrelSserrRhoForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_dynrel_sserr_rho');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_dynrel_sserr_trt_rho (case 22) routing -------
        function testNativeSupportedDiffDynrelSserrTrtRhoIsHmc(testCase)
            % Concurrent two-facet subject-level dynamic difference with a
            % per-subject residual correlation is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_dynrel_sserr_trt_rho'),{'hmc'});
        end

        function testResolveDiffDynrelSserrTrtRhoForcedHmc(testCase)
            spec = localSpec(3,'diff_dynrel_sserr_trt_rho');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffDynrelSserrTrtRhoForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_dynrel_sserr_trt_rho');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_sserr (case 8) registry + routing ------------
        function testNativeSupportedDiffSserrIsHmc(testCase)
            % Static subject-level two-event difference (three 2-D blocks +
            % per-subject rho, paired+unpaired) is HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_sserr'),{'hmc'});
        end

        function testResolveDiffSserrForcedHmc(testCase)
            spec = localSpec(3,'diff_sserr');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffSserrForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_sserr');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_trt (case 10, non-concurrent) registry/routing
        function testNativeSupportedDiffTrtIsHmc(testCase)
            % Static two-facet (test-retest) two-event difference, non-concurrent:
            % the KDIM=0 reduction of diff_dynrel_trt (case 15), HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_trt'),{'hmc'});
        end

        function testResolveDiffTrtForcedHmc(testCase)
            spec = localSpec(3,'diff_trt');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffTrtForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_trt');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native HMC: diff_trt_rescor (case 10, concurrent) registry/routing
        function testNativeSupportedDiffTrtRescorIsHmc(testCase)
            % Static concurrent two-facet two-event difference (paired bivariate
            % residual + rescor): the KDIM=0 reduction of diff_dynrel_trt_rescor
            % (case 18), HMC-only.
            testCase.verifyEqual(psyrat_native_supported('diff_trt_rescor'),{'hmc'});
        end

        function testResolveDiffTrtRescorForcedHmc(testCase)
            spec = localSpec(3,'diff_trt_rescor');
            av = localMockAvail(false,true,false);
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveDiffTrtRescorForcedFitlmeUnsupportedErrors(testCase)
            spec = localSpec(2,'diff_trt_rescor');
            av = localMockAvail(true,true,true);
            testCase.verifyError(@() psyrat_resolve_engine(2,spec,av), ...
                'PsyRAT:engineUnsupportedDesign');
        end

        % --- native: trt (case 5) registry + dual-engine routing -----------
        function testNativeSupportedTrtIsHmcAndFitlme(testCase)
            % Plain test-retest is the ONE design with both native engines; HMC
            % is listed first to document the cascade preference.
            testCase.verifyEqual(psyrat_native_supported('trt'),{'hmc','fitlme'});
        end

        function testResolveTrtForcedHmc(testCase)
            spec = localSpec(3,'trt');
            av = localMockAvail(false,true,false); % only hmc installed
            testCase.verifyEqual(psyrat_resolve_engine(3,spec,av),'hmc');
        end

        function testResolveTrtForcedFitlme(testCase)
            % trt SUPPORTS fitlme (unlike the location-scale designs), so a forced
            % fitlme request resolves rather than erroring.
            spec = localSpec(2,'trt');
            av = localMockAvail(false,false,true); % only fitlme installed
            testCase.verifyEqual(psyrat_resolve_engine(2,spec,av),'fitlme');
        end

        function testResolveTrtAutoPrefersHmcOverFitlme(testCase)
            spec = localSpec(0,'trt');
            av = localMockAvail(false,true,true); % no cmdstan, hmc+fitlme present
            testCase.verifyEqual(psyrat_resolve_engine(0,spec,av),'hmc');
        end

        function testResolveTrtAutoFitlmeWhenNoHmc(testCase)
            % With CmdStan and HMC both unavailable, trt still has fitlme as the
            % last-resort fallback.
            spec = localSpec(0,'trt');
            av = localMockAvail(false,false,true); % only fitlme installed
            testCase.verifyEqual(psyrat_resolve_engine(0,spec,av),'fitlme');
        end

        % --- duck-typed fit shim -------------------------------------------
        function testNativeFitExtractAndBlock(testCase)
            draws = struct('mu',(1:5)','sig_u',(6:10)');
            conv = {'name','n_eff','r_hat'; 'x',100,1.0};
            fit = PsyRATNativeFit(draws,conv,'fitlme',struct());
            fit.block(); % no-op, must not error
            testCase.verifyEqual(fit.extract('pars','mu').mu,(1:5)');
            testCase.verifyEqual(fit.extract('pars','sig_u').sig_u,(6:10)');
            testCase.verifyEqual(fit.convstats(0),conv);
            testCase.verifyEqual(fit.engine,'fitlme');
        end

        function testNativeFitMissingParErrors(testCase)
            fit = PsyRATNativeFit(struct('mu',(1:3)'),[],'fitlme',struct());
            testCase.verifyError(@() fit.extract('pars','sig_u'), ...
                'PsyRATNativeFit:missingPar');
        end

        % --- native fitlme one-facet engine --------------------------------
        function testFitlmeOneFacetMatchesDirectReml(testCase)
            testCase.assumeTrue(localHasFitlme(), ...
                'Statistics and Machine Learning Toolbox (fitlme) not available.');
            [data,Tbl] = localSimOneFacet(60,30,2.0,0.7,3.0,1.0,11);
            spec = localSpec(2,'sing'); spec.opts.bootstrap_reps = 100;
            fit = psyrat_run_native('fitlme',spec, ...
                'data',data,'seed',12345,'chains',4,'iter',200,'warmup',200);

            % the engine's stored point estimate must equal a direct REML fit
            % exactly (deterministic given the data), validating the extraction
            lme = fitlme(Tbl,'meas ~ 1 + (1|id) + (1|trl)');
            [psi,mse] = covarianceParameters(lme);
            fe = fixedEffects(lme);
            p = fit.meta.point;
            testCase.verifyEqual(p.sig_u,  sqrt(psi{1}(1,1)),'AbsTol',1e-9);
            testCase.verifyEqual(p.sig_trl,sqrt(psi{2}(1,1)),'AbsTol',1e-9);
            testCase.verifyEqual(p.sig_e,  sqrt(mse),        'AbsTol',1e-9);
            testCase.verifyEqual(p.mu,     fe(1),            'AbsTol',1e-9);

            % draws are aligned Bx1 SD vectors named as the case extracts them
            for f = {'mu','sig_u','sig_trl','sig_e'}
                v = fit.extract('pars',f{1}).(f{1});
                testCase.verifySize(v,[100 1]);
                testCase.verifyTrue(all(isfinite(v)));
            end
            % SD components are non-negative
            testCase.verifyTrue(all(fit.extract('pars','sig_u').sig_u >= 0));

            % convergence cell passes psyrat_checkconv's thresholds
            c = fit.convstats(0);
            testCase.verifyLessThan(c{2,3},1.1);        % r_hat
            testCase.verifyGreaterThanOrEqual(c{2,2},40); % n_eff >= 10*nchains
            testCase.verifyEqual(fit.engine,'fitlme');
        end

        function testFitlmeReproducibleWithSeed(testCase)
            testCase.assumeTrue(localHasFitlme(), ...
                'Statistics and Machine Learning Toolbox (fitlme) not available.');
            data = localSimOneFacet(40,20,2.0,0.7,3.0,1.0,7);
            spec = localSpec(2,'sing'); spec.opts.bootstrap_reps = 50;
            f1 = psyrat_run_native('fitlme',spec,'data',data,'seed',999, ...
                'chains',4,'iter',100,'warmup',100);
            f2 = psyrat_run_native('fitlme',spec,'data',data,'seed',999, ...
                'chains',4,'iter',100,'warmup',100);
            testCase.verifyEqual(f1.extract('pars','sig_u').sig_u, ...
                f2.extract('pars','sig_u').sig_u);
        end

        function testFitlmeFacetEmitsSuffixedDraws(testCase)
            % cases 2-3 (group/event): the engine fits each stratum and emits
            % per-stratum suffixed draws (mu_G1, sig_u_G2, ...) matching the
            % multi-stratum case extraction.
            testCase.assumeTrue(localHasFitlme(), ...
                'Statistics and Machine Learning Toolbox (fitlme) not available.');
            data = struct();
            for i = 1:2
                d = localSimOneFacet(30,15,2.0,0.7,3.0,1.0,20+i);
                data.(sprintf('id_G%d',i))   = d.id;
                data.(sprintf('trl_G%d',i))  = d.trl;
                data.(sprintf('meas_G%d',i)) = d.meas;
            end
            spec = localSpec(2,'sing_facet');
            spec.design.suffix = 'G';
            spec.design.nstrata = 2;
            spec.opts.bootstrap_reps = 60;
            fit = psyrat_run_native('fitlme',spec,'data',data,'seed',5, ...
                'chains',4,'iter',100,'warmup',100);
            for i = 1:2
                for f = {'mu','sig_u','sig_trl','sig_e'}
                    nm = sprintf('%s_G%d',f{1},i);
                    v = fit.extract('pars',nm).(nm);
                    testCase.verifySize(v,[60 1]);
                    testCase.verifyTrue(all(isfinite(v)));
                end
            end
            c = fit.convstats(0);
            testCase.verifyLessThan(c{2,3},1.1);
            testCase.verifyEqual(fit.engine,'fitlme');
        end

        function testFitlmeFacetRequiresSpecFields(testCase)
            testCase.assumeTrue(localHasFitlme(), ...
                'Statistics and Machine Learning Toolbox (fitlme) not available.');
            data = struct('id_G1',1,'trl_G1',1,'meas_G1',1);
            spec = localSpec(2,'sing_facet'); % missing suffix/nstrata
            testCase.verifyError(@() psyrat_run_native('fitlme',spec, ...
                'data',data,'seed',1,'chains',4,'iter',50,'warmup',50), ...
                'PsyRAT:fitlmeFacetSpec');
        end

        function testFitlmeSplitsRecoversAndRequiresWeight(testCase)
            % nonparallel splits (case 23): emits mu/sig_u/sig_trl/sig_splitxid/
            % sig_err draws and requires a weight field.
            testCase.assumeTrue(localHasFitlme(), ...
                'Statistics and Machine Learning Toolbox (fitlme) not available.');
            data = localSimSplits(40,6,12);
            spec = localSpec(2,'splits'); spec.opts.bootstrap_reps = 40;
            ws = warning('off','PsyRAT:fitlmeSplitsIdentifiability');
            cleanup = onCleanup(@() warning(ws));
            fit = psyrat_run_native('fitlme',spec,'data',data,'seed',3, ...
                'chains',4,'iter',100,'warmup',100);
            for f = {'mu','sig_u','sig_trl','sig_splitxid','sig_err'}
                v = fit.extract('pars',f{1}).(f{1});
                testCase.verifySize(v,[40 1]);
                testCase.verifyTrue(all(isfinite(v)));
            end
            testCase.verifyEqual(fit.engine,'fitlme');
            % weight-required guard
            nodata = rmfield(data,'weight');
            testCase.verifyError(@() psyrat_run_native('fitlme',spec, ...
                'data',nodata,'seed',3,'chains',4,'iter',50,'warmup',50), ...
                'PsyRAT:fitlmeSplitsWeight');
        end

        function testFitlmeSplitsTrtRecoversAndRequiresWeight(testCase)
            % splits test-retest (case 24): emits pop_int + 8 SD components and
            % requires a weight field.
            testCase.assumeTrue(localHasFitlme(), ...
                'Statistics and Machine Learning Toolbox (fitlme) not available.');
            data = localSimSplitsTrt(15,2,4,4);
            spec = localSpec(2,'splits_trt'); spec.opts.bootstrap_reps = 30;
            ws = warning('off','PsyRAT:fitlmeSplitsIdentifiability');
            w2 = warning('off','PsyRAT:fitlmeTrtFewOccasions');
            cleanup = onCleanup(@() warning([ws w2]));
            fit = psyrat_run_native('fitlme',spec,'data',data,'seed',1, ...
                'chains',4,'iter',60,'warmup',60);
            nm = {'pop_int','sig_id','sig_occ','sig_trl','sig_trlxid', ...
                'sig_occxid','sig_trlxocc','sig_posxid','sig_err'};
            for k = 1:numel(nm)
                v = fit.extract('pars',nm{k}).(nm{k});
                testCase.verifySize(v,[30 1]);
                testCase.verifyTrue(all(isfinite(v)));
            end
            nodata = rmfield(data,'weight');
            testCase.verifyError(@() psyrat_run_native('fitlme',spec, ...
                'data',nodata,'seed',1,'chains',4,'iter',50,'warmup',50), ...
                'PsyRAT:fitlmeSplitsTrtWeight');
        end

        function testHmcHarnessRecoversAndConverges(testCase)
            % Phase 2 sampling harness (psyrat_hmc_sample): recovers a known
            % Gaussian's mean/SD, reports R-hat near 1 with a usable n_eff, and
            % is reproducible for a fixed seed.
            testCase.assumeTrue(exist('hmcSampler','file') > 0, ...
                'hmcSampler (Statistics and Machine Learning Toolbox) not available.');
            rng(1); N = 200; y = 3.0 + 2.0*randn(N,1);
            f = @(th) localGaussLp(th,y,N);
            [draws,conv] = psyrat_hmc_sample(f,[0;0],2,300,400,123,[1 2],{'mu','logsig'});
            testCase.verifySize(draws,[800 2]);
            testCase.verifyEqual(mean(draws(:,1)),3.0,'AbsTol',0.4);
            testCase.verifyEqual(mean(exp(draws(:,2))),2.0,'AbsTol',0.5);
            testCase.verifyLessThan(conv{2,3},1.1);          % r_hat
            testCase.verifyGreaterThan(conv{2,2},40);        % n_eff
            [d2,~] = psyrat_hmc_sample(f,[0;0],2,300,400,123,[1 2],{'mu','logsig'});
            testCase.verifyEqual(draws,d2);                  % reproducible
        end

        function testFitlmeUnsupportedDesignErrors(testCase)
            testCase.assumeTrue(localHasFitlme(), ...
                'Statistics and Machine Learning Toolbox (fitlme) not available.');
            data = localSimOneFacet(20,10,2.0,0.7,3.0,1.0,3);
            spec = localSpec(2,'sserr');
            testCase.verifyError(@() psyrat_run_native('fitlme',spec, ...
                'data',data,'seed',1,'chains',4,'iter',50,'warmup',50), ...
                'PsyRAT:fitlmeUnsupportedDesign');
        end

        function testFitlmeTrtRecoversAndGuardsWeight(testCase)
            % test-retest fitter: recovers the 7 components + intercept, matches
            % a direct REML fit, and refuses a weighted (splits) data struct.
            testCase.assumeTrue(localHasFitlme(), ...
                'Statistics and Machine Learning Toolbox (fitlme) not available.');
            data = localSimTrt(25,2,12,5);
            spec = localSpec(2,'trt'); spec.opts.bootstrap_reps = 40;
            fit = psyrat_run_native('fitlme',spec,'data',data,'seed',7, ...
                'chains',4,'iter',100,'warmup',100);
            names = {'pop_int','sig_id','sig_occ','sig_trl','sig_trlxid', ...
                'sig_occxid','sig_trlxocc','sig_err'};
            for k = 1:numel(names)
                v = fit.extract('pars',names{k}).(names{k});
                testCase.verifySize(v,[40 1]);
                testCase.verifyTrue(all(isfinite(v)));
            end
            % point estimate equals a direct REML fit (deterministic extraction)
            Tbl = table(data.meas(:),categorical(data.id(:)),categorical(data.occ(:)), ...
                categorical(data.trl(:)),categorical(data.trlxid(:)), ...
                categorical(data.occxid(:)),categorical(data.trlxocc(:)), ...
                'VariableNames',{'meas','id','occ','trl','trlxid','occxid','trlxocc'});
            lme = fitlme(Tbl,['meas ~ 1 + (1|id)+(1|occ)+(1|trl)+(1|trlxid)+' ...
                '(1|occxid)+(1|trlxocc)']);
            [psi,mse] = covarianceParameters(lme);
            testCase.verifyEqual(fit.meta.point.sig_id, sqrt(psi{1}(1,1)),'AbsTol',1e-9);
            testCase.verifyEqual(fit.meta.point.sig_err, sqrt(mse),'AbsTol',1e-9);
            % weight-guard: a splits data struct must be refused, not silently fit
            wdata = data; wdata.weight = ones(numel(data.meas),1);
            testCase.verifyError(@() psyrat_run_native('fitlme',spec, ...
                'data',wdata,'seed',7,'chains',4,'iter',50,'warmup',50), ...
                'PsyRAT:fitlmeTrtWeighted');
        end

        % --- registry <-> dispatch consistency -----------------------------
        function testRegistryDispatchConsistency(testCase)
            %Guard the registry <-> dispatch agreement: every design key
            %psyrat_native_supported recognizes must be dispatchable by the
            %engine(s) it lists, and the engines must not silently lack a case.
            %The same key string is spelled at the psyrat_run_stan case site, in
            %the registry, and in the psyrat_native_hmc / psyrat_native_fitlme
            %switches; this catches drift here rather than mid-run on a long fit.
            %
            %THIS TEST COVERS THE LAST TWO OF THOSE THREE. The CASE SITE leg is
            %pinned by TestStanModelSnapshots.testNativeKeyMatchesCaseSite, which
            %captures the key each analysis actually passes and checks it against
            %a golden map AND against the registry list below. Both are needed: a
            %key mistyped at the case site is not an error at run time, it simply
            %returns no engine, so without that test a design silently becomes
            %CmdStan-only with every suite still green.
            %Each dispatcher is probed with empty data and the key is asserted
            %RECOGNIZED - it fails with a data/field error, NOT the
            %"unsupported/not-implemented" guard - without running a real fit.
            keys = psyrat_native_supported();
            testCase.verifyNotEmpty(keys, ...
                'psyrat_native_supported() must enumerate the design keys');
            for i = 1:numel(keys)
                key = keys{i};
                engines = psyrat_native_supported(key);
                testCase.verifyNotEmpty(engines, sprintf( ...
                    'registry returned no engine for enumerated key ''%s''',key));
                if ismember('hmc',engines)
                    testCase.verifyFalse(localDispatchRejects( ...
                        @() psyrat_native_hmc(struct('key',key),struct(), ...
                            12345,4,50,50,struct()),'PsyRAT:hmcNotImplemented'), ...
                        sprintf(['registry lists hmc for ''%s'' but ' ...
                        'psyrat_native_hmc has no matching dispatch case'],key));
                end
                if ismember('fitlme',engines)
                    testCase.verifyFalse(localDispatchRejects( ...
                        @() psyrat_native_fitlme(struct('key',key),struct(), ...
                            12345,4,struct('bootstrap_reps',100)), ...
                        'PsyRAT:fitlmeUnsupportedDesign'), ...
                        sprintf(['registry lists fitlme for ''%s'' but ' ...
                        'psyrat_native_fitlme has no matching dispatch case'],key));
                end
            end
        end

        function testResolveAutoTooOldCmdstanErrors(testCase)
            %A present-but-below-minimum CmdStan on the auto path must error with
            %an upgrade message (matching the forced-CmdStan preflight) rather
            %than silently falling back to an approximate native engine, even when
            %a native engine is available.
            spec = localSpec(0,'sing');
            av = localMockAvail(false,true,true);   % cmdstan unavailable...
            av.cmdstan_status = 'tooold';           % ...because it is too old
            testCase.verifyError(@() psyrat_resolve_engine(0,spec,av), ...
                'PsyRAT:cmdstanVersion');
        end

    end
end

% =============================== helpers =====================================
function spec = localSpec(engine,key)
spec = struct('engine',engine, ...
    'design',struct('key',key), ...
    'opts',struct('bootstrap_reps',100,'hmc_gradient','auto'));
end

function av = localMockAvail(cmdstan,hmc,fitlme)
%Build a mock availability struct (ad irrelevant to resolution). cmdstan_status
%mirrors the production field: 'ok' when CmdStan is available, 'absent' otherwise
%(so the auto cascade's too-old gate is not triggered by these mocks; the gate
%test sets cmdstan_status = 'tooold' explicitly).
if cmdstan
    cmdstan_status = 'ok';
else
    cmdstan_status = 'absent';
end
av = struct('cmdstan',logical(cmdstan),'cmdstan_status',cmdstan_status, ...
    'hmc',logical(hmc),'fitlme',logical(fitlme),'ad',false, ...
    'native',logical(hmc)||logical(fitlme), ...
    'details',struct('cmdstan','cmdstan(mock)','hmc','hmc(mock)', ...
        'fitlme','fitlme(mock)','ad','ad(mock)'));
end

function tf = localDispatchRejects(fn,unsupportedId)
%Run fn (a native-engine call with deliberately empty data) and report whether it
%failed with the "engine does not implement this design" guard. A RECOGNIZED key
%reaches its builder and fails with some OTHER error (a missing data field), so
%this returns false; an UNRECOGNIZED key hits the dispatch otherwise branch
%(unsupportedId), so this returns true.
tf = false;
try
    fn();
catch err
    tf = strcmp(err.identifier,unsupportedId);
end
end

function tf = localHasFitlme()
tf = (exist('fitlme','file') > 0) && license('test','Statistics_Toolbox');
end

function [lp,g] = localGaussLp(th,y,N)
%Log-density + gradient of y ~ N(mu, exp(logsig)) over theta = [mu; logsig],
%used to exercise the HMC harness against a target with a known answer.
mu = th(1); ls = th(2); s2 = exp(2*ls); r = y - mu;
lp = -0.5*N*log(2*pi) - N*ls - 0.5*sum(r.^2)/s2;
g = [ sum(r)/s2; -N + sum(r.^2)/s2 ];
end

function data = localSimSplitsTrt(NSUB,NOCC,NSPL,seed)
%Simulate a splits test-retest dataset (person x occasion x split; split-mean
%score with residual var sig_err^2/n_i) and build the case-24 Stan data struct.
rng(seed,'twister');
[ID,OCC,SP] = ndgrid(1:NSUB,1:NOCC,1:NSPL);
id = ID(:); occ = OCC(:); spl = SP(:); N = numel(id);
[~,~,txi] = unique([spl id],'rows');
[~,~,oxi] = unique([occ id],'rows');
[~,~,txo] = unique([spl occ],'rows');
w = randi([6 24],N,1);
uid=2.0*randn(NSUB,1); uocc=0.5*randn(NOCC,1); usp=0.6*randn(NSPL,1);
meas = 1.0 + uid(id) + uocc(occ) + usp(spl) + (3.0./sqrt(w)).*randn(N,1);
data = struct('NOBS',N,'NSUB',NSUB,'NOCC',NOCC,'NTRL',NSPL,'id',id,'occ',occ, ...
    'trl',spl,'trlxid',txi,'occxid',oxi,'trlxocc',txo,'meas',meas,'weight',w);
end

function data = localSimSplits(NSUB,NSPL,seed)
%Simulate a nonparallel-splits dataset (one row per person x split; the score is
%a split mean of n_i items so its residual variance is sig_err^2 / n_i) and build
%the case-23 Stan data struct.
rng(seed,'twister');
[ID,SP] = ndgrid(1:NSUB,1:NSPL);
id = ID(:); spl = SP(:); N = numel(id);
sigP=2.0; sigS=0.5; sigSX=0.8; sigE=3.0; mu=1.0;
uid = sigP*randn(NSUB,1); us = sigS*randn(NSPL,1); usx = sigSX*randn(N,1);
w = randi([5 30],N,1);
meas = mu + uid(id) + us(spl) + usx + (sigE./sqrt(w)).*randn(N,1);
data = struct('NOBS',N,'NSUB',NSUB,'NTRL',NSPL,'id',id,'trl',spl, ...
    'meas',meas,'weight',w);
end

function data = localSimTrt(NSUB,NOCC,NTRL,seed)
%Simulate a balanced test-retest (subject x occasion x trial) dataset and build
%the case-5 Stan data struct (id/occ/trl + the three interaction-level indices).
rng(seed,'twister');
[ID,OCC,TRL] = ndgrid(1:NSUB,1:NOCC,1:NTRL);
id = ID(:); occ = OCC(:); trl = TRL(:);
N = numel(id);
[~,~,trlxid]  = unique([trl id],'rows');
[~,~,occxid]  = unique([occ id],'rows');
[~,~,trlxocc] = unique([trl occ],'rows');
sP=2.0; sO=0.5; sT=0.6; sTI=0.4; sOI=0.4; sTO=0.3; sE=1.5; mu=1.0;
uid=sP*randn(NSUB,1); uocc=sO*randn(NOCC,1); utrl=sT*randn(NTRL,1);
utrlxid=sTI*randn(max(trlxid),1); uoccxid=sOI*randn(max(occxid),1);
utrlxocc=sTO*randn(max(trlxocc),1);
meas = mu + uid(id) + uocc(occ) + utrl(trl) + utrlxid(trlxid) + ...
    uoccxid(occxid) + utrlxocc(trlxocc) + sE*randn(N,1);
data = struct('NOBS',N,'NSUB',NSUB,'NOCC',NOCC,'NTRL',NTRL, ...
    'id',id,'occ',occ,'trl',trl,'trlxid',trlxid,'occxid',occxid, ...
    'trlxocc',trlxocc,'meas',meas);
end

function [data,Tbl] = localSimOneFacet(NSUB,NTRL,sigU,sigT,sigE,mu,seed)
%Simulate a balanced one-facet (person x trial) dataset and return both the
%Stan-style data struct the case builds and the long table.
rng(seed,'twister');
u = sigU*randn(NSUB,1);
t = sigT*randn(NTRL,1);
id = repelem((1:NSUB)',NTRL);
trl = repmat((1:NTRL)',NSUB,1);
meas = mu + u(id) + t(trl) + sigE*randn(NSUB*NTRL,1);
data = struct('NOBS',numel(meas),'NSUB',NSUB,'NTRL',NTRL, ...
    'id',id,'trl',trl,'meas',meas);
Tbl = table(meas,categorical(id),categorical(trl), ...
    'VariableNames',{'meas','id','trl'});
end
