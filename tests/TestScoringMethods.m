classdef TestScoringMethods < PsyRATTestBase
    % Ground-truth (characterization) tests for the single-trial scoring
    % engine in psyrat_project_score.m. These lock the *current* behavior of
    % the methods, polarity handling, and boundary policy so that later changes
    % (new methods) cannot silently alter existing scores. They mirror the
    % ERP-012 synthetic pattern in TestImportWorkflow.m
    % but isolate the scoring methods rather than the import/lifecycle path.
    %
    % Convention reminder: scores are computed per channel and then averaged
    % across channels (S7/ERP-005, owner-kept). Single-channel cases below test
    % the per-method math directly; the multichannel cases lock the averaging.

    methods (Test)

        function testCentroidLatencyValue(testCase)
            % centroid_latency = amplitude-weighted mean of the in-window
            % sample times. With weight_mode='absolute' the weights are
            % abs(seg). Asymmetric waveform chosen so the centroid is not the
            % peak latency, exercising the weighting.
            sig = zeros(1,1,5);
            sig(1,1,:) = [0 1 4 2 0];           % weights at 0,10,20,30,40 ms
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','Centroid', ...
                'method','centroid_latency', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1, ...
                'method_params',struct('weight_mode','absolute'));
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);

            % sum(t.*w)/sum(w) = (10*1 + 20*4 + 30*2) / (1+4+2) = 150/7.
            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,150/7,'AbsTol',1e-12);
            testCase.verifyEqual(row.latency_ms,150/7,'AbsTol',1e-12);
        end

        function testTimeWindowMeanAmplitudeValue(testCase)
            % Mean amplitude over the window; latency is undefined (NaN).
            sig = zeros(1,1,5);
            sig(1,1,:) = [0 2 4 2 0];
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','Mean', ...
                'method','time_window_mean_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1);
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);

            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,mean([0 2 4 2 0]),'AbsTol',1e-12); % 1.6
            testCase.verifyTrue(isnan(row.latency_ms));
        end

        function testNegativePolarityAbsolutePeak(testCase)
            % Negative-polarity absolute peak returns the most negative sample.
            sig = zeros(1,1,5);
            sig(1,1,:) = [0 -1 -5 -2 0];        % trough at 20 ms
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','NegPeak', ...
                'method','absolute_peak_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','negative');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);
            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,-5,'AbsTol',1e-12);
            testCase.verifyEqual(row.latency_ms,20,'AbsTol',1e-12);
        end

        function testBothPolarityAbsolutePeak(testCase)
            % polarity='both' selects by largest magnitude and returns the
            % signed value. Here the negative deflection is largest.
            sig = zeros(1,1,5);
            sig(1,1,:) = [0 3 0 -6 0];          % +3 at 10 ms, -6 at 30 ms
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','BothPeak', ...
                'method','absolute_peak_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','both');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);
            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,-6,'AbsTol',1e-12);
            testCase.verifyEqual(row.latency_ms,30,'AbsTol',1e-12);
        end

        function testBothPolarityTieHonorsTieBreakPolicy(testCase)
            % Workstream J item 3: polarity='both' selects by largest absolute
            % deviation. On a magnitude tie it must honor tie_break_policy like
            % the +/- paths do (pre-fix it always took the first via max()).
            % [0 5 0 -5 0]: +5 at 10 ms and -5 at 30 ms tie at |5|.
            sig = zeros(1,1,5);
            sig(1,1,:) = [0 5 0 -5 0];
            p0 = testCase.makeProject(sig,(0:10:40),{'Cz'});

            % Default policy ('first') -> earliest tie: +5 at 10 ms.
            specFirst = psyrat_scoring_spec_create('name','BothTieFirst', ...
                'method','absolute_peak_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','both');
            p = psyrat_project_add_scoring_spec('project',p0,'spec',specFirst);
            [p,~] = psyrat_project_score('project',p,'spec_id',specFirst.spec_id);
            rowF = p.scoring_results(end,:);
            testCase.verifyEqual(rowF.value,5,'AbsTol',1e-12);
            testCase.verifyEqual(rowF.latency_ms,10,'AbsTol',1e-12);

            % Policy 'last' -> latest tie: -5 at 30 ms (this is the J.3 fix).
            specLast = psyrat_scoring_spec_create('name','BothTieLast', ...
                'method','absolute_peak_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','both','tie_break_policy','last');
            p = psyrat_project_add_scoring_spec('project',p0,'spec',specLast);
            [p,~] = psyrat_project_score('project',p,'spec_id',specLast.spec_id);
            rowL = p.scoring_results(end,:);
            testCase.verifyEqual(rowL.value,-5,'AbsTol',1e-12);
            testCase.verifyEqual(rowL.latency_ms,30,'AbsTol',1e-12);
        end

        function testEdgePeakEmitsRunLevelWarning(testCase)
            % Workstream J item 2: an absolute peak that lands on the window
            % boundary (here a monotonic ramp -> max at the last sample) is a
            % valid score but signals a possibly too-narrow window. The scorer
            % emits one run-level warning and still scores the trial (no
            % exclusion). A centered peak must not warn.
            sig = zeros(1,1,5);
            sig(1,1,:) = [1 2 3 4 5];           % monotonic -> peak at 40 ms (edge)
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','EdgePeak', ...
                'method','absolute_peak_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','positive');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = testCase.verifyWarning( ...
                @() psyrat_project_score('project',p,'spec_id',spec.spec_id), ...
                'psyrat_scoring:edgePeak');
            testCase.verifyEqual(rpt.n_failed,0);
            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,5,'AbsTol',1e-12);
            testCase.verifyEqual(row.latency_ms,40,'AbsTol',1e-12);

            % A centered peak (max at the middle sample) must NOT warn.
            sig2 = zeros(1,1,5);
            sig2(1,1,:) = [0 2 9 2 0];          % peak at 20 ms (interior)
            p2 = testCase.makeProject(sig2,(0:10:40),{'Cz'});
            spec2 = psyrat_scoring_spec_create('name','MidPeak', ...
                'method','absolute_peak_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','positive');
            p2 = psyrat_project_add_scoring_spec('project',p2,'spec',spec2);
            testCase.verifyWarningFree( ...
                @() psyrat_project_score('project',p2,'spec_id',spec2.spec_id));
        end

        function testAllNaNWindowPeakMethodsFailCleanly(testCase)
            % Workstream J item 1: a fully-NaN measurement window must fail the
            % trial cleanly (excluded), not error. All polarities route through
            % localTie, which returns no index for an all-NaN window, so pre-fix
            % the scorer indexed seg(NaN) and threw; the guard returns a clean
            % failure instead. 'both' is included because J.3 folded it through
            % the same localTie/guard path (its tie set is empty for all-NaN).
            sig = NaN(1,1,5);
            p0 = testCase.makeProject(sig,(0:10:40),{'Cz'});

            combos = { ...
                'absolute_peak_amplitude','positive'; ...
                'absolute_peak_amplitude','negative'; ...
                'absolute_peak_amplitude','both'; ...
                'peak_latency','positive'; ...
                'peak_latency','negative'; ...
                'peak_latency','both'};

            for k = 1:size(combos,1)
                spec = psyrat_scoring_spec_create('name','NaNwin', ...
                    'method',combos{k,1}, ...
                    'window_start_ms',0,'window_end_ms',40, ...
                    'channel_selector',1,'polarity',combos{k,2});
                p = psyrat_project_add_scoring_spec('project',p0,'spec',spec);

                % Must not error (pre-fix this threw at seg(NaN)).
                [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);

                ctx = sprintf('%s/%s',combos{k,1},combos{k,2});
                testCase.verifyEqual(rpt.n_failed,1, ...
                    sprintf('%s should fail on an all-NaN window',ctx));
                row = p.scoring_results(end,:);
                testCase.verifyEqual(char(row.status),'failed',ctx);
                testCase.verifyTrue(row.excluded_from_handoff,ctx);
                testCase.verifyNotEmpty(char(row.status_reason),ctx);
            end
        end

        function testBoundaryPolicyStrictErrorsOutOfBounds(testCase)
            % Strict boundary policy rejects a window that extends past the
            % available time axis (error raised before any trial is scored).
            sig = zeros(1,1,5);
            sig(1,1,:) = [0 1 2 1 0];
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','StrictOOB', ...
                'method','time_window_mean_amplitude', ...
                'window_start_ms',-50,'window_end_ms',200, ...   % below min, above max
                'channel_selector',1,'boundary_policy','strict');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            testCase.verifyError(@() psyrat_project_score('project',p,'spec_id',spec.spec_id), ...
                'psyrat_scoring:windowOutOfBounds');
        end

        function testBoundaryPolicyNonStrictClampsWindow(testCase)
            % A non-strict boundary policy skips the bounds check and scores
            % only the in-range samples. Window [-50 25] selects 0,10,20 ms.
            sig = zeros(1,1,5);
            sig(1,1,:) = [1 2 3 99 99];          % only first three are in range
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','LaxClamp', ...
                'method','time_window_mean_amplitude', ...
                'window_start_ms',-50,'window_end_ms',25, ...
                'channel_selector',1,'boundary_policy','lax');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);
            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,mean([1 2 3]),'AbsTol',1e-12); % 2
        end

        % --- New methods: fractional area latency -----------------------------

        function testFractionalAreaLatencyFiftyPercent(testCase)
            % Symmetric triangle -> 50% area latency is the center sample.
            % cumtrapz of [0 2 4 2 0] over 0:10:40 is [0 10 40 70 80];
            % half of total (80) is 40, reached exactly at 20 ms.
            sig = zeros(1,1,5);
            sig(1,1,:) = [0 2 4 2 0];
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','FAL50', ...
                'method','fractional_area_latency', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','positive');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);
            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,20,'AbsTol',1e-12);
            testCase.verifyEqual(row.latency_ms,20,'AbsTol',1e-12);
        end

        function testFractionalAreaLatencyCustomFractionInterpolates(testCase)
            % area_fraction = 0.25 -> target area 20 falls between cum(2)=10
            % and cum(3)=40, i.e. one third of the way from 10 ms to 20 ms.
            sig = zeros(1,1,5);
            sig(1,1,:) = [0 2 4 2 0];
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','FAL25', ...
                'method','fractional_area_latency', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','positive', ...
                'method_params',struct('area_fraction',0.25));
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);
            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,10 + 10/3,'AbsTol',1e-12); % 13.333...
        end

        function testFractionalAreaLatencyNegativePolarity(testCase)
            % A negative-going triangle scored with negative polarity yields
            % the same area magnitudes as the positive case, so the 50% area
            % latency is again the center (20 ms).
            sig = zeros(1,1,5);
            sig(1,1,:) = [0 -2 -4 -2 0];
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','FALneg', ...
                'method','fractional_area_latency', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','negative');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);
            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,20,'AbsTol',1e-12);
        end

        function testFractionalAreaLatencyZeroAreaFailsCleanly(testCase)
            % Positive polarity on an all-negative waveform leaves no positive
            % area; the trial must fail cleanly (excluded), not error.
            sig = zeros(1,1,5);
            sig(1,1,:) = [-1 -2 -3 -2 -1];
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','FALzero', ...
                'method','fractional_area_latency', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','positive');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,1);
            row = p.scoring_results(end,:);
            testCase.verifyEqual(char(row.status),'failed');
            testCase.verifyTrue(row.excluded_from_handoff);
        end

        % --- New methods: integrated area amplitude ---------------------------

        function testIntegratedAreaAmplitudeValue(testCase)
            % Integral of a constant 2 uV over a 40 ms window = 2*40 = 80 uV*ms.
            sig = zeros(1,1,5);
            sig(1,1,:) = [2 2 2 2 2];
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','Area', ...
                'method','integrated_area_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','positive');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);
            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,80,'AbsTol',1e-12);
            testCase.verifyTrue(isnan(row.latency_ms));
        end

        function testIntegratedAreaAmplitudeNegativePolaritySign(testCase)
            % Negative polarity preserves the sign of the area (-80 uV*ms).
            sig = zeros(1,1,5);
            sig(1,1,:) = [-2 -2 -2 -2 -2];
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            spec = psyrat_scoring_spec_create('name','AreaNeg', ...
                'method','integrated_area_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','negative');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);
            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,-80,'AbsTol',1e-12);
        end

        function testIntegratedAreaAmplitudeBothRectifies(testCase)
            % A square biphasic wave [2 -2 2 -2 2]: 'both' rectifies (full 80),
            % 'positive' keeps only the positive half-samples (40).
            sig = zeros(1,1,5);
            sig(1,1,:) = [2 -2 2 -2 2];
            p = testCase.makeProject(sig,(0:10:40),{'Cz'});

            specBoth = psyrat_scoring_spec_create('name','AreaBoth', ...
                'method','integrated_area_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','both');
            p = psyrat_project_add_scoring_spec('project',p,'spec',specBoth);
            [p,~] = psyrat_project_score('project',p,'spec_id',specBoth.spec_id);
            testCase.verifyEqual(p.scoring_results(end,:).value,80,'AbsTol',1e-12);

            specPos = psyrat_scoring_spec_create('name','AreaPos', ...
                'method','integrated_area_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',1,'polarity','positive');
            p = psyrat_project_add_scoring_spec('project',p,'spec',specPos);
            [p,~] = psyrat_project_score('project',p,'spec_id',specPos.spec_id);
            testCase.verifyEqual(p.scoring_results(end,:).value,40,'AbsTol',1e-12);
        end

        function testNewMethodsMultiChannelAveraging(testCase)
            % Per-channel-then-average (S7/ERP-005) holds for the new methods:
            % channel areas 80 and 160 uV*ms average to 120.
            sig = zeros(1,2,5);
            sig(1,1,:) = [2 2 2 2 2];   % area 80
            sig(1,2,:) = [4 4 4 4 4];   % area 160
            p = testCase.makeProject(sig,(0:10:40),{'Fz','Cz'});

            spec = psyrat_scoring_spec_create('name','AreaROI', ...
                'method','integrated_area_amplitude', ...
                'window_start_ms',0,'window_end_ms',40, ...
                'channel_selector',[1 2],'polarity','positive');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, rpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(rpt.n_failed,0);
            testCase.verifyEqual(p.scoring_results(end,:).value,120,'AbsTol',1e-12);
        end

        function testFractionalAreaLatencyDefaultsAndValidates(testCase)
            % area_fraction defaults to 0.5 when omitted and is rejected when
            % outside the open interval (0,1) or non-scalar.
            spec = psyrat_scoring_spec_create('name','FALdefault', ...
                'method','fractional_area_latency', ...
                'window_start_ms',0,'window_end_ms',40,'channel_selector',1);
            testCase.verifyEqual(spec.method_params.area_fraction,0.5);

            mk = @(f) psyrat_scoring_spec_create('name','x', ...
                'method','fractional_area_latency', ...
                'window_start_ms',0,'window_end_ms',40,'channel_selector',1, ...
                'method_params',struct('area_fraction',f));
            testCase.verifyError(@() mk(0),'psyrat_scoring:areaLatencyParams');
            testCase.verifyError(@() mk(1),'psyrat_scoring:areaLatencyParams');
            testCase.verifyError(@() mk(1.5),'psyrat_scoring:areaLatencyParams');
            testCase.verifyError(@() mk([0.2 0.5]),'psyrat_scoring:areaLatencyParams');
        end

    end

    methods (Access = private)
        function p = makeProject(~,sig,time_ms,labels)
            % Build a minimal scoring-ready project around a synthetic
            % trial x channel x sample signal array. Every trial is assigned
            % subject 's1' / event 'A'; callers that need richer metadata can
            % overwrite p.canonical.trial_table afterward.
            p = psyrat_project_new('name','Scoring Method Test');
            nTrials = size(sig,1);

            dt = median(diff(double(time_ms(:))));
            fs = 1000/dt;

            p.canonical.has_single_trial = true;
            p.canonical.single_trial.signals = sig;
            p.canonical.single_trial.time_ms = time_ms;
            p.canonical.single_trial.fs_hz = fs;
            p.canonical.single_trial.amplitude_unit = 'uV';
            p.canonical.single_trial.channel_labels = labels;
            p.canonical.single_trial.data_hash = psyrat_hash_struct(struct('sig',sig));
            p.canonical.trial_table = psyrat_build_trial_table(nTrials, ...
                'subject_id',repmat({'s1'},nTrials,1), ...
                'event_id',repmat({'A'},nTrials,1), ...
                'source_ref',repmat({'synthetic'},nTrials,1));
        end
    end
end
