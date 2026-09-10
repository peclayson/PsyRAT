classdef TestImportWorkflow < PsyRATTestBase

    methods (Test)
        function testProjectLifecycleSaveLoadAndValidate(testCase)
            p = psyrat_project_new('name','Unit Test Project');

            rpt = psyrat_project_validate(p);
            testCase.verifyTrue(rpt.ok);

            outfile = [tempname '.psyratproj'];
            psyrat_project_save('project',p,'file',outfile);
            testCase.verifyTrue(exist(outfile,'file')==2);

            p2 = psyrat_project_load('file',outfile);
            testCase.verifyEqual(p2.project_id,p.project_id);
            testCase.verifyEqual(p2.schema_name,'PsyRATImportProject');
        end

        function testLoadMigratesLegacyBaselineSpecField(testCase)
            % Regression (J.4 follow-up): projects saved by older builds carry
            % an obsolete baseline_policy field on their scoring specs. Load must
            % strip it, otherwise adding a new (field-free) spec fails struct-
            % array concatenation ("Subscripted assignment between dissimilar
            % structures").
            p = psyrat_project_new('name','Legacy Spec Migration');
            spec1 = psyrat_scoring_spec_create('name','Old', ...
                'method','time_window_mean_amplitude', ...
                'window_start_ms',0,'window_end_ms',400,'channel_selector',1);
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec1);

            % Simulate a legacy save: re-attach the removed field to the spec.
            p.scoring_specs(1).baseline_policy = struct('mode','window','window_ms',[-200 0]);

            outfile = [tempname '.psyratproj'];
            testCase.addTeardown(@() delete(outfile));
            psyrat_project_save('project',p,'file',outfile);
            pLoaded = psyrat_project_load('file',outfile);

            % Migration dropped the obsolete field on load.
            testCase.verifyFalse(isfield(pLoaded.scoring_specs,'baseline_policy'));

            % Adding a new spec to the migrated project must not error.
            spec2 = psyrat_scoring_spec_create('name','New', ...
                'method','time_window_mean_amplitude', ...
                'window_start_ms',0,'window_end_ms',400,'channel_selector',1);
            pAdded = psyrat_project_add_scoring_spec('project',pLoaded,'spec',spec2);
            testCase.verifyEqual(numel(pAdded.scoring_specs),2);
        end

        function testScoringExclusionsAndHandoffExport(testCase)
            p = psyrat_project_new('name','Score Export Test');

            % trial x channel x sample
            sig = zeros(3,1,7);
            sig(1,1,:) = [0 1 3 1 0 -1 0];
            sig(2,1,:) = [0 0 0 0 0 0 0]; % no local peak
            sig(3,1,:) = [0 -1 -3 -1 0 1 0];

            p.canonical.has_single_trial = true;
            p.canonical.single_trial.signals = sig;
            p.canonical.single_trial.time_ms = (-300:100:300);
            p.canonical.single_trial.fs_hz = 10;
            p.canonical.single_trial.amplitude_unit = 'uV';
            p.canonical.single_trial.channel_labels = {'Cz'};
            p.canonical.single_trial.data_hash = psyrat_hash_struct(struct('sig',sig));
            p.canonical.trial_table = psyrat_build_trial_table(3, ...
                'subject_id',{'s1';'s1';'s2'}, ...
                'event_id',{'A';'A';'A'}, ...
                'source_ref',{'synthetic';'synthetic';'synthetic'});

            spec = psyrat_scoring_spec_create('name','Local Peak', ...
                'method','local_peak_amplitude', ...
                'window_start_ms',-200, ...
                'window_end_ms',200, ...
                'channel_selector',1, ...
                'polarity','positive', ...
                'method_params',struct('neighborhood_samples',1));
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, runRpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(runRpt.n_trials,3);
            testCase.verifyEqual(runRpt.n_failed,1);

            [tbl, meta, p] = psyrat_project_export_proc_table('project',p,'spec_id',spec.spec_id, ...
                'validate_with_preflight',true,'diffest',1,'sserrvar',1);

            testCase.verifyEqual(height(tbl),2);
            testCase.verifyFalse(any(~isfinite(tbl.meas)));
            testCase.verifyEqual(meta.n_excluded_trials,1);
            testCase.verifyEqual(meta.n_included_trials,2);
            testCase.verifyEqual(p.ui_state.active_handoff_spec_id,spec.spec_id);
        end

        function testDefaultExportOfEarlierSpecAfterScoringAnotherSpec(testCase)
            % REGRESSION for finding G35. active_run_id is spec-agnostic
            % (scoring ANY spec overwrites it), and the export default used it
            % unconditionally, so a default export of spec A after scoring spec
            % B filtered A's non-stale rows by B's run id and died
            % psyrat_handoff:runNotFound. The default must instead fall back to
            % the selected spec's latest surviving run. This is the first test
            % to score two specs and export; before the fix the first export
            % below errored.
            p = psyrat_project_new('name','Two Spec Export Test');

            sig = zeros(3,1,7);
            sig(1,1,:) = [0 1 3 1 0 -1 0];
            sig(2,1,:) = [0 2 1 0 0 0 0];
            sig(3,1,:) = [0 -1 -3 -1 0 1 0];

            p.canonical.has_single_trial = true;
            p.canonical.single_trial.signals = sig;
            p.canonical.single_trial.time_ms = (-300:100:300);
            p.canonical.single_trial.fs_hz = 10;
            p.canonical.single_trial.amplitude_unit = 'uV';
            p.canonical.single_trial.channel_labels = {'Cz'};
            p.canonical.single_trial.data_hash = psyrat_hash_struct(struct('sig',sig));
            p.canonical.trial_table = psyrat_build_trial_table(3, ...
                'subject_id',{'s1';'s1';'s2'}, ...
                'event_id',{'A';'A';'A'}, ...
                'source_ref',{'synthetic';'synthetic';'synthetic'});

            specA = psyrat_scoring_spec_create('name','Early Window', ...
                'method','time_window_mean_amplitude', ...
                'window_start_ms',-200,'window_end_ms',0,'channel_selector',1);
            specB = psyrat_scoring_spec_create('name','Late Window', ...
                'method','time_window_mean_amplitude', ...
                'window_start_ms',0,'window_end_ms',200,'channel_selector',1);
            p = psyrat_project_add_scoring_spec('project',p,'spec',specA);
            p = psyrat_project_add_scoring_spec('project',p,'spec',specB);

            [p, rptA] = psyrat_project_score('project',p,'spec_id',specA.spec_id);
            [p, rptB] = psyrat_project_score('project',p,'spec_id',specB.spec_id);
            % the precondition that made the defect fire: the active run now
            % belongs to spec B
            testCase.assertEqual(p.ui_state.active_run_id, rptB.run_id);

            % default export of the EARLIER spec must resolve to that spec's
            % own latest surviving run, not the active (spec B) run
            [tblA, metaA, p] = psyrat_project_export_proc_table('project',p, ...
                'spec_id',specA.spec_id);
            testCase.verifyEqual(metaA.run_id, rptA.run_id, ...
                'A default export must use the selected spec''s latest run.');
            testCase.verifyEqual(height(tblA), 3);

            % the active-run path is untouched when it DOES match the spec
            [~, metaB] = psyrat_project_export_proc_table('project',p, ...
                'spec_id',specB.spec_id);
            testCase.verifyEqual(metaB.run_id, rptB.run_id);

            % an explicitly requested run is not rescued: asking for a run the
            % spec does not have must still error, since the caller named it
            testCase.verifyError(@() psyrat_project_export_proc_table( ...
                'project',p,'spec_id',specA.spec_id,'run_id',rptB.run_id), ...
                'psyrat_handoff:runNotFound');
        end

        function testMultiChannelScoresAreAveragedAcrossPerChannelScores(testCase)
            p = psyrat_project_new('name','MultiChannel Score Test');

            % 1 trial x 2 channels x 5 samples, time: 0,10,20,30,40 ms
            sig = zeros(1,2,5);
            sig(1,1,:) = [0 5 0 0 0]; % peak at 10 ms
            sig(1,2,:) = [0 0 0 4 0]; % peak at 30 ms

            p.canonical.has_single_trial = true;
            p.canonical.single_trial.signals = sig;
            p.canonical.single_trial.time_ms = (0:10:40);
            p.canonical.single_trial.fs_hz = 100;
            p.canonical.single_trial.amplitude_unit = 'uV';
            p.canonical.single_trial.channel_labels = {'Fz','Cz'};
            p.canonical.single_trial.data_hash = psyrat_hash_struct(struct('sig',sig));
            p.canonical.trial_table = psyrat_build_trial_table(1, ...
                'subject_id',{'s1'},'event_id',{'A'},'source_ref',{'synthetic'});

            spec = psyrat_scoring_spec_create('name','Peak Lat ROI', ...
                'method','peak_latency', ...
                'window_start_ms',0, ...
                'window_end_ms',40, ...
                'channel_selector',[1 2], ...
                'polarity','positive', ...
                'method_params',struct());
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, runRpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(runRpt.n_failed,0);

            row = p.scoring_results(end,:);
            % Per-channel latencies: 10 and 30 ms -> averaged score 20 ms.
            testCase.verifyEqual(row.value,20,'AbsTol',1e-12);
            testCase.verifyEqual(row.latency_ms,20,'AbsTol',1e-12);
            testCase.verifyFalse(row.excluded_from_handoff);
        end



        function testMultiChannelAmplitudeIsAverageOfPerChannelAmplitudes(testCase)
            p = psyrat_project_new('name','MultiChannel Amplitude Test');

            % 1 trial x 2 channels x 5 samples, time: 0,10,20,30,40 ms
            % Channel peaks are 2 and 6 uV; expected averaged score = 4 uV.
            sig = zeros(1,2,5);
            sig(1,1,:) = [0 2 0 0 0];
            sig(1,2,:) = [0 0 6 0 0];

            p.canonical.has_single_trial = true;
            p.canonical.single_trial.signals = sig;
            p.canonical.single_trial.time_ms = (0:10:40);
            p.canonical.single_trial.fs_hz = 100;
            p.canonical.single_trial.amplitude_unit = 'uV';
            p.canonical.single_trial.channel_labels = {'Fz','Cz'};
            p.canonical.single_trial.data_hash = psyrat_hash_struct(struct('sig',sig));
            p.canonical.trial_table = psyrat_build_trial_table(1, ...
                'subject_id',{'s1'},'event_id',{'A'},'source_ref',{'synthetic'});

            spec = psyrat_scoring_spec_create('name','Abs Peak ROI', ...
                'method','absolute_peak_amplitude', ...
                'window_start_ms',0, ...
                'window_end_ms',40, ...
                'channel_selector',[1 2], ...
                'polarity','positive', ...
                'method_params',struct());
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);

            [p, runRpt] = psyrat_project_score('project',p,'spec_id',spec.spec_id);
            testCase.verifyEqual(runRpt.n_failed,0);

            row = p.scoring_results(end,:);
            testCase.verifyEqual(row.value,4,'AbsTol',1e-12);
            testCase.verifyEqual(row.latency_ms,15,'AbsTol',1e-12);
            testCase.verifyFalse(row.excluded_from_handoff);
        end

        function testImportEpToolkitMatSingleTrial(testCase)
            EPdata = struct();
            EPdata.dataType = 'single_trial';
            EPdata.Fs = 1000;
            EPdata.timeNames = [0 1 2 3 4];
            EPdata.chanNames = {'Fz';'Cz'};
            EPdata.cellNames = {'cond1';'cond2'};
            EPdata.subNames = {'sub1'};
            % chan x time x wave x sub
            EPdata.data = zeros(2,5,2,1);
            EPdata.data(:,:,1,1) = [0 1 2 1 0; 0 0 1 0 0];
            EPdata.data(:,:,2,1) = [0 -1 -2 -1 0; 0 0 -1 0 0];

            infile = [tempname '.mat'];
            save(infile,'EPdata','-v7');

            p = psyrat_project_new('name','Import MAT');
            [p, ir] = psyrat_project_import_file('project',p,'file',infile);

            testCase.verifyTrue(ir.single_trial_available);
            testCase.verifyTrue(p.canonical.has_single_trial);
            testCase.verifyEqual(size(p.canonical.single_trial.signals,1),2);
            testCase.verifyEqual(height(p.canonical.trial_table),2);
            testCase.verifyTrue(p.canonical.has_average);
        end

        function testImportErplabErpIsReferenceOnlyAndNotScoreable(testCase)
            % GAP-2 (.erp half): ERPLAB .erp ERPsets store AVERAGED waveforms
            % (channel x sample x bin), not single trials, so they are NOT
            % SUPPORTED for generalizability-theory scoring. They must import as
            % reference-only (average available, single-trial unavailable, with
            % a warning) and the scoring pipeline must refuse them. This guards
            % against an .erp source ever being treated as scoreable. (Real
            % owner-supplied ERPsets are verified in gui_review_artifacts/
            % gap2_erp_drive.m; this is the hermetic synthetic guard.)
            ERP = struct();
            ERP.datatype = 'ERP';
            ERP.nchan = 2;
            ERP.nbin = 2;
            ERP.pnts = 5;
            ERP.srate = 1000;
            ERP.xmin = 0;
            ERP.xmax = 0.004;
            ERP.times = [0 1 2 3 4];
            % channel x sample x bin: two AVERAGED condition waveforms.
            ERP.bindata = zeros(2,5,2);
            ERP.bindata(:,:,1) = [0 1 2 1 0; 0 0 1 0 0];   % "correct" average
            ERP.bindata(:,:,2) = [0 -1 -2 -1 0; 0 0 -1 0 0]; % "error" average
            ERP.binerror = zeros(2,5,2);                     % per-bin SEM
            ERP.bindescr = {'correct','error'};
            ERP.chanlocs = struct('labels',{'Fz','Cz'});
            ERP.ntrials = struct('accepted',[10 8],'rejected',[0 0],'invalid',[0 0]);

            infile = [tempname '.erp'];
            save(infile,'ERP','-v7');

            det = psyrat_import_detect_format(infile);
            testCase.verifyEqual(det.adapter,'erplab');

            p = psyrat_project_new('name','Import ERPLAB Average');
            [p, ir] = psyrat_project_import_file('project',p,'file',infile);
            testCase.verifyFalse(ir.single_trial_available);
            testCase.verifyTrue(ir.average_available);
            testCase.verifyFalse(p.canonical.has_single_trial);
            testCase.verifyTrue(p.canonical.has_average);
            testCase.verifyNotEmpty(ir.warnings);   % reference-only warning

            % The GT scoring pipeline must refuse averaged-only data.
            spec = psyrat_scoring_spec_create('name','TWM', ...
                'method','time_window_mean_amplitude', ...
                'window_start_ms',0,'window_end_ms',4, ...
                'channel_selector',1,'polarity','negative');
            p = psyrat_project_add_scoring_spec('project',p,'spec',spec);
            testCase.verifyError(@() psyrat_project_score('project',p,'spec_id',spec.spec_id), ...
                'psyrat_scoring:singleTrialRequired');
        end

        function testImportRealEpToolkitMatExactScoring(testCase)
            % GAP-2 (EP .mat half): a REAL EP Toolkit single-trial export carried
            % with a literal .mat extension must (a) route to the ep_toolkit_mat
            % adapter purely by file content/extension, and (b) score
            % time-window-mean amplitude exactly against the raw epochs. The
            % committed fixture is a small slice of real recorded flanker data
            % (8 ch x 151 samp x 5 trials), produced by
            % psyrat_make_ep_real_subset_fixture. Because the data are real,
            % the fixture is deliberately excluded from public releases; this
            % test therefore ASSUMES (skips, not fails) when it is absent.
            fixture = fullfile(localTestDataDir(), 'ep_toolkit_real_subset.mat');
            testCase.assumeTrue(exist(fixture, 'file') == 2, ...
                sprintf(['Real-EP .mat regression fixture not present: %s\n', ...
                'It is development-only and not distributed in public releases. ', ...
                'On a development checkout, regenerate with ', ...
                'psyrat_make_ep_real_subset_fixture.'], fixture));

            % Detection must select the EP adapter from a .mat-extension file.
            det = psyrat_import_detect_format(fixture);
            testCase.verifyEqual(det.adapter, 'ep_toolkit_mat');
            testCase.verifyEqual(det.display_name, 'EP Toolkit MAT');

            % Run the real import chain the GUI uses.
            p = psyrat_project_new('name', 'Import Real EP MAT');
            [p, ir] = psyrat_project_import_file('project', p, 'file', fixture);
            testCase.verifyTrue(ir.single_trial_available);
            testCase.verifyTrue(p.canonical.has_single_trial);

            st = p.canonical.single_trial;
            testCase.verifyEqual(size(st.signals, 1), 5);          % 5 single trials
            testCase.verifyEqual(numel(st.channel_labels), 8);     % E1..E8 retained

            % Score time-window-mean amplitude 0-100 ms at E6 (the canonical
            % FCz scoring channel, preserved at index 6 in the subset).
            chan = 6;
            spec = psyrat_scoring_spec_create('name', 'TWM', ...
                'method', 'time_window_mean_amplitude', ...
                'window_start_ms', 0, 'window_end_ms', 100, ...
                'channel_selector', chan, 'polarity', 'negative');
            p = psyrat_project_add_scoring_spec('project', p, 'spec', spec);
            [p, runRpt] = psyrat_project_score('project', p, 'spec_id', spec.spec_id);
            testCase.verifyEqual(runRpt.n_failed, 0);

            % Exact gate: each scored value equals the manual mean of the
            % fixture's own real samples in the 0-100 ms window, matched by
            % trial_uid. (Observed difference is floating-point zero.)
            sr = p.scoring_results;
            sr = sr(strcmp(sr.spec_id, spec.spec_id) & ~sr.is_stale, :);
            tt = p.canonical.trial_table;
            [tf, loc] = ismember(tt.trial_uid, sr.trial_uid);
            scored = nan(height(tt), 1);
            scored(tf) = sr.value(loc(tf));

            tmask = st.time_ms >= 0 & st.time_ms <= 100;
            manual = squeeze(mean(double(st.signals(:, chan, tmask)), 3));
            testCase.verifyEqual(scored, manual, 'AbsTol', 1e-9);

            % Handoff table exports and passes preflight validation.
            [tbl, meta] = psyrat_project_export_proc_table('project', p, ...
                'spec_id', spec.spec_id, 'validate_with_preflight', true, ...
                'diffest', 1, 'sserrvar', 1);
            testCase.verifyEqual(height(tbl), 5);
            testCase.verifyFalse(any(~isfinite(tbl.meas)));
            testCase.verifyEqual(meta.n_included_trials, 5);
        end

        function testFormatDetectionByExtension(testCase)
            setfile = [tempname '.set']; fclose(fopen(setfile,'w'));
            erpfile = [tempname '.erp']; fclose(fopen(erpfile,'w'));
            matfile = [tempname '.mat']; save(matfile,'setfile');
            eptfile = [tempname '.ept']; save(eptfile,'setfile');

            d1 = psyrat_import_detect_format(setfile);
            d2 = psyrat_import_detect_format(erpfile);
            d3 = psyrat_import_detect_format(matfile);
            d4 = psyrat_import_detect_format(eptfile);

            testCase.verifyEqual(d1.adapter,'eeglab');
            testCase.verifyEqual(d2.adapter,'erplab');
            testCase.verifyTrue(any(strcmp(d3.adapter,{'ep_toolkit_mat','eeglab','erplab'})) || ~d3.ok);
            testCase.verifyEqual(d4.adapter,'ep_toolkit_mat');
        end

        function testSingleTrialScoringMethodsMatchGroundTruth(testCase)
            % Ground-truth checks for each single-trial scoring method on
            % synthetic single-channel waveforms (absolute_peak_amplitude is
            % already covered by the multi-channel tests).
            times = (0:10:40);  % 0 10 20 30 40 ms

            r = localScoreSingleChannel(times, [0 2 4 6 8], ...
                'time_window_mean_amplitude', 'positive', struct());
            testCase.verifyEqual(r.value, 4, 'AbsTol', 1e-12);          % window mean

            r = localScoreSingleChannel(times, [0 5 1 0 0], ...
                'peak_latency', 'positive', struct());
            testCase.verifyEqual(r.value, 10, 'AbsTol', 1e-12);         % peak latency
            testCase.verifyEqual(r.latency_ms, 10, 'AbsTol', 1e-12);

            r = localScoreSingleChannel(times, [0 1 5 1 0], ...
                'local_peak_amplitude', 'positive', ...
                struct('neighborhood_samples', 1));
            testCase.verifyEqual(r.value, 5, 'AbsTol', 1e-12);          % local maximum
            testCase.verifyEqual(r.latency_ms, 20, 'AbsTol', 1e-12);

            r = localScoreSingleChannel(times, [0 10 10 0 0], ...
                'centroid_latency', 'positive', ...
                struct('weight_mode', 'positive_only'));
            testCase.verifyEqual(r.value, 15, 'AbsTol', 1e-12);         % weighted centroid
            testCase.verifyEqual(r.latency_ms, 15, 'AbsTol', 1e-12);
        end
    end
end

function row = localScoreSingleChannel(times, vals, method, polarity, methodParams)
% Build a one-trial, one-channel project, score it with the given spec, and
% return the resulting scoring-results row.
sig = zeros(1, 1, numel(vals));
sig(1, 1, :) = vals;

p = psyrat_project_new('name', 'ScoringGroundTruth');
p.canonical.has_single_trial = true;
p.canonical.single_trial.signals = sig;
p.canonical.single_trial.time_ms = times(:)';
p.canonical.single_trial.fs_hz = 100;
p.canonical.single_trial.amplitude_unit = 'uV';
p.canonical.single_trial.channel_labels = {'Fz'};
p.canonical.single_trial.data_hash = psyrat_hash_struct(struct('sig', sig));
p.canonical.trial_table = psyrat_build_trial_table(1, ...
    'subject_id', {'s1'}, 'event_id', {'A'}, 'source_ref', {'synthetic'});

spec = psyrat_scoring_spec_create('name', method, 'method', method, ...
    'window_start_ms', times(1), 'window_end_ms', times(end), ...
    'channel_selector', 1, 'polarity', polarity, 'method_params', methodParams);
p = psyrat_project_add_scoring_spec('project', p, 'spec', spec);
[p, ~] = psyrat_project_score('project', p, 'spec_id', spec.spec_id);
row = p.scoring_results(end, :);
end

function d = localTestDataDir()
% Resolve the root-level test_data/ dir from this test file's location
% (tests/TestImportWorkflow.m -> tests -> project root -> test_data).
d = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'test_data');
end
