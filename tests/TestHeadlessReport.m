classdef TestHeadlessReport < PsyRATTestBase
    %Tests for psyrat_report (headless, render-free reporting) and its helpers
    %psyrat_write_table and psyrat_provenance_lines. These run in the standard
    %unit suite: no CmdStan, no figures, no save dialogs.

    methods (Test)

        function testReportReturnsSingTables(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();

            report = psyrat_report(psyrat_data);

            testCase.verifyEqual(report.analysis, 'sing');
            testCase.verifyTrue(isfield(report.tables, 'cutoff'));
            testCase.verifyTrue(isfield(report.tables, 'overall'));
            testCase.verifyTrue(istable(report.tables.cutoff));
            testCase.verifyTrue(istable(report.tables.overall));
            % The returned tables match what the GUI builders produce headlessly.
            expectedCutoff = psyrat_depcutofft('psyrat_data', psyrat_data, 'gui', 0);
            testCase.verifyEqual(report.tables.cutoff, expectedCutoff);
            % Provenance lines are present and include the seed and convergence.
            testCase.verifyFalse(isempty(report.provenance));
            joined = strjoin(report.provenance, '|');
            testCase.verifySubstring(joined, 'Seed:');
            testCase.verifySubstring(joined, 'Analysis:');
            testCase.verifySubstring(joined, 'Converged:');
        end

        function testReportOpensNoFigures(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            nBefore = numel(findall(0, 'Type', 'figure'));
            psyrat_report(psyrat_data);
            nAfter = numel(findall(0, 'Type', 'figure'));
            testCase.verifyEqual(nAfter, nBefore, ...
                'psyrat_report must not open any figures.');
        end

        function testReportWritesFilesWithProvenance(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            tdir = localTempDir(testCase);

            report = psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv');

            % One file per table, all present on disk.
            testCase.verifyEqual(numel(report.files), numel(fieldnames(report.tables)));
            for k = 1:numel(report.files)
                testCase.verifyEqual(exist(report.files{k}, 'file'), 2);
            end

            % The written file carries the provenance header above the table.
            txt = fileread(report.files{1});
            testCase.verifySubstring(txt, 'Seed:');
            testCase.verifySubstring(txt, 'Analysis:');
            testCase.verifySubstring(txt, 'Converged:');
            testCase.verifySubstring(txt, 'PsyRAT Toolbox');
        end

        function testReportWritesXlsx(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            tdir = localTempDir(testCase);
            report = psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.xlsx');
            testCase.verifyFalse(isempty(report.files));
            for k = 1:numel(report.files)
                testCase.verifyEqual(exist(report.files{k}, 'file'), 2);
                [~, ~, ext] = fileparts(report.files{k});
                testCase.verifyEqual(ext, '.xlsx');
            end
        end

        function testReportDoDTables(testCase)
            % DoD reads rel directly (no relsummary) and is render-free.
            psyrat_data = PsyRATTestDataFactory.makeDoDRelData();
            nBefore = numel(findall(0, 'Type', 'figure'));
            report = psyrat_report(psyrat_data);
            nAfter = numel(findall(0, 'Type', 'figure'));

            testCase.verifyEqual(nAfter, nBefore, 'DoD report must not open figures.');
            testCase.verifyEqual(report.analysis, 'dod');
            testCase.verifyTrue(isfield(report.tables, 'dod_observed'));
            testCase.verifyTrue(isfield(report.tables, 'dod_varcomp'));
            testCase.verifyTrue(istable(report.tables.dod_observed));
            % The fixture has two groups -> one contrast row per group, with a
            % prepended Group column identifying each.
            testCase.verifyEqual(height(report.tables.dod_observed), 2);
            testCase.verifyTrue(ismember('Group', ...
                report.tables.dod_observed.Properties.VariableNames));
            testCase.verifyEqual(report.tables.dod_observed.Group, ...
                {'Control'; 'Clinical'});
        end

        function testReportDoDHonorsExplicitCi(testCase)
            % G52: the DoD route dropped an explicit 'CI' on the floor -- the
            % two DoD builders carry their own .95 defaults and psyrat_report
            % never threaded ciperc into local_dod_tables, and with no
            % relsummary there is no stored width, so neither the provenance
            % line nor the mismatch note could disclose the drop. The widths
            % compared here are far apart (.50 vs the .95 default) on
            % purpose: at small draw counts MATLAB's quantile clamps extreme
            % quantiles to min/max, so ADJACENT widths can render identically
            % and the differential below could not fail.
            psyrat_data = PsyRATTestDataFactory.makeDoDRelData();

            wide = psyrat_report(psyrat_data);
            explicit95 = psyrat_report(psyrat_data, 'CI', .95);
            narrow = psyrat_report(psyrat_data, 'CI', .50);

            % default and explicit .95 must agree exactly, pinning that the
            % default path did not move
            testCase.verifyEqual(explicit95.tables.dod_varcomp, ...
                wide.tables.dod_varcomp, ...
                'explicit CI .95 must reproduce the default tables');

            wWide = wide.tables.dod_varcomp.High_CrI - ...
                wide.tables.dod_varcomp.Low_CrI;
            wNarrow = narrow.tables.dod_varcomp.High_CrI - ...
                narrow.tables.dod_varcomp.Low_CrI;
            % precondition: the fixture's draws vary, so the .95 intervals
            % have real width for the differential below to bite on
            testCase.assertGreaterThan(max(wWide), 0, ...
                'fixture precondition broken: all DoD varcomp draws identical');
            testCase.verifyTrue(all(wNarrow <= wWide + 1e-12), ...
                'a .50 interval must never be wider than its .95 interval');
            testCase.verifyLessThan(sum(wNarrow), sum(wWide), ...
                ['an explicit CI .50 must narrow the DoD varcomp intervals; ', ...
                'equality means the option is still being dropped']);

            % the OBSERVED table is equally CI-governed (adversarial wave:
            % a regression threading ciperc into only one builder must not
            % pass) -- its intervals are string-formatted, so pin by
            % equality/difference of the rendered tables
            testCase.verifyEqual(explicit95.tables.dod_observed, ...
                wide.tables.dod_observed, ...
                'explicit CI .95 must reproduce the default observed table');
            testCase.verifyNotEqual(narrow.tables.dod_observed, ...
                wide.tables.dod_observed, ...
                ['an explicit CI .50 must change the DoD observed table''s ', ...
                'rendered intervals; equality means only one builder got ciperc']);
        end

        function testReportDoDExportHeaderRecordsTheWidth(testCase)
            % G52 disclosure (adversarial wave): DoD has no relsummary, so
            % the provenance block cannot record an interval width -- yet
            % the width is now user-settable, and two exports at different
            % widths must not be byte-indistinguishable in their headers.
            % The report writes a 'Credible interval' header line with the
            % width it actually applied.
            psyrat_data = PsyRATTestDataFactory.makeDoDRelData();
            tdir = localTempDir(testCase);
            report = psyrat_report(psyrat_data, 'CI', .50, ...
                'outdir', tdir, 'format', '.csv');
            testCase.assertFalse(isempty(report.files));
            txt = fileread(report.files{1});
            testCase.verifySubstring(txt, 'Credible interval: 50%', ...
                'the DoD export header must record the honored width');
        end

        function testReportDoDWritesFiles(testCase)
            psyrat_data = PsyRATTestDataFactory.makeDoDRelData();
            tdir = localTempDir(testCase);
            report = psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv');
            testCase.verifyEqual(numel(report.files), 2);
            for k = 1:numel(report.files)
                testCase.verifyEqual(exist(report.files{k}, 'file'), 2);
            end
            txt = fileread(report.files{1});
            testCase.verifySubstring(txt, 'Seed:');
            testCase.verifySubstring(txt, 'Analysis:');
        end

        function testReportRequiresRelsummary(testCase)
            psyrat_data = struct();
            psyrat_data.rel.analysis = 'ic';
            testCase.verifyError(@() psyrat_report(psyrat_data), ...
                'psyrat_report:norelsummary');
        end

        function testReportUnsupportedAnalysisErrors(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            % An unrecognized analysis string (family '', no dynrel/splits/dodiff
            % substring): passes the relsummary guard then hits the inner switch's
            % otherwise. sserr is now supported (see the sserr tests below), so this
            % uses a synthetic unknown label to keep exercising the guard.
            psyrat_data.rel.analysis = 'ic_notarealvariant';
            testCase.verifyError(@() psyrat_report(psyrat_data), ...
                'psyrat_report:unsupported');
        end

        function testVariancetGuiZeroOpensNoFigure(testCase)
            % C3 step 1: psyrat_variancet('gui',0) must return the table without
            % opening a figure (so psyrat_report can build the sserr group-level
            % table headlessly).
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            nBefore = numel(findall(0, 'Type', 'figure'));
            tbl = psyrat_variancet('psyrat_data', psyrat_data, 'gui', 0);
            nAfter = numel(findall(0, 'Type', 'figure'));
            testCase.verifyEqual(nAfter, nBefore, ...
                'psyrat_variancet(''gui'',0) must not open a figure.');
            testCase.verifyTrue(istable(tbl));
        end

        function testReportSserrTables(testCase)
            % Subject-level error-variance family. Build a real sserr relsummary
            % through the pipeline (makeICDiffSSErrRelData + psyrat_relsummary), so
            % the per-event ssrel_table and the location-scale data are genuine,
            % then report headlessly.
            psyrat_data = localMakeSserrData();
            nBefore = numel(findall(0, 'Type', 'figure'));
            report = psyrat_report(psyrat_data);
            nAfter = numel(findall(0, 'Type', 'figure'));

            testCase.verifyEqual(nAfter, nBefore, ...
                'sserr report must not open figures.');
            testCase.verifyEqual(report.analysis, 'sserr');
            % Group-level SD/SEM/ICC table (render-free psyrat_variancet).
            testCase.verifyTrue(isfield(report.tables, 'variance'));
            testCase.verifyTrue(istable(report.tables.variance));
            % Per-subject coefficients flattened into one long table.
            testCase.verifyTrue(isfield(report.tables, 'sscoeffs'));
            sc = report.tables.sscoeffs;
            testCase.verifyTrue(istable(sc));
            testCase.verifyGreaterThan(height(sc), 0);
            expectedCols = {'Group','Event','id','ind2include','trls', ...
                'dep_pt','dep_ll','dep_ul','icc_pt','icc_ll','icc_ul', ...
                'sem_pt','sem_ll','sem_ul','ss_errvar'};
            testCase.verifyEqual(sc.Properties.VariableNames, expectedCols);
            % Two events in the fixture -> both strata stacked.
            testCase.verifyEqual(numel(unique(sc.Event)), 2);
            % Coefficient label note is surfaced for the export header.
            testCase.verifyTrue(isfield(report, 'sscoeff_note'));
            testCase.verifySubstring(report.sscoeff_note, 'coefficient');
            % No reliability-vs-n curve for sserr.
            testCase.verifyFalse(isfield(report.tables, 'curve'));
        end

        function testReportTrtSserrTables(testCase)
            % Subject-level CROSSED TEST-RETEST error-variance family
            % (trt_sserrvar, analysis 25). Build a real relsummary through the
            % pipeline (makeTRTSSErrRelData + psyrat_relsummary 'trt_sserr') then
            % report headlessly. It routes to the same 'sserr' report family as
            % the single-occasion ic_sserrvar, with the crossed facet SDs.
            psyrat_data = localMakeTrtSserrData();
            nBefore = numel(findall(0, 'Type', 'figure'));
            report = psyrat_report(psyrat_data, 'reltype', 3);
            nAfter = numel(findall(0, 'Type', 'figure'));

            testCase.verifyEqual(nAfter, nBefore, ...
                'trt_sserr report must not open figures.');
            testCase.verifyEqual(report.analysis, 'sserr');
            testCase.verifyEqual(report.rel_analysis, 'trt_sserrvar');
            % Group-level SD/SEM/ICC table (render-free psyrat_variancet).
            testCase.verifyTrue(isfield(report.tables, 'variance'));
            testCase.verifyTrue(istable(report.tables.variance));
            % Per-subject coefficients flattened into one long table.
            testCase.verifyTrue(isfield(report.tables, 'sscoeffs'));
            sc = report.tables.sscoeffs;
            testCase.verifyTrue(istable(sc));
            expectedCols = {'Group','Event','id','ind2include','trls', ...
                'dep_pt','dep_ll','dep_ul','icc_pt','icc_ll','icc_ul', ...
                'sem_pt','sem_ll','sem_ul','ss_errvar'};
            testCase.verifyEqual(sc.Properties.VariableNames, expectedCols);
            % Three participants, single measure (no events) in the fixture.
            testCase.verifyEqual(height(sc), 3);
            % Coefficient label note is surfaced for the export header.
            testCase.verifyTrue(isfield(report, 'sscoeff_note'));
            testCase.verifySubstring(report.sscoeff_note, 'coefficient');
            % No reliability-vs-n curve for sserr.
            testCase.verifyFalse(isfield(report.tables, 'curve'));
        end

        function testReportSserrWritesFilesWithProvenance(testCase)
            psyrat_data = localMakeSserrData();
            tdir = localTempDir(testCase);
            report = psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv');
            testCase.verifyEqual(numel(report.files), numel(fieldnames(report.tables)));
            for k = 1:numel(report.files)
                testCase.verifyEqual(exist(report.files{k}, 'file'), 2);
            end
            % The sscoeffs file records the coefficient-label note + provenance.
            idx = find(contains(report.files, 'sscoeffs'), 1);
            testCase.verifyNotEmpty(idx);
            txt = fileread(report.files{idx});
            testCase.verifySubstring(txt, 'Seed:');
            testCase.verifySubstring(txt, 'Analysis:');
            testCase.verifySubstring(txt, 'coefficient');
        end

        function testReportCurveTableSing(testCase)
            % One-facet reliability-vs-n curve, tabulated from the calc layer.
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            report = psyrat_report(psyrat_data);
            testCase.verifyTrue(isfield(report.tables, 'curve'));
            crv = report.tables.curve;
            testCase.verifyTrue(istable(crv));
            % 1 group x 1 event x 50 trials (default ntrials).
            testCase.verifyEqual(height(crv), 50);
            testCase.verifyEqual(crv.n_trials, (1:50)');
            testCase.verifyEqual(crv.Properties.VariableNames, ...
                {'Group','Event','n_trials', ...
                'Dependability_LowerLimit','Dependability_PointEstimate', ...
                'Dependability_UpperLimit'});
            % Reliability is a proper coefficient in [0,1] and non-decreasing in n.
            testCase.verifyGreaterThanOrEqual(min(crv.Dependability_PointEstimate), 0);
            testCase.verifyLessThanOrEqual(max(crv.Dependability_PointEstimate), 1);
            testCase.verifyGreaterThanOrEqual( ...
                crv.Dependability_PointEstimate(end), crv.Dependability_PointEstimate(1));
            testCase.verifyTrue(isfield(report, 'curve_note'));
        end

        function testReportCurveTableTrt(testCase)
            % Test-retest reliability-vs-n curve (psyrat_rel_trt path).
            psyrat_data = PsyRATTestDataFactory.makeTRTSummaryData();
            report = psyrat_report(psyrat_data, 'ntrials', 20);
            testCase.verifyTrue(isfield(report.tables, 'curve'));
            crv = report.tables.curve;
            testCase.verifyEqual(height(crv), 20);
            testCase.verifyEqual(crv.n_trials, (1:20)');
        end

        function testReportCurveTableTrtHonorsOccasionEstimand(testCase)
            % RC-10. The headless reliability-vs-n curve passed no nocc, so it
            % silently used n'_o = 1 while the summary tables in the same export
            % used the k the user selected. Both artifacts answer "how many
            % trials do I need?" under an identical "Dependability" header, so a
            % reader had no way to see they disagreed.
            psyrat_data = PsyRATTestDataFactory.makeTRTSummaryData();
            rs = psyrat_data.relsummary;
            % The point-estimate column is named for the selected coefficient.
            if rs.gcoeff == 2
                ptcol = 'Generalizability_PointEstimate';
            else
                ptcol = 'Dependability_PointEstimate';
            end

            % The fixture carries no nocc field at all, so this baseline is the
            % legacy fallback path (a file saved before the option existed).
            testCase.assertFalse(isfield(rs, 'nocc'));
            single = psyrat_report(psyrat_data, 'ntrials', 10);

            multi_data = psyrat_data;
            multi_data.relsummary.noccmode = 2;
            multi_data.relsummary.nocc = 4;
            multi = psyrat_report(multi_data, 'ntrials', 10);

            % Non-vacuity: the multi-occasion curve must actually MOVE. Without
            % this the test would pass on the broken code, because under the
            % shipped default (noccmode = 1) the fix is a strict no-op.
            testCase.verifyGreaterThan( ...
                max(abs(multi.tables.curve.(ptcol) - ...
                    single.tables.curve.(ptcol))), 1e-3, ...
                'n''_o = 4 must give a different curve than n''_o = 1.');

            % And it must equal the coefficient computed at n'_o = 4 directly --
            % the same call psyrat_relsummary makes for the summary tables. The
            % table stores values rounded to 4 dp, so compare like with like.
            d = rs.data.g(1).e(1);
            [~, pt4, ~] = psyrat_rel_trt('gcoeff', rs.gcoeff, ...
                'reltype', rs.reltype, ...
                'bp', d.sig_id.raw, 'bo', d.sig_occ.raw, 'bt', d.sig_trl.raw, ...
                'txp', d.sig_trlxid.raw, 'oxp', d.sig_occxid.raw, ...
                'txo', d.sig_trlxocc.raw, 'err', d.sig_err.raw, ...
                'obs', [1 10], 'nocc', 4, 'CI', 0.95);
            testCase.verifyEqual(multi.tables.curve.(ptcol), ...
                round(pt4(:), 4), 'AbsTol', 1e-12);

            % The other direction of the guard: an explicit n'_o = 1 must give
            % exactly the legacy fallback curve, so honoring the field cannot
            % change a single-occasion result.
            explicit_data = psyrat_data;
            explicit_data.relsummary.noccmode = 1;
            explicit_data.relsummary.nocc = 1;
            explicit = psyrat_report(explicit_data, 'ntrials', 10);
            testCase.verifyEqual(explicit.tables.curve.(ptcol), ...
                single.tables.curve.(ptcol), 'AbsTol', 1e-12);
        end

        function testWriteTableRejectsBadExtension(testCase)
            tbl = table([1;2], {'a';'b'}, 'VariableNames', {'x', 'y'});
            testCase.verifyError( ...
                @() psyrat_write_table(fullfile(tempdir, 'x.txt'), {}, tbl), ...
                'psyrat_write_table:ext');
        end

        function testProvenanceLinesDegradeGracefully(testCase)
            % An empty struct (no rel) still yields the standard lines as
            % "not recorded"/"not assessed" rather than erroring.
            lines = psyrat_provenance_lines(struct());
            testCase.verifyFalse(isempty(lines));
            joined = strjoin(lines, '|');
            testCase.verifySubstring(joined, 'Seed: not recorded');
            testCase.verifySubstring(joined, 'Converged: not assessed');
        end

        function testProvenanceFamilyLineGammaOnly(testCase)
            % The gamma family adds a "Likelihood family" provenance line; the
            % Gaussian default and a legacy file (no family field) add nothing,
            % so their headers stay byte-identical.
            gammaLines = psyrat_provenance_lines( ...
                struct('rel', struct('family', 'gamma')));
            gammaJoined = strjoin(gammaLines, '|');
            testCase.verifySubstring(gammaJoined, ...
                'Likelihood family: Gamma / scaled chi-square');

            gaussLines = psyrat_provenance_lines( ...
                struct('rel', struct('family', 'gaussian')));
            testCase.verifyFalse( ...
                any(contains(gaussLines, 'Likelihood family')), ...
                'Gaussian runs must not emit a family line.');

            legacyLines = psyrat_provenance_lines(struct());
            testCase.verifyFalse( ...
                any(contains(legacyLines, 'Likelihood family')), ...
                'Legacy runs (no family field) must not emit a family line.');

            % A string-typed empty family (e.g. a hand-built/legacy rel) must be
            % treated as absent, not emit a blank "Likelihood family:" line.
            emptyStrLines = psyrat_provenance_lines( ...
                struct('rel', struct('family', "")));
            testCase.verifyFalse( ...
                any(contains(emptyStrLines, 'Likelihood family')), ...
                'A string-typed empty family must not emit a family line.');
        end

        function testProvenanceDisclosesModularCutPosterior(testCase)
            % Pre-merge review, Increment A. Analyses 19/20 under gamma are a
            % modular CUT posterior (no copula feedback into the margins;
            % intervals are not full joint-posterior credible intervals). The
            % provenance must say so, quoting the stage's archived framework
            % token, and must aggregate the stage-2 diagnostic flag counts so
            % absence of trouble is affirmative rather than silence.
            cop = struct('sel', (1:4)', 'rho_e', 0.3*ones(4,1), ...
                'invalidity_flag', logical([1;0;0;0]), ...
                'calibration_flag', logical([0;1;1;0]), ...
                'labels', struct( ...
                'inference_framework', ...
                'modular_cut_two_stage_gamma_margins_gaussian_copula_v1', ...
                'copula_feedback_to_margins', false, ...
                'joint_posterior_status', ...
                'not_a_full_joint_bayesian_posterior'));
            rel = struct('analysis', 'ic_diff_dynrel_sserrvar_rescor', ...
                'family', 'gamma', 'out', struct('copula', {{cop}}));
            txt = strjoin(psyrat_provenance_lines(struct('rel', rel)), '|');
            testCase.verifySubstring(txt, 'modular CUT posterior');
            testCase.verifySubstring(txt, ...
                'modular_cut_two_stage_gamma_margins_gaussian_copula_v1');
            testCase.verifySubstring(txt, ...
                'Not a full joint Bayesian posterior');
            % The location-scale estimand statement every other gamma LS
            % design (12/13) emits must appear here too (2026-08-11 review).
            testCase.verifySubstring(txt, 'DIFFERENT ESTIMAND');
            testCase.verifySubstring(txt, ...
                '1 of 4 retained draws flagged numerically invalid');
            testCase.verifySubstring(txt, '2 flagged for PIT calibration');

            % The two-facet analysis string discloses identically.
            rel20 = rel;
            rel20.analysis = 'ic_diff_dynrel_sserrvar_trt_rescor';
            txt20 = strjoin(psyrat_provenance_lines(struct('rel', rel20)), '|');
            testCase.verifySubstring(txt20, 'modular CUT posterior');

            % A clean stage-2 reports "none flagged" affirmatively.
            cop2 = cop;
            cop2.invalidity_flag(:) = false;
            cop2.calibration_flag(:) = false;
            rel2 = rel;
            rel2.out.copula = {cop2};
            txt2 = strjoin(psyrat_provenance_lines(struct('rel', rel2)), '|');
            testCase.verifySubstring(txt2, ...
                'none of 4 retained draws flagged');

            % The affirmative all-clear must be EARNED: a stratum the counter
            % could not summarize (no sel, flags present but uncountable) must
            % downgrade the line, not vanish into "none flagged" (2026-08-11
            % review: absence of evidence was printed as evidence of absence).
            badcop = struct('rho_e', 0.3*ones(500,1), ...
                'invalidity_flag', true(500,1));
            relmix = rel2;
            relmix.out.copula = {cop2, badcop};
            txtmix = strjoin(psyrat_provenance_lines(struct('rel', relmix)), '|');
            testCase.verifyFalse(contains(txtmix, 'none of'), ...
                'A skipped stratum must forfeit the affirmative all-clear.');
            testCase.verifySubstring(txtmix, 'not fully recorded');

            % Same downgrade when a recorded stratum merely lacks a flag field.
            nofl = struct('sel', (1:4)', 'rho_e', 0.3*ones(4,1));
            relnf = rel2;
            relnf.out.copula = {nofl};
            txtnf = strjoin(psyrat_provenance_lines(struct('rel', relnf)), '|');
            testCase.verifyFalse(contains(txtnf, 'none of'), ...
                'Missing flag fields must not read as an affirmative zero.');
            testCase.verifySubstring(txtnf, 'not fully recorded');

            % Negative controls: a GAUSSIAN analysis-19 run and a gamma
            % non-copula run must not emit the disclosure - existing output
            % stays byte-identical outside the modular cut designs.
            gauss = psyrat_provenance_lines(struct('rel', struct( ...
                'analysis', 'ic_diff_dynrel_sserrvar_rescor', ...
                'family', 'gaussian')));
            testCase.verifyFalse( ...
                any(contains(gauss, 'CUT posterior')), ...
                'Gaussian analysis 19 must not emit the cut disclosure.');
            g12 = psyrat_provenance_lines(struct('rel', struct( ...
                'analysis', 'ic_diff_dynrel', 'family', 'gamma')));
            testCase.verifyFalse( ...
                any(contains(g12, 'CUT posterior')), ...
                'Gamma analysis 12 must not emit the cut disclosure.');

            % A gamma 19 result WITHOUT a stored copula stage (older file)
            % emits nothing rather than erroring.
            nostage = psyrat_provenance_lines(struct('rel', struct( ...
                'analysis', 'ic_diff_dynrel_sserrvar_rescor', ...
                'family', 'gamma')));
            testCase.verifyFalse(any(contains(nostage, 'CUT posterior')));

            % Labels absent from the stored stage (pre-release result): the
            % estimator class is still disclosed, from the analysis string.
            cop3 = rmfield(cop, 'labels');
            rel3 = rel;
            rel3.out.copula = {cop3};
            txt3 = strjoin(psyrat_provenance_lines(struct('rel', rel3)), '|');
            testCase.verifySubstring(txt3, 'modular CUT posterior');
        end

        function testProvenanceFlagsUnestimatedTrialMainEffect(testCase)
            % RC-15. When the one-facet trial main effect is absent or zero, the
            % absolute-error term contains no trial component, so dependability
            % and generalizability are numerically identical and a "Dependability"
            % heading is unearned. The number cannot show this, so the header must.
            % Emitted only in that case, so ordinary output stays byte-identical.
            estimated = psyrat_provenance_lines(struct('rel', struct( ...
                'analysis', 'ic', ...
                'out', struct('sig_trl', [0.4; 0.5; 0.6]))));
            testCase.verifyFalse( ...
                any(contains(estimated, 'Trial main effect')), ...
                'An estimated sigma_i must not add a caveat line.');

            % Explicit producer flag wins over inference.
            flagged = psyrat_provenance_lines(struct('rel', struct( ...
                'analysis', 'ic', ...
                'out', struct('sig_trl', [0.4; 0.5], ...
                'sig_trl_source', 'legacy_default'))));
            testCase.verifySubstring(strjoin(flagged, '|'), ...
                'Trial main effect (sigma_i): NOT estimated');

            % Backfilled zeros: the case the flag exists to make visible.
            zeroed = psyrat_provenance_lines(struct('rel', struct( ...
                'analysis', 'ic', 'out', struct('sig_trl', [0; 0; 0]))));
            testCase.verifySubstring(strjoin(zeroed, '|'), ...
                'Trial main effect (sigma_i): NOT estimated');

            % Absent entirely (pre-sigma_i saved file).
            missing = psyrat_provenance_lines(struct('rel', struct( ...
                'analysis', 'ic', 'out', struct('sig_u', [1; 2]))));
            testCase.verifySubstring(strjoin(missing, '|'), ...
                'Trial main effect (sigma_i): NOT estimated');

            % Scoped: designs that partition error differently would be
            % misdescribed by this line, so they must never emit it.
            trt = psyrat_provenance_lines(struct('rel', struct( ...
                'analysis', 'trt', 'out', struct('sig_u', [1; 2]))));
            testCase.verifyFalse(any(contains(trt, 'Trial main effect')), ...
                'Test-retest designs must not emit the one-facet sigma_i caveat.');

            % And a struct with no analysis at all must not error or emit.
            testCase.verifyFalse( ...
                any(contains(psyrat_provenance_lines(struct()), 'Trial main effect')));
        end

        function testProvenanceLabelsCopulaDispersionParameterization(testCase)
            % S11 (scope NARROWED 2026-07-31). The Gaussian-copula concurrent
            % gamma difference designs run one of two structurally different
            % models, selected by 'dispersion': per-person log-nu (1) has a
            % person-level funnel whose convergence is not established in this
            % repo, while global fixed nu (2) is the group-level default and was
            % checked live. Nothing user-facing said which one produced a table.
            % The label must name the PARAMETERIZATION, not the family, and must
            % not fire where 'dispersion' is meaningless -- rel.dispersion is
            % stamped on EVERY run (defaulting to 1), including Gaussian ones.
            key = 'Gamma copula dispersion parameterization';

            % (1) Global fixed nu: reported as checked, but explicitly NOT a
            % convergence guarantee (S11: "establishes tractability, not
            % robustness" -- one simulated dataset, one settings choice, and the
            % two-facet run still produced divergent transitions).
            globalLines = psyrat_provenance_lines(localCopulaData(2, ...
                {'is_global_nu', true}));
            globalJoined = strjoin(globalLines, '|');
            testCase.verifySubstring(globalJoined, ...
                [key ': GLOBAL fixed nu per event (dispersion = 2)']);
            testCase.verifySubstring(globalJoined, ...
                'not a convergence guarantee for this dataset');
            testCase.verifyFalse(contains(globalJoined, 'EXPERIMENTAL'), ...
                'The global-nu parameterization must not be labeled experimental.');

            % (2) Per-person log-nu: experimental, convergence not established
            % here, with the single-fit evidence reported as a single fit.
            perPerson = psyrat_provenance_lines(localCopulaData(1, {}));
            ppJoined = strjoin(perPerson, '|');
            testCase.verifySubstring(ppJoined, ...
                [key ': PER-PERSON log-nu (dispersion = 1) - EXPERIMENTAL']);
            testCase.verifySubstring(ppJoined, ...
                'has NOT been established in this repository');
            testCase.verifySubstring(ppJoined, ...
                'one in-repo subject-level fit reached R-hat of about 2.5');
            testCase.verifyFalse(contains(ppJoined, 'invalid'), ...
                'The label must not claim the estimates are invalid.');

            % (3) Analysis 8 (subject-level concurrent gamma difference) always
            % runs per-person nu, because global dispersion is a group-level
            % estimand. Say so, so the user knows there is no other option.
            case8 = strjoin(psyrat_provenance_lines(localCopulaData(1, ...
                {'is_diff_sserr', true})), '|');
            testCase.verifySubstring(case8, ...
                [key ': PER-PERSON log-nu (dispersion = 1) - EXPERIMENTAL']);
            testCase.verifySubstring(case8, ...
                'Subject-level reliability requires per-person nu');
            % ...and the group-level run must NOT carry that clause.
            testCase.verifyFalse( ...
                contains(ppJoined, 'Subject-level reliability requires'), ...
                'The group-level copula run must not claim it is subject-level.');

            % (4) Legacy saved file: no dispersion field at all. Best-effort
            % contract -- a "not recorded" line, never an error.
            legacyData = localCopulaData(1, {});
            legacyData.rel = rmfield(legacyData.rel, 'dispersion');
            legacyLines = psyrat_provenance_lines(legacyData);
            testCase.verifySubstring(strjoin(legacyLines, '|'), ...
                [key ': not recorded']);

            % (5) Scoping. A gamma design that is NOT the copula concurrent
            % difference family records a dispersion (the default 1) but must
            % emit no label: the option selects nothing there.
            plainSserr = localCopulaData(1, {'is_diff_sserr', true});
            plainSserr.rel.out.gamma_disp = rmfield( ...
                plainSserr.rel.out.gamma_disp, 'is_diff_copula');
            testCase.verifyFalse( ...
                any(contains(psyrat_provenance_lines(plainSserr), key)), ...
                'A non-copula gamma design must not gain a dispersion label.');

            % (6) ...and neither must a plain Gaussian run, even though the
            % estimator stamps dispersion = 1 on it.
            gauss = psyrat_provenance_lines(struct('rel', struct( ...
                'family', 'gaussian', 'analysis', 'ic', 'dispersion', 1)));
            testCase.verifyFalse(any(contains(gauss, key)), ...
                'Gaussian runs must not gain a dispersion label.');

            % (7) The degenerate struct must neither error nor emit.
            testCase.verifyFalse( ...
                any(contains(psyrat_provenance_lines(struct()), key)));
        end

        function testTrtCoefLabelNamesTheEstimand(testCase)
            % RC-12. A test-retest run reports ONE of three coefficients and the
            % selection is recorded nowhere else in an exported table, so the
            % label must name all three parts of the estimand: which facets are
            % generalized over (CE/CS/CES), the error type, and n'_o.
            mk = @(rt, gc) struct('relsummary', ...
                struct('reltype', rt, 'gcoeff', gc));

            testCase.verifyEqual(psyrat_trt_coeflabel(mk(1, 1)), ...
                ['Coefficient of Equivalence (CE), Dependability ' ...
                '(absolute error); Single-occasion TRT reliability, n''_o = 1']);
            testCase.verifyEqual(psyrat_trt_coeflabel(mk(2, 1)), ...
                ['Coefficient of Stability (CS), Dependability ' ...
                '(absolute error); Single-occasion TRT reliability, n''_o = 1']);
            testCase.verifyEqual(psyrat_trt_coeflabel(mk(3, 1)), ...
                ['Coefficient of Equivalence and Stability (CES), Dependability ' ...
                '(absolute error); Single-occasion TRT reliability, n''_o = 1']);

            % The error type is part of the estimand, not decoration: the same
            % components give different numbers under gcoeff 1 vs 2.
            testCase.verifySubstring(psyrat_trt_coeflabel(mk(2, 2)), ...
                'Generalizability (relative error)');

            % The occasion clause is delegated to psyrat_trt_occlabel, so it
            % tracks the multi-occasion selection rather than restating n'_o = 1.
            multi = mk(3, 1);
            multi.relsummary.noccmode = 2;
            multi.relsummary.nocc = 4;
            testCase.verifySubstring(psyrat_trt_coeflabel(multi), ...
                'Multi-occasion composite TRT reliability, n''_o = 4');

            % No reltype -> empty, NOT a guessed default. CE/CS/CES are three
            % different estimands and none is what a struct without the field
            % "must have been", so silence beats a mislabel. Unlike
            % psyrat_trt_occlabel, which can honestly default to n'_o = 1.
            testCase.verifyEmpty(psyrat_trt_coeflabel(struct()));
            testCase.verifyEmpty(psyrat_trt_coeflabel( ...
                struct('relsummary', struct('gcoeff', 1))));
            testCase.verifyEmpty(psyrat_trt_coeflabel( ...
                struct('relsummary', struct('reltype', 7, 'gcoeff', 1))));

            % A legacy struct with no gcoeff omits that clause rather than
            % guessing which error term produced the archived number.
            partial = psyrat_trt_coeflabel(struct('relsummary', ...
                struct('reltype', 2)));
            testCase.verifySubstring(partial, 'Coefficient of Stability (CS)');
            testCase.verifyFalse(contains(partial, 'error)'), ...
                'A struct with no gcoeff must not assert an error type.');
        end

        function testTrtCoefLabelReportsBothErrorTypesOnDifferenceScores(testCase)
            % A trt_diff export holds TWO coefficients under two independently
            % selected error types: psyrat_depoverallt prints the per-condition
            % rows at gcoeff and a "diff score" row at diffgcoeff in the same
            % table, and psyrat_relfigures reads the two from separate prefs.
            % Naming only gcoeff put the wrong error type beside half the table.
            split = struct('relsummary', struct( ...
                'reltype', 2, 'gcoeff', 1, 'diffgcoeff', 2));
            lbl = psyrat_trt_coeflabel(split);
            testCase.verifySubstring(lbl, 'Dependability (absolute error)');
            testCase.verifySubstring(lbl, ...
                'difference score: Generalizability (relative error)');

            % ...and the other way round, so the clause tracks the field rather
            % than always appending "Generalizability".
            flipped = struct('relsummary', struct( ...
                'reltype', 2, 'gcoeff', 2, 'diffgcoeff', 1));
            testCase.verifySubstring(psyrat_trt_coeflabel(flipped), ...
                'difference score: Dependability (absolute error)');

            % When the two agree - the shared default, and the common case -
            % the second clause is suppressed so output stays unchanged.
            agree = struct('relsummary', struct( ...
                'reltype', 2, 'gcoeff', 1, 'diffgcoeff', 1));
            testCase.verifyFalse( ...
                contains(psyrat_trt_coeflabel(agree), 'difference score'), ...
                'Matching selectors must not add a redundant clause.');

            % A non-difference design records no diffgcoeff at all.
            plain = struct('relsummary', struct('reltype', 2, 'gcoeff', 1));
            testCase.verifyFalse( ...
                contains(psyrat_trt_coeflabel(plain), 'difference score'));
        end

        function testProvenanceRecordsTrtCoefficientOnly(testCase)
            % RC-12. The coefficient line reaches every export through
            % psyrat_provenance_lines, and is gated so that designs recording no
            % reltype (one-facet, splits, dynamic) stay byte-identical.
            trt = psyrat_provenance_lines(struct( ...
                'rel', struct('analysis', 'trt'), ...
                'relsummary', struct('reltype', 2, 'gcoeff', 1)));
            testCase.verifySubstring(strjoin(trt, '|'), ...
                ['Coefficient: Coefficient of Stability (CS), Dependability ' ...
                '(absolute error); Single-occasion TRT reliability, n''_o = 1']);

            % One-facet summary fixture: no reltype, so no line at all.
            ic = psyrat_provenance_lines( ...
                PsyRATTestDataFactory.makeICSummaryData());
            testCase.verifyFalse(any(contains(ic, 'Coefficient:')), ...
                'One-facet output must not gain a coefficient line.');

            % And the degenerate structs must neither error nor emit.
            testCase.verifyFalse( ...
                any(contains(psyrat_provenance_lines(struct()), 'Coefficient:')));
        end

        function testReportSurfacesGammaFamilyEndToEnd(testCase)
            % End-to-end: a gamma REL propagates the family label through
            % psyrat_report into report.provenance and into the written file
            % header (no live fit; reuses the one-facet summary fixture).
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            psyrat_data.rel.family = 'gamma';
            tdir = localTempDir(testCase);

            report = psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv');

            joined = strjoin(report.provenance, '|');
            testCase.verifySubstring(joined, ...
                'Likelihood family: Gamma / scaled chi-square');

            % The label rides the provenance header of the exported table.
            txt = fileread(report.files{1});
            testCase.verifySubstring(txt, ...
                'Likelihood family: Gamma / scaled chi-square');
        end

        function testReportReturnsSplitsTables(testCase)
            % Nonparallel data-splits (ic_splits): reads rel directly (no
            % relsummary) and is render-free, like DoD.
            psyrat_data = PsyRATTestDataFactory.makeSplitsSummaryData();

            report = psyrat_report(psyrat_data);

            testCase.verifyEqual(report.analysis, 'splits');
            testCase.verifyTrue(isfield(report.tables, 'coefficients'));
            testCase.verifyTrue(isfield(report.tables, 'varcomp'));
            testCase.verifyTrue(istable(report.tables.coefficients));
            testCase.verifyTrue(istable(report.tables.varcomp));
            % One stratum, one coefficient for the single-occasion design.
            testCase.verifyEqual(height(report.tables.coefficients), 1);
            % The flattened coefficient row reflects psyrat_splits_summary
            % (the compute layer) -- the Dependability point estimate matches.
            summ = psyrat_splits_summary(psyrat_data.rel, 'CI', .95, 'obs', [], 'nocc', []);
            testCase.verifyEqual(report.tables.coefficients.Dependability(1), ...
                round(summ.strata(1).coef(1).D(2), 3));
            % The generalizability lower-bound caveat is surfaced.
            testCase.verifyFalse(isempty(report.gbias_note));
            testCase.verifySubstring(report.gbias_note, 'downward-biased');
            % Provenance lines are present.
            testCase.verifyFalse(isempty(report.provenance));
        end

        function testReportSplitsOpensNoFigures(testCase)
            psyrat_data = PsyRATTestDataFactory.makeSplitsSummaryData();
            nBefore = numel(findall(0, 'Type', 'figure'));
            psyrat_report(psyrat_data);
            nAfter = numel(findall(0, 'Type', 'figure'));
            testCase.verifyEqual(nAfter, nBefore, ...
                'Data-splits report must not open any figures.');
        end

        function testReportSplitsWritesFilesWithProvenance(testCase)
            psyrat_data = PsyRATTestDataFactory.makeSplitsSummaryData();
            tdir = localTempDir(testCase);

            report = psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv');

            % One file per table (coefficients + varcomp).
            testCase.verifyEqual(numel(report.files), 2);
            for k = 1:numel(report.files)
                testCase.verifyEqual(exist(report.files{k}, 'file'), 2);
            end
            % The header carries provenance, the toolbox identity, and the
            % generalizability lower-bound caveat above the table.
            txt = fileread(report.files{1});
            testCase.verifySubstring(txt, 'Seed:');
            testCase.verifySubstring(txt, 'PsyRAT Toolbox');
            testCase.verifySubstring(txt, 'downward-biased');
        end

        function testReportSplitsTRT(testCase)
            % Multi-occasion data-splits (trt_splits): three test-retest
            % reltypes (CE/CS/CES) per stratum.
            psyrat_data = PsyRATTestDataFactory.makeSplitsTRTSummaryData();

            report = psyrat_report(psyrat_data);

            testCase.verifyEqual(report.analysis, 'splits');
            testCase.verifyEqual(height(report.tables.coefficients), 3);
            testCase.verifyFalse(isempty(report.tables.varcomp));
        end

        function testReportReturnsDynrelTables(testCase)
            % Dynamic reliability (ic_diff_dynrel): reads rel directly (no
            % relsummary) and is render-free. Surface + variance components; no
            % per-participant table for this group-level variant.
            psyrat_data = PsyRATTestDataFactory.makeDynrelDiffSummaryData();

            report = psyrat_report(psyrat_data);

            testCase.verifyEqual(report.analysis, 'dynrel');
            testCase.verifyTrue(isfield(report.tables, 'surface'));
            testCase.verifyTrue(isfield(report.tables, 'varcomp'));
            testCase.verifyTrue(istable(report.tables.surface));
            testCase.verifyTrue(istable(report.tables.varcomp));
            % One stratum, default ngrid = 25 -> 25 surface rows.
            testCase.verifyEqual(height(report.tables.surface), 25);
            % The flattened surface reflects psyrat_dynrel_summary (compute layer).
            summ = psyrat_dynrel_summary(psyrat_data.rel, 'CI', .95, ...
                'marginal', 0, 'reltype', 3, 'nocc', [], 'ngrid', 25);
            testCase.verifyEqual(report.tables.surface.Generalizability(1), ...
                round(summ.strata(1).G.pt(1), 4));
            testCase.verifyEqual(report.tables.surface.Dependability(end), ...
                round(summ.strata(1).D.pt(end), 4));
            % The raw surface struct is returned for programmatic use.
            testCase.verifyTrue(isfield(report, 'strata'));
            % Estimand note is surfaced (typical person for the default marginal).
            testCase.verifySubstring(report.estimand_note, 'typical person');
            testCase.verifyFalse(isempty(report.provenance));
        end

        function testReportDynrelOpensNoFigures(testCase)
            psyrat_data = PsyRATTestDataFactory.makeDynrelDiffSummaryData();
            nBefore = numel(findall(0, 'Type', 'figure'));
            psyrat_report(psyrat_data);
            nAfter = numel(findall(0, 'Type', 'figure'));
            testCase.verifyEqual(nAfter, nBefore, ...
                'Dynamic-reliability report must not open any figures.');
        end

        function testReportDynrelSubjectLevelHasSsrel(testCase)
            % Subject-level variant (ic_diff_dynrel_sserrvar) produces a
            % per-participant table in addition to surface + variance components.
            psyrat_data = PsyRATTestDataFactory.makeDynrelSubjSummaryData();

            report = psyrat_report(psyrat_data);

            testCase.verifyEqual(report.analysis, 'dynrel');
            testCase.verifyTrue(isfield(report.tables, 'ssrel'));
            testCase.verifyTrue(istable(report.tables.ssrel));
            % 3 subjects in the fixture -> 3 per-participant rows.
            testCase.verifyEqual(height(report.tables.ssrel), 3);
            % Parity with the compute layer's per-participant table.
            summ = psyrat_dynrel_summary(psyrat_data.rel, 'CI', .95, 'ngrid', 25);
            testCase.verifyEqual(height(report.tables.ssrel), ...
                height(summ.strata(1).ssrel_table));
        end

        function testReportDynrelNondiffSubjectLevelHasSsrel(testCase)
            % NON-difference subject-level variants (ic_dynrel_sserrvar, analysis
            % 26; ic_dynrel_sserrvar_trt, analysis 27) also produce a
            % per-participant table alongside the surface + variance components.
            factories = {@() PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData(), ...
                @() PsyRATTestDataFactory.makeDynrelSubjNondiffTrtSummaryData()};
            for k = 1:numel(factories)
                psyrat_data = factories{k}();
                report = psyrat_report(psyrat_data);
                testCase.verifyEqual(report.analysis, 'dynrel');
                testCase.verifyTrue(isfield(report.tables, 'ssrel'));
                testCase.verifyTrue(istable(report.tables.ssrel));
                testCase.verifyEqual(height(report.tables.ssrel), 3);
                summ = psyrat_dynrel_summary(psyrat_data.rel, 'CI', .95, 'ngrid', 25);
                testCase.verifyEqual(height(report.tables.ssrel), ...
                    height(summ.strata(1).ssrel_table));
            end
        end

        function testReportDynrelNonSubjectOmitsSsrel(testCase)
            % Group-level variants have no per-participant table; the field is
            % omitted rather than written as an empty table.
            psyrat_data = PsyRATTestDataFactory.makeDynrelDiffSummaryData();
            report = psyrat_report(psyrat_data);
            testCase.verifyFalse(isfield(report.tables, 'ssrel'));
        end

        function testReportDynrelTrtReltypePassthrough(testCase)
            % Two-facet variant (ic_dynrel_trt) consumes reltype/nocc: a
            % different coefficient yields a different surface, and the estimand
            % note records the coefficient + occasion n'.
            psyrat_data = PsyRATTestDataFactory.makeDynrelTrtSummaryData();

            r1 = psyrat_report(psyrat_data, 'reltype', 1);
            r3 = psyrat_report(psyrat_data, 'reltype', 3);

            testCase.verifyFalse(isequal(r1.tables.surface.Generalizability, ...
                r3.tables.surface.Generalizability), ...
                'reltype must change the two-facet surface.');
            testCase.verifySubstring(r3.estimand_note, 'coefficient:');
            testCase.verifySubstring(r3.estimand_note, 'occasion');
        end

        function testReportDynrel2DSurface(testCase)
            % Two-dimension surface flattens to ngrid^2 rows with a z2 column.
            psyrat_data = PsyRATTestDataFactory.makeDynrelTrt2DSummaryData();
            report = psyrat_report(psyrat_data, 'ngrid', 10);
            testCase.verifyEqual(height(report.tables.surface), 100);
            testCase.verifyTrue(ismember('z2', ...
                report.tables.surface.Properties.VariableNames));
        end

        function testReportDynrelWritesFilesWithProvenance(testCase)
            psyrat_data = PsyRATTestDataFactory.makeDynrelDiffSummaryData();
            tdir = localTempDir(testCase);

            report = psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv');

            % One file per table (surface + varcomp), all present on disk.
            testCase.verifyEqual(numel(report.files), numel(fieldnames(report.tables)));
            for k = 1:numel(report.files)
                testCase.verifyEqual(exist(report.files{k}, 'file'), 2);
            end
            % The header carries provenance and the estimand note above the table.
            txt = fileread(report.files{1});
            testCase.verifySubstring(txt, 'Seed:');
            testCase.verifySubstring(txt, 'PsyRAT Toolbox');
            testCase.verifySubstring(txt, 'Estimand');
        end

        function testReportRefusesOverwriteWhenGeneric(testCase)
            % With no proc.savename the tables are named generically
            % ('psyrat_report_*'); a second run into the same outdir must refuse
            % to silently clobber the first.
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData(); % no savename
            tdir = localTempDir(testCase);
            psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv'); % first OK
            testCase.verifyError( ...
                @() psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv'), ...
                'psyrat_report:overwrite');
        end

        function testReportOverwritesWhenNamed(testCase)
            % With proc.savename set, re-running into the same outdir overwrites
            % on purpose (deliberate re-runs) rather than erroring.
            psyrat_data = PsyRATTestDataFactory.makeSplitsSummaryData(); % has savename
            tdir = localTempDir(testCase);
            psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv');
            report = psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv');
            testCase.verifyEqual(numel(report.files), 2);
        end

        function testReportDatasetHeaderFallsBackToRaw(testCase)
            % Scripted/in-memory runs that never set rel.filename still record a
            % dataset in the export header via the raw.filename fallback.
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            psyrat_data.raw = struct('filename', 'inmemory_source.csv');
            tdir = localTempDir(testCase);
            report = psyrat_report(psyrat_data, 'outdir', tdir, 'format', '.csv');
            txt = fileread(report.files{1});
            testCase.verifySubstring(txt, 'Dataset: inmemory_source.csv');
        end

        function testReportSurfaceRecordsOccasionForTwoFacet(testCase)
            % Two-facet dynrel surfaces carry the occasion n' as a column (the
            % G(z)/D(z) values depend on it); one-facet variants omit it.
            r2 = psyrat_report(PsyRATTestDataFactory.makeDynrelTrtSummaryData());
            testCase.verifyTrue(ismember('nOccasions', ...
                r2.tables.surface.Properties.VariableNames));
            r1 = psyrat_report(PsyRATTestDataFactory.makeDynrelDiffSummaryData());
            testCase.verifyFalse(ismember('nOccasions', ...
                r1.tables.surface.Properties.VariableNames));
        end

        function testReportUnsupportedDynrelVariantErrorsClearly(testCase)
            % A dynrel-family analysis string not in the dispatch list errors
            % with the clear 'unsupported' id rather than the misleading
            % 'norelsummary' one (dynrel never populates relsummary).
            psyrat_data = PsyRATTestDataFactory.makeDynrelDiffSummaryData();
            psyrat_data.rel.analysis = 'ic_dynrel_unknownvariant';
            testCase.verifyError(@() psyrat_report(psyrat_data), ...
                'psyrat_report:unsupported');
        end

        function testReportDefaultsCiToTheStoredWidth(testCase)
            % G41(i): a plain psyrat_report call must compute its intervals
            % at the width the summary was built at, so one export carries
            % ONE width (the TestCiCoverageConsistency one-width rule). The
            % stored width is doctored to .50 so "default follows stored" is
            % distinguishable from the old hardcoded .95 -- and far enough
            % from it that the fixture's FOUR drifting draws move the
            % quantiles (at n = 4 the extreme quantiles clamp to min/max, so
            % adjacent widths render identically; the anti-vacuity assert
            % below caught exactly that on a first draft, and .50
            % interpolates interior positions instead).
            pd = PsyRATTestDataFactory.makeICSummaryData();
            pd.relsummary.ciperc = 0.50;

            repDefault = psyrat_report(pd);
            repStored  = psyrat_report(pd, 'CI', 0.50);
            repWide    = psyrat_report(pd, 'CI', 0.95);

            % anti-vacuity: the width genuinely moves the curve intervals on
            % this fixture, so the equality below cannot hold trivially
            testCase.assertFalse( ...
                isequal(repStored.tables.curve, repWide.tables.curve), ...
                'fixture draws must make the interval width observable');
            testCase.verifyEqual(repDefault.tables.curve, ...
                repStored.tables.curve, ...
                'the default CI must be the stored relsummary.ciperc');
            testCase.verifyEmpty(repDefault.ci_note, ...
                'one width in the export: no disclosure note');
            testCase.verifyEmpty(repStored.ci_note, ...
                'an explicit CI equal to the stored width needs no note');
        end

        function testReportDisclosesAnExplicitMismatchedCi(testCase)
            % G41 disclosure: an explicit 'CI' differing from the stored
            % width still puts two widths in one export (the stored-width
            % tables cannot be recomputed), so the header must say which
            % intervals carry which width -- and the note must reach the
            % written file, not just the return struct.
            pd = PsyRATTestDataFactory.makeICSummaryData();   % stored .95
            rep = psyrat_report(pd, 'CI', 0.99);
            testCase.verifyTrue(contains(rep.ci_note, '99%'), ...
                'the note must name the report''s computed width');
            testCase.verifyTrue(contains(rep.ci_note, '95%'), ...
                'the note must name the stored width');

            tdir = localTempDir(testCase);
            rep = psyrat_report(pd, 'CI', 0.99, 'outdir', tdir, ...
                'format', '.csv');
            txt = fileread(rep.files{1});
            %pin the note's substance in the written file, not just a
            %comma-free prefix: both widths must survive the CSV write
            testCase.verifySubstring(txt, ...
                'reliability-vs-n curve intervals in this report use 99%');
            testCase.verifySubstring(txt, ...
                '95% width recorded above');
        end

        function testMismatchNoteIsHonestWhenCiGovernsNothing(testCase)
            %G41 disclosure honesty (adversarial wave, 2026-08-29): on the
            %sserr route the report's 'CI' option reaches no table at all
            %(.variance and .sscoeffs are built from summary-time and
            %precomputed values), so a mismatched explicit 'CI' must NOT
            %claim curve intervals exist -- it must say the option did not
            %apply. The .curve-family sibling above pins the other wording.
            pd = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
            [pd, ~] = psyrat_relsummary( ...
                'psyrat_data', pd, ...
                'analysis', 'sing_sserr', ...
                'depcutoff', 0.1, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'CI', 0.95);
            testCase.assertFalse(isfield(psyrat_report(pd).tables, 'curve'), ...
                'weld: the sserr report must carry no curve table');

            rep = psyrat_report(pd, 'CI', 0.50);
            testCase.verifyTrue(contains(rep.ci_note, 'does not affect'), ...
                'a CI that governs nothing must be disclosed as inert');
            testCase.verifyFalse(contains(rep.ci_note, 'curve'), ...
                'the inert-CI note must not assert curve intervals exist');
        end

        function testMismatchGateIsNanSafe(testCase)
            %G41 (adversarial wave): an explicit 'CI',NaN used to defeat the
            %disclosure gate (abs(NaN - stored) > tol is false), so the
            %header asserted the stored width with no note. The gate is now
            %spelled ~(<= tol), so NaN counts as a mismatch. Driven through
            %the sserr family, where the option governs no computation, so
            %the NaN cannot hard-error in a quantile call before the note is
            %reached -- on the curve families a NaN width fails loudly
            %downstream, which is also acceptable (the silent-lie state is
            %what the gate fix removes).
            pd = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
            [pd, ~] = psyrat_relsummary( ...
                'psyrat_data', pd, ...
                'analysis', 'sing_sserr', ...
                'depcutoff', 0.1, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'CI', 0.95);
            rep = psyrat_report(pd, 'CI', NaN);
            testCase.verifyNotEmpty(rep.ci_note, ...
                'an explicit NaN width must be disclosed, not swallowed');
        end

        function testProvenanceRecordsTheStoredCiWidth(testCase)
            % G41(ii): the shared provenance header records the width the
            % summary was built at, so a saved table's intervals can be
            % interpreted on their own. Emit-only-when-applicable: a struct
            % with no built summary contributes no line, so those exports
            % stay byte-identical.
            pd = PsyRATTestDataFactory.makeICSummaryData();
            joined = strjoin(psyrat_provenance_lines(pd), '|');
            testCase.verifySubstring(joined, 'Credible interval: 95%');

            joinedNone = strjoin(psyrat_provenance_lines(struct()), '|');
            testCase.verifyFalse(contains(joinedNone, 'Credible interval'), ...
                'no summary, no width line (emit-only-when-applicable)');

            %the interesting negative branches (adversarial wave): a summary
            %that RECORDS a degenerate width -- NaN, empty, non-scalar --
            %must contribute no line either; relsummary's own validator
            %accepts these, so the gate is the only guard.
            for bad = {NaN, [], [0.9 0.95]}
                pdBad = pd;
                pdBad.relsummary.ciperc = bad{1};
                testCase.verifyFalse(contains( ...
                    strjoin(psyrat_provenance_lines(pdBad), '|'), ...
                    'Credible interval'), ...
                    'a degenerate recorded width must emit no line');
            end
        end

    end

    methods (Access = private)
        function tdir = localTempDir(testCase)
            tdir = tempname;
            mkdir(tdir);
            testCase.addTeardown(@() rmdir(tdir, 's'));
        end
    end

end

function psyrat_data = localCopulaData(dispersion, extraFlags)
%Minimal Gaussian-copula concurrent gamma difference result: the dispersion
%provenance struct psyrat_gamma_extract_diff_copula* writes into
%rel.out.gamma_disp (s_v / disp_factor draws plus the is_diff_copula flag),
%carried on a gamma rel with the RESOLVED dispersion the estimator stamps as
%REL.dispersion. extraFlags is a name/value cell appended to gamma_disp, e.g.
%{'is_global_nu', true} for the dispersion=2 extractor's structural zeros or
%{'is_diff_sserr', true} for the subject-level (analysis 8) inheritor.
gd = struct('s_v', [0.30; 0.35], 'rho_pv', [-0.10; -0.12], ...
    'disp_factor', [1.05; 1.06], 'is_diff_copula', true);
for k = 1:2:numel(extraFlags)
    gd.(extraFlags{k}) = extraFlags{k+1};
end
if isfield(gd, 'is_global_nu') && gd.is_global_nu
    % The global-nu extractor supplies structural zeros, not estimates.
    gd.s_v = zeros(2,2);
    gd.disp_factor = ones(2,2);
    gd = rmfield(gd, 'rho_pv');
end
psyrat_data = struct('rel', struct( ...
    'family', 'gamma', 'analysis', 'ic_diff', 'dispersion', dispersion, ...
    'out', struct('gamma_disp', gd)));
end

function psyrat_data = localMakeSserrData()
%Build a genuine subject-level (sserr) result by running psyrat_relsummary on the
%difference sserr fixture, so the per-event ssrel_table and the location-scale
%data are real (matching the pipeline). analysis becomes 'ic_diff_sserrvar'. This
%mirrors TestSummaryAndPresentationFunctions' subject-level-diff setup.
psyrat_data = PsyRATTestDataFactory.makeICDiffSSErrRelData();
[psyrat_data, ~] = psyrat_relsummary( ...
    'psyrat_data', psyrat_data, ...
    'analysis', 'sing_diff_sserr', ...
    'depcutoff', 0.5, ...
    'meascutoff', 2, ...
    'depcentmeas', 1, ...
    'diffgcoeff', 1, ...
    'CI', 0.95);
end

function psyrat_data = localMakeTrtSserrData()
%Build a genuine subject-level crossed test-retest (trt_sserrvar) result by
%running psyrat_relsummary on the crossed test-retest sserr fixture, so the
%per-participant ssrel_table and the group-level variance data are real (matching
%the pipeline). analysis stays 'trt_sserrvar'.
psyrat_data = PsyRATTestDataFactory.makeTRTSSErrRelData();
[psyrat_data, ~] = psyrat_relsummary( ...
    'psyrat_data', psyrat_data, ...
    'analysis', 'trt_sserr', ...
    'gcoeff', 1, ...
    'reltype', 3, ...
    'depcutoff', 0.5, ...
    'meascutoff', 2, ...
    'CI', 0.95);
end
