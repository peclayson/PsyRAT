classdef TestExportCallbacks < PsyRATTestBase
    % Regression guard for the GUI "Save Table" export callbacks.
    %
    % The Save buttons on PsyRAT output figures fire local callbacks that call
    % uiputfile + writecell/writetable/fprintf. The rest of the suite exercises
    % psyrat_relsummary and psyrat_relfigures but NOT these save callbacks, so two
    % real bugs reached the 2026-06 GUI review undetected:
    %
    %   * F-CM-M4trtcut: the test-retest trial-cutoff save callback read
    %     relsummary.depcutoff, which the plain-'trt' relsummary branch never sets
    %     (it sets relcutoff), so "Save Table" threw "Unrecognized field name
    %     depcutoff" for both CSV and XLSX.
    %   * F-CM-M4xlsxhdr: the XLSX variance export wrote the compact table variable
    %     names (Between_StdDev) instead of the on-screen labels (Between Std Dev).
    %
    % Strategy (no CmdStan, no untracked fixtures): build synthetic post-summary
    % psyrat_data via PsyRATTestDataFactory, call the real table builder with
    % gui=1 to create the output figure + Save button, shadow uiputfile to a temp
    % path, fire the Save callback for both .csv and .xlsx, and assert the file is
    % written without error and its header matches the on-screen labels.
    %
    % Coverage note: full CSV/XLSX parity across every design (criterion,
    % subject-level, DoD, the all-component difference tables) is exercised by the
    % manual export sweep gui_review_artifacts/current_main_export_sweep.m, which
    % loads precomputed .psyrat fixtures. This committed test pins the two
    % regressions that the headless path can build deterministically from the
    % shared data factory.

    properties
        ShadowDir
    end

    methods (TestMethodSetup)
        function setupSaveShadow(testCase)
            % Put a temp uiputfile shadow on the path so the Save callbacks
            % auto-answer their dialog from the PSYRAT_FAKE_SAVE global instead
            % of blocking on a modal file picker.
            testCase.ShadowDir = tempname;
            mkdir(testCase.ShadowDir);
            fid = fopen(fullfile(testCase.ShadowDir,'uiputfile.m'),'w');
            fprintf(fid,'%s\n', ...
                'function [filename,pathname,filterindex] = uiputfile(varargin)', ...
                '% TEST SHADOW: auto-answers the save dialog from a global.', ...
                'global PSYRAT_FAKE_SAVE', ...
                '[p,n,e] = fileparts(PSYRAT_FAKE_SAVE);', ...
                'filename = [n e]; pathname = [p filesep]; filterindex = 1;', ...
                'end');
            fclose(fid);
            addpath(testCase.ShadowDir,'-begin');
            rehash;
            testCase.addTeardown(@() teardownSaveShadow(testCase));
        end
    end

    methods (Test)
        function testTrtTrialCutoffSaveWritesBothFormats(testCase)
            % Guards F-CM-M4trtcut: firing the test-retest trial-cutoff Save
            % button must not error, and the written header must carry the
            % reliability cutoff (proving relsummary.relcutoff is read).
            pd = localTrtData();
            fig = localBuildOutputFig(testCase, ...
                @() psyrat_trt_relcutofft('psyrat_data', pd, 'gui', 1));

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);   % throws if the bug returns
            testCase.verifyTrue(isfile(csvfile), 'TRT trial-cutoff CSV not written');
            csvtxt = fileread(csvfile);
            testCase.verifyTrue(contains(csvtxt,'Dependability Cutoff: 0.8000'), ...
                'CSV header should carry the reliability cutoff value');

            xlsxfile = [tempname '.xlsx'];
            localFireSave(testCase, fig, xlsxfile);
            testCase.verifyTrue(isfile(xlsxfile), 'TRT trial-cutoff XLSX not written');
            raw = readcell(xlsxfile);                % errors if not a real workbook
            testCase.verifyGreaterThan(size(raw,1), 1, 'XLSX workbook looks empty');

            localCleanFiles({csvfile, xlsxfile});
        end

        function testVarianceSaveHeadersMatchOnScreenLabels(testCase)
            % Guards F-CM-M4xlsxhdr: the one-facet variance Save must write the
            % spaced on-screen labels (Between Std Dev / Trial Std Dev) in BOTH
            % CSV and XLSX, not the compact table variable names (Between_StdDev).
            pd = localIcData();
            fig = localBuildOutputFig(testCase, ...
                @() psyrat_variancet('psyrat_data', pd, 'gui', 1));

            xlsxfile = [tempname '.xlsx'];
            localFireSave(testCase, fig, xlsxfile);
            raw = readcell(xlsxfile);
            hdr = localHeaderRow(raw, 'Label');
            testCase.verifyTrue(localHasCell(hdr,'Between Std Dev'), ...
                'XLSX header should read "Between Std Dev"');
            testCase.verifyTrue(localHasCell(hdr,'Trial Std Dev'), ...
                'XLSX header should read "Trial Std Dev"');
            testCase.verifyFalse(localHasCell(hdr,'Between_StdDev'), ...
                'XLSX header must not write the compact variable name');

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            csvtxt = fileread(csvfile);
            testCase.verifyTrue(contains(csvtxt,'Between Std Dev'), ...
                'CSV header should read "Between Std Dev"');
            testCase.verifyTrue(contains(csvtxt,'Trial Std Dev'), ...
                'CSV header should read "Trial Std Dev"');

            localCleanFiles({xlsxfile, csvfile});
        end

        function testDiffRowCsvExportKeepsDashCellsIntact(testCase)
            % G47: the overall-table CSV loop pushed the diff row's '---'
            % cells through %d, which printed each character's numeric code
            % (45 for '-') and recycled the formatspec, so the one on-screen
            % diff row landed as three garbled lines in every saved
            % ic_diff/trt_diff CSV. Fire the real Save callback on an ic_diff
            % fixture: the diff row must come out as ONE line matching its
            % on-screen cells, and an all-numeric event row must stay
            % byte-identical to the pre-fix rendering (positive control that
            % the numeric path did not move).
            pd = localIcDiffData();
            fig = localBuildOutputFig(testCase, ...
                @() psyrat_depoverallt('psyrat_data', pd, 'gui', 1));

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile, 'Save Table');
            lines = localCsvLines(csvfile);

            hdr = find(startsWith(lines, 'Label,'), 1);
            testCase.assertNotEmpty(hdr, 'table header row not found');
            body = lines(hdr+1:end);
            testCase.assertEqual(numel(body), 3, sprintf( ...
                ['expected two event rows plus one diff row after the ', ...
                'header; a garbled diff row fragments into extra lines:\n%s'], ...
                strjoin(body, newline)));

            % the diff row, cell for cell: label, goodn, dashed badn, the
            % draw-based coefficient string, then five dashed trial cells
            ds = pd.relsummary.group(1).diffscore;
            expected_diff = sprintf('diff score,%d,---,%s,---,---,---,---,---', ...
                pd.relsummary.group(1).event(2).goodn, ...
                sprintf(' %0.2f CI [%0.2f %0.2f]', ds.pt, ds.ll, ds.ul));
            testCase.verifyEqual(body{3}, expected_diff, ...
                'diff row must survive as one intact line with its dash cells');

            % positive control: an all-numeric event row rendered through the
            % OLD formatspec (sans trailing newline) -- the numeric path must
            % be byte-identical before and after the all-%s rewrite
            ev = pd.relsummary.group(1).event(1);
            expected_row1 = sprintf('%s,%d,%d,%s,%0.2f,%d,%0.2f,%d,%d', ...
                'loss', ev.goodn, numel(pd.relsummary.group(1).badids), ...
                sprintf(' %0.2f CI [%0.2f %0.2f]', ev.dep.m, ev.dep.ll, ev.dep.ul), ...
                ev.trlinfo.mean, ev.trlinfo.med, ev.trlinfo.std, ...
                ev.trlinfo.min, ev.trlinfo.max);
            testCase.verifyEqual(body{1}, expected_row1, ...
                'all-numeric event row must stay byte-identical to the pre-fix rendering');

            % second control: the gain row's median is FRACTIONAL (7.5), so
            % %d renders it in scientific notation -- both old and new code
            % do, and this pins that quirky-but-shipping rendering
            evB = pd.relsummary.group(1).event(2);
            testCase.assertNotEqual(evB.trlinfo.med, round(evB.trlinfo.med), ...
                'fixture precondition broken: the gain median must be fractional');
            expected_row2 = sprintf('%s,%d,%d,%s,%0.2f,%g,%0.2f,%d,%d', ...
                'gain', evB.goodn, numel(pd.relsummary.group(1).badids), ...
                sprintf(' %0.2f CI [%0.2f %0.2f]', evB.dep.m, evB.dep.ll, evB.dep.ul), ...
                evB.trlinfo.mean, evB.trlinfo.med, evB.trlinfo.std, ...
                evB.trlinfo.min, evB.trlinfo.max);
            testCase.verifyEqual(body{2}, expected_row2, ...
                ['fractional-median row: the median prints as a decimal (pre-beta ', ...
                'median fix, %g); every other cell is byte-identical to the pre-G47 rendering']);

            % G44-shaped row through the SAME writer: a zeroed sserr group
            % dashes the coefficient and trial cells while keeping both
            % counts numeric (the opposite dash pattern from the diff row).
            % Injected via the Save button's own callback args, since the
            % headless factory cannot cheaply stage a zeroed sserr group
            % with a gui figure.
            btns = findall(fig, 'Style', 'pushbutton');
            cb = [];
            for b = btns(:)'
                s = b.String; if iscell(s), s = strjoin(s, ' '); end
                if contains(string(s), 'Save Table'); cb = b.Callback; break; end
            end
            testCase.assertFalse(isempty(cb), 'Save Table button not found');
            tbl = cb{3};
            for col = [4 5 6 7 8 9]
                tbl.(col){1} = '---';
            end
            global PSYRAT_FAKE_SAVE
            g44file = [tempname '.csv'];
            PSYRAT_FAKE_SAVE = g44file;
            feval(cb{1}, [], [], cb{2}, tbl);
            glines = localCsvLines(g44file);
            ghdr = find(startsWith(glines, 'Label,'), 1);
            expected_g44 = sprintf('%s,%d,%d,---,---,---,---,---,---', ...
                'loss', ev.goodn, numel(pd.relsummary.group(1).badids));
            testCase.verifyEqual(glines{ghdr+1}, expected_g44, ...
                'a G44-shaped dashed row must survive as one intact CSV line');

            localCleanFiles({csvfile, g44file});
        end

        function testOverallCsvFractionalMedianPrintsAsDecimal(testCase)
            %P4 observation (2026-09-04), fixed in the pre-beta batch: the
            %Save Table CSV printed a half-integer median trial count as
            %2.450000e+01, because column 6's format was '%d', which sprintf
            %renders as %e for a non-integer. The median of an even-count
            %trial vector is routinely x.5. An integer median must keep its
            %integer rendering (byte-identical to before).
            pd = localIcDiffData();
            fig = localBuildOutputFig(testCase, ...
                @() psyrat_depoverallt('psyrat_data', pd, 'gui', 1));
            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile, 'Save Table');
            lines = localCsvLines(csvfile);
            hdr = find(startsWith(lines, 'Label,'), 1);
            body = lines(hdr+1:end);

            evB = pd.relsummary.group(1).event(2);
            testCase.assertNotEqual(evB.trlinfo.med, round(evB.trlinfo.med), ...
                'fixture precondition broken: the gain median must be fractional');
            cellsB = strsplit(body{2}, ',');
            testCase.verifyEqual(cellsB{6}, sprintf('%g', evB.trlinfo.med), ...
                'a fractional median must print as a plain decimal in the CSV');
            testCase.verifyFalse(contains(body{2}, 'e+'), ...
                'no scientific notation may appear in an exported CSV row');

            evA = pd.relsummary.group(1).event(1);
            testCase.assertEqual(evA.trlinfo.med, round(evA.trlinfo.med), ...
                'fixture precondition broken: the loss median must be an integer');
            cellsA = strsplit(body{1}, ',');
            testCase.verifyEqual(cellsA{6}, sprintf('%d', evA.trlinfo.med), ...
                'an integer median keeps its integer rendering');

            localCleanFiles({csvfile});
        end

        function testTrtOverallCsvFractionalMedianPrintsAsDecimal(testCase)
            %The test-retest writer carries the same nine-column format as
            %the one-facet writer and read its median from trlinfo.med with
            %'%d', so a half-integer median printed as %e there too. The
            %fixture's median is bumped to x.5 to reach that cell.
            pd = localTrtData();
            med = pd.relsummary.group(1).event(1).trlinfo.med + 0.5;
            pd.relsummary.group(1).event(1).trlinfo.med = med;
            testCase.assertNotEqual(med, round(med), 'precondition: fractional median');
            fig = localBuildOutputFig(testCase, ...
                @() psyrat_trt_reloverallt('psyrat_data', pd, 'gui', 1));
            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile, 'Save Table');
            lines = localCsvLines(csvfile);
            hdr = find(startsWith(lines, 'Label,'), 1);
            testCase.assertNotEmpty(hdr, 'table header row not found');
            cells = strsplit(lines{hdr+1}, ',');
            testCase.verifyEqual(cells{6}, sprintf('%g', med), ...
                'a fractional median must print as a plain decimal in the trt CSV');
            testCase.verifyFalse(contains(lines{hdr+1}, 'e+'), ...
                'no scientific notation may appear in an exported CSV row');
            localCleanFiles({csvfile});
        end

        function testTrtOverallCsvNumericPathUnchangedAndDashSafe(testCase)
            % G47 latent-hardening pin for the trt writer's identical CSV
            % loop. No reachable route feeds psyrat_trt_reloverallt a '---'
            % cell today (trt_diff overall tables route through
            % psyrat_depoverallt on both the GUI and report paths), so this
            % pins two things directly against the trt save callback: (a) an
            % all-numeric row stays byte-identical to the pre-fix rendering,
            % and (b) a dash-bearing table injected through the callback's
            % own argument signature survives as one intact line rather than
            % the three-line numeric-code garble.
            pd = localTrtData();
            fig = localBuildOutputFig(testCase, ...
                @() psyrat_trt_reloverallt('psyrat_data', pd, 'gui', 1));

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile, 'Save Table');
            lines = localCsvLines(csvfile);
            hdr = find(startsWith(lines, 'Label,'), 1);
            testCase.assertNotEmpty(hdr, 'table header row not found');
            ev = pd.relsummary.group(1).event(1);
            relstr = sprintf(' %0.2f CI [%0.2f %0.2f]', ...
                ev.rel.m, ev.rel.ll, ev.rel.ul);
            expected = sprintf('%s,%d,%d,%s,%0.2f,%d,%0.2f,%d,%d', ...
                'Measurement', ev.goodn, numel(pd.relsummary.group(1).badids), ...
                relstr, ev.trlinfo.mean, ev.trlinfo.med, ev.trlinfo.std, ...
                ev.trlinfo.min, ev.trlinfo.max);
            testCase.verifyEqual(lines{hdr+1}, expected, ...
                'trt numeric row must stay byte-identical to the pre-fix rendering');

            % (b) the Save button stores {fcn, psyrat_data, overalltable};
            % re-fire the same callback with the badn and trial-summary cells
            % dashed out, the way the depoverallt diff rows dash theirs
            btns = findall(fig, 'Style', 'pushbutton');
            cb = [];
            for b = btns(:)'
                s = b.String; if iscell(s), s = strjoin(s, ' '); end
                if contains(string(s), 'Save Table'); cb = b.Callback; break; end
            end
            testCase.assertFalse(isempty(cb), 'Save Table button not found');
            tbl = cb{3};
            for col = [3 5 6 7 8 9]
                tbl.(col){1} = '---';
            end
            global PSYRAT_FAKE_SAVE
            dashfile = [tempname '.csv'];
            PSYRAT_FAKE_SAVE = dashfile;
            feval(cb{1}, [], [], cb{2}, tbl);
            dlines = localCsvLines(dashfile);
            dhdr = find(startsWith(dlines, 'Label,'), 1);
            testCase.assertNotEmpty(dhdr, 'table header row not found in dashed CSV');
            expected_dash = sprintf('%s,%d,---,%s,---,---,---,---,---', ...
                'Measurement', ev.goodn, relstr);
            testCase.verifyEqual(dlines{dhdr+1}, expected_dash, ...
                'a dash-bearing row must survive as one intact CSV line');
            testCase.verifyEqual(numel(dlines) - dhdr, 1, ...
                'the dashed table has one row; extra lines mean the formatspec recycled');

            localCleanFiles({csvfile, dashfile});
        end

        function testTrtCutoffSaveCarriesProvenance(testCase)
            % C4: the test-retest trial-cutoff export header must carry the run
            % provenance (seed/engine/convergence) in BOTH CSV and XLSX, so a
            % saved table is self-describing. Exercises the filehead/fprintf
            % (dual-branch) header path.
            pd = localTrtData();
            fig = localBuildOutputFig(testCase, ...
                @() psyrat_trt_relcutofft('psyrat_data', pd, 'gui', 1));

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            localVerifyProvenance(testCase, fileread(csvfile), 'CSV');

            xlsxfile = [tempname '.xlsx'];
            localFireSave(testCase, fig, xlsxfile);
            localVerifyProvenance(testCase, localXlsxText(xlsxfile), 'XLSX');

            localCleanFiles({csvfile, xlsxfile});
        end

        function testVarianceSaveCarriesProvenance(testCase)
            % C4: the one-facet group-level variance export header carries the
            % run provenance (CSV path).
            pd = localIcData();
            fig = localBuildOutputFig(testCase, ...
                @() psyrat_variancet('psyrat_data', pd, 'gui', 1));

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            localVerifyProvenance(testCase, fileread(csvfile), 'CSV');

            localCleanFiles({csvfile});
        end

        function testDodVarcompSaveCarriesProvenance(testCase)
            % C4: the difference-of-differences variance-component export header
            % carries the run provenance. Exercises the 'head'-cell (single
            % builder, both formats) header path used by the DoD/viewer exporters.
            pd = PsyRATTestDataFactory.makeDoDRelData();
            pd = localAddSaveHeaderFields(pd);
            map = local_dod_map(pd.rel);
            fig = localBuildOutputFig(testCase, ...
                @() psyrat_dod_varcompt('psyrat_data', pd, 'groupidx', 1, ...
                'map', map, 'gui', 1));

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            localVerifyProvenance(testCase, fileread(csvfile), 'CSV');

            xlsxfile = [tempname '.xlsx'];
            localFireSave(testCase, fig, xlsxfile);
            localVerifyProvenance(testCase, localXlsxText(xlsxfile), 'XLSX');

            localCleanFiles({csvfile, xlsxfile});
        end

        function testCriterionSingSaveWritesHeaderAndKeepsAllColumns(testCase)
            % RC-28: the criterion table export used to build its header inside
            % the .xlsx branch only, so the .csv branch wrote a bare table with
            % no header at all -- no seed, engine, priors or convergence
            % verdict. Both formats must now carry the provenance block.
            %
            % This also guards the mechanism change: the CSV branch went from a
            % bare writetable to header lines plus 'WriteMode','append'. If the
            % append lost the column-name row or the data rows, the export would
            % still "work" and still pass a provenance-only assertion, so the
            % column names and the row count are asserted explicitly.
            pd = localCriterionSingData();
            fig = localBuildOutputFig(testCase, @() psyrat_criterionfigures( ...
                'psyrat_prefs', localCriterionPrefs(0.5, 1), ...
                'psyrat_data', pd, 'analysis', 'sing'));

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            testCase.assertTrue(isfile(csvfile), 'Criterion CSV not written');
            csvtxt = fileread(csvfile);
            localVerifyProvenance(testCase, csvtxt, 'CSV');
            testCase.verifyTrue(contains(csvtxt,'Criterion Cutoff: 0.5000'), ...
                'CSV header should carry the criterion cutoff');
            testCase.verifyTrue(contains(csvtxt,'PsyRAT Toolbox v'), ...
                'CSV header should carry the toolbox version');
            testCase.verifyTrue(contains(csvtxt, localScopeSnippet()), ...
                'CSV header should carry the criterion scope line');

            % the appended table must still have its column-name row ...
            lines = localCsvLines(csvfile);
            idx = find(startsWith(lines,'Label,'), 1);
            testCase.assertNotEmpty(idx, ...
                'CSV lost the table column-name row (WriteVariableNames)');
            for name = {'Label','n_Included','Criterion_Cutoff','Trials_Used','Mu', ...
                    'Dependability_PointEstimate','Dependability_LowerCI', ...
                    'Dependability_UpperCI'}
                testCase.verifyTrue(contains(lines{idx}, name{1}), ...
                    sprintf('CSV column-name row is missing %s', name{1}));
            end
            % ... and its data. One group x one event = exactly one data row.
            testCase.verifyEqual(numel(lines) - idx, 1, ...
                'CSV should hold exactly one criterion data row');

            % RC-36: n_Included is asserted by VALUE, not just by presence. A
            % column-name-only check would pass on an all-NaN column, which is
            % exactly the failure local_goodn's isfield guard can produce if the
            % field it reads is ever renamed. The fixture retains 3 of 4
            % participants (PsyRATTestDataFactory ev.goodn = 3), and the column
            % is second, so the value sits in field 2 of the single data row.
            hdr = strsplit(strtrim(lines{idx}), ',');
            row = strsplit(strtrim(lines{idx+1}), ',');
            ncol = find(strcmp(hdr,'n_Included'), 1);
            testCase.assertEqual(ncol, 2, ...
                'n_Included should be the second exported column');
            testCase.verifyEqual(str2double(row{ncol}), 3, ...
                'n_Included should export the retained participant count (goodn)');

            localCleanFiles({csvfile});
        end

        function testCriterionSingXlsxKeepsOneBlankSpacerRow(testCase)
            % RC-28: inserting the provenance block lengthened the header, and
            % the .xlsx branch offsets the table by length(filehead)+1 against a
            % single trailing blank entry. That pairing must still leave exactly
            % one blank spacer row above the table, or the header block would
            % overwrite the column names (or strand the table further down).
            pd = localCriterionSingData();
            fig = localBuildOutputFig(testCase, @() psyrat_criterionfigures( ...
                'psyrat_prefs', localCriterionPrefs(0.5, 1), ...
                'psyrat_data', pd, 'analysis', 'sing'));

            xlsxfile = [tempname '.xlsx'];
            localFireSave(testCase, fig, xlsxfile);
            testCase.assertTrue(isfile(xlsxfile), 'Criterion XLSX not written');
            localVerifyProvenance(testCase, localXlsxText(xlsxfile), 'XLSX');

            raw = readcell(xlsxfile);
            r = localHeaderRowIndex(raw, 'Label');
            testCase.assertNotEmpty(r, 'XLSX lost the table column-name row');
            testCase.assertGreaterThan(r, 2, ...
                'Table starts too high for a header block to fit above it');
            testCase.verifyTrue(localCellIsBlank(raw{r-1,1}), ...
                'Exactly one blank spacer row should sit above the table');
            testCase.verifyFalse(localCellIsBlank(raw{r-2,1}), ...
                'Only ONE spacer row should sit above the table');

            localCleanFiles({xlsxfile});
        end

        function testCriterionTrtSaveRecordsCoefficientAndScope(testCase)
            % RC-28: the two-facet criterion export is directly reltype-driven
            % (psyrat_criterionfigures overwrites relsummary.reltype on its
            % local copy), so the Coefficient line is the point of the fix here.
            %
            % reltype is set to 3 (CES) while the fixture's relsummary carries
            % the stale reltype_name 'ic'. Asserting CES therefore also pins
            % that the label is read from the NUMERIC field, which is the whole
            % reason psyrat_trt_coeflabel reads it that way.
            pd = localCriterionTrtData();
            fig = localBuildOutputFig(testCase, @() psyrat_criterionfigures( ...
                'psyrat_prefs', localCriterionPrefs(0.5, 3), ...
                'psyrat_data', pd, 'analysis', 'trt'));

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            csvtxt = fileread(csvfile);
            localVerifyProvenance(testCase, csvtxt, 'CSV');
            testCase.verifyTrue(contains(csvtxt, ...
                'Coefficient of Equivalence and Stability (CES)'), ...
                'Header should name the selected coefficient from numeric reltype');
            % gcoeff is pinned to 1 for criterion output regardless of the
            % user's selection, so the error type must read absolute.
            testCase.verifyTrue(contains(csvtxt,'Dependability (absolute error)'), ...
                'Criterion header must report absolute error');
            testCase.verifyFalse(contains(csvtxt, ...
                'Coefficient of Equivalence and Stability (CES), Generalizability'), ...
                'Criterion header must not report relative error for the coefficient');
            % the cut-scores hardcode nocc = 1, so the occasion clause must agree
            testCase.verifyTrue(contains(csvtxt,'n''_o = 1'), ...
                'Criterion header should record the single-occasion estimand');
            testCase.verifyTrue(contains(csvtxt, localScopeSnippet()), ...
                'CSV header should carry the criterion scope line');

            localCleanFiles({csvfile});
        end

        function testCriterionTrtScopeLineQualifiesDifferenceClause(testCase)
            % RC-28, the reason the scope line exists. psyrat_trt_coeflabel
            % appends a "difference score: Generalizability (relative error)"
            % clause whenever relsummary.diffgcoeff differs from gcoeff -- and
            % on the criterion route gcoeff is pinned to 1, so a user who picked
            % relative error for difference scores gets that clause on a table
            % whose every number is an absolute-error cut-score.
            %
            % The clause is deliberately NOT suppressed (it reports a real user
            % selection, and psyrat_trt_coeflabel is shared by every other
            % export). It is qualified instead, so this pins that the qualifier
            % is present whenever the clause is.
            pd = localCriterionTrtData();
            pd.relsummary.diffgcoeff = 2;   % relative error for difference scores
            fig = localBuildOutputFig(testCase, @() psyrat_criterionfigures( ...
                'psyrat_prefs', localCriterionPrefs(0.5, 3), ...
                'psyrat_data', pd, 'analysis', 'trt'));

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            csvtxt = fileread(csvfile);
            testCase.assertTrue(contains(csvtxt, ...
                'difference score: Generalizability (relative error)'), ...
                'Fixture did not reach the difference-score clause');
            testCase.verifyTrue(contains(csvtxt, localScopeSnippet()), ...
                'The relative-error clause must be qualified by the scope line');

            localCleanFiles({csvfile});
        end

        function testDynrelVarcompSaveRecordsNPrimeDefault(testCase)
            % The dynrel viewer exposes a D-study trials-n' override
            % ("Trials n'", blank = per-stratum median), and the appended
            % provenance can describe the curve as using the median n' (the
            % analysis-13 gamma subject-level S17 lines) -- but nothing in the
            % export recorded the n' ACTUALLY used, so under an override a reader
            % could not tell what the surface was evaluated at (pre-existing
            % gap adjudicated 2026-08-11 during the S17-FIX adversarial
            % review). The header must state the resolved per-stratum n' and
            % its source. Default path first: a blank field resolves to each
            % stratum's stored median count (the diff fixture stores
            % ntrials = 20). Both formats share one header
            % builder, so XLSX is asserted once here and the other dynrel tests
            % pin the CSV branch only.
            pd = localDynrelDiffData();
            fig = localDynrelViewerFig(testCase, pd);

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            testCase.assertTrue(isfile(csvfile), 'Dynrel variance CSV not written');
            csvtxt = fileread(csvfile);
            testCase.verifyTrue(contains(csvtxt, ...
                ['Reliability surface n'' (trials), stratum none: 20 ' ...
                '(default: per-stratum median)']), ...
                'CSV header should record the resolved default n'' per stratum');

            xlsxfile = [tempname '.xlsx'];
            localFireSave(testCase, fig, xlsxfile);
            testCase.verifyTrue(contains(localXlsxText(xlsxfile), ...
                ['Reliability surface n'' (trials), stratum none: 20 ' ...
                '(default: per-stratum median)']), ...
                'XLSX header should record the resolved default n'' per stratum');

            localCleanFiles({csvfile, xlsxfile});
        end

        function testDynrelVarcompSaveRecordsNPrimeOverride(testCase)
            % The heart of the adjudicated gap: under a user override the
            % appended provenance description ("the median n'") is not what the
            % surface used, so the header line must name the override value and
            % say it was one.
            pd = localDynrelDiffData();
            fig = localDynrelViewerFig(testCase, pd);
            localSetDynrelTrials(testCase, fig, '7');

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            csvtxt = fileread(csvfile);
            testCase.verifyTrue(contains(csvtxt, ...
                'Reliability surface n'' (trials), stratum none: 7 (user override)'), ...
                'CSV header should record the override n'' and its source');
            testCase.verifyFalse(contains(csvtxt, 'stratum none: 20'), ...
                'CSV header must not claim the default n'' under an override');

            localCleanFiles({csvfile});
        end

        function testDynrelVarcompSaveRejectsComplexNPrime(testCase)
            % Complex text in the "Trials n'" field must fall back to the
            % default, like any other invalid input.
            %
            % local_parse_pos documented "return default on empty/invalid" but
            % accepted complex literals: str2double('3-4i') returns 3-4i,
            % isfinite is true for it, and MATLAB's ordering comparison uses
            % only the REAL part, so 'isfinite(x) && x > 0' passed. Downstream
            % guards missed it too (psyrat_dynrel_summary's n' check is
            % isfinite-only, and max(1,round(2+3i)) stays complex), so the
            % complex n' reached the calculators and errored in quantile/prctile
            % ("Invalid data type. First argument must be an array of real
            % values") -- an unexplained crash instead of the documented
            % fallback. Fixed by adding isreal to the accept condition.
            %
            % Asserted on the export header rather than the rendered variance
            % table on purpose: that table is INVARIANT to n' (measured), so an
            % equality check on it cannot tell "fell back to the default" from
            % "accepted the value" and would pin only the absence of the crash.
            % The header states the resolved n', and the two tests above pin
            % both of its branches, so this is a discriminating assertion.
            pd = localDynrelDiffData();
            fig = localDynrelViewerFig(testCase, pd);
            localSetDynrelTrials(testCase, fig, '3-4i');

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);   % errored pre-fix
            testCase.assertTrue(isfile(csvfile), 'Dynrel variance CSV not written');
            csvtxt = fileread(csvfile);
            testCase.verifyTrue(contains(csvtxt, ...
                ['Reliability surface n'' (trials), stratum none: 20 ' ...
                '(default: per-stratum median)']), ...
                ['complex n'' text must fall back to the per-stratum median, ' ...
                'so the header must report the default and name it one']);
            testCase.verifyFalse(contains(csvtxt, 'user override'), ...
                'complex n'' text must not be recorded as a user override');

            localCleanFiles({csvfile});
        end

        function testDynrelVarcompSaveTracksEditedNPrimeAcrossSaves(testCase)
            % Regression for the buildSummary cache key: the key printed the
            % requested trials n' with num2str at its default precision, which
            % keeps about six significant digits, while psyrat_dynrel_summary
            % resolves the count with max(1,round(obs)). Two DIFFERENT requests
            % could therefore collide in the key while rounding to different n'
            % (20.499999 -> 20 but 20.5 -> 21; num2str prints both as '20.5'),
            % so the second save silently reused the first surface and its
            % header and table reported the stale n'. Drive the SAME viewer
            % instance through both saves -- the cache lives in the figure's
            % nested-function workspace, so recreating the viewer would hide
            % the bug -- and require the second export to carry the freshly
            % resolved n'.
            pd = localDynrelDiffData();
            fig = localDynrelViewerFig(testCase, pd);

            localSetDynrelTrials(testCase, fig, '20.499999');
            csv1 = [tempname '.csv'];
            localFireSave(testCase, fig, csv1);
            txt1 = fileread(csv1);
            testCase.assertTrue(contains(txt1, ...
                'Reliability surface n'' (trials), stratum none: 20 (user override)'), ...
                'Precondition: 20.499999 should resolve to n'' = 20');

            localSetDynrelTrials(testCase, fig, '20.5');
            csv2 = [tempname '.csv'];
            localFireSave(testCase, fig, csv2);
            txt2 = fileread(csv2);
            testCase.verifyTrue(contains(txt2, ...
                'Reliability surface n'' (trials), stratum none: 21 (user override)'), ...
                'Second save must recompute the surface: 20.5 rounds to n'' = 21');
            testCase.verifyFalse(contains(txt2, 'stratum none: 20 '), ...
                'Second save must not reuse the stale n'' = 20 surface from the cache');

            localCleanFiles({csv1, csv2});
        end

        function testDynrelCsvHeaderLinesStayOneSpreadsheetField(testCase)
            % B19: the CSV branch wrote the metadata header with bare
            % fprintf('%s\n', ...) and no quoting, so every comma-bearing
            % provenance line (the resolved n' line, Chains/Iterations, the
            % DoD cell mapping) fragmented across spreadsheet columns when
            % the file was opened the ordinary way. The header is the file's
            % provenance record -- the on-screen note explicitly directs the
            % reader to it -- so each line must survive a spreadsheet
            % round-trip as ONE field. Comma-free lines stay byte-identical;
            % the XLSX branch (writecell) was never affected.
            pd = localDynrelDiffData();
            fig = localDynrelViewerFig(testCase, pd);

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            lines = splitlines(fileread(csvfile));

            % the header block ends at the first of the two spacer lines the
            % export writes before the table
            hend = find(cellfun(@(s) isempty(strtrim(s)), lines), 1) - 1;
            testCase.assertTrue(~isempty(hend) && hend >= 1, ...
                'no header/table separator found in the CSV');
            head = lines(1:hend);

            % precondition: the fixture actually produces a comma-bearing
            % header line, or the quoting assertion could pass vacuously
            hascomma = any(cellfun(@(s) contains(localCsvUnquote(s), ','), head));
            testCase.assertTrue(hascomma, ...
                'fixture precondition broken: no comma-bearing header line');

            for i = 1:numel(head)
                testCase.verifyEqual(localCsvFieldCount(head{i}), 1, sprintf( ...
                    'header line %d fragments across spreadsheet columns: %s', ...
                    i, head{i}));
            end

            localCleanFiles({csvfile});
        end

        function testTableCsvHeaderLinesStayOneSpreadsheetField(testCase)
            % B19 residue (2026-08-17 ruling: fix every residual CSV writer):
            % same contract as the dynrel test above, for the table writers'
            % direct-fprintf and head-cell headers. The always-comma Chains
            % line ('Chains: 3, Iterations: 600') fragments on every export
            % pre-fix. Each spec fires the writer's real Save button through
            % the uiputfile shadow, then asserts the FIRST multi-field CSV
            % line is the data-header row (its known first-column marker) --
            % pre-fix a fragmented header line is multi-field and appears
            % ABOVE the marker row, so the marker assertion is the
            % discriminator; given it, every earlier line is one field by
            % construction. DoD fixtures mirror
            % testDodVarcompSaveCarriesProvenance.
            dodpd = localAddSaveHeaderFields(PsyRATTestDataFactory.makeDoDRelData());
            dodmap = local_dod_map(dodpd.rel);
            specs = { ...
                {'trt_relcutofft', @() psyrat_trt_relcutofft('psyrat_data', localTrtData(), 'gui', 1), 'Save Table', 'Label,'}; ...
                {'trt_reloverallt', @() psyrat_trt_reloverallt('psyrat_data', localTrtData(), 'gui', 1), 'Save Table', 'Label,'}; ...
                {'trt_reloverallt ids', @() psyrat_trt_reloverallt('psyrat_data', localTrtData(), 'gui', 1), 'Save IDs', 'Good IDs,'}; ...
                {'depcutofft', @() psyrat_depcutofft('psyrat_data', localIcData(), 'gui', 1), 'Save Table', 'Label,'}; ...
                {'depoverallt', @() psyrat_depoverallt('psyrat_data', localIcData(), 'gui', 1), 'Save Table', 'Label,'}; ...
                {'depoverallt ids', @() psyrat_depoverallt('psyrat_data', localIcData(), 'gui', 1), 'Save IDs', 'Good IDs,'}; ...
                {'variancet', @() psyrat_variancet('psyrat_data', localIcData(), 'gui', 1), 'Save Table', 'Label,'}; ...
                {'dod_varcompt', @() psyrat_dod_varcompt('psyrat_data', dodpd, 'groupidx', 1, 'map', dodmap, 'gui', 1), 'Save Table', 'Component,'}; ...
                {'dod_observedt', @() psyrat_dod_observedt('psyrat_data', dodpd, 'groupidx', 1, 'map', dodmap, 'gui', 1), 'Save Table', 'Label,'}; ...
                {'criterion table', @() psyrat_criterionfigures('psyrat_prefs', localCriterionPrefs(0.5, 1), 'psyrat_data', localCriterionSingData(), 'analysis', 'sing'), 'Save', 'Label,'}; ...
                };
            for s = 1:size(specs,1)
                what = specs{s}{1};
                fig = localBuildOutputFig(testCase, specs{s}{2});
                csvfile = [tempname '.csv'];
                localFireSave(testCase, fig, csvfile, specs{s}{3});
                localAssertCsvHeaderOneField(testCase, csvfile, specs{s}{4}, what);
                localCleanFiles({csvfile});
            end
        end

        function testSplitsCsvHeaderLinesStayOneSpreadsheetField(testCase)
            % B19 residue: both splits-viewer save buttons route through the
            % shared exportSplitsTable header block, so one driven button
            % covers the writer. Same extraction as the dynrel test: the
            % header ends at the first blank spacer line before the appended
            % table; the fixture's Chains line carries the comma that makes
            % this RED without the quoting (the summary-only factory data
            % lacks the nchains/niter provenance fields the header reads, so
            % add them the same way the other save tests do).
            pd = localAddSaveHeaderFields(PsyRATTestDataFactory.makeSplitsSummaryData());
            fig = localSplitsViewerFig(testCase, pd);
            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile, 'Save coefficients');
            lines = splitlines(fileread(csvfile));
            hend = find(cellfun(@(s) isempty(strtrim(s)), lines), 1) - 1;
            testCase.assertTrue(~isempty(hend) && hend >= 1, ...
                'no header/table separator found in the splits CSV');
            head = lines(1:hend);
            hascomma = any(cellfun(@(s) contains(localCsvUnquote(s), ','), head));
            testCase.assertTrue(hascomma, ...
                'fixture precondition broken: no comma-bearing splits header line');
            for i = 1:numel(head)
                testCase.verifyEqual(localCsvFieldCount(head{i}), 1, sprintf( ...
                    'splits header line %d fragments across spreadsheet columns: %s', ...
                    i, head{i}));
            end
            localCleanFiles({csvfile});
        end

        function testCsvHeaderQuotingCoversEverySite(testCase)
            % census pin for the B19-residue fix: every CSV header write site
            % routes its comma-capable lines through psyrat_quote_csv_header,
            % including the sibling save callbacks the marker test above does
            % not drive (variancet's sscoeffs/allvar/alldiff blocks and the
            % two Save IDs blocks share their file's count). A dropped or
            % added wrap changes a count and trips this pin. The dynrel
            % viewer deliberately keeps its shipped inline copy of the same
            % conditional (pinned behaviorally by its own test above), so it
            % stays at zero calls.
            specs = { ...
                {fullfile('subroutines','tables','psyrat_fprintf_provenance.m'), 1}; ...
                {fullfile('subroutines','tables','psyrat_depcutofft.m'), 3}; ...
                {fullfile('subroutines','tables','psyrat_trt_relcutofft.m'), 3}; ...
                {fullfile('subroutines','tables','psyrat_trt_reloverallt.m'), 7}; ...
                {fullfile('subroutines','tables','psyrat_depoverallt.m'), 6}; ...
                {fullfile('subroutines','tables','psyrat_variancet.m'), 11}; ...
                {fullfile('subroutines','tables','psyrat_dod_varcompt.m'), 1}; ...
                {fullfile('subroutines','tables','psyrat_dod_observedt.m'), 1}; ...
                {fullfile('subroutines','plotting','psyrat_criterionfigures.m'), 1}; ...
                {fullfile('subroutines','guis','psyrat_startview_splits.m'), 1}; ...
                {fullfile('subroutines','guis','psyrat_startview_dynrel.m'), 0}; ...
                };
            here = fileparts(mfilename('fullpath'));
            for s = 1:size(specs,1)
                src = fileread(fullfile(fileparts(here), specs{s}{1}));
                testCase.verifyEqual( ...
                    numel(strfind(src, 'psyrat_quote_csv_header(')), specs{s}{2}, ...
                    sprintf('%s: unexpected psyrat_quote_csv_header call count', ...
                    specs{s}{1}));
            end
        end

        function testSplitsSaveTracksEditedSplitsNPrimeAcrossSaves(testCase)
            % The data-splits viewer's buildSummary carries the same cache-key
            % defect the dynrel test above pins, so the fix needs its own
            % coverage here: psyrat_splits_summary resolves the split count
            % with max(1,round(obs)), while the key printed it at num2str's
            % default precision, so 20.499999 and 20.5 shared a key ('20.5')
            % despite resolving to 20 vs 21 and the second save served the
            % stale summary. Drive the SAME viewer instance through both saves
            % -- the cache lives in the figure's nested-function workspace, so
            % recreating the viewer would hide the bug. The resolved count is
            % asserted from the export's own n_splits column.
            pd = localSplitsSingleData();
            fig = localSplitsViewerFig(testCase, pd);

            localSetEditByLabel(testCase, fig, 'Splits n', '20.499999');
            csv1 = [tempname '.csv'];
            localFireSave(testCase, fig, csv1, 'Save variance table');
            testCase.assertEqual( ...
                localSplitsDesignValue(testCase, csv1, 'n_splits'), 20, ...
                'Precondition: 20.499999 should resolve to n'' = 20');

            localSetEditByLabel(testCase, fig, 'Splits n', '20.5');
            csv2 = [tempname '.csv'];
            localFireSave(testCase, fig, csv2, 'Save variance table');
            testCase.verifyEqual( ...
                localSplitsDesignValue(testCase, csv2, 'n_splits'), 21, ...
                'Second save must recompute the summary: 20.5 rounds to n'' = 21');

            localCleanFiles({csv1, csv2});
        end

        function testSplitsSaveTracksEditedOccasionsNPrimeAcrossSaves(testCase)
            % The data-splits key has a SECOND half the dynrel one does not:
            % the occasion n', resolved by local_resolve_nocc with the same
            % max(1,round(nocc)). Only the multi-occasion design (trt_splits)
            % exposes that selector, so it needs its own fixture and lane --
            % fixing the splits half alone would leave this one colliding. The
            % splits n' field is left blank so the occasion value is the only
            % thing differing between the two saves.
            pd = localSplitsTrtData();
            fig = localSplitsViewerFig(testCase, pd);

            localSetEditByLabel(testCase, fig, 'Occasions n', '3.499999');
            csv1 = [tempname '.csv'];
            localFireSave(testCase, fig, csv1, 'Save variance table');
            testCase.assertEqual( ...
                localSplitsDesignValue(testCase, csv1, 'n_occasions'), 3, ...
                'Precondition: 3.499999 should resolve to n''_o = 3');

            localSetEditByLabel(testCase, fig, 'Occasions n', '3.5');
            csv2 = [tempname '.csv'];
            localFireSave(testCase, fig, csv2, 'Save variance table');
            testCase.verifyEqual( ...
                localSplitsDesignValue(testCase, csv2, 'n_occasions'), 4, ...
                'Second save must recompute the summary: 3.5 rounds to n''_o = 4');

            localCleanFiles({csv1, csv2});
        end

        function testDynrelTrtVarcompSaveRecordsOccasions(testCase)
            % Two-facet variants: the surface additionally depends on the
            % occasion n' (viewer default 1 -- an owner decision distinct from
            % the fixture's observed count of 2), so the header line must carry
            % it alongside the trials n' or the export is still not
            % interpretable from its own header.
            pd = localDynrelTrtData();
            fig = localDynrelViewerFig(testCase, pd);

            csvfile = [tempname '.csv'];
            localFireSave(testCase, fig, csvfile);
            csvtxt = fileread(csvfile);
            testCase.verifyTrue(contains(csvtxt, ...
                ['Reliability surface n'' (trials), stratum none: 8 ' ...
                '(default: per-stratum median)']), ...
                'CSV header should record the trials n'' for the two-facet variant');
            testCase.verifyTrue(contains(csvtxt, '; n'' (occasions): 1'), ...
                'CSV header should record the occasion n'' the surface used');

            localCleanFiles({csvfile});
        end
    end
end

% -- local helpers -------------------------------------------------------------

function teardownSaveShadow(testCase)
clear global PSYRAT_FAKE_SAVE
if ~isempty(testCase.ShadowDir) && isfolder(testCase.ShadowDir)
    % the base-class teardown may already have restored the original path,
    % so only rmpath when the shadow dir is still present
    if any(strcmp(strsplit(path, pathsep), testCase.ShadowDir))
        rmpath(testCase.ShadowDir);
    end
    rmdir(testCase.ShadowDir,'s');
end
close(findall(0,'Type','figure','Tag','psyrat_output'));
end

function pd = localTrtData()
pd = PsyRATTestDataFactory.makeTRTSummaryData();
pd = localAddSaveHeaderFields(pd);
end

function pd = localIcData()
pd = PsyRATTestDataFactory.makeICSummaryDataWithTrial();
pd = localAddSaveHeaderFields(pd);
end

function pd = localIcDiffData()
% one-facet difference fixture (ic_diff): two event rows plus the diff-score
% row whose badn and trial-summary cells carry the '---' placeholder (G47)
pd = PsyRATTestDataFactory.makeICDiffSummaryData();
pd = localAddSaveHeaderFields(pd);
end

function pd = localAddSaveHeaderFields(pd)
% the save-callback file headers read these provenance fields, which the
% summary-only factory data does not include
pd.ver = 'test';
pd.rel.filename = 'export_callback_fixture.psyrat';
pd.rel.nchains = 2;
pd.rel.niter = 600;
% run-provenance fields read by psyrat_provenance_lines (C4): a known seed,
% engine, prior status, and convergence verdict so the exported header lines
% are deterministic to assert against.
pd.rel.seed = 12345;
pd.rel.engine = 'cmdstan';
pd.rel.priors = psyrat_default_priors;
pd.rel.out.conv.converged = 1;
pd.rel.out.conv.ndivergent = 0;
end

function pd = localDynrelDiffData()
% group-level dynamic difference fixture (ic_diff_dynrel): one stratum
% ('none'), stored per-stratum median count ntrials = 20, no per-participant
% table, so the viewer shows exactly one Save button (the variance table)
pd = PsyRATTestDataFactory.makeDynrelDiffSummaryData();
pd = localAddSaveHeaderFields(pd);
end

function pd = localDynrelTrtData()
% two-facet dynamic reliability fixture (ic_dynrel_trt): ntrials = 8 with an
% observed occasion count of 2; the viewer's occasion-n' selector defaults to
% 1 occasion, which is what the export header must report
pd = PsyRATTestDataFactory.makeDynrelTrtSummaryData();
pd = localAddSaveHeaderFields(pd);
end

function fig = localDynrelViewerFig(testCase, pd)
% launch the dynamic-reliability viewer and return its figure. The viewer
% figure is located by its window name (it carries Tag psyrat_gui_dynrelview,
% which the psyrat_output table figures do not share), so it is
% located by its window name; the base-class teardown closes it.
name = 'Dynamic Reliability: View Results';
close(findall(0,'Type','figure','Name',name));
psyrat_startview_dynrel('psyrat_prefs', psyrat_defaults(), 'psyrat_data', pd);
fig = findall(0,'Type','figure','Name',name);
testCase.assertNotEmpty(fig, 'Dynrel viewer figure not created');
fig = fig(1);
end

function localSetDynrelTrials(testCase, fig, str)
% type an override into the dynrel viewer's "Trials n'" edit box
localSetEditByLabel(testCase, fig, 'Trials n', str);
end

function localSetEditByLabel(testCase, fig, labelprefix, str)
% type into the viewer edit box that shares its layout row with the label
% starting with labelprefix. The viewers' uicontrols carry no tags, so the row
% pairing is the only thing that disambiguates their several empty-by-default
% edit fields (dynrel: n' vs cutoff; data-splits: credible interval vs splits
% n' vs occasions n'). findall's -regexp is not supported for the uicontrol
% String property (it warns and falls back to ordinary comparison), so scan
% the labels. Same approach as the like-named helper in
% TestViewerEditFieldFallback; each file keeps its own copy.
lbl = [];
for h = findall(fig, 'Style', 'text')'
    s = h.String; if iscell(s); s = strjoin(s, ' '); end
    if startsWith(string(s), labelprefix)
        lbl = h;
        break;
    end
end
testCase.assertNotEmpty(lbl, ...
    sprintf('"%s" label not found on the viewer', labelprefix));
p = get(lbl(1), 'Position');
eds = findall(fig, 'Style', 'edit');
match = arrayfun(@(h) abs(h.Position(2) - p(2)) < 5, eds);
testCase.assertEqual(nnz(match), 1, ...
    sprintf('"%s" edit box not uniquely identified by its row', labelprefix));
set(eds(match), 'String', str);
end

function pd = localSplitsSingleData()
% single-occasion nonparallel data-splits fixture (ic_splits): one stratum,
% 3 persons x 4 splits with varying items-per-split, so the observed median
% splits-per-person is 4 and the viewer shows no occasions-n' selector
pd = PsyRATTestDataFactory.makeSplitsSummaryData();
pd = localAddSaveHeaderFields(pd);
end

function pd = localSplitsTrtData()
% multi-occasion nonparallel data-splits fixture (trt_splits): 3 persons x
% 2 occasions x 4 splits, the design that exposes the occasions-n' selector
% alongside the splits-n' one
pd = PsyRATTestDataFactory.makeSplitsTRTSummaryData();
pd = localAddSaveHeaderFields(pd);
end

function fig = localSplitsViewerFig(testCase, pd)
% launch the data-splits viewer and return its figure. Like the dynrel viewer
% it is located by its window name rather than its Tag; the base-class
% teardown closes it.
name = 'Data-Splits Reliability: View Results';
close(findall(0,'Type','figure','Name',name));
psyrat_startview_splits('psyrat_prefs', psyrat_defaults(), 'psyrat_data', pd);
fig = findall(0,'Type','figure','Name',name);
testCase.assertNotEmpty(fig, 'Data-splits viewer figure not created');
fig = fig(1);
end

function v = localSplitsDesignValue(testCase, csvfile, colname)
% Read a design column (n_splits / n_occasions) out of the first data row of
% an exported data-splits variance CSV. The metadata header lines sit above
% the table's column-name row, which psyrat_splits_varcomp_concat starts with
% 'stratum'. Only the four leading design columns are read, so a comma inside
% one of the later free-text component names cannot shift the parsed field.
lines = localCsvLines(csvfile);
idx = find(startsWith(lines,'stratum,'), 1);
testCase.assertNotEmpty(idx, 'Splits CSV lost the table column-name row');
hdr = strsplit(strtrim(lines{idx}), ',');
c = find(strcmp(hdr, colname), 1);
testCase.assertNotEmpty(c, sprintf('Splits CSV has no %s column', colname));
testCase.assertLessThanOrEqual(c, 4, ...
    sprintf('%s should sit among the four leading design columns', colname));
row = strsplit(strtrim(lines{idx+1}), ',');
v = str2double(row{c});
end

function pd = localCriterionSingData()
pd = PsyRATTestDataFactory.makeICSummaryDataWithTrial();
pd = localAddSaveHeaderFields(pd);
end

function pd = localCriterionTrtData()
% carries relsummary.reltype_name = 'ic', which psyrat_criterionfigures leaves
% stale when it overwrites the numeric reltype on its local copy
pd = PsyRATTestDataFactory.makeTRTSummaryData();
pd = localAddSaveHeaderFields(pd);
end

function prefs = localCriterionPrefs(cutoff, reltype)
% criterion viewer preferences. gcoeff is deliberately set to 2 (relative):
% criterion output must ignore it and stay absolute, which is what the scope
% line asserts in the exported header.
prefs = psyrat_defaults();
prefs.view.plotcriterion = 0;
prefs.view.tablecriterion = 1;
prefs.view.criterioncutoff = cutoff;
prefs.view.reltype = reltype;
prefs.view.gcoeff = 2;
prefs.view.depcentmeas = 1;
prefs.view.relcentmeas = 1;
end

function s = localScopeSnippet()
% distinctive fragment of the criterion scope line, kept in one place so the
% four criterion tests cannot drift apart from each other
s = 'Scope: criterion cut-scores are always absolute-error';
end

function lines = localCsvLines(csvfile)
% non-empty lines of an exported CSV, so header/table positions can be indexed
lines = strsplit(fileread(csvfile), newline);
lines = lines(~cellfun(@(s) isempty(strtrim(s)), lines));
end

function r = localHeaderRowIndex(raw, firstlabel)
% row INDEX of the table header (localHeaderRow returns its contents instead)
r = [];
for k = 1:size(raw,1)
    c = raw{k,1};
    if (ischar(c) || isstring(c)) && strcmp(char(c), firstlabel)
        r = k;
        return;
    end
end
end

function tf = localCellIsBlank(c)
% readcell returns missing for an empty spreadsheet cell; a blank header entry
% written by writecell can also come back as an empty char
tf = (isa(c,'missing') && ismissing(c)) || ...
    ((ischar(c) || isstring(c)) && isempty(strtrim(char(c))));
end

function fig = localBuildOutputFig(testCase, builderThunk)
% run a table builder (gui=1) and return its freshly created output figure
close(findall(0,'Type','figure','Tag','psyrat_output'));
builderThunk();
fig = findobj(0,'Type','figure','Tag','psyrat_output');
testCase.assertNotEmpty(fig, 'Builder did not create a psyrat_output figure');
fig = fig(1);
end

function localFireSave(testCase, fig, outfile, labelfrag)
% Fire the Save button whose label contains labelfrag. The default 'Save' is
% unambiguous on the figures that carry a single Save button; viewers with
% more than one (the data-splits viewer saves coefficients OR the variance
% table) must name the one they mean, since the order findall returns the
% buttons in is not part of any contract.
if nargin < 4; labelfrag = 'Save'; end
global PSYRAT_FAKE_SAVE
PSYRAT_FAKE_SAVE = outfile;
btns = findall(fig,'Style','pushbutton');
fired = false;
for b = btns(:)'
    s = b.String; if iscell(s), s = strjoin(s,' '); end
    if contains(string(s),labelfrag)
        cb = b.Callback;
        if iscell(cb)
            % table-figure Save buttons store {fcn, extra args...}
            feval(cb{1}, b, [], cb{2:end});
        else
            % viewer nested-callback Save buttons store a bare function handle
            cb(b, []);
        end
        fired = true;
        break;
    end
end
testCase.assertTrue(fired, ...
    sprintf('No "%s" button found on the output figure', labelfrag));
end

function hdr = localHeaderRow(raw, firstlabel)
% return the sheet row whose first cell equals firstlabel (the table header)
hdr = {};
for r = 1:size(raw,1)
    c = raw{r,1};
    if (ischar(c) || isstring(c)) && strcmp(char(c), firstlabel)
        hdr = raw(r,:);
        return;
    end
end
end

function tf = localHasCell(row, label)
tf = any(cellfun(@(x) (ischar(x) || isstring(x)) && strcmp(char(x), label), row));
end

function localCleanFiles(files)
for k = 1:numel(files)
    if isfile(files{k}), delete(files{k}); end
end
end

function localVerifyProvenance(testCase, txt, fmt)
% assert the standard run-provenance lines (C4) appear in an exported header,
% using the deterministic fixture values from localAddSaveHeaderFields
testCase.verifyTrue(contains(txt,'Seed: 12345'), ...
    sprintf('%s header should record the estimation seed', fmt));
testCase.verifyTrue(contains(txt,'Estimation engine: cmdstan'), ...
    sprintf('%s header should record the estimation engine', fmt));
testCase.verifyTrue(contains(txt,'Converged: Yes'), ...
    sprintf('%s header should record the convergence verdict', fmt));
end

function txt = localXlsxText(xlsxfile)
% join all text cells of an xlsx into one searchable string. The header lines
% sit in column A with empty/numeric cells elsewhere (the table below), so mask
% to char/string cells before joining (string() of a numeric/missing cell would
% otherwise pollute or break the join).
raw = readcell(xlsxfile);
mask = cellfun(@(c) ischar(c) || isstring(c), raw);
txt = strjoin(cellstr(string(raw(mask))), newline);
end

function map = local_dod_map(rel)
% numeric contrast indices into rel.events for the four DoD events (mirrors
% psyrat_report's local_dod_tables)
events = cellstr(string(rel.events(:)));
mapnames = cellstr(string(rel.dod_map(:)));
[~, map] = ismember(mapnames, events);
map = map(:)';
end

function n = localCsvFieldCount(line)
% count CSV fields the way a spreadsheet parses them: commas OUTSIDE
% double-quoted regions split fields; '""' inside a quoted region is an
% escaped quote (it toggles the state twice, so the net effect is correct)
inq = false;
n = 1;
for k = 1:numel(line)
    c = line(k);
    if c == '"'
        inq = ~inq;
    elseif c == ',' && ~inq
        n = n + 1;
    end
end
end

function localAssertCsvHeaderOneField(testCase, csvfile, marker, what)
% B19-residue helper: the header block is every line ABOVE the data-header
% row, which is anchored by its known first-column marker (e.g. 'Label,');
% every header line must parse as ONE spreadsheet field. The precondition
% probes the raw payloads (unquoted), so it is satisfied both pre-fix (the
% raw Chains line carries the comma) and post-fix (the quoted line does),
% and the per-line verification -- not an aborting assertion -- is what
% goes RED pre-fix, on the fragmenting line itself.
lines = splitlines(fileread(csvfile));
markerrow = find(cellfun(@(s) strncmp(s, marker, numel(marker)), lines), 1);
testCase.assertTrue(~isempty(markerrow) && markerrow >= 2, sprintf( ...
    '%s: no data-header row starting with "%s" found in the CSV', what, marker));
head = lines(1:markerrow-1);
hascomma = any(cellfun(@(s) contains(localCsvUnquote(s), ','), head));
testCase.assertTrue(hascomma, sprintf( ...
    '%s: fixture precondition broken -- no comma-bearing header line', what));
for i = 1:numel(head)
    testCase.verifyEqual(localCsvFieldCount(head{i}), 1, sprintf( ...
        '%s: header line %d fragments across spreadsheet columns: %s', ...
        what, i, head{i}));
end
end

function s = localCsvUnquote(line)
% strip one layer of CSV quoting so the precondition probe sees the payload
% whether or not the writer quotes it
s = line;
if numel(s) >= 2 && s(1) == '"' && s(end) == '"'
    s = strrep(s(2:end-1), '""', '"');
end
end
