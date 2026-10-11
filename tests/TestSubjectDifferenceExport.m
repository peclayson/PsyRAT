classdef TestSubjectDifferenceExport < PsyRATTestBase
    %G107 (analysis 8's per-participant DIFFERENCE coefficients were missing
    %from both subject exports) and G109 (Save Subject Coefficients could
    %overwrite a file the dialog never named, wrote text IDs unescaped, and,
    %like the other CSV-only Save callbacks, truncated its target before the
    %new file was complete). Written RED-first.
    %
    %CALLBACK TESTS, NOT NATIVE ACCEPTANCE. Every GUI test here builds the
    %real window with gui = 1 and fires the real button's own Callback with
    %feval; uiputfile, questdlg and uiwait are recording doubles
    %(G107ExportFixtures.installShadows). That proves what the callback does
    %for a given dialog answer. It is not a click, it does not show the
    %native macOS save panel, and whether that panel ever returns a name
    %without the filter's extension is not observed here; that needs a
    %native click-through. The headless tests drive the real psyrat_report.
    %The write-failure tests add recording I/O doubles that forward to the
    %built-ins (G107ExportFixtures.installFopenWatch, installIoFaults) and an
    %errordlg recorder that counts dialogs (installErrordlgRecorder). The
    %helper's own failure branches are tested directly in
    %TestCsvExportStageSafety. Nothing here fits a model.
    %
    %REFERENCES. Expected per-participant posterior means and equal-tailed
    %limits come from G107ExportFixtures.diffReference, recomputed from the
    %raw draws and from trial counts tallied from rel.data. The stored
    %summary table is checked against that reference first, then each file
    %against both. Table-to-file parity alone is not used as the scientific
    %reference. The outside-MATLAB reopen (Python csv and openpyxl) runs
    %separately on the files G107_write_reopen_artifacts writes; the
    %in-suite reopen uses an independent RFC 4180 reader and readcell.
    %
    %SCHEMA. The difference-section label, the GUI column names and the
    %headless variable names live in one place in G107ExportFixtures
    %(guiDiffLabel, directionText, guiDiffColumns, headlessDiffColumns).
    %
    %BEFORE THE FIX. Tests marked CONTROL pass before and after. Every other
    %test failed before the fix for the reason in its comment.

    properties (TestParameter)
        %Coefficient x bound x interval width. All three bounds, both
        %coefficients, and one width (50%) whose limits interpolate between
        %draws instead of landing on the extreme draw.
        Coef = struct( ...
            'dependabilityPoint80', struct('gc', 1, 'meas', 2, 'ci', 0.80), ...
            'generalizabilityLower80', struct('gc', 2, 'meas', 1, 'ci', 0.80), ...
            'dependabilityUpper50', struct('gc', 1, 'meas', 3, 'ci', 0.50));
        %G109: the three confirmed extension-normalization sites
        Site = struct( ...
            'subjectCoefficients', 'subjectCoefficients', ...
            'viewAllDifference', 'viewAllDifference', ...
            'criterionTable', 'criterionTable');
    end

    properties (SetAccess = private)
        ShadowDir
        TmpDir
    end

    methods (TestMethodSetup)
        function setupDoublesAndFolder(testCase)
            testCase.ShadowDir = G107ExportFixtures.installShadows();
            for f = {'uiputfile', 'questdlg', 'uiwait'}
                testCase.assertTrue(startsWith(which(f{1}), testCase.ShadowDir), ...
                    sprintf('the %s double must be the resolved %s', f{1}, f{1}));
            end
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);
            G107ExportFixtures.resetDialogs({}, 'decline');
            delete(findall(groot, 'Type', 'figure'));
        end
    end

    methods (TestMethodTeardown)
        function teardownDoublesAndFolder(testCase)
            delete(findall(groot, 'Type', 'figure'));
            G107ExportFixtures.removeShadows(testCase.ShadowDir);
            if isfolder(testCase.TmpDir)
                rmdir(testCase.TmpDir, 's');
            end
        end
    end

    methods (Test)

        %% ================= G107, GUI (Save Subject Coefficients) =========

        function testGuiWritesOneDifferenceSectionPerGroup(testCase, Coef)
            %(a) GUI half, for each coefficient/bound/width. Before the fix:
            %the callback writes the two event sections only, so the section
            %"Difference Score (loss minus gain)" is not found.
            pd = localMixedSummary(testCase, G107ExportFixtures.baseRel(), Coef);
            R = G107ExportFixtures.diffReference(pd, 1, Coef.gc, Coef.ci, ...
                pd.relsummary.depcutoff, Coef.meas);
            D = pd.relsummary.group(1).diffscore.ssrel_table;
            testCase.assertTrue(all(D.trls1 ~= D.trls2), ...
                'weld: the two counts must differ for everyone, or trls = min() could pass for both');
            localVerifyStoredAgainstReference(testCase, D, R, 'stored difference table');

            recs = localGuiSave(testCase, pd, 'subject_coefficients.csv');
            S = localGuiSections(recs, pd);
            testCase.verifyEqual({S.label}, [G107ExportFixtures.eventSectionLabels(pd), ...
                G107ExportFixtures.diffSectionLabels(pd)], ...
                'event sections first and unchanged, then one difference section per group');
            dl = G107ExportFixtures.diffSectionLabels(pd);
            s = localSection(testCase, S, dl{1});
            localVerifyGuiDiffSection(testCase, s, pd, 1, R, dl{1});
            testCase.verifyFalse(localHasEmptyStatus([s.notes, s.trailing], ''), ...
                'a populated difference section must not carry an empty status');
        end

        function testGuiEventSectionsKeepValuesAndPerEventFlags(testCase)
            %(b) The event sections keep their labels, column row, values and
            %their OWN inclusion flags (CONTROL half), while the difference
            %section carries the difference flags (RED half). The cutoff is
            %staged between one participant's event-1 and difference values,
            %so the two flags provably differ for that participant.
            pd0 = G107ExportFixtures.baseRel();
            clean = G107ExportFixtures.summarize(pd0, 1, 0.8, 0.001, 2);
            E1 = clean.relsummary.group(1).event(1).ssrel_table;
            D0 = clean.relsummary.group(1).diffscore.ssrel_table;
            [~, ia, ib] = intersect(E1.id, D0.id, 'stable');
            [gap, j] = max(E1.dep_pt(ia) - D0.dep_pt(ib));
            testCase.assertGreaterThan(gap, 0, ...
                'weld: some participant''s event-1 value must exceed their difference value');
            cut = (E1.dep_pt(ia(j)) + D0.dep_pt(ib(j))) / 2;
            pd = G107ExportFixtures.summarize(pd0, 1, 0.8, cut, 2);
            ev = pd.relsummary.group(1).event;
            D = pd.relsummary.group(1).diffscore.ssrel_table;
            who = D0.id{ib(j)};
            testCase.assertTrue(ev(1).ssrel_table.ind2include(strcmp(ev(1).ssrel_table.id, who)) && ...
                ~D.ind2include(strcmp(D.id, who)), ...
                'weld: the staged participant must be event-included and difference-excluded');

            recs = localGuiSave(testCase, pd, 'event_sections.csv');
            S = localGuiSections(recs, pd);
            labels = G107ExportFixtures.eventSectionLabels(pd);
            coef = G107ExportFixtures.eventCoefLabel(pd);
            for e = 1:numel(labels)
                s = localSection(testCase, S, labels{e});
                localVerifyGuiEventSection(testCase, s, ev(e).ssrel_table, coef, labels{e});
            end
            R = G107ExportFixtures.diffReference(pd, 1, 1, 0.8, cut, 2);
            dl = G107ExportFixtures.diffSectionLabels(pd);
            s = localSection(testCase, S, dl{1});
            localVerifyGuiDiffSection(testCase, s, pd, 1, R, dl{1});
        end

        function testGuiGroupedSectionsFollowGroupOrderAndAlignment(testCase)
            %(c) Two groups sharing IDs, sizes 3 and 2, reversed input rows,
            %a permuted id2 mapping and different draws and counts in the
            %second group. Before the fix: no difference sections.
            c = struct('gc', 1, 'meas', 2, 'ci', 0.8);
            pd = localMixedSummary(testCase, G107ExportFixtures.groupedRel(), c);
            DA = pd.relsummary.group(1).diffscore.ssrel_table;
            DB = pd.relsummary.group(2).diffscore.ssrel_table;
            testCase.assertEqual([height(DA), height(DB)], [3, 2], 'weld: unequal group sizes');
            [~, ia, ib] = intersect(DA.id, DB.id);
            testCase.assertNumElements(ia, 2, 'weld: the groups share two IDs');
            testCase.assertTrue(all(abs(DA.dep_pt(ia) - DB.dep_pt(ib)) > 1e-6), ...
                'weld: a shared ID must carry different values in the two groups');
            refs = cell(1, 2);
            for g = 1:2
                refs{g} = G107ExportFixtures.diffReference(pd, g, c.gc, c.ci, ...
                    pd.relsummary.depcutoff, c.meas);
                localVerifyStoredAgainstReference(testCase, ...
                    pd.relsummary.group(g).diffscore.ssrel_table, refs{g}, ...
                    sprintf('stored group %d table', g));
            end

            recs = localGuiSave(testCase, pd, 'grouped.csv');
            S = localGuiSections(recs, pd);
            elabels = G107ExportFixtures.eventSectionLabels(pd);
            dlabels = G107ExportFixtures.diffSectionLabels(pd);
            testCase.verifyEqual({S.label}, [elabels, dlabels], ...
                'group-major event sections unchanged, then GrpAlpha''s then GrpBeta''s difference section');
            coef = G107ExportFixtures.eventCoefLabel(pd);
            k = 0;
            for g = 1:2
                for e = 1:2
                    k = k + 1;
                    s = localSection(testCase, S, elabels{k});
                    localVerifyGuiEventSection(testCase, s, ...
                        pd.relsummary.group(g).event(e).ssrel_table, coef, elabels{k});
                end
            end
            for g = 1:2
                s = localSection(testCase, S, dlabels{g});
                localVerifyGuiDiffSection(testCase, s, pd, g, refs{g}, dlabels{g});
            end
        end

        function testGuiReversedEventOrderKeepsContrastAndCountsAligned(testCase)
            %(c) Reversed event order. The coefficient is symmetric in the two
            %events, so every value must equal the base run while Event A,
            %Event B, the direction and trls1/trls2 follow the reversed
            %order. Before the fix: no difference section.
            base = G107ExportFixtures.summarize(G107ExportFixtures.baseRel(), 1, 0.8, 0.5, 2);
            pd = G107ExportFixtures.summarize(G107ExportFixtures.reversedEventRel(), 1, 0.8, 0.5, 2);
            Db = base.relsummary.group(1).diffscore.ssrel_table;
            Dr = pd.relsummary.group(1).diffscore.ssrel_table;
            testCase.assertEqual(G107ExportFixtures.eventNames(pd), {'gain', 'loss'}, ...
                'weld: the reversed fixture must list gain first');
            testCase.assertEqual(Dr.id, Db.id, 'weld: same roster');
            testCase.assertEqual(Dr.trls1, Db.trls2, 'weld: the counts swap with the events');
            testCase.assertEqual(Dr.dep_pt, Db.dep_pt, 'AbsTol', 1e-12, ...
                'weld: the difference coefficient is invariant to event order');
            R = G107ExportFixtures.diffReference(pd, 1, 1, 0.8, 0.5, 2);
            localVerifyStoredAgainstReference(testCase, Dr, R, 'stored reversed table');

            recs = localGuiSave(testCase, pd, 'reversed.csv');
            S = localGuiSections(recs, pd);
            dl = G107ExportFixtures.diffSectionLabels(pd);
            testCase.verifyEqual(dl{1}, 'Difference Score (gain minus loss)', ...
                'precondition on the proposed label');
            s = localSection(testCase, S, dl{1});
            localVerifyGuiDiffSection(testCase, s, pd, 1, R, dl{1});
        end

        function testGuiMissingOrMalformedDifferenceTableFailsBeforeAnyWrite(testCase)
            %(d) A missing or malformed difference table is refused before any
            %file opens, with a message naming the re-summary route (and the
            %group, when grouped). The planted target keeps its bytes and no
            %other file appears. The window is built from a VALID summary and
            %the malformed one is injected through the callback's own
            %arguments. Before the fix: the callback writes the event rows
            %over the planted target and says nothing.
            pd = G107ExportFixtures.summarize(G107ExportFixtures.baseRel(), 1, 0.8, 0.5, 2);
            figs = G107ExportFixtures.openWindows(testCase, ...
                @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
            button = G107ExportFixtures.findButton(testCase, figs, 'Save Subject');
            target = fullfile(testCase.TmpDir, 'planted_target.csv');
            G107ExportFixtures.plant(target, sprintf('PLANTED SENTINEL G107\n'));
            bytes0 = G107ExportFixtures.fileBytes(target);
            files0 = G107ExportFixtures.folderFiles(testCase.TmpDir);
            V = localMalformedVariants(pd);
            for k = 1:size(V, 1)
                G107ExportFixtures.resetDialogs({target}, 'decline');
                out = G107ExportFixtures.fire(button, V{k, 2});
                localVerifyRefusal(testCase, out.errorText, out.dialogText, '', ...
                    ['GUI, ' V{k, 1}]);
                testCase.verifyEqual(G107ExportFixtures.fileBytes(target), bytes0, ...
                    ['GUI, ' V{k, 1} ': the planted target must keep its bytes']);
                testCase.verifyEqual(G107ExportFixtures.folderFiles(testCase.TmpDir), files0, ...
                    ['GUI, ' V{k, 1} ': no file may be written']);
                localCloseDialogs(out);
            end

            %grouped: only the second group is malformed; every group is
            %validated before anything is written
            pdG = G107ExportFixtures.summarize(G107ExportFixtures.groupedRel(), 1, 0.8, 0.5, 2);
            figsG = G107ExportFixtures.openWindows(testCase, ...
                @() psyrat_variancet('psyrat_data', pdG, 'gui', 1));
            buttonG = G107ExportFixtures.findButton(testCase, figsG, 'Save Subject');
            bad = pdG;
            bad.relsummary.group(2).diffscore = ...
                rmfield(bad.relsummary.group(2).diffscore, 'ssrel_table');
            G107ExportFixtures.resetDialogs({target}, 'decline');
            out = G107ExportFixtures.fire(buttonG, bad);
            localVerifyRefusal(testCase, out.errorText, out.dialogText, 'GrpBeta', ...
                'GUI, second group malformed');
            testCase.verifyEqual(G107ExportFixtures.fileBytes(target), bytes0, ...
                'GUI grouped: the planted target must keep its bytes');
            testCase.verifyEqual(G107ExportFixtures.folderFiles(testCase.TmpDir), files0, ...
                'GUI grouped: no file may be written');
        end

        function testGuiCancelPrecedesThePreflight(testCase)
            %(d) CONTROL. Cancel stays harmless even when the summary would
            %fail the preflight: nothing written, nothing said, no error.
            pd = G107ExportFixtures.summarize(G107ExportFixtures.baseRel(), 1, 0.8, 0.5, 2);
            figs = G107ExportFixtures.openWindows(testCase, ...
                @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
            button = G107ExportFixtures.findButton(testCase, figs, 'Save Subject');
            bad = pd;
            bad.relsummary.group = rmfield(bad.relsummary.group, 'diffscore');
            files0 = G107ExportFixtures.folderFiles(testCase.TmpDir);
            G107ExportFixtures.resetDialogs({0}, 'decline');
            out = G107ExportFixtures.fire(button, bad);
            testCase.verifyEmpty(out.errorText, 'Cancel must not throw');
            testCase.verifyEmpty(strtrim(out.dialogText), 'Cancel must open no dialog');
            testCase.verifyEqual(G107ExportFixtures.folderFiles(testCase.TmpDir), files0, ...
                'Cancel must write nothing');
            st = G107ExportFixtures.dialogState();
            testCase.verifyNumElements(st.uiputfileCalls, 1, 'one save dialog, then nothing');
            testCase.verifyEmpty(st.questdlgCalls, 'Cancel must not lead to a question');
        end

        function testGuiTypedEmptyDifferenceTableWritesEmptyStatus(testCase)
            %(e) A typed zero-row table (injected; a disjoint-roster zero-row
            %summary is not claimed reachable) writes its labeled section with
            %an explicit empty status and no rows; a populated section carries
            %no such status. Before the fix: no difference section at all.
            pd = G107ExportFixtures.summarize(G107ExportFixtures.baseRel(), 1, 0.8, 0.5, 2);
            D = pd.relsummary.group(1).diffscore.ssrel_table;
            bad = pd;
            bad.relsummary.group(1).diffscore.ssrel_table = D([], :);
            testCase.assertEqual(bad.relsummary.group(1).diffscore.ssrel_table.Properties.VariableNames, ...
                D.Properties.VariableNames, 'weld: the injected table keeps the schema');
            recs = localGuiSave(testCase, pd, 'typed_empty.csv', bad);
            S = localGuiSections(recs, pd);
            dl = G107ExportFixtures.diffSectionLabels(pd);
            s = localSection(testCase, S, dl{1});
            testCase.verifyEmpty(s.data, 'a zero-row table writes no rows and fabricates none');
            testCase.verifyTrue(localHasEmptyStatus([s.notes, s.trailing], ''), ...
                'the difference section must state its empty status explicitly');
            labels = G107ExportFixtures.eventSectionLabels(pd);
            for e = 1:numel(labels)
                se = localSection(testCase, S, labels{e});
                testCase.verifyNumElements(se.data, 3, [labels{e} ': the event rows stay']);
            end

            %grouped: the first group is empty, the second is populated
            pdG = G107ExportFixtures.summarize(G107ExportFixtures.groupedRel(), 1, 0.8, 0.5, 2);
            DG = pdG.relsummary.group(1).diffscore.ssrel_table;
            badG = pdG;
            badG.relsummary.group(1).diffscore.ssrel_table = DG([], :);
            recs = localGuiSave(testCase, pdG, 'typed_empty_grouped.csv', badG);
            S = localGuiSections(recs, pdG);
            dl = G107ExportFixtures.diffSectionLabels(pdG);
            sA = localSection(testCase, S, dl{1});
            sB = localSection(testCase, S, dl{2});
            testCase.verifyEmpty(sA.data, 'GrpAlpha: no rows');
            testCase.verifyTrue(localHasEmptyStatus([sA.notes, sA.trailing], ''), ...
                'GrpAlpha: explicit empty status');
            RB = G107ExportFixtures.diffReference(pdG, 2, 1, 0.8, 0.5, 2);
            localVerifyGuiDiffSection(testCase, sB, pdG, 2, RB, dl{2});
            testCase.verifyFalse(localHasEmptyStatus([sB.notes, sB.trailing], ''), ...
                'GrpBeta: a populated section carries no empty status');
        end

        function testGuiZeroRetentionKeepsEveryRowWithFlagZero(testCase)
            %(f) All excluded is not empty: every row stays, each flag 0,
            %every value unchanged. Before the fix: no difference section.
            pd0 = G107ExportFixtures.baseRel();
            clean = G107ExportFixtures.summarize(pd0, 1, 0.8, 0.001, 2);
            top = max(clean.relsummary.group(1).diffscore.ssrel_table.dep_pt);
            cut = top + (1 - top) / 2;
            testCase.assertTrue(cut > top && cut < 1, 'weld: the cutoff clears everyone and stays in (0,1)');
            pd = G107ExportFixtures.summarize(pd0, 1, 0.8, cut, 2);
            D = pd.relsummary.group(1).diffscore.ssrel_table;
            testCase.assertEqual(height(D), 3, 'weld: three participants');
            testCase.assertFalse(any(D.ind2include), 'weld: nobody is retained');
            R = G107ExportFixtures.diffReference(pd, 1, 1, 0.8, cut, 2);
            localVerifyStoredAgainstReference(testCase, D, R, 'stored zero-retention table');

            recs = localGuiSave(testCase, pd, 'zero_retention.csv');
            S = localGuiSections(recs, pd);
            dl = G107ExportFixtures.diffSectionLabels(pd);
            s = localSection(testCase, S, dl{1});
            localVerifyGuiDiffSection(testCase, s, pd, 1, R, dl{1});
            testCase.verifyFalse(localHasEmptyStatus([s.notes, s.trailing], ''), ...
                'zero retention is not an empty table');
        end

        function testGuiAllNaNGroupKeepsRowsAndFlags(testCase)
            %(f) A group whose per-participant residual draws are all NaN (a
            %failed fit) keeps its rows, with flags 0 and NaN values written
            %as NaN, beside a normal group. Before the fix: no difference
            %sections.
            pd0 = G107ExportFixtures.groupedRel();
            pd0.rel.out.er_var_ss1{1} = nan(size(pd0.rel.out.er_var_ss1{1}));
            pd0.rel.out.er_var_ss2{1} = nan(size(pd0.rel.out.er_var_ss2{1}));
            pd = G107ExportFixtures.summarize(pd0, 1, 0.8, 0.5, 2);
            DA = pd.relsummary.group(1).diffscore.ssrel_table;
            testCase.assertEqual(height(DA), 3, 'weld: GrpAlpha keeps three rows');
            testCase.assertTrue(all(isnan(DA.dep_pt)), 'weld: GrpAlpha''s values are NaN');
            testCase.assertFalse(any(DA.ind2include), 'weld: NaN never clears the cutoff');
            refs = {G107ExportFixtures.diffReference(pd, 1, 1, 0.8, 0.5, 2), ...
                G107ExportFixtures.diffReference(pd, 2, 1, 0.8, 0.5, 2)};

            recs = localGuiSave(testCase, pd, 'all_nan.csv');
            S = localGuiSections(recs, pd);
            dl = G107ExportFixtures.diffSectionLabels(pd);
            for g = 1:2
                s = localSection(testCase, S, dl{g});
                localVerifyGuiDiffSection(testCase, s, pd, g, refs{g}, dl{g});
                testCase.verifyFalse(localHasEmptyStatus([s.notes, s.trailing], ''), ...
                    [dl{g} ': all-NaN is not an empty table']);
            end
        end

        function testGuiDifferenceProvenanceNotesAndErrorVarianceLabel(testCase)
            %(h) The file records the run provenance (CONTROL: already
            %written today), the admissibility note (injected at the summary
            %seam, with a comma so it must be quoted), a note that the
            %difference error variance is the relative error variance of the
            %averaged contrast at the participant's own counts, the derived
            %(not Rocha Table 6) ICC disclosure, and a difference error
            %variance column that is not called 'Error Variance'. The notes
            %must also say that the SEM uses the selected coefficient's error
            %and that the difference error variance is always relative. Before
            %the fix: no admissibility note, no notes, no section.
            c = struct('gc', 2, 'meas', 2, 'ci', 0.8);
            pd = localMixedSummary(testCase, G107ExportFixtures.baseRel(), c);
            pd.relsummary.admissibility_note = localInjectedNote();
            recs = localGuiSave(testCase, pd, 'provenance.csv');
            S = localGuiSections(recs, pd);
            first = numel(recs) + 1;
            if ~isempty(S)
                first = min([S.index]);
            end
            head = G107ExportFixtures.singles(recs(1:first - 1));
            prov = psyrat_provenance_lines(pd);
            for k = 1:numel(prov)
                testCase.verifyTrue(any(strcmp(head, prov{k})), sprintf( ...
                    'provenance line missing or split across fields: %s', prov{k}));
            end
            lines = G107ExportFixtures.singles(recs);
            testCase.verifyTrue(any(strcmp(lines, pd.relsummary.admissibility_note)), ...
                'the admissibility note must be written as one field');
            testCase.verifyTrue(localAnyLine(lines, {'relative error variance', 'own', 'always'}), ...
                'a note must say the difference error variance is always the relative error variance at the participant''s own counts');
            testCase.verifyTrue(localAnyLine(lines, {'SEM', 'absolute error', 'relative error'}), ...
                'a note must say the SEM uses the selected coefficient''s error (absolute for dependability, relative for generalizability)');
            testCase.verifyTrue(localAnyLine(lines, {'derived', 'table 6'}), ...
                'a note must keep the derived (not Rocha Table 6) ICC disclosure');

            labels = G107ExportFixtures.eventSectionLabels(pd);
            se = localSection(testCase, S, labels{1});
            testCase.verifyTrue(any(strcmp(se.header, 'Error Variance')), ...
                'CONTROL: the event column keeps its label');
            dl = G107ExportFixtures.diffSectionLabels(pd);
            sd = localSection(testCase, S, dl{1});
            testCase.verifyFalse(any(strcmp(sd.header, 'Error Variance')), ...
                'the difference error variance must not share the event label');
            [keys, names] = G107ExportFixtures.guiDiffColumns();
            testCase.verifyTrue(any(strcmp(sd.header, names{strcmp(keys, 'errvar')})), ...
                'the difference error variance carries its own label');
        end

        function testGuiSpecialIdsRoundTripInEverySection(testCase)
            %G109 escaping. IDs with a comma, a double quote and an embedded
            %line feed must come back exactly, one field each, in the event
            %sections and the difference section. Before the fix: the comma
            %ID splits into 14 fields and the line-feed ID into two records.
            c = struct('gc', 1, 'meas', 2, 'ci', 0.8);
            pd = localMixedSummary(testCase, G107ExportFixtures.specialIdRel(), c);
            [~, newIds] = G107ExportFixtures.specialIdMap();
            D = pd.relsummary.group(1).diffscore.ssrel_table;
            testCase.assertTrue(all(ismember(newIds, D.id)), ...
                'weld: every special ID survives the summary');
            R = G107ExportFixtures.diffReference(pd, 1, 1, 0.8, pd.relsummary.depcutoff, 2);
            localVerifyStoredAgainstReference(testCase, D, R, 'stored special-ID table');

            recs = localGuiSave(testCase, pd, 'special_ids.csv');
            S = localGuiSections(recs, pd);
            labels = G107ExportFixtures.eventSectionLabels(pd);
            coef = G107ExportFixtures.eventCoefLabel(pd);
            for e = 1:numel(labels)
                s = localSection(testCase, S, labels{e});
                localVerifyGuiEventSection(testCase, s, ...
                    pd.relsummary.group(1).event(e).ssrel_table, coef, labels{e});
            end
            dl = G107ExportFixtures.diffSectionLabels(pd);
            s = localSection(testCase, S, dl{1});
            localVerifyGuiDiffSection(testCase, s, pd, 1, R, dl{1});
        end

        function testGuiSpecialGroupNamesAreEscaped(testCase)
            %G109 escaping of the other text fields: group names with a comma
            %and a double quote reach the difference rows' Group field and
            %the section labels intact. The event labels are a CONTROL (B19
            %quoting already covers them). Before the fix: no difference
            %sections.
            c = struct('gc', 1, 'meas', 2, 'ci', 0.8);
            pd = localMixedSummary(testCase, G107ExportFixtures.specialGroupRel(), c);
            recs = localGuiSave(testCase, pd, 'special_groups.csv');
            S = localGuiSections(recs, pd);
            testCase.verifyEqual({S.label}, [G107ExportFixtures.eventSectionLabels(pd), ...
                G107ExportFixtures.diffSectionLabels(pd)], 'section labels round-trip');
            dl = G107ExportFixtures.diffSectionLabels(pd);
            for g = 1:2
                R = G107ExportFixtures.diffReference(pd, g, 1, 0.8, pd.relsummary.depcutoff, 2);
                s = localSection(testCase, S, dl{g});
                localVerifyGuiDiffSection(testCase, s, pd, g, R, dl{g});
            end
        end

        %% ================= G107, headless (psyrat_report) =================

        function testReportAddsDifferenceTableBesideUnchangedEventTable(testCase, Coef)
            %(a) headless half and the tables.sscoeffs compatibility half.
            %Before the fix: report.tables has no sscoeffs_diff.
            pd = localMixedSummary(testCase, G107ExportFixtures.baseRel(), Coef);
            R = G107ExportFixtures.diffReference(pd, 1, Coef.gc, Coef.ci, ...
                pd.relsummary.depcutoff, Coef.meas);
            report = psyrat_report(pd);
            testCase.verifyTrue(isequaln(report.tables.sscoeffs, ...
                G107ExportFixtures.flattenEventTables(pd)), ...
                'CONTROL: tables.sscoeffs must stay exactly the flattened event tables');
            testCase.assertTrue(isfield(report.tables, 'sscoeffs_diff'), ...
                'psyrat_report must return tables.sscoeffs_diff for analysis 8');
            testCase.verifyEqual(fieldnames(report.tables)', ...
                {'variance', 'sscoeffs', 'sscoeffs_diff'}, ...
                'sscoeffs_diff comes after sscoeffs (TestHeadlessReport finds the event file first)');
            localVerifyHeadlessDiffTable(testCase, report.tables.sscoeffs_diff, pd, {R}, 'tables.sscoeffs_diff');
        end

        function testReportDifferenceFileReopensInBothFormats(testCase)
            %(g), (h) The difference table is written as its own CSV and XLSX
            %file; reopened with the independent reader (CSV) and readcell
            %(XLSX) it holds the stored and raw-draw values, and its header
            %carries every provenance line, the admissibility note, and the
            %two notes. The event file is still written with its six rows.
            %Before the fix: no _sscoeffs_diff file.
            c = struct('gc', 2, 'meas', 1, 'ci', 0.8);
            pd = localMixedSummary(testCase, G107ExportFixtures.baseRel(), c);
            pd.relsummary.admissibility_note = localInjectedNote();
            pd.proc.savename = 'g107_report.psyrat';
            R = G107ExportFixtures.diffReference(pd, 1, c.gc, c.ci, pd.relsummary.depcutoff, c.meas);
            for fmt = {'.csv', '.xlsx'}
                outdir = fullfile(testCase.TmpDir, fmt{1}(2:end));
                report = psyrat_report(pd, 'outdir', outdir, 'format', fmt{1});
                want = fullfile(outdir, ['g107_report_sscoeffs_diff' fmt{1}]);
                testCase.verifyTrue(any(strcmp(report.files, want)), ...
                    [fmt{1} ': report.files must list the difference file']);
                testCase.assertTrue(isfile(want), [fmt{1} ': the difference file must exist']);
                names = G107ExportFixtures.headlessDiffColumns();
                [headText, vals] = localReadHeadless(testCase, want, names);
                valueAt = @(k, n) vals{k, find(strcmp(names, n), 1)};
                localVerifyHeadlessRows(testCase, size(vals, 1), valueAt, pd, {R}, 1e-12, ...
                    [fmt{1} ' difference file']);
                localVerifyHeaderText(testCase, headText, pd, [fmt{1} ' difference file']);
                evfile = fullfile(outdir, ['g107_report_sscoeffs' fmt{1}]);
                testCase.verifyTrue(isfile(evfile), [fmt{1} ': the event file is still written']);
                if isfile(evfile)
                    [~, evals] = localReadHeadless(testCase, evfile, ...
                        G107ExportFixtures.headlessEventColumns());
                    testCase.verifyEqual(size(evals, 1), 6, [fmt{1} ': six event rows']);
                end
            end
        end

        function testReportGroupedAndReversedEventOrderAlignment(testCase)
            %(c) headless: grouped rows in group order and stored order, each
            %against its own group's reference; reversed event order keeps
            %the values and swaps the contrast labels and counts. Before the
            %fix: no sscoeffs_diff.
            c = struct('gc', 1, 'meas', 2, 'ci', 0.8);
            pd = localMixedSummary(testCase, G107ExportFixtures.groupedRel(), c);
            refs = {G107ExportFixtures.diffReference(pd, 1, 1, 0.8, pd.relsummary.depcutoff, 2), ...
                G107ExportFixtures.diffReference(pd, 2, 1, 0.8, pd.relsummary.depcutoff, 2)};
            report = psyrat_report(pd);
            testCase.assertTrue(isfield(report.tables, 'sscoeffs_diff'), 'grouped: sscoeffs_diff missing');
            localVerifyHeadlessDiffTable(testCase, report.tables.sscoeffs_diff, pd, refs, 'grouped');

            base = G107ExportFixtures.summarize(G107ExportFixtures.baseRel(), 1, 0.8, 0.5, 2);
            rev = G107ExportFixtures.summarize(G107ExportFixtures.reversedEventRel(), 1, 0.8, 0.5, 2);
            rb = psyrat_report(base);
            rr = psyrat_report(rev);
            testCase.assertTrue(isfield(rr.tables, 'sscoeffs_diff') && isfield(rb.tables, 'sscoeffs_diff'), ...
                'reversed: sscoeffs_diff missing');
            Rr = G107ExportFixtures.diffReference(rev, 1, 1, 0.8, 0.5, 2);
            localVerifyHeadlessDiffTable(testCase, rr.tables.sscoeffs_diff, rev, {Rr}, 'reversed');
            Tb = rb.tables.sscoeffs_diff;
            Tr = rr.tables.sscoeffs_diff;
            testCase.verifyEqual(Tr.dep_pt, Tb.dep_pt, 'AbsTol', 1e-12, 'reversed: values invariant');
            testCase.verifyEqual([Tr.trls1, Tr.trls2], [Tb.trls2, Tb.trls1], 'reversed: counts swap');
            testCase.verifyEqual(string(Tr.EventA), repmat("gain", height(Tr), 1), 'reversed: EventA');
        end

        function testReportMissingOrMalformedDifferenceTableFailsBeforeAnyFile(testCase)
            %(d) headless: the same variants raise an error naming the
            %re-summary route (and the group) with no outdir, into a fresh
            %outdir (nothing created), and into a folder of planted targets
            %under an explicit savename (every byte kept). Before the fix:
            %the report succeeds and overwrites the planted files.
            pd = G107ExportFixtures.summarize(G107ExportFixtures.baseRel(), 1, 0.8, 0.5, 2);
            V = localMalformedVariants(pd);
            pdG = G107ExportFixtures.summarize(G107ExportFixtures.groupedRel(), 1, 0.8, 0.5, 2);
            badG = pdG;
            badG.relsummary.group(2).diffscore = ...
                rmfield(badG.relsummary.group(2).diffscore, 'ssrel_table');
            V(end+1, :) = {'second group malformed', badG};
            for k = 1:size(V, 1)
                what = ['report, ' V{k, 1}];
                gname = '';
                if strcmp(V{k, 1}, 'second group malformed')
                    gname = 'GrpBeta';
                end
                bad = V{k, 2};
                msg = localError(@() psyrat_report(bad));
                localVerifyRefusal(testCase, msg, '', gname, [what ' (no outdir)']);

                fresh = fullfile(testCase.TmpDir, sprintf('fresh_%d', k));
                msg = localError(@() psyrat_report(bad, 'outdir', fresh, 'format', '.csv'));
                localVerifyRefusal(testCase, msg, '', gname, [what ' (fresh outdir)']);
                testCase.verifyTrue(~isfolder(fresh) || ...
                    isempty(G107ExportFixtures.folderFiles(fresh)), [what ': no file may be written']);

                pdir = fullfile(testCase.TmpDir, sprintf('planted_%d', k));
                mkdir(pdir);
                bad.proc.savename = 'g107_planted.psyrat';
                planted = fullfile(pdir, strcat('g107_planted_', ...
                    {'variance', 'sscoeffs', 'sscoeffs_diff'}, '.csv'));
                bytes0 = cell(size(planted));
                for p = 1:numel(planted)
                    G107ExportFixtures.plant(planted{p}, sprintf('PLANTED %d\n', p));
                    bytes0{p} = G107ExportFixtures.fileBytes(planted{p});
                end
                msg = localError(@() psyrat_report(bad, 'outdir', pdir, 'format', '.csv'));
                localVerifyRefusal(testCase, msg, '', gname, [what ' (planted targets)']);
                for p = 1:numel(planted)
                    testCase.verifyEqual(G107ExportFixtures.fileBytes(planted{p}), bytes0{p}, ...
                        sprintf('%s: %s must keep its bytes', what, planted{p}));
                end
                testCase.verifyNumElements(G107ExportFixtures.folderFiles(pdir), 3, ...
                    [what ': no other file may appear']);
            end
        end

        function testReportTypedEmptyDifferenceTableWritesEmptyStatus(testCase)
            %(e) headless: a typed zero-row table yields a zero-row
            %sscoeffs_diff with the full schema and a file whose header states
            %the empty status (naming the group when grouped) and holds no
            %rows. Before the fix: no sscoeffs_diff.
            pd = G107ExportFixtures.summarize(G107ExportFixtures.baseRel(), 1, 0.8, 0.5, 2);
            D = pd.relsummary.group(1).diffscore.ssrel_table;
            bad = pd;
            bad.relsummary.group(1).diffscore.ssrel_table = D([], :);
            bad.proc.savename = 'g107_empty.psyrat';
            outdir = fullfile(testCase.TmpDir, 'empty');
            report = psyrat_report(bad, 'outdir', outdir, 'format', '.csv');
            names = G107ExportFixtures.headlessDiffColumns();
            testCase.assertTrue(isfield(report.tables, 'sscoeffs_diff'), 'sscoeffs_diff missing');
            testCase.verifyEqual(height(report.tables.sscoeffs_diff), 0, 'no rows fabricated');
            testCase.verifyEqual(report.tables.sscoeffs_diff.Properties.VariableNames, names, ...
                'the empty table keeps the schema');
            testCase.verifyEqual(height(report.tables.sscoeffs), 6, 'the event rows stay');
            f = fullfile(outdir, 'g107_empty_sscoeffs_diff.csv');
            testCase.assertTrue(isfile(f), 'the empty difference file is still written');
            [headText, vals] = localReadHeadless(testCase, f, names);
            testCase.verifyEqual(size(vals, 1), 0, 'the file holds no rows');
            testCase.verifyTrue(localHasEmptyStatus(strsplit(headText, newline), ''), ...
                'the header must state the empty status');

            pdG = G107ExportFixtures.summarize(G107ExportFixtures.groupedRel(), 1, 0.8, 0.5, 2);
            DG = pdG.relsummary.group(1).diffscore.ssrel_table;
            badG = pdG;
            badG.relsummary.group(1).diffscore.ssrel_table = DG([], :);
            badG.proc.savename = 'g107_empty_grouped.psyrat';
            outdirG = fullfile(testCase.TmpDir, 'empty_grouped');
            reportG = psyrat_report(badG, 'outdir', outdirG, 'format', '.csv');
            testCase.assertTrue(isfield(reportG.tables, 'sscoeffs_diff'), 'grouped: sscoeffs_diff missing');
            T = reportG.tables.sscoeffs_diff;
            testCase.verifyEqual(height(T), 2, 'grouped: only GrpBeta''s two rows');
            testCase.verifyTrue(all(string(T.Group) == "GrpBeta"), 'grouped: rows belong to GrpBeta');
            [headText, ~] = localReadHeadless(testCase, ...
                fullfile(outdirG, 'g107_empty_grouped_sscoeffs_diff.csv'), names);
            testCase.verifyTrue(localHasEmptyStatus(strsplit(headText, newline), 'GrpAlpha'), ...
                'grouped: the header must name the empty group');
        end

        function testReportTypedEmptyXlsxKeepsHeaderAboveColumnRow(testCase)
            %(e), (g) The typed zero-row difference table through the real
            %psyrat_report in XLSX, beside its CSV sibling. Row 1 must hold
            %the first header line ('Table Generated on ...'), the header
            %above the column row must hold the CSV header's lines in the same
            %order (the timestamp line aside), it must carry the provenance,
            %the notes and the empty status, and no row may follow the column
            %row. The grouped run (the first group empty, the second
            %populated) is the CONTROL through the ordinary nonempty write:
            %GrpBeta's rows follow the column row with their values, and the
            %header names GrpAlpha as empty. Before the fix: writetable
            %ignored 'Range' for the zero-row table and wrote the column names
            %over row 1, so the header lines fell below the column row, where
            %a reader takes them for data rows (observed by the independent
            %openpyxl reopen: 19 rows instead of 0).
            names = G107ExportFixtures.headlessDiffColumns();
            pd = G107ExportFixtures.summarize(G107ExportFixtures.baseRel(), 1, 0.8, 0.5, 2);
            bad = pd;
            bad.relsummary.group(1).diffscore.ssrel_table = ...
                pd.relsummary.group(1).diffscore.ssrel_table([], :);
            bad.proc.savename = 'g107_xempty.psyrat';
            pdG = G107ExportFixtures.summarize(G107ExportFixtures.groupedRel(), 1, 0.8, 0.5, 2);
            badG = pdG;
            badG.relsummary.group(1).diffscore.ssrel_table = ...
                pdG.relsummary.group(1).diffscore.ssrel_table([], :);
            badG.proc.savename = 'g107_xempty_grouped.psyrat';
            RB = G107ExportFixtures.diffReference(pdG, 2, 1, 0.8, 0.5, 2);
            %{what, data, group named by the empty status, data rows, refs}
            V = {'zero rows', bad, '', 0, {}; ...
                'CONTROL grouped, one group empty', badG, 'GrpAlpha', 2, {[], RB}};
            for k = 1:size(V, 1)
                what = V{k, 1};
                p = V{k, 2};
                base = erase(p.proc.savename, '.psyrat');
                dc = fullfile(testCase.TmpDir, sprintf('xempty_csv_%d', k));
                dx = fullfile(testCase.TmpDir, sprintf('xempty_xlsx_%d', k));
                psyrat_report(p, 'outdir', dc, 'format', '.csv');
                psyrat_report(p, 'outdir', dx, 'format', '.xlsx');
                fc = fullfile(dc, [base '_sscoeffs_diff.csv']);
                fx = fullfile(dx, [base '_sscoeffs_diff.xlsx']);
                testCase.assertTrue(isfile(fc) && isfile(fx), [what ': both difference files are written']);

                recs = G107ExportFixtures.parseCsv(fileread(fc));
                hc = find(cellfun(@(r) isequal(r, names), recs), 1);
                testCase.assertNotEmpty(hc, [what ': the CSV column row']);
                csvHead = localHeaderTexts(G107ExportFixtures.singles(recs(1:hc - 1)));

                raw = readcell(fx);
                hx = localColumnRow(raw, names);
                testCase.assertNotEmpty(hx, [what ': the XLSX column row']);
                testCase.verifyTrue(startsWith(localText(raw(1, 1)), 'Table Generated on'), ...
                    [what ': row 1 holds the first header line, not the column row']);
                xHead = localHeaderTexts(raw(1:hx - 1, 1)');
                testCase.verifyEqual(xHead, csvHead, [what ': the XLSX header above the ' ...
                    'column row holds the CSV header''s lines, in order']);
                localVerifyHeaderText(testCase, strjoin(xHead, newline), p, [what ', XLSX']);
                testCase.verifyTrue(localHasEmptyStatus(xHead, V{k, 3}), ...
                    [what ': the XLSX header states the empty status']);
                body = raw(hx + 1:end, :);
                filled = any(~cellfun(@(x) isempty(localText({x})), body), 2);
                testCase.verifyEqual(nnz(filled), V{k, 4}, ...
                    [what ': the rows below the column row are the table''s rows only']);
                if V{k, 4} > 0
                    [~, vals] = localReadHeadless(testCase, fx, names);
                    valueAt = @(r, n) vals{r, find(strcmp(names, n), 1)};
                    localVerifyHeadlessRows(testCase, size(vals, 1), valueAt, p, V{k, 5}, ...
                        1e-12, [what ', XLSX rows']);
                end
            end
        end

        function testWriteTableZeroRowXlsxLayout(testCase)
            %psyrat_write_table on its own. A zero-row table below a
            %three-line header: the header stays in rows 1-3, the column-name
            %row is row 4 with the variable names unchanged, and nothing
            %follows it. Before the fix: the names landed in row 1 over the
            %first header line. CONTROLS (expected to pass before and after;
            %not run against the old code): a zero-row table with no header
            %(the names in row 1 and nothing else), a two-row table below the
            %same header (the ordinary layout), and the CSV writer's zero-row
            %layout, which already put the names below the header.
            T = table({'r1'; 'r2'}, [1; 2], 'VariableNames', {'Label', 'Value'});
            head = {'Header line 1'; 'Header line 2, with a comma'; 'Header line 3'};
            %{what, header, table, data rows}
            C = { ...
                'zero rows below a header', head, T([], :), cell(0, 2); ...
                'CONTROL zero rows, no header', {}, T([], :), cell(0, 2); ...
                'CONTROL two rows below a header', head, T, {'r1', 1; 'r2', 2}};
            for k = 1:size(C, 1)
                what = C{k, 1};
                f = fullfile(testCase.TmpDir, sprintf('write_table_%d.xlsx', k));
                psyrat_write_table(f, C{k, 2}, C{k, 3});
                testCase.assertTrue(isfile(f), [what ': the workbook is written']);
                raw = readcell(f);
                nh = numel(C{k, 2});
                nd = size(C{k, 4}, 1);
                testCase.verifyEqual(size(raw, 1), nh + 1 + nd, ...
                    [what ': the header rows, one column-name row and the data rows, nothing else']);
                if size(raw, 1) ~= nh + 1 + nd || size(raw, 2) < 2
                    continue;
                end
                for r = 1:nh
                    testCase.verifyEqual(localText(raw(r, 1)), C{k, 2}{r}, ...
                        sprintf('%s: header line %d stays in row %d', what, r, r));
                end
                testCase.verifyEqual({localText(raw(nh + 1, 1)), localText(raw(nh + 1, 2))}, ...
                    {'Label', 'Value'}, [what ': the column-name row directly below the header']);
                for r = 1:nd
                    testCase.verifyEqual(localText(raw(nh + 1 + r, 1)), C{k, 4}{r, 1}, ...
                        sprintf('%s: data row %d label', what, r));
                    testCase.verifyEqual(G107ExportFixtures.parseNum(raw(nh + 1 + r, 2)), ...
                        C{k, 4}{r, 2}, sprintf('%s: data row %d value', what, r));
                end
            end

            f = fullfile(testCase.TmpDir, 'write_table_zero.csv');
            psyrat_write_table(f, head, T([], :));
            testCase.verifyEqual(G107ExportFixtures.parseCsv(fileread(f)), ...
                {{'Header line 1'}, {'Header line 2, with a comma'}, {'Header line 3'}, ...
                {'Label', 'Value'}}, 'CONTROL CSV: the header, then the column-name row, and no row after it');
        end

        function testReportZeroRetentionAndAllNaNKeepRows(testCase)
            %(f) headless: all-excluded and all-NaN keep every row with flag 0
            %and no empty status. Before the fix: no sscoeffs_diff.
            pd0 = G107ExportFixtures.baseRel();
            clean = G107ExportFixtures.summarize(pd0, 1, 0.8, 0.001, 2);
            top = max(clean.relsummary.group(1).diffscore.ssrel_table.dep_pt);
            cut = top + (1 - top) / 2;
            pd = G107ExportFixtures.summarize(pd0, 1, 0.8, cut, 2);
            testCase.assertFalse(any(pd.relsummary.group(1).diffscore.ssrel_table.ind2include), ...
                'weld: nobody is retained');
            pd.proc.savename = 'g107_zero.psyrat';
            outdir = fullfile(testCase.TmpDir, 'zero');
            report = psyrat_report(pd, 'outdir', outdir, 'format', '.csv');
            testCase.assertTrue(isfield(report.tables, 'sscoeffs_diff'), 'zero retention: sscoeffs_diff missing');
            R = G107ExportFixtures.diffReference(pd, 1, 1, 0.8, cut, 2);
            localVerifyHeadlessDiffTable(testCase, report.tables.sscoeffs_diff, pd, {R}, 'zero retention');
            [headText, ~] = localReadHeadless(testCase, ...
                fullfile(outdir, 'g107_zero_sscoeffs_diff.csv'), G107ExportFixtures.headlessDiffColumns());
            testCase.verifyFalse(localHasEmptyStatus(strsplit(headText, newline), ''), ...
                'zero retention is not an empty table');

            pdN0 = G107ExportFixtures.groupedRel();
            pdN0.rel.out.er_var_ss1{1} = nan(size(pdN0.rel.out.er_var_ss1{1}));
            pdN0.rel.out.er_var_ss2{1} = nan(size(pdN0.rel.out.er_var_ss2{1}));
            pdN = G107ExportFixtures.summarize(pdN0, 1, 0.8, 0.5, 2);
            testCase.assertTrue(all(isnan(pdN.relsummary.group(1).diffscore.ssrel_table.dep_pt)), ...
                'weld: GrpAlpha''s values are NaN');
            reportN = psyrat_report(pdN);
            testCase.assertTrue(isfield(reportN.tables, 'sscoeffs_diff'), 'all-NaN: sscoeffs_diff missing');
            refs = {G107ExportFixtures.diffReference(pdN, 1, 1, 0.8, 0.5, 2), ...
                G107ExportFixtures.diffReference(pdN, 2, 1, 0.8, 0.5, 2)};
            localVerifyHeadlessDiffTable(testCase, reportN.tables.sscoeffs_diff, pdN, refs, 'all-NaN group');
        end

        function testReportEventCompatibilityControlAndSeparateDifferenceAssertions(testCase)
            %(j) The sscoeffs assertions are TestHeadlessReport.m:222-227 on
            %the same fixture and settings (CONTROL; TestHeadlessReport is not
            %edited). The difference assertions are separate and failed before
            %the fix.
            psyrat_data = PsyRATTestDataFactory.makeICDiffSSErrRelData();
            [psyrat_data, ~] = psyrat_relsummary('psyrat_data', psyrat_data, ...
                'analysis', 'sing_diff_sserr', 'depcutoff', 0.5, 'meascutoff', 2, ...
                'depcentmeas', 1, 'diffgcoeff', 1, 'CI', 0.95);
            report = psyrat_report(psyrat_data);
            sc = report.tables.sscoeffs;
            expectedCols = {'Group','Event','id','ind2include','trls', ...
                'dep_pt','dep_ll','dep_ul','icc_pt','icc_ll','icc_ul', ...
                'sem_pt','sem_ll','sem_ul','ss_errvar'};
            testCase.verifyEqual(sc.Properties.VariableNames, expectedCols);
            testCase.verifyEqual(numel(unique(sc.Event)), 2);
            testCase.verifyEqual(height(sc), 6, 'CONTROL: six event rows');

            testCase.assertTrue(isfield(report.tables, 'sscoeffs_diff'), ...
                'the difference rows must be returned in their own table');
            sd = report.tables.sscoeffs_diff;
            testCase.verifyFalse(any(strcmp(sd.Properties.VariableNames, 'Event')), ...
                'difference rows are not dressed up as a pseudo-event');
            testCase.verifyEqual(height(sd), ...
                height(psyrat_data.relsummary.group(1).diffscore.ssrel_table), ...
                'one row per stored difference row');
            testCase.verifyTrue(all(strcmpi(string(sd.ScoreType), "difference")), ...
                'every row says it is a difference score');
        end

        function testReportSpecialIdsRoundTrip(testCase)
            %G109 escaping on the headless route: the special IDs come back
            %exactly from the event file and the difference file, CSV and
            %XLSX. The event-file half is a CONTROL (writetable quotes); the
            %difference half failed before the fix (no difference file).
            c = struct('gc', 1, 'meas', 2, 'ci', 0.8);
            pd = localMixedSummary(testCase, G107ExportFixtures.specialIdRel(), c);
            pd.proc.savename = 'g107_ids.psyrat';
            R = G107ExportFixtures.diffReference(pd, 1, 1, 0.8, pd.relsummary.depcutoff, 2);
            for fmt = {'.csv', '.xlsx'}
                outdir = fullfile(testCase.TmpDir, ['ids_' fmt{1}(2:end)]);
                psyrat_report(pd, 'outdir', outdir, 'format', fmt{1});
                evnames = G107ExportFixtures.headlessEventColumns();
                [~, evals] = localReadHeadless(testCase, ...
                    fullfile(outdir, ['g107_ids_sscoeffs' fmt{1}]), evnames);
                want = G107ExportFixtures.flattenEventTables(pd);
                got = cellfun(@localText, evals(:, strcmp(evnames, 'id')), 'UniformOutput', false);
                testCase.verifyEqual(got, cellstr(string(want.id)), ...
                    [fmt{1} ': CONTROL: event-file IDs round-trip']);
                f = fullfile(outdir, ['g107_ids_sscoeffs_diff' fmt{1}]);
                testCase.assertTrue(isfile(f), [fmt{1} ': the difference file must exist']);
                names = G107ExportFixtures.headlessDiffColumns();
                [~, vals] = localReadHeadless(testCase, f, names);
                valueAt = @(k, n) vals{k, find(strcmp(names, n), 1)};
                localVerifyHeadlessRows(testCase, size(vals, 1), valueAt, pd, {R}, 1e-12, ...
                    [fmt{1} ' special IDs']);
            end
        end

        function testReportRouteDispatchAndProvenance(testCase)
            %Deterministic dispatch/schema controls for the supported analysis
            %8 routes: CmdStan and native HMC Gaussian, gamma log-nu
            %(non-concurrent), gamma concurrent Gaussian copula (per-person
            %log-nu, experimental) and gamma location-scale. The same
            %synthetic draws are stamped with each route's provenance; the
            %difference table must exist with the schema and the stored
            %values, and its file header must carry every provenance line
            %(engine, family, parameterization, maturity). No route is
            %certified by this; the gamma extractors are not exercised. Before
            %the fix: no sscoeffs_diff.
            routes = localRoutes();
            for k = 1:size(routes, 1)
                name = routes{k, 1};
                pd = G107ExportFixtures.summarize(routes{k, 2}(G107ExportFixtures.baseRel()), ...
                    1, 0.8, 0.5, 2);
                prov = psyrat_provenance_lines(pd);
                testCase.assertTrue(any(contains(prov, routes{k, 3})), sprintf( ...
                    '%s: fixture precondition, provenance must carry "%s"', name, routes{k, 3}));
                pd.proc.savename = ['g107_route_' name '.psyrat'];
                outdir = fullfile(testCase.TmpDir, name);
                report = psyrat_report(pd, 'outdir', outdir, 'format', '.csv');
                testCase.verifyTrue(isfield(report.tables, 'sscoeffs_diff'), ...
                    [name ': sscoeffs_diff missing']);
                if ~isfield(report.tables, 'sscoeffs_diff')
                    continue;
                end
                R = G107ExportFixtures.diffReference(pd, 1, 1, 0.8, 0.5, 2);
                localVerifyHeadlessDiffTable(testCase, report.tables.sscoeffs_diff, pd, {R}, name);
                f = fullfile(outdir, ['g107_route_' name '_sscoeffs_diff.csv']);
                testCase.verifyTrue(isfile(f), [name ': difference file missing']);
                if isfile(f)
                    headText = localReadHeadless(testCase, f, G107ExportFixtures.headlessDiffColumns());
                    for j = 1:numel(prov)
                        testCase.verifyTrue(contains(headText, prov{j}), sprintf( ...
                            '%s: header must carry "%s"', name, prov{j}));
                    end
                end
            end

            %the experimental concurrent route through the GUI as well
            pdC = G107ExportFixtures.summarize(localGammaCopula(G107ExportFixtures.baseRel()), ...
                1, 0.8, 0.5, 2);
            recs = localGuiSave(testCase, pdC, 'gamma_copula.csv');
            S = localGuiSections(recs, pdC);
            first = numel(recs) + 1;
            if ~isempty(S)
                first = min([S.index]);
            end
            head = G107ExportFixtures.singles(recs(1:first - 1));
            prov = psyrat_provenance_lines(pdC);
            for j = 1:numel(prov)
                testCase.verifyTrue(any(strcmp(head, prov{j})), sprintf( ...
                    'GUI gamma copula: CONTROL, header line missing: %s', prov{j}));
            end
            dl = G107ExportFixtures.diffSectionLabels(pdC);
            s = localSection(testCase, S, dl{1});
            localVerifyGuiDiffSection(testCase, s, pdC, 1, ...
                G107ExportFixtures.diffReference(pdC, 1, 1, 0.8, 0.5, 2), 'GUI gamma copula');
        end

        %% ================= controls ======================================

        function testAnalysis6And25GuiExportsUnchanged(testCase)
            %(i) CONTROL. Analyses 6 (ic_sserrvar, two groups) and 25
            %(trt_sserrvar) have no difference score, so their Save Subject
            %Coefficients file must stay byte-identical to the shipped layout
            %after the timestamp line.
            specs = {localAnalysis6(), localAnalysis25()};
            for k = 1:numel(specs)
                pd = specs{k};
                figs = G107ExportFixtures.openWindows(testCase, ...
                    @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
                button = G107ExportFixtures.findButton(testCase, figs, 'Save Subject');
                target = fullfile(testCase.TmpDir, sprintf('unchanged_%d.csv', k));
                G107ExportFixtures.resetDialogs({target}, 'decline');
                out = G107ExportFixtures.fire(button);
                testCase.assertEmpty(out.errorText, sprintf('%s: the callback threw: %s', ...
                    pd.rel.analysis, out.errorText));
                testCase.assertTrue(isfile(target), [pd.rel.analysis ': no file written']);
                [first, rest] = localSplitAfterDate(fileread(target));
                testCase.verifyEqual(first, 'Table Generated on', [pd.rel.analysis ': first line']);
                testCase.verifyEqual(rest, localExpectedEventOnlyGuiRest(pd), ...
                    [pd.rel.analysis ': everything after the timestamp is unchanged']);
                delete(figs(isvalid(figs)));
            end
        end

        function testAnalysis6And25ReportsUnchanged(testCase)
            %(i) CONTROL, headless: the same two designs return and write
            %exactly the variance and sscoeffs tables, as today.
            specs = {localAnalysis6(), localAnalysis25()};
            for k = 1:numel(specs)
                pd = specs{k};
                report = psyrat_report(pd);
                testCase.verifyEqual(fieldnames(report.tables)', {'variance', 'sscoeffs'}, ...
                    [pd.rel.analysis ': no new table']);
                testCase.verifyTrue(isequaln(report.tables.sscoeffs, ...
                    G107ExportFixtures.flattenEventTables(pd)), [pd.rel.analysis ': sscoeffs unchanged']);
                testCase.verifyTrue(isequaln(report.tables.variance, ...
                    psyrat_variancet('psyrat_data', pd, 'gui', 0)), [pd.rel.analysis ': variance unchanged']);
                pd.proc.savename = sprintf('g107_control_%d.psyrat', k);
                outdir = fullfile(testCase.TmpDir, sprintf('control_%d', k));
                report = psyrat_report(pd, 'outdir', outdir, 'format', '.csv');
                testCase.verifyNumElements(report.files, 2, [pd.rel.analysis ': two files']);
                testCase.verifyFalse(any(contains(report.files, 'sscoeffs_diff')), ...
                    [pd.rel.analysis ': no difference file']);
            end
        end

        function testSaveIdsFieldQuotingControl(testCase)
            %G109 CONTROL. Save IDs already quotes identity fields correctly
            %(psyrat_depoverallt local_id_csv_field); on analysis 8 its lists
            %are the DIFFERENCE-based retention, and the special IDs must
            %round-trip one field each.
            c = struct('gc', 1, 'meas', 2, 'ci', 0.8);
            pd = localMixedSummary(testCase, G107ExportFixtures.specialIdRel(), c);
            g = pd.relsummary.group(1);
            testCase.assertTrue(~isempty(g.goodids) && ~isempty(g.badids), ...
                'weld: both lists are populated');
            figs = G107ExportFixtures.openWindows(testCase, ...
                @() psyrat_depoverallt('psyrat_data', pd, 'gui', 1));
            button = G107ExportFixtures.findButton(testCase, figs, 'Save IDs');
            target = fullfile(testCase.TmpDir, 'ids.csv');
            G107ExportFixtures.resetDialogs({target}, 'decline');
            out = G107ExportFixtures.fire(button);
            testCase.assertEmpty(out.errorText, ['Save IDs threw: ' out.errorText]);
            recs = G107ExportFixtures.parseCsv(fileread(target));
            h = find(cellfun(@(r) isequal(r, {'Good IDs', 'Bad IDs', 'Good Group', 'Bad Group'}), recs), 1);
            testCase.assertNotEmpty(h, 'Save IDs column row not found');
            rows = recs(h + 1:end);
            testCase.verifyTrue(all(cellfun(@numel, rows) == 4), 'every Save IDs row has four fields');
            rows = rows(cellfun(@numel, rows) == 4);
            good = cellfun(@(r) r{1}, rows, 'UniformOutput', false);
            bad = cellfun(@(r) r{2}, rows, 'UniformOutput', false);
            good = good(~cellfun(@isempty, good));
            bad = bad(~cellfun(@isempty, bad));
            testCase.verifyEqual(sort(good(:)), sort(cellstr(string(g.goodids(:)))), 'good IDs round-trip');
            testCase.verifyEqual(sort(bad(:)), sort(cellstr(string(g.badids(:)))), 'bad IDs round-trip');
        end

        function testLegacySavedResultNeedsOnlyAReSummary(testCase)
            %Existing files. A result saved WITHOUT a summary (the .psyrat
            %wrapper; the GUI loader is bypassed with load -mat) re-summarizes
            %to the identical difference table, which both exports then carry
            %(failed before the fix). A file saved with a STALE summary that predates
            %the difference table is refused with the re-summary route, and a
            %re-summary of the same file exports normally.
            pd = G107ExportFixtures.summarize(G107ExportFixtures.baseRel(), 1, 0.8, 0.5, 2);
            f = fullfile(testCase.TmpDir, 'legacy.psyrat');
            psyrat_data = rmfield(pd, 'relsummary'); %#ok<NASGU>
            save(f, 'psyrat_data', '-mat');
            z = load(f, '-mat');
            q = G107ExportFixtures.summarize(z.psyrat_data, 1, 0.8, 0.5, 2);
            testCase.assertTrue(isequaln(q.relsummary.group(1).diffscore.ssrel_table, ...
                pd.relsummary.group(1).diffscore.ssrel_table), ...
                'weld: a re-summary of the saved draws reproduces the table');
            R = G107ExportFixtures.diffReference(q, 1, 1, 0.8, 0.5, 2);
            recs = localGuiSave(testCase, q, 'legacy.csv');
            S = localGuiSections(recs, q);
            dl = G107ExportFixtures.diffSectionLabels(q);
            s = localSection(testCase, S, dl{1});
            localVerifyGuiDiffSection(testCase, s, q, 1, R, 'legacy GUI');
            report = psyrat_report(q);
            testCase.assertTrue(isfield(report.tables, 'sscoeffs_diff'), 'legacy: sscoeffs_diff missing');
            localVerifyHeadlessDiffTable(testCase, report.tables.sscoeffs_diff, q, {R}, 'legacy headless');

            stale = pd;
            stale.relsummary.group = rmfield(stale.relsummary.group, 'diffscore');
            f2 = fullfile(testCase.TmpDir, 'stale.psyrat');
            psyrat_data = stale; %#ok<NASGU>
            save(f2, 'psyrat_data', '-mat');
            z2 = load(f2, '-mat');
            msg = localError(@() psyrat_report(z2.psyrat_data));
            localVerifyRefusal(testCase, msg, '', '', 'stale saved summary');
            fixed = G107ExportFixtures.summarize(z2.psyrat_data, 1, 0.8, 0.5, 2);
            report = psyrat_report(fixed);
            testCase.verifyTrue(isfield(report.tables, 'sscoeffs_diff') && ...
                height(report.tables.sscoeffs_diff) == 3, 'a re-summary exports the three rows');
        end

        %% ================= G109, destination safety =======================

        function testDestinationCancelWritesNothing(testCase, Site)
            %CONTROL at each site: Cancel writes nothing, asks nothing, says
            %nothing.
            S = localSite(testCase, Site);
            files0 = G107ExportFixtures.folderFiles(testCase.TmpDir);
            G107ExportFixtures.resetDialogs({0}, 'decline');
            out = G107ExportFixtures.fire(S.button);
            testCase.verifyEmpty(out.errorText, [Site ': Cancel must not throw']);
            testCase.verifyEmpty(strtrim(out.dialogText), [Site ': Cancel must open no dialog']);
            testCase.verifyEqual(G107ExportFixtures.folderFiles(testCase.TmpDir), files0, ...
                [Site ': Cancel must write nothing']);
            st = G107ExportFixtures.dialogState();
            testCase.verifyNumElements(st.uiputfileCalls, 1, [Site ': one save dialog']);
            testCase.verifyEmpty(st.questdlgCalls, [Site ': no question after Cancel']);
        end

        function testDestinationFreeNormalizedNameWrites(testCase, Site)
            %CONTROL at each site: a typed name without the CSV extension
            %whose normalized name is FREE is written there, keeping the
            %extension behavior, with no second prompt.
            S = localSite(testCase, Site);
            typed = fullfile(testCase.TmpDir, ['fresh' S.typedExt]);
            final = [typed '.csv'];
            G107ExportFixtures.resetDialogs({typed}, 'decline');
            out = G107ExportFixtures.fire(S.button);
            testCase.verifyEmpty(out.errorText, [Site ': the save must not throw']);
            testCase.verifyTrue(isfile(final), [Site ': the normalized name is written']);
            testCase.verifyFalse(isfile(typed), [Site ': nothing is written under the typed name']);
            if isfile(final)
                testCase.verifyTrue(contains(fileread(final), S.marker), [Site ': real content']);
            end
            [~, fn, fe] = fileparts(final);
            testCase.verifyEqual(G107ExportFixtures.folderFiles(testCase.TmpDir), {[fn fe]}, ...
                [Site ': only the saved file remains (no temporary file is left behind)']);
            st = G107ExportFixtures.dialogState();
            testCase.verifyFalse(G107ExportFixtures.prompted(st), ...
                [Site ': a free final name needs no second prompt']);
        end

        function testDestinationCollisionDeclinedLeavesFileIntact(testCase, Site)
            %At each site: the typed name normalizes onto an EXISTING file the
            %dialog never named. Without consent to that exact file nothing
            %is written: the file keeps its bytes and no other file appears.
            %If no prompt was issued, a visible refusal must name the file.
            %Before the fix: the existing file is truncated and overwritten.
            S = localSite(testCase, Site);
            typed = fullfile(testCase.TmpDir, ['collision' S.typedExt]);
            final = [typed '.csv'];
            [~, fn, fe] = fileparts(final);
            finalName = [fn fe];
            G107ExportFixtures.plant(final, sprintf('PLANTED SENTINEL G109\n'));
            bytes0 = G107ExportFixtures.fileBytes(final);
            files0 = G107ExportFixtures.folderFiles(testCase.TmpDir);
            G107ExportFixtures.resetDialogs({typed}, 'decline');
            out = G107ExportFixtures.fire(S.button);
            testCase.verifyEmpty(out.errorText, [Site ': declining must not throw']);
            testCase.verifyEqual(G107ExportFixtures.fileBytes(final), bytes0, ...
                [Site ': the existing ' finalName ' must keep its bytes without consent']);
            testCase.verifyEqual(G107ExportFixtures.folderFiles(testCase.TmpDir), files0, ...
                [Site ': no other file may be written without consent']);
            st = G107ExportFixtures.dialogState();
            if G107ExportFixtures.prompted(st)
                testCase.verifyTrue(contains(G107ExportFixtures.promptText(st), finalName), ...
                    [Site ': the consent prompt must name the final file ' finalName]);
            else
                testCase.verifyTrue(contains(out.dialogText, finalName), ...
                    [Site ': with no prompt, a visible refusal must name ' finalName]);
            end
        end

        function testDestinationCollisionConsentNamesExactTarget(testCase, Site)
            %At each site, the consenting answer. If the callback asked, the
            %question named the exact final file and the file is replaced
            %with real content; if it did not ask, the file keeps its bytes
            %and a visible refusal names it. Before the fix: no prompt, yet
            %the file was replaced.
            S = localSite(testCase, Site);
            typed = fullfile(testCase.TmpDir, ['collision' S.typedExt]);
            final = [typed '.csv'];
            [~, fn, fe] = fileparts(final);
            finalName = [fn fe];
            G107ExportFixtures.plant(final, sprintf('PLANTED SENTINEL G109\n'));
            bytes0 = G107ExportFixtures.fileBytes(final);
            G107ExportFixtures.resetDialogs({typed, final}, 'consent');
            out = G107ExportFixtures.fire(S.button);
            testCase.verifyEmpty(out.errorText, [Site ': the save must not throw']);
            st = G107ExportFixtures.dialogState();
            if G107ExportFixtures.prompted(st)
                testCase.verifyTrue(contains(G107ExportFixtures.promptText(st), finalName), ...
                    [Site ': the consent prompt must name ' finalName ', not only the typed name']);
                txt = fileread(final);
                testCase.verifyTrue(contains(txt, S.marker), [Site ': replaced with real content']);
                testCase.verifyFalse(contains(txt, 'PLANTED SENTINEL'), [Site ': the old bytes are gone']);
            else
                testCase.verifyEqual(G107ExportFixtures.fileBytes(final), bytes0, ...
                    [Site ': with no prompt naming it, ' finalName ' must not be replaced']);
                testCase.verifyTrue(contains(out.dialogText, finalName), ...
                    [Site ': with no prompt, a visible refusal must name ' finalName]);
            end
            testCase.verifyFalse(isfile(typed), [Site ': nothing under the typed name']);
        end

        function testDestinationExactCsvNameUsesDialogConsent(testCase, Site)
            %CONTROL at each site: a name that needs no normalization was
            %confirmed by the save dialog itself (the native panel asks before
            %returning an existing name), so it is replaced with no second
            %prompt.
            S = localSite(testCase, Site);
            final = fullfile(testCase.TmpDir, 'existing.csv');
            G107ExportFixtures.plant(final, sprintf('PLANTED SENTINEL G109\n'));
            G107ExportFixtures.resetDialogs({final}, 'decline');
            out = G107ExportFixtures.fire(S.button);
            testCase.verifyEmpty(out.errorText, [Site ': the save must not throw']);
            txt = fileread(final);
            testCase.verifyTrue(contains(txt, S.marker), [Site ': replaced with real content']);
            testCase.verifyFalse(contains(txt, 'PLANTED SENTINEL'), [Site ': the old bytes are gone']);
            testCase.verifyEqual(G107ExportFixtures.folderFiles(testCase.TmpDir), {'existing.csv'}, ...
                [Site ': only the saved file remains (no temporary file is left behind)']);
            st = G107ExportFixtures.dialogState();
            testCase.verifyFalse(G107ExportFixtures.prompted(st), ...
                [Site ': one valid consent suffices; no duplicate prompt']);
        end

        %% ================= G109, write-failure cleanup ====================

        function testWriteFailureAfterOpenKeepsTargetAndCleansUp(testCase, Site)
            %At each CSV callback, a write that fails AFTER a file was opened.
            %The existing target the save dialog confirmed must keep its
            %bytes, no partial or temporary file may remain in the folder,
            %every file the callback opened must be closed, and the user must
            %be told which file was not saved. The failure is injected
            %through the callback's own arguments (a table the writer cannot
            %finish: an event table without ss_errvar, a View All table
            %without its last column, a criterion "table" that is not a
            %table). The recording fopen, in watch-only mode, proves a file
            %really was opened for writing before the failure. Before the
            %fix: the callback opened the target itself with 'w', so the
            %target was already truncated to a partial file, the fprintf
            %writers left it open, and the only message was MATLAB's.
            S = localSite(testCase, Site);
            target = fullfile(testCase.TmpDir, 'kept.csv');
            G107ExportFixtures.plant(target, sprintf('PLANTED SENTINEL G109\n'));
            bytes0 = G107ExportFixtures.fileBytes(target);
            files0 = G107ExportFixtures.folderFiles(testCase.TmpDir);
            args = localFailingArgs(S);
            w = G107ExportFixtures.installFopenWatch(testCase.TmpDir, false);
            testCase.addTeardown(@() G107ExportFixtures.removeFopenWatch(w));
            testCase.assertTrue(startsWith(which('fopen'), w), ...
                'the recording fopen must be the resolved fopen');
            fids0 = G107ExportFixtures.openFileIDs();
            G107ExportFixtures.resetDialogs({target}, 'decline');
            out = G107ExportFixtures.fire(S.button, args{:});
            opened = G107ExportFixtures.fopenWrites();
            leaked = setdiff(G107ExportFixtures.openFileIDs(), fids0);
            for f = leaked(:)'
                fclose(f);   % a leak must not outlive the test
            end
            G107ExportFixtures.removeFopenWatch(w);
            testCase.assertNotEmpty(opened, [Site ': precondition, a file was ' ...
                'opened for writing before the injected failure']);
            testCase.assertTrue(~isempty(out.errorText) || ~isempty(strtrim(out.dialogText)), ...
                [Site ': precondition, the injected failure must be reported']);
            testCase.verifyEqual(G107ExportFixtures.fileBytes(target), bytes0, ...
                [Site ': the confirmed target must keep its bytes when the write fails part-way']);
            testCase.verifyEqual(G107ExportFixtures.folderFiles(testCase.TmpDir), files0, ...
                [Site ': no partial or temporary file may remain']);
            testCase.verifyEmpty(leaked, [Site ': every file the callback opened must be closed']);
            testCase.verifyTrue(contains([out.errorText newline out.dialogText], 'kept.csv'), ...
                [Site ': the user must be told that kept.csv was not saved']);
        end

        function testStageAllocationFailureKeepsTargetAndSaysSo(testCase, Site)
            %At each CSV callback, the folder refuses the file the callback
            %opens for writing (the recording fopen in refusal mode answers
            %-1, as a folder that refuses new files does). The existing
            %target keeps its bytes, nothing is left behind or left open, no
            %MATLAB error escapes the callback, and a dialog names the file
            %that was not saved. Before the fix: the refused open returned -1
            %and the next fprintf threw MATLAB's "Invalid file identifier"
            %with no word about the file. (The target survived before the fix
            %only because it was never opened, so on the old code this test
            %fails on the report, not on the bytes.)
            S = localSite(testCase, Site);
            target = fullfile(testCase.TmpDir, 'kept.csv');
            G107ExportFixtures.plant(target, sprintf('PLANTED SENTINEL G109\n'));
            bytes0 = G107ExportFixtures.fileBytes(target);
            files0 = G107ExportFixtures.folderFiles(testCase.TmpDir);
            w = G107ExportFixtures.installFopenWatch(testCase.TmpDir, true);
            testCase.addTeardown(@() G107ExportFixtures.removeFopenWatch(w));
            testCase.assertTrue(startsWith(which('fopen'), w), ...
                'the refusing fopen must be the resolved fopen');
            fids0 = G107ExportFixtures.openFileIDs();
            G107ExportFixtures.resetDialogs({target}, 'decline');
            out = G107ExportFixtures.fire(S.button);
            refused = G107ExportFixtures.fopenWrites();
            leaked = setdiff(G107ExportFixtures.openFileIDs(), fids0);
            for f = leaked(:)'
                fclose(f);
            end
            G107ExportFixtures.removeFopenWatch(w);
            testCase.assertNotEmpty(refused, [Site ': precondition, the callback ' ...
                'tried to open a file for writing and was refused']);
            testCase.verifyEqual(G107ExportFixtures.fileBytes(target), bytes0, ...
                [Site ': the existing target must keep its bytes']);
            testCase.verifyEqual(G107ExportFixtures.folderFiles(testCase.TmpDir), files0, ...
                [Site ': nothing may be left behind']);
            testCase.verifyEmpty(leaked, [Site ': nothing may be left open']);
            testCase.verifyEmpty(out.errorText, ...
                [Site ': a refused file is reported in a dialog, not thrown']);
            testCase.verifyTrue(contains(out.dialogText, 'kept.csv'), ...
                [Site ': the dialog must name kept.csv']);
        end

        %% ============ G109, close and write-status failures ===============
        %All fault tests below and their controls fire the real Save
        %callback with the save dialog confirming an existing kept.csv, with
        %the recording fopen (watch only) plus recording fclose, ferror and
        %movefile doubles (G107ExportFixtures.installIoFaults) and the
        %errordlg recorder on the path. Each I/O double forwards to its
        %built-in and fails only a budgeted number of calls on files under
        %the test folder; the recorder counts dialogs. The shared checks live
        %in localFaultRun, localVerifyFaultRefused (exactly one dialog) and
        %localVerifyFaultSaved (none); the dialog texts are checked by
        %G107ExportFixtures.verifyWriteFailureText and
        %verifyReplaceFailureText.

        function testCloseFailureKeepsTargetAndCleansUp(testCase, Site)
            %A close that fails after the data were written (fclose returns
            %-1 after the built-in has really closed the file; a flush-time
            %disk-full looks like this). The confirmed target must keep its
            %bytes, nothing may be left behind or left open, and the user must
            %be told which file was not saved. At the criterion site the first
            %close the double meets is the writer's own header close before
            %writetable appends. The one dialog reports an incomplete write,
            %not a failed replacement. Before the fix: the close status was
            %ignored, the file was renamed over the target (the baseline
            %wrote the target itself), and nothing was said.
            R = localFaultRun(testCase, Site, struct('fclose', 1, 'ferror', 0, 'movefile', 0));
            testCase.verifyNumElements(R.io.fcloseFailed, 1, [Site ': the injected close ' ...
                'failure was delivered to a file the callback wrote (seam fired)']);
            localVerifyFaultRefused(testCase, R, Site);
            G107ExportFixtures.verifyWriteFailureText(testCase, ...
                G107ExportFixtures.errordlgText(R.dialogCalls), Site);
        end

        function testCloseStatusControlReplacesTarget(testCase, Site)
            %CONTROL for the close double: with a budget of 0 it closes and
            %reports the built-in's own status, and the save replaces the
            %confirmed target with real content, leaving nothing behind.
            R = localFaultRun(testCase, Site, struct('fclose', 0, 'ferror', 0, 'movefile', 0));
            testCase.verifyNotEmpty(R.io.fcloseSeen, [Site ': the close double saw ' ...
                'the callback close the file it wrote (seam fired)']);
            localVerifyFaultSaved(testCase, R, Site);
        end

        function testWriteErrorStatusKeepsTargetAndCleansUp(testCase, Site)
            %A write error the stream recorded without fprintf throwing
            %(ferror returns a message for the file the callback wrote). The
            %callback must consult it before replacing the target, and then
            %refuse as above, with one dialog reporting an incomplete write.
            %Before the fix: ferror was never consulted, so the double never
            %fired and the target was replaced.
            R = localFaultRun(testCase, Site, struct('fclose', 0, 'ferror', 1, 'movefile', 0));
            testCase.verifyNumElements(R.io.ferrorFailed, 1, [Site ': the callback ' ...
                'consulted ferror on the file it wrote and received the injected error']);
            localVerifyFaultRefused(testCase, R, Site);
            G107ExportFixtures.verifyWriteFailureText(testCase, ...
                G107ExportFixtures.errordlgText(R.dialogCalls), Site);
        end

        function testWriteErrorEmptyControlReplacesTarget(testCase, Site)
            %CONTROL for the ferror double: with a budget of 0 it forwards the
            %built-in's (empty) message and the save succeeds.
            R = localFaultRun(testCase, Site, struct('fclose', 0, 'ferror', 0, 'movefile', 0));
            localVerifyFaultSaved(testCase, R, Site);
        end

        function testRenameRefusalKeepsTargetAndReportsHonestly(testCase, Site)
            %The finished file could not be moved onto the target (movefile
            %returns status 0 and moves nothing). Here the double guarantees
            %that nothing touched kept.csv, so its bytes are checked; a real
            %operating-system failure part-way through a replace is not
            %observable from MATLAB, so the ONE dialog must say the new file
            %could not replace kept.csv and ask the user to check it, and must
            %say none of what the write-failure dialog says ("could not be
            %written completely", "did not replace any existing file") and
            %not that the existing file is unchanged (review F2: the earlier
            %check, the word "replace", also matched the write-failure text).
            %The temporary file must still be removed. Before the fix: the
            %dialog said the export stopped before the file was complete (it
            %had not) and that the existing file was unchanged; the baseline
            %wrote the target itself and never renamed.
            R = localFaultRun(testCase, Site, struct('fclose', 0, 'ferror', 0, 'movefile', 1));
            testCase.verifyNumElements(R.io.movefileFailed, 1, [Site ': the callback ' ...
                'tried to move its finished file onto the target and was refused (seam fired)']);
            localVerifyFaultRefused(testCase, R, Site);
            G107ExportFixtures.verifyReplaceFailureText(testCase, ...
                G107ExportFixtures.errordlgText(R.dialogCalls), 'kept.csv', Site);
        end

        function testMoveExceptionReportsReplaceFailureOnce(testCase, Site)
            %Review F1. The final move RAISES instead of returning status 0
            %(the movefile double in throw mode moves nothing and throws even
            %though its outputs were requested; MATLAB documents status 0 for
            %a failed move, so this is the undocumented case). It must be
            %reported exactly like the refused move above: kept.csv keeps its
            %bytes, the temporary file is removed, nothing stays open, and
            %ONE dialog says the new file could not replace kept.csv, with
            %none of the write-failure text. EXPECTED TO FAIL on the phase-5
            %helper: the exception leaves commit before its own report and
            %before it marks the file settled, so the caller's abort finds
            %the temporary file and reports that the new file "could not be
            %written completely, so it did not replace any existing file",
            %which nobody observed.
            R = localFaultRun(testCase, Site, struct('fclose', 0, 'ferror', 0, ...
                'movefile', 1, 'movefileThrows', true));
            testCase.verifyNumElements(R.io.movefileFailed, 1, [Site ': the callback ' ...
                'tried to move its finished file onto the target and the move raised (seam fired)']);
            localVerifyFaultRefused(testCase, R, Site);
            G107ExportFixtures.verifyReplaceFailureText(testCase, ...
                G107ExportFixtures.errordlgText(R.dialogCalls), 'kept.csv', Site);
        end

        function testRenameControlReplacesTarget(testCase, Site)
            %CONTROL for the movefile double: with a budget of 0 it forwards
            %the move and the save replaces the target. (The seam check here
            %holds for a writer that renames a finished file onto the target;
            %the baseline, which wrote the target in place, never renamed.)
            R = localFaultRun(testCase, Site, struct('fclose', 0, 'ferror', 0, 'movefile', 0));
            testCase.verifyNotEmpty(R.io.movefileSeen, [Site ': the movefile double ' ...
                'saw the finished file moved onto the target (seam fired)']);
            localVerifyFaultSaved(testCase, R, Site);
        end

        %% ============ G109, temporary-file names ==========================

        function testStageNameCollisionsAndStageNamedTargetSave(testCase, Site)
            %Review F3. The temporary file is psyrat_partial_<k>.csv whatever
            %the target, so every export in a folder competes for the same
            %names and a crash can leave one behind.
            %(1) Beside a nonempty psyrat_partial_1.csv, an EMPTY
            %psyrat_partial_2.csv and a FOLDER psyrat_partial_3.csv
            %(G107ExportFixtures.plantStageDecoys), a save onto the
            %dialog-confirmed kept.csv must replace kept.csv with real
            %content and leave every decoy byte- and tree-identical, with
            %nothing added. The zero-budget movefile double shows that the
            %file moved onto kept.csv carried a psyrat_partial_ name that is
            %none of the decoys, so the decoys really stood in the way.
            %(2) A save to a NEW file named psyrat_partial_1.csv, in a folder
            %of its own, must write exactly that file: the file moved onto
            %it must not be the target itself.
            %Guards in psyrat_export_stage 'open' this depends on: the
            %existence skip (else the EMPTY decoy is adopted and moved away),
            %the folder skip (else the folder refuses the open and the save
            %is refused), the existence skip with the append-mode byte guard
            %(the nonempty decoy), and the target-equality skip (else the new
            %file is its own temporary file, the move onto itself fails and
            %the cleanup deletes it). Passes on the phase-5 helper; it was not
            %run against a helper without those guards.
            S = localSite(testCase, Site);
            target = fullfile(testCase.TmpDir, 'kept.csv');
            G107ExportFixtures.plant(target, sprintf('PLANTED SENTINEL G109\n'));
            decoys = G107ExportFixtures.plantStageDecoys(testCase.TmpDir);
            tree0 = G107ExportFixtures.folderTree(testCase.TmpDir);
            [out, io] = localObservedSave(testCase, S, target);
            testCase.verifyEmpty(out.errorText, [Site ': the save beside the decoys must not throw']);
            testCase.verifyEmpty(strtrim(out.dialogText), [Site ': no refusal beside the decoys']);
            moves = localStageMoves(io);
            testCase.assertNumElements(moves, 1, [Site ': one temporary file was ' ...
                'moved onto kept.csv (seam fired)']);
            [where, n, e] = fileparts(moves{1});
            testCase.verifyEqual(where, testCase.TmpDir, [Site ': the temporary file sat beside the target']);
            testCase.verifyTrue(~isempty(regexp([n e], '^psyrat_partial_\d+\.csv$', 'once')) && ...
                ~any(strcmpi([n e], decoys)), sprintf(['%s: the temporary file (%s) carried a ' ...
                'psyrat_partial_ name and was none of the decoys'], Site, [n e]));
            txt = fileread(target);
            testCase.verifyTrue(contains(txt, S.marker), [Site ': kept.csv replaced with real content']);
            testCase.verifyFalse(contains(txt, 'PLANTED SENTINEL'), [Site ': the old bytes are gone']);
            tree1 = G107ExportFixtures.folderTree(testCase.TmpDir);
            testCase.verifyEqual(tree1.paths, tree0.paths, [Site ': nothing added or removed ' ...
                'beside the decoys, at any depth']);
            if isequal(tree1.paths, tree0.paths)
                other = ~strcmp(tree0.paths, 'kept.csv');
                testCase.verifyEqual(tree1.bytes(other), tree0.bytes(other), ...
                    [Site ': every decoy keeps its bytes']);
            end

            sub = fullfile(testCase.TmpDir, 'named');
            mkdir(sub);
            named = fullfile(sub, 'psyrat_partial_1.csv');
            [out, io] = localObservedSave(testCase, S, named);
            testCase.verifyEmpty(out.errorText, [Site ': the save to psyrat_partial_1.csv must not throw']);
            testCase.verifyEmpty(strtrim(out.dialogText), [Site ': no refusal for psyrat_partial_1.csv']);
            moves = localStageMoves(io);
            testCase.assertNumElements(moves, 1, [Site ': one temporary file was ' ...
                'moved onto psyrat_partial_1.csv (seam fired)']);
            testCase.verifyFalse(strcmpi(moves{1}, named), [Site ': the target ' ...
                'psyrat_partial_1.csv was not its own temporary file']);
            testCase.assertTrue(isfile(named), [Site ': psyrat_partial_1.csv was written']);
            testCase.verifyTrue(contains(fileread(named), S.marker), ...
                [Site ': psyrat_partial_1.csv holds real content']);
            treeN = G107ExportFixtures.folderTree(sub);
            testCase.verifyEqual(treeN.paths, {'psyrat_partial_1.csv'}, ...
                [Site ': only the saved file remains in its folder']);
        end

        %% ============ G109, folder targets and long names ================

        function testExactCsvFolderTargetIsRefused(testCase, Site)
            %The save dialog returns a name ending in .csv that is an
            %EXISTING FOLDER (reachable through the callback's own dialog
            %answer; whether the native macOS panel can return a folder is
            %not verified here). Nothing may be written into that folder or
            %beside it, no file may be left open, nothing may be thrown, and
            %a dialog must name the folder. Before the fix: the name skipped
            %the folder check, and the finished file was moved INTO the
            %folder under its temporary name with no word (the baseline threw
            %"Invalid file identifier").
            S = localSite(testCase, Site);
            target = fullfile(testCase.TmpDir, 'existing.csv');
            mkdir(target);
            G107ExportFixtures.plant(fullfile(target, 'inside.txt'), sprintf('FOLDER SENTINEL G109\n'));
            localVerifyFolderRefusal(testCase, S, target, target, 'existing.csv', Site);
        end

        function testAppendedCsvFolderTargetIsRefused(testCase, Site)
            %CONTROL: a typed name that normalizes onto an existing folder
            %(x.xlsx -> x.xlsx.csv; x.txt -> x.txt.csv at the criterion
            %site) is refused the same way; the appended-name branch already
            %checks for a folder.
            S = localSite(testCase, Site);
            typed = fullfile(testCase.TmpDir, ['x' S.typedExt]);
            final = [typed '.csv'];
            mkdir(final);
            G107ExportFixtures.plant(fullfile(final, 'inside.txt'), sprintf('FOLDER SENTINEL G109\n'));
            [~, fn, fe] = fileparts(final);
            localVerifyFolderRefusal(testCase, S, typed, final, [fn fe], Site);
        end

        function testTargetAtNameLengthLimitStillSaves(testCase, Site)
            %A file name at the longest length the test folder actually
            %accepts must still save: the callback may not need a longer name
            %for its temporary file than the target itself. The limit is
            %measured in the folder by bisection (G107ExportFixtures
            %.longestCreatableName), not assumed, and the test fails, not
            %skips, when no limit exists between 16 and 4096 characters. On
            %the configured platform (macOS, APFS) that limit is the file
            %system's 255-byte name limit. Before the fix: the temporary name
            %appended ".psyrat_partial_1" to the target name, so it could not
            %be created and the save was refused (the baseline wrote the
            %target directly and saved).
            S = localSite(testCase, Site);
            n = G107ExportFixtures.longestCreatableName(testCase.TmpDir, '.csv', 16, 4096);
            testCase.assertFalse(isnan(n), sprintf(['precondition: the test folder must ' ...
                'refuse some file name between 16 and 4096 characters (%s); on a file ' ...
                'system without such a limit this test cannot reach the long-name case'], ...
                testCase.TmpDir));
            name = G107ExportFixtures.nameOfLength(n, '.csv');
            testCase.assertTrue(G107ExportFixtures.canCreateName(testCase.TmpDir, name) && ...
                ~G107ExportFixtures.canCreateName(testCase.TmpDir, ...
                G107ExportFixtures.nameOfLength(n + 1, '.csv')), ...
                'precondition: the target name is exactly at the folder''s measured limit');
            target = fullfile(testCase.TmpDir, name);
            fids0 = G107ExportFixtures.openFileIDs();
            G107ExportFixtures.resetDialogs({target}, 'decline');
            out = G107ExportFixtures.fire(S.button);
            leaked = setdiff(G107ExportFixtures.openFileIDs(), fids0);
            for f = leaked(:)'
                fclose(f);
            end
            testCase.verifyEmpty(out.errorText, [Site ': the save must not throw']);
            testCase.verifyEmpty(strtrim(out.dialogText), sprintf(['%s: a %d-character ' ...
                'name the folder accepts must save without a refusal'], Site, n));
            testCase.verifyTrue(isfile(target), [Site ': the long-named file is written']);
            if isfile(target)
                testCase.verifyTrue(contains(fileread(target), S.marker), [Site ': real content']);
            end
            tree = G107ExportFixtures.folderTree(testCase.TmpDir);
            testCase.verifyEqual(tree.paths, {name}, [Site ': only the saved file remains']);
            testCase.verifyEmpty(leaked, [Site ': nothing may be left open']);
        end

        %% ============ G109, label escaping on the sibling paths ==========

        function testGuiLineBreakGroupLabelsStayOneField(testCase)
            %Group names holding a line feed and a carriage return (a built
            %struct; the preflight refuses such labels on every new fit). The
            %event-section labels ('<group> - <event>'), the difference-section
            %labels and the empty-status line that names the group must each
            %stay one field, and every section's rows must keep their fields
            %and values. Before the fix: all three went through a writer that
            %quotes commas and quotes only, so each split into two records.
            %(The event-label half was added with the fix; that it failed
            %before follows from that writer, psyrat_quote_csv_header, leaving
            %a line break unquoted. It was not run against the old code.)
            pd = G107ExportFixtures.summarize(G107ExportFixtures.lineBreakGroupRel(), 1, 0.8, 0.5, 2);
            DA = pd.relsummary.group(1).diffscore.ssrel_table;
            bad = pd;
            bad.relsummary.group(1).diffscore.ssrel_table = DA([], :);
            recs = localGuiSave(testCase, pd, 'line_breaks.csv', bad);
            el = G107ExportFixtures.eventSectionLabels(pd);
            dl = G107ExportFixtures.diffSectionLabels(pd);
            testCase.assertTrue(all(contains(el(1:2), newline)) && ...
                all(contains(el(3:4), char(13))), ...
                'weld: the event-section labels carry the LF and CR group names');
            S = G107ExportFixtures.guiSections(recs, [el, dl]);
            testCase.verifyEqual({S.label}, [el, dl], ...
                'each event-section and difference-section label (LF and CR) is one field, in order');
            coef = G107ExportFixtures.eventCoefLabel(pd);
            k = 0;
            for g = 1:2
                for e = 1:2
                    k = k + 1;
                    idx = find(strcmp({S.label}, el{k}), 1);
                    if ~isempty(idx)
                        localVerifyGuiEventSection(testCase, S(idx), ...
                            pd.relsummary.group(g).event(e).ssrel_table, coef, ...
                            sprintf('event section %d', k));
                    end
                end
            end
            idx = find(strcmp({S.label}, dl{1}), 1);
            if ~isempty(idx)
                testCase.verifyEmpty(S(idx).data, 'the empty group writes no rows');
                testCase.verifyTrue(localHasEmptyStatus([S(idx).notes, S(idx).trailing], ...
                    pd.rel.groups{1}), 'the empty status naming the LF group is one field');
            end
            idx = find(strcmp({S.label}, dl{2}), 1);
            if ~isempty(idx)
                R = G107ExportFixtures.diffReference(pd, 2, 1, 0.8, 0.5, 2);
                localVerifyGuiDiffSection(testCase, S(idx), pd, 2, R, 'CR group');
            end
        end

        function testViewAllCommaLabelStaysOneField(testCase)
            %The View All difference table's Label cell carries runtime event
            %names; an event name with a comma (a built struct; the preflight
            %refuses such labels on every new fit) must stay one field, so no
            %row shifts. Before the fix: the label was written bare and
            %'loss, early' split into two fields.
            pd = G107ExportFixtures.icDiffViewData();
            pd.rel.events = {'loss, early'; 'gain'};
            figs = G107ExportFixtures.openWindows(testCase, ...
                @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
            win = G107ExportFixtures.openChild(testCase, figs, 'psyrat_viewdiffcomp');
            button = G107ExportFixtures.findButton(testCase, win, 'Save Table');
            target = fullfile(testCase.TmpDir, 'view_all.csv');
            G107ExportFixtures.resetDialogs({target}, 'decline');
            out = G107ExportFixtures.fire(button);
            testCase.assertEmpty(out.errorText, ['View All Save threw: ' out.errorText]);
            testCase.assertTrue(isfile(target), 'View All Save wrote no file');
            recs = G107ExportFixtures.parseCsv(fileread(target));
            h = find(cellfun(@(r) numel(r) > 1 && strcmp(r{1}, 'Label'), recs), 1);
            testCase.assertNotEmpty(h, 'the View All column row was not found');
            header = recs{h};
            rows = recs(h + 1:end);
            expected = {'loss, early', 'gain', 'diff score'};
            testCase.verifyNumElements(rows, numel(expected), 'one record per table row');
            for k = 1:min(numel(rows), numel(expected))
                testCase.verifyNumElements(rows{k}, numel(header), sprintf( ...
                    'row %d keeps the column row''s %d fields', k, numel(header)));
                testCase.verifyEqual(rows{k}{1}, expected{k}, sprintf('row %d label', k));
            end
        end

    end
end

% -- fixtures and drivers -------------------------------------------------------

function pd = localMixedSummary(testCase, pd0, c)
%Summarize at a cutoff staged between the two smallest selected-bound values,
%so the difference flags are mixed and a flag mix-up is visible.
cut = G107ExportFixtures.mixedCutoff(testCase, pd0, c.gc, c.ci, c.meas);
pd = G107ExportFixtures.summarize(pd0, c.gc, c.ci, cut, c.meas);
flags = [];
for g = 1:numel(pd.relsummary.group)
    flags = [flags; pd.relsummary.group(g).diffscore.ssrel_table.ind2include(:)]; %#ok<AGROW>
end
testCase.assertTrue(any(flags) && any(~flags), ...
    'weld: the staged cutoff must include some participants and exclude others');
end

function [recs, target] = localGuiSave(testCase, pd, name, override)
%Build the real variance window from PD, fire Save Subject Coefficients with
%the dialog answering TMPDIR/NAME, and parse the file. OVERRIDE replaces the
%callback's summary (a malformed or typed-empty table injected at the seam).
if nargin < 4
    override = pd;
end
figs = G107ExportFixtures.openWindows(testCase, ...
    @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
button = G107ExportFixtures.findButton(testCase, figs, 'Save Subject');
target = fullfile(testCase.TmpDir, name);
G107ExportFixtures.resetDialogs({target}, 'decline');
out = G107ExportFixtures.fire(button, override);
testCase.assertEmpty(out.errorText, ['Save Subject Coefficients threw: ' out.errorText]);
testCase.assertTrue(isfile(target), 'Save Subject Coefficients wrote no file');
recs = G107ExportFixtures.parseCsv(fileread(target));
delete(figs(isvalid(figs)));
end

function S = localGuiSections(recs, pd)
S = G107ExportFixtures.guiSections(recs, ...
    [G107ExportFixtures.eventSectionLabels(pd), G107ExportFixtures.diffSectionLabels(pd)]);
end

function s = localSection(testCase, S, label)
idx = [];
if ~isempty(S)
    idx = find(strcmp({S.label}, label));
end
testCase.assertNumElements(idx, 1, sprintf('section "%s" must appear exactly once', label));
s = S(idx);
end

function V = localMalformedVariants(pd)
%Missing and malformed difference tables for group 1 (the injection points are
%stated: the summary seam, never the draws).
D = pd.relsummary.group(1).diffscore.ssrel_table;
V = cell(0, 2);
p = pd;
p.relsummary.group = rmfield(p.relsummary.group, 'diffscore');
V(end+1, :) = {'no diffscore field', p};
p = pd;
p.relsummary.group(1).diffscore = [];
V(end+1, :) = {'empty diffscore', p};
p = pd;
p.relsummary.group(1).diffscore = rmfield(p.relsummary.group(1).diffscore, 'ssrel_table');
V(end+1, :) = {'no ssrel_table', p};
p = pd;
p.relsummary.group(1).diffscore.ssrel_table = [];
V(end+1, :) = {'ssrel_table not a table', p};
p = pd;
p.relsummary.group(1).diffscore.ssrel_table = removevars(D, 'trls2');
V(end+1, :) = {'trls2 column missing', p};
p = pd;
p.relsummary.group(1).diffscore.ssrel_table = removevars(D, 'ind2include');
V(end+1, :) = {'ind2include column missing', p};
p = pd;
T = D;
T.dep_pt = num2cell(T.dep_pt);
p.relsummary.group(1).diffscore.ssrel_table = T;
V(end+1, :) = {'dep_pt not numeric', p};
end

function S = localSite(testCase, name)
%The window and Save button for one G109 site, the file content marker, and
%the typed extension the site normalizes by appending .csv.
S.name = name;
switch name
    case 'subjectCoefficients'
        pd = G107ExportFixtures.summarize(G107ExportFixtures.baseRel(), 1, 0.8, 0.5, 2);
        figs = G107ExportFixtures.openWindows(testCase, ...
            @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
        S.button = G107ExportFixtures.findButton(testCase, figs, 'Save Subject');
        S.marker = 'Subject-Specific Estimates';
        S.typedExt = '.xlsx';
    case 'viewAllDifference'
        pd = G107ExportFixtures.icDiffViewData();
        figs = G107ExportFixtures.openWindows(testCase, ...
            @() psyrat_variancet('psyrat_data', pd, 'gui', 1));
        win = G107ExportFixtures.openChild(testCase, figs, 'psyrat_viewdiffcomp');
        S.button = G107ExportFixtures.findButton(testCase, win, 'Save Table');
        S.marker = 'Label,Between-Person Std Dev';
        S.typedExt = '.xlsx';
    case 'criterionTable'
        pd = G107ExportFixtures.criterionSingData();
        figs = G107ExportFixtures.openWindows(testCase, @() psyrat_criterionfigures( ...
            'psyrat_prefs', G107ExportFixtures.criterionPrefs(0.5), ...
            'psyrat_data', pd, 'analysis', 'sing'));
        S.button = G107ExportFixtures.findButton(testCase, figs, 'Save');
        S.marker = 'Criterion Cutoff:';
        S.typedExt = '.txt';
    otherwise
        error('TestSubjectDifferenceExport:site', 'Unknown site %s.', name);
end
end

function args = localFailingArgs(S)
%Replacement callback arguments that make each site's writer fail part-way,
%after its output file is open (the site's own stored arguments, with one table
%the writer cannot finish).
cb = get(S.button, 'Callback');
switch S.name
    case 'subjectCoefficients'
        %the second event section reads ss_errvar after the first section and
        %the header are written
        pd = cb{2};
        T = pd.relsummary.group(1).event(2).ssrel_table;
        pd.relsummary.group(1).event(2).ssrel_table = removevars(T, 'ss_errvar');
        args = {pd};
    case 'viewAllDifference'
        %the row loop reads column 8 after the header is written
        vt = cb{3};
        args = {cb{2}, vt(:, 1:7)};
    case 'criterionTable'
        %writetable appends after the header lines are written
        args = {cb{2}, struct('not', 'a table'), cb{4}, cb{5}};
    otherwise
        error('TestSubjectDifferenceExport:site', 'Unknown site %s.', S.name);
end
end

function R = localFaultRun(testCase, site, fails)
%Fire one site's real Save callback onto an existing, dialog-confirmed
%kept.csv with the recording fopen (watch only) and the fclose, ferror and
%movefile doubles at the failure budgets in FAILS, and the errordlg recorder.
%Returns what the user saw (R.out, and R.dialogCalls: every errordlg call),
%the target's bytes and a recursive snapshot of the test folder before and
%after, the file IDs left open, and what each double saw and failed. The
%doubles are removed before any file is read back.
S = localSite(testCase, site);
R = struct();
R.marker = S.marker;
R.target = fullfile(testCase.TmpDir, 'kept.csv');
G107ExportFixtures.plant(R.target, sprintf('PLANTED SENTINEL G109\n'));
R.bytes0 = G107ExportFixtures.fileBytes(R.target);
R.tree0 = G107ExportFixtures.folderTree(testCase.TmpDir);
for f = {'fclose', 'ferror', 'movefile'}
    testCase.assertEqual(exist(f{1}, 'builtin'), 5, sprintf( ...
        'precondition: %s is a built-in the double can forward to', f{1}));
end
wDlg = G107ExportFixtures.installErrordlgRecorder();
testCase.addTeardown(@() G107ExportFixtures.removeErrordlgRecorder(wDlg));
testCase.assertTrue(startsWith(which('errordlg'), wDlg), ...
    'the recording errordlg must be the resolved errordlg');
wOpen = G107ExportFixtures.installFopenWatch(testCase.TmpDir, false);
testCase.addTeardown(@() G107ExportFixtures.removeFopenWatch(wOpen));
wIo = G107ExportFixtures.installIoFaults(testCase.TmpDir, fails);
testCase.addTeardown(@() G107ExportFixtures.removeIoFaults(wIo));
testCase.assertTrue(startsWith(which('fopen'), wOpen), ...
    'the recording fopen must be the resolved fopen');
for f = {'fclose', 'ferror', 'movefile'}
    testCase.assertTrue(startsWith(which(f{1}), wIo), sprintf( ...
        'the %s double must be the resolved %s', f{1}, f{1}));
end
fids0 = G107ExportFixtures.openFileIDs();
G107ExportFixtures.resetDialogs({R.target}, 'decline');
R.out = G107ExportFixtures.fire(S.button);
R.dialogCalls = G107ExportFixtures.errordlgCalls();
R.opened = G107ExportFixtures.fopenWrites();
R.io = G107ExportFixtures.ioFaultLog();
R.leaked = setdiff(G107ExportFixtures.openFileIDs(), fids0);
G107ExportFixtures.removeErrordlgRecorder(wDlg);
G107ExportFixtures.removeIoFaults(wIo);
G107ExportFixtures.removeFopenWatch(wOpen);
for f = R.leaked(:)'
    fclose(f);   % after the doubles are gone; a leak must not outlive the test
end
R.bytes = G107ExportFixtures.fileBytes(R.target);
R.tree = G107ExportFixtures.folderTree(testCase.TmpDir);
testCase.assertNotEmpty(R.opened, [site ': precondition, the callback opened a ' ...
    'file for writing (a write was attempted)']);
end

function [out, io] = localObservedSave(testCase, S, answer)
%Fire S's Save callback with the save dialog answering ANSWER while
%zero-budget fclose, ferror and movefile doubles record what the export did
%(nothing is failed). Returns what the user saw and the doubles' log; the
%doubles are removed before returning, so files can be read back.
w = G107ExportFixtures.installIoFaults(testCase.TmpDir, ...
    struct('fclose', 0, 'ferror', 0, 'movefile', 0));
testCase.addTeardown(@() G107ExportFixtures.removeIoFaults(w));
testCase.assertTrue(startsWith(which('movefile'), w), ...
    'the recording movefile must be the resolved movefile');
G107ExportFixtures.resetDialogs({answer}, 'decline');
out = G107ExportFixtures.fire(S.button);
io = G107ExportFixtures.ioFaultLog();
G107ExportFixtures.removeIoFaults(w);
end

function moves = localStageMoves(io)
%The watched moves whose source is a psyrat_partial_<k>.csv file (the
%temporary file moved onto the target). Filtering by name keeps the check
%independent of any move writetable might make internally on the criterion
%path.
moves = {};
for k = 1:numel(io.movefileSeen)
    [~, n, e] = fileparts(io.movefileSeen{k});
    if ~isempty(regexpi([n e], '^psyrat_partial_\d+\.csv$', 'once'))
        moves{end+1} = io.movefileSeen{k}; %#ok<AGROW>
    end
end
end

function localVerifyFaultRefused(testCase, R, site)
%A refused save: the confirmed target byte-identical, the folder unchanged
%recursively (no temporary or partial file anywhere), nothing left open, and
%exactly one dialog (one errordlg call; the recorder counts calls because a
%'modal' errordlg reuses an open window with the same title), which names the
%target.
testCase.verifyEqual(R.bytes, R.bytes0, [site ': the confirmed kept.csv keeps its bytes']);
testCase.verifyEqual(R.tree, R.tree0, [site ': nothing is left behind (no ' ...
    'temporary or partial file, at any depth)']);
testCase.verifyEmpty(R.leaked, [site ': every file the callback opened is closed']);
testCase.verifyTrue(contains(R.out.dialogText, 'kept.csv'), [site ': a dialog ' ...
    'tells the user that kept.csv was not saved']);
testCase.verifyNumElements(R.dialogCalls, 1, [site ': exactly one dialog reports ' ...
    'the failure (no second, contradictory report)']);
end

function localVerifyFaultSaved(testCase, R, site)
%A successful save under the doubles: no error, no refusal, the target replaced
%with real content, only the target left in the folder, nothing left open.
testCase.verifyEmpty(R.out.errorText, [site ': the save must not throw']);
testCase.verifyEmpty(strtrim(R.out.dialogText), [site ': no refusal']);
testCase.verifyEmpty(R.dialogCalls, [site ': no errordlg call']);
txt = char(R.bytes(:)');
testCase.verifyTrue(contains(txt, R.marker), [site ': replaced with real content']);
testCase.verifyFalse(contains(txt, 'PLANTED SENTINEL'), [site ': the old bytes are gone']);
testCase.verifyEqual(R.tree.paths, {'kept.csv'}, [site ': only the saved file remains']);
testCase.verifyEmpty(R.leaked, [site ': nothing is left open']);
end

function localVerifyFolderRefusal(testCase, S, answer, folder, name, site)
%Answer the save dialog with ANSWER, which names (or normalizes onto) the
%existing FOLDER. Nothing may change in the test folder at any depth, nothing
%may be left open or thrown, and a dialog must name the folder.
testCase.assertTrue(isfolder(folder), 'precondition: the target is a folder');
tree0 = G107ExportFixtures.folderTree(testCase.TmpDir);
fids0 = G107ExportFixtures.openFileIDs();
G107ExportFixtures.resetDialogs({answer}, 'decline');
out = G107ExportFixtures.fire(S.button);
leaked = setdiff(G107ExportFixtures.openFileIDs(), fids0);
for f = leaked(:)'
    fclose(f);
end
testCase.verifyEqual(G107ExportFixtures.folderTree(testCase.TmpDir), tree0, ...
    [site ': nothing is written into ' name ' or beside it']);
testCase.verifyEmpty(leaked, [site ': nothing may be left open']);
testCase.verifyEmpty(out.errorText, [site ': a folder target is refused in a dialog, not thrown']);
testCase.verifyTrue(contains(out.dialogText, name), [site ': the refusal names ' name]);
end

function R = localRoutes()
%{name, stamp, a provenance line only that route emits}
R = { ...
    'cmdstanGaussian', @(pd) pd, 'Estimation engine: cmdstan'; ...
    'hmcGaussian', @localHmc, 'Estimation engine: hmc'; ...
    'gammaLogNu', @localGammaLogNu, 'Gamma subject-level difference residual (log-nu, POOLED)'; ...
    'gammaCopulaConcurrent', @localGammaCopula, 'Gamma copula dispersion parameterization: PER-PERSON log-nu'; ...
    'gammaLocationScale', @localGammaLs, 'Gamma subject-level difference reliability (location-scale, POOLED)'};
end

function pd = localHmc(pd)
pd.rel.engine = 'hmc';
end

function pd = localGammaLogNu(pd)
%Non-concurrent gamma log-nu analysis 8 (the is_diff_sserr stamp).
pd.rel.family = 'gamma';
pd.rel.gammascale = 1;
pd.rel.dispersion = 1;
pd.rel.out.gamma_disp = struct('s_v', [0.30; 0.31; 0.32; 0.33; 0.34], ...
    'rho_pv', -0.10 * ones(5, 1), 'disp_factor', 1.05 * ones(5, 1), ...
    'is_diff_sserr', true);
end

function pd = localGammaCopula(pd)
%Concurrent gamma analysis 8 (Gaussian copula, per-person log-nu,
%experimental): is_diff_copula and is_diff_sserr, dispersion = 1.
pd = localGammaLogNu(pd);
pd.rel.diffwpcov = 2;
pd.rel.out.gamma_disp.is_diff_copula = true;
end

function pd = localGammaLs(pd)
%Gamma location-scale analysis 8 (is_ls, is_diff, is_diff_sserr; two
%per-event columns).
pd.rel.family = 'gamma';
pd.rel.gammascale = 2;
ssd = repmat([0.20, 0.25], 5, 1);
pd.rel.out.gamma_disp = struct('is_ls', true, 'is_diff', true, ...
    'is_diff_sserr', true, 's_sd', ssd, 'rho_ps', repmat([0.10, 0.05], 5, 1), ...
    'disp_factor', exp(2 .* ssd .^ 2));
end

function pd = localAnalysis6()
%ic_sserrvar, two groups, no events (analysis 6).
pd = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
pd = G107ExportFixtures.addHeaderFields(pd, 'g107_analysis6.psyrat');
[pd, ~] = psyrat_relsummary('psyrat_data', pd, 'analysis', 'sing_sserr', ...
    'depcutoff', 0.5, 'meascutoff', 2, 'depcentmeas', 1, 'gcoeff', 1, 'CI', 0.95);
end

function pd = localAnalysis25()
%trt_sserrvar (analysis 25), the TestHeadlessReport settings.
pd = PsyRATTestDataFactory.makeTRTSSErrRelData();
pd = G107ExportFixtures.addHeaderFields(pd, 'g107_analysis25.psyrat');
[pd, ~] = psyrat_relsummary('psyrat_data', pd, 'analysis', 'trt_sserr', ...
    'gcoeff', 1, 'reltype', 3, 'depcutoff', 0.5, 'meascutoff', 2, 'CI', 0.95);
end

function note = localInjectedNote()
%Injected at the summary seam (no real admissibility record is generated);
%the comma makes the note need quoting.
note = ['G107 injected admissibility note, not a real record: 1 of 5 draws ' ...
    'of the relative error variance was negative.'];
end

function localCloseDialogs(out)
%Close the dialogs one fired callback opened, so the next variant starts clean.
if ~isempty(out.dialogs)
    delete(out.dialogs(isvalid(out.dialogs)));
end
end

function msg = localError(fn)
msg = '';
try
    fn();
catch ME
    msg = ME.message;
end
end

function [first, rest] = localSplitAfterDate(txt)
nl = find(txt == newline, 2);
if numel(nl) < 2
    first = txt;
    rest = '';
    return;
end
first = txt(1:nl(1) - 1);
rest = txt(nl(2) + 1:end);
end

function rest = localExpectedEventOnlyGuiRest(pd)
%The shipped Save Subject Coefficients file for a design with no difference
%score, rebuilt from the callback's own format (psyrat_variancet.m,
%psyrat_savesscoeffs) and the shared header helpers, from the line after the
%timestamp to the end. The callback now quotes section labels with
%psyrat_csv_field (G109); this builder keeps the older header quoting for them,
%which writes the same bytes for a label without a line break, so the control
%also pins that the switch left ordinary labels byte-identical.
q = @psyrat_quote_csv_header;
L = {sprintf('PsyRAT v%s', pd.ver); ' '; ...
    q(sprintf('Dataset: %s', pd.rel.filename)); q(psyrat_iterline(pd.rel))};
prov = psyrat_provenance_lines(pd);
for k = 1:numel(prov)
    L{end+1, 1} = q(prov{k}); %#ok<AGROW>
end
L = [L; {sprintf('Reliability Cutoff: %0.4f', pd.relsummary.depcutoff); ...
    'Subject-Specific Estimates'; ' '; ' '}];
labels = G107ExportFixtures.eventSectionLabels(pd);
header = strjoin(G107ExportFixtures.guiEventHeader(G107ExportFixtures.eventCoefLabel(pd)), ',');
k = 0;
for g = 1:numel(pd.relsummary.group)
    for e = 1:numel(pd.relsummary.group(g).event)
        k = k + 1;
        L{end+1, 1} = q(labels{k}); %#ok<AGROW>
        L{end+1, 1} = header; %#ok<AGROW>
        T = pd.relsummary.group(g).event(e).ssrel_table;
        for i = 1:height(T)
            L{end+1, 1} = sprintf(['%s,%d,%d,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,' ...
                '%0.3f,%0.3f,%0.3f,%0.3f'], char(T.id{i}), T.ind2include(i), T.trls(i), ...
                T.dep_pt(i), T.dep_ll(i), T.dep_ul(i), T.icc_pt(i), T.icc_ll(i), ...
                T.icc_ul(i), T.sem_pt(i), T.sem_ll(i), T.sem_ul(i), T.ss_errvar(i)); %#ok<AGROW>
        end
    end
end
rest = [strjoin(L', newline) newline];
end

% -- verification --------------------------------------------------------------

function localVerifyStoredAgainstReference(testCase, D, R, what)
%The stored summary table against the raw-draw reference (1e-12): the roster,
%both counts, every summary and the flag.
testCase.verifyEqual(sort(cellstr(string(D.id))), sort(R.id), ...
    [what ': the roster is exactly the participants with both events']);
pairs = {'dep_pt', 'pt'; 'dep_ll', 'll'; 'dep_ul', 'ul'; 'icc_pt', 'icc_pt'; ...
    'icc_ll', 'icc_ll'; 'icc_ul', 'icc_ul'; 'sem_pt', 'sem_pt'; 'sem_ll', 'sem_ll'; ...
    'sem_ul', 'sem_ul'; 'ss_errvar', 'errvar'; 'trls1', 'trls1'; 'trls2', 'trls2'};
for j = 1:height(D)
    id = char(string(D.id(j)));
    k = find(strcmp(R.id, id));
    testCase.verifyNumElements(k, 1, sprintf('%s: reference row for %s', what, id));
    if numel(k) ~= 1
        continue;
    end
    for p = 1:size(pairs, 1)
        localVerifyNum(testCase, D.(pairs{p, 1})(j), R.(pairs{p, 2})(k), 1e-12, ...
            sprintf('%s, %s: %s against raw-draw arithmetic', what, id, pairs{p, 1}));
    end
    testCase.verifyEqual(logical(D.ind2include(j)), R.include(k), ...
        sprintf('%s, %s: flag against the reference bound and cutoff', what, id));
end
end

function localVerifyGuiEventSection(testCase, s, T, coef, what)
%An event section exactly as shipped: the column row, no stray records, one
%row per stored row in order, the event's own flag, its count, and every value
%at the printed precision (at least 3 decimals).
testCase.verifyEqual(s.header, G107ExportFixtures.guiEventHeader(coef), [what ': column row unchanged']);
testCase.verifyEmpty(s.trailing, [what ': no stray single-field record (a split ID leaves one)']);
testCase.verifyNumElements(s.data, height(T), [what ': one row per stored row']);
vars = {'dep_pt', 'dep_ll', 'dep_ul', 'icc_pt', 'icc_ll', 'icc_ul', ...
    'sem_pt', 'sem_ll', 'sem_ul', 'ss_errvar'};
for j = 1:min(numel(s.data), height(T))
    r = s.data{j};
    id = char(string(T.id(j)));
    tag = sprintf('%s row %d (%s)', what, j, id);
    testCase.verifyNumElements(r, 13, [tag ': 13 fields']);
    if numel(r) ~= 13
        continue;
    end
    testCase.verifyEqual(r{1}, id, [tag ': ID and order']);
    testCase.verifyEqual(G107ExportFixtures.parseFlag(r{2}), double(T.ind2include(j)), ...
        [tag ': the event''s own flag']);
    testCase.verifyEqual(G107ExportFixtures.parseNum(r{3}), T.trls(j), [tag ': count']);
    for v = 1:numel(vars)
        localVerifyPrinted(testCase, r{3 + v}, T.(vars{v})(j), 3, [tag ': ' vars{v}]);
    end
end
end

function localVerifyGuiDiffSection(testCase, s, pd, g, R, what)
%A difference section against the schema's column row, the stored table (ID,
%order, both counts, flag, every value at the printed precision) and the
%raw-draw reference (posterior mean, limits, counts, flag).
[keys, names] = G107ExportFixtures.guiDiffColumns();
testCase.verifyEqual(s.header, names, [what ': difference column row']);
if ~isequal(s.header, names)
    return;
end
col = @(key) find(strcmp(keys, key), 1);
rs = pd.relsummary;
D = rs.group(g).diffscore.ssrel_table;
gl = G107ExportFixtures.groupLabels(pd);
ev = G107ExportFixtures.eventNames(pd);
testCase.verifyEmpty(s.trailing, [what ': no stray single-field record after the column row']);
testCase.verifyNumElements(s.data, height(D), [what ': one row per stored row, none dropped or fabricated']);
pairs = {'pt', 'dep_pt'; 'll', 'dep_ll'; 'ul', 'dep_ul'; 'icc_pt', 'icc_pt'; ...
    'icc_ll', 'icc_ll'; 'icc_ul', 'icc_ul'; 'sem_pt', 'sem_pt'; 'sem_ll', 'sem_ll'; ...
    'sem_ul', 'sem_ul'; 'errvar', 'ss_errvar'};
for j = 1:min(numel(s.data), height(D))
    r = s.data{j};
    id = char(string(D.id(j)));
    tag = sprintf('%s row %d (%s)', what, j, id);
    testCase.verifyNumElements(r, numel(names), [tag ': field count']);
    if numel(r) ~= numel(names)
        continue;
    end
    testCase.verifyEqual(r{col('id')}, id, [tag ': ID and order follow the stored table']);
    testCase.verifyEqual(r{col('group')}, gl{g}, [tag ': group']);
    testCase.verifyEqual(lower(strtrim(r{col('scoretype')})), 'difference', [tag ': score type']);
    testCase.verifyEqual(r{col('eventA')}, ev{1}, [tag ': fitted event A']);
    testCase.verifyEqual(r{col('eventB')}, ev{2}, [tag ': fitted event B']);
    testCase.verifyEqual(lower(strtrim(r{col('direction')})), ...
        lower(G107ExportFixtures.directionText(ev{1}, ev{2})), [tag ': direction']);
    testCase.verifyEqual(lower(strtrim(r{col('coefficient')})), ...
        lower(G107ExportFixtures.coefLabel(rs.diffgcoeff)), [tag ': coefficient']);
    testCase.verifyEqual(G107ExportFixtures.parseFlag(r{col('include')}), ...
        double(D.ind2include(j)), [tag ': the DIFFERENCE flag']);
    localVerifyPrinted(testCase, r{col('cutoff')}, rs.depcutoff, 4, [tag ': cutoff']);
    testCase.verifyEqual(strtrim(r{col('bound')}), strtrim(rs.meascutoff), [tag ': bound']);
    testCase.verifyEqual(G107ExportFixtures.parseNum(r{col('trls1')}), D.trls1(j), [tag ': trls1']);
    testCase.verifyEqual(G107ExportFixtures.parseNum(r{col('trls2')}), D.trls2(j), [tag ': trls2']);
    for p = 1:size(pairs, 1)
        localVerifyPrinted(testCase, r{col(pairs{p, 1})}, D.(pairs{p, 2})(j), 3, ...
            [tag ': ' pairs{p, 1} ' against the stored table']);
    end
    kind = lower(r{col('intervalkind')});
    testCase.verifyTrue(contains(kind, 'equal-tailed') && contains(kind, 'credible'), ...
        [tag ': the interval is named an equal-tailed credible interval']);
    testCase.verifyTrue(localWidthOk(r{col('intervalwidth')}, rs.ciperc), ...
        [tag ': the interval width is the one the summary used']);
    k = find(strcmp(R.id, id));
    testCase.verifyNumElements(k, 1, [tag ': reference row']);
    if numel(k) == 1
        localVerifyPrinted(testCase, r{col('pt')}, R.pt(k), 3, [tag ': posterior mean against raw draws']);
        localVerifyPrinted(testCase, r{col('ll')}, R.ll(k), 3, [tag ': lower limit against raw draws']);
        localVerifyPrinted(testCase, r{col('ul')}, R.ul(k), 3, [tag ': upper limit against raw draws']);
        testCase.verifyEqual(G107ExportFixtures.parseNum(r{col('trls1')}), R.trls1(k), [tag ': trls1 against rel.data']);
        testCase.verifyEqual(G107ExportFixtures.parseNum(r{col('trls2')}), R.trls2(k), [tag ': trls2 against rel.data']);
        testCase.verifyEqual(G107ExportFixtures.parseFlag(r{col('include')}), double(R.include(k)), ...
            [tag ': flag against the reference bound and cutoff']);
    end
end
end

function localVerifyHeadlessDiffTable(testCase, T, pd, refs, what)
%The returned table: the schema, then every row exactly equal to the
%stored table (tolerance 0) and within 1e-12 of the raw-draw reference.
names = G107ExportFixtures.headlessDiffColumns();
testCase.verifyEqual(T.Properties.VariableNames, names, [what ': variables']);
if ~isequal(T.Properties.VariableNames, names)
    return;
end
valueAt = @(k, n) T.(n)(k);
localVerifyHeadlessRows(testCase, height(T), valueAt, pd, refs, 0, what);
end

function localVerifyHeadlessRows(testCase, n, valueAt, pd, refs, tol, what)
%Row-by-row check shared by the returned table and the reopened files.
%VALUEAT(k, name) returns row k's raw value for a headless variable name. Rows
%run group by group, each group in stored order.
rs = pd.relsummary;
gl = G107ExportFixtures.groupLabels(pd);
ev = G107ExportFixtures.eventNames(pd);
E = struct('g', {}, 'j', {});
for g = 1:numel(rs.group)
    for j = 1:height(rs.group(g).diffscore.ssrel_table)
        E(end+1) = struct('g', g, 'j', j); %#ok<AGROW>
    end
end
testCase.verifyEqual(n, numel(E), [what ': one row per stored difference row across groups']);
pairs = {'dep_pt', 'dep_pt'; 'dep_ll', 'dep_ll'; 'dep_ul', 'dep_ul'; 'icc_pt', 'icc_pt'; ...
    'icc_ll', 'icc_ll'; 'icc_ul', 'icc_ul'; 'sem_pt', 'sem_pt'; 'sem_ll', 'sem_ll'; ...
    'sem_ul', 'sem_ul'; 'diff_rel_errvar', 'ss_errvar'; 'trls1', 'trls1'; 'trls2', 'trls2'};
for k = 1:min(n, numel(E))
    g = E(k).g;
    D = rs.group(g).diffscore.ssrel_table;
    j = E(k).j;
    id = char(string(D.id(j)));
    tag = sprintf('%s row %d (%s)', what, k, id);
    testCase.verifyEqual(localText(valueAt(k,'Group')), gl{g}, [tag ': Group']);
    testCase.verifyEqual(lower(localText(valueAt(k,'ScoreType'))), 'difference', [tag ': ScoreType']);
    testCase.verifyEqual(localText(valueAt(k,'EventA')), ev{1}, [tag ': EventA']);
    testCase.verifyEqual(localText(valueAt(k,'EventB')), ev{2}, [tag ': EventB']);
    testCase.verifyEqual(lower(localText(valueAt(k,'Direction'))), ...
        lower(G107ExportFixtures.directionText(ev{1}, ev{2})), [tag ': Direction']);
    testCase.verifyEqual(lower(localText(valueAt(k,'Coefficient'))), ...
        lower(G107ExportFixtures.coefLabel(rs.diffgcoeff)), [tag ': Coefficient']);
    testCase.verifyEqual(localText(valueAt(k,'id')), id, [tag ': id and order']);
    testCase.verifyEqual(G107ExportFixtures.parseFlag(valueAt(k,'ind2include')), ...
        double(D.ind2include(j)), [tag ': the DIFFERENCE flag']);
    localVerifyNum(testCase, valueAt(k,'Cutoff'), rs.depcutoff, tol, [tag ': Cutoff']);
    testCase.verifyEqual(strtrim(localText(valueAt(k,'CutoffBound'))), strtrim(rs.meascutoff), ...
        [tag ': CutoffBound']);
    for p = 1:size(pairs, 1)
        localVerifyNum(testCase, valueAt(k,pairs{p, 1}), D.(pairs{p, 2})(j), tol, ...
            [tag ': ' pairs{p, 1}]);
    end
    kind = lower(localText(valueAt(k,'IntervalKind')));
    testCase.verifyTrue(contains(kind, 'equal-tailed') && contains(kind, 'credible'), ...
        [tag ': IntervalKind names an equal-tailed credible interval']);
    localVerifyNum(testCase, valueAt(k,'IntervalWidth'), rs.ciperc, tol, ...
        [tag ': IntervalWidth is the proportion the summary used']);
    R = refs{g};
    kk = find(strcmp(R.id, id));
    testCase.verifyNumElements(kk, 1, [tag ': reference row']);
    if numel(kk) == 1
        localVerifyNum(testCase, valueAt(k,'dep_pt'), R.pt(kk), 1e-12, [tag ': posterior mean against raw draws']);
        localVerifyNum(testCase, valueAt(k,'dep_ll'), R.ll(kk), 1e-12, [tag ': lower limit against raw draws']);
        localVerifyNum(testCase, valueAt(k,'dep_ul'), R.ul(kk), 1e-12, [tag ': upper limit against raw draws']);
        localVerifyNum(testCase, valueAt(k,'trls1'), R.trls1(kk), 0, [tag ': trls1 against rel.data']);
        localVerifyNum(testCase, valueAt(k,'trls2'), R.trls2(kk), 0, [tag ': trls2 against rel.data']);
        testCase.verifyEqual(G107ExportFixtures.parseFlag(valueAt(k,'ind2include')), ...
            double(R.include(kk)), [tag ': flag against the reference']);
    end
end
end

function localVerifyHeaderText(testCase, headText, pd, what)
%A headless difference file's header: every provenance line, the admissibility
%note when the summary carries one, and the two difference notes.
prov = psyrat_provenance_lines(pd);
for k = 1:numel(prov)
    testCase.verifyTrue(contains(headText, prov{k}), sprintf('%s: header must carry "%s"', what, prov{k}));
end
if isfield(pd.relsummary, 'admissibility_note') && ~isempty(pd.relsummary.admissibility_note)
    testCase.verifyTrue(contains(headText, pd.relsummary.admissibility_note), ...
        [what ': header must carry the admissibility note']);
end
lines = strsplit(headText, newline);
testCase.verifyTrue(localAnyLine(lines, {'relative error variance', 'own', 'always'}), ...
    [what ': header must say the difference error variance is always relative, at the participant''s own counts']);
testCase.verifyTrue(localAnyLine(lines, {'SEM', 'absolute error', 'relative error'}), ...
    [what ': header must say the SEM uses the selected coefficient''s error']);
testCase.verifyTrue(localAnyLine(lines, {'derived', 'table 6'}), ...
    [what ': header must keep the derived (not Rocha Table 6) ICC disclosure']);
end

function [headText, vals] = localReadHeadless(testCase, file, names)
%Reopen a headless table file independently of the writer: CSV through the
%RFC 4180 reader, XLSX through readcell. Returns the header block as text and
%the data rows as an N x numel(NAMES) cell of raw values.
[~, ~, ext] = fileparts(file);
if strcmpi(ext, '.csv')
    recs = G107ExportFixtures.parseCsv(fileread(file));
    h = find(cellfun(@(r) isequal(r, names), recs), 1);
    testCase.assertNotEmpty(h, sprintf('%s: the table''s column row was not found', file));
    head = recs(1:h - 1);
    headText = strjoin(cellfun(@(r) strjoin(r, ','), head, 'UniformOutput', false), newline);
    body = recs(h + 1:end);
    vals = cell(numel(body), numel(names));
    for i = 1:numel(body)
        testCase.verifyNumElements(body{i}, numel(names), sprintf('%s: data row %d field count', file, i));
        m = min(numel(body{i}), numel(names));
        vals(i, 1:m) = body{i}(1:m);
    end
else
    raw = readcell(file);
    h = [];
    for r = 1:size(raw, 1)
        if size(raw, 2) >= numel(names) && all(cellfun(@(x, y) (ischar(x) || isstring(x)) && ...
                strcmp(char(x), y), raw(r, 1:numel(names)), names))
            h = r;
            break;
        end
    end
    testCase.assertNotEmpty(h, sprintf('%s: the table''s column row was not found', file));
    head = raw(1:h - 1, :);
    mask = cellfun(@(x) ischar(x) || isstring(x), head);
    headText = strjoin(reshape(cellstr(string(head(mask))), 1, []), newline);
    vals = raw(h + 1:end, 1:numel(names));
end
end

function h = localColumnRow(raw, names)
%The first row of a readcell result whose leading cells are exactly NAMES.
h = [];
if size(raw, 2) < numel(names)
    return;
end
for r = 1:size(raw, 1)
    if all(cellfun(@(x, y) (ischar(x) || isstring(x)) && strcmp(char(x), y), ...
            raw(r, 1:numel(names)), names))
        h = r;
        return;
    end
end
end

function out = localHeaderTexts(cells)
%Header lines as text, in order, without blank lines and without the timestamp
%line (its time differs between two writes).
out = {};
for i = 1:numel(cells)
    s = localText(cells(i));
    if ~isempty(strtrim(s)) && ~startsWith(s, 'Table Generated on')
        out{end+1} = s; %#ok<AGROW>
    end
end
end

function localVerifyRefusal(testCase, errorText, dialogText, groupName, what)
%A refusal: an error or a dialog, naming the re-summary route (and the group).
msg = strtrim([errorText newline dialogText]);
testCase.verifyTrue(~isempty(errorText) || ~isempty(strtrim(dialogText)), ...
    [what ': must raise an error, never export silently or fall back to event rows']);
testCase.verifyTrue(contains(lower(msg), 'summar'), ...
    [what ': the message must name the re-summary route']);
if ~isempty(groupName)
    testCase.verifyTrue(contains(msg, groupName), [what ': the message must name ' groupName]);
end
end

function localVerifyPrinted(testCase, s, expected, minDecimals, msg)
%A printed number against EXPECTED within the precision it was printed at,
%which must be at least MINDECIMALS places. NaN must be written as NaN (or an
%empty field), never dropped.
v = G107ExportFixtures.parseNum(s);
if isnan(expected)
    testCase.verifyTrue(isnan(v), [msg ': expected NaN']);
    return;
end
d = G107ExportFixtures.decimals(s);
testCase.verifyGreaterThanOrEqual(d, minDecimals, [msg ': printed with too few decimals']);
testCase.verifyLessThanOrEqual(abs(v - expected), 0.5 * 10 ^ (-d) * (1 + 1e-9) + 1e-12, msg);
end

function localVerifyNum(testCase, x, expected, tol, msg)
v = G107ExportFixtures.parseNum(x);
if isnan(expected)
    testCase.verifyTrue(isnan(v), [msg ': expected NaN (kept, not dropped)']);
elseif tol == 0
    testCase.verifyEqual(v, double(expected), msg);
else
    testCase.verifyEqual(v, double(expected), 'AbsTol', tol, msg);
end
end

function tf = localWidthOk(x, ci)
%The interval width as a proportion (0.8), or as a percent only when the field
%says so ('80%'); a bare 80 under 'Interval Width' would be ambiguous.
s = strtrim(char(string(x)));
v = G107ExportFixtures.parseNum(s);
if contains(s, '%')
    tf = abs(v - 100 * ci) <= 1e-6;
else
    tf = abs(v - ci) <= 1e-9;
end
end

function tf = localHasEmptyStatus(lines, groupName)
%An explicit empty-status line (naming the group when one is given).
tf = false;
for i = 1:numel(lines)
    s = lines{i};
    if ~isempty(regexpi(s, 'empty|zero rows|zero-row|no rows|0 rows', 'once')) && ...
            (isempty(groupName) || contains(s, groupName))
        tf = true;
        return;
    end
end
end

function tf = localAnyLine(lines, needles)
%True when one line contains every needle (case-insensitive).
tf = false;
for i = 1:numel(lines)
    s = lower(lines{i});
    if all(cellfun(@(n) contains(s, lower(n)), needles))
        tf = true;
        return;
    end
end
end

function s = localText(x)
%Text from a CSV field, an XLSX cell or a table element ('' for missing).
if iscell(x)
    if isempty(x)
        s = '';
        return;
    end
    x = x{1};
end
if isa(x, 'missing') || (isstring(x) && any(ismissing(x)))
    s = '';
elseif isnumeric(x) || islogical(x)
    s = num2str(x);
else
    s = char(x);
end
end
