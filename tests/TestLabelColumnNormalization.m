classdef TestLabelColumnNormalization < PsyRATTestBase
%Offline guard for finding B12: a label column that is not cellstr must be
%converted at the input boundary, before it can reach the estimation core, and
%missing label values must be rejected rather than quietly filled in.
%
%WHY THESE ARE BOUNDARY TESTS RATHER THAN END-TO-END ONES. The failure B12
%describes happens at REL.out.elabels(:,end+1) = REL.diff_names and its
%siblings, which run only AFTER CmdStan has finished sampling. There is no
%CmdStan-free path to those lines, so the crash site itself is unreachable
%from this lane. The guard therefore sits where the conversion happens and is
%testable, and the one fragile-core call site is pinned by a source scan --
%the tactic B13/B14 established for lines the offline lane cannot execute.

    methods (Test)

        % ---------- conversion, by input type ----------

        function testStringLabelColumnsBecomeCellstr(testCase)
            % The originally observed defect: readtable(f,'TextType','string')
            % produces string label columns, which survive every layer and then
            % fail to assign into the {}-initialized cell arrays in the core.
            tbl = localStringLabelTable();

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyClass(out.id, 'cell');
            testCase.verifyClass(out.group, 'cell');
            testCase.verifyClass(out.event, 'cell');
            testCase.verifyClass(out.time, 'cell');

            % Values must survive the conversion unchanged, not merely the type.
            testCase.verifyEqual(out.id{1}, 's1');
            testCase.verifyEqual(out.event{2}, 'E2');
            testCase.verifyEqual(unique(out.group), {'g1'; 'g2'});
            testCase.verifyEqual(height(out), height(tbl));
        end

        function testStringLabelsCanThenAppendIntoACellArray(testCase)
            % The decisive assertion, and the reason this class exists: it
            % reproduces the exact operation that killed a completed four-chain
            % fit, rather than only checking a class name.
            %
            % Verified in both directions. On the RAW column this append throws
            % MATLAB:UnableToConvert, which is the production failure; after
            % normalization it succeeds.
            tbl = localStringLabelTable();

            rawnames = unique(tbl.event(:), 'stable');
            testCase.verifyError(@() localAppend({}, rawnames), ...
                'MATLAB:UnableToConvert', ...
                ['A raw string event column must still fail this append. If ' ...
                'it does not, MATLAB semantics changed and the second half ' ...
                'of this test no longer proves anything.']);

            out = psyrat_normalize_label_columns(tbl);
            appended = localAppend({}, unique(out.event(:), 'stable'));
            testCase.verifyClass(appended, 'cell');
            testCase.verifyEqual(appended(:), {'E1'; 'E2'});
        end

        function testCategoricalLabelColumnsBecomeCellstr(testCase)
            tbl = localStringLabelTable();
            tbl.id    = categorical(cellstr(tbl.id));
            tbl.group = categorical(cellstr(tbl.group));

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyClass(out.id, 'cell');
            testCase.verifyEqual(out.id{1}, 's1');
            testCase.verifyEqual(unique(out.group), {'g1'; 'g2'});
        end

        function testDatetimeLabelColumnsConvertAndKeepChronologicalOrder(testCase)
            % THE HIGH-STAKES CASE. An occasion column holding session dates
            % arrives as datetime through the ORDINARY file path, because bare
            % readtable auto-detects date-like columns -- no unusual options
            % needed. So this is not an exotic input.
            %
            % The ordering half is the subtle part. Occasion levels are built
            % with unique(), which sorts. MATLAB's default datetime text form is
            % day-first ('20-Jan-2024'), and sorting THAT text orders by
            % day-of-month: 20-Jan-2024 and 15-Jun-2024 come back reversed
            % relative to the dates. Converting to ISO 8601 makes the
            % lexicographic order agree with the chronological one.
            %
            % SCOPE OF THE CLAIM: this protects labels, not coefficients.
            % Occasion is an exchangeable random facet here (psyrat_rel_trt
            % takes only a scalar nocc; REL.time is read only through numel();
            % the Stan occasion effect is i.i.d.), so a permutation does not
            % move a variance component. It changes which session a reader sees
            % as "Occasion 1".
            tbl = table( ...
                [datetime(2024,1,20); datetime(2024,6,15)], [1.1; 2.1], ...
                'VariableNames', {'time','meas'});

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyClass(out.time, 'cell');
            testCase.verifyEqual(out.time, {'2024-01-20'; '2024-06-15'});

            % The property that actually matters: sorting the labels reproduces
            % chronological order. Under the default day-first form it does not.
            testCase.verifyEqual(unique(out.time), {'2024-01-20'; '2024-06-15'}, ...
                ['Sorting the converted occasion labels must reproduce ' ...
                'chronological order. If this fails the labels are no longer ' ...
                'ISO 8601 and unique() can silently reorder occasions.']);

            defaulttext = cellstr(tbl.time);
            testCase.verifyNotEqual(unique(defaulttext), defaulttext, ...
                ['Control: MATLAB''s default datetime text really does sort ' ...
                'out of chronological order for these two dates. If this ' ...
                'stops being true the test above proves less than it claims.']);
        end

        function testDatetimeWithATimeOfDayKeepsTheTime(testCase)
            % Two sessions on the SAME DAY must stay two occasions. Dropping the
            % time for readability would merge them into one level.
            tbl = table( ...
                [datetime(2024,1,20,9,30,0); datetime(2024,1,20,14,0,0)], ...
                [1.1; 2.1], 'VariableNames', {'time','meas'});

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyEqual(out.time, ...
                {'2024-01-20 09:30:00'; '2024-01-20 14:00:00'});
            testCase.verifyNumElements(unique(out.time), 2, ...
                'Two same-day sessions must remain two distinct occasions.');
        end

        function testDurationLabelColumnsBecomeCellstr(testCase)
            tbl = table(duration({'00:30:00'; '01:00:00'}), [1.1; 2.1], ...
                'VariableNames', {'time','meas'});

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyClass(out.time, 'cell');
            testCase.verifyEqual(out.time{1}, '00:30:00');
        end

        function testCharMatrixIsOneLabelPerRowNotPerCharacter(testCase)
            % A char matrix holds one label per ROW. Linearizing it first would
            % turn ['ab';'cd'] into four single-character labels, silently
            % inventing levels -- so the conversion must not use (:) here.
            tbl = table(char('ab','c'), [1.1; 2.1], ...
                'VariableNames', {'id','meas'});

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyNumElements(out.id, 2);
            testCase.verifyEqual(out.id, {'ab'; 'c'});
        end

        function testNumericLabelColumnsPassThroughUntouched(testCase)
            % DELIBERATE CARVE-OUT, and it is load-bearing rather than an
            % oversight. This is the DEFAULT call, which is what psyrat_loadfile
            % uses above its which* row filters. Those filters branch on the
            % type of the user's selection and compare a numeric column with
            % `dataout.group == spec`. Converting a numeric column to text
            % before that point turns a comparison that works today into an
            % error, so numeric must survive the default call unchanged.
            % psyrat_loadfile converts numerics at a SECOND call at the end of
            % the function, once the filters are done with them; the behavioral
            % guard for that split is
            % testLoadFileKeepsTheNumericWhichFilterWorking.
            tbl = table([1; 10], [1.1; 2.1], 'VariableNames', {'id','meas'});

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyClass(out.id, 'double');
            testCase.verifyEqual(out.id, tbl.id, ...
                ['Numeric label columns must pass through untouched, or the ' ...
                'numeric branch of psyrat_loadfile''s which* filters breaks.']);
        end

        function testNumericMissingPassesTheDefaultCallUntouched(testCase)
            % Follows from the carve-out above: since numeric is not converted
            % on the DEFAULT call, its missing values are not policed there
            % either -- the NaN survives as a NaN for the which* filters to
            % compare, and is caught at the converting call instead.
            %
            % DO NOT read this as "a NaN label is accepted". It is not, on any
            % path. Both psyrat_loadfile's second call and psyrat_computevarcomp
            % pass convertnumeric = true and raise labels:missingvalues; see
            % testNumericMissingIsRejectedWhenConverting for the helper and
            % testLoadFileRejectsAMissingNumericLabel for the loader. Until
            % 2026-08-14 the loader really did turn a NaN into a level named
            % literally 'NaN'; that is the behavior those tests replaced.
            tbl = table([1; NaN], [1.1; 2.1], 'VariableNames', {'id','meas'});

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyEqual(out.id, tbl.id);
        end

        % ---------- the engine opt-in (finding B12's numeric residual) ----------

        function testNumericLabelColumnsConvertWhenOptedIn(testCase)
            % The engine call passes convertnumeric = true, because by then the
            % which* filters have run and a numeric label column is a live
            % defect rather than a working path.
            tbl = table([1; 10], [1.1; 2.1], 'VariableNames', {'id','meas'});

            out = psyrat_normalize_label_columns(tbl, true);

            testCase.verifyTrue(iscellstr(out.id), ...
                'the engine call must convert numeric label columns');
            testCase.verifyEqual(out.id, {'1'; '10'}, ...
                ['the spelling must be cellstr(string(...)), i.e. UNPADDED. ' ...
                'cellstr(num2str(...)) would give {'' 1'';''10''}, whose ' ...
                'leading space changes sort order and breaks string matches ' ...
                'against unpadded labels produced elsewhere in the pipeline.']);
        end

        function testNumericLabelCollisionIsRefused(testCase)
            % A CORRECTION TO THIS FILE'S OWN EARLIER REASONING, kept as a test
            % so it cannot be quietly re-adopted. An earlier draft asserted that
            % cellstr(string(...)) preserves full precision where
            % cellstr(num2str(...)) does not, and used that as the deciding
            % argument for the unpadded spelling. IT IS FALSE. Measured
            % 2026-08-12: string() uses the same short format as num2str, so
            % string(pi) is "3.1416" and BOTH render 0.123456789 and
            % 0.123456111 as '0.12346'.
            %
            % So the round trip is lossy under either spelling, and two
            % distinct levels can silently merge into one -- pooling two
            % groups, events or occasions and changing every variance component
            % computed from them. That is refused rather than assumed away.
            % (The real, and only verified, argument for the unpadded spelling
            % is padding: ' 1' against '10' changes sort order and breaks
            % string matches. See testNumericLabelColumnsConvertWhenOptedIn.)
            tbl = table([0.123456789; 0.123456111], [1.1; 2.1], ...
                'VariableNames', {'id','meas'});

            testCase.verifyError(@() psyrat_normalize_label_columns(tbl, true), ...
                'labels:numericlabelcollision', ...
                ['two numeric labels that render as the same text must be ' ...
                'refused, not silently pooled into one level.']);
        end

        function testIntegerLabelsRoundTripWithoutCollision(testCase)
            % POSITIVE CONTROL for the collision guard. Without it, the guard
            % could be rejecting everything and the test above would still
            % pass. Integers are the overwhelmingly common case and must be
            % completely unaffected, including wide ones where the padded
            % spelling would have differed.
            tbl = table([1; 10; 1000000], [1.1; 2.1; 3.1], ...
                'VariableNames', {'id','meas'});

            out = psyrat_normalize_label_columns(tbl, true);

            testCase.verifyEqual(out.id, {'1'; '10'; '1000000'}, ...
                'integer labels must round-trip exactly and unpadded');
        end

        function testNumericMissingIsRejectedWhenConverting(testCase)
            % The missing check has to run BEFORE the numeric conversion, and
            % the reason is sharper here than for the other types: string(NaN)
            % is the literal "NaN", which preflight's missing-cell guard cannot
            % see. Since the engine call runs BEFORE preflight, converting
            % first would turn a caught missing label into an accepted level
            % named NaN.
            tbl = table([1; NaN], [1.1; 2.1], 'VariableNames', {'id','meas'});

            testCase.verifyError(@() psyrat_normalize_label_columns(tbl, true), ...
                'labels:missingvalues', ...
                ['a missing numeric label must be rejected on the converting ' ...
                'call, not silently become the literal label ''NaN''.']);
        end

        function testNumericGroupLabelsProduceNonEmptyMasks(testCase)
            % THE DEFECT ITSELF, at the observable that actually moved.
            %
            % Before this fix the engine turned a numeric unique() into a CHAR
            % MATRIX with bare num2str, then linear-indexed it while the data
            % column stayed numeric, so ismember compared label values against
            % character CODE POINTS and every group mask came back EMPTY
            % (measured: NG1 = NG2 = 0, for same-width labels as well as
            % mixed-width). Asserting on the label TEXT alone would not have
            % caught that, because the text looked plausible; the mask is the
            % observable that distinguishes fixed from broken.
            %
            % THIS TEST IS DIFFERENTIAL, and it has to be. An earlier draft ran
            % the mask arithmetic only on the CONVERTED column, where
            % ismember(cellstr,cellstr) returning non-empty is a MATLAB
            % tautology -- it would have passed with the fix reverted, which is
            % exactly the vacuity this file's own header warns about. So the
            % OLD arithmetic is reproduced here as the negative leg: it must
            % produce empty masks, which is what makes the positive leg mean
            % something.
            %
            % The engine itself needs CmdStan, so the mask line from case 2 is
            % reproduced rather than called.
            % Participants are NESTED in groups (1 and 2 in group 1, 3 in group
            % 10), which is what the case-2 design actually requires -- an
            % earlier draft split participant 3 across both groups, which no
            % real design permits and which made the fixture contradict the
            % arithmetic it claims to reproduce. Widths differ (1 against 10)
            % because that is the case the padded spelling renders unequally.
            tbl = table([1;1;2;2;3;3], [1;1;1;1;10;10], ...
                [1.1;1.2;2.1;2.2;3.1;3.2], 'VariableNames', {'id','group','meas'});

            % NEGATIVE LEG: the pre-fix arithmetic, on the raw numeric column.
            oldnames = unique(tbl.group(:));
            oldnames = num2str(oldnames);          % char MATRIX, as before
            oldtotal = 0;
            for ii = 1:size(oldnames,1)
                oldtotal = oldtotal + sum(ismember(tbl.group, oldnames(ii)));
            end
            testCase.assertEqual(oldtotal, 0, ...
                ['negative leg: the pre-fix arithmetic must still produce ' ...
                'EMPTY masks. If this stops holding, the defect this test ' ...
                'guards has changed shape and the positive leg below proves ' ...
                'nothing.']);

            % POSITIVE LEG: the same arithmetic after normalization.
            out = psyrat_normalize_label_columns(tbl, true);
            groupnames = unique(out.group(:));

            testCase.assertNumElements(groupnames, 2, ...
                'the fixture must yield exactly two group levels');

            counts = zeros(1,numel(groupnames));
            for ii = 1:numel(groupnames)
                counts(ii) = sum(ismember(out.group, groupnames(ii)));
                testCase.verifyGreaterThan(counts(ii), 0, ...
                    sprintf(['group mask %d is empty. This is finding B12''s ' ...
                    'numeric residual: the mask arithmetic compares labels ' ...
                    'against character code points unless the column was ' ...
                    'converted to cellstr first.'], ii));
            end

            % Exactly one mask per row, not merely the right total: a row
            % matching two masks and another matching none would sum correctly.
            % Order follows sorted unique on the converted labels: '1' then
            % '10'. Asserting the counts elementwise (not just their sum) is
            % what rules out a row matching two masks while another matches
            % none, which a total-only check would accept.
            testCase.verifyEqual(counts, [4 2], ...
                'each group must claim exactly its own rows');
        end

        % ---------- the deliberate carve-out ----------

        function testCellstrColumnsPassThroughCompletelyUntouched(testCase)
            % The overwhelmingly common case, and a strict no-op: every golden,
            % baseline and shipped result was produced through cellstr labels.
            tbl = localCellstrLabelTable();

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyEqual(out, tbl, ...
                ['A table whose label columns are already cellstr must come ' ...
                'back byte-identical. Anything else would move existing ' ...
                'labels, and with them shipped output.']);
        end

        function testAnEmptyStringInsideACellstrColumnIsNotRejected(testCase)
            % DELIBERATE, and the reason is redundancy rather than safety: an
            % empty char inside a cellstr column is ALREADY rejected by
            % psyrat_preflight_validate's table-wide 'varargin:missingcells'
            % (measured on a cellstr group column with '' entries). Checking it
            % here too would change which error id the user sees, not whether
            % the run stops.
            %
            % This test pins the SCOPE of the helper's check, not a claim that
            % empty labels are acceptable input. They are not:
            % psyrat_relsummary_groupdata treats '' as "no group" and skips
            % group subsetting, so an empty group label that slipped past
            % preflight would silently return every row for that group.
            tbl = localCellstrLabelTable();
            tbl.group{2} = '';

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyEqual(out, tbl);
        end

        function testNonLabelColumnsAreNeverTouched(testCase)
            % Only the four label columns are in scope. meas is numeric data;
            % dim1/dim2 are dynamic-reliability covariates; weight is the
            % splits n_i. Converting any of them would break the I/O contract.
            tbl = localCellstrLabelTable();
            tbl.dim1   = [0.1; 0.2; 0.3; 0.4];
            tbl.weight = [10; 10; 12; 12];

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyEqual(out.meas, tbl.meas);
            testCase.verifyEqual(out.dim1, tbl.dim1);
            testCase.verifyEqual(out.weight, tbl.weight);
        end

        function testOptionalLabelColumnsMayBeAbsent(testCase)
            % group, event and time are optional. A one-facet design carries
            % none of them, so the function must skip rather than dereference.
            tbl = table(["s1"; "s2"], [1.1; 2.1], 'VariableNames', {'id','meas'});

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyClass(out.id, 'cell');
            testCase.verifyEqual(out.Properties.VariableNames, {'id','meas'});
        end

        function testZeroRowTableIsHandled(testCase)
            tbl = table(strings(0,1), zeros(0,1), 'VariableNames', {'id','meas'});

            out = psyrat_normalize_label_columns(tbl);

            testCase.verifyClass(out.id, 'cell');
            testCase.verifyEmpty(out.id);
        end

        % ---------- missing values are rejected, not filled in ----------

        function testMissingLabelValuesAreRejectedForEveryConvertedType(testCase)
            % WHY AN ERROR RATHER THAN A CONVERSION. cellstr() maps a <missing>
            % string to '' and an <undefined> categorical to '<undefined>', so
            % converting first turns a missing label into an odd-looking level
            % instead of a reported problem.
            %
            % WHAT THIS IS NOT. It is not the only guard: psyrat_preflight_
            % validate already rejects any missing cell table-wide with
            % 'varargin:missingcells', and it catches both the raw <missing>
            % and the '' that conversion produces (measured, both cases). So
            % this does not stop a wrong number that would otherwise be
            % computed -- it produces an earlier error that names the column
            % and the row, where the preflight message locates nothing.
            % Stated plainly so nobody later cites this test as proof that
            % missing labels were previously unguarded.
            cases = { ...
                'string',      [  "s1";  missing]; ...
                'categorical', categorical({'s1'; ''}); ...
                'datetime',    [datetime(2024,1,20); NaT]};

            for k = 1:size(cases,1)
                tbl = table(cases{k,2}, [1.1; 2.1], ...
                    'VariableNames', {'id','meas'});
                testCase.verifyError( ...
                    @() psyrat_normalize_label_columns(tbl), ...
                    'labels:missingvalues', ...
                    sprintf(['A missing label in a %s column must be ' ...
                    'rejected, not converted.'], cases{k,1}));
            end
        end

        function testTheMissingValueErrorNamesTheColumnAndTheRow(testCase)
            % An error a user cannot act on is barely better than the crash it
            % replaced, so the message must locate the problem.
            tbl = table({'s1'; 's2'}, [1.1; 2.1], ["g1"; missing], ...
                'VariableNames', {'id','meas','group'});

            % try/catch rather than verifyError, because verifyError returns a
            % 'missing' placeholder rather than the MException, so the message
            % cannot be inspected through it.
            msg = '';
            try
                psyrat_normalize_label_columns(tbl);
            catch ME
                testCase.verifyEqual(ME.identifier, 'labels:missingvalues');
                msg = ME.message;
            end
            testCase.assertNotEmpty(msg, 'Expected an error to be thrown.');

            testCase.verifySubstring(msg, 'group', ...
                'The error must name the offending column.');
            testCase.verifySubstring(msg, 'row 2', ...
                'The error must locate the first offending row.');
        end

        function testAnUnconvertibleColumnTypeIsRejectedClearly(testCase)
            % A cell that is not cellstr has no text form. It must fail with a
            % named error rather than whatever cellstr() happens to throw.
            tbl = table({1; 2}, [1.1; 2.1], 'VariableNames', {'id','meas'});

            msg = '';
            try
                psyrat_normalize_label_columns(tbl);
            catch ME
                testCase.verifyEqual(ME.identifier, 'labels:unsupportedtype');
                msg = ME.message;
            end
            testCase.assertNotEmpty(msg, 'Expected an error to be thrown.');
            testCase.verifySubstring(msg, 'id', ...
                'The error must name the offending column.');
        end

        % ---------- the two production call sites ----------

        function testLoadFileNormalizesStringLabelColumns(testCase)
            % INTEGRATION at the production boundary. All three production entry
            % points reach the estimation core through psyrat_loadfile -- the GUI
            % (psyrat_startproc), the CLI (psyrat_run) and the warp fallback
            % (psyrat_computevarcompwarp) -- so normalizing there covers each.
            %
            % The in-memory 'dataraw' form is used because it is what the
            % affected CLI call psyrat_run('data', readtable(f,'TextType',
            % 'string'), ...) reduces to: psyrat_run stores the caller's table
            % verbatim and hands that same table to the loader. 'file' is
            % supplied only for the table's Description; it is never read.
            tbl = localStringLabelTable();

            out = psyrat_loadfile('file', 'in_memory.csv', 'dataraw', tbl, ...
                'idcol', 'id', 'meascol', 'meas', ...
                'groupcol', 'group', 'eventcol', 'event', 'timecol', 'time');

            testCase.verifyClass(out.id, 'cell');
            testCase.verifyClass(out.group, 'cell');
            testCase.verifyClass(out.event, 'cell');
            testCase.verifyClass(out.time, 'cell');
            testCase.verifyEqual(out.id{1}, 's1');
        end

        function testLoadFileNormalizesADatetimeOccasionColumnFromAPlainCsv(testCase)
            % THE REGRESSION TEST FOR THE ORDINARY PATH. No special options: a
            % plain CSV with dated sessions, read the way psyrat_readtable reads
            % it. Before the widened guard this came back as datetime and went
            % on to fail after sampling.
            file = localWriteCsv(table( ...
                {'s1'; 's1'; 's2'; 's2'}, [1.1; 1.2; 2.1; 2.2], ...
                {'2024-01-20'; '2024-06-15'; '2024-01-20'; '2024-06-15'}, ...
                'VariableNames', {'id', 'meas', 'time'}));
            testCase.addTeardown(@() delete(file));

            % Control: confirm the premise is still true on this MATLAB, so a
            % future readtable that stops auto-detecting cannot make this test
            % pass vacuously.
            testCase.assumeClass(readtable(file).time, 'datetime', ...
                'Premise: bare readtable auto-detects a date column as datetime.');

            out = psyrat_loadfile('file', file, 'idcol', 'id', ...
                'meascol', 'meas', 'timecol', 'time');

            testCase.verifyClass(out.time, 'cell');
            testCase.verifyEqual(unique(out.time), {'2024-01-20'; '2024-06-15'});
        end

        function testLoadFileNormalizesBeforeTheWhichFiltersCompare(testCase)
            % ORDERING REGRESSION TEST, and it caught a real defect: the first
            % version of this fix normalized at the very END of psyrat_loadfile,
            % roughly 175 lines AFTER the which* row filters.
            %
            % Those filters compare the column against the user's selection with
            % strcmpi, and strcmpi returns a scalar false for a datetime or
            % categorical column instead of matching element-wise (measured).
            % So selecting an occasion on a dated column failed with a
            % misleading 'times:occasionmismatch' -- reporting that the occasion
            % does not exist in data that plainly contains it.
            %
            % Both halves matter: the filter must MATCH, and it must actually
            % SUBSET. A normalization that ran late would fail the first.
            file = localWriteCsv(table( ...
                {'s1'; 's1'; 's2'; 's2'}, [1.1; 1.2; 2.1; 2.2], ...
                {'2024-01-20'; '2024-06-15'; '2024-01-20'; '2024-06-15'}, ...
                'VariableNames', {'id', 'meas', 'time'}));
            testCase.addTeardown(@() delete(file));

            testCase.assumeClass(readtable(file).time, 'datetime', ...
                'Premise: bare readtable auto-detects a date column as datetime.');

            out = psyrat_loadfile('file', file, 'idcol', 'id', ...
                'meascol', 'meas', 'timecol', 'time', ...
                'whichtimes', {'2024-01-20'});

            testCase.verifyEqual(height(out), 2, ...
                'The occasion filter must subset to the selected occasion.');
            testCase.verifyEqual(unique(out.time), {'2024-01-20'});
        end

        % ---------- the loader's numeric conversion (B12's loader half) ----------

        function testLoadFileRendersNumericLabelsUnpadded(testCase)
            % THE POINT OF THE 2026-08-12 RULING, and the only test in the suite
            % that can see it. psyrat_loadfile used to convert numeric label
            % columns with four hand-rolled cellstr(num2str(...)) blocks, and
            % num2str RIGHT-ALIGN PADS a column vector to a common width. Every
            % other label producer in the toolbox emits unpadded
            % cellstr(string(...)), and the mismatch broke four cross-source
            % string joins and put a leading space into unquoted CSV exports.
            %
            % MIXED WIDTH IS LOAD-BEARING HERE. Every pre-existing fixture that
            % reaches the loader with a numeric label column uses uniform-width
            % values (TestFileFunctions' [1;1;2], testdata.csv's 4-digit subjid),
            % where num2str pads by zero characters and the two spellings are
            % byte-identical. That is why the whole suite stayed green through
            % this change and why 10 has to be in this fixture.
            file = localWriteCsv(table([1; 1; 2; 2; 10; 10], ...
                [1.1; 1.2; 2.1; 2.2; 3.1; 3.2], ...
                'VariableNames', {'id', 'meas'}));
            testCase.addTeardown(@() delete(file));

            testCase.assumeClass(readtable(file).id, 'double', ...
                'Premise: the id column must arrive numeric for this to test anything.');

            out = psyrat_loadfile('file', file, 'idcol', 'id', 'meascol', 'meas');

            testCase.verifyClass(out.id, 'cell');
            testCase.verifyEqual(unique(out.id), {'1'; '10'; '2'}, ...
                ['Numeric labels must render UNPADDED. The padded spelling ' ...
                'would give {'' 1'';'' 2'';''10''}, which sorts differently ' ...
                'and does not string-match the labels the rest of the ' ...
                'pipeline produces.']);
            testCase.verifyFalse(any(startsWith(out.id, ' ')), ...
                'No label may carry a leading space.');
        end

        function testLoadFileKeepsTheNumericWhichFilterWorking(testCase)
            % THE CARVE-OUT THIS DESIGN RESTS ON, and nothing else covers it.
            %
            % psyrat_loadfile calls psyrat_normalize_label_columns TWICE: once
            % above the which* row filters with convertnumeric defaulted to
            % false, and once at the very end with it true. The split exists
            % because the filters compare a numeric column with
            % `dataout.group == spec` -- so hoisting the numeric conversion up
            % to the first call turns that comparison into an error on a cell.
            %
            % This is a genuine differential: every other which* test in the
            % suite selects with TEXT, so the numeric branch at
            % psyrat_loadfile's whichgroups block was previously unexecuted by
            % any test and a hoist would have landed green.
            file = localWriteCsv(table([1; 1; 2; 2], [1.1; 1.2; 2.1; 2.2], ...
                [1; 1; 2; 2], 'VariableNames', {'id', 'meas', 'group'}));
            testCase.addTeardown(@() delete(file));

            testCase.assumeClass(readtable(file).group, 'double', ...
                'Premise: the group column must arrive numeric.');

            out = psyrat_loadfile('file', file, 'idcol', 'id', ...
                'meascol', 'meas', 'groupcol', 'group', 'whichgroups', {2});

            testCase.verifyEqual(height(out), 2, ...
                ['A NUMERIC group selection must still subset. If it does ' ...
                'not, the numeric conversion has been hoisted above the ' ...
                'which* filters.']);
            testCase.verifyEqual(unique(out.group), {'2'});
        end

        function testLoadFileRejectsAMissingNumericLabel(testCase)
            % APPROVED RIDER, and a user-visible change to previously-accepted
            % input. A NaN in a numeric label column used to become a level
            % named literally 'NaN': num2str renders it, and preflight's
            % missing-cell guard cannot see a string that is already text. So
            % the run continued with a bogus participant. Routing the loader's
            % conversion through psyrat_normalize_label_columns puts the
            % missing check ahead of the conversion, where it can still see it.
            file = localWriteCsv(table([1; 1; NaN; 2], ...
                [1.1; 1.2; 2.1; 2.2], 'VariableNames', {'id', 'meas'}));
            testCase.addTeardown(@() delete(file));

            testCase.verifyError( ...
                @() psyrat_loadfile('file', file, 'idcol', 'id', 'meascol', 'meas'), ...
                'labels:missingvalues', ...
                'A missing numeric label must be refused, not renamed to ''NaN''.');
        end

        function testLoadFileRejectsCollidingNumericLabels(testCase)
            % THE SECOND APPROVED RIDER. MATLAB renders numbers with about five
            % significant digits by default -- for string() exactly as for
            % num2str -- so 0.123456789 and 0.123456111 both become '0.12346'.
            % unique() would then see ONE level where the data have two,
            % silently pooling two participants into a single stratum and
            % changing every variance component computed from them. The loader
            % used to do this without complaint.
            file = localWriteCsv(table( ...
                [0.123456789; 0.123456789; 0.123456111; 0.123456111], ...
                [1.1; 1.2; 2.1; 2.2], 'VariableNames', {'id', 'meas'}));
            testCase.addTeardown(@() delete(file));

            testCase.verifyError( ...
                @() psyrat_loadfile('file', file, 'idcol', 'id', 'meascol', 'meas'), ...
                'labels:numericlabelcollision', ...
                'Two numeric labels rendering to the same text must be refused.');
        end

        function testLoadFileHasNoHandRolledNumericLabelConversion(testCase)
            % SOURCE-SCANNING GUARD against re-introduction. The behavioral
            % tests above pin the OUTPUT, but a future edit could reasonably
            % re-add a cellstr(num2str(...)) block for a column they do not
            % cover -- the four that were deleted here covered id, group, event
            % and time, and only id is exercised above.
            src = fileread(localLoaderSourcePath());

            % Strip comment lines before matching, the same way
            % testEstimationCoreNormalizesLabelsBeforePreflight does. This file
            % DISCUSSES the padded spelling at length in the comment block above
            % the replacement call, so an unstripped scan would match the very
            % prose explaining why the idiom is gone.
            lines = regexp(src, '\r\n|\n|\r', 'split');
            lines = lines(cellfun(@isempty, regexp(lines, '^\s*%', 'once')));
            body  = strjoin(lines, newline);

            testCase.verifyEmpty(strfind(body, 'num2str'), ...
                ['psyrat_loadfile must not convert labels with num2str. It ' ...
                'pads to a common width, and the padded labels do not ' ...
                'string-match the unpadded ones every other producer emits. ' ...
                'Use psyrat_normalize_label_columns(dataout,true) instead.']);

            % Both call sites must survive, and with the right opt-in. The
            % early one must NOT convert numerics (the which* filters compare
            % them with ==); the late one must.
            calls = strfind(body, 'psyrat_normalize_label_columns(dataout');
            testCase.assertNumElements(calls, 2, ...
                ['psyrat_loadfile must call the normalizer twice: once above ' ...
                'the which* filters for the non-numeric types, and once at ' ...
                'the end for the numeric ones.']);
            testCase.verifyNotEmpty( ...
                strfind(body, 'psyrat_normalize_label_columns(dataout,true)'), ...
                'The late call must opt in to numeric conversion.');
        end

        function testEstimationCoreNormalizesLabelsBeforePreflight(testCase)
            % SOURCE-SCANNING GUARD for the one fragile-core call site.
            %
            % WHY A SCAN. psyrat_computevarcomp cannot be run to completion in
            % this lane -- it needs CmdStan -- so no offline test can observe the
            % call happening. That is the same situation B13/B14 documented, and
            % this is the same remedy.
            %
            % WHY 'BEFORE PREFLIGHT' AND NOT MERELY 'PRESENT'. The normalization
            % has to precede psyrat_preflight_validate and every label read that
            % follows it. A call placed after them would scan green while leaving
            % the offending type in play for exactly the code that breaks on it,
            % so ordering is the property worth pinning, not existence.
            src = fileread(localEngineSourcePath());

            % Strip FULL-LINE comments first, so a comment mentioning the call
            % cannot satisfy the scan. Known limit: '^\s*%' catches full-line
            % comments only -- a trailing comment survives it -- but in MATLAB a
            % line whose first non-blank character is % is always a comment, so
            % this can never over-strip.
            lines = regexp(src, '\r\n|\n|\r', 'split');
            lines = lines(cellfun(@isempty, regexp(lines, '^\s*%', 'once')));
            body  = strjoin(lines, newline);

            % Match WITHOUT the closing paren so the scan survives argument
            % changes. It previously matched '...(datatable)' exactly and broke
            % the moment the numeric opt-in was added, which is a scanner
            % brittleness bug rather than a real regression.
            normIdx = strfind(body, 'psyrat_normalize_label_columns(datatable');
            testCase.assertNotEmpty(normIdx, ...
                ['psyrat_computevarcomp must normalize its label columns. ' ...
                'Without this call a caller that hands the core a table ' ...
                'directly -- which is how the entire test suite and any ' ...
                'programmatic user reaches it, bypassing psyrat_loadfile -- ' ...
                'can still lose a completed fit to B12.']);

            % The engine call must OPT IN to numeric conversion. This is not
            % cosmetic: with the default (false) a numeric label column reaches
            % the core untouched, and the core no longer has the num2str
            % branches that used to (badly) handle it, so REL.groups would be
            % stored numeric and every downstream gnames{i} would error.
            % Whitespace-tolerant on purpose: an exact 'datatable,true' match
            % would fail on the equally valid 'datatable, true', which is the
            % spacing this very file uses elsewhere. Pin the property, not the
            % formatting.
            convIdx = regexp(body, ...
                'psyrat_normalize_label_columns\(\s*datatable\s*,\s*true\s*\)', 'once');
            testCase.verifyNotEmpty(convIdx, ...
                ['The engine must call psyrat_normalize_label_columns with ' ...
                'convertnumeric = true. The loader deliberately passes the ' ...
                'default instead, because its which* row filters compare a ' ...
                'numeric column with == before its own conversion runs.']);

            preIdx = strfind(body, 'psyrat_preflight_validate(');
            testCase.assertNotEmpty(preIdx, ...
                ['Expected a psyrat_preflight_validate call in the engine. It ' ...
                'was renamed or removed - fix this scanner rather than ' ...
                'deleting it.']);

            % And it is the CONVERTING call whose position matters. normIdx
            % above now matches any call, so ordering must be re-pinned here:
            % a one-argument call placed early would satisfy normIdx while the
            % converting call sat after preflight.
            testCase.verifyLessThan(convIdx, preIdx(1), ...
                ['The converting call must run BEFORE preflight, not merely ' ...
                'some call to the normalizer.']);

            testCase.verifyLessThan(normIdx(1), preIdx(1), ...
                ['The label normalization must run BEFORE preflight ' ...
                'validation and before any label column is read. Placed ' ...
                'after, it would leave the offending type in play for the ' ...
                'very code that fails on it.']);
        end

    end

end

function tbl = localStringLabelTable()
%Balanced 2 id x 2 event x 2 occasion design with STRING label columns, which
%is the shape psyrat_loadfile's completeness checks expect (every participant
%has a measurement for each event at each occasion).
tbl = table( ...
    ["s1"; "s1"; "s1"; "s1"; "s2"; "s2"; "s2"; "s2"], ...
    [1.10; 1.20; 1.30; 1.40; 2.10; 2.20; 2.30; 2.40], ...
    ["g1"; "g1"; "g1"; "g1"; "g2"; "g2"; "g2"; "g2"], ...
    ["E1"; "E2"; "E1"; "E2"; "E1"; "E2"; "E1"; "E2"], ...
    ["T1"; "T1"; "T2"; "T2"; "T1"; "T1"; "T2"; "T2"], ...
    'VariableNames', {'id','meas','group','event','time'});
end

function tbl = localCellstrLabelTable()
%The ordinary case: label columns already cellstr.
tbl = table( ...
    {'s1'; 's1'; 's2'; 's2'}, ...
    [1.10; 1.20; 2.10; 2.20], ...
    {'g1'; 'g1'; 'g2'; 'g2'}, ...
    {'E1'; 'E2'; 'E1'; 'E2'}, ...
    'VariableNames', {'id','meas','group','event'});
end

function out = localAppend(target, names)
%Reproduce the estimation core's append verbatim:
%   REL.out.elabels(:,end+1) = REL.diff_names;
%against a cell array initialized as {}. Wrapped in a function so it can be
%handed to verifyError.
target(:, end+1) = names;
out = target;
end

function file = localWriteCsv(tbl)
file = [tempname '.csv'];
writetable(tbl, file);
end

function p = localEngineSourcePath()
%Absolute path to the estimator core, derived from THIS file's location so the
%scan does not depend on pwd or on which copy of the toolbox is first on the
%MATLAB path (a .claude/worktrees copy shadowing the real one has bitten this
%repo before).
here = fileparts(mfilename('fullpath'));          % <root>/tests
p = fullfile(fileparts(here), 'subroutines', 'estimation', 'psyrat_computevarcomp.m');
end

function p = localLoaderSourcePath()
%Absolute path to the loader, derived the same way and for the same reason as
%localEngineSourcePath above: which('psyrat_loadfile') would resolve against
%the MATLAB path, and a .claude/worktrees copy shadowing the real toolbox has
%bitten this repo before.
here = fileparts(mfilename('fullpath'));          % <root>/tests
p = fullfile(fileparts(here), 'subroutines', 'filefunctions', 'psyrat_loadfile.m');
end
