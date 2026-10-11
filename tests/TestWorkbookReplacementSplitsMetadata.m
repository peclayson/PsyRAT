classdef TestWorkbookReplacementSplitsMetadata < PsyRATTestBase
    %The multi-occasion data-splits viewer's two Save buttons (G110), pinned on
    %the design metadata they export. Each stratum's design is written as JSON
    %twice: in a header line ('Design <stratum>: <json>') and in a column of
    %every table row. The JSON grows with the design, and XLSX keeps at most
    %32,767 characters in a cell, so a large design must be refused in XLSX
    %and kept whole in CSV.
    %
    %THE LARGE DESIGN is the one the probe of 2026-09-23 used: 60 participants,
    %2 occasions, 10 splits in every participant-occasion cell (1,200 rows),
    %built on the factory's multi-occasion result. Its header line is 33,671
    %characters. Only the design matters here: the coefficients come from the
    %factory's scalar SDs and are not an estimate, and nothing here checks a
    %coefficient. THE SMALL DESIGN is the factory's own 24 rows.
    %
    %EXPECTATIONS are production parity: the JSON the summary itself encodes
    %(psyrat_splits_summary with the viewer's default fields, then jsonencode)
    %and the table the viewer's own builders make, never a re-derivation. The
    %design's 36 top-level fields, in order, are the export contract and are
    %pinned separately (localDesignFields), so a change to them is seen even
    %when the summary and the export change together.
    %
    %Three checks for each button: the large design's XLSX save is refused with
    %one error and leaves the chosen file byte-identical; its CSV save keeps
    %every design JSON whole (RFC 4180 parse); the small design's XLSX save
    %succeeds and reads back exactly (positive control).
    %
    %No fit, no CmdStan, no desktop.

    properties (TestParameter)
        %the viewer's two Save buttons
        Button = struct('coefficients', 'coef', 'varianceTable', 'varcomp');
    end

    methods (Test)

        function testLongDesignXlsxIsRefusedWithOneError(testCase, Button)
            %The large design through the real Save button, as .xlsx over an
            %existing workbook. Required: the save is refused before writing,
            %naming the XLSX text limit; errordlg is called once; the chosen
            %workbook and the folder are unchanged; and the "Saved" line is not
            %printed.
            env = localEnv(testCase);
            pd = G110WorkbookHarness.withHeader(localLargeDesign());
            x = localExpected(pd, Button);
            testCase.assertNumElements(x.headerLines, 1, 'precondition: one stratum with design metadata');
            testCase.assertGreaterThan(max(cellfun(@numel, x.headerLines)), 32767, ...
                'precondition: the large design''s header line must pass the XLSX cell limit');
            fig = localOpen(testCase, pd);
            name = G110WorkbookHarness.TargetName;
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            original = G110WorkbookHarness.fileBytes(target);
            G110WorkbookHarness.armErrordlgRecorder(env);
            r = G110WorkbookHarness.fireSave(fig, x.label, target);
            calls = G110WorkbookHarness.disarmErrordlgRecorder(env);
            testCase.assertEqual(r.saveCalls, 1, 'harness: the button did not ask the save dialog once');
            testCase.verifyTrue(startsWith(r.errText, 'psyrat:replaceWorkbook:refused'), ...
                ['the XLSX save must be refused before writing: ' r.errText]);
            testCase.verifyTrue(contains(r.errText, '32,767') || contains(r.errText, '32767'), ...
                ['the refusal must name the XLSX text limit: ' r.errText]);
            testCase.verifyNumElements(calls, 1, sprintf(['one refusal must show one error dialog; ', ...
                'errordlg was called %d times'], numel(calls)));
            testCase.verifyTrue(isequal(G110WorkbookHarness.fileBytes(target), original), ...
                'the chosen workbook keeps its bytes');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, true);
            testCase.verifyEmpty(p, ['the chosen file and the folder must be unchanged: ' localJoin(p)]);
            testCase.verifyFalse(contains(r.out, 'Saved '), ...
                ['a refused save must not print its "Saved" line: ' r.out]);
        end

        function testCsvKeepsLongDesignMetadataWhole(testCase, Button)
            %The large design through the real Save button, as .csv. Required:
            %no error and no error dialog, the "Saved" line printed, and, after
            %an RFC 4180 parse, the design header line of each stratum as one
            %field equal to the summary's, the table's column names, one record
            %per table row with as many fields as the table has columns, and
            %each row's design JSON equal to its stratum's, every JSON decoding
            %to the 36 contract fields.
            env = localEnv(testCase);
            pd = G110WorkbookHarness.withHeader(localLargeDesign());
            x = localExpected(pd, Button);
            testCase.assertNumElements(x.headerLines, 1, 'precondition: one stratum with design metadata');
            testCase.assertGreaterThan(max(cellfun(@numel, values(x.json))), 32767, ...
                'precondition: the large design''s JSON must pass the XLSX cell limit');
            localVerifyContract(testCase, x);
            fig = localOpen(testCase, pd);
            target = fullfile(env.work, 'g110_splits.csv');
            G110WorkbookHarness.armErrordlgRecorder(env);
            r = G110WorkbookHarness.fireSave(fig, x.label, target);
            calls = G110WorkbookHarness.disarmErrordlgRecorder(env);
            testCase.assertEqual(r.saveCalls, 1, 'harness: the button did not ask the save dialog once');
            testCase.assertEmpty(r.errText, ['the CSV save failed: ' r.errText]);
            testCase.verifyEmpty(calls, 'a successful CSV save shows no error dialog');
            testCase.verifyTrue(contains(r.out, 'Saved '), 'the CSV save prints its "Saved" line');
            testCase.assertTrue(isfile(target), 'the CSV file was not written');
            txt = native2unicode(G110WorkbookHarness.fileBytes(target)', 'UTF-8');
            recs = localCsvRecords(txt);
            localVerifyExport(testCase, recs, x, 'CSV');
        end

        function testSmallDesignXlsxSavesExactly(testCase, Button)
            %Positive control: the factory's 24-row design, whose JSON fits a
            %cell, through the same button as .xlsx over an existing workbook.
            %Required: the save succeeds with no error dialog and prints its
            %"Saved" line; the planted workbook is replaced (no planted or
            %stale cell, one sheet) and nothing else changes; and the workbook
            %reads back with the same design header line, column names, rows
            %and design JSON as the CSV check requires.
            env = localEnv(testCase);
            pd = G110WorkbookHarness.withHeader(PsyRATTestDataFactory.makeSplitsTRTSummaryData());
            testCase.assertEqual(height(pd.rel.data), 24, 'precondition: the factory''s 24-row design');
            x = localExpected(pd, Button);
            testCase.assertNumElements(x.headerLines, 1, 'precondition: one stratum with design metadata');
            testCase.assertLessThan(max(cellfun(@numel, x.headerLines)), 32767, ...
                'precondition: the small design''s header line fits an XLSX cell');
            localVerifyContract(testCase, x);
            fig = localOpen(testCase, pd);
            name = G110WorkbookHarness.TargetName;
            target = fullfile(env.work, name);
            G110WorkbookHarness.plantWorkbook(target);
            G110WorkbookHarness.plantDecoys(env.work, name);
            before = G110WorkbookHarness.snapshot(env.work);
            G110WorkbookHarness.armErrordlgRecorder(env);
            r = G110WorkbookHarness.fireSave(fig, x.label, target);
            calls = G110WorkbookHarness.disarmErrordlgRecorder(env);
            testCase.assertEqual(r.saveCalls, 1, 'harness: the button did not ask the save dialog once');
            testCase.assertEmpty(r.errText, ['the XLSX save failed: ' r.errText]);
            testCase.verifyEmpty(calls, 'a successful XLSX save shows no error dialog');
            testCase.verifyTrue(contains(r.out, 'Saved '), 'the XLSX save prints its "Saved" line');
            rep = G110WorkbookHarness.staleReport(target);
            testCase.verifyFalse(rep.plant || rep.near || rep.wide || rep.long || rep.sheet, ...
                'the planted workbook was replaced: no planted or stale cell and no old sheet remain');
            testCase.verifyNumElements(rep.sheets, 1, 'the new workbook has one sheet');
            p = G110WorkbookHarness.compareSnapshots(before, G110WorkbookHarness.snapshot(env.work), name, false);
            testCase.verifyEmpty(p, ['nothing may be left and the decoys untouched: ' localJoin(p)]);
            recs = localRecordsOf(readcell(target));
            localVerifyExport(testCase, recs, x, 'XLSX');
        end
    end
end

% -- local helpers -------------------------------------------------------------

function env = localEnv(testCase)
env = G110WorkbookHarness.install();
testCase.addTeardown(@() G110WorkbookHarness.uninstall(env));
end

function fig = localOpen(testCase, pd)
%Open the data-splits viewer on PD through the harness, and close every window
%this test opened from here on at teardown. A viewer that does not open fails
%this test only.
figsBefore = findall(groot, 'Type', 'figure');
testCase.addTeardown(@() localCloseNewFigures(figsBefore));
fig = [];
try
    fig = G110WorkbookHarness.openRoute(G110WorkbookHarness.route('splits_coef'), pd);
catch ME
    testCase.assertFail(['harness: the data-splits viewer did not open: ' ME.message]);
end
end

function localCloseNewFigures(figsBefore)
figs = findall(groot, 'Type', 'figure');
opened = figs(~ismember(figs, figsBefore));
delete(opened(isvalid(opened)));
end

function pd = localLargeDesign()
%The factory's multi-occasion result with REL.data replaced by the probe's
%60-participant design: for participant p = 1..60, occasion o = 1..2 and split
%s = 1..10, one row with id = 1000 + p, time = o, weight = 5 + mod(7p + 3o + s, 8)
%items (5 to 12) and meas = the row number, in that loop order. Every
%participant has both occasions and 10 splits in each.
pd = PsyRATTestDataFactory.makeSplitsTRTSummaryData();
N = 60;
S = 10;
n = 2 * N * S;
id = zeros(n, 1);
time = zeros(n, 1);
weight = zeros(n, 1);
k = 0;
for p = 1:N
    for o = 1:2
        for s = 1:S
            k = k + 1;
            id(k) = 1000 + p;
            time(k) = o;
            weight(k) = 5 + mod(7 * p + 3 * o + s, 8);
        end
    end
end
meas = (1:n)';
pd.rel.data = table(id, meas, weight, time, 'VariableNames', {'id', 'meas', 'weight', 'time'});
pd.rel.filename = 'g110_splits_p60_s10_balanced.psyrat';
pd.proc = struct('savename', pd.rel.filename);
end

function x = localExpected(pd, key)
%What the button KEY must export for PD, from production code only: the summary
%the viewer builds with its default fields (95% interval, observed split and
%occasion counts), the table the viewer's own builder makes from it, and each
%stratum's design JSON and header line.
REL = pd.rel;
summ = psyrat_splits_summary(REL, 'CI', 0.95, 'obs', [], 'nocc', []);
x = struct();
switch key
    case 'coef'
        x.label = 'Save coefficients';
        [data, cols] = psyrat_splits_coef_rows(summ, psyrat_unitlabels(REL.splits));
        x.T = cell2table(data, 'VariableNames', matlab.lang.makeValidName(cols));
        x.designVar = x.T.Properties.VariableNames{find(strcmp(cols, 'Design metadata'), 1)};
    otherwise
        x.label = 'Save variance table';
        x.T = psyrat_splits_varcomp_concat(summ);
        x.designVar = 'design_metadata';
end
x.json = containers.Map('KeyType', 'char', 'ValueType', 'any');
x.headerLines = {};
x.designs = {};
for t = 1:numel(summ.strata)
    label = char(string(summ.strata(t).label));
    json = jsonencode(summ.strata(t).design);
    x.json(label) = json;
    x.headerLines{end + 1} = sprintf('Design %s: %s', label, json);
    x.designs{end + 1} = summ.strata(t).design;
end
end

function f = localDesignFields()
%The 36 top-level fields of a multi-occasion stratum's design metadata, in the
%order psyrat_splits_summary builds them and jsonencode writes them (probe of
%2026-09-23, and the source: local_trt_design, then the mode, projection and
%reference fields). A change here is a change of the export contract.
f = {'group', 'event', 'participant_count', 'person_ids', 'occasion_ids', ...
    'occasion_union', 'person_occasion_counts', 'cell_person_index', ...
    'cell_occasion_index', 'cell_split_counts', 'row_indices', 'row_person_index', ...
    'row_occasion_index', 'row_cell_index', 'row_split_index', 'row_items', ...
    'row_weights', 'factor_order', 'person_factors', 'mean_factors', ...
    'missing_occasions', 'irregular', 'aggregation', 'split_identity', ...
    'row_residual_assumption', 'fixed_facet_target', 'mode', 'score_rule', ...
    'projection_nsplit', 'projection_nocc', 'coefficient_factors', 'obs_override', ...
    'nocc_override', 'reference_nsplit', 'reference_nocc', 'pooled_harmonic_items'};
end

function localVerifyContract(testCase, x)
%The summary's own design for each stratum has the 36 contract fields, in order.
for t = 1:numel(x.designs)
    testCase.verifyEqual(fieldnames(x.designs{t})', localDesignFields(), sprintf(['the ', ...
        'summary''s design for stratum %d must have the 36 contract fields, in order'], t));
end
end

function localVerifyExport(testCase, recs, x, what)
%The exported design metadata against the summary, for either format. RECS
%holds the CSV records or the sheet's rows, each a row cell of field texts.
names = x.T.Properties.VariableNames;
dcol = find(strcmp(names, x.designVar), 1);
%the column-name record
r0 = [];
for r = 1:numel(recs)
    if isequal(recs{r}, names)
        r0 = r;
        break;
    end
end
testCase.assertNotEmpty(r0, sprintf('%s: no record holds exactly the table''s %d column names', ...
    what, numel(names)));
%the design header lines above it: one field each, equal to the summary's
heads = {};
for r = 1:r0 - 1
    rec = recs{r};
    if ~isempty(rec) && startsWith(rec{1}, 'Design ')
        testCase.verifyTrue(all(cellfun(@isempty, rec(2:end))), ...
            sprintf('%s: the design header line in record %d must be one field', what, r));
        heads{end + 1} = rec{1}; %#ok<AGROW>
    end
end
testCase.verifyEqual(numel(heads), numel(x.headerLines), ...
    sprintf('%s: one design header line per stratum', what));
for t = 1:min(numel(heads), numel(x.headerLines))
    testCase.verifyTrue(strcmp(heads{t}, x.headerLines{t}), sprintf(['%s: design header ', ...
        'line %d differs from the summary''s (%s)'], what, t, localFirstDiff(heads{t}, x.headerLines{t})));
    k = strfind(heads{t}, ': {');
    if ~isempty(k)
        localVerifyDecoded(testCase, heads{t}(k(1) + 2:end), sprintf('%s header line %d', what, t));
    end
end
%the body: one record per table row, as wide as the table, with its stratum's JSON
body = recs(r0 + 1:end);
while ~isempty(body) && all(cellfun(@isempty, body{end}))
    body(end) = [];
end
testCase.verifyEqual(numel(body), height(x.T), sprintf('%s: one record per table row', what));
for i = 1:min(numel(body), height(x.T))
    rec = body{i};
    testCase.verifyEqual(numel(rec), numel(names), sprintf('%s: row %d has %d fields; the table has %d', ...
        what, i, numel(rec), numel(names)));
    if numel(rec) < dcol || ~isKey(x.json, rec{1})
        testCase.verifyFail(sprintf('%s: row %d names no summary stratum or has no design field', what, i));
        continue;
    end
    want = x.json(rec{1});
    testCase.verifyTrue(strcmp(rec{dcol}, want), sprintf(['%s: row %d''s design metadata differs ', ...
        'from the summary''s (%s)'], what, i, localFirstDiff(rec{dcol}, want)));
    localVerifyDecoded(testCase, rec{dcol}, sprintf('%s row %d', what, i));
end
end

function localVerifyDecoded(testCase, json, where)
%One exported design JSON decodes to the 36 contract fields, in order.
try
    got = fieldnames(jsondecode(json))';
catch ME
    got = {['<jsondecode failed: ' ME.message '>']};
end
testCase.verifyEqual(got, localDesignFields(), ...
    sprintf('%s: the design metadata must decode to the 36 contract fields, in order', where));
end

function s = localFirstDiff(a, b)
%Where two long texts first differ, for a diagnostic that does not print them.
n = min(numel(a), numel(b));
k = find(a(1:n) ~= b(1:n), 1);
if isempty(k)
    s = sprintf('lengths %d and %d, equal over the first %d characters', numel(a), numel(b), n);
else
    s = sprintf('lengths %d and %d, first difference at character %d: "%s" vs "%s"', ...
        numel(a), numel(b), k, a(k:min(end, k + 30)), b(k:min(end, k + 30)));
end
end

function recs = localCsvRecords(txt)
%The records of a CSV text per RFC 4180, each a row cell of field texts. A field
%that opens with a double quote runs to the next quote that is not doubled; it
%may hold commas, line breaks and doubled quotes (one literal quote each), and
%only a comma, a line break or the end of the text may follow it. Records end at
%LF or CR LF. An unterminated quote, text after a closing quote, or a bare CR
%outside quotes raises, so a record cannot be split or joined silently. Quote
%and separator positions are found once, so a long field costs one pass.
recs = {};
n = numel(txt);
q = find(txt == '"');
stops = find(txt == ',' | txt == newline | txt == char(13));
qi = 1;
si = 1;
fields = {};
i = 1;
while i <= n
    if txt(i) == '"'
        parts = {};
        j = i + 1;
        while true
            while qi <= numel(q) && q(qi) < j
                qi = qi + 1;
            end
            if qi > numel(q)
                error('TestWorkbookReplacementSplitsMetadata:csv', ...
                    'Unterminated quoted field starting at character %d.', i);
            end
            k = q(qi);
            if k < n && txt(k + 1) == '"'
                parts{end + 1} = txt(j:k); %#ok<AGROW> keeps one quote of the pair
                j = k + 2;
            else
                parts{end + 1} = txt(j:k - 1); %#ok<AGROW>
                i = k + 1;
                break;
            end
        end
        fields{end + 1} = [parts{:}]; %#ok<AGROW>
    else
        while si <= numel(stops) && stops(si) < i
            si = si + 1;
        end
        if si > numel(stops)
            k = n + 1;
        else
            k = stops(si);
        end
        fields{end + 1} = txt(i:k - 1); %#ok<AGROW>
        i = k;
    end
    %what ends the field
    if i > n
        break;
    elseif txt(i) == ','
        i = i + 1;
        if i > n
            fields{end + 1} = ''; %#ok<AGROW> a comma at the very end leaves an empty last field
        end
    elseif txt(i) == newline
        recs{end + 1} = fields; %#ok<AGROW>
        fields = {};
        i = i + 1;
    elseif txt(i) == char(13) && i < n && txt(i + 1) == newline
        recs{end + 1} = fields; %#ok<AGROW>
        fields = {};
        i = i + 2;
    else
        error('TestWorkbookReplacementSplitsMetadata:csv', ...
            'Character %d at position %d follows a field; only a comma or a line break may.', ...
            double(txt(i)), i);
    end
end
if ~isempty(fields)
    recs{end + 1} = fields;
end
end

function recs = localRecordsOf(C)
%The rows of a readcell result as records of field texts, in the form
%localCsvRecords returns: text as it is, an empty cell as '', any other value
%as a non-empty placeholder (so a value beside a header line is seen).
recs = cell(1, size(C, 1));
for r = 1:size(C, 1)
    row = cell(1, size(C, 2));
    for c = 1:size(C, 2)
        v = C{r, c};
        if ischar(v)
            row{c} = v;
        elseif isstring(v) && isscalar(v) && ~ismissing(v)
            row{c} = char(v);
        elseif isa(v, 'missing') || (isnumeric(v) && isempty(v))
            row{c} = '';
        elseif isnumeric(v) && isscalar(v)
            row{c} = sprintf('%.17g', v);
        else
            row{c} = ['<' class(v) '>'];
        end
    end
    recs{r} = row;
end
end

function s = localJoin(c)
%One line per entry; '' for an empty list.
if isempty(c)
    s = '';
else
    s = strjoin(reshape(c, 1, []), newline);
end
end
