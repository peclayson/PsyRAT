classdef TestCsvHeaderLineBreaks < PsyRATTestBase
    %A header line in a CSV export must reach a spreadsheet as ONE record that
    %holds ONE field, even when its text carries a line break (the B19 quoting
    %rule extended to carriage returns and line feeds; support for G109/G110).
    %
    %THE DEFECT. Two writers decide whether to quote a header line with the
    %same test: the shared helper psyrat_quote_csv_header (the table writers,
    %here psyrat_depoverallt's Save Table) and the inline copy in the dynamic-
    %reliability viewer's CSV arm (psyrat_startview_dynrel, exportDynrelTable).
    %Both quote only on a comma or a double quote. A header line that carries
    %a CR or an LF and no comma or quote is written bare, so a reader ends the
    %record at the break: the rest of the line becomes a record of its own and
    %every row below it moves down. psyrat_csv_field (G109) already quotes CR
    %and LF for single fields; these whole-line header writers do not.
    %
    %ROUTE. The dataset name (rel.filename) reaches both writers as the header
    %line "Dataset: <name>". Each test builds a fit-free result from
    %PsyRATTestDataFactory, sets rel.filename, opens the real output window and
    %fires its real Save button. The save dialog is answered by a uiputfile
    %shadow, as in TestExportCallbacks. The XLSX branches are out of scope.
    %
    %ORACLE. The file is read by the RFC 4180 reader below, written for this
    %class and independent of the production quoting code. Inside a quoted
    %field, commas, doubled quotes and line breaks are text. Outside quotes,
    %CR, LF or CRLF ends a record, as spreadsheet readers and Python's csv
    %module treat them. The expected Dataset field is literal text built from
    %the case name. Every record after it (the rest of the header, the table's
    %column-name row and its body) is compared with the same export from a
    %plain-named control run, which shows that nothing below the line moved
    %without a second copy of the table's expected content.
    %
    %CASES. Five texts per route, 15 tests:
    %  plain          - no comma, quote or break; must stay unquoted and
    %                   byte-identical (control)
    %  commaQuote     - a comma and a double quote; the existing B19 rule
    %                   quotes it and doubles the quote (control)
    %  lineFeed, carriageReturn, crlf - a break with no comma or quote, which
    %                   exposes the predicate
    %Before the header repair the six line-break tests of the two Dataset
    %routes failed and their four controls passed. The line-break tests fail
    %by verification, not by an aborted assertion: their assertions are
    %fixture welds only.
    %
    %OBSERVED-DoD ROUTE (third method). psyrat_dod_observedt writes the
    %contrast label, built from the four event names, into the header lines
    %(through psyrat_quote_csv_header) AND into the Label cell of its one
    %data row, which its own local_csv_scalar quotes on a comma, a quote or
    %an LF but not on a bare CR. Here one event of the fit-free DoD fixture is
    %renamed to the case text in the three places the name lives, so trial
    %counts and draws are unchanged; the Label cell must equal the literal
    %contrast and every other cell must match the unrenamed control. Expected
    %before that writer's repair: carriageReturn fails, the other four pass.

    % Copyright (C) 2016-2026 Peter E. Clayson
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program (gpl.txt). If not, see
    %     <http://www.gnu.org/licenses/>.

    properties
        ShadowDir
    end

    properties (TestParameter)
        %name: the text under test, written into rel.filename by the Dataset
        %methods and used as a renamed event by the observed-DoD method.
        %kind: which extra checks apply (see localVerifyDatasetLine).
        headerCase = struct( ...
            'plain', struct('name', 'csvhdr_plain.psyrat', 'kind', 'plain'), ...
            'commaQuote', struct('name', 'csvhdr_comma,quote"x.psyrat', 'kind', 'commaQuote'), ...
            'lineFeed', struct('name', ['csvhdr_lf_head' char(10) 'csvhdr_lf_tail.psyrat'], 'kind', 'break'), ...
            'carriageReturn', struct('name', ['csvhdr_cr_head' char(13) 'csvhdr_cr_tail.psyrat'], 'kind', 'break'), ...
            'crlf', struct('name', ['csvhdr_crlf_head' char([13 10]) 'csvhdr_crlf_tail.psyrat'], 'kind', 'break'));
    end

    methods (TestMethodSetup)
        function installSaveShadow(testCase)
            %answer every save dialog from a global, so the real Save callbacks
            %run without a modal file picker
            testCase.ShadowDir = tempname;
            mkdir(testCase.ShadowDir);
            fid = fopen(fullfile(testCase.ShadowDir, 'uiputfile.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function [filename,pathname,filterindex] = uiputfile(varargin)', ...
                '% TEST SHADOW: answers the save dialog from a global.', ...
                'global PSYRAT_CSVHDR_SAVE', ...
                '[p,n,e] = fileparts(PSYRAT_CSVHDR_SAVE);', ...
                'filename = [n e]; pathname = [p filesep]; filterindex = 1;', ...
                'end');
            fclose(fid);
            addpath(testCase.ShadowDir, '-begin');
            rehash;
            testCase.addTeardown(@() localRemoveShadow(testCase.ShadowDir));
            testCase.assertTrue(startsWith(which('uiputfile'), testCase.ShadowDir), ...
                'the save-dialog shadow is not first on the path');
        end
    end

    methods (Test)
        function testLegacyTableCsvKeepsDatasetLineOneRecord(testCase, headerCase)
            %psyrat_depoverallt, Save Table, CSV: the Dataset line is written
            %through psyrat_quote_csv_header
            control = localSaveLegacyTable(testCase, 'csvhdr_control.psyrat');
            target = localSaveLegacyTable(testCase, headerCase.name);
            localVerifyDatasetLine(testCase, target, control, headerCase, ...
                'depoverallt Save Table (CSV)');
        end

        function testDynrelVarianceCsvKeepsDatasetLineOneRecord(testCase, headerCase)
            %dynamic-reliability viewer, Save variance table, CSV: the Dataset
            %line is written through the viewer's inline quoting test
            control = localSaveDynrelVariance(testCase, 'csvhdr_control.psyrat');
            target = localSaveDynrelVariance(testCase, headerCase.name);
            localVerifyDatasetLine(testCase, target, control, headerCase, ...
                'dynrel Save variance table (CSV)');
        end

        function testObservedDodCsvKeepsEventLabelOneField(testCase, headerCase)
            %psyrat_dod_observedt, Save Table, CSV: event ERP_C renamed to the
            %case text; the Label cell of the table row is written through the
            %writer's own local_csv_scalar
            control = localSaveObservedDod(testCase, '');
            target = localSaveObservedDod(testCase, headerCase.name);
            localVerifyObservedDodRow(testCase, target, control, headerCase);
        end
    end
end


function raw = localSaveLegacyTable(testCase, datasetName)
%open psyrat_depoverallt's output window for a fit-free one-facet result and
%save its table as CSV through the real Save Table callback
pd = localWithHeaderFields(PsyRATTestDataFactory.makeICSummaryDataWithTrial(), datasetName);
close(findall(0, 'Type', 'figure', 'Tag', 'psyrat_output'));
psyrat_depoverallt('psyrat_data', pd, 'gui', 1);
fig = findobj(0, 'Type', 'figure', 'Tag', 'psyrat_output');
testCase.assertNotEmpty(fig, 'psyrat_depoverallt did not open its output window');
raw = localFireSaveAndRead(testCase, fig(1), 'Save Table');
close(fig);
end


function raw = localSaveDynrelVariance(testCase, datasetName)
%open the dynamic-reliability viewer on a fit-free group-level difference
%result (one Save button: the variance table) and save that table as CSV
pd = localWithHeaderFields(PsyRATTestDataFactory.makeDynrelDiffSummaryData(), datasetName);
name = 'Dynamic Reliability: View Results';
close(findall(0, 'Type', 'figure', 'Name', name));
psyrat_startview_dynrel('psyrat_prefs', psyrat_defaults(), 'psyrat_data', pd);
fig = findall(0, 'Type', 'figure', 'Name', name);
testCase.assertNotEmpty(fig, 'the dynamic-reliability viewer did not open');
raw = localFireSaveAndRead(testCase, fig(1), 'Save variance table');
close(fig);
end


function raw = localSaveObservedDod(testCase, eventName)
%open psyrat_dod_observedt's window for the fit-free DoD fixture (group 1,
%contrast (ERP_A - ERP_B) - (ERP_C - ERP_D)) and save its table as CSV through
%the real Save Table callback. A non-empty eventName renames ERP_C in the three
%places the name lives (rel.events, the data's event column, rel.dod_map), so
%the trial counts and draws are those of the unrenamed control ('').
pd = localWithHeaderFields(PsyRATTestDataFactory.makeDoDRelData(), 'csvhdr_dod.psyrat');
if ~isempty(eventName)
    old = 'ERP_C';
    pd.rel.events(strcmp(pd.rel.events, old)) = {eventName};
    pd.rel.data.event(strcmp(pd.rel.data.event, old)) = {eventName};
    pd.rel.dod_map(strcmp(pd.rel.dod_map, old)) = {eventName};
    testCase.assertTrue(any(strcmp(pd.rel.data.event, eventName)), ...
        'fixture weld: the renamed event has no data rows');
end
[~, map] = ismember(pd.rel.dod_map, pd.rel.events);
testCase.assertEqual(map(:)', [1 2 3 4], ...
    'fixture weld: the DoD map no longer resolves to events 1-4');
close(findall(0, 'Type', 'figure', 'Tag', 'psyrat_output'));
psyrat_dod_observedt('psyrat_data', pd, 'groupidx', 1, 'map', map(:)', 'gui', 1);
fig = findobj(0, 'Type', 'figure', 'Tag', 'psyrat_output');
testCase.assertNotEmpty(fig, 'psyrat_dod_observedt did not open its output window');
raw = localFireSaveAndRead(testCase, fig(1), 'Save Table');
close(fig);
end


function localVerifyObservedDodRow(testCase, raw, controlRaw, hc)
%(1) every header record is one field and the Contrast line holds the label as
%written; (2) the column-name row matches the control and the Label cell is
%the literal contrast built from the renamed event; (3) every cell after Label
%matches the unrenamed control and no record was added anywhere
what = 'dod_observedt Save Table (CSV)';
expected = ['(ERP_A - ERP_B) - (' hc.name ' - ERP_D)'];

%weld: a line-break name must reach the writer unaltered
if strcmp(hc.kind, 'break')
    testCase.assertTrue(contains(raw, hc.name), sprintf( ...
        '%s: the event name with its line break did not reach the file', what));
end

[recs, open] = localParseCsv(raw);
[crecs, copen] = localParseCsv(controlRaw);
testCase.assertFalse(copen, sprintf('%s: the control export ends inside a quoted field', what));
testCase.verifyFalse(open, sprintf('%s: the export ends inside a quoted field', what));
kc = localFindColumnRow(crecs);
k = localFindColumnRow(recs);
testCase.assertNotEmpty(kc, sprintf('%s: the control export has no column-name row', what));
testCase.assertNotEmpty(k, sprintf('%s: the export has no column-name row', what));
testCase.assertTrue(kc < numel(crecs) && k < numel(recs), sprintf( ...
    '%s: no data row below the column-name row', what));
crow = crecs{kc+1};
testCase.assertEqual(crow{1}, '(ERP_A - ERP_B) - (ERP_C - ERP_D)', ...
    'fixture weld: the control contrast is not the unrenamed fixture''s');
testCase.assertEqual(numel(crow), numel(crecs{kc}), ...
    'fixture weld: the control data row and column names differ in width');

%(1) the header: same length as the control, one field per record, and the
%Contrast line holding the label as written
testCase.verifyEqual(k, kc, sprintf( ...
    '%s: the header has %d records, the control %d', what, k - 1, kc - 1));
for i = 1:k-1
    testCase.verifyEqual(numel(recs{i}), 1, sprintf( ...
        '%s: header record %d splits into %d fields', what, i, numel(recs{i})));
end
j = find(cellfun(@(r) startsWith(r{1}, 'Contrast: '), recs(1:k-1)), 1);
testCase.assertNotEmpty(j, sprintf('%s: no Contrast header line', what));
testCase.verifyEqual(recs{j}{1}, ['Contrast: ' expected], sprintf( ...
    '%s: the Contrast header line does not hold the label as written', what));

%(2) the column names and the Label cell
testCase.verifyEqual(recs{k}, crecs{kc}, sprintf( ...
    '%s: the column-name row differs from the control', what));
row = recs{k+1};
testCase.verifyEqual(row{1}, expected, sprintf( ...
    '%s: the Label cell does not hold the contrast as written', what));

%(3) the counts, coefficients and mean trials, and the record count
testCase.verifyEqual(row(2:end), crow(2:end), sprintf( ...
    '%s: the cells after Label differ from the unrenamed control', what));
testCase.verifyEqual(numel(recs), numel(crecs), sprintf( ...
    '%s: %d records against the control''s %d', what, numel(recs), numel(crecs)));
end


function k = localFindColumnRow(recs)
%index of the table's column-name row: the first record whose first field is
%exactly 'Label'
k = find(cellfun(@(r) strcmp(r{1}, 'Label'), recs), 1);
end


function pd = localWithHeaderFields(pd, datasetName)
%the Save headers read these fields, which the summary-only factory data leaves
%out; the values match TestExportCallbacks' fixtures so the header is fixed
pd.ver = 'test';
pd.rel.filename = datasetName;
pd.rel.nchains = 2;
pd.rel.niter = 600;
pd.rel.seed = 12345;
pd.rel.engine = 'cmdstan';
pd.rel.priors = psyrat_default_priors;
pd.rel.out.conv.converged = 1;
pd.rel.out.conv.ndivergent = 0;
end


function raw = localFireSaveAndRead(testCase, fig, labelfrag)
%fire the Save button whose label contains labelfrag with the dialog answered
%by a fresh .csv path, then return the file's bytes as written
global PSYRAT_CSVHDR_SAVE %#ok<GVMIS>
outfile = [tempname '.csv'];
testCase.addTeardown(@() localDelete(outfile));
PSYRAT_CSVHDR_SAVE = outfile;
btn = [];
for b = findall(fig, 'Style', 'pushbutton')'
    s = b.String;
    if iscell(s); s = strjoin(s, ' '); end
    if contains(string(s), labelfrag)
        btn = b;
        break;
    end
end
testCase.assertNotEmpty(btn, sprintf('no "%s" button on the window', labelfrag));
cb = btn.Callback;
if iscell(cb)
    %table-window buttons store {fcn, extra args...}
    feval(cb{1}, btn, [], cb{2:end});
else
    %viewer buttons store a nested function handle
    cb(btn, []);
end
testCase.assertTrue(isfile(outfile), sprintf('the "%s" callback wrote no file', labelfrag));
raw = fileread(outfile);
end


function localVerifyDatasetLine(testCase, raw, controlRaw, hc, what)
%(1) the Dataset line is one record holding one field equal to the name as
%written; (2) every record after it matches the plain-named control; (3) the
%existing byte-level behavior for plain and comma/quote lines
LF = char(10);
expected = ['Dataset: ' hc.name];

%weld: a line-break name must reach the writer unaltered (quoted or not, the
%file then contains it verbatim), or the test would prove nothing
if strcmp(hc.kind, 'break')
    testCase.assertTrue(contains(raw, hc.name), sprintf( ...
        '%s: the dataset name with its line break did not reach the file', what));
end

[recs, open] = localParseCsv(raw);
[crecs, copen] = localParseCsv(controlRaw);
testCase.assertFalse(copen, sprintf('%s: the control export ends inside a quoted field', what));
testCase.verifyFalse(open, sprintf('%s: the export ends inside a quoted field', what));
kc = localFindDatasetRecord(crecs);
k = localFindDatasetRecord(recs);
testCase.assertNotEmpty(kc, sprintf('%s: the control export has no Dataset record', what));
testCase.assertNotEmpty(k, sprintf('%s: the export has no Dataset record', what));

%(1) one record, one field, the name as written
testCase.verifyEqual(numel(recs{k}), 1, sprintf( ...
    '%s: the Dataset line splits into %d fields', what, numel(recs{k})));
testCase.verifyEqual(recs{k}{1}, expected, sprintf( ...
    '%s: the Dataset field does not hold the dataset name as written', what));

%(2) nothing below it moved. The weld makes sure the comparison reaches the
%table (a multi-field row) and is not a comparison of empty tails.
tail = recs(k+1:end);
ctail = crecs(kc+1:end);
testCase.assertTrue(any(cellfun(@numel, ctail) > 1), sprintf( ...
    '%s: fixture weld: the control has no table row below the Dataset line', what));
testCase.verifyEqual(tail, ctail, sprintf( ...
    ['%s: the records below the Dataset line differ from the plain-named ' ...
    'control (%d records vs %d)'], what, numel(tail), numel(ctail)));

%(3) the rule that already holds must keep holding
switch hc.kind
    case 'plain'
        testCase.verifyTrue(contains(raw, [LF expected LF]), sprintf( ...
            '%s: a plain Dataset line must be written unquoted and byte-identical', what));
    case 'commaQuote'
        testCase.verifyTrue(contains(raw, ...
            [LF '"Dataset: csvhdr_comma,quote""x.psyrat"' LF]), sprintf( ...
            '%s: a comma/quote line must be quoted with the quote doubled', what));
end
end


function k = localFindDatasetRecord(recs)
%index of the first record whose first field starts with 'Dataset: '
k = find(cellfun(@(r) ~isempty(r) && startsWith(r{1}, 'Dataset: '), recs), 1);
end


function [recs, unterminated] = localParseCsv(txt)
%RFC 4180 reader: one cell row of char fields per record. A double quote opens
%a quoted field only at the start of a field; inside it, a doubled quote is one
%quote and commas, CR and LF are text. Outside quotes a comma ends a field and
%CR, LF or CRLF ends a record. unterminated is true when the text ends inside a
%quoted field.
CR = char(13);
LF = char(10);
Q = '"';
recs = {};
rec = {};
fld = '';
inq = false;
n = numel(txt);
i = 1;
while i <= n
    c = txt(i);
    if inq
        if c == Q
            if i < n && txt(i+1) == Q
                fld(end+1) = Q; %#ok<AGROW> doubled quote: one literal quote
                i = i + 1;
            else
                inq = false;    %closing quote
            end
        else
            fld(end+1) = c; %#ok<AGROW>
        end
    elseif c == Q && isempty(fld)
        inq = true;
    elseif c == ','
        rec{end+1} = fld; %#ok<AGROW>
        fld = '';
    elseif c == CR || c == LF
        if c == CR && i < n && txt(i+1) == LF
            i = i + 1;          %CRLF is one record terminator
        end
        rec{end+1} = fld; %#ok<AGROW>
        recs{end+1} = rec; %#ok<AGROW>
        rec = {};
        fld = '';
    else
        fld(end+1) = c; %#ok<AGROW>
    end
    i = i + 1;
end
if ~isempty(fld) || ~isempty(rec)
    rec{end+1} = fld;
    recs{end+1} = rec;
end
unterminated = inq;
end


function localDelete(f)
if isfile(f)
    delete(f);
end
end


function localRemoveShadow(d)
%the base-class teardown may already have restored the path, so rmpath only
%when the shadow folder is still on it
clear global PSYRAT_CSVHDR_SAVE
if any(strcmp(strsplit(path, pathsep), d))
    rmpath(d);
end
if isfolder(d)
    rmdir(d, 's');
end
end
