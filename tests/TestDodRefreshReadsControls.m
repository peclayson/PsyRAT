classdef TestDodRefreshReadsControls < matlab.unittest.TestCase
    %TESTDODREFRESHREADSCONTROLS The DoD Reliability Setup screen's Refresh
    %button must rebuild the screen from the CURRENT selections, with the
    %four trial counts reset to the selected group's defaults.
    %
    % Finding (2026-09-04, Session C; ruled 2026-09-05): psyrat_dod_refresh
    % rebuilt the screen from the preferences captured when the button was
    % built and read none of the controls, so a group change followed by
    % Refresh came back on the original group, while chapter 9 promised that
    % Refresh is how the trial defaults follow a group change. The fix reuses
    % the reader the Generate button already uses (local_read_main_inputs)
    % and clears the stored trial counts so local_dod_initprefs recomputes
    % the new group's defaults.
    %
    % Headless, on the real screen, the way TestPlotLabelInterpreter drives
    % its Generate button. The expected defaults are derived from the
    % fixture's data table here (mean trials per participant for each mapped
    % event in the selected group, rounded), independently of the builder,
    % and cross-checked against a hand computation from the fixture's
    % baseCounts so the helper cannot be trivially wrong.

    methods (Test)
        function testRefreshFollowsAGroupChange(testCase)
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>
            pd = PsyRATTestDataFactory.makeDoDRelData();
            if ~isfield(pd, 'ver'); pd.ver = '0-test'; end
            if ~isfield(pd, 'proc') || ~isfield(pd.proc, 'measheader')
                pd.proc.measheader = 'meas';
            end
            groups = cellstr(string(pd.rel.groups(:)));
            testCase.assertEqual(numel(groups), 2, 'fixture precondition: two groups');

            psyrat_startview_dodiff('psyrat_prefs', psyrat_defaults, 'psyrat_data', pd);
            [fig, grp, map] = localScreen(testCase, pd);

            %control: the screen opens on group 1 with group 1's defaults
            testCase.assertEqual(grp.Value, 1, 'control: the screen opens on group 1');
            exp1 = localExpectedDefaults(pd, groups{1}, map);
            testCase.assertEqual(exp1, [10 13 9 12], ...
                'helper cross-check: group 1 hand computation from baseCounts(1,:)');
            testCase.assertEqual(localTrialFields(testCase, fig, pd), exp1, ...
                'control: the fresh screen shows group 1''s trial defaults');

            %the user's action: pick group 2, press Refresh
            grp.Value = 2;
            btn = localButton(testCase, fig, 'Refresh');
            cb = get(btn, 'Callback');
            feval(cb{1}, btn, [], cb{2:end});

            [fig2, grp2, map2] = localScreen(testCase, pd);
            testCase.verifyEqual(grp2.Value, 2, ...
                'Refresh must rebuild on the group selected in the popup');
            testCase.verifyEqual(map2, map, 'the ERP mapping must survive a refresh');
            exp2 = localExpectedDefaults(pd, groups{2}, map);
            testCase.assertEqual(exp2, [8 8 10 7], ...
                'helper cross-check: group 2 hand computation from baseCounts(2,:)');
            testCase.verifyEqual(localTrialFields(testCase, fig2, pd), exp2, ...
                'the four trial counts must follow the group change (the new group''s defaults)');
        end
    end
end

function [fig, grp, map] = localScreen(testCase, pd)
%The one open setup screen, its group popup, and the ERP mapping its four
%event popups show (top to bottom).
figs = findall(groot, 'Type', 'figure', 'Name', 'DoD Reliability Setup');
testCase.assertEqual(numel(figs), 1, 'exactly one DoD Reliability Setup screen must be open');
fig = figs(1);
groups = cellstr(string(pd.rel.groups(:)));
events = cellstr(string(pd.rel.events(:)));
pops = findall(fig, 'Style', 'popupmenu');
grp = gobjects(0); erp = gobjects(1, 0);   % graphics arrays: a [] would coerce handles to doubles
for k = 1:numel(pops)
    s = cellstr(pops(k).String);
    if isequal(s(:), groups(:)); grp = pops(k); end
    if isequal(s(:), events(:)); erp(end+1) = pops(k); end %#ok<AGROW>
end
testCase.assertNotEmpty(grp, 'group popup not found');
testCase.assertEqual(numel(erp), 4, 'four ERP popups expected');
[~, order] = sort(arrayfun(@(h) h.Position(2), erp), 'descend');   % ERP_1 is the top row
erp = erp(order);
map = arrayfun(@(h) h.Value, erp);
end

function counts = localTrialFields(testCase, fig, pd)
%The four trial-count edit fields, read by the row they share with an ERP
%popup (ERP_1 at the top), as numbers.
events = cellstr(string(pd.rel.events(:)));
pops = findall(fig, 'Style', 'popupmenu');
erp = pops(arrayfun(@(h) isequal(cellstr(h.String), events(:)), pops));
[~, order] = sort(arrayfun(@(h) h.Position(2), erp), 'descend');
erp = erp(order);
edits = findall(fig, 'Style', 'edit');
counts = zeros(1, 4);
for k = 1:4
    row = erp(k).Position(2);
    onrow = edits(arrayfun(@(h) abs(h.Position(2) - row) < 1e-6, edits));
    testCase.assertEqual(numel(onrow), 1, sprintf('ERP row %d must carry exactly one edit field', k));
    counts(k) = str2double(onrow.String);
end
end

function d = localExpectedDefaults(pd, groupname, map)
%Group defaults as the screen defines them: for each mapped event, the mean
%number of trials per participant in that group, rounded, floored at zero.
tbl = pd.rel.data;
events = cellstr(string(pd.rel.events(:)));
d = zeros(1, 4);
for k = 1:4
    sub = tbl(string(tbl.group) == string(groupname) & string(tbl.event) == string(events{map(k)}), :);
    [~, ~, ic] = unique(sub.id);
    d(k) = max(round(mean(accumarray(ic, 1))), 0);
end
end

function btn = localButton(testCase, fig, label)
btns = findall(fig, 'Style', 'pushbutton');
btn = [];
for k = 1:numel(btns)
    s = btns(k).String; if iscell(s); s = strjoin(s, ' '); end
    if strcmp(strtrim(s), label); btn = btns(k); end
end
testCase.assertNotEmpty(btn, sprintf('button "%s" not found', label));
end
