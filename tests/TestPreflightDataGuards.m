classdef TestPreflightDataGuards < PsyRATTestBase
    %Audit 2026-09-05: data-integrity guards in psyrat_preflight_validate.
    %Each pins one condition the loader and engine accepted silently and that
    %changed the estimand or the exported numbers: +/-Inf scores (ismissing
    %flags NaN only), labels carrying a CSV delimiter or quote (the CSV data
    %rows are written unquoted), a participant under two groups (persons are
    %numbered within group, so one person became two), and an occasion column
    %with one level (its presence routes the run to the two-facet design).
    %Two warnings: labels differing only in case or surrounding spaces, and
    %exact duplicate rows. Every failing fixture has a clean twin that
    %asserts the guard stays silent, so no assertion can pass vacuously.
    %
    %RED against the pre-fix code: every "raises" test (the ids do not exist).
    %The clean twins pass on both builds.
    %
    %The Statistics-Toolbox preflight for the gamma concurrent designs cannot
    %be driven RED on a machine that has the toolbox; its positive control
    %(no such error with the toolbox present) is pinned here and the negative
    %arm is documented as unpinnable in this environment.

    methods (Test)

        function testInfScoreIsRejected(testCase)
            tbl = localCleanTable();
            tbl.meas(3) = Inf;
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(localHasId(report.errors, 'meas:notfinite'), ...
                'an Inf score must be rejected');
            testCase.verifyTrue(contains(localMessage(report.errors, 'meas:notfinite'), 'row 3'), ...
                'the message must name the first offending row');
            tbl.meas(3) = -Inf;
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(localHasId(report.errors, 'meas:notfinite'), '-Inf too');
        end

        function testReservedCharactersInLabelsAreRejected(testCase)
            tbl = localCleanTable();
            tbl.event{1} = 'Go, correct';
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(localHasId(report.errors, 'labels:reservedcharacters'));
            testCase.verifyTrue(contains(localMessage(report.errors, 'labels:reservedcharacters'), '''event'''));

            tbl = localCleanTable();
            tbl.id{2} = 's1"';
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(localHasId(report.errors, 'labels:reservedcharacters'), 'a quote too');

            tbl = localCleanTable();
            tbl.group{4} = sprintf('B%sx', newline);
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(localHasId(report.errors, 'labels:reservedcharacters'), 'a line break too');
        end

        function testParticipantUnderTwoGroupsWarns(testCase)
            tbl = localCleanTable();
            %s1 has rows in both groups: a WARNING, not a refusal (owner ruling
            %2026-09-07: ids numbered within group are a legitimate convention,
            %the core keys persons by (group, id), and the toolbox's own dev
            %fixture labels a within-person task as group), so the run proceeds
            %and the message names the count and an example id
            tbl.group{1} = 'B';
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(report.ok);
            testCase.verifyFalse(localHasId(report.errors, 'preflight:groupWithinParticipant'));
            testCase.verifyTrue(localHasId(report.warnings, 'preflight:groupWithinParticipant'));
            testCase.verifyTrue(contains(localMessage(report.warnings, 'preflight:groupWithinParticipant'), '''s1'''));
        end

        function testSingleOccasionLevelIsRejected(testCase)
            tbl = localCleanTable();
            tbl.time = repmat({'t1'}, height(tbl), 1);
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(localHasId(report.errors, 'preflight:singleOccasionLevel'));
            %two levels pass this guard
            tbl.time(1:2:end) = {'t2'};
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyFalse(localHasId(report.errors, 'preflight:singleOccasionLevel'));
        end

        function testCaseAndWhitespaceVariantsWarn(testCase)
            tbl = localCleanTable();
            tbl.event{1} = 'a';       % 'A' elsewhere
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(localHasId(report.warnings, 'preflight:labelCaseVariants'));
            testCase.verifyTrue(report.ok, 'a variant is a warning, not an error');

            tbl = localCleanTable();
            tbl.event{1} = 'A ';
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(localHasId(report.warnings, 'preflight:labelCaseVariants'), 'trailing space too');
        end

        function testDuplicateRowsWarn(testCase)
            tbl = localCleanTable();
            tbl = [tbl; tbl(1, :)];
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(localHasId(report.warnings, 'preflight:duplicateRows'));
            testCase.verifyTrue(contains(localMessage(report.warnings, 'preflight:duplicateRows'), '1 row(s)'));
            testCase.verifyTrue(report.ok, 'duplicates are a warning, not an error');
        end

        function testCleanTableRaisesNoneOfTheNewGuards(testCase)
            %positive control for every guard above
            tbl = localCleanTable();
            report = psyrat_preflight_validate('datatable', tbl);
            testCase.verifyTrue(report.ok);
            ids = {'meas:notfinite', 'labels:reservedcharacters', ...
                'preflight:singleOccasionLevel', 'preflight:statisticsToolboxRequired'};
            for k = 1:numel(ids)
                testCase.verifyFalse(localHasId(report.errors, ids{k}), ids{k});
            end
            wids = {'preflight:groupWithinParticipant', 'preflight:labelCaseVariants', ...
                'preflight:duplicateRows'};
            for k = 1:numel(wids)
                testCase.verifyFalse(localHasId(report.warnings, wids{k}), wids{k});
            end
        end

        function testToolboxPreflightStaysSilentWithTheToolbox(testCase)
            %Positive control only: with the Statistics and Machine Learning
            %Toolbox present the gamma concurrent request must not be refused
            %on toolbox grounds. The negative arm needs a toolbox-less MATLAB.
            testCase.assumeTrue(license('test', 'Statistics_Toolbox') == 1 && ...
                exist('gaminv', 'file') > 0, ...
                'positive control needs the Statistics and Machine Learning Toolbox');
            tbl = localCleanTable();
            tbl.meas = abs(tbl.meas) + 1;     % strictly positive for the gamma family
            report = psyrat_preflight_validate('datatable', tbl, 'family', 'gamma', ...
                'diffest', 2, 'diffrescor', 2);
            testCase.verifyFalse(localHasId(report.errors, 'preflight:statisticsToolboxRequired'));
        end

    end
end

function tbl = localCleanTable()
%two participants x two groups... no: two groups, each with one participant,
%two events, three trials per cell, all labels distinct in every spelling
ids = {}; groups = {}; events = {}; meas = [];
for g = 1:2
    gname = char('A' + g - 1);
    id = sprintf('s%d', g);
    for e = 1:2
        ename = char('A' + e - 1);
        for t = 1:3
            ids{end+1, 1} = id; %#ok<AGROW>
            groups{end+1, 1} = gname; %#ok<AGROW>
            events{end+1, 1} = ename; %#ok<AGROW>
            meas(end+1, 1) = 0.1 * numel(meas) + 0.01 * t; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, groups, events, 'VariableNames', {'id', 'meas', 'group', 'event'});
end

function tf = localHasId(items, id)
tf = ~isempty(items) && any(strcmp({items.id}, id));
end

function msg = localMessage(items, id)
%'' when the id is absent, so a contains() check on it fails the test cleanly
%instead of erroring out of the test method
msg = '';
if isempty(items)
    return;
end
k = find(strcmp({items.id}, id), 1, 'first');
if ~isempty(k)
    msg = items(k).message;
end
end
