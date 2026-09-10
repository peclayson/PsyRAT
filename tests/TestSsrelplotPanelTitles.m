classdef TestSsrelplotPanelTitles < matlab.unittest.TestCase
    %TESTSSRELPLOTPANELTITLES psyrat_ssrelplot must title every panel on
    %events-only and groups-only subject-level runs.
    %
    % Finding (2026-09-03, manual capture; fixed in the pre-beta batch): the
    % loader spells an absent facet's name '' (psyrat_ssrelplot.m sets gnames
    % or enames to {''} when rel.groups or rel.events is 'none'), but the two
    % middle title branches tested strcmpi(name,'none') with the negation on
    % the wrong side of the name they print, so neither could fire and every
    % panel fell through to title(''). A reader of an events-only run could
    % not tell which panel was which event. The groups x events case printed
    % 'GroupA: inc_cor' correctly and is kept here as the positive control.
    %
    % The fixture is the factory's single (group, event) summary replicated
    % across the requested strata, exactly as TestPlotLabelInterpreter does;
    % only the names matter. An absent facet is spelled 'none', the way
    % psyrat_loadfile spells it, so the fixture reaches the same branch the
    % loader's output reaches.

    properties
        FigsBefore
    end

    methods (TestMethodSetup)
        function recordOpenFigures(testCase)
            testCase.FigsBefore = findall(groot, 'Type', 'figure');
        end
    end

    methods (TestMethodTeardown)
        function closeNewFigures(testCase)
            figs = setdiff(findall(groot, 'Type', 'figure'), testCase.FigsBefore);
            if ~isempty(figs); delete(figs); end
        end
    end

    methods (Test)
        function testEventsOnlyDependabilityPanelsCarryEventNames(testCase)
            events = {'inc_cor'; 'inc_err'};
            pd = localExpandFacets(PsyRATTestDataFactory.makeICSummaryData(), ...
                'none', events);
            titles = localDrawAndCollectTitles(testCase, pd, 'dep');
            testCase.verifyEqual(sort(titles), sort(events), ...
                'events-only run: every panel must be titled with its event name');
        end

        function testGroupsOnlyDependabilityPanelsCarryGroupNames(testCase)
            groups = {'GroupA'; 'GroupB'};
            pd = localExpandFacets(PsyRATTestDataFactory.makeICSummaryData(), ...
                groups, 'none');
            titles = localDrawAndCollectTitles(testCase, pd, 'dep');
            testCase.verifyEqual(sort(titles), sort(groups), ...
                'groups-only run: every panel must be titled with its group name');
        end

        function testEventsOnlyIccPanelsCarryEventNames(testCase)
            %The ICC branch carries its own copy of the title sites.
            events = {'inc_cor'; 'inc_err'};
            pd = localExpandFacets(PsyRATTestDataFactory.makeICSummaryData(), ...
                'none', events);
            titles = localDrawAndCollectTitles(testCase, pd, 'icc');
            testCase.verifyEqual(sort(titles), sort(events), ...
                'events-only ICC run: every panel must be titled with its event name');
        end

        function testGroupsByEventsPanelsCarryBothNamesPositiveControl(testCase)
            groups = {'GroupA'; 'GroupB'};
            events = {'inc_cor'; 'inc_err'};
            pd = localExpandFacets(PsyRATTestDataFactory.makeICSummaryData(), ...
                groups, events);
            [gi, ei] = ndgrid(1:numel(groups), 1:numel(events));
            expected = arrayfun(@(g, e) sprintf('%s: %s', groups{g}, events{e}), ...
                gi(:), ei(:), 'UniformOutput', false);
            titles = localDrawAndCollectTitles(testCase, pd, 'dep');
            testCase.verifyEqual(sort(titles), sort(expected), ...
                'groups x events run (positive control): titles read "group: event"');
        end
    end
end

function titles = localDrawAndCollectTitles(testCase, pd, stat)
%Draw the caterpillar for PD and return the title of every axes it drew.
before = findall(groot, 'Type', 'figure');
psyrat_ssrelplot('psyrat_data', pd, 'stat', stat);
figs = setdiff(findall(groot, 'Type', 'figure'), before);
testCase.assertNotEmpty(figs, 'psyrat_ssrelplot drew no figure');
axs = findall(figs, 'Type', 'axes');
testCase.assertNotEmpty(axs, 'psyrat_ssrelplot drew no axes');
titles = cell(numel(axs), 1);
for k = 1:numel(axs)
    t = axs(k).Title.String;
    if iscell(t); t = strjoin(t, ' '); end
    titles{k} = char(t);
end
end

function pd = localExpandFacets(pd, groups, events)
%Replicate the fixture's single (group, event) cell across the requested
%strata. GROUPS and EVENTS are cellstr name lists, or the char 'none' for an
%absent facet (the loader's spelling); an absent facet keeps one unnamed
%stratum, which is what the loader produces.
if ischar(groups)
    gnames = {''};
    pd.rel.groups = groups;
else
    gnames = groups(:);
    pd.rel.groups = groups(:);
end
if ischar(events)
    enames = {''};
    pd.rel.events = events;
else
    enames = events(:);
    pd.rel.events = events(:);
end
g1 = pd.relsummary.group(1);
e1 = g1.event(1);
d1 = pd.relsummary.data.g(1).e(1);
for g = 1:numel(gnames)
    pd.relsummary.group(g) = g1;
    pd.relsummary.group(g).name = gnames{g};
    for e = 1:numel(enames)
        pd.relsummary.group(g).event(e) = e1;
        pd.relsummary.group(g).event(e).name = enames{e};
        pd.relsummary.data.g(g).e(e) = d1;
    end
end
end
