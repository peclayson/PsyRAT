classdef TestSelectDialogEmptySave < PsyRATTestBase
    %Sweep follow-up F1 (ruled 2026-09-16): the Save buttons of Specify Which
    %Groups / Events / Occasions to Process accepted an empty multi-select
    %listbox, stored cell(1,0), and psyrat_loadfile's ~isempty gate then
    %processed EVERY level while the Select button's tooltip read
    %"Selected: N/N" and the re-opened list showed nothing (observed live,
    %live GUI Item 11 A+). Select Measurements already refused the same state
    %with a modal errordlg; the three now do the same, with their own titles,
    %and leave the selection screen open so the user can choose or click Back.
    %
    %ROUTE. The Specify Inputs screen is entered through psyrat_startproc with
    %a uigetfile shadow (the drive route of TestStartprocAnalyzeDrive); with no
    %home figure open the entry loads psyrat_defaults and auto-detects the
    %fixture's group / event / time columns. Each Select button's own Callback
    %opens its dialog; the dialog's listbox Value is emptied (a Cmd-click does
    %that on a Min 0 / Max n list) and its Save button's own Callback fired.
    %
    %OBSERVABLES. The refusal dialog by title and its WindowStyle; the
    %selection dialog still open by Tag; no rebuilt Specify Inputs. The
    %positive control saves one level and reads the rebuilt screen's tooltip.

    properties (Access = private)
        ShadowDir
        TmpDir
    end

    methods (TestMethodSetup)
        function setupShadowAndFixture(testCase)
            testCase.ShadowDir = tempname;
            mkdir(testCase.ShadowDir);
            fid = fopen(fullfile(testCase.ShadowDir, 'uigetfile.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function [filename,pathname,filterindex] = uigetfile(varargin)', ...
                '% TEST SHADOW: answers the open dialog from a global.', ...
                'global PSYRAT_TEST_UIGETFILE', ...
                '[p,n,e] = fileparts(PSYRAT_TEST_UIGETFILE);', ...
                'filename = [n e]; pathname = [p filesep]; filterindex = 1;', ...
                'end');
            fclose(fid);
            addpath(testCase.ShadowDir, '-begin');
            rehash;
            testCase.assertTrue(startsWith(which('uigetfile'), testCase.ShadowDir), ...
                'the uigetfile shadow must be the resolved uigetfile');
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);
            csvfile = fullfile(testCase.TmpDir, 'select.csv');
            fid = fopen(csvfile, 'w');
            fprintf(fid, 'id,meas,group,event,time\n');
            row = 0;
            for g = 1:2
                for s = 1:4
                    for e = 1:2
                        for t = 1:2
                            for r = 1:2
                                row = row + 1;
                                fprintf(fid, 'g%ds%d,%0.3f,%s,e%d,t%d\n', g, s, ...
                                    mod(0.37 * row, 5) - 2.5, char('A' + g - 1), e, t);
                            end
                        end
                    end
                end
            end
            fclose(fid);
            global PSYRAT_TEST_UIGETFILE %#ok<GVMIS>
            PSYRAT_TEST_UIGETFILE = csvfile;
            close all force;
        end
    end

    methods (Test)

        function testEmptyGroupsSaveIsRefused(testCase)
            localRunFacet(testCase, 'Select Groups', 'psyrat_gui_whichgroups', 'No Groups Selected', 'No groups were selected.');
        end

        function testEmptyEventsSaveIsRefused(testCase)
            localRunFacet(testCase, 'Select Events', 'psyrat_gui_whichevents', 'No Events Selected', 'No events were selected.');
        end

        function testEmptyOccasionsSaveIsRefused(testCase)
            localRunFacet(testCase, 'Select Occasions', 'psyrat_gui_whichtimes', 'No Occasions Selected', 'No occasions were selected.');
        end

        function testOneLevelSaveStillCloses(testCase)
            %Positive control on every facet: a one-level selection is saved,
            %the dialog closes, the screen is rebuilt, and its Select button's
            %tooltip reports 1/2.
            facets = {'Select Groups', 'psyrat_gui_whichgroups'; ...
                'Select Events', 'psyrat_gui_whichevents'; ...
                'Select Occasions', 'psyrat_gui_whichtimes'};
            for k = 1:size(facets, 1)
                setup = localOpenSetup(testCase);
                dlg = localOpenDialog(testCase, setup, facets{k, 1}, facets{k, 2});
                lb = findall(dlg, 'Style', 'listbox');
                testCase.assertNumElements(lb, 1);
                testCase.assertEqual(numel(get(lb, 'String')), 2, 'the fixture has two levels');
                set(lb, 'Value', 1);
                localFireButton(testCase, dlg, 'Save');
                testCase.verifyEmpty(findall(groot, 'Type', 'figure', 'Tag', facets{k, 2}), ...
                    sprintf('%s: a one-level Save must close the dialog', facets{k, 1}));
                setup2 = findall(groot, 'Type', 'figure', 'Name', 'Specify Inputs');
                testCase.assertNumElements(setup2, 1, sprintf('%s: the screen must be rebuilt', facets{k, 1}));
                tt = localTooltip(testCase, setup2, facets{k, 1});
                testCase.verifyTrue(contains(tt, 'Selected: 1/2'), ...
                    sprintf('%s: the tooltip must report the one saved level, got "%s"', facets{k, 1}, tt));
                close all force;
            end
        end

    end

    methods (TestMethodTeardown)
        function teardownShadowAndFixture(testCase)
            close all force;
            clear global PSYRAT_TEST_UIGETFILE
            if any(strcmp(strsplit(path, pathsep), testCase.ShadowDir))
                rmpath(testCase.ShadowDir);
            end
            rehash;
            if isfolder(testCase.ShadowDir); rmdir(testCase.ShadowDir, 's'); end
            if isfolder(testCase.TmpDir); rmdir(testCase.TmpDir, 's'); end
        end
    end
end


function localRunFacet(testCase, button, dlgtag, title, firstline)
setup = localOpenSetup(testCase);
dlg = localOpenDialog(testCase, setup, button, dlgtag);
lb = findall(dlg, 'Style', 'listbox');
testCase.assertNumElements(lb, 1, 'the selection dialog holds one listbox');
testCase.assertGreaterThan(get(lb, 'Max') - get(lb, 'Min'), 1, 'the list is multi-select, so it can be emptied');
set(lb, 'Value', []);
testCase.assertTrue(isempty(get(lb, 'Value')), 'the listbox must be empty before Save');
localFireButton(testCase, dlg, 'Save');
d = findall(groot, 'Type', 'figure', 'Name', title);
testCase.verifyNumElements(d, 1, sprintf('%s: an empty Save must raise the "%s" dialog', button, title));
if ~isempty(d)
    testCase.verifyEqual(get(d(1), 'WindowStyle'), 'modal', sprintf('%s: the refusal must be modal', button));
    s = findall(d(1), '-property', 'String');
    txt = strjoin(cellfun(@(v) strjoin(cellstr(v), ' '), get(s, 'String'), 'UniformOutput', false), ' ');
    testCase.verifyTrue(contains(txt, firstline), sprintf('%s: the dialog must say "%s"', button, firstline));
    testCase.verifyTrue(contains(txt, 'How to fix:'), sprintf('%s: the dialog must carry a How-to-fix line', button));
    delete(d);
end
testCase.verifyNumElements(findall(groot, 'Type', 'figure', 'Tag', dlgtag), 1, ...
    sprintf('%s: the selection dialog must stay open after a refused Save', button));
testCase.verifyEmpty(findall(groot, 'Type', 'figure', 'Name', 'Specify Inputs'), ...
    sprintf('%s: a refused Save must not rebuild the setup screen', button));
end


function fig = localOpenSetup(testCase)
close all force;
psyrat_startproc();
fig = findall(groot, 'Type', 'figure', 'Name', 'Specify Inputs');
testCase.assertNumElements(fig, 1, 'expected exactly one Specify Inputs screen');
end


function dlg = localOpenDialog(testCase, setup, button, dlgtag)
localFireButton(testCase, setup, button);
dlg = findall(groot, 'Type', 'figure', 'Tag', dlgtag);
testCase.assertNumElements(dlg, 1, sprintf('%s must open its selection dialog', button));
end


function localFireButton(testCase, fig, label)
btn = findall(fig, 'Style', 'pushbutton', 'String', label);
testCase.assertNumElements(btn, 1, sprintf('button "%s" not found', label));
cb = btn.Callback;
if iscell(cb)
    cb{1}(btn, [], cb{2:end});
else
    cb(btn, []);
end
end


function tt = localTooltip(testCase, fig, label)
btn = findall(fig, 'Style', 'pushbutton', 'String', label);
testCase.assertNumElements(btn, 1, sprintf('button "%s" not found on the rebuilt screen', label));
tt = get(btn, 'Tooltip');
if iscell(tt); tt = strjoin(tt, ' '); end
end
