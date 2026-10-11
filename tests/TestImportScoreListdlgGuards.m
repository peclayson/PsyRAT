classdef TestImportScoreListdlgGuards < PsyRATTestBase
    %Sweep follow-up F2 (ruled 2026-09-16): the two listdlg pickers of the
    %import/scoring screen tested only the ok flag. MATLAB's listdlg copies
    %whatever the listbox Value is on OK (listdlg.m doOK), so an emptied
    %single-select list returns idx = [] with ok = 1, and the next line
    %(methods{idx} in onAddSpec; specs(idx).spec_id in localChooseSpec) threw
    %outside any try: a console error and an OK button that looked dead.
    %
    %ROUTE. The screen is built by psyrat_startimportscore (no figure needed
    %beforehand; it closes a home screen if one is open) and its state lives in
    %guidata, so a test can hand it a minimal project and fire the buttons'
    %own Callbacks. listdlg and inputdlg are shadowed from a temp folder that
    %answers from globals and counts its calls; the count is the weld that
    %proves the picker was consulted (a silent return with zero calls would be
    %a wrong route, not a guard). The spec form after a real method choice is
    %an inputdlg, shadowed to Cancel so the positive control returns.
    %
    %OBSERVABLES. An error escaping the callback is the defect; a 'Scoring
    %Error' modal from the fake project proves a chosen index was consumed
    %(the positive control); the ModalArg census in TestModalRejectionDialogs
    %stays at 14 because the guard adds no dialog.

    properties (Access = private)
        ShadowDir
        Fig
    end

    methods (TestMethodSetup)
        function setupShadowsAndScreen(testCase)
            testCase.ShadowDir = tempname;
            mkdir(testCase.ShadowDir);
            fid = fopen(fullfile(testCase.ShadowDir, 'listdlg.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function [selection, value] = listdlg(varargin)', ...
                '% TEST SHADOW: answers the picker from globals and counts calls.', ...
                'global PSYRAT_TEST_LISTDLG_IDX PSYRAT_TEST_LISTDLG_OK PSYRAT_TEST_LISTDLG_CALLS', ...
                'if isempty(PSYRAT_TEST_LISTDLG_CALLS), PSYRAT_TEST_LISTDLG_CALLS = 0; end', ...
                'PSYRAT_TEST_LISTDLG_CALLS = PSYRAT_TEST_LISTDLG_CALLS + 1;', ...
                'selection = PSYRAT_TEST_LISTDLG_IDX; value = PSYRAT_TEST_LISTDLG_OK;', ...
                'end');
            fclose(fid);
            fid = fopen(fullfile(testCase.ShadowDir, 'inputdlg.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function answer = inputdlg(varargin)', ...
                '% TEST SHADOW: the spec form, always cancelled; counts calls.', ...
                'global PSYRAT_TEST_INPUTDLG_CALLS', ...
                'if isempty(PSYRAT_TEST_INPUTDLG_CALLS), PSYRAT_TEST_INPUTDLG_CALLS = 0; end', ...
                'PSYRAT_TEST_INPUTDLG_CALLS = PSYRAT_TEST_INPUTDLG_CALLS + 1;', ...
                'answer = {};', ...
                'end');
            fclose(fid);
            addpath(testCase.ShadowDir, '-begin');
            rehash;
            testCase.assertTrue(startsWith(which('listdlg'), testCase.ShadowDir), ...
                'the listdlg shadow must be the resolved listdlg');
            testCase.assertTrue(startsWith(which('inputdlg'), testCase.ShadowDir), ...
                'the inputdlg shadow must be the resolved inputdlg');
            localResetGlobals();
            localCloseByTag('psyrat_import_gui');
            psyrat_startimportscore();
            testCase.Fig = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_import_gui');
            testCase.assertNumElements(testCase.Fig, 1, 'the import screen must be open');
        end
    end

    methods (Test)

        function testAddSpecEmptyIndexOnOkIsIgnored(testCase)
            %Site 1 (psyrat_startimportscore.m onAddSpec): OK with nothing
            %selected must return quietly; pre-fix methods{[]} threw.
            localSetProject(testCase, localProject(1));
            localAnswer([], 1);
            err = localFire(testCase, 'Add Scoring Spec');
            testCase.verifyEmpty(err, sprintf('Add Scoring Spec threw on an empty OK: %s', err));
            testCase.verifyEqual(localCalls('listdlg'), 1, 'the picker must have been consulted once');
            testCase.verifyEqual(localCalls('inputdlg'), 0, 'no form may open without a method');
        end

        function testAddSpecCancelIsIgnored(testCase)
            %Control: Cancel took the same silent route before the fix.
            localSetProject(testCase, localProject(1));
            localAnswer([], 0);
            err = localFire(testCase, 'Add Scoring Spec');
            testCase.verifyEmpty(err, sprintf('Add Scoring Spec threw on Cancel: %s', err));
            testCase.verifyEqual(localCalls('listdlg'), 1);
            testCase.verifyEqual(localCalls('inputdlg'), 0);
        end

        function testAddSpecSelectionReachesTheForm(testCase)
            %Positive control: a real choice reaches the spec form (shadowed to
            %Cancel), so the guard is not a blanket return.
            localSetProject(testCase, localProject(1));
            localAnswer(1, 1);
            err = localFire(testCase, 'Add Scoring Spec');
            testCase.verifyEmpty(err, sprintf('Add Scoring Spec threw on a real choice: %s', err));
            testCase.verifyEqual(localCalls('listdlg'), 1);
            testCase.verifyEqual(localCalls('inputdlg'), 1, 'the chosen method must open the form');
        end

        function testRunScoringEmptyIndexOnOkIsIgnored(testCase)
            %Site 2 (localChooseSpec via Run Scoring): OK with nothing selected
            %must return '' like Cancel; pre-fix specs([]).spec_id threw.
            localSetProject(testCase, localProject(2));
            localAnswer([], 1);
            err = localFire(testCase, 'Run Scoring');
            testCase.verifyEmpty(err, sprintf('Run Scoring threw on an empty OK: %s', err));
            testCase.verifyEqual(localCalls('listdlg'), 1, 'two specs must open the picker');
            testCase.verifyFalse(localDialogExists('Scoring Error'), ...
                'nothing may be scored when no spec was chosen');
        end

        function testRunScoringSelectionIsConsumed(testCase)
            %Positive control: a chosen index reaches the scoring call, which
            %fails on the fake project inside its try and reports modally.
            localSetProject(testCase, localProject(2));
            localAnswer(2, 1);
            err = localFire(testCase, 'Run Scoring');
            testCase.verifyEmpty(err, sprintf('Run Scoring threw outside its try: %s', err));
            testCase.verifyEqual(localCalls('listdlg'), 1);
            testCase.verifyTrue(localDialogExists('Scoring Error'), ...
                'the chosen spec must reach psyrat_project_score (which fails on the fake project)');
        end

        function testExportHandoffEmptyIndexOnOkIsIgnored(testCase)
            %Site 2 through its other caller.
            localSetProject(testCase, localProject(2));
            localAnswer([], 1);
            err = localFire(testCase, 'Export Handoff Table');
            testCase.verifyEmpty(err, sprintf('Export Handoff Table threw on an empty OK: %s', err));
            testCase.verifyEqual(localCalls('listdlg'), 1);
            testCase.verifyFalse(localDialogExists('Handoff Error'), ...
                'nothing may be exported when no spec was chosen');
        end

        function testSingleSpecSkipsThePicker(testCase)
            %Control for the picker weld: one spec never opens listdlg.
            localSetProject(testCase, localProject(1));
            localAnswer([], 1);
            err = localFire(testCase, 'Run Scoring');
            testCase.verifyEmpty(err, sprintf('Run Scoring threw with one spec: %s', err));
            testCase.verifyEqual(localCalls('listdlg'), 0, 'a single spec is used without a picker');
            testCase.verifyTrue(localDialogExists('Scoring Error'), ...
                'the single spec must reach the scoring call');
        end

    end

    methods (TestMethodTeardown)
        function teardownShadowsAndScreen(testCase)
            localCloseByName('Scoring Error');
            localCloseByName('Handoff Error');
            localCloseByName('Spec Error');
            localCloseByTag('psyrat_import_gui');
            clear global PSYRAT_TEST_LISTDLG_IDX PSYRAT_TEST_LISTDLG_OK ...
                PSYRAT_TEST_LISTDLG_CALLS PSYRAT_TEST_INPUTDLG_CALLS
            if any(strcmp(strsplit(path, pathsep), testCase.ShadowDir))
                rmpath(testCase.ShadowDir);
            end
            rehash;
            if isfolder(testCase.ShadowDir); rmdir(testCase.ShadowDir, 's'); end
        end
    end
end


function project = localProject(nspecs)
%A project holding nspecs scoring specs and nothing else. It is enough for
%the pickers; psyrat_project_score then fails on it inside the callers' try
%(no canonical field), which is the observable of the positive controls.
specs = struct('spec_id', {}, 'name', {}, 'method', {});
for k = 1:nspecs
    specs(k).spec_id = sprintf('spec_%d', k);
    specs(k).name = sprintf('Spec %d', k);
    specs(k).method = 'mean_amplitude';
end
project = struct();
project.scoring_specs = specs;
project.ui_state = struct('active_spec_id', specs(1).spec_id);
end


function localSetProject(testCase, project)
st = guidata(testCase.Fig);
st.project = project;
guidata(testCase.Fig, st);
end


function localAnswer(idx, ok)
global PSYRAT_TEST_LISTDLG_IDX PSYRAT_TEST_LISTDLG_OK %#ok<GVMIS>
PSYRAT_TEST_LISTDLG_IDX = idx;
PSYRAT_TEST_LISTDLG_OK = ok;
end


function localResetGlobals()
global PSYRAT_TEST_LISTDLG_CALLS PSYRAT_TEST_INPUTDLG_CALLS %#ok<GVMIS>
PSYRAT_TEST_LISTDLG_CALLS = 0;
PSYRAT_TEST_INPUTDLG_CALLS = 0;
end


function n = localCalls(which)
global PSYRAT_TEST_LISTDLG_CALLS PSYRAT_TEST_INPUTDLG_CALLS %#ok<GVMIS>
switch which
    case 'listdlg'
        n = PSYRAT_TEST_LISTDLG_CALLS;
    case 'inputdlg'
        n = PSYRAT_TEST_INPUTDLG_CALLS;
end
if isempty(n); n = 0; end
end


function err = localFire(testCase, label)
%Fire a button's own Callback and return the message of any escaping error
%('' when it returned). The buttons hold nested-function handles.
btn = findall(testCase.Fig, 'Style', 'pushbutton', 'String', label);
testCase.assertNumElements(btn, 1, sprintf('button "%s" not found', label));
err = '';
try
    cb = btn.Callback;
    if iscell(cb)
        cb{1}(btn, [], cb{2:end});
    else
        cb(btn, []);
    end
catch ME
    err = ME.message;
end
end


function tf = localDialogExists(name)
tf = ~isempty(findall(groot, 'Type', 'figure', 'Name', name));
end


function localCloseByName(name)
h = findall(groot, 'Type', 'figure', 'Name', name);
if ~isempty(h); delete(h); end
end


function localCloseByTag(tag)
h = findall(groot, 'Type', 'figure', 'Tag', tag);
if ~isempty(h); delete(h); end
end
