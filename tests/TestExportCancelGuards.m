classdef TestExportCancelGuards < PsyRATTestBase
    %Sweep follow-up F3 (ruled 2026-09-16): eight table Save callbacks ran
    %fileparts(fullfile(savepath,savename)) straight after uiputfile with no
    %cancel guard, and dispatched on the extension with no else. On Cancel
    %uiputfile returns (0,0); fullfile(0,0) is a three-byte char whose
    %extension is '' (probed on R2025b), so the callback fell through both
    %branches: no file, no message, a button that looked dead. A typed name
    %ending in neither .xlsx nor .csv took the same silent route. Seven sibling
    %Save sites already guarded (psyrat_criterionfigures.m, psyrat_savesscoeffs
    %and psyrat_savealldifftable in psyrat_variancet.m, the DoD tables, the
    %viewers); psyrat_startproc's save dialog is the messaged precedent.
    %
    %ROUTE. Each builder is run with gui = 1 and its Save button's own
    %Callback is fired through a uiputfile shadow that answers from a global:
    %0 is Cancel, a path is that name and folder. The eighth site,
    %psyrat_saveallvartable, lives on the window the trt variance table's
    %"View All Standard Deviations" button opens, so that button is fired
    %first. A modal errordlg does not block, so a dialog is observed by name.
    %
    %OBSERVABLES. Cancel: no error escapes, no dialog, no figure created (a
    %regression pin, silent both before and after because of the fullfile
    %behavior above). Unsupported extension: the 'Unsupported File Type'
    %dialog and no file (RED before the fix). The guard itself is pinned at
    %source: every site's guard sits between its uiputfile and its fileparts.

    properties (Access = private)
        ShadowDir
        TmpDir
    end

    methods (TestMethodSetup)
        function setupShadowAndTmp(testCase)
            testCase.ShadowDir = tempname;
            mkdir(testCase.ShadowDir);
            fid = fopen(fullfile(testCase.ShadowDir, 'uiputfile.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function [filename,pathname,filterindex] = uiputfile(varargin)', ...
                '% TEST SHADOW: answers the save dialog from a global; 0 is Cancel.', ...
                'global PSYRAT_TEST_UIPUTFILE', ...
                'if isempty(PSYRAT_TEST_UIPUTFILE) || isequal(PSYRAT_TEST_UIPUTFILE, 0)', ...
                '    filename = 0; pathname = 0; filterindex = 0; return;', ...
                'end', ...
                '[p,n,e] = fileparts(PSYRAT_TEST_UIPUTFILE);', ...
                'filename = [n e]; pathname = [p filesep]; filterindex = 1;', ...
                'end');
            fclose(fid);
            addpath(testCase.ShadowDir, '-begin');
            rehash;
            testCase.assertTrue(startsWith(which('uiputfile'), testCase.ShadowDir), ...
                'the uiputfile shadow must be the resolved uiputfile');
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);
            localCloseOutputs();
        end
    end

    methods (Test)

        function testCancelIsSilentAtEverySite(testCase)
            %Regression pin: Cancel writes nothing, says nothing, throws
            %nothing at all eight sites.
            global PSYRAT_TEST_UIPUTFILE %#ok<GVMIS>
            sites = localSites();
            for s = 1:size(sites, 1)
                PSYRAT_TEST_UIPUTFILE = 0;
                fig = localOpenSite(testCase, sites(s, :));
                nbefore = numel(findall(groot, 'Type', 'figure'));
                err = localFireSave(testCase, fig, sites{s, 3});
                testCase.verifyEmpty(err, sprintf('%s threw on Cancel: %s', sites{s, 1}, err));
                testCase.verifyFalse(localDialogExists('Unsupported File Type'), ...
                    sprintf('%s: Cancel must not be reported as an unsupported type', sites{s, 1}));
                testCase.verifyEqual(numel(findall(groot, 'Type', 'figure')), nbefore, ...
                    sprintf('%s: Cancel must create no window', sites{s, 1}));
                localCloseOutputs();
            end
        end

        function testUnsupportedExtensionRaisesTheDialogAtEverySite(testCase)
            %A name the filter did not supply (x.txt) must be refused with the
            %modal 'Unsupported File Type' and write nothing; before the fix it
            %fell through the xlsx/csv chain silently.
            global PSYRAT_TEST_UIPUTFILE %#ok<GVMIS>
            sites = localSites();
            for s = 1:size(sites, 1)
                target = fullfile(testCase.TmpDir, sprintf('site%d.txt', s));
                PSYRAT_TEST_UIPUTFILE = target;
                fig = localOpenSite(testCase, sites(s, :));
                err = localFireSave(testCase, fig, sites{s, 3});
                testCase.verifyEmpty(err, sprintf('%s threw on a .txt name: %s', sites{s, 1}, err));
                testCase.verifyTrue(localDialogExists('Unsupported File Type'), ...
                    sprintf('%s: a .txt name must raise the Unsupported File Type dialog', sites{s, 1}));
                testCase.verifyFalse(isfile(target), ...
                    sprintf('%s: nothing may be written for a .txt name', sites{s, 1}));
                d = findall(groot, 'Type', 'figure', 'Name', 'Unsupported File Type');
                if ~isempty(d)
                    testCase.verifyEqual(get(d(1), 'WindowStyle'), 'modal', ...
                        sprintf('%s: the refusal must be modal', sites{s, 1}));
                    delete(d);
                end
                localCloseOutputs();
            end
        end

        function testARealCsvNameStillSaves(testCase)
            %Positive control for the shadow and the guard: a real .csv name
            %writes the file at a table site and at the all-SDs site.
            global PSYRAT_TEST_UIPUTFILE %#ok<GVMIS>
            sites = localSites();
            for s = [3 5]
                target = fullfile(testCase.TmpDir, sprintf('site%d.csv', s));
                PSYRAT_TEST_UIPUTFILE = target;
                fig = localOpenSite(testCase, sites(s, :));
                err = localFireSave(testCase, fig, sites{s, 3});
                testCase.verifyEmpty(err, sprintf('%s threw on a .csv name: %s', sites{s, 1}, err));
                testCase.verifyTrue(isfile(target), sprintf('%s: the .csv must be written', sites{s, 1}));
                testCase.verifyFalse(localDialogExists('Unsupported File Type'));
                localCloseOutputs();
            end
        end

        function testEveryGuardSitsBetweenItsDialogAndItsFileparts(testCase)
            %Source pin: in each of the five table files every
            %fileparts(fullfile(savepath,savename)) is preceded, within its own
            %callback, by the cancel guard, and every xlsx/csv chain has the
            %refusal else. Counts are exact so a lost site cannot pass.
            root = testCase.projectRoot();
            %{file, fileparts(fullfile(...)) sites, refusal elses}; the variance
            %file's two already-guarded siblings use fileparts(savename), so
            %they are not among its two sites here
            files = {'psyrat_depoverallt.m', 2, 2; 'psyrat_depcutofft.m', 1, 1; ...
                'psyrat_variancet.m', 2, 2; 'psyrat_trt_relcutofft.m', 1, 1; ...
                'psyrat_trt_reloverallt.m', 2, 2; ...
                'psyrat_dod_observedt.m', 0, 1; 'psyrat_dod_varcompt.m', 0, 1};
            for k = 1:size(files, 1)
                src = fileread(fullfile(root, 'subroutines', 'tables', files{k, 1}));
                fp = strfind(src, '[~,~,ext] = fileparts(fullfile(savepath,savename));');
                testCase.assertEqual(numel(fp), files{k, 2}, ...
                    sprintf('%s: the fileparts sites moved', files{k, 1}));
                for i = 1:numel(fp)
                    before = src(max(1, fp(i) - 900):fp(i));
                    iu = strfind(before, 'uiputfile(');
                    testCase.assertNotEmpty(iu, sprintf('%s site %d: no uiputfile before fileparts', files{k, 1}, i));
                    between = before(iu(end):end);
                    testCase.verifyEqual(numel(strfind(between, 'if isequal(savename,0) || isequal(savepath,0)')), 1, ...
                        sprintf('%s site %d: the cancel guard must sit between uiputfile and fileparts', files{k, 1}, i));
                end
                testCase.verifyEqual(numel(strfind(src, '''Unsupported File Type'',''modal'')')), files{k, 3}, ...
                    sprintf('%s: every xlsx/csv chain must end in the modal refusal', files{k, 1}));
            end
        end

    end

    methods (TestMethodTeardown)
        function teardownShadowAndTmp(testCase)
            localCloseOutputs();
            delete(findall(groot, 'Type', 'figure', 'Name', 'Unsupported File Type'));
            clear global PSYRAT_TEST_UIPUTFILE
            if any(strcmp(strsplit(path, pathsep), testCase.ShadowDir))
                rmpath(testCase.ShadowDir);
            end
            rehash;
            if isfolder(testCase.ShadowDir); rmdir(testCase.ShadowDir, 's'); end
            if isfolder(testCase.TmpDir); rmdir(testCase.TmpDir, 's'); end
        end
    end
end


function sites = localSites()
%{name, builder thunk, Save-button label fragment, opener} for the eight
%unguarded sites plus the two DoD table windows, which already guarded Cancel
%but shared the else-less extension chain; the opener (empty or 'View All')
%fires a button on the builder's window first and the Save is then fired on
%the newest window.
sites = { ...
    'depoverallt Save Table',  @() psyrat_depoverallt('psyrat_data', localIcData(), 'gui', 1),  'Save Table', ''; ...
    'depoverallt Save IDs',    @() psyrat_depoverallt('psyrat_data', localIcData(), 'gui', 1),  'Save IDs',   ''; ...
    'depcutofft Save Table',   @() psyrat_depcutofft('psyrat_data', localIcData(), 'gui', 1),   'Save Table', ''; ...
    'variancet Save Table',    @() psyrat_variancet('psyrat_data', localIcData(), 'gui', 1),    'Save Table', ''; ...
    'variancet all SDs Save Table', @() psyrat_variancet('psyrat_data', localTrtData(), 'gui', 1), 'Save Table', 'View All'; ...
    'trt_relcutofft Save Table', @() psyrat_trt_relcutofft('psyrat_data', localTrtData(), 'gui', 1), 'Save Table', ''; ...
    'trt_reloverallt Save Table', @() psyrat_trt_reloverallt('psyrat_data', localTrtData(), 'gui', 1), 'Save Table', ''; ...
    'trt_reloverallt Save IDs', @() psyrat_trt_reloverallt('psyrat_data', localTrtData(), 'gui', 1), 'Save IDs', ''; ...
    'dod_varcompt Save Table',  @() psyrat_dod_varcompt('psyrat_data', localDodData(), 'groupidx', 1, 'map', localDodMap(localDodData().rel), 'gui', 1), 'Save Table', ''; ...
    'dod_observedt Save Table', @() psyrat_dod_observedt('psyrat_data', localDodData(), 'groupidx', 1, 'map', localDodMap(localDodData().rel), 'gui', 1), 'Save Table', ''; ...
    };
end


function fig = localOpenSite(testCase, site)
localCloseOutputs();
site{2}();
figs = findobj(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
testCase.assertNotEmpty(figs, sprintf('%s: the builder created no output window', site{1}));
fig = figs(1);
if ~isempty(site{4})
    localFireButton(testCase, fig, site{4});
    figs2 = findobj(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
    figs2 = figs2(~ismember(figs2, figs));
    testCase.assertNumElements(figs2, 1, sprintf('%s: the opener must open one window', site{1}));
    fig = figs2(1);
end
end


function err = localFireSave(testCase, fig, labelfrag)
err = '';
try
    localFireButton(testCase, fig, labelfrag);
catch ME
    err = ME.message;
end
end


function localFireButton(testCase, fig, labelfrag)
btns = findall(fig, 'Style', 'pushbutton');
for b = btns(:)'
    s = b.String; if iscell(s), s = strjoin(s, ' '); end
    if contains(string(s), labelfrag)
        cb = b.Callback;
        if iscell(cb)
            feval(cb{1}, b, [], cb{2:end});
        else
            cb(b, []);
        end
        return;
    end
end
testCase.assertFail(sprintf('no "%s" button on the window', labelfrag));
end


function tf = localDialogExists(name)
tf = ~isempty(findall(groot, 'Type', 'figure', 'Name', name));
end


function localCloseOutputs()
delete(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output'));
end


function pd = localDodData()
pd = localAddSaveHeaderFields(PsyRATTestDataFactory.makeDoDRelData());
end


function map = localDodMap(rel)
events = cellstr(string(rel.events(:)));
mapnames = cellstr(string(rel.dod_map(:)));
[~, map] = ismember(mapnames, events);
map = map(:)';
end


function pd = localTrtData()
pd = localAddSaveHeaderFields(PsyRATTestDataFactory.makeTRTSummaryData());
end


function pd = localIcData()
pd = localAddSaveHeaderFields(PsyRATTestDataFactory.makeICSummaryDataWithTrial());
end


function pd = localAddSaveHeaderFields(pd)
%the same provenance fields TestExportCallbacks adds: the writers' headers
%read them, and a real .csv save (the positive control) reaches the writers
pd.ver = 'test';
pd.rel.filename = 'export_cancel_fixture.psyrat';
pd.rel.nchains = 2;
pd.rel.niter = 600;
pd.rel.seed = 12345;
pd.rel.engine = 'cmdstan';
pd.rel.priors = psyrat_default_priors;
pd.rel.out.conv.converged = 1;
pd.rel.out.conv.ndivergent = 0;
end
